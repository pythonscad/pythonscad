#ifdef ENABLE_OCCT

#include "geometry/brep/ExactShape.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

#include <BRepAlgoAPI_BooleanOperation.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_GTransform.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <Standard_Failure.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Ax2.hxx>
#include <gp_GTrsf.hxx>
#include <gp_Mat.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <gp_XYZ.hxx>

// printutils.h cannot be included here: its struct Message clashes with OCCT's
// class Message (pulled in by the boolean headers). Logging goes through BRepGeometry.cc.
namespace BRepUtils {
void logWarning(const std::string& msg);
}

using ShapePtr = std::shared_ptr<const TopoDS_Shape>;

struct ExactShape::Impl {
  std::string what;                 // for diagnostics
  std::function<ShapePtr()> build;  // the recipe, released after evaluation
  std::once_flag once;
  ShapePtr result;
};

// Small factory so the helpers below can create nodes through the private constructor.
struct ExactShapeFactory {
  static std::shared_ptr<const ExactShape> make(std::string what, std::function<ShapePtr()> build)
  {
    auto impl = std::make_unique<ExactShape::Impl>();
    impl->what = std::move(what);
    impl->build = std::move(build);
    return std::shared_ptr<const ExactShape>(new ExactShape(std::move(impl)));
  }
};

namespace {

ShapePtr guarded(const std::string& what, const std::function<TopoDS_Shape()>& f)
{
  try {
    TopoDS_Shape s = f();
    if (s.IsNull()) {
      BRepUtils::logWarning("Exact geometry: " + what + " produced no shape");
      return nullptr;
    }
    return std::make_shared<const TopoDS_Shape>(std::move(s));
  } catch (const Standard_Failure& e) {
    BRepUtils::logWarning("Exact geometry: " + what + " failed (" +
                          (e.GetMessageString() ? e.GetMessageString() : "unknown error") + ")");
  } catch (const std::exception& e) {
    BRepUtils::logWarning("Exact geometry: " + what + " failed (" + e.what() + ")");
  }
  return nullptr;
}

/// Apply an OpenSCAD affine transform. Similarities (rotation, translation, uniform
/// scale, mirror) keep analytic surfaces; anything else (non-uniform scale, shear)
/// converts them to B-splines, which is still exact.
TopoDS_Shape applyTransform(const TopoDS_Shape& shape, const Transform3d& mat)
{
  const Eigen::Matrix3d lin = mat.linear();
  const Eigen::Vector3d tr = mat.translation();
  const Eigen::Matrix3d ltl = lin.transpose() * lin;
  const double s2 = ltl.trace() / 3.0;
  const bool similarity =
    s2 > 0 && (ltl - s2 * Eigen::Matrix3d::Identity()).norm() <= 1e-9 * std::max(1.0, s2);

  if (similarity) {
    gp_Trsf trsf;
    trsf.SetValues(lin(0, 0), lin(0, 1), lin(0, 2), tr[0],  //
                   lin(1, 0), lin(1, 1), lin(1, 2), tr[1],  //
                   lin(2, 0), lin(2, 1), lin(2, 2), tr[2]);
    BRepBuilderAPI_Transform op(shape, trsf, /*copy*/ true);
    return op.IsDone() ? op.Shape() : TopoDS_Shape();
  }
  if (std::abs(lin.determinant()) < 1e-12) return {};
  gp_GTrsf gtrsf;
  gtrsf.SetVectorialPart(gp_Mat(lin(0, 0), lin(0, 1), lin(0, 2),  //
                                lin(1, 0), lin(1, 1), lin(1, 2),  //
                                lin(2, 0), lin(2, 1), lin(2, 2)));
  gtrsf.SetTranslationPart(gp_XYZ(tr[0], tr[1], tr[2]));
  BRepBuilderAPI_GTransform op(shape, gtrsf, /*copy*/ true);
  return op.IsDone() ? op.Shape() : TopoDS_Shape();
}

/// Run one OCCT boolean. Non-destructive, because the inputs are shared with
/// other geometries and caches. Returns a null shape if OCCT fails or the
/// result is not a valid solid, so the caller falls back to the mesh.
TopoDS_Shape runBoolean(BRepAlgoAPI_BooleanOperation&& op, const TopoDS_Shape& argument,
                        const std::vector<ShapePtr>& tools)
{
  TopTools_ListOfShape args, tls;
  args.Append(argument);
  for (const auto& t : tools) tls.Append(*t);
  op.SetArguments(args);
  op.SetTools(tls);
  op.SetNonDestructive(Standard_True);
  op.Build();
  if (!op.IsDone() || op.HasErrors()) return {};
  const TopoDS_Shape raw = op.Shape();
  if (!BRepCheck_Analyzer(raw).IsValid()) return {};

  // Merge faces that lie on the same surface (e.g. two touching cubes), so the
  // STEP output is as clean as if the solid had been modelled in one piece.
  // Edges are deliberately not unified: in OCCT 7.6 that can corrupt solids
  // (e.g. a sphere poking through all faces of a cube). Keep the raw result if
  // the cleanup does not yield a valid solid.
  ShapeUpgrade_UnifySameDomain unify(raw, /*UnifyEdges*/ false, /*UnifyFaces*/ true,
                                     /*ConcatBSplines*/ false);
  unify.Build();
  const TopoDS_Shape clean = unify.Shape();
  return (!clean.IsNull() && BRepCheck_Analyzer(clean).IsValid()) ? clean : raw;
}

std::shared_ptr<const ExactShape> makeNode(std::string what, std::function<ShapePtr()> build)
{
  return ExactShapeFactory::make(std::move(what), std::move(build));
}

}  // namespace

ExactShape::ExactShape(std::unique_ptr<Impl> impl) : impl_(std::move(impl))
{
}
ExactShape::~ExactShape() = default;

const TopoDS_Shape *ExactShape::shape() const
{
  std::call_once(impl_->once, [this] {
    if (impl_->build) impl_->result = impl_->build();
    impl_->build = nullptr;  // drop references to children, they are not needed anymore
  });
  return impl_->result.get();
}

std::shared_ptr<const ExactShape> ExactShape::sphere(double r)
{
  return makeNode("sphere",
                  [r] { return guarded("sphere", [r] { return BRepPrimAPI_MakeSphere(r).Shape(); }); });
}

std::shared_ptr<const ExactShape> ExactShape::box(const Vector3d& p1, const Vector3d& p2)
{
  return makeNode("cube", [p1, p2] {
    return guarded("cube", [&] {
      return BRepPrimAPI_MakeBox(gp_Pnt(p1[0], p1[1], p1[2]), gp_Pnt(p2[0], p2[1], p2[2])).Shape();
    });
  });
}

std::shared_ptr<const ExactShape> ExactShape::cylinder(double r1, double r2, double z1, double z2,
                                                       double angle)
{
  return makeNode("cylinder", [=] {
    return guarded("cylinder", [&] {
      const gp_Ax2 axes(gp_Pnt(0, 0, z1), gp::DZ(), gp::DX());
      const double h = z2 - z1;
      const double rad = std::min(angle, 360.0) * M_PI / 180.0;
      if (r1 == r2) return BRepPrimAPI_MakeCylinder(axes, r1, h, rad).Shape();
      return BRepPrimAPI_MakeCone(axes, r1, r2, h, rad).Shape();
    });
  });
}

std::shared_ptr<const ExactShape> ExactShape::transformed(std::shared_ptr<const ExactShape> child,
                                                          const Transform3d& mat)
{
  if (!child) return nullptr;
  if (mat.matrix() == Eigen::Matrix4d::Identity()) return child;
  return makeNode("transform", [child = std::move(child), mat]() -> ShapePtr {
    const TopoDS_Shape *s = child->shape();
    if (!s) return nullptr;
    return guarded("transform", [&] { return applyTransform(*s, mat); });
  });
}

std::shared_ptr<const ExactShape> ExactShape::boolean(
  OpenSCADOperator op, std::vector<std::shared_ptr<const ExactShape>> children)
{
  if (children.empty()) return nullptr;
  if (std::any_of(children.begin(), children.end(), [](const auto& c) { return !c; })) return nullptr;
  if (children.size() == 1) return children.front();

  const char *what = nullptr;
  switch (op) {
  case OpenSCADOperator::UNION:        what = "union"; break;
  case OpenSCADOperator::DIFFERENCE:   what = "difference"; break;
  case OpenSCADOperator::INTERSECTION: what = "intersection"; break;
  default:                             return nullptr;
  }

  return makeNode(what, [op, what, children = std::move(children)]() -> ShapePtr {
    std::vector<ShapePtr> shapes;
    for (const auto& c : children) {
      const TopoDS_Shape *s = c->shape();
      if (!s) return nullptr;
      // non-owning alias; the children are kept alive by this lambda while it runs
      shapes.emplace_back(ShapePtr(c, s));
    }
    return guarded(what, [&]() -> TopoDS_Shape {
      const TopoDS_Shape& first = *shapes.front();
      const std::vector<ShapePtr> rest(shapes.begin() + 1, shapes.end());
      switch (op) {
      case OpenSCADOperator::UNION:      return runBoolean(BRepAlgoAPI_Fuse(), first, rest);
      case OpenSCADOperator::DIFFERENCE: return runBoolean(BRepAlgoAPI_Cut(), first, rest);
      default:                           {
        // OCCT's Common treats several tools as one union, so intersect pairwise.
        TopoDS_Shape acc = first;
        for (const auto& t : rest) {
          acc = runBoolean(BRepAlgoAPI_Common(), acc, {t});
          if (acc.IsNull()) break;
        }
        return acc;
      }
      }
    });
  });
}

#endif  // ENABLE_OCCT
