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
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <GC_MakeCircle.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <ShapeFix_Face.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Circ.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <Bnd_Box.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Poly_Triangulation.hxx>
#include <STEPControl_Reader.hxx>
#include <ShapeUpgrade_ShapeDivideClosed.hxx>
#include <TopAbs.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
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

/// Make a solid (or compound of solids) point outwards: some sweeps come out
/// inside-out depending on the direction of the profile normal and the sweep.
TopoDS_Shape orientOutwards(TopoDS_Shape shape)
{
  GProp_GProps props;
  BRepGProp::VolumeProperties(shape, props);
  if (props.Mass() < 0) shape.Reverse();
  return shape;
}

/// Valid solid with positive volume, or a null shape
TopoDS_Shape checkedSolid(const TopoDS_Shape& shape)
{
  if (shape.IsNull()) return {};
  TopoDS_Shape oriented = orientOutwards(shape);
  if (!BRepCheck_Analyzer(oriented).IsValid()) return {};
  GProp_GProps props;
  BRepGProp::VolumeProperties(oriented, props);
  return props.Mass() > 0 ? oriented : TopoDS_Shape();
}

double signedArea(const std::vector<Vector2d>& o)
{
  double a = 0;
  for (size_t i = 0, n = o.size(); i < n; ++i) {
    const auto& p = o[i];
    const auto& q = o[(i + 1) % n];
    a += p[0] * q[1] - q[0] * p[1];
  }
  return a / 2;
}

bool pointInOutline(const Vector2d& pt, const std::vector<Vector2d>& o)
{
  bool inside = false;
  for (size_t i = 0, j = o.size() - 1; i < o.size(); j = i++) {
    if ((o[i][1] > pt[1]) != (o[j][1] > pt[1]) &&
        pt[0] < (o[j][0] - o[i][0]) * (pt[1] - o[i][1]) / (o[j][1] - o[i][1]) + o[i][0]) {
      inside = !inside;
    }
  }
  return inside;
}

TopoDS_Wire polygonWire(const std::vector<Vector2d>& outline)
{
  BRepBuilderAPI_MakePolygon poly;
  for (const auto& v : outline) poly.Add(gp_Pnt(v[0], v[1], 0));
  poly.Close();
  if (!poly.IsDone()) return {};
  return poly.Wire();
}

/// All faces of shape, or the shape itself if it is a face
std::vector<TopoDS_Face> facesOf(const TopoDS_Shape& shape)
{
  std::vector<TopoDS_Face> faces;
  for (TopExp_Explorer ex(shape, TopAbs_FACE); ex.More(); ex.Next())
    faces.push_back(TopoDS::Face(ex.Current()));
  return faces;
}

TopoDS_Shape compoundOf(const std::vector<TopoDS_Shape>& shapes)
{
  if (shapes.size() == 1) return shapes.front();
  TopoDS_Compound compound;
  BRep_Builder builder;
  builder.MakeCompound(compound);
  for (const auto& s : shapes) builder.Add(compound, s);
  return compound;
}

gp_Trsf toTrsf(const Transform3d& mat)
{
  const Eigen::Matrix3d lin = mat.linear();
  const Eigen::Vector3d tr = mat.translation();
  gp_Trsf trsf;
  trsf.SetValues(lin(0, 0), lin(0, 1), lin(0, 2), tr[0], lin(1, 0), lin(1, 1), lin(1, 2), tr[1],
                 lin(2, 0), lin(2, 1), lin(2, 2), tr[2]);
  return trsf;
}

/// Ruled loft from wire to the same wire scaled by s about origin and moved by height
TopoDS_Shape loftWire(const TopoDS_Wire& wire, const gp_Pnt& origin, double s, const gp_Vec& height)
{
  BRepOffsetAPI_ThruSections loft(/*isSolid*/ Standard_True, /*ruled*/ Standard_True);
  loft.AddWire(wire);
  if (s == 0) {
    loft.AddVertex(BRepBuilderAPI_MakeVertex(origin.Translated(height)).Vertex());
  } else {
    gp_Trsf scaling;
    scaling.SetScale(origin, s);
    gp_Trsf move;
    move.SetTranslation(height);
    loft.AddWire(TopoDS::Wire(BRepBuilderAPI_Transform(wire, move * scaling, true).Shape()));
  }
  loft.Build();
  return loft.IsDone() ? loft.Shape() : TopoDS_Shape();
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

std::shared_ptr<const ExactShape> ExactShape::circle2d(double r, double angle)
{
  return makeNode("circle", [r, angle] {
    return guarded("circle", [&]() -> TopoDS_Shape {
      const gp_Circ circ(gp_Ax2(gp::Origin(), gp::DZ(), gp::DX()), r);
      if (angle >= 360) {
        return BRepBuilderAPI_MakeFace(BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(circ)).Wire())
          .Face();
      }
      const double a = angle * M_PI / 180;
      const gp_Pnt c = gp::Origin(), p0(r, 0, 0), p1(r * std::cos(a), r * std::sin(a), 0);
      BRepBuilderAPI_MakeWire wire;
      wire.Add(BRepBuilderAPI_MakeEdge(c, p0).Edge());
      wire.Add(BRepBuilderAPI_MakeEdge(GC_MakeArcOfCircle(circ, 0, a, true).Value()).Edge());
      wire.Add(BRepBuilderAPI_MakeEdge(p1, c).Edge());
      return BRepBuilderAPI_MakeFace(wire.Wire(), /*OnlyPlane*/ true).Face();
    });
  });
}

std::shared_ptr<const ExactShape> ExactShape::rect2d(const Vector2d& p1, const Vector2d& p2)
{
  return polygons2d({{p1, {p2[0], p1[1]}, p2, {p1[0], p2[1]}}});
}

std::shared_ptr<const ExactShape> ExactShape::polygons2d(std::vector<std::vector<Vector2d>> outlines)
{
  return makeNode("polygon", [outlines = std::move(outlines)] {
    return guarded("polygon", [&]() -> TopoDS_Shape {
      // Material outlines (counter-clockwise) each become a face; every hole goes
      // to the smallest material outline that contains it.
      std::vector<size_t> material, holes;
      std::vector<double> area(outlines.size());
      for (size_t i = 0; i < outlines.size(); ++i) {
        if (outlines[i].size() < 3) continue;
        area[i] = signedArea(outlines[i]);
        if (area[i] > 0) material.push_back(i);
        else if (area[i] < 0) holes.push_back(i);
      }
      std::vector<std::vector<size_t>> holesOf(outlines.size());
      for (size_t h : holes) {
        size_t best = outlines.size();
        for (size_t m : material) {
          if (area[m] <= -area[h]) continue;
          if (!pointInOutline(outlines[h][0], outlines[m])) continue;
          if (best == outlines.size() || area[m] < area[best]) best = m;
        }
        if (best != outlines.size()) holesOf[best].push_back(h);
      }
      std::vector<TopoDS_Shape> faces;
      for (size_t m : material) {
        const TopoDS_Wire outer = polygonWire(outlines[m]);
        if (outer.IsNull()) return {};
        BRepBuilderAPI_MakeFace face(outer, /*OnlyPlane*/ true);
        for (size_t h : holesOf[m]) {
          const TopoDS_Wire inner = polygonWire(outlines[h]);
          if (inner.IsNull()) return {};
          face.Add(inner);
        }
        if (!face.IsDone()) return {};
        ShapeFix_Face fix(face.Face());
        fix.FixOrientation();
        fix.Perform();
        faces.push_back(fix.Face());
      }
      if (faces.empty()) return {};
      return compoundOf(faces);
    });
  });
}

std::shared_ptr<const ExactShape> ExactShape::prism(std::shared_ptr<const ExactShape> profile,
                                                    const Transform3d& placement, const Vector3d& base,
                                                    const Vector3d& height, double scale)
{
  if (!profile) return nullptr;
  return makeNode(
    "linear_extrude", [profile = std::move(profile), placement, base, height, scale]() -> ShapePtr {
      const TopoDS_Shape *p = profile->shape();
      if (!p) return nullptr;
      return guarded("linear_extrude", [&]() -> TopoDS_Shape {
        gp_Trsf move;
        move.SetTranslation(gp_Vec(base[0], base[1], base[2]));
        const gp_Trsf trsf = move * toTrsf(placement);
        const gp_Vec h(height[0], height[1], height[2]);
        std::vector<TopoDS_Shape> solids;
        for (const TopoDS_Face& face : facesOf(*p)) {
          const TopoDS_Face placed = TopoDS::Face(BRepBuilderAPI_Transform(face, trsf, true).Shape());
          if (scale == 1) {
            BRepPrimAPI_MakePrism prism(placed, h);
            if (!prism.IsDone()) return {};
            solids.push_back(prism.Shape());
            continue;
          }
          // Scaled: loft the outer boundary and cut lofts of the holes out of it.
          const gp_Pnt origin = gp::Origin().Transformed(trsf);
          const TopoDS_Wire outer = BRepTools::OuterWire(placed);
          TopoDS_Shape solid = checkedSolid(loftWire(outer, origin, scale, h));
          if (solid.IsNull()) return {};
          for (TopExp_Explorer w(placed, TopAbs_WIRE); w.More(); w.Next()) {
            const TopoDS_Wire wire = TopoDS::Wire(w.Current());
            if (wire.IsSame(outer)) continue;
            const TopoDS_Shape hole = checkedSolid(loftWire(wire, origin, scale, h));
            if (hole.IsNull()) return {};
            solid = runBoolean(BRepAlgoAPI_Cut(), solid, {std::make_shared<const TopoDS_Shape>(hole)});
            if (solid.IsNull()) return {};
          }
          solids.push_back(solid);
        }
        if (solids.empty()) return {};
        return checkedSolid(compoundOf(solids));
      });
    });
}

std::shared_ptr<const ExactShape> ExactShape::revolve(std::shared_ptr<const ExactShape> profile,
                                                      double start, double angle)
{
  if (!profile) return nullptr;
  return makeNode("rotate_extrude", [profile = std::move(profile), start, angle]() -> ShapePtr {
    const TopoDS_Shape *p = profile->shape();
    if (!p) return nullptr;
    return guarded("rotate_extrude", [&]() -> TopoDS_Shape {
      // profile x stays the radius, profile y becomes z
      gp_Trsf upright;
      upright.SetRotation(gp::OX(), M_PI / 2);
      const TopoDS_Shape placed = BRepBuilderAPI_Transform(*p, upright, true).Shape();
      const double a = std::max(-360.0, std::min(360.0, angle)) * M_PI / 180;
      BRepPrimAPI_MakeRevol revol(placed, gp::OZ(), a);
      if (!revol.IsDone()) return {};
      gp_Trsf turn;
      turn.SetRotation(gp::OZ(), start * M_PI / 180);
      return checkedSolid(BRepBuilderAPI_Transform(revol.Shape(), turn, true).Shape());
    });
  });
}

std::shared_ptr<const ExactShape> ExactShape::fromShape(const TopoDS_Shape& shape)
{
  auto stored = std::make_shared<const TopoDS_Shape>(shape);
  return makeNode("import", [stored] { return stored; });
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

TopoDS_Shape ExactShapeUtils::splitClosedFaces(const TopoDS_Shape& shape)
{
  try {
    ShapeUpgrade_ShapeDivideClosed divide(shape);
    divide.SetNbSplitPoints(5);
    divide.Perform();
    return divide.Result();
  } catch (const Standard_Failure&) {
    return shape;
  }
}

bool ExactShapeUtils::readStep(const std::string& filename, TopoDS_Shape& shape, std::string& error)
{
  try {
    STEPControl_Reader reader;
    if (reader.ReadFile(filename.c_str()) != IFSelect_RetDone) {
      error = "cannot read file";
      return false;
    }
    if (reader.TransferRoots() == 0) {
      error = "no transferable geometry";
      return false;
    }
    shape = reader.OneShape();
    if (shape.IsNull()) {
      error = "empty shape";
      return false;
    }
    return true;
  } catch (const Standard_Failure& e) {
    error = e.GetMessageString() ? e.GetMessageString() : "unknown error";
  }
  return false;
}

double ExactShapeUtils::boundingDiagonal(const TopoDS_Shape& shape)
{
  Bnd_Box box;
  BRepBndLib::Add(shape, box);
  if (box.IsVoid()) return 0;
  double x0, y0, z0, x1, y1, z1;
  box.Get(x0, y0, z0, x1, y1, z1);
  return std::sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) + (z1 - z0) * (z1 - z0));
}

void ExactShapeUtils::solidInfo(const TopoDS_Shape& shape, int& solids, double& volume)
{
  solids = 0;
  for (TopExp_Explorer ex(shape, TopAbs_SOLID); ex.More(); ex.Next()) ++solids;
  GProp_GProps props;
  BRepGProp::VolumeProperties(shape, props);
  volume = props.Mass();
}

bool ExactShapeUtils::triangulate(const TopoDS_Shape& shape, double linearDeflection,
                                  double angularDeflection, Triangulation& result)
{
  try {
    BRepMesh_IncrementalMesh mesher(shape, linearDeflection, /*isRelative*/ Standard_False,
                                    angularDeflection, /*isInParallel*/ Standard_True);
    if (!mesher.IsDone()) return false;
  } catch (const Standard_Failure&) {
    return false;
  }

  for (TopExp_Explorer ex(shape, TopAbs_FACE); ex.More(); ex.Next()) {
    const TopoDS_Face& face = TopoDS::Face(ex.Current());
    TopLoc_Location location;
    const Handle(Poly_Triangulation) tri = BRep_Tool::Triangulation(face, location);
    if (tri.IsNull()) return false;  // a face BRepMesh could not handle: no closed mesh
    const gp_Trsf trsf = location.Transformation();
    const bool reversed = face.Orientation() == TopAbs_REVERSED;
    const int base = static_cast<int>(result.vertices.size());
    for (int i = 1; i <= tri->NbNodes(); ++i) {
      const gp_Pnt p = tri->Node(i).Transformed(trsf);
      result.vertices.emplace_back(p.X(), p.Y(), p.Z());
    }
    for (int i = 1; i <= tri->NbTriangles(); ++i) {
      int a, b, c;
      tri->Triangle(i).Get(a, b, c);
      if (reversed) std::swap(b, c);
      result.triangles.push_back({base + a - 1, base + b - 1, base + c - 1});
    }
  }
  return !result.triangles.empty();
}

#endif  // ENABLE_OCCT
