#ifdef ENABLE_OCCT

#include "geometry/brep/BRepGeometry.h"

#include <atomic>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <sstream>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include <boost/functional/hash.hpp>

#include <BRepBuilderAPI_GTransform.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRep_Builder.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Interface_Static.hxx>
#include <STEPControl_Writer.hxx>
#include <Standard_Failure.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_GTrsf.hxx>
#include <gp_Mat.hxx>
#include <gp_Trsf.hxx>
#include <gp_XYZ.hxx>

#include "geometry/Geometry.h"
#include "io/export.h"
#include "utils/printutils.h"

namespace {

/// Apply an OpenSCAD affine transform to an OCCT shape.
/// Similarities (rotation, translation, uniform scale, mirror) keep the analytic
/// surfaces; anything else (non-uniform scale, shear) converts them to B-splines,
/// which is still exact, just no longer e.g. a SPHERICAL_SURFACE.
/// Returns nullptr if OCCT cannot represent the transform (e.g. scale 0).
std::shared_ptr<const TopoDS_Shape> transformShape(const TopoDS_Shape& shape, const Transform3d& mat)
{
  const Eigen::Matrix3d lin = mat.linear();
  const Eigen::Vector3d tr = mat.translation();

  const Eigen::Matrix3d ltl = lin.transpose() * lin;
  const double s2 = ltl.trace() / 3.0;
  const bool similarity =
    s2 > 0 && (ltl - s2 * Eigen::Matrix3d::Identity()).norm() <= 1e-9 * std::max(1.0, s2);

  try {
    if (similarity) {
      gp_Trsf trsf;
      trsf.SetValues(lin(0, 0), lin(0, 1), lin(0, 2), tr[0],  //
                     lin(1, 0), lin(1, 1), lin(1, 2), tr[1],  //
                     lin(2, 0), lin(2, 1), lin(2, 2), tr[2]);
      BRepBuilderAPI_Transform op(shape, trsf, /*copy*/ true);
      if (!op.IsDone()) return nullptr;
      return std::make_shared<const TopoDS_Shape>(op.Shape());
    }
    if (std::abs(lin.determinant()) < 1e-12) return nullptr;
    gp_GTrsf gtrsf;
    gtrsf.SetVectorialPart(gp_Mat(lin(0, 0), lin(0, 1), lin(0, 2),  //
                                  lin(1, 0), lin(1, 1), lin(1, 2),  //
                                  lin(2, 0), lin(2, 1), lin(2, 2)));
    gtrsf.SetTranslationPart(gp_XYZ(tr[0], tr[1], tr[2]));
    BRepBuilderAPI_GTransform op(shape, gtrsf, /*copy*/ true);
    if (!op.IsDone()) return nullptr;
    return std::make_shared<const TopoDS_Shape>(op.Shape());
  } catch (const Standard_Failure& e) {
    LOG(message_group::Warning, "B-Rep transform failed (%1$s), keeping mesh only",
        e.GetMessageString() ? e.GetMessageString() : "unknown error");
    return nullptr;
  }
}

/// Collect the exact shapes of geom. Fails if any part is not a consistent BRepGeometry.
bool collectShapes(const std::shared_ptr<const Geometry>& geom,
                   std::vector<const TopoDS_Shape *>& shapes, bool& anyExact)
{
  if (!geom || geom->isEmpty()) return true;
  if (const auto list = std::dynamic_pointer_cast<const GeometryList>(geom)) {
    bool ok = true;
    for (const auto& item : list->flatten()) {
      ok &= collectShapes(item.second, shapes, anyExact);
    }
    return ok;
  }
  if (const auto brep = std::dynamic_pointer_cast<const BRepGeometry>(geom)) {
    if (brep->hasShape() && brep->isConsistent()) {
      anyExact = true;
      shapes.push_back(&brep->shape());
      return true;
    }
  }
  return false;
}

/// OCCT's STEP writer only writes to a file name (stream output appeared in 7.7),
/// so go through a temporary file to stay compatible with older distro packages.
bool writeStepToStream(const TopoDS_Shape& shape, std::ostream& output, const ExportInfo& exportInfo)
{
  static std::atomic<unsigned> counter{0};
  std::error_code ec;
  const auto tmpdir = std::filesystem::temp_directory_path(ec);
  if (ec) return false;
  const auto tmpfile =
    tmpdir / ("pythonscad-" + std::to_string(std::hash<std::string>{}(exportInfo.sourceFilePath)) + "-" +
              std::to_string(counter++) + ".step");

  bool ok = false;
  try {
    STEPControl_Writer writer;
    Interface_Static::SetCVal("write.step.schema", "AP214IS");
    Interface_Static::SetCVal("write.step.unit", "MM");
    if (!exportInfo.title.empty()) {
      Interface_Static::SetCVal("write.step.product.name", exportInfo.title.c_str());
    }
    if (writer.Transfer(shape, STEPControl_AsIs) == IFSelect_RetDone &&
        writer.Write(tmpfile.string().c_str()) == IFSelect_RetDone) {
      std::ifstream in(tmpfile, std::ios::binary);
      output << in.rdbuf();
      ok = in.good() || in.eof();
    }
  } catch (const Standard_Failure& e) {
    LOG(message_group::Error, "STEP export failed: %1$s",
        e.GetMessageString() ? e.GetMessageString() : "unknown error");
    ok = false;
  }
  std::filesystem::remove(tmpfile, ec);
  return ok;
}

}  // namespace

BRepGeometry::BRepGeometry(PolySet mesh, std::shared_ptr<const TopoDS_Shape> shape)
  : PolySet(std::move(mesh)), shape_(std::move(shape))
{
  fingerprint_ = meshFingerprint();
}

std::unique_ptr<Geometry> BRepGeometry::copy() const
{
  return std::make_unique<BRepGeometry>(*this);
}

size_t BRepGeometry::memsize() const
{
  // The shape is shared between copies and small for primitives; count a nominal size.
  return PolySet::memsize() + sizeof(*this) - sizeof(PolySet) + (shape_ ? 1024 : 0);
}

std::string BRepGeometry::dump() const
{
  std::ostringstream out;
  out << "BRepGeometry:\n exact shape: " << (hasShape() ? "yes" : "no") << "\n" << PolySet::dump();
  return out.str();
}

void BRepGeometry::transform(const Transform3d& mat)
{
  // Only follow the transform with the shape if it still describes this mesh.
  const bool wasConsistent = isConsistent();
  PolySet::transform(mat);
  if (shape_ && wasConsistent) shape_ = transformShape(*shape_, mat);
  else shape_.reset();
  fingerprint_ = meshFingerprint();
}

bool BRepGeometry::isConsistent() const
{
  return meshFingerprint() == fingerprint_;
}

size_t BRepGeometry::meshFingerprint() const
{
  size_t seed = vertices.size();
  for (const auto& v : vertices) {
    boost::hash_combine(seed, v[0]);
    boost::hash_combine(seed, v[1]);
    boost::hash_combine(seed, v[2]);
  }
  boost::hash_combine(seed, indices.size());
  for (const auto& face : indices) {
    for (const auto idx : face) boost::hash_combine(seed, idx);
  }
  return seed;
}

namespace BRepUtils {

std::unique_ptr<const Geometry> createSphere(std::unique_ptr<PolySet> mesh, double r)
{
  std::shared_ptr<const TopoDS_Shape> shape;
  try {
    shape = std::make_shared<const TopoDS_Shape>(BRepPrimAPI_MakeSphere(r).Shape());
  } catch (const Standard_Failure& e) {
    LOG(message_group::Warning, "Could not create exact sphere (%1$s), using mesh only",
        e.GetMessageString() ? e.GetMessageString() : "unknown error");
    return mesh;
  }
  return std::make_unique<BRepGeometry>(std::move(*mesh), std::move(shape));
}

bool exportStepExact(const std::shared_ptr<const Geometry>& geom, std::ostream& output,
                     const ExportInfo& exportInfo)
{
  std::vector<const TopoDS_Shape *> shapes;
  bool anyExact = false;
  if (!collectShapes(geom, shapes, anyExact) || shapes.empty()) {
    if (anyExact) {
      LOG(message_group::Warning,
          "STEP export: model mixes exact and mesh geometry, exporting everything as mesh");
    }
    return false;
  }

  if (shapes.size() == 1) return writeStepToStream(*shapes.front(), output, exportInfo);

  TopoDS_Compound compound;
  BRep_Builder builder;
  builder.MakeCompound(compound);
  for (const auto *s : shapes) builder.Add(compound, *s);
  return writeStepToStream(compound, output, exportInfo);
}

}  // namespace BRepUtils

#endif  // ENABLE_OCCT
