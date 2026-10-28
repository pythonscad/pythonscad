#ifdef ENABLE_OCCT

#include "geometry/brep/BRepGeometry.h"

#include <atomic>
#include <filesystem>
#include <fstream>
#include <functional>
#include <sstream>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include <boost/functional/hash.hpp>

#include <BRep_Builder.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Interface_Static.hxx>
#include <STEPControl_Writer.hxx>
#include <Standard_Failure.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Shape.hxx>

#include "geometry/Geometry.h"
#include "io/export.h"
#include "utils/printutils.h"
#ifdef ENABLE_MANIFOLD
#include "geometry/brep/BRepManifoldGeometry.h"
#endif

namespace {

/// Collect the exact shapes of geom (evaluating them). Fails if any non-empty part
/// has no exact solid, or if OpenCASCADE could not compute it.
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
  if (const auto exact = BRepUtils::exactOf(geom.get())) {
    anyExact = true;
    if (const TopoDS_Shape *s = exact->shape()) {
      shapes.push_back(s);
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

BRepGeometry::BRepGeometry(PolySet mesh, std::shared_ptr<const ExactShape> exact)
  : PolySet(std::move(mesh)), exact_(std::move(exact))
{
  fingerprint_ = meshFingerprint();
}

std::unique_ptr<Geometry> BRepGeometry::copy() const
{
  return std::make_unique<BRepGeometry>(*this);
}

size_t BRepGeometry::memsize() const
{
  return PolySet::memsize() + sizeof(*this) - sizeof(PolySet);
}

std::string BRepGeometry::dump() const
{
  std::ostringstream out;
  out << "BRepGeometry:\n exact solid: " << (exactShape() ? "yes" : "no") << "\n" << PolySet::dump();
  return out.str();
}

void BRepGeometry::transform(const Transform3d& mat)
{
  // Only follow the transform with the solid if it still describes this mesh.
  auto exact = exactShape();
  PolySet::transform(mat);
  exact_ = ExactShape::transformed(std::move(exact), mat);
  fingerprint_ = meshFingerprint();
}

std::shared_ptr<const ExactShape> BRepGeometry::exactShape() const
{
  if (!exact_ || meshFingerprint() != fingerprint_) return nullptr;
  return exact_;
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

void logWarning(const std::string& msg)
{
  LOG(message_group::Warning, "%1$s", msg);
}

std::unique_ptr<const Geometry> attachExact(std::unique_ptr<PolySet> mesh,
                                            std::shared_ptr<const ExactShape> exact)
{
  if (!mesh || !exact || mesh->isEmpty()) return mesh;
  return std::make_unique<BRepGeometry>(std::move(*mesh), std::move(exact));
}

std::shared_ptr<const ExactShape> exactOf(const Geometry *geom)
{
  const auto *carrier = dynamic_cast<const ExactCarrier *>(geom);
  return carrier ? carrier->exactShape() : nullptr;
}

namespace {

std::shared_ptr<const Geometry> wrap(const std::shared_ptr<const Geometry>& geom,
                                     std::shared_ptr<const ExactShape> exact)
{
  if (!geom || !exact || dynamic_cast<const ExactCarrier *>(geom.get())) return geom;
#ifdef ENABLE_MANIFOLD
  if (const auto *mani = dynamic_cast<const ManifoldGeometry *>(geom.get())) {
    return std::make_shared<const BRepManifoldGeometry>(*mani, std::move(exact));
  }
#endif
  if (const auto *ps = dynamic_cast<const PolySet *>(geom.get())) {
    return std::make_shared<const BRepGeometry>(*ps, std::move(exact));
  }
  return geom;  // e.g. CGAL Nef polyhedra: no exact companion (yet)
}

}  // namespace

std::shared_ptr<const Geometry> carryExact(const Geometry& src,
                                           std::shared_ptr<const Geometry> converted)
{
  return wrap(converted, exactOf(&src));
}

std::shared_ptr<const PolySet> carryExact(const Geometry& src, std::shared_ptr<const PolySet> converted)
{
  return std::static_pointer_cast<const PolySet>(wrap(converted, exactOf(&src)));
}

std::shared_ptr<const Geometry> attachBoolean(OpenSCADOperator op, const Geometry::Geometries& children,
                                              std::shared_ptr<const Geometry> result)
{
  if (!result || result->isEmpty()) return result;
  if (op != OpenSCADOperator::UNION && op != OpenSCADOperator::DIFFERENCE &&
      op != OpenSCADOperator::INTERSECTION) {
    return result;
  }

  // Same child semantics as the mesh booleans: union and difference ignore empty
  // children (except an empty first child of a difference, which empties the
  // result), intersection with an empty child is empty.
  std::vector<std::shared_ptr<const ExactShape>> exactChildren;
  bool first = true;
  for (const auto& item : children) {
    const Geometry *g = item.second.get();
    const bool empty = g == nullptr || g->isEmpty();
    if (empty) {
      if (op == OpenSCADOperator::INTERSECTION) return result;
      if (op == OpenSCADOperator::DIFFERENCE && first) return result;
      first = false;
      continue;
    }
    first = false;
    auto exact = exactOf(g);
    if (!exact) return result;  // a mesh-only child: the result is mesh-only too
    exactChildren.push_back(std::move(exact));
  }
  return wrap(result, ExactShape::boolean(op, std::move(exactChildren)));
}

bool exportStepExact(const std::shared_ptr<const Geometry>& geom, std::ostream& output,
                     const ExportInfo& exportInfo)
{
  std::vector<const TopoDS_Shape *> shapes;
  bool anyExact = false;
  if (!collectShapes(geom, shapes, anyExact) || shapes.empty()) {
    if (anyExact) {
      LOG(message_group::Warning,
          "STEP export: model is not completely exact (mesh-only parts or failed exact computation), "
          "exporting everything as mesh");
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
