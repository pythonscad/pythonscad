#ifdef ENABLE_OCCT

#include "geometry/brep/BRepGeometry.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <array>
#include <cmath>
#include <filesystem>
#include <map>
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
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>
#include <Standard_Failure.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Shape.hxx>

#include "Feature.h"
#include "core/AST.h"
#include "core/CurveDiscretizer.h"
#include "geometry/Geometry.h"
#include "geometry/PolySetBuilder.h"
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

double volumeOf(const TopoDS_Shape& shape)
{
  GProp_GProps props;
  BRepGProp::VolumeProperties(shape, props);
  return props.Mass();
}

/// Read the written file back and compare its volume with the solid it was made
/// from. OCCT's STEP writer can silently drop parts of a solid in degenerate
/// configurations (see ExactShape::genericAxes()); a mesh export is better than
/// an exact file that is wrong.
enum class Roundtrip { Ok, Mismatch, Unreadable };

/// Read the written file back and compare its volume with the solid it was made
/// from. OCCT's STEP translation can silently keep the wrong side of a face in
/// degenerate configurations (see ExactShapeUtils::splitClosedFaces()).
Roundtrip checkRoundtrip(const TopoDS_Shape& shape, const std::string& file, double& actual)
{
  STEPControl_Reader reader;
  if (reader.ReadFile(file.c_str()) != IFSelect_RetDone || reader.TransferRoots() == 0) {
    return Roundtrip::Unreadable;
  }
  const double expected = volumeOf(shape);
  actual = volumeOf(reader.OneShape());
  return std::abs(actual - expected) <= 1e-6 * std::max(1.0, std::abs(expected)) ? Roundtrip::Ok
                                                                                 : Roundtrip::Mismatch;
}

bool writeStepFile(const TopoDS_Shape& shape, const std::string& file, const ExportInfo& exportInfo)
{
  STEPControl_Writer writer;
  Interface_Static::SetCVal("write.step.schema", "AP214IS");
  Interface_Static::SetCVal("write.step.unit", "MM");
  if (!exportInfo.title.empty()) {
    Interface_Static::SetCVal("write.step.product.name", exportInfo.title.c_str());
  }
  return writer.Transfer(shape, STEPControl_AsIs) == IFSelect_RetDone &&
         writer.Write(file.c_str()) == IFSelect_RetDone;
}

/// Write shape as STEP and verify the file by reading it back. If the check fails,
/// retry once with closed faces split; if that fails too, return false so the
/// caller exports the mesh instead of an exact file that is wrong.
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
    const double expected = volumeOf(shape);
    double actual = expected;
    Roundtrip result = Roundtrip::Unreadable;
    for (int attempt = 0; attempt < 2 && result != Roundtrip::Ok; ++attempt) {
      const TopoDS_Shape candidate = attempt == 0 ? shape : ExactShapeUtils::splitClosedFaces(shape);
      if (!writeStepFile(candidate, tmpfile.string(), exportInfo)) break;
      result = checkRoundtrip(candidate, tmpfile.string(), actual);
    }
    if (result == Roundtrip::Ok) {
      std::ifstream in(tmpfile, std::ios::binary);
      output << in.rdbuf();
      ok = in.good() || in.eof();
    } else if (result == Roundtrip::Mismatch) {
      LOG(
        message_group::Warning,
        "STEP export: exact file failed the self check (volume %1$.6g instead of %2$.6g), exporting as "
        "mesh",
        actual, expected);
    } else {
      LOG(message_group::Warning,
          "STEP export: exact file could not be written or read back, exporting as mesh");
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

bool enabled()
{
  return Feature::ExperimentalOcct.is_enabled();
}

std::unique_ptr<const Geometry> attachExact(std::unique_ptr<PolySet> mesh,
                                            std::shared_ptr<const ExactShape> exact)
{
  if (!enabled() || !mesh || !exact || mesh->isEmpty()) return mesh;
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
  if (!enabled() || !result || result->isEmpty()) return result;
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

std::unique_ptr<Geometry> importStep(const std::string& filename, const CurveDiscretizer& discretizer,
                                     const Location& loc)
{
  if (!enabled()) return nullptr;  // silently use the built-in importer
  TopoDS_Shape shape;
  std::string error;
  if (!ExactShapeUtils::readStep(filename, shape, error)) {
    LOG(message_group::Warning, loc, "",
        "OpenCASCADE could not import STEP file '%1$s' (%2$s), trying the built-in importer", filename,
        error);
    return nullptr;
  }

  // Mesh as fine as OpenSCAD would make a circle of the model's size: the angle
  // between segments (curvature) and the chord error that goes with it.
  const double diag = ExactShapeUtils::boundingDiagonal(shape);
  const double r = std::max(diag / 2, 1e-6);
  const int segments = std::max(3, discretizer.getCircularSegmentCount(r).value_or(30));
  const double angular = 2 * M_PI / segments;
  const double linear = std::max(r * (1 - std::cos(M_PI / segments)), 1e-4 * r);

  ExactShapeUtils::Triangulation tri;
  if (!ExactShapeUtils::triangulate(shape, linear, angular, tri)) {
    LOG(message_group::Warning, loc, "",
        "OpenCASCADE could not triangulate STEP file '%1$s', trying the built-in importer", filename);
    return nullptr;
  }

  // Faces are meshed separately; weld the vertices they share along their edges.
  const double grid = 1e-9 * std::max(1.0, diag);
  std::map<std::array<long long, 3>, int> welded;
  std::vector<int> remap(tri.vertices.size());
  std::vector<Vector3d> vertices;
  for (size_t i = 0; i < tri.vertices.size(); ++i) {
    const auto& v = tri.vertices[i];
    const std::array<long long, 3> key{std::llround(v[0] / grid), std::llround(v[1] / grid),
                                       std::llround(v[2] / grid)};
    auto [it, inserted] = welded.emplace(key, static_cast<int>(vertices.size()));
    if (inserted) vertices.push_back(v);
    remap[i] = it->second;
  }
  PolySetBuilder builder(static_cast<int>(vertices.size()), static_cast<int>(tri.triangles.size()));
  for (const auto& t : tri.triangles) {
    const int a = remap[t[0]], b = remap[t[1]], c = remap[t[2]];
    if (a == b || b == c || a == c) continue;  // degenerate after welding
    builder.appendPolygon(std::vector<Vector3d>{vertices[a], vertices[b], vertices[c]});
  }
  std::unique_ptr<PolySet> mesh = builder.build();
  mesh->setTriangular(true);

  int solids = 0;
  double volume = 0;
  ExactShapeUtils::solidInfo(shape, solids, volume);
  if (solids == 0 || volume <= 0) {
    LOG(message_group::Warning, loc, "",
        "STEP file '%1$s' contains no solids, imported as surface mesh without exact geometry",
        filename);
    return mesh;
  }
  return std::make_unique<BRepGeometry>(std::move(*mesh), ExactShape::fromShape(shape));
}

bool exportStepExact(const std::shared_ptr<const Geometry>& geom, std::ostream& output,
                     const ExportInfo& exportInfo)
{
  if (!enabled()) return false;
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
