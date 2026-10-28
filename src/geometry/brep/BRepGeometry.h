#pragma once

#ifdef ENABLE_OCCT

#include <cstddef>
#include <memory>
#include <ostream>
#include <string>

#include "geometry/PolySet.h"

class TopoDS_Shape;
struct ExportInfo;

/**
 * A PolySet that additionally carries an exact B-Rep shape (OpenCASCADE).
 *
 * The PolySet part is the classic OpenSCAD tessellation (honouring $fn/$fa/$fs),
 * so every existing consumer (renderer, Manifold/CGAL booleans, STL/3MF export,
 * DRC checks, ...) keeps working unchanged and produces bit-identical results.
 *
 * The exact shape is only used where it is explicitly asked for, e.g. by the
 * STEP exporter. It survives copy() and transform(); any operation that builds
 * a new geometry (union, difference, extrude, ...) naturally yields a plain
 * PolySet again and the exact information is dropped.
 *
 * Code that edits the vertices of a copy directly (bypassing transform()) would
 * make mesh and shape disagree. isConsistent() detects that via a fingerprint
 * of the mesh, so consumers can fall back to the mesh in that case.
 */
class BRepGeometry : public PolySet
{
public:
  BRepGeometry(PolySet mesh, std::shared_ptr<const TopoDS_Shape> shape);

  [[nodiscard]] std::unique_ptr<Geometry> copy() const override;
  [[nodiscard]] size_t memsize() const override;
  [[nodiscard]] std::string dump() const override;
  // resize() is inherited from PolySet and routed through transform()
  void transform(const Transform3d& mat) override;

  /// True if an exact shape is present (a degenerate transform, e.g. scale 0, drops it)
  [[nodiscard]] bool hasShape() const { return shape_ != nullptr; }
  /// Exact shape; only valid if hasShape()
  [[nodiscard]] const TopoDS_Shape& shape() const { return *shape_; }
  /// True if the exact shape can be trusted to describe the same object as the mesh
  [[nodiscard]] bool isConsistent() const;

private:
  [[nodiscard]] size_t meshFingerprint() const;

  std::shared_ptr<const TopoDS_Shape> shape_;
  size_t fingerprint_ = 0;
};

namespace BRepUtils {

/// Wrap an already tessellated sphere mesh together with an exact sphere of radius r at the origin.
std::unique_ptr<const Geometry> createSphere(std::unique_ptr<PolySet> mesh, double r);

/**
 * Write geom as exact STEP if it consists only of consistent BRepGeometry objects
 * (a single one, or a GeometryList of them). Returns false without writing anything
 * otherwise, so the caller can fall back to the mesh based exporter.
 */
bool exportStepExact(const std::shared_ptr<const Geometry>& geom, std::ostream& output,
                     const ExportInfo& exportInfo);

}  // namespace BRepUtils

#endif  // ENABLE_OCCT
