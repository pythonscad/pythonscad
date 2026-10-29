#pragma once

#ifdef ENABLE_OCCT

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "core/enums.h"
#include "geometry/linalg.h"

class TopoDS_Shape;

/**
 * An exact (B-Rep) solid that is only computed when somebody asks for it.
 *
 * An ExactShape is a node in a small recipe tree: a primitive, a transformed
 * child, or a boolean of children. Building the tree is cheap and touches no
 * OpenCASCADE code; shape() evaluates the recipe on first use (thread safe),
 * caches the result and releases the recipe. Rendering and mesh exports never
 * call shape(), so they pay nothing for exact geometry.
 *
 * Nodes are immutable and shared between geometries (cache, copies), which is
 * why all factories return shared_ptr<const ExactShape>.
 */
class ExactShape
{
public:
  ~ExactShape();
  ExactShape(const ExactShape&) = delete;
  ExactShape& operator=(const ExactShape&) = delete;

  /// Sphere of radius r around the origin
  static std::shared_ptr<const ExactShape> sphere(double r);
  /// Axis aligned box between the corners p1 and p2
  static std::shared_ptr<const ExactShape> box(const Vector3d& p1, const Vector3d& p2);
  /// Cylinder/cone along +z from z1 to z2, radius r1 at z1 and r2 at z2 (either may be 0),
  /// restricted to the sector 0..angle degrees around +z (angle in (0, 360])
  static std::shared_ptr<const ExactShape> cylinder(double r1, double r2, double z1, double z2,
                                                    double angle);
  /// An already computed solid, e.g. read from a STEP file
  static std::shared_ptr<const ExactShape> fromShape(const TopoDS_Shape& shape);
  /// child transformed by mat
  static std::shared_ptr<const ExactShape> transformed(std::shared_ptr<const ExactShape> child,
                                                       const Transform3d& mat);
  /// UNION, DIFFERENCE (first minus all others) or INTERSECTION of children; nullptr for other ops
  static std::shared_ptr<const ExactShape> boolean(
    OpenSCADOperator op, std::vector<std::shared_ptr<const ExactShape>> children);

  /// The exact solid, computed on first call. nullptr if OpenCASCADE could not build it
  /// (a warning is logged once in that case).
  [[nodiscard]] const TopoDS_Shape *shape() const;

private:
  friend struct ExactShapeFactory;
  struct Impl;
  explicit ExactShape(std::unique_ptr<Impl> impl);
  std::unique_ptr<Impl> impl_;
};

namespace ExactShapeUtils {
/**
 * Split closed (periodic) faces such as a full sphere or cylinder into several
 * non-closed faces. Used as a fallback for STEP export: on a closed surface a
 * single boundary loop is ambiguous, and when the loop runs through a pole
 * (e.g. a sphere centered on the edge of a plate) OCCT's STEP translation keeps
 * the wrong side of it. Returns shape unchanged if splitting fails.
 */
TopoDS_Shape splitClosedFaces(const TopoDS_Shape& shape);

/**
 * Read a STEP file with OpenCASCADE's reader (units converted to mm).
 * Returns false and sets error if the file cannot be read or holds no geometry.
 */
bool readStep(const std::string& filename, TopoDS_Shape& shape, std::string& error);

/// Length of the diagonal of the bounding box of shape (0 if empty)
double boundingDiagonal(const TopoDS_Shape& shape);

/// Number of solids in shape and their total volume
void solidInfo(const TopoDS_Shape& shape, int& solids, double& volume);

/// A triangle mesh of a shape. Triangles are oriented counter-clockwise seen from
/// outside, like OpenSCAD's PolySets.
struct Triangulation {
  std::vector<Vector3d> vertices;
  std::vector<std::array<int, 3>> triangles;
};

/**
 * Triangulate shape with BRepMesh. linearDeflection is the maximum chord distance
 * in mm, angularDeflection the maximum angle between adjacent segments in radians.
 */
bool triangulate(const TopoDS_Shape& shape, double linearDeflection, double angularDeflection,
                 Triangulation& result);
}  // namespace ExactShapeUtils

#endif  // ENABLE_OCCT
