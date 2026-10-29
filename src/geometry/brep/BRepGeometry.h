#pragma once

#ifdef ENABLE_OCCT

#include <cstddef>
#include <memory>
#include <ostream>
#include <string>

#include "core/enums.h"
#include "geometry/Geometry.h"
#include "geometry/PolySet.h"
#include "geometry/brep/ExactShape.h"

struct ExportInfo;
class CurveDiscretizer;
class Location;
class LinearExtrudeNode;
class RotateExtrudeNode;

/**
 * Interface of geometries that carry an exact solid next to their mesh.
 *
 * The mesh stays the source of truth for everything mesh based (rendering,
 * mesh booleans, STL/3MF, DRC, ...). The exact solid is only read where it is
 * explicitly asked for, e.g. by the STEP exporter.
 */
class ExactCarrier
{
public:
  virtual ~ExactCarrier() = default;
  /// The exact solid, or nullptr if it cannot be trusted to describe the mesh anymore.
  [[nodiscard]] virtual std::shared_ptr<const ExactShape> exactShape() const = 0;
};

/**
 * A PolySet that additionally carries an exact solid.
 *
 * The PolySet part is the classic OpenSCAD tessellation (honouring $fn/$fa/$fs),
 * so every existing consumer keeps working unchanged and produces bit-identical
 * results. The exact solid follows copy() and transform(). Code that edits the
 * vertices directly (bypassing transform()) would make mesh and solid disagree;
 * a fingerprint of the mesh detects that and the exact solid is then ignored.
 */
class BRepGeometry : public PolySet, public ExactCarrier
{
public:
  BRepGeometry(PolySet mesh, std::shared_ptr<const ExactShape> exact);

  [[nodiscard]] std::unique_ptr<Geometry> copy() const override;
  [[nodiscard]] size_t memsize() const override;
  [[nodiscard]] std::string dump() const override;
  // resize() is inherited from PolySet and routed through transform()
  void transform(const Transform3d& mat) override;

  [[nodiscard]] std::shared_ptr<const ExactShape> exactShape() const override;

private:
  [[nodiscard]] size_t meshFingerprint() const;

  std::shared_ptr<const ExactShape> exact_;
  size_t fingerprint_ = 0;
};

/**
 * A Polygon2d that additionally carries an exact profile: planar face(s) in the
 * plane z = 0, in the same local frame as untransformedOutlines(). A pending 3D
 * placement (getTransform3d()) applies to both alike. Like BRepGeometry, a
 * fingerprint of the outlines detects changes that bypass BRepUtils::transform2d();
 * the exact profile is then ignored and the outlines are used instead.
 */
class ExactPolygon2d : public Polygon2d, public ExactCarrier
{
public:
  ExactPolygon2d(Polygon2d polygon, std::shared_ptr<const ExactShape> exact);

  [[nodiscard]] std::unique_ptr<Geometry> copy() const override;
  [[nodiscard]] size_t memsize() const override;
  [[nodiscard]] std::shared_ptr<const ExactShape> exactShape() const override;

  /// Replace the exact profile after the outlines were changed on purpose
  void setExact(std::shared_ptr<const ExactShape> exact);
  [[nodiscard]] size_t outlineFingerprint() const;

private:
  std::shared_ptr<const ExactShape> exact_;
  size_t fingerprint_ = 0;
};

namespace BRepUtils {

/// True if exact geometry is switched on (experimental feature "occt", off by default).
/// When off, no exact solids are created or imported and STEP export is faceted,
/// exactly as in a build without OpenCASCADE.
bool enabled();

/// Attach an exact solid to a freshly built primitive mesh.
std::unique_ptr<const Geometry> attachExact(std::unique_ptr<PolySet> mesh,
                                            std::shared_ptr<const ExactShape> exact);

/// The exact solid of geom, if it has a trustworthy one.
std::shared_ptr<const ExactShape> exactOf(const Geometry *geom);

/**
 * geom converted to another representation of the same solid (e.g. Manifold ->
 * PolySet, or a tessellation): if src carries an exact solid, return converted
 * with that solid attached, else converted unchanged.
 */
std::shared_ptr<const Geometry> carryExact(const Geometry& src,
                                           std::shared_ptr<const Geometry> converted);
std::shared_ptr<const PolySet> carryExact(const Geometry& src, std::shared_ptr<const PolySet> converted);

/**
 * Attach the exact boolean of the children to result, the mesh result of the
 * same operation. Returns result unchanged unless all relevant children are
 * exact. Handles UNION, DIFFERENCE and INTERSECTION.
 */
std::shared_ptr<const Geometry> attachBoolean(OpenSCADOperator op, const Geometry::Geometries& children,
                                              std::shared_ptr<const Geometry> result);

// --- 2D -----------------------------------------------------------------------

/// Attach an exact profile to a freshly built 2D primitive.
std::unique_ptr<const Geometry> attachExact2D(std::unique_ptr<Polygon2d> polygon,
                                              std::shared_ptr<const ExactShape> exact);

/// converted is src reshaped without changing the region (e.g. sanitized, copied):
/// attach src's exact profile to it if both share the same frame.
std::unique_ptr<Polygon2d> carryExact2D(const Polygon2d& src, std::unique_ptr<Polygon2d> converted);

/// polygon.transform3d(mat), keeping an exact profile in step with the outlines.
void transform2d(Polygon2d& polygon, const Transform3d& mat);

/// Attach the exact 2D boolean of the children to result (UNION, DIFFERENCE, INTERSECTION).
std::unique_ptr<Polygon2d> attachBoolean2D(OpenSCADOperator op,
                                           const std::vector<std::shared_ptr<const Polygon2d>>& children,
                                           std::unique_ptr<Polygon2d> result);

// --- extrusions -------------------------------------------------------------------

/// Attach the exact linear extrusion of poly to mesh, the mesh extrusion of it,
/// if the extrusion is one OpenCASCADE can represent exactly.
std::unique_ptr<Geometry> linearExtrude(const LinearExtrudeNode& node, const Polygon2d& poly,
                                        std::unique_ptr<Geometry> mesh);
/// Same for rotate_extrude.
std::unique_ptr<Geometry> rotateExtrude(const RotateExtrudeNode& node, const Polygon2d& poly,
                                        std::unique_ptr<Geometry> mesh);

/**
 * Import a STEP file with OpenCASCADE. Solids keep their exact shape (so a STEP
 * export writes them exactly again) and get a triangle mesh for everything mesh
 * based; its fineness follows $fn/$fa/$fs like a circle of the model's size.
 * Returns nullptr if OpenCASCADE cannot read the file, so the caller can fall
 * back to the built-in importer.
 */
std::unique_ptr<Geometry> importStep(const std::string& filename, const CurveDiscretizer& discretizer,
                                     const Location& loc);

/**
 * Write geom as exact STEP if it consists only of exact geometry (a single
 * object, or a GeometryList of them). Returns false without writing anything
 * otherwise, so the caller can fall back to the mesh based exporter.
 */
bool exportStepExact(const std::shared_ptr<const Geometry>& geom, std::ostream& output,
                     const ExportInfo& exportInfo);

}  // namespace BRepUtils

#endif  // ENABLE_OCCT
