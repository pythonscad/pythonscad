#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include "geometry/Geometry.h"
#include "geometry/PolySet.h"
#include "geometry/linalg.h"

// Width / space check between facing surfaces, the 3D counterpart of the DRC
// rules INTERNAL and EXTERNAL.
//
//   Internal: minimum wall thickness, measured through the material
//   External: minimum gap, measured through the air
//
// Both share one kernel: External runs the Internal kernel on the solid with
// all normals flipped, i.e. on its complement.
//
// Two triangles A and B "face" each other when their normals are at least
// 'min_angle_deg' apart and each lies behind the other one's plane. Both are
// clipped exactly (linear half-spaces) to the part that is closer than the
// limit to the other face's plane, optionally to the other face's projection
// prism, and the exact distance of the remaining convex polygons is computed
// (vertex-face and edge-edge). No sampling is involved.
namespace FacingCheck {

enum class Mode { Internal, External };

struct Options {
  // Minimum angle between the two face normals. 180 = exactly opposite. Must
  // be > 90, otherwise every convex box edge reports distance 0. With 120,
  // wedges sharper than 60 degrees count as thin.
  double min_angle_deg = 120.0;

  // Projection cone: 0 = strictly projecting (like PROJECTING in DRC decks),
  // values in between allow oblique measurement up to that angle, >= 90 means
  // any direction.
  double alpha_deg = 90.0;

  // Reject pairs whose connecting segment runs through the wrong region
  // (air for Internal, material for External), e.g. edge to edge across a
  // step. Costs nothing on clean designs.
  bool occlusion = true;

  // Keep the list of violating pairs. Without it only the summary
  // (count, minimum, per-triangle minimum) is produced in O(N) memory.
  bool store_pairs = true;

  // Keep the convex point set of each violation (needed for errorGeometry()).
  bool build_hulls = true;

  // Worker threads, 0 = hardware concurrency. Single threaded in WASM builds
  // without pthreads.
  unsigned threads = 0;
};

struct Violation {
  int tri_a = -1, tri_b = -1;  // triangle indices into Result::mesh
  Vector3d p, q;               // closest points on tri_a and tri_b
  double distance = 0;
  // Convex region of this violation:
  // Result::hull_points[hull_begin, hull_begin + hull_count)
  int hull_begin = 0, hull_count = 0;
};

struct Result {
  // The triangulated mesh all triangle indices refer to. For checkBetween()
  // it holds a's triangles first, then b's.
  std::shared_ptr<const PolySet> mesh;
  size_t split = 0;  // checkBetween(): number of triangles belonging to a

  std::vector<Violation> violations;     // sorted by (tri_a, tri_b)
  std::vector<Vector3d> hull_points;     // shared point pool for all hulls
  std::vector<double> tri_min_distance;  // per triangle, +inf if clean
  size_t violation_count = 0;            // counted even without store_pairs
  double min_distance;                   // +inf if clean

  bool clean() const { return violation_count == 0; }
};

// ps: closed, outward oriented solid. Non-triangular faces are tessellated.
Result check(const PolySet& ps, Mode mode, double distance, const Options& opt = {});

// External check between two solids only (like EXTERNAL layer1 layer2):
// gaps inside a or inside b are ignored.
Result checkBetween(const PolySet& a, const PolySet& b, double distance, const Options& opt = {});

// Error solid: union of the violation hulls, clipped to the material
// (Internal) or to the air (External) of 'bodies'. Each violating triangle
// contributes the hull to its nearest partner only (all pairs between two
// curved surfaces would overlap heavily and make the union explode). At most
// 'max_hulls' hulls are used, the smallest distances first. With grow > 0
// every vertex is moved outwards by about 'grow' along its vertex normal
// (cheap, not an exact offset) so the solid does not z-fight with the part. Returns nullptr if the
// result is clean or no error volume remains (e.g. only zero-volume edge contacts).
std::shared_ptr<const Geometry> errorGeometry(const Result& res, Mode mode,
                                              const std::vector<std::shared_ptr<const Geometry>>& bodies,
                                              size_t max_hulls = 20000, double grow = 0);

}  // namespace FacingCheck
