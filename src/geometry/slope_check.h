#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include "geometry/Geometry.h"
#include "geometry/PolySet.h"
#include "geometry/linalg.h"

// Slope check: the angle of every face against a direction, the 3D
// counterpart of an ANGLE rule. Overhang (3D printing) and draft
// (demolding, e.g. sand casting) are the same local test with different
// windows:
//
//   beta = asin(n . u)   n: face normal, u: pull / build direction
//
//   beta =  0   face parallel to u (vertical wall)
//   beta = 90   face looks along u, beta = -90 face looks against u
//
//   overhang(45): u = build direction (up), beta >= -45 (base skipped)
//   draft(2):     u = pull direction of the face's mold half, beta >= 2
//
// Draft additionally needs a global test: a face can have the right angle
// and still be hidden behind other material along u (C-profile lying on its
// side). With 'undercut' every face piece is checked for occluders in its
// pull direction; the shadow is computed exactly by clipping the piece
// against the projection prisms of the candidate occluders.
namespace SlopeCheck {

enum class Parting {
  None,   // one-sided: every face is checked against +dir
  Free,   // two mold halves, each face goes to the half its normal points to
  Plane,  // two mold halves split at the plane dir . x = parting_pos
};

struct Options {
  Vector3d dir = Vector3d(0, 0, 1);  // normalized internally
  double min_deg = -90;
  double max_deg = 90;
  Parting parting = Parting::None;
  double parting_pos = 0;
  bool undercut = false;   // global visibility test along the pull direction
  bool skip_base = false;  // ignore faces on the lowest level along dir (build plate)
  unsigned threads = 0;    // 0 = hardware concurrency
};

enum Flags { kAngle = 1, kUndercut = 2 };

// A violating face piece (a whole triangle, or the part on one side of the
// parting plane).
struct Violation {
  int tri = -1;
  int flags = 0;
  double beta_deg = 0;
  Vector3d pull;  // pull direction the piece was checked against
  // Points in Result::hull_points: the piece polygon itself, and for
  // undercuts one convex column per nearest occluder (the trapped region
  // between the piece and that occluder).
  int piece_begin = 0, piece_count = 0;
  std::vector<std::pair<int, int>> columns;  // (begin, count)
};

struct Result {
  std::shared_ptr<const PolySet> mesh;  // triangulated, indices refer to it
  std::vector<Violation> violations;    // sorted by tri
  std::vector<Vector3d> hull_points;
  size_t angle_count = 0;     // triangles with an angle violation
  size_t undercut_count = 0;  // triangles (partly) hidden along their pull direction
  size_t count = 0;           // distinct violating triangles
  double worst_deg = 0;       // beta of the worst angle violation (valid if angle_count > 0)
  bool clean() const { return count == 0; }
};

Result check(const PolySet& ps, const Options& opt);

// Error solid: angle violations as a skin of 'thickness' under the violating
// faces, undercuts as the trapped air between a face and its occluders.
// 'grow' lifts the result off the part (see FacingCheck::errorGeometry).
std::shared_ptr<const Geometry> errorGeometry(const Result& res,
                                              const std::shared_ptr<const Geometry>& body,
                                              double thickness, double grow, size_t max_hulls = 20000);

}  // namespace SlopeCheck
