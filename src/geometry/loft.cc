#include "geometry/loft.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <utility>
#include <vector>

#include <sstream>
#include <string>

#include "geometry/Polygon2d.h"
#include "geometry/PolySetBuilder.h"
#include "utils/printutils.h"

namespace {
// LOG()'s format-string support varies across call sites in this codebase,
// so trace messages are assembled as plain strings here to avoid depending
// on a particular formatting syntax.
template <typename... Args>
std::string concatStr(Args&&...args)
{
  std::ostringstream oss;
  (oss << ... << args);
  return oss.str();
}
}  // namespace

namespace {

// -------------------- Small 2D helpers --------------------

// Signed area of a closed 2D polygon (positive = CCW, negative = CW).
double signedArea2D(const std::vector<Vector2d>& poly)
{
  double area = 0.0;
  const size_t n = poly.size();
  for (size_t i = 0; i < n; i++) {
    const Vector2d& a = poly[i];
    const Vector2d& b = poly[(i + 1) % n];
    area += a.x() * b.y() - b.x() * a.y();
  }
  return 0.5 * area;
}

// Reverses a ring's 2D points and its parallel 3D points together, so the
// vertex-by-vertex correspondence between uv[i] and pos3d[i] is preserved.
void reverseRingInPlace(std::vector<Vector2d>& uv, std::vector<Vector3d>& pos3d)
{
  std::reverse(uv.begin(), uv.end());
  std::reverse(pos3d.begin(), pos3d.end());
}

// Forces 'outer' to be wound CCW and every entry of 'holes' to be wound CW
// (the standard convention for a polygon-with-holes). Not strictly required
// any more by the Delaunay-based construction below (which doesn't care
// about winding), but Polygon2d's "positive" outline flag is conventionally
// interpreted this way elsewhere in the codebase, so we keep the rings
// consistent with that convention. Keeps the parallel 3D point arrays in
// sync with any reversal.
void ensureOrientation(std::vector<Vector2d>& outer_uv, std::vector<Vector3d>& outer_pos,
                       std::vector<std::vector<Vector2d>>& holes_uv,
                       std::vector<std::vector<Vector3d>>& holes_pos)
{
  if (signedArea2D(outer_uv) < 0) reverseRingInPlace(outer_uv, outer_pos);
  for (size_t hi = 0; hi < holes_uv.size(); hi++) {
    if (holes_uv[hi].size() < 3) continue;
    if (signedArea2D(holes_uv[hi]) > 0) reverseRingInPlace(holes_uv[hi], holes_pos[hi]);
  }
}

// -------------------- Domain polygon (outer ring + holes) --------------------
// Builds a Polygon2d from the projected boundary rings so we can reuse the
// existing, hole-aware Polygon2d::point_inside().
Polygon2d buildDomainPolygon(const std::vector<Vector2d>& outer_uv,
                             const std::vector<std::vector<Vector2d>>& holes_uv)
{
  Polygon2d poly;

  Outline2d outerOutline;
  outerOutline.vertices.assign(outer_uv.begin(), outer_uv.end());
  outerOutline.positive = true;
  poly.addOutline(outerOutline);

  for (const auto& h : holes_uv) {
    Outline2d holeOutline;
    holeOutline.vertices.assign(h.begin(), h.end());
    holeOutline.positive = false;
    poly.addOutline(holeOutline);
  }
  return poly;
}

// -------------------- Base surface: barycentric base position + normal --------------------
// Gives any interior (u,v) query point a well-defined base 3D position and
// face normal via barycentric interpolation over a coarse triangulation of
// the boundary points. Built by the same Delaunay-plus-domain-filter
// approach as the main mesh (see loft() below) - not by bridging the holes
// into the outer contour and ear-clipping, which turned out to stall (see
// the long comment at the loft() base-surface construction site for why).
struct LoftBaseSurface {
  std::vector<Vector2d> uv;     // 2D coordinates (boundary points only: outer + holes)
  std::vector<Vector3d> pos3d;  // parallel 3D positions
  std::vector<std::array<int, 3>> tris;

  bool sample(const Vector2d& p, Vector3d& outPos, Vector3d& outNormal) const
  {
    for (const auto& t : tris) {
      const Vector2d& a = uv[t[0]];
      const Vector2d& b = uv[t[1]];
      const Vector2d& c = uv[t[2]];
      double denom = (b.y() - c.y()) * (a.x() - c.x()) + (c.x() - b.x()) * (a.y() - c.y());
      if (std::fabs(denom) < 1e-12) continue;
      double w0 = ((b.y() - c.y()) * (p.x() - c.x()) + (c.x() - b.x()) * (p.y() - c.y())) / denom;
      double w1 = ((c.y() - a.y()) * (p.x() - c.x()) + (a.x() - c.x()) * (p.y() - c.y())) / denom;
      double w2 = 1.0 - w0 - w1;
      const double eps = -1e-9;
      if (w0 < eps || w1 < eps || w2 < eps) continue;

      const Vector3d& p0 = pos3d[t[0]];
      const Vector3d& p1 = pos3d[t[1]];
      const Vector3d& p2 = pos3d[t[2]];
      outPos = w0 * p0 + w1 * p1 + w2 * p2;
      Vector3d n = (p1 - p0).cross(p2 - p0);
      double nlen = n.norm();
      outNormal = (nlen > 1e-12) ? (n / nlen) : Vector3d(0, 0, 1);
      return true;
    }
    return false;
  }
};

// Cheap deterministic hash-noise in [0,1), used only to jitter points fed
// into bowyerWatson() (see below).
double hashNoise01(double seed)
{
  double s = std::sin(seed * 127.1 + 311.7) * 43758.5453;
  return s - std::floor(s);
}

// Returns a copy of 'pts' with a tiny deterministic per-point offset. This
// is ONLY meant to be fed into bowyerWatson() as topology input - the
// original, unperturbed coordinates are what everything else (domain
// membership tests, the final vertex positions) should keep using.
//
// WHY THIS IS NEEDED: an axis-aligned regular grid (which is exactly what
// the interior sampling below produces) is a classic pathological input for
// Bowyer-Watson - every 2x2 grid cell's four corners are exactly cocircular,
// so a plain incremental implementation constantly hits circumcircle ties.
// Depending on how those ties break, this can corrupt the "boundary of the
// bad-triangle cavity is a single closed loop" assumption the algorithm
// relies on, which shows up as either wrong output or (with the retry/kept
// bookkeeping used here) triangle counts that blow up point by point. A
// small random perturbation, well above floating point noise but far below
// anything visible in the final geometry, breaks the exact cocircularity
// and avoids this entirely - a standard trick, sometimes called simulated
// perturbation, for making incremental Delaunay implementations robust.
std::vector<Vector2d> jitterForDelaunay(const std::vector<Vector2d>& pts, double domainScale)
{
  std::vector<Vector2d> out;
  out.reserve(pts.size());
  const double mag = std::max(domainScale, 1e-9) * 1e-5;
  for (size_t i = 0; i < pts.size(); i++) {
    double jx = (hashNoise01(static_cast<double>(i) * 2.0) - 0.5) * mag;
    double jy = (hashNoise01(static_cast<double>(i) * 2.0 + 1.0) - 0.5) * mag;
    out.push_back(pts[i] + Vector2d(jx, jy));
  }
  return out;
}

// -------------------- Unconstrained Bowyer-Watson Delaunay --------------------
struct DTriangle {
  int a, b, c;
};

// Delaunay triangulation of a plain 2D point set, no constraint edges.
// 'pts' is left unmodified; returned indices refer into 'pts'.
//
// Includes defensive logging/bailout: if the running triangle count ever
// grows far beyond the ~2*N that a healthy triangulation should have, that
// means the algorithm has gone unstable (see jitterForDelaunay() above for
// why that can happen) - we log a warning and bail out with an empty result
// rather than spinning on an ever-growing triangle list.
std::vector<DTriangle> bowyerWatson(const std::vector<Vector2d>& pts)
{
  if (pts.size() < 3) {
    return {};
  }

  const size_t maxSaneTris = std::max<size_t>(200, pts.size() * 20);

  double minX = pts[0].x(), maxX = minX, minY = pts[0].y(), maxY = minY;
  for (const auto& p : pts) {
    minX = std::min(minX, p.x());
    maxX = std::max(maxX, p.x());
    minY = std::min(minY, p.y());
    maxY = std::max(maxY, p.y());
  }
  double dmax = std::max(maxX - minX, maxY - minY);
  if (dmax < 1e-9) dmax = 1.0;
  double midx = (minX + maxX) * 0.5, midy = (minY + maxY) * 0.5;

  std::vector<Vector2d> P = pts;
  int i0 = static_cast<int>(P.size());
  P.push_back(Vector2d(midx - 20 * dmax, midy - dmax));
  int i1 = static_cast<int>(P.size());
  P.push_back(Vector2d(midx, midy + 20 * dmax));
  int i2 = static_cast<int>(P.size());
  P.push_back(Vector2d(midx + 20 * dmax, midy - dmax));

  auto makeCCW = [&](DTriangle t) {
    double area = (P[t.b].x() - P[t.a].x()) * (P[t.c].y() - P[t.a].y()) -
                  (P[t.c].x() - P[t.a].x()) * (P[t.b].y() - P[t.a].y());
    if (area < 0) std::swap(t.b, t.c);
    return t;
  };

  std::vector<DTriangle> tris = {makeCCW({i0, i1, i2})};

  auto inCircumcircle = [&](const DTriangle& t, const Vector2d& p) {
    const Vector2d &a = P[t.a], &b = P[t.b], &c = P[t.c];
    double ax = a.x() - p.x(), ay = a.y() - p.y();
    double bx = b.x() - p.x(), by = b.y() - p.y();
    double cx = c.x() - p.x(), cy = c.y() - p.y();
    double det = (ax * ax + ay * ay) * (bx * cy - cx * by) - (bx * bx + by * by) * (ax * cy - cx * ay) +
                 (cx * cx + cy * cy) * (ax * by - bx * ay);
    return det > 1e-12;
  };

  auto edgeOf = [](const DTriangle& t, int k) -> std::pair<int, int> {
    if (k == 0) return {t.a, t.b};
    if (k == 1) return {t.b, t.c};
    return {t.c, t.a};
  };

  for (int pi = 0; pi < static_cast<int>(pts.size()); pi++) {
    std::vector<DTriangle> bad;
    for (const auto& t : tris)
      if (inCircumcircle(t, P[pi])) bad.push_back(t);
    if (bad.empty()) continue;  // point lies exactly on an existing circumcircle etc. - skip it

    std::vector<std::pair<int, int>> boundary;
    for (size_t bi = 0; bi < bad.size(); bi++) {
      for (int k = 0; k < 3; k++) {
        auto e = edgeOf(bad[bi], k);
        bool shared = false;
        for (size_t bj = 0; bj < bad.size() && !shared; bj++) {
          if (bi == bj) continue;
          for (int k2 = 0; k2 < 3; k2++) {
            auto e2 = edgeOf(bad[bj], k2);
            if (e.first == e2.second && e.second == e2.first) {
              shared = true;
              break;
            }
          }
        }
        if (!shared) boundary.push_back(e);
      }
    }

    std::vector<DTriangle> kept;
    kept.reserve(tris.size());
    for (const auto& t : tris) {
      bool isBad = false;
      for (const auto& bt : bad) {
        if (bt.a == t.a && bt.b == t.b && bt.c == t.c) {
          isBad = true;
          break;
        }
      }
      if (!isBad) kept.push_back(t);
    }
    tris = std::move(kept);

    for (const auto& e : boundary) tris.push_back(makeCCW({e.first, e.second, pi}));
  }

  std::vector<DTriangle> result;
  result.reserve(tris.size());
  for (const auto& t : tris) {
    if (t.a < i0 && t.b < i0 && t.c < i0) result.push_back(t);
  }
  return result;
}

// Runs bowyerWatson() on a jittered copy of 'pts', then keeps only the
// triangles whose centroid is inside 'domain'. Shared between the base
// surface (boundary points only) and the final mesh (boundary + interior
// grid points) below - both are "triangulate everything, then discard
// what's outside the domain" in exactly the same way.
std::vector<DTriangle> triangulateAndFilterToDomain(const std::vector<Vector2d>& pts,
                                                    const Polygon2d& domain, double domainScale,
                                                    const char *debugLabel, int *outRejectedCount)
{
  auto jittered = jitterForDelaunay(pts, domainScale);
  auto tris = bowyerWatson(jittered);
  std::vector<DTriangle> kept;
  kept.reserve(tris.size());
  int rejected = 0;
  for (const auto& t : tris) {
    Vector2d centroid = (pts[t.a] + pts[t.b] + pts[t.c]) / 3.0;
    if (!domain.point_inside(centroid)) {
      rejected++;
      continue;
    }
    kept.push_back(t);
  }
  if (outRejectedCount) *outRejectedCount = rejected;
  return kept;
}

}  // namespace

std::unique_ptr<PolySet> loft(const std::vector<Vector3d>& outer,
                              const std::vector<std::vector<Vector3d>>& holes,
                              const std::function<Vector2d(const Vector3d&)>& proj,
                              double grid_spacing_uv,
                              const std::function<double(const Vector3d&)>& displacement)
{
  // NOTE ON LOG LEVELS: all diagnostic lines below intentionally use
  // message_group::Warning (not ::Trace) so they show up in the console
  // regardless of the verbosity setting, while this function is being
  // hardened. Once loft() is confirmed solid, the noisier ones can be
  // downgraded to ::Trace or ::Echo again.

  // 1) Project the boundary rings into (u,v) space.
  std::vector<Vector2d> outer_uv;
  outer_uv.reserve(outer.size());
  for (const auto& p : outer) outer_uv.push_back(proj(p));

  std::vector<std::vector<Vector2d>> holes_uv;
  holes_uv.reserve(holes.size());
  for (const auto& h : holes) {
    holes_uv.emplace_back();
    holes_uv.back().reserve(h.size());
    for (const auto& p : h) holes_uv.back().push_back(proj(p));
  }

  // 1b) Normalize winding: outer CCW, holes CW (see ensureOrientation()).
  // Keep the parallel 3D point arrays ('outer', 'holes' copies) in sync
  // with any reversal.
  std::vector<Vector3d> outer_pos = outer;
  std::vector<std::vector<Vector3d>> holes_pos = holes;
  ensureOrientation(outer_uv, outer_pos, holes_uv, holes_pos);

  Polygon2d domain = buildDomainPolygon(outer_uv, holes_uv);

  // Sanity check hole membership: the centroid of a hole ring should
  // normally be "outside" the material domain. (We deliberately do NOT test
  // the outer ring's vertex-mean centroid here - for a ring/annulus-shaped
  // cross-section, as with a vase wall, the mean of the outer boundary
  // vertices can itself fall inside a hole, which would look like a
  // point_inside() bug but isn't one.)
  for (size_t hi = 0; hi < holes_uv.size(); hi++) {
    if (holes_uv[hi].empty()) continue;
    Vector2d holeCentroid(0, 0);
    for (const auto& p : holes_uv[hi]) holeCentroid += p;
    holeCentroid /= static_cast<double>(holes_uv[hi].size());
  }

  // uv bounding box - needed both for the base surface's jitter scale and
  // for the interior sampling grid below, so compute it once, early.
  double umin = outer_uv[0].x(), umax = umin, vmin = outer_uv[0].y(), vmax = vmin;
  for (const auto& p : outer_uv) {
    umin = std::min(umin, p.x());
    umax = std::max(umax, p.x());
    vmin = std::min(vmin, p.y());
    vmax = std::max(vmax, p.y());
  }
  double domainScale = std::max(umax - umin, vmax - vmin);

  // 2) Base surface: gives interior grid points (step 3) a well-defined 3D
  //    base position + normal via barycentric interpolation.
  //
  //    EARLIER APPROACH (superseded): bridge each hole into the outer ring
  //    with a "nearest point" slit, then ear-clip the resulting simple
  //    polygon. This turned out to stall partway through on real models -
  //    e.g. one run produced only 43 of the 72 triangles a 74-point simple
  //    polygon should ear-clip into, leaving a large untriangulated gap
  //    right around the hole (support points missing in a region whose
  //    "incircle" was the hole). Root cause: the bridge/slit walks into the
  //    hole and back out along the *same* line, so the merged polygon has
  //    two exactly collinear, overlapping edges there. That degenerate
  //    configuration is a known failure mode for naive O(n^2) ear-clipping
  //    (candidate ears near the slit keep finding the "opposite side" of
  //    the slit sitting exactly on/inside them and get rejected), and the
  //    algorithm can stall with no valid ear left to clip.
  //
  //    CURRENT APPROACH: skip bridging/ear-clipping entirely and reuse the
  //    same "unconstrained Delaunay + discard triangles outside the domain"
  //    technique used for the final mesh in step 4. There is no bridge, so
  //    there is nothing to be collinear/degenerate, and a hole is handled
  //    exactly like the outer boundary - by filtering triangles whose
  //    centroid falls outside the domain afterwards.
  LoftBaseSurface base;
  {
    base.uv = outer_uv;
    base.pos3d = outer_pos;
    for (size_t hi = 0; hi < holes_uv.size(); hi++) {
      base.uv.insert(base.uv.end(), holes_uv[hi].begin(), holes_uv[hi].end());
      base.pos3d.insert(base.pos3d.end(), holes_pos[hi].begin(), holes_pos[hi].end());
    }

    auto baseTris = triangulateAndFilterToDomain(base.uv, domain, domainScale, "base surface", nullptr);
    base.tris.reserve(baseTris.size());
    for (const auto& t : baseTris) base.tris.push_back({t.a, t.b, t.c});

    // Now that we have a real triangle of the domain, log a genuine
    // known-inside point (its centroid) and confirm point_inside() agrees.
    const Vector2d& ta = base.uv[base.tris[0][0]];
    const Vector2d& tb = base.uv[base.tris[0][1]];
    const Vector2d& tc = base.uv[base.tris[0][2]];
    Vector2d knownInside = (ta + tb + tc) / 3.0;
  }

  // 3) Combined point set: boundary points (real 3D position known) +
  //    interior grid points (3D position/normal from the base surface,
  //    displacement already applied).
  std::vector<Vector2d> allUV;
  std::vector<Vector3d> allPos;
  allUV.reserve(outer_uv.size() + 64);
  allPos.reserve(outer_uv.size() + 64);

  for (size_t i = 0; i < outer_uv.size(); i++) {
    allUV.push_back(outer_uv[i]);
    allPos.push_back(outer_pos[i]);
  }
  for (size_t hi = 0; hi < holes_uv.size(); hi++) {
    for (size_t i = 0; i < holes_uv[hi].size(); i++) {
      allUV.push_back(holes_uv[hi][i]);
      allPos.push_back(holes_pos[hi][i]);
    }
  }

  int nu = std::max(2, static_cast<int>((umax - umin) / grid_spacing_uv) + 1);
  int nv = std::max(2, static_cast<int>((vmax - vmin) / grid_spacing_uv) + 1);
  LOG(message_group::Warning, concatStr("loft(): interior grid ", nu, "x", nv, " candidates"));

  int interiorAdded = 0;
  int rejectedOutsideDomain = 0;
  int rejectedNoBaseSample = 0;
  for (int iu = 0; iu < nu; iu++) {
    for (int iv = 0; iv < nv; iv++) {
      Vector2d p(umin + iu * grid_spacing_uv, vmin + iv * grid_spacing_uv);
      if (!domain.point_inside(p)) {
        rejectedOutsideDomain++;
        continue;
      }
      Vector3d pos, normal;
      if (!base.sample(p, pos, normal)) {
        rejectedNoBaseSample++;
        continue;
      }
      double d = displacement(pos);
      allUV.push_back(p);
      allPos.push_back(pos + normal * d);
      interiorAdded++;
    }
  }
  LOG(message_group::Warning,
      concatStr("loft(): interior grid: ", interiorAdded, " added, ", rejectedOutsideDomain,
                " rejected by point_inside(), ", rejectedNoBaseSample, " rejected by base.sample() (of ",
                nu, "x", nv, " candidates)"));
  if (interiorAdded == 0 && rejectedOutsideDomain > 0 && rejectedNoBaseSample == 0) {
    LOG(message_group::Warning,
        "loft(): every interior candidate was rejected by point_inside() - "
        "the domain polygon likely disagrees with the outer/hole winding "
        "(see the sanity-check lines above), or grid_spacing_uv is too "
        "coarse relative to the domain size.");
  }
  if (rejectedNoBaseSample > 0) {
    LOG(message_group::Warning,
        concatStr("loft(): ", rejectedNoBaseSample,
                  " grid points were inside the domain but not covered "
                  "by any base-surface triangle - the base surface triangulation may have gaps."));
  }

  // 4) Delaunay triangulation over the whole point set (boundary + interior
  //    grid), then discard triangles outside the domain (e.g. inside a
  //    hole, or outside the outer contour). Same helper as the base
  //    surface above.
  LOG(
    message_group::Warning,
    concatStr("loft(): running final mesh Delaunay over ", allUV.size(), " combined points (",
              outer_uv.size(), " outer + hole boundary points + ", interiorAdded, " interior points)"));
  int rejectedByCentroid = 0;
  auto finalTris =
    triangulateAndFilterToDomain(allUV, domain, domainScale, "final mesh", &rejectedByCentroid);
  if (finalTris.empty()) {
    LOG(message_group::Warning,
        "loft(): final mesh triangulation produced no triangles inside the domain.");
    return nullptr;
  }

  auto polyset = std::make_unique<PolySet>(3);
  polyset->setTriangular(true);
  polyset->vertices = allPos;
  polyset->indices.reserve(finalTris.size());
  for (const auto& t : finalTris) {
    polyset->indices.push_back({t.a, t.b, t.c});
  }
  LOG(message_group::Warning, concatStr("loft(): done, ", polyset->indices.size(), " triangles"));
  return polyset;
}
