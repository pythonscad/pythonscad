#include "geometry/loft.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

#include "geometry/Polygon2d.h"

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
// (the standard convention for a polygon-with-holes), keeping the parallel
// 3D point arrays in sync with any reversal.
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

// -------------------- Distance to the domain boundary --------------------
// Shortest distance from p to a single closed polyline (ring): the minimum
// distance to any of its edges.
double distanceToRing(const Vector2d& p, const std::vector<Vector2d>& ring)
{
  double best = std::numeric_limits<double>::max();
  const size_t n = ring.size();
  for (size_t i = 0; i < n; i++) {
    const Vector2d& a = ring[i];
    const Vector2d& b = ring[(i + 1) % n];
    Vector2d ab = b - a;
    double len2 = ab.squaredNorm();
    double t = (len2 > 1e-18) ? std::clamp((p - a).dot(ab) / len2, 0.0, 1.0) : 0.0;
    Vector2d closest = a + ab * t;
    best = std::min(best, (p - closest).norm());
  }
  return best;
}

// Shortest distance from p to the domain boundary, i.e. to the nearest
// point on the outer ring or on any hole ring.
double distanceToBoundary(const Vector2d& p, const std::vector<Vector2d>& outer_uv,
                          const std::vector<std::vector<Vector2d>>& holes_uv)
{
  double best = distanceToRing(p, outer_uv);
  for (const auto& h : holes_uv) best = std::min(best, distanceToRing(p, h));
  return best;
}

double smoothstep01(double x)
{
  x = std::clamp(x, 0.0, 1.0);
  return x * x * (3.0 - 2.0 * x);
}

// -------------------- Boundary-edge encroachment test --------------------
// A point p "encroaches" on a segment (a,b) if it lies inside that
// segment's diametral circle - the (smallest possible) circle that passes
// through a and b, centered at their midpoint. This is a standard,
// precisely justified test from constrained-Delaunay mesh generation
// (Ruppert/Chew): if a segment's diametral circle contains no other point,
// that segment is *guaranteed* to survive as an edge of the (unconstrained)
// Delaunay triangulation. Checking it reduces to a single sign: p sees a
// and b at an obtuse angle, i.e. (a-p).(b-p) < 0.
//
// This is used below to keep interior grid points from being placed close
// enough to any boundary-ring edge to threaten it - unlike a flat distance
// threshold, it scales automatically with each edge's actual local length,
// so it doesn't over-exclude where the ring is coarse and under-exclude
// where it's fine.
bool pointEncroachesEdge(const Vector2d& p, const Vector2d& a, const Vector2d& b)
{
  return (a - p).dot(b - p) < 0.0;
}

bool pointEncroachesBoundary(const Vector2d& p, const std::vector<Vector2d>& outer_uv,
                             const std::vector<std::vector<Vector2d>>& holes_uv)
{
  auto checkRing = [&](const std::vector<Vector2d>& ring) {
    const size_t n = ring.size();
    for (size_t i = 0; i < n; i++) {
      if (pointEncroachesEdge(p, ring[i], ring[(i + 1) % n])) return true;
    }
    return false;
  };
  if (checkRing(outer_uv)) return true;
  for (const auto& h : holes_uv) {
    if (checkRing(h)) return true;
  }
  return false;
}

// -------------------- Base surface: barycentric base position + normal --------------------
// Gives any interior (u,v) query point a well-defined base 3D position and
// face normal via barycentric interpolation over a coarse triangulation of
// the boundary points.
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
// membership tests, the final vertex positions) keeps using.
//
// WHY THIS IS NEEDED: an axis-aligned regular grid (which is exactly what
// the interior sampling below produces) is a classic pathological input for
// Bowyer-Watson - every 2x2 grid cell's four corners are exactly cocircular,
// which can corrupt the incremental algorithm's "bad-triangle cavity is a
// single closed loop" assumption. A small perturbation, well above floating
// point noise but far below anything visible in the final geometry, avoids
// this - a standard trick for making incremental Delaunay implementations
// robust against structured input.
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
std::vector<DTriangle> bowyerWatson(const std::vector<Vector2d>& pts)
{
  if (pts.size() < 3) return {};

  // Defensive bound: a healthy triangulation of N points has ~2N triangles.
  // If jitterForDelaunay() somehow still leaves a numerically degenerate
  // configuration, bail out rather than growing without limit.
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

    if (tris.size() > maxSaneTris) return {};
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
                                                    const Polygon2d& domain, double domainScale)
{
  auto jittered = jitterForDelaunay(pts, domainScale);
  auto tris = bowyerWatson(jittered);

  std::vector<DTriangle> kept;
  kept.reserve(tris.size());
  for (const auto& t : tris) {
    Vector2d centroid = (pts[t.a] + pts[t.b] + pts[t.c]) / 3.0;
    if (!domain.point_inside(centroid)) continue;
    kept.push_back(t);
  }
  return kept;
}

}  // namespace

std::unique_ptr<PolySet> loft(const std::vector<Vector3d>& outer,
                              const std::vector<std::vector<Vector3d>>& holes,
                              const std::function<Vector2d(const Vector3d&)>& proj,
                              double grid_spacing_uv,
                              const std::function<double(const Vector3d&)>& displacement)
{
  if (outer.size() < 3 || grid_spacing_uv <= 0.0) return nullptr;

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

  // uv bounding box - needed for the base surface's jitter scale and the
  // interior sampling grid below.
  double umin = outer_uv[0].x(), umax = umin, vmin = outer_uv[0].y(), vmax = vmin;
  for (const auto& p : outer_uv) {
    umin = std::min(umin, p.x());
    umax = std::max(umax, p.x());
    vmin = std::min(vmin, p.y());
    vmax = std::max(vmax, p.y());
  }
  double domainScale = std::max(umax - umin, vmax - vmin);

  // 2) Base surface: gives interior grid points (step 3) a well-defined 3D
  //    base position + normal via barycentric interpolation. Built by
  //    triangulating just the boundary points (outer + holes) with the same
  //    unconstrained-Delaunay-plus-domain-filter technique used for the
  //    final mesh in step 4, then discarding triangles outside the domain.
  LoftBaseSurface base;
  {
    base.uv = outer_uv;
    base.pos3d = outer_pos;
    for (size_t hi = 0; hi < holes_uv.size(); hi++) {
      base.uv.insert(base.uv.end(), holes_uv[hi].begin(), holes_uv[hi].end());
      base.pos3d.insert(base.pos3d.end(), holes_pos[hi].begin(), holes_pos[hi].end());
    }

    auto baseTris = triangulateAndFilterToDomain(base.uv, domain, domainScale);
    base.tris.reserve(baseTris.size());
    for (const auto& t : baseTris) base.tris.push_back({t.a, t.b, t.c});
    if (base.tris.empty()) return nullptr;
  }

  // 3) Combined point set: boundary points (real 3D position known, NEVER
  //    displaced - see below) + interior grid points (3D position/normal
  //    from the base surface, displacement applied with a smooth falloff
  //    towards the boundary).
  std::vector<Vector2d> allUV;
  std::vector<Vector3d> allPos;
  allUV.reserve(outer_uv.size() + 64);
  allPos.reserve(outer_uv.size() + 64);

  // Boundary points (outer ring + hole rings) keep their exact input
  // position. displacement() is intentionally never called for these - a
  // texture/bump function must not be able to move a point that lies on the
  // outer or inner (hole) contour, since those contours are what neighboring
  // loft() calls (e.g. the next ring up/down a vase wall) rely on to line up
  // exactly.
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

  // Interior grid points fade the displacement out smoothly as they
  // approach the boundary, reaching exactly zero at 'blendDist' or closer.
  // This isn't just cosmetic: without it, the row of interior points
  // immediately next to the boundary gets the full displacement while the
  // boundary itself (see above) gets none, producing a visible step right
  // at the rim - which looks like the boundary itself got bumped even
  // though, strictly, it never moved.
  const double blendDist = std::max(2.0 * grid_spacing_uv, 1e-9);

  for (int iu = 0; iu < nu; iu++) {
    for (int iv = 0; iv < nv; iv++) {
      Vector2d p(umin + iu * grid_spacing_uv, vmin + iv * grid_spacing_uv);
      if (!domain.point_inside(p)) continue;

      // Protects the *topology* of the final mesh, not just the
      // displacement: step 4's triangulation is an UNCONSTRAINED Delaunay
      // triangulation, which has no obligation to keep a boundary ring's
      // own short edges (e.g. between two adjacent points of 'outer' or a
      // hole) as edges of the result. If an interior point lies inside a
      // boundary edge's diametral circle, it can "steal" that edge (the
      // empty-circumcircle rule prefers connecting to it instead),
      // silently dropping the straight boundary edge and replacing it with
      // thinner triangles that bulge slightly off the true boundary line.
      // That breaks watertightness against any neighboring surface built
      // from the same ring (e.g. the next segment of a vase wall, or a cap
      // over the same hole). See pointEncroachesBoundary() above for the
      // exact (and tight - it scales with each edge's own local length,
      // not a fixed distance) test used here.
      if (pointEncroachesBoundary(p, outer_uv, holes_uv)) continue;

      Vector3d pos, normal;
      if (!base.sample(p, pos, normal)) continue;

      double falloff = smoothstep01(distanceToBoundary(p, outer_uv, holes_uv) / blendDist);
      double d = displacement(pos) * falloff;
      allUV.push_back(p);
      allPos.push_back(pos + normal * d);
    }
  }

  // 4) Delaunay triangulation over the whole point set (boundary + interior
  //    grid), then discard triangles outside the domain (e.g. inside a
  //    hole, or outside the outer contour). Same helper as the base
  //    surface above.
  auto finalTris = triangulateAndFilterToDomain(allUV, domain, domainScale);
  if (finalTris.empty()) return nullptr;

  auto polyset = std::make_unique<PolySet>(3);
  polyset->setTriangular(true);
  polyset->vertices = allPos;
  polyset->indices.reserve(finalTris.size());
  for (const auto& t : finalTris) {
    polyset->indices.push_back({t.a, t.b, t.c});
  }
  return polyset;
}
