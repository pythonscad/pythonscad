#include "geometry/loft.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

#include <Eigen/Eigenvalues>

#include "geometry/Polygon2d.h"
#include "utils/printutils.h"

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

// Reverses a ring's 2D points, its parallel 3D points, and (if present) its
// parallel tangent vectors together, so the vertex-by-vertex correspondence
// between uv[i], pos3d[i] and normal[i] is preserved.
void reverseRingInPlace(std::vector<Vector2d>& uv, std::vector<Vector3d>& pos3d,
                        std::vector<Vector3d>& normal)
{
  std::reverse(uv.begin(), uv.end());
  std::reverse(pos3d.begin(), pos3d.end());
  if (!normal.empty()) std::reverse(normal.begin(), normal.end());
}

// Forces 'outer' to be wound CCW and every entry of 'holes' to be wound CW
// (the standard convention for a polygon-with-holes), keeping the parallel
// 3D point arrays (and, when present, tangent arrays) in sync with any
// reversal.
void ensureOrientation(std::vector<Vector2d>& outer_uv, std::vector<Vector3d>& outer_pos,
                       std::vector<Vector3d>& outer_normal, std::vector<std::vector<Vector2d>>& holes_uv,
                       std::vector<std::vector<Vector3d>>& holes_pos,
                       std::vector<std::vector<Vector3d>>& holes_normal)
{
  if (signedArea2D(outer_uv) < 0) reverseRingInPlace(outer_uv, outer_pos, outer_normal);
  for (size_t hi = 0; hi < holes_uv.size(); hi++) {
    if (holes_uv[hi].size() < 3) continue;
    if (signedArea2D(holes_uv[hi]) > 0) {
      static std::vector<Vector3d> emptyNormal;
      std::vector<Vector3d>& hn = (hi < holes_normal.size()) ? holes_normal[hi] : emptyNormal;
      reverseRingInPlace(holes_uv[hi], holes_pos[hi], hn);
    }
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

// -------------------- Cubic (Hermite-style) triangle patch --------------------
//
// Boundary points (outer ring + hole rings) keep their EXACT input 3D
// position - this never changes, whether or not tangents are supplied (see
// the "boundary points keep their exact input position" comment further
// down). What tangents change is only how the INTERIOR of each boundary
// triangle is curved between those fixed corner points.
//
// This is deliberately NOT classic PN-Triangles: PN-Triangles take N as a
// *surface normal* and project the edge vector onto the tangent plane
// (i.e. they DISCARD the component of (Pj-Pi) along N). Here, N is the
// *tangent/leaving direction* of the surface at a boundary ring (e.g. the
// axis a circular port's wall leaves that port along) - exactly the
// direction PN-Triangles would throw away. So instead we build ordinary
// cubic Hermite control points directly from that tangent.
//
// A ring-id per corner point tells us whether an edge runs between two
// points of the SAME ring (in which case the tangent must NOT bend it -
// use the original straight/linear control points) or between two
// DIFFERENT rings (outer <-> hole), in which case the two tangents pull
// the edge's near-corner control points forward/backward respectively.
//
// The ring-id ordering convention (smaller ring id = "outer side", larger
// = "hole side") together with a caller-side global sign-normalization
// pass (see normalizeTangentSigns() below) is what keeps "leaving A" vs
// "arriving at B" consistent, independent of the arbitrary vertex order
// Bowyer-Watson happens to store a triangle's corners in.

void edgeControlPoints(const Vector3d& Pi, const Vector3d& Pj, const Vector3d& Ti, const Vector3d& Tj,
                       int ringI, int ringJ, Vector3d& outNearI, Vector3d& outNearJ)
{
  if (ringI == ringJ) {
    // Same ring: keep the original straight/linear control points so the
    // ring's own shape is never bent by a tangent meant for a *different*
    // boundary.
    outNearI = Pi + (Pj - Pi) / 3.0;
    outNearJ = Pj + (Pi - Pj) / 3.0;
    return;
  }
  const double L = (Pj - Pi).norm();
  if (ringI < ringJ) {
    outNearI = Pi + Ti * (L / 3.0);
    outNearJ = Pj - Tj * (L / 3.0);
  } else {
    outNearI = Pi - Ti * (L / 3.0);
    outNearJ = Pj + Tj * (L / 3.0);
  }

  // -------- DEBUG: nur die ersten paar Kanten protokollieren --------
  // Zeigt die eigentliche Rohdatenlage, mit der evalCubicTriangle() rechnet:
  // Kantenlaenge, beide Tangenten und die daraus abgeleiteten
  // Kontrollpunkte. Wenn Ti/Tj hier (fast) auf der Verbindungsgeraden
  // Pj-Pi liegen, bleibt die Flaeche zwangslaeufig fast gerade - das ist
  // dann kein Bug in evalCubicTriangle(), sondern eine Tangente, die
  // schon nach der Vorzeichen-Normalisierung praktisch axial zeigt.
  static int debugEdgeCount = 0;
  if (debugEdgeCount < 10) {
    debugEdgeCount++;
    Vector3d chordDir = (L > 1e-12) ? (Pj - Pi) / L : Vector3d(0, 0, 0);
    LOG(message_group::Warning,
        "loft debug edge #%1$d: ring %2$d->%3$d  L=%4$f  chordDir=(%5$f,%6$f,%7$f)", debugEdgeCount,
        ringI, ringJ, L, chordDir.x(), chordDir.y(), chordDir.z());
    LOG(message_group::Warning, "loft debug edge #%1$d: Ti=(%2$f,%3$f,%4$f) dot(chordDir)=%5$f",
        debugEdgeCount, Ti.x(), Ti.y(), Ti.z(), Ti.dot(chordDir));
    LOG(message_group::Warning, "loft debug edge #%1$d: Tj=(%2$f,%3$f,%4$f) dot(chordDir)=%5$f",
        debugEdgeCount, Tj.x(), Tj.y(), Tj.z(), Tj.dot(chordDir));
  }
}

Vector3d evalCubicTriangle(const Vector3d& P1, const Vector3d& P2, const Vector3d& P3,
                           const Vector3d& N1, const Vector3d& N2, const Vector3d& N3, int ring1,
                           int ring2, int ring3, double u, double v, double w)
{
  const Vector3d b300 = P1, b030 = P2, b003 = P3;
  Vector3d b210, b120, b021, b012, b102, b201;

  edgeControlPoints(P1, P2, N1, N2, ring1, ring2, b210, b120);
  edgeControlPoints(P2, P3, N2, N3, ring2, ring3, b021, b012);
  edgeControlPoints(P3, P1, N3, N1, ring3, ring1, b102, b201);

  const Vector3d E = (b210 + b120 + b021 + b012 + b102 + b201) / 6.0;
  const Vector3d V = (P1 + P2 + P3) / 3.0;
  const Vector3d b111 = E + (E - V) * 0.5;

  const double u2 = u * u, v2 = v * v, w2 = w * w, u3 = u2 * u, v3 = v2 * v, w3 = w2 * w;
  return b300 * u3 + b030 * v3 + b003 * w3 + 3 * b210 * u2 * v + 3 * b120 * u * v2 + 3 * b021 * v2 * w +
         3 * b012 * v * w2 + 3 * b102 * w2 * u + 3 * b201 * w * u2 + 6 * b111 * u * v * w;
}

// -------------------- Base surface: barycentric/cubic base position + normal --------------------
// Gives any interior (u,v) query point a well-defined base 3D position and
// face normal, via a coarse triangulation of the boundary points. If
// per-point tangents were supplied, the interior is curved with
// evalCubicTriangle(); otherwise it falls back to the original flat
// barycentric blend, unchanged.
struct LoftBaseSurface {
  std::vector<Vector2d> uv;      // 2D coordinates (boundary points only: outer + holes)
  std::vector<Vector3d> pos3d;   // parallel 3D positions
  std::vector<Vector3d> normal;  // parallel tangent/leaving directions; may be all-zero (unused)
  std::vector<int> ringId;       // parallel ring id: 0 = outer, i+1 = holes[i]
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

      const bool haveTangents =
        !normal.empty() && (normal[t[0]].squaredNorm() > 1e-12 || normal[t[1]].squaredNorm() > 1e-12 ||
                            normal[t[2]].squaredNorm() > 1e-12);

      if (!haveTangents) {
        outPos = w0 * p0 + w1 * p1 + w2 * p2;
      } else {
        outPos = evalCubicTriangle(p0, p1, p2, normal[t[0]], normal[t[1]], normal[t[2]], ringId[t[0]],
                                   ringId[t[1]], ringId[t[2]], w0, w1, w2);
      }

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

// Centroid of a list of 3D points (used only to establish a consistent
// global sign convention for tangents - see normalizeTangentSigns()).
Vector3d centroid3D(const std::vector<Vector3d>& pts)
{
  Vector3d c(0, 0, 0);
  for (const auto& p : pts) c += p;
  return pts.empty() ? c : c / static_cast<double>(pts.size());
}

// Makes sure every supplied tangent points, on balance, from the outer ring
// towards the (first) hole ring - i.e. all point "outer -> hole", not some
// forward and some backward. Without this, a tangent supplied by the
// caller in an arbitrary orientation (e.g. the plane normal of a 2D shape,
// which has no inherent "forward" vs "backward") can flip the curvature
// direction unpredictably from one boundary point to the next, producing
// self-intersecting surfaces instead of a simple bend.
//
// This is a coarse, whole-ring heuristic (single flip per ring based on the
// ring's average tangent vs. the outer->hole axis), not a per-point fix -
// it assumes each ring's tangent field is already roughly consistent with
// itself (true for a tangent derived from a single rigid 2D shape
// transform, as produced by python_loft_ring_from_shape()).
void normalizeTangentSigns(const std::vector<Vector3d>& outer_pos, std::vector<Vector3d>& outer_normal,
                           const std::vector<std::vector<Vector3d>>& holes_pos,
                           std::vector<std::vector<Vector3d>>& holes_normal)
{
  if (holes_pos.empty()) return;

  const Vector3d outerC = centroid3D(outer_pos);
  const Vector3d holeC = centroid3D(holes_pos[0]);
  Vector3d axis = holeC - outerC;
  if (axis.squaredNorm() < 1e-12) return;
  axis.normalize();

  auto flipRingIfNeeded = [&](std::vector<Vector3d>& ringNormal, const char *label) {
    if (ringNormal.empty()) return;
    Vector3d avg(0, 0, 0);
    for (const auto& n : ringNormal) avg += n;
    Vector3d avgBefore = avg;
    bool flipped = avg.dot(axis) < 0;
    if (flipped) {
      for (auto& n : ringNormal) n = -n;
    }
    // -------- DEBUG --------
    // 'dot' hier nahe 0 bedeutet: die Tangente steht (fast) senkrecht zur
    // outer->hole-Achse - dann ist die Vorzeichenwahl praktisch eine
    // Muenzwurf-Entscheidung und KEIN verlaesslicher Indikator, auf welcher
    // Seite die Woelbung landet. Ein 'flipped=1' bei einer Konfiguration,
    // die eigentlich keine ~180 Grad-Wende braucht, ist ein Hinweis, dass
    // hier die axis-Heuristik die falsche Seite waehlt.
    LOG(message_group::Warning,
        "loft debug normalizeTangentSigns[%1$s]: avgTangent=(%2$f,%3$f,%4$f) axis=(%5$f,%6$f,%7$f) "
        "dot=%8$f flipped=%9$d",
        label, avgBefore.x(), avgBefore.y(), avgBefore.z(), axis.x(), axis.y(), axis.z(),
        avgBefore.dot(axis), flipped ? 1 : 0);
  };

  flipRingIfNeeded(outer_normal, "outer");
  for (size_t hi = 0; hi < holes_normal.size(); hi++) {
    flipRingIfNeeded(holes_normal[hi], "hole");
  }
}

// -------------------- Automatic projection --------------------
//
// Builds a proj() perpendicular to the given axis (unit vector): u,v are
// the point's coordinates in the plane orthogonal to 'axis'. Same
// construction as the user's earlier hand-written make_axial_proj(), just
// in C++ and derived automatically instead of from two named points.
std::function<Vector2d(const Vector3d&)> makeAxialProjection(Vector3d axis)
{
  axis.normalize();
  Vector3d ref = (std::fabs(axis.y()) < 0.9) ? Vector3d(0, 1, 0) : Vector3d(1, 0, 0);
  Vector3d u_axis = axis.cross(ref).normalized();
  Vector3d v_axis = axis.cross(u_axis);  // already unit: axis, u_axis orthonormal
  return
    [u_axis, v_axis](const Vector3d& p) -> Vector2d { return Vector2d(p.dot(u_axis), p.dot(v_axis)); };
}

// Fits the best plane through 'pts' (least-squares, via the covariance
// matrix's eigenvectors) and returns a proj() onto it: u,v are the
// coordinates along the two directions of largest spread; the direction of
// smallest spread (the fitted plane's own normal) is dropped. Correct
// choice when outer and holes are all roughly coplanar (a flat panel with
// cutouts) - looking straight through the panel's own thickness.
std::function<Vector2d(const Vector3d&)> makePlaneProjection(const std::vector<Vector3d>& pts)
{
  Vector3d c = centroid3D(pts);
  Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
  for (const auto& p : pts) {
    Vector3d d = p - c;
    cov += d * d.transpose();
  }
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(cov);
  // Eigen returns eigenvalues/-vectors sorted ASCENDING: col(0) = smallest
  // variance (the panel's own normal, to be dropped), col(2)/col(1) = the
  // two largest-variance (in-plane) directions we keep as u/v.
  Vector3d u_axis = solver.eigenvectors().col(2);
  Vector3d v_axis = solver.eigenvectors().col(1);
  return [c, u_axis, v_axis](const Vector3d& p) -> Vector2d {
    Vector3d d = p - c;
    return Vector2d(d.dot(u_axis), d.dot(v_axis));
  };
}

// -------------------- Degenerate-hole detection & local re-embedding --------------------
//
// A hole ring's projected shape can become degenerate (near-zero extent in
// one direction) if its own boundary curve's natural spread doesn't align
// with whatever directions the *shared* projection happens to keep - e.g.
// a small hole drilled radially through a tube wall, when the tube itself
// is (correctly, for the tube's own sake) projected top-down, dropping
// exactly the axial coordinate the hole's own extent depends on. This is
// not fixable by picking a "smarter" shared projection: a rotationally
// symmetric outer/tube-partner ring and an arbitrarily oriented hole
// fundamentally need different 2D treatment - see the discussion that led
// here.
//
// The fix: only the TOPOLOGY actually matters for the Delaunay+domain-
// filter approach below - which points bound which loop, and that loops
// don't overlap - not that every ring share one single global projection
// formula. So each hole is projected with the shared 'effectiveProj' as
// usual, and only if THAT comes out degenerate, it is instead flattened
// with its OWN local best-fit plane (same PCA technique as
// makePlaneProjection(), just scoped to this one hole's own points) and
// then placed at the position its centroid *would* have under the shared
// projection - a single point is never degenerate, so this placement is
// always well-defined even when the hole's full extent isn't.

// True if 'uv's spread is too "flat" in one direction to form a sane
// simple polygon: computes the 2x2 covariance matrix's eigenvalues
// directly (no need to pull in Eigen for a 2x2) and checks their ratio.
bool isDegenerate2D(const std::vector<Vector2d>& uv)
{
  if (uv.size() < 3) return true;
  Vector2d c(0, 0);
  for (const auto& p : uv) c += p;
  c /= static_cast<double>(uv.size());

  double sxx = 0, syy = 0, sxy = 0;
  for (const auto& p : uv) {
    Vector2d d = p - c;
    sxx += d.x() * d.x();
    syy += d.y() * d.y();
    sxy += d.x() * d.y();
  }
  const double tr = sxx + syy;
  const double det = sxx * syy - sxy * sxy;
  const double disc = std::sqrt(std::max(0.0, tr * tr - 4 * det));
  const double lambdaMax = (tr + disc) * 0.5;
  const double lambdaMin = (tr - disc) * 0.5;
  if (lambdaMax < 1e-18) return true;
  const double aspectThreshold = 0.02;  // smaller/larger spread ratio
  return (lambdaMin / lambdaMax) < aspectThreshold;
}

// Locally flattens 'pts' (a single hole ring's own 3D boundary points) via
// its own best-fit plane, then places the result so its centroid sits at
// 'placeAt' - wherever the *shared* projection maps this hole's centroid
// to (see the long comment above). makePlaneProjection() already centers
// its output on the input points' own centroid, so 'placeAt' is simply
// added on top - no extra centroid bookkeeping needed here.
std::vector<Vector2d> localHoleReembed(const std::vector<Vector3d>& pts, const Vector2d& placeAt)
{
  auto localProj = makePlaneProjection(pts);
  std::vector<Vector2d> out;
  out.reserve(pts.size());
  for (const auto& p : pts) out.push_back(localProj(p) + placeAt);
  return out;
}

// Picks a proj() automatically when the caller doesn't supply one (see the
// long comment in loft.h). Decides between the "tube" (axial) and "panel"
// (best-fit-plane) case by comparing how far each hole's centroid sits
// from the outer ring's own centroid ('D') against the outer ring's own
// characteristic size ('R'): a hole clearly offset from the outer ring
// (D large relative to R) means the caller is connecting two separate
// rings in space (a tube/port), so the axial projection is used; holes
// close to the outer ring's own centroid (D small relative to R) mean
// they're most likely cutouts within the same flat panel, so the
// best-fit-plane projection is used instead.
//
// The axial direction itself is taken from the SINGLE FARTHEST hole only -
// not an average over all of them. A tube/port call typically has exactly
// one hole that's the genuine "other end" of the tube (large D), plus
// possibly several small, much closer holes (perforations/mounting holes
// through the wall, each individually handled by the degenerate-hole
// fallback above regardless of their own orientation). Averaging every
// hole's direction into the axis would let those close, off-axis holes tilt
// the axis away from the true tube direction - exactly the "correct near
// the holes, wrong on the far side" symptom this caused before switching to
// max-distance-only.
std::function<Vector2d(const Vector3d&)> computeAutoProj(const std::vector<Vector3d>& outer,
                                                         const std::vector<std::vector<Vector3d>>& holes)
{
  const Vector3d outerCentroid = centroid3D(outer);

  double R = 0.0;
  for (const auto& p : outer) R = std::max(R, (p - outerCentroid).norm());
  if (R < 1e-9) R = 1.0;

  double maxD = 0.0;
  Vector3d axis(0, 0, 0);
  for (const auto& h : holes) {
    if (h.empty()) continue;
    Vector3d holeCentroid = centroid3D(h);
    Vector3d toHole = holeCentroid - outerCentroid;
    double d = toHole.norm();
    if (d > maxD) {
      maxD = d;
      axis = (d > 1e-9) ? (toHole / d) : Vector3d(0, 0, 0);
    }
  }

  const bool tubeLike = !holes.empty() && (maxD > 0.5 * R) && axis.squaredNorm() > 1e-9;

  // -------- DEBUG --------
  LOG(message_group::Warning,
      "loft debug auto proj: outerR=%1$f maxHoleCentroidDist=%2$f axis=(%3$f,%4$f,%5$f) -> %6$s", R,
      maxD, axis.x(), axis.y(), axis.z(), tubeLike ? "axial (tube)" : "best-fit plane (panel)");

  if (tubeLike) {
    return makeAxialProjection(axis);
  }

  std::vector<Vector3d> allPts = outer;
  for (const auto& h : holes) allPts.insert(allPts.end(), h.begin(), h.end());
  return makePlaneProjection(allPts);
}

}  // namespace

std::unique_ptr<PolySet> loft(const std::vector<Vector3d>& outer,
                              const std::vector<std::vector<Vector3d>>& holes,
                              const std::function<Vector2d(const Vector3d&)>& proj,
                              double grid_spacing_uv,
                              const std::function<double(const Vector3d&)>& displacement,
                              const std::vector<Vector3d>& outer_normal,
                              const std::vector<std::vector<Vector3d>>& holes_normal)
{
  if (outer.size() < 3 || grid_spacing_uv <= 0.0) return nullptr;

  // 0) No proj() supplied -> pick one automatically (see loft.h / the long
  // comment on computeAutoProj() above). An explicitly supplied proj is
  // always used as-is and never overridden. 'autoProj' is a named local so
  // that, when used, 'effectiveProj' binds to a real object with function
  // scope lifetime rather than a dangling temporary.
  std::function<Vector2d(const Vector3d&)> autoProj;
  const std::function<Vector2d(const Vector3d&)>& effectiveProj =
    proj ? proj : (autoProj = computeAutoProj(outer, holes));

  // 1) Project the boundary rings into (u,v) space.
  std::vector<Vector2d> outer_uv;
  outer_uv.reserve(outer.size());
  for (const auto& p : outer) outer_uv.push_back(effectiveProj(p));

  std::vector<std::vector<Vector2d>> holes_uv;
  holes_uv.reserve(holes.size());
  for (const auto& h : holes) {
    holes_uv.emplace_back();
    holes_uv.back().reserve(h.size());
    for (const auto& p : h) holes_uv.back().push_back(effectiveProj(p));
  }

  // 1a) Degenerate-hole fallback: a hole whose shape under the SHARED
  // projection is too "flat" (see isDegenerate2D()) is re-flattened with
  // its own local best-fit plane and placed at its centroid's position
  // under the shared projection instead - see the long comment above
  // localHoleReembed(). Holes that already project fine (the normal case,
  // e.g. a genuine tube-partner ring) are left completely untouched.
  for (size_t hi = 0; hi < holes_uv.size(); hi++) {
    if (holes[hi].size() < 3 || !isDegenerate2D(holes_uv[hi])) continue;
    Vector2d placeAt = effectiveProj(centroid3D(holes[hi]));
    holes_uv[hi] = localHoleReembed(holes[hi], placeAt);
    LOG(message_group::Warning,
        "loft debug hole reembed: holes[%1$d] was degenerate under the shared projection - "
        "flattened locally, placed at (%2$f,%3$f)",
        (int)hi, placeAt.x(), placeAt.y());
  }

  // -------- DEBUG --------
  // Rohe UV-Werte DIREKT nach proj(), noch vor jeder Umorientierung/
  // Normalisierung. Wenn diese Zahlen sich zwischen zwei loft()-Aufrufen
  // im selben Skript (gleiche 3D-Punkte, mathematisch fast identische
  // proj()) unterscheiden, liefert proj() selbst unterschiedliche Werte -
  // z.B. weil eine falsche/alte proj()-Closure aufgerufen wird (Python-
  // Bindung/Cache-Verwechslung). Sind sie identisch, liegt die
  // Abweichung NICHT an proj(), sondern irgendwo danach (Cache-Kollision
  // via LoftNode::toString(), siehe Kommentar dort).
  for (size_t i = 0; i < std::min<size_t>(3, outer_uv.size()); i++) {
    LOG(message_group::Warning,
        "loft debug raw uv: outer_uv[%1$d]=(%2$f,%3$f) from outer[%4$d]=(%5$f,%6$f,%7$f)", (int)i,
        outer_uv[i].x(), outer_uv[i].y(), (int)i, outer[i].x(), outer[i].y(), outer[i].z());
  }
  for (size_t hi = 0; hi < holes_uv.size(); hi++) {
    for (size_t i = 0; i < std::min<size_t>(3, holes_uv[hi].size()); i++) {
      LOG(message_group::Warning,
          "loft debug raw uv: holes_uv[%1$d][%2$d]=(%3$f,%4$f) from holes[%5$d][%6$d]=(%7$f,%8$f,%9$f)",
          (int)hi, (int)i, holes_uv[hi][i].x(), holes_uv[hi][i].y(), (int)hi, (int)i, holes[hi][i].x(),
          holes[hi][i].y(), holes[hi][i].z());
    }
  }

  // 1b) Normalize winding: outer CCW, holes CW (see ensureOrientation()).
  // Keep the parallel 3D point arrays ('outer_pos', 'holes_pos' copies) and
  // the (optional) tangent arrays in sync with any reversal.
  std::vector<Vector3d> outer_pos = outer;
  std::vector<std::vector<Vector3d>> holes_pos = holes;
  std::vector<Vector3d> outer_normal_copy = outer_normal;
  std::vector<std::vector<Vector3d>> holes_normal_copy = holes_normal;
  ensureOrientation(outer_uv, outer_pos, outer_normal_copy, holes_uv, holes_pos, holes_normal_copy);

  // 1c) Make sure all supplied tangents point "outer -> hole" on balance,
  // so the cubic patch below bends consistently instead of flipping
  // direction from one boundary point to the next. No-op when no tangents
  // were supplied (both arrays stay empty).
  normalizeTangentSigns(outer_pos, outer_normal_copy, holes_pos, holes_normal_copy);

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
  //    base position + normal via barycentric (or, when tangents were
  //    supplied, cubic) interpolation. Built by triangulating just the
  //    boundary points (outer + holes) with the same unconstrained-
  //    Delaunay-plus-domain-filter technique used for the final mesh in
  //    step 4, then discarding triangles outside the domain.
  LoftBaseSurface base;
  {
    base.uv = outer_uv;
    base.pos3d = outer_pos;
    base.normal = outer_normal_copy.empty() ? std::vector<Vector3d>(outer_pos.size(), Vector3d(0, 0, 0))
                                            : outer_normal_copy;
    base.ringId = std::vector<int>(outer_pos.size(), 0);
    for (size_t hi = 0; hi < holes_uv.size(); hi++) {
      base.uv.insert(base.uv.end(), holes_uv[hi].begin(), holes_uv[hi].end());
      base.pos3d.insert(base.pos3d.end(), holes_pos[hi].begin(), holes_pos[hi].end());
      const bool haveHoleNormal = hi < holes_normal_copy.size() && !holes_normal_copy[hi].empty();
      const std::vector<Vector3d> zeros(holes_pos[hi].size(), Vector3d(0, 0, 0));
      const std::vector<Vector3d>& hn = haveHoleNormal ? holes_normal_copy[hi] : zeros;
      base.normal.insert(base.normal.end(), hn.begin(), hn.end());
      base.ringId.insert(base.ringId.end(), holes_pos[hi].size(), static_cast<int>(hi) + 1);
    }

    auto baseTris = triangulateAndFilterToDomain(base.uv, domain, domainScale);
    base.tris.reserve(baseTris.size());
    for (const auto& t : baseTris) base.tris.push_back({t.a, t.b, t.c});
    if (base.tris.empty()) return nullptr;

    // -------- DEBUG --------
    // "cross" = Dreiecke, die outer und mindestens ein Loch verbinden -
    // NUR diese werden von evalCubicTriangle() gekruemmt. Ist crossRing
    // hier 0 (oder sehr klein), gibt es schlicht keine Flaeche, an der
    // ueberhaupt eine Woelbung sichtbar werden koennte - dann liegt das
    // Problem VOR der Kruemmungsberechnung, in der UV-Domain/Projektion.
    int crossRing = 0, sameRing = 0;
    for (const auto& t : base.tris) {
      int r0 = base.ringId[t[0]], r1 = base.ringId[t[1]], r2 = base.ringId[t[2]];
      if (r0 == r1 && r1 == r2) sameRing++;
      else crossRing++;
    }
    LOG(message_group::Warning,
        "loft debug base surface: outer_pts=%1$d hole_pts_total=%2$d tris=%3$d crossRing=%4$d "
        "sameRing=%5$d",
        (int)outer_uv.size(), (int)(base.uv.size() - outer_uv.size()), (int)base.tris.size(), crossRing,
        sameRing);
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
  // exactly. The same reasoning is why boundary points are exempt from
  // curvature above (they are stored raw in base.pos3d and only ever used
  // as fixed corner points of the cubic patch, never re-evaluated by it).
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
  int debugGridCandidates = 0, debugGridOutsideDomain = 0, debugGridEncroached = 0,
      debugGridSampleFailed = 0, debugGridAccepted = 0;

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
      debugGridCandidates++;
      Vector2d p(umin + iu * grid_spacing_uv, vmin + iv * grid_spacing_uv);
      if (!domain.point_inside(p)) {
        debugGridOutsideDomain++;
        continue;
      }

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
      if (pointEncroachesBoundary(p, outer_uv, holes_uv)) {
        debugGridEncroached++;
        continue;
      }

      Vector3d pos, normal;
      if (!base.sample(p, pos, normal)) {
        debugGridSampleFailed++;
        continue;
      }

      double falloff = smoothstep01(distanceToBoundary(p, outer_uv, holes_uv) / blendDist);
      double d = displacement(pos) * falloff;
      allUV.push_back(p);
      allPos.push_back(pos + normal * d);
      debugGridAccepted++;
    }
  }

  // -------- DEBUG --------
  // Zeigt, wie viele der nu*nv Rasterpunkte tatsaechlich Flaeche wurden.
  // sampleFailed > 0 heisst: ein Punkt lag laut domain.point_inside() IN
  // der Domain, aber KEIN Dreieck von base.tris hat ihn eingeschlossen -
  // typischerweise ein Randfall zwischen zwei UV-Domains (Rundungsfehler
  // an der Grenze der coarsen Basistriangulierung, siehe LoftBaseSurface::sample()).
  // Viele "outsideDomain" bei einer erwartet fast-vollen Flaeche deutet
  // wieder auf das alte "Loch ausserhalb / proj() passt nicht"-Muster hin.
  LOG(message_group::Warning,
      "loft debug grid: nu=%1$d nv=%2$d candidates=%3$d outsideDomain=%4$d encroached=%5$d "
      "sampleFailed=%6$d accepted=%7$d",
      nu, nv, debugGridCandidates, debugGridOutsideDomain, debugGridEncroached, debugGridSampleFailed,
      debugGridAccepted);

  // 4) Delaunay triangulation over the whole point set (boundary + interior
  //    grid), then discard triangles outside the domain (e.g. inside a
  //    hole, or outside the outer contour). Same helper as the base
  //    surface above.
  auto finalTris = triangulateAndFilterToDomain(allUV, domain, domainScale);
  LOG(message_group::Warning, "loft debug final mesh: points=%1$d triangles=%2$d", (int)allPos.size(),
      (int)finalTris.size());
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
