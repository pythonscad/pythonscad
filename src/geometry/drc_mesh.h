#pragma once

// Mesh building blocks shared by the design rule checks (facing_check,
// slope_check): prepared triangles, a BVH with box / segment / ray queries,
// small convex polygon clipping, exact polygon distance, threading helpers
// and the cheap display offset of error solids.

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <memory>
#include <thread>
#include <utility>
#include <vector>

#include "geometry/PolySet.h"
#include "geometry/PolySetUtils.h"
#include "geometry/linalg.h"
#ifdef ENABLE_MANIFOLD
#include <manifold/manifold.h>
#endif

namespace DrcMesh {

constexpr double kInf = std::numeric_limits<double>::infinity();

inline double clamp01(double v)
{
  return v < 0 ? 0 : (v > 1 ? 1 : v);
}

// ------------------------------------------------------- prepared triangles

struct Tri {
  Vector3d v[3];  // CCW about n (already flipped for External)
  int vi[3];      // mesh vertex indices matching v[]
  Vector3d n;     // unit normal, pointing away from the checked region
  double off;     // n . v[0]
  BoundingBox box;
  bool valid;
};

inline BoundingBox expanded(const BoundingBox& b, double d)
{
  const Vector3d e(d, d, d);
  return BoundingBox(b.min() - e, b.max() + e);
}

// ----------------------------------------------------------------------- BVH

class Bvh
{
public:
  explicit Bvh(const std::vector<Tri>& tris) : tris_(tris)
  {
    idx_.reserve(tris.size());
    for (int i = 0; i < (int)tris.size(); ++i)
      if (tris[i].valid) idx_.push_back(i);
    if (!idx_.empty()) {
      nodes_.reserve(2 * idx_.size() / kLeaf + 2);
      nodes_.emplace_back();
      build(0, 0, (int)idx_.size());
    }
  }

  template <class F>
  void query(const BoundingBox& q, F&& f) const
  {
    if (nodes_.empty()) return;
    int stack[64];
    int sp = 0;
    stack[sp++] = 0;
    while (sp) {
      const Node& n = nodes_[stack[--sp]];
      if (!n.box.intersects(q)) continue;
      if (n.count) {
        for (int k = n.start; k < n.start + n.count; ++k)
          if (tris_[idx_[k]].box.intersects(q)) f(idx_[k]);
      } else {
        stack[sp++] = n.left;
        stack[sp++] = n.left + 1;
      }
    }
  }

  // Does the open segment a + t (b - a), t in (t0, t1), cross any triangle
  // other than skip0 / skip1?
  bool segmentHits(const Vector3d& a, const Vector3d& b, double t0, double t1, int skip0,
                   int skip1) const
  {
    if (nodes_.empty()) return false;
    const Vector3d d = b - a;
    const Vector3d inv = d.cwiseInverse();
    int stack[64];
    int sp = 0;
    stack[sp++] = 0;
    while (sp) {
      const Node& n = nodes_[stack[--sp]];
      if (!slab(n.box, a, inv, t0, t1)) continue;
      if (n.count) {
        for (int k = n.start; k < n.start + n.count; ++k) {
          const int ti = idx_[k];
          if (ti == skip0 || ti == skip1) continue;
          if (segTri(a, d, tris_[ti], t0, t1)) return true;
        }
      } else {
        stack[sp++] = n.left;
        stack[sp++] = n.left + 1;
      }
    }
    return false;
  }

  // Is p inside the (original, unflipped) solid? Majority vote of three ray
  // parities along skewed directions, robust against rays through edges.
  bool inside(const Vector3d& p) const
  {
    static const Vector3d dirs[3] = {
      {0.5773, 0.5774, 0.5775}, {-0.7071, 0.1234, 0.6963}, {0.2113, -0.9530, 0.2173}};
    int votes = 0;
    for (const Vector3d& d : dirs) votes += rayCrossings(p, d) & 1;
    return votes >= 2;
  }

private:
  static constexpr int kLeaf = 4;
  struct Node {
    BoundingBox box;
    int left = 0;  // children at left, left + 1
    int start = 0;
    int count = 0;  // > 0: leaf
  };

  void build(int slot, int start, int end)
  {
    BoundingBox b, cb;
    for (int k = start; k < end; ++k) {
      const Tri& t = tris_[idx_[k]];
      b.extend(t.box);
      cb.extend(t.box.center());
    }
    nodes_[slot].box = b;
    if (end - start <= kLeaf) {
      nodes_[slot].start = start;
      nodes_[slot].count = end - start;
      return;
    }
    int axis;
    cb.sizes().maxCoeff(&axis);
    const int mid = (start + end) / 2;
    std::nth_element(idx_.begin() + start, idx_.begin() + mid, idx_.begin() + end, [&](int a, int c) {
      return tris_[a].box.center()[axis] < tris_[c].box.center()[axis];
    });
    const int left = (int)nodes_.size();
    nodes_.emplace_back();
    nodes_.emplace_back();
    nodes_[slot].left = left;
    build(left, start, mid);
    build(left + 1, mid, end);
  }

  int rayCrossings(const Vector3d& o, const Vector3d& dir) const
  {
    const Vector3d inv = dir.cwiseInverse();
    int stack[64];
    int sp = 0, hits = 0;
    stack[sp++] = 0;
    while (sp) {
      const Node& n = nodes_[stack[--sp]];
      if (!slab(n.box, o, inv, 0, kInf)) continue;
      if (n.count) {
        for (int k = n.start; k < n.start + n.count; ++k)
          if (segTri(o, dir, tris_[idx_[k]], 0, kInf)) ++hits;
      } else {
        stack[sp++] = n.left;
        stack[sp++] = n.left + 1;
      }
    }
    return hits;
  }

  static bool slab(const BoundingBox& b, const Vector3d& o, const Vector3d& inv, double t0, double t1)
  {
    for (int k = 0; k < 3; ++k) {
      double ta = (b.min()[k] - o[k]) * inv[k], tb = (b.max()[k] - o[k]) * inv[k];
      if (std::isnan(ta) || std::isnan(tb)) {  // parallel and on a slab boundary
        if (o[k] < b.min()[k] || o[k] > b.max()[k]) return false;
        continue;
      }
      if (ta > tb) std::swap(ta, tb);
      t0 = std::max(t0, ta);
      t1 = std::min(t1, tb);
      if (t0 > t1) return false;
    }
    return true;
  }

  // Moller-Trumbore, hit parameter must lie in (t0, t1)
  static bool segTri(const Vector3d& o, const Vector3d& d, const Tri& t, double t0, double t1)
  {
    const Vector3d e1 = t.v[1] - t.v[0], e2 = t.v[2] - t.v[0];
    const Vector3d p = d.cross(e2);
    const double det = e1.dot(p);
    if (std::fabs(det) < 1e-300) return false;
    const double inv = 1.0 / det;
    const Vector3d s = o - t.v[0];
    const double u = s.dot(p) * inv;
    if (u < 0 || u > 1) return false;
    const Vector3d qv = s.cross(e1);
    const double v = d.dot(qv) * inv;
    if (v < 0 || u + v > 1) return false;
    const double tt = e2.dot(qv) * inv;
    return tt > t0 && tt < t1;
  }

  const std::vector<Tri>& tris_;
  std::vector<int> idx_;
  std::vector<Node> nodes_;
};

// ------------------------------------------------------------ small polygons

struct Poly {
  static constexpr int kMax = 24;
  Vector3d p[kMax];
  int n = 0;
};

// Keep the part with h . x + c >= -eps (Sutherland-Hodgman, one plane).
inline void clip(Poly& poly, const Vector3d& h, double c, double eps)
{
  if (poly.n == 0) return;
  double dv[Poly::kMax];
  bool all_in = true;
  for (int i = 0; i < poly.n; ++i) {
    dv[i] = h.dot(poly.p[i]) + c;
    if (dv[i] < -eps) all_in = false;
  }
  if (all_in) return;
  Poly out;
  for (int i = 0; i < poly.n; ++i) {
    const int j = (i + 1) % poly.n;
    const bool in_i = dv[i] >= -eps, in_j = dv[j] >= -eps;
    if (in_i && out.n < Poly::kMax) out.p[out.n++] = poly.p[i];
    if (in_i != in_j && out.n < Poly::kMax) {
      const double t = dv[i] / (dv[i] - dv[j]);
      out.p[out.n++] = poly.p[i] + (poly.p[j] - poly.p[i]) * t;
    }
  }
  poly = out;
}

// Closest points between segments p1-q1 and p2-q2 (Ericson, RTCD 5.1.9).
inline double segSeg(const Vector3d& p1, const Vector3d& q1, const Vector3d& p2, const Vector3d& q2,
                     Vector3d& c1, Vector3d& c2)
{
  const Vector3d d1 = q1 - p1, d2 = q2 - p2, r = p1 - p2;
  const double a = d1.dot(d1), e = d2.dot(d2), f = d2.dot(r);
  const double tiny = 1e-300;
  double s, t;
  if (a <= tiny && e <= tiny) {
    s = t = 0;
  } else if (a <= tiny) {
    s = 0;
    t = clamp01(f / e);
  } else {
    const double c = d1.dot(r);
    if (e <= tiny) {
      t = 0;
      s = clamp01(-c / a);
    } else {
      const double b = d1.dot(d2), denom = a * e - b * b;
      s = denom > 0 ? clamp01((b * f - c * e) / denom) : 0;
      t = (b * s + f) / e;
      if (t < 0) {
        t = 0;
        s = clamp01(-c / a);
      } else if (t > 1) {
        t = 1;
        s = clamp01((b - c) / a);
      }
    }
  }
  c1 = p1 + d1 * s;
  c2 = p2 + d2 * t;
  return (c1 - c2).squaredNorm();
}

// Distance of x to the interior of convex polygon P in the plane with normal
// n; edges are covered by the segment-segment tests.
inline bool pointFace(const Vector3d& x, const Poly& P, const Vector3d& n, double eps, double& d2,
                      Vector3d& foot)
{
  if (P.n < 3) return false;
  const double h = n.dot(x - P.p[0]);
  const Vector3d y = x - n * h;
  for (int i = 0; i < P.n; ++i) {
    const Vector3d& a = P.p[i];
    const Vector3d& b = P.p[(i + 1) % P.n];
    if (n.cross(b - a).dot(y - a) < -eps * (b - a).norm()) return false;
  }
  d2 = h * h;
  foot = y;
  return true;
}

// Exact distance between two convex planar polygons: vertex-face + edge-edge.
inline double polyDist(const Poly& A, const Vector3d& na, const Poly& B, const Vector3d& nb, double eps,
                       Vector3d& pa, Vector3d& pb)
{
  double best = kInf, d2;
  Vector3d c1, c2;
  for (int i = 0; i < A.n; ++i)
    if (pointFace(A.p[i], B, nb, eps, d2, c2) && d2 < best) {
      best = d2;
      pa = A.p[i];
      pb = c2;
    }
  for (int j = 0; j < B.n; ++j)
    if (pointFace(B.p[j], A, na, eps, d2, c1) && d2 < best) {
      best = d2;
      pa = c1;
      pb = B.p[j];
    }
  const int ea = A.n < 2 ? 1 : A.n, eb = B.n < 2 ? 1 : B.n;
  for (int i = 0; i < ea; ++i) {
    for (int j = 0; j < eb; ++j) {
      d2 = segSeg(A.p[i], A.p[(i + 1) % A.n], B.p[j], B.p[(j + 1) % B.n], c1, c2);
      if (d2 < best) {
        best = d2;
        pa = c1;
        pb = c2;
      }
    }
  }
  return std::sqrt(best);
}

inline void atomicMin(std::atomic<double>& a, double v)
{
  double cur = a.load(std::memory_order_relaxed);
  while (v < cur && !a.compare_exchange_weak(cur, v, std::memory_order_relaxed)) {
  }
}

inline unsigned workerCount(unsigned requested)
{
#if defined(__EMSCRIPTEN__) && !defined(__EMSCRIPTEN_PTHREADS__)
  (void)requested;
  return 1;
#else
  if (requested) return requested;
  return std::max(1u, std::thread::hardware_concurrency());
#endif
}

inline std::shared_ptr<const PolySet> triangulated(const PolySet& ps)
{
  if (ps.isTriangular()) return std::make_shared<PolySet>(ps);
  return PolySetUtils::tessellate_faces(ps);
}

#ifdef ENABLE_MANIFOLD
// Move every vertex outwards so that each incident face moves by about
// 'grow' (see drc_mesh.cc). Returns m unchanged if the result is invalid.
manifold::Manifold growAlongNormals(const manifold::Manifold& m, double grow);
#endif

}  // namespace DrcMesh
