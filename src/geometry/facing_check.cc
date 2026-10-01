/*
 *  PythonSCAD - facing surface width / space check (3D INTERNAL / EXTERNAL)
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */

#include "geometry/facing_check.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <thread>
#include <utility>

#include "geometry/PolySetUtils.h"
#ifdef ENABLE_MANIFOLD
#include "geometry/manifold/ManifoldGeometry.h"
#include "geometry/manifold/manifoldutils.h"
#include <manifold/manifold.h>
#endif

namespace FacingCheck {
namespace {

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
void clip(Poly& poly, const Vector3d& h, double c, double eps)
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
double segSeg(const Vector3d& p1, const Vector3d& q1, const Vector3d& p2, const Vector3d& q2,
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
bool pointFace(const Vector3d& x, const Poly& P, const Vector3d& n, double eps, double& d2,
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
double polyDist(const Poly& A, const Vector3d& na, const Poly& B, const Vector3d& nb, double eps,
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

// ------------------------------------------------------------------ kernel

struct Kernel {
  const std::vector<Tri>& tris;
  const Bvh& bvh;
  const std::vector<int>& vstart;  // vertex -> incident triangles (CSR)
  const std::vector<int>& vtris;
  double d, cos_min, tan_alpha, eps;
  bool projecting, flip, occlusion;

  struct Hit {
    Poly PA, PB;
    Vector3d p, q;
    double dist;
  };

  // Clip P to the region that faces triangle dst: behind dst's plane, closer
  // than d to it and, when projecting, inside the cone around dst: lateral
  // distance to dst <= depth * tan(alpha). That region has rounded corners;
  // it is approximated from outside by the tilted edge planes plus three
  // tangent planes per corner (overshoot < 8 %). Pairs are verified exactly
  // afterwards (withinCone).
  void clipTo(Poly& P, const Tri& dst) const
  {
    clip(P, -dst.n, dst.off, eps);     // behind: n . x - off <= 0
    clip(P, dst.n, d - dst.off, eps);  // near:   n . x - off >= -d
    if (!projecting) return;
    Vector3d out[3];  // outward edge normals in dst's plane
    for (int i = 0; i < 3; ++i) {
      const Vector3d& a = dst.v[i];
      const Vector3d& b = dst.v[(i + 1) % 3];
      out[i] = -dst.n.cross(b - a).normalized();
      // dot(out, x - a) <= tan(alpha) * depth, depth = -dot(n, x - a)
      const Vector3d h = -out[i] - dst.n * tan_alpha;
      clip(P, h, -h.dot(a), eps);
    }
    if (tan_alpha <= 0) return;
    for (int i = 0; i < 3; ++i) {  // corner i between edges i-1 and i
      const Vector3d& v = dst.v[i];
      const Vector3d& o1 = out[(i + 2) % 3];
      const Vector3d& o2 = out[i];
      const double phi = std::atan2(o1.cross(o2).dot(dst.n), o1.dot(o2));
      const Vector3d w = dst.n.cross(o1);
      for (int j = 1; j <= 3; ++j) {
        const double t = phi * j / 4;
        const Vector3d u = o1 * std::cos(t) + w * std::sin(t);
        const Vector3d h = -u - dst.n * tan_alpha;
        clip(P, h, -h.dot(v), eps);
      }
    }
  }

  // Distance from the projection of x onto T's plane to triangle T.
  double lateral(const Tri& T, const Vector3d& x) const
  {
    const Vector3d y = x - T.n * (T.n.dot(x) - T.off);
    bool inside = true;
    double best = kInf;
    for (int k = 0; k < 3; ++k) {
      const Vector3d& a = T.v[k];
      const Vector3d& b = T.v[(k + 1) % 3];
      if (T.n.cross(b - a).dot(y - a) < 0) inside = false;
      Vector3d c1, c2;
      best = std::min(best, segSeg(y, y, a, b, c1, c2));
    }
    return inside ? 0.0 : std::sqrt(best);
  }

  // Exact cone test for the measured pair: q seen from A and p seen from B.
  bool withinCone(const Tri& A, const Tri& B, const Vector3d& p, const Vector3d& q) const
  {
    const double tol = 1e3 * eps;
    const double da = A.off - A.n.dot(q), db = B.off - B.n.dot(p);  // depths
    return lateral(A, q) <= tan_alpha * da + tol && lateral(B, p) <= tan_alpha * db + tol;
  }

  bool strictlyInside(const Tri& T, const Vector3d& x) const
  {
    for (int k = 0; k < 3; ++k) {
      const Vector3d e = T.v[(k + 1) % 3] - T.v[k];
      if (e.cross(x - T.v[k]).dot(T.n) <= 1e3 * eps * e.norm()) return false;
    }
    return true;
  }

  // Does direction u, starting at x on the border of triangle ti, enter the
  // checked region right away? Sufficient condition: u points behind every
  // face incident to the vertex / edge x lies on. False when unsure.
  bool enters(int ti, const Vector3d& x, const Vector3d& u) const
  {
    const Tri& T = tris[ti];
    const double tol = 1e3 * eps, tu = 1e-9 * u.norm();
    auto behind_all = [&](int v0, int v1) {  // v1 < 0: whole vertex star
      for (int k = vstart[v0]; k < vstart[v0 + 1]; ++k) {
        const Tri& F = tris[vtris[k]];
        if (!F.valid) continue;
        if (v1 >= 0 && F.vi[0] != v1 && F.vi[1] != v1 && F.vi[2] != v1) continue;
        if (!(F.n.dot(u) < -tu)) return false;
      }
      return true;
    };
    for (int k = 0; k < 3; ++k)
      if ((x - T.v[k]).squaredNorm() <= tol * tol) return behind_all(T.vi[k], -1);
    int best = -1;
    double bd = kInf;
    for (int k = 0; k < 3; ++k) {
      Vector3d c1, c2;
      const double dd = segSeg(x, x, T.v[k], T.v[(k + 1) % 3], c1, c2);
      if (dd < bd) bd = dd, best = k;
    }
    if (bd > tol * tol) return false;
    return behind_all(T.vi[best], T.vi[(best + 1) % 3]);
  }

  bool test(int ia, int ib, Hit& h) const
  {
    const Tri& A = tris[ia];
    const Tri& B = tris[ib];
    if (A.n.dot(B.n) > cos_min) return false;

    // each triangle must reach into the slab [-d, 0] behind the other's plane
    auto reaches = [&](const Tri& P, const Tri& Q) {
      double lo = kInf, hi = -kInf;
      for (int k = 0; k < 3; ++k) {
        const double sd = P.n.dot(Q.v[k]) - P.off;
        lo = std::min(lo, sd);
        hi = std::max(hi, sd);
      }
      return lo <= eps && hi >= -d - eps;
    };
    if (!reaches(A, B) || !reaches(B, A)) return false;

    Poly& PA = h.PA;
    Poly& PB = h.PB;
    PA.n = PB.n = 3;
    for (int k = 0; k < 3; ++k) {
      PA.p[k] = A.v[k];
      PB.p[k] = B.v[k];
    }
    clipTo(PA, B);
    if (PA.n == 0) return false;
    clipTo(PB, A);
    if (PB.n == 0) return false;

    const Vector3d& p = h.p;
    const Vector3d& q = h.q;
    h.dist = polyDist(PA, A.n, PB, B.n, eps, h.p, h.q);
    if (!(h.dist < d)) return false;
    if (projecting && !withinCone(A, B, p, q)) return false;

    if (occlusion && h.dist > 10 * eps) {
      // 1) the connecting segment must not cross the surface
      const double tol = 1e-7;
      if (bvh.segmentHits(p, q, tol, 1 - tol, ia, ib)) return false;
      // 2) starting inside a face (or into the region at a vertex / edge) it
      //    stays in the checked region. Only rim to rim chords can lie in the
      //    wrong region entirely: decide those with a point-in-solid test,
      //    nudged off the boundary towards the centroid of both faces.
      if (!strictlyInside(A, p) && !strictlyInside(B, q) && !enters(ia, p, q - p) &&
          !enters(ib, q, p - q)) {
        Vector3d c = Vector3d::Zero();
        for (int k = 0; k < PA.n; ++k) c += PA.p[k];
        for (int k = 0; k < PB.n; ++k) c += PB.p[k];
        c /= double(PA.n + PB.n);
        const Vector3d m = (p + q) * 0.5;
        if (bvh.inside(m + (c - m) * 0.01) == flip) return false;
      }
    }
    return true;
  }
};

void atomicMin(std::atomic<double>& a, double v)
{
  double cur = a.load(std::memory_order_relaxed);
  while (v < cur && !a.compare_exchange_weak(cur, v, std::memory_order_relaxed)) {
  }
}

unsigned workerCount(unsigned requested)
{
#if defined(__EMSCRIPTEN__) && !defined(__EMSCRIPTEN_PTHREADS__)
  (void)requested;
  return 1;
#else
  if (requested) return requested;
  return std::max(1u, std::thread::hardware_concurrency());
#endif
}

std::shared_ptr<const PolySet> triangulated(const PolySet& ps)
{
  if (ps.isTriangular()) return std::make_shared<PolySet>(ps);
  return PolySetUtils::tessellate_faces(ps);
}

// split > 0: triangles [0, split) and [split, N) are two solids, only pairs
// between them are tested.
Result run(std::shared_ptr<const PolySet> mesh, double d, const Options& opt, bool flip, int split)
{
  Result res;
  res.mesh = mesh;
  res.split = split;
  res.min_distance = kInf;
  const int N = (int)mesh->indices.size();
  res.tri_min_distance.assign(N, kInf);
  if (!(d > 0) || N == 0) return res;

  BoundingBox all;
  for (const Vector3d& v : mesh->vertices) all.extend(v);
  const double eps = 1e-9 * std::max(all.sizes().norm(), d);

  // prepare triangles (flip: reverse winding and normal -> complement)
  std::vector<Tri> tris(N);
  for (int i = 0; i < N; ++i) {
    const auto& f = mesh->indices[i];
    Tri& T = tris[i];
    T.valid = f.size() == 3;
    if (!T.valid) continue;
    T.vi[0] = f[0];
    T.vi[1] = flip ? f[2] : f[1];
    T.vi[2] = flip ? f[1] : f[2];
    for (int k = 0; k < 3; ++k) {
      T.v[k] = mesh->vertices[T.vi[k]];
      T.box.extend(T.v[k]);
    }
    const Vector3d n = (T.v[1] - T.v[0]).cross(T.v[2] - T.v[0]);
    const double l = n.norm();
    T.valid = l > eps * eps;
    T.n = T.valid ? Vector3d(n / l) : Vector3d::Zero();
    T.off = T.n.dot(T.v[0]);
  }

  const Bvh bvh(tris);

  std::vector<int> vstart(mesh->vertices.size() + 1, 0), vtris;
  for (const Tri& T : tris)
    if (T.valid)
      for (int k = 0; k < 3; ++k) ++vstart[T.vi[k] + 1];
  for (size_t i = 1; i < vstart.size(); ++i) vstart[i] += vstart[i - 1];
  vtris.resize(vstart.back());
  {
    std::vector<int> fill(vstart.begin(), vstart.end() - 1);
    for (int i = 0; i < N; ++i)
      if (tris[i].valid)
        for (int k = 0; k < 3; ++k) vtris[fill[tris[i].vi[k]]++] = i;
  }

  const bool projecting = opt.alpha_deg < 90.0;
  const Kernel K{tris,
                 bvh,
                 vstart,
                 vtris,
                 d,
                 std::cos(opt.min_angle_deg * M_PI / 180.0),
                 projecting ? std::tan(std::max(0.0, opt.alpha_deg) * M_PI / 180.0) : 0.0,
                 eps,
                 projecting,
                 flip,
                 opt.occlusion};

  const bool pairs = opt.store_pairs, hulls = opt.store_pairs && opt.build_hulls;
  const unsigned nt = workerCount(opt.threads);
  const int chunk = 256;
  std::atomic<int> next{0};
  std::atomic<size_t> count{0};
  std::vector<std::atomic<double>> tmin(N);
  for (auto& a : tmin) a.store(kInf, std::memory_order_relaxed);
  std::vector<std::vector<Violation>> per(nt);
  std::vector<std::vector<Vector3d>> pool(nt);

  auto worker = [&](unsigned w) {
    Kernel::Hit h;
    size_t local = 0;
    for (;;) {
      const int s = next.fetch_add(chunk);
      if (s >= N) break;
      const int e = std::min(N, s + chunk);
      for (int ia = s; ia < e; ++ia) {
        if (!tris[ia].valid) continue;
        bvh.query(expanded(tris[ia].box, d), [&](int ib) {
          if (ib <= ia) return;  // each unordered pair once
          if (split > 0 && (ia < split) == (ib < split)) return;
          if (!K.test(ia, ib, h)) return;
          ++local;
          atomicMin(tmin[ia], h.dist);
          atomicMin(tmin[ib], h.dist);
          if (!pairs) return;
          Violation v;
          v.tri_a = ia;
          v.tri_b = ib;
          v.p = h.p;
          v.q = h.q;
          v.distance = h.dist;
          if (hulls) {
            v.hull_begin = (int)pool[w].size();
            v.hull_count = h.PA.n + h.PB.n;
            for (int k = 0; k < h.PA.n; ++k) pool[w].push_back(h.PA.p[k]);
            for (int k = 0; k < h.PB.n; ++k) pool[w].push_back(h.PB.p[k]);
          }
          per[w].push_back(v);
        });
      }
    }
    count += local;
  };

  if (nt == 1) {
    worker(0);
  } else {
    std::vector<std::thread> threads;
    for (unsigned w = 0; w < nt; ++w) threads.emplace_back(worker, w);
    for (auto& t : threads) t.join();
  }

  res.violation_count = count;
  for (int i = 0; i < N; ++i) {
    res.tri_min_distance[i] = tmin[i].load();
    res.min_distance = std::min(res.min_distance, res.tri_min_distance[i]);
  }
  if (pairs) {
    size_t total = 0, points = 0;
    for (unsigned w = 0; w < nt; ++w) total += per[w].size(), points += pool[w].size();
    res.violations.reserve(total);
    res.hull_points.reserve(points);
    for (unsigned w = 0; w < nt; ++w) {
      const int base = (int)res.hull_points.size();
      res.hull_points.insert(res.hull_points.end(), pool[w].begin(), pool[w].end());
      for (Violation& v : per[w]) {
        v.hull_begin += base;
        res.violations.push_back(v);
      }
      std::vector<Violation>().swap(per[w]);  // release early, keeps the peak low
      std::vector<Vector3d>().swap(pool[w]);
    }
    // deterministic order, independent of thread scheduling
    std::sort(res.violations.begin(), res.violations.end(), [](const Violation& a, const Violation& b) {
      return a.tri_a != b.tri_a ? a.tri_a < b.tri_a : a.tri_b < b.tri_b;
    });
  }
  return res;
}

}  // namespace

Result check(const PolySet& ps, Mode mode, double distance, const Options& opt)
{
  return run(triangulated(ps), distance, opt, mode == Mode::External, 0);
}

Result checkBetween(const PolySet& a, const PolySet& b, double distance, const Options& opt)
{
  const auto ta = triangulated(a), tb = triangulated(b);
  auto merged = std::make_shared<PolySet>(3);
  merged->vertices = ta->vertices;
  merged->vertices.insert(merged->vertices.end(), tb->vertices.begin(), tb->vertices.end());
  merged->indices = ta->indices;
  const int nv = (int)ta->vertices.size();
  for (auto f : tb->indices) {
    for (auto& k : f) k += nv;
    merged->indices.push_back(f);
  }
  merged->setTriangular(true);
  const int split = (int)ta->indices.size();
  if (split == 0 || tb->indices.empty()) {
    Result r;
    r.mesh = merged;
    r.split = split;
    r.min_distance = kInf;
    r.tri_min_distance.assign(merged->indices.size(), kInf);
    return r;
  }
  return run(merged, distance, opt, true, split);
}

std::shared_ptr<const Geometry> errorGeometry(const Result& res, Mode mode,
                                              const std::vector<std::shared_ptr<const Geometry>>& bodies,
                                              size_t max_hulls)
{
#ifdef ENABLE_MANIFOLD
  if (res.violations.empty()) return nullptr;
  // most severe pairs first
  std::vector<const Violation *> order;
  order.reserve(res.violations.size());
  for (const auto& v : res.violations)
    if (v.hull_count >= 4) order.push_back(&v);
  std::stable_sort(order.begin(), order.end(),
                   [](const Violation *a, const Violation *b) { return a->distance < b->distance; });
  if (order.size() > max_hulls) order.resize(max_hulls);

  std::vector<manifold::Manifold> parts;
  parts.reserve(order.size());
  std::vector<manifold::vec3> pts;
  for (const Violation *v : order) {
    pts.clear();
    for (int k = 0; k < v->hull_count; ++k) {
      const Vector3d& p = res.hull_points[v->hull_begin + k];
      pts.emplace_back(p[0], p[1], p[2]);
    }
    manifold::Manifold h = manifold::Manifold::Hull(pts);
    if (!h.IsEmpty()) parts.push_back(std::move(h));  // zero-volume contacts drop out
  }
  if (parts.empty()) return nullptr;
  manifold::Manifold err = manifold::Manifold::BatchBoolean(parts, manifold::OpType::Add);

  // clip to the material (Internal) or to the air (External)
  std::vector<manifold::Manifold> solids;
  for (const auto& g : bodies) {
    if (!g) continue;
    auto mg = ManifoldUtils::createManifoldFromGeometry(g);
    if (mg && !mg->isEmpty()) solids.push_back(mg->getManifold());
  }
  if (!solids.empty()) {
    const manifold::Manifold body = manifold::Manifold::BatchBoolean(solids, manifold::OpType::Add);
    err = mode == Mode::Internal ? (err ^ body) : (err - body);
  }
  if (err.IsEmpty()) return nullptr;
  return std::make_shared<ManifoldGeometry>(err);
#else
  (void)res;
  (void)mode;
  (void)bodies;
  (void)max_hulls;
  return nullptr;
#endif
}

}  // namespace FacingCheck
