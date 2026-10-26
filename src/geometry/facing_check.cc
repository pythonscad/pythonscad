/*
 *  PythonSCAD - facing surface width / space check (3D INTERNAL / EXTERNAL)
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */

#include "geometry/facing_check.h"
#include "geometry/drc_mesh.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <thread>
#include <utility>

#include "geometry/PolySetUtils.h"
#include <Eigen/QR>
#ifdef ENABLE_MANIFOLD
#include "geometry/manifold/ManifoldGeometry.h"
#include "geometry/manifold/manifoldutils.h"
#include <manifold/manifold.h>
#endif

namespace FacingCheck {
namespace {

using namespace DrcMesh;

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

  static double area(const Poly& P, const Vector3d& n)
  {
    Vector3d a = Vector3d::Zero();
    for (int i = 0; i < P.n; ++i) a += P.p[i].cross(P.p[(i + 1) % P.n]);
    return 0.5 * a.dot(n);
  }

  // Shrink the hull input of a violating pair to the column between the two
  // faces: the larger clipped face is cut down to the projection of the
  // smaller one. Hulls of neighbouring pairs then tile instead of overlapping,
  // which keeps the union of all hulls cheap. Oblique pairs without overlap
  // in projection keep the full clipped faces.
  void hullInput(const Tri& A, const Tri& B, Poly& PA, Poly& PB) const
  {
    const bool a_big = std::fabs(area(PA, A.n)) >= std::fabs(area(PB, B.n));
    Poly& L = a_big ? PA : PB;
    const Poly& Sm = a_big ? PB : PA;
    const Tri& TL = a_big ? A : B;
    if (Sm.n < 3) return;
    Poly proj;
    proj.n = Sm.n;
    for (int i = 0; i < Sm.n; ++i) proj.p[i] = Sm.p[i] - TL.n * (TL.n.dot(Sm.p[i]) - TL.off);
    const double ar = area(proj, TL.n);
    if (std::fabs(ar) <= eps * eps) return;
    const double sign = ar > 0 ? 1.0 : -1.0;
    Poly cut = L;
    for (int i = 0; i < proj.n && cut.n; ++i) {
      const Vector3d& a = proj.p[i];
      const Vector3d m = TL.n.cross(proj.p[(i + 1) % proj.n] - a) * sign;  // inward
      clip(cut, m, -m.dot(a), eps);
    }
    if (cut.n >= 3) L = cut;
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
            K.hullInput(tris[ia], tris[ib], h.PA, h.PB);
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
                                              size_t max_hulls, double grow)
{
#ifdef ENABLE_MANIFOLD
  if (res.violations.empty()) return nullptr;
  // Every violating triangle contributes the hull to its nearest partner
  // only. Between two curved surfaces each facet faces dozens of facets on
  // the other side; the hulls of all those pairs cross the same thin web and
  // their union explodes in complexity, while the nearest-partner hulls
  // already cover the thin region.
  const size_t ntri = res.tri_min_distance.size();
  std::vector<int> best(ntri, -1);
  for (int i = 0; i < (int)res.violations.size(); ++i) {
    const Violation& v = res.violations[i];
    if (v.hull_count < 4) continue;
    for (int t : {v.tri_a, v.tri_b})
      if (best[t] < 0 || v.distance < res.violations[best[t]].distance) best[t] = i;
  }
  std::sort(best.begin(), best.end());
  best.erase(std::unique(best.begin(), best.end()), best.end());
  std::vector<const Violation *> order;
  order.reserve(best.size());
  for (int i : best)
    if (i >= 0) order.push_back(&res.violations[i]);
  // most severe first when capped
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
  if (grow > 0) err = growAlongNormals(err, grow);
  return std::make_shared<ManifoldGeometry>(err);
#else
  (void)res;
  (void)mode;
  (void)bodies;
  (void)max_hulls;
  (void)grow;
  return nullptr;
#endif
}

}  // namespace FacingCheck
