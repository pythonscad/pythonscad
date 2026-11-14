/*
 *  PythonSCAD - slope check: overhang and draft / undercut
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */

#include "geometry/slope_check.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <thread>
#include <utility>

#include "geometry/drc_mesh.h"
#ifdef ENABLE_MANIFOLD
#include "geometry/manifold/ManifoldGeometry.h"
#include "geometry/manifold/manifoldutils.h"
#include <manifold/manifold.h>
#endif

namespace SlopeCheck {
namespace {

using namespace DrcMesh;

// Area of a planar polygon projected onto the plane perpendicular to u.
double projectedArea(const Poly& P, const Vector3d& u)
{
  Vector3d a = Vector3d::Zero();
  for (int i = 0; i < P.n; ++i) a += P.p[i].cross(P.p[(i + 1) % P.n]);
  return 0.5 * std::fabs(a.dot(u));
}

Vector3d centroid(const Poly& P)
{
  Vector3d c = Vector3d::Zero();
  for (int i = 0; i < P.n; ++i) c += P.p[i];
  return c / double(P.n);
}

struct Piece {
  Poly P;
  Vector3d u;
};

// Occluder of a piece: triangle C covering the part Q of the piece when
// looking along u, at depth t (measured at the centroid of Q).
struct Occluder {
  int tri;
  Poly Q;
  Vector3d q;
  double t;
};

struct Kernel {
  const std::vector<Tri>& tris;
  const Bvh& bvh;
  const Options& opt;
  Vector3d d;
  double eps, eps_area, tol, reach, base_level;
  double sin_min, sin_max;

  // inward normals of the prism of triangle C along u
  void prism(const Tri& C, const Vector3d& u, Vector3d m[3], double c[3]) const
  {
    const Vector3d cc = (C.v[0] + C.v[1] + C.v[2]) / 3.0;
    for (int k = 0; k < 3; ++k) {
      const Vector3d& a = C.v[k];
      const Vector3d& b = C.v[(k + 1) % 3];
      m[k] = u.cross(b - a);
      if (m[k].dot(cc - a) < 0) m[k] = -m[k];
      const double l = m[k].norm();
      if (l > 0) m[k] /= l;
      c[k] = -m[k].dot(a);
    }
  }

  // depth of C's plane above x along u, NaN if C is parallel to u
  static double depth(const Tri& C, const Vector3d& x, const Vector3d& u)
  {
    const double nu = C.n.dot(u);
    if (std::fabs(nu) < 1e-9) return std::nan("");
    return (C.off - C.n.dot(x)) / nu;
  }

  // nearest occluders of piece pc (from triangle ti) along its pull direction
  std::vector<Occluder> occluders(int ti, const Piece& pc) const
  {
    std::vector<Occluder> front;
    BoundingBox box;
    for (int k = 0; k < pc.P.n; ++k) {
      box.extend(pc.P.p[k]);
      box.extend(pc.P.p[k] + pc.u * reach);
    }
    bvh.query(box, [&](int ci) {
      if (ci == ti) return;
      const Tri& C = tris[ci];
      if (std::fabs(C.n.dot(pc.u)) < 1e-9) return;  // parallel to u: no shadow area
      Vector3d m[3];
      double c[3];
      prism(C, pc.u, m, c);
      Poly Q = pc.P;
      for (int k = 0; k < 3 && Q.n; ++k) clip(Q, m[k], c[k], eps);
      if (Q.n < 3 || projectedArea(Q, pc.u) < eps_area) return;
      const Vector3d q = centroid(Q);
      const double t = depth(C, q, pc.u);
      if (!(t > tol)) return;  // behind the piece
      front.push_back({ci, Q, q, t});
    });
    // keep only occluders that are nearest somewhere: drop k if another
    // occluder covers the centroid of k's shadow at a smaller depth
    std::vector<Occluder> nearest;
    for (size_t k = 0; k < front.size(); ++k) {
      bool hidden = false;
      for (size_t j = 0; j < front.size() && !hidden; ++j) {
        if (j == k) continue;
        const Tri& Cj = tris[front[j].tri];
        Vector3d m[3];
        double c[3];
        prism(Cj, pc.u, m, c);
        bool inside = true;
        for (int e = 0; e < 3; ++e) inside = inside && m[e].dot(front[k].q) + c[e] >= -eps;
        if (!inside) continue;
        const double tj = depth(Cj, front[k].q, pc.u);
        hidden = tj > tol && tj < front[k].t - tol;
      }
      if (!hidden) nearest.push_back(front[k]);
    }
    return nearest;
  }

  void pieces(const Tri& T, std::vector<Piece>& out) const
  {
    out.clear();
    Piece whole;
    whole.P.n = 3;
    for (int k = 0; k < 3; ++k) whole.P.p[k] = T.v[k];
    switch (opt.parting) {
    case Parting::None:
      whole.u = d;
      out.push_back(whole);
      break;
    case Parting::Free:
      whole.u = T.n.dot(d) >= -1e-12 ? d : Vector3d(-d);
      out.push_back(whole);
      break;
    case Parting::Plane: {
      Piece up = whole, dn = whole;
      clip(up.P, d, -opt.parting_pos, 0.0);  // d . x >= pos
      clip(dn.P, -d, opt.parting_pos, 0.0);  // d . x <= pos
      up.u = d;
      dn.u = -d;
      for (Piece *p : {&up, &dn})
        if (p->P.n >= 3 && projectedArea(p->P, T.n) > eps_area) out.push_back(*p);
      break;
    }
    }
  }

  bool onBase(const Poly& P) const
  {
    if (!opt.skip_base) return false;
    for (int k = 0; k < P.n; ++k)
      if (d.dot(P.p[k]) > base_level + tol) return false;
    return true;
  }
};

}  // namespace

Result check(const PolySet& ps, const Options& opt_in)
{
  Options opt = opt_in;
  Result res;
  res.mesh = triangulated(ps);
  const PolySet& mesh = *res.mesh;
  const int N = (int)mesh.indices.size();
  if (N == 0 || opt.dir.norm() == 0) return res;

  BoundingBox all;
  for (const Vector3d& v : mesh.vertices) all.extend(v);
  const double diag = std::max(all.sizes().norm(), 1e-12);

  std::vector<Tri> tris(N);
  for (int i = 0; i < N; ++i) {
    const auto& f = mesh.indices[i];
    Tri& T = tris[i];
    T.valid = f.size() == 3;
    if (!T.valid) continue;
    for (int k = 0; k < 3; ++k) {
      T.vi[k] = f[k];
      T.v[k] = mesh.vertices[f[k]];
      T.box.extend(T.v[k]);
    }
    const Vector3d n = (T.v[1] - T.v[0]).cross(T.v[2] - T.v[0]);
    const double l = n.norm();
    T.valid = l > 1e-18 * diag * diag;
    T.n = T.valid ? Vector3d(n / l) : Vector3d::Zero();
    T.off = T.n.dot(T.v[0]);
  }
  const Bvh bvh(tris);

  const Vector3d d = opt.dir.normalized();
  double base = kInf;
  for (const Vector3d& v : mesh.vertices) base = std::min(base, d.dot(v));

  const Kernel K{tris,
                 bvh,
                 opt,
                 d,
                 1e-9 * diag,
                 1e-12 * diag * diag,
                 1e-7 * diag,
                 2 * diag,
                 base,
                 std::sin(opt.min_deg * M_PI / 180.0),
                 std::sin(opt.max_deg * M_PI / 180.0)};

  const unsigned nt = workerCount(opt.threads);
  const int chunk = 256;
  std::atomic<int> next{0};
  std::vector<std::vector<Violation>> per(nt);
  std::vector<std::vector<Vector3d>> pool(nt);

  auto worker = [&](unsigned w) {
    std::vector<Piece> pcs;
    for (;;) {
      const int s = next.fetch_add(chunk);
      if (s >= N) break;
      const int e = std::min(N, s + chunk);
      for (int ti = s; ti < e; ++ti) {
        const Tri& T = tris[ti];
        if (!T.valid) continue;
        K.pieces(T, pcs);
        for (const Piece& pc : pcs) {
          if (K.onBase(pc.P)) continue;
          const double s_beta = std::clamp(T.n.dot(pc.u), -1.0, 1.0);
          int flags = 0;
          if (s_beta < K.sin_min - 1e-12 || s_beta > K.sin_max + 1e-12) flags |= kAngle;
          std::vector<Occluder> occ;
          // only faces looking towards their pull direction can be hidden;
          // faces looking away are angle violations (local undercuts) and
          // their column would run into their own material
          if (opt.undercut && s_beta > 1e-9) {
            occ = K.occluders(ti, pc);
            if (!occ.empty()) flags |= kUndercut;
          }
          if (!flags) continue;
          Violation v;
          v.tri = ti;
          v.flags = flags;
          v.beta_deg = std::asin(s_beta) * 180.0 / M_PI;
          if (std::fabs(v.beta_deg) < 1e-9) v.beta_deg = 0;  // vertical walls: no -1e-14
          v.pull = pc.u;
          auto& pts = pool[w];
          v.piece_begin = (int)pts.size();
          v.piece_count = pc.P.n;
          for (int k = 0; k < pc.P.n; ++k) pts.push_back(pc.P.p[k]);
          for (const Occluder& o : occ) {
            const Tri& C = tris[o.tri];
            const int b = (int)pts.size();
            for (int k = 0; k < o.Q.n; ++k) {
              const Vector3d& x = o.Q.p[k];
              const double t = K.depth(C, x, pc.u);
              pts.push_back(x);
              pts.push_back(x + pc.u * (std::isnan(t) ? o.t : std::max(t, 0.0)));
            }
            v.columns.emplace_back(b, (int)pts.size() - b);
          }
          per[w].push_back(std::move(v));
        }
      }
    }
  };
  if (nt == 1) {
    worker(0);
  } else {
    std::vector<std::thread> threads;
    for (unsigned w = 0; w < nt; ++w) threads.emplace_back(worker, w);
    for (auto& t : threads) t.join();
  }

  for (unsigned w = 0; w < nt; ++w) {
    const int base_idx = (int)res.hull_points.size();
    res.hull_points.insert(res.hull_points.end(), pool[w].begin(), pool[w].end());
    for (Violation& v : per[w]) {
      v.piece_begin += base_idx;
      for (auto& c : v.columns) c.first += base_idx;
      res.violations.push_back(std::move(v));
    }
  }
  std::stable_sort(res.violations.begin(), res.violations.end(),
                   [](const Violation& a, const Violation& b) { return a.tri < b.tri; });

  // per triangle counts and the worst angle
  double worst_excess = -1;
  for (size_t i = 0; i < res.violations.size();) {
    const int tri = res.violations[i].tri;
    int flags = 0;
    for (; i < res.violations.size() && res.violations[i].tri == tri; ++i) {
      const Violation& v = res.violations[i];
      flags |= v.flags;
      if (v.flags & kAngle) {
        const double excess = std::max(opt.min_deg - v.beta_deg, v.beta_deg - opt.max_deg);
        if (excess > worst_excess) {
          worst_excess = excess;
          res.worst_deg = v.beta_deg;
        }
      }
    }
    res.count++;
    if (flags & kAngle) res.angle_count++;
    if (flags & kUndercut) res.undercut_count++;
  }
  return res;
}

std::shared_ptr<const Geometry> errorGeometry(const Result& res,
                                              const std::shared_ptr<const Geometry>& body,
                                              double thickness, double grow, size_t max_hulls)
{
#ifdef ENABLE_MANIFOLD
  if (res.violations.empty()) return nullptr;
  const PolySet& mesh = *res.mesh;
  std::vector<manifold::Manifold> skin, trapped;
  std::vector<manifold::vec3> pts;
  size_t used = 0;
  for (const Violation& v : res.violations) {
    if (used >= max_hulls) break;
    const auto& f = mesh.indices[v.tri];
    const Vector3d n = (mesh.vertices[f[1]] - mesh.vertices[f[0]])
                         .cross(mesh.vertices[f[2]] - mesh.vertices[f[0]])
                         .normalized();
    if ((v.flags & kAngle) && thickness > 0) {
      pts.clear();
      for (int k = 0; k < v.piece_count; ++k) {
        const Vector3d& p = res.hull_points[v.piece_begin + k];
        const Vector3d q = p - n * thickness;
        pts.emplace_back(p[0], p[1], p[2]);
        pts.emplace_back(q[0], q[1], q[2]);
      }
      manifold::Manifold h = manifold::Manifold::Hull(pts);
      if (!h.IsEmpty()) skin.push_back(std::move(h)), ++used;
    }
    for (const auto& [b, cnt] : v.columns) {
      pts.clear();
      for (int k = 0; k < cnt; ++k) {
        const Vector3d& p = res.hull_points[b + k];
        pts.emplace_back(p[0], p[1], p[2]);
      }
      manifold::Manifold h = manifold::Manifold::Hull(pts);
      if (!h.IsEmpty()) trapped.push_back(std::move(h)), ++used;
    }
  }

  manifold::Manifold solid;
  if (body) {
    auto mg = ManifoldUtils::createManifoldFromGeometry(body);
    if (mg && !mg->isEmpty()) solid = mg->getManifold();
  }
  manifold::Manifold err;
  if (!skin.empty()) {
    manifold::Manifold s = manifold::Manifold::BatchBoolean(skin, manifold::OpType::Add);
    err = solid.IsEmpty() ? s : (s ^ solid);
  }
  if (!trapped.empty()) {
    manifold::Manifold t = manifold::Manifold::BatchBoolean(trapped, manifold::OpType::Add);
    if (!solid.IsEmpty()) t = t - solid;
    err = err.IsEmpty() ? t : (err + t);
  }
  if (err.IsEmpty()) return nullptr;
  if (grow > 0) err = growAlongNormals(err, grow);
  return std::make_shared<ManifoldGeometry>(err);
#else
  (void)res;
  (void)body;
  (void)thickness;
  (void)grow;
  (void)max_hulls;
  return nullptr;
#endif
}

}  // namespace SlopeCheck
