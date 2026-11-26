/*
 *  PythonSCAD - shared mesh helpers for design rule checks
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */

#include "geometry/drc_mesh.h"

#include <Eigen/QR>

namespace DrcMesh {

#ifdef ENABLE_MANIFOLD
// Move every vertex outwards so that each incident face moves by 'grow':
// per vertex the minimum-norm least-squares solution of n . f_i = grow over
// the distinct incident face normals f_i (exact at box corners and edges),
// capped at 4 * grow. Not an exact offset, but O(N) and enough to lift the
// error solid off the checked surface so both can be shown without
// z-fighting.
manifold::Manifold growAlongNormals(const manifold::Manifold& m, double grow)
{
  manifold::MeshGL64 mesh = m.GetMeshGL64();
  const size_t np = mesh.numProp, nv = mesh.vertProperties.size() / np;
  std::vector<size_t> rep(nv);  // seam duplicates share one displacement
  for (size_t i = 0; i < nv; ++i) rep[i] = i;
  for (size_t i = 0; i < mesh.mergeFromVert.size(); ++i)
    rep[mesh.mergeFromVert[i]] = mesh.mergeToVert[i];
  auto pos = [&](size_t v) {
    return Vector3d(mesh.vertProperties[v * np], mesh.vertProperties[v * np + 1],
                    mesh.vertProperties[v * np + 2]);
  };
  // distinct unit face normals per vertex
  std::vector<std::vector<Vector3d>> normals(nv);
  const size_t nt = mesh.triVerts.size() / 3;
  for (size_t t = 0; t < nt; ++t) {
    const size_t i0 = mesh.triVerts[3 * t], i1 = mesh.triVerts[3 * t + 1], i2 = mesh.triVerts[3 * t + 2];
    const Vector3d n = (pos(i1) - pos(i0)).cross(pos(i2) - pos(i0));
    const double l = n.norm();
    if (!(l > 0)) continue;
    const Vector3d u = n / l;
    for (size_t v : {i0, i1, i2}) {
      auto& list = normals[rep[v]];
      bool dup = false;
      for (const Vector3d& w : list) dup = dup || w.dot(u) > 0.9999;
      if (!dup) list.push_back(u);
    }
  }
  std::vector<Vector3d> disp(nv, Vector3d::Zero());
  for (size_t v = 0; v < nv; ++v) {
    const auto& list = normals[v];
    if (list.empty()) continue;
    Eigen::MatrixXd F(list.size(), 3);
    for (size_t i = 0; i < list.size(); ++i) F.row(i) = list[i].transpose();
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(F);
    cod.setThreshold(1e-3);
    Vector3d d = cod.solve(Eigen::VectorXd::Constant(list.size(), grow));
    if (d.norm() > 4 * grow) d *= 4 * grow / d.norm();
    disp[v] = d;
  }
  for (size_t v = 0; v < nv; ++v)
    for (int k = 0; k < 3; ++k) mesh.vertProperties[v * np + k] += disp[rep[v]][k];
  manifold::Manifold out(mesh);
  return out.Status() == manifold::Manifold::Error::NoError ? out : m;
}
#endif

}  // namespace DrcMesh
