/*
 *  OpenSCAD (www.openscad.org)
 *  Copyright (C) 2009-2011 Clifford Wolf <clifford@clifford.at> and
 *                          Marius Kintel <marius@kintel.net>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  As a special exception, you have permission to link this program
 *  with the CGAL library and distribute executables, as long as you
 *  follow the requirements of the GNU GPL in regard to all of the
 *  software in the executable aside from CGAL.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include "OversampleNode.h"
#include "module.h"
#include "ModuleInstantiation.h"
#include "Children.h"
#include "Parameters.h"
#include "src/utils/printutils.h"
#include "io/fileutils.h"
#include "Builtins.h"
#include "handle_dep.h"
#include "src/geometry/PolySetBuilder.h"

#include <cmath>
#include <sstream>

#include <src/geometry/PolySetUtils.h>
#include <src/core/Tree.h>
#include <src/geometry/GeometryEvaluator.h>
#include <boost/functional/hash.hpp>
#include <src/utils/hash.h>

void ov_add_poly(PolySetBuilder& builder, Vector3d p)
{
  builder.addVertex(builder.vertexIndex(p));
}

std::unique_ptr<const Geometry> ov_dynamic(const std::shared_ptr<const PolySet>& ps, double limit)
{
  // one round: walk all edges and test to sploit them
  auto ps_work = PolySetUtils::tessellate_faces(*ps);

  bool done = true;
  while (done) {
    done = false;
    PolySet ps_new(3);
    std::unordered_map<uint64_t, int> edges;
    ps_new.vertices = ps_work->vertices;
    ps_new.vertices.reserve(ps_work->vertices.size() * 2);
    ps_new.indices.reserve(ps_work->indices.size() * 4);

    // decide which vertices to split
    for (const auto& tri : ps_work->indices) {
      int ind_old = tri[2];
      for (int ind : tri) {
        double dist2 = (ps_work->vertices[ind_old] - ps_work->vertices[ind]).squaredNorm();
        if (dist2 > limit * limit) {
          uint64_t key = ((uint64_t)std::min(ind, ind_old) << 32) | std::max(ind, ind_old);
          auto it = edges.find(key);
          if (it == edges.end()) {
            Vector3d pmid = (ps_work->vertices[ind_old] + ps_work->vertices[ind]) / 2;

            int newind = ps_new.vertices.size();
            ps_new.vertices.push_back(pmid);
            edges[key] = newind;
            done = true;
          }
        }
        ind_old = ind;
      }
    }
    // build new triangles from split vertices
    for (const auto& tri : ps_work->indices) {
      int ind_old = tri[2];
      int tri_mid[3];
      for (int i = 0; i < 3; i++) {
        int ind = tri[i];
        uint64_t key = ((uint64_t)std::min(ind, ind_old) << 32) | std::max(ind, ind_old);
        auto it = edges.find(key);
        tri_mid[(i + 2) % 3] = (it != edges.end()) ? it->second : -1;
        ind_old = ind;
      }
      switch (((tri_mid[0] >= 0) ? 1 : 0) | ((tri_mid[1] >= 0) ? 2 : 0) | ((tri_mid[2] >= 0) ? 4 : 0)) {
      case 0: ps_new.indices.push_back(tri); break;
      case 1:
        ps_new.indices.push_back({tri[0], tri_mid[0], tri[2]});
        ps_new.indices.push_back({tri[2], tri_mid[0], tri[1]});
        break;
      case 2:
        ps_new.indices.push_back({tri[1], tri_mid[1], tri[0]});
        ps_new.indices.push_back({tri[0], tri_mid[1], tri[2]});
        break;
      case 3:
        ps_new.indices.push_back({tri[0], tri_mid[0], tri[2]});
        ps_new.indices.push_back({tri[2], tri_mid[0], tri_mid[1]});
        ps_new.indices.push_back({tri_mid[1], tri_mid[0], tri[1]});
        break;
      case 4:
        ps_new.indices.push_back({tri[2], tri_mid[2], tri[1]});
        ps_new.indices.push_back({tri[1], tri_mid[2], tri[0]});
        break;
      case 5:
        ps_new.indices.push_back({tri[2], tri_mid[2], tri[1]});
        ps_new.indices.push_back({tri[1], tri_mid[2], tri_mid[0]});
        ps_new.indices.push_back({tri_mid[0], tri_mid[2], tri[0]});
        break;
      case 6:
        ps_new.indices.push_back({tri[1], tri_mid[1], tri[0]});
        ps_new.indices.push_back({tri[0], tri_mid[1], tri_mid[2]});
        ps_new.indices.push_back({tri_mid[2], tri_mid[1], tri[2]});
        break;
      case 7:
        ps_new.indices.push_back({tri[0], tri_mid[0], tri_mid[2]});
        ps_new.indices.push_back({tri[1], tri_mid[1], tri_mid[0]});
        ps_new.indices.push_back({tri[2], tri_mid[2], tri_mid[1]});
        ps_new.indices.push_back({tri_mid[0], tri_mid[1], tri_mid[2]});
        break;
      }
    }
    if (done) *ps_work = ps_new;

    // now check which edges can be flipped
    std::unordered_map<EdgeKey, EdgeVal, boost::hash<EdgeKey>> edgeDb;
    edgeDb = createEdgeDb(ps_work->indices);
    for (auto it = edgeDb.begin(); it != edgeDb.end(); it++) {
      const EdgeKey& key = it->first;
      const EdgeVal& val = it->second;

      int left = -1, right = -1, afar = -1, bfar = -1, tmp;
      for (int i = 0; i < 3; i++) {
        tmp = ps_work->indices[val.facea][i];
        if (tmp != key.ind1 && tmp != key.ind2) {
          afar = ps_work->indices[val.facea][i];
          left = ps_work->indices[val.facea][(i + 1) % 3];
          right = ps_work->indices[val.facea][(i + 2) % 3];
        }
        tmp = ps_work->indices[val.faceb][i];
        if (tmp != key.ind1 && tmp != key.ind2) {
          bfar = ps_work->indices[val.faceb][i];
        }
      }
      const Vector3d& pleft = ps_work->vertices[left];
      const Vector3d& pright = ps_work->vertices[right];
      const Vector3d& pafar = ps_work->vertices[afar];
      const Vector3d& pbfar = ps_work->vertices[bfar];
      const Vector3d phor = pright - pleft;
      const Vector3d pver = pafar - pbfar;

      while (1) {
        if (pver.squaredNorm() > phor.squaredNorm()) break;
        EdgeKey opp1k(left, bfar);
        auto it = edgeDb.find(opp1k);
        if (it == edgeDb.end()) break;
        EdgeVal& opp1 = it->second;

        EdgeKey opp2k(right, afar);
        it = edgeDb.find(opp2k);
        if (it == edgeDb.end()) break;
        EdgeVal& opp2 = it->second;

        Vector3d n = pver.cross(phor);
        if (fabs(n.dot(pafar) - n.dot(pbfar)) > 1e-3) break;  // same level

        if ((pleft - pbfar).cross(pver).dot(n) <= 1e-3) break;   // triangle flip
        if ((pbfar - pright).cross(pver).dot(n) <= 1e-3) break;  // triangle flip

        EdgeKey keyt(bfar, afar);

        for (int i = 0; i < 3; i++) {
          if (ps_work->indices[val.facea][i] == right) ps_work->indices[val.facea][i] = bfar;
          if (ps_work->indices[val.faceb][i] == left) ps_work->indices[val.faceb][i] = afar;
        }
        if (opp1.facea == val.faceb) opp1.facea = val.facea;
        if (opp1.faceb == val.faceb) opp1.faceb = val.facea;
        if (opp2.facea == val.facea) opp2.facea = val.faceb;
        if (opp2.faceb == val.facea) opp2.faceb = val.faceb;

        break;
      }
    }
  }
  return std::make_unique<PolySet>(*ps_work);
}
// TODO check if triangles are changed sign
std::unique_ptr<const Geometry> ov_static(const std::shared_ptr<const PolySet>& ps, int n)
{
  // tesselate object
  auto ps_tess = PolySetUtils::tessellate_faces(*ps);
  std::vector<Vector3d> pt_dir;
  std::unordered_map<Vector3d, int, boost::hash<Vector3d>> pointIntMap;
  PolySetBuilder builder_ov(0, 0, 3, true);
  for (size_t i = 0; i < ps_tess->indices.size(); i++) {
    auto& pol = ps_tess->indices[i];
    Vector3d p1 = ps_tess->vertices[pol[0]];
    Vector3d p2 = ps_tess->vertices[pol[1]];
    Vector3d p3 = ps_tess->vertices[pol[2]];
    Vector3d p21 = (p2 - p1) / n;
    Vector3d p31 = (p3 - p1) / n;
    Vector3d botlast, botcur, toplast, topcur;
    double r = 1.0;
    Vector3d center(0, 0, 0);
    for (int j = 0; j < n; j++) {
      botcur = p1 + p31 * j;
      topcur = p1 + p31 * (j + 1);

      for (int k = 0; k < n - j; k++) {
        if (k != 0) {
          toplast = topcur;
          topcur = topcur + p21;
          builder_ov.beginPolygon(3);
          ov_add_poly(builder_ov, botcur);
          ov_add_poly(builder_ov, topcur);
          ov_add_poly(builder_ov, toplast);
        }
        botlast = botcur;
        botcur = botlast + p21;
        builder_ov.beginPolygon(3);
        ov_add_poly(builder_ov, botlast);
        ov_add_poly(builder_ov, botcur);
        ov_add_poly(builder_ov, topcur);
      }
    }
  }
  return builder_ov.build();
}

std::unique_ptr<const Geometry> OversampleNode::createGeometry() const
{
  if (this->children.size() == 0) {
    return std::unique_ptr<PolySet>();
  }
  std::shared_ptr<AbstractNode> child = this->children[0];
  Tree tree(child, "");
  GeometryEvaluator geomevaluator(tree);
  std::shared_ptr<const Geometry> geom = geomevaluator.evaluateGeometry(*tree.root(), true);
  std::shared_ptr<const PolySet> ps = PolySetUtils::getGeometryAsPolySet(geom);
  if (ps == nullptr) return std::unique_ptr<PolySet>();
  if (this->method == "dynamic") return ov_dynamic(ps, this->n);
  return ov_static(ps, this->n);
}
