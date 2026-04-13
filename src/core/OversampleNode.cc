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

#include "SurfaceNode.h"
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
#include "lodepng/lodepng.h"

const char *projectionNames[] = {"none",        "triplanar", "cubic",   "spherical",
                                 "cylindrical", "planarx",   "planary", "planarz"};

double OversampleNode::tcoord(std::shared_ptr<img_data_t> tex, double x, double y) const
{
  double u = x / texturewidth;
  double v = y / textureheight;
  // get texture coorindate
  if (u > 0) u -= (int)u;
  else {
    u = -u - (int)(-u);
    u = 1 - u;
  }
  if (v > 0) v -= (int)v;
  else {
    v = -v - (int)(-v);
    v = 1 - v;
  }

  int tx = (tex->width - 1) * u;
  int ty = (tex->height - 1) * v;
  Vector3f pixel = (*tex)[ty * tex->width + tx];
  double depth = (pixel[0] + pixel[1] + pixel[2]) * texturedepth / (3.0 * 256.0);
  return depth;
}

bool is_png(std::vector<uint8_t>& png)
{
  const size_t pngHeaderLength = 8;
  const uint8_t pngHeader[pngHeaderLength] = {0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a};
  return (png.size() >= pngHeaderLength && std::memcmp(png.data(), pngHeader, pngHeaderLength) == 0);
}

void convert_image(std::shared_ptr<img_data_t> data, std::vector<uint8_t>& img, unsigned int width,
                   unsigned int height)
{
  data->width = width;
  data->height = height;
  data->resize((size_t)width * height);
  for (unsigned int y = 0; y < height; ++y) {
    for (unsigned int x = 0; x < width; ++x) {
      long idx = 4l * (y * width + x);
      (*data)[x + (width * (height - 1 - y))] = Vector3f(img[idx], img[idx + 1], img[idx + 2]);
    }
  }
}

std::shared_ptr<img_data_t> load_png(const std::string& filename)
{
  std::shared_ptr<img_data_t> data = std::make_shared<img_data_t>();
  std::vector<uint8_t> png;
  int ret_val = 0;
  try {
    ret_val = lodepng::load_file(png, filename);
  } catch (std::bad_alloc& ba) {
    LOG(message_group::Warning, "bad_alloc caught for '%1$s'.", ba.what());
    return data;
  }

  if (ret_val == 78) {
    LOG(message_group::Warning, "The file '%1$s' couldn't be opened.", filename);
    return data;
  }

  if (!is_png(png)) {
  }

  unsigned int width, height;
  std::vector<uint8_t> img;
  auto error = lodepng::decode(img, width, height, png);
  if (error) {
    LOG(message_group::Warning, "Can't read PNG image '%1$s'", filename);
    return data;
  }

  convert_image(data, img, width, height);

  return data;
}

std::unique_ptr<const Geometry> OversampleNode::createGeometry_sub(
  const std::shared_ptr<const PolySet>& ps) const
{
  auto ps_work = PolySetUtils::tessellate_faces(*ps);

  // calculate orginal norm vectors
  std::vector<Vector3d> normals;
  normals.reserve(ps_work->indices.size());
  std::vector<int> orig_id;
  orig_id.reserve(ps_work->indices.size());

  for (const auto& ind : ps_work->indices) {
    Vector3d p1 = ps->vertices[ind[0]];
    Vector3d p2 = ps->vertices[ind[1]];
    Vector3d p3 = ps->vertices[ind[2]];
    normals.push_back((p2 - p1).cross(p3 - p1).normalized());
  }
  for (int i = 0; i < ps_work->indices.size(); i++) orig_id.push_back(i);

  bool done = true;
  while (done) {
    done = false;
    PolySet ps_new(3);
    std::unordered_map<uint64_t, int> edges;
    ps_new.vertices = ps_work->vertices;
    ps_new.vertices.reserve(ps_work->vertices.size() * 2);
    ps_new.indices.reserve(ps_work->indices.size() * 4);
    std::vector<int> new_id;
    new_id.reserve(ps_work->indices.size() * 4);

    // decide which vertices to split
    for (const auto& tri : ps_work->indices) {
      int ind_old = tri[2];
      for (int ind : tri) {
        double dist2 = (ps_work->vertices[ind_old] - ps_work->vertices[ind]).squaredNorm();
        if (dist2 > size * size) {
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
    for (int i = 0; i < ps_work->indices.size(); i++) {
      const auto& tri = ps_work->indices[i];
      int ind_old = tri[2];
      int tri_mid[3];
      for (int i = 0; i < 3; i++) {
        int ind = tri[i];
        uint64_t key = ((uint64_t)std::min(ind, ind_old) << 32) | std::max(ind, ind_old);
        auto it = edges.find(key);
        tri_mid[(i + 2) % 3] = (it != edges.end()) ? it->second : -1;
        ind_old = ind;
      }
      int a = 0;
      switch (((tri_mid[0] >= 0) ? 1 : 0) | ((tri_mid[1] >= 0) ? 2 : 0) | ((tri_mid[2] >= 0) ? 4 : 0)) {
      case 0:
        ps_new.indices.push_back(tri);
        a = 1;
        break;
      case 1:
        ps_new.indices.push_back({tri[0], tri_mid[0], tri[2]});
        ps_new.indices.push_back({tri[2], tri_mid[0], tri[1]});
        a = 2;
        break;
      case 2:
        ps_new.indices.push_back({tri[1], tri_mid[1], tri[0]});
        ps_new.indices.push_back({tri[0], tri_mid[1], tri[2]});
        a = 2;
        break;
      case 3:
        ps_new.indices.push_back({tri[0], tri_mid[0], tri[2]});
        ps_new.indices.push_back({tri[2], tri_mid[0], tri_mid[1]});
        ps_new.indices.push_back({tri_mid[1], tri_mid[0], tri[1]});
        a = 3;
        break;
      case 4:
        ps_new.indices.push_back({tri[2], tri_mid[2], tri[1]});
        ps_new.indices.push_back({tri[1], tri_mid[2], tri[0]});
        a = 2;
        break;
      case 5:
        ps_new.indices.push_back({tri[2], tri_mid[2], tri[1]});
        ps_new.indices.push_back({tri[1], tri_mid[2], tri_mid[0]});
        ps_new.indices.push_back({tri_mid[0], tri_mid[2], tri[0]});
        a = 3;
        break;
      case 6:
        ps_new.indices.push_back({tri[1], tri_mid[1], tri[0]});
        ps_new.indices.push_back({tri[0], tri_mid[1], tri_mid[2]});
        ps_new.indices.push_back({tri_mid[2], tri_mid[1], tri[2]});
        a = 3;
        break;
      case 7:
        ps_new.indices.push_back({tri[0], tri_mid[0], tri_mid[2]});
        ps_new.indices.push_back({tri[1], tri_mid[1], tri_mid[0]});
        ps_new.indices.push_back({tri[2], tri_mid[2], tri_mid[1]});
        ps_new.indices.push_back({tri_mid[0], tri_mid[1], tri_mid[2]});
        a = 4;
        break;
      }
      for (int j = 0; j < a; j++) new_id.push_back(orig_id[i]);
    }
    if (done) {
      *ps_work = ps_new;
      orig_id = new_id;
    }

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

  if (textureprojection != PROJECTION_NONE && texturefilename.size() > 0) {
    std::shared_ptr<img_data_t> texture = load_png(this->texturefilename);
    // now apply texture to all vertices
    // create  vertex-to-triangle mapping (beta)
    std::vector<int> vert2tri;
    vert2tri.reserve(ps_work->vertices.size());
    for (int i = 0; i < ps_work->vertices.size(); i++) vert2tri.push_back(0);
    for (int i = 0; i < ps_work->indices.size(); i++) {
      for (int j = 0; j < 3; j++) {
        vert2tri[ps_work->indices[i][j]] = i;
      }
    }
    switch (textureprojection) {
    case PROJECTION_NONE: break;
    case TRIPLANAR:       {
      Vector3d vx(1, 0, 0);
      Vector3d vy(0, 1, 0);
      Vector3d vz(0, 0, 1);
      for (int i = 0; i < ps_work->vertices.size(); i++) {
        Vector3d& pt = ps_work->vertices[i];
        Vector3d& n = normals[orig_id[vert2tri[i]]];

        // triplanar texturing
        pt = pt + vx * tcoord(texture, pt[1], pt[2]) * n[0] + vy * tcoord(texture, pt[0], pt[2]) * n[1] +
             vz * tcoord(texture, pt[0], pt[1]) * n[2]

          ;
      }
    } break;
    case CUBIC:     break;
    case SPHERICAL: break;
    case CYLINDRIC: break;
    }
  }
  return std::make_unique<PolySet>(*ps_work);
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
  return createGeometry_sub(ps);
}
