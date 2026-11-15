// Portions of this file are Copyright 2021 Google LLC, and licensed under GPL2+. See COPYING.
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/exceptions.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <list>
#include <memory>
#include <vector>

#include "geometry/PolySetBuilder.h"
#include "geometry/Polygon2d.h"
#include "geometry/cgal/cgalutils.h"
#include "utils/printutils.h"

namespace CGALUtils {

namespace Polygon2DCGAL {

struct FaceInfo {
  FaceInfo() = default;
  int nesting_level{42};
  [[nodiscard]] bool in_domain() const { return nesting_level % 2 == 1; }
};

struct IdInfo {
  int id = -1;
};

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb = CGAL::Triangulation_vertex_base_with_info_2<IdInfo, K>;
using Fbb = CGAL::Triangulation_face_base_with_info_2<FaceInfo, K>;
using Fb = CGAL::Constrained_triangulation_face_base_2<K, Fbb>;
using TDS = CGAL::Triangulation_data_structure_2<Vb, Fb>;
using Itag = CGAL::Exact_predicates_tag;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>;
using Point = CDT::Point;
using Polygon_2 = CGAL::Polygon_2<K>;

void mark_domains(CDT& ct, CDT::Face_handle start, int index, std::list<CDT::Edge>& border)
{
  if (start->info().nesting_level != -1) return;

  std::list<CDT::Face_handle> queue;
  queue.push_back(start);

  while (!queue.empty()) {
    auto fh = queue.front();
    queue.pop_front();
    if (fh->info().nesting_level == -1) {
      fh->info().nesting_level = index;
      for (int i = 0; i < 3; ++i) {
        CDT::Edge e(fh, i);
        auto n = fh->neighbor(i);
        if (n->info().nesting_level == -1) {
          if (ct.is_constrained(e)) border.push_back(e);
          else queue.push_back(n);
        }
      }
    }
  }
}

// Explore set of facets connected with non constrained edges,
// and attribute to each such set a nesting level.
// We start from facets incident to the infinite vertex, with a nesting
// level of 0. Then we recursively consider the non-explored facets incident
// to constrained edges bounding the former set and increase the nesting level by 1.
// Facets in the domain are those with an odd nesting level.
void mark_domains(CDT& cdt)
{
  for (CDT::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it) {
    it->info().nesting_level = -1;
  }

  int index = 0;
  std::list<CDT::Edge> border;
  mark_domains(cdt, cdt.infinite_face(), index++, border);
  while (!border.empty()) {
    CDT::Edge e = border.front();
    border.pop_front();
    CDT::Face_handle n = e.first->neighbor(e.second);
    if (n->info().nesting_level == -1) {
      mark_domains(cdt, n, e.first->info().nesting_level + 1, border);
    }
  }
}

}  // namespace Polygon2DCGAL

template <typename SurfaceMesh>
void triangulateFaces(SurfaceMesh& mesh)
{
  CGAL::Polygon_mesh_processing::triangulate_faces(mesh);
}

template void triangulateFaces(CGAL::Surface_mesh<CGAL::Point_3<CGAL::Epick>>& polyhedron);

std::unique_ptr<PolySet> createTriangulatedPolySetFromPolygon2d(const Polygon2d& polygon2d, bool in3d)
{
  auto polyset = std::make_unique<PolySet>(2);
  polyset->setTriangular(true);

  const auto& outlines = in3d ? polygon2d.untransformedOutlines() : polygon2d.outlines();
  std::vector<Color4f> region_colors;
  for (const auto& outline : outlines) {
    if (std::find(region_colors.begin(), region_colors.end(), outline.color) == region_colors.end()) {
      region_colors.push_back(outline.color);
      polyset->colors.push_back(outline.color.unpaintedIfFullyTransparent());
    }
  }

  // Each color has its own contours and holes. Combining the constraints of
  // adjacent colored regions would turn their shared boundary into a hole.
  for (size_t color_index = 0; color_index < region_colors.size(); ++color_index) {
    Polygon2DCGAL::CDT cdt;
    size_t vertex_index = 0;
    try {
      for (const auto& outline : outlines) {
        if (outline.vertices.empty()) continue;
        Polygon2DCGAL::CDT::Vertex_handle prev;
        for (size_t i = 0; i <= outline.vertices.size(); i++) {
          const size_t idx = i % outline.vertices.size();
          const auto& v = outline.vertices[idx];
          auto curr = cdt.insert({v[0], v[1]});
          if (curr != prev) {
            if (i < outline.vertices.size()) {
              // Keep the original vertex order for extrusion callers, including
              // vertices belonging to other color regions.
              if (color_index == 0) polyset->vertices.emplace_back(v[0], v[1], 0.0);
              curr->info().id = vertex_index++;
            }
            if (prev != nullptr && outline.color == region_colors[color_index]) {
              cdt.insert_constraint(prev, curr);
            }
            prev = curr;
          }
        }
      }
    } catch (const CGAL::Precondition_exception& e) {
      LOG("CGAL error in Polygon2d::tesselate(): %1$s", e.what());
      return nullptr;
    }

    mark_domains(cdt);
    for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
      if (fit->info().in_domain()) {
        // Intersecting constraints would add vertices without source indices.
        assert(fit->vertex(0)->info().id != -1);
        assert(fit->vertex(1)->info().id != -1);
        assert(fit->vertex(2)->info().id != -1);
        polyset->indices.push_back(
          {fit->vertex(0)->info().id, fit->vertex(1)->info().id, fit->vertex(2)->info().id});
        polyset->color_indices.push_back(color_index);
      }
    }
  }
  return polyset;
}

}  // namespace CGALUtils
