// Portions of this file are Copyright 2021 Google LLC, and licensed under GPL2+. See COPYING.
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/exceptions.h>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cmath>
#include <iterator>
#include <map>
#include <unordered_map>
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
  std::unordered_map<Color4f, size_t> region_lookup;
  std::vector<std::vector<size_t>> regions;
  for (size_t i = 0; i < outlines.size(); ++i) {
    const auto [it, inserted] = region_lookup.emplace(outlines[i].color, regions.size());
    if (inserted) {
      regions.emplace_back();
      polyset->colors.push_back(outlines[i].color.unpaintedIfFullyTransparent());
    }
    regions[it->second].push_back(i);
  }

  // Keep the source vertex order independently of the per-color triangulations.
  // Repeated coordinates retain their last source index, as in the original CDT.
  std::vector<std::vector<size_t>> source_vertices(outlines.size());
  std::map<Polygon2DCGAL::Point, size_t, Polygon2DCGAL::K::Less_xy_2> source_indices;
  for (size_t i = 0; i < outlines.size(); ++i) {
    auto& indices = source_vertices[i];
    for (const auto& point : outlines[i].vertices) {
      if (!std::isfinite(point[0]) || !std::isfinite(point[1])) {
        LOG("CGAL error in Polygon2d::tesselate(): non-finite vertex");
        return nullptr;
      }
      if (!indices.empty() && polyset->vertices[indices.back()].head<2>() == point) continue;
      const auto index = polyset->vertices.size();
      polyset->vertices.emplace_back(point[0], point[1], 0.0);
      indices.push_back(index);
      if (regions.size() > 1) source_indices[{point[0], point[1]}] = index;
    }
  }

  namespace bg = boost::geometry;
  namespace bgi = boost::geometry::index;
  using SpatialPoint = bg::model::d2::point_xy<double>;
  using SpatialBox = bg::model::box<SpatialPoint>;
  using SpatialVertex = std::pair<SpatialPoint, size_t>;
  std::vector<SpatialVertex> spatial_vertices;
  spatial_vertices.reserve(source_indices.size());
  for (const auto& [point, index] : source_indices) {
    spatial_vertices.emplace_back(SpatialPoint(point.x(), point.y()), index);
  }
  const bgi::rtree<SpatialVertex, bgi::quadratic<16>> point_index(spatial_vertices);
  const auto point_at = [&](size_t index) {
    const auto& v = polyset->vertices[index];
    return Polygon2DCGAL::Point(v[0], v[1]);
  };
  const auto less = Polygon2DCGAL::K().less_xy_2_object();

  // Each color has its own contours and holes. Unrelated interior vertices must
  // not be reinserted for each color, but shared boundary split points must be.
  for (size_t color_index = 0; color_index < regions.size(); ++color_index) {
    Polygon2DCGAL::CDT cdt;
    try {
      for (const auto outline_index : regions[color_index]) {
        const auto& indices = source_vertices[outline_index];
        if (regions.size() == 1) {
          // Preserve the inexpensive original path for ordinary uncolored models.
          Polygon2DCGAL::CDT::Vertex_handle first, prev;
          for (const auto index : indices) {
            auto curr = cdt.insert(point_at(index));
            curr->info().id = static_cast<int>(index);
            if (first == nullptr) first = curr;
            if (prev != nullptr && prev != curr) cdt.insert_constraint(prev, curr);
            prev = curr;
          }
          if (prev != nullptr && prev != first) cdt.insert_constraint(prev, first);
          continue;
        }
        for (size_t i = 0; i < indices.size(); ++i) {
          const auto a = point_at(indices[i]);
          const auto b = point_at(indices[(i + 1) % indices.size()]);
          if (a == b) continue;
          const Polygon2DCGAL::K::Segment_2 segment(a, b);
          const SpatialBox box(SpatialPoint(std::min(a.x(), b.x()), std::min(a.y(), b.y())),
                               SpatialPoint(std::max(a.x(), b.x()), std::max(a.y(), b.y())));
          std::vector<SpatialVertex> boundary;
          point_index.query(bgi::intersects(box), std::back_inserter(boundary));
          boundary.erase(
            std::remove_if(boundary.begin(), boundary.end(),
                           [&](const auto& v) { return !segment.has_on(point_at(v.second)); }),
            boundary.end());
          const bool forward = less(a, b);
          std::sort(boundary.begin(), boundary.end(), [&](const auto& left, const auto& right) {
            return forward ? less(point_at(left.second), point_at(right.second))
                           : less(point_at(right.second), point_at(left.second));
          });
          Polygon2DCGAL::CDT::Vertex_handle prev;
          for (const auto& [_, index] : boundary) {
            auto curr = cdt.insert(point_at(index));
            curr->info().id = static_cast<int>(index);
            if (prev != nullptr && prev != curr) cdt.insert_constraint(prev, curr);
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
