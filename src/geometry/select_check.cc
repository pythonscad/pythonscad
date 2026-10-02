#include "geometry/select_check.h"

#include <algorithm>
#include <memory>
#include <set>
#include <vector>
#include <map>

#include "geometry/Geometry.h"
#include "geometry/PolySet.h"
#include "geometry/linalg.h"

namespace SelectCheck {

// Forward declarations of spatial test functions
namespace {

// Check if bounding boxes don't overlap
bool boundingBoxesSeparated(const PolySet& part, const PolySet& B)
{
  const auto& part_bb = part.getBoundingBox();
  const auto& B_bb = B.getBoundingBox();

  // Check each axis
  if (part_bb.min().x() > B_bb.max().x() || part_bb.max().x() < B_bb.min().x()) return true;
  if (part_bb.min().y() > B_bb.max().y() || part_bb.max().y() < B_bb.min().y()) return true;
  if (part_bb.min().z() > B_bb.max().z() || part_bb.max().z() < B_bb.min().z()) return true;

  return false;
}

// Test if a single point is inside a closed mesh using ray casting
// Returns: 1 if inside, -1 if outside, 0 if on boundary (uncertain)
int pointInMesh(const Vector3d& p, const PolySet& mesh)
{
  // Use ray casting along Z axis
  // Count intersections with mesh triangles
  // Odd count = inside, even count = outside

  int intersection_count = 0;
  const auto& vertices = mesh.vertices;
  const auto& indices = mesh.indices;  // PolygonIndices = std::vector<IndexedFace>

  // Defensive checks for empty mesh
  if (vertices.empty() || indices.empty()) {
    return -1;  // Empty mesh = outside
  }

  // Cast a ray from point in +Z direction
  Vector3d ray_start = p;
  Vector3d ray_dir = Vector3d(0, 0, 1);

  for (const auto& face : indices) {
    if (face.size() < 3) continue;

    // Triangle fan triangulation from first vertex (fan apex = vertex 0)
    // This creates proper triangles: (0,1,2), (0,2,3), (0,3,4), etc.
    int v0_idx = face[0];
    if (v0_idx < 0 || v0_idx >= (int)vertices.size()) {
      continue;  // Skip face if apex is invalid
    }

    // Test ray-triangle intersection for each triangle in fan
    for (size_t i = 1; i + 1 < face.size(); ++i) {
      int v1_idx = face[i];
      int v2_idx = face[i + 1];

      if (v1_idx < 0 || v1_idx >= (int)vertices.size() || v2_idx < 0 || v2_idx >= (int)vertices.size()) {
        continue;
      }

      const Vector3d& v0 = vertices[v0_idx];
      const Vector3d& v1 = vertices[v1_idx];
      const Vector3d& v2 = vertices[v2_idx];

      // Möller–Trumbore ray-triangle intersection
      const double EPSILON = 1e-8;
      Vector3d edge1 = v1 - v0;
      Vector3d edge2 = v2 - v0;
      Vector3d h = ray_dir.cross(edge2);
      double a = edge1.dot(h);

      if (std::abs(a) < EPSILON) continue;  // Ray parallel to triangle

      double f = 1.0 / a;
      Vector3d s = ray_start - v0;
      double u = f * s.dot(h);

      if (u < 0.0 || u > 1.0) continue;

      Vector3d q = s.cross(edge1);
      double v = f * ray_dir.dot(q);

      if (v < 0.0 || u + v > 1.0) continue;

      double t = f * edge2.dot(q);

      if (t > EPSILON) {  // Intersection in front of ray
        intersection_count++;
      }
    }
  }

  // Odd = inside, even = outside
  return (intersection_count % 2 == 1) ? 1 : -1;
}

// Check if all vertices of part are inside B
bool allVerticesInside(const PolySet& part, const PolySet& B)
{
  for (const auto& vertex : part.vertices) {
    if (pointInMesh(vertex, B) != 1) {
      return false;
    }
  }
  return true;
}

// Check if any vertex of part is inside B
bool anyVertexInside(const PolySet& part, const PolySet& B)
{
  for (const auto& vertex : part.vertices) {
    if (pointInMesh(vertex, B) == 1) {
      return true;
    }
  }
  return false;
}

// Check if part is completely inside B
bool insideCheck(const PolySet& part, const PolySet& B)
{
  // Quick check: if bounding boxes don't overlap, part can't be inside B
  if (boundingBoxesSeparated(part, B)) {
    return false;
  }

  // All vertices of part must be inside B
  return allVerticesInside(part, B);
}

// Check if part is completely outside B
bool outsideCheck(const PolySet& part, const PolySet& B)
{
  // Quick check: if bounding boxes don't overlap, part is outside B
  if (boundingBoxesSeparated(part, B)) {
    return true;
  }

  // No vertex of part should be inside B, and no vertex of B inside part
  if (anyVertexInside(part, B)) {
    return false;
  }

  // Also check if any vertex of B is inside part
  for (const auto& vertex : B.vertices) {
    if (pointInMesh(vertex, part) == 1) {
      return false;
    }
  }

  return true;
}

// Check if part straddles B (crosses its boundary)
bool straddleCheck(const PolySet& part, const PolySet& B)
{
  // Straddle = neither completely inside nor completely outside
  // i.e., some vertices inside and some outside

  bool has_inside = false;
  bool has_outside = false;

  for (const auto& vertex : part.vertices) {
    int test = pointInMesh(vertex, B);
    if (test == 1) {
      has_inside = true;
    } else if (test == -1) {
      has_outside = true;
    }
    // If we have both, we're definitely straddling
    if (has_inside && has_outside) {
      return true;
    }
  }

  return has_inside && has_outside;
}

// Helper function to separate a PolySet into connected components
// Each component is identified by grouping faces that share vertices
std::vector<std::shared_ptr<PolySet>> separateParts(const PolySet& ps)
{
  std::vector<std::shared_ptr<PolySet>> parts;

  if (ps.indices.empty()) {
    return parts;  // Empty polyhedron
  }

  // Build vertex-to-faces map for efficient lookup O(1) instead of O(n)
  std::map<int, std::vector<size_t>> vertex_to_faces;
  for (size_t face_idx = 0; face_idx < ps.indices.size(); ++face_idx) {
    for (int vertex_idx : ps.indices[face_idx]) {
      vertex_to_faces[vertex_idx].push_back(face_idx);
    }
  }

  // Track which faces have been assigned to a component
  std::vector<bool> assigned(ps.indices.size(), false);

  // Process each face as the start of a potential component
  for (size_t i = 0; i < ps.indices.size(); ++i) {
    if (assigned[i]) continue;

    // Start a new component with this face
    auto component = std::make_shared<PolySet>(ps.getDimension());
    // Deep copy vertices to ensure they outlive the original polyset
    // This prevents use-after-free if ps is temporary
    component->vertices = ps.vertices;

    // BFS to find all connected faces
    std::vector<size_t> queue;
    queue.push_back(i);
    assigned[i] = true;

    while (!queue.empty()) {
      size_t curr_idx = queue.back();
      queue.pop_back();

      const auto& curr_face = ps.indices[curr_idx];
      component->indices.push_back(curr_face);

      // Find neighboring faces that share vertices (using pre-built map)
      std::set<size_t> neighbors;  // Use set to avoid duplicates
      for (int vertex_idx : curr_face) {
        auto it = vertex_to_faces.find(vertex_idx);
        if (it != vertex_to_faces.end()) {
          for (size_t neighbor_idx : it->second) {
            if (!assigned[neighbor_idx]) {
              neighbors.insert(neighbor_idx);
            }
          }
        }
      }

      // Add neighbors to queue
      for (size_t neighbor_idx : neighbors) {
        assigned[neighbor_idx] = true;
        queue.push_back(neighbor_idx);
      }
    }

    if (!component->indices.empty()) {
      parts.push_back(component);
    }
  }

  return parts;
}

}  // namespace

// Separate A into individual parts and test each against B
Result check(const PolySet& A, const PolySet& B, const Options& opt)
{
  Result result;

  // Separate A into connected components (individual parts)
  auto parts = separateParts(A);

  if (parts.empty()) {
    // Empty input, return clean result
    return result;
  }

  // Test each part against B
  // Collect parts that MATCH the relation (not violations)
  auto matched_polysets = std::make_shared<PolySet>(A.getDimension());
  matched_polysets->vertices = A.vertices;  // Use original vertex list

  size_t matched_count = 0;

  for (size_t i = 0; i < parts.size(); ++i) {
    const auto& part = *parts[i];
    bool matches = false;

    switch (opt.relation) {
    case Relation::Inside:      matches = insideCheck(part, B); break;
    case Relation::NotInside:   matches = !insideCheck(part, B); break;
    case Relation::Outside:     matches = outsideCheck(part, B); break;
    case Relation::NotOutside:  matches = !outsideCheck(part, B); break;
    case Relation::Straddle:    matches = straddleCheck(part, B); break;
    case Relation::NotStraddle: matches = !straddleCheck(part, B); break;
    }

    // Collect information about this part
    PartInfo info;
    info.part_index = (int)i;
    info.matches = matches;
    result.parts.push_back(info);

    // If part matches the relation, add its faces to the result solid
    if (matches) {
      matched_count++;
      // Append this part's faces to the matched solid
      matched_polysets->indices.insert(matched_polysets->indices.end(), part.indices.begin(),
                                       part.indices.end());
    } else {
      // Count non-matching parts (for DRC violation reporting)
      result.count++;
    }
  }

  // Set result geometry: matched parts (not violations)
  if (matched_count > 0) {
    result.error_solid = matched_polysets;  // Named 'error_solid' for DRC compatibility,
                                            // but contains matched parts for normal use
  }

  return result;
}

}  // namespace SelectCheck
