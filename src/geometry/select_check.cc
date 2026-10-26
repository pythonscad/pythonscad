#include "geometry/select_check.h"

#include <algorithm>
#include <map>
#include <memory>
#include <set>
#include <vector>

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

// Test if a single point is inside a closed mesh
// Returns: 1 if inside, -1 if outside, 0 if on boundary (uncertain)
// Uses the built-in PolySet::point_inside() which is well-tested
int pointInMesh(const Vector3d& p, const PolySet& mesh)
{
  if (mesh.vertices.empty() || mesh.indices.empty()) {
    return -1;  // Empty mesh = outside
  }

  // Use the existing PolySet::point_inside() function
  bool is_inside = mesh.point_inside(p);
  return is_inside ? 1 : -1;
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

// A connected component extracted from a PolySet.
// - test_mesh holds its own local vertex array with indices remapped to be
//   local to that array, so it is a valid standalone closed mesh suitable
//   for bounding-box and point_inside() queries.
// - original_faces holds the same faces but with the original (global)
//   vertex indices into the source PolySet, so matched parts can be
//   reassembled into output geometry that shares the source vertex array.
struct Part {
  std::shared_ptr<PolySet> test_mesh;
  PolygonIndices original_faces;
};

// Helper function to separate a PolySet into connected components
// Each component is identified by grouping faces that share vertices
std::vector<Part> separateParts(const PolySet& ps)
{
  std::vector<Part> parts;

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
    PolygonIndices original_faces;  // faces with original (global) vertex indices

    // BFS to find all connected faces
    std::vector<size_t> queue;
    queue.push_back(i);
    assigned[i] = true;
    std::set<int> used_vertex_indices;  // Track which vertices are used

    while (!queue.empty()) {
      size_t curr_idx = queue.back();
      queue.pop_back();

      const auto& curr_face = ps.indices[curr_idx];
      original_faces.push_back(curr_face);

      // Track vertices used by this face
      for (int vertex_idx : curr_face) {
        used_vertex_indices.insert(vertex_idx);
      }

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

    if (original_faces.empty()) continue;

    // Build the local test mesh: extract only the vertices used by this
    // component and remap indices to be local to this component.
    auto test_mesh = std::make_shared<PolySet>(ps.getDimension());

    std::map<int, int> old_to_new_idx;  // Mapping from old vertex index to new
    int new_idx = 0;
    for (int old_idx : used_vertex_indices) {
      old_to_new_idx[old_idx] = new_idx;
      test_mesh->vertices.push_back(ps.vertices[old_idx]);
      new_idx++;
    }

    test_mesh->indices = original_faces;  // copy, then remap in place
    for (auto& face : test_mesh->indices) {
      for (int& v_idx : face) {
        v_idx = old_to_new_idx[v_idx];
      }
    }

    Part part;
    part.test_mesh = test_mesh;
    part.original_faces = std::move(original_faces);
    parts.push_back(std::move(part));
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
  matched_polysets->vertices = A.vertices;  // Use original vertex list, since
                                            // original_faces below reference
                                            // indices into this array.

  size_t matched_count = 0;

  for (size_t i = 0; i < parts.size(); ++i) {
    const auto& part = parts[i];
    const PolySet& test_mesh = *part.test_mesh;
    bool matches = false;

    switch (opt.relation) {
    case Relation::Inside:      matches = insideCheck(test_mesh, B); break;
    case Relation::NotInside:   matches = !insideCheck(test_mesh, B); break;
    case Relation::Outside:     matches = outsideCheck(test_mesh, B); break;
    case Relation::NotOutside:  matches = !outsideCheck(test_mesh, B); break;
    case Relation::Straddle:    matches = straddleCheck(test_mesh, B); break;
    case Relation::NotStraddle: matches = !straddleCheck(test_mesh, B); break;
    }

    // Collect information about this part
    PartInfo info;
    info.part_index = (int)i;
    info.matches = matches;
    result.parts.push_back(info);

    // If part matches the relation, add its faces to the result solid.
    // Use original_faces (global indices into A.vertices), not the
    // test_mesh's locally-remapped indices.
    if (matches) {
      matched_count++;
      matched_polysets->indices.insert(matched_polysets->indices.end(), part.original_faces.begin(),
                                       part.original_faces.end());
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
