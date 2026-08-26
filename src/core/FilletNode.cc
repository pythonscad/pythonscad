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

#include "FilletNode.h"
#include "module.h"
#include "ModuleInstantiation.h"
#include "Children.h"
#include "Parameters.h"
#include "src/utils/printutils.h"
#include "core/FilletDiagnostics.h"
#include "io/fileutils.h"
#include "Builtins.h"
#include "handle_dep.h"
#include "src/geometry/PolySetBuilder.h"
#include "geometry/ClipperUtils.h"

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <iomanip>
#include <iterator>
#include <map>
#include <set>
#include <sstream>

#include <src/geometry/PolySetUtils.h>
#include <src/core/Tree.h>
#include <src/geometry/GeometryEvaluator.h>
#include <boost/functional/hash.hpp>
#include <src/utils/hash.h>
#include <src/geometry/PolySetUtils.h>

struct SearchReplace {
  int pol;
  int search;
  int other = -1;
  IndexedFace replace;
};

typedef std::vector<int> intList;

bool list_included(const std::vector<int>& list, int needle)
{
  for (size_t i = 0; i < list.size(); i++) {
    if (list[i] == needle) return true;
  }
  return false;
}
std::shared_ptr<const PolySet> childToPolySet(std::shared_ptr<AbstractNode> child)
{
  Tree tree(child, "");
  GeometryEvaluator geomevaluator(tree);
  std::shared_ptr<const Geometry> geom = geomevaluator.evaluateGeometry(*tree.root(), true);
  std::shared_ptr<const PolySet> ps;
  return PolySetUtils::getGeometryAsPolySet(geom);
}

// Credit: inphase Ryan Colyer
Vector3d Bezier(double t, Vector3d a, Vector3d b, Vector3d c)
{
  return (a * (1 - t) + b * t) * (1 - t) + (b * (1 - t) + c * t) * t;  // TODO improve
}

void bezier_patch(PolySetBuilder& builder, Vector3d center, Vector3d dir[3], int concave_1,
                  int concave_2, int concave_3, int N)
{
  if ((dir[1].cross(dir[0])).dot(dir[2]) < 0) {
    Vector3d tmp = dir[0];
    dir[0] = dir[1];
    dir[1] = tmp;
  }
  Vector3d xdir = dir[0].normalized();
  Vector3d ydir = dir[1].normalized();
  Vector3d zdir = dir[2].normalized();

  // zdir shall look upwards

  Matrix3d mat;
  mat << xdir[0], ydir[0], zdir[0], xdir[1], ydir[1], zdir[1], xdir[2], ydir[2], zdir[2];

  xdir = Vector3d(1, 0, 0) * dir[0].norm();
  ydir = Vector3d(0, 1, 0) * dir[1].norm();
  zdir = Vector3d(0, 0, 1) * dir[2].norm();

  // now use matrices to transform the vectors into std orientation
  //
  //  N = floor(N/2)*2 + 1;
  Vector3d pt;
  std::vector<Vector3d> points_xz;
  std::vector<Vector3d> points_yz;
  for (int i = 0; i < N; i++) {
    double t = (double)i / (double)(N - 1);
    points_xz.push_back(
      Bezier(t, xdir, xdir + zdir, zdir + 2 * (concave_1 + concave_2) * (xdir + ydir)));
    points_yz.push_back(
      Bezier(t, ydir, ydir + zdir, zdir + 2 * (concave_1 + concave_2) * (xdir + ydir)));
  }

  std::vector<int> points;
  for (int i = 0; i < N; i++) {
    if (i == N - 1) {
      pt = zdir + 2 * (concave_1 + concave_2) * (xdir + ydir);
      pt = mat * pt;
      points.push_back(builder.vertexIndex(pt + center));
    } else {
      int M = N - i;
      for (int j = 0; j < M; j++) {
        int k;
        if (concave_1 == 1 || concave_3 == 1) k = j;
        else k = M - 1 - j;
        double t2 = (double)k / (double)(M - 1);
        pt = Bezier(t2, points_xz[i], Vector3d(points_xz[i][0], points_yz[i][1], points_xz[i][2]),
                    points_yz[i]);
        pt = mat * pt;
        points.push_back(builder.vertexIndex(center + pt));
      }
    }
  }
  // total points = N*(N-1)/2
  int off = 0;
  for (int i = 0; i < N - 1; i++) {  // Zeile i, i-1
    int off_new = off + (N - i);
    for (int j = 0; j < N - i - 1; j++) {
      builder.appendPolygon({points[off + j], points[off + j + 1], points[off_new + j]});
      if (j < N - i - 2) {
        builder.appendPolygon({points[off + j + 1], points[off_new + j + 1], points[off_new + j]});
      }
    }
    off = off_new;
  }
}

void debug_pt(const char *msg, Vector3d pt)
{
  printf("%s %g/%g/%g\n", msg, pt[0], pt[1], pt[2]);
}

namespace {

struct FilletEdgePair {
  int facea = -1;
  int faceb = -1;
};

// Reject an edge collapse if it leaves an effectively zero-area face.
constexpr double COLLAPSED_FACE_MIN_RELATIVE_AREA = 1e-12;

bool validateFilletEdgePairs(const std::vector<IndexedFace>& indices, EdgeKey& failedEdge)
{
  std::unordered_map<EdgeKey, FilletEdgePair, boost::hash<EdgeKey>> edge_db;

  for (size_t faceIndex = 0; faceIndex < indices.size(); faceIndex++) {
    const auto& face = indices[faceIndex];
    int n = face.size();
    for (int pos = 0; pos < n; pos++) {
      int ind1 = face[pos];
      int ind2 = face[(pos + 1) % n];
      if (ind1 == ind2) continue;

      EdgeKey edge(ind1, ind2);
      auto edgeIt = edge_db.emplace(edge, FilletEdgePair{}).first;
      FilletEdgePair& value = edgeIt->second;
      if (ind2 > ind1) {
        if (value.facea != -1) {
          failedEdge = edge;
          return false;
        }
        value.facea = faceIndex;
      } else {
        if (value.faceb != -1) {
          failedEdge = edge;
          return false;
        }
        value.faceb = faceIndex;
      }
    }
  }

  for (const auto& edge : edge_db) {
    if (edge.second.facea == -1 || edge.second.faceb == -1) {
      failedEdge = edge.first;
      return false;
    }
  }
  return true;
}

bool validateCollapsedMesh(const std::vector<IndexedFace>& indices,
                           const std::vector<Vector3d>& vertices,
                           const std::vector<Vector4d>& referenceNormals)
{
  EdgeKey failedEdge;
  if (!validateFilletEdgePairs(indices, failedEdge) || indices.size() != referenceNormals.size()) {
    return false;
  }

  for (size_t faceIndex = 0; faceIndex < indices.size(); faceIndex++) {
    const auto& face = indices[faceIndex];
    if (face.size() < 3) return false;
    for (size_t i = 0; i < face.size(); i++) {
      for (size_t j = i + 1; j < face.size(); j++) {
        if (face[i] == face[j]) return false;
      }
    }

    Vector3d area = Vector3d::Zero();
    double maxEdgeSquared = 0;
    for (size_t i = 0; i < face.size(); i++) {
      const Vector3d& current = vertices[face[i]];
      const Vector3d& next = vertices[face[(i + 1) % face.size()]];
      area += current.cross(next);
      maxEdgeSquared = std::max(maxEdgeSquared, (next - current).squaredNorm());
    }
    if (!area.allFinite() || maxEdgeSquared == 0 ||
        area.squaredNorm() <= maxEdgeSquared * maxEdgeSquared * COLLAPSED_FACE_MIN_RELATIVE_AREA *
                                COLLAPSED_FACE_MIN_RELATIVE_AREA ||
        area.dot(referenceNormals[faceIndex].head<3>()) <= 0) {
      return false;
    }
  }
  return true;
}

using IntPoint = Clipper2Lib::Point64;
using IntPath = Clipper2Lib::Path64;
using WideInt = __int128;

struct FaceProjector {
  int axis1;
  int axis2;
  double scale;

  IntPoint operator()(const Vector3d& point) const
  {
    return {static_cast<int64_t>(std::llround(point[axis1] * scale)),
            static_cast<int64_t>(std::llround(point[axis2] * scale))};
  }
};

FaceProjector makeFaceProjector(const Vector3d& normal, const BoundingBox& bounds)
{
  int drop = 0;
  if (std::abs(normal[1]) > std::abs(normal[drop])) drop = 1;
  if (std::abs(normal[2]) > std::abs(normal[drop])) drop = 2;
  const int axis1 = (drop + 1) % 3;
  const int axis2 = (drop + 2) % 3;
  const int scaleBits = ClipperUtils::scaleBitsFromBounds(bounds, 50);
  return {axis1, axis2, std::ldexp(1.0, scaleBits)};
}

bool samePoint(const IntPoint& a, const IntPoint& b)
{
  return a.x == b.x && a.y == b.y;
}

WideInt orient2d(const IntPoint& a, const IntPoint& b, const IntPoint& c)
{
  return static_cast<WideInt>(b.x - a.x) * static_cast<WideInt>(c.y - a.y) -
         static_cast<WideInt>(b.y - a.y) * static_cast<WideInt>(c.x - a.x);
}

bool pointOnSegment(const IntPoint& point, const IntPoint& a, const IntPoint& b)
{
  return orient2d(a, b, point) == 0 && point.x >= std::min(a.x, b.x) && point.x <= std::max(a.x, b.x) &&
         point.y >= std::min(a.y, b.y) && point.y <= std::max(a.y, b.y);
}

bool segmentsIntersect(const IntPoint& a, const IntPoint& b, const IntPoint& c, const IntPoint& d)
{
  const WideInt ab_c = orient2d(a, b, c);
  const WideInt ab_d = orient2d(a, b, d);
  const WideInt cd_a = orient2d(c, d, a);
  const WideInt cd_b = orient2d(c, d, b);
  if (((ab_c > 0 && ab_d < 0) || (ab_c < 0 && ab_d > 0)) &&
      ((cd_a > 0 && cd_b < 0) || (cd_a < 0 && cd_b > 0))) {
    return true;
  }
  return (ab_c == 0 && pointOnSegment(c, a, b)) || (ab_d == 0 && pointOnSegment(d, a, b)) ||
         (cd_a == 0 && pointOnSegment(a, c, d)) || (cd_b == 0 && pointOnSegment(b, c, d));
}

IntPath projectFace(const IndexedFace& face, const std::vector<Vector3d>& vertices,
                    const FaceProjector& projector)
{
  IntPath result;
  result.reserve(face.size());
  for (int index : face) {
    const IntPoint point = projector(vertices[index]);
    if (result.empty() || !samePoint(result.back(), point)) result.push_back(point);
  }
  if (result.size() > 1 && samePoint(result.front(), result.back())) result.pop_back();
  return result;
}

WideInt pathArea(const IntPath& path)
{
  WideInt area = 0;
  for (size_t i = 0; i < path.size(); i++) {
    const auto& a = path[i];
    const auto& b = path[(i + 1) % path.size()];
    area += static_cast<WideInt>(a.x) * b.y - static_cast<WideInt>(a.y) * b.x;
  }
  return area;
}

bool pathSelfIntersects(const IntPath& path)
{
  const size_t count = path.size();
  for (size_t i = 0; i < count; i++) {
    const size_t inext = (i + 1) % count;
    for (size_t j = i + 1; j < count; j++) {
      const size_t jnext = (j + 1) % count;
      if (i == j || inext == j || jnext == i) continue;
      if (segmentsIntersect(path[i], path[inext], path[j], path[jnext])) return true;
    }
  }
  return false;
}

bool pathsIntersect(const IntPath& first, const IntPath& second)
{
  for (size_t i = 0; i < first.size(); i++) {
    for (size_t j = 0; j < second.size(); j++) {
      if (segmentsIntersect(first[i], first[(i + 1) % first.size()], second[j],
                            second[(j + 1) % second.size()])) {
        return true;
      }
    }
  }
  return false;
}

enum class PointLocation { Outside, Inside, Boundary };

PointLocation pointInPath(const IntPoint& point, const IntPath& path)
{
  bool inside = false;
  for (size_t i = 0; i < path.size(); i++) {
    const IntPoint& first = path[i];
    const IntPoint& second = path[(i + 1) % path.size()];
    if (pointOnSegment(point, first, second)) return PointLocation::Boundary;
    if ((first.y > point.y) == (second.y > point.y)) continue;
    const long double intersectionX = static_cast<long double>(second.x - first.x) *
                                        static_cast<long double>(point.y - first.y) /
                                        static_cast<long double>(second.y - first.y) +
                                      first.x;
    if (intersectionX > point.x) inside = !inside;
  }
  return inside ? PointLocation::Inside : PointLocation::Outside;
}

bool validateModifiedFaces(const std::vector<IndexedFace>& original,
                           const std::vector<IndexedFace>& modified, const std::vector<int>& faceParents,
                           const std::vector<Vector4d>& faceNormals,
                           const std::vector<Vector3d>& vertices, int& failedFace, std::string& reason)
{
  BoundingBox bounds;
  for (const auto& vertex : vertices) bounds.extend(vertex);

  for (size_t outer = 0; outer < modified.size(); outer++) {
    if (faceParents[outer] != -1) continue;
    const FaceProjector projector = makeFaceProjector(faceNormals[outer].head<3>(), bounds);
    std::vector<size_t> group{outer};
    for (size_t hole = 0; hole < modified.size(); hole++) {
      if (faceParents[hole] == static_cast<int>(outer)) group.push_back(hole);
    }

    std::vector<IntPath> modifiedPaths;
    std::vector<WideInt> modifiedAreas;
    for (size_t face : group) {
      modifiedPaths.push_back(projectFace(modified[face], vertices, projector));
      modifiedAreas.push_back(pathArea(modifiedPaths.back()));
      const WideInt originalArea = pathArea(projectFace(original[face], vertices, projector));
      if (originalArea != 0 && modifiedAreas.back() != 0 &&
          (originalArea > 0) != (modifiedAreas.back() > 0)) {
        failedFace = static_cast<int>(face);
        reason = "a face contour turns inside out";
        return false;
      }
      // A fillet may consume a face exactly. The replacement strips close the
      // mesh, so a zero-area residual contour is valid and needs no tessellation.
      if (modifiedAreas.back() != 0 && pathSelfIntersects(modifiedPaths.back())) {
        failedFace = static_cast<int>(face);
        reason = "a face contour intersects itself";
        return false;
      }
    }

    if (modifiedAreas.front() == 0) {
      for (size_t i = 1; i < modifiedAreas.size(); i++) {
        if (modifiedAreas[i] == 0) continue;
        failedFace = static_cast<int>(group[i]);
        reason = "a hole remains after its outer face is consumed";
        return false;
      }
      continue;
    }

    for (size_t i = 0; i < modifiedPaths.size(); i++) {
      for (size_t j = i + 1; j < modifiedPaths.size(); j++) {
        if (modifiedAreas[i] == 0 || modifiedAreas[j] == 0) continue;
        const bool boundariesIntersect = pathsIntersect(modifiedPaths[i], modifiedPaths[j]);
        const PointLocation secondInFirst = pointInPath(modifiedPaths[j].front(), modifiedPaths[i]);
        const PointLocation firstInSecond = pointInPath(modifiedPaths[i].front(), modifiedPaths[j]);
        const bool invalidContainment =
          i == 0 ? secondInFirst != PointLocation::Inside || firstInSecond != PointLocation::Outside
                 : secondInFirst != PointLocation::Outside || firstInSecond != PointLocation::Outside;
        if (boundariesIntersect || invalidContainment) {
          failedFace = static_cast<int>(group[j]);
          reason = "rounded outer and hole contours overlap";
          return false;
        }
      }
    }
  }
  return true;
}

std::string radiusFailure(double radius, const std::string& reason)
{
  std::ostringstream message;
  message << "fillet(): radius " << std::setprecision(8) << radius << " is too large (" << reason
          << "); reduce the fillet radius";
  return message.str();
}

void reportFilletError(const std::string& message)
{
  FilletDiagnostics::setError(message);
  LOG(message_group::Error, "%1$s", message);
}

using FilletEdgeMap = std::unordered_map<EdgeKey, EdgeVal, boost::hash<EdgeKey>>;
using BoundaryKey = std::pair<int, int>;

BoundaryKey boundaryKey(int first, int second)
{
  return first < second ? BoundaryKey{first, second} : BoundaryKey{second, first};
}

const IndexedFace& cornerChain(const EdgeKey& key, const EdgeVal& edge, int vertex)
{
  return vertex == key.ind1 ? edge.bez1 : edge.bez2;
}

int cornerEndpointOnFace(const EdgeKey& key, const EdgeVal& edge, int vertex, int face)
{
  const auto& chain = cornerChain(key, edge, vertex);
  if (face == edge.facea) return chain.front();
  if (face == edge.faceb) return chain.back();
  return -1;
}

bool addBoundarySegment(std::map<int, std::set<int>>& adjacency,
                        std::map<BoundaryKey, std::pair<int, int>>& directions, int first, int second,
                        int desiredFirst, int desiredSecond)
{
  if (first == second) return true;
  const BoundaryKey key = boundaryKey(first, second);
  auto inserted = directions.emplace(key, std::pair<int, int>{desiredFirst, desiredSecond});
  if (!inserted.second && inserted.first->second != std::pair<int, int>{desiredFirst, desiredSecond}) {
    return false;
  }
  adjacency[first].insert(second);
  adjacency[second].insert(first);
  return true;
}

bool orderBoundaryLoop(const std::map<int, std::set<int>>& adjacency,
                       const std::map<BoundaryKey, std::pair<int, int>>& directions,
                       IndexedFace& boundary)
{
  if (adjacency.size() < 3) return false;
  for (const auto& vertex : adjacency) {
    if (vertex.second.size() != 2) return false;
  }

  const int start = adjacency.begin()->first;
  int previous = -1;
  int current = start;
  do {
    boundary.push_back(current);
    const auto& neighbors = adjacency.at(current);
    int next = *neighbors.begin();
    if (next == previous) next = *std::next(neighbors.begin());
    previous = current;
    current = next;
    if (boundary.size() > adjacency.size()) return false;
  } while (current != start);
  if (boundary.size() != adjacency.size()) return false;

  const auto firstDirection = directions.at(boundaryKey(boundary[0], boundary[1]));
  if (firstDirection != std::pair<int, int>{boundary[0], boundary[1]}) {
    std::reverse(boundary.begin(), boundary.end());
  }
  for (size_t i = 0; i < boundary.size(); i++) {
    const int first = boundary[i];
    const int second = boundary[(i + 1) % boundary.size()];
    if (directions.at(boundaryKey(first, second)) != std::pair<int, int>{first, second}) return false;
  }
  return true;
}

bool appendHighValenceCorner(PolySetBuilder& builder, int vertex, int resolution, double radius,
                             const std::vector<Vector3d>& originalVertices,
                             const std::vector<Vector3d>& plannedVertices,
                             const std::vector<IndexedFace>& faces, const std::vector<int>& faceParents,
                             const std::vector<int>& incidentFaces,
                             const std::vector<int>& incidentPositions,
                             const std::vector<int>& roundedNeighbors, const FilletEdgeMap& edgeDb,
                             std::string& reason)
{
  std::map<int, std::set<int>> adjacency;
  std::map<BoundaryKey, std::pair<int, int>> directions;
  struct OpenBoundary {
    int endpoint;
    int desiredFirst;
    int desiredSecond;
  };
  std::vector<OpenBoundary> openBoundaries;
  for (int other : roundedNeighbors) {
    const EdgeKey key(vertex, other);
    const auto edgeIt = edgeDb.find(key);
    if (edgeIt == edgeDb.end() || edgeIt->second.sel != 1) {
      reason = "the selected-edge fan is incomplete";
      return false;
    }
    const auto& chain = cornerChain(key, edgeIt->second, vertex);
    if (chain.size() < 2) {
      reason = "a rounded edge has no usable corner boundary";
      return false;
    }
    for (size_t i = 0; i + 1 < chain.size(); i++) {
      const bool atStart = vertex == key.ind1;
      const int desiredFirst = atStart ? chain[i + 1] : chain[i];
      const int desiredSecond = atStart ? chain[i] : chain[i + 1];
      if (!addBoundarySegment(adjacency, directions, chain[i], chain[i + 1], desiredFirst,
                              desiredSecond)) {
        reason = "corner boundary directions disagree";
        return false;
      }
    }
  }

  for (size_t i = 0; i < incidentFaces.size(); i++) {
    const int faceIndex = incidentFaces[i];
    const auto& face = faces[faceIndex];
    const int position = incidentPositions[i];
    const int faceSize = static_cast<int>(face.size());
    const int previous = face[(position + faceSize - 1) % faceSize];
    const int next = face[(position + 1) % faceSize];
    const EdgeKey incomingKey(previous, vertex);
    const EdgeKey outgoingKey(vertex, next);
    const auto incoming = edgeDb.find(incomingKey);
    const auto outgoing = edgeDb.find(outgoingKey);
    const bool incomingSelected = incoming != edgeDb.end() && incoming->second.sel == 1;
    const bool outgoingSelected = outgoing != edgeDb.end() && outgoing->second.sel == 1;
    if (!incomingSelected && !outgoingSelected) continue;

    const int incomingPoint =
      incomingSelected ? cornerEndpointOnFace(incomingKey, incoming->second, vertex, faceIndex) : -1;
    const int outgoingPoint =
      outgoingSelected ? cornerEndpointOnFace(outgoingKey, outgoing->second, vertex, faceIndex) : -1;
    if ((incomingSelected && incomingPoint < 0) || (outgoingSelected && outgoingPoint < 0)) {
      reason = "face and rounded-edge boundaries do not connect";
      return false;
    }
    if (incomingSelected && outgoingSelected) {
      if (!addBoundarySegment(adjacency, directions, incomingPoint, outgoingPoint, outgoingPoint,
                              incomingPoint)) {
        reason = "face and rounded-edge boundary directions disagree";
        return false;
      }
    } else if (incomingSelected) {
      openBoundaries.push_back({incomingPoint, vertex, incomingPoint});
    } else {
      openBoundaries.push_back({outgoingPoint, outgoingPoint, vertex});
    }
  }

  const Vector3d corner = originalVertices[vertex];

  std::set<int> unvisited;
  for (const auto& entry : adjacency) unvisited.insert(entry.first);
  while (!unvisited.empty()) {
    std::set<int> component;
    std::vector<int> pending{*unvisited.begin()};
    while (!pending.empty()) {
      const int current = pending.back();
      pending.pop_back();
      if (!component.insert(current).second) continue;
      unvisited.erase(current);
      for (int neighbor : adjacency.at(current)) pending.push_back(neighbor);
    }

    std::map<int, std::set<int>> componentAdjacency;
    std::map<BoundaryKey, std::pair<int, int>> componentDirections;
    for (int current : component) {
      for (int neighbor : adjacency.at(current)) {
        if (current > neighbor || component.count(neighbor) == 0) continue;
        const BoundaryKey key = boundaryKey(current, neighbor);
        const auto direction = directions.at(key);
        if (!addBoundarySegment(componentAdjacency, componentDirections, current, neighbor,
                                direction.first, direction.second)) {
          reason = "corner boundary directions disagree";
          return false;
        }
      }
    }

    size_t openCount = 0;
    for (const auto& open : openBoundaries) {
      if (component.count(open.endpoint) == 0) continue;
      openCount++;
      if (!addBoundarySegment(componentAdjacency, componentDirections, open.endpoint, vertex,
                              open.desiredFirst, open.desiredSecond)) {
        reason = "open corner-fan boundary directions disagree";
        return false;
      }
    }
    if (openCount != 0 && openCount != 2) {
      reason = "a selected-edge fan does not have exactly two ends";
      return false;
    }

    IndexedFace boundary;
    if (!orderBoundaryLoop(componentAdjacency, componentDirections, boundary)) {
      reason = "the corner opening is not one simple closed loop";
      return false;
    }

    const auto originalVertex = std::find(boundary.begin(), boundary.end(), vertex);
    if (originalVertex != boundary.end()) {
      std::rotate(boundary.begin(), originalVertex, boundary.end());
      for (size_t i = 1; i + 1 < boundary.size(); i++) {
        builder.appendPolygon({vertex, boundary[i], boundary[i + 1]});
      }
      continue;
    }

    Eigen::MatrixXd planes(incidentFaces.size(), 3);
    Eigen::VectorXd offsets(incidentFaces.size());
    for (size_t i = 0; i < incidentFaces.size(); i++) {
      Vector3d normal = calcTriangleNormal(originalVertices, faces[incidentFaces[i]]).head<3>();
      if (faceParents[incidentFaces[i]] != -1) normal = -normal;
      normal.normalize();
      planes.row(i) = normal;
      offsets[i] = normal.dot(corner) - std::abs(radius);
    }
    const Vector3d ballCenter = planes.colPivHouseholderQr().solve(offsets);
    if (!ballCenter.allFinite()) {
      reason = "the incident face planes do not define a finite corner";
      return false;
    }
    const Vector3d centerDirection = corner - ballCenter;
    if (centerDirection.norm() < 1e-12) {
      reason = "the incident face planes define a degenerate corner";
      return false;
    }
    const Vector3d patchCenter = ballCenter + centerDirection.normalized() * std::abs(radius);
    if (!patchCenter.allFinite()) {
      reason = "the corner patch center is not finite";
      return false;
    }

    const int rings = std::max(2, resolution);
    IndexedFace previousRing = boundary;
    for (int ring = 1; ring < rings - 1; ring++) {
      const double t = static_cast<double>(ring) / static_cast<double>(rings - 1);
      IndexedFace currentRing;
      currentRing.reserve(boundary.size());
      for (int boundaryIndex : boundary) {
        const Vector3d point = plannedVertices[boundaryIndex];
        Vector3d tangent = corner - point;
        if (tangent.norm() < 1e-12) tangent = patchCenter - point;
        const double controlLength = std::min(std::abs(radius), (patchCenter - point).norm());
        const Vector3d control = point + tangent.normalized() * controlLength;
        const Vector3d innerPoint = Bezier(t, point, control, patchCenter);
        if (!innerPoint.allFinite()) {
          reason = "the corner patch contains a non-finite point";
          return false;
        }
        currentRing.push_back(builder.vertexIndex(innerPoint));
      }
      for (size_t i = 0; i < boundary.size(); i++) {
        const size_t next = (i + 1) % boundary.size();
        builder.appendPolygon({previousRing[i], previousRing[next], currentRing[next], currentRing[i]});
      }
      previousRing = std::move(currentRing);
    }

    const int centerIndex = builder.vertexIndex(patchCenter);
    for (size_t i = 0; i < boundary.size(); i++) {
      const size_t next = (i + 1) % boundary.size();
      builder.appendPolygon({previousRing[i], previousRing[next], centerIndex});
    }
  }
  return true;
}

}  // namespace

std::unique_ptr<const Geometry> createFilletInt(std::shared_ptr<const PolySet> ps,
                                                std::vector<bool> corner_selected, double r_, int bn,
                                                double minang, double min_edge_len)
{
  double cos_minang = cos(minang * 3.1415 / 180.0);
  std::vector<Vector4d> normals, newnormals;
  std::vector<int> faceParents;
  normals = calcTriangleNormals(ps->vertices, ps->indices);
  std::vector<IndexedFace> merged =
    mergeTriangles(ps->indices, normals, newnormals, faceParents, ps->vertices);

  EdgeKey failedEdge;
  if (!validateFilletEdgePairs(merged, failedEdge)) {
    std::ostringstream message;
    message << "fillet() cannot process the selected object: edge " << failedEdge.ind1 << "-"
            << failedEdge.ind2
            << " is not shared by exactly two oppositely-oriented faces. This can happen when the "
               "fillet radius collides with nearby geometry; try reducing the fillet radius or applying "
               "fillet() before unioning adjacent solids.";
    reportFilletError(message.str());
    return PolySet::createEmpty();
  }

  if (bn < 2) bn = 2;
  // Create vertex2face db
  auto vertices_copy = ps->vertices;

  bool improved = false;
  std::unordered_map<EdgeKey, EdgeVal, boost::hash<EdgeKey>> edge_db;
  std::vector<intList> polinds, polposs;

  std::vector<std::vector<int>> corner_rounds;
  do {
    improved = false;  // fix short edges until happy
    std::vector<int> lockouts;

    polinds.clear();
    polposs.clear();
    intList empty;
    for (size_t i = 0; i < vertices_copy.size(); i++) {
      polinds.push_back(empty);
      polposs.push_back(empty);
    }
    for (size_t i = 0; i < merged.size(); i++) {
      for (size_t j = 0; j < merged[i].size(); j++) {
        int ind = merged[i][j];
        polinds[ind].push_back(i);
        polposs[ind].push_back(j);
      }
    }

    // create Edge DB
    int error;
    edge_db = createEdgeDb(merged, error);
    if (error)
      LOG(message_group::Warning,
          "Resulting fillet is not manifold anymore, further processing might be inaccurate");

    // which rounded edges in a corner coner_rounds[vert]=[other_verts]
    corner_rounds.clear();
    for (size_t i = 0; i < vertices_copy.size(); i++) corner_rounds.push_back(empty);

    for (auto& e : edge_db) {
      if (corner_selected[e.first.ind1] && corner_selected[e.first.ind2]) {
        assert(e.second.facea >= 0);
        assert(e.second.faceb >= 0);
        auto& facea = merged[e.second.facea];
        auto& faceb = merged[e.second.faceb];
        Vector3d fan = calcTriangleNormal(vertices_copy, facea).head<3>();
        Vector3d fbn = calcTriangleNormal(vertices_copy, faceb).head<3>();
        double d = fan.dot(fbn);
        e.second.sel = 0;
        if (d >= cos_minang) continue;  // dont create facets when the angle conner is too small

        e.second.sel = 1;
        corner_rounds[e.first.ind1].push_back(e.first.ind2);
        corner_rounds[e.first.ind2].push_back(e.first.ind1);
      }
    }
    /* TODO activate

        // eliminate  too short edges by extrapolating the neighboring edges
        for(auto &e: edge_db) {
          if(!e.second.sel) continue;
          Vector3d line = vertices_copy[e.first.ind1] - vertices_copy[e.first.ind2];
          if(line.norm() < 2*r_) {

            int a_prev=-1, a_next=-1;
            int b_prev=-1, b_next=-1;

            if(std::find(lockouts.begin(), lockouts.end(), e.first.ind1) != lockouts.end()) continue;
            if(std::find(lockouts.begin(), lockouts.end(), e.first.ind2) != lockouts.end()) continue;
            auto &facea = merged[e.second.facea];
            int na=facea.size();
            for(int i=0;i<na;i++) {
              if(facea[i] == e.first.ind1){
                a_prev = facea[(i+na-1)%na];
                a_next = facea[(i+na+2)%na];
              }
            }

            auto &faceb = merged[e.second.faceb];
            int nb=faceb.size();
            for(int i=0;i<nb;i++) {
              if(faceb[i] == e.first.ind2){
                b_prev = faceb[(i+nb-1)%nb];
                b_next = faceb[(i+nb+2)%nb];
              }
            }

            if(std::find(lockouts.begin(), lockouts.end(), a_prev) != lockouts.end()) continue;
            if(std::find(lockouts.begin(), lockouts.end(), a_next) != lockouts.end()) continue;
            if(std::find(lockouts.begin(), lockouts.end(), b_prev) != lockouts.end()) continue;
            if(std::find(lockouts.begin(), lockouts.end(), b_next) != lockouts.end()) continue;

            // is it safe to take the bigger face ?
            int commonfaceind=-1, faceind1=-1, faceind2=-1;
            EdgeKey ek1, ek2;
            if(nb > na) {
              commonfaceind=e.second.faceb; // TODO b hat die richtigen punkte
              ek1 = EdgeKey(b_prev, e.first.ind2);
              ek2 = EdgeKey(e.first.ind1, b_next);
            } else { // na  > nb)
              commonfaceind=e.second.facea; // TODO b hat die richtigen punkte
              ek1 = EdgeKey(a_prev, e.first.ind1);
              ek2 = EdgeKey(e.first.ind2, a_next);
            }

            if(edge_db.count(ek1)) {
              auto &ev1 = edge_db.at(ek1);
              if(ev1.facea == commonfaceind) faceind1= ev1.faceb;
              if(ev1.faceb == commonfaceind) faceind1= ev1.facea;
            }

            // find opposite of e.first.ind1, b_next)
            if(edge_db.count(ek2)) {
              auto &ev2 = edge_db.at(ek2);
              if(ev2.facea == commonfaceind) faceind2= ev2.faceb;
              if(ev2.faceb == commonfaceind) faceind2= ev2.facea;
            }

            Vector3d fn1 =calcTriangleNormal(vertices_copy, merged[commonfaceind]).head<3>();
            Vector3d fn2 =calcTriangleNormal(vertices_copy, merged[faceind1]).head<3>();
            Vector3d fn3 =calcTriangleNormal(vertices_copy, merged[faceind2]).head<3>();

            Vector3d fp1 = vertices_copy[merged[commonfaceind][0]];
            Vector3d fp2 = vertices_copy[merged[faceind1][0]];
            Vector3d fp3 = vertices_copy[merged[faceind2][0]];
            Vector3d ptcut;
            if(cut_face_face_face(fp1, fn1, fp2, fn2, fp3, fn3, ptcut,nullptr)) {
              printf("Error during cutting\n");
              e.second.sel=0;
              continue;
            }
            //
            // change is going to happen
            vertices_copy[e.first.ind1]=ptcut;
            lockouts.push_back(ek1.ind1);
            lockouts.push_back(ek1.ind2);
            lockouts.push_back(ek2.ind1);
            lockouts.push_back(ek2.ind2);

            for(int j=0;j<merged.size();j++) {
              auto &tri = merged[j];
              int n = tri.size();
              int dupind=-1;
              for(int i=0;i<n;i++)
              {
                if(tri[i] == e.first.ind2){
                  tri[i] =e.first.ind1;
                  if(tri[(i+1)%n] == e.first.ind1 || tri[(i+n-1)%n] == e.first.ind1) {
                    dupind=i;
                  }
                }
              }
              if(dupind != -1) {
                IndexedFace tri_new;
                for(int i=0;i<dupind;i++) tri_new.push_back(tri[i]);
                for(int i=dupind+1;i<n;i++) tri_new.push_back(tri[i]);
                tri = tri_new;
                n--;
              }
              if(n < 3) {
                    merged.erase(merged.begin()+j);
                    j--;
              }
            }
            improved=true;
          // TODO lockout
          }

        }
    */
    // Boolean operations with a tiny overlap can leave sliver edges around the
    // seam. Collapse those before constructing fillet patches.
    if (min_edge_len > 0) {
      auto candidate = edge_db.end();
      for (auto edge = edge_db.begin(); edge != edge_db.end(); ++edge) {
        const int ind1 = edge->first.ind1;
        const int ind2 = edge->first.ind2;
        const double edge_len = (vertices_copy[ind1] - vertices_copy[ind2]).norm();
        if (edge_len >= min_edge_len) continue;
        if (candidate == edge_db.end() || ind1 < candidate->first.ind1 ||
            (ind1 == candidate->first.ind1 && ind2 < candidate->first.ind2)) {
          candidate = edge;
        }
      }

      if (candidate != edge_db.end()) {
        const auto vertices_before_collapse = vertices_copy;
        const auto merged_before_collapse = merged;
        const auto face_parents_before_collapse = faceParents;
        const auto normals_before_collapse = newnormals;
        const auto selected_before_collapse = corner_selected;
        auto& e = *candidate;
        int keep = e.first.ind1;
        int drop = e.first.ind2;

        const bool keep_selected =
          keep < static_cast<int>(corner_selected.size()) && corner_selected[keep];
        const bool drop_selected =
          drop < static_cast<int>(corner_selected.size()) && corner_selected[drop];
        if (!keep_selected && drop_selected) std::swap(keep, drop);
        if (keep_selected && drop_selected) {
          vertices_copy[keep] = 0.5 * (vertices_copy[keep] + vertices_copy[drop]);
        }
        if (drop < static_cast<int>(corner_selected.size())) corner_selected[drop] = false;

        for (size_t face_index = 0; face_index < merged.size(); face_index++) {
          auto& face = merged[face_index];
          const int size = static_cast<int>(face.size());
          int duplicate = -1;
          for (int i = 0; i < size; i++) {
            if (face[i] != drop) continue;
            face[i] = keep;
            if (face[(i + 1) % size] == keep || face[(i + size - 1) % size] == keep) duplicate = i;
          }
          if (duplicate != -1) face.erase(face.begin() + duplicate);
          if (face.size() >= 3) continue;

          merged.erase(merged.begin() + static_cast<long>(face_index));
          faceParents.erase(faceParents.begin() + static_cast<long>(face_index));
          newnormals.erase(newnormals.begin() + static_cast<long>(face_index));
          for (int& parent : faceParents) {
            if (parent == static_cast<int>(face_index)) parent = -1;
            else if (parent > static_cast<int>(face_index)) parent--;
          }
          face_index--;
        }
        if (validateCollapsedMesh(merged, vertices_copy, newnormals)) {
          newnormals = calcTriangleNormals(vertices_copy, merged);
          improved = true;
        } else {
          vertices_copy = vertices_before_collapse;
          merged = merged_before_collapse;
          faceParents = face_parents_before_collapse;
          newnormals = normals_before_collapse;
          corner_selected = selected_before_collapse;
        }
      }
    }
  } while (improved == true);

  const auto closesRoundedCorner = [&](int vertex) {
    return corner_rounds[vertex].size() >= 3 && corner_rounds[vertex].size() == polinds[vertex].size();
  };

  for (const auto& edge : edge_db) {
    if (edge.second.sel != 1) continue;
    const double edgeLength = (vertices_copy[edge.first.ind2] - vertices_copy[edge.first.ind1]).norm();
    const int trimmedEnds =
      (closesRoundedCorner(edge.first.ind1) ? 1 : 0) + (closesRoundedCorner(edge.first.ind2) ? 1 : 0);
    if (trimmedEnds == 0 || edgeLength > trimmedEnds * std::abs(r_) + 1e-9) continue;

    std::ostringstream message;
    message << "fillet(): radius " << std::setprecision(8) << r_ << " is too large for edge "
            << edge.first.ind1 << "-" << edge.first.ind2 << " (length " << edgeLength
            << ", maximum radius below " << edgeLength / trimmedEnds << "); reduce the fillet radius";
    reportFilletError(message.str());
    return PolySet::createEmpty();
  }

  // start builder with existing vertices to have VertexIndex available
  //
  PolySetBuilder builder;
  for (size_t i = 0; i < vertices_copy.size(); i++) {
    builder.vertexIndex(vertices_copy[i]);  // allocate all vertices in the right order
  }

  SearchReplace s;
  std::vector<SearchReplace> sp;
  std::unordered_map<EdgeKey, IndexedFace, boost::hash<EdgeKey>> taperedMiddles;

  // plan fillets of all edges now
  for (auto& e : edge_db) {
    if (e.second.sel == 1) {
      Vector3d p1 = vertices_copy[e.first.ind1];  // both ends of the selected edge
      Vector3d p2 = vertices_copy[e.first.ind2];
      Vector3d p1org = p1, p2org = p2;
      Vector3d dir = p2 - p1;
      if (closesRoundedCorner(e.first.ind1)) p1 += dir.normalized() * r_;
      if (closesRoundedCorner(e.first.ind2)) p2 -= dir.normalized() * r_;
      dir = dir.normalized();  // TODO
      auto& facea = merged[e.second.facea];
      auto& faceb = merged[e.second.faceb];
      //

      int facean = facea.size();
      int facebn = faceb.size();
      double fanf = (faceParents[e.second.facea] != -1) ? -1 : 1;  // is the edge part of a hole
      double fbnf = (faceParents[e.second.faceb] != -1) ? -1 : 1;
      Vector3d fan = calcTriangleNormal(vertices_copy, facea).head<3>();
      Vector3d fbn = calcTriangleNormal(vertices_copy, faceb).head<3>();

      // A 1st side of the edge
      // B 2nd face of the edge
      int indposao, indposbo, indposai, indposbi;
      Vector3d unit;

      indposao = facea[(e.second.posa + facean - 1) % facean];  // o away from edge
      indposai = facea[(e.second.posa + 1) % facean];           // i on edge

      indposbo = faceb[(e.second.posb + 2) % facebn];
      indposbi = faceb[e.second.posb];

      Vector3d e_fa1 = (vertices_copy[indposao] - vertices_copy[facea[e.second.posa]]).normalized() *
                       fanf;  // Facea neben ind1
      Vector3d e_fa1p =
        (vertices_copy[indposai] - vertices_copy[facea[e.second.posa]]) * fanf;  // Face1 nahe  richtung
                                                                                 //
      Vector3d e_fb1 =
        (vertices_copy[indposbo] - vertices_copy[faceb[(e.second.posb + 1) % facebn]]).normalized() *
        fbnf;  // Faceb neben ind1
      Vector3d e_fb1p =
        (vertices_copy[indposbi] - vertices_copy[faceb[(e.second.posb + 1) % facebn]]) * fbnf;

      if (corner_rounds[e.first.ind1].size() == 2) {
        double a = (e_fb1.cross(e_fa1)).dot(dir);
        double b = (fan.cross(fbn)).dot(e_fa1p) * fanf * fbnf;
        if (list_included(corner_rounds[e.first.ind1], indposao)) {
          double ang = (dir).dot(e_fa1.normalized());
          e_fa1 += dir * fanf;
          if (a * b < 0) e_fa1 = -e_fa1 * fanf;
          e_fa1 /= sqrt(1 - ang * ang);
        }

        if (list_included(corner_rounds[e.first.ind1], indposbo)) {
          double ang = (dir).dot(e_fb1.normalized());
          e_fb1 += dir * fbnf;
          if (a * b < 0) e_fb1 = -e_fb1 * fbnf;
          e_fb1 /= sqrt(1 - ang * ang);
        }
      }

      if (closesRoundedCorner(e.first.ind1)) {
        if ((fbn.cross(fan)).dot(e_fa1p) < 0 || (fbn.cross(fan)).dot(e_fb1p) < 0) {
          if ((e_fa1p.cross(e_fa1)).dot(fan) * fanf < 0) {
            e_fa1 = -e_fa1 * fanf - 2 * (p2org - p1org).normalized();
          }
          if ((e_fb1p.cross(e_fb1)).dot(fbn) * fbnf > 0) {
            e_fb1 = -e_fb1 * fbnf - 2 * (p2org - p1org).normalized();
          }
        }
        if ((fbn.cross(fan)).dot(e_fb1p) > 0) {
          if ((e_fa1p.cross(e_fa1)).dot(fan) * fanf > 0 && (e_fa1p.cross(e_fb1)).dot(fbn) * fbnf > 0) {
            e_fb1 = -e_fb1 * fbnf - 2 * dir;
          }
          if ((e_fb1p.cross(e_fb1)).dot(fbn) * fbnf < 0 && (e_fb1p.cross(e_fa1)).dot(fan) * fanf < 0) {
            e_fa1 = -e_fa1 * fanf - 2 * dir;
          }
        }
      }
      e_fa1 *= r_;
      e_fb1 *= r_;

      indposao = facea[(e.second.posa + 2) % facean];
      indposai = facea[e.second.posa];

      indposbo = faceb[(e.second.posb + facebn - 1) % facebn];
      indposbi = faceb[(e.second.posb + 1) % facebn];

      Vector3d e_fa2 =
        (vertices_copy[indposao] - vertices_copy[facea[(e.second.posa + 1) % facean]]).normalized() *
        fanf;  // Face1 entfernte richtung
      Vector3d e_fa2p = (vertices_copy[indposai] - vertices_copy[facea[(e.second.posa + 1) % facean]]) *
                        fanf;  // Face1 entfernte richtung
                               //
      Vector3d e_fb2 =
        (vertices_copy[indposbo] - vertices_copy[faceb[(e.second.posb + 0) % facebn]]).normalized() *
        fbnf;  // Face2 entfernte Rcithung
      Vector3d e_fb2p = (vertices_copy[indposbi] - vertices_copy[faceb[(e.second.posb + 0) % facebn]]) *
                        fbnf;  // Face2 entfernte Rcithung

      //
      if (corner_rounds[e.first.ind2].size() == 2) {
        double a = (e_fb2.cross(e_fa2)).dot(dir);
        double b = (fan.cross(fbn)).dot(e_fa2p) * fanf * fbnf;
        if (list_included(corner_rounds[e.first.ind2], indposao)) {
          double ang = (dir).dot(e_fa2.normalized());
          e_fa2 -= dir * fanf;
          if (a * b > 0) e_fa2 = -e_fa2 * fanf;
          e_fa2 /= sqrt(1 - ang * ang);
        }

        if (list_included(corner_rounds[e.first.ind2], indposbo)) {
          double ang = (dir).dot(e_fb2.normalized());
          e_fb2 -= dir * fbnf;
          if (a * b > 0) e_fb2 = -e_fb2 * fbnf;
          e_fb2 /= sqrt(1 - ang * ang);
        }
      }

      if (closesRoundedCorner(e.first.ind2)) {
        if (-(fbn.cross(fan)).dot(e_fa2p) < 0 || -(fbn.cross(fan)).dot(e_fb2p) < 0) {
          if (-(e_fa2p.cross(e_fa2)).dot(fan) * fanf < 0) {
            e_fa2 = -e_fa2 * fanf + 2 * dir;
          }
          if (-(e_fb2p.cross(e_fb2)).dot(fbn) * fbnf > 0) {
            e_fb2 = -e_fb2 * fbnf + 2 * dir;
          }
        }
        if (/* -(fbn.cross(fan)).dot(e_fa2p) > 0 || */ -(fbn.cross(fan)).dot(e_fb2p) > 0) {
          if (-(e_fb2p.cross(e_fb2)).dot(fbn) * fbnf < 0 && (e_fb2p.cross(e_fa2)).dot(fan) * fanf > 0) {
            e_fa2 = -e_fa2 * fanf + 2 * dir;  // laengs links
          }
          if (-(e_fa2p.cross(e_fa2)).dot(fan) * fanf > 0 && (e_fa2p.cross(e_fb2)).dot(fbn) * fbnf < 0) {
            e_fb2 = -e_fb2 * fbnf + 2 * dir;
          }
        }
      }

      e_fa2 *= r_;
      e_fb2 *= r_;

      // Calculate bezier patches
      const bool taperStart = polinds[e.first.ind1].size() > 3 && !closesRoundedCorner(e.first.ind1);
      const bool taperEnd = polinds[e.first.ind2].size() > 3 && !closesRoundedCorner(e.first.ind2);
      IndexedFace taperedMiddle;
      if (taperStart && taperEnd) taperedMiddle.reserve(bn);
      for (int i = 0; i < bn; i++) {
        double f = (double)i / (double)(bn - 1);  // from 0 to 1
        const Vector3d roundedStart = p1 + e_fa1 - 2 * f * e_fa1 + f * f * (e_fa1 + e_fb1);
        const Vector3d roundedEnd = p2 + e_fa2 - 2 * f * e_fa2 + f * f * (e_fa2 + e_fb2);
        const Vector3d startPoint = taperStart ? vertices_copy[e.first.ind1] : roundedStart;
        const Vector3d endPoint = taperEnd ? vertices_copy[e.first.ind2] : roundedEnd;
        e.second.bez1.push_back(builder.vertexIndex(startPoint));
        e.second.bez2.push_back(builder.vertexIndex(endPoint));
        if (taperStart && taperEnd) {
          taperedMiddle.push_back(builder.vertexIndex((roundedStart + roundedEnd) / 2.0));
        }
      }
      if (!taperedMiddle.empty()) taperedMiddles.emplace(e.first, std::move(taperedMiddle));

      const auto middle = taperedMiddles.find(e.first);
      const int middleFaceA = middle == taperedMiddles.end() ? -1 : middle->second.front();
      const int middleFaceB = middle == taperedMiddles.end() ? -1 : middle->second.back();
      const auto railReplacement = [&](int faceIndex, int search, int other, int endpoint,
                                       int middlePoint) {
        IndexedFace replacement{endpoint};
        if (middlePoint < 0) return replacement;
        const auto& face = merged[faceIndex];
        const auto position = std::find(face.begin(), face.end(), search);
        if (position == face.end()) return replacement;
        const size_t index = std::distance(face.begin(), position);
        const bool outgoing = face[(index + 1) % face.size()] == other;
        if (outgoing) replacement.push_back(middlePoint);
        else replacement.insert(replacement.begin(), middlePoint);
        return replacement;
      };
      s.pol = e.second.facea;  // laengsseite1
      s.search = e.first.ind1;
      s.other = e.first.ind2;
      s.replace = railReplacement(s.pol, s.search, s.other, e.second.bez1[0], middleFaceA);
      sp.push_back(s);
      s.pol = e.second.facea;  // laengsseite1
      s.search = e.first.ind2;
      s.other = e.first.ind1;
      s.replace = railReplacement(s.pol, s.search, s.other, e.second.bez2[0], middleFaceA);
      sp.push_back(s);

      s.pol = e.second.faceb;  // laengsseite2
      s.search = e.first.ind2;
      s.other = e.first.ind1;
      s.replace = railReplacement(s.pol, s.search, s.other, e.second.bez2[bn - 1], middleFaceB);
      sp.push_back(s);
      s.pol = e.second.faceb;  // laengsseite2
      s.search = e.first.ind1;
      s.other = e.first.ind2;
      s.replace = railReplacement(s.pol, s.search, s.other, e.second.bez1[bn - 1], middleFaceB);
      sp.push_back(s);

      // stirnseite 1
      if (corner_rounds[e.first.ind1].size() == 1 && polinds[e.first.ind1].size() == 3) {
        for (size_t i = 0; i < polinds[e.first.ind1].size(); i++) {
          int faceid = polinds[e.first.ind1][i];
          if (faceid == e.second.facea) continue;
          if (faceid == e.second.faceb) continue;
          s.pol = faceid;  // stirnseite1
          s.search = e.first.ind1;
          s.other = -1;
          s.replace = {e.second.bez1};
          std::reverse(s.replace.begin(), s.replace.end());
          sp.push_back(s);
        }
      }

      // stirnseite2
      if (corner_rounds[e.first.ind2].size() == 1 && polinds[e.first.ind2].size() == 3) {
        for (size_t i = 0; i < polinds[e.first.ind2].size(); i++) {
          int faceid = polinds[e.first.ind2][i];
          if (faceid == e.second.facea) continue;
          if (faceid == e.second.faceb) continue;
          s.pol = faceid;  // stirnseite2
          s.search = e.first.ind2;
          s.other = -1;
          s.replace = {e.second.bez2};
          sp.push_back(s);
        }
      }
      //     printf("\nNum=%d\n",debug);
      //     printf("P : %g/%g/%g EA: %g/%g/%g EB %g/%g/%g\n",p1[0], p1[1], p1[2], e_fa1[0], e_fa1[1],
      //     e_fa1[2], e_fb1[0], e_fb1[1], e_fb1[2]); printf("P : %g/%g/%g EA: %g/%g/%g EB
      //     %g/%g/%g\n",p2[0], p2[1], p2[2], e_fa2[0], e_fa2[1], e_fa2[2], e_fb2[0], e_fb2[1],
      //     e_fb2[2]);
    }
  }
  // copy modified faces
  std::vector<IndexedFace> newfaces;
  for (size_t i = 0; i < merged.size(); i++) {
    const IndexedFace& face = merged[i];
    const int faceSize = static_cast<int>(face.size());
    IndexedFace newface;
    for (int position = 0; position < faceSize; position++) {
      const int vertex = face[position];
      const int previous = face[(position + faceSize - 1) % faceSize];
      const int next = face[(position + 1) % faceSize];
      IndexedFace incoming;
      IndexedFace outgoing;
      IndexedFace extra;
      for (const auto& replacement : sp) {
        if (replacement.pol != static_cast<int>(i) || replacement.search != vertex) continue;
        if (replacement.other == previous) {
          incoming.insert(incoming.end(), replacement.replace.begin(), replacement.replace.end());
        } else if (replacement.other == next) {
          outgoing.insert(outgoing.end(), replacement.replace.begin(), replacement.replace.end());
        } else {
          extra.insert(extra.end(), replacement.replace.begin(), replacement.replace.end());
        }
      }
      if (incoming.empty() && outgoing.empty() && extra.empty()) {
        newface.push_back(vertex);
        continue;
      }

      IndexedFace points;
      const bool highValenceOpenFan =
        polinds[vertex].size() > 3 && extra.empty() && incoming.empty() != outgoing.empty();
      if (highValenceOpenFan && outgoing.empty()) {
        points.insert(points.end(), incoming.begin(), incoming.end());
        points.push_back(vertex);
      } else if (highValenceOpenFan) {
        points.push_back(vertex);
        points.insert(points.end(), outgoing.begin(), outgoing.end());
      } else {
        points.insert(points.end(), incoming.begin(), incoming.end());
        points.insert(points.end(), extra.begin(), extra.end());
        points.insert(points.end(), outgoing.begin(), outgoing.end());
      }
      for (int point : points) {
        if (newface.empty() || newface.back() != point) newface.push_back(point);
      }
    }
    if (newface.size() > 1 && newface.front() == newface.back()) newface.pop_back();
    newfaces.push_back(newface);
  }
  std::vector<Vector3d> vertices;
  builder.copyVertices(vertices);

  int failedFace = -1;
  std::string failureReason;
  if (!validateModifiedFaces(merged, newfaces, faceParents, newnormals, vertices, failedFace,
                             failureReason)) {
    reportFilletError(radiusFailure(r_, failureReason));
    return PolySet::createEmpty();
  }

  std::vector<Vector3f> verticesFloat;
  for (const auto& v : vertices) verticesFloat.push_back(v.cast<float>());

  for (size_t i = 0; i < newfaces.size(); i++) {
    // tessellate first with holes // search all holes
    if (faceParents[i] != -1) continue;
    std::vector<IndexedFace> faces;
    faces.push_back(newfaces[i]);
    for (size_t j = 0; j < newfaces.size(); j++)
      if ((size_t)faceParents[j] == i) faces.push_back(newfaces[j]);
    //    if(faces.size() >1 ) continue;
    std::vector<IndexedTriangle> triangles;
    Vector3f norm(newnormals[i][0], newnormals[i][1], newnormals[i][2]);
    GeometryUtils::tessellatePolygonWithHoles(verticesFloat, faces, triangles, &norm);
    for (const auto& t : triangles) {
      builder.appendPolygon({t[0], t[1], t[2]});
    }
  }

  // add Rounded edges
  const auto appendStripSection = [&](const IndexedFace& first, const IndexedFace& second) {
    for (int i = 0; i < bn - 1; i++) {
      const int start = first[i];
      const int startNext = first[i + 1];
      const int endNext = second[i + 1];
      const int end = second[i];
      if (start == startNext && end == endNext) continue;
      if (start == startNext) {
        builder.appendPolygon({start, endNext, end});
      } else if (end == endNext) {
        builder.appendPolygon({start, startNext, end});
      } else {
        builder.appendPolygon({start, startNext, endNext, end});
      }
    }
  };
  for (auto& e : edge_db) {
    if (e.second.sel == 1) {
      const auto middle = taperedMiddles.find(e.first);
      if (middle == taperedMiddles.end()) {
        appendStripSection(e.second.bez1, e.second.bez2);
      } else {
        appendStripSection(e.second.bez1, middle->second);
        appendStripSection(middle->second, e.second.bez2);
      }
    }
  }
  // add missing 3 corner patches
  //
  for (size_t i = 0; i < vertices_copy.size(); i++) {
    if (polinds[i].size() > 3 && corner_rounds[i].size() == polinds[i].size()) {
      std::string reason;
      if (!appendHighValenceCorner(builder, static_cast<int>(i), bn, r_, vertices_copy, vertices, merged,
                                   faceParents, polinds[i], polposs[i], corner_rounds[i], edge_db,
                                   reason)) {
        std::ostringstream message;
        message << "fillet(): cannot build a valid radius " << std::setprecision(8) << r_
                << " fillet at vertex " << i << " with " << corner_rounds[i].size() << " rounded edges ("
                << reason << ")";
        reportFilletError(message.str());
        return PolySet::createEmpty();
      }
    } else if (corner_rounds[i].size() == 3) {
      // now get the right ordering of corner_rounds[i]
      IndexedFace face[3];
      Vector3d facenorm[3];
      for (int j = 0; j < 3; j++) {
        face[j] = merged[polinds[i][j]];
        facenorm[j] = calcTriangleNormal(vertices_copy, face[j]).head<3>();
        if (faceParents[polinds[i][j]] != -1) facenorm[j] = -facenorm[j];
      }

      int facebeg[3];
      int faceend[3];
      for (int j = 0; j < 3; j++) {
        facebeg[j] = face[j][(polposs[i][j] + face[j].size() - 1) % face[j].size()];
        faceend[j] = face[j][(polposs[i][j] + 1) % face[j].size()];
      }

      std::vector<int> angle;
      std::vector<Vector3d> dir;
      Vector3d x;
      if (faceend[0] == facebeg[1]) {  // 0,1,2
        x = vertices_copy[faceend[1]] - vertices_copy[i];
        dir.push_back(x.normalized() * r_);
        angle.push_back(
          (facenorm[1].cross(facenorm[2])).dot(vertices_copy[faceend[1]] - vertices_copy[i]) > 0 ? 1
                                                                                                 : -1);
        x = vertices_copy[faceend[2]] - vertices_copy[i];
        dir.push_back(x.normalized() * r_);
        angle.push_back(
          (facenorm[2].cross(facenorm[0])).dot(vertices_copy[faceend[2]] - vertices_copy[i]) > 0 ? 1
                                                                                                 : -1);
        x = vertices_copy[faceend[0]] - vertices_copy[i];
        dir.push_back(x.normalized() * r_);
        angle.push_back(
          (facenorm[0].cross(facenorm[1])).dot(vertices_copy[faceend[0]] - vertices_copy[i]) > 0 ? 1
                                                                                                 : -1);
      } else if (faceend[0] == facebeg[2]) {
        x = vertices_copy[faceend[2]] - vertices_copy[i];
        dir.push_back(x.normalized() * r_);
        angle.push_back(
          (facenorm[2].cross(facenorm[1])).dot(vertices_copy[faceend[2]] - vertices_copy[i]) > 0 ? 1
                                                                                                 : -1);
        x = vertices_copy[faceend[0]] - vertices_copy[i];
        dir.push_back(x.normalized() * r_);
        angle.push_back(
          (facenorm[0].cross(facenorm[2])).dot(vertices_copy[faceend[0]] - vertices_copy[i]) > 0 ? 1
                                                                                                 : -1);
        x = vertices_copy[faceend[1]] - vertices_copy[i];
        dir.push_back(x.normalized() * r_);
        angle.push_back(
          (facenorm[1].cross(facenorm[0])).dot(vertices_copy[faceend[1]] - vertices_copy[i]) > 0 ? 1
                                                                                                 : -1);
      } else assert(0);
      int conc1 = -1, conc2 = -1, conc3 = -1;
      int dirshift = -1;
      Vector3d pdir[3];
      if (angle[0] == -1 && angle[1] == -1 && angle[2] == -1) {
        dirshift = 0;
        conc1 = 0;
        conc2 = 0;
        conc3 = 1;
      }
      if (angle[0] == -1 && angle[1] == -1 && angle[2] == 1) {
        dirshift = 0;
        conc1 = 0;
        conc2 = 1;
        conc3 = 0;
      }
      if (angle[0] == -1 && angle[1] == 1 && angle[2] == -1) {
        dirshift = 2;
        conc1 = 0;
        conc2 = 1;
        conc3 = 0;
      }
      if (angle[0] == -1 && angle[1] == 1 && angle[2] == 1) {
        dirshift = 1;
        conc1 = 1;
        conc2 = 0;
        conc3 = 0;
      }
      if (angle[0] == 1 && angle[1] == -1 && angle[2] == -1) {
        dirshift = 1;
        conc1 = 0;
        conc2 = 1;
        conc3 = 0;
      }
      if (angle[0] == 1 && angle[1] == -1 && angle[2] == 1) {
        dirshift = 2;
        conc1 = 1;
        conc2 = 0;
        conc3 = 0;
      }
      if (angle[0] == 1 && angle[1] == 1 && angle[2] == -1) {
        dirshift = 0;
        conc1 = 1;
        conc2 = 0;
        conc3 = 0;
      }
      if (angle[0] == 1 && angle[1] == 1 && angle[2] == 1) {
        dirshift = 0;
        conc1 = 0;
        conc2 = 0;
        conc3 = 0;
      }
      if (dirshift != -1) {
        for (int i = 0; i < 3; i++) {
          pdir[i] = -dir[(i + dirshift) % 3];
        }
        bezier_patch(builder, ps->vertices[i] - pdir[0] - pdir[1] - pdir[2], pdir, conc1, conc2, conc3,
                     bn);
      }
    }
  }
  //
  auto result = builder.build();

  return result;
}

std::unique_ptr<const Geometry> FilletNode::createGeometry() const
{
  std::shared_ptr<const PolySet> ps;
  std::shared_ptr<PolySet> ps_merged;
  std::vector<bool> corner_selected;
  if (this->children.size() >= 1) {
    ps = childToPolySet(this->children[0]);
    if (ps == nullptr) return std::unique_ptr<PolySet>();

    ps_merged = std::make_shared<PolySet>(*ps);

    std::vector<Vector4d> normals = calcTriangleNormals(ps->vertices, ps->indices);
    std::vector<int> faceParents;
    std::vector<Vector4d> newnormals;
    ps_merged->indices = mergeTriangles(ps->indices, normals, newnormals, faceParents, ps->vertices);

  } else return std::unique_ptr<PolySet>();
  if (this->children.size() >= 2) {
    std::shared_ptr<const PolySet> sel = childToPolySet(this->children[1]);
    if (sel != nullptr) {
      auto sel_tess = PolySetUtils::tessellate_faces(*sel);
      for (size_t i = 0; i < ps_merged->vertices.size(); i++) {
        corner_selected.push_back(sel_tess->point_inside(ps_merged->vertices[i]));
      }
    }

  } else {
    for (size_t i = 0; i < ps_merged->vertices.size(); i++) corner_selected.push_back(true);
  }
  return createFilletInt(ps_merged, corner_selected, this->r, this->fn, this->minang, 0);
}
