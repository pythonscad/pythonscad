#include <catch2/catch_all.hpp>

#include "geometry/Polygon2d.h"
#include "geometry/linalg.h"

namespace {

Polygon2d unitSquare()
{
  Outline2d outline;
  outline.vertices = {
    Vector2d(0, 0),
    Vector2d(1, 0),
    Vector2d(1, 1),
    Vector2d(0, 1),
  };
  outline.positive = true;
  return Polygon2d(outline);
}

}  // namespace

TEST_CASE("Polygon2d::point_location classifies interior, exterior, edge, and vertex", "[polygon2d]")
{
  const Polygon2d square = unitSquare();

  CHECK(square.point_location(Vector2d(0.5, 0.5)) == PointLocation2d::Inside);
  CHECK(square.point_location(Vector2d(2.0, 2.0)) == PointLocation2d::Outside);
  CHECK(square.point_location(Vector2d(0.5, 0.0)) == PointLocation2d::OnEdge);
  CHECK(square.point_location(Vector2d(0.0, 0.0)) == PointLocation2d::OnVertex);
}

TEST_CASE("PointLocation2d ordering supports threshold comparisons", "[polygon2d]")
{
  const Polygon2d square = unitSquare();

  // >= OnVertex: inside or on boundary (replaces the old point_inside).
  CHECK(square.point_location(Vector2d(0.5, 0.5)) >= PointLocation2d::OnVertex);
  CHECK(square.point_location(Vector2d(0.5, 0.0)) >= PointLocation2d::OnVertex);
  CHECK(square.point_location(Vector2d(0.0, 0.0)) >= PointLocation2d::OnVertex);
  CHECK_FALSE(square.point_location(Vector2d(2.0, 2.0)) >= PointLocation2d::OnVertex);

  // >= OnEdge: on edge or strictly inside (excludes vertices).
  CHECK(square.point_location(Vector2d(0.5, 0.5)) >= PointLocation2d::OnEdge);
  CHECK(square.point_location(Vector2d(0.5, 0.0)) >= PointLocation2d::OnEdge);
  CHECK_FALSE(square.point_location(Vector2d(0.0, 0.0)) >= PointLocation2d::OnEdge);
}
