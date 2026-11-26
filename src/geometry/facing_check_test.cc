#include "geometry/facing_check.h"

#include <catch2/catch_all.hpp>
#include <cmath>
#include <memory>
#include <utility>
#include <vector>

#include "geometry/PolySet.h"
#include "geometry/linalg.h"

using FacingCheck::Mode;
using FacingCheck::Options;

namespace {

// Extrude a CCW polygon (XY) from z0 to z1 into ps, outward oriented.
void extrude(PolySet& ps, const std::vector<std::pair<double, double>>& poly, double z0, double z1)
{
  const int n = (int)poly.size(), b = (int)ps.vertices.size();
  for (const auto& p : poly) ps.vertices.emplace_back(p.first, p.second, z0);
  for (const auto& p : poly) ps.vertices.emplace_back(p.first, p.second, z1);
  for (int i = 1; i + 1 < n; ++i) {
    ps.indices.push_back({b, b + i + 1, b + i});
    ps.indices.push_back({b + n, b + n + i, b + n + i + 1});
  }
  for (int i = 0; i < n; ++i) {
    const int j = (i + 1) % n;
    ps.indices.push_back({b + i, b + j, b + n + j});
    ps.indices.push_back({b + i, b + n + j, b + n + i});
  }
  ps.setTriangular(true);
}

void box(PolySet& ps, double x0, double y0, double z0, double x1, double y1, double z1)
{
  extrude(ps, {{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}}, z0, z1);
}

// Hollow cylinder along z, n segments.
PolySet tube(double R, double r, double h, int n)
{
  PolySet ps(3);
  for (int k = 0; k < 2; ++k) {
    for (int i = 0; i < n; ++i) {
      const double a = 2 * M_PI * i / n, rad = k == 0 ? R : r;
      ps.vertices.emplace_back(rad * std::cos(a), rad * std::sin(a), 0);
      ps.vertices.emplace_back(rad * std::cos(a), rad * std::sin(a), h);
    }
  }
  auto O0 = [&](int i) { return 2 * (i % n); };
  auto O1 = [&](int i) { return 2 * (i % n) + 1; };
  auto I0 = [&](int i) { return 2 * n + 2 * (i % n); };
  auto I1 = [&](int i) { return 2 * n + 2 * (i % n) + 1; };
  for (int i = 0; i < n; ++i) {
    const int j = i + 1;
    ps.indices.push_back({O0(i), O0(j), O1(j)});
    ps.indices.push_back({O0(i), O1(j), O1(i)});
    ps.indices.push_back({I0(i), I1(j), I0(j)});
    ps.indices.push_back({I0(i), I1(i), I1(j)});
    ps.indices.push_back({O1(i), O1(j), I1(j)});
    ps.indices.push_back({O1(i), I1(j), I1(i)});
    ps.indices.push_back({O0(i), I0(j), O0(j)});
    ps.indices.push_back({O0(i), I0(i), I0(j)});
  }
  ps.setTriangular(true);
  return ps;
}

Options with(double alpha, double angle = 120, bool occlusion = true)
{
  Options o;
  o.alpha_deg = alpha;
  o.min_angle_deg = angle;
  o.occlusion = occlusion;
  return o;
}

}  // namespace

TEST_CASE("FacingCheck plate thickness")
{
  PolySet ps(3);
  box(ps, 0, 0, 0, 10, 10, 0.5);
  const auto r = FacingCheck::check(ps, Mode::Internal, 0.8);
  REQUIRE_FALSE(r.clean());
  CHECK(r.min_distance == Catch::Approx(0.5).epsilon(1e-12));
  CHECK(FacingCheck::check(ps, Mode::Internal, 0.3).clean());
  CHECK(FacingCheck::check(ps, Mode::External, 0.8).clean());
}

TEST_CASE("FacingCheck gaps within one solid and between two")
{
  PolySet a(3), b(3);
  box(a, 0, 0, 0, 1, 1, 1);
  box(a, 1.3, 0, 0, 2.3, 1, 1);
  box(b, 2.7, 0, 0, 3.7, 1, 1);
  CHECK(FacingCheck::check(a, Mode::External, 0.5).min_distance == Catch::Approx(0.3));
  const auto r = FacingCheck::checkBetween(a, b, 0.5);
  REQUIRE_FALSE(r.clean());
  CHECK(r.min_distance == Catch::Approx(0.4));
  for (const auto& v : r.violations) {
    CHECK(v.tri_a < (int)r.split);
    CHECK(v.tri_b >= (int)r.split);
  }
  CHECK(FacingCheck::checkBetween(a, b, 0.35).clean());
  CHECK(FacingCheck::check(a, Mode::Internal, 0.8).clean());
}

TEST_CASE("FacingCheck tube: occlusion and oblique measurement")
{
  const PolySet t = tube(5, 4.5, 10, 64);
  const double apothem = 4.5 * std::cos(M_PI / 64);

  const auto wall = FacingCheck::check(t, Mode::Internal, 0.8);
  CHECK(wall.min_distance > 0.49);
  CHECK(wall.min_distance <= 0.5);

  // no connecting segment may run through the hole
  const auto big = FacingCheck::check(t, Mode::Internal, 20);
  for (const auto& v : big.violations) {
    const Vector3d m = (v.p + v.q) * 0.5;
    CHECK(std::hypot(m[0], m[1]) > apothem - 1e-6);
  }

  CHECK(FacingCheck::check(t, Mode::External, 20, with(0)).min_distance == Catch::Approx(2 * apothem));
  CHECK(FacingCheck::check(t, Mode::External, 20, with(90, 175)).min_distance ==
        Catch::Approx(2 * apothem));
  CHECK(FacingCheck::check(t, Mode::External, 20, with(90, 120)).min_distance < 7.8);
}

TEST_CASE("FacingCheck wedge error region")
{
  PolySet w(3);
  extrude(w, {{-1, 0}, {1, 0}, {0, 10}}, 0, 5);
  const auto r = FacingCheck::check(w, Mode::Internal, 1.0);
  REQUIRE_FALSE(r.clean());
  double ymin = 1e9;
  for (const auto& p : r.hull_points) ymin = std::min(ymin, p[1]);
  CHECK(ymin == Catch::Approx(5.0).margin(0.1));  // width 1.0 at y = 5
  CHECK(FacingCheck::check(w, Mode::Internal, 1.0, with(90, 170)).clean());
}

TEST_CASE("FacingCheck projection cone")
{
  PolySet ps(3);
  box(ps, 0, 0, 0, 10, 10, 1);
  box(ps, 10.2, 0, 1.3, 20, 10, 2.3);
  CHECK(FacingCheck::check(ps, Mode::External, 0.5, with(90)).min_distance ==
        Catch::Approx(std::sqrt(0.13)));
  CHECK_FALSE(FacingCheck::check(ps, Mode::External, 0.5, with(45)).clean());
  CHECK(FacingCheck::check(ps, Mode::External, 0.5, with(30)).clean());
  CHECK(FacingCheck::check(ps, Mode::External, 0.5, with(0)).clean());
}

TEST_CASE("FacingCheck projection cone at acute triangle corners")
{
  // Same offset plates, but triangulated so that the corner facing the gap
  // is an acute 45 degree triangle corner (as Manifold does for cubes).
  PolySet ps(3);
  extrude(ps, {{10, 0}, {10, 10}, {0, 10}, {0, 0}}, 0, 1);
  extrude(ps, {{10.2, 0}, {20, 0}, {20, 10}, {10.2, 10}}, 1.3, 2.3);
  CHECK(FacingCheck::check(ps, Mode::External, 0.5, with(90)).min_distance ==
        Catch::Approx(std::sqrt(0.13)));
  CHECK_FALSE(FacingCheck::check(ps, Mode::External, 0.5, with(45)).clean());  // needs 33.7 deg
  CHECK(FacingCheck::check(ps, Mode::External, 0.5, with(30)).clean());
  CHECK(FacingCheck::check(ps, Mode::External, 0.5, with(0)).clean());
}

TEST_CASE("FacingCheck step needs occlusion in oblique mode")
{
  PolySet ps(3);
  box(ps, 0, 0, 0, 5, 5, 10);
  box(ps, 6, 0, 8, 11, 5, 12.5);
  CHECK(FacingCheck::check(ps, Mode::Internal, 4.4, with(90)).clean());
  CHECK(FacingCheck::check(ps, Mode::Internal, 4.4, with(90, 120, false)).min_distance ==
        Catch::Approx(std::sqrt(5.0)));
  CHECK(FacingCheck::check(ps, Mode::Internal, 4.4, with(0, 120, false)).clean());
}

TEST_CASE("FacingCheck result is independent of thread count")
{
  const PolySet t = tube(5, 4.5, 10, 96);
  Options o1, o4;
  o1.threads = 1;
  o4.threads = 4;
  const auto a = FacingCheck::check(t, Mode::Internal, 0.8, o1);
  const auto b = FacingCheck::check(t, Mode::Internal, 0.8, o4);
  REQUIRE(a.violations.size() == b.violations.size());
  for (size_t i = 0; i < a.violations.size(); ++i) {
    CHECK(a.violations[i].tri_a == b.violations[i].tri_a);
    CHECK(a.violations[i].tri_b == b.violations[i].tri_b);
    CHECK(a.violations[i].distance == b.violations[i].distance);
  }
}

TEST_CASE("FacingCheck summary mode")
{
  const PolySet t = tube(5, 4.5, 10, 64);
  Options o;
  o.store_pairs = false;
  const auto s = FacingCheck::check(t, Mode::Internal, 0.8, o);
  const auto full = FacingCheck::check(t, Mode::Internal, 0.8);
  CHECK(s.violations.empty());
  CHECK(s.violation_count == full.violation_count);
  CHECK(s.min_distance == full.min_distance);
  CHECK(s.tri_min_distance == full.tri_min_distance);
}

TEST_CASE("FacingCheck error solid grow lifts every face")
{
  auto ps = std::make_shared<PolySet>(3);
  box(*ps, 0, 0, 0, 10, 10, 0.5);
  const auto r = FacingCheck::check(*ps, Mode::Internal, 0.8);
  const std::vector<std::shared_ptr<const Geometry>> bodies{ps};
  const auto exact = FacingCheck::errorGeometry(r, Mode::Internal, bodies, 20000, 0.0);
  const auto grown = FacingCheck::errorGeometry(r, Mode::Internal, bodies, 20000, 0.05);
  REQUIRE(exact);
  REQUIRE(grown);
  const BoundingBox be = exact->getBoundingBox(), bg = grown->getBoundingBox();
  for (int k = 0; k < 3; ++k) {
    CHECK(be.min()[k] == Catch::Approx(0).margin(1e-9));
    CHECK(bg.min()[k] <= -0.05 + 1e-9);  // every face moved out by at least grow
    CHECK(bg.max()[k] >= be.max()[k] + 0.05 - 1e-9);
    CHECK(bg.min()[k] >= -4 * 0.05 - 1e-9);  // capped
  }
}
