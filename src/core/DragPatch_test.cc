#include "core/DragPatch.h"

#include <catch2/catch_all.hpp>
#include <string>
#include <vector>

namespace {

DragMod mod(const std::string& name, int index, std::vector<int> arrinfo, double value)
{
  DragMod m;
  m.name = name;
  m.index = index;
  m.arrinfo = std::move(arrinfo);
  m.value = value;
  return m;
}

DragResult cubeResize(double x, double y, double z)
{
  DragResult r;
  r.modname = "cube";
  r.mods = {mod("size", 0, {0}, x), mod("size", 0, {1}, y), mod("size", 0, {2}, z)};
  return r;
}

std::string patched(const std::string& src, const DragResult& r)
{
  std::string out;
  REQUIRE(applyDragResult(src, r, out));
  return out;
}

}  // namespace

TEST_CASE("formatDragValue removes floating point noise", "[core][drag]")
{
  CHECK(formatDragValue(12.300000000000001) == "12.3");
  CHECK(formatDragValue(7.0) == "7");
  CHECK(formatDragValue(5.55555) == "5.5556");
  CHECK(formatDragValue(-0.00001) == "0");
  CHECK(formatDragValue(1234567.5) == "1234567.5");
}

TEST_CASE("applyDragResult writes all values of a call", "[core][drag]")
{
  const auto r = cubeResize(12.300000000000001, 7, 5.55555);

  SECTION("double and single quoted values, positional")
  {
    CHECK(patched("c = cube([\"10\", '20', \"30\"])\n", r) == "c = cube([\"12.3\", '7', \"5.5556\"])\n");
  }

  SECTION("named parameter with whitespace")
  {
    CHECK(patched("x = cube( size = [ \"1\" , \"2\" ,'3' ] , center=True)\n", r) ==
          "x = cube( size = [ \"12.3\" , \"7\" ,'5.5556' ] , center=True)\n");
  }

  SECTION("marked() values")
  {
    CHECK(patched("cube([marked(10), marked( 20 ), marked(value=30)])", r) ==
          "cube([marked(12.3), marked( 7 ), marked(value=5.5556)])");
  }

  SECTION("module qualified marked() and mixed markers")
  {
    CHECK(patched("cube([pythonscad.marked(1), \"2\", marked(3)])", r) ==
          "cube([pythonscad.marked(12.3), \"7\", marked(5.5556)])");
  }
}

TEST_CASE("applyDragResult finds the right call", "[core][drag]")
{
  const auto r = cubeResize(4, 5, 6);

  SECTION("longer identifiers and unmarked calls are skipped")
  {
    CHECK(patched("a = cuboid(3)\nb = mycube([\"1\",\"1\",\"1\"])\nc = cube([1,2,3])\n"
                  "d = cube([\"1\",\"2\",\"3\"])\n",
                  r) ==
          "a = cuboid(3)\nb = mycube([\"1\",\"1\",\"1\"])\nc = cube([1,2,3])\n"
          "d = cube([\"4\",\"5\",\"6\"])\n");
  }

  SECTION("python comments inside the call are ignored")
  {
    CHECK(patched("d = cube([\"1\",\"2\",\"3\"]) # cube(\n", r) ==
          "d = cube([\"4\",\"5\",\"6\"]) # cube(\n");
  }

  SECTION("method call syntax")
  {
    DragResult t;
    t.modname = "rotate";
    t.mods = {mod("angle", 0, {2}, 45)};
    CHECK(patched("obj = cube(5).rotate([0, 0, marked(30)])\n", t) ==
          "obj = cube(5).rotate([0, 0, marked(45)])\n");
  }
}

TEST_CASE("applyDragResult refuses unmarked or unknown values", "[core][drag]")
{
  const auto r = cubeResize(4, 5, 6);
  std::string out;
  CHECK_FALSE(applyDragResult("cube([1,2,3])", r, out));
  CHECK_FALSE(applyDragResult("cube([\"1\",\"2\"])", r, out));    // index out of range
  CHECK_FALSE(applyDragResult("cube([\"1\",\"2\",3])", r, out));  // z not marked
  CHECK_FALSE(applyDragResult("cube([marker(1),\"2\",\"3\"])", r, out));
  CHECK_FALSE(applyDragResult("sphere(\"1\")", r, out));
  CHECK_FALSE(applyDragResult("cube([\"1\",\"2\",\"3\"]", r, out));  // unterminated

  DragResult empty;
  CHECK_FALSE(applyDragResult("cube([\"1\",\"2\",\"3\"])", empty, out));
}
