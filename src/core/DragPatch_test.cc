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

TEST_CASE("applyDragResult handles scalar parameters", "[core][drag]")
{
  SECTION("cube with scalar size")
  {
    DragResult r;
    r.modname = "cube";
    r.mods = {mod("size", 0, {}, 12)};
    CHECK(patched("cube(marked(10), center=True)", r) == "cube(marked(12), center=True)");
    CHECK(patched("cube(\"10\")", r) == "cube(\"12\")");
  }

  SECTION("sphere r and d")
  {
    DragResult r;
    r.modname = "sphere";
    r.mods = {mod("r", 0, {}, 7.5)};
    CHECK(patched("sphere(marked(5))", r) == "sphere(marked(7.5))");
    CHECK(patched("sphere(r=\"5\", fn=40)", r) == "sphere(r=\"7.5\", fn=40)");
    r.mods = {mod("d", 1, {}, 15)};
    CHECK(patched("sphere(d=marked(10))", r) == "sphere(d=marked(15))");
  }

  SECTION("cylinder: named r1/r2 are not confused with r")
  {
    DragResult r;
    r.modname = "cylinder";
    r.mods = {mod("h", 0, {}, 20), mod("r2", 2, {}, 3)};
    CHECK(patched("cylinder(marked(10), r1=5, r2=marked(5), r=1)", r) ==
          "cylinder(marked(20), r1=5, r2=marked(3), r=1)");
    r.mods = {mod("angle", 8, {}, 270)};
    CHECK(patched("cylinder(h=10, r=5, angle=\"90\")", r) == "cylinder(h=10, r=5, angle=\"270\")");
  }
}

TEST_CASE("applyDragResult patches polyhedron points", "[core][drag]")
{
  DragResult r;
  r.modname = "polyhedron";
  r.mods = {mod("points", 0, {2, 0}, 1.5), mod("points", 0, {2, 2}, -2)};

  SECTION("literal point list")
  {
    CHECK(patched("polyhedron([[0,0,0],[1,0,0],[marked(0),1,\"0\"],[0,0,1]], faces)", r) ==
          "polyhedron([[0,0,0],[1,0,0],[marked(1.5),1,\"-2\"],[0,0,1]], faces)");
  }

  SECTION("point list in a variable, multi line")
  {
    const std::string src =
      "pts = [\n"
      "  [0, 0, 0],\n"
      "  [1, 0, 0],\n"
      "  [marked(0), 1, marked(0)],  # apex\n"
      "]\n"
      "show(polyhedron(pts, [[0,1,2]]))\n";
    const std::string expected =
      "pts = [\n"
      "  [0, 0, 0],\n"
      "  [1, 0, 0],\n"
      "  [marked(1.5), 1, marked(-2)],  # apex\n"
      "]\n"
      "show(polyhedron(pts, [[0,1,2]]))\n";
    CHECK(patched(src, r) == expected);
  }

  SECTION("single point in a variable, latest assignment before the call wins")
  {
    const std::string src =
      "p = [0, 0, 0]\n"
      "p = [marked(0), 1, \"0\"]\n"
      "polyhedron(points=[[0,0,0],[1,0,0],p], faces=f)\n"
      "p = [9, 9, 9]\n";
    const std::string expected =
      "p = [0, 0, 0]\n"
      "p = [marked(1.5), 1, \"-2\"]\n"
      "polyhedron(points=[[0,0,0],[1,0,0],p], faces=f)\n"
      "p = [9, 9, 9]\n";
    CHECK(patched(src, r) == expected);
  }

  SECTION("unmarked point is refused")
  {
    std::string out;
    CHECK_FALSE(applyDragResult("polyhedron([[0,0,0],[1,0,0],[0,1,0]], f)", r, out));
  }
}

TEST_CASE("applyDragResult patches transforms in function and method form", "[core][drag]")
{
  DragResult r;
  r.objFirst = true;

  SECTION("translate")
  {
    r.modname = "translate";
    r.mods = {mod("v", 1, {0}, 12), mod("v", 1, {2}, -3)};
    CHECK(patched("c = cube(5).translate([marked(10), 0, \"0\"])", r) ==
          "c = cube(5).translate([marked(12), 0, \"-3\"])");
    CHECK(patched("c = translate(cube(5), [marked(10), 0, \"0\"])", r) ==
          "c = translate(cube(5), [marked(12), 0, \"-3\"])");
    CHECK(patched("c = pythonscad.translate(cube(5), [marked(10), 0, \"0\"])", r) ==
          "c = pythonscad.translate(cube(5), [marked(12), 0, \"-3\"])");
    CHECK(patched("c = cube(5).translate(v=[marked(10), 0, \"0\"])", r) ==
          "c = cube(5).translate(v=[marked(12), 0, \"-3\"])");
  }

  SECTION("rotate euler angle and axis-angle form")
  {
    r.modname = "rotate";
    r.mods = {mod("a", 1, {2}, 135)};
    CHECK(patched("cube(5).rotate([0, 0, marked(90)])", r) == "cube(5).rotate([0, 0, marked(135)])");
    CHECK(patched("rotate(cube(5), [0, 0, \"90\"])", r) == "rotate(cube(5), [0, 0, \"135\"])");
    r.mods = {mod("a", 1, {}, -45)};
    CHECK(patched("cube(5).rotate(marked(30), [1, 1, 0])", r) ==
          "cube(5).rotate(marked(-45), [1, 1, 0])");
    CHECK(patched("cube(5).rotate(a=\"30\", v=[1, 1, 0])", r) ==
          "cube(5).rotate(a=\"-45\", v=[1, 1, 0])");
  }

  SECTION("scale vector and scalar")
  {
    r.modname = "scale";
    r.mods = {mod("v", 1, {1}, 2.5)};
    CHECK(patched("scale(cube(5), [1, marked(2), 1])", r) == "scale(cube(5), [1, marked(2.5), 1])");
    r.mods = {mod("v", 1, {}, 0.75)};
    CHECK(patched("cube(5).scale(marked(1))", r) == "cube(5).scale(marked(0.75))");
  }
}
