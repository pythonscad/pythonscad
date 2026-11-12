#include "core/PreviewObjectNavigator.h"

#include <catch2/catch_all.hpp>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

#include "core/CSGNode.h"
#include "core/CSGTreeEvaluator.h"
#include "core/CsgOpNode.h"
#include "core/ModuleInstantiation.h"
#include "core/RenderNode.h"
#include "core/TransformNode.h"
#include "core/Tree.h"
#include "core/primitives.h"
#include "geometry/GeometryEvaluator.h"

namespace {

std::shared_ptr<ModuleInstantiation> source(const std::string& name, int line,
                                            const std::string& file = "/tmp/preview-navigation.py")
{
  return std::make_shared<ModuleInstantiation>(
    name, AssignmentList(), Location(line, 1, line, 1, std::make_shared<fs::path>(file)));
}

std::shared_ptr<ModuleInstantiation> noSource(const std::string& name)
{
  return std::make_shared<ModuleInstantiation>(name, AssignmentList(), Location::NONE);
}

// The representative (first) node of the group a step returned.
std::shared_ptr<const AbstractNode> front(const PreviewObjectNavigator::Group *group)
{
  return group && !group->nodes.empty() ? group->nodes.front() : nullptr;
}

}  // namespace

TEST_CASE("Preview objects are listed in pre-order without modifiers", "[preview-navigation]")
{
  // union (line 5)
  //   cube (line 2)
  //   translate (line 3)  <- modifier, skipped
  //     cube (line 3)
  //   group without location  <- skipped, but its child is listed
  //     cube (line 4)
  //   echo-only group (line 6)  <- no geometry, skipped
  auto unionNode = std::make_shared<CsgOpNode>(source("union", 5), OpenSCADOperator::UNION);
  auto first = std::make_shared<CubeNode>(source("cube", 2));
  auto transform = std::make_shared<TransformNode>(source("translate", 3), "translate");
  auto second = std::make_shared<CubeNode>(source("cube", 3));
  transform->children.push_back(second);
  auto synthetic = std::make_shared<GroupNode>(noSource("group"));
  auto third = std::make_shared<CubeNode>(source("cube", 4));
  synthetic->children.push_back(third);
  auto echoOnly = std::make_shared<GroupNode>(source("echo", 6));
  echoOnly->children.push_back(std::make_shared<GroupNode>(source("echo", 6)));
  unionNode->children = {first, transform, synthetic, echoOnly};

  PreviewObjectNavigator navigation;
  CHECK(navigation.step(1) == nullptr);  // nothing before the first preview
  navigation.reset(unionNode);
  REQUIRE(navigation.size() == 4);

  SECTION("+1 after a preview starts at the first object, then clamps at the end")
  {
    CHECK(front(navigation.step(1)) == unionNode);
    CHECK(navigation.currentPosition() == 0);
    CHECK(front(navigation.step(1)) == first);
    CHECK(front(navigation.step(1)) == second);
    CHECK(front(navigation.step(1)) == third);
    CHECK(front(navigation.step(1)) == third);
    CHECK(navigation.currentPosition() == 3);
    CHECK(front(navigation.step(-1)) == second);
    CHECK(front(navigation.step(0)) == second);
  }

  SECTION("-1 after a preview starts at the last object, then clamps at the start")
  {
    CHECK(front(navigation.step(-1)) == third);
    CHECK(front(navigation.step(-1)) == second);
    CHECK(front(navigation.step(-1)) == first);
    CHECK(front(navigation.step(-1)) == unionNode);
    CHECK(front(navigation.step(-1)) == unionNode);
    CHECK(front(navigation.step(1)) == first);
  }

  SECTION("rewind and reset restart from an end, clear empties the list")
  {
    CHECK(navigation.step(0) == nullptr);
    navigation.step(1);
    navigation.step(1);
    navigation.rewind();
    CHECK(navigation.current() == nullptr);
    CHECK(front(navigation.step(-1)) == third);
    navigation.reset(unionNode);
    CHECK(front(navigation.step(1)) == unionNode);
    navigation.clear();
    CHECK(navigation.size() == 0);
    CHECK(navigation.step(1) == nullptr);
    CHECK(navigation.step(-1) == nullptr);
  }
}

TEST_CASE("Objects created by one line form one step, also in a loop", "[preview-navigation]")
{
  // for i in range(100):        # line 1
  //     parts.append(cube(1))   # line 2
  //     parts.append(sphere(1)) # line 3   (modelled as cubes on line 3)
  // show(union(parts))          # line 4
  auto root = std::make_shared<CsgOpNode>(source("union", 4), OpenSCADOperator::UNION);
  std::vector<std::shared_ptr<CubeNode>> line2, line3;
  for (int i = 0; i < 100; ++i) {
    line2.push_back(std::make_shared<CubeNode>(source("cube", 2)));
    line3.push_back(std::make_shared<CubeNode>(source("sphere", 3)));
    root->children.push_back(line2.back());
    root->children.push_back(line3.back());
  }
  // Same line number in another file is a different line.
  auto other = std::make_shared<CubeNode>(source("cube", 2, "/tmp/module.py"));
  root->children.push_back(other);

  PreviewObjectNavigator navigation;
  navigation.reset(root);
  REQUIRE(navigation.size() == 4);  // lines 4, 2, 3 and module.py:2

  const auto *unionGroup = navigation.step(1);
  REQUIRE(unionGroup);
  CHECK(unionGroup->nodes.size() == 1);
  CHECK(unionGroup->highlightIndices == std::vector<int>{root->index()});

  const auto *cubes = navigation.step(1);
  REQUIRE(cubes);
  REQUIRE(cubes->nodes.size() == 100);
  CHECK(cubes->nodes.front() == line2.front());
  CHECK(cubes->nodes.back() == line2.back());
  CHECK(cubes->highlightIndices.size() == 100);

  const auto *spheres = navigation.step(1);
  REQUIRE(spheres);
  CHECK(spheres->nodes.size() == 100);
  CHECK(spheres->nodes.front() == line3.front());

  CHECK(front(navigation.step(1)) == other);
  CHECK(front(navigation.step(-1)) == line3.front());
  CHECK(front(navigation.step(-1)) == line2.front());
}

TEST_CASE("Nested objects of one line are highlighted through the outermost one", "[preview-navigation]")
{
  // show(a | b | c) on line 5: (a | b) | c gives two unions on the same line.
  auto a = std::make_shared<CubeNode>(source("cube", 2));
  auto b = std::make_shared<CubeNode>(source("sphere", 3));
  auto c = std::make_shared<CubeNode>(source("cylinder", 4));
  auto inner = std::make_shared<CsgOpNode>(source("union", 5), OpenSCADOperator::UNION);
  inner->children = {a, b};
  auto outer = std::make_shared<CsgOpNode>(source("union", 5), OpenSCADOperator::UNION);
  outer->children = {inner, c};

  PreviewObjectNavigator navigation;
  navigation.reset(outer);
  REQUIRE(navigation.size() == 4);  // line 5, then 2, 3, 4
  const auto *unions = navigation.step(1);
  REQUIRE(unions);
  CHECK(unions->nodes.size() == 2);
  CHECK(unions->highlightIndices == std::vector<int>{outer->index()});
  CHECK(front(navigation.step(1)) == a);
  CHECK(front(navigation.step(1)) == b);
  CHECK(front(navigation.step(1)) == c);
}

TEST_CASE("A deep chain of operations does not overflow the stack", "[preview-navigation]")
{
  // r = cube(); for ...: r = r | cube()   -> a left-deep chain of unions
  std::shared_ptr<AbstractNode> chain = std::make_shared<CubeNode>(source("cube", 1));
  for (int i = 0; i < 5000; ++i) {
    auto op = std::make_shared<CsgOpNode>(source("union", 2), OpenSCADOperator::UNION);
    op->children = {chain, std::make_shared<CubeNode>(source("cube", 2))};
    chain = op;
  }
  PreviewObjectNavigator navigation;
  navigation.reset(chain);
  REQUIRE(navigation.size() == 2);  // line 2 (unions + their cubes), line 1
  const auto *line2 = navigation.step(1);
  REQUIRE(line2);
  CHECK(line2->nodes.size() == 10000);
  CHECK(line2->highlightIndices == std::vector<int>{chain->index()});
}

TEST_CASE("Preview navigation skips nodes without a source location", "[preview-navigation]")
{
  auto root = std::make_shared<GroupNode>(noSource("group"));
  auto located = std::make_shared<CubeNode>(source("cube", 1));
  auto unlocated = std::make_shared<CubeNode>(noSource("cube"));
  root->children = {unlocated, located};
  PreviewObjectNavigator navigation;
  navigation.reset(root);
  REQUIRE(navigation.size() == 1);
  CHECK(front(navigation.step(1)) == located);
  CHECK(front(navigation.step(-1)) == located);
}

TEST_CASE("Preview navigation treats render() as a single object", "[preview-navigation]")
{
  auto render = std::make_shared<RenderNode>(source("render", 1));
  render->children = {std::make_shared<CubeNode>(source("cube", 2))};
  PreviewObjectNavigator navigation;
  navigation.reset(render);
  REQUIRE(navigation.size() == 1);
  CHECK(front(navigation.step(1)) == render);
}

TEST_CASE("Preview highlight marks nodes without changing their source tags", "[preview-navigation]")
{
  auto transformSource = source("translate", 1);
  auto cubeSource = source("cube", 2);
  auto root = std::make_shared<CsgOpNode>(source("union", 3), OpenSCADOperator::UNION);
  auto transform = std::make_shared<TransformNode>(transformSource, "translate");
  transform->matrix.translate(Vector3d(12, -3, 4));
  auto first = std::make_shared<CubeNode>(cubeSource);
  auto second = std::make_shared<CubeNode>(source("cube", 4));
  transform->children = {first};
  root->children = {transform, second};
  Tree tree(root);
  GeometryEvaluator geometry(tree);
  const auto original = tree.getString(*root, "");

  CSGTreeEvaluator normal(tree, &geometry);
  normal.buildCSGTree(*root);
  CHECK(normal.getHighlightNodes().empty());

  CSGTreeEvaluator selected(tree, &geometry, {first->index()});
  selected.buildCSGTree(*root);
  REQUIRE(selected.getHighlightNodes().size() == 1);
  const auto leaf = std::dynamic_pointer_cast<CSGLeaf>(selected.getHighlightNodes().front());
  REQUIRE(leaf);
  CHECK(leaf->index == first->index());
  CHECK(leaf->matrix.translation().isApprox(Vector3d(12, -3, 4)));
  CHECK_FALSE(cubeSource->isHighlight());
  CHECK_FALSE(transformSource->isHighlight());
  CHECK(tree.getString(*root, "") == original);

  // A whole group: both cubes light up at once. As with `#` on both operands,
  // the evaluator merges them into one highlighted union term.
  CSGTreeEvaluator both(tree, &geometry, {first->index(), second->index()});
  both.buildCSGTree(*root);
  REQUIRE(both.getHighlightNodes().size() == 1);
  const auto merged = std::dynamic_pointer_cast<CSGOperation>(both.getHighlightNodes().front());
  REQUIRE(merged);
  CHECK(merged->getType() == OpenSCADOperator::UNION);
  const auto left = std::dynamic_pointer_cast<CSGLeaf>(merged->left());
  const auto right = std::dynamic_pointer_cast<CSGLeaf>(merged->right());
  REQUIRE(left);
  REQUIRE(right);
  CHECK(std::unordered_set<int>{left->index, right->index} ==
        std::unordered_set<int>{first->index(), second->index()});

  CSGTreeEvaluator wholeRoot(tree, &geometry, {root->index()});
  wholeRoot.buildCSGTree(*root);
  CHECK(wholeRoot.getHighlightNodes().size() == 1);
}

TEST_CASE("Preview highlight does not merge with the design's own # objects", "[preview-navigation]")
{
  // body = union([sphere(), sphere()], +cube())   # `+` is `#`
  // show(difference(body, cutter))
  // Selecting the spheres must highlight only them. With the `#` cube taking
  // part, the union rule would merge spheres and cube into one highlighted
  // operand and the difference rule would then highlight the whole design.
  auto first = std::make_shared<CubeNode>(source("sphere", 9));
  auto second = std::make_shared<CubeNode>(source("sphere", 9));
  auto spheres = std::make_shared<CsgOpNode>(source("union", 8), OpenSCADOperator::UNION);
  spheres->children = {first, second};
  auto tagged = source("union", 10);
  tagged->tag_highlight = true;
  auto marked = std::make_shared<CsgOpNode>(tagged, OpenSCADOperator::UNION);
  marked->children = {std::make_shared<CubeNode>(source("cube", 10))};
  auto body = std::make_shared<CsgOpNode>(source("union", 8), OpenSCADOperator::UNION);
  body->children = {spheres, marked};
  auto cutter = std::make_shared<CubeNode>(source("cube", 13));
  auto root = std::make_shared<CsgOpNode>(source("difference", 21), OpenSCADOperator::DIFFERENCE);
  root->children = {body, cutter};
  Tree tree(root);
  GeometryEvaluator geometry(tree);

  CSGTreeEvaluator selected(tree, &geometry, {first->index(), second->index()});
  selected.buildCSGTree(*root);
  REQUIRE(selected.getHighlightNodes().size() == 1);
  const auto term = std::dynamic_pointer_cast<CSGOperation>(selected.getHighlightNodes().front());
  REQUIRE(term);
  CHECK(term->getType() == OpenSCADOperator::UNION);
  const auto left = std::dynamic_pointer_cast<CSGLeaf>(term->left());
  const auto right = std::dynamic_pointer_cast<CSGLeaf>(term->right());
  REQUIRE(left);
  REQUIRE(right);
  CHECK(std::unordered_set<int>{left->index, right->index} ==
        std::unordered_set<int>{first->index(), second->index()});

  // Without a selection the design's `#` works as before.
  CSGTreeEvaluator normal(tree, &geometry);
  normal.buildCSGTree(*root);
  REQUIRE(normal.getHighlightNodes().size() == 1);
  CHECK(std::dynamic_pointer_cast<CSGLeaf>(normal.getHighlightNodes().front()));
}
