#pragma once

#include <cstddef>
#include <list>
#include <map>
#include <memory>
#include <unordered_set>
#include <utility>
#include <vector>

#include "core/BaseVisitable.h"
#include "core/CSGNode.h"
#include "core/ModuleInstantiation.h"
#include "core/NodeVisitor.h"
#include "core/enums.h"
#include "core/node.h"
#include "geometry/Geometry.h"

class CSGNode;
class GeometryEvaluator;
class Tree;

class CSGTreeEvaluator : public NodeVisitor
{
public:
  /*!
    previewHighlightIndices: indices of nodes to render as if they carried the
    `#` (highlight) modifier, without touching their ModuleInstantiation. Used
    for the temporary selection of preview object navigation, which highlights
    all objects created by one source line; empty disables it.

    When the set is not empty it is the only source of highlighting: the
    design's own `#` modifiers are ignored. The `#` rules merge highlighted
    operands (two highlighted union operands become one highlighted union, a
    highlighted first operand of a difference highlights the whole difference),
    so mixing both would let a selected object merge with an unrelated `#`
    object and spread to its parents, e.g. highlight a whole body instead of
    the selected spheres. The caller draws the design's own highlights
    separately.
   */
  CSGTreeEvaluator(const Tree& tree, GeometryEvaluator *geomevaluator = nullptr,
                   std::unordered_set<int> previewHighlightIndices = {})
    : tree(tree),
      geomevaluator(geomevaluator),
      previewHighlightIndices(std::move(previewHighlightIndices))
  {
  }

  Response visit(State& state, const AbstractNode& node) override;
  Response visit(State& state, const AbstractIntersectionNode& node) override;
  Response visit(State& state, const AbstractPolyNode& node) override;
  Response visit(State& state, const ListNode& node) override;
  Response visit(State& state, const CsgOpNode& node) override;
  Response visit(State& state, const TransformNode& node) override;
  Response visit(State& state, const ColorNode& node) override;
  Response visit(State& state, const RenderNode& node) override;
  Response visit(State& state, const CgalAdvNode& node) override;

  std::shared_ptr<CSGNode> buildCSGTree(const AbstractNode& node);

  [[nodiscard]] const std::shared_ptr<CSGNode>& getRootNode() const { return this->rootNode; }
  [[nodiscard]] const std::vector<std::shared_ptr<CSGNode>>& getHighlightNodes() const
  {
    return this->highlightNodes;
  }
  [[nodiscard]] const std::vector<std::shared_ptr<CSGNode>>& getBackgroundNodes() const
  {
    return this->backgroundNodes;
  }

private:
  [[nodiscard]] bool isHighlighted(const AbstractNode& node) const
  {
    if (!previewHighlightIndices.empty()) return previewHighlightIndices.count(node.index()) > 0;
    return node.modinst && node.modinst->isHighlight();
  }
  void addToParent(const State& state, const AbstractNode& node);
  void applyToChildren(State& state, const AbstractNode& node, OpenSCADOperator op);
  std::shared_ptr<CSGNode> evaluateCSGNodeFromGeometry(
    State& state, const std::shared_ptr<const Geometry>& geom,
    const std::shared_ptr<const ModuleInstantiation>& modinst, const AbstractNode& node);
  void applyBackgroundAndHighlight(State& state, const AbstractNode& node);

  using ChildList = std::list<std::shared_ptr<const AbstractNode>>;
  std::map<int, ChildList> visitedchildren;

protected:
  const Tree& tree;
  GeometryEvaluator *geomevaluator;
  const std::unordered_set<int> previewHighlightIndices;
  std::shared_ptr<CSGNode> rootNode;
  std::vector<std::shared_ptr<CSGNode>> highlightNodes;
  std::vector<std::shared_ptr<CSGNode>> backgroundNodes;
  std::map<int, std::shared_ptr<CSGNode>> stored_term;  // The term evaluated from each node index
};
