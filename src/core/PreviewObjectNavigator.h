#pragma once

#include <algorithm>
#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "core/CgalAdvNode.h"
#include "core/ColorNode.h"
#include "core/ModuleInstantiation.h"
#include "core/RenderNode.h"
#include "core/TransformNode.h"
#include "core/node.h"

/*!
  The selectable objects of one evaluated preview tree, grouped by the source
  line that created them, used to step through the design (Alt+Down / Alt+Up in
  the GUI).

  Why groups: one source line can create many objects, e.g. a cube() call in a
  loop that runs hundreds of times, or `a | b | c` creating two unions. Stepping
  through every single instance would make navigation unusable, and the
  instances cannot be told apart in the editor anyway, because the editor can
  only select the line. So one step selects one line and highlights everything
  that line created.

  Objects are visited in pre-order (a parent before its children, children in
  order). A group is placed at the position of its first object; later objects
  of the same line join that group even when other lines come in between (two
  lines in a loop body give two groups, not 2 * iterations). Transform and
  color nodes, empty/echo-only groups and nodes without a source location are
  not objects. Nodes that turn into a single CSG leaf (AbstractPolyNode,
  RenderNode, CgalAdvNode) are not descended into.

  Navigation never changes the tree, the AST or the source text.
 */
class PreviewObjectNavigator
{
public:
  struct Group {
    /// All objects created by the line, in tree order. front() represents the
    /// group (source selection, status bar).
    std::vector<std::shared_ptr<const AbstractNode>> nodes;
    /// Node indices to highlight: the members that are not inside another
    /// member of the same group. Highlighting a member's descendants as well
    /// would only draw the same geometry twice (e.g. the inner union of
    /// `a | b | c`). Free of duplicates (a shared node has a single index).
    std::vector<int> highlightIndices;
  };

  void reset(const std::shared_ptr<const AbstractNode>& root)
  {
    clear();
    if (!root) return;

    // Pass 1: pre-order list, used to find the nodes that produce geometry.
    std::vector<std::shared_ptr<const AbstractNode>> ordered;
    {
      std::vector<std::shared_ptr<const AbstractNode>> pending{root};
      while (!pending.empty()) {
        auto node = pending.back();
        pending.pop_back();
        ordered.push_back(node);
        // These operations produce one CSG leaf; their input geometry is not an
        // independently visible object in the compiled preview.
        if (isGeometryLeaf(*node)) continue;
        for (auto child = node->children.rbegin(); child != node->children.rend(); ++child) {
          if (*child) pending.push_back(*child);
        }
      }
    }
    // Skip echo/assert-only branches and empty structural groups.
    std::unordered_set<int> geometric;
    for (auto node = ordered.rbegin(); node != ordered.rend(); ++node) {
      if (isGeometryLeaf(**node) || std::any_of((*node)->children.begin(), (*node)->children.end(),
                                                [&geometric](const auto& child) {
                                                  return child && geometric.count(child->index()) > 0;
                                                })) {
        geometric.insert((*node)->index());
      }
    }

    // Pass 2: depth-first again, with enter/exit events, to know for each object
    // whether an ancestor on the current path belongs to the same line. The
    // traversal is iterative: Python designs built in a loop (`r = r | x`) can be
    // thousands of levels deep.
    std::map<LineKey, size_t> groupOfLine;
    std::map<LineKey, int> openMembers;  // members of a line on the current path
    std::vector<std::unordered_set<int>> highlighted;
    struct Event {
      std::shared_ptr<const AbstractNode> node;
      bool exit;
      const LineKey *key;  // set for objects on exit, to close them
    };
    std::vector<Event> pending{{root, false, nullptr}};
    while (!pending.empty()) {
      Event event = pending.back();
      pending.pop_back();
      if (event.exit) {
        if (event.key) --openMembers[*event.key];
        continue;
      }
      const auto& node = event.node;
      if (!geometric.count(node->index())) continue;
      const LineKey *openedKey = nullptr;
      if (isObject(*node)) {
        const auto& location = node->modinst->location();
        auto [it, inserted] =
          groupOfLine.emplace(LineKey{location.fileName(), location.firstLine()}, groups.size());
        if (inserted) {
          groups.emplace_back();
          highlighted.emplace_back();
        }
        auto& group = groups[it->second];
        group.nodes.push_back(node);
        auto& open = openMembers[it->first];
        if (open == 0 && highlighted[it->second].insert(node->index()).second) {
          group.highlightIndices.push_back(node->index());
        }
        ++open;
        openedKey = &it->first;  // std::map keys are stable
      }
      pending.push_back({node, true, openedKey});
      if (isGeometryLeaf(*node)) continue;
      for (auto child = node->children.rbegin(); child != node->children.rend(); ++child) {
        if (*child) pending.push_back({*child, false, nullptr});
      }
    }
  }

  void clear()
  {
    groups.clear();
    rewind();
  }

  /// Forget the current position; the next step starts at an end again.
  void rewind() { position = -1; }
  /// Number of groups (source lines), i.e. the number of navigation steps.
  [[nodiscard]] size_t size() const { return groups.size(); }
  [[nodiscard]] int currentPosition() const { return position; }
  [[nodiscard]] const Group *current() const
  {
    if (position < 0 || position >= static_cast<int>(groups.size())) return nullptr;
    return &groups[position];
  }

  /*!
    Move to the next (direction > 0) or previous (direction < 0) group and
    return it. Right after reset()/rewind(), +1 selects the first group and -1
    the last one. The position is clamped at both ends. direction == 0 returns
    the current group (or nullptr if there is none). The pointer is valid until
    the next reset()/clear().
   */
  const Group *step(int direction)
  {
    if (groups.empty()) return nullptr;
    if (position < 0) {
      if (direction == 0) return nullptr;
      position = direction < 0 ? static_cast<int>(groups.size()) - 1 : 0;
    } else if (direction > 0 && position + 1 < static_cast<int>(groups.size())) {
      ++position;
    } else if (direction < 0 && position > 0) {
      --position;
    }
    return &groups[position];
  }

private:
  using LineKey = std::pair<std::string, int>;  // file, first line

  static bool isObject(const AbstractNode& node)
  {
    return !isModifier(node) && node.modinst && !node.modinst->location().isNone();
  }

  static bool isModifier(const AbstractNode& node)
  {
    return dynamic_cast<const TransformNode *>(&node) || dynamic_cast<const ColorNode *>(&node);
  }

  static bool isGeometryLeaf(const AbstractNode& node)
  {
    return dynamic_cast<const AbstractPolyNode *>(&node) || dynamic_cast<const RenderNode *>(&node) ||
           dynamic_cast<const CgalAdvNode *>(&node);
  }

  std::vector<Group> groups;
  int position{-1};
};
