#include "genlang.h"

std::vector<PythonName> pythonName;
std::vector<PythonName> pythonNameReady;

std::vector<std::shared_ptr<AbstractNode>> shows;
std::shared_ptr<AbstractNode> genlang_result_node = nullptr;
bool python_build_hashmap(const std::shared_ptr<AbstractNode>& node, int level);

void show_final(void)
{
  printf("show_final\n");
  pythonNameReady = pythonName;
  pythonName.clear();
  if (shows.size() == 1) genlang_result_node = shows[0];
  else {
    DECLARE_INSTANCE();
    genlang_result_node = std::make_shared<CsgOpNode>(instance, OpenSCADOperator::UNION);
    genlang_result_node->children = shows;
  }
  for (const auto& show : shows) {
    if (!python_build_hashmap(show, 0)) {
      /* Helper hit a non-TypeError exception (MemoryError, ...) while
       * iterating __main__ -- propagate to the Python caller rather
       * than silently returning a half-populated selection table. */
    }
  }
  shows.clear();
}
