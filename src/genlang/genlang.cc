#include "genlang.h"

#include "Feature.h"

std::vector<std::string> mapping_name;
std::vector<std::string> mapping_code;
std::vector<int> mapping_level;

std::vector<std::shared_ptr<AbstractNode>> shows;
std::shared_ptr<AbstractNode> genlang_result_node = nullptr;
void show_final(void)
{
  mapping_name.clear();
  mapping_code.clear();
  mapping_level.clear();
  if (shows.size() == 1) genlang_result_node = shows[0];
  else {
    DECLARE_INSTANCE();
    if (Feature::ExperimentalPythonSeparateObjects.is_enabled()) {
      /* A ListNode sitting at the root is unpacked into a GeometryList by the
       * geometry evaluator rather than being unioned, which is what lets the
       * 3MF and AMF writers emit one object per part. Every other format is
       * handed a fused solid by the exporter, so nothing downstream breaks. */
      genlang_result_node = std::make_shared<ListNode>(instance);
    } else {
      genlang_result_node = std::make_shared<CsgOpNode>(instance, OpenSCADOperator::UNION);
    }
    genlang_result_node->children = shows;
  }
  shows.clear();
}
