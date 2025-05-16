#pragma once
#include "src/core/node.h"
#include "src/core/function.h"
#include "src/geometry/Polygon2d.h"
#include "src/core/Selection.h"

void initJs(double time);
void finishJs();
std::string evaluateJs(const std::string &code);

extern std::shared_ptr<AbstractNode> js_result_node;
