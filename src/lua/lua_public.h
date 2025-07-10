#pragma once
#include "src/core/node.h"
#include "src/core/function.h"
#include "src/geometry/Polygon2d.h"
#include "src/core/Selection.h"

void initLua(double time);
void finishLua();
std::string evaluateLua(const std::string &code);

extern std::shared_ptr<AbstractNode> lua_result_node;
