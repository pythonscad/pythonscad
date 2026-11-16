/* This is a collection of functions to share among all new languages */
#pragma once

#include "genlang/language.h"

#include <memory>
#include "core/node.h"
#include "core/function.h"
#include "core/ScopeContext.h"
#include "core/UserModule.h"
#include "core/CsgOpNode.h"

// Source location of the call that is creating a node right now. Each language runtime
// (e.g. Python) registers a provider that inspects its own call stack; without one, or
// outside of any script frame, Location::NONE is returned.
using GenlangLocationProvider = Location (*)();
extern GenlangLocationProvider genlang_location_provider;
Location genlang_current_location();

#define DECLARE_INSTANCE()     \
  std::string instance_name;   \
  AssignmentList inst_asslist; \
  auto instance =              \
    std::make_shared<ModuleInstantiation>(instance_name, inst_asslist, genlang_current_location());

extern std::vector<std::string> mapping_name;
extern std::vector<std::string> mapping_code;
extern std::vector<int> mapping_level;

void show_final(void);  // this is called when the new language terminates
extern std::vector<std::shared_ptr<AbstractNode>> shows;
extern std::shared_ptr<AbstractNode> genlang_result_node;
