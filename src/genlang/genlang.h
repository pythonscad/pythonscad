/* This is a collection of functions to share among all new languages */
#pragma once

#include "genlang/language.h"

#include <memory>
#include "core/node.h"
#include "core/function.h"
#include "core/ScopeContext.h"
#include "core/UserModule.h"
#include "core/CsgOpNode.h"

#define DECLARE_INSTANCE()                                                                              \
  std::string instance_name;                                                                            \
  AssignmentList inst_asslist;                                                                          \
  ModuleInstantiation *instance = new ModuleInstantiation(instance_name, inst_asslist, Location::NONE); \
  modinsts_list.push_back(instance);

typedef struct {
  std::string name;
  std::string code;
  int level;
} PythonName;
extern std::vector<PythonName> pythonName;
extern std::vector<PythonName> pythonNameReady;

void show_final(void);  // this is called when the new language terminates
extern std::vector<std::shared_ptr<AbstractNode>> shows;
extern std::shared_ptr<AbstractNode> genlang_result_node;
