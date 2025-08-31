/* This is a colletion of function to share among all new languages */

#include <memory>
#include "core/node.h"
#include "core/function.h"
#include "core/ScopeContext.h"
#include "core/UserModule.h"
#include "core/CsgOpNode.h"

enum { LANG_NONE = -1, LANG_SCAD = 0, LANG_PYTHON = 1, LANG_JS = 2, LANG_LUA = 3 };

#define DECLARE_INSTANCE                                                                                \
  std::string instance_name;                                                                            \
  AssignmentList inst_asslist;                                                                          \
  ModuleInstantiation *instance = new ModuleInstantiation(instance_name, inst_asslist, Location::NONE); \
  modinsts_list.push_back(instance);

extern std::vector<std::string> mapping_name;
extern std::vector<std::string> mapping_code;
extern std::vector<int> mapping_level;

class LanguageDesc
{
public:
  LanguageDesc(const char *suffix, int (*check_valid)(const char *, const char *))
  {
    this->suffix = suffix;
    this->check_validfile = check_valid;
  }
  const char *suffix;
  int (*check_validfile)(const char *filename, const char *content);
};

extern int language;
extern LanguageDesc *languageDesc[4];
void show_final(void);  // this is called when the new language terminates
extern std::vector<std::shared_ptr<AbstractNode>> shows;
extern std::shared_ptr<AbstractNode> genlang_result_node;
