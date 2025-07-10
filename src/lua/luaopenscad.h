#include <memory>
#include "node.h"
#include <geometry/Polygon2d.h>
#include "src/core/function.h"

#pragma GCC diagnostic ignored "-Wwrite-strings"

#define DECLARE_INSTANCE	std::string instance_name; \
	AssignmentList inst_asslist;\
	ModuleInstantiation *instance = new ModuleInstantiation(instance_name,inst_asslist, Location::NONE);

//extern js_State *js_interp;
extern std::shared_ptr<AbstractNode> lua_result_node;
void registerLuaFunctions(void);
void LuaOpenSCADObjectFromNode(const std::shared_ptr<AbstractNode> &node);
std::shared_ptr<AbstractNode> LuaOpenSCADObjectToNode(void *data);
void initLua(double time);
std::string evaluateLua(const std::string & code);
void finishLua(void);


