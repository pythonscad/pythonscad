#include <memory>
#include "node.h"
#include <geometry/Polygon2d.h>
#include "src/core/function.h"
extern "C" {
#include<lua.h>
#include<lauxlib.h>
#include<lualib.h>
}

#pragma GCC diagnostic ignored "-Wwrite-strings"

typedef struct
{
  public:	
  int type_id; // 0 = AbstractNode		
  std::shared_ptr<AbstractNode> node;
} LuaOpenSCADObject;

extern lua_State  *L;
extern std::shared_ptr<AbstractNode> lua_result_node;
void registerLuaFunctions(void);
int LuaOpenSCADObjectFromNode(const std::shared_ptr<AbstractNode> &node);
std::shared_ptr<AbstractNode> LuaOpenSCADObjectToNode(lua_State *lua, int argnum);
void initLua(double time);
std::string evaluateLua(const std::string & code);
void finishLua(void);


