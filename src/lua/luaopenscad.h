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

int LuArg_ParseTupleAndKeywords(lua_State *L, const char *fmt, char **kwlist , ...);
void lua_get_fnas(double& fn, double& fa, double& fs)  ;

typedef struct
{
  public:	
  int type_id; // 0 = AbstractNode		
  std::shared_ptr<AbstractNode> node;
} LuaOpenSCADObject;

extern lua_State  *L;
void registerLuaFunctions(void);
int LuaOpenSCADObjectFromNode(const std::shared_ptr<AbstractNode> &node);
std::shared_ptr<AbstractNode> LuaOpenSCADObjectToNode(lua_State *lua, int argnum);
void initLua(double time);
std::string evaluateLua(const std::string & code);
void finishLua(void);


