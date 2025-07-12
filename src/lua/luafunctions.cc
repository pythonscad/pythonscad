/*
 *  OpenSCAD (www.openscad.org)
 *  Copyright (C) 2009-2011 Clifford Wolf <clifford@clifford.at> and
 *                          Marius Kintel <marius@kintel.net>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  As a special exception, you have permission to link this program
 *  with the CGAL library and distribute executables, as long as you
 *  follow the requirements of the GNU GPL in regard to all of the
 *  software in the executable aside from CGAL.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */


#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include "linalg.h"
#include "export.h"
#include "GeometryUtils.h"
#include "SourceFile.h"
#include "BuiltinContext.h"
#include <PolySetBuilder.h>
extern bool parse(SourceFile *& file, const std::string& text, const std::string& filename, const std::string& mainFile, int debug);

#include "GeometryUtils.h"
#include "primitives.h"
#include "TransformNode.h"
#include "RotateExtrudeNode.h"
#include "LinearExtrudeNode.h"
#include "PathExtrudeNode.h"
#include "PullNode.h"
#include "WrapNode.h"
#include "OversampleNode.h"
#include "FilletNode.h"
#include "CgalAdvNode.h"
#include "CsgOpNode.h"
#include "ColorNode.h"
#include "Expression.h"
#include "RoofNode.h"
#include "RenderNode.h"
#include "SurfaceNode.h"
#include "TextNode.h"
#include "OffsetNode.h"
#include <hash.h>
#include <PolySetUtils.h>
#include "ProjectionNode.h"
#include "ImportNode.h"
#include <Tree.h>
#include <GeometryEvaluator.h>

#include "degree_trig.h"
#include "printutils.h"
#include "io/fileutils.h"
#include "handle_dep.h"
#include "luaopenscad.h"
#include "genlang/genlang.h"
#include "stdarg.h"

//using namespace boost::assign; // bring 'operator+=()' into scope

#define TAG_NODE "node"


// Colors extracted from https://drafts.csswg.org/css-color/ on 2015-08-02
// CSS Color Module Level 4 - Editor’s Draft, 29 May 2015
extern std::unordered_map<std::string, Color4f> webcolors;
extern boost::optional<Color4f> parse_hex_color(const std::string& hex);

static int lua_cube(lua_State *L)
{
  DECLARE_INSTANCE

  char *kwlist[] = {"x", "y", "z",NULL};
  void *size = NULL;
  double x=10,y=10,z=10;

  void *center = NULL;

  if (!LuArg_ParseTupleAndKeywords(L, "ddd", kwlist, &x,&y,&z)){
    printf("Error during cube\n");	  
    return 0;
  }       
	 
  auto node = std::make_shared<CubeNode>(instance);
 
  node->dim[0]=x;	  
  node->dim[1]=y;	  
  node->dim[2]=z;	  
  for(int i=0;i<3;i++) { // TODO fix
    node->center[i]=0;
  }

  return LuaOpenSCADObjectFromNode( node);
}

static int lua_sphere(lua_State *L)
{
  DECLARE_INSTANCE
  auto node = std::make_shared<SphereNode>(instance);

  char *kwlist[] = {"r",  NULL};
  double r = NAN;
  double d = NAN;
  double fn = NAN, fa = NAN, fs = NAN;

  double vr = 1;

  if (!LuArg_ParseTupleAndKeywords(L, "dd", kwlist, &r)) {
    printf( "Error during parsing sphere(r|d)");
    return 1;
  } 
    if(r <= 0) {
      printf("Parameter r must be positive");
      return 0;
    }       
    vr = r;
    if(!std::isnan(d)) {
      printf("Cant specify r and d at the same time for sphere");
      return 1;
    }
  if (!std::isnan(d)) {
    if(d <= 0) {
      printf( "Parameter d must be positive");
      return 1;
    }       
    vr = d / 2.0;
  }

  lua_get_fnas(node->fn, node->fa, node->fs);
  if (!std::isnan(fn)) node->fn = fn;
  if (!std::isnan(fa)) node->fa = fa;
  if (!std::isnan(fs)) node->fs = fs;

  node->r = vr;

  return LuaOpenSCADObjectFromNode( node);
}


static int lua_show(lua_State *J)
{
  // TODO type and argnum checking	
  auto child = LuaOpenSCADObjectToNode(J, 1);	
  if(child == nullptr) return 0;
  shows.push_back(child);
  return 0;	
}



#define JS_ADDFUNCTION(name) \
  js_newcfunction(js_interp, js_##name, #name, 1); js_setglobal(js_interp, #name);

void registerLuaFunctions(void) {
  lua_register(L,"cube", lua_cube);
  lua_register(L,"sphere", lua_sphere);
  lua_register(L,"show", lua_show);
#if 0	
  JS_ADDFUNCTION(print)	
  JS_ADDFUNCTION(cube)
  JS_ADDFUNCTION(sphere)
  JS_ADDFUNCTION(cylinder)
  JS_ADDFUNCTION(difference)
  JS_ADDFUNCTION(intersection)
  JS_ADDFUNCTION(translate)
  JS_ADDFUNCTION(rotate)
  JS_ADDFUNCTION(output)
#endif  
}

