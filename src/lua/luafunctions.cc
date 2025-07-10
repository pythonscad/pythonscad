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

//using namespace boost::assign; // bring 'operator+=()' into scope

#define TAG_NODE "node"


// Colors extracted from https://drafts.csswg.org/css-color/ on 2015-08-02
// CSS Color Module Level 4 - Editor’s Draft, 29 May 2015
extern std::unordered_map<std::string, Color4f> webcolors;
extern boost::optional<Color4f> parse_hex_color(const std::string& hex);

#if 0
static void js_cube(js_State *J)
{
  DECLARE_INSTANCE
  std::shared_ptr<CubeNode>  node = std::make_shared<CubeNode>(instance);
  for (int i=0;i<3;i++) node->center[i]=0;
  node->dim[0]=1;
  node->dim[1]=1;
  node->dim[2]=1;
  if(js_isarray(J,1)) {
    if(js_hasindex(J, 1, 0)) {
      if(js_isnumber(J,-1)) node->dim[0] = js_tonumber(J, -1);
      js_pop(J,1);
    }
    if(js_hasindex(J, 1, 1)) {
      if(js_isnumber(J,-1)) node->dim[1] = js_tonumber(J, -1);
      js_pop(J,1);
    }
    if(js_hasindex(J, 1, 2)) {
      if(js_isnumber(J,-1)) node->dim[2] = js_tonumber(J, -1);
      js_pop(J,1);
    }
  } else if(js_isnumber(J,1)) {
    double size = js_tonumber(J, 1);
    node->dim[0] = size;
    node->dim[1] = size;
    node->dim[2] = size;
  }
  if(js_isobject(J,2)) {
    if(js_hasproperty(J, 2, "center")) {
      if(js_isboolean(J,-1)){
        for(int i=0;i<3;i++) 
          node->center[i] = js_toboolean(J, -1)?1:0;
      }
      js_pop(J,1);
    }
  }
//  char *kwlist[] = {"size", "center", NULL};
//  PyObject *size = NULL;

//  PyObject *center = NULL;



  JsOpenSCADObjectFromNode(node);
}


void get_fnas(js_State *J, int objindex, double &fn, double &fa, double &fs) {
  fn=0;
  fa=12;
  fs=2;
  if(js_isobject(J,objindex)) {
    if(js_hasproperty(J, objindex, "fn")) {
      if(js_isnumber(J,-1))
        fn = js_tonumber(J, -1);
      js_pop(J,1);
    }
    if(js_hasproperty(J, objindex, "fa")) {
      if(js_isnumber(J,-1))
        fa = js_tonumber(J, -1);
      js_pop(J,1);
    }
    if(js_hasproperty(J, objindex, "fs")) {
      if(js_isnumber(J,-1))
        fs = js_tonumber(J, -1);
      js_pop(J,1);
    }
  }  
}


static void js_sphere(js_State *J)
{
  DECLARE_INSTANCE
  auto node = std::make_shared<SphereNode>(instance);

  node->r=1;
  if(js_isnumber(J,1)) node->r = js_tonumber(J, 1);
  get_fnas(J,2,node->fn, node->fa, node->fs);

//  js_retrieve_pyname(node);
  JsOpenSCADObjectFromNode(node);
}
static void js_cylinder(js_State *J)
{
  DECLARE_INSTANCE
  auto node = std::make_shared<CylinderNode>(instance);
  node->r1=1;
  node->r2=1;
  node->h=1;

  if(js_isnumber(J,1)) node->h = js_tonumber(J, 1);
  if(js_isnumber(J,2)) { node->r1 = js_tonumber(J, 2);  node->r2 = node->r1; }
  if(js_isobject(J,3)) {
    if(js_hasproperty(J, 3, "center")) {
      if(js_isboolean(J,-1))
        node->center = js_toboolean(J, -1)?1:0;
      js_pop(J,1);
    }
    if(js_hasproperty(J, 3, "angle")) {
      if(js_isnumber(J,-1))
        node->angle = js_tonumber(J, -1)?1:0;
      js_pop(J,1);
    }
    if(js_hasproperty(J, 3, "r2")) {
      if(js_isnumber(J,-1))
        node->r2 = js_tonumber(J, -1)?1:0;
      js_pop(J,1);
    }
  }
  get_fnas(J,3,node->fn, node->fa, node->fs);

  JsOpenSCADObjectFromNode(node);
}


void js_rotate(js_State *J)
{
  DECLARE_INSTANCE
  auto node = std::make_shared<TransformNode>(instance, "rotate");

  if(js_isuserdata(J,1,TAG_NODE)){
    void *data = js_touserdata(J, 1, TAG_NODE);
    std::shared_ptr<AbstractNode> child = JsOpenSCADObjectToNode(data);
    node->children.push_back(child);
  }
  Matrix3d M;
  Vector3d vec3(0,0,0);	
  if(js_isarray(J,2)) {
    if(js_hasindex(J, 2, 0)) {
      if(js_isnumber(J,-1)) vec3[0] = js_tonumber(J, -1);
      js_pop(J,1);
    }
    if(js_hasindex(J, 2, 1)) {
      if(js_isnumber(J,-1)) vec3[1] = js_tonumber(J, -1);
      js_pop(J,1);
    }
    if(js_hasindex(J, 2, 2)) {
      if(js_isnumber(J,-1)) vec3[2] = js_tonumber(J, -1);
      js_pop(J,1);
    }

    double sx = 0, sy = 0, sz = 0;
    double cx = 1, cy = 1, cz = 1;
    double a = 0.0;
    bool ok = true;
    if(vec3[2] != 0) {
      a = vec3[2];
      sz = sin_degrees(a);
      cz = cos_degrees(a);
    }
    if(vec3[1] != 0) {
      a = vec3[1];
      sy = sin_degrees(a);
      cy = cos_degrees(a);
    }
    if(vec3[0] != 0) {
      a = vec3[0];
      sx = sin_degrees(a);
      cx = cos_degrees(a);
    }

    M << cy * cz,  cz *sx *sy - cx * sz,   cx *cz *sy + sx * sz,
        cy *sz,  cx *cz + sx * sy * sz,  -cz * sx + cx * sy * sz,
      -sy,       cy *sx,                  cx *cy;
    node->matrix.rotate(M);
  }
  if(js_isnumber(J,2)) {
    double angle=js_tonumber(J,2);	  
    if(js_isarray(J,3)) {
      if(js_hasindex(J, 3, 0)) {
        if(js_isnumber(J,-1)) vec3[0] = js_tonumber(J, -1);
        js_pop(J,1);
      }
      if(js_hasindex(J, 3, 1)) {
        if(js_isnumber(J,-1)) vec3[1] = js_tonumber(J, -1);
        js_pop(J,1);
      }
      if(js_hasindex(J, 3, 2)) {
        if(js_isnumber(J,-1)) vec3[2] = js_tonumber(J, -1);
        js_pop(J,1);
      }
    }  
    M = angle_axis_degrees(angle, vec3);
    node->matrix.rotate(M);

  }

//  node->setPyName(child->getPyName());

  JsOpenSCADObjectFromNode(node);
}

void js_translate(js_State *J)
{
  DECLARE_INSTANCE
  auto node = std::make_shared<TransformNode>(instance, "translate");
  if(js_isuserdata(J,1,TAG_NODE)){
    void *data = js_touserdata(J, 1, TAG_NODE);
    std::shared_ptr<AbstractNode> child = JsOpenSCADObjectToNode(data);
    node->children.push_back(child);
  }

  Vector3d vec3(0,0,0);	
  if(js_isarray(J,2)) {
    if(js_hasindex(J, 2, 0)) {
      if(js_isnumber(J,-1)) vec3[0] = js_tonumber(J, -1);
      js_pop(J,1);
    }
    if(js_hasindex(J, 2, 1)) {
      if(js_isnumber(J,-1)) vec3[1] = js_tonumber(J, -1);
      js_pop(J,1);
    }
    if(js_hasindex(J, 2, 2)) {
      if(js_isnumber(J,-1)) vec3[2] = js_tonumber(J, -1);
      js_pop(J,1);
    }
  }

  node->matrix.translate(vec3);
  JsOpenSCADObjectFromNode(node);
}



static void js_output(js_State *J)
{
//  PyObject *child_dict;
  if(js_isuserdata(J,1,TAG_NODE)){
    void *data = js_touserdata(J, 1, TAG_NODE);
    js_result_node = JsOpenSCADObjectToNode(data);
  }
	

}
void js_csg_sub(js_State *J, OpenSCADOperator mode)
{
  DECLARE_INSTANCE
  int i;

  auto node = std::make_shared<CsgOpNode>(instance, mode);
  node->r=0;
  node->fn=1;
  int n=1;
  while(1) {
    if(!js_isuserdata(J,n,TAG_NODE)) break;
    void *data = js_touserdata(J, n, TAG_NODE);
    std::shared_ptr<AbstractNode> child = JsOpenSCADObjectToNode(data);
    node->children.push_back(child);
    n++;
  }
  JsOpenSCADObjectFromNode(node);
}

static void js_union(js_State *J)
{
  return js_csg_sub(J, OpenSCADOperator::UNION);
}

static void js_difference(js_State *J)
{
  return js_csg_sub(J, OpenSCADOperator::DIFFERENCE);
}

static void js_intersection(js_State *J)
{
  return js_csg_sub(J, OpenSCADOperator::INTERSECTION);
}

#endif

#define JS_ADDFUNCTION(name) \
  js_newcfunction(js_interp, js_##name, #name, 1); js_setglobal(js_interp, #name);

void registerLuaFunctions(void) {
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

