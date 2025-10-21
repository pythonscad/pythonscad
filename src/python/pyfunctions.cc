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
#include <Python.h>
#include "pyfunctions.h"
#include "python/pyopenscad.h"
#include "core/primitives.h"
#include "core/CsgOpNode.h"
#include "core/ColorNode.h"
#include "core/ColorUtil.h"
#include "SourceFile.h"
#include "BuiltinContext.h"
#include <PolySetBuilder.h>
#include "io/export.h"
extern bool parse(SourceFile *& file, const std::string& text, const std::string& filename,
                  const std::string& mainFile, int debug);

#include <python/pydata.h>
#ifdef ENABLE_LIBFIVE
#include "python/FrepNode.h"
#endif
#include "GeometryUtils.h"
#include "core/TransformNode.h"
#include "core/LinearExtrudeNode.h"
#include "core/RotateExtrudeNode.h"
#include "core/PathExtrudeNode.h"
#include "core/PullNode.h"
#include "core/WrapNode.h"
#include "core/OversampleNode.h"
#include "core/DebugNode.h"
#include "core/RepairNode.h"
#include "core/FilletNode.h"
#include "core/SkinNode.h"
#include "core/ConcatNode.h"
#include "core/CgalAdvNode.h"
#include "Expression.h"
#include "core/RoofNode.h"
#include "core/RenderNode.h"
#include "core/SurfaceNode.h"
#include "core/SheetNode.h"
#include "core/TextNode.h"
#include "core/OffsetNode.h"
#include <hash.h>
#include "geometry/PolySetUtils.h"
#include "core/ProjectionNode.h"
#include "core/ImportNode.h"
#include "core/Tree.h"
#include "geometry/PolySet.h"
#include "geometry/GeometryEvaluator.h"
#include "utils/degree_trig.h"
#include "printutils.h"
#include "io/fileutils.h"
#include "handle_dep.h"
#include <fstream>
#include <ostream>
#include <boost/functional/hash.hpp>
#include <ScopeContext.h>
#include "PlatformUtils.h"
#include <iostream>
#include <filesystem>

// using namespace boost::assign; // bring 'operator+=()' into scope

// Colors extracted from https://drafts.csswg.org/css-color/ on 2015-08-02
// CSS Color Module Level 4 - Editor’s Draft, 29 May 2015
extern std::unordered_map<std::string, Color4f> webcolors;

PyObject *python_edge(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE
  auto node = std::make_shared<EdgeNode>(instance);

  char *kwlist[] = {"size", "center", NULL};
  double size = 1;

  PyObject *center = NULL;

  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "|dO", kwlist, &size, &center)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing edge(size)");
    return NULL;
  }

  if (size < 0) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Edge Length must be positive");
    return NULL;
  }
  node->size = size;
  if (center == Py_FALSE || center == NULL)
    ;
  else if (center == Py_TRUE) {
    node->center = true;
  }
  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

PyObject *python_marked(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE
  auto node = std::make_shared<EdgeNode>(instance);

  char *kwlist[] = {"value", NULL};
  double value = 0.0;

  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "d", kwlist, &value)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing marked(value)");
    return NULL;
  }
  return PyDataObjectFromValue(&PyDataType, value);
}

PyObject *python_cube(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE
  auto node = std::make_shared<CubeNode>(instance);

  char *kwlist[] = {"size", "center", NULL};
  PyObject *size = NULL;

  PyObject *center = NULL;

  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "|OO", kwlist, &size, &center)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing cube(size)");
    return NULL;
  }

  if (size != NULL) {
    int flags = 0;
    if (python_vectorval(size, 3, 3, &(node->dim[0]), &(node->dim[1]), &(node->dim[2]), nullptr,
                         &(node->dragflags))) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid Cube dimensions");
      return NULL;
    }
  }
  if (node->dim[0] <= 0 || node->dim[1] <= 0 || node->dim[2] <= 0) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Cube Dimensions must be positive");
    return NULL;
  }
  for (int i = 0; i < 3; i++) node->center[i] = 1;
  if (center == Py_FALSE || center == NULL)
    ;
  else if (center == Py_TRUE) {
    for (int i = 0; i < 3; i++) node->center[i] = 0;
  } else if (PyUnicode_Check(center)) {
    PyObject *centerval = pf.PyUnicode_AsEncodedString(center, "utf-8", "~");
    const char *centerstr = PyBytes_AS_STRING(centerval);
    if (centerstr == nullptr) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "Cannot parse center code");
      return NULL;
    }
    if (strlen(centerstr) != 3) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "Center code must be exactly 3 characters");
      return NULL;
    }
    for (int i = 0; i < 3; i++) {
      switch (centerstr[i]) {
      case '|': node->center[i] = 0; break;
      case '0': node->center[i] = 0; break;
      case ' ': node->center[i] = 0; break;
      case '_': node->center[i] = 0; break;

      case '>': node->center[i] = -1; break;
      case ']': node->center[i] = -1; break;
      case ')': node->center[i] = -1; break;
      case '+': node->center[i] = -1; break;

      case '<': node->center[i] = 1; break;
      case '[': node->center[i] = 1; break;
      case '(': node->center[i] = 1; break;
      case '-': node->center[i] = 1; break;

      default:
        pf.PyErr_SetString(pf.PyExc_TypeError, "Center code chars not recognized, must be + - or 0");
        return NULL;
      }
    }
  } else {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Unknown Value for center parameter");
    return NULL;
  }
  python_retrieve_pyname(node);
  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

int sphereCalcIndInt(PyObject *func, Vector3d& dir)
{
  dir.normalize();
  PyObject *dir_p = pf.PyList_New(3);
  for (int i = 0; i < 3; i++) pf.PyList_SetItem(dir_p, i, pf.PyFloat_FromDouble(dir[i]));
  PyObject *args = pf.PyTuple_Pack(1, dir_p);
  PyObject *len_p = pf.PyObject_CallObject(func, args);
  double len = 0;
  if (len_p == nullptr) {
    std::string errorstr;
    python_catch_error(errorstr);
    pf.PyErr_SetString(pf.PyExc_TypeError, errorstr.c_str());
    LOG(message_group::Error, errorstr.c_str());
    return 1;
  }
  python_numberval(len_p, &len);
  dir *= len;
  return 0;
}

int sphereCalcInd(PolySetBuilder& builder, std::vector<Vector3d>& vertices, PyObject *func, Vector3d dir)
{
  std::string errorstr;
  if (sphereCalcIndInt(func, dir)) return -1;  // TODO fix
  unsigned int ind = builder.vertexIndex(dir);
  if (ind == vertices.size()) vertices.push_back(dir);
  return ind;
}

class SphereEdgeDb
{
public:
  SphereEdgeDb(int a, int b)
  {
    ind1 = (a < b) ? a : b;
    ind2 = (a > b) ? a : b;
  }
  int ind1, ind2;
  int operator==(const SphereEdgeDb ref)
  {
    if (this->ind1 == ref.ind1 && this->ind2 == ref.ind2) return 1;
    return 0;
  }
};

unsigned int hash_value(const SphereEdgeDb& r)
{
  unsigned int i;
  i = r.ind1 | (r.ind2 << 16);
  return i;
}

int operator==(const SphereEdgeDb& t1, const SphereEdgeDb& t2)
{
  if (t1.ind1 == t2.ind1 && t1.ind2 == t2.ind2) return 1;
  return 0;
}

int sphereCalcSplitInd(PolySetBuilder& builder, std::vector<Vector3d>& vertices,
                       std::unordered_map<SphereEdgeDb, int, boost::hash<SphereEdgeDb>>& edges,
                       PyObject *func, int ind1, int ind2)
{
  SphereEdgeDb edge(ind1, ind2);
  if (edges.count(edge) > 0) {
    return edges[edge];
  }
  int result = sphereCalcInd(builder, vertices, func, vertices[ind1] + vertices[ind2]);
  if (result != -1) edges[edge] = result;
  return result;
}

std::unique_ptr<const Geometry> sphereCreateFuncGeometry(void *funcptr, double fs, int n)
{
  PyObject *func = (PyObject *)funcptr;
  std::unordered_map<SphereEdgeDb, int, boost::hash<SphereEdgeDb>> edges;

  PolySetBuilder builder;
  std::vector<Vector3d> vertices;

  int topind, botind, leftind, rightind, frontind, backind;
  leftind = sphereCalcInd(builder, vertices, func, Vector3d(-1, 0, 0));
  if (leftind < 0) return builder.build();
  rightind = sphereCalcInd(builder, vertices, func, Vector3d(1, 0, 0));
  frontind = sphereCalcInd(builder, vertices, func, Vector3d(0, -1, 0));
  backind = sphereCalcInd(builder, vertices, func, Vector3d(0, 1, 0));
  botind = sphereCalcInd(builder, vertices, func, Vector3d(0, 0, -1));
  topind = sphereCalcInd(builder, vertices, func, Vector3d(0, 0, 1));
  if (rightind < 0 || frontind < 0 || backind < 0 || botind < 0 || topind < 0) return builder.build();

  std::vector<IndexedTriangle> triangles;
  std::vector<IndexedTriangle> tri_new;
  tri_new.push_back(IndexedTriangle(leftind, frontind, topind));
  tri_new.push_back(IndexedTriangle(frontind, rightind, topind));
  tri_new.push_back(IndexedTriangle(rightind, backind, topind));
  tri_new.push_back(IndexedTriangle(backind, leftind, topind));
  tri_new.push_back(IndexedTriangle(leftind, botind, frontind));
  tri_new.push_back(IndexedTriangle(frontind, botind, rightind));
  tri_new.push_back(IndexedTriangle(rightind, botind, backind));
  tri_new.push_back(IndexedTriangle(backind, botind, leftind));

  int round = 0;
  unsigned int i1, i2, imid;
  Vector3d p1, p2, p3, pmin, pmax, pmid, pmid_test, dir1, dir2;
  double dist, ang, ang_test;
  do {
    triangles = tri_new;
    if (round == n) break;
    tri_new.clear();
    std::vector<int> midinds;
    for (const IndexedTriangle& tri : triangles) {
      int zeroang = 0;
      unsigned int midind = -1;
      for (int i = 0; i < 3; i++) {
        i1 = tri[i];
        i2 = tri[(i + 1) % 3];
        SphereEdgeDb edge(i1, i2);
        if (edges.count(edge) > 0) continue;
        dist = (vertices[i1] - vertices[i2]).norm();
        if (dist < fs) continue;
        p1 = vertices[i1];
        p2 = vertices[i2];
        pmin = p1;
        pmax = p2;
        pmid = (pmin + pmax) / 2;
        if (sphereCalcIndInt(func, pmid)) return builder.build();
        dir1 = (pmid - p1).normalized();
        dir2 = (p2 - pmid).normalized();
        ang = acos(dir1.dot(dir2));
        //	printf("ang=%g\n",ang*180/3.14);
        if (ang < 0.001) {
          zeroang++;
          continue;
        }
        do {
          pmid_test = (pmin + pmid) / 2;
          if (sphereCalcIndInt(func, pmid_test)) return builder.build();
          dir1 = (pmid_test - p1).normalized();
          dir2 = (p2 - pmid_test).normalized();
          ang_test = acos(dir1.dot(dir2));
          if (ang_test > ang) {
            pmax = pmid;
            ang = ang_test;
            pmid = pmid_test;
            if ((pmax - pmin).norm() > fs) continue;
          }

          pmid_test = (pmax + pmid) / 2;
          if (sphereCalcIndInt(func, pmid_test)) return builder.build();
          dir1 = (pmid_test - p1).normalized();
          dir2 = (p2 - pmid_test).normalized();
          ang_test = acos(dir1.dot(dir2));
          if (ang_test > ang) {
            pmin = pmid;
            ang = ang_test;
            pmid = pmid_test;
            if ((pmax - pmin).norm() > fs) continue;
          }
        } while (0);

        imid = builder.vertexIndex(pmid);
        if (imid == vertices.size()) vertices.push_back(pmid);
        edges[edge] = imid;
      }
      if (zeroang == 3) {
        p1 = vertices[tri[0]];
        p2 = vertices[tri[1]];
        p3 = vertices[tri[2]];
        pmid = (p1 + p2 + p3) / 3.0;
        if (sphereCalcIndInt(func, pmid)) return builder.build();
        Vector4d norm = calcTriangleNormal(vertices, {tri[0], tri[1], tri[2]});
        if (fabs(pmid.dot(norm.head<3>()) - norm[3]) > 1e-3) {
          midind = builder.vertexIndex(pmid);
          if (midind == vertices.size()) vertices.push_back(pmid);
        }
      }
      midinds.push_back(midind);
    }
    // create new triangles from split edges
    int ind = 0;
    for (const IndexedTriangle& tri : triangles) {
      int splitind[3];
      for (int i = 0; i < 3; i++) {
        SphereEdgeDb e(tri[i], tri[(i + 1) % 3]);
        splitind[i] = edges.count(e) > 0 ? edges[e] : -1;
      }

      if (midinds[ind] != -1) {
        for (int i = 0; i < 3; i++) {
          if (splitind[i] == -1) {
            tri_new.push_back(IndexedTriangle(tri[i], tri[(i + 1) % 3], midinds[ind]));
          } else {
            tri_new.push_back(IndexedTriangle(tri[i], splitind[i], midinds[ind]));
            tri_new.push_back(IndexedTriangle(splitind[i], tri[(i + 1) % 3], midinds[ind]));
          }
        }
        ind++;
        continue;
      }

      int bucket =
        ((splitind[0] != -1) ? 1 : 0) | ((splitind[1] != -1) ? 2 : 0) | ((splitind[2] != -1) ? 4 : 0);
      switch (bucket) {
      case 0: tri_new.push_back(IndexedTriangle(tri[0], tri[1], tri[2])); break;
      case 1:
        tri_new.push_back(IndexedTriangle(tri[0], splitind[0], tri[2]));
        tri_new.push_back(IndexedTriangle(tri[2], splitind[0], tri[1]));
        break;
      case 2:
        tri_new.push_back(IndexedTriangle(tri[1], splitind[1], tri[0]));
        tri_new.push_back(IndexedTriangle(tri[0], splitind[1], tri[2]));
        break;
      case 3:
        tri_new.push_back(IndexedTriangle(tri[0], splitind[0], tri[2]));
        tri_new.push_back(IndexedTriangle(splitind[0], splitind[1], tri[2]));
        tri_new.push_back(IndexedTriangle(splitind[0], tri[1], splitind[1]));
        break;
      case 4:
        tri_new.push_back(IndexedTriangle(tri[2], splitind[2], tri[1]));
        tri_new.push_back(IndexedTriangle(tri[1], splitind[2], tri[0]));
        break;
      case 5:
        tri_new.push_back(IndexedTriangle(tri[0], splitind[0], splitind[2]));
        tri_new.push_back(IndexedTriangle(splitind[0], tri[2], splitind[2]));
        tri_new.push_back(IndexedTriangle(splitind[0], tri[1], tri[2]));
        break;
      case 6:
        tri_new.push_back(IndexedTriangle(tri[0], tri[1], splitind[2]));
        tri_new.push_back(IndexedTriangle(tri[1], splitind[1], splitind[2]));
        tri_new.push_back(IndexedTriangle(splitind[1], tri[2], splitind[2]));
        break;
      case 7:
        tri_new.push_back(IndexedTriangle(splitind[2], tri[0], splitind[0]));
        tri_new.push_back(IndexedTriangle(splitind[0], tri[1], splitind[1]));
        tri_new.push_back(IndexedTriangle(splitind[1], tri[2], splitind[2]));
        tri_new.push_back(IndexedTriangle(splitind[0], splitind[1], splitind[2]));
        break;
      }
      ind++;
    }

    round++;
  } while (tri_new.size() != triangles.size());
  for (const IndexedTriangle& tri : tri_new) {
    builder.appendPolygon({tri[0], tri[1], tri[2]});
  }
  auto ps = builder.build();

  int done = 0;
  round = 0;
  do {
    done = 0;
    auto edge_db = createEdgeDb(ps->indices);
    for (size_t i = 0; i < ps->indices.size(); i++) {
      auto& tri = ps->indices[i];
      if (tri[0] == tri[1] || tri[0] == tri[2] || tri[1] == tri[2]) continue;
      for (int j = 0; j < 3; j++) {
        int i1 = tri[j];
        int i2 = tri[(j + 1) % 3];
        double l1 = (ps->vertices[i1] - ps->vertices[i2]).norm();
        EdgeKey ek(tri[j], tri[(j + 1) % 3]);
        if (edge_db.count(ek) != 0) {
          auto ev = edge_db.at(ek);
          int face_o, pos_o;
          if (i2 > i1) {
            face_o = ev.faceb;
            pos_o = ev.posb;
          } else {
            face_o = ev.facea;
            pos_o = ev.posa;
          }
          if (face_o == -1 || pos_o == -1) continue;
          auto& tri_oth = ps->indices[face_o];
          double l2 = (ps->vertices[tri[(j + 2) % 3]] - ps->vertices[tri_oth[(pos_o + 2) % 3]]).norm();
          if (l2 < l1) {
            Vector3d norm = calcTriangleNormal(ps->vertices, tri).head<3>();
            Vector3d norm_oth = calcTriangleNormal(ps->vertices, tri_oth).head<3>();

            auto tri_ = tri;
            auto tri_oth_ = tri_oth;

            tri_[(j + 1) % 3] = tri_oth[(pos_o + 2) % 3];
            for (int k = 0; k < 3; k++)
              if (tri_oth[k] == i1) tri_oth_[k] = tri[(j + 2) % 3];
            // reorganize

            Vector3d norm_ = calcTriangleNormal(ps->vertices, tri_).head<3>();
            Vector3d norm_oth_ = calcTriangleNormal(ps->vertices, tri_oth_).head<3>();

            if (norm.dot(norm_) > 0 && norm_oth.dot(norm_oth_) > 0) {
              tri = tri_;
              tri_oth = tri_oth_;

              for (int k = 0; k < 3; k++) {
                edge_db.erase(EdgeKey(tri[k], tri[(k + 1) % 3]));
                edge_db.erase(EdgeKey(tri_oth[k], tri_oth[(k + 1) % 3]));
              }
              done++;
              break;  // dont proceed with
            }
          }
        }
      }
    }
    for (size_t i = 0; i < ps->indices.size(); i++) {
      auto& tri = ps->indices[i];
      if (tri[0] == tri[1] && tri[0] == tri[2]) {
        ps->indices.erase(ps->indices.begin() + i);
        i--;
      }
    }
    round++;
  } while (done > 0);  //  && round < 3);

  return ps;
}

PyObject *python_sphere(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE
  auto node = std::make_shared<SphereNode>(instance);

  char *kwlist[] = {"r", "d", "fn", "fa", "fs", NULL};
  double r = NAN;
  PyObject *rp = nullptr;
  double d = NAN;
  double fn = NAN, fa = NAN, fs = NAN;

  double vr = 1;

  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "|Odddd", kwlist, &rp, &d, &fn, &fa, &fs)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing sphere(r|d)");
    return NULL;
  }
  if (rp != nullptr) {
    if (python_numberval(rp, &r, &(node->dragflags), 1))
      if (rp->ob_type == pf.PyFunction_Type) node->r_func = rp;
  }
  if (!isnan(r)) {
    if (r <= 0) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "Parameter r must be positive");
      return NULL;
    }
    vr = r;
    if (!isnan(d)) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "Cant specify r and d at the same time for sphere");
      return NULL;
    }
  }
  if (!isnan(d)) {
    if (d <= 0) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "Parameter d must be positive");
      return NULL;
    }
    vr = d / 2.0;
  }

  get_fnas(node->fn, node->fa, node->fs);
  if (!isnan(fn)) node->fn = fn;
  if (!isnan(fa)) node->fa = fa;
  if (!isnan(fs)) node->fs = fs;

  node->r = vr;

  python_retrieve_pyname(node);
  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

PyObject *python_cylinder(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE
  auto node = std::make_shared<CylinderNode>(instance);

  char *kwlist[] = {"h", "r1", "r2", "center", "r", "d", "d1", "d2", "angle", "fn", "fa", "fs", NULL};
  PyObject *h_ = nullptr;
  PyObject *r_ = nullptr;
  double r1 = NAN;
  double r2 = NAN;
  double d = NAN;
  double d1 = NAN;
  double d2 = NAN;
  double angle = NAN;

  double fn = NAN, fa = NAN, fs = NAN;

  PyObject *center = NULL;
  double vr1 = 1, vr2 = 1, vh = 1;

  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "|OddOOddddddd", kwlist, &h_, &r1, &r2, &center, &r_,
                                      &d, &d1, &d2, &angle, &fn, &fa, &fs)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing cylinder(h,r|r1+r2|d1+d2)");
    return NULL;
  }
  double r = NAN;
  double h = NAN;

  python_numberval(r_, &r, &(node->dragflags), 1);
  python_numberval(h_, &h, &(node->dragflags), 2);

  if (h <= 0) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Cylinder height must be positive");
    return NULL;
  }
  vh = h;

  if (!isnan(d) && d <= 0) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Cylinder d must be positive");
    return NULL;
  }
  if (!isnan(r1) && r1 < 0) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Cylinder r1 must not be negative");
    return NULL;
  }
  if (!isnan(r2) && r2 < 0) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Cylinder r2 must not be negative");
    return NULL;
  }
  if (!isnan(d1) && d1 < 0) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Cylinder d1 must not be negative");
    return NULL;
  }
  if (!isnan(d2) && d2 < 0) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Cylinder d2 must not be negative");
    return NULL;
  }

  if (!isnan(r1) && !isnan(r2)) {
    vr1 = r1;
    vr2 = r2;
  } else if (!isnan(r1) && isnan(r2)) {
    vr1 = r1;
    vr2 = r1;
  } else if (!isnan(d1) && !isnan(d2)) {
    vr1 = d1 / 2.0;
    vr2 = d2 / 2.0;
  } else if (!isnan(r)) {
    vr1 = r;
    vr2 = r;
  } else if (!isnan(d)) {
    vr1 = d / 2.0;
    vr2 = d / 2.0;
  }

  if (!isnan(angle)) node->angle = angle;
  get_fnas(node->fn, node->fa, node->fs);
  if (!isnan(fn)) node->fn = fn;
  if (!isnan(fa)) node->fa = fa;
  if (!isnan(fs)) node->fs = fs;

  node->r1 = vr1;
  node->r2 = vr2;
  node->h = vh;

  if (center == Py_TRUE) node->center = 1;
  else if (center == Py_FALSE || center == NULL) node->center = 0;
  else {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Unknown Value for center parameter");
    return NULL;
  }

  python_retrieve_pyname(node);
  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

PyObject *python_polyhedron(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE
  unsigned int i, j, pointIndex;
  auto node = std::make_shared<PolyhedronNode>(instance);

  char *kwlist[] = {"points", "faces", "convexity", "triangles", "colors", NULL};
  PyObject *points = NULL;
  PyObject *faces = NULL;
  int convexity = 2;
  PyObject *triangles = NULL;
  PyObject *colors = NULL;

  PyObject *element;
  Vector3d point;

  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!|iO!", kwlist, pf.PyList_Type, &points,
                                      pf.PyList_Type, &faces, &convexity, pf.PyList_Type, &triangles)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing polyhedron(points, faces)");
    return NULL;
  }

  if (points != NULL && PyList_CHECK(points)) {
    if (pf.PyList_Size(points) == 0) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "There must at least be one point in the polyhedron");
      return NULL;
    }
    for (i = 0; i < pf.PyList_Size(points); i++) {
      element = pf.PyList_GetItem(points, i);
      if (PyList_CHECK(element) && pf.PyList_Size(element) == 3) {
        point[0] = pf.PyFloat_AsDouble(pf.PyList_GetItem(element, 0));
        point[1] = pf.PyFloat_AsDouble(pf.PyList_GetItem(element, 1));
        point[2] = pf.PyFloat_AsDouble(pf.PyList_GetItem(element, 2));
        node->points.push_back(point);
      } else {
        pf.PyErr_SetString(pf.PyExc_TypeError, "Coordinate must exactly contain 3 numbers");
        return NULL;
      }
    }
  } else {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Polyhedron Points must be a list of coordinates");
    return NULL;
  }

  if (triangles != NULL) {
    faces = triangles;
    //	LOG(message_group::Deprecated, inst->location(), parameters.documentRoot(),
    //"polyhedron(triangles=[]) will be removed in future releases. Use polyhedron(faces=[]) instead.");
  }

  if (faces != NULL && PyList_CHECK(faces)) {
    if (pf.PyList_Size(faces) == 0) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "must specify at least 1 face");
      return NULL;
    }
    for (i = 0; i < pf.PyList_Size(faces); i++) {
      element = pf.PyList_GetItem(faces, i);
      if (PyList_CHECK(element)) {
        IndexedFace face;
        for (j = 0; j < pf.PyList_Size(element); j++) {
          pointIndex = pf.PyLong_AsLong(pf.PyList_GetItem(element, j));
          if (pointIndex < 0 || pointIndex >= node->points.size()) {
            pf.PyErr_SetString(pf.PyExc_TypeError, "Polyhedron Point Index out of range");
            return NULL;
          }
          face.push_back(pointIndex);
        }
        if (face.size() >= 3) {
          node->faces.push_back(std::move(face));
        } else {
          pf.PyErr_SetString(pf.PyExc_TypeError, "Polyhedron Face must sepcify at least 3 indices");
          return NULL;
        }

      } else {
        pf.PyErr_SetString(pf.PyExc_TypeError, "Polyhedron Face must be a list of indices");
        return NULL;
      }
    }
  } else {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Polyhedron faces must be a list of indices");
    return NULL;
  }

  if (colors != NULL && PyList_CHECK(colors)) {
    if ((size_t)pf.PyList_Size(colors) != node->faces.size()) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "when specified must match number of faces");
      return NULL;
    }
    for (i = 0; i < pf.PyList_Size(colors); i++) {
      element = pf.PyList_GetItem(colors, i);
      if (PyList_CHECK(element) && pf.PyList_Size(element) == 3) {
        Vector4f color(0, 0, 0, 1.0);
        for (j = 0; j < 3; j++) {
          color[j] = pf.PyFloat_AsDouble(PyList_GetItem(element, j));
        }
        int colind = -1;
        int ind = 0;
        for (auto& c : node->colors) {
          if (c == color) {
            colind = ind;
            break;
          }
          ind++;
        }
        if (colind == -1) {
          colind = node->colors.size();
          node->colors.push_back(color);
        }
        node->color_indices.push_back(colind);
      } else {
        pf.PyErr_SetString(pf.PyExc_TypeError, "Face Color must be a list with 3 values");
        return NULL;
      }
    }
  }

  node->convexity = convexity;
  if (node->convexity < 1) node->convexity = 1;

  python_retrieve_pyname(node);
  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

#ifdef ENABLE_LIBFIVE
PyObject *python_frep(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE
  auto node = std::make_shared<FrepNode>(instance);
  PyObject *expression = NULL;
  PyObject *bmin = NULL, *bmax = NULL;
  double res = 10;

  char *kwlist[] = {"exp", "min", "max", "res", NULL};

  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "OOO|d", kwlist, &expression, &bmin, &bmax, &res))
    return NULL;

  python_vectorval(bmin, 3, 3, &(node->x1), &(node->y1), &(node->z1));
  python_vectorval(bmax, 3, 3, &(node->x2), &(node->y2), &(node->z2));
  node->res = res;

  if (expression->ob_type == &PyDataType) {
    node->expression = expression;
  } else if (expression->ob_type == pf.PyFunction_Type) {
    node->expression = expression;
  } else {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Unknown frep expression type\n");
    return NULL;
  }

  python_retrieve_pyname(node);
  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

PyObject *python_ifrep(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE
  PyObject *object = NULL;
  PyObject *dummydict;

  char *kwlist[] = {"obj", nullptr};
  std::shared_ptr<AbstractNode> child;

  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O!", kwlist, &PyOpenSCADType, &object)) return NULL;

  child = PyOpenSCADObjectToNodeMulti(object, &dummydict);
  LeafNode *node = (LeafNode *)child.get();
  const std::shared_ptr<const Geometry> geom = node->createGeometry();
  const std::shared_ptr<const PolySet> ps = std::dynamic_pointer_cast<const PolySet>(geom);

  return ifrep(ps);
}

#endif

PyObject *python_square(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE
  auto node = std::make_shared<SquareNode>(instance);

  char *kwlist[] = {"dim", "center", NULL};
  PyObject *dim = NULL;

  PyObject *center = NULL;
  double z = NAN;

  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "|OO", kwlist, &dim, &center)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing square(dim)");
    return NULL;
  }
  if (dim != NULL) {
    if (python_vectorval(dim, 2, 2, &(node->x), &(node->y), &z)) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid Square dimensions");
      return NULL;
    }
  }
  if (center == Py_TRUE) node->center = 1;
  else if (center == Py_FALSE || center == NULL) node->center = 0;
  else {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Unknown Value for center parameter");
    return NULL;
  }

  python_retrieve_pyname(node);
  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

PyObject *python_circle(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE
  auto node = std::make_shared<CircleNode>(instance);

  char *kwlist[] = {"r", "d", "angle", "fn", "fa", "fs", NULL};
  double r = NAN;
  double d = NAN;
  double angle = NAN;
  double fn = NAN, fa = NAN, fs = NAN;

  double vr = 1;

  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "|dddddd", kwlist, &r, &d, &angle, &fn, &fa, &fs)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing circle(r|d)");
    return NULL;
  }

  if (!isnan(r)) {
    if (r <= 0) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "Parameter r must be positive");
      return NULL;
    }
    vr = r;
    if (!isnan(d)) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "Cant specify r and d at the same time for circle");
      return NULL;
    }
  }
  if (!isnan(d)) {
    if (d <= 0) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "Parameter d must be positive");
      return NULL;
    }
    vr = d / 2.0;
  }

  if (!isnan(angle)) node->angle = angle;

  get_fnas(node->fn, node->fa, node->fs);
  if (!isnan(fn)) node->fn = fn;
  if (!isnan(fa)) node->fa = fa;
  if (!isnan(fs)) node->fs = fs;

  node->r = vr;

  python_retrieve_pyname(node);
  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

PyObject *python_polygon(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE
  unsigned int i, j, pointIndex;
  auto node = std::make_shared<PolygonNode>(instance);

  char *kwlist[] = {"points", "paths", "convexity", NULL};
  PyObject *points = NULL;
  PyObject *paths = NULL;
  int convexity = 2;

  PyObject *element;
  Vector3d point;

  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O!|O!i", kwlist, pf.PyList_Type, &points,
                                      pf.PyList_Type, &paths, &convexity)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing polygon(points,paths)");
    return NULL;
  }

  if (points != NULL && PyList_CHECK(points)) {
    if (pf.PyList_Size(points) == 0) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "There must at least be one point in the polygon");
      return NULL;
    }
    for (i = 0; i < pf.PyList_Size(points); i++) {
      element = pf.PyList_GetItem(points, i);
      point[2] = 0;  // default no radius
      if (python_vectorval(element, 2, 3, &point[0], &point[1], &point[2])) {
        pf.PyErr_SetString(pf.PyExc_TypeError, "Coordinate must contain 2 or 3 numbers");
        return NULL;
      }
      node->points.push_back(point);
    }
  } else {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Polygon points must be a list of coordinates");
    return NULL;
  }

  if (paths != NULL && PyList_CHECK(paths)) {
    if (pf.PyList_Size(paths) == 0) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "must specify at least 1 path when specified");
      return NULL;
    }
    for (i = 0; i < pf.PyList_Size(paths); i++) {
      element = pf.PyList_GetItem(paths, i);
      if (PyList_CHECK(element)) {
        std::vector<size_t> path;
        for (j = 0; j < pf.PyList_Size(element); j++) {
          pointIndex = pf.PyLong_AsLong(pf.PyList_GetItem(element, j));
          if (pointIndex < 0 || pointIndex >= node->points.size()) {
            pf.PyErr_SetString(pf.PyExc_TypeError, "Polyhedron Point Index out of range");
            return NULL;
          }
          path.push_back(pointIndex);
        }
        node->paths.push_back(std::move(path));
      } else {
        pf.PyErr_SetString(pf.PyExc_TypeError, "Polygon path must be a list of indices");
        return NULL;
      }
    }
  }

  node->convexity = convexity;
  if (node->convexity < 1) node->convexity = 1;

  python_retrieve_pyname(node);
  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

PyObject *python_spline(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE
  unsigned int i;
  auto node = std::make_shared<SplineNode>(instance);

  char *kwlist[] = {"points", "fn", "fa", "fs", NULL};
  PyObject *points = NULL;
  double fn = 0, fa = 0, fs = 0;

  PyObject *element;
  Vector2d point;

  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O!|ddd", kwlist, pf.PyList_Type, &points, &fn, &fa,
                                      &fs)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing spline(points)");
    return NULL;
  }

  if (points != NULL && PyList_CHECK(points)) {
    if (pf.PyList_Size(points) == 0) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "There must at least be one point in the polygon");
      return NULL;
    }
    for (i = 0; i < pf.PyList_Size(points); i++) {
      element = pf.PyList_GetItem(points, i);
      if (PyList_CHECK(element) && pf.PyList_Size(element) == 2) {
        point[0] = pf.PyFloat_AsDouble(pf.PyList_GetItem(element, 0));
        point[1] = pf.PyFloat_AsDouble(pf.PyList_GetItem(element, 1));
        node->points.push_back(point);
      } else {
        pf.PyErr_SetString(pf.PyExc_TypeError, "Coordinate must exactly contain 2 numbers");
        return NULL;
      }
    }
  } else {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Polygon points must be a list of coordinates");
    return NULL;
  }
  node->fn = fn;
  node->fa = fa;
  node->fs = fs;

  python_retrieve_pyname(node);
  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

int python_tomatrix(PyObject *pyt, Matrix4d& mat)
{
  if (pyt == nullptr) return 1;
  PyObject *row, *cell;
  double val;
  if (!PyList_CHECK(pyt)) return 1;
  if (pf.PyList_Size(pyt) != 4) return 1;
  for (int i = 0; i < 4; i++) {
    row = pf.PyList_GetItem(pyt, i);
    if (!PyList_CHECK(row)) return 1;
    if (pf.PyList_Size(row) != 4) return 1;
    for (int j = 0; j < 4; j++) {
      cell = pf.PyList_GetItem(row, j);
      if (python_numberval(cell, &val)) return 1;
      mat(i, j) = val;
    }
  }
  return 0;
}

int python_tovector(PyObject *pyt, Vector3d& vec)
{
  if (pyt == nullptr) return 1;
  PyObject *cell;
  double val;
  if (!PyList_CHECK(pyt)) return 1;
  if (pf.PyList_Size(pyt) != 3) return 1;
  for (int i = 0; i < 3; i++) {
    cell = pf.PyList_GetItem(pyt, i);
    if (python_numberval(cell, &val)) return 1;
    vec[i] = val;
  }
  return 0;
}

PyObject *python_frommatrix(const Matrix4d& mat)
{
  PyObject *pyo = pf.PyList_New(4);
  PyObject *row;
  for (int i = 0; i < 4; i++) {
    row = pf.PyList_New(4);
    for (int j = 0; j < 4; j++) pf.PyList_SetItem(row, j, pf.PyFloat_FromDouble(mat(i, j)));
    pf.PyList_SetItem(pyo, i, row);
    //      Py_XDECREF(row);
  }
  return pyo;
}

PyObject *python_fromvector(const Vector3d vec)
{
  PyObject *res = pf.PyList_New(3);
  for (int i = 0; i < 3; i++) pf.PyList_SetItem(res, i, pf.PyFloat_FromDouble(vec[i]));
  return res;
}

PyObject *python_number_scale(PyObject *pynum, Vector3d scalevec)
{
  Matrix4d mat;
  if (!python_tomatrix(pynum, mat)) {
    Transform3d matrix = Transform3d::Identity();
    matrix.scale(scalevec);
    Vector3d n;
    for (int i = 0; i < 3; i++) {
      n = Vector3d(mat(0, i), mat(1, i), mat(2, i));
      n = matrix * n;
      for (int j = 0; j < 3; j++) mat(j, i) = n[j];
    }
    return python_frommatrix(mat);
  }
  Vector3d vec;
  if (!python_tovector(pynum, vec)) {
    for (int i = 0; i < 3; i++) vec[i] *= scalevec[i];
    return python_fromvector(vec);
  }
  return nullptr;
}

PyObject *python_scale_sub(PyObject *obj, Vector3d scalevec)
{
  PyObject *mat = python_number_scale(obj, scalevec);
  if (mat != nullptr) return mat;

  DECLARE_INSTANCE
  std::shared_ptr<AbstractNode> child;
  auto node = std::make_shared<TransformNode>(instance, "scale");
  PyObject *child_dict;
  child = PyOpenSCADObjectToNodeMulti(obj, &child_dict);
  if (child == NULL) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid type for Object in scale");
    return NULL;
  }
  node->matrix.scale(scalevec);
  node->setPyName(child->getPyName());
  node->children.push_back(child);
  PyObject *pyresult = PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
  if (child_dict != nullptr) {
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (pf.PyDict_Next(child_dict, &pos, &key, &value)) {
      PyObject *value1 = python_number_scale(value, scalevec, 4);
      if (value1 != nullptr) pf.PyDict_SetItem(((PyOpenSCADObject *)pyresult)->dict, key, value1);
      else pf.PyDict_SetItem(((PyOpenSCADObject *)pyresult)->dict, key, value);
    }
  }
  return pyresult;
}

PyObject *python_scale_core(PyObject *obj, PyObject *val_v)
{
  double x = 1, y = 1, z = 1;
  if (python_vectorval(val_v, 2, 3, &x, &y, &z)) {
    pf.PyErr_SetString(pf.PyExc_TypeError,
                       "Invalid vector specifiaction in scale, use 1 to 3 ordinates.");
    return NULL;
  }
  Vector3d scalevec(x, y, z);

  if (OpenSCAD::rangeCheck) {
    if (scalevec[0] == 0 || scalevec[1] == 0 || scalevec[2] == 0 || !std::isfinite(scalevec[0]) ||
        !std::isfinite(scalevec[1]) || !std::isfinite(scalevec[2])) {
      //      LOG(message_group::Warning, instance->location(), parameters.documentRoot(), "scale(%1$s)",
      //      parameters["v"].toEchoStringNoThrow());
    }
  }

  return python_scale_sub(obj, scalevec);
}

PyObject *python_scale(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"obj", "v", NULL};
  PyObject *obj = NULL;
  PyObject *val_v = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "OO", kwlist, &obj, &val_v)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing scale(object, scale)");
    return NULL;
  }
  return python_scale_core(obj, val_v);
}

PyObject *python_oo_scale(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"v", NULL};
  PyObject *val_v = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &val_v)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing scale(object, scale)");
    return NULL;
  }
  return python_scale_core(obj, val_v);
}

PyObject *python_explode(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"obj", "v", NULL};
  PyObject *obj = NULL;
  PyObject *val_v = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "OO", kwlist, &obj, &val_v)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing explode(object, list)");
    return NULL;
  }
  return python_nb_sub_vec3(obj, val_v, 3);
}

PyObject *python_oo_explode(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"v", NULL};
  PyObject *val_v = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &val_v)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing explode(object, list)");
    return NULL;
  }
  return python_nb_sub_vec3(obj, val_v, 3);
}

PyObject *python_number_rot(PyObject *mat, Matrix3d rotvec, int vecs)
{
  Transform3d matrix = Transform3d::Identity();
  matrix.rotate(rotvec);
  Matrix4d raw;
  if (python_tomatrix(mat, raw)) return nullptr;
  Vector3d n;
  for (int i = 0; i < vecs; i++) {
    n = Vector3d(raw(0, i), raw(1, i), raw(2, i));
    n = matrix * n;
    for (int j = 0; j < 3; j++) raw(j, i) = n[j];
  }
  return python_frommatrix(raw);
}

PyObject *python_rotate_sub(PyObject *obj, Vector3d vec3, double angle, PyObject *ref, int dragflags)
{
  Matrix3d M;
  if (isnan(angle)) {
    double sx = 0, sy = 0, sz = 0;
    double cx = 1, cy = 1, cz = 1;
    double a = 0.0;
    if (vec3[2] != 0) {
      a = vec3[2];
      sz = sin_degrees(a);
      cz = cos_degrees(a);
    }
    if (vec3[1] != 0) {
      a = vec3[1];
      sy = sin_degrees(a);
      cy = cos_degrees(a);
    }
    if (vec3[0] != 0) {
      a = vec3[0];
      sx = sin_degrees(a);
      cx = cos_degrees(a);
    }

    M << cy * cz, cz * sx * sy - cx * sz, cx * cz * sy + sx * sz, cy * sz, cx * cz + sx * sy * sz,
      -cz * sx + cx * sy * sz, -sy, cy * sx, cx * cy;
  } else {
    M = angle_axis_degrees(angle, vec3);
  }
  PyObject *mat = python_number_rot(obj, M, 4);
  if (mat != nullptr) return mat;

  DECLARE_INSTANCE
  auto node = std::make_shared<TransformNode>(instance, "rotate");
  node->dragflags = dragflags;

  PyObject *child_dict;
  std::shared_ptr<AbstractNode> child = PyOpenSCADObjectToNodeMulti(obj, &child_dict);
  PyTypeObject *type = PyOpenSCADObjectType(obj);
  if (child == NULL) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid type for Object in rotate");
    return NULL;
  }
  node->matrix.rotate(M);
  node->setPyName(child->getPyName());

  PyObject *pyresult;
  if (ref == nullptr) {
    node->children.push_back(child);
    pyresult = PyOpenSCADObjectFromNode(type, node);
  } else {
    Vector3d vec3;
    python_vectorval(ref, 1, 3, &(vec3[0]), &(vec3[1]), &(vec3[2]), nullptr, &dragflags);

    std::shared_ptr<TransformNode> prenode, postnode;
    {
      DECLARE_INSTANCE
      prenode = std::make_shared<TransformNode>(instance, "translate");
      prenode->matrix.translate(-vec3);
    }
    {
      DECLARE_INSTANCE
      postnode = std::make_shared<TransformNode>(instance, "translate");
      postnode->matrix.translate(vec3);
    }
    prenode->children.push_back(child);
    node->children.push_back(prenode);
    postnode->children.push_back(node);
    pyresult = PyOpenSCADObjectFromNode(type, postnode);
  }
  if (child_dict != nullptr) {
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (pf.PyDict_Next(child_dict, &pos, &key, &value)) {
      PyObject *value1 = python_number_rot(value, M, 4);
      if (value1 != nullptr) pf.PyDict_SetItem(((PyOpenSCADObject *)pyresult)->dict, key, value1);
      else pf.PyDict_SetItem(((PyOpenSCADObject *)pyresult)->dict, key, value);
    }
  }
  return pyresult;
}

PyObject *python_rotate_core(PyObject *obj, PyObject *val_a, PyObject *val_v, PyObject *ref)
{
  Vector3d vec3(0, 0, 0);
  double angle;
  int dragflags = 0;
  if (val_a != nullptr && PyList_CHECK(val_a) && val_v == nullptr) {
    python_vectorval(val_a, 1, 3, &(vec3[0]), &(vec3[1]), &(vec3[2]), nullptr, &dragflags);
    return python_rotate_sub(obj, vec3, NAN, dragflags);
  } else if (val_a != nullptr && val_v != nullptr && !python_numberval(val_a, &angle) &&
             PyList_CHECK(val_v) && pf.PyList_Size(val_v) == 3) {
    vec3[0] = pf.PyFloat_AsDouble(pf.PyList_GetItem(val_v, 0));
    vec3[1] = pf.PyFloat_AsDouble(pf.PyList_GetItem(val_v, 1));
    vec3[2] = pf.PyFloat_AsDouble(pf.PyList_GetItem(val_v, 2));
    return python_rotate_sub(obj, vec3, angle, dragflags);
  }
  pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid arguments to rotate()");
  return nullptr;
}

PyObject *python_rotate(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"obj", "a", "v", "ref", nullptr};
  PyObject *val_a = nullptr;
  PyObject *val_v = nullptr;
  PyObject *obj = nullptr;
  PyObject *ref = nullptr;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "OO|O", kwlist, &obj, &val_a, &val_v), &ref) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing rotate(object, vec3)");
    return NULL;
  }
  return python_rotate_core(obj, val_a, val_v, ref);
}

PyObject *python_oo_rotate(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"a", "v", "ref", nullptr};
  PyObject *val_a = nullptr;
  PyObject *val_v = nullptr;
  PyObject *ref = nullptr;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O|O", kwlist, &val_a, &val_v)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing rotate(object, vec3)");
    return NULL;
  }
  return python_rotate_core(obj, val_a, val_v, ref);
}

PyObject *python_matrix_mirror(PyObject *mat, Matrix4d m)
{
  Matrix4d raw;
  if (python_tomatrix(mat, raw)) return nullptr;
  Vector4d n;
  for (int i = 0; i < 4; i++) {
    n = Vector4d(raw(0, i), raw(1, i), raw(2, i), 0);
    n = m * n;
    for (int j = 0; j < 3; j++) raw(j, i) = n[j];
  }
  return python_frommatrix(raw);
}

PyObject *python_mirror_sub(PyObject *obj, Matrix4d& m)
{
  PyObject *mat = python_matrix_mirror(obj, m);
  if (mat != nullptr) return mat;

  DECLARE_INSTANCE
  auto node = std::make_shared<TransformNode>(instance, "mirror");
  node->matrix = m;
  PyObject *child_dict;
  PyTypeObject *type = PyOpenSCADObjectType(obj);
  std::shared_ptr<AbstractNode> child = PyOpenSCADObjectToNodeMulti(obj, &child_dict);
  if (child == NULL) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid type for Object in mirror");
    return NULL;
  }
  node->children.push_back(child);
  node->setPyName(child->getPyName());
  PyObject *pyresult = PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
  if (child_dict != nullptr) {
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (pf.PyDict_Next(child_dict, &pos, &key, &value)) {
      PyObject *value1 = python_number_mirror(value, m, 4);
      if (value1 != nullptr) pf.PyDict_SetItem(((PyOpenSCADObject *)pyresult)->dict, key, value1);
      else pf.PyDict_SetItem(((PyOpenSCADObject *)pyresult)->dict, key, value);
    }
  }
  return pyresult;
}

PyObject *python_mirror_core(PyObject *obj, PyObject *val_v)
{
  Vector3d mirrorvec;
  double x = 1.0, y = 0.0, z = 0.0;
  if (python_vectorval(val_v, 2, 3, &x, &y, &z)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid vector specifiaction in mirror");
    return NULL;
  }
  // x /= sqrt(x*x + y*y + z*z)
  // y /= sqrt(x*x + y*y + z*z)
  // z /= sqrt(x*x + y*y + z*z)
  Matrix4d m = Matrix4d::Identity();
  if (x != 0.0 || y != 0.0 || z != 0.0) {
    // skip using sqrt to normalize the vector since each element of matrix contributes it with two
    // multiplied terms instead just divide directly within each matrix element simplified calculation
    // leads to less float errors
    double a = x * x + y * y + z * z;

    m << 1 - 2 * x * x / a, -2 * y * x / a, -2 * z * x / a, 0, -2 * x * y / a, 1 - 2 * y * y / a,
      -2 * z * y / a, 0, -2 * x * z / a, -2 * y * z / a, 1 - 2 * z * z / a, 0, 0, 0, 0, 1;
  }
  return python_mirror_sub(obj, m);
}

PyObject *python_mirror(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"obj", "v", NULL};

  PyObject *obj = NULL;
  PyObject *val_v = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "OO", kwlist, &obj, &val_v)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing mirror(object, vec3)");
    return NULL;
  }
  return python_mirror_core(obj, val_v);
}

PyObject *python_oo_mirror(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"v", NULL};

  PyObject *val_v = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &val_v)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing mirror(object, vec3)");
    return NULL;
  }
  return python_mirror_core(obj, val_v);
}

PyObject *python_number_trans(PyObject *pynum, Vector3d transvec, int vecs)
{
  Matrix4d mat;
  if (!python_tomatrix(pynum, mat)) {
    for (int i = 0; i < 3; i++) mat(i, 3) += transvec[i];
    return python_frommatrix(mat);
  }

  Vector3d vec;
  if (!python_tovector(pynum, vec)) {
    return python_fromvector(vec + transvec);
  }

  return nullptr;
}

PyObject *python_translate_sub(PyObject *obj, Vector3d translatevec, int dragflags)
{
  PyObject *child_dict;
  PyObject *mat = python_number_trans(obj, translatevec, 4);
  if (mat != nullptr) return mat;

  DECLARE_INSTANCE
  auto node = std::make_shared<TransformNode>(instance, "translate");
  std::shared_ptr<AbstractNode> child;
  PyTypeObject *type = PyOpenSCADObjectType(obj);
  child = PyOpenSCADObjectToNodeMulti(obj, &child_dict);
  node->setPyName(child->getPyName());
  node->dragflags = dragflags;
  if (child == NULL) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid type for Object in translate");
    return NULL;
  }
  node->matrix.translate(translatevec);

  node->children.push_back(child);
  PyObject *pyresult = PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
  if (child_dict != nullptr) {  // TODO dies ueberall
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (pf.PyDict_Next(child_dict, &pos, &key, &value)) {
      PyObject *value1 = python_number_trans(value, translatevec, 4);
      if (value1 != nullptr) pf.PyDict_SetItem(((PyOpenSCADObject *)pyresult)->dict, key, value1);
      else pf.PyDict_SetItem(((PyOpenSCADObject *)pyresult)->dict, key, value);
    }
  }
  return pyresult;
}

PyObject *python_nb_sub_vec3(PyObject *arg1, PyObject *arg2, int mode);
PyObject *python_translate_core(PyObject *obj, PyObject *v) { return python_nb_sub_vec3(obj, v, 0); }

PyObject *python_translate(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"obj", "v", NULL};
  PyObject *obj = NULL;
  PyObject *v = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O|O", kwlist, &obj, &v)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing translate(object,vec3)");
    return NULL;
  }
  return python_translate_core(obj, v);
}

PyObject *python_oo_translate(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"v", NULL};
  PyObject *v = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "|O", kwlist, &v)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing translate(object,vec3)");
    return NULL;
  }
  return python_translate_core(obj, v);
}

PyObject *python_dir_sub_core(PyObject *obj, double arg, int mode)
{
  if (mode < 6) {
    Vector3d trans;
    switch (mode) {
    case 0: trans = Vector3d(arg, 0, 0); break;
    case 1: trans = Vector3d(-arg, 0, 0); break;
    case 2: trans = Vector3d(0, -arg, 0); break;
    case 3: trans = Vector3d(0, arg, 0); break;
    case 4: trans = Vector3d(0, 0, -arg); break;
    case 5: trans = Vector3d(0, 0, arg); break;
    }
    return python_translate_sub(obj, trans, 0);
  } else {
    Vector3d rot;
    switch (mode) {
    case 6: rot = Vector3d(arg, 0, 0); break;
    case 7: rot = Vector3d(0, arg, 0); break;
    case 8: rot = Vector3d(0, 0, arg); break;
    }
    return python_rotate_sub(obj, rot, NAN, nullptr, 0);
  }
}

PyObject *python_dir_sub(PyObject *self, PyObject *args, PyObject *kwargs, int mode)
{
  char *kwlist[] = {"obj", "v", NULL};
  PyObject *obj = NULL;
  double arg;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "Od", kwlist, &obj, &arg)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing translate(object,vec3)");
    return NULL;
  }
  return python_dir_sub_core(obj, arg, mode);
}

PyObject *python_oo_dir_sub(PyObject *obj, PyObject *args, PyObject *kwargs, int mode)
{
  char *kwlist[] = {"v", NULL};
  double arg;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "d", kwlist, &arg)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing translate(object,vec3)");
    return NULL;
  }
  return python_dir_sub_core(obj, arg, mode);
}

PyObject *python_right(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_dir_sub(self, args, kwargs, 0);
}
PyObject *python_oo_right(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_oo_dir_sub(self, args, kwargs, 0);
}
PyObject *python_left(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_dir_sub(self, args, kwargs, 1);
}
PyObject *python_oo_left(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_oo_dir_sub(self, args, kwargs, 1);
}
PyObject *python_front(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_dir_sub(self, args, kwargs, 2);
}
PyObject *python_oo_front(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_oo_dir_sub(self, args, kwargs, 2);
}
PyObject *python_back(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_dir_sub(self, args, kwargs, 3);
}
PyObject *python_oo_back(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_oo_dir_sub(self, args, kwargs, 3);
}
PyObject *python_down(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_dir_sub(self, args, kwargs, 4);
}
PyObject *python_oo_down(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_oo_dir_sub(self, args, kwargs, 4);
}
PyObject *python_up(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_dir_sub(self, args, kwargs, 5);
}
PyObject *python_oo_up(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_oo_dir_sub(self, args, kwargs, 5);
}
PyObject *python_rotx(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_dir_sub(self, args, kwargs, 6);
}
PyObject *python_oo_rotx(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_oo_dir_sub(self, args, kwargs, 6);
}
PyObject *python_roty(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_dir_sub(self, args, kwargs, 7);
}
PyObject *python_oo_roty(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_oo_dir_sub(self, args, kwargs, 7);
}
PyObject *python_rotz(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_dir_sub(self, args, kwargs, 8);
}
PyObject *python_oo_rotz(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_oo_dir_sub(self, args, kwargs, 8);
}

PyObject *python_math_sub1(PyObject *self, PyObject *args, PyObject *kwargs, int mode)
{
  char *kwlist[] = {"value", NULL};
  double arg;
  double result = 0;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "d", kwlist, &arg)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing math function");
    return NULL;
  }
  switch (mode) {
  case 0: result = sin(arg * G_PI / 180.0); break;
  case 1: result = cos(arg * G_PI / 180.0); break;
  case 2: result = tan(arg * G_PI / 180.0); break;
  case 3: result = asin(arg) * 180.0 / G_PI; break;
  case 4: result = acos(arg) * 180.0 / G_PI; break;
  case 5: result = atan(arg) * 180.0 / G_PI; break;
  }
  return pf.PyFloat_FromDouble(result);
}

PyObject *python_math_sub2(PyObject *self, PyObject *args, PyObject *kwargs, int mode)
{
  int dragflags = 0;
  char *kwlist[] = {"vec1", "vec2", NULL};
  PyObject *obj1 = nullptr;
  PyObject *obj2 = nullptr;
  Vector3d vec31(0, 0, 0);
  Vector3d vec32(0, 0, 0);
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "OO", kwlist, &obj1, &obj2)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing norm(vec3)");
    return NULL;
  }
  python_vectorval(obj1, 1, 3, &(vec31[0]), &(vec31[1]), &(vec31[2]), nullptr, &dragflags);
  python_vectorval(obj2, 1, 3, &(vec32[0]), &(vec32[1]), &(vec32[2]), nullptr, &dragflags);

  switch (mode) {
  case 0: return pf.PyFloat_FromDouble(vec31.dot(vec32)); break;
  case 1:
    Vector3d res = vec31.cross(vec32);
    return python_fromvector(vec31.cross(vec32));
    break;
  }
  return Py_NONE;
}

PyObject *python_sin(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_math_sub1(self, args, kwargs, 0);
}
PyObject *python_cos(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_math_sub1(self, args, kwargs, 1);
}
PyObject *python_tan(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_math_sub1(self, args, kwargs, 2);
}
PyObject *python_asin(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_math_sub1(self, args, kwargs, 3);
}
PyObject *python_acos(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_math_sub1(self, args, kwargs, 4);
}
PyObject *python_atan(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_math_sub1(self, args, kwargs, 5);
}

PyObject *python_dot(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_math_sub2(self, args, kwargs, 0);
}
PyObject *python_cross(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_math_sub2(self, args, kwargs, 1);
}

PyObject *python_norm(PyObject *self, PyObject *args, PyObject *kwargs)
{
  int dragflags = 0;
  char *kwlist[] = {"vec", NULL};
  double result = 0;
  PyObject *obj = nullptr;
  Vector3d vec3(0, 0, 0);
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &obj)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing norm(vec3)");
    return NULL;
  }
  python_vectorval(obj, 1, 3, &(vec3[0]), &(vec3[1]), &(vec3[2]), nullptr, &dragflags);

  result = sqrt(vec3[0] * vec3[0] + vec3[1] * vec3[1] + vec3[2] * vec3[2]);
  return pf.PyFloat_FromDouble(result);
}

PyObject *python_multmatrix_sub(PyObject *pyobj, PyObject *pymat, int div)
{
  Matrix4d mat;
  if (!python_tomatrix(pymat, mat)) {
    double w = mat(3, 3);
    if (w != 1.0) mat = mat / w;
  } else {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Matrix vector should be 4x4 array");
    return NULL;
  }
  if (div) {
    auto tmp = mat.inverse().eval();
    mat = tmp;
  }

  Matrix4d objmat;
  if (!python_tomatrix(pyobj, objmat)) {
    objmat = mat * objmat;
    return python_frommatrix(objmat);
  }

  DECLARE_INSTANCE
  auto node = std::make_shared<TransformNode>(instance, "multmatrix");
  std::shared_ptr<AbstractNode> child;
  PyObject *child_dict;
  PyTypeObject *type = PyOpenSCADObjectType(pyobj);
  child = PyOpenSCADObjectToNodeMulti(pyobj, &child_dict);
  node->setPyName(child->getPyName());
  if (child == NULL) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid type for Object in multmatrix");
    return NULL;
  }

  node->matrix = mat;
  node->children.push_back(child);
  PyObject *pyresult = PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
  if (child_dict != nullptr) {
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (pf.PyDict_Next(child_dict, &pos, &key, &value)) {
      Matrix4d raw;
      if (python_tomatrix(value, raw)) return nullptr;
      PyObject *value1 = python_frommatrix(node->matrix * raw);
      if (value1 != nullptr) pf.PyDict_SetItem(((PyOpenSCADObject *)pyresult)->dict, key, value1);
      else pf.PyDict_SetItem(((PyOpenSCADObject *)pyresult)->dict, key, value);
    }
  }
  return pyresult;
}

PyObject *python_multmatrix(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"obj", "m", NULL};
  PyObject *obj = NULL;
  PyObject *mat = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "OO!", kwlist, &obj, pf.PyList_Type, &mat)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing multmatrix(object, vec16)");
    return NULL;
  }
  return python_multmatrix_sub(obj, mat, 0);
}

PyObject *python_oo_multmatrix(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"m", NULL};
  PyObject *mat = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O!", kwlist, pf.PyList_Type, &mat)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing multmatrix(object, vec16)");
    return NULL;
  }
  return python_multmatrix_sub(obj, mat, 0);
}

PyObject *python_divmatrix(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"obj", "m", NULL};
  PyObject *obj = NULL;
  PyObject *mat = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "OO!", kwlist, &obj, pf.PyList_Type, &mat)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing divmatrix(object, vec16)");
    return NULL;
  }
  return python_multmatrix_sub(obj, mat, 1);
}

PyObject *python_oo_divmatrix(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"m", NULL};
  PyObject *mat = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O!", kwlist, pf.PyList_Type, &mat)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing divmatrix(object, vec16)");
    return NULL;
  }
  return python_multmatrix_sub(obj, mat, 1);
}

PyObject *python_pull_core(PyObject *obj, PyObject *anchor, PyObject *dir)
{
  DECLARE_INSTANCE
  auto node = std::make_shared<PullNode>(instance);
  PyObject *dummydict;
  PyTypeObject *type = PyOpenSCADObjectType(obj);
  std::shared_ptr<AbstractNode> child = PyOpenSCADObjectToNodeMulti(obj, &dummydict);
  if (child == NULL) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid type for  Object in translate\n");
    return NULL;
  }

  double x = 0, y = 0, z = 0;
  if (python_vectorval(anchor, 3, 3, &x, &y, &z)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid vector specifiaction in anchor\n");
    return NULL;
  }
  node->anchor = Vector3d(x, y, z);

  if (python_vectorval(dir, 3, 3, &x, &y, &z)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid vector specifiaction in dir\n");
    return NULL;
  }
  node->dir = Vector3d(x, y, z);

  node->children.push_back(child);
  return PyOpenSCADObjectFromNode(type, node);
}

PyObject *python_pull(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"obj", "src", "dst", NULL};
  PyObject *obj = NULL;
  PyObject *anchor = NULL;
  PyObject *dir = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "OOO|", kwlist, &obj, &anchor, &dir)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error during parsing\n");
    return NULL;
  }
  return python_pull_core(obj, anchor, dir);
}

PyObject *python_oo_pull(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"src", "dst", NULL};
  PyObject *anchor = NULL;
  PyObject *dir = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "OO|", kwlist, &anchor, &dir)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error during parsing\n");
    return NULL;
  }
  return python_pull_core(obj, anchor, dir);
}

PyObject *python_wrap_core(PyObject *obj, PyObject *target, double r, double d, double fn, double fa,
                           double fs)
{
  DECLARE_INSTANCE
  auto node = std::make_shared<WrapNode>(instance);

  PyObject *dummydict;
  PyTypeObject *type = PyOpenSCADObjectType(obj);
  std::shared_ptr<AbstractNode> child = PyOpenSCADObjectToNodeMulti(obj, &dummydict);
  if (child == NULL) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid type for  Object in Wrap\n");
    return NULL;
  }

  if (PyFloat_CHECK(target)) {
    printf("is float\n");
    node->r = pf.PyFloat_AsDouble(target);
    node->shape = nullptr;
  } else if (Py_TYPE(target) == &PyOpenSCADType) {
    std::shared_ptr<AbstractNode> abstr = ((PyOpenSCADObject *)target)->node;
    node->shape = abstr;
  } else {
    pf.PyErr_SetString(pf.PyExc_TypeError,
                       "warpign object must bei either Polygon or cylidner radius\n");
    return NULL;
  }

  get_fnas(node->fn, node->fa, node->fs);
  if (!isnan(fn)) node->fn = fn;
  if (!isnan(fa)) node->fa = fa;
  if (!isnan(fs)) node->fs = fs;
  node->children.push_back(child);
  return PyOpenSCADObjectFromNode(type, node);
}

PyObject *python_wrap(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"obj", "target", "fn", "fa", "fs", NULL};
  PyObject *obj = NULL, *target = NULL;
  double fn, fa, fs;
  double r = NAN, d = NAN;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "OO|ddd", kwlist, &obj, &target, &fn, &fa, &fs)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error during parsing wrap\n");
    return NULL;
  }
  return python_wrap_core(obj, target, r, d, fn, fa, fs);
}

PyObject *python_oo_wrap(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"target", "fn", "fa", "fs", NULL};
  double fn = NAN, fa = NAN, fs = NAN;
  PyObject *target = NULL;
  double r = NAN, d = NAN;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O|ddd", kwlist, &target, &fn, &fa, &fs)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error during parsing\n");
    return NULL;
  }
  std::vector<double> xsteps;
  xsteps.push_back(4);
  xsteps.push_back(6);

  return python_wrap_core(obj, target, r, d, fn, fa, fs);
}

PyObject *python_show_core(PyObject *obj)
{
  python_result_obj = obj;
  PyObject *child_dict = nullptr;
  std::shared_ptr<AbstractNode> child = PyOpenSCADObjectToNodeMulti(obj, &child_dict);
  if (child == NULL) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid type for Object in show");
    return NULL;
  }
  if (child == void_node) {
    return nullptr;
  }

  if (child == full_node) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Cannot display infinite space");
    return nullptr;
  }

  PyObject *key, *value;
  Py_ssize_t pos = 0;
  python_build_hashmap(child, 0);
  std::string varname = child->getPyName();
  if (child_dict != nullptr) {
    while (pf.PyDict_Next(child_dict, &pos, &key, &value)) {
      Matrix4d raw;
      if (python_tomatrix(value, raw)) continue;
      PyObject *value1 = pf.PyUnicode_AsEncodedString(key, "utf-8", "~");
      const char *value_str = PyBytes_AS_STRING(value1);
      SelectedObject sel;
      sel.pt.clear();
      sel.pt.push_back(Vector3d(raw(0, 3), raw(1, 3), raw(2, 3)));
      sel.pt.push_back(Vector3d(raw(0, 0), raw(1, 0), raw(2, 0)));
      sel.pt.push_back(Vector3d(raw(0, 1), raw(1, 1), raw(2, 1)));
      sel.pt.push_back(Vector3d(raw(0, 2), raw(1, 2), raw(2, 2)));
      sel.type = SelectionType::SELECTION_HANDLE;
      sel.name = varname + "." + value_str;
      python_result_handle.push_back(sel);
    }
  }
  shows.push_back(child);
  Py_INCREF(obj);
  return obj;
}

PyObject *python_show(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyObject *obj = NULL;
  char *kwlist[] = {"obj", NULL};
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &obj)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing output(object)");
    return NULL;
  }
  return python_show_core(obj);
}

PyObject *python_oo_show(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {NULL};
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "", kwlist)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing output(object)");
    return NULL;
  }
  return python_show_core(obj);
}

PyObject *python_output(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  LOG(message_group::Deprecated, "output is deprecated, please use show() instead");
  return python_show(obj, args, kwargs);
}

PyObject *python_oo_output(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  LOG(message_group::Deprecated, "output is deprecated, please use show() instead");
  return python_oo_show(obj, args, kwargs);
}

void Export3mfPartInfo::writeProps(void *obj) const
{
  if (this->props == nullptr) return;
  PyObject *prop = (PyObject *)this->props;
  if (!PyDict_Check(prop)) return;
  PyObject *key, *value;
  Py_ssize_t pos = 0;
  while (pf.PyDict_Next(prop, &pos, &key, &value)) {
    PyObject *key1 = pf.PyUnicode_AsEncodedString(key, "utf-8", "~");
    const char *key_str = PyBytes_AS_STRING(key1);
    if (key_str == nullptr) continue;
    if (PyFloat_CHECK(value)) {
      writePropsFloat(obj, key_str, pf.PyFloat_AsDouble(value));
    }
    if (PyLong_Check(value)) {
      writePropsLong(obj, key_str, pf.PyLong_AsLong(value));
    }
    if (PyUnicode_Check(value)) {
      PyObject *val1 = pf.PyUnicode_AsEncodedString(value, "utf-8", "~");
      const char *val_str = PyBytes_AS_STRING(val1);
      writePropsString(obj, key_str, val_str);
    }
  }
}

void python_export_obj_att(std::ostream& output)
{
  PyObject *child_dict = nullptr;
  if (python_result_obj == nullptr) return;
  std::shared_ptr<AbstractNode> child = PyOpenSCADObjectToNodeMulti(python_result_obj, &child_dict);
  if (child_dict == nullptr) return;
  if (!PyDict_Check(child_dict)) return;
  PyObject *key, *value;
  Py_ssize_t pos = 0;
  while (pf.PyDict_Next(child_dict, &pos, &key, &value)) {
    PyObject *key1 = pf.PyUnicode_AsEncodedString(key, "utf-8", "~");
    const char *key_str = PyBytes_AS_STRING(key1);
    if (key_str == nullptr) continue;

    if (PyLong_Check(value)) output << "# " << key_str << " = " << pf.PyLong_AsLong(value) << "\n";

    if (PyFloat_CHECK(value)) output << "# " << key_str << " = " << pf.PyFloat_AsDouble(value) << "\n";

    if (PyUnicode_Check(value)) {
      auto valuestr = std::string(pf.PyUnicode_AsUTF8(value));
      output << "# " << key_str << " = \"" << valuestr << "\"\n";
    }
  }
}

PyObject *python_export_core(PyObject *obj, char *file)
{
  std::string filename;
  if (python_scriptpath.string().size() > 0)
    filename = lookup_file(file, python_scriptpath.parent_path().u8string(), ".");  // TODO problem hbier
  else filename = file;
  const auto path = fs::path(filename);
  std::string suffix = path.has_extension() ? path.extension().generic_string().substr(1) : "";
  boost::algorithm::to_lower(suffix);
  python_result_obj = obj;

  FileFormat exportFileFormat = FileFormat::BINARY_STL;
  if (!fileformat::fromIdentifier(suffix, exportFileFormat)) {
    LOG("Invalid suffix %1$s. Defaulting to binary STL.", suffix);
  }

  std::vector<Export3mfPartInfo> export3mfPartInfos;
  std::vector<std::string> names;

  PyObject *child_dict;
  std::shared_ptr<AbstractNode> child = PyOpenSCADObjectToNodeMulti(obj, &child_dict);
  if (child != nullptr) {
    Tree tree(child, "parent");
    GeometryEvaluator geomevaluator(tree);
    Export3mfPartInfo info(geomevaluator.evaluateGeometry(*tree.root(), false), "OpenSCAD Model",
                           nullptr);
    export3mfPartInfos.push_back(info);
  } else if (PyDict_Check(obj)) {
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (pf.PyDict_Next(obj, &pos, &key, &value)) {
      PyObject *value1 = pf.PyUnicode_AsEncodedString(key, "utf-8", "~");
      const char *value_str = PyBytes_AS_STRING(value1);
      if (value_str == nullptr) continue;
      std::shared_ptr<AbstractNode> child = PyOpenSCADObjectToNodeMulti(value, &child_dict);
      if (child == nullptr) continue;

      void *prop = nullptr;
      if (child_dict != nullptr && PyDict_Check(child_dict)) {
        PyObject *key = pf.PyUnicode_FromStringAndSize("props_3mf", 9);
        prop = pf.PyDict_GetItem(child_dict, key);
      }
      Tree tree(child, "parent");
      GeometryEvaluator geomevaluator(tree);
      Export3mfPartInfo info(geomevaluator.evaluateGeometry(*tree.root(), false), value_str, prop);
      export3mfPartInfos.push_back(info);
    }
  }
  if (export3mfPartInfos.size() == 0) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Object not recognized");
    return NULL;
  }

  Export3mfOptions options3mf;
  options3mf.decimalPrecision = 6;
  options3mf.color = "#f9d72c";
  ExportInfo exportInfo = {.format = exportFileFormat,
                           .sourceFilePath = file,
                           .options3mf = std::make_shared<Export3mfOptions>(options3mf)};

  if (exportFileFormat == FileFormat::_3MF) {
    std::ofstream fstream(file, std::ios::out | std::ios::trunc | std::ios::binary);
    if (!fstream.is_open()) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "Can't write export file");
      return nullptr;
    }
    export_3mf(export3mfPartInfos, fstream, exportInfo);
  } else {
    if (export3mfPartInfos.size() > 1) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "This Format can at most export one object");
      return nullptr;
    }
    exportFileByName(export3mfPartInfos[0].geom, file, exportInfo);
  }
  return Py_NONE;
}

PyObject *python_export(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyObject *obj = NULL;
  char *file = nullptr;
  char *kwlist[] = {"obj", "file", NULL};
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "Os|O", kwlist, &obj, &file)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing output(object)");
    return NULL;
  }
  return python_export_core(obj, file);
}

PyObject *python_oo_export(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"file", NULL};
  char *file = nullptr;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "s", kwlist, &file)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing output(object)");
    return NULL;
  }
  return python_export_core(obj, file);
}

PyObject *python_find_face_core(PyObject *obj, PyObject *vec_p)
{
  Vector3d vec;
  double dummy;
  PyObject *child_dict;
  if (python_vectorval(vec_p, 3, 3, &vec[0], &vec[1], &vec[2], &dummy)) return Py_NONE;
  std::shared_ptr<AbstractNode> child = PyOpenSCADObjectToNodeMulti(obj, &child_dict);

  if (child == nullptr) return Py_NONE;

  Tree tree(child, "");
  GeometryEvaluator geomevaluator(tree);
  std::shared_ptr<const Geometry> geom = geomevaluator.evaluateGeometry(*tree.root(), true);
  std::shared_ptr<const PolySet> ps = PolySetUtils::getGeometryAsPolySet(geom);

  double dmax = -1;
  Vector4d vectormax;
  for (auto ind : ps->indices) {
    Vector4d norm = calcTriangleNormal(ps->vertices, ind);
    double d = norm.head<3>().dot(vec);
    if (d > dmax) {
      dmax = d;
      vectormax = norm;
    }
  }
  PyObject *coord = pf.PyList_New(4);
  for (int i = 0; i < 4; i++) pf.PyList_SetItem(coord, i, pf.PyFloat_FromDouble(vectormax[i]));
  return coord;
}

PyObject *python_find_face(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"obj", "vec", NULL};
  PyObject *obj = nullptr;
  PyObject *vec = nullptr;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "OO", kwlist, &obj, &vec)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing find_face(object)");
    return NULL;
  }
  return python_find_face_core(obj, vec);
}

PyObject *python_oo_find_face(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"vec", NULL};
  PyObject *vec = nullptr;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &vec)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing find_face(object)");
    return NULL;
  }
  return python_find_face_core(obj, vec);
}

PyObject *python_sitonto_core(PyObject *pyobj, PyObject *vecx_p, PyObject *vecy_p, PyObject *vecz_p)
{
  Vector4d vecx, vecy, vecz;
  Vector3d cut;
  if (python_vectorval(vecx_p, 3, 3, &vecx[0], &vecx[1], &vecx[2], &vecx[3])) return Py_NONE;
  if (python_vectorval(vecy_p, 3, 3, &vecy[0], &vecy[1], &vecy[2], &vecy[3])) return Py_NONE;
  if (python_vectorval(vecz_p, 3, 3, &vecz[0], &vecz[1], &vecz[2], &vecz[3])) return Py_NONE;

  if (cut_face_face_face(vecx.head<3>() * vecx[3], vecx.head<3>(), vecy.head<3>() * vecy[3],
                         vecy.head<3>(), vecz.head<3>() * vecz[3], vecz.head<3>(), cut))
    return Py_NONE;
  DECLARE_INSTANCE
  auto node = std::make_shared<TransformNode>(instance, "sitontonode");
  std::shared_ptr<AbstractNode> child;
  PyObject *dummydict;
  PyTypeObject *type = PyOpenSCADObjectType(pyobj);
  child = PyOpenSCADObjectToNodeMulti(pyobj, &dummydict);
  node->setPyName(child->getPyName());
  if (child == NULL) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid type for Object in multmatrix");
    return NULL;
  }
  Matrix4d mat;
  Vector3d vecx_n = vecz.head<3>().cross(vecy.head<3>()).normalized();
  Vector3d vecy_n = vecx.head<3>().cross(vecz.head<3>()).normalized();
  if (vecx.head<3>().cross(vecy.head<3>()).dot(vecz.head<3>()) < 0) {
    vecx_n = -vecx_n;
    vecy_n = -vecy_n;
  }
  mat << vecx_n[0], vecy_n[0], vecz[0], cut[0], vecx_n[1], vecy_n[1], vecz[1], cut[1], vecx_n[2],
    vecy_n[2], vecz[2], cut[2], 0, 0, 0, 1;

  node->matrix = mat;
  node->children.push_back(child);
  return PyOpenSCADObjectFromNode(type, node);
}

PyObject *python_sitonto(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"obj", "vecz", "vecx", "vecy", NULL};
  PyObject *obj = nullptr;
  PyObject *vecx = nullptr;
  PyObject *vecy = nullptr;
  PyObject *vecz = nullptr;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "OOOO", kwlist, &obj, &vecz, &vecx, &vecy)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing sitonto(object)");
    return NULL;
  }
  return python_sitonto_core(obj, vecx, vecy, vecz);
}

PyObject *python_oo_sitonto(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"vecz", " vecx", "vecy", NULL};
  PyObject *vecx = nullptr;
  PyObject *vecy = nullptr;
  PyObject *vecz = nullptr;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "OOO", kwlist, &vecz, &vecx, &vecy)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing sitonto(object)");
    return NULL;
  }
  return python_sitonto_core(obj, vecx, vecy, vecz);
}

PyObject *python__getitem__(PyObject *obj, PyObject *key)
{
  PyOpenSCADObject *self = (PyOpenSCADObject *)obj;
  if (self->dict == nullptr) {
    return nullptr;
  }
  PyObject *result = pf.PyDict_GetItem(self->dict, key);
  if (result == NULL) {
    PyObject *keyname = pf.PyUnicode_AsEncodedString(key, "utf-8", "~");
    if (keyname == nullptr) return nullptr;
    std::string keystr = PyBytes_AS_STRING(keyname);
    result = Py_NONE;
    if (keystr == "matrix") {
      PyObject *dummy_dict;
      std::shared_ptr<AbstractNode> node = PyOpenSCADObjectToNode(obj, &dummy_dict);
      std::shared_ptr<const TransformNode> trans = std::dynamic_pointer_cast<const TransformNode>(node);
      Matrix4d matrix = Matrix4d::Identity();
      if (trans != nullptr) matrix = trans->matrix.matrix();
      result = python_frommatrix(matrix);
    } else if (keystr == "size") {
      return python_size_core(obj);
    } else if (keystr == "position") {
      return python_position_core(obj);
    } else if (keystr == "bbox") {
      return python_bbox_core(obj);
    }
  } else Py_INCREF(result);
  return result;
}

int python__setitem__(PyObject *dict, PyObject *key, PyObject *v)
{
  PyOpenSCADObject *self = (PyOpenSCADObject *)dict;
  if (self->dict == NULL) {
    return 0;
  }
  Py_INCREF(v);
  pf.PyDict_SetItem(self->dict, key, v);
  return 0;
}

PyObject *python_color_core(PyObject *obj, PyObject *color, double alpha)
{
  PyObject *child_dict;
  std::shared_ptr<AbstractNode> child;
  PyTypeObject *type = PyOpenSCADObjectType(obj);
  child = PyOpenSCADObjectToNodeMulti(obj, &child_dict);
  if (child == NULL) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid type for Object in color");
    return NULL;
  }
  DECLARE_INSTANCE
  auto node = std::make_shared<ColorNode>(instance);

  Vector4d col(0, 0, 0, alpha);
  if (!python_vectorval(color, 3, 4, &col[0], &col[1], &col[2], &col[3])) {
    node->color.setRgba(float(col[0]), float(col[1]), float(col[2]), float(col[3]));
  } else if (PyUnicode_Check(color)) {
    PyObject *value = pf.PyUnicode_AsEncodedString(color, "utf-8", "~");
    char *colorname = PyBytes_AS_STRING(value);
    const auto color = OpenSCAD::parse_color(colorname);
    if (color) {
      node->color = *color;
      node->color.setAlpha(alpha);
    } else {
      pf.PyErr_SetString(pf.PyExc_TypeError, "Cannot parse color");
      return NULL;
    }
  } else {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Unknown color representation");
    return nullptr;
  }

  node->textureind = -1;
  node->children.push_back(child);

  PyObject *pyresult = PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
  if (child_dict != nullptr) {
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (pf.PyDict_Next(child_dict, &pos, &key, &value)) {
      pf.PyDict_SetItem(((PyOpenSCADObject *)pyresult)->dict, key, value);
    }
  }
  return pyresult;
}

PyObject *python_color(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"obj", "c", "alpha", NULL};
  PyObject *obj = NULL;
  PyObject *color = NULL;
  double alpha = 1.0;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O|Odi", kwlist, &obj, &color, &alpha)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error during parsing color");
    return NULL;
  }
  return python_color_core(obj, color, alpha);
}

PyObject *python_oo_color(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"c", "alpha", NULL};
  PyObject *color = NULL;
  double alpha = 1.0;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "|Odi", kwlist, &color, &alpha)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error during parsing color");
    return NULL;
  }
  return python_color_core(obj, color, alpha);
}

typedef std::vector<int> intList;

PyObject *python_mesh_core(PyObject *obj, bool tessellate)
{
  PyObject *dummydict;
  std::shared_ptr<AbstractNode> child = PyOpenSCADObjectToNodeMulti(obj, &dummydict);
  if (child == NULL) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid type for  Object in mesh \n");
    return NULL;
  }
  Tree tree(child, "");
  GeometryEvaluator geomevaluator(tree);
  std::shared_ptr<const Geometry> geom = geomevaluator.evaluateGeometry(*tree.root(), true);
  std::shared_ptr<const PolySet> ps = PolySetUtils::getGeometryAsPolySet(geom);

  if (ps != nullptr) {
    if (tessellate == true) {
      ps = PolySetUtils::tessellate_faces(*ps);
    }
    // Now create Python Point array
    PyObject *ptarr = pf.PyList_New(ps->vertices.size());
    for (unsigned int i = 0; i < ps->vertices.size(); i++) {
      PyObject *coord = pf.PyList_New(3);
      for (int j = 0; j < 3; j++) pf.PyList_SetItem(coord, j, pf.PyFloat_FromDouble(ps->vertices[i][j]));
      pf.PyList_SetItem(ptarr, i, coord);
      Py_XINCREF(coord);
    }
    Py_XINCREF(ptarr);
    // Now create Python Point array
    PyObject *polarr = pf.PyList_New(ps->indices.size());
    for (unsigned int i = 0; i < ps->indices.size(); i++) {
      PyObject *face = pf.PyList_New(ps->indices[i].size());
      for (unsigned int j = 0; j < ps->indices[i].size(); j++)
        pf.PyList_SetItem(face, j, pf.PyLong_FromLong(ps->indices[i][j]));
      pf.PyList_SetItem(polarr, i, face);
      Py_XINCREF(face);
    }
    Py_XINCREF(polarr);

    PyObject *result = pf.PyTuple_New(2);
    pf.PyTuple_SetItem(result, 0, ptarr);
    pf.PyTuple_SetItem(result, 1, polarr);

    return result;
  }
  if (auto polygon2d = std::dynamic_pointer_cast<const Polygon2d>(geom)) {
    const std::vector<Outline2d> outlines = polygon2d->outlines();
    PyObject *pyth_outlines = pf.PyList_New(outlines.size());
    for (unsigned int i = 0; i < outlines.size(); i++) {
      const Outline2d& outline = outlines[i];
      PyObject *pyth_outline = pf.PyList_New(outline.vertices.size());
      for (unsigned int j = 0; j < outline.vertices.size(); j++) {
        Vector2d pt = outline.vertices[j];
        PyObject *pyth_pt = pf.PyList_New(2);
        for (int k = 0; k < 2; k++) pf.PyList_SetItem(pyth_pt, k, pf.PyFloat_FromDouble(pt[k]));
        pf.PyList_SetItem(pyth_outline, j, pyth_pt);
      }
      pf.PyList_SetItem(pyth_outlines, i, pyth_outline);
    }
    return pyth_outlines;
  }
  return Py_NONE;
}

PyObject *python_mesh(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"obj", "triangulate", NULL};
  PyObject *obj = NULL;
  PyObject *tess = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O|O", kwlist, &obj, &tess)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error during parsing\n");
    return NULL;
  }
  return python_mesh_core(obj, tess == Py_TRUE);
}

PyObject *python_oo_mesh(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"triangulate", NULL};
  PyObject *tess = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "|O", kwlist, &tess)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error during parsing\n");
    return NULL;
  }
  return python_mesh_core(obj, tess == Py_TRUE);
}

PyObject *python_bbox_core(PyObject *obj)
{
  // Get position and size attributes from the object
  PyObject *position_key = PyUnicode_FromString("position");
  PyObject *size_key = PyUnicode_FromString("size");

  PyObject *position = python__getitem__(obj, position_key);
  PyObject *size = python__getitem__(obj, size_key);

  Py_DECREF(position_key);
  Py_DECREF(size_key);

  if (position == Py_NONE || size == Py_NONE) {
    if (position != Py_NONE) Py_DECREF(position);
    if (size != Py_NONE) Py_DECREF(size);
    return Py_NONE;
  }

  // Create cube with the size, not centered (starts at origin [0,0,0])
  PyObject *cube_args = pf.PyTuple_New(0);
  PyObject *cube_kwargs = pf.PyDict_New();
  pf.PyDict_SetItemString(cube_kwargs, "size", size);
  pf.PyDict_SetItemString(cube_kwargs, "center", Py_FALSE);

  PyObject *cube = python_cube(NULL, cube_args, cube_kwargs);
  if (cube == NULL) {
    Py_DECREF(position);
    Py_DECREF(size);
    Py_DECREF(cube_args);
    Py_DECREF(cube_kwargs);
    return NULL;
  }

  // Translate cube to the object's position
  PyObject *bbox_box = python_translate_core(cube, position);

  // Clean up
  Py_DECREF(position);
  Py_DECREF(size);
  Py_DECREF(cube_args);
  Py_DECREF(cube_kwargs);
  Py_DECREF(cube);

  return bbox_box;
}

PyObject *python_bbox(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"obj", NULL};
  PyObject *obj = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &obj)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error during parsing\n");
    return NULL;
  }
  return python_bbox_core(obj);
}

PyObject *python_oo_bbox(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {NULL};
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "", kwlist)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error during parsing\n");
    return NULL;
  }
  return python_bbox_core(obj);
}

PyObject *python_size_core(PyObject *obj)
{
  PyObject *dummydict;
  std::shared_ptr<AbstractNode> child = PyOpenSCADObjectToNodeMulti(obj, &dummydict);
  if (child == NULL) {
    return Py_NONE;
  }
  Tree tree(child, "");
  GeometryEvaluator geomevaluator(tree);
  std::shared_ptr<const Geometry> geom = geomevaluator.evaluateGeometry(*tree.root(), true);
  std::shared_ptr<const PolySet> ps = PolySetUtils::getGeometryAsPolySet(geom);

  if (ps != nullptr && ps->vertices.size() > 0) {
    Vector3d pmin = ps->vertices[0];
    Vector3d pmax = pmin;
    for (const auto& pt : ps->vertices) {
      for (int i = 0; i < 3; i++) {
        if (pt[i] > pmax[i]) pmax[i] = pt[i];
        if (pt[i] < pmin[i]) pmin[i] = pt[i];
      }
    }
    // Calculate size directly
    Vector3d size = pmax - pmin;
    PyObject *size_list = python_fromvector(size);
    Py_INCREF(size_list);
    return size_list;
  }
  return Py_NONE;
}

PyObject *python_position_core(PyObject *obj)
{
  PyObject *dummydict;
  std::shared_ptr<AbstractNode> child = PyOpenSCADObjectToNodeMulti(obj, &dummydict);
  if (child == NULL) {
    return Py_NONE;
  }
  Tree tree(child, "");
  GeometryEvaluator geomevaluator(tree);
  std::shared_ptr<const Geometry> geom = geomevaluator.evaluateGeometry(*tree.root(), true);
  std::shared_ptr<const PolySet> ps = PolySetUtils::getGeometryAsPolySet(geom);

  if (ps != nullptr && ps->vertices.size() > 0) {
    Vector3d pmin = ps->vertices[0];
    for (const auto& pt : ps->vertices) {
      for (int i = 0; i < 3; i++) {
        if (pt[i] < pmin[i]) pmin[i] = pt[i];
      }
    }
    PyObject *position_list = python_fromvector(pmin);
    Py_INCREF(position_list);
    return position_list;
  }
  return Py_NONE;
}

PyObject *python_size(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"obj", NULL};
  PyObject *obj = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &obj)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error during parsing size(obj)\n");
    return NULL;
  }
  return python_size_core(obj);
}

PyObject *python_position(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"obj", NULL};
  PyObject *obj = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &obj)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error during parsing position(obj)\n");
    return NULL;
  }
  return python_position_core(obj);
}

PyObject *python_separate_core(PyObject *obj)
{
  PyObject *dummydict;
  PyTypeObject *type = PyOpenSCADObjectType(obj);
  std::shared_ptr<AbstractNode> child = PyOpenSCADObjectToNodeMulti(obj, &dummydict);
  if (child == NULL) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid type for  Object in separate \n");
    return NULL;
  }
  Tree tree(child, "");
  GeometryEvaluator geomevaluator(tree);
  std::shared_ptr<const Geometry> geom = geomevaluator.evaluateGeometry(*tree.root(), true);
  std::shared_ptr<const PolySet> ps = PolySetUtils::getGeometryAsPolySet(geom);

  if (ps != nullptr) {
    // setup databases
    intList empty_list;
    std::vector<intList> pt2tri;

    std::vector<int> vert_db;
    for (size_t i = 0; i < ps->vertices.size(); i++) {
      vert_db.push_back(-1);
      pt2tri.push_back(empty_list);
    }

    std::vector<int> tri_db;
    for (size_t i = 0; i < ps->indices.size(); i++) {
      tri_db.push_back(-1);
      for (auto ind : ps->indices[i]) pt2tri[ind].push_back(i);
    }

    // now sort for objects
    int obj_num = 0;
    for (size_t i = 0; i < vert_db.size(); i++) {
      if (vert_db[i] != -1) continue;
      std::vector<int> vert_todo;
      vert_todo.push_back(i);
      while (vert_todo.size() > 0) {
        int vert_ind = vert_todo[vert_todo.size() - 1];
        vert_todo.pop_back();
        if (vert_db[vert_ind] != -1) continue;
        vert_db[vert_ind] = obj_num;
        for (int tri_ind : pt2tri[vert_ind]) {
          if (tri_db[tri_ind] != -1) continue;
          tri_db[tri_ind] = obj_num;
          for (int vert1_ind : ps->indices[tri_ind]) {
            if (vert_db[vert1_ind] != -1) continue;
            vert_todo.push_back(vert1_ind);
          }
        }
      }
      obj_num++;
    }

    PyObject *objects = pf.PyList_New(obj_num);
    for (int i = 0; i < obj_num; i++) {
      // create a polyhedron for each
      DECLARE_INSTANCE
      auto node = std::make_shared<PolyhedronNode>(instance);
      node->convexity = 2;
      std::vector<int> vert_map;
      for (size_t j = 0; j < ps->vertices.size(); j++) {
        if (vert_db[j] == i) {
          vert_map.push_back(node->points.size());
          node->points.push_back(ps->vertices[j]);
        } else vert_map.push_back(-1);
      }
      for (size_t j = 0; j < ps->indices.size(); j++) {
        if (tri_db[j] == i) {
          IndexedFace face_map;
          for (auto ind : ps->indices[j]) {
            face_map.push_back(vert_map[ind]);
          }
          std::reverse(face_map.begin(), face_map.end());
          node->faces.push_back(face_map);
        }
      }
      pf.PyList_SetItem(objects, i, PyOpenSCADObjectFromNode(&PyOpenSCADType, node));
    }
    return objects;
  }
  return Py_NONE;
}

PyObject *python_separate(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"obj", NULL};
  PyObject *obj = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &obj)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error during parsing\n");
    return NULL;
  }
  return python_separate_core(obj);
}

PyObject *python_oo_separate(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {NULL};
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "", kwlist)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error during parsing\n");
    return NULL;
  }
  return python_separate_core(obj);
}

PyObject *python_edges_core(PyObject *obj)
{
  PyObject *dummydict;
  PyTypeObject *type = PyOpenSCADObjectType(obj);
  std::shared_ptr<AbstractNode> child = PyOpenSCADObjectToNodeMulti(obj, &dummydict);
  if (child == NULL) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid type for  Object in faces \n");
    return NULL;
  }

  Tree tree(child, "");
  GeometryEvaluator geomevaluator(tree);
  std::shared_ptr<const Geometry> geom = geomevaluator.evaluateGeometry(*tree.root(), true);
  const std::shared_ptr<const Polygon2d> poly = std::dynamic_pointer_cast<const Polygon2d>(geom);
  if (poly == nullptr) return Py_NONE;
  int edgenum = 0;
  Vector3d zdir(0, 0, 0);
  Transform3d trans = poly->getTransform3d();
  for (auto ol : poly->untransformedOutlines()) {
    int n = ol.vertices.size();
    for (int i = 0; i < n; i++) {
      Vector3d p1 = trans * Vector3d(ol.vertices[i][0], ol.vertices[i][1], 0);
      Vector3d p2 = trans * Vector3d(ol.vertices[(i + 1) % n][0], ol.vertices[(i + 1) % n][1], 0);
      Vector3d p3 = trans * Vector3d(ol.vertices[(i + 2) % n][0], ol.vertices[(i + 2) % n][1], 0);
      zdir += (p1 - p2).cross(p2 - p3);
    }
    edgenum += n;
  }
  zdir.normalize();
  PyObject *pyth_edges = pf.PyList_New(edgenum);
  int ind = 0;

  for (auto ol : poly->untransformedOutlines()) {
    int n = ol.vertices.size();
    for (int i = 0; i < n; i++) {
      Vector3d p1 = trans * Vector3d(ol.vertices[i][0], ol.vertices[i][1], 0);
      Vector3d p2 = trans * Vector3d(ol.vertices[(i + 1) % n][0], ol.vertices[(i + 1) % n][1], 0);
      Vector3d pt = (p1 + p2) / 2.0;
      Vector3d xdir = (p2 - p1).normalized();
      Vector3d ydir = xdir.cross(zdir).normalized();

      Matrix4d mat;
      mat << xdir[0], ydir[0], zdir[0], pt[0], xdir[1], ydir[1], zdir[1], pt[1], xdir[2], ydir[2],
        zdir[2], pt[2], 0, 0, 0, 1;

      DECLARE_INSTANCE
      auto edge = std::make_shared<EdgeNode>(instance);
      edge->size = (p2 - p1).norm();
      edge->center = true;
      {
        DECLARE_INSTANCE
        auto mult = std::make_shared<TransformNode>(instance, "multmatrix");
        mult->matrix = mat;
        mult->children.push_back(edge);

        PyObject *pyth_edge = PyOpenSCADObjectFromNode(&PyOpenSCADType, mult);
        pf.PyList_SetItem(pyth_edges, ind++, pyth_edge);
      }
    }
  }
  return pyth_edges;
}

PyObject *python_edges(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"obj", NULL};
  PyObject *obj = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &obj)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error during parsing\n");
    return NULL;
  }
  return python_edges_core(obj);
}

PyObject *python_oo_edges(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {NULL};
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "", kwlist)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error during parsing\n");
    return NULL;
  }
  return python_edges_core(obj);
}

PyObject *python_faces_core(PyObject *obj, bool tessellate)
{
  PyObject *dummydict;
  PyTypeObject *type = PyOpenSCADObjectType(obj);
  std::shared_ptr<AbstractNode> child = PyOpenSCADObjectToNodeMulti(obj, &dummydict);
  if (child == NULL) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid type for  Object in faces \n");
    return NULL;
  }

  Tree tree(child, "");
  GeometryEvaluator geomevaluator(tree);
  std::shared_ptr<const Geometry> geom = geomevaluator.evaluateGeometry(*tree.root(), true);
  std::shared_ptr<const PolySet> ps = PolySetUtils::getGeometryAsPolySet(geom);

  if (ps != nullptr) {
    PolygonIndices inds;
    std::vector<int> face_parents;
    if (tessellate == true) {
      ps = PolySetUtils::tessellate_faces(*ps);
      inds = ps->indices;
      for (size_t i = 0; i < inds.size(); i++) face_parents.push_back(-1);
    } else {
      std::vector<Vector4d> normals, new_normals;
      normals = calcTriangleNormals(ps->vertices, ps->indices);
      inds = mergeTriangles(ps->indices, normals, new_normals, face_parents, ps->vertices);
    }
    int resultlen = 0, resultiter = 0;
    for (size_t i = 0; i < face_parents.size(); i++)
      if (face_parents[i] == -1) resultlen++;

    PyObject *pyth_faces = pf.PyList_New(resultlen);

    for (size_t j = 0; j < inds.size(); j++) {
      if (face_parents[j] != -1) continue;
      auto& face = inds[j];
      if (face.size() < 3) continue;
      Vector3d zdir = calcTriangleNormal(ps->vertices, face).head<3>().normalized();
      // calc center of face
      Vector3d ptmin, ptmax;
      for (size_t i = 0; i < face.size(); i++) {
        Vector3d pt = ps->vertices[face[i]];
        for (int k = 0; k < 3; k++) {
          if (i == 0 || pt[k] < ptmin[k]) ptmin[k] = pt[k];
          if (i == 0 || pt[k] > ptmax[k]) ptmax[k] = pt[k];
        }
      }
      Vector3d pt =
        Vector3d((ptmin[0] + ptmax[0]) / 2.0, (ptmin[1] + ptmax[1]) / 2.0, (ptmin[2] + ptmax[2]) / 2.0);
      Vector3d xdir = (ps->vertices[face[1]] - ps->vertices[face[0]]).normalized();
      Vector3d ydir = zdir.cross(xdir);

      Matrix4d mat;
      mat << xdir[0], ydir[0], zdir[0], pt[0], xdir[1], ydir[1], zdir[1], pt[1], xdir[2], ydir[2],
        zdir[2], pt[2], 0, 0, 0, 1;

      Matrix4d invmat = mat.inverse();

      DECLARE_INSTANCE
      auto poly = std::make_shared<PolygonNode>(instance);
      std::vector<size_t> path;
      for (size_t i = 0; i < face.size(); i++) {
        Vector3d pt = ps->vertices[face[i]];
        Vector4d pt4(pt[0], pt[1], pt[2], 1);
        pt4 = invmat * pt4;
        path.push_back(poly->points.size());
        Vector3d pt3 = pt4.head<3>();
        pt3[2] = 0;  // no radius
        poly->points.push_back(pt3);
      }

      // xdir should be parallel to the longest edge
      int n = poly->points.size();
      Vector3d maxline(0, 0, 0);
      for (int i = 0; i < n; i++) {
        Vector3d line = poly->points[(i + 1) % n] - poly->points[i];
        if (line.norm() > maxline.norm()) maxline = line;
      }
      double arcshift = atan2(maxline[1], maxline[0]);
      Vector3d newpt(0, 0, 0);
      for (auto& pt : poly->points) {
        newpt[0] = pt[0] * cos(-arcshift) - pt[1] * sin(-arcshift);
        newpt[1] = pt[0] * sin(-arcshift) + pt[1] * cos(-arcshift);
        pt = newpt;
      }

      poly->paths.push_back(path);

      // check if there are holes
      for (size_t k = 0; k < inds.size(); k++) {
        if ((size_t)face_parents[k] == j) {
          auto& hole = inds[k];

          std::vector<size_t> path;
          for (size_t i = 0; i < hole.size(); i++) {
            Vector3d pt = ps->vertices[hole[i]];
            Vector4d pt4(pt[0], pt[1], pt[2], 1);
            pt4 = invmat * pt4;
            path.push_back(poly->points.size());
            Vector3d pt3 = pt4.head<3>();
            pt3[2] = 0;  // no radius
            poly->points.push_back(pt3);
          }
          poly->paths.push_back(path);
        }
      }
      {
        DECLARE_INSTANCE
        auto mult = std::make_shared<TransformNode>(instance, "multmatrix");
        // mat um arcshift drehen
        Vector3d xvec(mat(0, 0), mat(1, 0), mat(2, 0));
        Vector3d yvec(mat(0, 1), mat(1, 1), mat(2, 1));

        PyObject *pyth_face = PyOpenSCADObjectFromNode(&PyOpenSCADType, mult);
        pf.PyList_SetItem(pyth_faces, resultiter++, pyth_face);
      }
    }
    return pyth_faces;
  }
  return Py_NONE;
}

PyObject *python_faces(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"obj", "triangulate", NULL};
  PyObject *obj = NULL;
  PyObject *tess = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O|O", kwlist, &obj, &tess)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error during parsing\n");
    return NULL;
  }
  return python_faces_core(obj, tess == Py_TRUE);
}

PyObject *python_oo_faces(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"triangulate", NULL};
  PyObject *tess = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "|O", kwlist, &tess)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error during parsing\n");
    return NULL;
  }
  return python_faces_core(obj, tess == Py_TRUE);
}

PyObject *python_oversample_core(PyObject *obj, int n, PyObject *round)
{
  PyObject *dummydict;
  PyTypeObject *type = PyOpenSCADObjectType(obj);
  std::shared_ptr<AbstractNode> child = PyOpenSCADObjectToNodeMulti(obj, &dummydict);
  if (child == NULL) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid type for  Object in oversample \n");
    return NULL;
  }

  DECLARE_INSTANCE
  auto node = std::make_shared<OversampleNode>(instance);
  node->children.push_back(child);
  node->n = n;
  node->round = 0;
  if (round == Py_TRUE) node->round = 1;

  return PyOpenSCADObjectFromNode(type, node);
}

PyObject *python_oversample(PyObject *self, PyObject *args, PyObject *kwargs)
{
  int n = 2;
  char *kwlist[] = {"obj", "n", "round", NULL};
  PyObject *obj = NULL;
  PyObject *round = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "Oi|O", kwlist, &obj, &n, &round)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error during parsing\n");
    return NULL;
  }
  return python_oversample_core(obj, n, round);
}

PyObject *python_oo_oversample(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  int n = 2;
  char *kwlist[] = {"n", "round", NULL};
  PyObject *round = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "i|O", kwlist, &n, &round)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error during parsing\n");
    return NULL;
  }
  return python_oversample_core(obj, n, round);
}

PyObject *python_debug_core(PyObject *obj, PyObject *faces)
{
  PyObject *dummydict;
  PyTypeObject *type = PyOpenSCADObjectType(obj);
  std::shared_ptr<AbstractNode> child = PyOpenSCADObjectToNodeMulti(obj, &dummydict);
  if (child == NULL) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid type for  Object in debug \n");
    return NULL;
  }

  DECLARE_INSTANCE
  auto node = std::make_shared<DebugNode>(instance);
  node->children.push_back(child);
  if (faces != nullptr) {
    std::vector<int> intfaces = python_intlistval(faces);
    node->faces = intfaces;
  }
  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

PyObject *python_debug(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"obj", "faces", NULL};
  PyObject *obj = NULL;
  PyObject *faces = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O|O", kwlist, &obj, &faces)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error duing parsing\n");
    return NULL;
  }
  return python_debug_core(obj, faces);
}

PyObject *python_oo_debug(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"faces", NULL};
  PyObject *faces = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "|O", kwlist, &faces)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error duing parsing\n");
    return NULL;
  }
  return python_debug_core(self, faces);
}

PyObject *python_repair_core(PyObject *obj, PyObject *color)
{
  PyObject *dummydict;
  PyTypeObject *type = PyOpenSCADObjectType(obj);
  std::shared_ptr<AbstractNode> child = PyOpenSCADObjectToNode(obj, &dummydict);
  if (child == NULL) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid type for  Object in repair \n");
    return NULL;
  }

  DECLARE_INSTANCE
  auto node = std::make_shared<RepairNode>(instance);
  node->children.push_back(child);
  if (color != nullptr) {
    Vector4d col(0, 0, 0, 1.0);
    if (!python_vectorval(color, 3, 4, &col[0], &col[1], &col[2], &col[3])) {
      node->color.setRgba(float(col[0]), float(col[1]), float(col[2]), float(col[3]));
    } else if (PyUnicode_Check(color)) {
      PyObject *value = PyUnicode_AsEncodedString(color, "utf-8", "~");
      char *colorname = PyBytes_AS_STRING(value);
      const auto color = OpenSCAD::parse_color(colorname);
      if (color) {
        node->color = *color;
        node->color.setAlpha(1.0);
      } else {
        pf.PyErr_SetString(pf.PyExc_TypeError, "Cannot parse color");
        return NULL;
      }
    } else {
      pf.PyErr_SetString(pf.PyExc_TypeError, "Unknown color representation");
      return nullptr;
    }
  }
  return PyOpenSCADObjectFromNode(type, node);
}

PyObject *python_repair(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"obj", "color", NULL};
  PyObject *obj = NULL;
  PyObject *color = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O|O", kwlist, &obj, &color)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error duing parsing\n");
    return NULL;
  }
  return python_repair_core(obj, color);
}

PyObject *python_oo_repair(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"color", NULL};
  PyObject *color = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "|O", kwlist, &color)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error duing parsing\n");
    return NULL;
  }
  return python_repair_core(self, color);
}

PyObject *python_fillet_core(PyObject *obj, double r, int fn, PyObject *sel, double minang)
{
  PyObject *dummydict;
  DECLARE_INSTANCE
  auto node = std::make_shared<FilletNode>(instance);
  PyTypeObject *type = &PyOpenSCADType;
  node->r = r;
  node->fn = fn;
  node->minang = minang;
  if (obj != nullptr) node->children.push_back(PyOpenSCADObjectToNodeMulti(obj, &dummydict));
  else {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid type for  Object in fillet \n");
    return NULL;
  }

  if (sel != nullptr) node->children.push_back(PyOpenSCADObjectToNodeMulti(sel, &dummydict));

  return PyOpenSCADObjectFromNode(type, node);
}

PyObject *python_fillet(PyObject *self, PyObject *args, PyObject *kwargs)
{
  double r = 1.0;
  double fn = NAN;
  double minang = 30;
  char *kwlist[] = {"obj", "r", "sel", "n", "minang", NULL};
  PyObject *obj = NULL;
  PyObject *sel = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "Od|Odd", kwlist, &obj, &r, &sel, &fn, &minang)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error during parsing\n");
    return NULL;
  }
  double dummy;
  if (isnan(fn)) get_fnas(fn, dummy, dummy);
  return python_fillet_core(obj, r, fn, sel, minang);
}

PyObject *python_oo_fillet(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  double r = 1.0;
  double fn = NAN;
  double minang = 30;
  PyObject *sel = nullptr;
  char *kwlist[] = {"r", "sel", "fn", "minang", NULL};
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "d|Odd", kwlist, &r, &sel, &fn, &minang)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error during parsing\n");
    return NULL;
  }
  double dummy;
  if (isnan(fn)) get_fnas(fn, dummy, dummy);
  return python_fillet_core(obj, r, fn, sel, minang);
}

PyObject *rotate_extrude_core(PyObject *obj, int convexity, double scale, double angle, PyObject *twist,
                              PyObject *origin, PyObject *offset, PyObject *vp, char *method, double fn,
                              double fa, double fs)
{
  DECLARE_INSTANCE
  std::shared_ptr<AbstractNode> child;
  auto node = std::make_shared<RotateExtrudeNode>(instance);
  PyTypeObject *type = &PyOpenSCADType;
  node->profile_func = NULL;
  node->twist_func = NULL;
  if (obj->ob_type == pf.PyFunction_Type) {
    Py_XINCREF(obj);  // TODO there to decref it ?
    node->profile_func = obj;
    auto dummy_node = std::make_shared<SquareNode>(instance);
    node->children.push_back(dummy_node);
  } else {
    PyObject *dummydict;
    type = PyOpenSCADObjectType(obj);
    child = PyOpenSCADObjectToNodeMulti(obj, &dummydict);
    if (child == NULL) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid type for  Object in rotate_extrude\n");
      return NULL;
    }
    node->children.push_back(child);
  }

  get_fnas(node->fn, node->fa, node->fs);
  if (!isnan(fn)) node->fn = fn;
  if (!isnan(fa)) node->fa = fa;
  if (!isnan(fs)) node->fs = fs;

  node->convexity = convexity;
  node->scale = scale;
  node->angle = angle;
  if (twist != NULL) {
    if (twist->ob_type == pf.PyFunction_Type) {
      Py_XINCREF(twist);  // TODO there to decref it ?
      node->twist_func = twist;
    } else node->twist = pf.PyFloat_AsDouble(twist);
  }

  if (origin != NULL && PyList_CHECK(origin) && pf.PyList_Size(origin) == 2) {
    node->origin_x = pf.PyFloat_AsDouble(pf.PyList_GetItem(origin, 0));
    node->origin_y = pf.PyFloat_AsDouble(pf.PyList_GetItem(origin, 1));
  }
  if (offset != NULL && PyList_CHECK(offset) && pf.PyList_Size(offset) == 2) {
    node->offset_x = pf.PyFloat_AsDouble(pf.PyList_GetItem(offset, 0));
    node->offset_y = pf.PyFloat_AsDouble(pf.PyList_GetItem(offset, 1));
  }
  double dummy;
  Vector3d v(0, 0, 0);
  if (vp != nullptr && !python_vectorval(vp, 3, 3, &v[0], &v[1], &v[2], &dummy)) {
  }
  node->v = v;
  if (method != nullptr) node->method = method;
  else node->method = "centered";

  if (node->convexity <= 0) node->convexity = 2;
  if (node->scale <= 0) node->scale = 1;
  if (node->angle <= -360) node->angle = 360;

  return PyOpenSCADObjectFromNode(type, node);
}

PyObject *python_rotate_extrude(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyObject *obj = NULL;
  int convexity = 1;
  double scale = 1.0;
  double angle = 360.0;
  PyObject *twist = NULL;
  PyObject *v = NULL;
  char *method = NULL;
  PyObject *origin = NULL;
  PyObject *offset = NULL;
  double fn = NAN, fa = NAN, fs = NAN;
  get_fnas(fn, fa, fs);
  char *kwlist[] = {"obj", "convexity", "scale", "angle", "twist", "origin", "offset",
                    "v",   "method",    "fn",    "fa",    "fs",    NULL};
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O|iddOOOOsddd", kwlist, &obj, &convexity, &scale,
                                      &angle, &twist, &origin, &offset, &v, &method, &fn, &fa, &fs)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing rotate_extrude(object,...)");
    return NULL;
  }
  return rotate_extrude_core(obj, convexity, scale, angle, twist, origin, offset, v, method, fn, fa, fs);
}

PyObject *python_oo_rotate_extrude(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  int convexity = 1;
  double scale = 1.0;
  double angle = 360.0;
  PyObject *twist = NULL;
  PyObject *origin = NULL;
  PyObject *offset = NULL;
  double fn = NAN, fa = NAN, fs = NAN;
  get_fnas(fn, fa, fs);
  PyObject *v = NULL;
  char *method = NULL;
  char *kwlist[] = {"convexity", "scale",  "angle", "twist", "origin", "offset",
                    "v",         "method", "fn",    "fa",    "fs",     NULL};
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "|iddOOOOsddd", kwlist, &convexity, &scale, &angle,
                                      &twist, &origin, &offset, &v, &method, &fn, &fa, &fs)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error during parsing\n");
    return NULL;
  }
  return rotate_extrude_core(obj, convexity, scale, angle, twist, origin, offset, v, method, fn, fa, fs);
}

PyObject *linear_extrude_core(PyObject *obj, PyObject *height, int convexity, PyObject *origin,
                              PyObject *scale, PyObject *center, int slices, int segments,
                              PyObject *twist, double fn, double fa, double fs)
{
  DECLARE_INSTANCE
  std::shared_ptr<AbstractNode> child;
  auto node = std::make_shared<LinearExtrudeNode>(instance);
  PyTypeObject *type = &PyOpenSCADType;

  node->profile_func = NULL;
  node->twist_func = NULL;
  get_fnas(node->fn, node->fa, node->fs);
  if (obj->ob_type == pf.PyFunction_Type) {
    Py_XINCREF(obj);  // TODO there to decref it ?
    node->profile_func = obj;
    node->fn = 2;
    auto dummy_node = std::make_shared<SquareNode>(instance);
    node->children.push_back(dummy_node);
  } else {
    PyObject *dummydict;
    child = PyOpenSCADObjectToNodeMulti(obj, &dummydict);
    if (child == NULL) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid type for  Object in linear_extrude\n");
      return NULL;
    }
    node->children.push_back(child);
  }

  if (!isnan(fn)) node->fn = fn;
  if (!isnan(fa)) node->fa = fa;
  if (!isnan(fs)) node->fs = fs;

  Vector3d height_vec(0, 0, 0);
  double dummy;
  if (!python_numberval(height, &height_vec[2])) {
    node->height = height_vec;
  } else if (!python_vectorval(height, 3, 3, &height_vec[0], &height_vec[1], &height_vec[2], &dummy)) {
    node->height = height_vec;
    node->has_heightvector = true;
  } else {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Height must be either a number or a vector\n");
    return NULL;
  }

  node->convexity = convexity;

  node->origin_x = 0.0;
  node->origin_y = 0.0;
  if (origin != NULL && PyList_CHECK(origin) && pf.PyList_Size(origin) == 2) {
    node->origin_x = pf.PyFloat_AsDouble(pf.PyList_GetItem(origin, 0));
    node->origin_y = pf.PyFloat_AsDouble(pf.PyList_GetItem(origin, 1));
  }

  node->scale_x = 1.0;
  node->scale_y = 1.0;
  if (scale != NULL && PyList_CHECK(scale) && pf.PyList_Size(scale) == 2) {
    node->scale_x = pf.PyFloat_AsDouble(pf.PyList_GetItem(scale, 0));
    node->scale_y = pf.PyFloat_AsDouble(pf.PyList_GetItem(scale, 1));
  }

  if (center == Py_TRUE) node->center = 1;
  else if (center == Py_FALSE || center == NULL) node->center = 0;
  else {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Unknown Value for center parameter");
    return NULL;
  }

  node->slices = slices;
  node->has_slices = slices != 1 ? 1 : 0;

  node->segments = segments;
  node->has_segments = segments != 1 ? 1 : 0;

  if (twist != NULL) {
    if (twist->ob_type == pf.PyFunction_Type) {
      Py_XINCREF(twist);  // TODO there to decref it ?
      node->twist_func = twist;
    } else node->twist = pf.PyFloat_AsDouble(twist);
    node->has_twist = 1;
  } else node->has_twist = 0;
  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

PyObject *python_linear_extrude(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyObject *obj = NULL;
  PyObject *height = NULL;
  int convexity = 1;
  PyObject *origin = NULL;
  PyObject *scale = NULL;
  PyObject *center = NULL;
  int slices = 1;
  int segments = 0;
  PyObject *twist = NULL;
  double fn = NAN, fa = NAN, fs = NAN;

  char *kwlist[] = {"obj",      "height", "convexity", "origin", "scale", "center", "slices",
                    "segments", "twist",  "fn",        "fa",     "fs",    NULL};
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O|OiOOOiiOddd", kwlist, &obj, &height, &convexity,
                                      &origin, &scale, &center, &slices, &segments, &twist, &fn, &fs,
                                      &fs)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error during parsing\n");
    return NULL;
  }

  return linear_extrude_core(obj, height, convexity, origin, scale, center, slices, segments, twist, fn,
                             fa, fs);
}

PyObject *python_oo_linear_extrude(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  PyObject *height = NULL;
  int convexity = 1;
  PyObject *origin = NULL;
  PyObject *scale = NULL;
  PyObject *center = NULL;
  int slices = 1;
  int segments = 0;
  PyObject *twist = NULL;
  double fn = NAN, fa = NAN, fs = NAN;

  char *kwlist[] = {"height",   "convexity", "origin", "scale", "center", "slices",
                    "segments", "twist",     "fn",     "fa",    "fs",     NULL};
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "|OiOOOiiOddd", kwlist, &height, &convexity, &origin,
                                      &scale, &center, &slices, &segments, &twist, &fn, &fs, &fs)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error during parsing\n");
    return NULL;
  }

  return linear_extrude_core(obj, height, convexity, origin, scale, center, slices, segments, twist, fn,
                             fa, fs);
}

PyObject *path_extrude_core(PyObject *obj, PyObject *path, PyObject *xdir, int convexity,
                            PyObject *origin, PyObject *scale, PyObject *twist, PyObject *closed,
                            PyObject *allow_intersect, double fn, double fa, double fs)
{
  DECLARE_INSTANCE
  std::shared_ptr<AbstractNode> child;
  auto node = std::make_shared<PathExtrudeNode>(instance);
  PyTypeObject *type = &PyOpenSCADType;
  node->profile_func = NULL;
  node->twist_func = NULL;
  if (obj->ob_type == pf.PyFunction_Type) {
    Py_XINCREF(obj);  // TODO there to decref it ?
    node->profile_func = obj;
    auto dummy_node = std::make_shared<SquareNode>(instance);
    node->children.push_back(dummy_node);
  } else {
    PyObject *dummydict;
    type = PyOpenSCADObjectType(obj);
    child = PyOpenSCADObjectToNodeMulti(obj, &dummydict);
    if (child == NULL) {
      PyErr_SetString(PyExc_TypeError, "Invalid type for  Object in path_extrude\n");
      return NULL;
    }
    node->children.push_back(child);
  }
  if (path != NULL && PyList_Check(path)) {
    int n = PyList_Size(path);
    for (int i = 0; i < n; i++) {
      PyObject *point = PyList_GetItem(path, i);
      double x, y, z, w = 0;
      if (python_vectorval(point, 3, 4, &x, &y, &z, &w)) {
        PyErr_SetString(PyExc_TypeError, "Cannot parse vector in path_extrude path\n");
        return NULL;
      }
      Vector4d pt3d(x, y, z, w);
      if (i > 0 && node->path[i - 1] == pt3d) continue;  //  prevent double pts
      node->path.push_back(pt3d);
    }
  }
  node->xdir_x = 1;
  node->xdir_y = 0;
  node->xdir_z = 0;
  node->closed = false;
  if (closed == Py_True) node->closed = true;
  if (allow_intersect == Py_True) node->allow_intersect = true;
  if (xdir != NULL) {
    if (python_vectorval(xdir, 3, 3, &(node->xdir_x), &(node->xdir_y), &(node->xdir_z))) {
      PyErr_SetString(PyExc_TypeError, "error in path_extrude xdir parameter\n");
      return NULL;
    }
  }
  if (fabs(node->xdir_x) < 0.001 && fabs(node->xdir_y) < 0.001 && fabs(node->xdir_z) < 0.001) {
    PyErr_SetString(PyExc_TypeError, "error in path_extrude xdir parameter has zero size\n");
    return NULL;
  }

  get_fnas(node->fn, node->fa, node->fs);
  if (fn != -1) node->fn = fn;
  if (fa != -1) node->fa = fa;
  if (fs != -1) node->fs = fs;

  node->convexity = convexity;

  node->origin_x = 0.0;
  node->origin_y = 0.0;
  if (origin != NULL) {
    double dummy;
    if (python_vectorval(origin, 2, 2, &(node->origin_x), &(node->origin_y), &dummy)) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "error in path_extrude origin parameter\n");
      return NULL;
    }
  }

  node->scale_x = 1.0;
  node->scale_y = 1.0;
  if (scale != NULL) {
    double dummy;
    if (python_vectorval(scale, 2, 2, &(node->scale_x), &(node->scale_y), &dummy)) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "error in path_extrude scale parameter\n");
      return NULL;
    }
  }

  if (scale != NULL && PyList_Check(scale) && pf.PyList_Size(scale) == 2) {
    node->scale_x = pf.PyFloat_AsDouble(pf.PyList_GetItem(scale, 0));
    node->scale_y = pf.PyFloat_AsDouble(pf.PyList_GetItem(scale, 1));
  }
  if (twist != NULL) {
    if (twist->ob_type == pf.PyFunction_Type) {
      Py_XINCREF(twist);  // TODO there to decref it ?
      node->twist_func = twist;
    } else node->twist = pf.PyFloat_AsDouble(twist);
    node->has_twist = 1;
  } else node->has_twist = 0;

  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

PyObject *python_path_extrude(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyObject *obj = NULL;
  int convexity = 1;
  PyObject *origin = NULL;
  PyObject *scale = NULL;
  PyObject *path = NULL;
  PyObject *xdir = NULL;
  PyObject *closed = NULL;
  PyObject *allow_intersect = NULL;
  PyObject *twist = NULL;
  double fn = -1, fa = -1, fs = -1;

  char *kwlist[] = {"obj",   "path",   "xdir", "convexity", "origin", "scale",
                    "twist", "closed", "fn",   "fa",        "fs",     NULL};
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "OO!|O!iOOOOOddd", kwlist, &obj, pf.PyList_Type,
                                      &path, pf.PyList_Type, &xdir, &convexity, &origin, &scale, &twist,
                                      &closed, &allow_intersect, &fn, &fs, &fs)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error during parsing\n");
    return NULL;
  }

  return path_extrude_core(obj, path, xdir, convexity, origin, scale, twist, closed, allow_intersect, fn,
                           fa, fs);
}

PyObject *python_concat(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE
  int i;

  auto node = std::make_shared<ConcatNode>(instance);
  PyObject *obj;
  PyObject *obj1;
  PyObject *child_dict = nullptr;
  std::shared_ptr<AbstractNode> child;
  PyTypeObject *type = &PyOpenSCADType;
  // dont do union in any circumstance
  for (i = 0; i < pf.PyTuple_Size(args); i++) {
    obj = pf.PyTuple_GetItem(args, i);
    if (PyObject_IsInstance(obj, reinterpret_cast<PyObject *>(&PyOpenSCADType))) {
      type = PyOpenSCADObjectType(obj);
      node->children.push_back(((PyOpenSCADObject *)obj)->node);
    } else if (PyList_CHECK(obj)) {
      for (int j = 0; j < pf.PyList_Size(obj); j++) {
        obj1 = pf.PyList_GetItem(obj, j);
        if (PyObject_IsInstance(obj1, reinterpret_cast<PyObject *>(&PyOpenSCADType))) {
          type = PyOpenSCADObjectType(obj1);
          node->children.push_back(((PyOpenSCADObject *)obj1)->node);
        } else {
          pf.PyErr_SetString(pf.PyExc_TypeError, "Error during concat. arguments must be solids");
          return nullptr;
        }
      }
    } else {
      pf.PyErr_SetString(pf.PyExc_TypeError, "Error during concat. arguments must be solids");
      return nullptr;
    }
  }

  PyObject *pyresult = PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
  if (child_dict != nullptr) {
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (pf.PyDict_Next(child_dict, &pos, &key, &value)) {
      pf.PyDict_SetItem(((PyOpenSCADObject *)pyresult)->dict, key, value);
    }
  }
  return pyresult;
}

PyObject *python_skin(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE
  int i;

  auto node = std::make_shared<SkinNode>(instance);
  PyTypeObject *type = &PyOpenSCADType;
  PyObject *obj;
  PyObject *child_dict = nullptr;
  PyObject *dummy_dict = nullptr;
  std::shared_ptr<AbstractNode> child;
  if (kwargs != nullptr) {
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (pf.PyDict_Next(kwargs, &pos, &key, &value)) {
      PyObject *value1 = pf.PyUnicode_AsEncodedString(key, "utf-8", "~");
      const char *value_str = PyBytes_AS_STRING(value1);
      double tmp;
      if (value_str == nullptr) {
        pf.PyErr_SetString(pf.PyExc_TypeError, "Unkown parameter name in CSG.");
        return nullptr;
      } else if (strcmp(value_str, "convexity") == 0) {
        python_numberval(value, &tmp);
        node->convexity = (int)tmp;
      } else if (strcmp(value_str, "align_angle") == 0) {
        python_numberval(value, &tmp);
        node->align_angle = tmp;
        node->has_align_angle = true;
      } else if (strcmp(value_str, "segments") == 0) {
        python_numberval(value, &tmp);
        node->has_segments = true;
        node->segments = (int)tmp;
      } else if (strcmp(value_str, "interpolate") == 0) {
        python_numberval(value, &tmp);
        node->has_interpolate = true;
        node->interpolate = tmp;
      } else {
        pf.PyErr_SetString(pf.PyExc_TypeError, "Unkown parameter name in skin.");
        return nullptr;
      }
    }
  }
  for (i = 0; i < pf.PyTuple_Size(args); i++) {
    obj = pf.PyTuple_GetItem(args, i);
    if (i == 0) child = PyOpenSCADObjectToNodeMulti(obj, &child_dict);
    else child = PyOpenSCADObjectToNodeMulti(obj, &dummy_dict);
    if (child != NULL) {
      node->children.push_back(child);
    } else {
      pf.PyErr_SetString(pf.PyExc_TypeError, "Error during skin. arguments must be solids or arrays.");
      return nullptr;
    }
  }

  PyObject *pyresult = PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
  if (child_dict != nullptr) {
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (pf.PyDict_Next(child_dict, &pos, &key, &value)) {
      pf.PyDict_SetItem(((PyOpenSCADObject *)pyresult)->dict, key, value);
    }
  }
  return pyresult;
}

PyObject *python_oo_path_extrude(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  int convexity = 1;
  PyObject *origin = NULL;
  PyObject *scale = NULL;
  PyObject *path = NULL;
  PyObject *xdir = NULL;
  PyObject *closed = NULL;
  PyObject *allow_intersect = NULL;
  PyObject *twist = NULL;
  double fn = -1, fa = -1, fs = -1;

  char *kwlist[] = {"path",   "xdir", "convexity", "origin", "scale", "twist",
                    "closed", "fn",   "fa",        "fs",     NULL};
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O!|O!iOOOOOddd", kwlist, pf.PyList_Type, &path,
                                      pf.PyList_Type, &xdir, &convexity, &origin, &scale, &twist,
                                      &closed, &allow_intersect, &fn, &fs, &fs)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "error during parsing\n");
    return NULL;
  }

  return path_extrude_core(obj, path, xdir, convexity, origin, scale, twist, closed, allow_intersect, fn,
                           fa, fs);
}

PyObject *python_csg_core(std::shared_ptr<CsgOpNode>& node,
                          const std::vector<std::shared_ptr<AbstractNode>>& childs)
{
  PyTypeObject *type = &PyOpenSCADType;
  for (size_t i = 0; i < childs.size(); i++) {
    const auto& child = childs[i];
    if (child.get() == void_node.get()) {
      if (node->type == OpenSCADOperator::DIFFERENCE && i == 0)
        return PyOpenSCADObjectFromNode(type, void_node);
      if (node->type == OpenSCADOperator::INTERSECTION) return PyOpenSCADObjectFromNode(type, void_node);
    } else if (child.get() == full_node.get()) {
      if (node->type == OpenSCADOperator::UNION) return PyOpenSCADObjectFromNode(type, full_node);
      if (node->type == OpenSCADOperator::DIFFERENCE) {
        if (i == 0) return PyOpenSCADObjectFromNode(type, full_node);  // eigentlich negativ
        else return PyOpenSCADObjectFromNode(type, void_node);
      }
    } else node->children.push_back(child);
  }
  return PyOpenSCADObjectFromNode(type, node);
}

PyObject *python_csg_sub(PyObject *self, PyObject *args, PyObject *kwargs, OpenSCADOperator mode)
{
  DECLARE_INSTANCE
  int i;
  auto node = std::make_shared<CsgOpNode>(instance, mode);
  node->r = 0;
  node->fn = 1;
  PyObject *obj;
  std::vector<PyObject *> child_dict;
  std::shared_ptr<AbstractNode> child;
  if (kwargs != nullptr) {
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (pf.PyDict_Next(kwargs, &pos, &key, &value)) {
      PyObject *value1 = pf.PyUnicode_AsEncodedString(key, "utf-8", "~");
      const char *value_str = PyBytes_AS_STRING(value1);
      if (value_str == nullptr) {
        pf.PyErr_SetString(pf.PyExc_TypeError, "Unkown parameter name in CSG.");
        return nullptr;
      } else if (strcmp(value_str, "r") == 0) {
        python_numberval(value, &(node->r));
      } else if (strcmp(value_str, "fn") == 0) {
        double fn;
        python_numberval(value, &fn);
        node->fn = (int)fn;
      } else {
        pf.PyErr_SetString(pf.PyExc_TypeError, "Unkown parameter name in CSG.");
        return nullptr;
      }
    }
  }
  std::vector<std::shared_ptr<AbstractNode>> child_solid;
  for (i = 0; i < pf.PyTuple_Size(args); i++) {
    obj = pf.PyTuple_GetItem(args, i);
    PyObject *dict = nullptr;
    child = PyOpenSCADObjectToNodeMulti(obj, &dict);
    child_dict.push_back(dict);
    if (child != NULL) {
      if (child.get() == void_node.get() && mode == OpenSCADOperator::UNION) {
      } else if (child.get() == void_node.get() && i > 0 && mode == OpenSCADOperator::DIFFERENCE) {
      } else if (child.get() == full_node.get() && mode == OpenSCADOperator::INTERSECTION) {
      } else {
        node->children.push_back(child);
      }
    } else {
      switch (mode) {
      case OpenSCADOperator::UNION:
        pf.PyErr_SetString(pf.PyExc_TypeError,
                           "Error during parsing union. arguments must be solids or arrays.");
        return nullptr;
        break;
      case OpenSCADOperator::DIFFERENCE:
        pf.PyErr_SetString(pf.PyExc_TypeError,
                           "Error during parsing difference. arguments must be solids or arrays.");
        return nullptr;
        break;
      case OpenSCADOperator::INTERSECTION:
        pf.PyErr_SetString(pf.PyExc_TypeError,
                           "Error during parsing intersection. arguments must be solids or arrays.");
        return nullptr;
        break;
      case OpenSCADOperator::MINKOWSKI: break;
      case OpenSCADOperator::HULL:      break;
      case OpenSCADOperator::FILL:      break;
      case OpenSCADOperator::RESIZE:    break;
      case OpenSCADOperator::OFFSET:    break;
      }
      return NULL;
    }
  }
  PyObject *pyresult = python_csg_core(node, child_solid);

  for (int i = child_dict.size() - 1; i >= 0; i--)  // merge from back  to give 1st child most priority
  {
    auto& dict = child_dict[i];
    if (dict == nullptr) continue;
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (pf.PyDict_Next(dict, &pos, &key, &value)) {
      PyObject *value1 = pf.PyUnicode_AsEncodedString(key, "utf-8", "~");
      const char *value_str = PyBytes_AS_STRING(value1);
      pf.PyDict_SetItem(((PyOpenSCADObject *)pyresult)->dict, key, value);
    }
  }
  return pyresult;
}

PyObject *python_union(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_csg_sub(self, args, kwargs, OpenSCADOperator::UNION);
}

PyObject *python_difference(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_csg_sub(self, args, kwargs, OpenSCADOperator::DIFFERENCE);
}

PyObject *python_intersection(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_csg_sub(self, args, kwargs, OpenSCADOperator::INTERSECTION);
}

PyObject *python_oo_csg_sub(PyObject *self, PyObject *args, PyObject *kwargs, OpenSCADOperator mode)
{
  DECLARE_INSTANCE
  int i;

  auto node = std::make_shared<CsgOpNode>(instance, mode);
  node->r = 0;
  node->fn = 1;

  PyObject *obj;
  std::vector<PyObject *> child_dict;
  std::shared_ptr<AbstractNode> child;
  PyObject *dict;

  dict = nullptr;
  child = PyOpenSCADObjectToNodeMulti(self, &dict);
  if (child != NULL) {
    node->children.push_back(child);
    child_dict.push_back(dict);
  }

  if (kwargs != nullptr) {
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (pf.PyDict_Next(kwargs, &pos, &key, &value)) {
      PyObject *value1 = pf.PyUnicode_AsEncodedString(key, "utf-8", "~");
      const char *value_str = PyBytes_AS_STRING(value1);
      if (value_str == nullptr) {
        pf.PyErr_SetString(pf.PyExc_TypeError, "Unkown parameter name in CSG.");
        return nullptr;
      } else if (strcmp(value_str, "r") == 0) {
        python_numberval(value, &(node->r));
      } else if (strcmp(value_str, "fn") == 0) {
        double fn;
        python_numberval(value, &fn);
        node->fn = (int)fn;
      } else {
        pf.PyErr_SetString(pf.PyExc_TypeError, "Unkown parameter name in CSG.");
        return nullptr;
      }
    }
  }
  std::vector<std::shared_ptr<AbstractNode>> child_solid;
  for (i = 0; i < pf.PyTuple_Size(args); i++) {
    obj = pf.PyTuple_GetItem(args, i);
    child = PyOpenSCADObjectToNodeMulti(obj, &dict);
    if (child != NULL) {
      if (child.get() == void_node.get() && mode == OpenSCADOperator::UNION) {
      } else if (child.get() == void_node.get() && i > 0 && mode == OpenSCADOperator::DIFFERENCE) {
      } else if (child.get() == full_node.get() && mode == OpenSCADOperator::INTERSECTION) {
      } else {
        node->children.push_back(child);
      }
    } else {
      switch (mode) {
      case OpenSCADOperator::UNION:
        pf.PyErr_SetString(pf.PyExc_TypeError,
                           "Error during parsing union. arguments must be solids or arrays.");
        break;
      case OpenSCADOperator::DIFFERENCE:
        pf.PyErr_SetString(pf.PyExc_TypeError,
                           "Error during parsing difference. arguments must be solids or arrays.");
        break;
      case OpenSCADOperator::INTERSECTION:
        pf.PyErr_SetString(pf.PyExc_TypeError,
                           "Error during parsing intersection. arguments must be solids or arrays.");
        break;
      case OpenSCADOperator::MINKOWSKI: break;
      case OpenSCADOperator::HULL:      break;
      case OpenSCADOperator::FILL:      break;
      case OpenSCADOperator::RESIZE:    break;
      case OpenSCADOperator::OFFSET:    break;
      }
      return NULL;
    }
  }

  PyObject *pyresult = python_csg_core(node, child_solid);
  for (int i = child_dict.size() - 1; i >= 0; i--)  // merge from back  to give 1st child most priority
  {
    auto& dict = child_dict[i];
    if (dict == nullptr) continue;
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (pf.PyDict_Next(dict, &pos, &key, &value)) {
      pf.PyDict_SetItem(((PyOpenSCADObject *)pyresult)->dict, key, value);
    }
  }
  return pyresult;
}

PyObject *python_oo_union(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_oo_csg_sub(self, args, kwargs, OpenSCADOperator::UNION);
}

PyObject *python_oo_difference(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_oo_csg_sub(self, args, kwargs, OpenSCADOperator::DIFFERENCE);
}

PyObject *python_oo_intersection(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_oo_csg_sub(self, args, kwargs, OpenSCADOperator::INTERSECTION);
}

PyObject *python_nb_sub(PyObject *arg1, PyObject *arg2, OpenSCADOperator mode)
{
  DECLARE_INSTANCE
  std::vector<std::shared_ptr<AbstractNode>> child;
  std::vector<PyObject *> child_dict;

  if (arg1 == Py_NONE && mode == OpenSCADOperator::UNION) return arg2;
  if (arg2 == Py_NONE && mode == OpenSCADOperator::UNION) return arg1;
  if (arg2 == Py_NONE && mode == OpenSCADOperator::DIFFERENCE) return arg1;

  for (int i = 0; i < 2; i++) {
    PyObject *dict;
    dict = nullptr;
    auto solid = PyOpenSCADObjectToNodeMulti(i == 1 ? arg2 : arg1, &dict);
    child_dict.push_back(dict);
    if (solid != nullptr) child.push_back(solid);
    else {
      PyErr_SetString(PyExc_TypeError, "invalid argument left to operator");
      return NULL;
    }
  }
  auto node = std::make_shared<CsgOpNode>(instance, mode);
  PyObject *pyresult = python_csg_core(node, child);

  python_retrieve_pyname(node);
  for (int i = 1; i >= 0; i--) {
    if (child_dict[i] != nullptr) {
      std::string name = child[i]->getPyName();
      PyObject *key, *value;
      Py_ssize_t pos = 0;
      while (pf.PyDict_Next(child_dict[i], &pos, &key, &value)) {
        if (name.size() > 0) {
          PyObject *key1 = pf.PyUnicode_AsEncodedString(key, "utf-8", "~");
          const char *key_str = PyBytes_AS_STRING(key1);
          std::string handle_name = name + "_" + key_str;
          PyObject *key_mod =
            pf.PyUnicode_FromStringAndSize(handle_name.c_str(), strlen(handle_name.c_str()));
          pf.PyDict_SetItem(((PyOpenSCADObject *)pyresult)->dict, key_mod, value);
        } else if (i == 0) {
          pf.PyDict_SetItem(((PyOpenSCADObject *)pyresult)->dict, key, value);
        }
      }
    }
  }
  return pyresult;
}

PyObject *python_nb_sub_vec3(PyObject *arg1, PyObject *arg2,
                             int mode)  // 0: translate, 1: scale, 2: translateneg, 3=translate-exp
{
  DECLARE_INSTANCE
  std::shared_ptr<AbstractNode> child;
  PyObject *child_dict;

  PyTypeObject *type = PyOpenSCADObjectType(arg1);
  child = PyOpenSCADObjectToNodeMulti(arg1, &child_dict);
  if (arg2 == nullptr) return PyOpenSCADObjectFromNode(&PyOpenSCADType, child);
  std::vector<Vector3d> vecs;
  int dragflags = 0;
  if (mode == 3) {
    if (!PyList_CHECK(arg2)) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "explode arg must be a list");
      return NULL;
    }
    int n = pf.PyList_Size(arg2);
    if (pf.PyList_Size(arg2) > 3) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "explode arg list can have maximal 3 directions");
      return NULL;
    }
    double dmy;
    std::vector<float> vals[3];
    for (int i = 0; i < 3; i++) vals[i].push_back(0.0);
    for (int i = 0; i < n; i++) {
      vals[i].clear();
      auto *item = pf.PyList_GetItem(arg2, i);  // TODO fix here
      if (!python_numberval(item, &dmy, &dragflags, 1 << i)) vals[i].push_back(dmy);
      else if (PyList_CHECK(item)) {
        int m = pf.PyList_Size(item);
        for (int j = 0; j < m; j++) {
          auto *item1 = pf.PyList_GetItem(item, j);
          if (!python_numberval(item1, &dmy)) vals[i].push_back(dmy);
        }
      } else {
        pf.PyErr_SetString(pf.PyExc_TypeError, "Unknown explode spec");
        return NULL;
      }
    }
    for (auto z : vals[2])
      for (auto y : vals[1])
        for (auto x : vals[0]) vecs.push_back(Vector3d(x, y, z));
  } else vecs = python_vectors(arg2, 2, 3, &dragflags);

  if (mode == 0 && vecs.size() == 1) {  // translate on numbers
    PyObject *mat = python_number_trans(arg1, vecs[0], 4);
    if (mat != nullptr) return mat;
  }

  if (vecs.size() > 0) {
    if (child == NULL) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "invalid argument left to operator");
      return NULL;
    }
    std::vector<std::shared_ptr<TransformNode>> nodes;
    for (size_t j = 0; j < vecs.size(); j++) {
      std::shared_ptr<TransformNode> node;
      switch (mode) {
      case 0:
      case 3:
        node = std::make_shared<TransformNode>(instance, "translate");
        node->matrix.translate(vecs[j]);
        break;
      case 1:
        node = std::make_shared<TransformNode>(instance, "scale");
        node->matrix.scale(vecs[j]);
        break;
      case 2:
        node = std::make_shared<TransformNode>(instance, "translate");
        node->matrix.translate(-vecs[j]);
        break;
      }
      node->children.push_back(child);
      nodes.push_back(node);
    }
    if (nodes.size() == 1) {
      nodes[0]->dragflags = dragflags;
      PyObject *pyresult = PyOpenSCADObjectFromNode(&PyOpenSCADType, nodes[0]);
      if (child_dict != nullptr) {
        PyObject *key, *value;
        Py_ssize_t pos = 0;
        while (pf.PyDict_Next(child_dict, &pos, &key, &value)) {
          PyObject *value1 = python_number_trans(value, vecs[0], 4);
          if (value1 != nullptr) pf.PyDict_SetItem(((PyOpenSCADObject *)pyresult)->dict, key, value1);
          else pf.PyDict_SetItem(((PyOpenSCADObject *)pyresult)->dict, key, value);
        }
      }
      return pyresult;
    } else {
      auto node = std::make_shared<CsgOpNode>(instance, OpenSCADOperator::UNION);
      DECLARE_INSTANCE
      for (auto x : nodes) node->children.push_back(x->clone());
      return PyOpenSCADObjectFromNode(type, node);
    }
  }
  pf.PyErr_SetString(pf.PyExc_TypeError, "invalid argument right to operator");
  return NULL;
}

PyObject *python_nb_add(PyObject *arg1, PyObject *arg2)
{
  return python_nb_sub_vec3(arg1, arg2, 0);
}  // translate

PyObject *python_nb_xor(PyObject *arg1, PyObject *arg2)
{
  PyObject *dummy_dict;
  if (PyObject_IsInstance(arg2, reinterpret_cast<PyObject *>(&PyOpenSCADType))) {
    auto node1 = PyOpenSCADObjectToNode(arg1, &dummy_dict);
    auto node2 = PyOpenSCADObjectToNode(arg2, &dummy_dict);
    if (node1 == nullptr || node2 == nullptr) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing hull. arguments must be solids.");
      return nullptr;
    }
    DECLARE_INSTANCE
    std::shared_ptr<AbstractNode> child;
    auto node = std::make_shared<CgalAdvNode>(instance, CgalAdvType::HULL);
    node->children.push_back(node1);
    node->children.push_back(node2);
    return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
  }
  return python_nb_sub_vec3(arg1, arg2, 3);
}

PyObject *python_nb_remainder(PyObject *arg1, PyObject *arg2)
{
  if (PyObject_IsInstance(arg2, reinterpret_cast<PyObject *>(&PyOpenSCADType))) {
    PyObject *dummy_dict;
    auto node1 = PyOpenSCADObjectToNode(arg1, &dummy_dict);
    auto node2 = PyOpenSCADObjectToNode(arg2, &dummy_dict);
    if (node1 == nullptr || node2 == nullptr) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing hull. arguments must be solids.");
      return nullptr;
    }
    DECLARE_INSTANCE
    std::shared_ptr<AbstractNode> child;
    auto node = std::make_shared<CgalAdvNode>(instance, CgalAdvType::MINKOWSKI);
    node->children.push_back(node1);
    node->children.push_back(node2);
    return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
  }
  Vector3d vec3(0, 0, 0);
  if (!python_vectorval(arg2, 1, 3, &(vec3[0]), &(vec3[1]), &(vec3[2]), nullptr)) {
    return python_rotate_sub(arg1, vec3, NAN, nullptr, 0);
  }

  pf.PyErr_SetString(pf.PyExc_TypeError, "Unknown types for % oprator");
  return nullptr;
}

PyObject *python_nb_mul(PyObject *arg1, PyObject *arg2)
{
  return python_nb_sub_vec3(arg1, arg2, 1);
}  // scale
PyObject *python_nb_or(PyObject *arg1, PyObject *arg2)
{
  return python_nb_sub(arg1, arg2, OpenSCADOperator::UNION);
}
PyObject *python_nb_subtract(PyObject *arg1, PyObject *arg2)
{
  double dmy;
  if (PyList_CHECK(arg2) && pf.PyList_Size(arg2) > 0) {
    PyObject *sub = pf.PyList_GetItem(arg2, 0);
    if (!python_numberval(sub, &dmy) || PyList_CHECK(sub)) {
      return python_nb_sub_vec3(arg1, arg2, 2);
    }
  }
  return python_nb_sub(arg1, arg2, OpenSCADOperator::DIFFERENCE);  // if its solid
}
PyObject *python_nb_and(PyObject *arg1, PyObject *arg2)
{
  return python_nb_sub(arg1, arg2, OpenSCADOperator::INTERSECTION);
}

PyObject *python_nb_matmult(PyObject *arg1, PyObject *arg2)
{
  return python_multmatrix_sub(arg1, arg2, 0);
}

PyObject *python_csg_adv_sub(PyObject *self, PyObject *args, PyObject *kwargs, CgalAdvType mode)
{
  DECLARE_INSTANCE
  std::shared_ptr<AbstractNode> child;
  PyTypeObject *type = &PyOpenSCADType;
  int i;
  PyObject *dummydict;

  auto node = std::make_shared<CgalAdvNode>(instance, mode);
  PyObject *obj;
  for (i = 0; i < pf.PyTuple_Size(args); i++) {
    obj = pf.PyTuple_GetItem(args, i);
    child = PyOpenSCADObjectToNodeMulti(obj, &dummydict);
    if (child != NULL) {
      node->children.push_back(child);
    } else {
      switch (mode) {
      case CgalAdvType::HULL:
        pf.PyErr_SetString(pf.PyExc_TypeError,
                           "Error during parsing hull. arguments must be solids or arrays.");
        break;
      case CgalAdvType::FILL:
        pf.PyErr_SetString(pf.PyExc_TypeError,
                           "Error during parsing fill. arguments must be solids or arrays.");
        break;
      case CgalAdvType::RESIZE:    break;
      case CgalAdvType::MINKOWSKI: break;
      }
      return NULL;
    }
  }

  return PyOpenSCADObjectFromNode(type, node);
}

PyObject *python_minkowski(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE
  std::shared_ptr<AbstractNode> child;
  int convexity = 2;

  auto node = std::make_shared<CgalAdvNode>(instance, CgalAdvType::MINKOWSKI);
  char *kwlist[] = {"obj1", "obj2", "convexity", NULL};
  PyObject *obj1, *obj2;
  PyObject *dummydict;

  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "OO|i", kwlist, &obj1, &obj2, &convexity)) {
    pf.PyErr_SetString(pf.PyExc_TypeError,
                       "Error during parsing minkowski(object1, object2[, convexity])");
    return NULL;
  }
  PyTypeObject *type = PyOpenSCADObjectType(obj1);
  child = PyOpenSCADObjectToNodeMulti(obj1, &dummydict);
  node->children.push_back(child);

  child = PyOpenSCADObjectToNodeMulti(obj2, &dummydict);
  node->children.push_back(child);

  node->convexity = convexity;

  return PyOpenSCADObjectFromNode(type, node);
}

PyObject *python_hull(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_csg_adv_sub(self, args, kwargs, CgalAdvType::HULL);
}

PyObject *python_fill(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_csg_adv_sub(self, args, kwargs, CgalAdvType::FILL);
}

PyObject *python_resize_core(PyObject *obj, PyObject *newsize, PyObject *autosize, int convexity)
{
  DECLARE_INSTANCE
  std::shared_ptr<AbstractNode> child;

  auto node = std::make_shared<CgalAdvNode>(instance, CgalAdvType::RESIZE);
  PyObject *dummydict;
  PyTypeObject *type = PyOpenSCADObjectType(obj);
  child = PyOpenSCADObjectToNodeMulti(obj, &dummydict);
  if (child == NULL) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid type for Object in resize");
    return NULL;
  }

  if (newsize != NULL) {
    double x, y, z;
    if (python_vectorval(newsize, 3, 3, &x, &y, &z)) {
      pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid resize dimensions");
      return NULL;
    }
    node->newsize[0] = x;
    node->newsize[1] = y;
    node->newsize[2] = z;
  }

  /* TODO what is that ?
     const auto& autosize = parameters["auto"];
     node->autosize << false, false, false;
     if (autosize.type() == Value::Type::VECTOR) {
     const auto& va = autosize.toVector();
     if (va.size() >= 1) node->autosize[0] = va[0].toBool();
     if (va.size() >= 2) node->autosize[1] = va[1].toBool();
     if (va.size() >= 3) node->autosize[2] = va[2].toBool();
     } else if (autosize.type() == Value::Type::BOOL) {
     node->autosize << autosize.toBool(), autosize.toBool(), autosize.toBool();
     }
   */

  node->children.push_back(child);
  node->convexity = convexity;

  return PyOpenSCADObjectFromNode(type, node);
}

PyObject *python_resize(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"obj", "newsize", "auto", "convexity", NULL};
  PyObject *obj;
  PyObject *newsize = NULL;
  PyObject *autosize = NULL;
  int convexity = 2;

  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O|O!O!i", kwlist, &obj, pf.PyList_Type, &newsize,
                                      pf.PyList_Type, &autosize, &convexity)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing resize(object,vec3)");
    return NULL;
  }
  return python_resize_core(obj, newsize, autosize, convexity);
}

PyObject *python_oo_resize(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"newsize", "auto", "convexity", NULL};
  PyObject *newsize = NULL;
  PyObject *autosize = NULL;
  int convexity = 2;

  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "|O!O!i", kwlist, pf.PyList_Type, &newsize,
                                      pf.PyList_Type, &autosize, &convexity)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing resize(object,vec3)");
    return NULL;
  }
  return python_resize_core(obj, newsize, autosize, convexity);
}

#if defined(ENABLE_EXPERIMENTAL) && defined(ENABLE_CGAL)
PyObject *python_roof_core(PyObject *obj, const char *method, int convexity, double fn, double fa,
                           double fs)
{
  DECLARE_INSTANCE
  std::shared_ptr<AbstractNode> child;
  auto node = std::make_shared<RoofNode>(instance);
  PyObject *dummydict;
  PyTypeObject *type = PyOpenSCADObjectType(obj);
  child = PyOpenSCADObjectToNodeMulti(obj, &dummydict);
  if (child == NULL) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid type for Object in roof");
    return NULL;
  }

  get_fnas(node->fn, node->fa, node->fs);
  if (!isnan(fn)) node->fn = fn;
  if (!isnan(fa)) node->fa = fa;
  if (!isnan(fs)) node->fs = fs;

  node->fa = std::max(node->fa, 0.01);
  node->fs = std::max(node->fs, 0.01);
  if (node->fn > 0) {
    node->fa = 360.0 / node->fn;
    node->fs = 0.0;
  }

  if (method == NULL) {
    node->method = "voronoi";
  } else {
    node->method = method;
    // method can only be one of...
    if (node->method != "voronoi" && node->method != "straight") {
      //      LOG(message_group::Warning, inst->location(), parameters.documentRoot(),
      //          "Unknown roof method '" + node->method + "'. Using 'voronoi'.");
      node->method = "voronoi";
    }
  }

  double tmp_convexity = convexity;
  node->convexity = static_cast<int>(tmp_convexity);
  if (node->convexity <= 0) node->convexity = 1;

  node->children.push_back(child);
  return PyOpenSCADObjectFromNode(type, node);
}

PyObject *python_roof(PyObject *self, PyObject *args, PyObject *kwargs)
{
  double fn = NAN, fa = NAN, fs = NAN;
  char *kwlist[] = {"obj", "method", "convexity", "fn", "fa", "fs", NULL};
  PyObject *obj = NULL;
  const char *method = NULL;
  int convexity = 2;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O|sdddd", kwlist, &obj, &method, convexity, &fn,
                                      &fa, &fs)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing roof(object)");
    return NULL;
  }
  return python_roof_core(obj, method, convexity, fn, fa, fs);
}

PyObject *python_oo_roof(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  double fn = NAN, fa = NAN, fs = NAN;
  char *kwlist[] = {"method", "convexity", "fn", "fa", "fs", NULL};
  const char *method = NULL;
  int convexity = 2;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "|sdddd", kwlist, &method, convexity, &fn, &fa,
                                      &fs)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing roof(object)");
    return NULL;
  }
  return python_roof_core(obj, method, convexity, fn, fa, fs);
}
#endif

PyObject *python_render_core(PyObject *obj, int convexity)
{
  DECLARE_INSTANCE
  auto node = std::make_shared<RenderNode>(instance);

  PyObject *dummydict;
  PyTypeObject *type = PyOpenSCADObjectType(obj);
  std::shared_ptr<AbstractNode> child = PyOpenSCADObjectToNode(obj, &dummydict);
  node->convexity = convexity;
  node->children.push_back(child);
  return PyOpenSCADObjectFromNode(type, node);
}

PyObject *python_render(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"obj", "convexity", NULL};
  PyObject *obj = NULL;
  long convexity = 2;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O!|i", kwlist, &PyOpenSCADType, &obj, &convexity)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing render(object)");
    return NULL;
  }
  return python_render_core(obj, convexity);
}

PyObject *python_oo_render(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"convexity", NULL};
  long convexity = 2;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "|i", kwlist, &convexity)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing render(object)");
    return NULL;
  }
  return python_render_core(obj, convexity);
}

PyObject *python_surface_core(const char *file, PyObject *center, PyObject *invert, PyObject *color,
                              int convexity)
{
  DECLARE_INSTANCE
  std::shared_ptr<AbstractNode> child;

  auto node = std::make_shared<SurfaceNode>(instance);

  std::string fileval = file == NULL ? "" : file;

  std::string filename = lookup_file(fileval, python_scriptpath.parent_path().u8string(),
                                     instance->location().filePath().parent_path().string());
  node->filename = filename;
  handle_dep(fs::path(filename).generic_string());

  if (center == Py_TRUE) node->center = 1;
  else if (center == Py_FALSE || center == NULL) node->center = 0;
  else {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Unknown Value for center parameter");
    return NULL;
  }
  node->convexity = 2;
  if (invert == Py_TRUE) node->invert = 1;
  else if (center == Py_FALSE || center == NULL) node->center = 0;
  else {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Unknown Value for invert parameter");
    return NULL;
  }

  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

PyObject *python_surface(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"file", "center", "convexity", "invert", "color", NULL};
  const char *file = NULL;
  PyObject *center = NULL;
  PyObject *invert = NULL;
  PyObject *color = NULL;
  long convexity = 2;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "s|OlO", kwlist, &file, &center, &convexity)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing surface(object)");
    return NULL;
  }

  return python_surface_core(file, center, invert, color, convexity);
}

int sheetCalcIndInt(PyObject *func, double i, double j, Vector3d& pos)
{
  PyObject *args = pf.PyTuple_Pack(2, pf.PyFloat_FromDouble(i), pf.PyFloat_FromDouble(j));
  PyObject *pos_p = PyObject_CallObject(func, args);
  if (pos_p == nullptr) {
    std::string errorstr;
    python_catch_error(errorstr);
    pf.PyErr_SetString(pf.PyExc_TypeError, errorstr.c_str());
    LOG(message_group::Error, errorstr.c_str());
    return 1;
  }
  return python_vectorval(pos_p, 3, 3, &pos[0], &pos[1], &pos[2], nullptr, nullptr);
}

int sheetCalcInd(PolySetBuilder& builder, std::vector<Vector3d>& vertices, std::vector<double>& istore,
                 std::vector<double>& jstore, PyObject *func, double i, double j)
{
  std::string errorstr;
  Vector3d pos;
  if (sheetCalcIndInt(func, i, j, pos)) return -1;
  //  printf("pos %g/%g/%g\n", pos[0], pos[1], pos[2]);
  unsigned int ind = builder.vertexIndex(pos);
  if (ind == vertices.size()) {
    vertices.push_back(pos);
    istore.push_back(i);
    jstore.push_back(j);
  }
  return ind;
}

std::unique_ptr<const Geometry> sheetCreateFuncGeometry(void *funcptr, double imin, double imax,
                                                        double jmin, double jmax, double fs, bool ispan,
                                                        bool jspan)
{
  PyObject *func = (PyObject *)funcptr;
  std::unordered_map<SphereEdgeDb, int, boost::hash<SphereEdgeDb>> edges;

  PolySetBuilder builder;
  std::vector<Vector3d> vertices;
  std::vector<double> istore;
  std::vector<double> jstore;
  //
  int ind11, ind21, ind12, ind22;

  ind11 = sheetCalcInd(builder, vertices, istore, jstore, func, imin, jmin);
  ind12 = sheetCalcInd(builder, vertices, istore, jstore, func, imin, jmax);
  ind21 = sheetCalcInd(builder, vertices, istore, jstore, func, imax, jmin);
  ind22 = sheetCalcInd(builder, vertices, istore, jstore, func, imax, jmax);
  if (ind11 < 0 || ind12 < 0 || ind21 < 0 || ind22 < 0) return builder.build();

  std::vector<IndexedTriangle> triangles;
  std::vector<IndexedTriangle> tri_new;
  tri_new.push_back(IndexedTriangle(ind11, ind21, ind22));
  tri_new.push_back(IndexedTriangle(ind11, ind22, ind12));
  int round = 0;
  unsigned int i1, i2, imid;
  Vector3d p1, p2, p3, pmin, pmax, pmid, pmid_test, dir1, dir2;
  double dist, ang, ang_test;
  do {
    triangles = tri_new;
    if (round >= 15) break;  // emergency stop for non-continous models
    tri_new.clear();
    std::vector<int> midinds;
    for (const IndexedTriangle& tri : triangles) {
      int zeroang = 0;
      unsigned int midind = -1;
      for (int i = 0; i < 3; i++) {
        i1 = tri[i];
        i2 = tri[(i + 1) % 3];
        SphereEdgeDb edge(i1, i2);
        if (edges.count(edge) > 0) continue;
        dist = (vertices[i1] - vertices[i2]).norm();
        if (dist < fs && round > 2) continue;
        p1 = vertices[i1];
        p2 = vertices[i2];
        double i_mid = (istore[i1] + istore[i2]) / 2.0;
        double j_mid = (jstore[i1] + jstore[i2]) / 2.0;
        if (sheetCalcIndInt(func, i_mid, j_mid, pmid)) return builder.build();
        dir1 = (pmid - p1).normalized();
        dir2 = (p2 - pmid).normalized();
        ang = acos(dir1.dot(dir2));
        imid = builder.vertexIndex(pmid);
        if (imid == vertices.size()) {
          vertices.push_back(pmid);
          istore.push_back(i_mid);
          jstore.push_back(j_mid);
        }
        if (ang < 0.001 && round > 2) {
          zeroang++;
          continue;
        }
        edges[edge] = imid;
      }
      if (zeroang == 3) {
        p1 = vertices[tri[0]];
        p2 = vertices[tri[1]];
        p3 = vertices[tri[2]];
        double i_mid = (istore[tri[0]] + istore[tri[1]] + istore[tri[2]]) / 3.0;
        double j_mid = (jstore[tri[0]] + jstore[tri[1]] + jstore[tri[2]]) / 3.0;

        if (sheetCalcIndInt(func, i_mid, j_mid, pmid)) return builder.build();
        Vector4d norm = calcTriangleNormal(vertices, {tri[0], tri[1], tri[2]});
        if (fabs(pmid.dot(norm.head<3>()) - norm[3]) > 1e-3) {
          midind = builder.vertexIndex(pmid);
          if (midind == vertices.size()) {
            vertices.push_back(pmid);
            istore.push_back(i_mid);
            jstore.push_back(j_mid);
          }
        }
      }
      midinds.push_back(midind);
    }
    // create new triangles from split edges
    int ind = 0;
    for (const IndexedTriangle& tri : triangles) {
      int splitind[3];
      for (int i = 0; i < 3; i++) {
        SphereEdgeDb e(tri[i], tri[(i + 1) % 3]);
        splitind[i] = edges.count(e) > 0 ? edges[e] : -1;
      }

      if (midinds[ind] != -1) {
        for (int i = 0; i < 3; i++) {
          if (splitind[i] == -1) {
            tri_new.push_back(IndexedTriangle(tri[i], tri[(i + 1) % 3], midinds[ind]));
          } else {
            tri_new.push_back(IndexedTriangle(tri[i], splitind[i], midinds[ind]));
            tri_new.push_back(IndexedTriangle(splitind[i], tri[(i + 1) % 3], midinds[ind]));
          }
        }
        ind++;
        continue;
      }

      int bucket =
        ((splitind[0] != -1) ? 1 : 0) | ((splitind[1] != -1) ? 2 : 0) | ((splitind[2] != -1) ? 4 : 0);
      switch (bucket) {
      case 0: tri_new.push_back(IndexedTriangle(tri[0], tri[1], tri[2])); break;
      case 1:
        tri_new.push_back(IndexedTriangle(tri[0], splitind[0], tri[2]));
        tri_new.push_back(IndexedTriangle(tri[2], splitind[0], tri[1]));
        break;
      case 2:
        tri_new.push_back(IndexedTriangle(tri[1], splitind[1], tri[0]));
        tri_new.push_back(IndexedTriangle(tri[0], splitind[1], tri[2]));
        break;
      case 3:
        tri_new.push_back(IndexedTriangle(tri[0], splitind[0], tri[2]));
        tri_new.push_back(IndexedTriangle(splitind[0], splitind[1], tri[2]));
        tri_new.push_back(IndexedTriangle(splitind[0], tri[1], splitind[1]));
        break;
      case 4:
        tri_new.push_back(IndexedTriangle(tri[2], splitind[2], tri[1]));
        tri_new.push_back(IndexedTriangle(tri[1], splitind[2], tri[0]));
        break;
      case 5:
        tri_new.push_back(IndexedTriangle(tri[0], splitind[0], splitind[2]));
        tri_new.push_back(IndexedTriangle(splitind[0], tri[2], splitind[2]));
        tri_new.push_back(IndexedTriangle(splitind[0], tri[1], tri[2]));
        break;
      case 6:
        tri_new.push_back(IndexedTriangle(tri[0], tri[1], splitind[2]));
        tri_new.push_back(IndexedTriangle(tri[1], splitind[1], splitind[2]));
        tri_new.push_back(IndexedTriangle(splitind[1], tri[2], splitind[2]));
        break;
      case 7:
        tri_new.push_back(IndexedTriangle(splitind[2], tri[0], splitind[0]));
        tri_new.push_back(IndexedTriangle(splitind[0], tri[1], splitind[1]));
        tri_new.push_back(IndexedTriangle(splitind[1], tri[2], splitind[2]));
        tri_new.push_back(IndexedTriangle(splitind[0], splitind[1], splitind[2]));
        break;
      }
      ind++;
    }

    round++;
  } while (tri_new.size() != triangles.size());

  // create point mapping
  std::vector<int> mapping;
  for (int i = 0; i < istore.size(); i++) mapping.push_back(i);

  // now sewing ispan together
  if (ispan) {
    // collect low and high points
    intList minList, maxList;
    for (int i = 0; i < istore.size(); i++) {
      if (istore[i] == imin) minList.push_back(i);
      if (istore[i] == imax) maxList.push_back(i);
    }

    // now sort the list for acending jstore
    std::sort(minList.begin(), minList.end(),
              [jstore](const int& a, const int& b) { return jstore[a] < jstore[b]; });
    std::sort(maxList.begin(), maxList.end(),
              [jstore](const int& a, const int& b) { return jstore[a] < jstore[b]; });

    // pair up min with max
    int minptr = 0, minlen = minList.size();
    int maxptr = 0, maxlen = minList.size();
    while (minptr < minlen && maxptr < maxlen) {  // if one hit the end, no matches meaningful
      double diff = jstore[minList[minptr]] - jstore[maxList[maxptr]];
      if (fabs(diff) < 1e-6) {
        int tgt = maxList[maxptr++];
        int src = mapping[minList[minptr++]];
        vertices[src] = (vertices[tgt] + vertices[src]) / 2.0;
        mapping[tgt] = src;
      } else if (diff > 0) maxptr++;
      else minptr++;
    }
  }

  if (jspan) {
    // collect low and high points
    intList minList, maxList;
    for (int i = 0; i < jstore.size(); i++) {
      if (jstore[i] == jmin) minList.push_back(i);
      if (jstore[i] == jmax) maxList.push_back(i);
    }

    // now sort the list for acending jstore
    std::sort(minList.begin(), minList.end(),
              [istore](const int& a, const int& b) { return istore[a] < istore[b]; });
    std::sort(maxList.begin(), maxList.end(),
              [istore](const int& a, const int& b) { return istore[a] < istore[b]; });

    // pair up min with max
    int minptr = 0, minlen = minList.size();
    int maxptr = 0, maxlen = minList.size();
    while (minptr < minlen && maxptr < maxlen) {  // if one hit the end, no matches meaningful
      double diff = istore[minList[minptr]] - istore[maxList[maxptr]];
      if (fabs(diff) < 1e-6) {
        int tgt = maxList[maxptr++];
        int src = mapping[minList[minptr++]];
        vertices[src] = (vertices[tgt] + vertices[src]) / 2.0;
        mapping[tgt] = src;
      } else if (diff > 0) maxptr++;
      else minptr++;
    }
  }

  for (const IndexedTriangle& tri : tri_new) {
    builder.appendPolygon({mapping[tri[0]], mapping[tri[1]], mapping[tri[2]]});
  }
  return builder.build();
}

PyObject *python_sheet_core(PyObject *func, double imin, double imax, double jmin, double jmax,
                            double fs, PyObject *ispan, PyObject *jspan)
{
  DECLARE_INSTANCE
  auto node = std::make_shared<SheetNode>(instance);
  // TODO check type of func
  node->func = (void *)func;
  node->imin = imin;
  node->imax = imax;
  node->jmin = jmin;
  node->jmax = jmax;
  node->fs = fs;
  node->ispan = (ispan == Py_TRUE) ? true : false;
  node->jspan = (jspan == Py_TRUE) ? true : false;

  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

PyObject *python_sheet(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"func", "imin", "imax", "jmin", "jmax", "fs", "iclose", "jclose", NULL};
  PyObject *func = NULL;
  double imin, imax, jmin, jmax;
  PyObject *ispan = nullptr, *jspan = nullptr;
  double dum1, dum2, fs;
  get_fnas(dum1, dum2, fs);
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "Odddd|dOO", kwlist, &func, &imin, &imax, &jmin,
                                      &jmax, &fs, &ispan, &jspan)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing sheet(func, imin, imax, jmin, jmax)");
    return NULL;
  }
  if (func->ob_type != &PyFunction_Type) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "must specify a function");
    return NULL;
  }

  return python_sheet_core(func, imin, imax, jmin, jmax, fs, ispan, jspan);
}

PyObject *python_text(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE
  auto node = std::make_shared<TextNode>(instance);

  char *kwlist[] = {"text",   "size",   "font", "spacing", "direction", "language", "script",
                    "halign", "valign", "fn",   "fa",      "fs",        NULL};

  double size = 1.0, spacing = 1.0;
  double fn = NAN, fa = NAN, fs = NAN;

  get_fnas(fn, fa, fs);

  const char *text = "", *font = NULL, *direction = "ltr", *language = "en", *script = "latin",
             *valign = "baseline", *halign = "left";

  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "s|dsdsssssddd", kwlist, &text, &size, &font,
                                      &spacing, &direction, &language, &script, &halign, &valign, &fn,
                                      &fa, &fs)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing text(string, ...))");
    return NULL;
  }

  node->params.set_fn(fn);
  node->params.set_fa(fa);
  node->params.set_fs(fs);
  node->params.set_size(size);
  if (text != NULL) node->params.set_text(text);
  node->params.set_spacing(spacing);
  if (font != NULL) node->params.set_font(font);
  if (direction != NULL) node->params.set_direction(direction);
  if (language != NULL) node->params.set_language(language);
  if (script != NULL) node->params.set_script(script);
  if (valign != NULL) node->params.set_halign(halign);
  if (halign != NULL) node->params.set_valign(valign);
  node->params.set_loc(instance->location());

  /*
     node->params.set_documentPath(session->documentRoot());
     }
   */
  node->params.detect_properties();

  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

PyObject *python_textmetrics(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE
  auto node = std::make_shared<TextNode>(instance);

  char *kwlist[] = {"text",     "size",   "font",   "spacing", "direction",
                    "language", "script", "halign", "valign",  NULL};

  double size = 1.0, spacing = 1.0;

  const char *text = "", *font = NULL, *direction = "ltr", *language = "en", *script = "latin",
             *valign = "baseline", *halign = "left";

  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "s|dsdsssss", kwlist, &text, &size, &font, &spacing,
                                      &direction, &language, &script, &valign, &halign)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing textmetrics");
    return NULL;
  }

  FreetypeRenderer::Params ftparams;

  ftparams.set_size(size);
  if (text != NULL) ftparams.set_text(text);
  ftparams.set_spacing(spacing);
  if (font != NULL) ftparams.set_font(font);
  if (direction != NULL) ftparams.set_direction(direction);
  if (language != NULL) ftparams.set_language(language);
  if (script != NULL) ftparams.set_script(script);
  if (valign != NULL) ftparams.set_halign(halign);
  if (halign != NULL) ftparams.set_valign(valign);
  ftparams.set_loc(instance->location());

  FreetypeRenderer::TextMetrics metrics(ftparams);
  if (!metrics.ok) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid Metric");
    return NULL;
  }
  PyObject *offset = pf.PyList_New(2);
  pf.PyList_SetItem(offset, 0, pf.PyFloat_FromDouble(metrics.x_offset));
  pf.PyList_SetItem(offset, 1, pf.PyFloat_FromDouble(metrics.y_offset));

  PyObject *advance = pf.PyList_New(2);
  pf.PyList_SetItem(advance, 0, pf.PyFloat_FromDouble(metrics.advance_x));
  pf.PyList_SetItem(advance, 1, pf.PyFloat_FromDouble(metrics.advance_y));

  PyObject *position = pf.PyList_New(2);
  pf.PyList_SetItem(position, 0, pf.PyFloat_FromDouble(metrics.bbox_x));
  pf.PyList_SetItem(position, 1, pf.PyFloat_FromDouble(metrics.bbox_y));

  PyObject *dims = pf.PyList_New(2);
  pf.PyList_SetItem(dims, 0, pf.PyFloat_FromDouble(metrics.bbox_w));
  pf.PyList_SetItem(dims, 1, pf.PyFloat_FromDouble(metrics.bbox_h));

  PyObject *dict;
  dict = pf.PyDict_New();
  pf.PyDict_SetItemString(dict, "ascent", pf.PyFloat_FromDouble(metrics.ascent));
  pf.PyDict_SetItemString(dict, "descent", pf.PyFloat_FromDouble(metrics.descent));
  pf.PyDict_SetItemString(dict, "offset", offset);
  pf.PyDict_SetItemString(dict, "advance", advance);
  pf.PyDict_SetItemString(dict, "position", position);
  pf.PyDict_SetItemString(dict, "size", dims);
  return (PyObject *)dict;
}

PyObject *python_osversion(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {NULL};
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "", kwlist)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing version()");
    return NULL;
  }

  PyObject *version = pf.PyList_New(3);
  pf.PyList_SetItem(version, 0, pf.PyFloat_FromDouble(OPENSCAD_YEAR));
  pf.PyList_SetItem(version, 1, pf.PyFloat_FromDouble(OPENSCAD_MONTH));
#ifdef OPENSCAD_DAY
  pf.PyList_SetItem(version, 2, pf.PyFloat_FromDouble(OPENSCAD_DAY));
#else
  pf.PyList_SetItem(version, 2, pf.PyFloat_FromDouble(0));
#endif

  return version;
}

PyObject *python_osversion_num(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {NULL};
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "", kwlist)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing version_num()");
    return NULL;
  }

  double version = OPENSCAD_YEAR * 10000 + OPENSCAD_MONTH * 100;
#ifdef OPENSCAD_DAY
  version += OPENSCAD_DAY;
#endif
  return pf.PyFloat_FromDouble(version);
}

PyObject *python_offset_core(PyObject *obj, double r, double delta, PyObject *chamfer, double fn,
                             double fa, double fs)
{
  DECLARE_INSTANCE
  auto node = std::make_shared<OffsetNode>(instance);

  PyObject *dummydict;
  PyTypeObject *type = PyOpenSCADObjectType(obj);
  std::shared_ptr<AbstractNode> child = PyOpenSCADObjectToNodeMulti(obj, &dummydict);
  if (child == NULL) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid type for Object in offset");
    return NULL;
  }

  get_fnas(node->fn, node->fa, node->fs);
  if (!isnan(fn)) node->fn = fn;
  if (!isnan(fa)) node->fa = fa;
  if (!isnan(fs)) node->fs = fs;

  node->delta = 1;
  node->chamfer = false;
  node->join_type = Clipper2Lib::JoinType::Round;
  if (!isnan(r)) {
    node->delta = r;
  } else if (!isnan(delta)) {
    node->delta = delta;
    node->join_type = Clipper2Lib::JoinType::Miter;
    if (chamfer == Py_TRUE) {
      node->chamfer = true;
      node->join_type = Clipper2Lib::JoinType::Square;
    } else if (chamfer == Py_FALSE || chamfer == NULL) node->chamfer = 0;
    else {
      pf.PyErr_SetString(pf.PyExc_TypeError, "Unknown Value for chamfer parameter");
      return NULL;
    }
  }
  node->children.push_back(child);
  return PyOpenSCADObjectFromNode(type, node);
}

PyObject *python_offset(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"obj", "r", "delta", "chamfer", "fn", "fa", "fs", NULL};
  PyObject *obj = NULL;
  double r = NAN, delta = NAN;
  PyObject *chamfer = NULL;
  double fn = NAN, fa = NAN, fs = NAN;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "Od|dOddd", kwlist, &obj, &r, &delta, &chamfer, &fn,
                                      &fa, &fs)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing offset(object,r)");
    return NULL;
  }
  return python_offset_core(obj, r, delta, chamfer, fn, fa, fs);
}

PyObject *python_oo_offset(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"r", "delta", "chamfer", "fn", "fa", "fs", NULL};
  double r = NAN, delta = NAN;
  PyObject *chamfer = NULL;
  double fn = NAN, fa = NAN, fs = NAN;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "d|dOddd", kwlist, &r, &delta, &chamfer, &fn, &fa,
                                      &fs)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing offset(object,r)");
    return NULL;
  }
  return python_offset_core(obj, r, delta, chamfer, fn, fa, fs);
}

PyObject *python_projection_core(PyObject *obj, PyObject *cut, int convexity)
{
  DECLARE_INSTANCE
  auto node = std::make_shared<ProjectionNode>(instance);
  PyObject *dummydict;
  PyTypeObject *type = PyOpenSCADObjectType(obj);
  std::shared_ptr<AbstractNode> child = PyOpenSCADObjectToNodeMulti(obj, &dummydict);
  if (child == NULL) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid type for Object in projection");
    return NULL;
  }
  node->convexity = convexity;
  node->cut_mode = 0;
  if (cut == Py_TRUE) node->cut_mode = 1;
  else if (cut == Py_FALSE) node->cut_mode = 0;
  else {
    pf.PyErr_SetString(pf.PyExc_TypeError, "cut can be either True or false");
    return NULL;
  }

  node->children.push_back(child);
  return PyOpenSCADObjectFromNode(type, node);
}

PyObject *python_projection(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"obj", "cut", "convexity", NULL};
  PyObject *obj = NULL;
  PyObject *cutmode = Py_FALSE;
  long convexity = 2;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O|Ol", kwlist, &obj, &cutmode, &convexity)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing projection(object)");
    return NULL;
  }
  return python_projection_core(obj, cutmode, convexity);
}

PyObject *python_oo_projection(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"cut", "convexity", NULL};
  PyObject *cutmode = Py_FALSE;
  long convexity = 2;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "|Ol", kwlist, &cutmode, &convexity)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing projection(object)");
    return NULL;
  }
  return python_projection_core(obj, cutmode, convexity);
}

PyObject *python_group(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE
  std::shared_ptr<AbstractNode> child;

  auto node = std::make_shared<GroupNode>(instance);

  char *kwlist[] = {"obj", NULL};
  PyObject *obj = NULL;
  PyObject *dummydict;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O!", kwlist, &PyOpenSCADType, &obj)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing group(group)");
    return NULL;
  }
  PyTypeObject *type = PyOpenSCADObjectType(obj);
  child = PyOpenSCADObjectToNode(obj, &dummydict);

  node->children.push_back(child);
  return PyOpenSCADObjectFromNode(type, node);
}

PyObject *python_align_core(PyObject *obj, PyObject *pyrefmat, PyObject *pydstmat)
{
  if (obj->ob_type != &PyOpenSCADType) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Must specify Object as 1st parameter");
    return nullptr;
  }
  PyObject *child_dict = nullptr;
  std::shared_ptr<AbstractNode> dstnode = PyOpenSCADObjectToNode(obj, &child_dict);
  if (dstnode == nullptr) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid align object");
    return Py_NONE;
  }
  DECLARE_INSTANCE
  auto multmatnode = std::make_shared<TransformNode>(instance, "align");
  multmatnode->children.push_back(dstnode);
  Matrix4d mat;
  Matrix4d MT = Matrix4d::Identity();

  if (!python_tomatrix(pyrefmat, mat)) MT = MT * mat;
  if (!python_tomatrix(pydstmat, mat)) MT = MT * mat.inverse();

  multmatnode->matrix = MT;
  multmatnode->setPyName(dstnode->getPyName());

  PyObject *pyresult = PyOpenSCADObjectFromNode(&PyOpenSCADType, multmatnode);
  if (child_dict != nullptr) {
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (pf.PyDict_Next(child_dict, &pos, &key, &value)) {
      //       PyObject* value1 = pf.PyUnicode_AsEncodedString(key, "utf-8", "~");
      //       const char *value_str =  PyBytes_AS_STRING(value1);
      if (!python_tomatrix(value, mat)) {
        mat = MT * mat;
        pf.PyDict_SetItem(((PyOpenSCADObject *)pyresult)->dict, key, python_frommatrix(mat));
      } else pf.PyDict_SetItem(((PyOpenSCADObject *)pyresult)->dict, key, value);
    }
  }
  return pyresult;
}

PyObject *python_align(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"obj", "refmat", "objmat", NULL};
  PyObject *obj = NULL;
  PyObject *pyrefmat = NULL;
  PyObject *pyobjmat = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "OO|O", kwlist, &obj, &pyrefmat, &pyobjmat)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during align");
    return NULL;
  }
  return python_align_core(obj, pyrefmat, pyobjmat);
}

PyObject *python_oo_align(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"refmat", "objmat", NULL};
  PyObject *pyrefmat = NULL;
  PyObject *pyobjmat = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O|O", kwlist, &pyrefmat, &pyobjmat)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during align");
    return NULL;
  }
  return python_align_core(obj, pyrefmat, pyobjmat);
}

PyObject *do_import_python(PyObject *self, PyObject *args, PyObject *kwargs, ImportType type)
{
  DECLARE_INSTANCE
  char *kwlist[] = {"file",   "layer",    "convexity", "origin", "scale", "width",
                    "height", "filename", "center",    "dpi",    "id",    NULL};
  double fn = NAN, fa = NAN, fs = NAN;

  std::string filename;
  const char *v = NULL, *layer = NULL, *id = NULL;
  PyObject *center = NULL;
  int convexity = 2;
  double scale = 1.0, width = 1, height = 1, dpi = 1.0;
  PyObject *origin = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "s|slO!dddsfOddd", kwlist, &v, &layer, &convexity,
                                      pf.PyList_Type, origin, &scale, &width, &height, &center, &dpi,
                                      &id, &fn, &fa, &fs

                                      )) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing osimport(filename)");
    return NULL;
  }
  filename = lookup_file(v == NULL ? "" : v, python_scriptpath.parent_path().u8string(),
                         instance->location().filePath().parent_path().string());
  if (!filename.empty()) handle_dep(filename);
  ImportType actualtype = type;
  if (actualtype == ImportType::UNKNOWN) {
    std::string extraw = fs::path(filename).extension().generic_string();
    std::string ext = boost::algorithm::to_lower_copy(extraw);
    if (ext == ".stl") actualtype = ImportType::STL;
    else if (ext == ".off") actualtype = ImportType::OFF;
    else if (ext == ".obj") actualtype = ImportType::OBJ;
    else if (ext == ".dxf") actualtype = ImportType::DXF;
    else if (ext == ".nef3") actualtype = ImportType::NEF3;
    else if (ext == ".3mf") actualtype = ImportType::_3MF;
    else if (ext == ".amf") actualtype = ImportType::AMF;
    else if (ext == ".svg") actualtype = ImportType::SVG;
    else if (ext == ".stp") actualtype = ImportType::STEP;
    else if (ext == ".step") actualtype = ImportType::STEP;
  }

  auto node = std::make_shared<ImportNode>(instance, actualtype);

  get_fnas(node->fn, node->fa, node->fs);
  if (!isnan(fn)) node->fn = fn;
  if (!isnan(fa)) node->fa = fa;
  if (!isnan(fs)) node->fs = fs;

  node->filename = filename;

  if (layer != NULL) node->layer = layer;
  if (id != NULL) node->id = id;
  node->convexity = convexity;
  if (node->convexity <= 0) node->convexity = 1;

  if (origin != NULL && PyList_CHECK(origin) && pf.PyList_Size(origin) == 2) {
    node->origin_x = pf.PyFloat_AsDouble(pf.PyList_GetItem(origin, 0));
    node->origin_y = pf.PyFloat_AsDouble(pf.PyList_GetItem(origin, 1));
  }

  node->center = 0;
  if (center == Py_TRUE) node->center = 1;

  node->scale = scale;
  if (node->scale <= 0) node->scale = 1;

  node->dpi = ImportNode::SVG_DEFAULT_DPI;
  double val = dpi;
  if (val < 0.001) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Invalid dpi value giving");
    return NULL;
  } else {
    node->dpi = val;
  }

  node->width = width;
  node->height = height;
  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

PyObject *python_oo_clone(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyObject *dict;
  PyObject *obj = NULL;
  char *kwlist[] = {"obj", NULL};

  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &obj)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during clone");
    return NULL;
  }
  std::shared_ptr<AbstractNode> node = PyOpenSCADObjectToNodeMulti(obj, &dict);
  if (node.use_count() > 1) ((PyOpenSCADObject *)self)->node = node->clone();
  else ((PyOpenSCADObject *)self)->node = node;

  if (dict != nullptr) {
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (pf.PyDict_Next(dict, &pos, &key, &value)) {
      pf.PyDict_SetItem(((PyOpenSCADObject *)self)->dict, key, value);
    }
  }
  return Py_NONE;
}

PyObject *python_import(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return do_import_python(self, args, kwargs, ImportType::UNKNOWN);
}

#ifndef OPENSCAD_NOGUI
std::vector<std::string> nimport_downloaded;

extern int curl_download(std::string url, std::string path);
PyObject *python_nimport(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static bool called_already = false;
  char *kwlist[] = {"url", NULL};
  const char *c_url = nullptr;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "s", kwlist, &c_url)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing nimport(filename)");
    return NULL;
  }
  if (c_url == nullptr) return Py_NONE;

  std::string url = c_url;
  std::string filename, path, importcode;
  filename = url.substr(url.find_last_of("/") + 1);
  importcode = "from " + filename.substr(0, filename.find_last_of(".")) + " import *";

  path = PlatformUtils::userLibraryPath() + "/" + filename;
  bool do_download = false;
  if (std::find(nimport_downloaded.begin(), nimport_downloaded.end(), url) == nimport_downloaded.end()) {
    do_download = true;
    nimport_downloaded.push_back(url);
  }

  std::ifstream f(path.c_str());
  if (!f.good()) {
    do_download = true;
  }

  if (do_download) {
    curl_download(url, path);
  }

  pf.PyRun_SimpleString(importcode.c_str());
  return Py_NONE;
}
#endif

void python_str_sub(std::ostringstream& stream, const std::shared_ptr<AbstractNode>& node, int ident)
{
  for (int i = 0; i < ident; i++) stream << "  ";
  stream << node->toString();
  switch (node->children.size()) {
  case 0: stream << ";\n"; break;
  case 1:
    stream << "\n";
    python_str_sub(stream, node->children[0], ident + 1);
    break;
  default:
    stream << "{\n";
    for (const auto& child : node->children) {
      python_str_sub(stream, child, ident + 1);
    }
    for (int i = 0; i < ident; i++) stream << "  ";
    stream << "}\n";
  }
}

PyObject *python_str(PyObject *self)
{
  std::ostringstream stream;
  PyObject *dummydict;
  std::shared_ptr<AbstractNode> node = PyOpenSCADObjectToNode(self, &dummydict);
  if (node != nullptr) stream << "OpenSCAD (" << (int)node->index() << ")";
  else stream << "Invalid OpenSCAD Object";

  return pf.PyUnicode_FromStringAndSize(stream.str().c_str(), stream.str().size());
}

PyObject *python_add_parameter(PyObject *self, PyObject *args, PyObject *kwargs, ImportType type)
{
  char *kwlist[] = {"name", "default", NULL};
  char *name = NULL;
  PyObject *value = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "sO", kwlist, &name, &value)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing add_parameter(name,defval)");
    return NULL;
  }
  bool found = false;
  std::shared_ptr<Literal> lit;
  if (value == Py_TRUE) {
    lit = std::make_shared<Literal>(true, Location::NONE);
    found = true;
  } else if (value == Py_FALSE) {
    lit = std::make_shared<Literal>(false, Location::NONE);
    found = true;
  } else if (PyFloat_CHECK(value)) {
    lit = std::make_shared<Literal>(pf.PyFloat_AsDouble(value), Location::NONE);
    found = true;
  } else if (PyLong_Check(value)) {
    lit = std::make_shared<Literal>(pf.PyLong_AsLong(value) * 1.0, Location::NONE);
    found = true;
  } else if (PyUnicode_Check(value)) {
    PyObject *value1 = pf.PyUnicode_AsEncodedString(value, "utf-8", "~");
    const char *value_str = PyBytes_AS_STRING(value1);
    lit = std::make_shared<Literal>(value_str, Location::NONE);
    found = true;
  }

  if (found) {
    AnnotationList annotationList;
    //    annotationList.push_back(Annotation("Parameter",std::make_shared<Literal>("Parameter")));
    //    annotationList.push_back(Annotation("Description",std::make_shared<Literal>("Description")));
    //    annotationList.push_back(Annotation("Group",std::make_shared<Literal>("Group")));
    auto assignment = std::make_shared<Assignment>(name, lit);
    //    assignment->addAnnotations(&annotationList);
    customizer_parameters.push_back(assignment);
    PyObject *value_effective = value;
    for (unsigned int i = 0; i < customizer_parameters_finished.size(); i++) {
      if (customizer_parameters_finished[i]->getName() == name) {
        auto expr = customizer_parameters_finished[i]->getExpr();
        const auto& lit = std::dynamic_pointer_cast<Literal>(expr);
        if (lit != nullptr) {
          if (lit->isDouble()) value_effective = pf.PyFloat_FromDouble(lit->toDouble());
          if (lit->isString()) value_effective = pf.PyUnicode_FromString(lit->toString().c_str());
        }
      }
    }
    PyObject *maindict = pf.PyModule_GetDict(pythonMainModule.get());
    pf.PyDict_SetItemString(maindict, name, value_effective);
  }
  return Py_NONE;
}

PyObject *python_scad(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE
  char *kwlist[] = {"code", NULL};
  const char *code = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "s", kwlist, &code)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing scad(code)");
    return NULL;
  }

  SourceFile *parsed_file = NULL;
  if (!parse(parsed_file, code, "python", "python", false)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error in SCAD code");
    return Py_NONE;
  }
  parsed_file->handleDependencies(true);

  EvaluationSession session{"python"};
  ContextHandle<BuiltinContext> builtin_context{Context::create<BuiltinContext>(&session)};
  std::shared_ptr<const FileContext> file_context;
  std::shared_ptr<AbstractNode> resultnode = parsed_file->instantiate(*builtin_context, &file_context);
  resultnode = resultnode->clone();  // instmod will go out of scope
  delete parsed_file;
  parsed_file = nullptr;
  return PyOpenSCADObjectFromNode(&PyOpenSCADType, resultnode);
}

PyObject *python_osuse_include(int mode, PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE
  auto empty = std::make_shared<CubeNode>(instance);
  char *kwlist[] = {"file", NULL};
  const char *file = NULL;
  std::ostringstream stream;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "s", kwlist, &file)) {
    if (mode) pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing osinclude(path)");
    else pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing osuse(path)");
    return NULL;
  }
  const std::string filename = lookup_file(file, python_scriptpath.parent_path().u8string(), ".");
  stream << "include <" << filename << ">\n";

  SourceFile *source;
  if (!parse(source, stream.str(), "python", "python", false)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error in SCAD code");
    return Py_NONE;
  }
  if (mode == 0) source->scope->moduleInstantiations.clear();
  source->handleDependencies(true);

  EvaluationSession *session = new EvaluationSession("python");
  ContextHandle<BuiltinContext> builtin_context{Context::create<BuiltinContext>(session)};

  std::shared_ptr<const FileContext> osinclude_context;
  std::shared_ptr<AbstractNode> resultnode =
    source->instantiate(*builtin_context, &osinclude_context);  // TODO keine globakle var, kollision!

  auto scope = source->scope;
  PyOpenSCADObject *result = (PyOpenSCADObject *)PyOpenSCADObjectFromNode(&PyOpenSCADType, empty);

  for (auto mod : source->scope->modules) {  // copy modules
    std::shared_ptr<UserModule> usmod = mod.second;
    InstantiableModule m;
    //    m.defining_context=osinclude_context;
    //    m.module=mod.second.get();
    //    boost::optional<InstantiableModule> res(m);
    pf.PyDict_SetItemString(result->dict, mod.first.c_str(),
                            PyDataObjectFromModule(&PyDataType, filename, mod.first));
  }

  for (auto fun : source->scope->functions) {           // copy functions
    std::shared_ptr<UserFunction> usfunc = fun.second;  // install lambda functions ?
                                                        //    printf("%s\n",fun.first.c_str());
                                                        //    InstantiableModule m;
                                                        //    m.defining_context=osinclude_context;
                                                        //    m.module=mod.second.get();
                                                        //    boost::optional<InstantiableModule> res(m);
    //    pf.PyDict_SetItemString(result->dict, mod.first.c_str(),PyDataObjectFromModule(&PyDataType, res
    //    ));
  }

  for (auto ass : source->scope->assignments) {  // copy assignments
                                                 //    printf("Var %s\n",ass->getName().c_str());
    const std::shared_ptr<Expression> expr = ass->getExpr();
    Value val = expr->evaluate(osinclude_context);
    if (val.isDefined()) {
      PyObject *res = python_fromopenscad(std::move(val));
      pf.PyDict_SetItemString(result->dict, ass->getName().c_str(), res);
    }
  }
  // SourceFileCache::instance()->clear();
  return (PyObject *)result;
}

PyObject *python_osuse(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_osuse_include(0, self, args, kwargs);
}

PyObject *python_osinclude(PyObject *self, PyObject *args, PyObject *kwargs)
{
  LOG(message_group::Deprecated, "osinclude  is deprecated, please use osuse() instead");
  return python_osuse_include(1, self, args, kwargs);
}

PyObject *python_debug_modifier(PyObject *arg, int mode)
{
  DECLARE_INSTANCE
  PyObject *dummydict;
  PyTypeObject *type = PyOpenSCADObjectType(arg);
  auto child = PyOpenSCADObjectToNode(arg, &dummydict);
  switch (mode) {
  case 0: instance->tag_highlight = true; break;   // #
  case 1: instance->tag_background = true; break;  // %
  case 2: instance->tag_root = true; break;        // !
  }
  auto node = std::make_shared<CsgOpNode>(instance, OpenSCADOperator::UNION);
  node->children.push_back(child);
  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);  // TODO 1st loswerden
}

PyObject *python_debug_modifier_func(PyObject *self, PyObject *args, PyObject *kwargs, int mode)
{
  char *kwlist[] = {"obj", NULL};
  PyObject *obj = NULL;
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "O!", kwlist, &PyOpenSCADType, &obj)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing group(group)");
    return NULL;
  }
  return python_debug_modifier(obj, mode);
}

PyObject *python_debug_modifier_func_oo(PyObject *obj, PyObject *args, PyObject *kwargs, int mode)
{
  char *kwlist[] = {NULL};
  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "", kwlist)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing group(group)");
    return NULL;
  }
  return python_debug_modifier(obj, mode);
}
PyObject *python_highlight(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_debug_modifier_func(self, args, kwargs, 0);
}
PyObject *python_oo_highlight(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_debug_modifier_func_oo(self, args, kwargs, 0);
}
PyObject *python_background(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_debug_modifier_func(self, args, kwargs, 1);
}
PyObject *python_oo_background(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_debug_modifier_func_oo(self, args, kwargs, 1);
}
PyObject *python_only(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_debug_modifier_func(self, args, kwargs, 2);
}
PyObject *python_oo_only(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_debug_modifier_func_oo(self, args, kwargs, 2);
}

PyObject *python_nb_invert(PyObject *arg) { return python_debug_modifier(arg, 0); }
PyObject *python_nb_neg(PyObject *arg) { return python_debug_modifier(arg, 1); }
PyObject *python_nb_pos(PyObject *arg) { return python_debug_modifier(arg, 2); }

#ifndef OPENSCAD_NOGUI
extern void add_menuitem_trampoline(const char *menuname, const char *itemname, const char *callback);
PyObject *python_add_menuitem(PyObject *self, PyObject *args, PyObject *kwargs, int mode)
{
  char *kwlist[] = {"menuname", "itemname", "callback", NULL};
  const char *menuname = nullptr, *itemname = nullptr, *callback = nullptr;

  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "sss", kwlist, &menuname, &itemname, &callback)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing add_menuitem");
    return NULL;
  }
  add_menuitem_trampoline(menuname, itemname, callback);
  return Py_NONE;
}
#endif

PyObject *python_model(PyObject *self, PyObject *args, PyObject *kwargs, int mode)
{
  char *kwlist[] = {NULL};

  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "", kwlist)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing model");
    return NULL;
  }
  if (genlang_result_node == nullptr) return Py_NONE;
  return PyOpenSCADObjectFromNode(&PyOpenSCADType, genlang_result_node);
}

PyObject *python_modelpath(PyObject *self, PyObject *args, PyObject *kwargs, int mode)
{
  char *kwlist[] = {NULL};

  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "", kwlist)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing model");
    return NULL;
  }
  return pf.PyUnicode_FromString(python_scriptpath.u8string().c_str());
}

PyObject *python_oo_dict(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyObject *dict = ((PyOpenSCADObject *)self)->dict;
  Py_INCREF(dict);
  return dict;
}

int PyDict_SetDefaultRef(PyObject *d, PyObject *key, PyObject *default_value, PyObject **result)
{
  PyDict_SetDefault(d, key, default_value);
  return 0;
}

int type_add_method(PyTypeObject *type, PyMethodDef *meth)  // from typeobject.c
{
  PyObject *descr;
  int isdescr = 1;
  if (meth->ml_flags & METH_CLASS) {
    if (meth->ml_flags & METH_STATIC) {
      //      pf.PyErr_SetString(PyExc_ValueError, "method cannot be both class and static");
      return -1;
    }
    descr = PyDescr_NewClassMethod(type, meth);
  } else if (meth->ml_flags & METH_STATIC) {
    PyObject *cfunc = PyCFunction_NewEx(meth, (PyObject *)type, NULL);
    if (cfunc == NULL) {
      return -1;
    }
    descr = PyStaticMethod_New(cfunc);
    isdescr = 0;  // PyStaticMethod is not PyDescrObject
    Py_DECREF(cfunc);
  } else {
    descr = PyDescr_NewMethod(type, meth);
  }
  if (descr == NULL) {
    return -1;
  }

  PyObject *name;
  if (isdescr) {
    name = PyDescr_NAME(descr);
  } else {
    name = PyUnicode_FromString(meth->ml_name);
    if (name == NULL) {
      Py_DECREF(descr);
      return -1;
    }
  }

  int err;
  PyObject *dict = type->tp_dict;
  if (!(meth->ml_flags & METH_COEXIST)) {
    err = PyDict_SetDefaultRef(dict, name, descr, NULL) < 0;
  } else {
    err = 0;  // TODO (pf.PyDict_SetItem(dict, name, descr) < 0);
  }
  if (!isdescr) {
    Py_DECREF(name);
  }
  Py_DECREF(descr);
  if (err) {
    return -1;  // return here
  }
  return 0;
}

std::vector<PyObject *> python_member_callables;
std::vector<std::string> python_member_names;
int python_member_callind;

PyObject *python_member_trampoline(PyObject *self, PyObject *args, PyObject *kwargs)
{
  int n = pf.PyTuple_Size(args);
  PyObject *newargs = pf.PyTuple_New(n + 1);
  pf.PyTuple_SetItem(newargs, 0, self);
  for (int i = 0; i < n; i++) pf.PyTuple_SetItem(newargs, i + 1, pf.PyTuple_GetItem(args, i));

  return PyObject_Call(python_member_callables[python_member_callind], newargs, kwargs);
}

#define PYTHON_MAX_USERMEMBERS 20

PyObject *python_member_trampoline_0(PyObject *self, PyObject *args, PyObject *kwargs)
{
  python_member_callind = 0;
  return python_member_trampoline(self, args, kwargs);
}
PyObject *python_member_trampoline_1(PyObject *self, PyObject *args, PyObject *kwargs)
{
  python_member_callind = 1;
  return python_member_trampoline(self, args, kwargs);
}
PyObject *python_member_trampoline_2(PyObject *self, PyObject *args, PyObject *kwargs)
{
  python_member_callind = 2;
  return python_member_trampoline(self, args, kwargs);
}
PyObject *python_member_trampoline_3(PyObject *self, PyObject *args, PyObject *kwargs)
{
  python_member_callind = 3;
  return python_member_trampoline(self, args, kwargs);
}
PyObject *python_member_trampoline_4(PyObject *self, PyObject *args, PyObject *kwargs)
{
  python_member_callind = 4;
  return python_member_trampoline(self, args, kwargs);
}
PyObject *python_member_trampoline_5(PyObject *self, PyObject *args, PyObject *kwargs)
{
  python_member_callind = 5;
  return python_member_trampoline(self, args, kwargs);
}
PyObject *python_member_trampoline_6(PyObject *self, PyObject *args, PyObject *kwargs)
{
  python_member_callind = 6;
  return python_member_trampoline(self, args, kwargs);
}
PyObject *python_member_trampoline_7(PyObject *self, PyObject *args, PyObject *kwargs)
{
  python_member_callind = 7;
  return python_member_trampoline(self, args, kwargs);
}
PyObject *python_member_trampoline_8(PyObject *self, PyObject *args, PyObject *kwargs)
{
  python_member_callind = 8;
  return python_member_trampoline(self, args, kwargs);
}
PyObject *python_member_trampoline_9(PyObject *self, PyObject *args, PyObject *kwargs)
{
  python_member_callind = 9;
  return python_member_trampoline(self, args, kwargs);
}
PyObject *python_member_trampoline_10(PyObject *self, PyObject *args, PyObject *kwargs)
{
  python_member_callind = 10;
  return python_member_trampoline(self, args, kwargs);
}
PyObject *python_member_trampoline_11(PyObject *self, PyObject *args, PyObject *kwargs)
{
  python_member_callind = 11;
  return python_member_trampoline(self, args, kwargs);
}
PyObject *python_member_trampoline_12(PyObject *self, PyObject *args, PyObject *kwargs)
{
  python_member_callind = 12;
  return python_member_trampoline(self, args, kwargs);
}
PyObject *python_member_trampoline_13(PyObject *self, PyObject *args, PyObject *kwargs)
{
  python_member_callind = 13;
  return python_member_trampoline(self, args, kwargs);
}
PyObject *python_member_trampoline_14(PyObject *self, PyObject *args, PyObject *kwargs)
{
  python_member_callind = 14;
  return python_member_trampoline(self, args, kwargs);
}
PyObject *python_member_trampoline_15(PyObject *self, PyObject *args, PyObject *kwargs)
{
  python_member_callind = 15;
  return python_member_trampoline(self, args, kwargs);
}
PyObject *python_member_trampoline_16(PyObject *self, PyObject *args, PyObject *kwargs)
{
  python_member_callind = 16;
  return python_member_trampoline(self, args, kwargs);
}
PyObject *python_member_trampoline_17(PyObject *self, PyObject *args, PyObject *kwargs)
{
  python_member_callind = 17;
  return python_member_trampoline(self, args, kwargs);
}
PyObject *python_member_trampoline_18(PyObject *self, PyObject *args, PyObject *kwargs)
{
  python_member_callind = 18;
  return python_member_trampoline(self, args, kwargs);
}
PyObject *python_member_trampoline_19(PyObject *self, PyObject *args, PyObject *kwargs)
{
  python_member_callind = 19;
  return python_member_trampoline(self, args, kwargs);
}

PyObject *python_memberfunction(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"membername", "memberfunc", "docstring", NULL};
  char *membername = nullptr;
  PyObject *memberfunc = nullptr;
  char *memberdoc = nullptr;

  if (!pf.PyArg_ParseTupleAndKeywords(args, kwargs, "sO|s", kwlist, &membername, &memberfunc)) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Error during parsing member");
    return NULL;
  }
  std::string member_name = membername;
  int curind = std::find(python_member_names.begin(), python_member_names.end(), member_name) -
               python_member_names.begin();

  if (memberdoc == nullptr) {
    memberdoc = "Added by member function";
  }

  if (curind >= PYTHON_MAX_USERMEMBERS) {
    pf.PyErr_SetString(pf.PyExc_TypeError, "Maximum user member amount reached");
    return NULL;
  }

  PyCFunction next_trampoline;
  switch (curind) {
  case 0:  next_trampoline = (PyCFunction)python_member_trampoline_0; break;
  case 1:  next_trampoline = (PyCFunction)python_member_trampoline_1; break;
  case 2:  next_trampoline = (PyCFunction)python_member_trampoline_2; break;
  case 3:  next_trampoline = (PyCFunction)python_member_trampoline_3; break;
  case 4:  next_trampoline = (PyCFunction)python_member_trampoline_4; break;
  case 5:  next_trampoline = (PyCFunction)python_member_trampoline_5; break;
  case 6:  next_trampoline = (PyCFunction)python_member_trampoline_6; break;
  case 7:  next_trampoline = (PyCFunction)python_member_trampoline_7; break;
  case 8:  next_trampoline = (PyCFunction)python_member_trampoline_8; break;
  case 9:  next_trampoline = (PyCFunction)python_member_trampoline_9; break;
  case 10: next_trampoline = (PyCFunction)python_member_trampoline_10; break;
  case 11: next_trampoline = (PyCFunction)python_member_trampoline_11; break;
  case 12: next_trampoline = (PyCFunction)python_member_trampoline_12; break;
  case 13: next_trampoline = (PyCFunction)python_member_trampoline_13; break;
  case 14: next_trampoline = (PyCFunction)python_member_trampoline_14; break;
  case 15: next_trampoline = (PyCFunction)python_member_trampoline_15; break;
  case 16: next_trampoline = (PyCFunction)python_member_trampoline_16; break;
  case 17: next_trampoline = (PyCFunction)python_member_trampoline_17; break;
  case 18: next_trampoline = (PyCFunction)python_member_trampoline_18; break;
  case 19: next_trampoline = (PyCFunction)python_member_trampoline_19; break;
  default: next_trampoline = nullptr;
  }

  PyMethodDef *meth = (PyMethodDef *)malloc(sizeof(PyMethodDef));  // never freed
  meth->ml_name = strdup(membername);
  meth->ml_meth = next_trampoline;
  meth->ml_flags = METH_VARARGS | METH_KEYWORDS;
  meth->ml_doc = memberdoc;
  if (type_add_method(&PyOpenSCADType, meth) < 0) return Py_NONE;

  Py_INCREF(memberfunc);  // needed because pythons garbage collector eats it when not used.
  if (curind < python_member_names.size()) {
    python_member_callables[curind] = memberfunc;
  } else {
    python_member_names.push_back(member_name);
    python_member_callables.push_back(memberfunc);
  }

  return Py_NONE;
}

PyMethodDef PyOpenSCADFunctions[] = {
  {"edge", (PyCFunction)python_edge, METH_VARARGS | METH_KEYWORDS, "Create Edge."},
  {"square", (PyCFunction)python_square, METH_VARARGS | METH_KEYWORDS, "Create Square."},
  {"circle", (PyCFunction)python_circle, METH_VARARGS | METH_KEYWORDS, "Create Circle."},
  {"polygon", (PyCFunction)python_polygon, METH_VARARGS | METH_KEYWORDS, "Create Polygon."},
  {"spline", (PyCFunction)python_spline, METH_VARARGS | METH_KEYWORDS, "Create Spline."},
  {"text", (PyCFunction)python_text, METH_VARARGS | METH_KEYWORDS, "Create Text."},
  {"textmetrics", (PyCFunction)python_textmetrics, METH_VARARGS | METH_KEYWORDS, "Get textmetrics."},

  {"cube", (PyCFunction)python_cube, METH_VARARGS | METH_KEYWORDS, "Create Cube."},
  {"cylinder", (PyCFunction)python_cylinder, METH_VARARGS | METH_KEYWORDS, "Create Cylinder."},
  {"sphere", (PyCFunction)python_sphere, METH_VARARGS | METH_KEYWORDS, "Create Sphere."},
  {"polyhedron", (PyCFunction)python_polyhedron, METH_VARARGS | METH_KEYWORDS, "Create Polyhedron."},
#ifdef ENABLE_LIBFIVE
  {"frep", (PyCFunction)python_frep, METH_VARARGS | METH_KEYWORDS, "Create F-Rep."},
  {"ifrep", (PyCFunction)python_ifrep, METH_VARARGS | METH_KEYWORDS, "Create Inverse F-Rep."},
#endif

  {"translate", (PyCFunction)python_translate, METH_VARARGS | METH_KEYWORDS, "Move  Object."},
  {"right", (PyCFunction)python_right, METH_VARARGS | METH_KEYWORDS, "Move  Object."},
  {"left", (PyCFunction)python_left, METH_VARARGS | METH_KEYWORDS, "Move Left Object."},
  {"back", (PyCFunction)python_back, METH_VARARGS | METH_KEYWORDS, "Move Back Object."},
  {"front", (PyCFunction)python_front, METH_VARARGS | METH_KEYWORDS, "Move Front Object."},
  {"up", (PyCFunction)python_up, METH_VARARGS | METH_KEYWORDS, "Move Up Object."},
  {"down", (PyCFunction)python_down, METH_VARARGS | METH_KEYWORDS, "Move Down Object."},
  {"rotx", (PyCFunction)python_rotx, METH_VARARGS | METH_KEYWORDS, "Rotate X Object."},
  {"roty", (PyCFunction)python_roty, METH_VARARGS | METH_KEYWORDS, "Rotate Y Object."},
  {"rotz", (PyCFunction)python_rotz, METH_VARARGS | METH_KEYWORDS, "Rotate Z Object."},
  {"rotate", (PyCFunction)python_rotate, METH_VARARGS | METH_KEYWORDS, "Rotate Object."},
  {"scale", (PyCFunction)python_scale, METH_VARARGS | METH_KEYWORDS, "Scale Object."},
  {"mirror", (PyCFunction)python_mirror, METH_VARARGS | METH_KEYWORDS, "Mirror Object."},
  {"multmatrix", (PyCFunction)python_multmatrix, METH_VARARGS | METH_KEYWORDS, "Multmatrix Object."},
  {"divmatrix", (PyCFunction)python_divmatrix, METH_VARARGS | METH_KEYWORDS, "Divmatrix Object."},
  {"offset", (PyCFunction)python_offset, METH_VARARGS | METH_KEYWORDS, "Offset Object."},
#if defined(ENABLE_EXPERIMENTAL) && defined(ENABLE_CGAL)
  {"roof", (PyCFunction)python_roof, METH_VARARGS | METH_KEYWORDS, "Roof Object."},
#endif
  {"pull", (PyCFunction)python_pull, METH_VARARGS | METH_KEYWORDS, "Pull apart Object."},
  {"wrap", (PyCFunction)python_wrap, METH_VARARGS | METH_KEYWORDS, "Wrap Object around cylidner."},
  {"color", (PyCFunction)python_color, METH_VARARGS | METH_KEYWORDS, "Color Object."},
  {"output", (PyCFunction)python_output, METH_VARARGS | METH_KEYWORDS, "Output the result."},
  {"show", (PyCFunction)python_show, METH_VARARGS | METH_KEYWORDS, "Show the result."},
  {"separate", (PyCFunction)python_separate, METH_VARARGS | METH_KEYWORDS, "Split into separate parts."},
  {"export", (PyCFunction)python_export, METH_VARARGS | METH_KEYWORDS, "Export the result."},
  {"find_face", (PyCFunction)python_find_face, METH_VARARGS | METH_KEYWORDS, "find_face."},
  {"sitonto", (PyCFunction)python_sitonto, METH_VARARGS | METH_KEYWORDS, "sitonto"},

  {"linear_extrude", (PyCFunction)python_linear_extrude, METH_VARARGS | METH_KEYWORDS,
   "Linear_extrude Object."},
  {"rotate_extrude", (PyCFunction)python_rotate_extrude, METH_VARARGS | METH_KEYWORDS,
   "Rotate_extrude Object."},
  {"path_extrude", (PyCFunction)python_path_extrude, METH_VARARGS | METH_KEYWORDS,
   "Path_extrude Object."},
  {"skin", (PyCFunction)python_skin, METH_VARARGS | METH_KEYWORDS, "Path_extrude Object."},

  {"union", (PyCFunction)python_union, METH_VARARGS | METH_KEYWORDS, "Union Object."},
  {"difference", (PyCFunction)python_difference, METH_VARARGS | METH_KEYWORDS, "Difference Object."},
  {"intersection", (PyCFunction)python_intersection, METH_VARARGS | METH_KEYWORDS,
   "Intersection Object."},
  {"hull", (PyCFunction)python_hull, METH_VARARGS | METH_KEYWORDS, "Hull Object."},
  {"minkowski", (PyCFunction)python_minkowski, METH_VARARGS | METH_KEYWORDS, "Minkowski Object."},
  {"fill", (PyCFunction)python_fill, METH_VARARGS | METH_KEYWORDS, "Fill Object."},
  {"resize", (PyCFunction)python_resize, METH_VARARGS | METH_KEYWORDS, "Resize Object."},
  {"concat", (PyCFunction)python_concat, METH_VARARGS | METH_KEYWORDS, "Concatenate Object."},

  {"highlight", (PyCFunction)python_highlight, METH_VARARGS | METH_KEYWORDS, "Highlight Object."},
  {"background", (PyCFunction)python_background, METH_VARARGS | METH_KEYWORDS, "Background Object."},
  {"only", (PyCFunction)python_only, METH_VARARGS | METH_KEYWORDS, "Only Object."},

  {"projection", (PyCFunction)python_projection, METH_VARARGS | METH_KEYWORDS, "Projection Object."},
  {"surface", (PyCFunction)python_surface, METH_VARARGS | METH_KEYWORDS, "Surface Object."},
  {"sheet", (PyCFunction)python_sheet, METH_VARARGS | METH_KEYWORDS, "Sheet Object."},
  {"mesh", (PyCFunction)python_mesh, METH_VARARGS | METH_KEYWORDS, "exports mesh."},
  {"bbox", (PyCFunction)python_bbox, METH_VARARGS | METH_KEYWORDS, "caluculate bbox of object."},
  {"size", (PyCFunction)python_size, METH_VARARGS | METH_KEYWORDS, "get size dimensions of object."},
  {"position", (PyCFunction)python_position, METH_VARARGS | METH_KEYWORDS,
   "get position (minimum coordinates) of object."},
  {"faces", (PyCFunction)python_faces, METH_VARARGS | METH_KEYWORDS, "exports a list of faces."},
  {"edges", (PyCFunction)python_edges, METH_VARARGS | METH_KEYWORDS,
   "exports a list of edges from a face."},
  {"explode", (PyCFunction)python_explode, METH_VARARGS | METH_KEYWORDS,
   "explode a solid with a vector"},
  {"oversample", (PyCFunction)python_oversample, METH_VARARGS | METH_KEYWORDS, "oversample."},
  {"debug", (PyCFunction)python_debug, METH_VARARGS | METH_KEYWORDS, "debug a face."},
  {"fillet", (PyCFunction)python_fillet, METH_VARARGS | METH_KEYWORDS, "fillet."},

  {"group", (PyCFunction)python_group, METH_VARARGS | METH_KEYWORDS, "Group Object."},
  {"render", (PyCFunction)python_render, METH_VARARGS | METH_KEYWORDS, "Render Object."},
  {"osimport", (PyCFunction)python_import, METH_VARARGS | METH_KEYWORDS, "Import Object."},
  {"osuse", (PyCFunction)python_osuse, METH_VARARGS | METH_KEYWORDS, "Use OpenSCAD Library."},
  {"osinclude", (PyCFunction)python_osinclude, METH_VARARGS | METH_KEYWORDS,
   "Include OpenSCAD Library."},
  {"version", (PyCFunction)python_osversion, METH_VARARGS | METH_KEYWORDS, "Output openscad Version."},
  {"version_num", (PyCFunction)python_osversion_num, METH_VARARGS | METH_KEYWORDS,
   "Output openscad Version."},
  {"add_parameter", (PyCFunction)python_add_parameter, METH_VARARGS | METH_KEYWORDS,
   "Add Parameter for Customizer."},
  {"scad", (PyCFunction)python_scad, METH_VARARGS | METH_KEYWORDS, "Source OpenSCAD code."},
  {"align", (PyCFunction)python_align, METH_VARARGS | METH_KEYWORDS, "Align Object to another."},
#ifndef OPENSCAD_NOGUI
  {"add_menuitem", (PyCFunction)python_add_menuitem, METH_VARARGS | METH_KEYWORDS,
   "Add Menuitem to the the openscad window."},
  {"nimport", (PyCFunction)python_nimport, METH_VARARGS | METH_KEYWORDS, "Import Networked Object."},
#endif
  {"model", (PyCFunction)python_model, METH_VARARGS | METH_KEYWORDS, "Yield Model"},
  {"modelpath", (PyCFunction)python_modelpath, METH_VARARGS | METH_KEYWORDS,
   "Returns absolute Path to script"},
  {"memberfunction", (PyCFunction)python_memberfunction, METH_VARARGS | METH_KEYWORDS,
   "Registers additional openscad memberfunction functions"},
  {"marked", (PyCFunction)python_marked, METH_VARARGS | METH_KEYWORDS, "Create a marked value."},
  {"Sin", (PyCFunction)python_sin, METH_VARARGS | METH_KEYWORDS, "Calculate sin."},
  {"Cos", (PyCFunction)python_cos, METH_VARARGS | METH_KEYWORDS, "Calculate cos."},
  {"Tan", (PyCFunction)python_tan, METH_VARARGS | METH_KEYWORDS, "Calculate tan."},
  {"Asin", (PyCFunction)python_asin, METH_VARARGS | METH_KEYWORDS, "Calculate asin."},
  {"Acos", (PyCFunction)python_acos, METH_VARARGS | METH_KEYWORDS, "Calculate acos."},
  {"Atan", (PyCFunction)python_atan, METH_VARARGS | METH_KEYWORDS, "Calculate atan."},
  {"norm", (PyCFunction)python_norm, METH_VARARGS | METH_KEYWORDS, "Calculate vector size."},
  {"dot", (PyCFunction)python_dot, METH_VARARGS | METH_KEYWORDS, "Calculate dot product."},
  {"cross", (PyCFunction)python_cross, METH_VARARGS | METH_KEYWORDS, "Calculate cross product."},
  {NULL, NULL, 0, NULL}};

#define OO_METHOD_ENTRY(name, desc) \
  {#name, (PyCFunction)python_oo_##name, METH_VARARGS | METH_KEYWORDS, desc},

PyMethodDef PyOpenSCADMethods[] = {
  OO_METHOD_ENTRY(translate, "Move Object") OO_METHOD_ENTRY(rotate, "Rotate Object") OO_METHOD_ENTRY(
    right, "Right Object") OO_METHOD_ENTRY(left, "Left Object") OO_METHOD_ENTRY(back, "Back Object")
    OO_METHOD_ENTRY(front, "Front Object") OO_METHOD_ENTRY(up, "Up Object") OO_METHOD_ENTRY(
      down, "Lower Object")

      OO_METHOD_ENTRY(union, "Union Object") OO_METHOD_ENTRY(
        difference, "Difference Object") OO_METHOD_ENTRY(intersection, "Intersection Object")

        OO_METHOD_ENTRY(rotx, "Rotx Object") OO_METHOD_ENTRY(roty, "Roty Object") OO_METHOD_ENTRY(
          rotz, "Rotz Object")

          OO_METHOD_ENTRY(scale, "Scale Object") OO_METHOD_ENTRY(mirror, "Mirror Object")
            OO_METHOD_ENTRY(multmatrix, "Multmatrix Object") OO_METHOD_ENTRY(
              divmatrix, "Divmatrix Object") OO_METHOD_ENTRY(offset, "Offset Object")
#if defined(ENABLE_EXPERIMENTAL) && defined(ENABLE_CGAL)
              OO_METHOD_ENTRY(roof, "Roof Object")
#endif
                OO_METHOD_ENTRY(color, "Color Object") OO_METHOD_ENTRY(
                  separate, "Split into separate Objects") OO_METHOD_ENTRY(export, "Export Object")
                  OO_METHOD_ENTRY(find_face, "Find Face") OO_METHOD_ENTRY(sitonto, "Sit onto")

                    OO_METHOD_ENTRY(linear_extrude, "Linear_extrude Object")
                      OO_METHOD_ENTRY(rotate_extrude, "Rotate_extrude Object") OO_METHOD_ENTRY(
                        path_extrude, "Path_extrude Object") OO_METHOD_ENTRY(resize, "Resize Object")

                        OO_METHOD_ENTRY(explode, "Explode a solid with a vector") OO_METHOD_ENTRY(
                          mesh, "Mesh Object") OO_METHOD_ENTRY(bbox, "Evaluate Bound Box of object")
                          OO_METHOD_ENTRY(faces, "Create Faces list") OO_METHOD_ENTRY(
                            edges, "Create Edges list") OO_METHOD_ENTRY(oversample, "Oversample Object")
                            OO_METHOD_ENTRY(debug, "Debug Object Faces") OO_METHOD_ENTRY(
                              repair, "Make solid watertight") OO_METHOD_ENTRY(fillet, "Fillet Object")
                              OO_METHOD_ENTRY(align, "Align Object to another")

                                OO_METHOD_ENTRY(highlight, "Highlight Object") OO_METHOD_ENTRY(
                                  background, "Background Object") OO_METHOD_ENTRY(only, "Only Object")
                                  OO_METHOD_ENTRY(show, "Show Object")
                                    OO_METHOD_ENTRY(projection, "Projection Object")
                                      OO_METHOD_ENTRY(pull, "Pull Obejct apart")
                                        OO_METHOD_ENTRY(wrap, "Wrap Object around Cylinder")
                                          OO_METHOD_ENTRY(render, "Render Object")
                                            OO_METHOD_ENTRY(clone, "Clone Object") OO_METHOD_ENTRY(
                                              dict, "return all dictionary"){NULL, NULL, 0, NULL}};

PyNumberMethods PyOpenSCADNumbers = {
  python_nb_add,        // binaryfunc nb_add
  python_nb_subtract,   // binaryfunc nb_subtract
  python_nb_mul,        // binaryfunc nb_multiply
  python_nb_remainder,  // binaryfunc nb_remainder
  0,                    // binaryfunc nb_divmod
  0,                    // ternaryfunc nb_power
  python_nb_neg,        // unaryfunc nb_negative
  python_nb_pos,        // unaryfunc nb_positive
  0,                    // unaryfunc nb_absolute
  0,                    // inquiry nb_bool
  python_nb_invert,     // unaryfunc nb_invert
  0,                    // binaryfunc nb_lshift
  0,                    // binaryfunc nb_rshift
  python_nb_and,        // binaryfunc nb_and
  python_nb_xor,        // binaryfunc nb_xor
  python_nb_or,         // binaryfunc nb_or
  0,                    // unaryfunc nb_int
  0,                    // void *nb_reserved
  0,                    // unaryfunc nb_float

  0,  // binaryfunc nb_inplace_add
  0,  // binaryfunc nb_inplace_subtract
  0,  // binaryfunc nb_inplace_multiply
  0,  // binaryfunc nb_inplace_remainder
  0,  // ternaryfunc nb_inplace_power
  0,  // binaryfunc nb_inplace_lshift
  0,  // binaryfunc nb_inplace_rshift
  0,  // binaryfunc nb_inplace_and
  0,  // binaryfunc nb_inplace_xor
  0,  // binaryfunc nb_inplace_or

  0,  // binaryfunc nb_floor_divide
  0,  // binaryfunc nb_true_divide
  0,  // binaryfunc nb_inplace_floor_divide
  0,  // binaryfunc nb_inplace_true_divide

  0,  // unaryfunc nb_index

  python_nb_matmult,  // binaryfunc nb_matrix_multiply
  0                   // binaryfunc nb_inplace_matrix_multiply
};

PyMappingMethods PyOpenSCADMapping = {0, python__getitem__, python__setitem__};
