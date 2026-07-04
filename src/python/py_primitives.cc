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

#include "linalg.h"
#include "GeometryUtils.h"
#include "genlang/genlang.h"
#include <Python.h>
#include "pyopenscad.h"
#include "pyfunctions.h"
#include "pyconversion.h"
#include "primitives.h"
#include "pydata.h"
// #include "Geometry.h"
#include "PolySet.h"
#include "PolySetBuilder.h"
#include "GeometryEvaluator.h"
#ifdef ENABLE_LIBFIVE
#include "FrepNode.h"
#endif
#include "core/FreetypeRenderer.h"
#include "core/TextNode.h"

PyObject *python_edge(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE();
  auto node = std::make_shared<EdgeNode>(instance);

  char *kwlist[] = {"size", "center", NULL};
  double size = 1;

  PyObject *center = NULL;

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|dO", kwlist, &size, &center)) {
    PyErr_SetString(PyExc_TypeError, "Error during parsing edge(size)");
    return NULL;
  }

  if (size < 0) {
    PyErr_SetString(PyExc_TypeError, "Edge Length must be positive");
    return NULL;
  }
  node->size = size;
  if (center == Py_False || center == NULL)
    ;
  else if (center == Py_True) {
    node->center = true;
  }
  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

PyObject *python_marked(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE();
  auto node = std::make_shared<EdgeNode>(instance);

  char *kwlist[] = {"value", NULL};
  double value = 0.0;

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "d", kwlist, &value)) {
    PyErr_SetString(PyExc_TypeError, "Error during parsing marked(value)");
    return NULL;
  }
  return PyDataObjectFromValue(&PyDataType, value);
}

PyObject *python_cube(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE();
  auto node = std::make_shared<CubeNode>(instance);

  char *kwlist[] = {"size", "center", NULL};
  PyObject *size = NULL;

  PyObject *center = NULL;

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|OO", kwlist, &size, &center)) {
    PyErr_SetString(PyExc_TypeError, "Error during parsing cube(size)");
    return NULL;
  }

  if (size != NULL) {
    if (python_vectorval(size, 3, 3, &(node->dim[0]), &(node->dim[1]), &(node->dim[2]), nullptr,
                         &(node->dragflags))) {
      PyErr_SetString(PyExc_TypeError, "Invalid Cube dimensions");
      return NULL;
    }
  }
  if (node->dim[0] <= 0 || node->dim[1] <= 0 || node->dim[2] <= 0) {
    PyErr_SetString(PyExc_TypeError, "Cube Dimensions must be positive");
    return NULL;
  }
  for (int i = 0; i < 3; i++) node->center[i] = 1;
  if (center == Py_False || center == NULL)
    ;
  else if (center == Py_True) {
    for (int i = 0; i < 3; i++) node->center[i] = 0;
  } else if (PyUnicode_Check(center)) {
    std::string centerstr;
    if (!python_pyobject_to_utf8(center, centerstr, "cube() center")) {
      return NULL;
    }
    if (centerstr.size() != 3) {
      PyErr_SetString(PyExc_TypeError, "Center code must be exactly 3 characters");
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
        PyErr_SetString(PyExc_TypeError, "Center code chars not recognized, must be + - or 0");
        return NULL;
      }
    }
  } else {
    PyErr_SetString(PyExc_TypeError, "Unknown Value for center parameter");
    return NULL;
  }
  python_retrieve_pyname(node);
  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

int sphereCalcIndInt(PyObject *func, Vector3d& dir)
{
  dir.normalize();
  PyObject *dir_p = PyList_New(3);
  for (int i = 0; i < 3; i++) PyList_SetItem(dir_p, i, PyFloat_FromDouble(dir[i]));
  PyObject *args = PyTuple_Pack(1, dir_p);
  PyObject *len_p = PyObject_CallObject(func, args);
  double len = 0;
  if (len_p == nullptr) {
    std::string errorstr;
    python_catch_error(errorstr);
    PyErr_SetString(PyExc_TypeError, errorstr.c_str());
    LOG(message_group::Error, errorstr.c_str());
    return 1;
  }
  python_numberval(len_p, &len, nullptr, 0);
  dir *= len;
  return 0;
}

int sphereCalcInd(PolySetBuilder& builder, std::vector<Vector3d>& vertices, PyObject *func, Vector3d dir)
{
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
    int error;
    auto edge_db = createEdgeDb(ps->indices, error);
    if (error)
      LOG(message_group::Warning,
          "Resulting sphere is not manifold anymore, further processing might be inaccurate");
    for (size_t i = 0; i < ps->indices.size(); i++) {
      auto& tri = ps->indices[i];
      if (tri[0] == tri[1] || tri[0] == tri[2] || tri[1] == tri[2]) continue;
      for (int j = 0; j < 3; j++) {
        int vi1 = tri[j];
        int vi2 = tri[(j + 1) % 3];
        double l1 = (ps->vertices[vi1] - ps->vertices[vi2]).norm();
        EdgeKey ek(tri[j], tri[(j + 1) % 3]);
        if (edge_db.count(ek) != 0) {
          auto ev = edge_db.at(ek);
          int face_o, pos_o;
          if (vi2 > vi1) {
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
              if (tri_oth[k] == vi1) tri_oth_[k] = tri[(j + 2) % 3];
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
  DECLARE_INSTANCE();

  char *kwlist[] = {"r", "d", NULL};
  double r = NAN;
  PyObject *rp = nullptr;
  double d = NAN;

  double vr = 1;

  auto discretizer = CreateCurveDiscretizer(kwargs);
  auto node = std::make_shared<SphereNode>(instance, discretizer);
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|Od", kwlist, &rp, &d)) {
    PyErr_SetString(PyExc_TypeError, "Error during parsing sphere(r|d)");
    return NULL;
  }
  if (rp != nullptr) {
    if (python_numberval(rp, &r, &(node->dragflags), 1))
      if (rp->ob_type == &PyFunction_Type) node->r_func = rp;
  }
  if (!isnan(r)) {
    if (r <= 0) {
      PyErr_SetString(PyExc_TypeError, "Parameter r must be positive");
      return NULL;
    }
    vr = r;
    if (!isnan(d)) {
      PyErr_SetString(PyExc_TypeError, "Cant specify r and d at the same time for sphere");
      return NULL;
    }
  }
  if (!isnan(d)) {
    if (d <= 0) {
      PyErr_SetString(PyExc_TypeError, "Parameter d must be positive");
      return NULL;
    }
    vr = d / 2.0;
  }

  node->r = vr;

  python_retrieve_pyname(node);
  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

PyObject *python_cylinder(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE();

  char *kwlist[] = {"h", "r1", "r2", "center", "r", "d", "d1", "d2", "angle", NULL};
  PyObject *h_ = nullptr;
  PyObject *r_ = nullptr;
  double r1 = NAN;
  double r2 = NAN;
  double d = NAN;
  double d1 = NAN;
  double d2 = NAN;
  double angle = NAN;

  PyObject *center = NULL;
  double vr1 = 1, vr2 = 1, vh = 1;

  auto discretizer = CreateCurveDiscretizer(kwargs);
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|OddOOdddd", kwlist, &h_, &r1, &r2, &center, &r_, &d,
                                   &d1, &d2, &angle)) {
    PyErr_SetString(PyExc_TypeError, "Error during parsing cylinder(h,r|r1+r2|d1+d2)");
    return NULL;
  }
  double r = NAN;
  double h = NAN;

  auto node = std::make_shared<CylinderNode>(instance, discretizer);
  python_numberval(r_, &r, &(node->dragflags), 1);
  python_numberval(h_, &h, &(node->dragflags), 2);

  if (h <= 0) {
    PyErr_SetString(PyExc_TypeError, "Cylinder height must be positive");
    return NULL;
  }
  vh = h;

  if (!isnan(d) && d <= 0) {
    PyErr_SetString(PyExc_TypeError, "Cylinder d must be positive");
    return NULL;
  }
  if (!isnan(r1) && r1 < 0) {
    PyErr_SetString(PyExc_TypeError, "Cylinder r1 must not be negative");
    return NULL;
  }
  if (!isnan(r2) && r2 < 0) {
    PyErr_SetString(PyExc_TypeError, "Cylinder r2 must not be negative");
    return NULL;
  }
  if (!isnan(d1) && d1 < 0) {
    PyErr_SetString(PyExc_TypeError, "Cylinder d1 must not be negative");
    return NULL;
  }
  if (!isnan(d2) && d2 < 0) {
    PyErr_SetString(PyExc_TypeError, "Cylinder d2 must not be negative");
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

  node->r1 = vr1;
  node->r2 = vr2;
  node->h = vh;

  if (center == Py_True) node->center = 1;
  else if (center == Py_False || center == NULL) node->center = 0;
  else {
    PyErr_SetString(PyExc_TypeError, "Unknown Value for center parameter");
    return NULL;
  }

  python_retrieve_pyname(node);
  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

PyObject *python_polyhedron(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE();
  auto node = std::make_shared<PolyhedronNode>(instance);

  char *kwlist[] = {"points", "faces", "convexity", "triangles", "colors", NULL};
  PyObject *points = NULL;
  PyObject *faces = NULL;
  int convexity = 2;
  PyObject *triangles = NULL;
  PyObject *colors = NULL;

  PyObject *element;
  Vector3d point;

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!|iO!O!", kwlist, &PyList_Type, &points,
                                   &PyList_Type, &faces, &convexity, &PyList_Type, &triangles,
                                   &PyList_Type, &colors)) {
    PyErr_SetString(PyExc_TypeError,
                    "Error during parsing polyhedron(points, faces, convexity, triangles, colors)");
    return NULL;
  }

  if (points != NULL && PyList_Check(points)) {
    if (PyList_Size(points) == 0) {
      PyErr_SetString(PyExc_TypeError, "There must at least be one point in the polyhedron");
      return NULL;
    }
    for (Py_ssize_t i = 0; i < PyList_Size(points); i++) {
      element = PyList_GetItem(points, i);
      if (PyList_Check(element) && PyList_Size(element) == 3) {
        point[0] = PyFloat_AsDouble(PyList_GetItem(element, 0));
        point[1] = PyFloat_AsDouble(PyList_GetItem(element, 1));
        point[2] = PyFloat_AsDouble(PyList_GetItem(element, 2));
        node->points.push_back(point);
      } else {
        PyErr_SetString(PyExc_TypeError, "Coordinate must exactly contain 3 numbers");
        return NULL;
      }
    }
  } else {
    PyErr_SetString(PyExc_TypeError, "Polyhedron Points must be a list of coordinates");
    return NULL;
  }

  if (triangles != NULL) {
    faces = triangles;
    //	LOG(message_group::Deprecated, inst->location(), parameters.documentRoot(),
    //"polyhedron(triangles=[]) will be removed in future releases. Use polyhedron(faces=[]) instead.");
  }

  if (faces != NULL && PyList_Check(faces)) {
    if (PyList_Size(faces) == 0) {
      PyErr_SetString(PyExc_TypeError, "must specify at least 1 face");
      return NULL;
    }
    for (Py_ssize_t i = 0; i < PyList_Size(faces); i++) {
      element = PyList_GetItem(faces, i);
      if (PyList_Check(element)) {
        IndexedFace face;
        for (Py_ssize_t j = 0; j < PyList_Size(element); j++) {
          long pointIndex = PyLong_AsLong(PyList_GetItem(element, j));
          if (pointIndex < 0 || pointIndex >= static_cast<long>(node->points.size())) {
            PyErr_SetString(PyExc_TypeError, "Polyhedron Point Index out of range");
            return NULL;
          }
          face.push_back(pointIndex);
        }
        if (face.size() >= 3) {
          node->faces.push_back(std::move(face));
        } else {
          PyErr_SetString(PyExc_TypeError, "Polyhedron Face must sepcify at least 3 indices");
          return NULL;
        }

      } else {
        PyErr_SetString(PyExc_TypeError, "Polyhedron Face must be a list of indices");
        return NULL;
      }
    }
  } else {
    PyErr_SetString(PyExc_TypeError, "Polyhedron faces must be a list of indices");
    return NULL;
  }

  if (colors != NULL && PyList_Check(colors)) {
    if ((size_t)PyList_Size(colors) != node->faces.size()) {
      PyErr_SetString(PyExc_TypeError, "when specified must match number of faces");
      return NULL;
    }
    for (Py_ssize_t i = 0; i < PyList_Size(colors); i++) {
      element = PyList_GetItem(colors, i);
      if (PyList_Check(element) && PyList_Size(element) == 3) {
        Vector4f color(0, 0, 0, 1.0);
        for (int j = 0; j < 3; j++) {
          color[j] = PyFloat_AsDouble(PyList_GetItem(element, j));
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
        PyErr_SetString(PyExc_TypeError, "Face Color must be a list with 3 values");
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
  DECLARE_INSTANCE();
  auto node = std::make_shared<FrepNode>(instance);
  PyObject *expression = NULL;
  PyObject *bmin = NULL, *bmax = NULL;
  double res = 10;

  char *kwlist[] = {"exp", "min", "max", "res", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOO|d", kwlist, &expression, &bmin, &bmax, &res))
    return NULL;

  python_vectorval(bmin, 3, 3, &(node->x1), &(node->y1), &(node->z1));
  python_vectorval(bmax, 3, 3, &(node->x2), &(node->y2), &(node->z2));
  node->res = res;

  if (expression->ob_type == &PyDataType) {
    node->expression = expression;
  } else if (expression->ob_type == &PyFunction_Type) {
    node->expression = expression;
  } else {
    PyErr_SetString(PyExc_TypeError, "Unknown frep expression type\n");
    return NULL;
  }

  python_retrieve_pyname(node);
  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

PyObject *python_ifrep(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE();
  PyObject *object = NULL;
  PyObject *dummydict = nullptr;

  char *kwlist[] = {"obj", nullptr};
  std::shared_ptr<AbstractNode> child;

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!", kwlist, &PyOpenSCADType, &object)) return NULL;

  child = PyOpenSCADObjectToNodeMulti(object, &dummydict);
  auto dummydict_owner = py_owned(dummydict);
  // Two failure modes the original C-style cast quietly missed:
  //   1) child == nullptr (PyOpenSCADObjectToNodeMulti failed,
  //      possibly with a Python exception already set);
  //   2) child is non-null but not actually a LeafNode -- e.g. the
  //      caller passed a CSG operation like `union(a, b)` whose
  //      runtime type is CsgOpNode. The C-style cast would happily
  //      pretend the bytes are a LeafNode and the subsequent
  //      `createGeometry()` call would walk uninitialized vtable
  //      slots, i.e. undefined behaviour.
  // Validate both with a dynamic_cast (AbstractNode is polymorphic
  // via BaseVisitable) and surface a Python-visible error instead.
  if (child == nullptr) return propagate_or_typeerror("Invalid type for Object in ifrep");
  LeafNode *node = dynamic_cast<LeafNode *>(child.get());
  if (node == nullptr) {
    PyErr_SetString(PyExc_TypeError,
                    "ifrep() expects a primitive object (cube, sphere, polyhedron, ...); "
                    "for CSG combinations use frep() with an explicit expression");
    return nullptr;
  }
  const std::shared_ptr<const Geometry> geom = node->createGeometry();
  const std::shared_ptr<const PolySet> ps = std::dynamic_pointer_cast<const PolySet>(geom);

  return ifrep(ps);
}

#endif

PyObject *python_square(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE();
  auto node = std::make_shared<SquareNode>(instance);

  char *kwlist[] = {"dim", "center", NULL};
  PyObject *dim = NULL;

  PyObject *center = NULL;
  double z = NAN;

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|OO", kwlist, &dim, &center)) {
    PyErr_SetString(PyExc_TypeError, "Error during parsing square(dim)");
    return NULL;
  }
  if (dim != NULL) {
    if (python_vectorval(dim, 2, 2, &(node->x), &(node->y), &z)) {
      PyErr_SetString(PyExc_TypeError, "Invalid Square dimensions");
      return NULL;
    }
  }
  if (center == Py_True) node->center = 1;
  else if (center == Py_False || center == NULL) node->center = 0;
  else {
    PyErr_SetString(PyExc_TypeError, "Unknown Value for center parameter");
    return NULL;
  }

  python_retrieve_pyname(node);
  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

PyObject *python_circle(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE();

  char *kwlist[] = {"r", "d", "angle", NULL};
  double r = NAN;
  double d = NAN;
  double angle = NAN;
  double vr = 1;

  auto discretizer = CreateCurveDiscretizer(kwargs);
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|ddd", kwlist, &r, &d, &angle)) {
    PyErr_SetString(PyExc_TypeError, "Error during parsing circle(r|d)");
    return NULL;
  }

  if (!isnan(r)) {
    if (r <= 0) {
      PyErr_SetString(PyExc_TypeError, "Parameter r must be positive");
      return NULL;
    }
    vr = r;
    if (!isnan(d)) {
      PyErr_SetString(PyExc_TypeError, "Cant specify r and d at the same time for circle");
      return NULL;
    }
  }
  if (!isnan(d)) {
    if (d <= 0) {
      PyErr_SetString(PyExc_TypeError, "Parameter d must be positive");
      return NULL;
    }
    vr = d / 2.0;
  }

  auto node = std::make_shared<CircleNode>(instance, discretizer);

  if (!isnan(angle)) node->angle = angle;

  node->r = vr;

  python_retrieve_pyname(node);
  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

PyObject *python_polygon(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE();
  auto node = std::make_shared<PolygonNode>(instance, CreateCurveDiscretizer(kwargs));

  char *kwlist[] = {"points", "paths", "convexity", NULL};
  PyObject *pypoints = NULL;
  PyObject *pypaths = NULL;
  int convexity = 2;

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!|O!i", kwlist, &PyList_Type, &pypoints, &PyList_Type,
                                   &pypaths, &convexity)) {
    PyErr_SetString(PyExc_TypeError, "Error during parsing polygon(points,paths)");
    return NULL;
  }

  std::vector<Vector3d> points = python_to2dvarpointlist(pypoints);
  if (points.size() == 0) {
    return NULL;
  }

  node->points = points;

  std::vector<std::vector<size_t>> paths = python_to2dintlist(pypaths);
  //  if (paths.size() == 0) { TODO fehlercode ?
  //    return NULL;
  //  }
  node->paths = paths;

  node->convexity = convexity;
  if (node->convexity < 1) node->convexity = 1;

  python_retrieve_pyname(node);
  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

PyObject *python_polyline(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE();
  auto node = std::make_shared<PolylineNode>(instance, CreateCurveDiscretizer(kwargs));

  char *kwlist[] = {"points", NULL};
  PyObject *pypoints = NULL;

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!", kwlist, &PyList_Type, &pypoints)) {
    PyErr_SetString(PyExc_TypeError, "Error during parsing polyline(points,paths)");
    return NULL;
  }

  std::vector<Vector3d> points = python_to2dvarpointlist(pypoints);
  if (points.size() == 0) {
    return NULL;
  }

  node->points = points;

  python_retrieve_pyname(node);
  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

PyObject *python_spline(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE();
  auto node = std::make_shared<SplineNode>(instance);

  char *kwlist[] = {"points", "fn", "fa", "fs", NULL};
  PyObject *points = NULL;
  double fn = 0, fa = 0, fs = 0;

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!|ddd", kwlist, &PyList_Type, &points, &fn, &fa,
                                   &fs)) {
    PyErr_SetString(PyExc_TypeError, "Error during parsing spline(points)");
    return NULL;
  }

  if (points != NULL && PyList_Check(points)) {
    if (PyList_Size(points) == 0) {
      PyErr_SetString(PyExc_TypeError, "There must at least be one point in the polygon");
      return NULL;
    }
    for (Py_ssize_t i = 0; i < PyList_Size(points); i++) {
      PyObject *element = PyList_GetItem(points, i);
      if (PyList_Check(element) && PyList_Size(element) == 2) {
        Vector2d point;
        point[0] = PyFloat_AsDouble(PyList_GetItem(element, 0));
        point[1] = PyFloat_AsDouble(PyList_GetItem(element, 1));
        node->points.push_back(point);
      } else {
        PyErr_SetString(PyExc_TypeError, "Coordinate must exactly contain 2 numbers");
        return NULL;
      }
    }
  } else {
    PyErr_SetString(PyExc_TypeError, "Polygon points must be a list of coordinates");
    return NULL;
  }
  node->fn = fn;
  node->fa = fa;
  node->fs = fs;

  python_retrieve_pyname(node);
  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

std::optional<std::string> to_optional_string(const char *ptr)
{
  if (ptr != nullptr) {
    return std::string(ptr);
  }
  return {};
}

PyObject *python_text(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE();
  char *kwlist[] = {"text",     "size",   "font",   "spacing", "direction",
                    "language", "script", "halign", "valign",  NULL};

  double size = 1.0, spacing = 1.0;

  const char *text = "", *font = NULL, *direction = "ltr", *language = "en", *script = "latin",
             *valign = "baseline", *halign = "left";

  auto discretizer = CreateCurveDiscretizer(kwargs);
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s|dsdsssss", kwlist, &text, &size, &font, &spacing,
                                   &direction, &language, &script, &halign, &valign)) {
    PyErr_SetString(PyExc_TypeError, "Error during parsing text(string, ...))");
    return NULL;
  }

  auto node = std::make_shared<TextNode>(
    instance, FreetypeRenderer::Params(FreetypeRenderer::Params::ParamsOptions{
                .curve_discretizer = std::make_shared<CurveDiscretizer>(discretizer),
                .size = size,
                .spacing = spacing,
                .text = to_optional_string(text),
                .font = to_optional_string(font),
                .direction = to_optional_string(direction),
                .language = to_optional_string(language),
                .script = to_optional_string(script),
                .halign = to_optional_string(halign),
                .valign = to_optional_string(valign),
                .loc = instance->location(),
              }));

  node->params.detect_properties();

  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

PyObject *python_textmetrics(PyObject *self, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE();
  char *kwlist[] = {"text",     "size",   "font",   "spacing", "direction",
                    "language", "script", "halign", "valign",  NULL};

  double size = 1.0, spacing = 1.0;

  const char *text = "", *font = NULL, *direction = "ltr", *language = "en", *script = "latin",
             *valign = "baseline", *halign = "left";

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s|dsdsssss", kwlist, &text, &size, &font, &spacing,
                                   &direction, &language, &script, &valign, &halign)) {
    PyErr_SetString(PyExc_TypeError, "Error during parsing textmetrics");
    return NULL;
  }

  FreetypeRenderer::Params ftparams(FreetypeRenderer::Params::ParamsOptions{
    .curve_discretizer = {},
    .size = size,
    .spacing = spacing,
    .text = to_optional_string(text),
    .font = to_optional_string(font),
    .direction = to_optional_string(direction),
    .language = to_optional_string(language),
    .script = to_optional_string(script),
    .halign = to_optional_string(halign),
    .valign = to_optional_string(valign),
    .loc = instance->location(),
  });

  FreetypeRenderer::TextMetrics metrics(ftparams);
  if (!metrics.ok) {
    PyErr_SetString(PyExc_TypeError, "Invalid Metric");
    return NULL;
  }
  PyObject *offset = PyList_New(2);
  PyList_SetItem(offset, 0, PyFloat_FromDouble(metrics.x_offset));
  PyList_SetItem(offset, 1, PyFloat_FromDouble(metrics.y_offset));

  PyObject *advance = PyList_New(2);
  PyList_SetItem(advance, 0, PyFloat_FromDouble(metrics.advance_x));
  PyList_SetItem(advance, 1, PyFloat_FromDouble(metrics.advance_y));

  PyObject *position = PyList_New(2);
  PyList_SetItem(position, 0, PyFloat_FromDouble(metrics.bbox_x));
  PyList_SetItem(position, 1, PyFloat_FromDouble(metrics.bbox_y));

  PyObject *dims = PyList_New(2);
  PyList_SetItem(dims, 0, PyFloat_FromDouble(metrics.bbox_w));
  PyList_SetItem(dims, 1, PyFloat_FromDouble(metrics.bbox_h));

  PyObject *dict;
  dict = PyDict_New();
  PyDict_SetItemString(dict, "ascent", PyFloat_FromDouble(metrics.ascent));
  PyDict_SetItemString(dict, "descent", PyFloat_FromDouble(metrics.descent));
  PyDict_SetItemString(dict, "offset", offset);
  PyDict_SetItemString(dict, "advance", advance);
  PyDict_SetItemString(dict, "position", position);
  PyDict_SetItemString(dict, "size", dims);
  return (PyObject *)dict;
}

int sheetCalcIndInt(PyObject *func, double i, double j, Vector3d& pos)
{
  PyObject *args = PyTuple_Pack(2, PyFloat_FromDouble(i), PyFloat_FromDouble(j));
  PyObject *pos_p = PyObject_CallObject(func, args);
  if (pos_p == nullptr) {
    std::string errorstr;
    python_catch_error(errorstr);
    PyErr_SetString(PyExc_TypeError, errorstr.c_str());
    LOG(message_group::Error, errorstr.c_str());
    return 1;
  }
  return python_vectorval(pos_p, 3, 3, &pos[0], &pos[1], &pos[2], nullptr, nullptr);
}

int sheetCalcInd(PolySetBuilder& builder, std::vector<Vector3d>& vertices, std::vector<double>& istore,
                 std::vector<double>& jstore, PyObject *func, double i, double j)
{
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
  double dist, ang;
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
  DECLARE_INSTANCE();
  auto node = std::make_shared<SheetNode>(instance);
  // TODO check type of func
  node->func = (void *)func;
  node->imin = imin;
  node->imax = imax;
  node->jmin = jmin;
  node->jmax = jmax;
  node->fs = fs;
  node->ispan = (ispan == Py_True) ? true : false;
  node->jspan = (jspan == Py_True) ? true : false;

  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

PyObject *python_sheet(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"func", "imin", "imax", "jmin", "jmax", "fs", "iclose", "jclose", NULL};
  PyObject *func = NULL;
  double imin, imax, jmin, jmax;
  PyObject *ispan = nullptr, *jspan = nullptr;
  double fs = 1.0;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Odddd|dOO", kwlist, &func, &imin, &imax, &jmin, &jmax,
                                   &fs, &ispan, &jspan)) {
    PyErr_SetString(PyExc_TypeError, "Error during parsing sheet(func, imin, imax, jmin, jmax)");
    return NULL;
  }
  if (func->ob_type != &PyFunction_Type) {
    PyErr_SetString(PyExc_TypeError, "must specify a function");
    return NULL;
  }

  return python_sheet_core(func, imin, imax, jmin, jmax, fs, ispan, jspan);
}

// organic_mesh.cpp
// -----------------------------------------------------------------------
// Eigenstaendige C++-Funktionen (keine CPython-Bindung noetig):
//
//   PolySet organic_resample(const std::vector<Vector3d>& points,
//                             double max_mesh_size);
//
// Ablauf:
//   1) delaunay2d_bowyer_watson()   - klassische 2D-Delaunay-Triangulierung
//                                     (Bowyer-Watson, inkrementell)
//   2) project_points_to_plane()    - PCA-Projektion der 3D-Punkte auf ihre
//                                     Hauptebene, damit (1) angewendet werden
//                                     kann (Topologie kommt aus der Ebene,
//                                     die Positionen bleiben die echten 3D-
//                                     Koordinaten)
//   3) compute_vertex_normals()     - flaechengewichtete Eckennormalen
//   4) compute_pn_patch()/eval_pn() - PN-Triangle (Vlachos et al. 2001):
//                                     kubisches Bezier-Dreieck, dessen
//                                     Randkurven nur von Position+Normale
//                                     der beiden Eckpunkte abhaengen ->
//                                     Normalensprung an gemeinsamen Kanten
//                                     verschwindet.
//   5) tessellate()                 - baryzentrisches Resampling aller
//                                     Patches mit global einheitlicher
//                                     Unterteilungstiefe, so gewaehlt, dass
//                                     keine Kante > max_mesh_size lang ist.
//                                     Kantenpunkte werden dedupliziert, das
//                                     Ergebnis ist wasserdicht.
//
// Abhaengigkeit: nur Eigen (Vector2d/Vector3d), keine weiteren Libraries.
//
// Hinweis zur Integration in PythonSCAD/OpenSCAD:
//   Falls ihr die echte OpenSCAD-Klasse PolySet (src/geometry/PolySet.h)
//   statt der hier definierten Miniversion verwenden wollt, ersetzt die
//   Struct-Definition unten durch ein #include "PolySet.h" und passt
//   build_polyset() an eure PolySet-Konstruktion an (append_vertex /
//   append_poly bzw. PolySetBuilder, je nach Version). Die eigentliche
//   Geometrie-Berechnung (Schritte 1-5) bleibt unveraendert.
// -----------------------------------------------------------------------

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <vector>
#include <array>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cassert>

using Vector2d = Eigen::Vector2d;
using Vector3d = Eigen::Vector3d;

// =========================================================================
// 1) 2D-Delaunay-Triangulierung (Bowyer-Watson, inkrementell)
// =========================================================================
namespace detail {

struct Tri2 {
  int a, b, c;
};

// Umkreis-Test: liegt Punkt p im Umkreis des Dreiecks (a,b,c)?
static bool in_circumcircle(const Vector2d& a, const Vector2d& b, const Vector2d& c, const Vector2d& p)
{
  const double ax = a.x() - p.x(), ay = a.y() - p.y();
  const double bx = b.x() - p.x(), by = b.y() - p.y();
  const double cx = c.x() - p.x(), cy = c.y() - p.y();

  const double det = (ax * ax + ay * ay) * (bx * cy - cx * by) -
                     (bx * bx + by * by) * (ax * cy - cx * ay) +
                     (cx * cx + cy * cy) * (ax * by - bx * ay);

  // Orientierung von (a,b,c) beruecksichtigen (Vorzeichenkonvention)
  const double orient = (b.x() - a.x()) * (c.y() - a.y()) - (c.x() - a.x()) * (b.y() - a.y());
  return orient > 0 ? (det > 0) : (det < 0);
}

struct EdgeKey {
  int a, b;  // a < b
  bool operator==(const EdgeKey& o) const { return a == o.a && b == o.b; }
};
struct EdgeKeyHash {
  size_t operator()(const EdgeKey& e) const
  {
    return (static_cast<size_t>(e.a) << 32) ^ static_cast<size_t>(e.b);
  }
};

// Liefert Dreiecke als Indizes in 'points' (points enthaelt NUR die
// Eingabepunkte; das Super-Dreieck wird intern separat verwaltet und am
// Ende herausgefiltert).
inline std::vector<Tri2> delaunay2d_bowyer_watson(const std::vector<Vector2d>& points)
{
  const int n = static_cast<int>(points.size());
  std::vector<Vector2d> pts = points;

  // Super-Dreieck, das alle Punkte weit umschliesst
  double minx = 1e300, miny = 1e300, maxx = -1e300, maxy = -1e300;
  for (const auto& p : points) {
    minx = std::min(minx, p.x());
    maxx = std::max(maxx, p.x());
    miny = std::min(miny, p.y());
    maxy = std::max(maxy, p.y());
  }
  const double dx = maxx - minx, dy = maxy - miny;
  const double delta = std::max(dx, dy) * 10.0 + 1.0;
  const double midx = (minx + maxx) * 0.5, midy = (miny + maxy) * 0.5;

  const int s0 = n, s1 = n + 1, s2 = n + 2;
  pts.push_back(Vector2d(midx - 20 * delta, midy - delta));
  pts.push_back(Vector2d(midx, midy + 20 * delta));
  pts.push_back(Vector2d(midx + 20 * delta, midy - delta));

  std::vector<Tri2> tris;
  tris.push_back({s0, s1, s2});

  for (int pi = 0; pi < n; ++pi) {
    const Vector2d& p = pts[pi];

    // 1) "schlechte" Dreiecke finden (Umkreis enthaelt p)
    std::vector<int> bad;
    bad.reserve(8);
    for (int t = 0; t < static_cast<int>(tris.size()); ++t) {
      const Tri2& tr = tris[t];
      if (in_circumcircle(pts[tr.a], pts[tr.b], pts[tr.c], p)) bad.push_back(t);
    }

    // 2) Rand-Polygon = Kanten, die genau EINEM schlechten Dreieck angehoeren
    std::unordered_map<EdgeKey, int, EdgeKeyHash> edge_count;
    auto add_edge = [&](int a, int b) {
      EdgeKey k{std::min(a, b), std::max(a, b)};
      edge_count[k]++;
    };
    for (int t : bad) {
      add_edge(tris[t].a, tris[t].b);
      add_edge(tris[t].b, tris[t].c);
      add_edge(tris[t].c, tris[t].a);
    }

    std::vector<std::array<int, 2>> boundary;
    for (int t : bad) {
      const Tri2& tr = tris[t];
      int es[3][2] = {{tr.a, tr.b}, {tr.b, tr.c}, {tr.c, tr.a}};
      for (auto& e : es) {
        EdgeKey k{std::min(e[0], e[1]), std::max(e[0], e[1])};
        if (edge_count[k] == 1) boundary.push_back({e[0], e[1]});
      }
    }

    // 3) schlechte Dreiecke entfernen (rueckwaerts, damit Indizes stabil bleiben)
    std::sort(bad.rbegin(), bad.rend());
    for (int t : bad) tris.erase(tris.begin() + t);

    // 4) neue Dreiecke aus Rand-Polygon + neuem Punkt
    for (auto& e : boundary) tris.push_back({e[0], e[1], pi});
  }

  // Dreiecke entfernen, die eine Ecke des Super-Dreiecks benutzen
  std::vector<Tri2> result;
  result.reserve(tris.size());
  for (const auto& t : tris) {
    if (t.a >= n || t.b >= n || t.c >= n) continue;
    result.push_back(t);
  }
  return result;
}

}  // namespace detail

// =========================================================================
// 2) PCA-Projektion der 3D-Punkte auf ihre Hauptebene (fuer die Topologie)
// =========================================================================
inline std::vector<Vector2d> project_points_to_plane(const std::vector<Vector3d>& pts)
{
  Vector3d centroid = Vector3d::Zero();
  for (const auto& p : pts) centroid += p;
  centroid /= static_cast<double>(pts.size());

  Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
  for (const auto& p : pts) {
    Vector3d d = p - centroid;
    cov += d * d.transpose();
  }

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(cov);
  // Eigenwerte aufsteigend sortiert -> die zwei groessten sind die letzten Spalten
  Vector3d axis1 = solver.eigenvectors().col(2);
  Vector3d axis2 = solver.eigenvectors().col(1);

  std::vector<Vector2d> uv;
  uv.reserve(pts.size());
  for (const auto& p : pts) {
    Vector3d d = p - centroid;
    uv.emplace_back(d.dot(axis1), d.dot(axis2));
  }
  return uv;
}

// =========================================================================
// 3) Eckennormalen (flaechengewichtet)
// =========================================================================
inline std::vector<Vector3d> compute_vertex_normals(const std::vector<Vector3d>& pts,
                                                    const std::vector<std::array<int, 3>>& tris)
{
  std::vector<Vector3d> normals(pts.size(), Vector3d::Zero());
  for (const auto& t : tris) {
    const Vector3d &a = pts[t[0]], &b = pts[t[1]], &c = pts[t[2]];
    Vector3d n = (b - a).cross(c - a);  // Laenge ~ 2x Dreiecksflaeche
    normals[t[0]] += n;
    normals[t[1]] += n;
    normals[t[2]] += n;
  }
  for (auto& n : normals) {
    double len = n.norm();
    if (len > 1e-12) n /= len;
  }
  return normals;
}

// =========================================================================
// 4) PN-Triangle: Kontrollpunkte + Auswertung
// =========================================================================
struct PNPatch {
  Vector3d b300, b030, b003, b210, b120, b021, b012, b102, b201, b111;
};

// Max. Anteil der Kantenlaenge, um den ein Kanten-Kontrollpunkt von der
// Verbindungsgeraden abweichen darf. Schuetzt vor Ausreissern durch
// instabile Normalen an duennen/entarteten Randdreiecken (haeufig bei
// spaerlich abgetasteten Punktwolken).
static constexpr double kMaxOffsetFraction = 0.5;

inline PNPatch compute_pn_patch(const Vector3d& p1, const Vector3d& p2, const Vector3d& p3,
                                const Vector3d& n1, const Vector3d& n2, const Vector3d& n3)
{
  PNPatch cp;
  cp.b300 = p1;
  cp.b030 = p2;
  cp.b003 = p3;

  // Liefert den Kanten-Kontrollpunkt b_ij = (2*Pi+Pj - w*Ni)/3, wobei die
  // Korrektur -w*Ni/3 auf kMaxOffsetFraction * |Pj-Pi| begrenzt wird.
  auto edge_ctrl = [](const Vector3d& pi, const Vector3d& pj, const Vector3d& ni) {
    double w = (pj - pi).dot(ni);
    Vector3d correction = -w * ni / 3.0;
    double max_len = kMaxOffsetFraction * (pj - pi).norm();
    double len = correction.norm();
    if (len > max_len && len > 1e-12) correction *= (max_len / len);
    return (2 * pi + pj) / 3.0 + correction;
  };

  cp.b210 = edge_ctrl(p1, p2, n1);
  cp.b120 = edge_ctrl(p2, p1, n2);
  cp.b021 = edge_ctrl(p2, p3, n2);
  cp.b012 = edge_ctrl(p3, p2, n3);
  cp.b102 = edge_ctrl(p3, p1, n3);
  cp.b201 = edge_ctrl(p1, p3, n1);

  Vector3d E = (cp.b210 + cp.b120 + cp.b021 + cp.b012 + cp.b102 + cp.b201) / 6.0;
  Vector3d V = (p1 + p2 + p3) / 3.0;
  cp.b111 = E + (E - V) * 0.5;
  return cp;
}

inline Vector3d eval_pn(const PNPatch& cp, double u, double v, double w)
{
  const double u2 = u * u, v2 = v * v, w2 = w * w;
  const double u3 = u2 * u, v3 = v2 * v, w3 = w2 * w;
  return cp.b300 * u3 + cp.b030 * v3 + cp.b003 * w3 + cp.b210 * (3 * u2 * v) + cp.b120 * (3 * u * v2) +
         cp.b021 * (3 * v2 * w) + cp.b012 * (3 * v * w2) + cp.b102 * (3 * u * w2) +
         cp.b201 * (3 * u2 * w) + cp.b111 * (6 * u * v * w);
}

// =========================================================================
// 5) Subdiv-Tiefe aus max_mesh_size ableiten + baryzentrisches Resampling
// =========================================================================
namespace detail {

inline int subdiv_from_max_edge(const std::vector<Vector3d>& pts,
                                const std::vector<std::array<int, 3>>& tris, double max_mesh_size)
{
  double max_edge = 0.0;
  for (const auto& t : tris) {
    max_edge = std::max(max_edge, (pts[t[0]] - pts[t[1]]).norm());
    max_edge = std::max(max_edge, (pts[t[1]] - pts[t[2]]).norm());
    max_edge = std::max(max_edge, (pts[t[2]] - pts[t[0]]).norm());
  }
  if (max_mesh_size <= 0.0 || max_edge <= max_mesh_size) return 1;
  int subdiv = static_cast<int>(std::ceil(max_edge / max_mesh_size));
  return std::max(1, std::min(subdiv, 64));  // Obergrenze gegen Explosion
}

// baryzentrisches Gitter (u,v,w) fuer gegebene Unterteilungstiefe n,
// inkl. Kantenmetadaten fuer die Dedup-Logik.
struct Grid {
  std::vector<std::array<double, 3>> coords;  // (u,v,w)
  std::vector<std::array<int, 3>> faces;      // lokale Indizes
  std::vector<int> edge_type;                 // 0=innen,1=ab,2=bc,3=ca
  std::vector<int> edge_param;
};

inline Grid build_grid(int n)
{
  Grid g;
  std::vector<std::vector<int>> index(n + 1, std::vector<int>(n + 1, -1));

  for (int i = 0; i <= n; ++i) {
    for (int j = 0; j <= n - i; ++j) {
      int k = n - i - j;
      index[i][j] = static_cast<int>(g.coords.size());
      g.coords.push_back({double(i) / n, double(j) / n, double(k) / n});
      if (k == 0) {
        g.edge_type.push_back(1);
        g.edge_param.push_back(i);
      } else if (i == 0) {
        g.edge_type.push_back(2);
        g.edge_param.push_back(j);
      } else if (j == 0) {
        g.edge_type.push_back(3);
        g.edge_param.push_back(i);
      } else {
        g.edge_type.push_back(0);
        g.edge_param.push_back(0);
      }
    }
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n - i; ++j) {
      int k = n - i - j;
      int a = index[i][j], b = index[i + 1][j], c = index[i][j + 1];
      g.faces.push_back({a, b, c});
      if (k - 1 > 0) {
        int d = index[i + 1][j + 1];
        g.faces.push_back({b, d, c});
      }
    }
  }
  return g;
}

struct EdgeCacheKey {
  int lo, hi, p;
  bool operator==(const EdgeCacheKey& o) const { return lo == o.lo && hi == o.hi && p == o.p; }
};
struct EdgeCacheKeyHash {
  size_t operator()(const EdgeCacheKey& k) const
  {
    return (static_cast<size_t>(k.lo) << 42) ^ (static_cast<size_t>(k.hi) << 21) ^
           static_cast<size_t>(k.p);
  }
};

}  // namespace detail

// =========================================================================
// Hauptfunktion
// =========================================================================
PolySet organic_resample(const std::vector<Vector3d>& points, double max_mesh_size)
{
  PolySet out(3);
  if (points.size() < 3) return out;

  // --- 1+2) Delaunay-Topologie ueber PCA-Projektion ---
  std::vector<Vector2d> uv = project_points_to_plane(points);
  std::vector<detail::Tri2> tris2d = detail::delaunay2d_bowyer_watson(uv);

  std::vector<std::array<int, 3>> tris;
  tris.reserve(tris2d.size());
  for (const auto& t : tris2d) tris.push_back({t.a, t.b, t.c});

  // --- 3) Eckennormalen ---
  std::vector<Vector3d> normals = compute_vertex_normals(points, tris);

  // --- 4+5) subdiv-Tiefe bestimmen, Gitter einmalig bauen ---
  int subdiv = detail::subdiv_from_max_edge(points, tris, max_mesh_size);
  detail::Grid grid = detail::build_grid(subdiv);

  std::unordered_map<detail::EdgeCacheKey, int, detail::EdgeCacheKeyHash> edge_cache;
  edge_cache.reserve(tris.size() * subdiv * 3);

  std::vector<int> local_index(grid.coords.size());

  for (const auto& t : tris) {
    const int ia = t[0], ib = t[1], ic = t[2];
    PNPatch cp =
      compute_pn_patch(points[ia], points[ib], points[ic], normals[ia], normals[ib], normals[ic]);

    for (size_t gi = 0; gi < grid.coords.size(); ++gi) {
      const auto& uvw = grid.coords[gi];
      const int et = grid.edge_type[gi];

      if (et != 0) {
        int lo, hi, p;
        const int param = grid.edge_param[gi];
        // Wichtig: die lokalen Parameter i/j messen je Kantentyp aus
        // unterschiedlicher Richtung (siehe build_grid). Damit der
        // globale Schluessel (lo,hi,p) unabhaengig vom Dreieck
        // konsistent "Abstand ab lo" bedeutet, muss die Richtung
        // pro Kantentyp passend gespiegelt werden.
        if (et == 1) {
          lo = std::min(ia, ib);
          hi = std::max(ia, ib);
          p = (ia < ib) ? subdiv - param : param;
        } else if (et == 2) {
          lo = std::min(ib, ic);
          hi = std::max(ib, ic);
          p = (ib < ic) ? subdiv - param : param;
        } else {
          lo = std::min(ic, ia);
          hi = std::max(ic, ia);
          p = (ic < ia) ? param : subdiv - param;
        }

        detail::EdgeCacheKey key{lo, hi, p};
        auto it = edge_cache.find(key);
        if (it != edge_cache.end()) {
          local_index[gi] = it->second;
          continue;
        }
        Vector3d pos = eval_pn(cp, uvw[0], uvw[1], uvw[2]);
        int idx = static_cast<int>(out.vertices.size());
        out.vertices.push_back(pos);
        edge_cache.emplace(key, idx);
        local_index[gi] = idx;
      } else {
        Vector3d pos = eval_pn(cp, uvw[0], uvw[1], uvw[2]);
        int idx = static_cast<int>(out.vertices.size());
        out.vertices.push_back(pos);
        local_index[gi] = idx;
      }
    }

    for (const auto& f : grid.faces) {
      out.indices.push_back({local_index[f[0]], local_index[f[1]], local_index[f[2]]});
    }
  }

  return out;
}

PyObject *python_organic(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  DECLARE_INSTANCE();
  auto node = std::make_shared<OrganicNode>(instance);
  char *kwlist[] = {"pts", "r", nullptr};
  double d = NAN;
  PyObject *points = nullptr;
  Vector3d point;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Od", kwlist, &points, &d)) {
    PyErr_SetString(PyExc_TypeError, "Error during parsing oragnic(ptlist,d)");
    return NULL;
  }
  if (points != NULL && PyList_Check(points)) {
    if (PyList_Size(points) == 0) {
      PyErr_SetString(PyExc_TypeError, "There must at least be one point in the polyhedron");
      return NULL;
    }
    for (Py_ssize_t i = 0; i < PyList_Size(points); i++) {
      PyObject *element = PyList_GetItem(points, i);
      if (PyList_Check(element) && PyList_Size(element) == 3) {
        point[0] = PyFloat_AsDouble(PyList_GetItem(element, 0));
        point[1] = PyFloat_AsDouble(PyList_GetItem(element, 1));
        point[2] = PyFloat_AsDouble(PyList_GetItem(element, 2));
        node->points.push_back(point);
      } else {
        PyErr_SetString(PyExc_TypeError, "Coordinate must exactly contain 3 numbers");
        return NULL;
      }
    }
  } else {
    PyErr_SetString(PyExc_TypeError, "Polyhedron Points must be a list of coordinates");
    return NULL;
  }
  node->d = d;
  return PyOpenSCADObjectFromNode(&PyOpenSCADType, node);
}

// -------------------------------------------------------------------------
// Optionaler Selbsttest: mit  g++ -DORGANIC_MESH_TEST_MAIN -I/usr/include/eigen3
//                             organic_mesh.cpp -o test && ./test
// -------------------------------------------------------------------------
#ifdef ORGANIC_MESH_TEST_MAIN
#include <cstdio>
#include <random>

int main()
{
  std::mt19937 rng(0);
  std::uniform_real_distribution<double> ang_dist(0.0, 2 * M_PI);
  std::uniform_real_distribution<double> rad_dist(0.0, 10.0);

  std::vector<Vector3d> pts_in;
  for (int i = 0; i < 40; ++i) {
    double ang = ang_dist(rng);
    double rad = rad_dist(rng);
    double x = rad * std::cos(ang);
    double y = rad * std::sin(ang);
    double z = 2.5 * std::exp(-(rad * rad) / 40.0);
    pts_in.emplace_back(x, y, z);
  }

  PolySet ps = organic_resample(pts_in, /*max_mesh_size=*/1.0);

  std::printf("Eingang: %zu Punkte\n", pts_in.size());
  std::printf("Ausgang: %zu Vertices, %zu Dreiecke\n", ps.vertices.size(), ps.triangles.size());

  // Grobe Sanity-Checks
  double max_edge = 0.0;
  for (const auto& t : ps.triangles) {
    max_edge = std::max(max_edge, (ps.vertices[t[0]] - ps.vertices[t[1]]).norm());
    max_edge = std::max(max_edge, (ps.vertices[t[1]] - ps.vertices[t[2]]).norm());
    max_edge = std::max(max_edge, (ps.vertices[t[2]] - ps.vertices[t[0]]).norm());
  }
  std::printf("Groesste Kantenlaenge im Ergebnis: %.4f (Ziel <= 1.0)\n", max_edge);
  return 0;
}
#endif
