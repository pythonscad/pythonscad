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

// organic_mesh.cpp
// -----------------------------------------------------------------------
// Eigenstaendige C++-Funktionen (keine CPython-Bindung noetig):
//
//   PolySet organic_resample(const std::vector<Vector3d>& points,
//                             double max_mesh_size);
//
// Immer fuer GESCHLOSSENE Koerper gedacht (kein Hoehenfeld-Modus).
// Auch konkave Formen werden unterstuetzt, da die Oberflaeche ueber eine
// echte 3D-Delaunay-Tetraedrisierung + Alpha-Shape-Filterung rekonstruiert
// wird (nicht nur ueber die konvexe Huelle).
//
// Ablauf:
//   1) delaunay3d_bowyer_watson()   - 3D-Delaunay-Tetraedrisierung
//                                     (Bowyer-Watson, inkrementell)
//   2) alpha_shape_boundary()       - Alpha-Shape-Filterung: nur "kompakte"
//                                     Tetraeder (kleiner Umkugelradius)
//                                     behalten, davon die Aussenflaechen
//                                     extrahieren. Dadurch koennen auch
//                                     konkave Einbuchtungen abgebildet
//                                     werden (die reine konvexe Huelle
//                                     wuerde sie "wegbuegeln").
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
//                                     Kantenpunkte UND Original-Ecken werden
//                                     dedupliziert, das Ergebnis ist
//                                     wasserdicht.
//
// Abhaengigkeit: nur Eigen (Vector3d), keine weiteren Libraries.
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
#include <Eigen/Dense>
#include <vector>
#include <array>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cassert>
#include <limits>

using Vector3d = Eigen::Vector3d;

// =========================================================================
// 1) 3D-Delaunay-Tetraedrisierung (Bowyer-Watson, inkrementell)
// =========================================================================
namespace detail {

struct Tet4 {
  int a, b, c, d;
};

inline double signed_volume6(const Vector3d& a, const Vector3d& b, const Vector3d& c, const Vector3d& d)
{
  return (b - a).dot((c - a).cross(d - a));
}

inline bool circumsphere(const Vector3d& a, const Vector3d& b, const Vector3d& c, const Vector3d& d,
                         Vector3d& center, double& radius)
{
  Eigen::Matrix3d M;
  Eigen::Vector3d rhs;
  M.row(0) = (b - a).transpose();
  M.row(1) = (c - a).transpose();
  M.row(2) = (d - a).transpose();
  rhs(0) = 0.5 * (b.squaredNorm() - a.squaredNorm());
  rhs(1) = 0.5 * (c.squaredNorm() - a.squaredNorm());
  rhs(2) = 0.5 * (d.squaredNorm() - a.squaredNorm());

  double det = M.determinant();
  if (std::fabs(det) < 1e-14) return false;

  center = M.fullPivLu().solve(rhs);
  radius = (center - a).norm();
  return true;
}

struct FaceKey {
  int a, b, c;
  bool operator==(const FaceKey& o) const { return a == o.a && b == o.b && c == o.c; }
};
struct FaceKeyHash {
  size_t operator()(const FaceKey& f) const
  {
    return (static_cast<size_t>(f.a) << 40) ^ (static_cast<size_t>(f.b) << 20) ^
           static_cast<size_t>(f.c);
  }
};
inline FaceKey make_face_key(int a, int b, int c)
{
  int arr[3] = {a, b, c};
  std::sort(arr, arr + 3);
  return {arr[0], arr[1], arr[2]};
}

inline std::vector<Tet4> delaunay3d_bowyer_watson(const std::vector<Vector3d>& points)
{
  const int n = static_cast<int>(points.size());
  std::vector<Vector3d> pts = points;

  // Mini-Jitter GEGEN NUMERISCHE ENTARTUNG: liegen viele Punkte (fast)
  // exakt auf einer gemeinsamen Kugel (z.B. Kugel-Sampling), haben viele
  // verschiedene Tetraeder exakt denselben Umkugelradius - die
  // Bowyer-Watson-Kombinatorik wird dann mehrdeutig und kann nicht-
  // mannigfaltige Facetten erzeugen. Ein winziges, deterministisches
  // Jitter (nur fuer die Topologie-Berechnung hier, NICHT fuer die
  // spaeter zurueckgegebene Geometrie) bricht solche Ties robust auf.
  {
    Vector3d lo0(1e300, 1e300, 1e300), hi0(-1e300, -1e300, -1e300);
    for (const auto& p : points) {
      lo0 = lo0.cwiseMin(p);
      hi0 = hi0.cwiseMax(p);
    }
    double diag = (hi0 - lo0).norm();
    if (diag < 1e-12) diag = 1.0;
    const double jitter_scale = diag * 1e-7;
    for (int i = 0; i < n; ++i) {
      uint32_t h = static_cast<uint32_t>(i * 2654435761u);
      auto rnd = [&](int k) {
        h ^= h << 13;
        h ^= h >> 17;
        h ^= h << 5;
        h += k * 2246822519u;
        return (static_cast<double>(h) / 4294967295.0) * 2.0 - 1.0;
      };
      pts[i] += Vector3d(rnd(1), rnd(2), rnd(3)) * jitter_scale;
    }
  }

  Vector3d lo(1e300, 1e300, 1e300), hi(-1e300, -1e300, -1e300);
  for (const auto& p : points) {
    lo = lo.cwiseMin(p);
    hi = hi.cwiseMax(p);
  }
  Vector3d center = (lo + hi) * 0.5;
  double extent = (hi - lo).norm() * 10.0 + 1.0;

  const int s0 = n, s1 = n + 1, s2 = n + 2, s3 = n + 3;
  pts.push_back(center + Vector3d(-extent, -extent, -extent));
  pts.push_back(center + Vector3d(extent, -extent, -extent));
  pts.push_back(center + Vector3d(0, extent, -extent));
  pts.push_back(center + Vector3d(0, 0, extent));

  std::vector<Tet4> tets;
  tets.push_back({s0, s1, s2, s3});
  if (signed_volume6(pts[s0], pts[s1], pts[s2], pts[s3]) < 0) std::swap(tets[0].a, tets[0].b);

  for (int pi = 0; pi < n; ++pi) {
    const Vector3d& p = pts[pi];

    std::vector<int> bad;
    bad.reserve(16);
    for (int t = 0; t < static_cast<int>(tets.size()); ++t) {
      const Tet4& tt = tets[t];
      Vector3d c;
      double r;
      if (!circumsphere(pts[tt.a], pts[tt.b], pts[tt.c], pts[tt.d], c, r)) continue;
      if ((p - c).norm() < r - 1e-9) bad.push_back(t);
    }
    if (bad.empty()) continue;

    std::unordered_map<FaceKey, int, FaceKeyHash> face_count;
    std::unordered_map<FaceKey, std::array<int, 3>, FaceKeyHash> face_orient;
    auto add_face = [&](int a, int b, int c) {
      FaceKey k = make_face_key(a, b, c);
      face_count[k]++;
      face_orient[k] = {a, b, c};
    };
    for (int t : bad) {
      const Tet4& tt = tets[t];
      add_face(tt.b, tt.c, tt.d);
      add_face(tt.a, tt.c, tt.d);
      add_face(tt.a, tt.b, tt.d);
      add_face(tt.a, tt.b, tt.c);
    }

    std::sort(bad.rbegin(), bad.rend());
    for (int t : bad) tets.erase(tets.begin() + t);

    for (auto& kv : face_count) {
      if (kv.second != 1) continue;
      auto& f = face_orient[kv.first];
      int a = f[0], b = f[1], c = f[2];
      if (signed_volume6(pts[a], pts[b], pts[c], pts[pi]) < 0) std::swap(a, b);
      tets.push_back({a, b, c, pi});
    }
  }

  std::vector<Tet4> result;
  result.reserve(tets.size());
  for (const auto& t : tets) {
    if (t.a >= n || t.b >= n || t.c >= n || t.d >= n) continue;
    result.push_back(t);
  }
  return result;
}

}  // namespace detail

// =========================================================================
// 2) Alpha-Shape: aus der Delaunay-Tetraedrisierung nur "kompakte"
//    Tetraeder behalten und deren Aussenflaechen extrahieren. Das erlaubt
//    (im Gegensatz zur reinen konvexen Huelle) auch konkave Formen.
// =========================================================================
namespace detail {

struct Face3 {
  int a, b, c;
};

inline std::vector<Face3> alpha_shape_boundary(const std::vector<Vector3d>& points,
                                               const std::vector<Tet4>& tets, double alpha)
{
  std::vector<char> keep(tets.size(), 0);
  for (size_t i = 0; i < tets.size(); ++i) {
    const Tet4& t = tets[i];
    Vector3d c;
    double r;
    if (!circumsphere(points[t.a], points[t.b], points[t.c], points[t.d], c, r)) continue;
    if (r <= alpha) keep[i] = 1;
  }

  std::unordered_map<FaceKey, std::vector<std::pair<int, int>>, FaceKeyHash> face_owners;
  for (size_t i = 0; i < tets.size(); ++i) {
    if (!keep[i]) continue;
    const Tet4& t = tets[i];
    int verts[4] = {t.a, t.b, t.c, t.d};
    for (int excl = 0; excl < 4; ++excl) {
      int f[3];
      int k = 0;
      for (int v = 0; v < 4; ++v)
        if (v != excl) f[k++] = verts[v];
      FaceKey key = make_face_key(f[0], f[1], f[2]);
      face_owners[key].push_back({static_cast<int>(i), verts[excl]});
    }
  }

  std::vector<Face3> boundary;
  boundary.reserve(face_owners.size() / 2);
  for (auto& kv : face_owners) {
    if (kv.second.size() != 1) continue;
    int tet_idx = kv.second[0].first;
    int excluded = kv.second[0].second;
    const Tet4& t = tets[tet_idx];
    int verts[4] = {t.a, t.b, t.c, t.d};
    int f[3];
    int k = 0;
    for (int v = 0; v < 4; ++v)
      if (verts[v] != excluded) f[k++] = verts[v];

    Vector3d nrm = (points[f[1]] - points[f[0]]).cross(points[f[2]] - points[f[0]]);
    Vector3d toExcluded = points[excluded] - points[f[0]];
    if (nrm.dot(toExcluded) > 0) std::swap(f[1], f[2]);

    boundary.push_back({f[0], f[1], f[2]});
  }
  return boundary;
}

// Automatische Alpha-Wahl: die reine Abschaetzung ueber den naechsten-
// Nachbar-Abstand versagt bei Punktwolken, die (wie eine Kugel) auf einer
// gemeinsamen Oberflaeche liegen - dort ist der Umkugelradius jedes
// Delaunay-Tetraeders durch die Gesamtkruemmung bestimmt, nicht durch den
// lokalen Punktabstand, und kann um Groessenordnungen groesser sein.
//
// Stattdessen wird das KLEINSTE Alpha gesucht, bei dem noch jeder
// Eingabepunkt auf der resultierenden Oberflaeche erscheint UND das
// Ergebnis mannigfaltig ist (jede Kante gehoert zu genau 2 Dreiecken).
// Reine Punkt-Abdeckung allein reicht nicht: bei duenn abgetasteten,
// fast linienfoermigen "Spitzen" (z.B. nur 3 Punkte nahe einer Geraden)
// kann das kleinste abdeckende Alpha eine beidseitige, quasi nulldicke
// "Lasche" erzeugen (zwei fast deckungsgleiche, entgegengesetzt
// orientierte Dreiecke). Der naechstgroessere kritische Alpha-Wert
// loest das in der Praxis fast immer auf, da genau dort die Lasche durch
// ein einzelnes, "dickeres" Tetraeder ersetzt wird.
//
// (Die kritischen Alpha-Werte sind exakt die Umkugelradien selbst,
// dazwischen aendert sich die Topologie nicht - daher reicht es, diese
// der Groesse nach zu durchsuchen.)
inline bool is_manifold(const std::vector<Face3>& faces)
{
  std::unordered_map<int64_t, int> edge_count;
  auto key = [](int a, int b) {
    int lo = std::min(a, b), hi = std::max(a, b);
    return (static_cast<int64_t>(lo) << 32) | static_cast<int64_t>(hi);
  };
  for (const auto& f : faces) {
    edge_count[key(f.a, f.b)]++;
    edge_count[key(f.b, f.c)]++;
    edge_count[key(f.c, f.a)]++;
  }
  for (auto& kv : edge_count)
    if (kv.second != 2) return false;
  return true;
}

inline double estimate_alpha(const std::vector<Vector3d>& points, const std::vector<Tet4>& tets)
{
  std::vector<double> radii;
  radii.reserve(tets.size());
  for (const auto& t : tets) {
    Vector3d c;
    double r;
    if (circumsphere(points[t.a], points[t.b], points[t.c], points[t.d], c, r)) radii.push_back(r);
  }
  if (radii.empty()) return 1.0;
  std::sort(radii.begin(), radii.end());

  // Erste Runde: kleinstes Alpha mit voller Punkt-Abdeckung UND
  // Mannigfaltigkeit.
  for (double r : radii) {
    auto faces = alpha_shape_boundary(points, tets, r);
    std::vector<char> covered(points.size(), 0);
    for (const auto& f : faces) {
      covered[f.a] = 1;
      covered[f.b] = 1;
      covered[f.c] = 1;
    }
    bool all_covered = true;
    for (char c : covered)
      if (!c) {
        all_covered = false;
        break;
      }
    if (all_covered && is_manifold(faces)) return r;
  }
  // Fallback: kleinstes Alpha mit voller Abdeckung, auch wenn nicht
  // perfekt mannigfaltig (seltener Grenzfall bei extrem entarteten
  // Punktwolken) - besser als komplett leer zurueckzugeben.
  for (double r : radii) {
    auto faces = alpha_shape_boundary(points, tets, r);
    std::vector<char> covered(points.size(), 0);
    for (const auto& f : faces) {
      covered[f.a] = 1;
      covered[f.b] = 1;
      covered[f.c] = 1;
    }
    bool all_covered = true;
    for (char c : covered)
      if (!c) {
        all_covered = false;
        break;
      }
    if (all_covered) return r;
  }
  return radii.back();  // Fallback: entspricht der reinen konvexen Huelle
}

}  // namespace detail

// =========================================================================
// 3) Loop-Subdivision: echte G1-glatte Oberflaeche ohne sichtbare Kanten
// =========================================================================
//
// PN-Triangles (frueherer Ansatz) garantieren nur, dass die Randkurve
// zwischen zwei Nachbardreiecken IDENTISCH ist (keine Risse) - die
// Tangentialebene QUER zur Kante wird dabei aber nicht angeglichen. Das
// ergibt einen festen "Knick" entlang der Original-Kanten, der auch bei
// feinerer Unterteilung nicht verschwindet.
//
// Loop-Subdivision (Charles Loop, 1987) loest das strukturell: jede Kante
// bekommt einen neuen Punkt als gewichtetes Mittel aus den beiden
// Kanten-Endpunkten UND den beiden gegenueberliegenden Dreieckspitzen;
// bestehende Punkte werden leicht in Richtung ihrer Nachbarn verschoben.
// Das Ergebnis ist echte G1-Glattheit (tatsaechlich C2 fast ueberall) -
// keine Knicke mehr, unabhaengig vom Winkel der Ursprungskanten.
//
// Trade-off: Loop-Subdivision ist ein APPROXIMIERENDES Verfahren - die
// Original-Punkte werden dabei leicht verschoben (nicht mehr exakt
// interpoliert wie beim vorherigen PN-Triangle-Ansatz).
// =========================================================================
namespace detail {

struct EdgeKey2 {
  int lo, hi;
  bool operator==(const EdgeKey2& o) const { return lo == o.lo && hi == o.hi; }
};
struct EdgeKey2Hash {
  size_t operator()(const EdgeKey2& e) const
  {
    return (static_cast<size_t>(e.lo) << 32) ^ static_cast<size_t>(e.hi);
  }
};
inline EdgeKey2 make_edge_key2(int a, int b)
{
  return {std::min(a, b), std::max(a, b)};
}

struct LoopMesh {
  std::vector<Vector3d> verts;
  PolygonIndices tris;
};

inline double max_edge_length(const LoopMesh& m)
{
  double e = 0.0;
  for (const auto& t : m.tris) {
    e = std::max(e, (m.verts[t[0]] - m.verts[t[1]]).norm());
    e = std::max(e, (m.verts[t[1]] - m.verts[t[2]]).norm());
    e = std::max(e, (m.verts[t[2]] - m.verts[t[0]]).norm());
  }
  return e;
}

inline LoopMesh loop_subdivide_once(const LoopMesh& in)
{
  const int n = static_cast<int>(in.verts.size());

  // Kante -> Liste von (gegenueberliegender Eckpunkt). Bei einer
  // geschlossenen, mannigfaltigen Flaeche hat jede Kante genau 2
  // Eintraege; bei offenen/nicht ganz sauberen Raendern (Sicherheitsnetz
  // fuer seltene Alpha-Shape-Randfaelle) koennen es auch 1 oder >2 sein.
  std::unordered_map<EdgeKey2, std::vector<int>, EdgeKey2Hash> edge_opposite;
  std::unordered_map<int, std::vector<int>> neighbors;  // alle Nachbarn (fuer Valenz-Regel)
  std::unordered_map<int, std::vector<int>>
    boundary_neighbors;  // nur ueber Rand-Kanten erreichte Nachbarn

  auto add_edge = [&](int a, int b, int opp) { edge_opposite[make_edge_key2(a, b)].push_back(opp); };
  for (const auto& t : in.tris) {
    add_edge(t[0], t[1], t[2]);
    add_edge(t[1], t[2], t[0]);
    add_edge(t[2], t[0], t[1]);
  }
  for (auto& kv : edge_opposite) {
    int a = kv.first.lo, b = kv.first.hi;
    neighbors[a].push_back(b);
    neighbors[b].push_back(a);
    if (kv.second.size() == 1) {
      boundary_neighbors[a].push_back(b);
      boundary_neighbors[b].push_back(a);
    }
  }

  LoopMesh out;
  out.verts.resize(n);  // Platz fuer die repositionierten Original-Punkte

  // --- Original-Punkte reposition (Loop-Glaettungsregel) ---
  for (int v = 0; v < n; ++v) {
    auto bit = boundary_neighbors.find(v);
    if (bit != boundary_neighbors.end() && bit->second.size() == 2) {
      // Randpunkt: 3/4*alt + 1/8*(beide Randnachbarn)
      out.verts[v] = 0.75 * in.verts[v] + 0.125 * (in.verts[bit->second[0]] + in.verts[bit->second[1]]);
    } else {
      auto nit = neighbors.find(v);
      if (nit == neighbors.end() || nit->second.empty()) {
        out.verts[v] = in.verts[v];  // isolierter Punkt (sollte nicht vorkommen)
        continue;
      }
      const auto& nb = nit->second;
      int valence = static_cast<int>(nb.size());
      double beta = (valence == 3) ? (3.0 / 16.0) : (3.0 / (8.0 * valence));
      Vector3d sum = Vector3d::Zero();
      for (int u : nb) sum += in.verts[u];
      out.verts[v] = (1.0 - valence * beta) * in.verts[v] + beta * sum;
    }
  }

  // --- Neue Kantenpunkte ---
  std::unordered_map<EdgeKey2, int, EdgeKey2Hash> edge_point_index;
  for (auto& kv : edge_opposite) {
    int a = kv.first.lo, b = kv.first.hi;
    Vector3d pos;
    if (kv.second.size() == 2) {
      int o1 = kv.second[0], o2 = kv.second[1];
      pos = 0.375 * (in.verts[a] + in.verts[b]) + 0.125 * (in.verts[o1] + in.verts[o2]);
    } else {
      // Randkante oder nicht-mannigfaltiger Sonderfall: einfacher Mittelpunkt
      pos = 0.5 * (in.verts[a] + in.verts[b]);
    }
    int idx = static_cast<int>(out.verts.size());
    out.verts.push_back(pos);
    edge_point_index[kv.first] = idx;
  }

  // --- Neue Topologie: jedes alte Dreieck -> 4 neue ---
  out.tris.reserve(in.tris.size() * 4);
  for (const auto& t : in.tris) {
    int a = t[0], b = t[1], c = t[2];
    int mab = edge_point_index[make_edge_key2(a, b)];
    int mbc = edge_point_index[make_edge_key2(b, c)];
    int mca = edge_point_index[make_edge_key2(c, a)];
    out.tris.push_back({a, mab, mca});
    out.tris.push_back({b, mbc, mab});
    out.tris.push_back({c, mca, mbc});
    out.tris.push_back({mab, mbc, mca});
  }

  return out;
}

inline LoopMesh loop_subdivide_to_target(const std::vector<Vector3d>& points, const PolygonIndices& tris,
                                         double max_mesh_size)
{
  LoopMesh mesh{points, tris};
  if (max_mesh_size <= 0.0) return mesh;

  double cur_max = max_edge_length(mesh);
  int levels = 0;
  double predicted = cur_max;
  // Jede Stufe halbiert die Kantenlaenge ungefaehr (typisches Verhalten
  // von Loop-Subdivision auf einem einigermassen regulaeren Netz).
  while (predicted > max_mesh_size && levels < 8) {
    predicted *= 0.5;
    levels++;
  }
  for (int i = 0; i < levels; ++i) mesh = loop_subdivide_once(mesh);
  return mesh;
}

}  // namespace detail

// =========================================================================
// Hauptfunktion
// =========================================================================
PolySet organic_resample(const std::vector<Vector3d>& points, double max_mesh_size, double alpha = -1.0)
{
  PolySet out(3);
  if (points.size() < 4) return out;

  // --- 1) 3D-Delaunay-Tetraedrisierung ---
  std::vector<detail::Tet4> tets = detail::delaunay3d_bowyer_watson(points);

  // --- 2) Alpha-Shape: Aussenflaeche extrahieren (erlaubt Konkavitaeten) ---
  double eff_alpha = (alpha > 0.0) ? alpha : detail::estimate_alpha(points, tets);
  std::vector<detail::Face3> faces = detail::alpha_shape_boundary(points, tets, eff_alpha);

  PolygonIndices tris;
  tris.reserve(faces.size());
  for (const auto& f : faces) tris.push_back({f.a, f.b, f.c});
  if (tris.empty()) return out;

  // --- 3) Loop-Subdivision: glatt, ohne sichtbare Kanten ---
  detail::LoopMesh smooth = detail::loop_subdivide_to_target(points, tris, max_mesh_size);

  out.vertices = std::move(smooth.verts);
  out.indices = std::move(smooth.tris);
  return out;
}

// -------------------------------------------------------------------------
// Optionaler Selbsttest: mit  g++ -DORGANIC_MESH_TEST_MAIN -I/usr/include/eigen3
//                             organic_mesh.cpp -o test && ./test
// -------------------------------------------------------------------------
#ifdef ORGANIC_MESH_TEST_MAIN
#include <cstdio>

int main()
{
  std::vector<Vector3d> octa = {Vector3d(1, 0, 0),  Vector3d(-1, 0, 0), Vector3d(0, 1, 0),
                                Vector3d(0, -1, 0), Vector3d(0, 0, 1),  Vector3d(0, 0, -1)};
  PolySet ps = organic_resample(octa, 0.3);
  std::printf("Oktaeder: %zu Vertices, %zu Dreiecke\n", ps.vertices.size(), ps.triangles.size());
  return 0;
}
#endif
