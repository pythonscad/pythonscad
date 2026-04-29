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

#include <math.h>
#include "linalg.h"
#include "GeometryUtils.h"
#include <Python.h>
#include "pyfunctions.h"
#include "python/pyconversion.h"

PyObject *python_math_sub1(PyObject *self, PyObject *args, PyObject *kwargs, int mode)
{
  char *kwlist[] = {"value", NULL};
  double arg;
  double result = 0;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "d", kwlist, &arg)) {
    PyErr_SetString(PyExc_TypeError, "Error during parsing math function");
    return NULL;
  }
  switch (mode) {
  case 0: result = sin(arg * M_PI / 180.0); break;
  case 1: result = cos(arg * M_PI / 180.0); break;
  case 2: result = tan(arg * M_PI / 180.0); break;
  case 3: result = asin(arg) * 180.0 / M_PI; break;
  case 4: result = acos(arg) * 180.0 / M_PI; break;
  case 5: result = atan(arg) * 180.0 / M_PI; break;
  }
  return PyFloat_FromDouble(result);
}

PyObject *python_math_sub2(PyObject *self, PyObject *args, PyObject *kwargs, int mode)
{
  int dragflags = 0;
  char *kwlist[] = {"vec1", "vec2", NULL};
  PyObject *obj1 = nullptr;
  PyObject *obj2 = nullptr;
  Vector3d vec31(0, 0, 0);
  Vector3d vec32(0, 0, 0);
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO", kwlist, &obj1, &obj2)) {
    PyErr_SetString(PyExc_TypeError, "Error during parsing norm(vec3)");
    return NULL;
  }
  python_vectorval(obj1, 1, 3, &(vec31[0]), &(vec31[1]), &(vec31[2]), nullptr, &dragflags);
  python_vectorval(obj2, 1, 3, &(vec32[0]), &(vec32[1]), &(vec32[2]), nullptr, &dragflags);

  switch (mode) {
  case 0: return PyFloat_FromDouble(vec31.dot(vec32)); break;
  case 1: return python_fromvector(vec31.cross(vec32)); break;
  }
  Py_RETURN_NONE;
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
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &obj)) {
    PyErr_SetString(PyExc_TypeError, "Error during parsing norm(vec3)");
    return NULL;
  }
  python_vectorval(obj, 1, 3, &(vec3[0]), &(vec3[1]), &(vec3[2]), nullptr, &dragflags);

  result = sqrt(vec3[0] * vec3[0] + vec3[1] * vec3[1] + vec3[2] * vec3[2]);
  return PyFloat_FromDouble(result);
}
