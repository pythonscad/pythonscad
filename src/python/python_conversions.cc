#include "python_conversions.h"
#include <src/core/primitives.h>
#include <src/core/matrix.h>
#include <cstring>

// Python to C++ Conversions
int python_tomatrix(PyObject *pyt, Matrix4d& mat)
{
  if (!PyList_Check(pyt)) return -1;
  if (PyList_Size(pyt) != 16) return -2;

  for (int i = 0; i < 16; i++) {
    PyObject *obj = PyList_GetItem(pyt, i);
    if (!PyNumber_Check(obj)) return -3;
    mat.data[i] = PyFloat_AsDouble(obj);
    if (PyErr_Occurred()) return -4;
  }
  return 0;
}

int python_tovector(PyObject *pyt, Vector3d& vec)
{
  if (!PyList_Check(pyt)) return -1;
  int size = PyList_Size(pyt);
  if (size < 2 || size > 3) return -2;

  for (int i = 0; i < size; i++) {
    PyObject *obj = PyList_GetItem(pyt, i);
    if (!PyNumber_Check(obj)) return -3;
    vec[i] = PyFloat_AsDouble(obj);
    if (PyErr_Occurred()) return -4;
  }
  if (size == 2) vec[2] = 0.0;
  return 0;
}

int python_numberval(PyObject *number, double *result, int *flags, int flagor)
{
  if (!PyNumber_Check(number)) return -1;

  *result = PyFloat_AsDouble(number);
  if (PyErr_Occurred()) {
    PyErr_Clear();
    if (PyLong_Check(number)) {
      *result = (double)PyLong_AsLong(number);
      if (PyErr_Occurred()) return -2;
    } else {
      return -3;
    }
  }

  if (flags) *flags |= flagor;
  return 0;
}

int python_vectorval(PyObject *vec, int minarg, int maxarg, double *x, double *y, double *z, double *w,
                     int *flags)
{
  if (!PyList_Check(vec) && !PyTuple_Check(vec)) return -1;

  int size = PyList_Check(vec) ? PyList_Size(vec) : PyTuple_Size(vec);
  if (size < minarg || size > maxarg) return -2;

  for (int i = 0; i < size; i++) {
    PyObject *obj = PyList_Check(vec) ? PyList_GetItem(vec, i) : PyTuple_GetItem(vec, i);
    if (!PyNumber_Check(obj)) return -3;

    double val = PyFloat_AsDouble(obj);
    if (PyErr_Occurred()) {
      PyErr_Clear();
      if (PyLong_Check(obj)) {
        val = (double)PyLong_AsLong(obj);
        if (PyErr_Occurred()) return -4;
      } else {
        return -5;
      }
    }

    if (i == 0) *x = val;
    else if (i == 1) *y = val;
    else if (i == 2) *z = val;
    else if (i == 3 && w) *w = val;
  }

  return 0;
}

// C++ to Python Conversions
PyObject *python_frommatrix(const Matrix4d& mat)
{
  PyObject *list = PyList_New(16);
  if (!list) return NULL;

  for (int i = 0; i < 16; i++) {
    PyObject *val = PyFloat_FromDouble(mat.data[i]);
    if (!val) {
      Py_DECREF(list);
      return NULL;
    }
    PyList_SetItem(list, i, val);
  }
  return list;
}

PyObject *python_fromvector(const Vector3d vec)
{
  PyObject *list = PyList_New(3);
  if (!list) return NULL;

  for (int i = 0; i < 3; i++) {
    PyObject *val = PyFloat_FromDouble(vec[i]);
    if (!val) {
      Py_DECREF(list);
      return NULL;
    }
    PyList_SetItem(list, i, val);
  }
  return list;
}

PyObject *python_frompointsflexible(const std::vector<Vector3d>& points)
{
  PyObject *list = PyList_New(points.size());
  if (!list) return NULL;

  for (size_t i = 0; i < points.size(); i++) {
    PyObject *point;
    if (points[i][2] == 0.0) {
      point = PyList_New(2);
      if (!point) {
        Py_DECREF(list);
        return NULL;
      }
      PyList_SetItem(point, 0, PyFloat_FromDouble(points[i][0]));
      PyList_SetItem(point, 1, PyFloat_FromDouble(points[i][1]));
    } else {
      point = PyList_New(3);
      if (!point) {
        Py_DECREF(list);
        return NULL;
      }
      PyList_SetItem(point, 0, PyFloat_FromDouble(points[i][0]));
      PyList_SetItem(point, 1, PyFloat_FromDouble(points[i][1]));
      PyList_SetItem(point, 2, PyFloat_FromDouble(points[i][2]));
    }
    PyList_SetItem(list, i, point);
  }
  return list;
}

PyObject *python_frompaths(const std::vector<std::vector<size_t>>& paths)
{
  PyObject *list = PyList_New(paths.size());
  if (!list) return NULL;

  for (size_t i = 0; i < paths.size(); i++) {
    PyObject *path = PyList_New(paths[i].size());
    if (!path) {
      Py_DECREF(list);
      return NULL;
    }
    for (size_t j = 0; j < paths[i].size(); j++) {
      PyList_SetItem(path, j, PyLong_FromLong(paths[i][j]));
    }
    PyList_SetItem(list, i, path);
  }
  return list;
}

PyObject *python_fromfaces(const std::vector<IndexedFace>& faces)
{
  PyObject *list = PyList_New(faces.size());
  if (!list) return NULL;

  for (size_t i = 0; i < faces.size(); i++) {
    PyObject *face = PyList_New(faces[i].indices.size());
    if (!face) {
      Py_DECREF(list);
      return NULL;
    }
    for (size_t j = 0; j < faces[i].indices.size(); j++) {
      PyList_SetItem(face, j, PyLong_FromLong(faces[i].indices[j]));
    }
    PyList_SetItem(list, i, face);
  }
  return list;
}

PyObject *python_fromopenscad(const Value& val)
{
  if (val.type() == Value::NUMBER) {
    return PyFloat_FromDouble(val.toDouble());
  } else if (val.type() == Value::STRING) {
    return PyUnicode_FromString(val.toString().c_str());
  } else if (val.type() == Value::VECTOR) {
    PyObject *list = PyList_New(val.vector_size());
    if (!list) return NULL;
    for (size_t i = 0; i < val.vector_size(); i++) {
      PyList_SetItem(list, i, python_fromopenscad(val.vector_get(i)));
    }
    return list;
  }
  Py_RETURN_NONE;
}

Value python_convertresult(PyObject *arg, int& error)
{
  error = 0;

  if (PyFloat_Check(arg)) {
    return Value(PyFloat_AsDouble(arg));
  } else if (PyLong_Check(arg)) {
    return Value((double)PyLong_AsLong(arg));
  } else if (PyUnicode_Check(arg)) {
    return Value(std::string(PyUnicode_AsUTF8(arg)));
  } else if (PyList_Check(arg)) {
    std::vector<Value> vec;
    for (Py_ssize_t i = 0; i < PyList_Size(arg); i++) {
      int err = 0;
      vec.push_back(python_convertresult(PyList_GetItem(arg, i), err));
      if (err) {
        error = err;
        return Value();
      }
    }
    return Value(vec);
  }

  error = -1;
  return Value();
}
