#pragma once

#include <Python.h>
#include <src/core/matrix.h>
#include <src/core/primitives.h>
#include <vector>

// Python to C++ conversions
int python_tomatrix(PyObject *pyt, Matrix4d& mat);
int python_tovector(PyObject *pyt, Vector3d& vec);
int python_numberval(PyObject *number, double *result, int *flags = nullptr, int flagor = 0);
int python_vectorval(PyObject *vec, int minarg, int maxarg, double *x, double *y, double *z,
                     double *w = NULL, int *flags = nullptr);

// C++ to Python conversions
PyObject *python_frommatrix(const Matrix4d& mat);
PyObject *python_fromvector(const Vector3d vec);
PyObject *python_frompointsflexible(const std::vector<Vector3d>& points);
PyObject *python_frompaths(const std::vector<std::vector<size_t>>& paths);
PyObject *python_fromfaces(const std::vector<IndexedFace>& faces);
PyObject *python_fromopenscad(const Value& val);

// OpenSCAD Value conversions
Value python_convertresult(PyObject *arg, int& error);
