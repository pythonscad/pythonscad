/*
 *  PythonSCAD - internal() / external(): width and spacing checks between
 *  facing surfaces, the 3D counterpart of the DRC rules INTERNAL / EXTERNAL.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */

#include "genlang/genlang.h"
#include <Python.h>
#include <memory>
#include "pyopenscad.h"
#include <Tree.h>
#include <GeometryEvaluator.h>
#include <PolySetUtils.h>
#include <FacingCheckNode.h>
#include "geometry/facing_check.h"
#include "pyfunctions.h"

namespace {

std::shared_ptr<const PolySet> evaluateToPolySet(const std::shared_ptr<AbstractNode>& node)
{
  Tree tree(node, "");
  GeometryEvaluator geomevaluator(tree);
  std::shared_ptr<const Geometry> geom = geomevaluator.evaluateGeometry(*tree.root(), true);
  return geom ? PolySetUtils::getGeometryAsPolySet(geom) : nullptr;
}

// report=True: evaluate right away and return {"count": n, "min": d or None}
PyObject *facingReport(const std::shared_ptr<AbstractNode>& child,
                       const std::shared_ptr<AbstractNode>& other, FacingCheck::Mode mode, double d,
                       const FacingCheck::Options& opt_in)
{
  FacingCheck::Options opt = opt_in;
  opt.store_pairs = false;  // summary only, O(N) memory
  const auto ps = evaluateToPolySet(child);
  FacingCheck::Result res;
  if (ps && other) {
    const auto ps2 = evaluateToPolySet(other);
    if (ps2) res = FacingCheck::checkBetween(*ps, *ps2, d, opt);
  } else if (ps) {
    res = FacingCheck::check(*ps, mode, d, opt);
  }
  PyObject *dict = PyDict_New();
  auto dict_owner = py_owned(dict);
  PyObject *count = PyLong_FromSize_t(res.violation_count);
  PyDict_SetItemString(dict, "count", count);
  Py_DECREF(count);
  if (res.clean()) {
    PyDict_SetItemString(dict, "min", Py_None);
  } else {
    PyObject *mn = PyFloat_FromDouble(res.min_distance);
    PyDict_SetItemString(dict, "min", mn);
    Py_DECREF(mn);
  }
  return dict_owner.release();
}

PyObject *python_facing_core(PyObject *obj, PyObject *other, double d, double angle, double alpha,
                             int occlusion, PyObject *grow_obj, int report, FacingCheck::Mode mode)
{
  const char *fname = mode == FacingCheck::Mode::Internal ? "internal" : "external";
  if (!(d > 0)) {
    PyErr_Format(PyExc_ValueError, "%s(): d must be > 0", fname);
    return nullptr;
  }
  if (!(angle > 90 && angle <= 180)) {
    PyErr_Format(PyExc_ValueError, "%s(): angle must be in (90, 180]", fname);
    return nullptr;
  }
  if (!(alpha >= 0)) {
    PyErr_Format(PyExc_ValueError, "%s(): alpha must be >= 0", fname);
    return nullptr;
  }

  PyObject *dummydict = nullptr;
  PyTypeObject *type = PyOpenSCADObjectType(obj);
  std::shared_ptr<AbstractNode> child = PyOpenSCADObjectToNodeMulti(obj, &dummydict);
  auto dummydict_owner = py_owned(dummydict);
  if (child == nullptr) return propagate_or_typeerror("Invalid type for Object in internal/external\n");

  std::shared_ptr<AbstractNode> child2;
  if (other != nullptr && other != Py_None) {
    if (mode == FacingCheck::Mode::Internal) {
      PyErr_SetString(PyExc_TypeError, "internal(): other is only supported by external()");
      return nullptr;
    }
    PyObject *dummydict2 = nullptr;
    child2 = PyOpenSCADObjectToNodeMulti(other, &dummydict2);
    auto dummydict2_owner = py_owned(dummydict2);
    if (child2 == nullptr) return propagate_or_typeerror("Invalid type for other in external\n");
  }

  double grow = -1;  // automatic
  if (grow_obj != nullptr && grow_obj != Py_None) {
    grow = PyFloat_AsDouble(grow_obj);
    if (PyErr_Occurred()) return nullptr;
    if (grow < 0) {
      PyErr_Format(PyExc_ValueError, "%s(): grow must be >= 0", fname);
      return nullptr;
    }
  }

  FacingCheck::Options opt;
  opt.min_angle_deg = angle;
  opt.alpha_deg = alpha;
  opt.occlusion = occlusion != 0;
  if (report) return facingReport(child, child2, mode, d, opt);

  DECLARE_INSTANCE();
  auto node = std::make_shared<FacingCheckNode>(instance);
  node->mode = mode;
  node->distance = d;
  node->min_angle = angle;
  node->alpha = alpha;
  node->occlusion = occlusion != 0;
  node->grow = grow;
  node->children.push_back(child);
  if (child2) node->children.push_back(child2);
  return PyOpenSCADObjectFromNode(type, node);
}

PyObject *python_facing(PyObject *args, PyObject *kwargs, FacingCheck::Mode mode)
{
  const bool ext = mode == FacingCheck::Mode::External;
  char *kwlist_int[] = {"obj", "d", "angle", "alpha", "occlusion", "grow", "report", NULL};
  char *kwlist_ext[] = {"obj", "d", "other", "angle", "alpha", "occlusion", "grow", "report", NULL};
  PyObject *obj = nullptr, *other = nullptr;
  double d = 0, angle = 120, alpha = 90;
  int occlusion = 1, report = 0;
  PyObject *grow = nullptr;
  const bool ok = ext ? PyArg_ParseTupleAndKeywords(args, kwargs, "Od|OddpOp", kwlist_ext, &obj, &d,
                                                    &other, &angle, &alpha, &occlusion, &grow, &report)
                      : PyArg_ParseTupleAndKeywords(args, kwargs, "Od|ddpOp", kwlist_int, &obj, &d,
                                                    &angle, &alpha, &occlusion, &grow, &report);
  if (!ok) {
    PyErr_SetString(PyExc_TypeError, "error during parsing\n");
    return nullptr;
  }
  return python_facing_core(obj, other, d, angle, alpha, occlusion, grow, report, mode);
}

PyObject *python_oo_facing(PyObject *self, PyObject *args, PyObject *kwargs, FacingCheck::Mode mode)
{
  const bool ext = mode == FacingCheck::Mode::External;
  char *kwlist_int[] = {"d", "angle", "alpha", "occlusion", "grow", "report", NULL};
  char *kwlist_ext[] = {"d", "other", "angle", "alpha", "occlusion", "grow", "report", NULL};
  PyObject *other = nullptr;
  double d = 0, angle = 120, alpha = 90;
  int occlusion = 1, report = 0;
  PyObject *grow = nullptr;
  const bool ok = ext ? PyArg_ParseTupleAndKeywords(args, kwargs, "d|OddpOp", kwlist_ext, &d, &other,
                                                    &angle, &alpha, &occlusion, &grow, &report)
                      : PyArg_ParseTupleAndKeywords(args, kwargs, "d|ddpOp", kwlist_int, &d, &angle,
                                                    &alpha, &occlusion, &grow, &report);
  if (!ok) {
    PyErr_SetString(PyExc_TypeError, "error during parsing\n");
    return nullptr;
  }
  return python_facing_core(self, other, d, angle, alpha, occlusion, grow, report, mode);
}

}  // namespace

PyObject *python_internal(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_facing(args, kwargs, FacingCheck::Mode::Internal);
}

PyObject *python_external(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_facing(args, kwargs, FacingCheck::Mode::External);
}

PyObject *python_oo_internal(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_oo_facing(self, args, kwargs, FacingCheck::Mode::Internal);
}

PyObject *python_oo_external(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return python_oo_facing(self, args, kwargs, FacingCheck::Mode::External);
}
