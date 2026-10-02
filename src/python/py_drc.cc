/*
 *  PythonSCAD - Design Rule Checks (DRC)
 *  Unified check() function for internal/external/slope/overhang/draft checks.
 *  All parameter interpretation and type selection happens in the Python layer.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */

#include "genlang/genlang.h"
#include <Python.h>
#include <memory>
#include <cstring>
#include "pyopenscad.h"
#include <Tree.h>
#include <GeometryEvaluator.h>
#include <PolySetUtils.h>
#include <CheckNode.h>
#include "geometry/facing_check.h"
#include "geometry/slope_check.h"
#include "geometry/select_check.h"
#include "pyfunctions.h"

namespace {

std::shared_ptr<const PolySet> evaluateToPolySet(const std::shared_ptr<AbstractNode>& node)
{
  Tree tree(node, "");
  GeometryEvaluator geomevaluator(tree);
  std::shared_ptr<const Geometry> geom = geomevaluator.evaluateGeometry(*tree.root(), true);
  return geom ? PolySetUtils::getGeometryAsPolySet(geom) : nullptr;
}

// Unified check() function: type_str is one of "internal", "external", "slope", "overhang", "draft"
PyObject *python_check(PyObject *obj, const char *type_str, PyObject *other, double distance,
                       double min_deg, double max_deg, double angle_param, double alpha, int occlusion,
                       const SlopeCheck::Options& slope_opt, PyObject *grow_obj, int report)
{
  const bool is_facing =
    std::strcmp(type_str, "internal") == 0 || std::strcmp(type_str, "external") == 0;
  const bool is_internal = std::strcmp(type_str, "internal") == 0;
  const bool is_external = std::strcmp(type_str, "external") == 0;

  PyObject *dummydict = nullptr;
  PyTypeObject *type = PyOpenSCADObjectType(obj);
  std::shared_ptr<AbstractNode> child = PyOpenSCADObjectToNodeMulti(obj, &dummydict);
  auto dummydict_owner = py_owned(dummydict);
  if (child == nullptr) {
    return propagate_or_typeerror("Invalid type for Object in check()\n");
  }

  std::shared_ptr<AbstractNode> child2;
  if (other != nullptr && other != Py_None) {
    PyObject *dummydict2 = nullptr;
    child2 = PyOpenSCADObjectToNodeMulti(other, &dummydict2);
    auto dummydict2_owner = py_owned(dummydict2);
    if (child2 == nullptr) {
      return propagate_or_typeerror("Invalid type for other in check()\n");
    }
  }

  double grow = -1;  // automatic
  if (grow_obj != nullptr && grow_obj != Py_None) {
    grow = PyFloat_AsDouble(grow_obj);
    if (PyErr_Occurred()) return nullptr;
  }

  // Facing checks report
  if (is_facing && report) {
    FacingCheck::Options opt;
    opt.min_angle_deg = alpha;
    opt.alpha_deg = angle_param;
    opt.occlusion = occlusion != 0;
    opt.store_pairs = false;  // summary only, O(N) memory
    const auto ps = evaluateToPolySet(child);
    FacingCheck::Result res;
    if (ps && child2) {
      const auto ps2 = evaluateToPolySet(child2);
      if (ps2) res = FacingCheck::checkBetween(*ps, *ps2, distance, opt);
    } else if (ps) {
      FacingCheck::Mode mode = is_internal ? FacingCheck::Mode::Internal : FacingCheck::Mode::External;
      res = FacingCheck::check(*ps, mode, distance, opt);
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

  // Slope checks report
  if (!is_facing && report) {
    const auto ps = evaluateToPolySet(child);
    SlopeCheck::Result res;
    if (ps) res = SlopeCheck::check(*ps, slope_opt);
    const bool is_overhang = std::strcmp(type_str, "overhang") == 0;
    PyObject *dict = PyDict_New();
    auto dict_owner = py_owned(dict);
    auto set_size = [&](const char *key, size_t v) {
      PyObject *o = PyLong_FromSize_t(v);
      PyDict_SetItemString(dict, key, o);
      Py_DECREF(o);
    };
    set_size("count", res.count);
    set_size("angle", res.angle_count);
    set_size("undercut", res.undercut_count);
    if (res.angle_count) {
      PyObject *o = PyFloat_FromDouble(is_overhang ? -res.worst_deg : res.worst_deg);
      PyDict_SetItemString(dict, "worst", o);
      Py_DECREF(o);
    } else {
      PyDict_SetItemString(dict, "worst", Py_None);
    }
    return dict_owner.release();
  }

  // Node creation for deferred evaluation
  DECLARE_INSTANCE();
  auto node = std::make_shared<CheckNode>(instance);

  if (is_internal) {
    node->type = CheckNode::Type::Internal;
    node->distance = distance;
    node->min_angle = angle_param;
    node->alpha = alpha;
    node->occlusion = occlusion != 0;
  } else if (is_external) {
    node->type = CheckNode::Type::External;
    node->distance = distance;
    node->min_angle = angle_param;
    node->alpha = alpha;
    node->occlusion = occlusion != 0;
  } else {
    // Slope, Overhang, Draft all use CheckNode::Type::Slope
    node->type = CheckNode::Type::Slope;
    node->slope = slope_opt;
  }

  node->grow = grow;
  node->children.push_back(child);
  if (child2) node->children.push_back(child2);
  return PyOpenSCADObjectFromNode(type, node);
}

bool parseDir(PyObject *obj, const char *fname, Vector3d& dir)
{
  if (obj == nullptr || obj == Py_None) return true;
  double x, y, z;
  if (python_vectorval(obj, 3, 3, &x, &y, &z)) {
    PyErr_Format(PyExc_TypeError, "%s(): dir must be a 3D vector", fname);
    return false;
  }
  dir = Vector3d(x, y, z);
  if (dir.norm() == 0) {
    PyErr_Format(PyExc_ValueError, "%s(): dir must not be zero", fname);
    return false;
  }
  return true;
}

bool parseOptionalNumber(PyObject *obj, const char *fname, const char *arg, double& out, bool& given)
{
  given = obj != nullptr && obj != Py_None;
  if (!given) return true;
  out = PyFloat_AsDouble(obj);
  if (PyErr_Occurred()) {
    PyErr_Clear();
    PyErr_Format(PyExc_TypeError, "%s(): %s must be a number", fname, arg);
    return false;
  }
  return true;
}

// overhang: report the overhang angle from the vertical instead of beta
PyObject *python_slope_core(PyObject *obj, const SlopeCheck::Options& opt, PyObject *grow_obj,
                            int report)
{
  double grow = -1;
  bool grow_given;
  if (!parseOptionalNumber(grow_obj, "slope", "grow", grow, grow_given)) return nullptr;
  if (grow_given && grow < 0) {
    PyErr_SetString(PyExc_ValueError, "slope(): grow must be >= 0");
    return nullptr;
  }
  PyObject *dummydict = nullptr;
  PyTypeObject *type = PyOpenSCADObjectType(obj);
  std::shared_ptr<AbstractNode> child = PyOpenSCADObjectToNodeMulti(obj, &dummydict);
  auto dummydict_owner = py_owned(dummydict);
  if (child == nullptr) return propagate_or_typeerror("Invalid type for Object in slope\n");

  if (report) {
    const auto ps = evaluateToPolySet(child);
    SlopeCheck::Result res;
    if (ps) res = SlopeCheck::check(*ps, opt);
    PyObject *dict = PyDict_New();
    auto dict_owner = py_owned(dict);
    auto set_size = [&](const char *key, size_t v) {
      PyObject *o = PyLong_FromSize_t(v);
      PyDict_SetItemString(dict, key, o);
      Py_DECREF(o);
    };
    set_size("count", res.count);
    set_size("angle", res.angle_count);
    set_size("undercut", res.undercut_count);
    if (res.angle_count) {
      PyObject *o = PyFloat_FromDouble(res.worst_deg);
      PyDict_SetItemString(dict, "worst", o);
      Py_DECREF(o);
    } else {
      PyDict_SetItemString(dict, "worst", Py_None);
    }
    return dict_owner.release();
  }

  DECLARE_INSTANCE();
  auto node = std::make_shared<CheckNode>(instance);
  node->type = CheckNode::Type::Slope;
  node->slope = opt;
  node->grow = grow;
  node->children.push_back(child);
  return PyOpenSCADObjectFromNode(type, node);
}

}  // namespace

// Object methods (OO_METHOD_ENTRY calls these)
PyObject *python_oo_internal(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"d", "angle", "alpha", "occlusion", "grow", "report", NULL};
  PyObject *other = nullptr;
  double d = 0, angle = 120, alpha = 90;
  int occlusion = 1, report = 0;
  PyObject *grow = nullptr;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "d|ddpOp", kwlist, &d, &angle, &alpha, &occlusion,
                                   &grow, &report)) {
    PyErr_SetString(PyExc_TypeError, "error during parsing\n");
    return nullptr;
  }
  return python_check(self, "internal", other, d, -90, 90, angle, alpha, occlusion,
                      SlopeCheck::Options(), grow, report);
}

PyObject *python_oo_external(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"d", "other", "angle", "alpha", "occlusion", "grow", "report", NULL};
  PyObject *other = nullptr;
  double d = 0, angle = 120, alpha = 90;
  int occlusion = 1, report = 0;
  PyObject *grow = nullptr;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "d|OddpOp", kwlist, &d, &other, &angle, &alpha,
                                   &occlusion, &grow, &report)) {
    PyErr_SetString(PyExc_TypeError, "error during parsing\n");
    return nullptr;
  }
  return python_check(self, "external", other, d, -90, 90, angle, alpha, occlusion,
                      SlopeCheck::Options(), grow, report);
}

PyObject *python_oo_slope(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"dir", "min", "max", "undercut", "grow", "report", NULL};
  PyObject *dir = nullptr, *mn = nullptr, *mx = nullptr, *grow = nullptr;
  int undercut = 0, report = 0;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|OOOpOp", kwlist, &dir, &mn, &mx, &undercut, &grow,
                                   &report)) {
    PyErr_SetString(PyExc_TypeError, "error during parsing\n");
    return nullptr;
  }
  SlopeCheck::Options opt;
  bool given;
  if (!parseDir(dir, "slope", opt.dir)) return nullptr;
  if (!parseOptionalNumber(mn, "slope", "min", opt.min_deg, given)) return nullptr;
  if (!parseOptionalNumber(mx, "slope", "max", opt.max_deg, given)) return nullptr;
  if (opt.min_deg > opt.max_deg) {
    PyErr_SetString(PyExc_ValueError, "slope(): min must not be larger than max");
    return nullptr;
  }
  opt.undercut = undercut != 0;
  return python_slope_core(self, opt, grow, report);
}

PyObject *python_oo_overhang(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"angle", "dir", "grow", "report", NULL};
  PyObject *dir = nullptr, *grow = nullptr;
  double angle = 45;
  int report = 0;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|dOOp", kwlist, &angle, &dir, &grow, &report)) {
    PyErr_SetString(PyExc_TypeError, "error during parsing\n");
    return nullptr;
  }
  SlopeCheck::Options opt;
  if (!parseDir(dir, "overhang", opt.dir)) return nullptr;
  if (!(angle >= 0 && angle <= 90)) {
    PyErr_SetString(PyExc_ValueError, "overhang(): angle must be in [0, 90]");
    return nullptr;
  }
  opt.min_deg = -angle;  // dir = build direction: beta >= -angle
  opt.skip_base = true;
  return python_slope_core(self, opt, grow, report);
}

PyObject *python_oo_draft(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"angle", "dir", "parting", "undercut", "grow", "report", NULL};
  PyObject *dir = nullptr, *parting = nullptr, *grow = nullptr;
  double angle = 2;
  int undercut = 1, report = 0;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|dOOpOp", kwlist, &angle, &dir, &parting, &undercut,
                                   &grow, &report)) {
    PyErr_SetString(PyExc_TypeError, "error during parsing\n");
    return nullptr;
  }
  SlopeCheck::Options opt;
  bool given;
  if (!parseDir(dir, "draft", opt.dir)) return nullptr;
  if (!(angle >= -90 && angle <= 90)) {
    PyErr_SetString(PyExc_ValueError, "draft(): angle must be in [-90, 90]");
    return nullptr;
  }
  opt.min_deg = angle;
  if (!parseOptionalNumber(parting, "draft", "parting", opt.parting_pos, given)) return nullptr;
  opt.parting = given ? SlopeCheck::Parting::Plane : SlopeCheck::Parting::Free;
  opt.undercut = undercut != 0;
  return python_slope_core(self, opt, grow, report);
}

// Object-oriented select check: obj.select(other, relation)
// relation: "inside", "not_inside", "outside", "not_outside", "straddle", "not_straddle"
PyObject *python_oo_select(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"other", "relation", "report", NULL};
  PyObject *other_obj = nullptr;
  const char *relation_str = "inside";
  int report = 0;

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|sp", kwlist, &other_obj, &relation_str, &report)) {
    PyErr_SetString(PyExc_TypeError, "error during parsing select()\n");
    return nullptr;
  }

  // Parse relation string
  SelectCheck::Relation relation = SelectCheck::Relation::Inside;
  if (std::strcmp(relation_str, "inside") == 0) {
    relation = SelectCheck::Relation::Inside;
  } else if (std::strcmp(relation_str, "not_inside") == 0) {
    relation = SelectCheck::Relation::NotInside;
  } else if (std::strcmp(relation_str, "outside") == 0) {
    relation = SelectCheck::Relation::Outside;
  } else if (std::strcmp(relation_str, "not_outside") == 0) {
    relation = SelectCheck::Relation::NotOutside;
  } else if (std::strcmp(relation_str, "straddle") == 0) {
    relation = SelectCheck::Relation::Straddle;
  } else if (std::strcmp(relation_str, "not_straddle") == 0) {
    relation = SelectCheck::Relation::NotStraddle;
  } else {
    PyErr_SetString(PyExc_ValueError,
                    "select(): invalid relation (must be inside, not_inside, outside, not_outside, "
                    "straddle, or not_straddle)\n");
    return nullptr;
  }

  // Convert self and other to nodes
  PyTypeObject *type = PyOpenSCADObjectType(self);
  PyObject *dummydict = nullptr;
  std::shared_ptr<AbstractNode> self_node = PyOpenSCADObjectToNodeMulti(self, &dummydict);
  auto dummydict_owner = py_owned(dummydict);
  if (self_node == nullptr) {
    return propagate_or_typeerror("Invalid type for Object in select()\n");
  }

  dummydict = nullptr;
  std::shared_ptr<AbstractNode> other_node = PyOpenSCADObjectToNodeMulti(other_obj, &dummydict);
  auto dummydict_owner2 = py_owned(dummydict);
  if (other_node == nullptr) {
    return propagate_or_typeerror("Invalid type for other in select()\n");
  }

  // For report mode, evaluate directly and return stats
  if (report) {
    auto self_ps = evaluateToPolySet(self_node);
    auto other_ps = evaluateToPolySet(other_node);

    if (!self_ps || !other_ps) {
      return propagate_or_typeerror("Unable to evaluate objects to PolySet in select()\n");
    }

    SelectCheck::Options opt;
    opt.relation = relation;
    SelectCheck::Result res = SelectCheck::check(*self_ps, *other_ps, opt);

    PyObject *report_dict = PyDict_New();
    if (!report_dict) return nullptr;
    auto dict_owner = py_owned(report_dict);

    PyObject *count = PyLong_FromSize_t(res.count);
    PyDict_SetItemString(report_dict, "count", count);
    Py_DECREF(count);

    return dict_owner.release();
  }

  // For normal mode, create a CheckNode and return it wrapped as PyOpenSCAD
  auto node = std::make_shared<CheckNode>(nullptr);
  node->type = CheckNode::Type::Select;
  node->select_relation = relation;
  node->children.push_back(self_node);
  node->children.push_back(other_node);
  return PyOpenSCADObjectFromNode(type, node);
}

// Generic check() entry point from Python.
// Parameters are passed exactly as the Python layer provides them.
// type_str: "internal", "external", "slope", "overhang", or "draft"
// Facing checks use: distance, alpha, angle_param, occlusion, other
// Slope checks use: check_options (optional, pre-configured SlopeCheck::Options)
PyObject *python_check(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *kwlist[] = {"obj",         "type_str", "distance",  "min_deg",       "max_deg",
                    "angle_param", "alpha",    "occlusion", "check_options", "other",
                    "grow",        "report",   NULL};
  PyObject *obj = nullptr;
  const char *type_str = nullptr;
  double distance = 0, min_deg = -90, max_deg = 90, angle_param = 120, alpha = 90;
  int occlusion = 1, report = 0;
  PyObject *check_options = nullptr, *other = nullptr, *grow_obj = nullptr;

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Os|dddddiOOOp", kwlist, &obj, &type_str, &distance,
                                   &min_deg, &max_deg, &angle_param, &alpha, &occlusion, &check_options,
                                   &other, &grow_obj, &report)) {
    PyErr_SetString(PyExc_TypeError, "error during parsing check()\n");
    return nullptr;
  }

  // Build SlopeCheck::Options from Python dict (optional for slope checks)
  SlopeCheck::Options slope_opt;
  if (check_options != nullptr && check_options != Py_None) {
    // Extract fields from dict
    PyObject *dir_obj = PyDict_GetItemString(check_options, "dir");
    if (dir_obj && dir_obj != Py_None) {
      double x, y, z;
      if (!python_vectorval(dir_obj, 3, 3, &x, &y, &z)) {
        slope_opt.dir = Vector3d(x, y, z);
      }
    }
    PyObject *min_deg_obj = PyDict_GetItemString(check_options, "min_deg");
    if (min_deg_obj) {
      slope_opt.min_deg = PyFloat_AsDouble(min_deg_obj);
    }
    PyObject *max_deg_obj = PyDict_GetItemString(check_options, "max_deg");
    if (max_deg_obj) {
      slope_opt.max_deg = PyFloat_AsDouble(max_deg_obj);
    }
    PyObject *parting_obj = PyDict_GetItemString(check_options, "parting");
    if (parting_obj) {
      long parting_val = PyLong_AsLong(parting_obj);
      if (parting_val == 0) slope_opt.parting = SlopeCheck::Parting::None;
      else if (parting_val == 1) slope_opt.parting = SlopeCheck::Parting::Free;
      else if (parting_val == 2) slope_opt.parting = SlopeCheck::Parting::Plane;
    }
    PyObject *parting_pos_obj = PyDict_GetItemString(check_options, "parting_pos");
    if (parting_pos_obj) {
      slope_opt.parting_pos = PyFloat_AsDouble(parting_pos_obj);
    }
    PyObject *undercut_obj = PyDict_GetItemString(check_options, "undercut");
    if (undercut_obj) {
      slope_opt.undercut = PyObject_IsTrue(undercut_obj) > 0;
    }
    PyObject *skip_base_obj = PyDict_GetItemString(check_options, "skip_base");
    if (skip_base_obj) {
      slope_opt.skip_base = PyObject_IsTrue(skip_base_obj) > 0;
    }
  }

  return python_check(obj, type_str, other, distance, min_deg, max_deg, angle_param, alpha, occlusion,
                      slope_opt, grow_obj, report);
}
