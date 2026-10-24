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

#include "core/TransformNode.h"

#include <algorithm>
#include <boost/assign/std/vector.hpp>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <memory>
#include <sstream>
#include <utility>
#include <vector>

#include "core/Builtins.h"
#include "core/Children.h"
#include "core/ModuleInstantiation.h"
#include "core/Parameters.h"
#include "core/Value.h"
#include "core/module.h"
#include "geometry/linalg.h"
#include "utils/degree_trig.h"
#include "utils/printutils.h"
using namespace boost::assign;  // bring 'operator+=()' into scope

enum class transform_type_e { SCALE, ROTATE, MIRROR, TRANSLATE, MULTMATRIX };

std::shared_ptr<AbstractNode> builtin_scale(const std::shared_ptr<const ModuleInstantiation>& inst,
                                            Arguments arguments, const Children& children)
{
  auto node = std::make_shared<TransformNode>(inst, "scale");

  Parameters parameters = Parameters::parse(std::move(arguments), inst->location(), {"v"});

  Vector3d scalevec(1, 1, 1);
  if (!parameters["v"].getVec3(scalevec[0], scalevec[1], scalevec[2], 1.0)) {
    double num;
    if (parameters["v"].getDouble(num)) {
      scalevec.setConstant(num);
    } else {
      LOG(message_group::Warning, inst->location(), parameters.documentRoot(),
          "Unable to convert scale(%1$s) parameter to a number or a vec3 or vec2 of numbers",
          parameters["v"].toEchoStringNoThrow());
    }
  }
  if (OpenSCAD::rangeCheck) {
    if (scalevec[0] == 0 || scalevec[1] == 0 || scalevec[2] == 0 || !std::isfinite(scalevec[0]) ||
        !std::isfinite(scalevec[1]) || !std::isfinite(scalevec[2])) {
      LOG(message_group::Warning, inst->location(), parameters.documentRoot(), "scale(%1$s)",
          parameters["v"].toEchoStringNoThrow());
    }
  }
  node->matrix.scale(scalevec);

  return children.instantiate(node);
}

std::shared_ptr<AbstractNode> builtin_rotate(const std::shared_ptr<const ModuleInstantiation>& inst,
                                             Arguments arguments, const Children& children)
{
  auto node = std::make_shared<TransformNode>(inst, "rotate");

  Parameters parameters = Parameters::parse(std::move(arguments), inst->location(), {"a", "v"});

  const auto& val_a = parameters["a"];
  const auto& val_v = parameters["v"];
  if (val_a.type() == Value::Type::VECTOR) {
    double sx = 0, sy = 0, sz = 0;
    double cx = 1, cy = 1, cz = 1;
    double a = 0.0;
    bool ok = true;
    const auto& vec_a = val_a.toVector();
    switch (vec_a.size()) {
    default: ok &= false; [[fallthrough]];
    case 3:
      ok &= vec_a[2].getDouble(a);
      ok &= !std::isinf(a) && !std::isnan(a);
      sz = sin_degrees(a);
      cz = cos_degrees(a);
      [[fallthrough]];
    case 2:
      ok &= vec_a[1].getDouble(a);
      ok &= !std::isinf(a) && !std::isnan(a);
      sy = sin_degrees(a);
      cy = cos_degrees(a);
      [[fallthrough]];
    case 1:
      ok &= vec_a[0].getDouble(a);
      ok &= !std::isinf(a) && !std::isnan(a);
      sx = sin_degrees(a);
      cx = cos_degrees(a);
      break;
    case 0: break;
    }

    bool v_supplied = val_v.isDefined();
    if (ok) {
      if (v_supplied) {
        LOG(message_group::Warning, inst->location(), parameters.documentRoot(),
            "When parameter a is supplied as vector, v is ignored rotate(a=%1$s, v=%2$s)",
            val_a.toEchoStringNoThrow(), val_v.toEchoStringNoThrow());
      }
    } else {
      if (v_supplied) {
        LOG(message_group::Warning, inst->location(), parameters.documentRoot(),
            "Problem converting rotate(a=%1$s, v=%2$s) parameter", val_a.toEchoStringNoThrow(),
            val_v.toEchoStringNoThrow());
      } else {
        LOG(message_group::Warning, inst->location(), parameters.documentRoot(),
            "Problem converting rotate(a=%1$s) parameter", val_a.toEchoStringNoThrow());
      }
    }
    Matrix3d M;
    M << cy * cz, cz * sx * sy - cx * sz, cx * cz * sy + sx * sz, cy * sz, cx * cz + sx * sy * sz,
      -cz * sx + cx * sy * sz, -sy, cy * sx, cx * cy;
    node->matrix.rotate(M);
  } else {
    double a = 0.0;
    bool aConverted = val_a.getDouble(a);
    aConverted &= !std::isinf(a) && !std::isnan(a);

    Vector3d v(0, 0, 1);
    bool vConverted = val_v.getVec3(v[0], v[1], v[2], 0.0);
    node->matrix.rotate(angle_axis_degrees(aConverted ? a : 0, v));
    if (val_v.isDefined() && !vConverted) {
      if (aConverted) {
        LOG(message_group::Warning, inst->location(), parameters.documentRoot(),
            "Problem converting rotate(..., v=%1$s) parameter", val_v.toEchoStringNoThrow());
      } else {
        LOG(message_group::Warning, inst->location(), parameters.documentRoot(),
            "Problem converting rotate(a=%1$s, v=%2$s) parameter", val_a.toEchoStringNoThrow(),
            val_v.toEchoStringNoThrow());
      }
    } else if (!aConverted) {
      LOG(message_group::Warning, inst->location(), parameters.documentRoot(),
          "Problem converting rotate(a=%1$s) parameter", val_a.toEchoStringNoThrow());
    }
  }

  return children.instantiate(node);
}

std::shared_ptr<AbstractNode> builtin_mirror(const std::shared_ptr<const ModuleInstantiation>& inst,
                                             Arguments arguments, const Children& children)
{
  auto node = std::make_shared<TransformNode>(inst, "mirror");

  Parameters parameters = Parameters::parse(std::move(arguments), inst->location(), {"v"});

  double x = 1.0, y = 0.0, z = 0.0;
  if (!parameters["v"].getVec3(x, y, z, 0.0)) {
    LOG(message_group::Warning, inst->location(), parameters.documentRoot(),
        "Unable to convert mirror(%1$s) parameter to a vec3 or vec2 of numbers",
        parameters["v"].toEchoStringNoThrow());
  }

  // x /= sqrt(x*x + y*y + z*z)
  // y /= sqrt(x*x + y*y + z*z)
  // z /= sqrt(x*x + y*y + z*z)
  if (x != 0.0 || y != 0.0 || z != 0.0) {
    // skip using sqrt to normalize the vector since each element of matrix contributes it with two
    // multiplied terms instead just divide directly within each matrix element simplified calculation
    // leads to less float errors
    double a = x * x + y * y + z * z;

    Matrix4d m;
    m << 1 - 2 * x * x / a, -2 * y * x / a, -2 * z * x / a, 0, -2 * x * y / a, 1 - 2 * y * y / a,
      -2 * z * y / a, 0, -2 * x * z / a, -2 * y * z / a, 1 - 2 * z * z / a, 0, 0, 0, 0, 1;
    node->matrix = m;
  }

  return children.instantiate(node);
}

std::shared_ptr<AbstractNode> builtin_translate(const std::shared_ptr<const ModuleInstantiation>& inst,
                                                Arguments arguments, const Children& children)
{
  auto node = std::make_shared<TransformNode>(inst, "translate");

  Parameters parameters = Parameters::parse(std::move(arguments), inst->location(), {"v"});

  Vector3d translatevec(0, 0, 0);
  bool ok = parameters["v"].getVec3(translatevec[0], translatevec[1], translatevec[2], 0.0);
  ok &=
    std::isfinite(translatevec[0]) && std::isfinite(translatevec[1]) && std::isfinite(translatevec[2]);
  if (ok) {
    node->matrix.translate(translatevec);
  } else {
    LOG(message_group::Warning, inst->location(), parameters.documentRoot(),
        "Unable to convert translate(%1$s) parameter to a vec3 or vec2 of numbers",
        parameters["v"].toEchoStringNoThrow());
  }

  return children.instantiate(node);
}

std::shared_ptr<AbstractNode> builtin_multmatrix(const std::shared_ptr<const ModuleInstantiation>& inst,
                                                 Arguments arguments, const Children& children)
{
  auto node = std::make_shared<TransformNode>(inst, "multmatrix");

  Parameters parameters = Parameters::parse(std::move(arguments), inst->location(), {"m"});

  if (parameters["m"].type() == Value::Type::VECTOR) {
    Matrix4d rawmatrix{Matrix4d::Identity()};
    const auto& mat = parameters["m"].toVector();
    for (size_t row_i = 0; row_i < std::min(mat.size(), size_t(4)); ++row_i) {
      const auto& row = mat[row_i].toVector();
      for (size_t col_i = 0; col_i < std::min(row.size(), size_t(4)); ++col_i) {
        row[col_i].getDouble(rawmatrix(row_i, col_i));
      }
    }
    double w = rawmatrix(3, 3);
    if (w != 1.0) node->matrix = rawmatrix / w;
    else node->matrix = rawmatrix;
  }

  return children.instantiate(node);
}

std::string TransformNode::toString() const
{
  std::ostringstream stream;

  stream << "multmatrix([";
  for (int j = 0; j < 4; ++j) {
    stream << "[";
    for (int i = 0; i < 4; ++i) {
      stream << this->matrix(j, i);
      if (i != 3) stream << ", ";
    }
    stream << "]";
    if (j != 3) stream << ", ";
  }
  stream << "])";

  return stream.str();
}

TransformNode::TransformNode(std::shared_ptr<const ModuleInstantiation> mi, std::string verbose_name)
  : AbstractNode(std::move(mi)), matrix(Transform3d::Identity()), _name(std::move(verbose_name))
{
}

void TransformNode::dragBegin()
{
  matrix_ = matrix;
  drag_init = true;
  drag_axis = -1;
  drag_prev_angle = NAN;
  AbstractNode::dragBegin();
}

namespace {

constexpr double DRAG_TOL = 1e-3;
constexpr double DRAG_SNAP = 0.1;  // translate grid
constexpr double DRAG_SNAP_SCALE = 0.01;
constexpr double DRAG_SNAP_ANGLE = 1;  // degrees

DragMod transformMod(const std::string& name, std::vector<int> arrinfo, double value)
{
  DragMod mod;
  mod.name = name;
  mod.index = 1;  // translate(obj, v), rotate(obj, a, ...), scale(obj, v)
  mod.arrinfo = std::move(arrinfo);
  mod.value = value;
  return mod;
}

Matrix3d eulerMatrix(const Vector3d& e)
{
  return (Eigen::AngleAxisd(e[2] * M_PI / 180, Vector3d::UnitZ()) *
          Eigen::AngleAxisd(e[1] * M_PI / 180, Vector3d::UnitY()) *
          Eigen::AngleAxisd(e[0] * M_PI / 180, Vector3d::UnitX()))
    .toRotationMatrix();
}

/// Choose the angle closest to prev (360 degree periods), so a drag can turn
/// past +-180 degrees and even several times around without jumping.
double unwrapAngle(double angle, double prev)
{
  if (std::isnan(prev)) return angle;
  return angle + 360 * std::round((prev - angle) / 360);
}

}  // namespace

/*!
  Drag of a point which belongs to a child of this transform.

  translate: the point moves along the marked axes (line, plane or screen plane).
  scale:     like translate; the factor is new position / old position per axis.
  rotate:    the point moves on a circle around the rotation axis. The mouse ray
             is intersected with the plane of that circle and the angle is
             measured from the point's position at angle 0. This is absolute
             (no accumulation of errors); only the previous angle is kept to
             continue smoothly beyond +-180 degrees.
 */
std::shared_ptr<const Geometry> TransformNode::dragPoint(const Vector3d& pt, const Vector3d& newpt,
                                                         const Vector3d& dir, DragResult& result)
{
  if (!drag_init) dragBegin();
  const Transform3d invmat = matrix_.inverse();
  const Vector3d pt_tran = invmat * pt;
  const Vector3d newpt_tran = invmat * newpt;
  const Vector3d dir_tran = invmat.linear() * dir;

  // The innermost draggable node wins: give the children the first chance.
  auto result_geom = std::make_shared<PolySet>(3);
  bool child_claimed = false, found = false;
  for (auto& child : children) {
    const bool fresh = std::isnan(result.anchor[0]);
    auto child_geom = child->dragPoint(pt_tran, newpt_tran, dir_tran, result);
    if (fresh && !std::isnan(result.anchor[0])) {
      result.anchor = matrix * result.anchor;
      child_claimed = true;
    }
    auto ps = std::dynamic_pointer_cast<const PolySet>(child_geom);
    if (ps == nullptr) continue;
    for (const auto& v : ps->vertices) {
      if (!found && (v - pt_tran).norm() < DRAG_TOL) found = true;
      result_geom->vertices.push_back(matrix_ * v);
    }
  }
  if (child_claimed || !found || !(dragflags & 15)) return result_geom;

  std::vector<DragMod> mods;
  if (_name == "translate") {
    // in parent coordinates the translation simply adds
    const Vector3d target = dragConstrain(pt, newpt, dir, dragflags);
    const Vector3d v0 = matrix_.translation();
    Vector3d v = v0;
    for (int k = 0; k < 3; k++) {
      if (!(dragflags & (1 << k))) continue;
      v[k] = dragSnap(v0[k] + target[k] - pt[k], DRAG_SNAP);
      mods.push_back(transformMod("v", {k}, v[k]));
    }
    matrix = matrix_;
    matrix.translation() = v;
  } else if (_name == "scale") {
    const Vector3d s0 = matrix_.linear().diagonal();
    Vector3d s = s0;
    if (dragflags & 8) {  // scale(marked(2)): uniform factor along the line origin -> point
      Vector3d target;
      if (!dragRayLine(newpt, dir, Vector3d::Zero(), pt, target)) target = newpt;
      if (pt.norm() > DRAG_TOL) {
        const double f = dragSnap(s0[0] * target.dot(pt) / pt.dot(pt), DRAG_SNAP_SCALE);
        if (f != 0) {
          s = Vector3d(f, f, f);
          mods.push_back(transformMod("v", {}, f));
        }
      }
    } else {
      int mask = 0;  // axes which can change: marked and point not on the plane
      for (int k = 0; k < 3; k++)
        if ((dragflags & (1 << k)) && fabs(pt_tran[k]) > DRAG_TOL) mask |= 1 << k;
      const Vector3d target = dragConstrain(pt, newpt, dir, mask);
      for (int k = 0; k < 3; k++) {
        if (!(mask & (1 << k))) continue;
        const double f = dragSnap(target[k] / pt_tran[k], DRAG_SNAP_SCALE);
        if (f == 0) continue;
        s[k] = f;
        mods.push_back(transformMod("v", {k}, f));
      }
    }
    matrix = matrix_;
    matrix.linear() = s.asDiagonal();
  } else if (_name == "rotate") {
    Vector3d axis_parent;  // rotation axis in parent coordinates
    Vector3d zero_parent;  // grabbed point at angle 0 of the dragged axis
    Matrix3d after = Matrix3d::Identity();
    double angle0 = NAN;
    if (!std::isnan(rot_euler[0]) && (dragflags & 7)) {
      if (drag_axis < 0) {
        // several marked angles: use the axis we look along most, its circle is
        // seen most face-on and is the easiest to follow with the mouse
        double best = -1;
        for (int k = 0; k < 3; k++) {
          if (!(dragflags & (1 << k))) continue;
          Matrix3d a = Matrix3d::Identity();
          if (k == 0) a = eulerMatrix(Vector3d(0, rot_euler[1], rot_euler[2]));
          if (k == 1) a = eulerMatrix(Vector3d(0, 0, rot_euler[2]));
          const double score = fabs((a * Vector3d::Unit(k)).normalized().dot(dir.normalized()));
          if (score > best) {
            best = score;
            drag_axis = k;
          }
        }
      }
      const int k = drag_axis;
      Vector3d before_e(0, 0, 0), after_e(0, 0, 0);
      for (int j = 0; j < k; j++) before_e[j] = rot_euler[j];
      for (int j = k + 1; j < 3; j++) after_e[j] = rot_euler[j];
      after = eulerMatrix(after_e);
      axis_parent = after * Vector3d::Unit(k);
      zero_parent = after * (eulerMatrix(before_e) * pt_tran);
      angle0 = rot_euler[k];
    } else if (!std::isnan(rot_angle) && (dragflags & 8)) {
      axis_parent = rot_axis.normalized();
      zero_parent = pt_tran;
      angle0 = rot_angle;
    }
    if (!std::isnan(angle0)) {
      // The circle lies in the plane through the grabbed point, normal to the axis.
      Vector3d target;
      if (!dragRayPlane(newpt, dir, pt, axis_parent, target)) {
        target = newpt - axis_parent.normalized() * axis_parent.normalized().dot(newpt - pt);
      }
      double angle = dragAngleAbout(axis_parent, zero_parent, target);
      if (!std::isnan(angle)) {
        angle = unwrapAngle(angle, std::isnan(drag_prev_angle) ? angle0 : drag_prev_angle);
        drag_prev_angle = angle;
        angle = dragSnap(angle, DRAG_SNAP_ANGLE);
        Matrix3d M;
        if (!std::isnan(rot_euler[0]) && (dragflags & 7)) {
          Vector3d e = rot_euler;
          e[drag_axis] = angle;
          M = eulerMatrix(e);
          mods.push_back(transformMod("a", {drag_axis}, angle));
        } else {
          M = Eigen::AngleAxisd(angle * M_PI / 180, rot_axis.normalized()).toRotationMatrix();
          mods.push_back(transformMod("a", {}, angle));
        }
        matrix = matrix_;
        matrix.linear() = M;
      }
    }
  }
  if (mods.empty()) return result_geom;

  result.modname = _name;
  result.mods = mods;
  result.objFirst = true;
  result.anchor = matrix * pt_tran;
  return result_geom;
}

std::string TransformNode::name() const
{
  return "transform";
}

std::string TransformNode::verbose_name() const
{
  return _name;
}

void register_builtin_transform()
{
  Builtins::init("scale", new BuiltinModule(builtin_scale),
                 {
                   "scale([x, y, z])",
                 });

  Builtins::init("rotate", new BuiltinModule(builtin_rotate),
                 {
                   "rotate([x, y, z])",
                 });

  Builtins::init("mirror", new BuiltinModule(builtin_mirror),
                 {
                   "mirror([x, y, z])",
                 });

  Builtins::init("translate", new BuiltinModule(builtin_translate),
                 {
                   "translate([x, y, z])",
                 });

  Builtins::init("multmatrix", new BuiltinModule(builtin_multmatrix),
                 {
                   "multmatrix(matrix_4_by_4)",
                 });
}
