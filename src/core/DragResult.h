#pragma once

#include <cmath>
#include <string>
#include <vector>

#include "geometry/linalg.h"

/*! One parameter change caused by an interactive drag, e.g. cube size[1] = 12.5. */
class DragMod
{
public:
  std::string name;          ///< parameter name, e.g. "size"
  int index = 0;             ///< positional index of the parameter in the call
  std::vector<int> arrinfo;  ///< path into nested lists, e.g. {1} for size[1]
  double value = 0;          ///< new value
};

/*! Result of AbstractNode::dragPoint(): which call to patch and how. */
class DragResult
{
public:
  Vector3d anchor;            ///< new position of the grabbed point (for the marker)
  std::string modname;        ///< name of the call to patch, e.g. "cube"
  std::vector<DragMod> mods;  ///< parameters to change
  /// DragMod::index counts the object as first argument (translate(obj, v)).
  /// For the method form obj.translate(v) the patcher subtracts one.
  bool objFirst = false;
  /// Script line of the dragged call (AbstractNode::py_line), 0 if unknown.
  int line = 0;
};

/*
   Geometry helpers for dragging. The mouse defines a ray: it passes through
   newpt (the mouse position unprojected at the depth of the grabbed point)
   with direction dir (the viewing direction at that pixel). Instead of using
   newpt directly, every parameter intersects this ray with the line or plane
   along which it can actually move the grabbed point.
 */

/// Intersection of the ray p + t*d with the plane through q with normal n.
/// Returns false when the ray is (nearly) parallel to the plane.
inline bool dragRayPlane(const Vector3d& p, const Vector3d& d, const Vector3d& q, const Vector3d& n,
                         Vector3d& out)
{
  const double dn = d.norm(), nn = n.norm();
  if (dn == 0 || nn == 0) return false;
  const double den = d.dot(n);
  if (std::fabs(den) < 0.05 * dn * nn) return false;  // within ~3 degrees of edge-on
  out = p + d * ((q - p).dot(n) / den);
  return true;
}

/// Point on the line q + s*a which is closest to the ray p + t*d.
/// Returns false when ray and line are (nearly) parallel.
inline bool dragRayLine(const Vector3d& p, const Vector3d& d, const Vector3d& q, const Vector3d& a,
                        Vector3d& out)
{
  if (d.norm() == 0 || a.norm() == 0) return false;
  const Vector3d u = a.normalized(), v = d.normalized(), w = q - p;
  const double b = u.dot(v);
  const double den = 1 - b * b;
  if (den < 0.0025) return false;  // within ~3 degrees of parallel
  const double s = (b * v.dot(w) - u.dot(w)) / den;
  out = q + u * s;
  return true;
}

/// New position of a point pt that may only move along the coordinate axes in
/// mask (bit 0..2 = x,y,z): one axis -> along a line, two axes -> in a plane,
/// otherwise in the screen plane (newpt). Falls back to newpt when degenerate.
inline Vector3d dragConstrain(const Vector3d& pt, const Vector3d& newpt, const Vector3d& dir, int mask)
{
  mask &= 7;
  const int count = (mask & 1) + ((mask >> 1) & 1) + ((mask >> 2) & 1);
  Vector3d out = newpt;
  if (count == 1) {
    Vector3d axis(mask & 1 ? 1 : 0, mask & 2 ? 1 : 0, mask & 4 ? 1 : 0);
    if (!dragRayLine(newpt, dir, pt, axis, out)) out = newpt;
  } else if (count == 2) {
    Vector3d normal(mask & 1 ? 0 : 1, mask & 2 ? 0 : 1, mask & 4 ? 0 : 1);
    if (!dragRayPlane(newpt, dir, pt, normal, out)) out = newpt;
  }
  return out;
}

/// Signed angle in degrees rotating vector a into vector b around axis
/// (right hand rule). Both vectors are projected into the plane normal to axis.
/// Returns NAN if one of them lies on the axis.
inline double dragAngleAbout(const Vector3d& axis, const Vector3d& a, const Vector3d& b)
{
  const Vector3d n = axis.normalized();
  const Vector3d pa = a - n * n.dot(a), pb = b - n * n.dot(b);
  if (pa.norm() < 1e-9 || pb.norm() < 1e-9) return NAN;
  return std::atan2(n.dot(pa.cross(pb)), pa.dot(pb)) * 180 / M_PI;
}
