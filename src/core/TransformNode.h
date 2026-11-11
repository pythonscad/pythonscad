#pragma once

#include <string>

#include "core/ModuleInstantiation.h"
#include "core/node.h"
#include "geometry/linalg.h"
#include "geometry/PolySet.h"

class TransformNode : public AbstractNode
{
public:
  VISITABLE();
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  TransformNode(std::shared_ptr<const ModuleInstantiation> mi, std::string verbose_name);
  std::string toString() const override;
  std::string name() const override;
  std::string verbose_name() const override;
  Transform3d matrix;
  virtual std::shared_ptr<const Geometry> dragPoint(const Vector3d& pt, const Vector3d& delta,
                                                    const Vector3d& dir, DragResult& result) override;
  void dragBegin() override;
  Transform3d matrix_;  // matrix at drag start
  bool drag_init = false;
  // bit 0..2: component x,y,z of the vector is a marked value
  // bit 3: the scalar argument is a marked value (scale(2), rotate(a=45, v=...))
  int dragflags = 0;
  // rotate() only: the arguments the matrix was built from, needed to drag
  // along the real circle instead of linearly
  Vector3d rot_euler{NAN, NAN, NAN};  // rotate([x,y,z]) in degrees
  Vector3d rot_axis{NAN, NAN, NAN};   // rotate(a, v): axis v
  double rot_angle = NAN;             // rotate(a, v): angle a in degrees
  int drag_axis = -1;                 // euler axis chosen for this drag
  double drag_prev_angle = NAN;       // last angle, keeps the angle continuous past +-180

private:
  const std::string _name;
};
