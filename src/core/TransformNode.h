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
                                                    DragResult& result) override;
  void dragBegin() override;
  Transform3d matrix_;  // matrix at drag start
  bool drag_init = false;
  int dragflags = 0;  // bit 0..2: x,y,z are draggable

private:
  const std::string _name;
};
