#pragma once

#include <string>

#include "node.h"
#include "geometry/facing_check.h"

// internal(obj, d) / external(obj, d[, other]): error solid of a width or
// spacing check between facing surfaces (see geometry/facing_check.h).
// children[0] is the checked object; for External a second child restricts
// the check to gaps between the two.
class FacingCheckNode : public AbstractPolyNode
{
public:
  VISITABLE();
  FacingCheckNode(std::shared_ptr<const ModuleInstantiation> mi) : AbstractPolyNode(std::move(mi)) {}
  std::string toString() const override;
  std::string name() const override
  {
    return mode == FacingCheck::Mode::Internal ? "internal" : "external";
  }

  FacingCheck::Mode mode = FacingCheck::Mode::Internal;
  double distance = 0;
  double min_angle = 120;
  double alpha = 90;
  bool occlusion = true;
};
