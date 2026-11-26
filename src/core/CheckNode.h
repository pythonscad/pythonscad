#pragma once

#include <string>

#include "node.h"
#include "geometry/facing_check.h"
#include "geometry/slope_check.h"
#include "geometry/select_check.h"

// One node for all design rule checks. children[0] is the checked object;
// External may have a second child to restrict the check to the gaps
// between the two. The node's geometry is the error solid of the check.
class CheckNode : public AbstractPolyNode
{
public:
  enum class Type {
    Internal,  // min wall thickness            (geometry/facing_check.h)
    External,  // min gap                       (geometry/facing_check.h)
    Slope,     // face angle check (slope/overhang/draft) (geometry/slope_check.h)
    Select,    // spatial relationship filtering (geometry/select_check.h)
  };

  VISITABLE();
  CheckNode(std::shared_ptr<const ModuleInstantiation> mi) : AbstractPolyNode(std::move(mi)) {}
  std::string toString() const override;
  std::string name() const override;

  bool isFacing() const { return type == Type::Internal || type == Type::External; }
  FacingCheck::Mode facingMode() const
  {
    return type == Type::Internal ? FacingCheck::Mode::Internal : FacingCheck::Mode::External;
  }

  Type type = Type::Internal;
  double grow = -1;  // display offset of the error solid, < 0: automatic

  // internal / external
  double distance = 0;
  double min_angle = 120;
  double alpha = 90;
  bool occlusion = true;

  // slope / overhang / draft
  SlopeCheck::Options slope;

  // select
  SelectCheck::Relation select_relation = SelectCheck::Relation::Inside;
};
