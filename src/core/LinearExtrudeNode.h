#pragma once

#include <string>

#include "core/node.h"
#include "core/ModuleInstantiation.h"
#include "core/Value.h"
#include "geometry/linalg.h"
#ifdef ENABLE_PYTHON
#include <src/python/python_public.h>
#endif

class LinearExtrudeNode : public AbstractPolyNode
{
public:
  VISITABLE();
  LinearExtrudeNode(const ModuleInstantiation *mi) : AbstractPolyNode(mi) {
  }
  std::string toString() const override;
  std::string name() const override { return "linear_extrude"; }

  Vector3d height=Vector3d(0, 0, 1);
  double origin_x = 0.0, origin_y = 0.0;
  double fn = 0.0, fs = 0.0, fa = 0.0;
  double scale_x = 1.0, scale_y = 1.0;
  double twist = 0.0;
  unsigned int convexity = 1u;
  unsigned int slices = 1u, segments = 0u;
  bool has_twist = false, has_slices = false, has_segments = false, has_heightvector=false;
  bool center = false;
 #ifdef ENABLE_PYTHON
  void *profile_func;
  void *twist_func;
 #endif  
};
