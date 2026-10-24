#pragma once

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
};
