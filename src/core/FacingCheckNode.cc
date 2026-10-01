/*
 *  PythonSCAD
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */

#include "FacingCheckNode.h"

#include <sstream>

std::string FacingCheckNode::toString() const
{
  std::ostringstream stream;
  stream << name() << "(d = " << distance << ", angle = " << min_angle << ", alpha = " << alpha
         << ", occlusion = " << (occlusion ? "true" : "false") << ")";
  return stream.str();
}
