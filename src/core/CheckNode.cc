/*
 *  PythonSCAD
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */

#include "CheckNode.h"

#include <sstream>

std::string CheckNode::name() const
{
  switch (type) {
  case Type::Internal: return "internal";
  case Type::External: return "external";
  case Type::Slope:    return "slope";
  case Type::Select:   return "select";
  }
  return "check";
}

std::string CheckNode::toString() const
{
  std::ostringstream stream;
  stream << name() << "(";
  if (isFacing()) {
    stream << "d = " << distance << ", angle = " << min_angle << ", alpha = " << alpha
           << ", occlusion = " << (occlusion ? "true" : "false");
  } else if (type == Type::Slope) {
    stream << "dir = [" << slope.dir[0] << ", " << slope.dir[1] << ", " << slope.dir[2]
           << "], min = " << slope.min_deg << ", max = " << slope.max_deg << ", parting = ";
    switch (slope.parting) {
    case SlopeCheck::Parting::None:  stream << "none"; break;
    case SlopeCheck::Parting::Free:  stream << "free"; break;
    case SlopeCheck::Parting::Plane: stream << slope.parting_pos; break;
    }
    stream << ", undercut = " << (slope.undercut ? "true" : "false")
           << ", skip_base = " << (slope.skip_base ? "true" : "false");
  } else if (type == Type::Select) {
    stream << "relation = ";
    switch (select_relation) {
    case SelectCheck::Relation::Inside:      stream << "inside"; break;
    case SelectCheck::Relation::NotInside:   stream << "not_inside"; break;
    case SelectCheck::Relation::Outside:     stream << "outside"; break;
    case SelectCheck::Relation::NotOutside:  stream << "not_outside"; break;
    case SelectCheck::Relation::Straddle:    stream << "straddle"; break;
    case SelectCheck::Relation::NotStraddle: stream << "not_straddle"; break;
    }
  }
  stream << ", grow = " << grow << ")";
  return stream.str();
}
