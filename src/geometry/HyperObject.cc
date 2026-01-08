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

#include "geometry/HyperObject.h"
#include "geometry/Geometry.h"
#include "geometry/PolySetUtils.h"
#include "geometry/linalg.h"
#include "utils/printutils.h"
#include "geometry/Grid.h"
#include <algorithm>
#include <sstream>
#include <memory>
#include <Eigen/LU>
#include <cstddef>
#include <string>
#include <vector>

HyperObject::HyperObject(unsigned int dim, boost::tribool convex) : dim_(dim), convex_(convex) {}

std::unique_ptr<Geometry> HyperObject::copy() const { return std::make_unique<HyperObject>(*this); }

std::string HyperObject::dump() const
{
  std::ostringstream out;
  out << "HyperObject:" << "\n dimensions:" << dim_ << "\n convexity:" << this->convexity
      << "\n manifold: " << this->is_manifold_ << "\n num polygons: " << indices.size()
      << "\n polygons data:";
  for (const auto& polygon : indices) {
    out << "\n  polygon begin:";
    for (auto v : polygon) {
      out << "\n   vertex:" << this->vertices[v].transpose();
    }
  }
  out << "\nHyperObject end";
  return out.str();
}

BoundingBox HyperObject::getBoundingBox() const
{
  if (bbox_.isNull()) {
    for (const auto& v : vertices) {
      //      bbox_.extend(v);
    }
  }
  return bbox_;
}

size_t HyperObject::memsize() const
{
  size_t mem = 0;
  for (const auto& p : this->indices) mem += p.size() * sizeof(int);
  for (const auto& p : this->vertices) mem += p.size() * sizeof(Vector4d);
  mem += sizeof(HyperObject);
  return mem;
}
void HyperObject::transform(const Transform3d& mat)
{
  // If mirroring transform, flip faces to avoid the object to end up being inside-out
  bool mirrored = mat.matrix().determinant() < 0;

  for (auto& v : this->vertices) v = mat * v;

  if (mirrored)
    for (auto& p : this->indices) {
      std::reverse(p.begin(), p.end());
    }
  bbox_.setNull();
}

void HyperObject::setColor(const Color4f& c) {}

bool HyperObject::isConvex() const
{
  if (convex_ || this->isEmpty()) return true;
  if (!convex_) return false;
  const bool is_convex = true;
  convex_ = is_convex;
  return is_convex;
}

void HyperObject::resize(const Vector3d& newsize, const Eigen::Matrix<bool, 3, 1>& autosize)
{
  //  this->transform(GeometryUtils::getResizeTransform(this->getBoundingBox(), newsize, autosize));
}
