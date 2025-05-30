#include "geometry/Polygon2d.h"

#include <sstream>
#include <utility>
#include <cstddef>
#include <string>
#include <memory>

#include "geometry/Geometry.h"
#include "geometry/linalg.h"
#include "utils/printutils.h"
#include "Feature.h"
#ifdef ENABLE_MANIFOLD
#include "geometry/manifold/manifoldutils.h"
#endif
#include "geometry/cgal/cgalutils.h"
#include "Feature.h"
#include "geometry/PolySet.h"
#include "glview/RenderSettings.h"


Polygon2d::Polygon2d(Outline2d outline) : sanitized(true) {
  addOutline(std::move(outline));
}

std::unique_ptr<Geometry> Polygon2d::copy() const
{
  return std::make_unique<Polygon2d>(*this);
}

BoundingBox Outline2d::getBoundingBox() const {
  BoundingBox bbox;
  // Note: this uses ->outlines, and so automatically gets trans3d applied
  for (const auto& v : this->vertices) {
    bbox.extend(Vector3d(v[0], v[1], 0));
  }
  return bbox;
}

/*!
   Class for holding 2D geometry.

   This class will hold 2D geometry consisting of a number of closed
   polygons. Each polygon can contain holes and islands. Both polygons,
   holes and island contours may intersect each other.

   We can store sanitized vs. unsanitized polygons. Sanitized polygons
   will have opposite winding order for holes and is guaranteed to not
   have intersecting geometry. The winding order will be counter-clockwise
   for positive outlines and clockwise for holes. Sanitization is typically
   done by ClipperUtils, but if you create geometry which you know is sanitized,
   the flag can be set manually.
 */

size_t Polygon2d::memsize() const
{
  size_t mem = 0;
  for (const auto &o : this->theoutlines) {
    mem += o.vertices.size() * sizeof(Vector2d) + sizeof(Outline2d);
  }
  for (const auto &o : this->trans3dOutlines) {
    mem += o.vertices.size() * sizeof(Vector2d) + sizeof(Outline2d);
  }
  mem += sizeof(Polygon2d);
  return mem;
}

BoundingBox Polygon2d::getBoundingBox() const
{
  BoundingBox bbox;
  for (const auto& o : this->outlines()) {
    bbox.extend(o.getBoundingBox());
  }
  return bbox;
}

std::string Polygon2d::dump() const
{
  std::ostringstream out;
  for (const auto& o : this->theoutlines) {
    out << "contour:\n";
    for (const auto& v : o.vertices) {
      out << "  " << v.transpose();
    }
    out << "\n";
  }
  if (trans3dState != Transform3dState::NONE) {
    out << "with trans3d: [";
    for (int j=0;j<4;j++) {
      out << "[";
      for (int i=0;i<4;i++) {
        double v(this->trans3d(j, i));
        out << v;
        if (i != 3) out << ", ";
      }
      out << "]";
      if (j != 3) out << ", ";
    }
    out << "]\n";
  }
  return out.str();
}

bool Polygon2d::isEmpty() const
{
  return this->theoutlines.empty();
}

void Polygon2d::transform(const Transform2d& mat)
{
  if (mat.matrix().determinant() == 0) {
    LOG(message_group::Warning, "Scaling a 2D object with 0 - removing object");
    this->theoutlines.clear();
    trans3dState= Transform3dState::NONE;
    return;
  }
  if (trans3dState != Transform3dState::NONE) mergeTrans3d();
  for (auto& o : this->theoutlines) {
    for (auto& v : o.vertices) {
      v = mat * v;
    }
  }
}

void Polygon2d::resize(const Vector2d& newsize, const Eigen::Matrix<bool, 2, 1>& autosize)
{
  auto bbox = this->getBoundingBox();

  // Find largest dimension
  int maxdim = (newsize[1] && newsize[1] > newsize[0]) ? 1 : 0;

  // Default scale (scale with 1 if the new size is 0)
  Vector2d scale(newsize[0] > 0 ? newsize[0] / bbox.sizes()[0] : 1,
                 newsize[1] > 0 ? newsize[1] / bbox.sizes()[1] : 1);

  // Autoscale where applicable
  double autoscale = newsize[maxdim] > 0 ? newsize[maxdim] / bbox.sizes()[maxdim] : 1;
  Vector2d newscale(!autosize[0] || (newsize[0] > 0) ? scale[0] : autoscale,
                    !autosize[1] || (newsize[1] > 0) ? scale[1] : autoscale);

  Transform2d t;
  t.matrix() <<
    newscale[0], 0, 0,
    0, newscale[1], 0,
    0, 0, 1;

  this->transform(t);
}

bool Polygon2d::is_convex() const
{
  if (theoutlines.size() > 1) return false;
  if (theoutlines.empty()) return true;

  auto const& pts = theoutlines[0].vertices;
  int N = pts.size();

  // Check for a right turn. This assumes the polygon is simple.
  for (int i = 0; i < N; ++i) {
    const auto& d1 = pts[(i + 1) % N] - pts[i];
    const auto& d2 = pts[(i + 2) % N] - pts[(i + 1) % N];
    double zcross = d1[0] * d2[1] - d1[1] * d2[0];
    if (zcross < 0) return false;
  }
  return true;
}

double Polygon2d::area() const
{
  auto ps = tessellate();
  if (ps == nullptr) {
    return 0;
  }

  double area = 0.0;
  for (const auto& poly : ps->indices) {
    const auto& v1 = ps->vertices[poly[0]];
    const auto& v2 = ps->vertices[poly[1]];
    const auto& v3 = ps->vertices[poly[2]];
    area += 0.5 * (
      v1.x() * (v2.y() - v3.y())
      + v2.x() * (v3.y() - v1.y())
      + v3.x() * (v1.y() - v2.y()));
  }
  return area;
}

/*!
   Triangulates this polygon2d and returns a 2D-in-3D PolySet.

   This is used for various purposes:
   * Geometry evaluation for roof, linear_extrude, rotate_extrude
   * Rendering (both preview and render mode)
   * Polygon area calculation
   *
   * One use-case is special: For geometry construction in Manifold mode, we require this function to
   * guarantee that vertices and their order are untouched (apart from adding a zero 3rd dimension)
   *
 */
std::unique_ptr<PolySet> Polygon2d::tessellate(bool in3d) const
{
  PRINTDB("Polygon2d::tessellate(): %d outlines", this->outlines().size());
  std::unique_ptr<PolySet> res;
#if defined(ENABLE_MANIFOLD) && defined(USE_MANIFOLD_TRIANGULATOR)
  if (RenderSettings::inst()->backend3D == RenderBackend3D::ManifoldBackend) {
    res = ManifoldUtils::createTriangulatedPolySetFromPolygon2d(*this, in3d);
  }
  else
#endif
    res = CGALUtils::createTriangulatedPolySetFromPolygon2d(*this, in3d);
  if (in3d)
    res->transform(this->getTransform3d());
  return res;
}

void Polygon2d::transform3d(const Transform3d &mat)
{
  // Check whether it can be a 2d transform, and avoid the 3d overhead
  if (trans3dState == Transform3dState::NONE
    && mat(2,0) == 0 && mat(2,1) == 0 && mat(2,2) == 1 && mat(2,3) == 0
    && mat(0,2) == 0 && mat(1,2) == 0 && mat(3,2) == 0
  ) {
    Transform2d t;
    t.matrix() <<
      mat(0,0), mat(0,1), mat(0,3),
      mat(1,0), mat(1,1), mat(1,3),
      mat(3,0), mat(3,1), mat(3,3);
    if (t.matrix().determinant() == 0) {
      LOG(message_group::Warning,"Scaling a 2D object with 0 - removing object");
      this->theoutlines.clear();
      trans3dState= Transform3dState::NONE;
      return;
    }
    transform(t);
    // A 2D transformation may flip the winding order of a polygon.
    // If that happens with a sanitized polygon, we need to reverse
    // the winding order for it to be correct.
    if (sanitized && t.matrix().determinant() < 0)
      for (auto &o : this->theoutlines)
        std::reverse(o.vertices.begin(), o.vertices.end());
  }
  else {
    if (mat.matrix().determinant() == 0) {
      LOG(message_group::Warning,"Scaling a 2D object with 0 - removing object");
      this->theoutlines.clear();
      trans3dState= Transform3dState::NONE;
      return;
    }
    trans3d= (trans3dState == Transform3dState::NONE)? mat : mat * trans3d;
    trans3dState= Transform3dState::PENDING;
  }
}


// This returns the outlines after applying any Transform3d that might be Transform3dState::PENDING.
// If there is no Transform3d, this returns the outlines vector.
// If there is a Transform3dState::CACHED Transform3d, this uses the cache.
// Else it creates and returns the cache
const Polygon2d::Outlines2d &Polygon2d::transformedOutlines() const {
  if (trans3dState == Transform3dState::NONE) return theoutlines;
  if (trans3dState != Transform3dState::CACHED) {
    // Need to remove const from the cache object.  It maintains proper const semantics to the public API though.
    Polygon2d::Outlines2d &cache= const_cast<Polygon2d::Outlines2d&>(trans3dOutlines);
    cache= theoutlines;
    applyTrans3dToOutlines(cache);
    const_cast<Polygon2d*>(this)->trans3dState= Transform3dState::CACHED;
  }
  return trans3dOutlines;
}

// This flattens the 3D transform into the 2D transform that it would have been
// originally.
void Polygon2d::mergeTrans3d() {
  if (trans3dState == Transform3dState::CACHED)
    theoutlines.swap(trans3dOutlines);
  else if (trans3dState == Transform3dState::PENDING)
    applyTrans3dToOutlines(theoutlines);
  trans3dOutlines.clear();
  trans3dState= Transform3dState::NONE;
}

void Polygon2d::applyTrans3dToOutlines(Polygon2d::Outlines2d &outlines) const {
  Transform2d t;
  t.matrix() <<
    trans3d(0,0), trans3d(0,1), trans3d(0,3),
    trans3d(1,0), trans3d(1,1), trans3d(1,3),
    trans3d(3,0), trans3d(3,1), trans3d(3,3);
  for (auto &o : outlines) {
    for (auto &v : o.vertices) {
      v = t * v;
    }
  }
  // A 2D transformation may flip the winding order of a polygon.
  // If that happens with a sanitized polygon, we need to reverse
  // the winding order for it to be correct.
  if (sanitized && t.matrix().determinant() < 0)
    for (auto &o : outlines)
      std::reverse(o.vertices.begin(), o.vertices.end());
}

void Polygon2d::reverse(void) {
  for (auto &o : theoutlines)
    std::reverse(o.vertices.begin(), o.vertices.end());
  for (auto &o : trans3dOutlines)
    std::reverse(o.vertices.begin(), o.vertices.end());
}
