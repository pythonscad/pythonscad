#pragma once

#include "geometry/Geometry.h"
#include "geometry/linalg.h"
#include "geometry/GeometryUtils.h"
#include "geometry/Polygon2d.h"
#include "geometry/Curve.h"
#include "geometry/Surface.h"
#include "utils/boost-utils.h"

#include <cstdint>
#include <memory>
#include <cstddef>
#include <string>
#include <vector>

class HyperObject : public Geometry
{
public:
  VISITABLE_GEOMETRY();
  PolygonIndices indices;
  std::vector<Vector4d> vertices;

  HyperObject(unsigned int dim, boost::tribool convex = unknown);

  size_t memsize() const override;
  BoundingBox getBoundingBox() const override;
  std::string dump() const override;
  unsigned int getDimension() const override { return dim_; }
  bool isEmpty() const override { return indices.empty(); }
  std::unique_ptr<Geometry> copy() const override;

  void quantizeVertices(std::vector<Vector3d> *pPointsOut = nullptr);
  size_t numFacets() const override { return indices.size(); }
  void transform(const Transform3d& mat) override;
  void resize(const Vector3d& newsize, const Eigen::Matrix<bool, 3, 1>& autosize) override;
  void setColor(const Color4f& c) override;

  bool isConvex() const;
  boost::tribool convexValue() const { return convex_; }

  bool isManifold() const { return is_manifold_; }
  void setManifold(bool manifold) { is_manifold_ = manifold; }

  bool isTriangular() const { return is_triangular_; }
  void setTriangular(bool triangular) { is_triangular_ = triangular; }

  static std::unique_ptr<HyperObject> createEmpty() { return std::make_unique<HyperObject>(4); }

private:
  bool is_triangular_ = false;
  unsigned int dim_;
  mutable boost::tribool convex_;
  mutable BoundingBox bbox_;

  bool is_manifold_ = false;  // false means "unknown"
};
