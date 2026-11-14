#pragma once

#if defined(ENABLE_OCCT) && defined(ENABLE_MANIFOLD)

#include <cstddef>
#include <memory>

#include "geometry/brep/BRepGeometry.h"
#include "geometry/manifold/ManifoldGeometry.h"

/**
 * A ManifoldGeometry (typically the result of a boolean) that additionally
 * carries the exact solid of the same operation. See BRepGeometry.
 * ManifoldGeometry can only be changed through its methods, all of which go
 * through transform(), so no fingerprint is needed here.
 */
class BRepManifoldGeometry : public ManifoldGeometry, public ExactCarrier
{
public:
  BRepManifoldGeometry(const ManifoldGeometry& mesh, std::shared_ptr<const ExactShape> exact)
    : ManifoldGeometry(mesh), exact_(std::move(exact))
  {
  }

  [[nodiscard]] std::unique_ptr<Geometry> copy() const override
  {
    return std::make_unique<BRepManifoldGeometry>(*this);
  }
  [[nodiscard]] size_t memsize() const override { return ManifoldGeometry::memsize() + sizeof(exact_); }
  // resize() is inherited and routed through transform()
  void transform(const Transform3d& mat) override
  {
    ManifoldGeometry::transform(mat);
    exact_ = ExactShape::transformed(exact_, mat);
  }

  [[nodiscard]] std::shared_ptr<const ExactShape> exactShape() const override { return exact_; }

private:
  std::shared_ptr<const ExactShape> exact_;
};

#endif
