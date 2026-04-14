#include "node.h"
#include "Value.h"
#include <src/geometry/linalg.h>
#include "src/geometry/PolySet.h"
#include "utils/png_util.h"

enum textureProjections {
  PROJECTION_NONE,
  TRIPLANAR,
  CUBIC,
  SPHERICAL,
  CYLINDRIC,
  PLANARX,
  PLANARY,
  PLANARZ,
  TEXTUREPROJECTION_NUM
};
extern const char *projectionNames[];

class BaseProjection
{
public:
  BaseProjection(const img_data_t& tex, double w, double h, double dep)
  {
    texture = tex;
    width = w;
    height = h;
    depth = dep;
  }
  img_data_t texture;
  double width, height, depth;
  double tcoord(double x, double y) const;
  virtual Vector3d calcDisplacement(const Vector3d& pt, const Vector3d& n) = 0;
};

class TriPlanarProjection : public BaseProjection
{
public:
  TriPlanarProjection(const img_data_t& tex, double width, double height, double depth)
    : BaseProjection(tex, width, height, depth)
  {
  }
  virtual Vector3d calcDisplacement(const Vector3d& pt, const Vector3d& n);
};

class CubicProjection : public BaseProjection
{
public:
  CubicProjection(const img_data_t& tex, double width, double height, double depth)
    : BaseProjection(tex, width, height, depth)
  {
  }
  virtual Vector3d calcDisplacement(const Vector3d& pt, const Vector3d& n);
};

class SphericalProjection : public BaseProjection
{
public:
  SphericalProjection(const img_data_t& tex, double width, double height, double depth, Vector3d c)
    : BaseProjection(tex, width, height, depth)
  {
    center = c;
  }
  virtual Vector3d calcDisplacement(const Vector3d& pt, const Vector3d& n);
  Vector3d center;
};

class CylindricProjection : public BaseProjection
{
public:
  CylindricProjection(const img_data_t& tex, double width, double height, double depth, Vector3d c)
    : BaseProjection(tex, width, height, depth)
  {
    center = c;
  }
  virtual Vector3d calcDisplacement(const Vector3d& pt, const Vector3d& n);
  Vector3d center;
};

class PlanarXProjection : public BaseProjection
{
public:
  PlanarXProjection(const img_data_t& tex, double width, double height, double depth)
    : BaseProjection(tex, width, height, depth)
  {
  }
  virtual Vector3d calcDisplacement(const Vector3d& pt, const Vector3d& n);
};

class PlanarYProjection : public BaseProjection
{
public:
  PlanarYProjection(const img_data_t& tex, double width, double height, double depth)
    : BaseProjection(tex, width, height, depth)
  {
  }
  virtual Vector3d calcDisplacement(const Vector3d& pt, const Vector3d& n);
};

class PlanarZProjection : public BaseProjection
{
public:
  PlanarZProjection(const img_data_t& tex, double width, double height, double depth)
    : BaseProjection(tex, width, height, depth)
  {
  }
  virtual Vector3d calcDisplacement(const Vector3d& pt, const Vector3d& n);
};

class OversampleNode : public LeafNode
{
public:
  OversampleNode(const ModuleInstantiation *mi) : LeafNode(mi) {}
  std::string toString() const override
  {
    std::ostringstream stream;
    stream << "oversample( size = " << size;
    if (texturefilename.size() > 0) {
      stream << ", texture = " << texturefilename
             << ", projection = " << projectionNames[textureprojection]
             << ", texturewidth = " << texturewidth << ", textureheight = " << textureheight
             << ", texturedepth = " << texturedepth;
    }
    stream << ")";
    return stream.str();
  }
  std::string name() const override { return "oversample"; }
  std::unique_ptr<const Geometry> createGeometry() const override;

  double size;  // how fine is the oversampling in units

  std::string texturefilename;
  int textureprojection;
  double texturewidth;
  double textureheight;
  double texturedepth = 0.5;

  std::unique_ptr<const Geometry> createGeometry_sub(const std::shared_ptr<const PolySet>& ps) const;
};
