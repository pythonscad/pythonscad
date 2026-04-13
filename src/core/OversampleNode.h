#pragma once

#include "node.h"
#include "Value.h"
#include <src/geometry/linalg.h>
#include "src/geometry/PolySet.h"

class OversampleNode : public LeafNode
{
public:
  OversampleNode(const ModuleInstantiation *mi) : LeafNode(mi) {}
  std::string toString() const override
  {
    std::ostringstream stream;
    stream << "oversample( size = " << size;
    if (texturefilename.size() > 0) {
      stream << ", texture = " << texturefilename << ", texturewidth = " << texturewidth
             << ", textureheight = " << textureheight << ", texturedepth = " << texturedepth;
    }
    stream << ")";
    return stream.str();
  }
  std::string name() const override { return "oversample"; }
  std::unique_ptr<const Geometry> createGeometry() const override;

  double size;  // how fine is the oversampling in units

  std::string texturefilename;
  std::string textureprojection;
  double texturewidth;
  double textureheight;
  double texturedepth = 0.5;

  std::unique_ptr<const Geometry> createGeometry_sub(const std::shared_ptr<const PolySet>& ps) const;
  double tcoord(std::shared_ptr<img_data_t> tex, double u, double v) const;
};
