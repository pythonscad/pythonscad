/*
 *  OpenSCAD (www.openscad.org)
 *  Copyright (C) 2009-2011 Clifford Wolf <clifford@clifford.at> and
 *                          Marius Kintel <marius@kintel.net>
 *  Copyright (C) 2021      Konstantin Podsvirov <konstantin@podsvirov.pro>
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

#include "io/export.h"
 
#include <ostream>
#include <memory>

#include "geometry/Geometry.h"
#include "geometry/PolySetUtils.h"
#include "geometry/PolySet.h"
#ifdef ENABLE_PYTHON
#include "python/python_public.h"
#endif

void export_obj(const std::shared_ptr<const Geometry>& geom, std::ostream& output)
{
  // FIXME: In lazy union mode, should we export multiple objects?
  
  std::shared_ptr<const PolySet> out = PolySetUtils::getGeometryAsPolySet(geom);
  if (!out->isTriangular()) {
    // While the OBJ format allows for faces to have more than 3
    // vertices, this seems to confuse a number of applications
    // we care about, so for now this will just always tesselate
    // faces to be composed of triangles only.
    //
    // See: https://github.com/openscad/openscad/issues/5993
    out = PolySetUtils::tessellate_faces(*out);
  }
  if (Feature::ExperimentalPredictibleOutput.is_enabled()) {
    out = createSortedPolySet(*out);
  }

  output << "# OpenSCAD obj exporter\n";

#ifdef ENABLE_PYTHON
  python_export_obj_att(output);
#endif  

  for (const auto &v : out->vertices) {
    output << "v " <<v[0] << " " << v[1] << " " << v[2] << "\n";
  }

  for (const auto& poly : out->indices) {
    output << "f ";
    for (const auto idx : poly) {
      output << " " << idx + 1;
    }
    output << "\n";
  }
}
