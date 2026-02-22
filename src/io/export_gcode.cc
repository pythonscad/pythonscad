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
#include "io/export.h"

#include <cassert>
#include <clocale>
#include <cmath>
#include <memory>
#include <ostream>

#include "geometry/Geometry.h"
#include "geometry/linalg.h"
#include "geometry/Polygon2d.h"
#include "geometry/PolySet.h"

void output_gcode_pars(std::ostream& output, int gnum, double x, double y, double feed, double power) {
  static int gnum_cached=-1;
  static double x_cached=NAN;
  static double y_cached=NAN;
  static double feed_cached=NAN;
  static double power_cached=NAN;

  if(gnum == -2) {
    // reinitialize the cache
    gnum_cached=-1;
    x_cached=NAN;
    y_cached=NAN;
    feed_cached=NAN;
    power_cached=NAN;

    return;
  }

  if(gnum != -1 && gnum != gnum_cached) {
    output << "G" << gnum;
    gnum_cached = gnum;
  }
  if(!std::isnan(x) && x != x_cached) {
    output << "X" << x;
    x_cached = x;
  }
  if(!std::isnan(y) && y != y_cached) {
    output << "Y" << y;
    y_cached = y;
  }
  if(!std::isnan(feed) && feed != feed_cached) {
    output << "F" << feed;
    feed_cached = feed;
  }
  if(!std::isnan(power) && power != power_cached) {
    output << "S" << power;
    power_cached = power;
  }

  output << "\r\n";
}

static double color_to_parm(const Color4f color, const double max, const int pos, const int dynamic)
{
  int r,g,b,a;
  double parm;

  color.getRgba(r,g,b,a);
  double power = double(uint(r)) * 4.0;
  double feed  = double(uint(g)<<8) + double(b);
  switch (pos) {
    case 0: // power
      parm = (power <= max) ? power : max;
      break;
    case 1: // feed/speed
      parm = (feed <= max) ? feed : max;
      break;
    default:
      fprintf(stderr, "Internal Error: invalid colar param position.\n");
      return -1;
  }

  return (parm);
}

static void append_gcode(const Polygon2d& poly, std::ostream& output, const ExportInfo& exportInfo)
{
  auto  options = exportInfo.optionsGcode;

  for (const auto& o : poly.outlines()) {
    const Eigen::Vector2d& p0 = o.vertices[0];
    const double laserpower = color_to_parm(o.color, options->laserpower, 0, -1.0);
    const double feedrate   = color_to_parm(o.color, options->feedrate, 1, -1.0);

    output_gcode_pars(output, 0, p0.x(), p0.y(), NAN, NAN);
    output_gcode_pars(output, -1, NAN, NAN, NAN, laserpower);
    int n=o.vertices.size();
    for (unsigned int idx = 1; idx <=  n; ++idx) {
      const Eigen::Vector2d& p = o.vertices[idx%n];
      output_gcode_pars(output, 1, p.x(), p.y(), feedrate, laserpower);
    }
    output_gcode_pars(output, -1, NAN, NAN, NAN, 0);
  }
}

static void append_gcode(const std::shared_ptr<const Geometry>& geom, std::ostream& output,
                       const ExportInfo& exportInfo)
{
  if (const auto geomlist = std::dynamic_pointer_cast<const GeometryList>(geom)) {
    for (const auto& item : geomlist->getChildren()) {
      append_gcode(item.second, output, exportInfo);
    }
  } else if (const auto poly = std::dynamic_pointer_cast<const Polygon2d>(geom)) {
    append_gcode(*poly, output, exportInfo);
  } else if (std::dynamic_pointer_cast<const PolySet>(geom)) {  // NOLINT(bugprone-branch-clone)
    assert(false && "Unsupported file format");
  } else {  // NOLINT(bugprone-branch-clone)
    assert(false && "Export as SVG for this geometry type is not supported");
  }
}

void export_gcode(const std::shared_ptr<const Geometry>& geom, std::ostream& output,
                const ExportInfo& exportInfo)
{
  setlocale(LC_NUMERIC, "C");  // Ensure radix is . (not ,) in output
  BoundingBox bbox = geom->getBoundingBox();

  // reset the cached parameters
  output_gcode_pars(output, -2, NAN, NAN, NAN, NAN);

  auto  options = exportInfo.optionsGcode;

  // read in the configuration file
  auto  configfile = options->configfile;
  std::cout << "***** ConfigFile: " << configfile << std::endl;

  output << options->initCode << "\r\n";
  if(options->lasermode == 1) {
    output	<< "M4 S0\r\n";
  } else {
    output	<< "M3 S0\r\n";
  }
  append_gcode(geom, output, exportInfo);
  output	<< "M5 S0\r\n";
  output << options->exitCode;
  setlocale(LC_NUMERIC, "");  // Set default locale
}
