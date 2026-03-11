/*
 *  OpenSCAD (www.openscad.org)
 *  Copyright (C) 2009-2019 Clifford Wolf <clifford@clifford.at> and
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
// #include "io/import.h"
//
// #include <clipper2/clipper.h>
//
// #include <Eigen/Core>
// #include <Eigen/Geometry>
// #include <exception>
#include <memory>
#include <string>
#include <vector>
//
// #include "core/AST.h"
#include "core/CurveDiscretizer.h"
// #include "geometry/ClipperUtils.h"
#include "geometry/Polygon2d.h"
// #include "libsvg/libsvg.h"
// #include "libsvg/shape.h"
// #include "libsvg/svgpage.h"
// #include "libsvg/util.h"
// #include "utils/printutils.h"
// #include "src/core/ColorUtil.h"
#include <iostream>
#include <fstream>

#include <libcdr-0.1/libcdr/libcdr.h>
#include <librevenge/librevenge.h>

class CDRReader : public librevenge::RVNGDrawingInterface
{
public:
  void startDocument(const librevenge::RVNGPropertyList& propList) override
  {
    printf("start\n");
    //        std::cout << "Start document\n";
  }

  void endDocument() override
  {
    printf("Emnd\n");
    //       std::cout << "End document\n";
  }

  void drawPath(const librevenge::RVNGPropertyList& propList) override
  {
    printf("path\n");
    //        std::cout << "Vector path found\n";
  }

  void drawGraphicObject(const librevenge::RVNGPropertyList& propList) override
  {
    printf("graph\n");
    //       std::cout << "Graphic object found\n";
  }
  void setDocumentMetaData(const librevenge::RVNGPropertyList&) override {}
  void defineEmbeddedFont(const librevenge::RVNGPropertyList& propList) override {}
  void startPage(const librevenge::RVNGPropertyList& propList) override {}
  void startMasterPage(const librevenge::RVNGPropertyList& propList) override {}
  void setStyle(const librevenge::RVNGPropertyList& propList) override {}
  void startLayer(const librevenge::RVNGPropertyList& propList) override {}
  void endLayer() override {}
  void startEmbeddedGraphics(const librevenge::RVNGPropertyList& propList) override {}
  void endEmbeddedGraphics() override {}
  void openGroup(const librevenge::RVNGPropertyList& propList) override {}
  void closeGroup() override {}
  void drawRectangle(const librevenge::RVNGPropertyList& propList) override { printf("rect\n"); }
  void drawEllipse(const librevenge::RVNGPropertyList& propList) override { printf("ellips\n"); }
  void drawPolygon(const librevenge::RVNGPropertyList& propList) override { printf("polygon\n"); }
  void drawPolyline(const librevenge::RVNGPropertyList& propList) override { printf("polyline\n"); }
  void drawConnector(const librevenge::RVNGPropertyList& propList) override { printf("rect\n"); }
  void startTextObject(const librevenge::RVNGPropertyList& propList) override {}
  void endTextObject() override {}
  void startTableObject(const librevenge::RVNGPropertyList& propList) override {}
  void openTableRow(const librevenge::RVNGPropertyList& propList) override {}
  void closeTableRow() override {}
  void openTableCell(const librevenge::RVNGPropertyList& propList) override {}
  void closeTableCell() override {}
  void insertCoveredTableCell(const librevenge::RVNGPropertyList& propList) override {}
  void endTableObject() override {}
  void insertTab() override {}
  void insertSpace() override {}
  void insertText(const librevenge::RVNGString& text) override {}
  void insertLineBreak() override {}
  void insertField(const librevenge::RVNGPropertyList& propList) override {}
  void openOrderedListLevel(const librevenge::RVNGPropertyList& propList) override {}
  void openUnorderedListLevel(const librevenge::RVNGPropertyList& propList) override {}
  void closeOrderedListLevel() override {}
  void closeUnorderedListLevel() override {}
  void openListElement(const librevenge::RVNGPropertyList& propList) override {}
  void closeListElement() override {}
  void defineParagraphStyle(const librevenge::RVNGPropertyList& propList) override {}
  void openParagraph(const librevenge::RVNGPropertyList& propList) override {}
  void closeParagraph() override {}
  void defineCharacterStyle(const librevenge::RVNGPropertyList& propList) override {}
  void openSpan(const librevenge::RVNGPropertyList& propList) override {}
  void closeSpan() override {}
  void openLink(const librevenge::RVNGPropertyList& propList) override {}
  void closeLink() override {}
  void endPage() override {}
  void endMasterPage() override {}
};

std::unique_ptr<Polygon2d> import_cdr(CurveDiscretizer discretizer, const std::string& filename,
                                      const Location& loc)
{
  Polygon2d result;
  printf("import cdr\n");

  librevenge::RVNGFileStream input(filename.c_str());
  CDRReader collector;

  if (!libcdr::CDRDocument::parse(&input, &collector)) {
    std::cerr << "Parse failed\n";
    return nullptr;
  }
  printf("suicceeded\n");

  return std::make_unique<Polygon2d>(result);
}
