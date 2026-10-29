/*Copyright(c) 2018, slugdev
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met :
1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.
3. All advertising materials mentioning features or use of this software
must display the following acknowledgement :
This product includes software developed by slugdev.
4. Neither the name of the slugdev nor the
names of its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY SLUGDEV ''AS IS'' AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.IN NO EVENT SHALL SLUGDEV BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <sstream>
#include <math.h>
#include <cmath>
#include <cstdio>
#include "src/geometry/GeometryUtils.h"

/**
 * Format a number as an ISO 10303-21 REAL: always with a decimal point and an
 * upper case exponent ("10." instead of "10", "1.E-05" instead of "1e-05"),
 * locale independent and with 15 significant digits so coordinates survive
 * the round trip. Plain ostream output violates the first and truncates to 6.
 */
inline std::string stepReal(double v)
{
  if (!std::isfinite(v) || v == 0.0) return "0.";
  char buf[40];
  std::snprintf(buf, sizeof(buf), "%.15G", v);
  std::string s(buf);
  std::replace(s.begin(), s.end(), ',', '.');  // in case a locale sneaked in
  const auto e = s.find('E');
  std::string mantissa = s.substr(0, e);
  const std::string exponent = e == std::string::npos ? "" : s.substr(e);
  if (mantissa.find('.') == std::string::npos) mantissa += '.';
  return mantissa + exponent;
}

class StepKernel
{
public:
  class Entity
  {
  public:
    Entity(std::vector<Entity *>& ent_list)
    {
      ent_list.push_back(this);
      id = int(ent_list.size());
    }
    virtual ~Entity() {}

    std::vector<std::string> tokenize(const std::string& str, const std::string& delimiters = ",")
    {
      std::vector<std::string> tokens;
      // Skip delimiters at beginning.
      std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);

      // Find first non-delimiter.
      std::string::size_type pos = str.find_first_of(delimiters, lastPos);

      while (std::string::npos != pos || std::string::npos != lastPos) {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));

        // Skip delimiters.
        lastPos = str.find_first_not_of(delimiters, pos);

        // Find next non-delimiter.
        pos = str.find_first_of(delimiters, lastPos);
      }
      return tokens;
    }

    virtual void serialize(std::ostream& stream_in) = 0;
    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args) = 0;
    int id;
    std::string label;
  };

  class Direction : public Entity
  {
  public:
    Direction(std::vector<Entity *>& ent_list) : Entity(ent_list)
    {
      pt[0] = 0;
      pt[1] = 0;
      pt[2] = 0;
    }
    Direction(std::vector<Entity *>& ent_list, Vector3d pt_in) : Entity(ent_list) { pt = pt_in; }

    virtual ~Direction() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = DIRECTION('" << label << "', (" << stepReal(pt[0]) << ", "
                << stepReal(pt[1]) << ", " << stepReal(pt[2]) << "));\n";
    }

    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args)
    {
      auto st = args.find_first_of('(');
      auto en = args.find_last_of(')');
      auto arg_str = args.substr(st + 1, en - st - 1);
      std::replace(arg_str.begin(), arg_str.end(), ',', ' ');
      std::stringstream ss(arg_str);
      ss >> pt[0] >> pt[1] >> pt[2];
    }
    Vector3d pt;
  };

  class Point : public Entity
  {
  public:
    Point(std::vector<Entity *>& ent_list) : Entity(ent_list)
    {
      pt[0] = 0;
      pt[1] = 0;
      pt[2] = 0;
    }
    Point(std::vector<Entity *>& ent_list, Vector3d pt_in) : Entity(ent_list) { pt = pt_in; }

    virtual ~Point() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = CARTESIAN_POINT('" << label << "', (" << stepReal(pt[0]) << ","
                << stepReal(pt[1]) << "," << stepReal(pt[2]) << "));\n";
    }
    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args)
    {
      auto st = args.find_first_of('(');
      auto en = args.find_last_of(')');
      auto arg_str = args.substr(st + 1, en - st - 1);
      std::replace(arg_str.begin(), arg_str.end(), ',', ' ');
      std::stringstream ss(arg_str);
      ss >> pt[0] >> pt[1] >> pt[2];
    }
    Vector3d pt;
  };

  class Axis2Placement : public Entity
  {
  public:
    Axis2Placement(std::vector<Entity *>& ent_list) : Entity(ent_list)
    {
      dir1 = 0;
      dir2 = 0;
      point = 0;
    }

    Axis2Placement(std::vector<Entity *>& ent_list, Direction *dir1_in, Direction *dir2_in,
                   Point *point_in)
      : Entity(ent_list)
    {
      dir1 = dir1_in;
      dir2 = dir2_in;
      point = point_in;
    }

    virtual ~Axis2Placement() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = AXIS2_PLACEMENT_3D('" << label << "',#" << point->id << ",#"
                << dir1->id << ",#" << dir2->id << ");\n";
    }

    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args)
    {
      auto st = args.find_first_of(',');
      auto arg_str = args.substr(st + 1);
      std::replace(arg_str.begin(), arg_str.end(), ',', ' ');
      std::replace(arg_str.begin(), arg_str.end(), '#', ' ');
      std::stringstream ss(arg_str);
      int d1_id, d2_id, p_id;
      ss >> p_id >> d1_id >> d2_id;

      dir1 = dynamic_cast<Direction *>(ent_map[d1_id]);
      dir2 = dynamic_cast<Direction *>(ent_map[d2_id]);
      point = dynamic_cast<Point *>(ent_map[p_id]);
    }

    Direction *dir1;
    Direction *dir2;
    Point *point;
  };

  class SurfaceType : public Entity
  {
  public:
    SurfaceType(std::vector<Entity *>& ent_list) : Entity(ent_list) {}
  };

  class Plane : public SurfaceType
  {
  public:
    Plane(std::vector<Entity *>& ent_list) : SurfaceType(ent_list) { axis = 0; }

    Plane(std::vector<Entity *>& ent_list, Axis2Placement *axis_in) : SurfaceType(ent_list)
    {
      axis = axis_in;
    }
    virtual ~Plane() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = PLANE('" << label << "',#" << axis->id << ");\n";
    }
    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args)
    {
      auto st = args.find_first_of(',');
      auto arg_str = args.substr(st + 1);
      std::replace(arg_str.begin(), arg_str.end(), ',', ' ');
      std::replace(arg_str.begin(), arg_str.end(), '#', ' ');
      std::stringstream ss(arg_str);
      int p_id;
      ss >> p_id;

      axis = dynamic_cast<Axis2Placement *>(ent_map[p_id]);
    }

    Axis2Placement *axis;
  };

  class CylindricalSurface : public SurfaceType
  {
  public:
    CylindricalSurface(std::vector<Entity *>& ent_list) : SurfaceType(ent_list)
    {
      axis = 0;
      r = 0;
    }

    CylindricalSurface(std::vector<Entity *>& ent_list, std::string name, Axis2Placement *axis_in,
                       double r_in)
      : SurfaceType(ent_list)
    {
      axis = axis_in;
      r = r_in;
    }
    virtual ~CylindricalSurface() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = CYLINDRICAL_SURFACE('" << label << "',#" << axis->id << ","
                << stepReal(r) << ");\n";
    }
    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args)
    {
      auto st = args.find_first_of(',');
      auto arg_str = args.substr(st + 1);
      std::replace(arg_str.begin(), arg_str.end(), ',', ' ');
      std::replace(arg_str.begin(), arg_str.end(), '#', ' ');
      std::stringstream ss(arg_str);
      int p_id;
      ss >> p_id >> r;
      axis = dynamic_cast<Axis2Placement *>(ent_map[p_id]);
    }
    std::string name;
    double r;
    Axis2Placement *axis;
  };

  class RoundType : public Entity
  {
  public:
    RoundType(std::vector<Entity *>& ent_list) : Entity(ent_list) {}
  };

  class Circle : public RoundType
  {
  public:
    Circle(std::vector<Entity *>& ent_list) : RoundType(ent_list)
    {
      axis = 0;
      r = 0;
    }

    Circle(std::vector<Entity *>& ent_list, std::string name, Axis2Placement *axis_in, double r_in)
      : RoundType(ent_list)
    {
      axis = axis_in;
      r = r_in;
    }
    virtual ~Circle() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = CIRCLE('" << label << "',#" << axis->id << ");\n";
      stream_in << "#" << id << " = CIRCLE('" << label << "',#" << axis->id << "," << stepReal(r)
                << ");\n";
    }
    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args)
    {
      auto st = args.find_first_of(',');
      auto arg_str = args.substr(st + 1);
      std::replace(arg_str.begin(), arg_str.end(), ',', ' ');
      std::replace(arg_str.begin(), arg_str.end(), '#', ' ');
      std::stringstream ss(arg_str);
      int p_id;
      ss >> p_id >> r;
      axis = dynamic_cast<Axis2Placement *>(ent_map[p_id]);
    }
    std::string name;
    double r;
    Axis2Placement *axis;
  };

  class OrientedEdge;

  class EdgeLoop : public Entity
  {
  public:
    EdgeLoop(std::vector<Entity *>& ent_list) : Entity(ent_list) {}
    EdgeLoop(std::vector<Entity *>& ent_list, std::vector<OrientedEdge *>& edges_in) : Entity(ent_list)
    {
      faces = edges_in;
    }
    virtual ~EdgeLoop() {}

    virtual void serialize(std::ostream& stream_in)
    {
      // #17 = ADVANCED_FACE('', (#18), #32, .T.);
      stream_in << "#" << id << " = EDGE_LOOP('" << label << "', (";
      for (size_t i = 0; i < faces.size(); i++) {
        stream_in << "#" << faces[i]->id;
        if (i != faces.size() - 1) stream_in << ",";
      }
      stream_in << "));\n";
    }
    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args)
    {
      auto st = args.find_first_of('(');
      auto en = args.find_last_of(')');
      auto arg_str = args.substr(st + 1, en - st - 1);
      std::replace(arg_str.begin(), arg_str.end(), '#', ' ');
      auto vals = tokenize(arg_str);
      for (auto v : vals) {
        int id = std::atoi(v.c_str());
        faces.push_back(dynamic_cast<OrientedEdge *>(ent_map[id]));
      }
      // axis = dynamic_cast<Axis2Placement*>(ent_map[p_id]);
    }
    std::vector<OrientedEdge *> faces;
  };

  class FaceBound : public Entity
  {
  public:
    FaceBound(std::vector<Entity *>& ent_list) : Entity(ent_list)
    {
      edgeLoop = 0;
      dir = true;
    }
    FaceBound(std::vector<Entity *>& ent_list, EdgeLoop *edge_loop_in, bool dir_in) : Entity(ent_list)
    {
      edgeLoop = edge_loop_in;
      dir = dir_in;
    }
    virtual ~FaceBound() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = FACE_BOUND('" << label << "', #" << edgeLoop->id << ","
                << (dir ? ".T." : ".F.") << ");\n";
    }
    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args)
    {
      auto st = args.find_first_of(',');
      auto arg_str = args.substr(st + 1);
      std::replace(arg_str.begin(), arg_str.end(), ',', ' ');
      std::replace(arg_str.begin(), arg_str.end(), '#', ' ');
      std::stringstream ss(arg_str);
      int p_id;
      std::string tf;
      ss >> p_id >> tf;

      edgeLoop = dynamic_cast<EdgeLoop *>(ent_map[p_id]);
      dir = (tf == ".T.");
    }
    EdgeLoop *edgeLoop;
    bool dir;
  };

  class Face : public Entity
  {
  public:
    Face(std::vector<Entity *>& ent_list) : Entity(ent_list)
    {
      dir = true;
      surface = 0;
    }
    Face(std::vector<Entity *>& ent_list, std::vector<FaceBound *> face_bounds_in,
         SurfaceType *surface_in, bool dir_in)
      : Entity(ent_list)
    {
      faceBounds = face_bounds_in;
      dir = dir_in;
      surface = surface_in;
    }
    virtual ~Face() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = ADVANCED_FACE('" << label << "', (";
      for (size_t i = 0; i < faceBounds.size(); i++) {
        stream_in << "#" << faceBounds[i]->id;
        if (i != faceBounds.size() - 1) stream_in << ",";
      }
      stream_in << "),#" << surface->id << "," << (dir ? ".T." : ".F.") << ");\n";
    }
    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args)
    {
      auto st = args.find_first_of('(');
      auto en = args.find_last_of(')');
      auto arg_str = args.substr(st + 1, en - st - 1);
      std::replace(arg_str.begin(), arg_str.end(), '#', ' ');
      auto vals = tokenize(arg_str);
      for (auto v : vals) {
        int id = std::atoi(v.c_str());
        faceBounds.push_back(dynamic_cast<FaceBound *>(ent_map[id]));
      }
      auto remaining = args.substr(en + 1);
      std::replace(remaining.begin(), remaining.end(), '#', ' ');
      std::replace(remaining.begin(), remaining.end(), ',', ' ');
      std::stringstream ss(remaining);
      int p_id;
      std::string tf;
      ss >> p_id >> tf;

      surface = dynamic_cast<SurfaceType *>(ent_map[p_id]);
      dir = (tf == ".T.");
    }

    std::vector<FaceBound *> faceBounds;
    bool dir;
    SurfaceType *surface;
  };

  class Shell : public Entity
  {
  public:
    Shell(std::vector<Entity *>& ent_list) : Entity(ent_list) { isOpen = true; }
    Shell(std::vector<Entity *>& ent_list, std::vector<Face *>& faces_in) : Entity(ent_list)
    {
      faces = faces_in;
      isOpen = true;
    }
    virtual ~Shell() {}

    virtual void serialize(std::ostream& stream_in)
    {
      if (isOpen) stream_in << "#" << id << " = OPEN_SHELL('" << label << "',(";
      else stream_in << "#" << id << " = CLOSED_SHELL('" << label << "',(";

      for (size_t i = 0; i < faces.size(); i++) {
        stream_in << "#" << faces[i]->id;
        if (i != faces.size() - 1) stream_in << ",";
      }
      stream_in << "));\n";
    }
    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args)
    {
      auto st = args.find_first_of('(');
      auto en = args.find_last_of(')');
      auto arg_str = args.substr(st + 1, en - st - 1);
      std::replace(arg_str.begin(), arg_str.end(), '#', ' ');
      auto vals = tokenize(arg_str);
      for (auto v : vals) {
        int id = std::atoi(v.c_str());
        faces.push_back(dynamic_cast<Face *>(ent_map[id]));
      }
    }

    std::vector<Face *> faces;
    bool isOpen;
  };

  class ShellModel : public Entity
  {
  public:
    ShellModel(std::vector<Entity *>& ent_list) : Entity(ent_list) {}
    ShellModel(std::vector<Entity *>& ent_list, std::vector<Shell *> shells_in) : Entity(ent_list)
    {
      shells = shells_in;
    }
    virtual ~ShellModel() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = SHELL_BASED_SURFACE_MODEL('" << label << "', (";
      for (size_t i = 0; i < shells.size(); i++) {
        stream_in << "#" << shells[i]->id;
        if (i != shells.size() - 1) stream_in << ",";
      }
      stream_in << "));\n";
    }
    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args)
    {
      auto st = args.find_first_of('(');
      auto en = args.find_last_of(')');
      auto arg_str = args.substr(st + 1, en - st - 1);
      std::replace(arg_str.begin(), arg_str.end(), '#', ' ');
      auto vals = tokenize(arg_str);
      for (auto v : vals) {
        int id = std::atoi(v.c_str());
        shells.push_back(dynamic_cast<Shell *>(ent_map[id]));
      }
    }

    std::vector<Shell *> shells;
  };

  class ManifoldShape : public Entity
  {
  public:
    ManifoldShape(std::vector<Entity *>& ent_list) : Entity(ent_list)
    {
      axis = 0;
      shellModel = 0;
    }
    ManifoldShape(std::vector<Entity *>& ent_list, Axis2Placement *axis_in, ShellModel *shell_model_in)
      : Entity(ent_list)
    {
      axis = axis_in;
      shellModel = shell_model_in;
    }
    virtual ~ManifoldShape() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = MANIFOLD_SURFACE_SHAPE_REPRESENTATION('" << label << "', (#"
                << axis->id << ", #" << shellModel->id << "));\n";
    }
    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args)
    {
      auto st = args.find_first_of('(');
      auto en = args.find_last_of(')');
      auto arg_str = args.substr(st + 1, en - st - 1);
      std::replace(arg_str.begin(), arg_str.end(), ',', ' ');
      std::replace(arg_str.begin(), arg_str.end(), '#', ' ');
      std::stringstream ss(arg_str);
      int p1_id, p2_id;
      ss >> p1_id >> p2_id;

      axis = dynamic_cast<Axis2Placement *>(ent_map[p1_id]);
      shellModel = dynamic_cast<ShellModel *>(ent_map[p2_id]);

      if (!axis && !shellModel) {
        axis = dynamic_cast<Axis2Placement *>(ent_map[p2_id]);
        shellModel = dynamic_cast<ShellModel *>(ent_map[p1_id]);
      }
    }

    Axis2Placement *axis;
    ShellModel *shellModel;
  };

  class ManifoldSolid : public Entity
  {
  public:
    ManifoldSolid(std::vector<Entity *>& ent_list) : Entity(ent_list)
    {
      axis = 0;
      shell = 0;
    }
    ManifoldSolid(std::vector<Entity *>& ent_list, Axis2Placement *axis_in, Shell *shell_in)
      : Entity(ent_list)
    {
      axis = axis_in;
      shell = shell_in;
    }
    virtual ~ManifoldSolid() {}

    virtual void serialize(std::ostream& stream_in)
    {
      if (axis == 0)
        stream_in << "#" << id << " = MANIFOLD_SOLID_BREP('" << label << "',#" << shell->id << ");\n";
      else
        stream_in << "#" << id << " = MANIFOLD_SURFACE_SHAPE_REPRESENTATION('" << label << "', (#"
                  << axis->id << ", #" << shell->id << "));\n";
    }
    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args)
    {
      auto st = args.find_first_of('(');
      auto en = args.find_last_of(')');
      auto arg_str = args.substr(st + 1, en - st - 1);
      //			std::replace(arg_str.begin(), arg_str.end(), ',', ' ');
      std::replace(arg_str.begin(), arg_str.end(), '#', ' ');
      auto toks = tokenize(arg_str);
      int p1_id = -1, p2_id = -1;
      if (toks.size() >= 1) {
        std::stringstream ss(toks[0]);
        ss >> p1_id;
      }
      if (toks.size() >= 2) {
        std::stringstream ss(toks[1]);
        ss >> p2_id;
      }

      axis = dynamic_cast<Axis2Placement *>(ent_map[p1_id]);
      shell = dynamic_cast<Shell *>(ent_map[p2_id]);

      if (!axis) axis = dynamic_cast<Axis2Placement *>(ent_map[p1_id]);
      if (!shell) {
        shell = dynamic_cast<Shell *>(ent_map[p2_id]);
        printf("t %p\n", ent_map[p2_id]);
      }
    }

    Axis2Placement *axis;
    Shell *shell;
  };

  class Vertex : public Entity
  {
  public:
    Vertex(std::vector<Entity *>& ent_list) : Entity(ent_list) { point = 0; }
    Vertex(std::vector<Entity *>& ent_list, Point *point_in) : Entity(ent_list) { point = point_in; }
    virtual ~Vertex() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = VERTEX_POINT('" << label << "', #" << point->id << ");\n";
    }
    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args)
    {
      auto st = args.find_first_of(',');
      auto arg_str = args.substr(st + 1);
      std::replace(arg_str.begin(), arg_str.end(), ',', ' ');
      std::replace(arg_str.begin(), arg_str.end(), '#', ' ');
      std::stringstream ss(arg_str);
      int p1_id;
      ss >> p1_id;

      point = dynamic_cast<Point *>(ent_map[p1_id]);
    }

    Point *point;
  };

  class Line;

  class SurfaceCurve : public RoundType
  {
  public:
    SurfaceCurve(std::vector<Entity *>& ent_list) : RoundType(ent_list) { line = 0; }
    SurfaceCurve(std::vector<Entity *>& ent_list, Line *surface_curve_in) : RoundType(ent_list)
    {
      line = surface_curve_in;
    }
    virtual ~SurfaceCurve() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = SURFACE_CURVE('" << label << "', #" << line->id << ");\n";
    }
    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args)
    {
      auto st = args.find_first_of(',');
      auto arg_str = args.substr(st + 1);
      std::replace(arg_str.begin(), arg_str.end(), ',', ' ');
      std::replace(arg_str.begin(), arg_str.end(), '#', ' ');
      std::stringstream ss(arg_str);
      int p1_id;
      ss >> p1_id;

      line = dynamic_cast<Line *>(ent_map[p1_id]);
    }

    Line *line;
  };

  class EdgeCurve : public Entity
  {
  public:
    EdgeCurve(std::vector<Entity *>& ent_list) : Entity(ent_list)
    {
      vert1 = 0;
      vert2 = 0;
      dir = true;
    }
    EdgeCurve(std::vector<Entity *>& ent_list, Vertex *vert1_in, Vertex *vert2_in, RoundType *round_in,
              bool dir_in)
      : Entity(ent_list)
    {
      vert1 = vert1_in;
      vert2 = vert2_in;
      round = round_in;
      dir = dir_in;
    }
    virtual ~EdgeCurve() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = EDGE_CURVE('', #" << vert1->id << ", #" << vert2->id << ",#"
                << round->id << "," << (dir ? ".T." : ".F.") << ");\n";
    }

    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args)
    {
      auto st = args.find_first_of(',');
      auto arg_str = args.substr(st + 1);
      std::replace(arg_str.begin(), arg_str.end(), ',', ' ');
      std::replace(arg_str.begin(), arg_str.end(), '#', ' ');
      std::stringstream ss(arg_str);
      int p1_id, p2_id, p3_id;
      std::string tf;
      ss >> p1_id >> p2_id >> p3_id >> tf;

      vert1 = dynamic_cast<Vertex *>(ent_map[p1_id]);
      vert2 = dynamic_cast<Vertex *>(ent_map[p2_id]);
      round = dynamic_cast<RoundType *>(ent_map[p3_id]);
      dir = (tf == ".T.");
    }

    Vertex *vert1;
    Vertex *vert2;
    RoundType *round;
    bool dir;
  };

  class OrientedEdge : public Entity
  {
  public:
    OrientedEdge(std::vector<Entity *>& ent_list) : Entity(ent_list)
    {
      edge = 0;
      dir = 0;
    }
    OrientedEdge(std::vector<Entity *>& ent_list, EdgeCurve *edge_curve_in, bool dir_in)
      : Entity(ent_list)
    {
      edge = edge_curve_in;
      dir = dir_in;
    }
    virtual ~OrientedEdge() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = ORIENTED_EDGE('" << label << "',*,*,#" << edge->id << ","
                << (dir ? ".T." : ".F.") << ");\n";
    }

    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args)
    {
      auto st = args.find_first_of(',');
      auto arg_str = args.substr(st + 1);
      std::replace(arg_str.begin(), arg_str.end(), ',', ' ');
      std::replace(arg_str.begin(), arg_str.end(), '#', ' ');
      std::stringstream ss(arg_str);
      int p1_id;
      std::string s1, s2, tf;
      ss >> s1 >> s2 >> p1_id >> tf;

      edge = dynamic_cast<EdgeCurve *>(ent_map[p1_id]);
      dir = (tf == ".T.");
    }

    bool dir;
    EdgeCurve *edge;
  };

  class Vector : public Entity
  {
  public:
    Vector(std::vector<Entity *>& ent_list) : Entity(ent_list)
    {
      dir = 0;
      length = 0;
    }
    Vector(std::vector<Entity *>& ent_list, Direction *dir_in, double len_in) : Entity(ent_list)
    {
      dir = dir_in;
      length = len_in;
    }
    virtual ~Vector() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = VECTOR('" << label << "',#" << dir->id << "," << stepReal(length)
                << ");\n";
    }

    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args)
    {
      auto st = args.find_first_of(',');
      auto arg_str = args.substr(st + 1);
      std::replace(arg_str.begin(), arg_str.end(), ',', ' ');
      std::replace(arg_str.begin(), arg_str.end(), '#', ' ');
      std::stringstream ss(arg_str);
      int p1_id;
      ss >> p1_id >> length;

      dir = dynamic_cast<Direction *>(ent_map[p1_id]);
    }

    double length;
    Direction *dir;
  };
  // ---------------------------------------------------------------------------
  // Product structure (ISO 10303-41, AP214 automotive_design)
  //
  //   APPLICATION_CONTEXT <- APPLICATION_PROTOCOL_DEFINITION
  //        ^         ^
  //        |         +-- PRODUCT_CONTEXT <- PRODUCT <- PRODUCT_RELATED_PRODUCT_CATEGORY
  //        |                                   ^
  //        |                     PRODUCT_DEFINITION_FORMATION
  //        |                                   ^
  //        +-- PRODUCT_DEFINITION_CONTEXT <- PRODUCT_DEFINITION <- PRODUCT_DEFINITION_SHAPE
  // ---------------------------------------------------------------------------

  // STEP string literal: quote and double embedded apostrophes
  static std::string stepString(const std::string& s)
  {
    std::string out = "'";
    for (char c : s) {
      if (c == '\'') out += "''";
      else out += c;
    }
    return out + "'";
  }

  class ApplicationContext : public Entity
  {
  public:
    ApplicationContext(
      std::vector<Entity *>& ent_list,
      const std::string& application_in = "core data for automotive mechanical design processes")
      : Entity(ent_list), application(application_in)
    {
    }
    virtual ~ApplicationContext() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = APPLICATION_CONTEXT(" << stepString(application) << ");\n";
    }
    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args) {}

    std::string application;
  };

  class ApplicationProtocolDefinition : public Entity
  {
  public:
    ApplicationProtocolDefinition(std::vector<Entity *>& ent_list, ApplicationContext *app_ctx_in,
                                  const std::string& schema_in = "automotive_design", int year_in = 2000)
      : Entity(ent_list), app_ctx(app_ctx_in), schema(schema_in), year(year_in)
    {
    }
    virtual ~ApplicationProtocolDefinition() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = APPLICATION_PROTOCOL_DEFINITION('international standard',"
                << stepString(schema) << "," << year << ",#" << app_ctx->id << ");\n";
    }
    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args) {}

    ApplicationContext *app_ctx;
    std::string schema;
    int year;
  };

  class ProductContext : public Entity
  {
  public:
    ProductContext(std::vector<Entity *>& ent_list, ApplicationContext *app_ctx_in,
                   const std::string& discipline_in = "mechanical")
      : Entity(ent_list), app_ctx(app_ctx_in), discipline(discipline_in)
    {
    }
    virtual ~ProductContext() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = PRODUCT_CONTEXT('',#" << app_ctx->id << "," << stepString(discipline)
                << ");\n";
    }
    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args) {}

    ApplicationContext *app_ctx;
    std::string discipline;
  };

  class Product : public Entity
  {
  public:
    Product(std::vector<Entity *>& ent_list, const std::string& product_id_in,
            const std::string& name_in, std::vector<ProductContext *> contexts_in,
            const std::string& description_in = "")
      : Entity(ent_list),
        product_id(product_id_in),
        name(name_in),
        description(description_in),
        contexts(std::move(contexts_in))
    {
    }
    virtual ~Product() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = PRODUCT(" << stepString(product_id) << "," << stepString(name) << ","
                << stepString(description) << ",(";
      for (size_t i = 0; i < contexts.size(); i++) stream_in << (i ? "," : "") << "#" << contexts[i]->id;
      stream_in << "));\n";
    }
    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args) {}

    std::string product_id;
    std::string name;
    std::string description;
    std::vector<ProductContext *> contexts;
  };

  class ProductRelatedProductCategory : public Entity
  {
  public:
    ProductRelatedProductCategory(std::vector<Entity *>& ent_list, std::vector<Product *> products_in,
                                  const std::string& name_in = "part")
      : Entity(ent_list), name(name_in), products(std::move(products_in))
    {
    }
    virtual ~ProductRelatedProductCategory() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = PRODUCT_RELATED_PRODUCT_CATEGORY(" << stepString(name) << ",$,(";
      for (size_t i = 0; i < products.size(); i++) stream_in << (i ? "," : "") << "#" << products[i]->id;
      stream_in << "));\n";
    }
    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args) {}

    std::string name;
    std::vector<Product *> products;
  };

  class ProductDefinitionFormation : public Entity
  {
  public:
    ProductDefinitionFormation(std::vector<Entity *>& ent_list, Product *product_in,
                               const std::string& formation_id_in = "")
      : Entity(ent_list), product(product_in), formation_id(formation_id_in)
    {
    }
    virtual ~ProductDefinitionFormation() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = PRODUCT_DEFINITION_FORMATION(" << stepString(formation_id) << ",'',#"
                << product->id << ");\n";
    }
    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args) {}

    Product *product;
    std::string formation_id;
  };

  class ProductDefinitionContext : public Entity
  {
  public:
    ProductDefinitionContext(std::vector<Entity *>& ent_list, ApplicationContext *app_ctx_in,
                             const std::string& life_cycle_stage_in = "design")
      : Entity(ent_list), app_ctx(app_ctx_in), life_cycle_stage(life_cycle_stage_in)
    {
    }
    virtual ~ProductDefinitionContext() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = PRODUCT_DEFINITION_CONTEXT('part definition',#" << app_ctx->id << ","
                << stepString(life_cycle_stage) << ");\n";
    }
    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args) {}

    ApplicationContext *app_ctx;
    std::string life_cycle_stage;
  };

  class ProductDefinition : public Entity
  {
  public:
    ProductDefinition(std::vector<Entity *>& ent_list, ProductDefinitionFormation *formation_in,
                      ProductDefinitionContext *context_in, const std::string& def_id_in = "design")
      : Entity(ent_list), def_id(def_id_in), formation(formation_in), context(context_in)
    {
    }
    virtual ~ProductDefinition() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = PRODUCT_DEFINITION(" << stepString(def_id) << ",'',#"
                << formation->id << ",#" << context->id << ");\n";
    }

    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args) {}

    std::string def_id;
    ProductDefinitionFormation *formation;
    ProductDefinitionContext *context;
  };
  class ProductDefinitionShape : public Entity
  {
  public:
    ProductDefinitionShape(std::vector<Entity *>& ent_list, ProductDefinition *prod_in)
      : Entity(ent_list)
    {
      prod = prod_in;
    }
    virtual ~ProductDefinitionShape() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = PRODUCT_DEFINITION_SHAPE('', '', #" << prod->id << ");\n";
    }

    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args) {}

    ProductDefinition *prod;
  };

  // ---------------------------------------------------------------------------
  // Units and representation context (ISO 10303-41 / -43)
  //
  // All unit entities are complex (multi-leaf) instances and must be written
  // as "#id = ( A() B() C() );" with the partial entities in alphabetical
  // order (ISO 10303-21, 11.2.5.5).
  // ---------------------------------------------------------------------------

  // Common base of all SI based named units: carries the SI_UNIT attributes.
  class NamedUnit : public Entity
  {
  public:
    // prefix: SI prefix enum without dots ("MILLI", "CENTI", ...) or "" for none ($)
    // si_name: SI unit name enum without dots ("METRE", "RADIAN", ...)
    NamedUnit(std::vector<Entity *>& ent_list, const std::string& prefix_in,
              const std::string& si_name_in)
      : Entity(ent_list), prefix(prefix_in), si_name(si_name_in)
    {
    }
    virtual ~NamedUnit() {}

    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args) {}

    std::string prefix;
    std::string si_name;

  protected:
    std::string siUnitArgs() const
    {
      return "SI_UNIT(" + (prefix.empty() ? std::string("$") : "." + prefix + ".") + ",." + si_name +
             ".)";
    }
  };

  class LengthUnit : public NamedUnit
  {
  public:
    LengthUnit(std::vector<Entity *>& ent_list, const std::string& prefix_in = "MILLI")
      : NamedUnit(ent_list, prefix_in, "METRE")
    {
    }
    virtual ~LengthUnit() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = ( LENGTH_UNIT() NAMED_UNIT(*) " << siUnitArgs() << " );\n";
    }
  };

  class PlaneAngleUnit : public NamedUnit
  {
  public:
    PlaneAngleUnit(std::vector<Entity *>& ent_list) : NamedUnit(ent_list, "", "RADIAN") {}
    virtual ~PlaneAngleUnit() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = ( NAMED_UNIT(*) PLANE_ANGLE_UNIT() " << siUnitArgs() << " );\n";
    }
  };

  class SolidAngleUnit : public NamedUnit
  {
  public:
    SolidAngleUnit(std::vector<Entity *>& ent_list) : NamedUnit(ent_list, "", "STERADIAN") {}
    virtual ~SolidAngleUnit() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = ( NAMED_UNIT(*) " << siUnitArgs() << " SOLID_ANGLE_UNIT() );\n";
    }
  };

  // Complex instance of
  //   GEOMETRIC_REPRESENTATION_CONTEXT + GLOBAL_UNIT_ASSIGNED_CONTEXT + REPRESENTATION_CONTEXT
  // Referenced as context_of_items by every shape representation.
  class GeometricRepresentationContext : public Entity
  {
  public:
    GeometricRepresentationContext(std::vector<Entity *>& ent_list, std::vector<NamedUnit *> units_in,
                                   int dimension_in = 3)
      : Entity(ent_list), units(std::move(units_in)), dimension(dimension_in)
    {
    }
    virtual ~GeometricRepresentationContext() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = ( GEOMETRIC_REPRESENTATION_CONTEXT(" << dimension << ") "
                << "GLOBAL_UNIT_ASSIGNED_CONTEXT((";
      for (size_t i = 0; i < units.size(); i++) stream_in << (i ? "," : "") << "#" << units[i]->id;
      stream_in << ")) REPRESENTATION_CONTEXT('','" << dimension << "D') );\n";
    }

    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args) {}

    std::vector<NamedUnit *> units;
    int dimension;
  };

  class ShapeRepresentation : public Entity
  {
  public:
    ShapeRepresentation(std::vector<Entity *>& ent_list, const char *name_in,
                        std::vector<Entity *> items_in, GeometricRepresentationContext *context_in)
      : Entity(ent_list), items(std::move(items_in)), context(context_in)
    {
      name = strdup(name_in);
    }
    virtual ~ShapeRepresentation() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = SHAPE_REPRESENTATION(" << stepString(name) << ",(";
      for (size_t i = 0; i < items.size(); i++) stream_in << (i ? "," : "") << "#" << items[i]->id;
      stream_in << "),#" << context->id << ");\n";
    }

    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args) {}

    char *name;
    std::vector<Entity *> items;
    GeometricRepresentationContext *context;
  };
  class ShapeDefinition_Representation : public Entity
  {
  public:
    ShapeDefinition_Representation(std::vector<Entity *>& ent_list,
                                   ProductDefinitionShape *prod_shape_in, ShapeRepresentation *repr_in)
      : Entity(ent_list)
    {
      repr = repr_in;
      prod_shape = prod_shape_in;
    }
    virtual ~ShapeDefinition_Representation() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = SHAPE_DEFINITION_REPRESENTATION(#" << prod_shape->id << " ,#"
                << repr->id << ");\n";
    }

    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args) {}
    ProductDefinitionShape *prod_shape;
    ShapeRepresentation *repr;

    double length;
    Direction *dir;
  };

  class AdvancesBrepRepresentation : public Entity
  {
  public:
    AdvancesBrepRepresentation(std::vector<Entity *>& ent_list, const char *name_in,
                               ManifoldSolid *solid_in, GeometricRepresentationContext *context_in)
      : Entity(ent_list)
    {
      name = strdup(name_in);
      solid = solid_in;
      context = context_in;
    }
    virtual ~AdvancesBrepRepresentation() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = ADVANCED_BREP_SHAPE_REPRESENTATION(" << stepString(name) << ",(#"
                << solid->id << "),#" << context->id << ");\n";
    }

    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args) {}
    char *name;
    ManifoldSolid *solid;
    GeometricRepresentationContext *context;
    double length;
  };

  class ShapeRepresentationRelationShip : public Entity
  {
  public:
    ShapeRepresentationRelationShip(std::vector<Entity *>& ent_list, ShapeRepresentation *shape_repr_in,
                                    AdvancesBrepRepresentation *adv_brep_in)
      : Entity(ent_list)
    {
      shape_repr = shape_repr_in;
      adv_brep = adv_brep_in;
    }
    virtual ~ShapeRepresentationRelationShip() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = SHAPE_REPRESENTATION_RELATIONSHIP('', '', #" << shape_repr->id
                << ", #" << adv_brep->id << ");\n";
    }

    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args) {}
    ShapeRepresentation *shape_repr;
    AdvancesBrepRepresentation *adv_brep;
  };

  class Line : public RoundType
  {
  public:
    Line(std::vector<Entity *>& ent_list) : RoundType(ent_list)
    {
      vector = 0;
      point = 0;
    }
    Line(std::vector<Entity *>& ent_list, Point *point_in, Vector *vec_in) : RoundType(ent_list)
    {
      vector = vec_in;
      point = point_in;
    }
    virtual ~Line() {}

    virtual void serialize(std::ostream& stream_in)
    {
      stream_in << "#" << id << " = LINE('" << label << "',#" << point->id << ", #" << vector->id
                << ");\n";
    }

    virtual void parse_args(std::map<int, Entity *>& ent_map, std::string args)
    {
      auto st = args.find_first_of(',');
      auto arg_str = args.substr(st + 1);
      std::replace(arg_str.begin(), arg_str.end(), ',', ' ');
      std::replace(arg_str.begin(), arg_str.end(), '#', ' ');
      std::stringstream ss(arg_str);
      int p1_id, p2_id;
      ss >> p1_id >> p2_id;

      point = dynamic_cast<Point *>(ent_map[p1_id]);
      vector = dynamic_cast<Vector *>(ent_map[p2_id]);
    }

    Point *point;
    Vector *vector;
  };

public:
  StepKernel();
  virtual ~StepKernel();

  StepKernel::EdgeCurve *create_line_edge_curve(StepKernel::Vertex *vert1, StepKernel::Vertex *vert2,
                                                bool dir);
  StepKernel::EdgeCurve *create_arc_edge_curve(StepKernel::Vertex *vert1, StepKernel::Vertex *vert2,
                                               bool dir);

  void build_tri_body(const char *name, std::vector<Vector3d> tris, std::vector<IndexedFace> faces,
                      const std::vector<int>& faceParents, double tol);
  EdgeCurve *get_line_from_map(Vector3d p0, Vector3d p1,
                               std::map<std::tuple<double, double, double, double, double, double>,
                                        StepKernel::EdgeCurve *>& edge_map,
                               StepKernel::Vertex *vert1, StepKernel::Vertex *vert2, bool& edge_dir,
                               int& merge_cnt);
  EdgeCurve *get_arc_from_map(Vector3d p0, Vector3d p1,
                              std::map<std::tuple<double, double, double, double, double, double>,
                                       StepKernel::EdgeCurve *>& edge_map,
                              StepKernel::Vertex *vert1, StepKernel::Vertex *vert2, bool& edge_dir,
                              int& merge_cnt);
  std::string read_line(std::ifstream& stp_file, bool skip_all_space);
  void read_step(std::string file_name);
  std::vector<Entity *> entities;
};
