// Exercise theme changes on existing renderers and geometry, without script reevaluation.
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

#include "Feature.h"
#include "core/CSGNode.h"
#include "core/CsgOpNode.h"
#include "core/Tree.h"
#include "core/ProjectionNode.h"
#include "core/TransformNode.h"
#include "glview/RenderSettings.h"
#include "geometry/GeometryEvaluator.h"
#include "geometry/ClipperUtils.h"
#include "core/primitives.h"
#include "geometry/PolySet.h"
#include "geometry/PolySetBuilder.h"
#include "geometry/PolySetUtils.h"
#include "glview/ColorMap.h"
#include "glview/OffscreenView.h"
#include "glview/PolySetRenderer.h"
#ifndef USE_POLYSET_FOR_CGAL
#include "glview/cgal/CGALRenderer.h"
#endif
#include "glview/preview/OpenCSGRenderer.h"
#include "glview/preview/ThrownTogetherRenderer.h"
#include "platform/PlatformUtils.h"
#include "io/export.h"
#ifdef ENABLE_MANIFOLD
#include "geometry/manifold/ManifoldGeometry.h"
#include "geometry/manifold/manifoldutils.h"
#endif
#ifdef ENABLE_CGAL
#include "geometry/cgal/cgalutils.h"
#include "geometry/cgal/CGALNefGeometry.h"
#endif

namespace {
constexpr int size = 120;
const Color4f explicit_color(220, 40, 90, 128);

void require(bool condition, const std::string& message)
{
  if (!condition) throw std::runtime_error(message);
}

// Three separate triangles, sampled at their centers: default, explicit, cut surface.
std::shared_ptr<PolySet> triangles()
{
  auto ps = std::make_shared<PolySet>(3);
  for (int i = 0; i < 3; ++i) {
    double x = -0.8 + i * 0.8;
    ps->vertices.emplace_back(x - 0.3, -0.5, 0);
    ps->vertices.emplace_back(x + 0.3, -0.5, 0);
    ps->vertices.emplace_back(x, 0.6, 0);
    ps->indices.push_back({i * 3, i * 3 + 1, i * 3 + 2});
  }
  // Put an explicit color first to catch palette[0] leaking into uncolored faces.
  ps->colors = {explicit_color, Color4f::defaultCutout()};
  ps->color_indices = {-1, 0, 1};
  ps->setTriangular(true);
  return ps;
}

std::shared_ptr<CSGProducts> products(const std::shared_ptr<const PolySet>& ps)
{
  auto p = std::make_shared<CSGProducts>();
  p->products[0].intersections.emplace_back(
    std::make_shared<CSGLeaf>(ps, Transform3d::Identity(), Color4f(), "theme test", 1));
  return p;
}

void expectPixel(int x, const Color4f& expected, const std::string& context)
{
  std::array<unsigned char, 4> pixel{};
  glReadPixels(x, size / 2, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, pixel.data());
  const auto rgba = expected.toVector4f();
  GLint alpha_bits = 0;
  glGetIntegerv(GL_ALPHA_BITS, &alpha_bits);
  for (int i = 0; i < (alpha_bits ? 4 : 3); ++i) {
    require(std::abs(int(pixel[i]) - int(std::lround(rgba[i] * 255))) <= 1,
            context + ": channel " + std::to_string(i) + " got " + std::to_string(pixel[i]) +
              " expected " + std::to_string(std::lround(rgba[i] * 255)));
  }
}

void draw(OffscreenView& view, Renderer& renderer)
{
  view.fbo->bind();
  glViewport(0, 0, size, size);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-1.2, 1.2, -1.2, 1.2, -10, 10);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glDisable(GL_LIGHTING);
  glDisable(GL_BLEND);
  glDisable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);
  glClearColor(0, 0, 0, 0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
  renderer.prepare(view.edge_shader.get());
  renderer.draw(false, view.edge_shader.get());
  glFinish();
}

// A marker in the existing VBO detects rebuilds even if the driver reuses buffer names.
void checkUnchangedScheme(OffscreenView& view, Renderer& renderer, const ColorScheme& scheme)
{
  GLint buffer = 0;
  glGetIntegerv(GL_VERTEX_ARRAY_BUFFER_BINDING, &buffer);
  require(buffer != 0, "Renderer did not leave a vertex buffer to inspect");
  glBindBuffer(GL_ARRAY_BUFFER, buffer);
  float original = 0;
  glGetBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(original), &original);
  const float marker = 12345.0f;
  glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(marker), &marker);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  // Also exercise an equivalent scheme stored at a different address.
  const auto copy = scheme;
  renderer.setColorScheme(copy);
  renderer.prepare(view.edge_shader.get());
  require(glIsBuffer(buffer), "Unchanged scheme deleted a cached vertex buffer");
  glBindBuffer(GL_ARRAY_BUFFER, buffer);
  float actual = 0;
  glGetBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(actual), &actual);
  require(actual == marker, "Unchanged scheme rebuilt a cached vertex buffer");
  glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(original), &original);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  renderer.setColorScheme(scheme);
}

void checkRenderer(OffscreenView& view, Renderer& renderer, const std::string& name,
                   const std::string& root)
{
  // Switching back also catches one-shot invalidation and accumulated draw states.
  for (const auto *theme : {"tomorrow-night", "nord-dark", "tomorrow-night"}) {
    RenderColorScheme scheme(std::filesystem::path(root) / "color-schemes/render" /
                             (std::string(theme) + ".json"));
    require(scheme.valid(), "Cannot load test theme");
    auto& cs = scheme.colorScheme();
    renderer.setColorScheme(cs);
    draw(view, renderer);
    checkUnchangedScheme(view, renderer, cs);
    expectPixel(20, ColorMap::getColor(cs, RenderColor::OPENCSG_FACE_FRONT_COLOR), name + " default");
    expectPixel(60, explicit_color, name + " explicit");
    expectPixel(100, ColorMap::getColor(cs, RenderColor::OPENCSG_FACE_BACK_COLOR), name + " cutout");
  }
}

void checkPolygonRenderer(OffscreenView& view, Renderer& renderer, const std::string& root)
{
  for (const auto *theme : {"Cornfield", "tomorrow-night", "nord-dark", "Cornfield"}) {
    RenderColorScheme scheme =
      std::string(theme) == "Cornfield"
        ? RenderColorScheme()
        : RenderColorScheme(std::filesystem::path(root) / "color-schemes/render" /
                            (std::string(theme) + ".json"));
    require(scheme.valid(), "Cannot load polygon test theme");
    auto& cs = scheme.colorScheme();
    renderer.setColorScheme(cs);
    draw(view, renderer);
    checkUnchangedScheme(view, renderer, cs);
    expectPixel(60, ColorMap::getColor(cs, RenderColor::OPENCSG_FACE_FRONT_COLOR), "2D default");
  }
}

#ifdef ENABLE_CGAL
// Nef rendering is replaced by PolySetRenderer in this build configuration.
#ifndef USE_POLYSET_FOR_CGAL
void checkNefScheme(OffscreenView& view)
{
  CubeNode cube(nullptr);
  std::fill(std::begin(cube.center), std::end(cube.center), 0);
  auto geometry = cube.createGeometry();
  const auto& mesh = dynamic_cast<const PolySet&>(*geometry);
  std::shared_ptr<const Geometry> nef = CGALUtils::createNefPolyhedronFromPolySet(mesh);
  require(nef && !nef->isEmpty(), "Could not create Nef color test geometry");
  CGALRenderer renderer(nef);
  glPushAttrib(GL_LIGHTING_BIT);
  const GLfloat ambient[] = {1, 1, 1, 1};
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);
  glDisable(GL_LIGHT0);
  glDisable(GL_LIGHT1);
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);
  auto scheme = ColorMap::instance().defaultColorScheme();
  for (const auto color : {Color4f(210, 40, 80), Color4f(30, 180, 220), Color4f(210, 40, 80)}) {
    // Preview, edge, and 2D colors deliberately stay unchanged.
    scheme[RenderColor::CGAL_FACE_FRONT_COLOR] = color;
    scheme[RenderColor::CGAL_FACE_BACK_COLOR] = color;
    renderer.setColorScheme(scheme);
    draw(view, renderer);
    expectPixel(60, color, "Nef face-only scheme change");
    checkUnchangedScheme(view, renderer, scheme);
  }
  glPopAttrib();
}

#endif

void checkCGALRegionVertices()
{
  const Color4f red(255, 0, 0), blue(0, 0, 255);
  const auto rectangle = [](double x, double y, double width, double height, const Color4f& color) {
    Outline2d outline;
    outline.vertices = {{x, y}, {x + width, y}, {x + width, y + height}, {x, y + height}};
    outline.color = color;
    return outline;
  };
  Polygon2d overlap;
  overlap.addOutline(rectangle(0, 0, 10, 10, red));
  overlap.addOutline(rectangle(2, 2, 1, 1, blue));
  auto minimal = CGALUtils::createTriangulatedPolySetFromPolygon2d(overlap, false);
  require(minimal && minimal->indices.size() == 4,
          "CGAL inserted unrelated interior vertices into a color region");

  Polygon2d adjacent;
  adjacent.addOutline(rectangle(0, 0, 1, 1, red));
  Outline2d neighbor;
  neighbor.vertices = {{1, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 0.5}};
  neighbor.color = blue;
  adjacent.addOutline(neighbor);
  auto shared = CGALUtils::createTriangulatedPolySetFromPolygon2d(adjacent, false);
  require(shared != nullptr, "CGAL lost adjacent color regions");
  size_t left_faces = 0;
  bool has_split = false;
  for (size_t i = 0; i < shared->indices.size(); ++i) {
    if (shared->colors.at(shared->color_indices.at(i)) != red) continue;
    ++left_faces;
    const auto& face = shared->indices[i];
    has_split |= std::find(face.begin(), face.end(), 8) != face.end();
  }
  require(left_faces == 3 && has_split, "CGAL lost a shared-boundary split point");

  Polygon2d many_colors;
  for (int i = 0; i < 500; ++i) {
    many_colors.addOutline(rectangle(i * 3, 0, 1, 1, Color4f(float(i + 1) / 501, 0.25f, 0.5f, 1.0f)));
  }
  auto many = CGALUtils::createTriangulatedPolySetFromPolygon2d(many_colors, false);
  require(
    many && many->vertices.size() == 2000 && many->indices.size() == 1000 && many->colors.size() == 500,
    "CGAL did not preserve 500 independent color regions");
}
#endif

#ifdef ENABLE_MANIFOLD
std::shared_ptr<const Polygon2d> projectedCutout(bool recessed)
{
  auto outer = std::make_shared<CubeNode>(std::make_shared<ModuleInstantiation>("cube"));
  auto inner = std::make_shared<CubeNode>(std::make_shared<ModuleInstantiation>("cube"));
  for (int i = 0; i < 3; ++i) {
    outer->dim[i] = 2;
    outer->center[i] = inner->center[i] = 0;
  }
  inner->dim[0] = inner->dim[1] = recessed ? 1 : 4;
  inner->dim[2] = 2;
  auto translate =
    std::make_shared<TransformNode>(std::make_shared<ModuleInstantiation>("translate"), "translate");
  translate->matrix = Transform3d::Identity();
  translate->matrix.translation().z() = 1;
  translate->children = {inner};
  auto difference = std::make_shared<CsgOpNode>(std::make_shared<ModuleInstantiation>("difference"),
                                                OpenSCADOperator::DIFFERENCE);
  difference->children = {outer, translate};
  auto projection =
    std::make_shared<ProjectionNode>(std::make_shared<ModuleInstantiation>("projection"));
  projection->detail_mode = true;
  projection->children = {difference};
  Tree tree(projection);
  GeometryEvaluator evaluator(tree);
  const auto previous_backend = RenderSettings::inst()->backend3D;
  RenderSettings::inst()->backend3D = RenderBackend3D::ManifoldBackend;
  auto polygon =
    std::dynamic_pointer_cast<const Polygon2d>(evaluator.evaluateGeometry(*projection, false));
  RenderSettings::inst()->backend3D = previous_backend;
  require(polygon && !polygon->isEmpty(), "Detail projection lost geometry");
  bool has_cutout = false;
  for (const auto& outline : polygon->outlines()) has_cutout |= outline.color.isDefaultCutout();
  require(has_cutout, "Detail projection lost the upward cut-surface role");
  return polygon;
}

void checkPolygonSamples(OffscreenView& view, const std::shared_ptr<const Polygon2d>& polygon,
                         const Color4f& side, const Color4f& center, const std::string& context)
{
  PolySetRenderer polyset_renderer(polygon);
  std::vector<Renderer *> renderers{&polyset_renderer};
#ifndef USE_POLYSET_FOR_CGAL
  CGALRenderer cgal_renderer(polygon);
  renderers.push_back(&cgal_renderer);
#endif
  for (Renderer *renderer : renderers) {
    RenderColorScheme cornfield;
    renderer->setColorScheme(cornfield.colorScheme());
    draw(view, *renderer);
    expectPixel(20, side, context + " side");
    expectPixel(60, center, context + " center");
  }
}

void checkProjectedCutout(OffscreenView& view)
{
  auto polygon = projectedCutout(true);
  auto cut_only = projectedCutout(false);
  auto offset = std::shared_ptr<const Polygon2d>(
    ClipperUtils::applyOffset(*cut_only, 0.1, Clipper2Lib::JoinType::Miter, 2, 0.01));
  require(offset && !offset->isEmpty(), "Offset lost the projected face");
  for (const auto& outline : offset->outlines())
    require(outline.color.isDefaultCutout(), "Offset lost the projected cut-surface role");

  // A colored ring plus a separate explicit island exercises per-color holes.
  auto ring = std::make_shared<Polygon2d>();
  Outline2d outer, hole, island;
  outer.vertices = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
  outer.color = Color4f::defaultCutout();
  hole.vertices = {{-0.4, -0.4}, {-0.4, 0.4}, {0.4, 0.4}, {0.4, -0.4}};
  hole.positive = false;
  hole.color = Color4f::defaultCutout();
  island.vertices = {{2, -0.5}, {3, -0.5}, {3, 0.5}, {2, 0.5}};
  island.color = explicit_color;
  ring->addOutline(outer);
  ring->addOutline(hole);
  ring->addOutline(island);
#ifdef ENABLE_CGAL
  auto triangulated = CGALUtils::createTriangulatedPolySetFromPolygon2d(*ring, false);
  require(triangulated != nullptr, "CGAL did not tessellate colored regions");
  size_t vertex = 0;
  for (const auto& outline : ring->outlines()) {
    for (const auto& point : outline.vertices) {
      require(triangulated->vertices.at(vertex++).head<2>() == point,
              "CGAL changed source vertex order");
    }
  }
  require(vertex == triangulated->vertices.size(), "CGAL changed source vertex count");
  double area = 0;
  for (const auto& face : triangulated->indices) {
    const auto ab = triangulated->vertices[face[1]] - triangulated->vertices[face[0]];
    const auto ac = triangulated->vertices[face[2]] - triangulated->vertices[face[0]];
    area += ab.cross(ac).norm() / 2;
  }
  require(std::abs(area - 4.36) < 1e-9, "CGAL changed colored-region area or filled a hole");
  bool has_explicit = false;
  for (const auto& color : triangulated->colors) has_explicit |= color == explicit_color;
  require(has_explicit, "CGAL lost explicit RGBA while retaining symbolic colors");
#endif

  RenderColorScheme cornfield;
  const auto material =
    ColorMap::getColor(cornfield.colorScheme(), RenderColor::OPENCSG_FACE_FRONT_COLOR);
  const auto cutout = ColorMap::getColor(cornfield.colorScheme(), RenderColor::OPENCSG_FACE_BACK_COLOR);
  const auto previous_backend = RenderSettings::inst()->backend3D;
  std::vector<RenderBackend3D> triangulators = {RenderBackend3D::ManifoldBackend};
#ifdef ENABLE_CGAL
  triangulators.push_back(RenderBackend3D::CGALBackend);
#endif
  for (auto backend : triangulators) {
    // Select CGAL explicitly even in builds using Manifold's triangulator by default.
    RenderSettings::inst()->backend3D = backend;
    checkPolygonSamples(view, polygon, material, cutout, "Projected cut surface");
    checkPolygonSamples(view, offset, cutout, cutout, "Offset projected cut surface");
    checkPolygonSamples(view, ring, cutout, Color4f(0, 0, 0, 0), "Colored hole");
  }
  RenderSettings::inst()->backend3D = previous_backend;
}

#endif

void checkExports(const std::string& root)
{
  RenderColorScheme scheme(std::filesystem::path(root) / "color-schemes/render/tomorrow-night.json");
  require(scheme.valid(), "Cannot load export theme");
  ExportInfo info{
    .format = FileFormat::OFF,
    .defaultColor = ColorMap::getColor(scheme.colorScheme(), RenderColor::CGAL_FACE_FRONT_COLOR),
    .colorScheme = &scheme.colorScheme()};
  auto ps = triangles();
  // Use a symbolic material palette entry here as produced by mesh conversion.
  ps->colors.push_back(Color4f());
  ps->color_indices[0] = 2;
  std::ostringstream off;
  export_off(ps, off, info);
  require(off.str().find("129 162 190") != std::string::npos,
          "OFF export did not resolve material color");
  require(off.str().find("222 147 95") != std::string::npos,
          "OFF export did not resolve cut-surface color");
  std::ostringstream wrl;
  export_wrl(ps, wrl, info);
  require(wrl.str().find("0.505882 0.635294 0.745098") != std::string::npos,
          "WRL export did not resolve material color");
  require(wrl.str().find("0.870588 0.576471 0.372549") != std::string::npos,
          "WRL export did not resolve cut-surface color");
}

void checkUnpaintedExports(const std::string& root)
{
  RenderColorScheme scheme(std::filesystem::path(root) / "color-schemes/render/tomorrow-night.json");
  ExportInfo info{
    .format = FileFormat::OFF,
    .defaultColor = ColorMap::getColor(scheme.colorScheme(), RenderColor::CGAL_FACE_FRONT_COLOR),
    .colorScheme = &scheme.colorScheme()};
  auto ps = triangles();
  std::ostringstream off;
  export_off(ps, off, info);
  // OFF export always sorts vertices lexicographically, so the first triangle becomes 0 2 1.
  require(off.str().find("3 0 2 1 129 162 190\n") != std::string::npos,
          "OFF negative color index lost the export default");
  std::ostringstream wrl;
  export_wrl(ps, wrl, info);
  require(wrl.str().find("0.505882 0.635294 0.745098, # default colour") != std::string::npos,
          "WRL negative color index lost the export default");
  ps->colors.clear();
  ps->color_indices.clear();
  std::ostringstream unpainted_wrl;
  export_wrl(ps, unpainted_wrl, info);
  require(unpainted_wrl.str().find("diffuseColor 0.505882 0.635294 0.745098") != std::string::npos,
          "WRL empty palette lost the export default");

  SquareNode square(nullptr);
  auto polygon = std::shared_ptr<const Geometry>(square.createGeometry());
  for (int mode : {0, 1}) {
    info.optionsGcode = std::make_shared<const ExportGcodeOptions>(ExportGcodeOptions{
      .feedrate = 4321, .laserpower = 123, .lasermode = mode, .initCode = "", .exitCode = ""});
    std::ostringstream gcode;
    export_gcode(polygon, gcode, info);
    require(
      gcode.str().find("F4321") != std::string::npos && gcode.str().find("S123") != std::string::npos,
      "Unpainted G-code did not use configured power/feed");
  }
}

void checkGeometry()
{
  auto ps = triangles();
  // Inconsistent palettes must not turn concat into an out-of-bounds read.
  auto invalid = std::make_shared<PolyhedronNode>(std::make_shared<ModuleInstantiation>("polyhedron"));
  invalid->points = ps->vertices;
  invalid->faces = ps->indices;
  invalid->colors = ps->colors;
  invalid->color_indices = ps->color_indices;
  invalid->color_indices[0] = std::numeric_limits<int32_t>::max();
  auto valid = std::make_shared<PolyhedronNode>(std::make_shared<ModuleInstantiation>("polyhedron"));
  valid->points = ps->vertices;
  valid->faces = ps->indices;
  valid->colors = ps->colors;
  valid->color_indices = ps->color_indices;
  auto concat = std::make_shared<CsgOpNode>(std::make_shared<ModuleInstantiation>("concat"),
                                            OpenSCADOperator::CONCAT);
  concat->children = {invalid, valid};
  Tree tree(concat);
  GeometryEvaluator evaluator(tree);
  auto joined = std::dynamic_pointer_cast<const PolySet>(evaluator.evaluateGeometry(*concat, false));
  require(joined && joined->indices.size() == 6, "Concat lost geometry with invalid palette index");
  require(joined->color_indices.size() == 6 && joined->color_indices[0] == -1,
          "Concat did not leave an invalid palette index unpainted");
  require(joined->colors.at(joined->color_indices[1]) == explicit_color,
          "Concat lost valid explicit colors alongside an invalid palette index");
  auto tessellated = PolySetUtils::tessellate_faces(*ps);
  require(tessellated->colors == ps->colors, "Tessellation lost symbolic colors");
  auto red = triangles();
  red->setColor(Color4f(255, 0, 0));
  auto blue = triangles();
  blue->setColor(Color4f(0, 0, 255));
  PolySetBuilder append_builder;
  append_builder.appendPolySet(*red);
  append_builder.appendPolySet(*blue);
  auto appended = append_builder.build();
  require(appended->color_indices.size() == appended->indices.size(),
          "Appending meshes duplicated face color indices");
  for (size_t i = 0; i < appended->indices.size(); ++i) {
    const auto expected = i < red->indices.size() ? Color4f(255, 0, 0) : Color4f(0, 0, 255);
    require(appended->colors.at(appended->color_indices[i]) == expected,
            "Appending meshes lost explicit face colors");
  }
  PolySetBuilder builder;
  builder.beginPolygon(3);
  builder.addVertex(Vector3d(0, 0, 0));
  builder.addVertex(Vector3d(1, 0, 0));
  builder.addVertex(Vector3d(0, 1, 0));
  builder.endPolygon(Color4f::defaultCutout());
  auto built = builder.build();
  require(!built->colors.empty() && built->colors[0].isDefaultCutout(),
          "PolySetBuilder lost cut-surface role");
#ifdef ENABLE_MANIFOLD
  auto cube = manifold::Manifold::Cube({1, 1, 1});
  ManifoldGeometry unpainted(cube, {static_cast<uint32_t>(cube.OriginalID())});
  auto mesh = unpainted.toPolySet();
  for (const auto& color : mesh->colors) require(color == Color4f(), "Manifold baked in a theme");
#ifdef ENABLE_CGAL
  auto nef = CGALUtils::createNefPolyhedronFromPolySet(*mesh);
  require(nef && !nef->isEmpty(), "Could not convert mesh to CGAL");
  auto cgal_mesh = CGALUtils::createPolySetFromNefPolyhedron3(*nef->p3);
  require(cgal_mesh != nullptr, "Could not convert CGAL to mesh");
  for (const auto& color : cgal_mesh->colors) require(!color.isValid(), "CGAL baked in a theme");
#endif
  auto converted = ManifoldUtils::createManifoldFromPolySet(*mesh);
  auto cutter = manifold::Manifold::Cube({0.5, 0.5, 2}).Translate({0.25, 0.25, -0.5}).AsOriginal();
  ManifoldGeometry negative(cutter, {static_cast<uint32_t>(cutter.OriginalID())});
  auto difference = (*converted - negative).toPolySet();
  bool has_cutout = false;
  for (const auto& color : difference->colors) {
    require(!color.isValid(), "Boolean baked in a theme");
    has_cutout |= color.isDefaultCutout();
  }
  require(has_cutout, "Boolean lost cut-surface role after mesh round trip");
#endif
}
}  // namespace

int main(int argc, char **argv)
{
  try {
    require(argc == 2 || argc == 3, "Expected source directory and optional --indexing");
    if (argc == 3) {
      require(std::string(argv[2]) == "--indexing", "Unexpected argument");
      Feature::enable_feature("vertex-object-renderers-indexing");
    }
    const std::string root = argv[1];
    PlatformUtils::registerApplicationPath(root);
    for (const auto *name : {"tomorrow", "tomorrow-night", "clearsky", "nord-dark", "nord-light"}) {
      RenderColorScheme scheme(std::filesystem::path(root) / "color-schemes/render" /
                               (std::string(name) + ".json"));
      require(scheme.valid(), "Cannot load theme");
      const auto& cs = scheme.colorScheme();
      require(ColorMap::getColor(cs, RenderColor::OPENCSG_FACE_FRONT_COLOR) ==
                ColorMap::getColor(cs, RenderColor::CGAL_FACE_FRONT_COLOR),
              std::string(name) + ": F5/F6 material defaults differ");
      require(ColorMap::getColor(cs, RenderColor::OPENCSG_FACE_BACK_COLOR) ==
                ColorMap::getColor(cs, RenderColor::CGAL_FACE_BACK_COLOR),
              std::string(name) + ": F5/F6 cutout defaults differ");
    }
    auto& colors = ColorMap::instance();
    require(&colors.colorSchemeOrDefault("__deleted_custom_theme__") == &colors.defaultColorScheme(),
            "A deleted custom theme did not fall back to Cornfield");
    require(&colors.colorSchemeOrDefault("Tomorrow Night") == colors.findColorScheme("Tomorrow Night"),
            "An available custom theme was replaced by the fallback");
#ifdef ENABLE_CGAL
    checkCGALRegionVertices();
#endif
    checkGeometry();
    checkExports(root);
    checkUnpaintedExports(root);
    OffscreenView view(size, size);
#if defined(ENABLE_CGAL) && !defined(USE_POLYSET_FOR_CGAL)
    checkNefScheme(view);
#endif
#ifdef ENABLE_MANIFOLD
    checkProjectedCutout(view);
#endif
    auto ps = triangles();
    PolySetRenderer manifold_renderer(ps);
    checkRenderer(view, manifold_renderer, "PolySet", root);
    auto invalid_palette = triangles();
    invalid_palette->color_indices[0] = std::numeric_limits<int32_t>::max();
    PolySetRenderer invalid_renderer(invalid_palette);
    checkRenderer(view, invalid_renderer, "PolySet invalid palette", root);
    SquareNode square(nullptr);
    square.center = true;
    auto square_geometry = square.createGeometry();
    auto polygon = std::make_shared<Polygon2d>(dynamic_cast<const Polygon2d&>(*square_geometry));
    for (const auto& outline : polygon->outlines())
      require(outline.color == Color4f(), "Square baked in yellow");
#ifdef ENABLE_CGAL
    auto triangulated = CGALUtils::createTriangulatedPolySetFromPolygon2d(*polygon, false);
    for (const auto& color : triangulated->colors)
      require(!color.isValid(), "Polygon tessellation baked in yellow");
#endif
    PolySetRenderer polygon_renderer(polygon);
    checkPolygonRenderer(view, polygon_renderer, root);
#ifndef USE_POLYSET_FOR_CGAL
    CGALRenderer cgal_polygon_renderer(polygon);
    checkPolygonRenderer(view, cgal_polygon_renderer, root);
    CGALRenderer cgal_renderer(ps);
    checkRenderer(view, cgal_renderer, "CGAL PolySet", root);
#endif
#ifdef ENABLE_OPENCSG
    auto p = products(ps);
    OpenCSGRenderer preview(p, nullptr, nullptr);
    checkRenderer(view, preview, "OpenCSG", root);
    ThrownTogetherRenderer thrown(p, nullptr, nullptr);
    checkRenderer(view, thrown, "ThrownTogether", root);
#endif
    std::cout << "Color scheme rendering checks passed\n";
    return 0;
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
