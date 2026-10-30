// Exercise theme changes on existing renderers and geometry, without script reevaluation.
#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <sstream>
#include <string>

#include "Feature.h"
#include "core/CSGNode.h"
#include "core/primitives.h"
#include "geometry/PolySet.h"
#include "geometry/PolySetBuilder.h"
#include "geometry/PolySetUtils.h"
#include "glview/ColorMap.h"
#include "glview/OffscreenView.h"
#include "glview/PolySetRenderer.h"
#include "glview/cgal/CGALRenderer.h"
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
    expectPixel(20, ColorMap::getColor(cs, RenderColor::OPENCSG_FACE_FRONT_COLOR), name + " default");
    expectPixel(60, explicit_color, name + " explicit");
    expectPixel(100, ColorMap::getColor(cs, RenderColor::OPENCSG_FACE_BACK_COLOR), name + " cutout");
  }
}

void checkPolygonRenderer(OffscreenView& view, Renderer& renderer, const std::string& root)
{
  for (const auto *theme : {"tomorrow-night", "nord-dark", "tomorrow-night"}) {
    RenderColorScheme scheme(std::filesystem::path(root) / "color-schemes/render" /
                             (std::string(theme) + ".json"));
    require(scheme.valid(), "Cannot load polygon test theme");
    auto& cs = scheme.colorScheme();
    renderer.setColorScheme(cs);
    draw(view, renderer);
    expectPixel(60, ColorMap::getColor(cs, RenderColor::CGAL_FACE_2D_COLOR), "2D default");
  }
}

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

void checkGeometry()
{
  auto ps = triangles();
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
    checkGeometry();
    checkExports(root);
    OffscreenView view(size, size);
    auto ps = triangles();
    PolySetRenderer manifold_renderer(ps);
    checkRenderer(view, manifold_renderer, "PolySet", root);
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
    CGALRenderer cgal_polygon_renderer(polygon);
    checkPolygonRenderer(view, cgal_polygon_renderer, root);
    CGALRenderer cgal_renderer(ps);
    checkRenderer(view, cgal_renderer, "CGAL PolySet", root);
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
