#include "core/PatchNode.h"

#include <Python.h>

#include <cstdint>
#include <sstream>

#include "geometry/PolySet.h"
#include "geometry/patch.h"
#include "python/pyconversion.h"
#include "python/pyopenscad.h"
#include "utils/printutils.h"

// ---------------------------------------------------------------------
// Refcounting: proj_func/displacement_func werden hier (nicht am Aufrufort
// in py_primitives.cc) ge-incref-t/decref-t, damit JEDER Weg, wie ein
// PatchNode entsteht oder vergeht - auch ueber den generischen
// Copy-Constructor-Klon-Pfad in node_clone.cc (NodeCloneFunc(PatchNode),
// std::make_shared<PatchNode>(*node)) - automatisch korrekt behandelt wird.
// (Siehe die SheetNode-Debugging-Historie: fehlendes INCREF fuehrte zu
// Use-after-free, fehlendes INCREF beim Klonen zu Double-Free.)
// ---------------------------------------------------------------------

PatchNode::PatchNode(const PatchNode& other) : LeafNode(other)
{
  outer = other.outer;
  holes = other.holes;
  outer_normal = other.outer_normal;
  holes_normal = other.holes_normal;
  use_tangents = other.use_tangents;
  grid_spacing_uv = other.grid_spacing_uv;
  proj_func_hash = other.proj_func_hash;
  displacement_func_hash = other.displacement_func_hash;

  proj_func = other.proj_func;
  displacement_func = other.displacement_func;

  PyGILState_STATE gstate = PyGILState_Ensure();
  Py_XINCREF(static_cast<PyObject *>(proj_func));
  Py_XINCREF(static_cast<PyObject *>(displacement_func));
  PyGILState_Release(gstate);
}

PatchNode::~PatchNode()
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  Py_XDECREF(static_cast<PyObject *>(proj_func));
  Py_XDECREF(static_cast<PyObject *>(displacement_func));
  PyGILState_Release(gstate);
}

namespace {

uint64_t fnv1a_mix(uint64_t h, const char *buf, size_t len)
{
  for (size_t i = 0; i < len; i++) {
    h ^= static_cast<unsigned char>(buf[i]);
    h *= 1099511628211ULL;
  }
  return h;
}

uint64_t hashPoints(uint64_t h, const std::vector<Vector3d>& pts)
{
  for (const auto& p : pts) {
    h = fnv1a_mix(h, reinterpret_cast<const char *>(p.data()), sizeof(double) * 3);
  }
  return h;
}

// -------------------- GIL-sichere Python-Callback-Aufrufe --------------------
// Konvention: proj(p) / displacement(p) bekommen jeweils EIN Argument -
// eine 3-elementige Python-Liste [x,y,z] - wie in den zuvor besprochenen
// Beispielen (proj_xy(p), lambda p: bump(p[0], p[1], p[2])).

PyObject *pointToPyList(const Vector3d& p)
{
  PyObject *px = PyFloat_FromDouble(p.x());
  PyObject *py = PyFloat_FromDouble(p.y());
  PyObject *pz = PyFloat_FromDouble(p.z());
  PyObject *list = PyList_New(3);
  PyList_SET_ITEM(list, 0, px);  // steals references
  PyList_SET_ITEM(list, 1, py);
  PyList_SET_ITEM(list, 2, pz);
  return list;
}

bool callProjFunc(PyObject *func, const Vector3d& p, Vector2d& out)
{
  PyGILState_STATE gstate = PyGILState_Ensure();

  PyObject *point = pointToPyList(p);
  PyObject *args = PyTuple_Pack(1, point);
  Py_DECREF(point);

  PyObject *result = PyObject_CallObject(func, args);
  Py_DECREF(args);

  bool ok = false;
  if (result != nullptr) {
    double x = 0, y = 0;
    if (python_vectorval(result, 2, 2, &x, &y, nullptr, nullptr, nullptr) == 0) {
      out = Vector2d(x, y);
      ok = true;
    } else {
      LOG(message_group::Error, "patch(): proj() muss einen 2er-Vektor [u,v] zurueckgeben.");
    }
    Py_DECREF(result);
  } else {
    std::string errorstr;
    python_catch_error(errorstr);
    LOG(message_group::Error, "patch(): Fehler beim Aufruf von proj(): %1$s", errorstr.c_str());
  }

  PyGILState_Release(gstate);
  return ok;
}

bool callDisplacementFunc(PyObject *func, const Vector3d& p, double& out)
{
  PyGILState_STATE gstate = PyGILState_Ensure();

  PyObject *point = pointToPyList(p);
  PyObject *args = PyTuple_Pack(1, point);
  Py_DECREF(point);

  PyObject *result = PyObject_CallObject(func, args);
  Py_DECREF(args);

  bool ok = false;
  if (result != nullptr) {
    if (PyFloat_Check(result) || PyLong_Check(result)) {
      out = PyFloat_AsDouble(result);
      ok = true;
    } else {
      LOG(message_group::Error, "patch(): displacement() muss eine Zahl zurueckgeben.");
    }
    Py_DECREF(result);
  } else {
    std::string errorstr;
    python_catch_error(errorstr);
    LOG(message_group::Error, "patch(): Fehler beim Aufruf von displacement(): %1$s", errorstr.c_str());
  }

  PyGILState_Release(gstate);
  return ok;
}

}  // namespace

std::string PatchNode::toString() const
{
  uint64_t h = 1469598103934665603ULL;
  h = hashPoints(h, outer);
  for (const auto& hole : holes) h = hashPoints(h, hole);
  // Tangenten mit in den Hash aufnehmen: 'outer'/'holes' (reine 3D-
  // Positionen) koennen fuer zwei PatchNodes identisch sein, waehrend
  // outer_normal/holes_normal (und damit das Ergebnis) sich unterscheiden
  // - siehe Kommentar bei 'use_tangents' in PatchNode.h.
  h = hashPoints(h, outer_normal);
  for (const auto& hn : holes_normal) h = hashPoints(h, hn);

  std::ostringstream stream;
  stream << this->name() << "(outer_n=" << outer.size() << ", holes_n=" << holes.size()
         << ", points_hash=" << std::hex << h << std::dec << ", proj=" << proj_func_hash
         << ", displacement=" << displacement_func_hash << ", grid_spacing_uv=" << grid_spacing_uv
         << ", use_tangents=" << (use_tangents ? 1 : 0) << ")";
  return stream.str();
}

std::unique_ptr<const Geometry> PatchNode::createGeometry() const
{
  // -------- DEBUG --------
  // toString() ist (laut eigenem Kommentar am Ende dieser Datei) genau
  // die Grundlage fuer den Geometrie-Cache-Schluessel. Wenn zwei
  // verschiedene patch()-Aufrufe im selben Skript hier denselben String
  // ausgeben (z.B. weil proj_func_hash fuer zwei inhaltlich verschiedene,
  // aber gleichnamige "proj"-Funktionen kollidiert), erklaert das eine
  // Cache-Verwechslung zwischen den beiden Aufrufen vollstaendig.

  auto *proj_py = static_cast<PyObject *>(proj_func);
  auto *disp_py = static_cast<PyObject *>(displacement_func);

  bool failed = false;

  // Kein proj() angegeben (proj_py == nullptr, z.B. patch() ohne proj-
  // Argument aufgerufen) -> ein leeres std::function weitergeben, damit
  // patch() selbst automatisch eine Projektion waehlt (siehe
  // computeAutoProj() in geometry/patch.cc). NICHT hier schon irgendeine
  // Fallback-Funktion aufrufen - ein leeres std::function ist das
  // vereinbarte Signal fuer "bitte automatisch bestimmen".
  std::function<Vector2d(const Vector3d&)> projFn;
  if (proj_py != nullptr) {
    projFn = [proj_py, &failed](const Vector3d& p) -> Vector2d {
      Vector2d out(0, 0);
      if (!callProjFunc(proj_py, p, out)) failed = true;
      return out;
    };
  }
  auto dispFn = [&](const Vector3d& p) -> double {
    double out = 0.0;
    if (disp_py != nullptr) {
      if (!callDisplacementFunc(disp_py, p, out)) failed = true;
    }
    return out;
  };

  auto result = patch(outer, holes, projFn, grid_spacing_uv, dispFn, outer_normal, holes_normal);

  if (failed || result == nullptr) {
    LOG(message_group::Error, "patch(): Geometrieerzeugung fehlgeschlagen.");
    return std::make_unique<PolySet>(3);  // leeres, aber gueltiges Ergebnis statt Crash
  }
  return result;
}
