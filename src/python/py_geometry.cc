#include "python/py_geometry.h"

#include <Python.h>

#include "core/FilletDiagnostics.h"
#include "geometry/GeometryEvaluator.h"

std::shared_ptr<const Geometry> python_evaluate_geometry_checked(GeometryEvaluator& evaluator,
                                                                 const AbstractNode& node, bool allownef)
{
  auto geometry = evaluator.evaluateGeometry(node, allownef);
  if (auto error = FilletDiagnostics::takeError()) {
    PyErr_SetString(PyExc_ValueError, error->c_str());
    return nullptr;
  }
  return geometry;
}
