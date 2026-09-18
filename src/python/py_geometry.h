#pragma once

#include <memory>

class AbstractNode;
class Geometry;
class GeometryEvaluator;

std::shared_ptr<const Geometry> python_evaluate_geometry_checked(GeometryEvaluator& evaluator,
                                                                 const AbstractNode& node,
                                                                 bool allownef);
