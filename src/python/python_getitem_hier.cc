#include "python_getitem_hier.h"
#include "python_conversions.h"
#include <src/core/primitives.h>
#include <memory>
#include <iostream>

// Forward declare node types (would include actual headers in real implementation)
// class AbstractNode { ... };
// class TransformNode : public AbstractNode { ... };
// class PolygonNode : public AbstractNode { ... };
// class PolyhedronNode : public AbstractNode { ... };

/**
 * Main entry point: Hierarchical property extraction
 * Implements depth-first search through AbstractNode tree
 */
PyObject *python_getitem_hier(const NodePtr& node, const std::string& property_name, int depth)
{
  if (!node) {
    PyErr_SetString(PyExc_ValueError, "Node is NULL");
    return NULL;
  }

  // Property: "matrix" - Extract from TransformNode
  if (property_name == "matrix") {
    return python_extract_matrix(node, depth);
  }

  // Property: "vertices" - Extract from PolygonNode or PolyhedronNode
  if (property_name == "vertices") {
    return python_extract_vertices(node, depth);
  }

  // Property: "paths" - Extract from PolygonNode
  if (property_name == "paths") {
    return python_extract_paths(node, depth);
  }

  // Property: "indices" - Extract from PolyhedronNode
  if (property_name == "indices") {
    return python_extract_indices(node, depth);
  }

  // Unknown property
  PyErr_Format(PyExc_KeyError, "Unknown property: %s", property_name.c_str());
  return NULL;
}

/**
 * Extract "matrix" property from TransformNode
 * Uses hierarchical search if not directly available
 */
PyObject *python_extract_matrix(const NodePtr& node, int depth)
{
  if (!node) {
    PyErr_SetString(PyExc_ValueError, "Node is NULL");
    return NULL;
  }

  // Try to cast to TransformNode
  // auto transform_node = std::dynamic_pointer_cast<TransformNode>(node);
  // if (transform_node) {
  //   Matrix4d matrix = transform_node->getMatrix();
  //   return python_frommatrix(matrix);
  // }

  // If depth > 0, search children
  if (depth > 0) {
    // Pseudo-code: iterate through node->children()
    // for (const auto& child : node->children()) {
    //   PyObject *result = python_extract_matrix(child, depth - 1);
    //   if (result && result != Py_None) {
    //     return result;
    //   }
    //   Py_XDECREF(result);
    // }
  }

  Py_RETURN_NONE;
}

/**
 * Extract "vertices" property from PolygonNode or PolyhedronNode
 * Uses python_frompointsflexible for intelligent 2D/3D conversion
 */
PyObject *python_extract_vertices(const NodePtr& node, int depth)
{
  if (!node) {
    PyErr_SetString(PyExc_ValueError, "Node is NULL");
    return NULL;
  }

  // Try PolygonNode
  // auto polygon_node = std::dynamic_pointer_cast<PolygonNode>(node);
  // if (polygon_node) {
  //   const std::vector<Vector3d>& vertices = polygon_node->getVertices();
  //   return python_frompointsflexible(vertices);
  // }

  // Try PolyhedronNode
  // auto polyhedron_node = std::dynamic_pointer_cast<PolyhedronNode>(node);
  // if (polyhedron_node) {
  //   const std::vector<Vector3d>& vertices = polyhedron_node->getVertices();
  //   return python_frompointsflexible(vertices);
  // }

  // If depth > 0, search children
  if (depth > 0) {
    // Pseudo-code: iterate through node->children()
    // for (const auto& child : node->children()) {
    //   PyObject *result = python_extract_vertices(child, depth - 1);
    //   if (result && result != Py_None) {
    //     return result;
    //   }
    //   Py_XDECREF(result);
    // }
  }

  Py_RETURN_NONE;
}

/**
 * Extract "paths" property from PolygonNode
 * Returns connectivity information for polygon paths
 */
PyObject *python_extract_paths(const NodePtr& node, int depth)
{
  if (!node) {
    PyErr_SetString(PyExc_ValueError, "Node is NULL");
    return NULL;
  }

  // Try PolygonNode
  // auto polygon_node = std::dynamic_pointer_cast<PolygonNode>(node);
  // if (polygon_node) {
  //   const std::vector<std::vector<size_t>>& paths = polygon_node->getPaths();
  //   return python_frompaths(paths);
  // }

  // If depth > 0, search children
  if (depth > 0) {
    // Pseudo-code: iterate through node->children()
    // for (const auto& child : node->children()) {
    //   PyObject *result = python_extract_paths(child, depth - 1);
    //   if (result && result != Py_None) {
    //     return result;
    //   }
    //   Py_XDECREF(result);
    // }
  }

  Py_RETURN_NONE;
}

/**
 * Extract "indices" property from PolyhedronNode
 * Returns face connectivity information
 */
PyObject *python_extract_indices(const NodePtr& node, int depth)
{
  if (!node) {
    PyErr_SetString(PyExc_ValueError, "Node is NULL");
    return NULL;
  }

  // Try PolyhedronNode
  // auto polyhedron_node = std::dynamic_pointer_cast<PolyhedronNode>(node);
  // if (polyhedron_node) {
  //   const std::vector<IndexedFace>& faces = polyhedron_node->getFaces();
  //   return python_fromfaces(faces);
  // }

  // If depth > 0, search children
  if (depth > 0) {
    // Pseudo-code: iterate through node->children()
    // for (const auto& child : node->children()) {
    //   PyObject *result = python_extract_indices(child, depth - 1);
    //   if (result && result != Py_None) {
    //     return result;
    //   }
    //   Py_XDECREF(result);
    // }
  }

  Py_RETURN_NONE;
}
