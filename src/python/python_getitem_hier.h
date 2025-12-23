#pragma once

#include <Python.h>
#include <memory>
#include <string>

// Forward declarations
class AbstractNode;
typedef std::shared_ptr<AbstractNode> NodePtr;

/**
 * Hierarchical property extraction from geometry nodes
 * Supports depth-based search through AbstractNode tree
 *
 * Supported properties:
 * - "matrix" (TransformNode) - Returns 4x4 transformation matrix
 * - "vertices" (PolygonNode, PolyhedronNode) - Returns vertex coordinates
 * - "paths" (PolygonNode) - Returns path connectivity indices
 * - "indices" (PolyhedronNode) - Returns face connectivity indices
 */

/**
 * Extract property from node with hierarchical depth search
 * @param node The AbstractNode to extract from
 * @param property_name Name of property ("matrix", "vertices", "paths", "indices")
 * @param depth Maximum depth to search (default 2)
 * @return PyObject* containing the property value, or NULL on error
 */
PyObject *python_getitem_hier(const NodePtr& node, const std::string& property_name, int depth = 2);

/**
 * Helper function: Extract "matrix" property from TransformNode
 * @param node The AbstractNode (must be TransformNode)
 * @param depth Search depth for hierarchical access
 * @return PyObject* containing 4x4 matrix as nested Python list
 */
PyObject *python_extract_matrix(const NodePtr& node, int depth);

/**
 * Helper function: Extract "vertices" property from geometry nodes
 * @param node The AbstractNode (PolygonNode or PolyhedronNode)
 * @param depth Search depth for hierarchical access
 * @return PyObject* containing vertices with flexible 2D/3D conversion
 */
PyObject *python_extract_vertices(const NodePtr& node, int depth);

/**
 * Helper function: Extract "paths" property from PolygonNode
 * @param node The AbstractNode (must be PolygonNode)
 * @param depth Search depth for hierarchical access
 * @return PyObject* containing path connectivity indices
 */
PyObject *python_extract_paths(const NodePtr& node, int depth);

/**
 * Helper function: Extract "indices" property from PolyhedronNode
 * @param node The AbstractNode (must be PolyhedronNode)
 * @param depth Search depth for hierarchical access
 * @return PyObject* containing face connectivity indices
 */
PyObject *python_extract_indices(const NodePtr& node, int depth);
