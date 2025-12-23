# python_getitem_hier() - Hierarchical Property Extraction

## Overview
The `python_getitem_hier()` function implements hierarchical property extraction from geometry nodes using depth-based search through the AbstractNode tree.

## Features

### Supported Properties
1. **"matrix"** - Extract 4x4 transformation matrix from TransformNode
2. **"vertices"** - Extract vertex coordinates from PolygonNode or PolyhedronNode
3. **"paths"** - Extract path connectivity from PolygonNode
4. **"indices"** - Extract face indices from PolyhedronNode

### Intelligent Conversion
- Uses `python_frompointsflexible()` for flexible 2D/3D vertex conversion
- When z==0: returns [x, y] (2D points)
- When z!=0: returns [x, y, z] (3D points)

### Hierarchical Search
- Depth-based search through AbstractNode children
- Default depth: 2 levels
- Returns first match found
- Fallback to parent/sibling nodes

## Implementation Architecture

### Function Hierarchy
```
python_getitem_hier()
├── python_extract_matrix()     [handles "matrix" property]
├── python_extract_vertices()   [handles "vertices" property]
├── python_extract_paths()      [handles "paths" property]
└── python_extract_indices()    [handles "indices" property]
```

### Hierarchical Search Pattern
Each helper function follows this pattern:
```cpp
1. Try direct_cast to specific node type
2. If successful, extract property and return
3. If depth > 0:
   - Iterate through node->children()
   - Recursively search each child with depth-1
   - Return first non-None result
4. If not found, return None
```

## Usage Example

### Python Code
```python
# Get transformation matrix from geometry
matrix = geometry_node['matrix']  # Returns [[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]]

# Get vertices from polygon
vertices = polygon_node['vertices']  # Returns [[x1,y1], [x2,y2], [x3,y3,z3]]

# Get path connectivity
paths = polygon_node['paths']  # Returns [[0,1,2], [2,3,0]]

# Get face indices from polyhedron
indices = polyhedron_node['indices']  # Returns [[0,1,2], [1,2,3], ...]
```

### C++ Integration
```cpp
// In PyModule init code:
NodePtr geometry = create_geometry();

// Extract property with hierarchical search (depth=2)
PyObject *matrix = python_getitem_hier(geometry, "matrix", 2);
PyObject *vertices = python_getitem_hier(geometry, "vertices", 2);
```

## Files Created/Modified

### New Files
- `src/python/python_getitem_hier.h` - Header with function declarations
- `src/python/python_getitem_hier.cc` - Implementation with hierarchical search

### Dependencies
- `python_conversions.h/cc` - Uses conversion helpers
  - `python_frommatrix()`
  - `python_frompointsflexible()`
  - `python_frompaths()`
  - `python_fromfaces()`

### Integration Points
1. Update `src/python/pyopenscad.cc` - Integrate `python__getitem__()` with new hierarchy support
2. Update CMakeLists.txt - Add `python_getitem_hier.cc`
3. Update setup.py - Add `python_getitem_hier.cc`

## Technical Details

### Node Type Detection
Uses C++ `dynamic_pointer_cast<>()` for safe type conversion:
```cpp
auto transform_node = std::dynamic_pointer_cast<TransformNode>(node);
if (transform_node) {
  // Handle TransformNode
}
```

### Memory Management
- Python objects use Py_INCREF/Py_DECREF for reference counting
- Smart pointers (std::shared_ptr) for C++ nodes
- Py_XDECREF for safe NULL checking

### Error Handling
- Returns NULL with PyErr_SetString on error
- Returns Py_RETURN_NONE for not-found (graceful fallback)
- Format errors include property name for debugging

## Pseudo-code Implementation

The actual implementation uses:
1. **dynamic_pointer_cast** - Safe C++ RTTI casting
2. **Getter methods** - node->getMatrix(), node->getVertices(), etc.
3. **Recursive search** - Through node->children() with depth tracking
4. **Python conversion** - Via python_conversions module helpers

## Future Enhancements

1. **Caching** - Memoize property lookups
2. **Custom properties** - Allow user-defined property extractors
3. **Deep search** - Configurable depth limits
4. **Multiple matches** - Return array of matches instead of first
5. **Performance optimization** - Indexed property lookups

## Testing Checklist

- [ ] TransformNode.matrix extraction works
- [ ] PolygonNode.vertices extraction works
- [ ] PolygonNode.paths extraction works
- [ ] PolyhedronNode.vertices extraction works
- [ ] PolyhedronNode.indices extraction works
- [ ] Hierarchical search finds properties at depth 2
- [ ] Flexible 2D/3D conversion works correctly
- [ ] Error handling for invalid properties
- [ ] Error handling for NULL nodes
- [ ] Memory management is correct (no leaks)
