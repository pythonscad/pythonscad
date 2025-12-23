#!/bin/bash

echo "=== C++ CODE QUALITY CHECK ==="
echo ""

echo "1. Checking for undefined types and forward declarations..."
echo ""

# Check python_conversions.cc
echo "python_conversions.cc:"
grep -E "Matrix4d|Vector3d|IndexedFace|Value|PyObject" src/python/python_conversions.cc | head -5
echo ""

# Check python_getitem_hier.cc
echo "python_getitem_hier.cc:"
grep -E "AbstractNode|TransformNode|PolygonNode|PolyhedronNode" src/python/python_getitem_hier.cc | head -5
echo ""

echo "2. Checking function return types..."
echo ""
grep -E "^(int|PyObject|Value|void)" src/python/python_conversions.cc | sort | uniq -c
echo ""

echo "3. Checking for memory management (Py_INCREF/DECREF/Py_RETURN)..."
echo ""
echo "python_conversions.cc:"
echo "  Py_INCREF: $(grep -c 'Py_INCREF' src/python/python_conversions.cc)"
echo "  Py_DECREF: $(grep -c 'Py_DECREF' src/python/python_conversions.cc)"
echo "  Py_XDECREF: $(grep -c 'Py_XDECREF' src/python/python_conversions.cc)"
echo "  Py_RETURN: $(grep -c 'Py_RETURN' src/python/python_conversions.cc)"
echo ""

echo "python_getitem_hier.cc:"
echo "  Py_XDECREF: $(grep -c 'Py_XDECREF' src/python/python_getitem_hier.cc)"
echo "  Py_RETURN: $(grep -c 'Py_RETURN' src/python/python_getitem_hier.cc)"
echo ""

echo "4. Checking for NULL checks..."
echo ""
echo "python_conversions.cc NULL checks: $(grep -c '!=' src/python/python_conversions.cc | head -1)"
echo "python_getitem_hier.cc NULL checks: $(grep -c '!node' src/python/python_getitem_hier.cc)"
echo ""

echo "5. Checking error handling patterns..."
echo ""
echo "python_conversions.cc:"
grep "PyErr_" src/python/python_conversions.cc | sort | uniq -c
echo ""
echo "python_getitem_hier.cc:"
grep "PyErr_" src/python/python_getitem_hier.cc | sort | uniq -c
echo ""

echo "6. Verifying function declarations match definitions..."
echo ""

# Extract function names from header
echo "Functions declared in python_conversions.h:"
grep "^[a-zA-Z]" src/python/python_conversions.h | grep -v "//" | head -15
echo ""

echo "Functions defined in python_conversions.cc:"
grep "^[a-zA-Z].*{$" src/python/python_conversions.cc | sed 's/ {$//'
echo ""

echo "7. Checking for common C++ issues..."
echo ""

# Check for vector usage
echo "Vector usage in python_conversions.cc: $(grep -c 'std::vector' src/python/python_conversions.cc)"
echo "Vector usage in python_getitem_hier.cc: $(grep -c 'std::vector' src/python/python_getitem_hier.cc)"
echo ""

# Check for shared_ptr usage
echo "shared_ptr usage: $(grep -c 'shared_ptr' src/python/python_getitem_hier.cc)"
echo ""

# Check for casts
echo "dynamic_pointer_cast usage: $(grep -c 'dynamic_pointer_cast' src/python/python_getitem_hier.cc)"
echo ""
