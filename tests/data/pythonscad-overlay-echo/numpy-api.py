"""Check NumPy inputs and PythonSCAD-only vector helpers."""

import atexit as _atexit
import os as _os
import tempfile as _tempfile

import numpy as np
import openscad as _openscad
import pythonscad as _pythonscad
from pythonscad import (
    HAS_NUMPY,
    Matrix4x4,
    Vector3,
    cube,
    linear_extrude,
    polygon,
    polyhedron,
)


assert not hasattr(_openscad, "Vector3"), "Vector3 must remain PythonSCAD-only"
assert _pythonscad.Vector3 is Vector3
assert HAS_NUMPY, "PythonSCAD vector helpers should use NumPy when it is installed"
assert isinstance(Vector3([1, 2, 3]), np.ndarray)
assert isinstance(Matrix4x4(), np.ndarray)

_tmpdir_ctx = _tempfile.TemporaryDirectory(prefix="numpy_api_")
_atexit.register(_tmpdir_ctx.cleanup)
_tmpdir = _tmpdir_ctx.name
_counter = 0


def _stl(obj):
    """Export obj to a fresh STL file and return its bytes."""
    global _counter
    _counter += 1
    path = _os.path.join(_tmpdir, f"obj_{_counter}.stl")
    obj.export(path)
    with open(path, "rb") as fh:
        data = fh.read()
    _os.remove(path)
    return data


def _assert_same(name, obj_list, obj_numpy):
    list_stl = _stl(obj_list)
    numpy_stl = _stl(obj_numpy)
    assert list_stl, f"{name}: empty export"
    assert list_stl == numpy_stl, f"{name}: NumPy result differs from list result"


_assert_same(
    "translate",
    cube(10).translate([1.0, 2.0, 3.0]),
    cube(10).translate(np.array([1.0, 2.0, 3.0])),
)
_assert_same(
    "Vector3",
    cube(10).translate([1.0, 2.0, 3.0]),
    cube(10).translate(Vector3([1.0, 2.0, 3.0])),
)
_assert_same(
    "scale",
    cube(10).scale([1.0, 2.0, 3.0]),
    cube(10).scale(np.array([1.0, 2.0, 3.0])),
)
_assert_same(
    "rotate",
    cube(10).rotate([0.0, 0.0, 45.0]),
    cube(10).rotate(np.array([0.0, 0.0, 45.0])),
)

_matrix = [
    [1.0, 0.0, 0.0, 5.0],
    [0.0, 1.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0, 1.0],
]
_assert_same(
    "multmatrix",
    cube(10).multmatrix(_matrix),
    cube(10).multmatrix(np.array(_matrix, dtype=float)),
)
_assert_same(
    "Matrix4x4",
    cube(10).multmatrix(_matrix),
    cube(10).multmatrix(Matrix4x4(_matrix)),
)

_polygon = [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]]
_assert_same(
    "polygon",
    linear_extrude(polygon(_polygon), height=5),
    linear_extrude(polygon(np.array(_polygon, dtype=float)), height=5),
)

_points = [
    [0.0, 0.0, 0.0],
    [10.0, 0.0, 0.0],
    [10.0, 10.0, 0.0],
    [0.0, 10.0, 0.0],
    [5.0, 5.0, 10.0],
]
_faces = [[0, 1, 2, 3], [0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]]
_assert_same(
    "polyhedron",
    polyhedron(_points, _faces),
    polyhedron(np.array(_points, dtype=float), _faces),
)

print("OK")
