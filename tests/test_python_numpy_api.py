#!/usr/bin/env python3
"""Integration test: the PythonSCAD C API accepts NumPy arrays.

Runs a script through the ``pythonscad`` binary that drives a representative
set of list-accepting APIs (translate/scale/rotate/multmatrix, polygon,
polyhedron, linear_extrude) once with plain Python lists and once with NumPy
arrays, exports each result to STL, and asserts the two exports are
byte-for-byte identical. This proves NumPy inputs flow through the same code
paths and produce the same geometry as list inputs.

If NumPy is not importable inside the bundled interpreter the script prints a
line starting with ``SKIP`` and exits 0; the CTest registration turns that
into a skipped result via SKIP_REGULAR_EXPRESSION.
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile


TIMEOUT_SECONDS = 120

# This program runs *inside* the pythonscad interpreter.
PYTHONSCAD_PROGRAM = r'''
import atexit as _atexit
import os as _os
import sys as _sys
import tempfile as _tempfile

try:
    import numpy as np
except ImportError:
    print("SKIP: numpy not available in the pythonscad interpreter")
    _sys.exit(0)

from pythonscad import (
    cube,
    linear_extrude,
    polygon,
    polyhedron,
)

# TemporaryDirectory + atexit so the scratch dir is always removed when the
# interpreter exits, including on failure, rather than leaking every run.
_tmpdir_ctx = _tempfile.TemporaryDirectory(prefix="numpy_api_")
_atexit.register(_tmpdir_ctx.cleanup)
_tmpdir = _tmpdir_ctx.name
_counter = [0]


def _stl(obj):
    """Export obj to a fresh STL file and return its bytes."""
    _counter[0] += 1
    path = _os.path.join(_tmpdir, "obj_%d.stl" % _counter[0])
    obj.export(path)
    with open(path, "rb") as fh:
        data = fh.read()
    _os.remove(path)
    return data


def _assert_same(name, obj_list, obj_numpy):
    a = _stl(obj_list)
    b = _stl(obj_numpy)
    assert len(a) > 0, name + ": empty export"
    assert a == b, name + ": numpy result differs from list result"
    print("OK", name)


# translate
_assert_same(
    "translate",
    cube(10).translate([1.0, 2.0, 3.0]),
    cube(10).translate(np.array([1.0, 2.0, 3.0])),
)

# scale
_assert_same(
    "scale",
    cube(10).scale([1.0, 2.0, 3.0]),
    cube(10).scale(np.array([1.0, 2.0, 3.0])),
)

# rotate
_assert_same(
    "rotate",
    cube(10).rotate([0.0, 0.0, 45.0]),
    cube(10).rotate(np.array([0.0, 0.0, 45.0])),
)

# multmatrix
_M = [
    [1.0, 0.0, 0.0, 5.0],
    [0.0, 1.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0, 1.0],
]
_assert_same(
    "multmatrix",
    cube(10).multmatrix(_M),
    cube(10).multmatrix(np.array(_M, dtype=float)),
)

# polygon (2D coordinates) -> extrude so it can be exported as STL
_poly = [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]]
_assert_same(
    "polygon",
    linear_extrude(polygon(_poly), height=5),
    linear_extrude(polygon(np.array(_poly, dtype=float)), height=5),
)

# polyhedron (3D points + face indices)
_pts = [
    [0.0, 0.0, 0.0],
    [10.0, 0.0, 0.0],
    [10.0, 10.0, 0.0],
    [0.0, 10.0, 0.0],
    [5.0, 5.0, 10.0],
]
_faces = [[0, 1, 2, 3], [0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]]
_assert_same(
    "polyhedron",
    polyhedron(_pts, _faces),
    polyhedron(np.array(_pts, dtype=float), _faces),
)

print("NUMPY_API_OK")

# Provide a valid top-level object so the -o batch export exits cleanly.
from pythonscad import cube as _cube, show as _show
_show(_cube(1))
'''


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: test_python_numpy_api.py <pythonscad-binary>", file=sys.stderr)
        return 2

    pythonscad = sys.argv[1]
    if not os.path.isfile(pythonscad):
        print(f"binary not found: {pythonscad}", file=sys.stderr)
        return 2

    with tempfile.NamedTemporaryFile(
        "w", suffix=".py", delete=False, encoding="utf-8"
    ) as fh:
        fh.write(PYTHONSCAD_PROGRAM)
        script_path = fh.name

    # Pass an output file so the binary runs headless and exits instead of
    # launching the GUI. The script does its own exports; this top-level
    # output is unused (the script shows nothing) but forces batch mode.
    out_stl = script_path + ".out.stl"

    try:
        proc = subprocess.run(
            [pythonscad, "--trust-python", script_path, "-o", out_stl],
            capture_output=True,
            text=True,
            timeout=TIMEOUT_SECONDS,
        )
    finally:
        os.remove(script_path)
        if os.path.exists(out_stl):
            os.remove(out_stl)

    print("===== stdout =====")
    print(proc.stdout)
    print("===== stderr =====", file=sys.stderr)
    print(proc.stderr, file=sys.stderr)

    # pythonscad routes Python print() output to stderr, so search both.
    combined = (proc.stdout or "") + (proc.stderr or "")

    if "SKIP" in combined:
        # Surface the SKIP marker on stdout for SKIP_REGULAR_EXPRESSION.
        print("SKIP: numpy not available")
        return 0

    if proc.returncode != 0:
        print(f"FAIL: pythonscad exited with {proc.returncode}", file=sys.stderr)
        return 1

    if "NUMPY_API_OK" not in combined:
        print("FAIL: numpy API script did not complete", file=sys.stderr)
        return 1

    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
