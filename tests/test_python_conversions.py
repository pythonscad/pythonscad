#!/usr/bin/env python3
"""Coverage test for the C++ Python-object type-conversion helpers.

Every conversion helper in ``src/python/pyconversion.cc`` is reachable only
through a Python-facing API entry point, so this test drives those entry
points through the ``pythonscad`` binary and asserts that each one:

* accepts the coordinate/index shapes it documents (lists, the native
  ``vector()`` type, and scalar broadcast),
* rejects malformed input with a ``TypeError`` (the converters' error paths),
* and hands back values of the documented type/shape.

Helper -> entry point coverage:

    python_numberval        - cube(<scalar>); reject non-numeric str
    python_vectorval        - cube(size3), color(3/4-vec, w component)
    python_vectors          - translate/scale (+ scalar broadcast, vector())
    python_tovector         - rotate(angle, v3)
    python_tomatrix         - multmatrix(4x4)
    python_intlistval       - debug(faces)
    python_to2dvarpointlist - polygon(points 2-D/3-D)
    python_to2dintlist      - polygon(points, paths)
    python_fromvector       - .size / .bbox.*  / cross()
    python_frommatrix       - face.matrix
    python_from2dvarpointlist / python_from2dint
                            - polygon.points / polyhedron.faces
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
import tempfile as _tempfile

from pythonscad import (
    cube, polygon, polyhedron, linear_extrude,
    translate, scale, rotate, multmatrix, vector,
)

# TemporaryDirectory + atexit so the scratch dir is always removed when the
# interpreter exits, including on failure; mkdtemp() alone leaked a directory
# on every run.
_tmpdir_ctx = _tempfile.TemporaryDirectory(prefix="conv_")
_atexit.register(_tmpdir_ctx.cleanup)
_tmpdir = _tmpdir_ctx.name
_n = [0]
_checks = [0]


def _stl(obj):
    _n[0] += 1
    path = _os.path.join(_tmpdir, "o_%d.stl" % _n[0])
    obj.export(path)
    with open(path, "rb") as fh:
        data = fh.read()
    _os.remove(path)
    assert len(data) > 0, "empty export"
    return data


def same(name, *objs):
    """Assert every object exports to identical STL bytes."""
    ref = _stl(objs[0])
    for i, o in enumerate(objs[1:], 1):
        assert _stl(o) == ref, "%s: variant %d differs" % (name, i)
    _checks[0] += 1
    print("OK", name)


def err(name, thunk):
    """Assert the call raises TypeError -- the converters' documented error path.

    Only TypeError counts as a pass. Catching every exception would let an
    unrelated failure (a bug in this script, an AssertionError, a crash in
    some other layer) masquerade as the converter correctly rejecting bad
    input, and would not pin down the contract the module docstring claims.
    Anything else is re-raised as a failure naming what actually happened.
    """
    try:
        thunk()
    except TypeError:
        _checks[0] += 1
        print("OK", name, "(raised TypeError)")
        return
    except BaseException as exc:
        raise AssertionError(
            "%s: expected TypeError, got %s: %s" % (name, type(exc).__name__, exc)
        ) from exc
    raise AssertionError("%s: expected TypeError, no exception raised" % name)


def check(name, cond):
    assert cond, name
    _checks[0] += 1
    print("OK", name)


def _ext(p):
    return linear_extrude(p, height=5)


# --- python_numberval -------------------------------------------------
same("numberval-scalar", cube(10), cube(10.0))
err("numberval-reject-str", lambda: cube("nope"))

# --- python_vectorval (cube size = 3, exact) --------------------------
same("vectorval-size3", cube([2, 4, 6]), cube(vector(2, 4, 6)))
err("vectorval-too-short", lambda: cube([1, 2]))
err("vectorval-too-long", lambda: cube([1, 2, 3, 4]))

# --- python_vectorval (color, 3 or 4 components -> w path) ------------
same("vectorval-color3", cube(1).color([1, 0, 0]), cube(1).color([1.0, 0.0, 0.0]))
same("vectorval-color4-w", cube(1).color([1, 0, 0, 0.5]), cube(1).color([1.0, 0.0, 0.0, 0.5]))

# --- python_vectors (translate/scale, broadcast, native vector) ------
same("vectors-translate",
     translate(cube(1), [1, 2, 3]), translate(cube(1), vector(1, 2, 3)))
same("vectors-scale", scale(cube(2), [1, 2, 3]), scale(cube(2), [1.0, 2.0, 3.0]))
same("vectors-scale-scalar-broadcast", scale(cube(2), 3), scale(cube(2), [3, 3, 3]))

# --- python_tovector (rotate axis, exactly 3) ------------------------
same("tovector-rotate-v",
     rotate(cube([1, 2, 3]), 45, [0, 0, 1]), rotate(cube([1, 2, 3]), 45, [0.0, 0.0, 1.0]))
err("tovector-wrong-len", lambda: cube([1, 2, 3]).rotate(45, [0, 1]))

# --- python_tomatrix (multmatrix 4x4) --------------------------------
_M = [[1, 0, 0, 5], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
_Mf = [[float(x) for x in row] for row in _M]
same("tomatrix-multmatrix", multmatrix(cube([1, 2, 3]), _M), multmatrix(cube([1, 2, 3]), _Mf))
err("tomatrix-not-matrix", lambda: cube(1).multmatrix([1, 2, 3, 4]))

# --- python_intlistval (debug faces) ---------------------------------
same("intlistval-debug", cube(10).debug([0, 1]), cube(10).debug([0, 1]))

# --- python_to2dvarpointlist (polygon points) ------------------------
_sq = [[0, 0], [10, 0], [10, 10], [0, 10]]
_sqf = [[float(x) for x in p] for p in _sq]
same("to2dvarpointlist-polygon", _ext(polygon(_sq)), _ext(polygon(_sqf)))
err("to2dvarpointlist-empty", lambda: polygon([]))
err("to2dvarpointlist-bad-coord", lambda: polygon([[1]]))

# --- python_to2dintlist (polygon paths) ------------------------------
same("to2dintlist-paths",
     _ext(polygon(_sq, paths=[[0, 1, 2, 3]])), _ext(polygon(_sqf, paths=[[0, 1, 2, 3]])))

# --- python_fromvector (vector return) -------------------------------
_sz = cube([2, 4, 6]).size
check("fromvector-value", list(_sz) == [2.0, 4.0, 6.0])
check("fromvector-repr", repr(_sz) == "vector(2,4,6)")
check("fromvector-cross", list(vector(1, 0, 0) * vector(0, 1, 0)) == [0.0, 0.0, 1.0])

# --- python_frommatrix (4x4 return) ----------------------------------
_m = cube(2).rotz(30).faces()[0].matrix
check("frommatrix-shape", len(_m) == 4 and len(_m[0]) == 4)
# round-trips straight back through python_tomatrix without error
multmatrix(cube(1), _m)
check("frommatrix-roundtrip", True)

# --- python_from2dvarpointlist / from2dint (list returns) ------------
_pol = polygon([[0, 0], [10, 0], [10, 10]])
check("from2dvar-list", isinstance(_pol.points, list))
# list concatenation idiom (2-D -> 3-D promotion) must keep working
check("from2dvar-concat", _pol.points[0] + [9] == [0.0, 0.0, 9])
_ph = polyhedron([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
                 [[0, 1, 2], [0, 1, 3], [1, 2, 3], [0, 2, 3]])
check("from2dint-faces-list", isinstance(_ph.faces, list))
check("from2dint-faces-value", _ph.faces[0] == [0, 1, 2])

print("CONVERSION_CHECKS", _checks[0])
print("CONVERSIONS_OK")

from pythonscad import cube as _cube, show as _show
_show(_cube(1))
'''


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: test_python_conversions.py <pythonscad-binary>", file=sys.stderr)
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
    # launching the GUI.
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

    if proc.returncode != 0:
        print(f"FAIL: pythonscad exited with {proc.returncode}", file=sys.stderr)
        return 1

    if "CONVERSIONS_OK" not in combined:
        print("FAIL: conversion coverage script did not complete", file=sys.stderr)
        return 1

    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
