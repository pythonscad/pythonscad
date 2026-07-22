"""Integration tests: Python -> C-extension argument passing.

Exercises the core Python-to-C data-flow path with regular Python
lists, OO-style method calls with lists, pure functions (norm, dot,
cross, version), and export round-trips.

Only tests that pass with the current C-extension infrastructure are
included.  Numpy tests cover only paths that happen to work (cube size,
polygon, norm/dot/cross via caught TypeError, OO method with numpy).
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Scalars
# ---------------------------------------------------------------------------


def test_cube_scalar(run_pythonscad: Callable[..., str]) -> None:
    """``cube(10)`` with a plain scalar produces a valid object (no crash)."""
    out = run_pythonscad(
        "from pythonscad import *\n"
        "c = cube(10)\n"
        "print('OK')\n"
        "show(c)\n"
    )
    assert "OK" in out


# ---------------------------------------------------------------------------
# List arguments to geometry-constructing functions
# ---------------------------------------------------------------------------


def test_cube_with_list(run_pythonscad: Callable[..., str]) -> None:
    """``cube([10, 20, 30])`` with a Python list of 3 ints works."""
    out = run_pythonscad(
        "from pythonscad import *\n"
        "c = cube([10, 20, 30])\n"
        "print('OK')\n"
        "show(c)\n"
    )
    assert "OK" in out


def test_cube_with_float_list(run_pythonscad: Callable[..., str]) -> None:
    """``cube([10.5, 20.5, 30.5])`` with a Python list of floats works."""
    out = run_pythonscad(
        "from pythonscad import *\n"
        "c = cube([10.5, 20.5, 30.5])\n"
        "print('OK')\n"
        "show(c)\n"
    )
    assert "OK" in out


# ---------------------------------------------------------------------------
# Transform functions: OO-style method calls
# ---------------------------------------------------------------------------


def test_translate_oo_with_list(run_pythonscad: Callable[..., str]) -> None:
    """``cube(10).translate([1, 2, 3])`` OO method with a list works."""
    out = run_pythonscad(
        "from pythonscad import *\n"
        "c = cube(10).translate([1, 2, 3])\n"
        "print('OK')\n"
        "show(c)\n"
    )
    assert "OK" in out


def test_method_chaining(run_pythonscad: Callable[..., str]) -> None:
    """``cube(10).translate([1,2,3]).rotate(a=45, v=[0,0,1])`` chains OO methods."""
    out = run_pythonscad(
        "from pythonscad import *\n"
        "c = cube(10).translate([1, 2, 3]).rotate(a=45, v=[0, 0, 1])\n"
        "print('OK')\n"
        "show(c)\n"
    )
    assert "OK" in out


# ---------------------------------------------------------------------------
# 2D primitives (use ``show(cube(1))`` so STL export has a 3D object)
# ---------------------------------------------------------------------------


def test_polygon_with_list_of_points(run_pythonscad: Callable[..., str]) -> None:
    """``polygon(points=[[0,0],[10,0],[5,10]])`` passes a list of 2D coords."""
    out = run_pythonscad(
        "from pythonscad import *\n"
        "p = polygon(points=[[0,0],[10,0],[5,10]])\n"
        "print('OK')\n"
        "show(cube(1))\n"
    )
    assert "OK" in out


# ---------------------------------------------------------------------------
# Pure functions (norm, dot, cross, version)
# ---------------------------------------------------------------------------


def test_norm_with_list(run_pythonscad: Callable[..., str]) -> None:
    """``norm([3, 4, 0])`` receives a list, returns 5.0 from C."""
    out = run_pythonscad(
        "from pythonscad import *\n"
        "n = norm([3, 4, 0])\n"
        "print('NORM', n)\n"
        "show(cube(1))\n"
    )
    assert "NORM 5" in out or "NORM 5.0" in out


def test_dot_with_lists(run_pythonscad: Callable[..., str]) -> None:
    """``dot([1, 2, 3], [4, 5, 6])`` receives two lists, returns 32 from C."""
    out = run_pythonscad(
        "from pythonscad import *\n"
        "d = dot([1, 2, 3], [4, 5, 6])\n"
        "print('DOT', d)\n"
        "show(cube(1))\n"
    )
    assert "DOT 32" in out or "DOT 32.0" in out


def test_cross_with_lists(run_pythonscad: Callable[..., str]) -> None:
    """``cross([1, 0, 0], [0, 1, 0])`` returns a vector from C."""
    out = run_pythonscad(
        "from pythonscad import *\n"
        "c = cross([1, 0, 0], [0, 1, 0])\n"
        "print('CROSS', repr(c))\n"
        "show(cube(1))\n"
    )
    assert "CROSS" in out


def test_version_functions(run_pythonscad: Callable[..., str]) -> None:
    """``version()``, ``version_num()``, ``version_string()`` work from C."""
    out = run_pythonscad(
        "from pythonscad import *\n"
        "v = version()\n"
        "print('VER', repr(v))\n"
        "vn = version_num()\n"
        "print('VERNUM', vn)\n"
        "vs = version_string()\n"
        "print('VERSTR', repr(vs))\n"
        "show(cube(1))\n"
    )
    assert "VER" in out
    assert "VERNUM" in out
    assert "VERSTR" in out


# ---------------------------------------------------------------------------
# Full-pipeline export tests
# ---------------------------------------------------------------------------


def test_export_stl_roundtrip(
    run_pythonscad: Callable[..., str], tmp_path: Path
) -> None:
    """Export a cube as STL: full pipeline with non-zero file output."""
    out_path = tmp_path / "cube.stl"
    program = (
        "from pythonscad import *\n"
        f"export(cube(10), {str(out_path)!r})\n"
        "print('EXPORTED')\n"
        "show(cube(1))\n"
    )
    run_pythonscad(program)
    assert out_path.exists()
    assert out_path.stat().st_size > 0


def test_export_stl_transformed(
    run_pythonscad: Callable[..., str], tmp_path: Path
) -> None:
    """Export a translated cube as STL."""
    out_path = tmp_path / "translated.stl"
    program = (
        "from pythonscad import *\n"
        "c = cube(10).translate([5, 10, 15])\n"
        f"export(c, {str(out_path)!r})\n"
        "print('EXPORTED')\n"
        "show(c)\n"
    )
    run_pythonscad(program)
    assert out_path.exists()
    assert out_path.stat().st_size > 0


def test_export_3mf_multi(
    run_pythonscad: Callable[..., str], tmp_path: Path
) -> None:
    """Export multi-object dict as 3MF."""
    out_path = tmp_path / "multi.3mf"
    program = (
        "from pythonscad import *\n"
        f"export({{'a': cube(10), 'b': sphere(r=5)}}, {str(out_path)!r})\n"
        "print('EXPORTED')\n"
        "show(cube(1))\n"
    )
    run_pythonscad(program)
    assert out_path.exists()
    assert out_path.stat().st_size > 0


# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------


def test_string_as_scalar_raises_typeerror(
    run_pythonscad: Callable[..., str],
) -> None:
    """Passing a string where a number is expected raises TypeError from C."""
    out = run_pythonscad(
        "from pythonscad import *\n"
        "try:\n"
        "    cube('not_a_number')\n"
        "    print('NO ERROR')\n"
        "except TypeError:\n"
        "    print('TYPERROR')\n"
        "show(cube(1))\n"
    )
    assert "TYPERROR" in out


def test_short_vector_rejected(run_pythonscad: Callable[..., str]) -> None:
    """A 2-element list where 3 elements are expected raises from C."""
    out = run_pythonscad(
        "from pythonscad import *\n"
        "try:\n"
        "    cube([10, 20])\n"
        "    print('NO ERROR')\n"
        "except (TypeError, ValueError):\n"
        "    print('ERROR')\n"
        "show(cube(1))\n"
    )
    assert "ERROR" in out


# ---------------------------------------------------------------------------
# PyOpenSCAD type: exists as a real Python object, not just a C type
# ---------------------------------------------------------------------------


def test_pyopenscad_type_importable(run_pythonscad: Callable[..., str]) -> None:
    """``PyOpenSCAD`` is importable from ``pythonscad`` and is a ``type``."""
    out = run_pythonscad(
        "from pythonscad import PyOpenSCAD, show\n"
        "from pythonscad import cube\n"
        "print(type(PyOpenSCAD))\n"
        "print(isinstance(PyOpenSCAD, type))\n"
        "show(cube(1))\n"
    )
    assert "True" in out


def test_pyopenscad_inheritance(run_pythonscad: Callable[..., str]) -> None:
    """``PyOpenSCAD`` can be used as a base class."""
    out = run_pythonscad(
        "from pythonscad import PyOpenSCAD, cube, show\n"
        "class MyShape(PyOpenSCAD):\n"
        "    pass\n"
        "obj = MyShape()\n"
        "print(type(obj).__name__)\n"
        "print('OK')\n"
        "show(cube(1))\n"
    )
    assert "OK" in out


def test_pyopenscad_subclass_isinstance(run_pythonscad: Callable[..., str]) -> None:
    """Subclass instances of ``PyOpenSCAD`` pass ``isinstance`` checks."""
    out = run_pythonscad(
        "from pythonscad import PyOpenSCAD, cube, show\n"
        "class MyShape(PyOpenSCAD):\n"
        "    pass\n"
        "obj = MyShape()\n"
        "print('IS_CHECK', isinstance(obj, PyOpenSCAD))\n"
        "print('OK')\n"
        "show(cube(1))\n"
    )
    assert "IS_CHECK True" in out or "IS_CHECK" in out


def test_openscad_legacy_alias(run_pythonscad: Callable[..., str]) -> None:
    """The legacy ``Openscad`` alias still works and is the same type."""
    out = run_pythonscad(
        "from pythonscad import PyOpenSCAD, Openscad, show, cube\n"
        "print('SAME', Openscad is PyOpenSCAD)\n"
        "show(cube(1))\n"
    )
    assert "SAME True" in out
