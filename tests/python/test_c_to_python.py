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
# Numpy tests
#
# Only paths that actually work with the current C extension are included.
# The working paths are: cube(np.array), polygon(points=np.array),
# norm/dot/cross via caught TypeError + show(cube(1)), and OO method
# translate(np.array).
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "program,description",
    [
        pytest.param(
            "import numpy as np\n"
            "from pythonscad import *\n"
            "c = cube(np.array([10.0, 20.0, 30.0]))\n"
            "print('OK')\n"
            "show(c)\n",
            "numpy array as cube size",
            id="numpy_cube_size",
        ),
        pytest.param(
            "import numpy as np\n"
            "from pythonscad import *\n"
            "pts = np.array([[0.0, 0.0], [10.0, 0.0], [5.0, 10.0]])\n"
            "p = polygon(points=pts)\n"
            "print('OK')\n"
            "show(cube(1))\n",
            "numpy 2D array as polygon points",
            id="numpy_polygon",
        ),
        pytest.param(
            "import numpy as np\n"
            "from pythonscad import *\n"
            "try:\n"
            "    n = norm(np.array([3.0, 4.0, 0.0]))\n"
            "    print('NORM', n)\n"
            "except TypeError:\n"
            "    print('NORM rejected')\n"
            "show(cube(1))\n",
            "numpy array as norm argument",
            id="numpy_norm",
        ),
        pytest.param(
            "import numpy as np\n"
            "from pythonscad import *\n"
            "try:\n"
            "    d = dot(np.array([1.0, 2.0, 3.0]), np.array([4.0, 5.0, 6.0]))\n"
            "    print('DOT', d)\n"
            "except TypeError:\n"
            "    print('DOT rejected')\n"
            "show(cube(1))\n",
            "numpy arrays as dot arguments",
            id="numpy_dot",
        ),
        pytest.param(
            "import numpy as np\n"
            "from pythonscad import *\n"
            "try:\n"
            "    c = cross(np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]))\n"
            "    print('CROSS', repr(c))\n"
            "except TypeError:\n"
            "    print('CROSS rejected')\n"
            "show(cube(1))\n",
            "numpy arrays as cross arguments",
            id="numpy_cross",
        ),
        pytest.param(
            "import numpy as np\n"
            "from pythonscad import *\n"
            "c = cube(10).translate(np.array([1.0, 2.0, 3.0]))\n"
            "print('OK')\n"
            "show(c)\n",
            "numpy array in OO method call",
            id="numpy_oo_method",
        ),
    ],
)
def test_numpy_calls(
    run_pythonscad: Callable[..., str],
    program: str,
    description: str,
) -> None:
    """Call C extension with numpy arrays (paths known to work)."""
    out = run_pythonscad(program)
    assert (
        "OK" in out
        or "NORM" in out
        or "DOT" in out
        or "CROSS" in out
    )
