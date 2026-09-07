"""Integration tests that drive a built pythonscad binary.

These prove the pure-Python overlay cooperates with the real ``_openscad``
extension end-to-end. They skip when no binary is provided (see the
``run_pythonscad`` fixture in conftest).
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path


def test_basic_cube_renders(run_pythonscad: Callable[..., str]) -> None:
    """A simple ``cube(10)`` can be rendered and exported without errors."""
    run_pythonscad("from pythonscad import *\nshow(cube(10))\n")


def test_rounded_cube_renders(run_pythonscad: Callable[..., str]) -> None:
    """The pure-Python ``rounded_cube`` helper works end-to-end with real geometry.

    ``rounded_cube`` is composed of ``minkowski(cube(...), sphere(...))``,
    exercising the full geometry pipeline.
    """
    run_pythonscad("from pythonscad import *\nshow(rounded_cube(20, r=2))\n")


def test_basic_vector_usage(run_pythonscad: Callable[..., str]) -> None:
    """Basic Python list comparison works inside the embedded interpreter."""
    out = run_pythonscad(
        "from pythonscad import *\n"
        "c = cube(10)\n"
        "print('MARK', [1, 2, 3] == [1, 2, 3])\n"
        "show(c)\n"
    )
    assert "MARK True" in out


def test_multitool_exporter_writes_per_part_files(
    run_pythonscad: Callable[..., str], tmp_path: Path
) -> None:
    """``MultiToolExporter`` writes one STL per part when run inside the real interpreter.

    Verifies ``parts()`` returns the correct names and that the output files
    exist on disk.
    """
    prefix = str(tmp_path / "part-")
    program = (
        "from pythonscad import *\n"
        f"exp = MultiToolExporter({prefix!r}, '.stl')\n"
        "exp.append(('base', cube(10)))\n"
        "exp.append(('pin', cylinder(h=20, r=2)))\n"
        "exp.export()\n"
        "print('PARTS', [n for n, _ in exp.parts()])\n"
        "show(cube(1))\n"
    )
    out = run_pythonscad(program)
    assert "PARTS ['base', 'pin']" in out
    assert (tmp_path / "part-base.stl").exists()
    assert (tmp_path / "part-pin.stl").exists()
