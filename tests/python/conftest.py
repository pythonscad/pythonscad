"""pytest configuration for the PythonSCAD Python test-suite.

The Python side of PythonSCAD comes in two testable flavours, and this
file provides the shared machinery for both:

* **Overlay unit tests** for the pure-Python code that lives *inside* the
  ``openscad`` / ``pythonscad`` packages and the top-level helper modules
  (``pyutil``, ``pytexture``, ``pymachineconfig`` ...). These modules do
  ``from pythonscad import *`` / ``import _openscad`` at import time, so they
  cannot normally be imported without a build. The :func:`overlay` fixture
  installs a stub ``_openscad`` (every geometry name resolves to a
  :class:`~unittest.mock.MagicMock`) plus any extra library stubs the module
  needs, then imports it fresh. That makes the module's own Python logic --
  validation, path handling, arithmetic -- unit-testable in isolation while
  the geometry calls become inspectable mocks.

* **Integration tests** that drive a built ``pythonscad`` binary in a
  subprocess via :func:`run_pythonscad`. The binary path comes from
  ``--pythonscad-binary`` or ``PYTHONSCAD_BINARY`` (CMake sets the latter);
  without it those tests skip, so ``pytest tests/python`` still works from a
  source checkout with no build.
"""

from __future__ import annotations

import importlib
import os
import subprocess
import sys
import types
from collections.abc import Callable
from pathlib import Path
from unittest.mock import MagicMock

import pytest


_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir, os.pardir))
PYLIB_DIR = os.path.join(_REPO_ROOT, "libraries", "python")

#: Integration scripts run a full render; give them room on slow CI workers.
RUN_TIMEOUT_SECONDS = 120

#: Geometry / API names the overlay modules pull in via ``from pythonscad
#: import *`` (or ``from _openscad import *``) and then call as bare globals.
#: The stub exposes them as mocks so importing the overlay succeeds and the
#: calls can be asserted. Any name not listed here is still served lazily by
#: the stub's ``__getattr__``; the explicit list is what ``import *`` copies.
# The C extension exports these names via ``PyOpenSCADFunctions[]``.
# Keep this list in sync with ``src/python/pyfunctions.cc``.
_OPENSCAD_API_NAMES: list[str] = [
    # Primitives
    "edge", "square", "circle", "polygon", "polyline", "spline",
    "text", "textmetrics",
    "cube", "cylinder", "sphere", "polyhedron",
    # Transforms
    "translate", "right", "left", "back", "front", "up", "down",
    "rotx", "roty", "rotz", "rotate",
    "scale", "mirror", "multmatrix", "divmatrix", "offset",
    "color",
    # CSG / modelling ops
    "union", "difference", "intersection", "hull", "minkowski", "fill",
    "resize", "concat",
    # Analysis
    "bbox", "size", "position", "faces", "edges", "children",
    "inside", "mesh", "separate",
    # Extrusions
    "linear_extrude", "rotate_extrude", "path_extrude", "skin",
    # Modifiers
    "highlight", "background", "only",
    "projection", "surface",
    "explode", "oversample",
    "debug", "repair", "fillet", "render", "pull", "wrap",
    "group", "organic", "sheet",
    # IO
    "output", "show", "export",
    "osimport", "osuse", "osinclude",
    # System / info
    "version", "version_num", "version_string",
    "machineconfig", "rendervars",
    "model", "modelpath",
    "memberfunction", "marked",
    "scad", "align", "add_parameter",
    # Math functions (capitalised for OpenSCAD compat)
    "Sin", "Cos", "Tan", "Asin", "Acos", "Atan",
    "norm", "dot", "cross", "vector",
    # object protocol types
    "ChildIterator", "ChildRef", "Openscad",
]


def make_openscad_stub() -> types.ModuleType:
    """Build a stand-in ``_openscad`` module whose every symbol is a mock.

    ``import *`` copies the names in :data:`_OPENSCAD_API_NAMES`; anything
    else (e.g. the explicit ``ChildIterator`` import in ``openscad`` or a
    name only some overlay uses) is served lazily and memoised by
    ``__getattr__`` so repeated access returns the *same* mock and tests can
    assert on it.
    """
    mod = types.ModuleType("_openscad")
    for name in _OPENSCAD_API_NAMES:
        setattr(mod, name, MagicMock(name=f"_openscad.{name}"))
    mod.__all__ = list(_OPENSCAD_API_NAMES)  # type: ignore[attr-defined]

    def _getattr(name: str) -> MagicMock:
        mock = MagicMock(name=f"_openscad.{name}")
        setattr(mod, name, mock)
        return mock

    mod.__getattr__ = _getattr  # type: ignore[method-assign]  # PEP 562
    return mod


# Modules that must be re-imported from scratch whenever the stubbed
# ``_openscad`` changes, so a previous import can't leak a real/other stub.
_OVERLAY_MODULES: tuple[str, ...] = (
    "openscad",
    "pythonscad",
    "pyutil",
    "pytexture",
    "pymachineconfig",
    "pylibfive",
    "pylaser",
    "pybuild123d",
    "jupyterdisplay",
)


@pytest.fixture
def overlay(monkeypatch: pytest.MonkeyPatch) -> Callable[..., types.ModuleType]:
    """Import a PythonSCAD overlay/helper module against a stubbed ``_openscad``.

    Returns a callable ``import_overlay(module_name, extra_stubs=None)`` that
    installs the stub (plus any ``extra_stubs`` such as ``libfive``), imports
    the module fresh and returns it. The stubbed ``_openscad`` module is also
    returned via the fixture's ``.stub`` attribute for call assertions.
    """
    monkeypatch.syspath_prepend(PYLIB_DIR)
    stub = make_openscad_stub()
    monkeypatch.setitem(sys.modules, "_openscad", stub)
    for name in _OVERLAY_MODULES:
        monkeypatch.delitem(sys.modules, name, raising=False)

    def import_overlay(
        module_name: str,
        extra_stubs: dict[str, types.ModuleType] | None = None,
    ) -> types.ModuleType:
        for stub_name, stub_mod in (extra_stubs or {}).items():
            monkeypatch.setitem(sys.modules, stub_name, stub_mod)
        monkeypatch.delitem(sys.modules, module_name, raising=False)
        return importlib.import_module(module_name)

    import_overlay.stub = stub  # type: ignore[attr-defined]
    return import_overlay  # type: ignore[return-value]


# --------------------------------------------------------------------------
# Integration: drive a built pythonscad binary (skips without one).
# --------------------------------------------------------------------------
def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--pythonscad-binary",
        action="store",
        default=None,
        help=(
            "Path to the built pythonscad binary used by the integration "
            "tests. Defaults to the PYTHONSCAD_BINARY environment variable."
        ),
    )


@pytest.fixture(scope="session")
def pythonscad_binary(request: pytest.FixtureRequest) -> str:
    """Absolute path to the pythonscad binary, or skip if unavailable."""
    path = request.config.getoption("--pythonscad-binary") or os.environ.get(
        "PYTHONSCAD_BINARY"
    )
    if not path:
        pytest.skip(
            "no pythonscad binary given "
            "(pass --pythonscad-binary or set PYTHONSCAD_BINARY)"
        )
    if not os.path.isfile(path):
        pytest.skip(f"pythonscad binary not found: {path}")
    return os.path.abspath(path)


@pytest.fixture
def run_pythonscad(
    pythonscad_binary: str, tmp_path: Path
) -> Callable[..., str]:
    """Run a Python program inside the pythonscad interpreter.

    Returns the program's combined stdout+stderr (pythonscad routes the
    embedded interpreter's ``print()`` to stderr). An output file is always
    passed so the binary runs headless and exits instead of opening the GUI;
    the program is expected to ``show()`` something for that export to work.
    """

    def _run(program: str, timeout: int = RUN_TIMEOUT_SECONDS) -> str:
        script = tmp_path / "program.py"
        script.write_text(program, encoding="utf-8")
        out_stl = tmp_path / "out.stl"
        proc = subprocess.run(
            [pythonscad_binary, "--trust-python", str(script), "-o", str(out_stl)],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        combined = (proc.stdout or "") + (proc.stderr or "")
        for line in combined.splitlines():
            if line.startswith("SKIP:"):
                pytest.skip(line[len("SKIP:"):].strip())
        assert proc.returncode == 0, (
            f"pythonscad exited with {proc.returncode}\n--- output ---\n{combined}"
        )
        return combined

    return _run
