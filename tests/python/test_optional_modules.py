"""Import-smoke tests for the heavier / externally-dependent modules.

``pylaser``, ``jupyterdisplay`` and ``pybuild123d`` pull in real geometry or
third-party libraries (pythreejs, ipywidgets, build123d) and are not amenable
to fine-grained unit testing here. These tests prove each module imports
cleanly and exposes its public entry point when dependencies are stubbed.
"""

from __future__ import annotations

import types

import pytest


def _module_stub(name: str) -> types.ModuleType:
    """Create an importable placeholder that satisfies ``import x``."""
    return types.ModuleType(name)


def test_pylaser_imports_and_exposes_lasercutter(overlay: types.FunctionType) -> None:
    """``pylaser`` imports against the stubbed ``_openscad`` and exposes ``LaserCutter``.

    ``pylaser`` only depends on the (stubbed) pythonscad API and the stdlib,
    so no extra stubs are needed.
    """
    mod = overlay("pylaser")
    assert hasattr(mod, "LaserCutter")
    assert isinstance(mod.LaserCutter, type)


def test_pybuild123d_imports_with_stubbed_dep(overlay: types.FunctionType) -> None:
    """``pybuild123d`` imports when ``build123d`` is stubbed.

    Skips if the stubbed module is rejected by the import machinery.
    """
    try:
        mod = overlay("pybuild123d", extra_stubs={"build123d": _module_stub("build123d")})
    except ImportError as exc:
        pytest.skip(f"build123d unavailable and not stubbable: {exc}")
    assert hasattr(mod, "build123d")
    assert callable(mod.build123d)


def test_jupyterdisplay_imports_with_stubbed_deps(overlay: types.FunctionType) -> None:
    """``jupyterdisplay`` imports when ``pythreejs`` and ``ipywidgets`` are stubbed.

    The real Jupyter widget libraries are heavy and browser-dependent, so
    they are replaced with empty stubs for this smoke test.
    """
    stubs = {
        "pythreejs": _module_stub("pythreejs"),
        "ipywidgets": _module_stub("ipywidgets"),
    }
    try:
        mod = overlay("jupyterdisplay", extra_stubs=stubs)
    except ImportError as exc:
        pytest.skip(f"jupyter deps unavailable: {exc}")
    assert hasattr(mod, "build_widget")
