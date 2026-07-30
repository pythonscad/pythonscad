"""Unit tests for ``pyutil`` (the ``loft`` helpers).

``loft_func`` is pure trigonometry and is tested directly. ``loft`` /
``loft_prepare`` need objects with a ``.mesh()`` method, which are supplied as
lightweight fakes; the module is imported against a stubbed ``_openscad``.
"""

from __future__ import annotations

import math
import types

import pytest


@pytest.fixture
def pu(overlay: types.FunctionType) -> types.ModuleType:
    """Import the ``pyutil`` overlay module against a stubbed ``_openscad``."""
    return overlay("pyutil")


class FakeSolid:
    """Minimal stand-in for a PythonSCAD solid: only ``mesh()`` is used."""

    def __init__(self, points2d: list[list[float]]) -> None:
        self._points = points2d

    def mesh(self) -> tuple[list[list[float]], None]:
        # loft_prepare only reads mesh()[0]; the rest is unused here.
        return (self._points, None)


def test_loft_func_interpolates_endpoints(pu: types.ModuleType) -> None:
    """``loft_func`` linearly interpolates magnitude between the two shapes
    across height: at height 0 the result matches the first shape's radius;
    at ``loft_height`` it matches the second; halfway gives the average.
    """
    data = [[0.0], [1.0], [3.0]]
    assert pu.loft_func(data, 10, 0, 0)[0] == pytest.approx([1.0, 0.0])
    assert pu.loft_func(data, 10, 10, 0)[0] == pytest.approx([3.0, 0.0])
    assert pu.loft_func(data, 10, 5, 0)[0] == pytest.approx([2.0, 0.0])


def test_loft_func_applies_rotation(pu: types.ModuleType) -> None:
    """At maximum height the full rotation angle is applied to each point."""
    data = [[0.0], [1.0], [1.0]]
    x, y = pu.loft_func(data, 10, 10, math.pi / 2)[0]
    assert x == pytest.approx(0.0, abs=1e-9)
    assert y == pytest.approx(1.0)


def test_loft_func_point_count_matches_angles(pu: types.ModuleType) -> None:
    """``loft_func`` returns one result per input angle."""
    data = [[0.0, 1.0, 2.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]]
    assert len(pu.loft_func(data, 10, 5, 0)) == 3


def test_loft_returns_height_function_over_two_shapes(pu: types.ModuleType) -> None:
    """``loft`` returns a callable that produces interpolated cross-sections
    at any height; the point count is consistent across all heights.
    """
    square: list[list[float]] = [[-5, -5], [5, -5], [5, 5], [-5, 5]]
    fn = pu.loft(FakeSolid(square), FakeSolid(square), loft_height=10, n=8)
    assert callable(fn)
    bottom = fn(0)
    top = fn(10)
    assert len(bottom) == len(top)
    assert all(len(p) == 2 for p in bottom)


def test_loft_prepare_builds_angle_and_two_magnitude_rows(pu: types.ModuleType) -> None:
    """``loft_prepare`` returns three parallel lists (angles, magnitudes for
    shape 1, magnitudes for shape 2) with the angle list sorted ascending.
    """
    square: list[list[float]] = [[-5, -5], [5, -5], [5, 5], [-5, 5]]
    data = pu.loft_prepare(FakeSolid(square), FakeSolid(square), n=4, rot=0)
    assert len(data) == 3  # [angles, mags_shape1, mags_shape2]
    assert len(data[1]) == len(data[0])
    assert len(data[2]) == len(data[0])
    assert data[0] == sorted(data[0])
