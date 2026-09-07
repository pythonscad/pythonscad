"""Unit tests for the pure-arithmetic helpers in ``pylibfive``.

The module imports ``libfive``; only a handful of primitives (``min``/``max``/
``sqrt``) are needed by the numeric helpers, so a tiny ``libfive`` stub backed
by real math is injected. The functions that build libfive CSG trees are left
to callers that actually have libfive.
"""

from __future__ import annotations

import math
import types

import pytest


@pytest.fixture
def lv(overlay: types.FunctionType) -> types.ModuleType:
    """Import ``pylibfive`` with a minimal ``libfive`` stub backed by ``math``.

    The stub provides ``min``, ``max``, ``sqrt`` and no-op ``x``/``y``/``z``
    identity functions so the vector/matrix helpers can be exercised without a
    real libfive installation.
    """
    stub = types.ModuleType("libfive")
    stub.min = min  # type: ignore[attr-defined]
    stub.max = max  # type: ignore[attr-defined]
    stub.sqrt = math.sqrt  # type: ignore[attr-defined]
    stub.x = lambda: 0.0  # type: ignore[attr-defined]
    stub.y = lambda: 0.0  # type: ignore[attr-defined]
    stub.z = lambda: 0.0  # type: ignore[attr-defined]
    return overlay("pylibfive", extra_stubs={"libfive": stub})


def test_dot_and_scalar(lv: types.ModuleType) -> None:
    """``lv_dot`` and ``lv_scalar`` both compute the dot product of two 3-vectors."""
    assert lv.lv_dot([1, 2, 3], [4, 5, 6]) == 32
    assert lv.lv_scalar([1, 2, 3], [4, 5, 6]) == 32


def test_vecsub(lv: types.ModuleType) -> None:
    """``lv_vecsub`` subtracts two 3-vectors component-wise."""
    assert lv.lv_vecsub([1, 2, 3], [1, 1, 1]) == (0, 1, 2)


def test_lerp_scalar_and_vector(lv: types.ModuleType) -> None:
    """``lv_lerp`` linearly interpolates between scalars or component-wise between vectors."""
    assert lv.lv_lerp(0, 10, 0.5) == 5
    assert lv.lv_lerp([0, 0, 0], [10, 20, 30], 0.5) == (5, 10, 15)


def test_vec_scale(lv: types.ModuleType) -> None:
    """``lv_vec_scale`` multiplies each component of a 3-vector by a scalar."""
    assert lv.lv_vec_scale([1, 2, 3], 2) == (2, 4, 6)


def test_len_and_length_use_real_sqrt(lv: types.ModuleType) -> None:
    """``lv_len`` / ``lv_length`` compute Euclidean norm using ``math.sqrt``."""
    assert lv.lv_len([3, 4, 0]) == pytest.approx(5.0)
    assert lv.lv_length([3, 4, 0]) == pytest.approx(5.0)


def test_vec_unit_normalises(lv: types.ModuleType) -> None:
    """``lv_vec_unit`` returns a unit vector (length 1) in the same direction."""
    x, y, z = lv.lv_vec_unit([3, 4, 0])
    assert (x, y, z) == pytest.approx((0.6, 0.8, 0.0))


@pytest.mark.parametrize("value,expected", [(-1, 0), (2, 2), (5, 3)])
def test_clamp(lv: types.ModuleType, value: int, expected: int) -> None:
    """``lv_clamp`` restricts a value to the closed interval ``[lo, hi]``."""
    assert lv.lv_clamp(value, 0, 3) == expected


def test_lv_rotyz_has_a_known_typo_bug(lv: types.ModuleType) -> None:
    """``lv_rotyz`` currently references an undefined variable ``sp``.

    This documents an existing bug where the body references the undefined
    name ``sp`` (``p[2]*c-p[1]*sp[1]*c+p[2]*s``). If this test ever starts
    passing, the typo was fixed and the test should be replaced with a real
    rotation-about-YZ assertion.
    """
    with pytest.raises(NameError):
        lv.lv_rotyz([1, 2, 3], 90)
