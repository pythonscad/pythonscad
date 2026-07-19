"""Unit tests for the geometry-free helpers in ``pytexture``.

``find_face`` / ``face_bbox`` / ``solid_bbox`` only read ``.matrix`` and
``.mesh()``; they are tested with lightweight fakes. ``load_texture`` /
``add_texture`` drive real geometry (surface/offset/roof/boolean ops) and are
left to the integration layer.
"""

from __future__ import annotations

import types

import pytest


@pytest.fixture
def pt(overlay: types.FunctionType) -> types.ModuleType:
    """Import the ``pytexture`` overlay module against a stubbed ``_openscad``."""
    return overlay("pytexture")


class FakeFace:
    """Minimal stand-in for a PythonSCAD face: provides ``.matrix`` and ``.mesh()``."""

    def __init__(self, matrix: object = None, points: object = None) -> None:
        self.matrix = matrix
        self._points = points

    def mesh(self) -> tuple[object, None]:
        return (self._points, None)


def _matrix_with_normal(nx: float, ny: float, nz: float) -> list[list[float]]:
    """Build a 4x4 matrix whose third column (rows 0..2) is the given face
    normal, mimicking the structure of a real OpenSCAD transformation matrix
    returned by ``face.matrix``.
    """
    return [
        [1, 0, nx, 0],
        [0, 1, ny, 0],
        [0, 0, nz, 0],
        [0, 0, 0, 1],
    ]


def test_find_face_picks_the_one_aligned_with_normal(pt: types.ModuleType) -> None:
    """``find_face`` returns the face whose normal most closely matches the given vector."""
    up = FakeFace(matrix=_matrix_with_normal(0, 0, 1))
    right = FakeFace(matrix=_matrix_with_normal(1, 0, 0))
    faces = [up, right]
    assert pt.find_face(faces, [0, 0, 1]) is up
    assert pt.find_face(faces, [1, 0, 0]) is right


def test_face_bbox_returns_min_max_xy(pt: types.ModuleType) -> None:
    """``face_bbox`` returns ``(min_x, min_y, max_x, max_y)`` of a 2D face's vertices."""
    face = FakeFace(points=[[0, 0], [10, 0], [10, 5], [0, 5]])
    assert pt.face_bbox(face) == (0, 0, 10, 5)


def test_solid_bbox_returns_min_max_xyz(pt: types.ModuleType) -> None:
    """``solid_bbox`` returns ``(min_x, min_y, min_z, max_x, max_y, max_z)`` of a solid's vertices."""
    solid = FakeFace(points=[[-1, -2, -3], [4, 5, 6], [0, 0, 0]])
    assert pt.solid_bbox(solid) == (-1, -2, -3, 4, 5, 6)
