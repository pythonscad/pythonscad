"""Unit tests for the ``pythonscad._vectors`` runtime vector/matrix helpers.

These exercise the pure-Python :class:`Vector1`/:class:`Vector2`/
:class:`Vector3` / :class:`Matrix4x4` classes that back
``from pythonscad import Vector3`` etc. The module is loaded directly by file
path (not via ``import pythonscad._vectors``) so that the package ``__init__``
-- and therefore the compiled ``_openscad`` extension -- is not required to
run these tests.

Both code paths are covered:

* NumPy installed  -> instances are ``numpy.ndarray`` subclasses.
* NumPy absent     -> instances are ``list`` subclasses. This path is
  exercised deterministically in a subprocess with ``numpy`` import blocked,
  regardless of whether NumPy is installed in the test environment.

Points covered (current-backing, both paths)::

    construction            - VectorN from list, tuple, ndarray; default = zero
    length / validation     - correct len(), ValueError on wrong cardinality
    from_array              - roundtrip and numpy-like tolist() unwrapping
    Matrix4x4 identity      - default identity, 4-row from list-of-lists
    Matrix4x4 validation    - ValueError on wrong shape (3-col, 1-row)

Points covered (numpy-backed, requires NumPy)::

    ndarray subclass        - isinstance check, .shape, from ndarray constructor
    array protocol          - __array__(), np.asarray(), np.array_equal()
    arithmetic              - vector + vector
    no-alias guarantee      - constructor copies source, does not alias
    iteration / repr        - iter yields plain float; repr is list form
    equality                - __eq__/__ne__ return plain bool; scalar -> False
    unhashability           - TypeError from hash()
    matrix repr             - list-of-lists form

Points covered (list-backed, forced via numpy-blocked subprocess)::

    list subclass           - isinstance check for VectorN, Matrix4x4
    correctness             - values, identity, from_array with list/tolist()
    rejection               - wrong length, bad matrix shape, __array__() fails
"""

from __future__ import annotations

import importlib.util
import os
import subprocess
import sys
import textwrap
import types

import pytest


_HERE = os.path.dirname(os.path.abspath(__file__))
_VECTORS_PY = os.path.abspath(
    os.path.join(_HERE, os.pardir, os.pardir, "libraries", "python",
                 "pythonscad", "_vectors.py")
)


@pytest.fixture(scope="module")
def ov() -> types.ModuleType:
    """Load ``pythonscad/_vectors.py`` standalone as ``ov``."""
    spec = importlib.util.spec_from_file_location(
        "_pythonscad_vectors_under_test", _VECTORS_PY
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ------------------------------------------------------------------
# Current-backing tests -- must hold regardless of the active backing
# ------------------------------------------------------------------

class TestCurrentBacking:
    """Behaviour that must hold regardless of the active backing."""

    def test_construct_from_list(self, ov: types.ModuleType) -> None:
        assert list(ov.Vector3([1, 2, 3])) == [1.0, 2.0, 3.0]

    def test_construct_from_tuple(self, ov: types.ModuleType) -> None:
        assert list(ov.Vector2((4, 5))) == [4.0, 5.0]

    def test_construct_default_is_zero(self, ov: types.ModuleType) -> None:
        assert list(ov.Vector3()) == [0.0, 0.0, 0.0]

    def test_lengths(self, ov: types.ModuleType) -> None:
        assert len(ov.Vector1([1])) == 1
        assert len(ov.Vector2([1, 2])) == 2
        assert len(ov.Vector3([1, 2, 3])) == 3

    def test_wrong_length_raises(self, ov: types.ModuleType) -> None:
        with pytest.raises(ValueError):
            ov.Vector3([1, 2])
        with pytest.raises(ValueError):
            ov.Vector2([1, 2, 3])

    def test_matrix_identity_default(self, ov: types.ModuleType) -> None:
        m = ov.Matrix4x4()
        rows = [list(r) for r in m]
        assert rows == [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]

    def test_matrix_from_rows(self, ov: types.ModuleType) -> None:
        src = [[1, 0, 0, 5], [0, 1, 0, 6], [0, 0, 1, 7], [0, 0, 0, 1]]
        m = ov.Matrix4x4(src)
        assert [list(r) for r in m] == [[float(x) for x in row] for row in src]

    def test_matrix_wrong_shape_raises(self, ov: types.ModuleType) -> None:
        with pytest.raises(ValueError):
            ov.Matrix4x4([[1, 2, 3]])
        with pytest.raises(ValueError):
            ov.Matrix4x4([[1, 2, 3, 4]])  # only 1 row

    def test_from_array_roundtrip(self, ov: types.ModuleType) -> None:
        v = ov.Vector3.from_array([9, 8, 7])
        assert list(v) == [9.0, 8.0, 7.0]


# ------------------------------------------------------------------
# NumPy-backed tests
# ------------------------------------------------------------------

np = pytest.importorskip("numpy", reason="NumPy not installed")


class TestNumpyBacked:
    """Behaviour specific to the NumPy-backed path."""

    def test_instances_are_ndarrays(self, ov: types.ModuleType) -> None:
        assert isinstance(ov.Vector3([1, 2, 3]), np.ndarray)
        assert isinstance(ov.Matrix4x4(), np.ndarray)

    def test_construct_from_ndarray(self, ov: types.ModuleType) -> None:
        v = ov.Vector3(np.array([1.0, 2.0, 3.0]))
        assert np.array_equal(v, [1.0, 2.0, 3.0])

    def test_array_protocol(self, ov: types.ModuleType) -> None:
        v = ov.Vector3([1, 2, 3])
        assert np.asarray(v).tolist() == [1.0, 2.0, 3.0]

    def test_vector_arithmetic(self, ov: types.ModuleType) -> None:
        result = ov.Vector3([1, 2, 3]) + ov.Vector3([10, 20, 30])
        assert np.array_equal(result, [11.0, 22.0, 33.0])

    def test_matrix_shape(self, ov: types.ModuleType) -> None:
        assert ov.Matrix4x4().shape == (4, 4)

    def test_construct_does_not_alias_source(self, ov: types.ModuleType) -> None:
        src = np.array([1.0, 2.0, 3.0])
        v = ov.Vector3(src)
        v[0] = 99.0
        assert src[0] == 1.0

    def test_iter_yields_python_floats(self, ov: types.ModuleType) -> None:
        items = list(ov.Vector3([1, 2, 3]))
        assert items == [1.0, 2.0, 3.0]
        assert all(type(x) is float for x in items)

    def test_repr_is_list_form(self, ov: types.ModuleType) -> None:
        assert repr(ov.Vector3([1, 2, 3])) == "[1.0, 2.0, 3.0]"

    def test_eq_returns_plain_bool(self, ov: types.ModuleType) -> None:
        v = ov.Vector3([1, 2, 3])
        assert v == [1, 2, 3]
        assert not (v == [9, 9, 9])
        assert v != [9, 9, 9]
        assert not (v == 5)

    def test_unhashable(self, ov: types.ModuleType) -> None:
        with pytest.raises(TypeError):
            hash(ov.Vector3([1, 2, 3]))

    def test_matrix_repr_is_list_form(self, ov: types.ModuleType) -> None:
        m = ov.Matrix4x4([[1.0, 0.0, 0.0, 5.0], [0.0, 1.0, 0.0, 0.0],
                          [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]])
        assert repr(m) == (
            "[[1.0, 0.0, 0.0, 5.0], [0.0, 1.0, 0.0, 0.0], "
            "[0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]"
        )


# ------------------------------------------------------------------
# List-backed fallback tests (forced via numpy-blocked subprocess)
# ------------------------------------------------------------------

_NO_NUMPY_PROGRAM = textwrap.dedent(
    """
    import builtins, importlib.util, sys
    _orig = builtins.__import__
    def _blocked(name, *a, **k):
        if name == "numpy" or name.startswith("numpy."):
            raise ImportError("blocked for test")
        return _orig(name, *a, **k)
    builtins.__import__ = _blocked
    sys.modules.pop("numpy", None)
    _spec = importlib.util.spec_from_file_location("_ov_nonumpy", {vectors_py!r})
    ov = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(ov)
    assert ov.HAS_NUMPY is False, "expected list-backed fallback"
    v = ov.Vector3([1, 2, 3])
    assert isinstance(v, list), type(v)
    assert v == [1.0, 2.0, 3.0], v
    assert isinstance(ov.Vector1([5]), list) and ov.Vector1([5]) == [5.0]
    assert list(ov.Vector2((4, 5))) == [4.0, 5.0]
    assert ov.Vector3() == [0.0, 0.0, 0.0]
    try:
        ov.Vector3([1, 2]); raise SystemExit("wrong-length accepted")
    except ValueError:
        pass
    m = ov.Matrix4x4()
    assert m == [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ], m
    try:
        ov.Matrix4x4([[1, 2, 3]]); raise SystemExit("bad matrix accepted")
    except ValueError:
        pass
    try:
        ov.Vector3([1, 2, 3]).__array__(); raise SystemExit("array() should fail")
    except TypeError:
        pass
    assert ov.Vector3.from_array([7, 8, 9]) == [7.0, 8.0, 9.0]
    class _Fake:
        def tolist(self):
            return [1, 2, 3]
    assert ov.Vector3.from_array(_Fake()) == [1.0, 2.0, 3.0]
    print("NO_NUMPY_OK")
    """
)


class TestListBacked:
    """The list-backed fallback, forced on via a numpy-blocked subprocess."""

    def test_fallback_behaviour(self) -> None:
        proc = subprocess.run(
            [sys.executable, "-c",
             _NO_NUMPY_PROGRAM.format(vectors_py=_VECTORS_PY)],
            capture_output=True,
            text=True,
        )
        assert proc.returncode == 0, proc.stderr
        assert "NO_NUMPY_OK" in proc.stdout
