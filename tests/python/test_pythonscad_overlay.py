"""Unit tests for the pure-Python logic in ``pythonscad/__init__.py``.

Covers the parts that do not depend on real geometry -- ``MultiToolExporter``
validation / uniqueness / part computation, ``_normalize_filename_key`` and
``rounded_cube`` argument checking -- by importing the overlay against a
stubbed ``_openscad`` so the geometry calls are inspectable mocks.
"""

from __future__ import annotations

import sys
import types

import pytest


@pytest.fixture
def ps(overlay: types.FunctionType) -> types.ModuleType:
    """Import the ``pythonscad`` overlay module against a stubbed ``_openscad``."""
    return overlay("pythonscad")


# --------------------------------------------------------------------------
# _normalize_filename_key
# --------------------------------------------------------------------------
def test_normalize_collapses_path_aliases(ps: types.ModuleType) -> None:
    """``_normalize_filename_key`` resolves ``..`` components like
    ``os.path.normpath`` so that ``"a/../b.stl"`` and ``"b.stl"`` produce the
    same key and are rejected as duplicates.
    """
    assert ps._normalize_filename_key("a/../b.stl") == ps._normalize_filename_key("b.stl")


def test_normalize_case_sensitivity_matches_platform(ps: types.ModuleType) -> None:
    """Case folding in ``_normalize_filename_key`` matches the platform's
    filesystem semantics -- case-insensitive on macOS/Windows, case-sensitive on Linux.
    """
    same = ps._normalize_filename_key("a.stl") == ps._normalize_filename_key("A.stl")
    if sys.platform in ("darwin", "win32"):
        assert same
    else:
        assert not same


# --------------------------------------------------------------------------
# MultiToolExporter -- validation
# --------------------------------------------------------------------------
def test_constructor_validates_initial_items(ps: types.ModuleType) -> None:
    """``MultiToolExporter(items=...)`` accepts valid ``(name, object)`` tuples."""
    exp = ps.MultiToolExporter("out/", ".stl", items=[("a", object()), ("b", object())])
    assert [name for name, _ in exp] == ["a", "b"]


@pytest.mark.parametrize(
    "bad",
    [
        ["a", object()],          # list, not tuple
        ("a",),                    # wrong arity
        ("a", object(), 1),        # wrong arity
        (1, object()),             # non-str name
    ],
)
def test_append_type_errors(ps: types.ModuleType, bad: object) -> None:
    """``append`` raises ``TypeError`` for non-tuple, wrong-arity or non-str-name items."""
    exp = ps.MultiToolExporter("out/", ".stl")
    with pytest.raises(TypeError):
        exp.append(bad)


def test_append_empty_name_is_value_error(ps: types.ModuleType) -> None:
    """``append`` raises ``ValueError`` when the name is an empty string."""
    exp = ps.MultiToolExporter("out/", ".stl")
    with pytest.raises(ValueError):
        exp.append(("", object()))


def test_extend_is_atomic(ps: types.ModuleType) -> None:
    """``extend`` validates all items before appending any of them; on error the exporter is unchanged."""
    exp = ps.MultiToolExporter("out/", ".stl")
    with pytest.raises(TypeError):
        exp.extend([("ok", object()), ("bad-not-tuple",)])
    assert len(exp) == 0


def test_iadd_validates(ps: types.ModuleType) -> None:
    """``+=`` validates all items before extending, leaving the exporter unchanged on error."""
    exp = ps.MultiToolExporter("out/", ".stl")
    with pytest.raises(TypeError):
        exp += [("bad",)]
    assert len(exp) == 0


def test_setitem_validates(ps: types.ModuleType) -> None:
    """Item assignment (``exp[i] = ...``) validates the replacement value."""
    exp = ps.MultiToolExporter("out/", ".stl", items=[("a", object())])
    with pytest.raises(TypeError):
        exp[0] = ("bad", object(), object())


# --------------------------------------------------------------------------
# MultiToolExporter -- part computation and export
# --------------------------------------------------------------------------
def test_parts_last_wins_difference(ps: types.ModuleType) -> None:
    """``parts()`` returns items in order; the last item is returned as-is
    while earlier items wrap their geometry in a ``difference()`` node
    that subtracts all later items ("last wins" overlap policy).
    """
    a, b, c = object(), object(), object()
    exp = ps.MultiToolExporter("out/", ".stl", items=[("a", a), ("b", b), ("c", c)])
    names = [n for n, _ in exp.parts()]
    assert names == ["a", "b", "c"]
    assert exp.parts()[-1][1] is c
    ps.difference.assert_called()


def test_export_per_file_writes_one_per_part(ps: types.ModuleType, tmp_path: pytest.TempPathFactory) -> None:
    """``export()`` writes one file per item using ``f"{prefix}{name}{suffix}"``."""
    exp = ps.MultiToolExporter(f"{tmp_path}/flag-", ".stl",
                               items=[("base", object()), ("overlay", object())])
    exp.export()
    assert ps.export.call_count == 2
    written = [call.args[1] for call in ps.export.call_args_list]
    assert written == [f"{tmp_path}/flag-base.stl", f"{tmp_path}/flag-overlay.stl"]


def test_export_rejects_colliding_output_paths(ps: types.ModuleType) -> None:
    """``export()`` raises ``ValueError`` when two items would write to the same file."""
    exp = ps.MultiToolExporter("out/", ".stl", items=[("a", object()), ("a", object())])
    with pytest.raises(ValueError):
        exp.export()


def test_export_single_file_requires_3mf(ps: types.ModuleType) -> None:
    """``single_file`` export only accepts ``.3mf`` extension, otherwise raises ``ValueError``."""
    exp = ps.MultiToolExporter("out/", ".stl", items=[("a", object())])
    with pytest.raises(ValueError):
        exp.export(single_file="out/model.stl")


def test_export_single_file_3mf_calls_multi_object_export(ps: types.ModuleType, tmp_path: pytest.TempPathFactory) -> None:
    """``export(single_file="out.3mf")`` calls the multi-object ``export(dict, path)`` form."""
    target = f"{tmp_path}/model.3mf"
    exp = ps.MultiToolExporter("out/", ".stl",
                               items=[("a", object()), ("b", object())])
    exp.export(single_file=target)
    ps.export.assert_called_once()
    payload, filename = ps.export.call_args.args
    assert filename == target
    assert isinstance(payload, dict) and set(payload) == {"a", "b"}


# --------------------------------------------------------------------------
# rounded_cube -- argument validation and geometry composition
# --------------------------------------------------------------------------
def test_rounded_cube_requires_exactly_one_of_r_d(ps: types.ModuleType) -> None:
    """``rounded_cube`` raises ``TypeError`` unless exactly one of ``r`` or ``d`` is given."""
    with pytest.raises(TypeError):
        ps.rounded_cube(20)
    with pytest.raises(TypeError):
        ps.rounded_cube(20, r=2, d=4)


def test_rounded_cube_positive_radius(ps: types.ModuleType) -> None:
    """``rounded_cube`` raises ``ValueError`` when ``r`` is zero or negative."""
    with pytest.raises(ValueError):
        ps.rounded_cube(20, r=0)


def test_rounded_cube_dimension_must_exceed_2r(ps: types.ModuleType) -> None:
    """Each outer dimension must be strictly greater than ``2 * radius``."""
    with pytest.raises(ValueError):
        ps.rounded_cube(3, r=2)
    with pytest.raises(ValueError):
        ps.rounded_cube([30, 20, 3], r=2)


def test_rounded_cube_bad_size_type(ps: types.ModuleType) -> None:
    """``rounded_cube`` raises ``TypeError`` for non-numeric or wrong-length size."""
    with pytest.raises(TypeError):
        ps.rounded_cube("nope", r=2)
    with pytest.raises(TypeError):
        ps.rounded_cube([1, 2], r=0.1)


def test_rounded_cube_builds_minkowski_of_shrunken_cube_and_sphere(ps: types.ModuleType) -> None:
    """``rounded_cube(size, r)`` composes ``minkowski(cube(size - 2r), sphere(r))`` then translates."""
    ps.rounded_cube(20, r=2)
    ps.cube.assert_called_once_with(16)
    ps.sphere.assert_called_once_with(r=2, fn=None, fa=None, fs=None)
    ps.minkowski.assert_called_once()
    ps.minkowski.return_value.translate.assert_called_once_with([2, 2, 2])


def test_rounded_cube_accepts_diameter(ps: types.ModuleType) -> None:
    """``rounded_cube(size, d=4)`` is equivalent to ``rounded_cube(size, r=2)``."""
    ps.rounded_cube([30, 20, 10], d=4)
    ps.cube.assert_called_once_with([26, 16, 6])
