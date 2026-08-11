"""Unit tests for the pure-Python bits of ``openscad/__init__.py``.

The ``_deprecated`` helper exercised against a stubbed ``_openscad``.
"""

from __future__ import annotations

import types
import warnings

import pytest


@pytest.fixture
def osc(overlay: types.FunctionType) -> types.ModuleType:
    """Import the ``openscad`` overlay module against a stubbed ``_openscad``."""
    return overlay("openscad")


def test_binds_underlying_extension_name(osc: types.ModuleType) -> None:
    """Verify ``import _openscad`` binds the C extension name at module scope.

    This binding is required by the documented per-symbol deprecation recipe
    in ``doc/python-modules.md``, which users copy-paste as::

        foo = _deprecated("foo", replacement="foo")(_openscad.foo)
    """
    assert hasattr(osc, "_openscad")


def test_deprecated_emits_warning_and_forwards(osc: types.ModuleType) -> None:
    """Calling a ``_deprecated``-wrapped function emits a DeprecationWarning.

    The warning message includes the old name and, when provided, the
    replacement hint. The underlying function is still called and its
    return value forwarded unchanged.
    """
    calls: list[tuple[int, int]] = []

    @osc._deprecated("old_name", replacement="new_name")
    def impl(a: int, b: int) -> int:
        calls.append((a, b))
        return a + b

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = impl(2, 3)

    assert result == 5
    assert calls == [(2, 3)]
    assert len(caught) == 1
    assert issubclass(caught[0].category, DeprecationWarning)
    msg = str(caught[0].message)
    assert "old_name is deprecated" in msg
    assert "new_name" in msg


def test_deprecated_preserves_wrapped_metadata(osc: types.ModuleType) -> None:
    """The ``_deprecated`` decorator preserves ``__name__`` and ``__doc__``."""
    @osc._deprecated("foo")
    def foo() -> None:
        """docstring stays."""

    assert foo.__name__ == "foo"
    assert foo.__doc__ == "docstring stays."
