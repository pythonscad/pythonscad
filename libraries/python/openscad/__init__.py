"""OpenSCAD-compatible Python API.

This package re-exports the C extension :mod:`_openscad` 1:1 and is the
home for any pure-Python reimplementations or compatibility shims that
should match upstream OpenSCAD's API surface.

PythonSCAD-only additions live in the :mod:`pythonscad` package, which
re-exports :mod:`openscad` and adds extra functionality on top.

The three-module layout is:

* :mod:`_openscad`  - C extension (low level, implementation detail).
* :mod:`openscad`   - this module: ``_openscad`` + OpenSCAD-compatible
  pure-Python additions/overrides.
* :mod:`pythonscad` - ``openscad`` + PythonSCAD-only extensions.

Switching a script between ``from openscad import *`` and
``from pythonscad import *`` requires no other code changes.
"""

import functools as _functools
import warnings as _warnings

# `import _openscad` (in addition to the star-import below) binds the name
# `_openscad` at module scope so the documented per-symbol deprecation recipe
# in `doc/python-modules.md` works as a literal copy-paste:
#     foo = _deprecated("foo", replacement="foo")(_openscad.foo)
import _openscad  # noqa: F401
from _openscad import *  # noqa: F401,F403
from _openscad import (  # noqa: F401
    ChildIterator,
    ChildRef,
    Openscad,
)

# Vector/matrix helper types. These live in the ``openscad._vectors``
# submodule (so they ship inside this package in the pip wheel) and are
# NumPy-backed when NumPy is installed, list-backed otherwise. Re-exported
# here so scripts can write ``from openscad import Vector3`` /
# ``from pythonscad import *``.
from ._vectors import (  # noqa: F401
    HAS_NUMPY,
    Matrix4x4,
    Vector1,
    Vector2,
    Vector3,
)


def _deprecated(name, replacement=None):
    """Decorator wrapping a callable so that calling it emits a
    :class:`DeprecationWarning` once per call site.

    Intended for future per-symbol deprecation of legacy aliases that
    should be pruned from the OpenSCAD-compatible API. Not currently
    in use; keeping the helper available so deprecations can be added
    without touching :mod:`_openscad` or :mod:`pythonscad`.

    See :doc:`doc/python-modules` for the documented recipe and how to
    keep ``pythonscad`` users from seeing the warning.
    """

    def _wrap(fn):
        @_functools.wraps(fn)
        def _inner(*args, **kwargs):
            msg = "openscad." + name + " is deprecated"
            if replacement:
                msg += "; use pythonscad." + replacement
            _warnings.warn(msg, DeprecationWarning, stacklevel=2)
            return fn(*args, **kwargs)

        return _inner

    return _wrap
