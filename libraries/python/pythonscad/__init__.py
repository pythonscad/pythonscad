"""PythonSCAD public API.

Strict superset of :mod:`openscad`. Currently a 1:1 re-export; future
PythonSCAD-only additions (new classes, helper APIs) will land here as
either symbols in this package or as submodules.

This package re-exports :mod:`openscad` rather than :mod:`_openscad`
directly, so any pure-Python implementation or override added in
:mod:`openscad` is automatically picked up.

The three-module layout is:

* :mod:`_openscad`  - C extension (low level, implementation detail).
* :mod:`openscad`   - ``_openscad`` + OpenSCAD-compatible pure-Python
  additions/overrides.
* :mod:`pythonscad` - this module: ``openscad`` + PythonSCAD-only
  extensions.

Switching a script between ``from openscad import *`` and
``from pythonscad import *`` requires no other code changes.
"""

from openscad import *  # noqa: F401,F403
from openscad import (  # noqa: F401
    ChildIterator,
    ChildRef,
    Openscad,
)

import os as _os
import typing as _typing


class MultiToolExporter(list[tuple[_typing.Any, str]]):
    """List-based helper for exporting multi-tool / multi-color 3D models.

    Each item in the list is a ``(object, name)`` tuple, where ``object`` is a
    PythonSCAD geometry and ``name`` is a label used to build the output
    filename. The exporter is designed for workflows where a single model is
    split into several parts (for example, one part per filament color on a
    multi-tool 3D printer).

    Semantics
    ---------
    For each index ``i`` in the list, the exporter produces the geometry
    obtained by taking ``self[i][0]`` and subtracting all later objects
    ``self[i+1:][...]`` from it. Overlapping regions are therefore assigned to
    exactly one part: **later entries "win" over earlier ones**, so each part
    only keeps the volume not claimed by any subsequent part. The last entry
    is emitted as-is (no degenerate one-child ``difference`` node) and
    therefore "wins" everything that overlaps with it.

    Attributes:
        prefix: String prepended to each output filename. Typically a path
            and/or base name, e.g. ``"out/model-"``.
        suffix: String appended to each output filename, typically the file
            extension, e.g. ``".stl"`` or ``".3mf"``.
        mkdir: If ``True``, the directory portion of each output filename is
            created (with :func:`os.makedirs`) before exporting. Defaults to
            ``False``. Filenames without a directory component are exported
            as-is, no error is raised.

    Validation
    ----------
    Items are validated on insertion (constructor argument, :meth:`append`,
    :meth:`extend`, :meth:`insert`, ``self[i] = ...``). A :class:`TypeError`
    is raised if the item is not a 2-tuple of ``(object, str)``, and a
    :class:`ValueError` is raised if the name is empty.

    Duplicate names are detected at :meth:`export` time and raise
    :class:`ValueError` (rather than silently letting later entries clobber
    files written by earlier ones).

    Example:
        >>> # Append base/background parts first; later entries "win" overlap.
        >>> exporter = MultiToolExporter("out/flag-", ".stl", mkdir=True)
        >>> exporter.append((base_geometry, "base"))
        >>> exporter.append((overlay_geometry, "overlay"))
        >>> exporter.export()  # writes out/flag-base.stl and out/flag-overlay.stl
    """

    def __init__(
        self,
        prefix: str,
        suffix: str,
        mkdir: bool = False,
        items: _typing.Iterable[tuple[_typing.Any, str]] = (),
    ):
        """Initialize a (possibly empty) MultiToolExporter.

        Args:
            prefix: String prepended to each output filename.
            suffix: String appended to each output filename, usually the file
                extension.
            mkdir: If ``True``, create the output directory for each file
                before exporting. Defaults to ``False``.
            items: Optional iterable of initial ``(object, name)`` tuples.
                Each item is validated as if it were appended.

        Raises:
            TypeError: If any item in ``items`` is not a 2-tuple of
                ``(object, str)``.
            ValueError: If any name in ``items`` is empty.
        """
        super().__init__()
        self.prefix = prefix
        self.suffix = suffix
        self.mkdir = mkdir
        for item in items:
            self.append(item)

    @staticmethod
    def _validate_item(item: _typing.Any) -> tuple[_typing.Any, str]:
        """Return ``item`` if it is a valid ``(object, str)`` 2-tuple.

        Raises:
            TypeError: If ``item`` is not a 2-tuple whose second element is
                a string. Lists and other sequences are rejected.
            ValueError: If the name is empty.
        """
        if not isinstance(item, tuple) or len(item) != 2:
            raise TypeError(
                f"MultiToolExporter items must be (object, name) 2-tuples, "
                f"got {type(item).__name__}: {item!r}"
            )
        _, name = item
        if not isinstance(name, str):
            raise TypeError(
                f"MultiToolExporter item name must be a str, "
                f"got {type(name).__name__}: {name!r}"
            )
        if not name:
            raise ValueError("MultiToolExporter item name must be a non-empty string")
        return item

    def append(self, item: tuple[_typing.Any, str]) -> None:
        """Append a validated ``(object, name)`` tuple."""
        super().append(self._validate_item(item))

    def extend(self, items: _typing.Iterable[tuple[_typing.Any, str]]) -> None:
        """Append each ``(object, name)`` tuple from ``items``.

        Atomic: every item is validated *before* anything is appended, so a
        bad item in the middle of the iterable leaves the exporter unchanged.
        """
        validated = [self._validate_item(item) for item in items]
        super().extend(validated)

    def insert(self, index: _typing.SupportsIndex, item: tuple[_typing.Any, str]) -> None:
        """Insert a validated ``(object, name)`` tuple at ``index``."""
        super().insert(index, self._validate_item(item))

    def __setitem__(self, index, value):
        """Replace one or more items, validating each new ``(object, name)``."""
        if isinstance(index, slice):
            super().__setitem__(index, [self._validate_item(v) for v in value])
        else:
            super().__setitem__(index, self._validate_item(value))

    def _filename(self, i: int) -> str:
        """Return the output filename for part ``i``."""
        return f"{self.prefix}{self[i][1]}{self.suffix}"

    def _part(self, i: int):
        """Return the geometry for part ``i``: ``self[i][0]`` minus all later parts.

        For the last entry, the object is returned as-is rather than wrapped
        in a one-child ``difference`` node.
        """
        rest = [obj for obj, _name in self[i:]]
        return rest[0] if len(rest) == 1 else difference(*rest)  # noqa: F405

    def _check_unique_names(self) -> None:
        """Raise :class:`ValueError` if any two items share the same name.

        Duplicate names would produce duplicate filenames in :meth:`export`,
        causing later parts to overwrite earlier ones. We refuse rather than
        silently lose data.
        """
        seen = set()
        for _, name in self:
            if name in seen:
                raise ValueError(
                    f"MultiToolExporter has duplicate name {name!r}; "
                    f"each item must have a unique name to avoid overwriting files"
                )
            seen.add(name)

    def export(self) -> None:
        """Export each part to a file via PythonSCAD.

        For each index ``i``, exports the difference of ``self[i][0]`` and
        all subsequent objects to ``f"{prefix}{name}{suffix}"``. The last
        entry is exported as-is without a degenerate ``difference`` node.

        If :attr:`mkdir` is ``True``, the parent directory of each output
        file is created beforehand (filenames without a directory component
        are skipped silently).

        Raises:
            ValueError: If two or more items share the same name.
        """
        self._check_unique_names()
        for i in range(len(self)):
            filename = self._filename(i)
            if self.mkdir:
                directory = _os.path.dirname(filename)
                if directory:
                    _os.makedirs(directory, exist_ok=True)
            export(self._part(i), filename)  # noqa: F405

    def show(self) -> None:
        """Display each part in the PythonSCAD preview.

        For each index ``i``, computes the difference of ``self[i][0]`` and
        all subsequent objects and passes the result to :func:`show`,
        producing a layered preview equivalent to what :meth:`export` would
        write to disk.
        """
        for i in range(len(self)):
            show(self._part(i))  # noqa: F405
