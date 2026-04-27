"""Type stubs for the `pythonscad` package.

The `pythonscad` package is a strict superset of `openscad` (which itself
re-exports `_openscad`). PythonSCAD-only additions are surfaced here.
"""
import typing as _typing

from openscad import *  # noqa: F401,F403
from openscad import (  # noqa: F401
    Color,
    Matrix4x4,
    PyLibFive,
    PyOpenSCAD,
    PyOpenSCADs,
    Vector2,
    Vector3,
)


class MultiToolExporter(list[tuple[_typing.Any, str]]):
    """List-based helper for exporting multi-tool / multi-color 3D models.

    Each item is an ``(object, name)`` 2-tuple. For each index ``i``,
    :meth:`export` writes the geometry obtained by subtracting all later
    objects from ``self[i][0]`` into the file ``f"{prefix}{name}{suffix}"``.
    The last entry is emitted as-is (no degenerate one-child ``difference``
    node). Names must be unique; duplicates raise :class:`ValueError` at
    export time.
    """

    prefix: str
    """String prepended to each output filename."""

    suffix: str
    """String appended to each output filename (typically the file extension)."""

    mkdir: bool
    """If True, create each output file's directory before exporting."""

    def __init__(
        self,
        prefix: str,
        suffix: str,
        mkdir: bool = ...,
        items: _typing.Iterable[tuple[_typing.Any, str]] = ...,
    ) -> None:
        """Initialize the exporter, optionally seeding it with ``items``."""
        ...

    def append(self, item: tuple[_typing.Any, str]) -> None:
        """Append a validated ``(object, name)`` tuple."""
        ...

    def extend(self, items: _typing.Iterable[tuple[_typing.Any, str]]) -> None:
        """Append each validated ``(object, name)`` tuple from ``items``."""
        ...

    def insert(self, index: _typing.SupportsIndex, item: tuple[_typing.Any, str]) -> None:
        """Insert a validated ``(object, name)`` tuple at ``index``."""
        ...

    def parts(self) -> list[tuple[str, _typing.Any]]:
        """Return computed ``(name, geometry)`` pairs in declaration order."""
        ...

    def export(self) -> None:
        """Export each part to a file via PythonSCAD."""
        ...

    def show(self) -> None:
        """Display each part in the PythonSCAD preview."""
        ...
