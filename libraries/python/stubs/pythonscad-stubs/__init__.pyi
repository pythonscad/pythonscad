"""Type stubs for the `pythonscad` package.

The `pythonscad` package is a strict superset of `openscad` (which itself
re-exports `_openscad`). PythonSCAD-only additions are surfaced here.
"""

# Convention (mirrors libraries/python/pythonscad/__init__.py): any
# import that is NOT part of the public `pythonscad` stub must be
# aliased with a leading underscore (e.g. ``import typing as _typing``).
# Type-checkers treat unaliased names in a stub as part of the public
# API surface, so leaking helpers here is just as bad as leaking them
# at runtime.
import typing as _typing

try:
    import numpy as _np
except ImportError:
    _np = _typing.Any
from openscad import *  # noqa: F401,F403
from openscad import (  # noqa: F401
    Color,
    PyLibFive,
    PyOpenSCAD,
    PyOpenSCADs,
)

HAS_NUMPY: bool

class Customizer(_typing.Mapping[str, _typing.Any]):
    """Read-only mapping of Customizer parameter names to effective values."""

    def __init__(self) -> None: ...
    def __getitem__(self, name: str) -> _typing.Any: ...
    def __iter__(self) -> _typing.Iterator[str]: ...
    def __len__(self) -> int: ...
    def add_parameter(
        self,
        name: str,
        default: _typing.Any,
        description: str | None = ...,
        range: _typing.Any = ...,
        step: float | None = ...,
        max_length: int | None = ...,
        options: list[_typing.Any] | dict[_typing.Any, str] | None = ...,
    ) -> _typing.Any: ...
    def group(self, title: str) -> _CustomizerGroup: ...

class _CustomizerGroup(_typing.Mapping[str, _typing.Any]):
    def __getitem__(self, name: str) -> _typing.Any: ...
    def __iter__(self) -> _typing.Iterator[str]: ...
    def __len__(self) -> int: ...
    def add_parameter(
        self,
        name: str,
        default: _typing.Any,
        description: str | None = ...,
        range: _typing.Any = ...,
        step: float | None = ...,
        max_length: int | None = ...,
        options: list[_typing.Any] | dict[_typing.Any, str] | None = ...,
    ) -> _typing.Any: ...

class _VectorBase(_np.ndarray[_typing.Any, _np.dtype[_np.float64]]):
    """Base class for NumPy-backed fixed-length PythonSCAD vectors."""

    def __init__(
        self, iterable: _typing.Iterable[float] | None = ...
    ) -> None: ...
    def __array__(
        self, dtype: _typing.Any = ..., copy: _typing.Any = ...
    ) -> _typing.Any: ...
    @classmethod
    def from_array(cls, array: _typing.Any) -> _typing.Self: ...

class Vector1(_VectorBase):
    """1D vector represented as [x]."""

class Vector2(_VectorBase):
    """2D vector represented as [x, y]."""

class Vector3(_VectorBase):
    """3D vector represented as [x, y, z]."""

class Matrix4x4(_np.ndarray[_typing.Any, _np.dtype[_np.float64]]):
    """NumPy-backed 4x4 transformation matrix helper."""

    def __init__(
        self,
        iterable: _typing.Iterable[_typing.Iterable[float]] | None = ...,
    ) -> None: ...
    def __array__(
        self, dtype: _typing.Any = ..., copy: _typing.Any = ...
    ) -> _typing.Any: ...
    @classmethod
    def from_array(cls, array: _typing.Any) -> "Matrix4x4": ...

_MultiToolExporterItem = tuple[str, _typing.Any] | tuple[str, _typing.Any, bool]


class MultiToolExporter(list[_MultiToolExporterItem]):
    """List-based helper for exporting multi-tool / multi-color 3D models.

    Each item is a ``(name, object)`` 2-tuple or ``(name, object, export)``
    3-tuple (matching :func:`dict.items` and the multi-object form of
    :func:`export`). The optional ``export`` flag defaults to ``True``;
    ``export=False`` entries are *cutters* that still subtract from earlier
    items but are omitted from :meth:`parts`, :meth:`show`, and all export
    paths. For each index ``i``, :meth:`export` writes the geometry obtained
    by subtracting every later item's object from ``self[i]``'s object into
    either per-part files named ``f"{prefix}{name}{suffix}"`` or one
    multi-object 3MF file when ``single_file`` is provided. The last list
    entry is emitted as-is (no degenerate one-child ``difference`` node).
    Output paths and single-file part names must be unique among exportable
    items; collisions raise :class:`ValueError` at export time.
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
        items: _typing.Iterable[_MultiToolExporterItem] = ...,
    ) -> None:
        """Initialize the exporter, optionally seeding it with ``items``."""
        ...

    def append(self, item: _MultiToolExporterItem) -> None:
        """Append a validated ``(name, object)`` or ``(name, object, export)`` tuple."""
        ...

    def extend(self, items: _typing.Iterable[_MultiToolExporterItem]) -> None:
        """Append each validated item tuple from ``items``."""
        ...

    def insert(
        self, index: _typing.SupportsIndex, item: _MultiToolExporterItem
    ) -> None:
        """Insert a validated item tuple at ``index``."""
        ...

    def __iadd__(  # type: ignore[override]
        self,
        other: _typing.Iterable[_MultiToolExporterItem],
    ) -> "MultiToolExporter":
        """Validate each item then in-place extend (``self += other``)."""
        ...

    def parts(self) -> list[tuple[str, _typing.Any]]:
        """Return computed ``(name, geometry)`` pairs for exportable items."""
        ...

    def export(self, single_file: str | None = ...) -> None:
        """Export each part separately, or all parts into one 3MF file."""
        ...

    def show(self) -> None:
        """Display each part in the PythonSCAD preview."""
        ...

@_typing.overload
def rounded_cube(
    size: float | Vector3,
    r: float,
    *,
    center: bool | None = ...,
    fn: float | None = ...,
    fa: float | None = ...,
    fs: float | None = ...,
) -> PyOpenSCAD: ...
@_typing.overload
def rounded_cube(
    size: float | Vector3,
    *,
    d: float,
    center: bool | None = ...,
    fn: float | None = ...,
    fa: float | None = ...,
    fs: float | None = ...,
) -> PyOpenSCAD: ...
def rounded_cube(
    size: float | Vector3,
    r: float | None = ...,
    *,
    d: float | None = ...,
    center: bool | None = ...,
    fn: float | None = ...,
    fa: float | None = ...,
    fs: float | None = ...,
) -> PyOpenSCAD:
    """Create a cube or box with uniformly rounded edges and corners.

    Specify exactly one of ``r`` (radius) or ``d`` (diameter). Set
    ``center=True`` to center the generated shape's bounding box on the origin;
    ``False`` or ``None`` leaves it in the positive octant. Optional ``fn``,
    ``fa``, and ``fs`` control rounding-sphere tessellation.
    """
    ...

def loft(
    shape1: PyOpenSCAD,
    shape2: PyOpenSCAD,
    height: float,
    n: int = ...,
    rot: float = ...,
) -> _typing.Callable[[float], list[list[float]]]:
    """Interpolate a cross-section between two 2D shapes.

    Returns a function ``f(h)`` that, for a height ``h`` between ``0`` and
    ``height``, linearly interpolates by radius between the sampled outlines
    of ``shape1`` (at height ``0``) and ``shape2`` (at height ``height``),
    optionally twisted by ``rot`` degrees over that span.
    """
    ...
