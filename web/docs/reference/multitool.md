# Multi-tool Export

PythonSCAD ships a small helper, `MultiToolExporter`, for the common
multi-tool / multi-color workflow where a single model is split into several
parts (typically one per filament color or print head) and each part is
exported to its own file.

The class lives in the `pythonscad` package (a strict superset of `openscad`),
so it is available under either of these imports:

=== "Python"

```python
from pythonscad import *
# or, equivalently for this class:
from pythonscad import MultiToolExporter
```

## MultiToolExporter

`MultiToolExporter` is a `list` subclass whose items are
`(object, name)` 2-tuples. The `name` is a label used to build the output
filename and must be a non-empty string and unique within the exporter.

**Filename layout:** for each item, the exporter writes to:

```text
{prefix}{name}{suffix}
```

so a typical use is `prefix="out/model-"` and `suffix=".stl"`.

**Cumulative-difference semantics:** for each index `i`, the geometry
exported is

```text
self[i][0] − self[i+1][0] − self[i+2][0] − ... − self[-1][0]
```

i.e. the first object's volume minus every later object's volume.
This guarantees that overlapping regions are claimed by exactly one part:
earlier entries "win", later parts only contain the volume not already
claimed. The last entry is exported as-is (no degenerate one-child
`difference` node).

**Constructor:**

=== "Python"

```python
MultiToolExporter(prefix, suffix, mkdir=False, items=())
```

**Parameters:**

| Parameter | Type                            | Default | Description                                                                            |
|-----------|---------------------------------|---------|----------------------------------------------------------------------------------------|
| `prefix`  | `str`                           | —       | Prepended to every output filename                                                     |
| `suffix`  | `str`                           | —       | Appended to every output filename (typically the extension, e.g. `".stl"`, `".3mf"`)   |
| `mkdir`   | `bool`                          | `False` | If `True`, create each output file's directory with `os.makedirs(..., exist_ok=True)`. See note below. |
| `items`   | iterable of `(object, name)`    | `()`    | Optional initial items, validated as if they were appended one at a time.              |

When `mkdir=True`, filenames without a directory component (e.g.
`"flag-"` rather than `"out/flag-"`) are exported as-is - no directory
is created and no error is raised.

**Methods:**

| Method            | Description                                                                                                              |
|-------------------|--------------------------------------------------------------------------------------------------------------------------|
| `append(item)`    | Append a single `(object, name)` 2-tuple. Validates the shape.                                                           |
| `extend(items)`   | Append each `(object, name)` from an iterable.                                                                           |
| `insert(i, item)` | Insert a single `(object, name)` 2-tuple at position `i`.                                                                |
| `export()`        | Write each part to its file. Raises `ValueError` if any two items share the same name.                                  |
| `show()`          | Render each part into the preview viewport (same cumulative-difference semantics as `export`).                           |

**Validation:**

* Every inserted item must be a 2-tuple of `(object, str)`. Anything else
  raises `TypeError`.
* The `name` must be a non-empty string. Empty names raise `ValueError`.
* At `export()` time, duplicate names raise `ValueError` rather than letting
  later parts silently overwrite the files written by earlier parts.

**Examples:**

A two-color flag (red star cut out of a blue background):

=== "Python"

```python
from pythonscad import *

background = cube([200, 100, 1]).color("blue")
star       = cylinder(r=20, h=2, fn=5).translate([100, 50, -0.5]).color("red")

exporter = MultiToolExporter("out/flag-", ".stl", mkdir=True)
exporter.append((star,       "red"))     # earlier wins
exporter.append((background, "blue"))    # blue gets background minus star
exporter.export()
# -> writes out/flag-red.stl and out/flag-blue.stl
```

Seeding from the constructor:

=== "Python"

```python
from pythonscad import *

red  = cube(10).color("red")
blue = cube(10).color("blue").right(5)

MultiToolExporter(
    prefix="out/cube-",
    suffix=".3mf",
    mkdir=True,
    items=[(red, "red"), (blue, "blue")],
).export()
```

Previewing the same split inside the GUI without writing files:

=== "Python"

```python
from pythonscad import *

exporter = MultiToolExporter("ignored-", ".stl")
exporter.append((red,  "red"))
exporter.append((blue, "blue"))
exporter.show()
```
