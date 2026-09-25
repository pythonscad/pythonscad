# I/O and Integration

## Choosing an import function

PythonSCAD has several ways to pull in external content. They are **not**
interchangeable — pick by what you are loading:

| Goal | Function | Returns |
|------|----------|---------|
| Local mesh / 2D drawing (STL, 3MF, SVG, …) | [`osimport`](#osimport) | Geometry object (or color→object dict for SVG `split_by_color`) |
| Local OpenSCAD library (`.scad` modules/functions/vars) | [`osuse`](#osuse) (prefer over deprecated [`osinclude`](#osinclude)) | Handle with attributes for modules, functions, and variables |
| Remote **Python** library over HTTP(S) | [`nimport`](#nimport) (GUI only) | `None` — side effect is `from <module> import *` into the current namespace |
| Local Python / PythonSCAD script | Ordinary Python `import` / `from … import …` | Whatever that module exports |

Inline OpenSCAD snippets (not a file) use [`scad`](#scad).

---

## osimport

Import geometry from a file. This is the PythonSCAD equivalent of OpenSCAD's `import()` (renamed because `import` is a Python keyword).

**Syntax:**

=== "Python"

    ```python
    osimport(file, layer=None, convexity=2, origin=None, scale=1,
             width=1, height=1, center=False, dpi=72, id=None, stroke=False,
             fn=None, fa=None, fs=None, split_by_color=False)
    ```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `file` | string | — | Path to the file to import |
| `layer` | string | `None` | Layer name (for DXF files) |
| `convexity` | int | `2` | Convexity hint |
| `origin` | `[x, y]` | `None` | Origin offset (for DXF) |
| `scale` | float | `1` | Scale factor (for DXF) |
| `width` | float | `1` | Width (for image-based imports) |
| `height` | float | `1` | Height (for image-based imports) |
| `center` | bool | `False` | Center the imported geometry |
| `dpi` | float | `72` | DPI for SVG imports |
| `id` | string | `None` | Element ID (for SVG) |
| `stroke` | bool | `False` | Include stroke paths (for SVG) |
| `fn`, `fa`, `fs` | float | global | Curve discretization; defaults to the global `fn`/`fa`/`fs` values |
| `split_by_color` | bool | `False` | SVG only: return `{hex_color: 2D object}` instead of one merged object; raises `ValueError` for non-SVG input or when no colors match the SVG selection. |

**Supported formats** (by file extension):

- **3D:** STL, OFF, OBJ, 3MF, STEP / STP; NEF3 when the build includes CGAL
- **2D:** DXF, SVG; CDR when the build includes optional CDR support (`ENABLE_CDR`)

`osimport` does **not** load `.scad` or `.py` files — use [`osuse`](#osuse) or
Python's `import` for those.

**Examples:**

=== "Python"

    ```python
    from pythonscad import *

    model = osimport("model.stl")
    model.show()

    drawing = osimport("design.svg", dpi=96)
    drawing.linear_extrude(height=2).show()
    ```

**Multi-color SVG to multi-object 3MF (e.g. for multi-toolhead printing):**

`split_by_color=True` returns one 2D object per SVG color, keyed by hex color
string. Extrude each separately and pass the resulting dict to `export()` to
produce a 3MF file with one separately named/colored object per SVG color —
ready to assign to different tool heads in a slicer:

=== "Python"

    ```python
    from pythonscad import *

    parts = osimport("watermelon.svg", split_by_color=True)
    # {"#4a7d2eff": <2D object>, "#c0392bff": <2D object>, ...}

    solids = {name: obj.linear_extrude(2) for name, obj in parts.items()}
    export(solids, "watermelon.3mf")
    ```

**OpenSCAD reference:** [import](https://en.wikibooks.org/wiki/OpenSCAD_User_Manual/Importing_Geometry#import)

---

## osuse

Load an OpenSCAD library file and **return a handle object** whose attributes
are the library's modules, functions, and top-level variables. PythonSCAD's
analog of OpenSCAD's `use <file.scad>`, with two semantic differences:

- The imported symbols are *not* injected into the global namespace — you
  must access them through the returned handle.
- The handle also exposes top-level variable assignments, which OpenSCAD's
  `use` does not import.

**Syntax:**

=== "Python"

    ```python
    lib = osuse(file)
    ```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `file` | string | Path to the `.scad` file |

**Returns:** an object exposing the imported library's contents:

- **Modules** become methods that return geometry objects (e.g. `lib.gear(...)`).
- **Functions** become methods that return values (e.g. `lib.calc_pitch(...)`).
- **Top-level variables** are accessible via attribute access for plain names
  (e.g. `lib.my_var`) and via item access for `$`-prefixed special variables
  (e.g. `lib["$fn"]`).

!!! note
    Calling a module from the imported library produces a geometry object —
    you still need to call `.show()` on it (or assign it to a variable and
    call `.show()` on that variable later) for it to appear in the output,
    just like any other PythonSCAD geometry.

**Examples:**

Calling a module from an imported library:

=== "Python"

    ```python
    from pythonscad import *

    mcad = osuse("MCAD/gears.scad")

    g = mcad.gear(
        number_of_teeth=20,
        circular_pitch=200,
        pressure_angle=20,
        clearance=0,
        verbose=False,
    )
    g.show()
    ```

Calling a function and reading variables from an imported library:

=== "Python"

    ```python
    from pythonscad import *

    lib = osuse("mylib.scad")

    width = lib.get_width()
    radius = lib.calculate_radius(diameter=10)

    print(lib.my_constant)
    print(lib["$fn"])
    ```

---

## osinclude

!!! warning "Deprecated"
    `osinclude` is deprecated — use [`osuse`](#osuse) for new code. Calling
    `osinclude` logs a deprecation message to the PythonSCAD console (it does
    not raise a Python `DeprecationWarning`).

PythonSCAD's analog of OpenSCAD's `include <file.scad>`. The returned handle
exposes the imported modules, functions, and top-level variables exactly the
same way as [`osuse`](#osuse) — see that section for usage examples.

`osuse` and `osinclude` differ in whether they evaluate the imported
file's top-level *module instantiations* — in OpenSCAD that category covers
calls like `cube()`, `echo()`, and `assert()`, which are syntactically
module calls. `osinclude` evaluates them, `osuse` suppresses them. Output
or errors from those top-level calls therefore only surface with
`osinclude`.

Top-level *variable assignments* (e.g. `x = 10;`) are always evaluated by
both functions — they are needed to populate the returned handle, so any
errors in their expressions will propagate from either call.

In both cases, any geometry produced by top-level module instantiations is
discarded — neither `osuse` nor `osinclude` exposes top-level geometry on
the returned handle. To use geometry from an OpenSCAD file, call its
modules explicitly via the handle (e.g. `lib.my_module().show()`).

**Syntax:**

=== "Python"

    ```python
    lib = osinclude(file)
    ```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `file` | string | Path to the `.scad` file |

**Examples:**

=== "Python"

    ```python
    from pythonscad import *

    lib = osinclude("config.scad")
    print(lib.default_radius)
    ```

---

## scad

Execute inline OpenSCAD code from within a Python script.

**Syntax:**

=== "Python"

    ```python
    scad(code)
    ```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `code` | string | OpenSCAD source code to execute |

**Examples:**

=== "Python"

    ```python
    from pythonscad import *

    result = scad("cube(10);")
    result.show()
    ```

---

## nimport

Download a **Python** module from a network URL into the user library directory,
then run `from <module_stem> import *` so its symbols appear in the current
namespace. This is **not** a geometry importer — it does not load STL, 3MF,
SVG, or other mesh/drawing files (use [`osimport`](#osimport) for those, after
downloading the file yourself if needed).

`nimport` is only available in **GUI** builds (it is omitted from headless /
`OPENSCAD_NOGUI` builds). It returns `None`; useful results come from the
imported module's exported names (for example a function or solid you then
call or `.show()`).

**Behavior:**

1. Take the last path segment of `url` as the filename (e.g. `mylib.py`).
2. Download it to the PythonSCAD user library path (skipped if this session
   already downloaded the same URL and the file is still present).
3. Execute `from <stem> import *` where `<stem>` is the filename without its
   final extension.

Preferences → Python can list default network-import URLs; new editor tabs
pre-fill matching `nimport("…")` lines from that list.

**Syntax:**

=== "Python"

    ```python
    nimport(url)
    ```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `url` | string | HTTP(S) URL of a `.py` module to download and star-import |

**Examples:**

=== "Python"

    ```python
    from pythonscad import *

    # Remote library that defines e.g. make_widget() / WIDGET_SIZE
    nimport("https://example.com/mylib.py")
    make_widget().show()
    ```

For a **local** PythonSCAD library on disk, prefer a normal Python import
(with the file on `sys.path` or next to your script) instead of `nimport`.
