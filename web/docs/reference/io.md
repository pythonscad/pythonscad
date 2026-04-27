# I/O and Integration

## osimport

Import geometry from a file. This is the PythonSCAD equivalent of OpenSCAD's `import()` (renamed because `import` is a Python keyword).

**Syntax:**

=== "Python"

```python
osimport(file, layer="", convexity=2, origin=None, scale=None,
         width=0, height=0, center=False, dpi=96, id="", stroke=False,
         fn=0, fa=0, fs=0)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `file` | string | — | Path to the file to import |
| `layer` | string | `""` | Layer name (for DXF files) |
| `convexity` | int | `2` | Convexity hint |
| `origin` | `[x, y]` | `None` | Origin offset (for DXF) |
| `scale` | float | `None` | Scale factor (for DXF) |
| `width` | float | `0` | Width (for image-based imports) |
| `height` | float | `0` | Height (for image-based imports) |
| `center` | bool | `False` | Center the imported geometry |
| `dpi` | float | `96` | DPI for SVG imports |
| `id` | string | `""` | Element ID (for SVG) |
| `stroke` | bool | `False` | Include stroke paths (for SVG) |
| `fn`, `fa`, `fs` | float | global | Curve discretization; defaults to the global `fn`/`fa`/`fs` values |

**Supported formats:** STL, OFF, AMF, 3MF (3D); DXF, SVG (2D)

**Examples:**

=== "Python"

```python
from openscad import *

model = osimport("model.stl")
model.show()

drawing = osimport("design.svg", dpi=96)
drawing.linear_extrude(height=2).show()
```

**OpenSCAD reference:** [import](https://en.wikibooks.org/wiki/OpenSCAD_User_Manual/Importing_Geometry#import)

---

## osuse

Load an OpenSCAD library file and **return a handle object** whose attributes
are the library's modules, functions, and top-level variables. Equivalent to
OpenSCAD's `use <file.scad>`, but unlike OpenSCAD's `use`, the imported
symbols are *not* injected into the global namespace — you must access them
through the returned object.

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
    you still need to call `.show()` on it (or assign it to a variable that
    you later `show()`) for it to appear in the output, just like any other
    PythonSCAD geometry.

**Examples:**

Calling a module from an imported library:

=== "Python"

```python
from openscad import *

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
from openscad import *

lib = osuse("mylib.scad")

width = lib.get_width()
radius = lib.calculate_radius(diameter=10)

print(lib.my_constant)
print(lib["$fn"])
```

---

## osinclude

!!! warning "Deprecated"
    `osinclude` is deprecated — use [`osuse`](#osuse) instead. Both functions
    return the same handle object; `osinclude` exists only for backwards
    compatibility with early PythonSCAD scripts and emits a deprecation
    warning when called.

Include an OpenSCAD file, executing its top-level code. Equivalent to
OpenSCAD's `include <file.scad>`. Returns the same kind of handle object as
[`osuse`](#osuse); see that section for how to access the imported modules,
functions, and variables.

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
from openscad import *

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
from openscad import *

result = scad("cube(10);")
result.show()
```

---

## nimport

Import a model from a network URL. This function is only available in GUI mode.

**Syntax:**

=== "Python"

```python
nimport(url)
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `url` | string | URL of the model to download and import |

**Examples:**

=== "Python"

```python
from openscad import *

model = nimport("https://example.com/model.stl")
model.show()
```
