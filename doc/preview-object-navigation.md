# Preview object navigation

After a preview (**F5**) you can step through the objects of the design with
the keyboard, one source line at a time. The line is selected in the editor
and every object it created is highlighted in the 3D view (like the `#`
modifier). This helps to find which code describes a given part of the design.
The source text is never changed, only the editor selection.

A line that creates many objects, e.g. a `cube()` call in a loop that runs a
hundred times, is a single step that highlights all hundred cubes. The status
bar shows the number of objects, e.g. `Preview line 2 of 4: 100 objects (cube,
...) (design.py:3)`.

| Shortcut                          | Action (Design menu)                                       |
| --------------------------------- | ---------------------------------------------------------- |
| **Alt+Down** (macOS: Option+Down) | Next Preview Object; right after a preview: the first line  |
| **Alt+Up** (macOS: Option+Up)     | Previous Preview Object; right after a preview: the last line |
| (menu only)                       | Clear Preview Object Selection                             |

Navigation stops at the first and last line. The order is the order of the
object tree: an operation comes before its operands, operands in order; a line
takes the place of the first object it created. Pure
modifiers (`translate`, `rotate`, `color`, ...) are not separate objects; the
object they modify is selected instead. A `render()` is one object.

```python
from openscad import *
a = cube(10)
b = sphere(6).translate(20, 0, 0)
c = cylinder(r=3, h=15).translate(-20, 0, 0)
show(a | b | c)
```

**Alt+Down** after F5 selects line 5 (`a | b | c` is `(a | b) | c`, two unions,
both highlighted as one step), then the cube (line 2), the sphere (line 3) and
the cylinder (line 4). **Alt+Up** after F5 starts with the cylinder (line 4).

While navigating, the `#` objects of the design itself stay highlighted, but
they are kept apart from the selection: selecting the spheres of a body that
also contains a `#` part highlights only the spheres. A selected first
operand of a `difference` is shown with the cut applied, like a `#` there.

Objects created in another file, e.g. in a function of an imported `.py`
module, select their line in that file; the file is opened in a tab if needed.
Objects created by helpers of the `openscad`/`pythonscad` packages (such as
`rounded_cube()`) point to the line of your script that called the helper.

The navigation is reset by a new preview or render, by editing the document,
and by closing the previewed tab. Unsaved (untitled) Python documents have no
file name, so their objects cannot be navigated; save the file first.

On Linux and Windows **Alt+Up/Down** in the editor also changes the number
under the cursor. When the cursor sits inside a number, that wins; otherwise
the keys navigate preview objects. (On macOS that number feature uses
Option+Shift+Up/Down, so there is no overlap.)

## Implementation notes

- Every node created from Python stores the location (absolute file, line) of
  the innermost user Python frame in `node->modinst->location()`. Python
  registers this via `genlang_location_provider` (`src/genlang/genlang.h`);
  `DECLARE_INSTANCE()` uses it.
- `evaluatePython(code, dry_run, filename)` compiles the code under its real
  file name, so frames and tracebacks name the user's file.
- `PreviewObjectNavigator` (`src/core/PreviewObjectNavigator.h`) groups the
  objects by (file, line); `CSGTreeEvaluator`'s `previewHighlightIndices`
  produces the temporary highlight of a group without touching the AST. Only
  the outermost members of a group are passed, so nested objects of the same
  line are not drawn twice.
