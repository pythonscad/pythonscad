# Extend Your Design with an Interactive Form

PythonSCAD's editor can show more than a plain calltip when you pause after
typing a function call like `cube(`. If a Python class defines a small set of
conventional methods, PythonSCAD offers a clickable indicator in the calltip
that opens a real PyQt6 dialog — a graphical form for setting parameters, a
canvas for placing points, or anything else a `QDialog` can show. Accepting
the dialog rewrites the call in the editor with the new arguments.

This is not a built-in feature of any single class — it is a small, generic
protocol. Any class that implements it gets the same editor integration, with
no changes to PythonSCAD itself. This page assumes you're already familiar
with [PyQt6 and Forms in PythonSCAD](pyqt6-and-forms.md) and
[.pythonscadrc](pythonscadrc-and-menuitems.md), since editor-trigger classes
normally live in one or the other.

Interactive forms require a Qt6 build of PythonSCAD. Qt5 builds do not bundle
PyQt and continue to use the standard, non-interactive editor calltips.

## Why This Exists

Some parameters are awkward to reason about as plain numbers — a 2D
polygon's point list is the clearest example. Rather than asking a user to
hand-type coordinates, a class can present a small drawing canvas instead,
and write the resulting point list back into the script for them. The same
mechanism works for anything with a fixed, named parameter set — including
primitives that already exist, like `cube`, once wrapped.

## The Protocol

A class becomes editor-aware by implementing three things:

```python
class mypolygon:
    @staticmethod
    def get_calltip():
        """Returns the text shown in the calltip popup."""
        return "mypolygon(points=[[x1,y1], [x2,y2], ...])"

    @staticmethod
    def on_editor_trigger(pos):
        """Called when the user clicks the calltip. `pos` is the character
        offset just after the opening parenthesis of the call in the editor."""
        ...

    def __new__(cls, points):
        """The actual geometry construction. Runs whenever a script calls
        mypolygon(...), whether or not the form was ever used."""
        return polygon([tuple(p) for p in points])
```

None of these three are required together — a class with only `__new__` is
just a normal callable, and PythonSCAD falls back to its regular built-in
calltip if `get_calltip`/`on_editor_trigger` are absent. The point-and-click
form is opt-in per class.

### `get_calltip()`

A `staticmethod` returning the string PythonSCAD shows in the popup when you
type the class name followed by `(`. Keep it short — it's a single-line hint,
the same as the built-in calltip for a native primitive.

### `on_editor_trigger(pos)`

A `staticmethod` called when the calltip is clicked. `pos` is the character
offset (into the current document) just after the call's opening parenthesis — pass
it straight through to the two helper functions below, which read and
replace the argument text between that parenthesis and its matching close.
This is where you build and run your `QDialog` (see
[PyQt6 and Forms in PythonSCAD](pyqt6-and-forms.md) for the basics of
getting a dialog on screen).

### `__new__(cls, ...)`

The class itself is the thing a script actually calls. Using `__new__`
instead of `__init__` lets the class return a completely different object —
typically the real geometry node — rather than an instance of the wrapper
class itself. Calling `mypolygon(points=[...])` in a script returns actual
solid/2D geometry, indistinguishable from a native primitive.

## Helper Functions

The editor helpers are exported by `pythonscad`. The generic argument parser
is a pure-Python helper in `pyscadforms`:

| Function | Purpose |
|---|---|
| `editor_get_call_args(pos)` | Returns the raw text between the parentheses of the call at `pos` — empty string if there is no matching close paren (a fresh, unfinished call). |
| `editor_replace_call_args(pos, new_args_text)` | Replaces the argument text of the call at `pos` with `new_args_text`, adding a closing parenthesis if none was found. |
| `pyscadforms.parse_arguments(cls, args_text)` | Parses an existing call's argument text against `cls.__new__`'s actual signature. Returns `(known, unknown)`. |

`pyscadforms.parse_arguments` is intentionally generic — it doesn't know anything about
what a `points` or `size` parameter *means*, only that `cls.__new__` accepts
one. This is what lets the same form both create a new call and edit an
existing one: on open, pre-fill from `parse_arguments`; on accept, write
everything back in keyword form via `editor_replace_call_args`.

## Worked Example

A minimal, complete class — no canvas, just two numeric fields, editable both
for new and existing calls:

```python
from PyQt6 import QtWidgets, sip
from pythonscad import *
from pyscadforms import parse_arguments


class mysphere:
    @staticmethod
    def get_calltip():
        return "mysphere(radius=10.0, twist=0.0)"

    @staticmethod
    def on_editor_trigger(pos):
        mw = sip.wrapinstance(mainwindow_ptr(), QtWidgets.QMainWindow)
        existing_args = editor_get_call_args(pos)
        known, unknown = parse_arguments(mysphere, existing_args)

        dialog = QtWidgets.QDialog(mw)
        dialog.setWindowTitle("mysphere")
        layout = QtWidgets.QFormLayout(dialog)

        radius_box = QtWidgets.QDoubleSpinBox()
        radius_box.setRange(0.1, 1000)
        radius_box.setValue(known.get("radius", 10.0))
        layout.addRow("radius", radius_box)

        twist_box = QtWidgets.QDoubleSpinBox()
        twist_box.setRange(-360, 360)
        twist_box.setValue(known.get("twist", 0.0))
        layout.addRow("twist", twist_box)

        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        layout.addRow(buttons)

        if dialog.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            parts = [f"radius={radius_box.value()}", f"twist={twist_box.value()}"]
            parts += [f"{k}={v!r}" for k, v in unknown.items()]
            editor_replace_call_args(pos, ", ".join(parts))

    def __new__(cls, radius=10.0, twist=0.0):
        solid = sphere(radius)
        if twist:
            solid = rotate([0, 0, twist])(solid)
        return solid
```

## Wrapping a Native Primitive

The same pattern works for primitives that already exist, such as `cube`.
Because a class definition simply rebinds the name in the current namespace,
you must capture the original callable **first**, or you lose access to it —
and, if `__new__` calls the name being defined instead of the saved original,
every call recurses into itself:

```python
from pythonscad import *
_native_cube = cube   # capture the original BEFORE the name is overwritten


class cube:
    @staticmethod
    def get_calltip():
        return "cube(size=[1,1,1], center=False)"

    def __new__(cls, size=1, center=False):
        return _native_cube(size, center)   # not "cube(...)" -- that would recurse

    @staticmethod
    def on_editor_trigger(pos):
        ...
```

Existing scripts that call `cube(10)` or `cube([1,2,3], center=True)` keep
working unmodified — `__new__` simply forwards everything to the saved
original.

## Where to Put This

Editor-trigger classes need to be defined somewhere PythonSCAD loads
automatically and keeps around for the life of the session — typically
`.pythonscadrc`, or a library imported from it. See
[.pythonscadrc and add_menuitem](pythonscadrc-and-menuitems.md) for why this
matters and what happens to definitions made elsewhere.

```python
# .pythonscadrc
import sys
sys.path.insert(0, "/path/to/your/library/directory")

from pythonscad import *
from your_forms_library import *
```

## Common Pitfalls

- **Forgetting `_native_x = x` before redefining a built-in.** Once the class
  statement runs, the original is gone unless you saved a reference to it
  first.
- **Calling the shadowed name instead of the saved original inside
  `__new__`.** This produces infinite recursion the first time a script
  actually calls the class.
- **Defining editor-trigger classes only inside a script, not
  `.pythonscadrc`.** Works only intermittently, depending on render order.
