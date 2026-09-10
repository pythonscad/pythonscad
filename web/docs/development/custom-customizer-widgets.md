# Custom PyQt6 Customizer Widgets

PythonSCAD ships PyQt6 bundled with the application. Beyond letting user
modules open their own PyQt6 dialogs (e.g. an interactive polygon-point
editor), you can register a custom PyQt6 widget as a first-class
Customizer parameter type. Once registered, `add_parameter()` can use it
just like the built-in slider/checkbox/text types: it shows up in the
Customizer sidebar, participates in parameter sets, and drives a live
preview re-render when changed.

## Quick example

```python
from PyQt6 import QtWidgets, QtCore

class ColorPicker(QtWidgets.QPushButton):
    valueChanged = QtCore.pyqtSignal()

    def __init__(self, name, default):
        super().__init__(default)
        self._value = default
        self.clicked.connect(self._pick)

    def _pick(self):
        set_modal_dialog_active(True)
        try:
            c = QtWidgets.QColorDialog.getColor()
        finally:
            set_modal_dialog_active(False)
        if c.isValid():
            self._value = c.name()
            self.setText(self._value)
            self.valueChanged.emit()

    def get_value(self):
        return self._value

    def set_value(self, value):
        self._value = value
        self.setText(value)


def color_picker_factory(name, default):
    return ColorPicker(name, default)


def setup():
    add_parameter_widget("color_picker", color_picker_factory)
```

```python
# In a model script:
linecolor = add_parameter("linecolor", "#ff00ff", type="color_picker")
cube(10).color(linecolor)
```

## Where registration belongs

`add_parameter_widget()` registers a widget *type*, not something tied to
a particular model. It belongs in your `.pythonscadrc`, inside a `setup()`
function, so it's available to every script you open — the same place
you'd register a menu item or an editor helper. Model scripts only call
`add_parameter(..., type="...")` to *use* an already-registered type.

## The factory function

```python
def factory(name, default):
    ...
    return widget
```

- `name`: the parameter name, as passed to `add_parameter()`.
- `default`: the default value, decoded from the value passed to
  `add_parameter()` (string, number, bool, or a list of those — whatever
  JSON can represent).
- Return value: any PyQt6 `QWidget` subclass.

## The widget contract

A widget registered this way must implement three things:

| Member | Direction | Purpose |
|---|---|---|
| `get_value(self)` | widget → parameter | Called when the Customizer needs to read the widget's current value (e.g. before saving a parameter set). |
| `set_value(self, value)` | parameter → widget | Called when the Customizer needs to push a stored value into the widget's display (e.g. after loading a parameter set, or after a script re-run rebuilds the widget). |
| `valueChanged = QtCore.pyqtSignal()` | widget → Customizer | Emit this whenever the user changes the value interactively, so the Customizer knows to save it and trigger a preview update. |

`get_value()`/`set_value()` exchange plain JSON-compatible values (str,
int, float, bool, list) — not PyQt6 types.

## Modal dialogs: `set_modal_dialog_active()`

If your widget's interaction opens a **modal** dialog (`QColorDialog`,
`QFileDialog`, a custom `QDialog.exec()`, or any blocking `getXxx()`
call), wrap that call:

```python
set_modal_dialog_active(True)
try:
    ...  # the blocking/modal call
finally:
    set_modal_dialog_active(False)
```

Background geometry evaluation acquires and releases the Python GIL from
worker threads. A modal Qt dialog runs a nested event loop and manages
its own GIL handoff internally; without this guard, a background render
that happens to be in flight while the dialog is open can deadlock.
Always pair the call with `try/finally` so the flag is cleared even if
the dialog call raises.

If your widget only reacts to non-modal events (e.g. a plain button that
never opens a blocking dialog), you don't need this at all.

## Why `valueChanged` should defer, not emit synchronously

A widget's Qt slot (e.g. a `clicked` handler) runs inside PyQt6's own
signal-dispatch machinery, which manages the GIL independently of
PythonSCAD's own lock/unlock bookkeeping. Emitting `valueChanged`
directly from inside that same call stack can re-enter PythonSCAD's
compile/render path (which itself needs the GIL) before PyQt6's dispatch
has cleanly returned. To avoid this, defer the signal to the next event
loop iteration instead of emitting it immediately:

```python
QtCore.QTimer.singleShot(0, self.valueChanged.emit)
```

Use this instead of emitting `valueChanged` directly at the end of a
slot that was itself invoked by a Qt signal (which is the common case:
a `clicked` handler, a `colorSelected` handler, etc.).

## Persistence behavior

Like built-in parameter types, a custom widget's value follows the
active parameter set:

- Editing the value while `<design default>` is selected creates a new
  set to hold the edit.
- Re-running the script (F5/F6) while `<design default>` is active
  resets every parameter — custom widgets included — to the value
  hard-coded in `add_parameter(...)`. This matches the existing
  behavior of built-in types; it is not specific to custom widgets.
- Saving and selecting a named set preserves the value across script
  re-runs.

## Limitations

- Values exchanged via `get_value()`/`set_value()` must be
  JSON-representable (str, number, bool, or a list of those).
- The factory is called once per Customizer rebuild (e.g. after F5/F6),
  not once per application lifetime — don't assume widget instance
  identity persists across script re-runs.
