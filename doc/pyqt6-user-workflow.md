# Using PyQt6 from a PythonSCAD script

The user-facing documentation already exists, in the development section of
the manual:

- [PyQt6 and Forms in PythonSCAD](../web/docs/development/pyqt6-and-forms.md)
  covers where the bundled bindings come from, why you wrap the existing
  application instead of creating one, and a minimal dialog.
- [Extend Your Design with an Interactive Form](../web/docs/development/extend-with-interactive-form.md)
  covers the `get_calltip` / `on_editor_trigger` / `__new__` protocol, the
  editor helper functions, and wrapping a builtin.
- [Custom PyQt6 Customizer Widgets](../web/docs/development/custom-customizer-widgets.md)
  covers registering a widget type for the Customizer sidebar.
- [.pythonscadrc and add_menuitem](../web/docs/development/pythonscadrc-and-menuitems.md)
  covers the startup file and menu entries.

`web/mkdocs.yml` links the first, third and fourth of those. The interactive
form page is reached only through links inside the other pages.

What follows is the same material condensed into the sequence you actually
perform, plus the parts the manual leaves implicit.

## What you need

A build that can import `PyQt6`. The Qt6 AppImage, the macOS app and the
Windows build ship it, so `from PyQt6 import QtWidgets` works in any script
with nothing installed. A Qt5 build does not, and the features below are
absent there.

Debian and Fedora packages depend on the distro's `python3-pyqt6`. A
locally built binary uses that same package when it was linked against the
distro Qt. Ubuntu 22.04 and Rocky Linux 9 do not ship PyQt6, so those
packages have none. You can check a Qt6 build with `pythonscad --info`:
the Qt version it prints should match
`python3 -c "from PyQt6.QtCore import QT_VERSION_STR; print(QT_VERSION_STR)"`.

## Where code has to live

PythonSCAD reads `~/.config/PythonSCAD/.pythonscadrc` once, at startup, in the
same Python namespace your scripts use, and then calls `setup()` if it
exists. Anything defined while that file loads stays available for the whole
session and survives rendering.

That persistence is the reason editor forms and custom widgets belong in this
file. A class defined only inside a model script is visible after that script
renders, and vanishes again the next time a different document renders. Menu
entries belong specifically inside `setup()`, because the menus already exist
at that point.

Editing `.pythonscadrc` does nothing until PythonSCAD is restarted. The file is
only read by the GUI. Running `pythonscad file.py` on the command line never
loads it, so menu entries and editor forms are GUI-only by construction.

One trap: rendering a script that itself runs `from pythonscad import *`
re-imports the builtins into that shared namespace and silently undoes any
builtin you redefined, such as a `cube` class with a form. Keep that import
in `.pythonscadrc` and out of the model scripts.

## Opening a dialog

```python
from PyQt6 import sip, QtWidgets
from pythonscad import mainwindow_ptr

mw = sip.wrapinstance(mainwindow_ptr(), QtWidgets.QMainWindow)

dialog = QtWidgets.QDialog(mw)
dialog.setWindowTitle("Example")
layout = QtWidgets.QVBoxLayout(dialog)
layout.addWidget(QtWidgets.QLabel("Hello from PyQt6"))

buttons = QtWidgets.QDialogButtonBox(
    QtWidgets.QDialogButtonBox.StandardButton.Ok
    | QtWidgets.QDialogButtonBox.StandardButton.Cancel
)
buttons.accepted.connect(dialog.accept)
buttons.rejected.connect(dialog.reject)
layout.addWidget(buttons)

dialog.exec()
```

`mainwindow_ptr()` returns the address of the window PythonSCAD is already
running. `sip.wrapinstance` turns that address into a usable PyQt object, so
the dialog is parented to the real window. Creating a new `QApplication` here
fails, because one already exists.

## Adding a menu entry

```python
from pythonscad import add_menuitem

def hello():
    print("hello from the menu")

def setup():
    add_menuitem("Tools", "Say hello", "hello()")
```

The third argument is source text, executed later in the same namespace, not
a function object.

## An editor form

A class becomes editable from the calltip when it defines three things.

```python
from PyQt6 import QtWidgets, sip
from pythonscad import (
    editor_get_call_args,
    editor_replace_call_args,
    mainwindow_ptr,
    sphere,
)

class mysphere:
    @staticmethod
    def get_calltip():
        return "mysphere(radius=10.0)"

    @staticmethod
    def on_editor_trigger(pos):
        mw = sip.wrapinstance(mainwindow_ptr(), QtWidgets.QMainWindow)
        dialog = QtWidgets.QDialog(mw)
        dialog.setWindowTitle("mysphere")
        layout = QtWidgets.QFormLayout(dialog)

        radius = QtWidgets.QDoubleSpinBox()
        radius.setRange(0.1, 1000)
        radius.setValue(10.0)
        layout.addRow("radius", radius)

        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok
        )
        buttons.accepted.connect(dialog.accept)
        layout.addRow(buttons)

        if dialog.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            editor_replace_call_args(pos, f"radius={radius.value()}")

    def __new__(cls, radius=10.0):
        return sphere(radius)
```

`get_calltip()` supplies the text shown when you type `mysphere(`. With it
present, the calltip becomes clickable. `on_editor_trigger(pos)` runs on that
click; `pos` marks the call in the buffer, and `editor_get_call_args(pos)`
and `editor_replace_call_args(pos, text)` read and rewrite the argument text.
`__new__` returns the actual geometry, so the class behaves like a builtin
when the script runs.

`pyscadforms.parse_arguments(cls, text)` parses existing arguments against
`__new__`'s signature, which is what lets a form edit a call that is already
written instead of only creating one. `libraries/python/pyscadforms.py`
contains two full implementations, an interactive `polygon` canvas and a
`cube` size dialog, that show the pattern end to end.

## A Customizer widget

Register a type once, in `setup()`:

```python
from PyQt6 import QtCore, QtWidgets
from pythonscad import add_parameter_widget, set_modal_dialog_active

class ColorButton(QtWidgets.QPushButton):
    valueChanged = QtCore.pyqtSignal()

    def __init__(self, name, default):
        super().__init__(default)
        self._value = default
        self.clicked.connect(self._pick)

    def _pick(self):
        set_modal_dialog_active(True)
        try:
            color = QtWidgets.QColorDialog.getColor()
        finally:
            set_modal_dialog_active(False)
        if color.isValid():
            self._value = color.name()
            self.setText(self._value)
            self.valueChanged.emit()

    def get_value(self):
        return self._value

    def set_value(self, value):
        self._value = value
        self.setText(value)

def setup():
    add_parameter_widget("color_button", lambda name, default: ColorButton(name, default))
```

A model script then uses it like any built-in parameter type:

```python
linecolor = add_parameter("linecolor", "#ff00ff", type="color_button")
cube(10).color(linecolor)
```

The widget must offer `get_value()`, `set_value(value)` and a `valueChanged`
signal, and the values must be JSON-serialisable. `set_modal_dialog_active`
is only needed around a blocking dialog; it prevents a background render from
deadlocking against the dialog's event loop. The manual recommends emitting
`valueChanged` through `QTimer.singleShot(0, ...)` rather than directly from
inside another Qt slot.

## Trying it

`doc/examples/pyqt6-form/` contains a small `.pythonscadrc` and a model
script that exercise the menu entry and the editor form. See the `README.md`
next to them for the exact steps.
