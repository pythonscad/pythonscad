# PyQt6 and Forms in PythonSCAD

PythonSCAD ships with a prebuilt copy of PyQt6, ready to import from any
script or from `.pythonscadrc` — no separate installation or build step is
needed. This page covers where it comes from, how to reach the running
application window from Python, and how to build a basic dialog.

## Where PyQt6 Comes From

PyQt6 is built for each supported platform against the exact Qt6 version
PythonSCAD itself links against, and shipped inside every installation under
`libraries/python/PyQt6`. It is available the moment PythonSCAD starts —
there is nothing to install, and no internet access is required at runtime.

```python
from PyQt6 import QtWidgets
```

This works out of the box in any script, and in `.pythonscadrc`.

## One Application, Not Two

PythonSCAD's C++ side already runs a `QApplication` and a main window. Qt
does not allow a second `QApplication` in the same process, so Python code
does not create one — it wraps the existing objects instead, using `sip`
(the same binding technology PyQt6 itself is built with).

Two functions are available for this:

```python
mainwindow_ptr()   # raw pointer to the running MainWindow
qapp_ptr()          # raw pointer to the running QApplication
```

Wrap them with `sip.wrapinstance` to get real, usable PyQt6 objects backed by
the actual running application:

```python
from PyQt6 import sip, QtWidgets

app = sip.wrapinstance(qapp_ptr(), QtWidgets.QApplication)
mw = sip.wrapinstance(mainwindow_ptr(), QtWidgets.QMainWindow)
```

`mw` behaves like any other `QMainWindow` from here on — you can read its
title, use it as a parent for your own dialogs, and so on.

## A Minimal Dialog

```python
from PyQt6 import sip, QtWidgets

mw = sip.wrapinstance(mainwindow_ptr(), QtWidgets.QMainWindow)

dialog = QtWidgets.QDialog(mw)
dialog.setWindowTitle("Example")
layout = QtWidgets.QVBoxLayout(dialog)

label = QtWidgets.QLabel("Hello from PyQt6")
layout.addWidget(label)

buttons = QtWidgets.QDialogButtonBox(
    QtWidgets.QDialogButtonBox.StandardButton.Ok
    | QtWidgets.QDialogButtonBox.StandardButton.Cancel
)
buttons.accepted.connect(dialog.accept)
buttons.rejected.connect(dialog.reject)
layout.addWidget(buttons)

if dialog.exec() == QtWidgets.QDialog.DialogCode.Accepted:
    print("Accepted")
```

You can trigger this from a menu item (see
[.pythonscadrc and add_menuitem](pythonscadrc-and-menuitems.md)), or bind it
to a class so it opens automatically from the editor — covered next in
[Extend your design with an Interactive Form](extend-with-interactive-form.md).

## What's Available

Since the bundled PyQt6 is a normal, complete build, the full range of Qt
Widgets is available: form layouts, spin boxes, tables, custom `QWidget`
subclasses with their own `paintEvent` and mouse handling, and so on.
Anything you can do in a standalone PyQt6 application, you can do here,
parented to PythonSCAD's own window.

## Summary

- PyQt6 is prebuilt and shipped with every installation — `from PyQt6 import
  ...` works immediately.
- Don't create a second `QApplication`; wrap the existing one via
  `mainwindow_ptr()` / `qapp_ptr()` and `sip.wrapinstance`.
- From there, it's ordinary PyQt6 — dialogs, layouts, widgets, custom
  painting.
