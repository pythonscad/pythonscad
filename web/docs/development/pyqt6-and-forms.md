# PyQt6 and Forms in PythonSCAD

The Qt6 AppImage, the macOS app and the Windows installer can
`import PyQt6` without a separate install. Debian and Fedora packages use
the distro's `python3-pyqt6` instead, which the package manager installs
with PythonSCAD. This page covers where it comes from, how to reach the
running application window from Python, and how to build a basic dialog.

The Ubuntu 22.04 package is Qt5 and has no PyQt. Rocky Linux 9's
repositories do not ship `python3-pyqt6`, so that package has no PyQt
either. Dialogs, forms and custom widgets are unavailable on those two
packages. The rest of PythonSCAD continues to work normally.

## Where PyQt6 Comes From

Portable Qt6 builds carry a PyQt6 that was built against the same Qt the
application links. On Debian and Fedora it is the distro package. It is
available as soon as PythonSCAD starts.

```python
from PyQt6 import QtWidgets
```

This works out of the box in any script, and in `.pythonscadrc`.

### Build-system lifecycle

Debian and Fedora declare a dependency on `python3-pyqt6`. The AppImage
copies the build runner's `python3-pyqt6` into the image. macOS links the
Qt frameworks inside the PyPI PyQt6 wheel and ships that wheel's bindings.
Windows compiles PyQt6 in the release workflow and caches the result. When
Windows support moves to a new Qt minor, add its reviewed PyQt version to
`scripts/ci/pyqt6-versions.json` first.

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

- PyQt6 is prebuilt and shipped with the Qt6 AppImage, macOS, and Windows
  packages — `from PyQt6 import ...` works immediately there. Qt5
  installations do not provide this feature.
- Don't create a second `QApplication`; wrap the existing one via
  `mainwindow_ptr()` / `qapp_ptr()` and `sip.wrapinstance`.
- From there, it's ordinary PyQt6 — dialogs, layouts, widgets, custom
  painting.
