"""Startup file that adds a menu entry and an editor form for `badge`.

PythonSCAD loads ~/.config/PythonSCAD/.pythonscadrc, so this file has to be
copied or symlinked there. It is named pythonscadrc.py in the repository so
it stays readable; the installed copy must be called .pythonscadrc.

Restart PythonSCAD after installing it. The file is read once, at startup.
"""

from PyQt6 import QtWidgets, sip
from pythonscad import (
    add_menuitem,
    cylinder,
    editor_get_call_args,
    editor_replace_call_args,
    mainwindow_ptr,
    show,
    text,
)
from pyscadforms import parse_arguments

_DEFAULT_LABEL = "PythonSCAD"
_DEFAULT_RADIUS = 30.0


def pyqt6_works():
    """Menu callback. Registered from setup() below.

    Opens a dialog parented to the running main window, which is the smallest
    possible proof that the startup file loaded and that PyQt6 can reach the
    application.
    """
    mw = sip.wrapinstance(mainwindow_ptr(), QtWidgets.QMainWindow)
    QtWidgets.QMessageBox.information(
        mw,
        "PyQt6 works",
        "This dialog was opened by .pythonscadrc through PyQt6.\n\n"
        "The startup file loaded, and Python can reach the running window.",
    )


class badge:
    """A disc with a raised label, editable from the editor calltip.

    __new__ returns real geometry, so the script renders normally whether or
    not the dialog was ever opened.
    """

    @staticmethod
    def get_calltip():
        return f'badge(label="{_DEFAULT_LABEL}", radius={_DEFAULT_RADIUS})'

    @staticmethod
    def on_editor_trigger(pos):
        mw = sip.wrapinstance(mainwindow_ptr(), QtWidgets.QMainWindow)
        known, unknown = parse_arguments(badge, editor_get_call_args(pos))

        dialog = QtWidgets.QDialog(mw)
        dialog.setWindowTitle("badge")
        layout = QtWidgets.QFormLayout(dialog)

        label = QtWidgets.QLineEdit(str(known.get("label", _DEFAULT_LABEL)))
        layout.addRow("label", label)

        radius = QtWidgets.QDoubleSpinBox()
        radius.setRange(5.0, 200.0)
        radius.setSingleStep(1.0)
        radius.setValue(float(known.get("radius", _DEFAULT_RADIUS)))
        layout.addRow("radius", radius)

        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        layout.addRow(buttons)

        if dialog.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            parts = [f"label={label.text()!r}", f"radius={radius.value()}"]
            parts += [f"{name}={value!r}" for name, value in unknown.items()]
            editor_replace_call_args(pos, ", ".join(parts))

    def __new__(cls, label=_DEFAULT_LABEL, radius=_DEFAULT_RADIUS):
        disc = cylinder(r=radius, h=2, center=True)
        glyph = text(label, size=radius / 3, halign="center", valign="center")
        show(disc + glyph.linear_extrude(3).translate([0, 0, -1]))


def setup():
    add_menuitem("Tools", "PyQt6 works?", "pyqt6_works()")
