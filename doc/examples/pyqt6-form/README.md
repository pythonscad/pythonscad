# Editor form demo

Shows the two things the PyQt6 integration adds for a user: a menu entry, and
a dialog that edits a function call directly in the editor.

## Install

PythonSCAD reads `.pythonscadrc` from `~/.config/PythonSCAD/` exactly once, at
startup.

```bash
mkdir -p ~/.config/PythonSCAD
ln -sf "$PWD/pythonscadrc.py" ~/.config/PythonSCAD/.pythonscadrc
```

Run that from this directory, then restart PythonSCAD. If you already have a
`.pythonscadrc`, append the contents of `pythonscadrc.py` to it instead of
replacing it.

This relies on the system PyQt6, which works for a locally built binary and
for the Debian package because both link against the distribution's Qt. Check
they agree before bothering with the rest:

```bash
pythonscad --info | grep "Qt version"
python3 -c "from PyQt6.QtCore import QT_VERSION_STR; print(QT_VERSION_STR)"
```

The two versions should be the same. The Qt6 AppImage, macOS and Windows
builds ship their own PyQt6 and need nothing installed.

## Try the menu entry

Open PythonSCAD. The Tools menu should contain "About the editor form demo".
Choosing it opens a message box parented to the main window. That exercises
`add_menuitem` and `sip.wrapinstance(mainwindow_ptr(), ...)`.

`.pythonscadrc` is loaded by the GUI only. `pythonscad file.py` on the command
line never reads it, so the menu entry and the calltip dialog exist only in
the application window. `badge.py` defines its own `badge` for exactly that
reason: the command line has no editor and no startup file, but the geometry
should still render.

```bash
pythonscad --trust-python -o badge.stl badge.py
```

## Try the editor form

Open `badge.py` from this directory in the GUI. The `badge` you get there is
the one from `.pythonscadrc`, which shadows the function defined in the
script.

1. Put the cursor inside the parentheses of `badge("PythonSCAD", 30)`.
2. The calltip shows `badge(label=..., radius=...)` with an arrow at the end.
3. Click the calltip. A dialog opens with the current label and radius filled
   in.
4. Change them and press Ok. The line in the editor is rewritten as keyword
   arguments, for example `badge(label='Hello', radius=40.0)`.
5. Press F5. You get a disc with the new label raised on top.

To see it create a call from scratch, delete the line, type `badge(` and
click the calltip. Ok writes a complete call where you were typing.

## What is happening

`pythonscadrc.py` defines `badge` with three hooks. `get_calltip()` provides
the calltip text. `on_editor_trigger(pos)` builds the dialog and writes the
result back through `editor_replace_call_args`. `__new__` returns ordinary
geometry, so the script renders the same way whether or not the dialog was
used. `pyscadforms.parse_arguments` is what lets the dialog show the values
already in the file instead of always starting from the defaults.

The background is in `doc/pyqt6-user-workflow.md`, and the build and
packaging side is in `doc/pyqt6-internals.md`.
