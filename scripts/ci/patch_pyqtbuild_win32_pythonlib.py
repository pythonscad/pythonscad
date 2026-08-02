#!/usr/bin/env python3
"""Patch PyQt-builder to link Python explicitly when building with MinGW.

MinGW/g++ ignores the ``#pragma comment(lib, ...)`` that Python.h embeds on
Windows to link pythonXY.lib automatically; that pragma is an MSVC-only
feature. PyQt-builder adds only the library search path (``-L...``) on win32,
never an explicit ``-lpythonX.Y``, because it assumes the pragma is honored.
This script adds the missing linker flag directly to the installed
PyQt-builder source.
"""
import os
import pyqtbuild

path = os.path.join(os.path.dirname(pyqtbuild.__file__), "builder.py")
text = open(path, encoding="utf-8").read()

old = """            pro_lines.extend(['win32 {',
                    '    LIBS += -L{}'.format(
                            self.qmake_quote(project.py_pylib_dir)),
                    '}'])"""

new = """            pro_lines.extend(['win32 {',
                    '    LIBS += -L{}'.format(
                            self.qmake_quote(project.py_pylib_dir)),
                    '    LIBS += -lpython{}.{}'.format(
                            project.py_major_version, project.py_minor_version),
                    '}'])"""

if old not in text:
    raise SystemExit(
        "Expected code block not found in pyqtbuild/builder.py; the "
        "PyQt-builder version may have changed and the patch must be updated."
    )

text = text.replace(old, new)
open(path, "w", encoding="utf-8").write(text)
print(f"Patched: {path}")
