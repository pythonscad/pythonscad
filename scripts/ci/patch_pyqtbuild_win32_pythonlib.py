#!/usr/bin/env python3
"""MinGW/g++ ignoriert das '#pragma comment(lib, ...)', das Python.h unter
Windows einbettet, um pythonXY.lib automatisch zu linken -- das ist reine
MSVC-Funktionalitaet. PyQt-builder fuegt deswegen auf win32 nur den
Bibliotheks-Suchpfad (-L...) hinzu, nie ein explizites -lpythonX.Y, und
verlaesst sich stillschweigend auf das (MSVC-only) Pragma. Dieses Skript
ergaenzt den fehlenden expliziten Link-Flag direkt im installierten
PyQt-builder-Quellcode.
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
        "Erwarteter Codeblock in pyqtbuild/builder.py nicht gefunden -- "
        "PyQt-builder-Version hat sich vermutlich geaendert, Patch muss "
        "angepasst werden."
    )

text = text.replace(old, new)
open(path, "w", encoding="utf-8").write(text)
print(f"Patched: {path}")
