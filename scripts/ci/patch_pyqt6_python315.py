#!/usr/bin/env python3
"""Make older PyQt6 sdists compile against Python 3.15.

Python 3.15 removed PyWeakref_GetObject(). Its replacement,
PyWeakref_GetRef(), is only part of the limited API from 3.13 on, but PyQt6
builds with Py_LIMITED_API=0x030a0000, so it is not declared either.
Dereferencing the weak reference by calling it is part of the 3.10 limited
API and behaves the same: it returns the referent, or None once it is dead.

PyWeakref_GetObject() returned a borrowed reference, so the new reference
from the call is released again to keep PyQt's ownership rules unchanged.
The replacement is guarded by PY_VERSION_HEX, so older Pythons keep the
original code path.

Usage: patch_pyqt6_python315.py <pyqt6-source-dir>
"""
from __future__ import annotations

import sys
from pathlib import Path

TARGET = "qpy/QtCore/qpycore_pyqtslot.cpp"
OLD = "        return PyWeakref_GetObject(mself_wr);\n"
NEW = """\
    {
#if PY_VERSION_HEX >= 0x030f0000
        // PyWeakref_GetObject() was removed in Python 3.15 (patched by
        // PythonSCAD CI, scripts/ci/patch_pyqt6_python315.py).
        PyObject *obj = PyObject_CallNoArgs(mself_wr);
        Py_XDECREF(obj);
        return obj;
#else
        return PyWeakref_GetObject(mself_wr);
#endif
    }
"""


def main() -> int:
    if len(sys.argv) != 2:
        print(__doc__.strip().splitlines()[-1], file=sys.stderr)
        return 1

    path = Path(sys.argv[1]) / TARGET
    source = path.read_text()
    if "patch_pyqt6_python315.py" in source:
        print(f"{path}: already patched")
        return 0
    if OLD not in source:
        # A newer PyQt6 no longer uses the removed API.
        print(f"{path}: PyWeakref_GetObject not found, nothing to patch")
        return 0
    path.write_text(source.replace(OLD, NEW, 1))
    print(f"{path}: replaced PyWeakref_GetObject for Python >= 3.15")
    return 0


if __name__ == "__main__":
    sys.exit(main())
