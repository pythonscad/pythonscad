#!/usr/bin/env python3
"""Regression test for the sys.modules cleanup on interpreter re-use.

Background
----------
When the embedded interpreter is re-used for another script run (GUI re-run,
or ``--animate`` on the command line, which calls ``initPython()`` once per
frame), ``initPython()`` drops user modules from ``sys.modules`` so edited
libraries are reloaded.  The old filter decided per entry from ``repr()``
path fragments ("usr/lib", "site-packages", ...).  On Windows and macOS the
bundled stdlib does not live below ``usr/lib``, so most of the stdlib was
dropped - including ``importlib`` but *not* the frozen
``importlib.machinery``.  The next ``import pkgutil`` (pulled in by e.g.
``PyQt6/__init__.py``) then failed with::

    AttributeError: module 'importlib' has no attribute 'machinery'

This test renders two animation frames and asserts, from inside the second
frame, that

  * ``import pkgutil`` works and ``importlib.machinery`` is consistent,
  * stdlib modules survived the cleanup (same ``json`` module object),
  * the user's own package next to the script *was* reloaded.

On Linux builds against a ``/usr/lib`` Python the old code happened to keep
the stdlib, so the first two checks only bite on Windows/macOS layouts; the
reload check guards against over-eager protection everywhere.

Usage:
    test_python_reinit_module_cleanup.py <path-to-pythonscad>
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path


TIMEOUT_SECONDS = 120
FRAMES = 2
PKG = "pythonscad_reinit_probe_pkg"


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: test_python_reinit_module_cleanup.py <pythonscad-binary>", file=sys.stderr)
        return 2

    pythonscad = sys.argv[1]
    if not os.path.isfile(pythonscad):
        print(f"binary not found: {pythonscad}", file=sys.stderr)
        return 2

    with tempfile.TemporaryDirectory(prefix="pythonscad-reinit-") as tmp:
        root = Path(tmp)
        frames_log = root / "frames.txt"
        imports_log = root / "imports.txt"

        pkg_dir = root / PKG
        pkg_dir.mkdir()
        (pkg_dir / "__init__.py").write_text(
            "from pathlib import Path\n"
            f"with Path({str(imports_log)!r}).open('a', encoding='utf-8') as _f:\n"
            "    _f.write('import\\n')\n"
            "from . import sub\n",
            encoding="utf-8",
        )
        (pkg_dir / "sub.py").write_text("VALUE = 1\n", encoding="utf-8")

        main_py = root / "main.py"
        main_py.write_text(
            "import sys\n"
            "from pathlib import Path\n"
            "import pkgutil\n"
            "import importlib, importlib.machinery\n"
            "import json\n"
            f"import {PKG}\n"
            "from openscad import cube, show\n"
            "consistent = importlib.machinery is sys.modules['importlib.machinery']\n"
            "prev = getattr(sys, '_pythonscad_reinit_probe', None)\n"
            "stdlib_kept = prev is None or prev == id(json)\n"
            "sys._pythonscad_reinit_probe = id(json)\n"
            f"with Path({str(frames_log)!r}).open('a', encoding='utf-8') as f:\n"
            "    f.write(f'consistent={consistent} stdlib_kept={stdlib_kept}\\n')\n"
            f"show(cube([{PKG}.sub.VALUE] * 3))\n",
            encoding="utf-8",
        )

        proc = subprocess.run(
            [
                pythonscad,
                "--trust-python",
                "--animate",
                str(FRAMES),
                "-o",
                str(root / "out.stl"),
                str(main_py),
            ],
            cwd=root,
            capture_output=True,
            text=True,
            timeout=TIMEOUT_SECONDS,
        )
        frames = frames_log.read_text(encoding="utf-8").splitlines() if frames_log.is_file() else []
        imports = imports_log.read_text(encoding="utf-8").splitlines() if imports_log.is_file() else []

    print("===== stdout =====")
    print(proc.stdout)
    print("===== stderr =====", file=sys.stderr)
    print(proc.stderr, file=sys.stderr)

    failures = []
    if proc.returncode != 0:
        failures.append(f"pythonscad exited with {proc.returncode}")
    if "AttributeError" in proc.stdout + proc.stderr:
        failures.append("an AttributeError was reported (importlib/sys.modules inconsistency?)")
    expected = ["consistent=True stdlib_kept=True"] * FRAMES
    if frames != expected:
        failures.append(f"frame log {frames!r}, expected {expected!r}")
    if len(imports) != FRAMES:
        failures.append(
            f"user package imported {len(imports)} time(s), expected {FRAMES} "
            "(it must be reloaded on every re-init)"
        )

    if failures:
        for f in failures:
            print(f"FAIL: {f}", file=sys.stderr)
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
