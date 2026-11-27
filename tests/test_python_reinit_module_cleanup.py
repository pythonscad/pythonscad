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
  * the user's own package next to the script *was* reloaded,
  * a user module that *shadows a stdlib name* (``mailbox.py``) is also
    reloaded — protection must check the module's origin, not just the name.

On Linux builds against a ``/usr/lib`` Python the old code happened to keep
the stdlib, so the first two checks only bite on Windows/macOS layouts; the
reload checks guard against over-eager protection everywhere.

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
# A top-level name that is also in sys.stdlib_module_names on Python >= 3.10.
STDLIB_SHADOW = "mailbox"


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
        diagnostics_log = root / "diagnostics.txt"
        imports_log = root / "imports.txt"
        shadow_log = root / "shadow.txt"

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

        # User file that deliberately shares a stdlib top-level name. It must
        # still be dropped on re-init so edits are picked up.
        (root / f"{STDLIB_SHADOW}.py").write_text(
            "from pathlib import Path\n"
            f"with Path({str(shadow_log)!r}).open('a', encoding='utf-8') as _f:\n"
            "    _f.write('shadow\\n')\n"
            "VALUE = 1\n",
            encoding="utf-8",
        )

        main_py = root / "main.py"
        main_py.write_text(
            "import sys\n"
            "from pathlib import Path\n"
            "import pkgutil\n"
            "import importlib, importlib.machinery\n"
            "import json\n"
            "# Exercise Windows prefix/file paths with opposite separators.\n"
            "if sys.platform == 'win32':\n"
            "    src, dst = ('\\\\', '/') if '\\\\' in json.__file__ else ('/', '\\\\')\n"
            "    sys.prefix = sys.prefix.replace(src, dst)\n"
            "    sys.base_prefix = sys.base_prefix.replace(src, dst)\n"
            f"import {PKG}\n"
            f"sys.path.insert(0, {str(root)!r})\n"
            "# Startup may already have imported the real stdlib module.\n"
            "# Evict it only once so the fixture can shadow it; subsequent\n"
            "# frames must rely on PythonSCAD's cleanup to reload the shadow.\n"
            "if not getattr(sys, '_pythonscad_shadow_probe_started', False):\n"
            f"    sys.modules.pop({STDLIB_SHADOW!r}, None)\n"
            "    sys._pythonscad_shadow_probe_started = True\n"
            f"import {STDLIB_SHADOW}\n"
            "from openscad import cube, show\n"
            "consistent = importlib.machinery is sys.modules['importlib.machinery']\n"
            "prev = getattr(sys, '_pythonscad_reinit_probe', None)\n"
            "stdlib_kept = prev is None or prev == id(json)\n"
            "sys._pythonscad_reinit_probe = id(json)\n"
            f"with Path({str(diagnostics_log)!r}).open('a', encoding='utf-8') as f:\n"
            "    f.write(f'json={json.__file__!r} prefix={sys.prefix!r} base_prefix={sys.base_prefix!r} stdlib_dir={getattr(sys, \"_stdlib_dir\", None)!r} path={sys.path!r}\\n')\n"
            f"with Path({str(frames_log)!r}).open('a', encoding='utf-8') as f:\n"
            "    f.write(f'consistent={consistent} stdlib_kept={stdlib_kept}\\n')\n"
            f"show(cube([{PKG}.sub.VALUE + {STDLIB_SHADOW}.VALUE] * 3))\n",
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
        diagnostics = diagnostics_log.read_text(encoding="utf-8") if diagnostics_log.is_file() else ""
        imports = imports_log.read_text(encoding="utf-8").splitlines() if imports_log.is_file() else []
        shadows = shadow_log.read_text(encoding="utf-8").splitlines() if shadow_log.is_file() else []

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
    if len(shadows) != FRAMES:
        failures.append(
            f"stdlib-shadowing user module ({STDLIB_SHADOW}.py) imported "
            f"{len(shadows)} time(s), expected {FRAMES} "
            "(protection must check origin, not just the name)"
        )

    if failures:
        print(f"Interpreter diagnostics:\n{diagnostics}", file=sys.stderr)
        for f in failures:
            print(f"FAIL: {f}", file=sys.stderr)
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
