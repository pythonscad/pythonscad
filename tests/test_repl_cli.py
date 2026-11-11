#!/usr/bin/env python3
"""CI smoke test for `pythonscad --repl`.

Pipes a tiny script into the binary and asserts that:
  * the basic embedded Python REPL is reachable (exit 0 within timeout),
  * the script ran (the OK marker appears on stdout),
  * NO `IPython is not installed` diagnostic is emitted -- `--repl` is
    the explicit, intentional path to the basic REPL and must never
    print the IPython fallback warning.

Usage:
    test_repl_cli.py <path-to-pythonscad>
"""
from __future__ import annotations

import os
import subprocess
import sys


# --repl drops the user into an empty `__main__` namespace; the user
# is responsible for importing whatever they need. The explicit
# `from pythonscad import cube` line below is the contract under test.
#
# Import `ctypes`/`sysconfig` on non-macOS: on Windows this exercises the
# sysconfig path that used to fail in packaged MSYS2 builds when
# `sys._is_mingw` was missing; on Linux it is cheap coverage. Keep it
# off Darwin — Homebrew Python 3.15 still SIGBUS-es on that import under
# --repl on macos-26 CI even after the sys.stdin REPL rewrite (#1138).
# Cube import alone remains the macOS contract under test.
_SCRIPT_LINES = []
if sys.platform != "darwin":
    _SCRIPT_LINES.append("import ctypes, sysconfig")
_SCRIPT_LINES.extend(
    [
        "from pythonscad import cube",
        "c = cube([1, 1, 1])",
        "print('OK', type(c).__name__)",
    ]
)
SCRIPT = "\n".join(_SCRIPT_LINES) + "\n"

TIMEOUT_SECONDS = 60

# Piped prints before an exception must still reach stdout: the --repl
# path returns without Py_Finalize, so the bootstrap must flush in a
# finally block (Copilot review on #1148).
FLUSH_ON_ERROR_SCRIPT = "print('MARKER')\nraise RuntimeError('boom')\n"

# User scripts that rebind `sys` in `__main__` must not break the
# bootstrap's post-exec flush (bootstrap runs in a function with its
# own locals).
SYS_REBIND_SCRIPT = "print('OK-rebind')\nsys = 123\n"


def _run_repl(pythonscad: str, script: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [pythonscad, "--repl"],
        input=script,
        capture_output=True,
        text=True,
        timeout=TIMEOUT_SECONDS,
    )


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: test_repl_cli.py <path-to-pythonscad>", file=sys.stderr)
        return 2

    pythonscad = sys.argv[1]
    if not os.path.isfile(pythonscad):
        print(f"binary not found: {pythonscad}", file=sys.stderr)
        return 2

    proc = _run_repl(pythonscad, SCRIPT)

    print("===== stdout =====")
    print(proc.stdout)
    print("===== stderr =====", file=sys.stderr)
    print(proc.stderr, file=sys.stderr)

    if proc.returncode != 0:
        # Negative codes are usually -signal (e.g. -11 = SIGSEGV, -10 = SIGBUS).
        detail = f" (signal {-proc.returncode})" if proc.returncode < 0 else ""
        print(
            f"FAIL: pythonscad --repl exited with {proc.returncode}{detail}",
            file=sys.stderr,
        )
        return 1

    if "OK" not in proc.stdout:
        print(
            "FAIL: piped script did not run (missing 'OK' marker on stdout)",
            file=sys.stderr,
        )
        return 1

    # Check BOTH streams: on Windows the `pythonscad.com` console shim
    # merges stderr into stdout, so a `proc.stderr`-only check would
    # silently miss a regression where --repl started advertising the
    # IPython fallback path on Windows.
    fallback_msg = "IPython is not installed"
    if fallback_msg in proc.stderr or fallback_msg in proc.stdout:
        print(
            "FAIL: --repl emitted the IPython fallback diagnostic; the "
            "explicit --repl path must never advertise IPython.",
            file=sys.stderr,
        )
        return 1

    hook_msg = "Failed calling sys.__interactivehook__"
    if hook_msg in proc.stderr or hook_msg in proc.stdout:
        print(
            "FAIL: --repl emitted the sys.__interactivehook__ startup warning.",
            file=sys.stderr,
        )
        return 1

    flush_proc = _run_repl(pythonscad, FLUSH_ON_ERROR_SCRIPT)
    print("===== flush-on-error stdout =====")
    print(flush_proc.stdout)
    print("===== flush-on-error stderr =====", file=sys.stderr)
    print(flush_proc.stderr, file=sys.stderr)
    if "MARKER" not in flush_proc.stdout:
        print(
            "FAIL: print before exception was not flushed to stdout",
            file=sys.stderr,
        )
        return 1
    if "RuntimeError" not in flush_proc.stderr and "RuntimeError" not in flush_proc.stdout:
        print(
            "FAIL: expected RuntimeError traceback from flush-on-error script",
            file=sys.stderr,
        )
        return 1

    rebind_proc = _run_repl(pythonscad, SYS_REBIND_SCRIPT)
    print("===== sys-rebind stdout =====")
    print(rebind_proc.stdout)
    if rebind_proc.returncode != 0 or "OK-rebind" not in rebind_proc.stdout:
        print(
            "FAIL: rebinding sys in __main__ broke --repl stdout flush",
            file=sys.stderr,
        )
        return 1

    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
