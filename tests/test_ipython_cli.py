#!/usr/bin/env python3
"""CI smoke test for `pythonscad --ipython`.

Pipes a tiny script into the binary and asserts that:
  * the process exits with code 0 within a reasonable timeout,
  * the script ran (the OK marker appears on stdout),
  * the output carries an IPython fingerprint (so we know the real
    IPython shell was used and not the fallback REPL).

When IPython is *not* installed in the test environment, the binary is
expected to print the fallback diagnostic to stderr and still execute
the piped script via the basic REPL. In that case the IPython
fingerprint check is skipped, but the OK marker check still runs.

Usage:
    test_ipython_cli.py <path-to-pythonscad>
"""
from __future__ import annotations

import os
import subprocess
import sys


# NB: no trailing `exit` line. Both the basic REPL and IPython exit
# cleanly when stdin reaches EOF, and an explicit `exit` (without parens)
# would either no-op or trip IPython's "Use exit()..." advisory message,
# both of which make the test less deterministic.
SCRIPT = (
    "from pythonscad import cube\n"
    "c = cube([1, 1, 1])\n"
    "print('OK', type(c).__name__)\n"
)

TIMEOUT_SECONDS = 60


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: test_ipython_cli.py <pythonscad-binary>", file=sys.stderr)
        return 2

    pythonscad = sys.argv[1]
    if not os.path.isfile(pythonscad):
        print(f"binary not found: {pythonscad}", file=sys.stderr)
        return 2

    proc = subprocess.run(
        [pythonscad, "--ipython"],
        input=SCRIPT,
        capture_output=True,
        text=True,
        timeout=TIMEOUT_SECONDS,
    )

    print("===== stdout =====")
    print(proc.stdout)
    print("===== stderr =====", file=sys.stderr)
    print(proc.stderr, file=sys.stderr)

    if proc.returncode != 0:
        print(
            f"FAIL: pythonscad --ipython exited with {proc.returncode}",
            file=sys.stderr,
        )
        return 1

    if "OK" not in proc.stdout:
        print(
            "FAIL: piped script did not run (missing 'OK' marker on stdout)",
            file=sys.stderr,
        )
        return 1

    fallback_msg = "IPython is not installed"
    is_fallback = fallback_msg in proc.stderr

    if is_fallback:
        print(
            "INFO: IPython not installed in this environment; the fallback "
            "REPL handled the script. This is acceptable for the smoke "
            "test, but Layer-1 coverage of the real IPython prompt requires "
            "IPython to be available.",
            file=sys.stderr,
        )
    else:
        if "In [" not in proc.stdout and "IPython" not in proc.stdout:
            print(
                "FAIL: real IPython appears to be installed but no IPython "
                "fingerprint (`In [`, `IPython`) was found on stdout. Either "
                "IPython failed to launch silently, or the bootstrap "
                "snippet regressed.",
                file=sys.stderr,
            )
            return 1

    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
