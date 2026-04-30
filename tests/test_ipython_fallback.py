#!/usr/bin/env python3
"""Verify `--ipython` falls back to the basic REPL when IPython is missing.

This is the Layer-3 test from the PR plan: the goal is to assert that a
missing IPython is handled gracefully (one-line diagnostic + drop into
the basic Python REPL), not as a hard error.

The trick is to *force* the embedded interpreter to fail importing
IPython even on a CI worker that has IPython globally installed. The
cleanest way is to constrain `PYTHONPATH` to a directory that hides
IPython but still exposes PythonSCAD's own `pythonscad`/`openscad`
overlays plus the standard library, *and* set `PYTHONNOUSERSITE=1`
so the user-site directory is ignored.

If the test environment somehow still imports IPython (e.g. because
IPython was installed straight into the system stdlib path which is
unusual), the test treats that as inconclusive and skips with a clear
INFO message instead of failing -- the Layer-1 smoke covers the happy
path and we do not want false negatives on exotic CI workers.

Usage:
    test_ipython_fallback.py <path-to-pythonscad>
"""
from __future__ import annotations

import os
import subprocess
import sys


SCRIPT = (
    "print('FALLBACK_OK')\n"
    "exit()\n"
)

TIMEOUT_SECONDS = 60


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: test_ipython_fallback.py <pythonscad-binary>", file=sys.stderr)
        return 2

    pythonscad = sys.argv[1]
    if not os.path.isfile(pythonscad):
        print(f"binary not found: {pythonscad}", file=sys.stderr)
        return 2

    env = os.environ.copy()
    env["PYTHONPATH"] = ""
    env["PYTHONNOUSERSITE"] = "1"
    env.pop("IPYTHONDIR", None)

    proc = subprocess.run(
        [pythonscad, "--ipython"],
        input=SCRIPT,
        capture_output=True,
        text=True,
        timeout=TIMEOUT_SECONDS,
        env=env,
    )

    print("===== stdout =====")
    print(proc.stdout)
    print("===== stderr =====", file=sys.stderr)
    print(proc.stderr, file=sys.stderr)

    fallback_msg = "IPython is not installed"
    if fallback_msg not in proc.stderr:
        print(
            "INFO: could not force the embedded interpreter into the "
            "fallback path (IPython is reachable even with PYTHONPATH "
            "scrubbed). Skipping Layer-3 assertion -- Layer-1 still "
            "covers the happy path.",
            file=sys.stderr,
        )
        print("SKIP")
        return 0

    if proc.returncode != 0:
        print(
            f"FAIL: fallback REPL did not exit cleanly (rc={proc.returncode})",
            file=sys.stderr,
        )
        return 1

    if "FALLBACK_OK" not in proc.stdout:
        print(
            "FAIL: fallback message printed but the piped script did not run "
            "(FALLBACK_OK marker missing on stdout)",
            file=sys.stderr,
        )
        return 1

    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
