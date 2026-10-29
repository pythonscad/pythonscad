#!/usr/bin/env python3
"""Number format of the faceted (mesh based) STEP export.

ISO 10303-21 REALs need a decimal point ("10." not "10") and an upper case
exponent ("1.E-05" not "1e-05"); strict STEP readers reject anything else.
The values must also keep enough digits to survive the round trip.
Runs with and without OpenCASCADE: mesh-only models always use this writer.
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
import tempfile

TIMEOUT_SECONDS = 120

# A tetrahedron with awkward coordinates: integers, many digits, tiny and negative values.
SCRIPT = """from pythonscad import *
show(polyhedron(
    [[0, 0, 0], [123.456789012, 0, 0], [0, 1e-5 + 40, -0.000012345], [-7, 3, 50]],
    [[0, 1, 2], [0, 3, 1], [0, 2, 3], [1, 3, 2]]))
"""

STEP_REAL = re.compile(r"^[+-]?\d+\.\d*(E[+-]?\d+)?$")
NUMERIC_ENTITY = re.compile(r"=\s*(CARTESIAN_POINT|DIRECTION)\s*\('[^']*',\s*\(([^)]*)\)\s*\)")
VECTOR_ENTITY = re.compile(r"=\s*VECTOR\s*\('[^']*',\s*#\d+\s*,\s*([^)]*)\)")


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: test_step_faceted.py <pythonscad-binary>", file=sys.stderr)
        return 2
    pythonscad = sys.argv[1]

    with tempfile.TemporaryDirectory() as workdir:
        script = os.path.join(workdir, "tetra.py")
        out = os.path.join(workdir, "tetra.step")
        with open(script, "w", encoding="utf-8") as f:
            f.write(SCRIPT)
        proc = subprocess.run(
            [pythonscad, "--trust-python", "-o", out, script],
            capture_output=True,
            text=True,
            timeout=TIMEOUT_SECONDS,
        )
        if proc.returncode != 0 or not os.path.isfile(out):
            print(f"FAIL export: exit {proc.returncode}\n{proc.stdout}{proc.stderr}", file=sys.stderr)
            return 1
        with open(out, encoding="utf-8") as f:
            text = f.read()

    numbers = []
    for _, values in NUMERIC_ENTITY.findall(text):
        numbers += [v.strip() for v in values.split(",")]
    numbers += [v.strip() for v in VECTOR_ENTITY.findall(text)]
    if not numbers:
        print("FAIL: no CARTESIAN_POINT/DIRECTION/VECTOR found", file=sys.stderr)
        return 1

    bad = sorted({n for n in numbers if not STEP_REAL.match(n)})
    failures = 0
    if bad:
        failures += 1
        print(f"FAIL: {len(bad)} values are not STEP REALs, e.g. {bad[:5]}", file=sys.stderr)
    else:
        print(f"ok   all {len(numbers)} values are STEP REALs")

    for expected in ("123.456789012", "40.00001", "-1.2345E-05"):
        if expected not in numbers:
            failures += 1
            print(f"FAIL: coordinate {expected} not written with full precision", file=sys.stderr)
        else:
            print(f"ok   precision kept for {expected}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
