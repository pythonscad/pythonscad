#!/usr/bin/env python3
"""Exact (B-Rep) STEP export of spheres, only registered with ENABLE_OCCT.

Checks that a sphere is written with an analytic SPHERICAL_SURFACE instead of
a faceted mesh, that transforms keep it exact, and that models mixing exact
and mesh geometry still fall back to the mesh-based STEP writer.
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
import tempfile

TIMEOUT_SECONDS = 120

CASES = [
    # name, script body, predicate(step_text) -> error message or None
    (
        "plain-sphere",
        "show(sphere(10))",
        lambda t: None
        if "SPHERICAL_SURFACE" in t and count(t, "ADVANCED_FACE") <= 2
        else "expected one analytic SPHERICAL_SURFACE",
    ),
    (
        "sphere-fn-independent",
        "show(sphere(10, fn=6))",
        lambda t: None
        if "SPHERICAL_SURFACE" in t and count(t, "ADVANCED_FACE") <= 2
        else "coarse $fn must not change the exact export",
    ),
    (
        "rotated-translated-sphere",
        "show(sphere(5).rotate([30, 40, 50]).translate([10, 20, 30]))",
        lambda t: None
        if "SPHERICAL_SURFACE" in t
        else "similarity transforms must keep the SPHERICAL_SURFACE",
    ),
    (
        "nonuniform-scale",
        "show(sphere(5).scale([2, 1, 1]))",
        lambda t: None
        if "SPHERICAL_SURFACE" not in t and "B_SPLINE_SURFACE" in t
        else "non-uniform scale should give an exact B-spline ellipsoid",
    ),
    (
        "mixed-falls-back-to-mesh",
        "show(sphere(5) | cube(4))",
        lambda t: None
        if "SPHERICAL_SURFACE" not in t and count(t, "ADVANCED_FACE") > 10
        else "boolean results are meshes and must use the faceted writer",
    ),
]


def count(text: str, entity: str) -> int:
    return len(re.findall(r"=\s*" + entity + r"\s*\(", text))


def run_case(pythonscad: str, workdir: str, name: str, body: str) -> tuple[str | None, str]:
    script = os.path.join(workdir, f"{name}.py")
    out = os.path.join(workdir, f"{name}.step")
    with open(script, "w", encoding="utf-8") as f:
        f.write("from pythonscad import *\n" + body + "\n")
    proc = subprocess.run(
        [pythonscad, "--trust-python", "-o", out, script],
        capture_output=True,
        text=True,
        timeout=TIMEOUT_SECONDS,
    )
    log = proc.stdout + proc.stderr
    if proc.returncode != 0:
        return None, log + f"\nexit code {proc.returncode}"
    if not os.path.isfile(out):
        return None, log + "\nno output file"
    with open(out, encoding="utf-8", errors="replace") as f:
        return f.read(), log


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: test_step_exact.py <pythonscad-binary>", file=sys.stderr)
        return 2
    pythonscad = sys.argv[1]
    if not os.path.isfile(pythonscad):
        print(f"binary not found: {pythonscad}", file=sys.stderr)
        return 2

    failures = 0
    with tempfile.TemporaryDirectory() as workdir:
        for name, body, check in CASES:
            text, log = run_case(pythonscad, workdir, name, body)
            error = "export failed" if text is None else check(text)
            if error:
                failures += 1
                print(f"FAIL {name}: {error}\n{log}", file=sys.stderr)
            else:
                print(f"ok   {name}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
