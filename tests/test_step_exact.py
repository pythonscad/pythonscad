#!/usr/bin/env python3
"""Exact (B-Rep) STEP export, only registered with ENABLE_OCCT.

Checks that primitives are written with analytic surfaces instead of a faceted
mesh, that transforms and booleans keep them exact, and that models mixing
exact and mesh geometry still fall back to the mesh-based STEP writer.
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
import tempfile

TIMEOUT_SECONDS = 120


def count(text: str, entity: str) -> int:
    # Also matches entities inside complex instances such as
    # "#22 = ( BOUNDED_SURFACE() B_SPLINE_SURFACE(...) RATIONAL_B_SPLINE_SURFACE(...) ... )"
    return len(re.findall(r"(?<![A-Z0-9_])" + entity + r"\s*\(", text))


def expect(**want):
    """Build a checker: entity=n means exactly n, entity=(lo, hi) a range (hi None = open)."""

    def check(text: str) -> str | None:
        errors = []
        for entity, n in want.items():
            have = count(text, entity)
            lo, hi = n if isinstance(n, tuple) else (n, n)
            if have < lo or (hi is not None and have > hi):
                errors.append(f"{entity}: expected {n}, got {have}")
        return "; ".join(errors) or None

    return check


EXACT = (1, None)  # at least one
NONE = 0

# Faceted fallback marker: the mesh writer emits lots of planar faces and no curved surfaces.
MESH = expect(SPHERICAL_SURFACE=NONE, CYLINDRICAL_SURFACE=NONE, ADVANCED_FACE=(10, None))

TETRA = "polyhedron([[0,0,0],[4,0,0],[0,4,0],[0,0,4]], [[0,1,2],[0,3,1],[0,2,3],[1,3,2]])"

CASES = [
    # name, script body, checker
    ("sphere", "show(sphere(10))", expect(SPHERICAL_SURFACE=1, MANIFOLD_SOLID_BREP=1)),
    ("sphere-fn-independent", "show(sphere(10, fn=6))", expect(SPHERICAL_SURFACE=1, ADVANCED_FACE=1)),
    (
        "rotated-translated-sphere",
        "show(sphere(5).rotate([30, 40, 50]).translate([10, 20, 30]))",
        expect(SPHERICAL_SURFACE=1),
    ),
    (
        "nonuniform-scale",
        "show(sphere(5).scale([2, 1, 1]))",
        expect(SPHERICAL_SURFACE=NONE, RATIONAL_B_SPLINE_SURFACE=EXACT),
    ),
    ("cube", "show(cube([10, 20, 30]))", expect(PLANE=6, ADVANCED_FACE=6, MANIFOLD_SOLID_BREP=1)),
    ("cube-center", "show(cube(10, center=True))", expect(PLANE=6, ADVANCED_FACE=6)),
    ("cylinder", "show(cylinder(h=10, r=3, fn=5))", expect(CYLINDRICAL_SURFACE=1, PLANE=2)),
    ("cone", "show(cylinder(h=10, r1=4, r2=2))", expect(CONICAL_SURFACE=1, PLANE=2)),
    ("pointed-cone", "show(cylinder(h=10, r1=4, r2=0))", expect(CONICAL_SURFACE=1, PLANE=1)),
    (
        "cylinder-sector",
        "show(cylinder(h=10, r=5, angle=90))",
        expect(CYLINDRICAL_SURFACE=1, PLANE=4, ADVANCED_FACE=5),
    ),
    (
        "difference",
        "show(sphere(10) - cube(10))",
        expect(SPHERICAL_SURFACE=EXACT, PLANE=3, MANIFOLD_SOLID_BREP=1),
    ),
    (
        "intersection",
        "show(sphere(7) & cube(10, center=True))",
        expect(SPHERICAL_SURFACE=EXACT, PLANE=6, MANIFOLD_SOLID_BREP=1),
    ),
    (
        # touching cubes fuse into one box: coplanar faces must be merged
        "union-merges-faces",
        "show(cube(10) | cube(10).translate([10, 0, 0]))",
        expect(PLANE=6, ADVANCED_FACE=6, MANIFOLD_SOLID_BREP=1),
    ),
    (
        "drilled-plate",
        "show(cube([30, 20, 5]) - cylinder(h=20, r=3, center=True).translate([15, 10, 0]))",
        expect(CYLINDRICAL_SURFACE=1, PLANE=6, MANIFOLD_SOLID_BREP=1),
    ),
    (
        "transformed-boolean",
        "show((cube(10) - sphere(4).translate([10, 10, 10])).rotate([0, 0, 45]).translate([5, 0, 0]))",
        expect(SPHERICAL_SURFACE=EXACT, PLANE=6),
    ),
    (
        # sphere pokes through all six cube faces; regression for an OCCT edge-unification bug
        "union-of-three",
        "show(union(sphere(6), cube(10, center=True), cylinder(h=30, r=2, center=True)))",
        expect(SPHERICAL_SURFACE=EXACT, CYLINDRICAL_SURFACE=EXACT, MANIFOLD_SOLID_BREP=1),
    ),
    ("mixed-falls-back-to-mesh", f"show(sphere(5) | {TETRA})", MESH),
    (
        "mixed-boolean-poisons-parent",
        f"show((sphere(5) | {TETRA}) - cube(2))",
        MESH,
    ),
]

# Same models, but written through the Python export() call (allownef=false path)
EXPORT_CASES = [
    ("export-call-sphere", "sphere(10)", expect(SPHERICAL_SURFACE=1)),
    ("export-call-difference", "sphere(10) - cube(10)", expect(SPHERICAL_SURFACE=EXACT, PLANE=3)),
]


def run(pythonscad: str, workdir: str, name: str, script_text: str, out: str) -> tuple[str | None, str]:
    script = os.path.join(workdir, f"{name}.py")
    with open(script, "w", encoding="utf-8") as f:
        f.write(script_text)
    target = os.path.join(workdir, f"{name}-render.step")
    proc = subprocess.run(
        [pythonscad, "--trust-python", "-o", target, script],
        capture_output=True,
        text=True,
        timeout=TIMEOUT_SECONDS,
    )
    log = proc.stdout + proc.stderr
    if proc.returncode != 0:
        return None, log + f"\nexit code {proc.returncode}"
    if not os.path.isfile(out):
        return None, log + f"\nno output file {out}"
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
        jobs = []
        for name, body, check in CASES:
            out = os.path.join(workdir, f"{name}-render.step")
            jobs.append((name, f"from pythonscad import *\n{body}\n", out, check))
        for name, obj, check in EXPORT_CASES:
            out = os.path.join(workdir, f"{name}.step")
            script = f"from pythonscad import *\nexport({obj}, {out!r})\nshow(cube(1))\n"
            jobs.append((name, script, out, check))

        for name, script, out, check in jobs:
            text, log = run(pythonscad, workdir, name, script, out)
            error = "export failed" if text is None else check(text)
            if error:
                failures += 1
                print(f"FAIL {name}: {error}\n{log}", file=sys.stderr)
            else:
                print(f"ok   {name}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
