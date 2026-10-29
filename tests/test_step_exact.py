#!/usr/bin/env python3
"""Exact (B-Rep) STEP export, only registered with ENABLE_OCCT; runs with --enable=occt.

Checks that primitives are written with analytic surfaces instead of a faceted
mesh, that transforms and booleans keep them exact, and that models mixing
exact and mesh geometry still fall back to the mesh-based STEP writer.
"""
from __future__ import annotations

import math
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
    (
        # Spheres centred on the edges/corner of a plate: OCCT's STEP translation keeps the
        # wrong side of the closed sphere face here; the exporter must detect that and retry
        # with split faces instead of falling back to the mesh or writing a wrong file.
        "spheres-on-plate-edges",
        "show(union(cube([60, 60, 5]), *[sphere(3).translate([20*i, 20*j, 5]) for i in range(3) for j in range(3)]))",
        expect(SPHERICAL_SURFACE=(9, None), MANIFOLD_SOLID_BREP=1),
    ),
    # --- linear_extrude / rotate_extrude ---
    ("extrude-circle", "show(linear_extrude(circle(5), height=10))", expect(CYLINDRICAL_SURFACE=1, PLANE=2)),
    (
        "extrude-ring-centered",
        "show(linear_extrude(circle(5) - circle(3), height=4, center=True))",
        expect(CYLINDRICAL_SURFACE=2, PLANE=2, MANIFOLD_SOLID_BREP=1),
    ),
    (
        "extrude-polygon",
        "show(linear_extrude(polygon([[0,0],[10,0],[10,10],[5,4],[0,10]]), height=2))",
        expect(PLANE=7, CYLINDRICAL_SURFACE=NONE, MANIFOLD_SOLID_BREP=1),
    ),
    ("extrude-sector", "show(linear_extrude(circle(5, angle=90), height=3))", expect(CYLINDRICAL_SURFACE=1, PLANE=4)),
    (
        "extrude-2d-boolean",
        "show(linear_extrude(square([20, 10]) - circle(2).translate([5, 5, 0]), height=4))",
        expect(CYLINDRICAL_SURFACE=1, PLANE=6),
    ),
    ("extrude-placed-in-3d", "show(linear_extrude(circle(3).rotate([90, 0, 0]), height=5))", expect(CYLINDRICAL_SURFACE=1)),
    (
        "extrude-ellipse",
        "show(linear_extrude(circle(3).scale([2, 1, 1]), height=5))",
        # a non-uniform scale turns the circle into an exact rational B-spline (caps included)
        expect(RATIONAL_B_SPLINE_CURVE=EXACT, MANIFOLD_SOLID_BREP=1),
    ),
    ("extrude-cone", "show(linear_extrude(circle(5), height=10, scale=[0, 0]))", expect(CONICAL_SURFACE=1, PLANE=1)),
    (
        "extrude-frustum",
        "show(linear_extrude(square(4, center=True), height=6, scale=[0.5, 0.5]))",
        # the loft stores its (planar) side faces as ruled B-spline surfaces
        expect(PLANE=2, ADVANCED_FACE=6, MANIFOLD_SOLID_BREP=1),
    ),
    (
        "extrude-then-boolean",
        "show(linear_extrude(circle(5), height=10) - sphere(4).translate([0, 0, 10]))",
        expect(CYLINDRICAL_SURFACE=1, SPHERICAL_SURFACE=EXACT, MANIFOLD_SOLID_BREP=1),
    ),
    ("extrude-twist-falls-back", "show(linear_extrude(square(4), height=5, twist=90))", MESH),
    ("extrude-nonuniform-scale-falls-back", "show(linear_extrude(circle(5), height=5, scale=[0.5, 1]))", MESH),
    ("rotate-torus", "show(rotate_extrude(circle(2).translate([5, 0, 0])))", expect(TOROIDAL_SURFACE=1)),
    (
        "rotate-partial",
        "show(rotate_extrude(circle(2).translate([5, 0, 0]), angle=180))",
        expect(TOROIDAL_SURFACE=1, PLANE=2, MANIFOLD_SOLID_BREP=1),
    ),
    (
        "rotate-half-disk-is-sphere",
        "show(rotate_extrude(circle(5, angle=180).rotate([0, 0, -90])))",
        expect(SPHERICAL_SURFACE=1, MANIFOLD_SOLID_BREP=1),
    ),
    ("mixed-falls-back-to-mesh", f"show(sphere(5) | {TETRA})", MESH),
    (
        "mixed-boolean-poisons-parent",
        f"show((sphere(5) | {TETRA}) - cube(2))",
        MESH,
    ),
]

# STEP import through OpenCASCADE: (name, model written to STEP first, script using it as {src}, checker)
IMPORT_CASES = [
    ("import-reexport", "sphere(10) - cube(10)", "show(osimport({src}))", expect(SPHERICAL_SURFACE=EXACT, PLANE=3)),
    (
        "import-boolean",
        "sphere(10) - cube(10)",
        "show(osimport({src}) - cube(6, center=True))",
        expect(SPHERICAL_SURFACE=EXACT, PLANE=(6, None), MANIFOLD_SOLID_BREP=1),
    ),
    # a faceted file from the mesh writer becomes an exact planar solid
    ("import-faceted", TETRA, "show(osimport({src}))", expect(PLANE=4, MANIFOLD_SOLID_BREP=1)),
]

# Imported STEP rendered as mesh: (name, model, fn, exact volume, max relative error, min triangles)
IMPORT_MESH_CASES = [
    ("import-mesh-coarse", "sphere(10) - cube(10)", 8, 4 / 3 * math.pi * 1000 * 7 / 8, 0.15, 50),
    ("import-mesh-fine", "sphere(10) - cube(10)", 64, 4 / 3 * math.pi * 1000 * 7 / 8, 0.005, 2000),
]

# Same models, but written through the Python export() call (allownef=false path)
EXPORT_CASES = [
    ("export-call-sphere", "sphere(10)", expect(SPHERICAL_SURFACE=1)),
    ("export-call-difference", "sphere(10) - cube(10)", expect(SPHERICAL_SURFACE=EXACT, PLANE=3)),
]


def stl_volume(text: str) -> tuple[float, int]:
    """Volume and triangle count of an ASCII STL."""
    v = [tuple(map(float, m)) for m in re.findall(r"vertex\s+(\S+)\s+(\S+)\s+(\S+)", text)]
    volume = 0.0
    for i in range(0, len(v) - 2, 3):
        a, b, c = v[i : i + 3]
        volume += (
            a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])
        ) / 6
    return volume, len(v) // 3


def run(
    pythonscad: str, workdir: str, name: str, script_text: str, out: str, occt: bool = True
) -> tuple[str | None, str]:
    script = os.path.join(workdir, f"{name}.py")
    with open(script, "w", encoding="utf-8") as f:
        f.write(script_text)
    target = os.path.join(workdir, f"{name}-render.step")
    proc = subprocess.run(
        [pythonscad, "--trust-python", *(["--enable=occt"] if occt else []), "-o", target, script],
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

        for name, model, body, check in IMPORT_CASES:
            src = os.path.join(workdir, f"{name}-source.step")
            text, log = run(pythonscad, workdir, f"{name}-source", f"from pythonscad import *\nexport({model}, {src!r})\nshow(cube(1))\n", src)
            if text is None:
                failures += 1
                print(f"FAIL {name}: writing the source file failed\n{log}", file=sys.stderr)
                continue
            out = os.path.join(workdir, f"{name}-render.step")
            jobs.append((name, "from pythonscad import *\n" + body.format(src=repr(src)) + "\n", out, check))

        for name, model, fn, exact, tolerance, min_triangles in IMPORT_MESH_CASES:
            src = os.path.join(workdir, f"{name}-source.step")
            text, log = run(pythonscad, workdir, f"{name}-source", f"from pythonscad import *\nexport({model}, {src!r})\nshow(cube(1))\n", src)
            stl = os.path.join(workdir, f"{name}.stl")
            script = f"from pythonscad import *\nexport(osimport({src!r}, fn={fn}), {stl!r})\nshow(cube(1))\n"
            mesh, log2 = run(pythonscad, workdir, name, script, stl) if text is not None else (None, log)
            if mesh is None:
                failures += 1
                print(f"FAIL {name}: export failed\n{log2}", file=sys.stderr)
                continue
            volume, triangles = stl_volume(mesh)
            error = abs(volume - exact) / exact
            if error > tolerance or triangles < min_triangles:
                failures += 1
                print(f"FAIL {name}: volume {volume:.2f} vs {exact:.2f} ({error:.2%}), {triangles} triangles", file=sys.stderr)
            else:
                print(f"ok   {name} (volume error {error:.2%}, {triangles} triangles)")

        # Without the experimental feature "occt" everything must behave as without OpenCASCADE.
        for name, body in (("flag-off-sphere", "show(sphere(10))"), ("flag-off-difference", "show(sphere(10) - cube(10))")):
            out = os.path.join(workdir, f"{name}-render.step")
            text, log = run(pythonscad, workdir, name, f"from pythonscad import *\n{body}\n", out, occt=False)
            error = "export failed" if text is None else MESH(text)
            if error:
                failures += 1
                print(f"FAIL {name}: {error}\n{log}", file=sys.stderr)
            else:
                print(f"ok   {name}")

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
