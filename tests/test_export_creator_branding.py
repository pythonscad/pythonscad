#!/usr/bin/env python3
"""Assert raw EXPORT_CREATOR branding on POV / PDF / 3MF before any golden normalize.

Regression goldens stay OpenSCAD-shaped and normalizers rewrite the live
PythonSCAD string, so a revert of ``EXPORT_CREATOR`` would still pass those
compares. This harness exports each format once and checks the unnormalized
bytes for ``PythonSCAD (https://pythonscad.org/)``.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import test_cmdline_tool as tct  # noqa: E402

TIMEOUT_SECONDS = 120


def _export(pythonscad: str, script: Path, out: Path) -> None:
    proc = subprocess.run(
        [
            pythonscad,
            "--enable=predictible-output",
            "--backend=manifold",
            "-o",
            str(out),
            str(script),
        ],
        capture_output=True,
        text=True,
        timeout=TIMEOUT_SECONDS,
        check=False,
    )
    if proc.returncode != 0 or not out.is_file():
        sys.stderr.write(proc.stdout or "")
        sys.stderr.write(proc.stderr or "")
        raise SystemExit(
            f"export failed ({proc.returncode}): {script.name} -> {out.name}"
        )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pythonscad")
    parser.add_argument("--pov-script", type=Path, required=True)
    parser.add_argument("--pdf-script", type=Path, required=True)
    parser.add_argument("--threemf-script", type=Path, default=None)
    args = parser.parse_args()

    if not os.path.isfile(args.pythonscad):
        print(f"binary not found: {args.pythonscad}", file=sys.stderr)
        return 2

    with tempfile.TemporaryDirectory(prefix="export-creator-") as tmp:
        tmpdir = Path(tmp)
        checks = [
            ("pov", args.pov_script, tmpdir / "out.pov"),
            ("pdf", args.pdf_script, tmpdir / "out.pdf"),
        ]
        if args.threemf_script is not None:
            checks.append(("3mf", args.threemf_script, tmpdir / "out.3mf"))

        for label, script, out in checks:
            print(f"exporting {label} from {script.name} ...")
            _export(args.pythonscad, script, out)
            tct.assert_raw_export_creator(str(out))
            print(f"OK {label}: raw EXPORT_CREATOR present")

    return 0


if __name__ == "__main__":
    sys.exit(main())
