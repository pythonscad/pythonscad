#!/usr/bin/env python3
#
# Regression test driver for fixtures that produce one or more output files
# via in-script PythonSCAD `export()` calls (including helpers like
# `MultiToolExporter`).
#
# The driver:
#
#   1. Creates a per-test scratch directory (`output/<testname>/<basename>/`).
#   2. Runs the PythonSCAD binary with the fixture as a script, with that
#      scratch directory as CWD, so any `export("foo.stl")` call lands there.
#      A throwaway `-o _cli_driver_dummy.echo` is supplied to force CLI /
#      headless mode; without it, PythonSCAD treats the script as "open in
#      GUI" and hits the single-instance lock if a desktop pythonscad is
#      already running. The dummy is not inspected.
#   3. Auto-discovers every produced file in the scratch directory matching
#      `*.<suffix>`.
#   4. Applies format-aware post-processing (header progname rewrite for
#      STL/SVG/OBJ, inner-XML extraction for 3MF), then compares each
#      produced file against
#      `tests/regression/<testname>/<basename>/<filename>` using
#      `test_cmdline_tool.compare_default()` -- a normalized text
#      comparison (line-ending normalization + unified diff), not a raw
#      bytes-equality check. For the ASCII / text-derived formats that the
#      post-processors normalize, this is effectively bytes-equality;
#      true binary outputs (binary STL, AMF, ...) would need a separate
#      bytes-equality branch added to `_post_process` / the comparison
#      step. The expected directory mirrors the actual directory layout
#      one-for-one, so missing or unexpected files are caught by simple
#      set diffing.
#
# When TEST_GENERATE=1 is set in the environment (or -g/--generate is passed),
# the produced files are copied into place as the new goldens instead of
# being compared. Any stale goldens with the same suffix that the fixture no
# longer writes are removed so the expected directory always mirrors the run
# output.
#
# Usage:
#   test_export_files.py
#       --pythonscad <pythonscad-binary>
#       --testname <test-group-name>
#       --basename <fixture-basename>
#       --suffix <stl|3mf|svg|obj|...>
#       [--regressiondir <dir>]
#       [--generate]
#       <fixture.py> [<extra-args-for-pythonscad>...]
#
# Exit codes match test_cmdline_tool.py: 0 pass, 1 failure, 2 invalid args.

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import test_cmdline_tool as tct  # noqa: E402  reuse helpers and globals


def _setup_tct_options():
    """Populate the globals that ``tct.compare_default`` reads from."""
    tct.options = tct.Options()
    tct.options.exclude_line_re = None
    tct.options.exclude_debug = False


def _post_process(filename, suffix):
    """Apply the same normalization test_cmdline_tool.py uses post-export."""
    if suffix in ("stl", "svg", "obj"):
        tct.post_process_progname(filename)
    elif suffix == "3mf":
        tct.post_process_3mf(filename)


def _run_pythonscad(pythonscad, fixture, extra_args, rundir):
    """Run the binary inside ``rundir`` so in-script ``export()`` writes there.

    A dummy ``-o`` is supplied to force CLI/headless mode; without it,
    PythonSCAD treats the script as "open this file in the GUI" and hits
    the single-instance lock if a desktop pythonscad is already running.
    The dummy output file is not inspected; the fixture's own in-script
    ``export()`` calls produce the artifacts we compare.
    """
    # Bare filename: resolves relative to pythonscad's CWD (= rundir).
    cli_dummy_name = "_cli_driver_dummy.echo"
    cmdline = [
        pythonscad,
        "--trust-python",
        "--enable=predictible-output",
        "--render",
        "-o", cli_dummy_name,
        os.path.abspath(fixture),
    ] + list(extra_args)

    fontdir = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "data/ttf"))
    env = os.environ.copy()
    env["OPENSCAD_FONT_PATH"] = fontdir

    print("export-files run cmdline:", " ".join(cmdline))
    print("export-files run cwd:", str(rundir))
    sys.stdout.flush()
    sys.stderr.flush()

    proc = subprocess.run(
        cmdline,
        cwd=str(rundir),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if proc.stdout:
        sys.stderr.write(proc.stdout.decode("utf-8", "replace"))
    if proc.stderr:
        sys.stderr.write(proc.stderr.decode("utf-8", "replace"))
    if proc.returncode != 0:
        print(
            f"Error: pythonscad failed with return code {proc.returncode}",
            file=sys.stderr,
        )
        return False
    return True


def _discover(directory, suffix):
    if not directory.is_dir():
        return []
    return sorted(p for p in directory.glob(f"*.{suffix}") if p.is_file())


def main():
    parser = argparse.ArgumentParser(
        description="Run a PythonSCAD fixture that writes files via in-script "
                    "export() and compare the produced set against goldens.")
    parser.add_argument("--pythonscad", required=True,
                        help="Path to the pythonscad executable.")
    parser.add_argument("--testname", required=True,
                        help="ctest group name; also the regression subdir.")
    parser.add_argument("--basename", required=True,
                        help="Fixture basename (without extension).")
    parser.add_argument("--suffix", required=True,
                        help="Output suffix without leading dot, e.g. 'stl'.")
    parser.add_argument(
        "--regressiondir",
        default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "regression"),
        help="Root regression directory (defaults to tests/regression).")
    parser.add_argument("-g", "--generate", action="store_true",
                        help="Generate goldens instead of comparing. Also "
                             "honored via TEST_GENERATE=1 in the environment.")
    parser.add_argument("fixture", help="Python fixture script to run.")
    parser.add_argument("extra_args", nargs=argparse.REMAINDER,
                        help="Extra arguments forwarded to pythonscad.")
    args = parser.parse_args()

    suffix = args.suffix.lstrip(".")
    if not suffix:
        print("Error: --suffix must be non-empty", file=sys.stderr)
        return 2

    generate = args.generate or bool(os.getenv("TEST_GENERATE"))

    rundir = Path("output") / args.testname / args.basename
    rundir.mkdir(parents=True, exist_ok=True)
    for stale in rundir.glob(f"*.{suffix}"):
        if stale.is_file():
            stale.unlink()

    expecteddir = Path(args.regressiondir) / args.testname / args.basename
    if generate:
        expecteddir.mkdir(parents=True, exist_ok=True)

    if not _run_pythonscad(
            args.pythonscad, args.fixture, args.extra_args, rundir):
        return 1

    actual_paths = _discover(rundir, suffix)
    actual_names = {p.name for p in actual_paths}

    if not actual_paths:
        print(
            f"Error: fixture produced no .{suffix} files in {rundir}",
            file=sys.stderr,
        )
        return 1

    if generate:
        # Drop stale goldens with this suffix that the fixture no longer
        # writes, so the expected directory always mirrors the run output.
        for stale in _discover(expecteddir, suffix):
            if stale.name not in actual_names:
                print(f"removing stale golden: {stale}", file=sys.stderr)
                stale.unlink()
        for produced in actual_paths:
            _post_process(str(produced), suffix)
            dst = expecteddir / produced.name
            shutil.copyfile(str(produced), str(dst))
            print(f"generated golden: {dst}", file=sys.stderr)
        return 0

    expected_paths = _discover(expecteddir, suffix)
    if not expected_paths:
        print(
            f"Error: no .{suffix} goldens in {expecteddir}; regenerate with "
            f"TEST_GENERATE=1.",
            file=sys.stderr,
        )
        return 1
    expected_names = {p.name for p in expected_paths}

    missing = expected_names - actual_names
    extra = actual_names - expected_names
    ok = True
    if missing:
        print(
            f"Error: fixture failed to produce expected file(s): "
            f"{sorted(missing)}",
            file=sys.stderr,
        )
        ok = False
    if extra:
        print(
            f"Error: fixture produced unexpected file(s) without goldens: "
            f"{sorted(extra)}",
            file=sys.stderr,
        )
        ok = False

    _setup_tct_options()
    for produced in actual_paths:
        if produced.name not in expected_names:
            continue
        _post_process(str(produced), suffix)
        expected = expecteddir / produced.name
        tct.expectedfilename = str(expected)
        if not tct.compare_default(str(produced)):
            ok = False

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
