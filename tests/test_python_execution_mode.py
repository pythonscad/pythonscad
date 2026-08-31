#!/usr/bin/env python3

import argparse
import subprocess
import tempfile
from pathlib import Path


def run(args, cwd):
    return subprocess.run(
        args,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pythonscad", required=True)
    args = parser.parse_args()

    with tempfile.TemporaryDirectory(prefix="pythonscad-python-mode-") as tmp:
        tmpdir = Path(tmp)
        script = tmpdir / "model.py"
        script.write_text("from pythonscad import *\nshow(cube(1))\n", encoding="utf-8")

        native_output = tmpdir / "native.csg"
        native = run(
            [args.pythonscad, "--python=native", "-o", str(native_output), str(script)],
            tmpdir,
        )
        assert native.returncode == 0, native.stdout
        assert native_output.exists() and native_output.stat().st_size > 0

        invalid = run([args.pythonscad, "--python=bogus", "-o", str(tmpdir / "bad.csg"), str(script)], tmpdir)
        assert invalid.returncode != 0
        assert "--python must be either 'sandboxed' or 'native'" in invalid.stdout

        deprecated = run(
            [args.pythonscad, "--trust-python", "-o", str(tmpdir / "deprecated.csg"), str(script)],
            tmpdir,
        )
        assert "--trust-python is deprecated" in deprecated.stdout
        deprecated_output = tmpdir / "deprecated.csg"
        if deprecated.returncode == 0:
            assert deprecated_output.exists() and deprecated_output.stat().st_size > 0
        else:
            assert "Sandboxed Python" in deprecated.stdout


if __name__ == "__main__":
    main()
