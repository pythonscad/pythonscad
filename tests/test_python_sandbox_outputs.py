#!/usr/bin/env python3

import argparse
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


FAKE_WASM_MODULE = r"""
const dirs = new Set(['/']);
const files = new Map();

function normalize(p) {
  return p.replace(/\/+/g, '/').replace(/\/$/, '') || '/';
}

function parent(p) {
  const normalized = normalize(p);
  if (normalized === '/') return '/';
  const index = normalized.lastIndexOf('/');
  return index <= 0 ? '/' : normalized.slice(0, index);
}

function mkdirp(p) {
  const normalized = normalize(p);
  if (normalized === '/') return;
  mkdirp(parent(normalized));
  dirs.add(normalized);
}

function children(dir) {
  const normalized = normalize(dir);
  const prefix = normalized === '/' ? '/' : `${normalized}/`;
  const names = new Set(['.', '..']);
  for (const candidate of dirs) {
    if (candidate === normalized || !candidate.startsWith(prefix)) continue;
    const rest = candidate.slice(prefix.length);
    names.add(rest.split('/')[0]);
  }
  for (const candidate of files.keys()) {
    if (!candidate.startsWith(prefix)) continue;
    const rest = candidate.slice(prefix.length);
    names.add(rest.split('/')[0]);
  }
  return Array.from(names);
}

async function OpenSCAD() {
  const FS = {
    mkdir: mkdirp,
    writeFile(name, data) {
      mkdirp(parent(name));
      files.set(normalize(name), Buffer.isBuffer(data) ? data : Buffer.from(String(data)));
    },
    readFile(name, options) {
      const data = files.get(normalize(name));
      if (!data) throw new Error(`missing file: ${name}`);
      return options && options.encoding === 'utf8' ? data.toString('utf8') : new Uint8Array(data);
    },
    readdir: children,
    stat(name) {
      const normalized = normalize(name);
      if (dirs.has(normalized)) return {mode: 0o040000};
      if (files.has(normalized)) return {mode: 0o100000};
      throw new Error(`missing path: ${name}`);
    },
    isDir(mode) {
      return (mode & 0o170000) === 0o040000;
    },
    isFile(mode) {
      return (mode & 0o170000) === 0o100000;
    },
  };
  return {
    FS,
    callMain(argv) {
      FS.writeFile('/work/output.csg', 'cube(size = [1, 1, 1], center = false);\n');
      FS.writeFile('/work/out/a.stl', 'solid a\nendsolid a\n');
      FS.writeFile('/work/out/nested/b.stl', 'solid b\nendsolid b\n');
    },
  };
}

module.exports = OpenSCAD;
"""


def run(args, cwd, env):
    return subprocess.run(
        args,
        cwd=cwd,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pythonscad", required=True)
    parser.add_argument("--runner", required=True)
    args = parser.parse_args()

    if shutil.which("node") is None:
        print("node was not found; skipping sandbox output test")
        return

    with tempfile.TemporaryDirectory(prefix="pythonscad-sandbox-outputs-") as tmp:
        tmpdir = Path(tmp)
        wasm_dir = tmpdir / "fake-wasm"
        wasm_dir.mkdir()
        (wasm_dir / "pythonscad.js").write_text(FAKE_WASM_MODULE, encoding="utf-8")
        (wasm_dir / "pythonscad.wasm").write_bytes(b"")
        (wasm_dir / "pythonscad.data").write_bytes(b"")

        script = tmpdir / "model.py"
        script.write_text("from pythonscad import *\nshow(cube(1))\n", encoding="utf-8")
        out_csg = tmpdir / "out.csg"
        sandbox_outputs = tmpdir / "sandbox-outputs"
        sandbox_outputs.mkdir()

        env = os.environ.copy()
        env["PYTHONSCAD_WASM_DIR"] = str(wasm_dir)
        env["PYTHONSCAD_SANDBOX_RUNNER"] = args.runner

        first = run(
            [
                args.pythonscad,
                "--python=sandboxed",
                "--sandbox-output-dir",
                str(sandbox_outputs),
                "-o",
                str(out_csg),
                str(script),
            ],
            tmpdir,
            env,
        )
        assert first.returncode == 0, first.stdout
        assert out_csg.read_text(encoding="utf-8").startswith("cube(")
        assert (sandbox_outputs / "a.stl").read_text(encoding="utf-8").startswith("solid a")
        assert (sandbox_outputs / "nested" / "b.stl").read_text(encoding="utf-8").startswith("solid b")

        second = run(
            [
                args.pythonscad,
                "--python=sandboxed",
                "--sandbox-output-dir",
                str(sandbox_outputs),
                "-o",
                str(tmpdir / "out2.csg"),
                str(script),
            ],
            tmpdir,
            env,
        )
        assert second.returncode != 0
        assert "Refusing to overwrite sandbox output file" in second.stdout


if __name__ == "__main__":
    main()
