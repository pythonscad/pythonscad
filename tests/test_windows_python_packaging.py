"""Host-side checks for Windows Python packaging, without a Windows toolchain."""

import os
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class WindowsPythonPackagingTests(unittest.TestCase):
    def test_import_library_uses_shipped_dll(self):
        with tempfile.TemporaryDirectory(prefix="python packaging ") as directory:
            root = Path(directory)
            tools = root / "tools"
            tools.mkdir()
            # Model the external tools' file/argument contract, including paths with spaces.
            (tools / "gendef").write_text(
                '#!/bin/sh\nset -eu\ntest -f "$1"\n'
                'printf "EXPORTS\\nPy_Initialize\\n" > "${1%.dll}.def"\n'
            )
            (tools / "dlltool").write_text(
                '#!/bin/sh\nset -eu\n'
                'test "$1" = --dllname\ntest "$3" = --def\n'
                'test "$5" = --output-lib\n'
                'printf "%s\\n" "$2" > "$6"\ncat "$4" >> "$6"\n'
            )
            for tool in tools.iterdir():
                tool.chmod(0o755)
            dll = root / "python314.dll"
            dll.write_bytes(b"DLL fixture")
            archive = root / "output directory" / "libpython3.14.dll.a"
            subprocess.run(
                ["sh", str(ROOT / "scripts/libpython_patch.sh"), str(dll), str(archive)],
                env={**os.environ, "PATH": str(tools) + os.pathsep + os.environ["PATH"]},
                cwd=root,
                check=True,
            )
            self.assertEqual(archive.read_text(), "python314.dll\nEXPORTS\nPy_Initialize\n")
            self.assertEqual(dll.read_bytes(), b"DLL fixture")

    def test_relative_pyqt_destination_survives_build_directory_change(self):
        script = (ROOT / "scripts/ci/build-pyqt6-bindings-windows.sh").read_text()
        setup = script.split('root="$(cd ', 1)[0]
        with tempfile.TemporaryDirectory(prefix="pyqt output ") as directory:
            root = Path(directory)
            build = root / "build"
            build.mkdir()
            source = root / "module.pyd"
            source.write_bytes(b"extension fixture")
            subprocess.run(
                ["bash", "-c", setup + '\ncd "$2"\ncp "$3" "$output_parent/QtCore.pyd"',
                 "pyqt-setup", "libraries/python/PyQt6", str(build), str(source)],
                cwd=root,
                check=True,
            )
            self.assertEqual(
                (root / "libraries/python/PyQt6/QtCore.pyd").read_bytes(),
                b"extension fixture",
            )

    def test_standalone_workflow_generates_archive_after_preparing_headers(self):
        workflow = (ROOT / ".github/workflows/build-pyqt6-bindings.yml").read_text()
        block = workflow.split('          "$GITHUB_WORKSPACE/scripts/prepare-msys2-python.sh"', 1)[1]
        block = block.split("          python -m venv", 1)[0]
        with tempfile.TemporaryDirectory(prefix="bindings workflow ") as directory:
            root = Path(directory)
            scripts = root / "scripts"
            scripts.mkdir()
            dll = root / "embed" / "python314.dll"
            dll.parent.mkdir()
            dll.write_bytes(b"DLL fixture")
            generator = scripts / "libpython_patch.sh"
            generator.write_text(
                '#!/bin/sh\nset -eu\ntest -f "$1"\n'
                'mkdir -p "$(dirname "$2")"\ncp "$1" "$2"\n'
            )
            generator.chmod(0o755)
            # Stop before copying to the runner's /ucrt64 tree.
            generation = block.split('          cp "$GITHUB_WORKSPACE', 1)[0]
            subprocess.run(
                ["bash", "-ec", ': ' + generation],
                env={**os.environ, "GITHUB_WORKSPACE": str(root), "PYTHON_VERSION": "3.14",
                     "DEST_PATH": str(dll.parent)},
                check=True,
            )
            archive = root / "python_mingw/ucrt64/lib/libpython3.14.dll.a"
            self.assertEqual(archive.read_bytes(), b"DLL fixture")


if __name__ == "__main__":
    unittest.main()
