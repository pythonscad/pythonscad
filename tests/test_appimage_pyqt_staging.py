"""Exercise AppImage PyQt staging without building an AppImage."""

import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest


class PyQtStagingTests(unittest.TestCase):
    def test_namespace_package_is_importable_from_embedded_search_path(self):
        script = (Path(__file__).resolve().parents[1] / "scripts/create_appimage.sh").read_text()
        staging = script.split('    PYQT_SRC=', 1)[1].split('    info "Bundled PyQt6', 1)[0]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "dist-packages" / "PyQt6"
            source.mkdir(parents=True)
            (source / "QtCore.py").write_text("MARKER = 'bundled'\n")
            destination = root / "app" / "usr" / "lib" / "python3"
            destination.mkdir(parents=True)
            # -S excludes any system-installed PyQt6, preserving namespace semantics.
            command = (
                f'python3() {{ {shlex.quote(sys.executable)} -S "$@"; }}\n'
                'die() { echo "$*" >&2; exit 1; }\n'
                'set -e\n    PYQT_SRC=' + staging
            )
            subprocess.run(
                ["bash", "-c", command],
                env={**os.environ, "PYTHONPATH": str(source.parent),
                     "PYTHON_LIB_DIR": str(destination)},
                check=True,
            )
            subprocess.run(
                [sys.executable, "-S", "-c",
                 "from PyQt6.QtCore import MARKER; assert MARKER == 'bundled'"],
                env={**os.environ, "PYTHONPATH": str(destination)},
                cwd=root,
                check=True,
            )


if __name__ == "__main__":
    unittest.main()
