"""Smoke test for freshly-built PyQt6 bindings, used by the
pyqt6-bindings CI workflow across Linux/macOS/Windows."""
import sys

if len(sys.argv) != 2:
    print("Usage: pyqt6_smoke_test.py <path-to-bindings-dir>")
    sys.exit(1)

sys.path.insert(0, sys.argv[1])

from PyQt6 import sip, QtWidgets  # noqa: E402
from PyQt6.QtCore import QT_VERSION_STR  # noqa: E402

app = QtWidgets.QApplication([])
widget = QtWidgets.QWidget()
print(f"OK: PyQt6 importiert, QApplication erzeugt, Qt {QT_VERSION_STR}")
