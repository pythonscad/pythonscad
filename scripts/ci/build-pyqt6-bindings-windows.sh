#!/usr/bin/env bash
# Compile PyQt6 against the MSYS2 UCRT64 Qt and the python.org embeddable
# runtime. Only the Qt modules the application imports are built.
#
# Usage: build-pyqt6-bindings-windows.sh <output-parent>
# The PyQt6 package is written to <output-parent>/PyQt6.
set -euo pipefail

output_parent="${1:?Usage: $0 <output-parent>}"
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

PYTHON_VERSION=$(python --version | cut -d' ' -f2 | cut -d'.' -f1,2)
import_lib="$root/python_mingw/ucrt64/lib/libpython${PYTHON_VERSION}.dll.a"
if [[ ! -f "$import_lib" ]]; then
  echo "ERROR: $import_lib is missing. Build PythonSCAD first so dlltool can generate it." >&2
  exit 1
fi
cp "$import_lib" "/ucrt64/lib/libpython${PYTHON_VERSION}.dll.a"

python -m venv /tmp/pyqt-build-venv
# shellcheck disable=SC1091
source /tmp/pyqt-build-venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r "$root/scripts/ci/pyqt6-build-requirements.txt"
python3 "$root/scripts/ci/patch_pyqtbuild_win32_pythonlib.py"

QMAKE=$(which qmake6)
if [[ "$QMAKE" != *.exe ]]; then
  QMAKE="${QMAKE}.exe"
fi
QT_VERSION=$("$QMAKE" -query QT_VERSION)
echo "Detected Qt version: $QT_VERSION"
MAJOR_MINOR=$(echo "$QT_VERSION" | cut -d. -f1,2)

workdir=$(mktemp -d)
cd "$workdir"
RESOLVED=$(python3 "$root/scripts/ci/resolve_pyqt6_version.py" "$MAJOR_MINOR")
read -r PYQT_VERSION PYQT_SDIST_URL <<< "$RESOLVED"
echo "Selected PyQt6 version: $PYQT_VERSION"
curl -fsSL -o PyQt6.tar.gz "$PYQT_SDIST_URL"
tar xzf PyQt6.tar.gz
extracted=$(tar tzf PyQt6.tar.gz | sed -n '1p' | cut -d/ -f1)
cd "$extracted"

python3 "$root/scripts/ci/patch_pyqt6_free_operators.py" sip/
python3 "$root/scripts/ci/patch_pyqt6_removed_qt_apis.py" sip/ "$QT_VERSION"
python3 "$root/scripts/ci/patch_pyqt6_python315.py" .

mkdir -p /tmp/pyqt-staging
win_staging=$(cygpath -w /tmp/pyqt-staging)
# QtCore, QtGui and QtWidgets are what forms and the customizer use.
# QtSvg covers color and icon helpers. Leaving the rest out keeps a
# cache miss closer to ten minutes than forty.
sip-build --qmake="$QMAKE" --confirm-license --verbose \
  --enable QtCore --enable QtGui --enable QtWidgets --enable QtSvg \
  --target-dir="$win_staging"
cd build
mingw32-make -j"$(nproc)"
mingw32-make install

found=$(find /tmp/pyqt-staging -maxdepth 12 -type f \( -iname "QtWidgets*.pyd" -o -iname "QtWidgets*.dll" \) -exec dirname {} \; | head -1)
if [[ -z "$found" ]]; then
  echo "ERROR: could not find installed PyQt6 modules" >&2
  find /tmp/pyqt-staging -maxdepth 10 >&2
  exit 1
fi

mkdir -p "$output_parent"
rm -rf "$output_parent/PyQt6"
cp -a "$found" "$output_parent/PyQt6"

sip_only=$(mktemp -d)
python -m pip install PyQt6-sip \
  --constraint="$root/scripts/ci/pyqt6-build-requirements.txt" \
  --target="$sip_only" --no-deps
sip_pyd=$(find "$sip_only" -maxdepth 6 \( -iname "sip*.so" -o -iname "sip*.pyd" \) | head -1)
if [[ -z "$sip_pyd" ]]; then
  echo "ERROR: PyQt6-sip extension not found" >&2
  exit 1
fi
cp "$sip_pyd" "$output_parent/PyQt6/sip.pyd"

qtcore=$(find "$output_parent/PyQt6" -maxdepth 2 -iname "QtCore*.pyd" | head -1)
expected_dll="python$(echo "$PYTHON_VERSION" | tr -d '.').dll"
if ! objdump -p "$qtcore" | grep -qi "$expected_dll"; then
  echo "ERROR: $qtcore does not link to $expected_dll" >&2
  objdump -p "$qtcore" | grep -i python >&2 || true
  exit 1
fi
echo "PyQt6 $PYQT_VERSION installed to $output_parent/PyQt6"
