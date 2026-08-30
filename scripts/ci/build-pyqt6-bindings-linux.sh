#!/usr/bin/env bash
# Build PyQt6 bindings against the system Qt6 already installed via
# get-dependencies.py, and drop them into libraries/python/PyQt6 so
# that CMake's install() picks them up automatically as part of the
# normal package build (dpkg-buildpackage / rpmbuild). No changes to
# debian/control or pythonscad.spec are needed: dh_shlibdeps and RPM's
# automatic ELF dependency scanning discover the Qt6 .so dependencies
# of the bindings themselves.
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$PROJECT_ROOT"

QMAKE=$(command -v qmake6 || command -v qmake)
if [ -z "$QMAKE" ]; then
  echo "ERROR: qmake6/qmake not found - is qt6-base-dev(-tools) installed?" >&2
  exit 1
fi

QT_VERSION=$("$QMAKE" -query QT_VERSION)
MAJOR_MINOR=$(echo "$QT_VERSION" | cut -d. -f1,2)
echo "Building PyQt6 bindings against system Qt $QT_VERSION"

# Build tooling not covered by pythonscad's own Build-Depends/BuildRequires
if command -v apt-get >/dev/null; then
  apt-get install -y python3-venv python3-pip curl >/dev/null
elif command -v dnf >/dev/null; then
  dnf install -y python3-pip >/dev/null
 command -v curl >/dev/null || dnf install -y curl-minimal >/dev/null
fi

PYTHON_BIN=$(command -v python3.12 || command -v python3.11 || command -v python3.10 || command -v python3)
echo "Using Python: $PYTHON_BIN ($($PYTHON_BIN --version))"

"$PYTHON_BIN" -m venv /tmp/pyqt-build-venv

# shellcheck disable=SC1091
source /tmp/pyqt-build-venv/bin/activate
pip install --upgrade pip
pip install -r "$PROJECT_ROOT/scripts/ci/pyqt6-build-requirements.txt"

mkdir -p /tmp/pyqt-src && cd /tmp/pyqt-src
read -r PYQT_VERSION PYQT_SDIST_URL <<< "$(python3 "$PROJECT_ROOT/scripts/ci/resolve_pyqt6_version.py" "$MAJOR_MINOR")"
echo "Selected PyQt6 version: $PYQT_VERSION"
curl -L -o PyQt6.tar.gz "$PYQT_SDIST_URL"
tar xzf PyQt6.tar.gz
EXTRACTED_DIR=$(tar tzf PyQt6.tar.gz | sed -n '1p' | cut -d/ -f1)
cd "$EXTRACTED_DIR"

python3 "$PROJECT_ROOT/scripts/ci/patch_pyqt6_free_operators.py" sip/
sip-build --qmake="$QMAKE" --confirm-license --verbose --target-dir=/tmp/pyqt-staging
cd build
make -j"$(nproc)"
make INSTALL_ROOT=/tmp/pyqt-install-root install

FOUND_DIR=$(find /tmp/pyqt-install-root -maxdepth 12 -type f -iname "QtWidgets*.so" -exec dirname {} \; | head -1)
if [ -z "$FOUND_DIR" ]; then
  echo "ERROR: could not find installed PyQt6 modules" >&2
  find /tmp/pyqt-install-root -maxdepth 10 >&2
  exit 1
fi

mkdir -p "$PROJECT_ROOT/libraries/python"
rm -rf "$PROJECT_ROOT/libraries/python/PyQt6"
cp -r "$FOUND_DIR" "$PROJECT_ROOT/libraries/python/PyQt6"

mkdir -p /tmp/pyqt-sip-only
pip install PyQt6-sip \
  --constraint="$PROJECT_ROOT/scripts/ci/pyqt6-build-requirements.txt" \
  --target=/tmp/pyqt-sip-only --no-deps
SIP_SO=$(find /tmp/pyqt-sip-only -maxdepth 6 -iname "sip*.so" | head -1)
if [ -z "$SIP_SO" ]; then
  echo "ERROR: PyQt6-sip .so not found" >&2
  exit 1
fi
cp "$SIP_SO" "$PROJECT_ROOT/libraries/python/PyQt6/"

echo "PyQt6 bindings placed in $PROJECT_ROOT/libraries/python/PyQt6 (Qt $QT_VERSION, PyQt $PYQT_VERSION)"
