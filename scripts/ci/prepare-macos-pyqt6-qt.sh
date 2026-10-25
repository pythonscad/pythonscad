#!/usr/bin/env bash
# Build macOS against the Qt frameworks shipped inside the PyQt6 wheel.
#
# The wheel contains runtime frameworks only (no headers, moc, or
# macdeployqt). This script installs the matching official Qt SDK for
# those tools and headers, then replaces the SDK's framework binaries
# and plugins with the copies from the wheel, so the link and the
# deployed app use the same Qt the bindings import. QScintilla is built
# from source against that prefix. Homebrew is not asked for Qt.
set -euo pipefail

if [[ "$(uname -s)" != "Darwin" ]]; then
  echo "prepare-macos-pyqt6-qt.sh runs on macOS only" >&2
  exit 1
fi

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
dest="$root/third_party/macos-pyqt6-qt"
versions_file="$root/scripts/ci/pyqt6-versions.json"

pyqt_version=$(python3 - "$versions_file" <<'PY'
import json, sys
data = json.load(open(sys.argv[1]))
def key(version):
    return tuple(int(part) for part in version.split("."))
print(max(data.values(), key=key))
PY
)

if [[ -f "$dest/pyqt-version" && "$(cat "$dest/pyqt-version")" == "$pyqt_version" && -f "$dest/prefix" ]]; then
  prefix=$(cat "$dest/prefix")
  if [[ -x "$prefix/bin/qmake" && -d "$root/libraries/python/PyQt6" ]]; then
    echo "macOS wheel Qt already prepared at $prefix (PyQt6 $pyqt_version)"
    exit 0
  fi
fi

rm -rf "$dest"
mkdir -p "$dest/wheels"
python3 -m venv "$dest/venv"
# shellcheck disable=SC1091
source "$dest/venv/bin/activate"
python -m pip install --upgrade pip aqtinstall
python -m pip download "PyQt6==${pyqt_version}" "PyQt6-sip==13.13.0" \
  --only-binary=:all: -d "$dest/wheels"

unpack="$dest/unpacked"
rm -rf "$unpack"
mkdir -p "$unpack"
for wheel in "$dest/wheels"/*.whl; do
  unzip -q -o "$wheel" -d "$unpack"
done

qt_plist="$unpack/PyQt6/Qt6/lib/QtCore.framework/Resources/Info.plist"
if [[ ! -f "$qt_plist" ]]; then
  echo "ERROR: PyQt6-Qt6 wheel did not contain QtCore.framework" >&2
  exit 1
fi
qt_version=$(python3 -c 'import plistlib,sys; print(plistlib.load(open(sys.argv[1],"rb"))["CFBundleVersion"])' "$qt_plist")
echo "PyQt6 $pyqt_version ships Qt $qt_version"

aqt install-qt mac desktop "$qt_version" clang_64 \
  -m qtmultimedia qt5compat \
  -O "$dest/aqt"
prefix=$(find "$dest/aqt" -type f -path '*/bin/qmake' | head -1)
prefix=$(dirname "$(dirname "$prefix")")
if [[ ! -x "$prefix/bin/qmake" ]]; then
  echo "ERROR: qmake not found after installing Qt $qt_version" >&2
  exit 1
fi

wheel_qt="$unpack/PyQt6/Qt6"
# Keep the SDK headers and cmake files. Replace the binaries and plugins
# with the wheel's, which are the ones the bindings were built for.
for framework in "$wheel_qt/lib/"*.framework; do
  name=$(basename "$framework" .framework)
  src="$framework/Versions/A/$name"
  dest_bin="$prefix/lib/$name.framework/Versions/A/$name"
  if [[ -f "$src" && -f "$dest_bin" ]]; then
    cp -f "$src" "$dest_bin"
  fi
done
for extra in "$wheel_qt/lib/"*.dylib; do
  [[ -e "$extra" ]] || continue
  cp -f "$extra" "$prefix/lib/"
done
if [[ -d "$wheel_qt/plugins" ]]; then
  mkdir -p "$prefix/plugins"
  cp -R "$wheel_qt/plugins/." "$prefix/plugins/"
fi

qsci_version=2.14.1
qsci_tar="$dest/QScintilla_src-${qsci_version}.tar.gz"
curl -fsSL -o "$qsci_tar" \
  "https://www.riverbankcomputing.com/static/Downloads/QScintilla/${qsci_version}/QScintilla_src-${qsci_version}.tar.gz"
tar -C "$dest" -xzf "$qsci_tar"
(
  cd "$dest/QScintilla_src-${qsci_version}/src"
  "$prefix/bin/qmake" qscintilla.pro
  make -j"$(sysctl -n hw.ncpu)"
  make install
)

if [[ ! -d "$unpack/PyQt6" ]]; then
  echo "ERROR: PyQt6 bindings wheel did not unpack" >&2
  exit 1
fi
# Drop the wheel's Qt6 runtime tree from the importable package. The app
# links and deploys those frameworks from $prefix instead.
rm -rf "$unpack/PyQt6/Qt6"
rm -rf "$root/libraries/python/PyQt6"
mkdir -p "$root/libraries/python"
cp -R "$unpack/PyQt6" "$root/libraries/python/PyQt6"

printf '%s\n' "$prefix" > "$dest/prefix"
printf '%s\n' "$pyqt_version" > "$dest/pyqt-version"
echo "macOS wheel Qt ready at $prefix"
