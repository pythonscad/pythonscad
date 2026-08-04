#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 <platform> <profile> <qmake> <python> <destination>" >&2
  exit 2
fi

PLATFORM=$1
PROFILE=$2
QMAKE=$3
PYTHON=$4
DESTINATION=$5
REPO_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
REPOSITORY=${GITHUB_REPOSITORY:-pythonscad/pythonscad}

QT_VERSION=$("$QMAKE" -query QT_VERSION)
QT_MAJOR_MINOR=$(cut -d. -f1,2 <<< "$QT_VERSION")
PYQT_VERSION=$("$PYTHON" "$REPO_ROOT/scripts/ci/resolve_pyqt6_version.py" \
  "$QT_MAJOR_MINOR" --version-only)
PYTHON_VERSION=$("$PYTHON" -c 'import platform; print(platform.python_version())')
SOABI=$("$PYTHON" -c 'import sysconfig; print(sysconfig.get_config_var("SOABI") or "")')

IDENTITY_ARGS=(
  --repo-root "$REPO_ROOT"
  --platform "$PLATFORM"
  --profile "$PROFILE"
  --qt-version "$QT_VERSION"
  --pyqt-version "$PYQT_VERSION"
  --python-version "$PYTHON_VERSION"
  --soabi "$SOABI"
)

RELEASE_TAG=$("$PYTHON" "$REPO_ROOT/scripts/ci/pyqt6_artifact.py" field \
  "${IDENTITY_ARGS[@]}" --field release_tag)
ASSET_NAME=$("$PYTHON" "$REPO_ROOT/scripts/ci/pyqt6_artifact.py" field \
  "${IDENTITY_ARGS[@]}" --field asset_name)

DOWNLOAD_DIR=$(mktemp -d)
trap 'rm -rf "$DOWNLOAD_DIR"' EXIT

echo "Fetching immutable PyQt6 artifact $RELEASE_TAG/$ASSET_NAME"
if ! gh release download "$RELEASE_TAG" --repo "$REPOSITORY" \
  --pattern "$ASSET_NAME" --dir "$DOWNLOAD_DIR"; then
  echo "No compatible PyQt6 bindings artifact exists for:" >&2
  echo "  platform=$PLATFORM Qt=$QT_VERSION PyQt=$PYQT_VERSION" >&2
  echo "  Python=$PYTHON_VERSION SOABI=$SOABI" >&2
  echo "Run the 'Build PyQt6 Bindings (Qt6)' workflow on the default branch." >&2
  exit 1
fi

rm -rf "$DESTINATION/PyQt6"
unzip -q "$DOWNLOAD_DIR/$ASSET_NAME" -d "$DESTINATION"

"$PYTHON" "$REPO_ROOT/scripts/ci/pyqt6_artifact.py" validate \
  "${IDENTITY_ARGS[@]}" --manifest "$DESTINATION/PYQT6_MANIFEST.json"
