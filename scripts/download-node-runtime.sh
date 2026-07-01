#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: download-node-runtime.sh --dest DIR [--platform PLATFORM]

Downloads the pinned Node.js runtime from packaging/node-runtime.json,
verifies its SHA256 checksum, and writes a normalized runtime tree:

  POSIX:   DIR/bin/node
  Windows: DIR/node.exe

PLATFORM defaults to the current host and is one of:
  darwin-arm64, darwin-x64, linux-x64, win-x64
EOF
}

die() {
  echo "error: $*" >&2
  exit 1
}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/.." && pwd)"
metadata_file="${project_root}/packaging/node-runtime.json"
dest=""
platform=""

while [ "$#" -gt 0 ]; do
  case "$1" in
    --dest)
      dest="${2:-}"
      shift 2
      ;;
    --platform)
      platform="${2:-}"
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      die "unknown argument: $1"
      ;;
  esac
done

[ -n "$dest" ] || die "--dest is required"
[ -f "$metadata_file" ] || die "metadata file not found: $metadata_file"

if [ -z "$platform" ]; then
  case "$(uname -s)" in
    Darwin)
      case "$(uname -m)" in
        arm64) platform="darwin-arm64" ;;
        x86_64) platform="darwin-x64" ;;
        *) die "unsupported macOS architecture: $(uname -m)" ;;
      esac
      ;;
    Linux)
      case "$(uname -m)" in
        x86_64|amd64) platform="linux-x64" ;;
        aarch64|arm64) platform="linux-arm64" ;;
        *) die "unsupported Linux architecture: $(uname -m)" ;;
      esac
      ;;
    MINGW*|MSYS*|CYGWIN*)
      case "$(uname -m)" in
        x86_64|amd64) platform="win-x64" ;;
        *) die "unsupported Windows architecture: $(uname -m)" ;;
      esac
      ;;
    *)
      die "unsupported host OS: $(uname -s)"
      ;;
  esac
fi

tmpdir="$(mktemp -d)"
cleanup() {
  rm -rf "$tmpdir"
}
trap cleanup EXIT

metadata_tmp="${tmpdir}/metadata.txt"
python3 - "$metadata_file" "$platform" > "$metadata_tmp" <<'PY'
import json
import pathlib
import sys

metadata_path = pathlib.Path(sys.argv[1])
platform = sys.argv[2]
data = json.loads(metadata_path.read_text())
entry = data["platforms"].get(platform)
if entry is None:
    raise SystemExit(f"unsupported platform: {platform}")

print(data["nodeVersion"])
print(entry["url"])
print(entry["sha256"])
print(entry["type"])
print(pathlib.PurePosixPath(entry["url"]).name)
PY

node_version="$(sed -n '1p' "$metadata_tmp")"
url="$(sed -n '2p' "$metadata_tmp")"
expected_sha256="$(sed -n '3p' "$metadata_tmp")"
package_type="$(sed -n '4p' "$metadata_tmp")"
archive_name="$(sed -n '5p' "$metadata_tmp")"

download_path="${tmpdir}/${archive_name}"
echo "Downloading Node.js ${node_version} (${platform})"
curl -fsSL "$url" -o "$download_path"

if command -v sha256sum >/dev/null 2>&1; then
  actual_sha256="$(sha256sum "$download_path" | awk '{print $1}')"
else
  actual_sha256="$(shasum -a 256 "$download_path" | awk '{print $1}')"
fi

if [ "$actual_sha256" != "$expected_sha256" ]; then
  die "checksum mismatch for ${archive_name}: expected ${expected_sha256}, got ${actual_sha256}"
fi

rm -rf "$dest"
mkdir -p "$dest"

case "$package_type" in
  tar)
    mkdir -p "${tmpdir}/extract"
    tar -xf "$download_path" -C "${tmpdir}/extract"
    node_dir="$(find "${tmpdir}/extract" -mindepth 1 -maxdepth 1 -type d | head -1)"
    [ -n "$node_dir" ] || die "extracted Node.js archive did not contain a top-level directory"
    mkdir -p "${dest}/bin"
    cp "${node_dir}/bin/node" "${dest}/bin/node"
    chmod +x "${dest}/bin/node"
    ;;
  exe)
    cp "$download_path" "${dest}/node.exe"
    chmod +x "${dest}/node.exe"
    ;;
  *)
    die "unsupported package type: ${package_type}"
    ;;
esac

echo "Node.js ${node_version} runtime installed to ${dest}"
