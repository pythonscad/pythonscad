#!/usr/bin/env bash
#
# Runs a command in the OpenSCAD Base Wasm Docker image for Emscripten builds.
# (mounts $PWD as readwrite and sets up ccache)
#
# Example usage:
#
#   ./scripts/wasm-base-docker-run.sh emcmake cmake -B build-node -DWASM_BUILD_TYPE=node -DCMAKE_BUILD_TYPE=Release -DEXPERIMENTAL=1
#   ./scripts/wasm-base-docker-run.sh cmake --build build-node -j
#   build-node/openscad.js -h
#
# If docker fails because of a platform mismatch, e.g. on Silicon Macs
# (build env currently only built for linux/amd64), register QEMU with Docker using:
#
#   docker run --privileged --rm tonistiigi/binfmt --install all
#
# Also see:
# - https://github.com/openscad/openscad-wasm
#   For a barebones setup
# - https://github.com/openscad/openscad-playground
#   For a full-fledged example
#
set -euo pipefail

CCACHE_DIR=${CCACHE_DIR:-$HOME/.ccache/}
mkdir -p "$CCACHE_DIR"

# Prefer the immutable public GHCR image for this exact toolchain input set.
# A source checkout with unpublished toolchain changes falls back to a local
# build, preserving the ability to develop and test upgrades before merging.
TOOLCHAIN_IMAGE=${PYTHONSCAD_WASM_TOOLCHAIN_IMAGE:-ghcr.io/pythonscad/wasm-python-base}
TOOLCHAIN_HASH=$(./scripts/wasm-toolchain-id.sh)
TOOLCHAIN_REF="${TOOLCHAIN_IMAGE}:toolchain-v1-${TOOLCHAIN_HASH}"
LOCAL_TOOLCHAIN_TAG=pythonscad-wasm-python-base:local

manifest_file=$(mktemp)
trap 'rm -f "$manifest_file"' EXIT
if docker buildx imagetools inspect "$TOOLCHAIN_REF" \
  --raw > "$manifest_file" 2>/dev/null; then
  digest="sha256:$(sha256sum "$manifest_file" | cut -d' ' -f1)"
  immutable_ref="${TOOLCHAIN_REF}@${digest}"
  echo "Pulling ${immutable_ref}..."
  docker pull --platform=linux/amd64 "$immutable_ref"
  docker tag "$immutable_ref" "$TOOLCHAIN_REF"
elif ! docker image inspect "$TOOLCHAIN_REF" &>/dev/null; then
  echo "No published or cached image exists for ${TOOLCHAIN_REF}."
  if ! docker pull --platform=linux/amd64 "$TOOLCHAIN_REF"; then
    echo "No published image exists; building the WASM toolchain locally..."
    docker build \
      --platform=linux/amd64 \
      -f docker/wasm/sysroot.dockerfile \
      --target wasm-python-base \
      -t "$TOOLCHAIN_REF" \
      .
  fi
else
  echo "Registry unavailable; using the cached immutable toolchain image." >&2
fi
docker tag "$TOOLCHAIN_REF" "$LOCAL_TOOLCHAIN_TAG"

echo "
  FROM ${LOCAL_TOOLCHAIN_TAG}
  RUN apt update && \
      apt install -y ccache && \
      apt clean
" | docker build \
  --platform=linux/amd64 \
  -t pythonscad-wasm-ccache:local \
  -f - .

docker run --rm -i \
  --platform=linux/amd64 \
  -w /src \
  -v "$PWD:/src:rw" \
  -v "$CCACHE_DIR:/root/.ccache:rw" \
  pythonscad-wasm-ccache:local \
  "$@"
