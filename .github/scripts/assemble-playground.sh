#!/usr/bin/env bash
# Assemble the /playground/ subtree from a downloaded WASM web bundle artifact.
set -euo pipefail

bundle_dir="${1:?WASM bundle directory required}"
playground_dir="${2:-playground}"
fonts_src="${3:-web/docs/fonts/IntelOneMono/IntelOneMono-Regular.woff2}"

mkdir -p "${playground_dir}/vendor" "${playground_dir}/fonts"

cp "${bundle_dir}/pythonscad.js" \
   "${bundle_dir}/pythonscad.wasm" \
   "${bundle_dir}/pythonscad.data" \
   "${bundle_dir}/notebook-geometry.mjs" \
   "${bundle_dir}/notebook-protocol.mjs" \
   "${bundle_dir}/notebook-kernel.mjs" \
   "${bundle_dir}/notebook-worker.mjs" \
   "${playground_dir}/"

cp "${bundle_dir}/vendor/three.module.min.js" \
   "${bundle_dir}/vendor/three.core.min.js" \
   "${playground_dir}/vendor/"

cp "${bundle_dir}/notebook.html" "${playground_dir}/index.html"
cp "${fonts_src}" "${playground_dir}/fonts/"
cp wasm-test/playground.htaccess "${playground_dir}/.htaccess"

ls -la "${playground_dir}" "${playground_dir}/vendor" "${playground_dir}/fonts"
