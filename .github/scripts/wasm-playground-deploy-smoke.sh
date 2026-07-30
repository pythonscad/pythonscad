#!/usr/bin/env bash
# Post-deploy HTTP smoke checks for the /playground/ subtree only.
set -euo pipefail

base_url="${1:?base URL required (e.g. https://test.pythonscad.org/)}"

playground="${base_url%/}/playground/"
wasm_url="${playground}pythonscad.wasm"
body="$(mktemp)"
trap 'rm -f "$body"' EXIT

echo "Checking playground landing page ${playground}"
curl -fsSL "$playground" -o "$body"
grep -q '<html' "$body"

echo "Checking WASM MIME type for ${wasm_url}"
content_type="$(
  curl -fsSI "$wasm_url" | tr -d '\r' | awk -F': ' 'tolower($1)=="content-type"{print tolower($2); exit}'
)"
case "$content_type" in
  application/wasm*)
    echo "WASM Content-Type OK: ${content_type}"
    ;;
  *)
    echo "Error: expected application/wasm, got '${content_type:-<missing>}'" >&2
    exit 1
    ;;
esac

echo "Playground smoke check passed."
