#!/usr/bin/env bash
# Post-deploy HTTP smoke checks for the mkdocs website root only.
# Playground checks live in the WASM workflow because the two deploys are independent.
set -euo pipefail

base_url="${1:?base URL required (e.g. https://test.pythonscad.org/)}"
expect_test_banner="${2:-false}"

homepage="${base_url%/}/"
body="$(mktemp)"
trap 'rm -f "$body"' EXIT

echo "Checking homepage ${homepage}"
curl -fsSL "$homepage" -o "$body"
grep -q '<html' "$body"

if [ "$expect_test_banner" = "true" ]; then
  grep -q 'test-site-banner' "$body"
  grep -q 'tracks' "$body"
  grep -q 'master' "$body"
  echo "Test-site banner present."
else
  if grep -q 'test-site-banner' "$body"; then
    echo "Error: production homepage unexpectedly contains the test-site banner." >&2
    exit 1
  fi
  echo "Production homepage has no test-site banner."
fi

echo "Homepage smoke check passed."
