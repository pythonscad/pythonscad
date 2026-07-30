#!/usr/bin/env bash
# Build the mkdocs site for production or test deployment targets.
set -euo pipefail

target="${1:?production or test required}"
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
config="${repo_root}/web/mkdocs.yml"
build_config="${config}"

cleanup() {
  rm -f "${repo_root}/web/.mkdocs-test.build.yml"
}
trap cleanup EXIT

case "$target" in
  production)
    ;;
  test)
    build_config="${repo_root}/web/.mkdocs-test.build.yml"
    sed -e 's|^site_url:.*|site_url: https://test.pythonscad.org/|' \
        -e 's|^  test_site: false|  test_site: true|' \
        "$config" > "$build_config"
    ;;
  *)
    echo "Error: unknown mkdocs deploy target '${target}' (expected production or test)." >&2
    exit 1
    ;;
esac

cd "${repo_root}/web"
uv run --project .. --no-sync mkdocs build --strict -f "$build_config"

# shellcheck disable=SC1091
source "${repo_root}/.github/scripts/web-deploy-common.sh"
validate_website_build "${repo_root}/web/site"
