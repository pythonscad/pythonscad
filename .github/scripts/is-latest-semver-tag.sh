#!/usr/bin/env bash
# Exit 0 if the given tag is the highest vMAJOR.MINOR.PATCH tag in the repo.
#
# Used by website/playground deploy jobs so a patch release on an older
# maintenance line (e.g. v1.2.1 while v1.3.0 exists) does not roll production
# back to that older tree.
#
# Requires a checkout with tags available (fetch-depth: 0 or fetch-tags: true).
set -euo pipefail

raw="${1:?usage: is-latest-semver-tag.sh <tag>}"
this="${raw#v}"

if ! [[ "$this" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
  echo "Not a plain semver tag: ${raw}" >&2
  exit 1
fi

latest="$(
  git tag -l 'v*' \
    | sed 's/^v//' \
    | grep -E '^[0-9]+\.[0-9]+\.[0-9]+$' \
    | sort -V \
    | tail -1
)"

if [ -z "$latest" ]; then
  echo "No semver tags found in this checkout" >&2
  exit 1
fi

if [ "$this" = "$latest" ]; then
  echo "${raw} is the latest semver tag (v${latest})"
  exit 0
fi

echo "Skipping: ${raw} is not the latest semver tag (latest is v${latest})"
exit 1
