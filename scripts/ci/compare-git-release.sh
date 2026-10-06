#!/usr/bin/env bash
# Compare a delayed upstream release tag with the submodule HEAD.
#
# Used by .github/workflows/update-submodule-manifold.yml. Tag comparison
# peels annotated tags to the underlying commit SHA (`git rev-list -n1`,
# equivalent to `tag^{commit}` but without caret/brace syntax that cmd.exe
# and Git-for-Windows wrappers mangle) so a submodule already at that
# commit is not treated as a new release.
#
# Usage:
#   compare-git-release.sh --git-dir DIR --tag TAG --current-commit SHA
#                          [--output FILE]
#
# Writes GitHub Actions output keys (new_release_available, latest_commit,
# release_tag) to FILE when --output is set, otherwise to stdout.

set -euo pipefail

usage() {
  echo "usage: $0 --git-dir DIR --tag TAG --current-commit SHA [--output FILE]" >&2
  exit 2
}

GIT_DIR=""
TAG=""
CURRENT=""
OUTPUT=""

while [ $# -gt 0 ]; do
  case "$1" in
    --git-dir)
      [ $# -ge 2 ] || usage
      GIT_DIR="$2"
      shift 2
      ;;
    --tag)
      [ $# -ge 2 ] || usage
      TAG="$2"
      shift 2
      ;;
    --current-commit)
      [ $# -ge 2 ] || usage
      CURRENT="$2"
      shift 2
      ;;
    --output)
      [ $# -ge 2 ] || usage
      OUTPUT="$2"
      shift 2
      ;;
    *)
      usage
      ;;
  esac
done

emit() {
  if [ -n "$OUTPUT" ]; then
    echo "$1" >> "$OUTPUT"
  else
    echo "$1"
  fi
}

if [ -z "$TAG" ] || [ "$TAG" = "null" ]; then
  echo "No delayed release found."
  emit "new_release_available=false"
  exit 0
fi

if [ -z "$GIT_DIR" ] || [ -z "$CURRENT" ]; then
  usage
fi

LATEST_COMMIT=$(git -C "$GIT_DIR" rev-list -n 1 "$TAG")
emit "latest_commit=${LATEST_COMMIT}"

if [ "$LATEST_COMMIT" != "$CURRENT" ]; then
  echo "New release detected: ${TAG} (${LATEST_COMMIT}) (current: ${CURRENT})"
  emit "new_release_available=true"
  emit "release_tag=${TAG}"
else
  echo "Already at ${TAG} (${CURRENT})."
  emit "new_release_available=false"
fi
