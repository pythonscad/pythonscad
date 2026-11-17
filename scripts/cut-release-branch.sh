#!/usr/bin/env bash
# Create and push a maintenance branch for an already-tagged minor release.
#
# Usage:
#   ./scripts/cut-release-branch.sh 1.2.0
#   ./scripts/cut-release-branch.sh v1.2.0 --dry-run
#
# Creates release/X.Y from tag vX.Y.Z so patch releases (X.Y.Z+1, ...) can
# ship while master continues toward the next minor. See VERSIONING.md.
set -euo pipefail

usage() {
  echo "usage: $0 <version> [--dry-run] [--remote <name>]" >&2
  exit 1
}

if [ $# -lt 1 ]; then
  usage
fi

raw="$1"
shift
dry_run=0
remote="origin"

while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run)
      dry_run=1
      shift
      ;;
    --remote)
      remote="${2:?--remote requires a value}"
      shift 2
      ;;
    *)
      usage
      ;;
  esac
done

version="${raw#v}"
if ! [[ "$version" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
  echo "error: version must be MAJOR.MINOR.PATCH (got: ${raw})" >&2
  exit 1
fi

tag="v${version}"
# shellcheck disable=SC2001
branch="release/$(echo "$version" | sed 's/\.[0-9]*$//')"

if ! git rev-parse -q --verify "refs/tags/${tag}" >/dev/null; then
  echo "error: tag ${tag} not found locally; fetch tags first" >&2
  exit 1
fi

if git show-ref --verify --quiet "refs/heads/${branch}" \
  || git show-ref --verify --quiet "refs/remotes/${remote}/${branch}"; then
  echo "error: branch ${branch} already exists" >&2
  exit 1
fi

echo "Creating ${branch} from ${tag}"
if [ "$dry_run" -eq 1 ]; then
  echo "(dry-run) git branch ${branch} ${tag}"
  echo "(dry-run) git push -u ${remote} ${branch}"
  exit 0
fi

git branch "${branch}" "${tag}"
git push -u "${remote}" "${branch}"
echo "Pushed ${branch}. Release Please will open patch release PRs against it."
