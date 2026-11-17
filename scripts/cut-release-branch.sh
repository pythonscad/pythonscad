#!/usr/bin/env bash
# Create and push a maintenance branch for an already-tagged minor release.
#
# Usage:
#   ./scripts/cut-release-branch.sh 1.2.0
#   ./scripts/cut-release-branch.sh v1.2.0 --dry-run
#
# Creates release/X.Y from tag vX.Y.Z so patch releases (X.Y.Z+1, ...) can
# ship while master continues toward the next minor. See VERSIONING.md.
#
# If the tag predates the maintenance-branch Release Please / deploy guards,
# the script commits those files from <remote>/master onto the new branch
# before pushing (a no-version-bump `ci:` commit).
set -euo pipefail

usage() {
  echo "usage: $0 <version> [--dry-run] [--remote <name>]" >&2
  exit 1
}

# Return 0 when treeish is missing the release/** Release Please setup.
needs_release_automation_bootstrap() {
  local treeish="$1"
  if ! git cat-file -e "${treeish}:.release-please-config.release.json" 2>/dev/null; then
    return 0
  fi
  if ! git show "${treeish}:.github/workflows/release-please.yml" 2>/dev/null \
    | grep -q 'release/\*\*'; then
    return 0
  fi
  return 1
}

# Files copied from master onto older tags so Release Please, deploy guards,
# and the PR helper test job all work on the new maintenance branch.
bootstrap_paths=(
  .github/workflows/release-please.yml
  .release-please-config.json
  .release-please-config.release.json
  .github/scripts/is-latest-semver-tag.sh
  .github/scripts/is-latest-semver-tag.test.py
  .github/workflows/deploy-website.yml
  .github/workflows/build-wasm.yml
  scripts/cut-release-branch.sh
  scripts/cut-release-branch.test.py
)

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

# Switching branches / committing bootstrap files must not drag along
# unrelated local edits (staged or unstaged tracked changes).
if ! git diff-index --quiet HEAD --; then
  echo "error: working tree has uncommitted changes; commit or stash first" >&2
  git diff-index --name-status HEAD -- >&2 || true
  exit 1
fi

bootstrap_ref="${remote}/master"
if ! git rev-parse -q --verify "refs/remotes/${bootstrap_ref}" >/dev/null \
  && ! git rev-parse -q --verify "refs/heads/master" >/dev/null; then
  echo "error: need ${remote}/master (or local master) to bootstrap older tags" >&2
  exit 1
fi
if ! git rev-parse -q --verify "refs/remotes/${bootstrap_ref}" >/dev/null; then
  bootstrap_ref="master"
fi

do_bootstrap=0
if needs_release_automation_bootstrap "${tag}"; then
  do_bootstrap=1
fi

if [ "$do_bootstrap" -eq 1 ]; then
  # Untracked files at bootstrap paths survive branch switches and would be
  # clobbered by `git checkout <ref> -- <path>`.
  collisions=0
  for path in "${bootstrap_paths[@]}"; do
    if [ -e "$path" ] || [ -L "$path" ]; then
      if ! git ls-files --error-unmatch -- "$path" >/dev/null 2>&1; then
        echo "error: untracked file would be overwritten by bootstrap: ${path}" >&2
        collisions=1
      fi
    fi
  done
  if [ "$collisions" -ne 0 ]; then
    exit 1
  fi
fi

echo "Creating ${branch} from ${tag}"
if [ "$do_bootstrap" -eq 1 ]; then
  echo "Tag ${tag} predates maintenance-branch automation; will bootstrap from ${bootstrap_ref}"
fi

if [ "$dry_run" -eq 1 ]; then
  echo "(dry-run) git branch ${branch} ${tag}"
  if [ "$do_bootstrap" -eq 1 ]; then
    echo "(dry-run) checkout release-please + deploy-guard files from ${bootstrap_ref}"
    echo "(dry-run) git commit -m 'ci: bootstrap release maintenance-branch automation'"
  fi
  echo "(dry-run) git push -u ${remote} ${branch}"
  exit 0
fi

start_ref="$(git symbolic-ref -q --short HEAD || git rev-parse HEAD)"

git branch "${branch}" "${tag}"
git checkout "${branch}"

if [ "$do_bootstrap" -eq 1 ]; then
  missing=0
  for path in "${bootstrap_paths[@]}"; do
    if ! git cat-file -e "${bootstrap_ref}:${path}" 2>/dev/null; then
      echo "error: ${bootstrap_ref} is missing ${path}; merge the maintenance-branch CI PR to master first" >&2
      missing=1
    fi
  done
  if [ "$missing" -ne 0 ]; then
    git checkout "${start_ref}"
    git branch -D "${branch}"
    exit 1
  fi
  git checkout "${bootstrap_ref}" -- "${bootstrap_paths[@]}"
  # --only commits just these paths even if something else got staged.
  # Options must come before `--` / pathspecs.
  git commit --only -m "$(cat <<'EOF'
ci: bootstrap release maintenance-branch automation

EOF
)" -- "${bootstrap_paths[@]}"
fi

git push -u "${remote}" "${branch}"
git checkout "${start_ref}"
echo "Pushed ${branch}. Release Please will open patch release PRs against it."
