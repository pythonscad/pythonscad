# PythonSCAD Versioning System

This document describes how PythonSCAD manages versions, releases, and changelogs.

## Overview

PythonSCAD uses **automated semantic versioning** based on [Conventional Commits](https://www.conventionalcommits.org/). Version numbers are determined automatically from commit messages, eliminating manual version decisions.

Development happens on `master`. Each published minor also gets a long-lived
`release/X.Y` maintenance branch so patch releases can ship after work on the
next minor has already started on `master`.

## Semantic Versioning

Version format: `MAJOR.MINOR.PATCH` (e.g., `1.2.0`, `1.2.3`)

- **MAJOR**: Breaking changes (incompatible API changes)
- **MINOR**: New features (backward-compatible)
- **PATCH**: Bug fixes (backward-compatible)

### Version bumps (post-1.0.0)

- Breaking changes bump **MAJOR** version: `1.2.3` → `2.0.0`
- New features bump **MINOR** version: `1.2.3` → `1.3.0`
- Bug fixes bump **PATCH** version: `1.2.3` → `1.2.4`

On a `release/X.Y` maintenance branch, Release Please is configured with
`always-bump-patch`, so every releasable commit stays on that minor line
(`1.2.0` → `1.2.1` → `1.2.2`), even if a `feat:` commit is cherry-picked there.

## How Versions Are Determined

Versions are calculated automatically from commit message prefixes:

| Commit Type | Example | Bump on `master` | Bump on `release/X.Y` |
|-------------|---------|------------------|------------------------|
| `feat:` | `feat: add Python 3.12 support` | `1.2.0` → `1.3.0` | patch only |
| `fix:` | `fix: resolve memory leak` | `1.2.0` → `1.2.1` | patch only |
| `feat!:` | `feat!: redesign coordinate system` | `1.2.0` → `2.0.0` | patch only |
| `fix!:` | `fix!: correct broken API behavior` | `1.2.0` → `2.0.0` | patch only |
| Other types | `docs:`, `style:`, `refactor:`, etc. | No version bump | No version bump |

**Breaking Change Footer:**
```
feat: redesign rendering engine

BREAKING CHANGE: The origin point has changed from bottom-left to center.
All existing models need to be adjusted.
```

This also triggers a MAJOR version bump on `master`.

## Branching Model

```text
v1.2.0 ── master ────── feat ────── fix ─────── … ──► v1.3.0
            │                         ▲
            │                         │ cherry-pick
            └─ release/1.2 ── fix ────┘──► v1.2.1
```

| Branch | Purpose | Release Please config |
|--------|---------|------------------------|
| `master` | Next minor / major. Feature and bugfix PRs land here first. | `.release-please-config.json` (default SemVer) |
| `release/X.Y` | Patch line for that minor. Created from `vX.Y.0` (or the latest `vX.Y.Z`). | `.release-please-config.release.json` (`always-bump-patch`) |

Rules:

- **Never merge** a `release/*` branch back into `master`. That would fight
  the next-minor Release PR (`VERSION.txt`, `CHANGELOG.md`, manifest).
- **Cherry-pick** the fix commit onto each maintenance line that still needs it.
- Prefer fixing on `master` first when the bug still exists there, then
  backport. Fix on the oldest line first only when `master` has rewritten the
  code and the patch would not apply.
- OpenSCAD upstream syncs stay on `master` as merge commits (to preserve
  upstream SHAs). Sync branches are not merged into `release/*`; cherry-pick
  individual upstream commits onto a maintenance line only when needed.
- Feature PRs continue to squash-merge into `master`. Backport PRs into
  `release/*` should preferably be cherry-picks of those squash commits
  (`git cherry-pick -x <sha>`).

### Cutting a maintenance branch

After a minor release tag exists (for example `v1.2.0`):

```bash
./scripts/cut-release-branch.sh 1.2.0
# or: git branch release/1.2 v1.2.0 && git push -u origin release/1.2
```

Prefer the helper script. If the tag already contains the maintenance-branch
Release Please workflow and `.release-please-config.release.json`, the new
branch inherits them unchanged. If you cut from an **older** tag that predates
that automation, `cut-release-branch.sh` commits those files (plus the
website/playground “latest tag only” deploy guards) from `origin/master`
onto the new branch before pushing, as a `ci:` commit that does not bump the
version.

Pushes to `release/1.2` then open patch Release PRs against that branch.

Protect `release/*` the same way as `master` (required PR + required checks).

### Support policy

- The **current** minor (`release/X.Y` matching the latest `vX.Y.*` on
  `master`'s release line) receives regular bugfix patches.
- **Older** minors receive security fixes only, and only while their CI
  runners and packaging targets still exist.
- Close the matching GitHub milestone when a line is retired.

A severe fix that must land on several old lines is cherry-picked onto each
still-supported `release/X.Y`, then each line's Release PR is merged
independently (`1.3.x`, `1.4.x`, …).

## Release Process

### Automated workflow on `master`

1. **Developer makes changes**:
   - Create feature branch: `git checkout -b feature/my-feature`
   - Commit with conventional format: `git commit -m "feat: add feature"`
   - Push and create PR (assign milestone, labels, assignee)

2. **PR is merged to master**:
   - Prefer squash merge for feature/fix PRs
   - Code is now in master

3. **Release Please bot activates**:
   - Analyzes all commits since last release on that branch
   - Calculates new version number
   - Creates or updates a **Release PR** targeting `master`

4. **Release PR contents**:
   - Updated `VERSION.txt`
   - Updated `CHANGELOG.md`
   - Summary of all changes
   - Proposed version number in PR title

5. **Maintainers review Release PR**:
   - Check if version bump is correct
   - Review changelog entries
   - Verify all features are ready for release

6. **Release PR is merged**:
   - Git tag is created (e.g., `v1.3.0`)
   - GitHub Release is published
   - Changelog is finalized

7. **Cut the maintenance branch** from the new tag:
   - `./scripts/cut-release-branch.sh 1.3.0`

8. **CI builds release artifacts**:
   - Binaries built for all platforms
   - Attached to GitHub Release

### Patch releases from `release/X.Y`

1. Land the fix on `master` (unless the code no longer exists there).
2. Cherry-pick onto `release/X.Y` and open a PR into that branch.
3. Release Please opens or updates a patch Release PR
   (title scoped as `chore(release/X.Y): …`).
4. Merge the patch Release PR → tag `vX.Y.Z` → the same release build
   workflows run from the `release` event.
5. Website and playground production deploys only when the new tag is the
   highest semver tag in the repo, so an older-line patch cannot roll the
   site back.

### Batching releases

Multiple features can be included in one minor release on `master`:

```
Day 1: Merge feat: add feature A  → Release PR updated to v1.3.0
Day 2: Merge fix: fix bug B        → still v1.3.0 (fix absorbed into minor)
Day 3: Merge feat: add feature C   → still v1.3.0
Day 4: Merge Release PR            → v1.3.0 released with A, B, and C
Day 5: cut release/1.3 from v1.3.0
```

After `v1.2.0`, features may continue merging to `master` toward `1.3.0`
while hotfixes go to `release/1.2` as `1.2.1`, `1.2.2`, …

## Version Storage

### VERSION.txt

The canonical version is stored in `VERSION.txt` in the repository root:

```
1.2.0
```

This file is:
- Read by CMake during build
- Updated by Release Please when creating releases
- Tracked in version control (per branch)

### Git Tags

Each release is tagged in git:

```bash
git tag -l 'v1.*'
# Output:
# v1.0.0
# v1.1.0
# v1.1.1
# v1.1.2
```

Tags follow the format: `v{MAJOR}.{MINOR}.{PATCH}`

### Release Please config

| File | Used when |
|------|-----------|
| `.release-please-config.json` | Pushes to `master` |
| `.release-please-config.release.json` | Pushes to `release/**` (`always-bump-patch`) |
| `.release-please-manifest.json` | Last released version on that branch |

Each branch carries its own manifest tip. A newer tag on `master` does not
become "latest" for `release/1.2`.

## Changelog

### CHANGELOG.md

Automatically generated and updated by Release Please. Format:

```markdown
# Changelog

## [1.3.0](https://github.com/pythonscad/pythonscad/compare/v1.2.0...v1.3.0) (2026-04-01)

### Features

* add GPU-accelerated rendering ([abc1234](https://github.com/pythonscad/pythonscad/commit/abc1234))
* support Python 3.12 ([def5678](https://github.com/pythonscad/pythonscad/commit/def5678))

### Bug Fixes

* resolve memory leak in CGAL renderer ([ghi9012](https://github.com/pythonscad/pythonscad/commit/ghi9012))
```

Each branch builds its changelog from commits on that branch, so the same
fix can appear in both `master`'s next minor notes and a `release/X.Y`
patch notes.

### Changelog Sections

Commits are grouped by type:

- **Features**: `feat:` commits
- **Bug Fixes**: `fix:` commits
- **Performance Improvements**: `perf:` commits
- **Documentation**: `docs:` commits (optional in changelog)

Other types (`style:`, `refactor:`, `test:`, `build:`, `ci:`, `chore:`) are hidden by default.

## Version in Code

### CMake

During build, CMake reads `VERSION.txt` and extracts version components:

```cmake
# From VERSION.txt: "1.2.0"
OPENSCAD_VERSION = "1.2.0"
OPENSCAD_MAJOR = 1
OPENSCAD_MINOR = 2
OPENSCAD_PATCH = 0
```

### C++ Code

Version is available in `src/version.h`:

```cpp
#include "version.h"

// Full version: "1.2.0"
std::string version = std::string(openscad_versionnumber);

// With git commit (if built from non-release): "1.2.0 (git a1b2c3d)"
std::string display_version = std::string(openscad_displayversionnumber);
```

### Runtime

Check version with:

```bash
./pythonscad --version
# Output: PythonSCAD 1.2.0
```

## Development Builds

### Between Releases

Commits merged to master but not yet released will still show the current version:

```
VERSION.txt contains: 1.2.0
Current commit: a1b2c3d (after v1.2.0 tag)
Display version: 1.2.0 (git a1b2c3d)
```

The git commit hash helps identify development builds.

### CI Builds

GitHub Actions can access the version for artifact naming:

```yaml
- name: Get version
  run: |
    VERSION=$(cat VERSION.txt)
    echo "VERSION=$VERSION" >> $GITHUB_ENV

- name: Create artifact
  run: |
    cp build/pythonscad pythonscad-${VERSION}-linux-x64
```

## Manual Overrides

### Emergency Releases

If you need to create a release without waiting for conventional commits:

1. Manually edit `VERSION.txt`
2. Commit with `chore(release): prepare version X.Y.Z`
3. Create git tag manually: `git tag vX.Y.Z`
4. Push tag: `git push origin vX.Y.Z`

⚠️ This should be rare! Prefer using the automated process.

### Skipping Release Please

If a commit shouldn't trigger a release, use types that don't bump versions:

```bash
git commit -m "docs: update README"        # No version bump
git commit -m "style: fix code formatting" # No version bump
git commit -m "chore: update dependencies" # No version bump
```

## FAQ

### Q: Can I choose my own version number?

**A:** No, versions are determined automatically from commit messages. This is by design to ensure consistency and remove subjective decisions.

### Q: What if the bot chooses the wrong version?

**A:** The bot follows strict rules. If the version seems wrong, check your commit message prefixes. You can amend commits before merging to master. On a maintenance branch, confirm the push used `.release-please-config.release.json`.

### Q: How do I create a major version bump?

**A:** Use `feat!:` or `fix!:`, or include `BREAKING CHANGE:` in the commit footer on `master`.

### Q: Can I see what version will be released before merging?

**A:** Yes! After merging to the target branch, the Release PR will show the proposed version. You can review it before merging the Release PR.

### Q: What if I merge a feature but want to delay the release?

**A:** That's fine! Your feature will be in `master`, and the Release PR will accumulate changes. Merge the Release PR when ready to publish. Leave the next-minor Release PR open while shipping patches from `release/X.Y`.

### Q: How do I ship a hotfix after the next minor has already started on master?

**A:** Cut `release/X.Y` from the last `vX.Y.*` tag (if it does not exist yet), cherry-pick the fix there, and merge that branch's Release PR. Do not merge the fix-only work into `master`'s Release PR if you need the patch number on the old line — `master` will take the same fix as part of the next minor.

### Q: Does Release Please support multiple release lines?

**A:** Yes. It versions each `target-branch` from commits and the manifest on that branch. It does not create maintenance branches or cherry-picks for you.

## References

- [Semantic Versioning](https://semver.org/)
- [Conventional Commits](https://www.conventionalcommits.org/)
- [Release Please](https://github.com/googleapis/release-please)
- [Release Please: multiple release branches](https://github.com/googleapis/release-please-action?tab=readme-ov-file#supporting-multiple-release-branches)
- [Contributing Guide](CONTRIBUTING.md)
