# Upstream Sync Process (OpenSCAD → PythonSCAD)

This repository is a long-lived fork of **OpenSCAD**. Changes from upstream
(`openscad/openscad`) are periodically synced into it using
**small, reviewable chunks**, so that:

- merge conflicts are handled by the right domain expert (build vs Python C-API),
- regressions are easier to pinpoint,
- Git/GitHub can correctly determine whether we are “behind” upstream
  (commit-graph ancestry).

> **Important:** To make GitHub stop showing “this branch is X commits behind”,
  **upstream commits themselves** must be integrated (their SHAs appear in
  PythonSCAD's history). Cherry-picking upstream commits must be avoided,
  because it creates new SHAs and can leave GitHub thinking we are still behind.

This document describes the recommended process for syncing upstream changes
into PythonSCAD.

---

## 1 Terminology

- **upstream**: the OpenSCAD repository remote (`https://github.com/openscad/openscad`).
- **origin**: our fork (PythonSCAD) GitHub repository.
- **sync branch**: a temporary branch created for each sync cycle.
- **last sync tag**: an annotated tag marking the last upstream commit we
  synced to.

---

## 2 One-time setup

### 2.1 Add the upstream remote

```bash
git remote add upstream https://github.com/openscad/openscad.git
git fetch upstream
```

### 2.2 (Optional, recommended) Enable rerere

This makes Git remember how recurring conflicts have been resolved and can
automatically reapply those resolutions in future syncs.

```bash
git config --global rerere.enabled true
```

---

## 3 “Last synced” convention (annotated tag)

The last synced upstream state is tracked by creating an **annotated tag** that
points to the upstream `master` commit that was synced up to.

### 3.1 Tag name format

Use:

- `upstream-sync/openscad-YYYY-MM-DD`

Example:

- `upstream-sync/openscad-2026-01-29`

---

## 4 Sync workflow (PR-driven, small chunks)

### Goal

- Identify which PRs were merged upstream since the last sync.
- Merge upstream changes into our fork in the
  **same order upstream landed them** (first-parent order), but still keep
  PR attribution.
- Allow domain experts to resolve conflicts where appropriate.

### Requirements

- GitHub CLI installed and authenticated: `gh auth login`
- `upstream` remote configured

---

### 4.1 Create `sync/openscad-YYYY-MM-DD` branch

```bash
git checkout -b sync/openscad-$(date +%Y-%m-%d) origin/master

git fetch upstream master
```
---

### 4.2 Generate a sync plan (recommended)

Generate a “merge plan” using:

```bash
python3 scripts/plan_openscad_sync.py
```

### 4.3 Apply the plan (merge commits one-by-one)

For each SHA listed by the plan (in order):

```bash
git merge --no-ff <SHA>
# resolve conflicts
# run build/tests
```

### 4.4 Push the branch to origin and open a PR

Push the sync branch to origin:

```bash
git push -u origin sync/openscad-$(date +%Y-%m-%d)
```

Open a PR into `master`.

---

### 4.6 Finish: tag the upstream tip we synced to

```bash
# ensure upstream refs are current
git fetch upstream

UP_TIP=$(git rev-parse upstream/master)

# create an annotated tag pointing at the upstream tip commit
TAG="upstream-sync/openscad-$(date +%Y-%m-%d)"
git tag -a "$TAG" "$UP_TIP" -m "Synced OpenSCAD up to $UP_TIP"

# push the tag to PythonSCAD
git push origin "$TAG"
```

---

## 5 Troubleshooting

### 5.1 If a PR’s mergeCommit is missing

Sometimes PR metadata can be odd depending on merge mode or repo settings.

Fallback options:

- Find the commit on upstream `master` that mentions the PR number:

  ```bash
  git log upstream/master --first-parent --grep "#<PR_NUMBER>" --oneline
  ```

- Or search via GitHub CLI:

  ```bash
  gh pr list -R openscad/openscad --state merged --search "<SHA>"
  ```

---
