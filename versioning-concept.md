# PythonSCAD Versioning Concept

## Overview

This document outlines a strategy for migrating PythonSCAD from date-based versioning to semantic versioning with automated version management, changelog generation, and frequent releases following the "release early, release often" paradigm.

## Current Situation

- **PythonSCAD**: Currently uses date-based versioning (YYYY.MM.DD)
- **OpenSCAD (upstream)**: Still uses date-based versioning, last release 4+ years ago
- **Goal**: Frequent releases with semantic versioning and automation
- **Platform**: GitHub with planned GitHub Actions CI/CD

## Semantic Versioning Strategy

### Version Format: `MAJOR.MINOR.PATCH`

- **MAJOR**: Breaking changes, incompatible API changes
- **MINOR**: New features, backward-compatible functionality
- **PATCH**: Bug fixes, backward-compatible changes

### Starting Version
Since PythonSCAD is a fork, consider starting at `1.0.0` to indicate:
- Stable fork status
- Distinct identity from OpenSCAD
- Ready for production use

## Automated Versioning Solutions

### Option 1: Conventional Commits + semantic-release

**Tools:**
- [Conventional Commits](https://www.conventionalcommits.org/) specification
- [semantic-release](https://github.com/semantic-release/semantic-release)
- GitHub Actions integration

**Workflow:**
```
feat: add new Python wrapper system    → MINOR bump
fix: resolve operator precedence issue → PATCH bump
feat!: change core API structure      → MAJOR bump
docs: update README                    → no version bump
```

**Benefits:**
- Fully automated versioning based on commit messages
- Automatic changelog generation
- Release notes generation
- GitHub releases creation
- No manual version management

**Implementation:**
```yaml
# .github/workflows/release.yml
name: Release
on:
  push:
    branches: [main]
jobs:
  release:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
        with:
          fetch-depth: 0
      - uses: actions/setup-node@v4
      - run: npx semantic-release
```

### Option 2: Release Please

**Tool:** [Release Please](https://github.com/googleapis/release-please) (Google)

**Features:**
- Conventional commits based
- Automatic PR creation for releases
- Multi-language support (supports CMake)
- Changelog automation
- Version bumping in source files

**Benefits:**
- More conservative approach (PR-based)
- Better visibility into upcoming releases
- Supports CMake projects natively
- Can update version in CMakeLists.txt automatically

### Option 3: Changesets

**Tool:** [Changesets](https://github.com/changesets/changesets)

**Approach:**
- Developers write changeset files for each PR
- Tool aggregates changes for releases
- Manual control over what goes into each release

**Benefits:**
- More control over release content
- Better for teams wanting to curate releases
- Excellent changelog quality

## Git Branching Strategy

### Recommended: GitHub Flow with Release Automation

```
main (protected)
├── feature/python-wrapper-enhancement
├── feature/new-3d-operations
├── hotfix/critical-bug-fix
└── upstream/openscad-sync
```

**Branch Types:**

1. **`main`** - Always deployable, protected branch
   - All commits trigger semantic analysis
   - Automated releases on merge
   - Requires PR for changes

2. **`feature/*`** - New features and enhancements
   - Merged via PR to main
   - Use conventional commit messages

3. **`hotfix/*`** - Critical bug fixes
   - Can be fast-tracked to main
   - Automatically trigger patch releases

4. **`upstream/openscad-sync`** - OpenSCAD synchronization
   - Regular merges from OpenSCAD upstream
   - Carefully managed to avoid conflicts

### Alternative: GitFlow with Automation

```
main (production releases)
├── develop (integration branch)
│   ├── feature/python-improvements
│   ├── feature/new-geometry-ops
│   └── release/1.2.0 (release preparation)
└── hotfix/1.1.1 (emergency fixes)
```

**Benefits:**
- More structured for complex projects
- Separate development and release branches
- Better for coordinated releases

## Recommended Toolchain

### Core Tools

1. **semantic-release** with GitHub Actions
   - Zero-configuration semantic versioning
   - Automatic changelog generation
   - GitHub releases with assets
   - CMake version file updates

2. **Conventional Commits**
   - Standardized commit message format
   - Enables automatic analysis
   - Clear communication of changes

3. **Commitizen** (optional)
   - CLI tool for writing conventional commits
   - Reduces human error in commit messages

### GitHub Actions Workflow

```yaml
# .github/workflows/release.yml
name: Release

on:
  push:
    branches: [main]

jobs:
  release:
    runs-on: ubuntu-latest
    permissions:
      contents: write
      issues: write
      pull-requests: write

    steps:
      - name: Checkout
        uses: actions/checkout@v4
        with:
          fetch-depth: 0
          token: ${{ secrets.GITHUB_TOKEN }}

      - name: Setup Node.js
        uses: actions/setup-node@v4
        with:
          node-version: 20

      - name: Install semantic-release
        run: |
          npm install -g semantic-release
          npm install -g @semantic-release/changelog
          npm install -g @semantic-release/git

      - name: Run semantic-release
        env:
          GITHUB_TOKEN: ${{ secrets.GITHUB_TOKEN }}
        run: semantic-release

  build:
    needs: release
    if: needs.release.outputs.new_release_published == 'true'
    strategy:
      matrix:
        os: [ubuntu-latest, windows-latest, macos-latest]
    runs-on: ${{ matrix.os }}

    steps:
      - name: Checkout
        uses: actions/checkout@v4

      - name: Build PythonSCAD
        run: |
          # Platform-specific build commands
          # Upload artifacts to GitHub release
```

### Configuration Files

**`.releaserc.json`:**
```json
{
  "branches": ["main"],
  "plugins": [
    "@semantic-release/commit-analyzer",
    "@semantic-release/release-notes-generator",
    "@semantic-release/changelog",
    [
      "@semantic-release/exec",
      {
        "prepareCmd": "cmake -DOPENSCAD_VERSION=${nextRelease.version} ."
      }
    ],
    "@semantic-release/github",
    [
      "@semantic-release/git",
      {
        "assets": ["CHANGELOG.md", "CMakeLists.txt"],
        "message": "chore(release): ${nextRelease.version} [skip ci]\n\n${nextRelease.notes}"
      }
    ]
  ]
}
```

## Integration with Existing Build System

### CMake Integration

Update `cmake/Modules/openscad_version.cmake`:

```cmake
# Check if OPENSCAD_VERSION is semantic version
if(OPENSCAD_VERSION MATCHES "^([0-9]+)\\.([0-9]+)\\.([0-9]+).*")
  set(OPENSCAD_MAJOR ${CMAKE_MATCH_1})
  set(OPENSCAD_MINOR ${CMAKE_MATCH_2})
  set(OPENSCAD_PATCH ${CMAKE_MATCH_3})
  set(OPENSCAD_SHORTVERSION "${OPENSCAD_MAJOR}.${OPENSCAD_MINOR}.${OPENSCAD_PATCH}")

  # Set legacy date fields for compatibility
  string(TIMESTAMP CURRENT_YEAR "%Y")
  string(TIMESTAMP CURRENT_MONTH "%m")
  string(TIMESTAMP CURRENT_DAY "%d")
  set(OPENSCAD_YEAR ${CURRENT_YEAR})
  set(OPENSCAD_MONTH ${CURRENT_MONTH})
  set(OPENSCAD_DAY ${CURRENT_DAY})
else()
  # Fallback to existing date-based logic
  if("${OPENSCAD_VERSION}" STREQUAL "")
    string(TIMESTAMP OPENSCAD_VERSION "%Y.%m.%d")
  endif()
  # ... existing logic
endif()
```

## Migration Plan

### Phase 1: Preparation (Week 1)
1. Implement semantic versioning detection in CMake
2. Set up GitHub Actions workflows
3. Configure semantic-release
4. Train team on conventional commits

### Phase 2: Transition (Week 2)
1. Create initial semantic version (1.0.0)
2. Update all documentation
3. Test automation on feature branches
4. Migrate to new branching strategy

### Phase 3: Full Automation (Week 3+)
1. Enable automated releases on main branch
2. Set up cross-platform build automation
3. Monitor and refine the process
4. Establish OpenSCAD upstream sync workflow

## Benefits of This Approach

1. **Automated Version Management**: No manual version bumping
2. **Consistent Releases**: Every merge to main can trigger a release
3. **Clear Communication**: Semantic versions convey meaning
4. **Rich Changelogs**: Automatically generated from commits
5. **Fast Iterations**: Support for "release early, release often"
6. **Upstream Compatibility**: Can still sync with OpenSCAD easily
7. **Professional Ecosystem**: Integrates with package managers and tools

## OpenSCAD Upstream Synchronization

### Strategy for Staying Current

1. **Regular Sync Branch**: `upstream/openscad-sync`
2. **Automated Monitoring**: GitHub Actions to check for OpenSCAD updates
3. **Careful Merging**: Review all upstream changes for conflicts
4. **Version Mapping**: Document which PythonSCAD versions include which OpenSCAD changes

### Example Sync Workflow

```bash
# Monthly or when OpenSCAD updates
git checkout -b upstream/openscad-sync
git pull upstream main  # OpenSCAD repository
# Resolve conflicts, test changes
git checkout main
git merge upstream/openscad-sync
# Commits trigger automatic release with appropriate version bump
```

## Conclusion

The recommended approach combines semantic-release with GitHub Actions for maximum automation while maintaining flexibility. This enables rapid, frequent releases while ensuring version numbers are meaningful and changelogs are comprehensive.

The key is starting simple with semantic-release and conventional commits, then expanding the automation as the project grows and requirements become clearer.
