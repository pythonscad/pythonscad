# AGENTS.md

Guidance for AI coding agents working in this repository.

This file follows the cross-agent [AGENTS.md](https://agents.md/) convention
(Agentic AI Foundation / Linux Foundation). The repo root also ships a
git-tracked symlink `CLAUDE.md` → `AGENTS.md` so Claude Code and other tools
that look for vendor-specific filenames pick up the same guidance.

## Project Overview

PythonSCAD is a fork of OpenSCAD that adds native Python language support for
3D modeling. It's a C++ application with integrated Python interpreter,
Qt-based GUI, and OpenGL rendering. The project maintains close sync with
upstream OpenSCAD while adding Python-specific features.

## Build Commands

### Standard Build Process

```bash
# Install dependencies (Linux/BSD; Qt6 matches the default USE_QT6=ON)
sudo ./scripts/get-dependencies.py --profile pythonscad-qt6

# Configure and build
mkdir build
cd build
cmake ..
make -j$(nproc)

# Run the application
./pythonscad
```

For a Qt5 build, install `--profile pythonscad-qt5` and configure with
`-DUSE_QT6=OFF`.

### Build Configuration Options

Key CMake options (pass with `-D` flag):

- `ENABLE_TESTS=ON/OFF` - Enable test suite (default: ON)
- `HEADLESS=ON/OFF` - Build without GUI
- `EXPERIMENTAL=ON/OFF` - Enable experimental features (default: ON)
- `ENABLE_PYTHON=ON/OFF` - Enable Python support (default: ON)
- `USE_QT6=ON/OFF` - Use Qt6 instead of Qt5 (default: ON)
- `ENABLE_CGAL=ON/OFF` - Enable CGAL geometry backend
- `ENABLE_MANIFOLD=ON/OFF` - Enable Manifold geometry backend

Example: `cmake -DHEADLESS=ON -DEXPERIMENTAL=ON ..`

### Testing

```bash
# From build directory
ctest -j8                      # Run default tests with 8 parallel jobs
ctest -R <regex>               # Run tests matching pattern (e.g., ctest -R dxf)
ctest -C Heavy                 # Run time-consuming tests
ctest -C Examples              # Test all examples
ctest -C All                   # Run all tests
```

Catch2 unit-test sources live next to the code under `src/` (for example
`src/core/*_test.cc`, `src/io/*_test.cc`, `src/utils/*_test.cc`), but the
`OpenSCADUnitTests` executable target is currently commented out in the root
`CMakeLists.txt` and is **not** produced by a standard build. Prefer `ctest`
for regression coverage. If you re-enable that target locally, quote Catch2
filter args so the shell does not treat `#` as a comment, e.g.
`./OpenSCADUnitTests -# '#vector_math_test'`.

### Code Formatting

```bash
# Format changed files
./scripts/beautify.sh

# Format all files
./scripts/beautify.sh --all

# Check formatting (used in CI)
./scripts/beautify.sh --check
```

## Development Workflow

### Commit Message Format

**CRITICAL:** PythonSCAD uses Conventional Commits for automated versioning. All commits MUST follow this format:

```text
<type>[optional scope]: <description>

[optional body]

[optional footer]
```

Valid types:

- `feat:` - New feature (bumps version)
- `fix:` - Bug fix (bumps version)
- `docs:` - Documentation only
- `style:` - Code style/formatting
- `refactor:` - Code refactoring
- `perf:` - Performance improvement
- `test:` - Test changes
- `build:` - Build system changes
- `ci:` - CI configuration
- `chore:` - Maintenance tasks

Breaking changes: Add `!` after type (e.g., `feat!:`) or include `BREAKING CHANGE:` in footer.

Examples:

```text
feat: add GPU-accelerated rendering
fix(python): resolve memory leak in pyopenscad module
docs: update build instructions for Ubuntu 24.04
```

### Issue and PR Conventions

**CRITICAL:** Every new issue and pull request MUST be properly categorized.
When creating an issue or PR (whether by hand or via `gh`), always set:

- **Labels** — apply at least one descriptive label that matches the change
  (e.g. `bug`, `enhancement`, `build`, `ci`, `docs`, `refactor`, `perf`,
  `test`, `dependencies`, `security`, `UI/UX`). Never leave an issue or PR
  unlabeled.
- **Assignee** — by default, assign the issue/PR to the person actually
  working on it (the author). Do not leave work items unassigned.
- **Milestone** — required for merge. Every PR must be assigned to a
  minor-version milestone (for example `Quadratum v1.2`, `Rete v1.3`,
  `Sphaera v1.4`). Unless the user explicitly names a different milestone,
  assign new issues and PRs to the **next minor-version milestone** after
  the latest tagged release. Example: if the latest tag is `v1.1.2`, use
  the milestone for `1.2` (currently `Quadratum v1.2`). Do not leave work
  on `Backlog` by default.
- **Issue type** (issues only) — set the correct GitHub issue type. The
  available types are:
  - `Bug` — an unexpected problem or behavior.
  - `Feature` — a request, idea, or new functionality.
  - `Task` — a specific piece of work that is neither a bug nor a feature.

Examples:

```bash
# Create a labeled, typed, assigned, milestoned issue
gh issue create \
  --title "fix(wasm): svg target not compiled with -fPIC" \
  --label bug --label build \
  --assignee @me \
  --milestone "Quadratum v1.2" \
  --body "..."
# Issue type is org-level; set it via the API or the web UI after creation.

# Create a labeled, assigned, milestoned PR
gh pr create \
  --title "..." \
  --label bug --label build \
  --assignee @me \
  --milestone "Quadratum v1.2" \
  --body "..."
```

GitHub issue *types* are an organization-level field separate from labels and
are not settable with `gh issue create` flags; set the type through the issue's
web UI or the GraphQL/REST API after the issue exists.

### Pre-commit Hooks

The project uses pre-commit hooks that run automatically:

- Code formatting (clang-format)
- Trailing whitespace removal
- YAML validation
- Commit message validation
- Private key detection

Install with:

```bash
pip install pre-commit
pre-commit install --hook-type commit-msg --hook-type pre-commit
```

Both `--hook-type` flags are required to validate both code and commit messages.

### PR Iteration Workflow

The repository is configured to auto-cancel superseded workflow runs on push
and to auto-trigger a fresh Copilot review on every push (including
force-pushes). No manual `gh run cancel` or `gh copilot-review` re-request is
needed in the normal case — just push.

Copilot's review typically starts and finishes within 2-10 minutes, well
inside the runtime of the CI matrix. The only fallback case: if every
required check has gone green and no Copilot review has landed at all since
the current HEAD was pushed, request one manually:

```bash
gh copilot-review "https://github.com/$OWNER/$REPO/pull/$PR" \
  --wait --wait-timeout 15min --wait-interval 30sec
```

This should be rare. A completed review with zero new comments (e.g.
"Copilot reviewed 5 out of 6 changed files in this pull request and
generated no new comments.") counts as a completed review — don't treat "no
comments" as "no review happened."

This applies whether the push is a force-push (amend) or a fast-forward (new
commit). The new push will start fresh runs of every required check, so the
previous in-flight ones add no signal.

## Architecture Overview

### Core Components

**Evaluation Pipeline:**

1. **Lexer/Parser** (`src/core/lexer.l`, `src/core/parser.y`) - Tokenizes and parses OpenSCAD/Python scripts into AST
2. **AST & Expression** (`src/core/AST.h`, `src/core/Expression.h`) - Abstract syntax tree representation
3. **Module/Function** (`src/core/module.h`, `src/core/function.h`) - Built-in and user-defined operations
4. **Context & Evaluation** (`src/core/Context.h`, `src/core/EvaluationSession.h`) - Variable scoping and script execution
5. **CSG Tree** (`src/core/CSGNode.h`, `src/core/CSGTreeEvaluator.cc`) - Constructive solid geometry tree construction
6. **Geometry Evaluation** (`src/geometry/GeometryEvaluator.cc`) - Converts CSG nodes to actual geometry
7. **Rendering** (`src/glview/`) - OpenGL visualization

**Key Architectural Concepts:**

- **Node System**: Operations are represented as nodes (e.g., `CsgOpNode`,
  `TransformNode`, `LinearExtrudeNode`). Each node type handles specific
  geometric operations.
- **Geometry Backends**: Supports multiple backends (CGAL for exact geometry,
  Manifold for fast mesh operations). Backend selection happens at geometry
  evaluation.
- **Value System**: `Value` class (`src/core/Value.h`) is a variant type
  supporting numbers, strings, vectors, ranges, functions, and Python objects.
- **Context System**: Variables and functions are looked up through nested
  `Context` objects that maintain scope hierarchy.
- **Caching**: Heavy use of caching (`GeometryCache`, `FontCache`) to avoid
  recomputation during interactive editing.

### Python Integration

Python support is implemented in `src/python/`:

- `pyfunctions.cc` - Python builtin functions exposed to scripts
- `pyopenscad.cc` - Python type for 3D objects (solids as first-class objects)
- `pydata.cc` - Conversion between Python and OpenSCAD Value types
- `FrepNode.cc` - Function representation nodes for Python functions

Python scripts use OpenSCAD as a library:

```python
from pythonscad import *
c = cube([10, 20, 30]).color("Tomato")
show(c)
```

#### Three-module layout

PythonSCAD ships three importable modules; new C-side functionality goes
into the C extension, new pure-Python OpenSCAD-compatible code goes into
the `openscad` overlay, and PythonSCAD-only additions go into the
`pythonscad` overlay:

- `_openscad` - C extension built from `src/python/`. Registered via
  `PyImport_AppendInittab("_openscad", ...)` in the embedded interpreter
  and shipped as `_openscad.so`/`_openscad.pyd` in the pip wheel.
- `openscad` - pure-Python overlay at
  `libraries/python/openscad/__init__.py`. Re-exports `_openscad` 1:1 and
  is the home for any pure-Python reimplementations or compatibility
  shims that should match upstream OpenSCAD's API.
- `pythonscad` - pure-Python overlay at
  `libraries/python/pythonscad/__init__.py`. Re-exports `openscad` and is
  the home for PythonSCAD-only features that should not bleed into the
  OpenSCAD-compatible surface.

Switching a script between `from openscad import *` and
`from pythonscad import *` requires no other code change. See
`doc/python-modules.md` for the per-symbol deprecation pattern that
becomes available because of this layout.

### Directory Structure

- `src/core/` - Language core (parser, AST, evaluation, modules, nodes)
- `src/python/` - Python interpreter integration
- `src/geometry/` - Geometry operations (CGAL backend in `cgal/`, Manifold in `manifold/`)
- `src/glview/` - OpenGL rendering (CSG preview, mesh rendering)
- `src/gui/` - Qt-based GUI (MainWindow, Editor, Console)
- `src/io/` - Import/Export (STL, DXF, SVG, PNG, etc.)
- `src/utils/` - Utility functions
- `tests/` - Test suite and test data
- `examples/` - Example models (used for testing and documentation)
- `scripts/` - Build and maintenance scripts

### Geometry Backends

PythonSCAD supports two geometry computation backends:

**CGAL** (`src/geometry/cgal/`):

- Exact arithmetic, robust for complex operations
- Used for: minkowski, hull, convex operations
- Primary types: `CGALNefGeometry`, `CGALPolyhedron`

**Manifold** (`src/geometry/manifold/`):

- Fast mesh-based operations
- Used for: boolean operations, extrusions (when enabled)
- Can be toggled with `USE_MANIFOLD_TRIANGULATOR` build option

### Rendering Pipeline

1. **CSG Tree Building**: Script evaluation creates tree of CSG operations
2. **Geometry Evaluation**: CSG tree converted to concrete geometry (polygons/meshes)
3. **OpenCSG Preview**: Fast approximate preview using depth buffer techniques
4. **Full Render**: Complete geometric evaluation for export/final view

## Adding New Features

### Test coverage requirement

**CRITICAL:** New functionality MUST ship with test coverage in the same PR
(or an inseparable follow-up commit on the same branch). Do not merge feature
or bug-fix work without a regression test, unit test, or other automated
check that would fail if the change were reverted. Pure docs/chore changes
are exempt when there is nothing meaningful to assert.

### Adding a New Module/Function

1. Implement in `src/core/builtin_functions.cc` (functions) or create new node in `src/core/`
2. If creating new node type:
   - Declare class inheriting from `AbstractNode`
   - Implement geometry evaluation in `src/geometry/GeometryEvaluator.cc`
3. Add Python bindings in `src/python/pyfunctions.cc` if needed
4. Add tests under the appropriate category directory (see [Adding Tests](#adding-tests))
5. Add examples in `examples/`
6. Run with `TEST_GENERATE=1` to create expected outputs: `TEST_GENERATE=1 ctest -R mytest`
7. Document in release notes

### Adding Tests

Tests are organized by **category directories** under `tests/data/`. Many suites
`file(GLOB ...)` those directories from `tests/CMakeLists.txt`, so dropping a
script into the right folder registers it automatically — prefer that over
adding one-off entries to CMake.

**Category directories (examples):**

| Directory | Typical use |
| --- | --- |
| `tests/data/scad/` | OpenSCAD (`.scad`) regression inputs (often further subdivided) |
| `tests/data/pythonscad/` | PythonSCAD geometry / render (`.py` → PNG) |
| `tests/data/pythonscad-echo/` | PythonSCAD echo / textual output |
| `tests/data/pythonscad-failing/` | Expected-failure PythonSCAD cases |
| `tests/data/pythonscad-laser/` | Laser / SVG-oriented PythonSCAD cases |
| `tests/data/pythonscad-3d-export/` | 3D export (e.g. 3MF) |
| `tests/data/pythonscad-multitool-export/` | Multi-tool export |
| `tests/data/pythonscad-export-{stl,stl-bin,3mf,obj,svg,pov}/` | Format-specific export |
| `tests/data/pythonscad-overlay/` | Overlay / import-surface render tests |
| `tests/data/pythonscad-overlay-echo/` | Overlay echo tests |

Example: a new PythonSCAD-specific geometry function belongs as
`tests/data/pythonscad/my-feature.py` (plus expected output under
`tests/regression/`), **not** as a hand-listed one-liner in
`tests/CMakeLists.txt`.

Only add or edit entries in `tests/CMakeLists.txt` when the test does not fit
an existing category (new runner args, a new category directory, special
flags, platform exclusions, etc.).

**Regression tests:**

1. Place the test file in the matching `tests/data/<category>/` directory
2. Reconfigure / rebuild so CMake picks up new GLOB'd files if needed
3. Generate expected output: `TEST_GENERATE=1 ctest -R testname`
4. Verify output in `tests/regression/*/testname-expected.*`
5. Commit the test file and expected output

**Unit tests:**

- Add to an appropriate `*_test.cc` next to the code under `src/` (e.g.
  `src/core/`, `src/io/`, `src/utils/`) or create a new one there
- Uses the Catch2 framework
- The dedicated `OpenSCADUnitTests` binary is not built by default (target is
  commented out in the root `CMakeLists.txt`); re-enable it locally if you need
  to run Catch2 filters directly

## Platform-Specific Notes

### macOS Build

```bash
# Install dependencies via Homebrew (Qt6 only on macOS)
./scripts/get-dependencies.py --yes --profile pythonscad-qt6

# Configure and build
mkdir build
cd build
cmake ..
make -j$(sysctl -n hw.ncpu)
```

### Windows Build

The release packaging path is a native MSYS2 UCRT64 build. Start commands from
the repo root with:

```powershell
C:/msys64/msys2_shell.cmd -defterm -here -no-start -ucrt64 -shell bash
```

Inside that shell, install/update the package set used by CI:

```bash
PACBOY_PACKAGES=$(python3 ./scripts/get-dependencies.py --distro msys2 \
  --profile pythonscad-qt6 --list | tr '\n' ' ')
pacboy -S --noconfirm --needed $PACBOY_PACKAGES nsis:p python-psutil:p uv:p
```

Configure and build the native package target:

```bash
cmake -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DEXPERIMENTAL=ON \
  -DSNAPSHOT=ON \
  -DUSE_QT6=ON \
  -DENABLE_PYTHON=ON \
  -DENABLE_LIBFIVE=ON \
  -DPORTABLE_BINARY=ON \
  -DENABLE_TESTS=OFF \
  -DPACKAGE_ARCH=windows-x86-64
cmake --build build -j4
```

Before packaging, refresh the bundled runtime Python dependencies and staging
tree, then generate ZIP and NSIS artifacts:

```bash
BUNDLE_PY_AUTO_INSTALL_PIP_LICENSES=1 \
  ./scripts/bundle-runtime-python.sh build/pythonscad-bundled-py --python python
cmake --install build --prefix build/staging
cd build
cpack -G ZIP
cpack -G NSIS
```

Smoke-test local artifacts from PowerShell with:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\release-smoke\test-windows.ps1 `
  --skip-download --artifact-dir build --keep-workdir
```

See `doc/win-build.md` for the older MSVC/vcpkg build notes. The historical MXE
cross-build path is `./scripts/mingw-x-build-dependencies.sh 64`.

### WebAssembly Build

PythonSCAD ships a Python-enabled WASM build on **Emscripten 6.0** with **CPython
3.14**. Two variants: `node` (NODERAWFS, for smoke testing) and `web` (MEMFS +
preloaded stdlib, for browser distribution).

Docker uses a content-addressed public toolchain image (no
`openscad/wasm-base` dependency):

**`ghcr.io/pythonscad/wasm-python-base:toolchain-v1-<input-sha256>`** — Emscripten
sysroot, third-party libraries, and cross-compiled CPython 3.14.

`scripts/wasm-base-docker-run.sh` pulls the matching image automatically. A
checkout with unpublished toolchain changes falls back to a local cold build.

```bash
# Unified build (sysroot + CPython, ~60 min first time)
docker build -f docker/wasm/sysroot.dockerfile --target wasm-python-base \
  -t pythonscad-wasm-python-base:local .

# Node variant (fast smoke test)
./scripts/wasm-base-docker-run.sh emcmake cmake -B build-wasm-node \
  -DWASM_BUILD_TYPE=node -DCMAKE_BUILD_TYPE=Release -DEXPERIMENTAL=1
./scripts/wasm-base-docker-run.sh cmake --build build-wasm-node -j$(nproc)
node build-wasm-node/pythonscad.js -o out.stl --trust-python script.py

# Web variant (browser distribution; MAIN_MODULE=2 for Chromium JSPI)
./scripts/wasm-base-docker-run.sh emcmake cmake -B build-wasm-web \
  -DWASM_BUILD_TYPE=web -DCMAKE_BUILD_TYPE=Release -DEXPERIMENTAL=1
./scripts/wasm-base-docker-run.sh cmake --build build-wasm-web -j$(nproc)

# Test in browser (copy pages via docker if build dir is docker-owned)
docker run --rm -v "$PWD:/src" -w /src pythonscad-wasm-python-base:local \
  bash -c 'mkdir -p build-wasm-web/vendor \
    && cp wasm-test/test.html wasm-test/notebook.html build-wasm-web/ \
    && cp wasm-test/vendor/three.module.min.js wasm-test/vendor/three.core.min.js build-wasm-web/vendor/'
python3 wasm-test/serve.py 8080 build-wasm-web/
# http://localhost:8080/test.html  or  /notebook.html
```

See `doc/wasm-build.md` for the full guide including architecture rationale,
smoke testing, JavaScript API (`EmsInitPython` / `EmsEvaluatePython`), and known gotchas.

## Common Debugging Workflows

### Finding Where Feature Is Implemented

1. Search for module/function name in `src/core/builtin_functions.cc`
2. Look for corresponding Node class (e.g., `LinearExtrudeNode` for `linear_extrude`)
3. Find geometry evaluation in `GeometryEvaluator::visit()` methods

### Debugging Geometry Issues

1. Enable debug output in geometry evaluator
2. Use `--render` flag to force full geometry computation
3. Check intermediate CSG tree with `--debug` flag
4. Compare CGAL vs Manifold backends by toggling in build

### Test Failures

- Check `build/Testing/Temporary/` for test logs
- Expected results: `tests/regression/*/`
- Actual results: `build/tests/output/*/`
- Regenerate expected: `TEST_GENERATE=1 ctest -R failing_test`

## Important Files

- `AGENTS.md` - Agent / contributor workflow guidance (this file)
- `CLAUDE.md` - Symlink to `AGENTS.md` for Claude Code compatibility
- `CMakeLists.txt` - Main build configuration
- `VERSION.txt` - Current version (auto-updated by release bot)
- `CHANGELOG.md` - Auto-generated from commit messages
- `.pre-commit-config.yaml` - Pre-commit hook configuration
- `.clang-format` - Code style rules
- `doc/testing.md` - Detailed testing documentation
- `CONTRIBUTING.md` - Contribution guidelines
- `VERSIONING.md` - Versioning and release process

## Versioning

PythonSCAD uses automated semantic versioning via Release Please bot:

- Version is determined from conventional commit messages
- `feat:` commits bump patch version before 1.0, and minor version after 1.0
- `fix:` commits bump patch version
- `feat!:` or `BREAKING CHANGE:` bump minor version before 1.0, and major
  version after 1.0
- Release PR is auto-created when commits are merged to master
- Merging Release PR creates git tag and GitHub release

See `VERSIONING.md` for complete details.
