# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PythonSCAD is a fork of OpenSCAD that adds support for Python as a native scripting language while maintaining compatibility with OpenSCAD scripts. It's a programmatic 3D modeling application for creating models suitable for 3D printing.

Key differences from OpenSCAD:
- Uses Python as the primary scripting language (in addition to OpenSCAD's domain-specific language)
- Provides a functional programming model vs OpenSCAD's descriptive language
- Solids are first-class objects that can be function parameters/return values
- Access to the full Python ecosystem and libraries
- Additional methods like fillets and vertex access

## Build System and Common Commands

### Dependencies Installation
```bash
# Install dependencies on Ubuntu/Debian
sudo ./scripts/uni-get-dependencies.sh

# Check dependency versions
./scripts/check-dependencies.sh
```

### Building
```bash
# Standard build process
mkdir build
cd build
cmake ..
make

# Build with specific options
cmake -DHEADLESS=ON -DEXPERIMENTAL=ON ..
make

# Install after building
sudo make install
```

### Key CMake Configuration Options
- `-DHEADLESS=ON/OFF` - Build without GUI frontend
- `-DEXPERIMENTAL=ON/OFF` - Enable experimental features
- `-DNULLGL=ON/OFF` - Build without OpenGL (implies HEADLESS=ON)
- `-DENABLE_TESTS=ON/OFF` - Run testsuite after building
- `-DENABLE_PYTHON=ON/OFF` - Enable Python interpreter (default ON)
- `-DENABLE_CGAL=ON/OFF` - Enable CGAL backend
- `-DENABLE_MANIFOLD=ON/OFF` - Enable Manifold backend
- `-DUSE_QT6=ON/OFF` - Build GUI with Qt6 (default ON for macOS, OFF elsewhere)

### Testing
```bash
# Run all tests
make test

# Run tests from build directory
ctest

# Run specific test categories (from tests directory)
ctest -R "render-"        # Run rendering tests
ctest -R "preview-"       # Run preview tests
ctest -R "3mf"           # Run 3MF format tests (requires lib3mf)
ctest -R "manifold"      # Run Manifold backend tests
```

### Test Structure
Tests are located in the `tests/` directory with the following structure:
- `tests/CMakeLists.txt` - Main test configuration
- `tests/testdata/` - Test input files (.scad scripts)
- Tests are categorized by type: render, preview, export, import
- Test data includes 2D/3D features, regression tests for issues
- 3MF tests are conditionally enabled based on lib3mf availability

### Platform-Specific Builds

#### macOS
```bash
# Install dependencies with Homebrew
./scripts/macosx-build-homebrew.sh

# Or build from source
source scripts/setenv-macos.sh
./scripts/macosx-build-dependencies.sh
```

#### Windows Cross-compilation
```bash
source ./scripts/setenv-mingw-xbuild.sh 64
./scripts/mingw-x-build-dependencies.sh 64
./scripts/release-common.sh mingw64
```

#### WebAssembly
```bash
./scripts/wasm-base-docker-run.sh emcmake cmake -B build-web -DCMAKE_BUILD_TYPE=Debug -DEXPERIMENTAL=1
./scripts/wasm-base-docker-run.sh cmake --build build-web -j2
```

## Architecture and Code Structure

### Core Components

**Source Code Layout (`src/`):**
- `openscad.cc` - Main command-line application entry point
- `openscad_gui.cc` - GUI application entry point
- `core/` - Core geometry and language processing
- `geometry/` - Geometric operations and data structures
- `glview/` - OpenGL rendering and visualization
- `gui/` - Qt-based user interface
- `io/` - File input/output (STL, 3MF, etc.)
- `python/` - Python interpreter integration
- `platform/` - Platform-specific code

**Key Directories:**
- `libraries/` - External library dependencies
- `examples/` - Example .scad files
- `tests/` - Comprehensive test suite
- `scripts/` - Build and utility scripts
- `submodules/` - Git submodules for dependencies

### Backend Architecture
PythonSCAD supports multiple geometry backends:
- **CGAL** - Computational Geometry Algorithms Library (default)
- **Manifold** - Fast, robust manifold-based operations
- Backend selection affects rendering, boolean operations, and mesh processing

### Python Integration
- Python interpreter is embedded and enabled by default
- Python scripts can import the `openscad` module to access 3D modeling functions
- Supports both .scad (OpenSCAD language) and .py files
- Python objects can represent 3D solids as first-class entities

### Test Architecture
- Tests use CMake's CTest framework
- Image comparison tests for visual regression testing
- Export/import tests for file format compatibility
- Backend-specific tests (CGAL vs Manifold)
- Conditional test execution based on available dependencies (e.g., lib3mf)

### Build Dependencies
Major dependencies include:
- Qt5/6 for GUI
- CGAL for computational geometry
- Manifold for alternative geometry backend
- Python 3.8+ for scripting support
- OpenCSG for constructive solid geometry
- Various format libraries (lib3mf for 3MF files)

### Key Build Features
- Supports headless builds for server environments
- WebAssembly compilation for browser usage
- Cross-platform support (Linux, macOS, Windows)
- Optional features controlled via CMake options
- Extensive test suite with multiple backend validation
