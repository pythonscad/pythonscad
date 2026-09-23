#!/usr/bin/env bash
# Install native build dependencies inside the manylinux_2_28 build container.
# Runs as cibuildwheel [tool.cibuildwheel.linux] before-all.
set -euo pipefail

echo "=== install-deps-manylinux.sh: installing pip wheel build deps ==="

# manylinux_2_28 images are AlmaLinux/RHEL 8 based.
if command -v dnf >/dev/null 2>&1; then
    PKG_MGR=dnf
else
    PKG_MGR=yum
fi

# EPEL + CRB/powertools provide CGAL and several -devel packages on EL8.
$PKG_MGR install -y epel-release
$PKG_MGR config-manager --set-enabled powertools 2>/dev/null \
    || $PKG_MGR config-manager --set-enabled crb 2>/dev/null \
    || true

$PKG_MGR install -y \
    bison \
    boost-devel \
    cairo-devel \
    CGAL-devel \
    cmake \
    curl \
    double-conversion-devel \
    eigen3-devel \
    flex \
    fontconfig-devel \
    freetype-devel \
    gcc \
    gcc-c++ \
    glib2-devel \
    gmp-devel \
    harfbuzz-devel \
    libzip-devel \
    libxml2-devel \
    make \
    mpfr-devel \
    patch \
    pkgconfig

# cibuildwheel exports CC/CXX to gcc-toolset-12 before running before-all,
# so install it before any source-built fallback dependencies invoke CMake.
$PKG_MGR install -y gcc-toolset-12
/opt/rh/gcc-toolset-12/root/usr/bin/g++ --version

if ! $PKG_MGR install -y lib3mf-devel; then
    echo "lib3mf-devel unavailable from EL8 repos; building lib3mf from source"
    LIB3MF_VERSION=2.4.1
    LIB3MF_SRC="/tmp/lib3mf-${LIB3MF_VERSION}"
    rm -rf "$LIB3MF_SRC"
    curl --fail --show-error --retry 3 --retry-delay 5 -L \
        "https://github.com/3MFConsortium/lib3mf/archive/v${LIB3MF_VERSION}.tar.gz" \
        -o "/tmp/lib3mf-${LIB3MF_VERSION}.tar.gz"
    tar -C /tmp -xzf "/tmp/lib3mf-${LIB3MF_VERSION}.tar.gz"
    cd "$LIB3MF_SRC"
    # The patch removes non-portable compiler/linker assumptions and applies to manylinux too.
    patch -p1 < /project/patches/lib3mf-macos.patch
    # Do NOT pass a *relative* -DCMAKE_INSTALL_LIBDIR=lib64: with the
    # lib3mf-macos.patch order (project() then GNUInstallDirs), CMake
    # absolutizes a pre-set relative LIBDIR against the source tree
    # (/tmp/lib3mf-*/lib64), so pkg-config never sees the install (#1025).
    # Use an absolute libdir so the .so/.pc land in manylinux's search path
    # (/usr/lib64/pkgconfig). The lib3mf.pc.in template then emits a broken
    # libdir=${exec_prefix}//usr/lib64 — rewrite the .pc after install.
    cmake -S . -B build \
        -DLIB3MF_TESTS=OFF \
        -DUSE_INCLUDED_ZLIB=OFF \
        -DUSE_INCLUDED_LIBZIP=OFF \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=/usr \
        -DCMAKE_INSTALL_LIBDIR=/usr/lib64
    cmake --build build -j"$(nproc)"
    cmake --install build

    mkdir -p /usr/lib64/pkgconfig
    PC_FILE=""
    for candidate in \
        /usr/lib64/pkgconfig/lib3mf.pc \
        /usr/lib/pkgconfig/lib3mf.pc \
        "$LIB3MF_SRC"/lib64/pkgconfig/lib3mf.pc \
        "$LIB3MF_SRC"/lib/pkgconfig/lib3mf.pc
    do
        if [[ -f "$candidate" ]]; then
            PC_FILE=$candidate
            break
        fi
    done
    if [[ -z "$PC_FILE" ]]; then
        echo "ERROR: source-built lib3mf.pc was not installed" >&2
        find /usr "$LIB3MF_SRC" -name 'lib3mf.pc' 2>/dev/null || true
        exit 1
    fi
    if [[ "$PC_FILE" != /usr/lib64/pkgconfig/lib3mf.pc ]]; then
        cp -a "$PC_FILE" /usr/lib64/pkgconfig/lib3mf.pc
        PC_FILE=/usr/lib64/pkgconfig/lib3mf.pc
    fi
    # Normalize paths after absolute CMAKE_INSTALL_LIBDIR (see comment above).
    sed -i \
        -e 's|^prefix=.*|prefix=/usr|' \
        -e 's|^exec_prefix=.*|exec_prefix=/usr|' \
        -e 's|^libdir=.*|libdir=/usr/lib64|' \
        -e 's|^includedir=.*|includedir=/usr/include|' \
        "$PC_FILE"

    if [[ ! -e /usr/lib64/lib3mf.so ]]; then
        echo "ERROR: source-built lib3mf.so was not installed under /usr/lib64" >&2
        find /usr "$LIB3MF_SRC" -name 'lib3mf.so*' 2>/dev/null || true
        exit 1
    fi

    # before-all env does not persist into the wheel compile; put the .pc on
    # the default manylinux search path and refresh the linker cache instead.
    export PKG_CONFIG_PATH="/usr/lib64/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
    ldconfig

    # Fail closed: wheel builds must not silently ship dummy 3MF stubs.
    if ! pkg-config --exists lib3mf; then
        echo "ERROR: source-built lib3mf is not visible to pkg-config" >&2
        echo "PC_FILE=${PC_FILE} PKG_CONFIG_PATH=${PKG_CONFIG_PATH:-}" >&2
        pkg-config --debug --exists lib3mf 2>&1 | tail -40 || true
        find /usr "$LIB3MF_SRC" -name 'lib3mf.pc' 2>/dev/null || true
        find /usr "$LIB3MF_SRC" -name 'lib3mf.so*' 2>/dev/null || true
        exit 1
    fi
    echo "lib3mf: pkg-config $(pkg-config --modversion lib3mf) libs=$(pkg-config --libs lib3mf) cflags=$(pkg-config --cflags lib3mf)"
    cd /project
fi

# Distro packages and the source-build fallback both must be discoverable
# before any wheel compile starts.
export PKG_CONFIG_PATH="/usr/lib64/pkgconfig:/usr/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
if ! pkg-config --exists lib3mf && ! pkg-config --exists lib3MF; then
    echo "ERROR: lib3mf is required for manylinux wheels but was not found via pkg-config" >&2
    exit 1
fi

# libfive tree.cpp uses std::optional without including <optional>; EL8 libstdc++
# does not pull it in transitively.
TREE_CPP="submodules/libfive/libfive/src/tree/tree.cpp"
if [[ -f "$TREE_CPP" ]] && ! grep -q '#include <optional>' "$TREE_CPP"; then
    sed -i '/#include <stack>/a #include <optional>' "$TREE_CPP"
    echo "Patched libfive tree.cpp: added #include <optional>"
fi

echo "=== install-deps-manylinux.sh: done ==="
