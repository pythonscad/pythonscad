#!/bin/sh
#
# Generate a MinGW import library from the Python DLL that will be shipped.
#
# The native Windows build links with MinGW but ships python.org's embeddable
# runtime, whose DLL is named python3XY.dll. MSYS2's own import library names
# libpython3.X.dll instead. Generating the import library from the shipped
# DLL keeps the linker and the runtime on the same name.
#
# Usage: libpython_patch.sh <python.dll> <output.dll.a>
#
# Requires gendef and dlltool (MSYS2 package mingw-w64-ucrt-x86_64-tools).

set -eu

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <python.dll> <output.dll.a>" >&2
    exit 1
fi

DLL=$1
OUT=$2

if [ ! -f "$DLL" ]; then
    echo "Error: Python DLL not found: $DLL" >&2
    exit 1
fi

if ! command -v gendef >/dev/null 2>&1 || ! command -v dlltool >/dev/null 2>&1; then
    echo "Error: gendef and dlltool are required (pacboy -S tools:p)" >&2
    exit 1
fi

DLL_ABS=$(cd "$(dirname "$DLL")" && pwd)/$(basename "$DLL")
WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT

cp "$DLL_ABS" "$WORKDIR/"
(
    cd "$WORKDIR"
    gendef "$(basename "$DLL_ABS")"
)

DEF=$(find "$WORKDIR" -maxdepth 1 -name '*.def' | head -1)
if [ -z "$DEF" ]; then
    echo "Error: gendef produced no .def for $(basename "$DLL_ABS")" >&2
    exit 1
fi

mkdir -p "$(dirname "$OUT")"
dlltool --dllname "$(basename "$DLL_ABS")" --def "$DEF" --output-lib "$OUT"
echo "Generated $OUT from $(basename "$DLL_ABS")"
