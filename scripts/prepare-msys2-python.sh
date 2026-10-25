#!/usr/bin/env bash
set -euo pipefail

# Copy MSYS2 Python headers into a local folder. The import library is
# generated separately from the embeddable python3XY.dll.
# Intended for native MSYS2 Windows builds (not MXE cross-compilation).
#
# Usage:
#   prepare-msys2-python.sh <callstack> <python_version_minor> <output_dir>
#
# Example:
#   prepare-msys2-python.sh ucrt 3.13 /path/to/source/python_mingw

callstack="${1:?Usage: $0 <callstack> <python_version_minor> <output_dir>}"
python_minor="${2:?Missing python version (major.minor)}"
output_dir="${3:?Missing output directory}"

system_lib="/${callstack}64/lib/libpython${python_minor}.dll.a"
system_include="/${callstack}64/include/python${python_minor}"

if [[ ! -f "$system_lib" ]]; then
  echo "ERROR: Python import lib not found: $system_lib" >&2
  exit 1
fi
if [[ ! -d "$system_include" ]]; then
  echo "ERROR: Python include dir not found: $system_include" >&2
  exit 1
fi

rm -rf "$output_dir"
mkdir -p "$output_dir/${callstack}64/include"

cp -r "$system_include" "$output_dir/${callstack}64/include/"
# The import library is generated later from the embeddable python3XY.dll.
touch "$output_dir/.headers-ready"

echo "MSYS2 Python headers copied to: $output_dir"
