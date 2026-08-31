#!/usr/bin/env python3
"""Remove free comparison-operator declarations from all PyQt6 .sip files.

These declarations assume that Qt provides a matching free
``bool operator==/!=/</<=/>/>=(...)``. Newer Qt 6 patch releases break ADL
resolution for several of them because they collide with QCborTag and
QCborKnownTags operators from unrelated headers. Removing the declarations
only disables Python-side comparisons (==, !=, <, <=, >, >=) for the affected
Qt type; all other functionality remains unchanged.
"""
import re
import sys
from pathlib import Path

# Match only unindented (column-zero) free operator declarations that end
# directly in ");". Member operators ("    bool operator==(...) const;") have
# leading whitespace and additional text before the semicolon, so the "^bool"
# anchor and the absence of " const" before ";" exclude them.
PATTERN = re.compile(r'^bool operator(==|!=|<=|>=|<|>)\([^)]*\);\s*$')


def patch_file(path: Path) -> int:
    lines = path.read_text().splitlines(keepends=True)
    kept = []
    removed = 0
    for line in lines:
        if PATTERN.match(line):
            removed += 1
            continue
        kept.append(line)
    if removed:
        path.write_text("".join(kept))
    return removed


def main():
    if len(sys.argv) != 2:
        print("Usage: patch_pyqt6_free_operators.py <sip-root-dir>", file=sys.stderr)
        sys.exit(1)

    root = Path(sys.argv[1])
    total = 0
    for sip_file in root.rglob("*.sip"):
        n = patch_file(sip_file)
        if n:
            print(f"{sip_file}: removed {n} operator declaration(s)")
            total += n
    print(f"Total removed: {total}")


if __name__ == "__main__":
    main()
