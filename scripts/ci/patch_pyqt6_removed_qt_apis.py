#!/usr/bin/env python3
"""Drop PyQt6 .sip declarations for Qt APIs that newer Qt releases removed.

We sometimes build an older PyQt6 sdist against a newer Qt 6 minor (PyPI has
no PyQt6 release for that minor yet, see pyqt6-versions.json). PyQt's
version timeline cannot express "removed in a Qt version it does not know",
so the generated C++ then refers to symbols that no longer exist.

Each entry below lists the first Qt version that removed the API and the
exact .sip lines to drop. The patch is a no-op for older Qt versions.

Usage: patch_pyqt6_removed_qt_apis.py <sip-root-dir> <qt-version>
"""
from __future__ import annotations

import sys
from pathlib import Path

# (first Qt version without the API, sip file relative to sip/, lines to drop)
REMOVALS = [
    # Qt 6.12 removed keypad navigation; QEvent::EnterEditFocus and
    # QEvent::LeaveEditFocus are gone. PyQt6 6.11 always compiles QtCore with
    # QT_KEYPAD_NAVIGATION and lists them unconditionally.
    ((6, 12), "QtCore/qcoreevent.sip", ["EnterEditFocus,", "LeaveEditFocus,"]),
]


def parse_version(text: str) -> tuple[int, int]:
    parts = text.strip().split(".")
    return int(parts[0]), int(parts[1])


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__.strip().splitlines()[-1], file=sys.stderr)
        return 1

    root = Path(sys.argv[1])
    qt = parse_version(sys.argv[2])
    for since, rel_path, entries in REMOVALS:
        if qt < since:
            continue
        path = root / rel_path
        lines = path.read_text().splitlines(keepends=True)
        kept = [line for line in lines if line.strip() not in entries]
        removed = len(lines) - len(kept)
        if removed == 0:
            # Nothing matched: a newer PyQt6 already handles it. Stay quiet but
            # visible so stale entries can be cleaned up.
            print(f"{path}: nothing to remove for Qt >= {since[0]}.{since[1]}")
            continue
        path.write_text("".join(kept))
        print(f"{path}: removed {removed} line(s) for Qt >= {since[0]}.{since[1]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
