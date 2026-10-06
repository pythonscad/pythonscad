#!/usr/bin/env python3
"""macOS Homebrew dependency installs must be Qt6-only."""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _formula_loop(script: str) -> list[str]:
    match = re.search(
        r"for formula in ([^;]+); do",
        script,
        re.MULTILINE,
    )
    if not match:
        raise AssertionError("could not find Homebrew formula loop")
    return match.group(1).split()


def main() -> int:
    root = _repo_root()
    script_path = root / "scripts" / "macosx-build-homebrew.sh"
    script = script_path.read_text(encoding="utf-8")

    assert "brew install qt5" not in script, "Homebrew script still installs qt5"
    assert "qt@5" not in script, "Homebrew script still references qt@5"
    assert "$USE_QT6" not in script, "Homebrew script still branches on USE_QT6"

    formulas = _formula_loop(script)
    assert "qt" in formulas, formulas
    assert "qscintilla2" in formulas, formulas
    assert "qt5" not in formulas, formulas

    qt5_profile = json.loads(
        (root / "scripts" / "deps" / "profiles" / "qt5.json").read_text(
            encoding="utf-8"
        )
    )
    assert "macos" not in qt5_profile.get("distros", {}), (
        "qt5 profile still lists macOS packages"
    )

    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
