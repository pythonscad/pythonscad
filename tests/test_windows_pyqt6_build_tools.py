#!/usr/bin/env python3
"""Regression test: the Windows PyQt6 build must have mingw32-make available.

The v1.2.1 release build failed after an hour in "Compile PyQt6 bindings" with
`mingw32-make: command not found`, because the msys2 dependency profile did not
install mingw-w64-ucrt-x86_64-make and the toolchain group no longer provided
it. This checks that the profile installs it ("make:p" via pacboy) and that
scripts/ci/build-pyqt6-bindings-windows.sh fails fast when it is missing.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _load_get_dependencies_module():
    script_path = _repo_root() / "scripts" / "get-dependencies.py"
    spec = importlib.util.spec_from_file_location("get_dependencies", script_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> int:
    gd = _load_get_dependencies_module()

    for profiles in (["pythonscad-qt6"], ["pythonscad"]):
        cfg = gd.load_multiple_profiles(profiles)
        msys2 = cfg["distros"].get("msys2")
        if msys2 is None:
            print(f"FAIL: no msys2 distro resolved for profiles {profiles}", file=sys.stderr)
            return 1
        packages = msys2.get("packages", [])
        if "make:p" not in packages:
            print(
                f"FAIL: msys2 packages for {profiles} do not include 'make:p' "
                f"(mingw-w64-ucrt-x86_64-make): {packages}",
                file=sys.stderr,
            )
            return 1

    script = (_repo_root() / "scripts" / "ci" / "build-pyqt6-bindings-windows.sh").read_text(
        encoding="utf-8"
    )
    check = "command -v mingw32-make"
    first_make_use = script.find("\nmingw32-make ")
    if check not in script or first_make_use == -1 or script.find(check) > first_make_use:
        print(
            "FAIL: build-pyqt6-bindings-windows.sh must check for mingw32-make "
            "before it is first used",
            file=sys.stderr,
        )
        return 1

    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
