#!/usr/bin/env python3
"""Regression test for scripts/get-dependencies.py profile merging.

Combining profiles (either via a single profile's `extends` list, or via
multiple `--profile` CLI flags) must never drop an earlier profile's
`pre_commands` / `post_commands`. A prior revision replaced the whole
command list whenever the *current* profile in the merge loop also defined
a (shorter, inherited) command list for the same distro, silently dropping
the macOS Qt tap/pin sequence when combining `qt6` + `pythonscad` — see
PR #1134 review discussion.
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

    # Commands that must all survive any way of combining the qt6 and
    # pythonscad profiles for the macos distro.
    required_pre_commands = {
        "brew tap openscad/homebrew-tap",
    }

    combinations = [
        ["qt6", "pythonscad"],
        ["pythonscad", "qt6"],
        ["pythonscad-qt6"],
    ]

    for profiles in combinations:
        cfg = gd.load_multiple_profiles(profiles)
        macos = cfg["distros"].get("macos")
        if macos is None:
            print(f"FAIL: no macos distro resolved for profiles {profiles}", file=sys.stderr)
            return 1

        pre_commands = set(macos.get("pre_commands", []))
        missing = required_pre_commands - pre_commands
        if missing:
            print(
                f"FAIL: profiles {profiles} dropped macos pre_commands: {sorted(missing)}\n"
                f"  resolved pre_commands: {macos.get('pre_commands')}",
                file=sys.stderr,
            )
            return 1

        post_commands = set(macos.get("post_commands", []))
        required_post_command = "bash scripts/ci/prepare-macos-pyqt6-qt.sh"
        if required_post_command not in post_commands:
            print(
                f"FAIL: profiles {profiles} dropped the macOS Qt preparation command: "
                f"{macos.get('post_commands')}",
                file=sys.stderr,
            )
            return 1

    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
