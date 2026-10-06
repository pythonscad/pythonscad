#!/usr/bin/env python3
"""Regression tests for scripts/ci/compare-git-release.sh.

Covers the delayed-release compare used by the Manifold submodule
update workflow: empty/null tags, matching peeled tag SHAs, and
mismatched SHAs (including annotated tags, which must not be
compared as tag names).
"""
from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def _run_git(repo: Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(repo), *args],
        text=True,
    ).strip()


def _init_repo() -> Path:
    repo = Path(tempfile.mkdtemp(prefix="compare-git-release-"))
    _run_git(repo, "init")
    _run_git(repo, "config", "user.name", "PythonSCAD CI")
    _run_git(repo, "config", "user.email", "ci@example.com")
    _run_git(repo, "config", "commit.gpgsign", "false")
    (repo / "README").write_text("one\n", encoding="utf-8")
    _run_git(repo, "add", "README")
    _run_git(repo, "commit", "-m", "first")
    return repo


def _compare(script: Path, repo: Path, tag: str, current: str) -> dict[str, str]:
    proc = subprocess.run(
        [
            "bash",
            str(script),
            "--git-dir",
            str(repo),
            "--tag",
            tag,
            "--current-commit",
            current,
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise AssertionError(
            f"compare-git-release.sh exited {proc.returncode}\n"
            f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )
    out: dict[str, str] = {}
    for line in proc.stdout.splitlines():
        if "=" in line and not line.startswith("New ") and not line.startswith("Already "):
            key, _, value = line.partition("=")
            out[key] = value
    return out


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: test_compare_git_release.py <compare-git-release.sh>", file=sys.stderr)
        return 2

    script = Path(sys.argv[1])
    if not script.is_file():
        print(f"script not found: {script}", file=sys.stderr)
        return 2

    if shutil.which("bash") is None or shutil.which("git") is None:
        print("SKIP: bash or git not on PATH")
        return 0

    repo = _init_repo()
    try:
        first = _run_git(repo, "rev-parse", "HEAD")
        _run_git(
            repo,
            "tag",
            "-a",
            "v1.0.0",
            "-m",
            "annotated release",
        )
        peeled = _run_git(repo, "rev-parse", "v1.0.0^{commit}")
        assert peeled == first, (peeled, first)

        empty = _compare(script, repo, "", first)
        assert empty.get("new_release_available") == "false", empty
        assert "latest_commit" not in empty, empty

        null = _compare(script, repo, "null", first)
        assert null.get("new_release_available") == "false", null

        matching = _compare(script, repo, "v1.0.0", first)
        assert matching.get("new_release_available") == "false", matching
        assert matching.get("latest_commit") == first, matching

        (repo / "README").write_text("two\n", encoding="utf-8")
        _run_git(repo, "add", "README")
        _run_git(repo, "commit", "-m", "second")
        second = _run_git(repo, "rev-parse", "HEAD")
        assert second != first

        mismatched = _compare(script, repo, "v1.0.0", second)
        assert mismatched.get("new_release_available") == "true", mismatched
        assert mismatched.get("latest_commit") == first, mismatched
        assert mismatched.get("release_tag") == "v1.0.0", mismatched
    finally:
        shutil.rmtree(repo, ignore_errors=True)

    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
