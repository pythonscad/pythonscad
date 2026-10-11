#!/usr/bin/env python3
"""Regression tests for scripts/ci/compare-git-release.sh.

Covers the delayed-release compare used by the Manifold submodule
update workflow: empty/null tags, matching peeled tag SHAs, and
mismatched SHAs (including annotated tags, which must not be
compared as tag names).
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def _is_wsl_bash_stub(path: Path) -> bool:
    """True for the Windows System32/SysWOW64 WSL launcher bash.exe.

    GitHub-hosted Windows runners ship that stub with no distro installed.
    CreateProcess also prefers System32 for bare ``bash``, so callers must
    spawn an absolute path to a real Bash (MSYS2 / Git for Windows).
    """
    try:
        resolved = path.resolve()
    except OSError:
        return False
    if resolved.name.lower() not in ("bash", "bash.exe"):
        return False
    return any(part.lower() in ("system32", "syswow64") for part in resolved.parts)


def _bash_works(path: str) -> bool:
    try:
        proc = subprocess.run(
            [path, "--version"],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return proc.returncode == 0 and "GNU bash" in (proc.stdout + proc.stderr)


def _find_bash() -> str | None:
    candidates: list[Path] = []
    for name in ("bash", "bash.exe"):
        found = shutil.which(name)
        if found:
            candidates.append(Path(found))

    if sys.platform == "win32":
        msystem_prefix = os.environ.get("MSYSTEM_PREFIX")
        if msystem_prefix:
            candidates.append(Path(msystem_prefix).parent / "usr" / "bin" / "bash.exe")
        program_files = Path(os.environ.get("ProgramFiles", r"C:\Program Files"))
        candidates.extend(
            [
                Path(r"C:\msys64\usr\bin\bash.exe"),
                program_files / "Git" / "bin" / "bash.exe",
                program_files / "Git" / "usr" / "bin" / "bash.exe",
            ]
        )
    else:
        candidates.extend([Path("/usr/bin/bash"), Path("/bin/bash")])

    seen: set[str] = set()
    for candidate in candidates:
        try:
            resolved = candidate.resolve(strict=True)
        except OSError:
            continue
        key = str(resolved).lower()
        if key in seen or _is_wsl_bash_stub(resolved):
            continue
        seen.add(key)
        path = str(resolved)
        if _bash_works(path):
            return path
    return None


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


def _compare(bash: str, script: Path, repo: Path, tag: str, current: str) -> dict[str, str]:
    proc = subprocess.run(
        [
            bash,
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

    bash = _find_bash()
    if bash is None or shutil.which("git") is None:
        print("SKIP: usable bash or git not on PATH")
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
        # rev-list peels annotated tags without ^{commit}, which cmd.exe
        # and Git-for-Windows wrappers strip to `v1.0.0^commit`.
        peeled = _run_git(repo, "rev-list", "-n", "1", "v1.0.0")
        assert peeled == first, (peeled, first)

        empty = _compare(bash, script, repo, "", first)
        assert empty.get("new_release_available") == "false", empty
        assert "latest_commit" not in empty, empty

        null = _compare(bash, script, repo, "null", first)
        assert null.get("new_release_available") == "false", null

        matching = _compare(bash, script, repo, "v1.0.0", first)
        assert matching.get("new_release_available") == "false", matching
        assert matching.get("latest_commit") == first, matching

        (repo / "README").write_text("two\n", encoding="utf-8")
        _run_git(repo, "add", "README")
        _run_git(repo, "commit", "-m", "second")
        second = _run_git(repo, "rev-parse", "HEAD")
        assert second != first

        mismatched = _compare(bash, script, repo, "v1.0.0", second)
        assert mismatched.get("new_release_available") == "true", mismatched
        assert mismatched.get("latest_commit") == first, mismatched
        assert mismatched.get("release_tag") == "v1.0.0", mismatched
    finally:
        shutil.rmtree(repo, ignore_errors=True)

    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
