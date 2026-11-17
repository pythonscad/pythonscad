#!/usr/bin/env python3
"""Integration tests for scripts/cut-release-branch.sh."""

from __future__ import annotations

import os
import pathlib
import shutil
import subprocess
import tempfile
import unittest


def _run(
    args: list[str],
    *,
    cwd: pathlib.Path,
    check: bool = True,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        args,
        cwd=cwd,
        check=check,
        capture_output=True,
        text=True,
        env=env,
    )


def _git(cwd: pathlib.Path, *args: str, check: bool = True) -> subprocess.CompletedProcess[str]:
    return _run(["git", *args], cwd=cwd, check=check)


class CutReleaseBranchTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.repo_root = pathlib.Path(__file__).resolve().parents[1]
        cls.script = cls.repo_root / "scripts" / "cut-release-branch.sh"

    def setUp(self) -> None:
        self._tmpdir = tempfile.TemporaryDirectory()
        self.root = pathlib.Path(self._tmpdir.name)
        self.bare = self.root / "remote.git"
        self.work = self.root / "work"
        _git(self.root, "init", "--bare", str(self.bare))
        _git(self.root, "clone", str(self.bare), str(self.work))
        _git(self.work, "config", "user.email", "test@example.com")
        _git(self.work, "config", "user.name", "Test")
        # Default branch name varies by git version; normalize to master.
        _git(self.work, "checkout", "-B", "master")

    def tearDown(self) -> None:
        self._tmpdir.cleanup()

    def _write(self, rel: str, content: str) -> None:
        path = self.work / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")

    def _commit_all(self, message: str) -> None:
        _git(self.work, "add", "-A")
        _git(self.work, "commit", "-m", message)

    def _seed_old_tag(self) -> None:
        """Create v1.1.2 without maintenance-branch automation."""
        self._write("VERSION.txt", "1.1.2\n")
        self._write(
            ".github/workflows/release-please.yml",
            "on:\n  push:\n    branches:\n      - master\n",
        )
        self._write(
            ".release-please-config.json",
            '{"packages":{".":{"release-type":"simple"}}}\n',
        )
        self._commit_all("chore: seed old release")
        _git(self.work, "tag", "v1.1.2")
        _git(self.work, "push", "origin", "HEAD:master")
        _git(self.work, "push", "origin", "v1.1.2")

    def _seed_master_automation(self) -> None:
        """Add the automation files that cut-release-branch bootstraps from master."""
        for rel in (
            ".github/workflows/release-please.yml",
            ".release-please-config.json",
            ".release-please-config.release.json",
            ".github/scripts/is-latest-semver-tag.sh",
            ".github/scripts/is-latest-semver-tag.test.py",
            ".github/workflows/deploy-website.yml",
            ".github/workflows/build-wasm.yml",
            "scripts/cut-release-branch.sh",
            "scripts/cut-release-branch.test.py",
        ):
            src = self.repo_root / rel
            dst = self.work / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
        # Make the copied shell scripts executable in the fixture tree.
        (self.work / ".github/scripts/is-latest-semver-tag.sh").chmod(0o755)
        (self.work / "scripts/cut-release-branch.sh").chmod(0o755)
        self._commit_all("ci: add maintenance-branch automation")
        _git(self.work, "push", "origin", "HEAD:master")

    def _run_script(self, *args: str, check: bool = True) -> subprocess.CompletedProcess[str]:
        env = os.environ.copy()
        env["GIT_AUTHOR_NAME"] = "Test"
        env["GIT_AUTHOR_EMAIL"] = "test@example.com"
        env["GIT_COMMITTER_NAME"] = "Test"
        env["GIT_COMMITTER_EMAIL"] = "test@example.com"
        return _run(
            [str(self.script), *args, "--remote", "origin"],
            cwd=self.work,
            check=check,
            env=env,
        )

    def test_rejects_dirty_worktree(self) -> None:
        self._seed_old_tag()
        self._seed_master_automation()
        (self.work / "VERSION.txt").write_text("dirty\n", encoding="utf-8")
        result = self._run_script("1.1.2", check=False)
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("uncommitted changes", result.stderr)
        # Must not create the branch.
        self.assertFalse((self.work / ".git" / "refs" / "heads" / "release" / "1.1").exists())

    def test_rejects_untracked_bootstrap_collision(self) -> None:
        self._seed_old_tag()
        self._seed_master_automation()
        # Start from the old tag tree (bootstrap paths absent) with an
        # untracked file that checkout-from-master would clobber.
        _git(self.work, "checkout", "-B", "from-old", "v1.1.2")
        collision = self.work / ".release-please-config.release.json"
        self.assertFalse(collision.exists())
        collision.write_text("{}\n", encoding="utf-8")
        result = self._run_script("1.1.2", check=False)
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("untracked file would be overwritten", result.stderr)
        self.assertTrue(collision.exists())
        self.assertEqual(collision.read_text(encoding="utf-8"), "{}\n")
        self.assertFalse((self.work / ".git" / "refs" / "heads" / "release" / "1.1").exists())

    def test_bootstraps_older_tag_and_restores_checkout(self) -> None:
        self._seed_old_tag()
        self._seed_master_automation()
        start = _git(self.work, "rev-parse", "--abbrev-ref", "HEAD").stdout.strip()
        self.assertEqual(start, "master")

        result = self._run_script("1.1.2")
        self.assertIn("predates maintenance-branch automation", result.stdout)

        # Back on the starting branch.
        self.assertEqual(
            _git(self.work, "rev-parse", "--abbrev-ref", "HEAD").stdout.strip(),
            "master",
        )

        # Branch exists locally and on the bare remote.
        _git(self.work, "rev-parse", "release/1.1")
        _git(self.work, "ls-remote", "--exit-code", "origin", "refs/heads/release/1.1")

        # Bootstrap commit brought in the release config and release/** trigger.
        tip = _git(self.work, "rev-parse", "release/1.1").stdout.strip()
        for rel in (
            ".release-please-config.release.json",
            "scripts/cut-release-branch.sh",
            "scripts/cut-release-branch.test.py",
            ".github/scripts/is-latest-semver-tag.test.py",
        ):
            self.assertEqual(
                _git(self.work, "cat-file", "-e", f"{tip}:{rel}", check=False).returncode,
                0,
                f"missing bootstrap path {rel}",
            )
        workflow = _git(self.work, "show", f"{tip}:.github/workflows/release-please.yml").stdout
        self.assertIn("release/**", workflow)
        log = _git(self.work, "log", "-1", "--format=%s", tip).stdout.strip()
        self.assertEqual(log, "ci: bootstrap release maintenance-branch automation")

    def test_compatible_tag_skips_bootstrap(self) -> None:
        self._seed_old_tag()
        self._seed_master_automation()
        # Retag current master (with automation) as v1.2.0.
        _git(self.work, "tag", "v1.2.0")
        _git(self.work, "push", "origin", "v1.2.0")

        result = self._run_script("1.2.0")
        self.assertNotIn("predates maintenance-branch automation", result.stdout)
        tip = _git(self.work, "rev-parse", "release/1.2").stdout.strip()
        # Tip should be exactly the tag (no bootstrap commit).
        tag = _git(self.work, "rev-parse", "v1.2.0").stdout.strip()
        self.assertEqual(tip, tag)

    def test_dry_run_does_not_create_branch(self) -> None:
        self._seed_old_tag()
        self._seed_master_automation()
        result = self._run_script("1.1.2", "--dry-run")
        self.assertIn("(dry-run)", result.stdout)
        self.assertNotEqual(
            _git(self.work, "show-ref", "--verify", "--quiet", "refs/heads/release/1.1", check=False).returncode,
            0,
        )


if __name__ == "__main__":
    unittest.main()
