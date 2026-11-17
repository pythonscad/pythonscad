#!/usr/bin/env python3
"""Unit tests for is-latest-semver-tag.sh."""

from __future__ import annotations

import pathlib
import subprocess
import tempfile
import unittest


class IsLatestSemverTagTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.repo_root = pathlib.Path(__file__).resolve().parents[2]
        cls.script = cls.repo_root / ".github/scripts/is-latest-semver-tag.sh"

    def _run(self, cwd: pathlib.Path, tag: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [str(self.script), tag],
            check=False,
            capture_output=True,
            text=True,
            cwd=cwd,
        )

    def test_latest_and_older_tags(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            repo = pathlib.Path(temp_dir)
            subprocess.run(["git", "init"], check=True, cwd=repo, capture_output=True)
            subprocess.run(
                ["git", "config", "user.email", "test@example.com"],
                check=True,
                cwd=repo,
                capture_output=True,
            )
            subprocess.run(
                ["git", "config", "user.name", "Test"],
                check=True,
                cwd=repo,
                capture_output=True,
            )
            (repo / "README").write_text("x\n", encoding="utf-8")
            subprocess.run(["git", "add", "README"], check=True, cwd=repo, capture_output=True)
            subprocess.run(
                ["git", "commit", "-m", "init"],
                check=True,
                cwd=repo,
                capture_output=True,
            )
            for tag in ("v1.0.0", "v1.1.0", "v1.1.2"):
                subprocess.run(
                    ["git", "tag", tag],
                    check=True,
                    cwd=repo,
                    capture_output=True,
                )

            latest = self._run(repo, "v1.1.2")
            self.assertEqual(latest.returncode, 0, latest.stdout + latest.stderr)
            self.assertIn("is the latest semver tag", latest.stdout)

            older = self._run(repo, "v1.0.0")
            self.assertEqual(older.returncode, 1, older.stdout + older.stderr)
            self.assertIn("is not the latest semver tag", older.stdout)

            bare = self._run(repo, "1.1.2")
            self.assertEqual(bare.returncode, 0, bare.stdout + bare.stderr)

    def test_rejects_non_semver(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            repo = pathlib.Path(temp_dir)
            subprocess.run(["git", "init"], check=True, cwd=repo, capture_output=True)
            result = self._run(repo, "v1.2.0-rc.1")
            self.assertEqual(result.returncode, 1, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
