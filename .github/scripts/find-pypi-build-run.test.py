#!/usr/bin/env python3
"""Unit tests for find-pypi-build-run.py."""

from __future__ import annotations

import importlib.util
import pathlib
import re
import unittest

SCRIPT = pathlib.Path(__file__).resolve().parent / "find-pypi-build-run.py"
spec = importlib.util.spec_from_file_location("find_pypi_build_run", SCRIPT)
assert spec is not None and spec.loader is not None
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

TAG = "v1.2.2"
SHA = "a" * 40
GOOD_ARTIFACTS = [
    {"name": "pypi-sdist", "expired": False},
    {"name": "pypi-wheels-ubuntu-24.04", "expired": False},
    {"name": "pypi-wheels-windows-2022", "expired": False},
]


def run(run_id: int, **overrides):
    base = {
        "id": run_id,
        "event": "release",
        "head_branch": TAG,
        "head_sha": SHA,
        "status": "completed",
        "conclusion": "success",
        "created_at": f"2026-11-25T10:00:{run_id % 60:02d}Z",
    }
    base.update(overrides)
    return base


class SelectRunTest(unittest.TestCase):
    def select(self, runs, artifacts):
        return mod.select_run(runs, lambda rid: artifacts.get(rid, []), TAG, SHA)

    def test_returns_matching_run(self) -> None:
        self.assertEqual(self.select([run(1)], {1: GOOD_ARTIFACTS}), 1)

    def test_prefers_newest_usable_run(self) -> None:
        runs = [run(1), run(2)]
        self.assertEqual(self.select(runs, {1: GOOD_ARTIFACTS, 2: GOOD_ARTIFACTS}), 2)

    def test_skips_newer_run_whose_artifacts_expired(self) -> None:
        expired = [dict(a, expired=True) for a in GOOD_ARTIFACTS]
        runs = [run(1), run(2)]
        self.assertEqual(self.select(runs, {1: GOOD_ARTIFACTS, 2: expired}), 1)

    def test_ignores_other_tag(self) -> None:
        self.assertIsNone(self.select([run(1, head_branch="v1.2.1")], {1: GOOD_ARTIFACTS}))

    def test_ignores_other_commit(self) -> None:
        self.assertIsNone(self.select([run(1, head_sha="b" * 40)], {1: GOOD_ARTIFACTS}))

    def test_ignores_failed_and_in_progress_runs(self) -> None:
        runs = [run(1, conclusion="failure"), run(2, status="in_progress", conclusion=None)]
        self.assertIsNone(self.select(runs, {1: GOOD_ARTIFACTS, 2: GOOD_ARTIFACTS}))

    def test_ignores_non_release_events(self) -> None:
        self.assertIsNone(self.select([run(1, event="workflow_dispatch")], {1: GOOD_ARTIFACTS}))

    def test_requires_sdist_and_wheels(self) -> None:
        only_sdist = [{"name": "pypi-sdist", "expired": False}]
        only_wheels = [{"name": "pypi-wheels-macos-15", "expired": False}]
        self.assertIsNone(self.select([run(1)], {1: only_sdist}))
        self.assertIsNone(self.select([run(1)], {1: only_wheels}))

    def test_no_runs(self) -> None:
        self.assertIsNone(self.select([], {}))


class WorkflowWiringTest(unittest.TestCase):
    """The workflow must build on `published` and upload on `released`."""

    @classmethod
    def setUpClass(cls) -> None:
        path = SCRIPT.parents[1] / "workflows" / "publish-to-pypi.yml"
        cls.text = path.read_text(encoding="utf-8")

    def job(self, name: str) -> str:
        start = self.text.index(f"\n  {name}:\n")
        # Jobs are separated by a blank line followed by a two-space-indented key.
        m = re.search(r"\n\n  [a-z0-9-]+:\n", self.text[start + 1 :])
        end = start + 1 + m.start() if m else len(self.text)
        return self.text[start:end]

    def test_release_types(self) -> None:
        self.assertIn("types: [published, released]", self.text)

    def test_production_upload_only_on_released_non_prerelease(self) -> None:
        job = self.job("publish-to-pypi")
        self.assertIn("github.event.action == 'released'", job)
        self.assertIn("!github.event.release.prerelease", job)
        self.assertIn("find-build-run", job)

    def test_testpypi_only_on_published(self) -> None:
        job = self.job("publish-to-test-pypi")
        self.assertIn("github.event.action == 'published'", job)

    def test_builds_skipped_on_released_when_run_found(self) -> None:
        for name in ("build-wheels", "build-sdist"):
            job = self.job(name)
            self.assertIn("needs.find-build-run.outputs.run_id == ''", job, name)


if __name__ == "__main__":
    unittest.main()
