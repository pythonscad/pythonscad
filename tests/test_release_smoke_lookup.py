"""Prevent release smoke tests from silently selecting packages from another commit."""

import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
COMMIT = "a" * 40
OTHER = "b" * 40


class LookupTests(unittest.TestCase):
    def invoke(self, command, runs, actual=COMMIT, gh_exit=0):
        with tempfile.TemporaryDirectory() as directory:
            gh = Path(directory) / "gh"
            gh.write_text('''#!/usr/bin/env python3
import json, os, sys
args = sys.argv[1:]
if int(os.environ['GH_EXIT']):
    sys.exit(int(os.environ['GH_EXIT']))
if args[:2] == ['run', 'list']:
    # Exact SHA filtering must happen at the server, before --limit is applied.
    assert args[args.index('--commit') + 1] == os.environ['EXPECTED_SHA'], args
    runs = [r for r in json.loads(os.environ['RUNS']) if r['headSha'] == os.environ['EXPECTED_SHA']]
    if '--branch' in args:
        runs = [r for r in runs if r['headBranch'] == args[args.index('--branch') + 1]]
    print(runs[0]['databaseId'] if runs else '')
elif args[:2] == ['run', 'view']:
    print(os.environ['ACTUAL_SHA'])
else:
    sys.exit('unexpected gh call: ' + repr(args))
''')
            gh.chmod(0o755)
            return subprocess.run(
                ["bash", "-euo", "pipefail", "-c",
                 'source "$1/scripts/release-smoke/common-smoke.sh"; ' + command,
                 "bash", str(ROOT)],
                env={**os.environ, "PATH": directory + os.pathsep + os.environ["PATH"],
                     "RUNS": json.dumps(runs), "EXPECTED_SHA": COMMIT,
                     "ACTUAL_SHA": actual, "GH_EXIT": str(gh_exit),
                     "RELEASE_SMOKE_COMMIT": COMMIT},
                capture_output=True, text=True, timeout=10,
            )

    def test_selects_exact_commit_even_with_many_newer_builds(self):
        runs = [{"databaseId": n, "headSha": OTHER, "headBranch": "master"} for n in range(100)]
        runs += [{"databaseId": 1234, "headSha": COMMIT, "headBranch": "master"}]
        result = self.invoke('rs_latest_successful_run owner/repo build.yml master', runs)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), "1234")

    def test_other_commit_or_branch_never_falls_back(self):
        for sha, branch in [(OTHER, "master"), (COMMIT, "another-branch")]:
            with self.subTest(sha=sha, branch=branch):
                result = self.invoke('rs_latest_successful_run owner/repo build.yml master',
                                     [{"databaseId": 1, "headSha": sha, "headBranch": branch}])
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("trigger fresh release builds", result.stderr)
                self.assertEqual(result.stdout, "")

    def test_explicit_run_from_wrong_commit_is_rejected(self):
        result = self.invoke('rs_validate_run_commit owner/repo 123 "$RELEASE_SMOKE_COMMIT"',
                             [], actual=OTHER)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("expected " + COMMIT, result.stderr)

    def test_github_failure_is_not_treated_as_a_build(self):
        result = self.invoke('rs_latest_successful_run owner/repo build.yml master', [], gh_exit=1)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")


if __name__ == "__main__":
    unittest.main()
