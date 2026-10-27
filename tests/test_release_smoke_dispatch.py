"""Exercise release build orchestration with completed GitHub run fixtures."""

import os
from pathlib import Path
import subprocess
import tempfile
import unittest


class DispatchTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        version = subprocess.run(
            ["bash", "-c", '[[ ${BASH_VERSINFO[0]} -ge 4 ]]'], check=False,
        )
        if version.returncode:
            raise unittest.SkipTest("release build tracking requires Bash 4 or newer")

    def run_dispatch(self, conclusion):
        root = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as directory:
            gh = Path(directory) / "gh"
            gh.write_text('''#!/usr/bin/env python3
import os
import sys
args = sys.argv[1:]
if args[:2] == ["workflow", "run"]:
    pass
elif args[:2] == ["run", "list"]:
    workflows = ["build-appimage.yml", "build-debian-packages.yml", "build-rpm-packages.yml",
                 "build-windows-native.yml", "build-macos-release.yml"]
    print(1000 + workflows.index(args[args.index("--workflow") + 1]))
elif args[:2] == ["run", "view"]:
    print("completed\\t" + os.environ["TEST_CONCLUSION"] + "\\thttps://example.test/run\\tBuild")
elif args[:2] == ["run", "cancel"]:
    print("CANCELLED " + args[2])
else:
    sys.exit("unexpected gh invocation: " + repr(args))
''')
            gh.chmod(0o755)
            result = subprocess.run(
                ["bash", str(root / "scripts/release-smoke/trigger-release-builds.sh"),
                 "--repo", "owner/repo", "--ref", "master", "--wait"],
                env={**os.environ, "PATH": directory + os.pathsep + os.environ["PATH"],
                     "TEST_CONCLUSION": conclusion},
                capture_output=True, text=True, timeout=15,
            )
            return result

    def test_all_successful_builds_return_success(self):
        result = self.run_dispatch("success")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("All dispatched workflows completed successfully", result.stdout)
        self.assertNotIn("CANCELLED", result.stdout)

    def test_failed_build_cancels_remaining_builds(self):
        result = self.run_dispatch("failure")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("failed with conclusion: failure", result.stderr)
        self.assertNotIn("CANCELLED 1000", result.stdout)
        for run_id in range(1001, 1005):
            self.assertIn(f"CANCELLED {run_id}", result.stdout)


if __name__ == "__main__":
    unittest.main()
