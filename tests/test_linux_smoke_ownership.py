"""Run package-container shell bodies with host-side command fixtures."""

import os
from pathlib import Path
import subprocess
import tempfile
import unittest


class OutputOwnershipTests(unittest.TestCase):
    def test_ownership_restored_on_success_and_failure(self):
        root = Path(__file__).resolve().parents[1]
        script = (root / "scripts/release-smoke/test-linux-packages.sh").read_text()
        functions = script.split("run_deb_package_test() {", 1)[1].split("\nfor package in", 1)[0]
        functions = "run_deb_package_test() {" + functions
        fixtures = r'''
set -euo pipefail
ensure_deb_smoke_image() { echo fixture-image; }
ensure_rpm_smoke_image() { echo fixture-image; }
apt-get() { :; }
dnf() { :; }
xvfb-run() {
  mkdir -p "$FIXTURE_OUTPUT/nested"
  touch "$FIXTURE_OUTPUT/nested/result.stl"
  return "$FIXTURE_STATUS"
}
chown() {
  printf '%s\n' "$@" > "$FIXTURE_LOG"
  command chown "$@"
}
export -f apt-get dnf xvfb-run chown
docker() {
  local body="${@: -1}"
  while (($#)); do
    if [[ "$1" == -e ]]; then export "$2"; shift; fi
    shift
  done
  body=${body//\/tmp\/pythonscad-smoke/$FIXTURE_OUTPUT}
  bash -c "$body"
}
SCRIPT_DIR=unused-fixture-directory
'''
        for package_type in ("deb", "rpm"):
            for status in (0, 7):
                with self.subTest(package_type=package_type, status=status), tempfile.TemporaryDirectory() as directory:
                    output = Path(directory) / "output"
                    output.mkdir()
                    log = Path(directory) / "ownership.log"
                    result = subprocess.run(
                        ["bash", "-c", fixtures + functions +
                         f'\nrun_{package_type}_package_test fixture.{package_type} fixture "$FIXTURE_OUTPUT"'],
                        env={**os.environ, "FIXTURE_OUTPUT": str(output),
                             "FIXTURE_LOG": str(log), "FIXTURE_STATUS": str(status)},
                        capture_output=True, text=True, timeout=10,
                    )
                    self.assertEqual(result.returncode, status, result.stderr)
                    self.assertTrue(log.is_file(), "container did not restore output ownership")
                    self.assertEqual(log.read_text().splitlines(),
                                     ["-R", f"{os.getuid()}:{os.getgid()}", str(output)])
                    self.assertTrue((output / "nested/result.stl").is_file())


if __name__ == "__main__":
    unittest.main()
