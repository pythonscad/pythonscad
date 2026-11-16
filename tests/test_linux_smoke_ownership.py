"""Run package-container shell bodies with host-side command fixtures."""

import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


@unittest.skipIf(os.name == "nt", "Linux container shell fixtures require Unix")
class OutputOwnershipTests(unittest.TestCase):
    def test_startup_check_can_avoid_opengl(self):
        root = Path(__file__).resolve().parents[1]
        fixture = r'''
set -e
source "$FIXTURE_COMMON"
rs_write_common_test_inputs() { :; }
rs_log() { :; }
rs_run_logged() { printf '%s\n' "${@: -1}"; return 77; }
rs_smoke_binary "$FIXTURE_EXE" fixture /tmp 0 "$FIXTURE_STARTUP"
'''
        for option in ("--info", "--version"):
            result = subprocess.run(
                ["bash", "-c", fixture], capture_output=True, text=True,
                env={**os.environ, "FIXTURE_COMMON": str(root / "scripts/release-smoke/common-smoke.sh"),
                     "FIXTURE_STARTUP": option, "FIXTURE_EXE": sys.executable}, timeout=10)
            self.assertEqual(result.returncode, 77, result.stderr)
            self.assertEqual(result.stdout.strip(), option)

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
  if [[ "${QT_QPA_PLATFORM:-}" == offscreen ]]; then
    [[ "$*" == *"--version"* ]] || return 99
  fi
  printf '%s\n' "${QT_QPA_PLATFORM:-}" > "$FIXTURE_OUTPUT/platform"
  mkdir -p "$FIXTURE_OUTPUT/nested"
  touch "$FIXTURE_OUTPUT/nested/result.stl"
  return "$FIXTURE_STATUS"
}
chown() {
  printf '%s\n' "$@" > "$FIXTURE_LOG"
  command chown "$@"
}
bash() {
  if [[ "${QT_QPA_PLATFORM:-}" == offscreen ]]; then
    xvfb-run "$@"
  else
    command bash "$@"
  fi
}
export -f apt-get dnf xvfb-run chown bash
docker() {
  local body="${@: -1}"
  while (($#)); do
    if [[ "$1" == -e ]]; then export "$2"; shift; fi
    shift
  done
  body=${body//\/tmp\/pythonscad-smoke/$FIXTURE_OUTPUT}
  command bash -c "$body"
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
                             "FIXTURE_LOG": str(log), "FIXTURE_STATUS": str(status),
                             "QT_QPA_PLATFORM": "offscreen"},
                        capture_output=True, text=True, timeout=10,
                    )
                    self.assertEqual(result.returncode, status, result.stderr)
                    self.assertEqual((output / "platform").read_text().strip(),
                                     "offscreen" if package_type == "rpm" else "xcb")
                    self.assertTrue(log.is_file(), "container did not restore output ownership")
                    self.assertEqual(log.read_text().splitlines(),
                                     ["-R", f"{os.getuid()}:{os.getgid()}", str(output)])
                    self.assertTrue((output / "nested/result.stl").is_file())


if __name__ == "__main__":
    unittest.main()
