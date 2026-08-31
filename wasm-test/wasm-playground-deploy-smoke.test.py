#!/usr/bin/env python3

import os
import pathlib
import subprocess
import tempfile
import unittest


class PlaygroundDeploySmokeTest(unittest.TestCase):
    def test_smoke_script_downloads_before_searching(self):
        repo_root = pathlib.Path(__file__).resolve().parents[1]
        script = repo_root / ".github/scripts/wasm-playground-deploy-smoke.sh"
        with tempfile.TemporaryDirectory() as temp_dir:
            fake_curl = pathlib.Path(temp_dir) / "curl"
            fake_curl.write_text(
                """#!/usr/bin/env python3
import pathlib
import sys

args = sys.argv[1:]
if any(arg.startswith("-") and "I" in arg[1:] for arg in args):
    print("HTTP/1.1 200 OK\\r\\nContent-Type: application/wasm\\r\\n")
    sys.exit()

body = b"<html>" + (b"x" * (2 * 1024 * 1024))
if "-o" in args:
    pathlib.Path(args[args.index("-o") + 1]).write_bytes(body)
    sys.exit()

try:
    sys.stdout.buffer.write(body)
    sys.stdout.buffer.flush()
except BrokenPipeError:
    pass
sys.exit(23)
"""
            )
            fake_curl.chmod(0o755)
            env = os.environ.copy()
            env["PATH"] = f"{temp_dir}:{env['PATH']}"
            result = subprocess.run(
                [script, "https://example.test/"],
                check=False,
                capture_output=True,
                text=True,
                env=env,
            )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("Playground smoke check passed.", result.stdout)


if __name__ == "__main__":
    unittest.main()
