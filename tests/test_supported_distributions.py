"""Exercise the package workflows against the November 2026 lifecycle audit."""

import contextlib
import io
import json
import os
from pathlib import Path
import re
import tempfile
import types
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]


class DistributionTests(unittest.TestCase):
    def matrix(self, workflow, month):
        source = (ROOT / '.github/workflows' / workflow).read_text()
        script = source.split("python3 - <<'EOF'\n", 1)[1].split('\n          EOF', 1)[0]
        script = '\n'.join(line[10:] for line in script.splitlines())
        clock = types.SimpleNamespace(
            datetime=types.SimpleNamespace(utcnow=lambda: types.SimpleNamespace(
                strftime=lambda _format: month)))
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'output'
            with patch.dict(os.environ, GITHUB_OUTPUT=str(output)), \
                    patch.dict('sys.modules', datetime=clock), \
                    contextlib.redirect_stdout(io.StringIO()):
                previous = Path.cwd()
                try:
                    os.chdir(ROOT)
                    exec(compile(script, workflow, 'exec'), {})
                finally:
                    os.chdir(previous)
            rows = json.loads(output.read_text().split('=', 1)[1])['include']
            return [dict(distro=x.get('family', x['distro']),
                         version=x.get('distro_version', x.get('version')),
                         container=x.get('docker_image', x.get('container'))) for x in rows]

    def test_active_targets_at_audit_date(self):
        expected = {
            'build-debian-packages.yml': {
                ('ubuntu', '22.04'), ('ubuntu', '24.04'), ('ubuntu', '26.04'),
                ('ubuntu', '26.10'), ('debian', '12'), ('debian', '13')},
            'build-rpm-packages.yml': {
                ('fedora', '43'), ('fedora', '44'), ('fedora', '45'),
                ('el', '9'), ('el', '10')},
        }
        for workflow, targets in expected.items():
            with self.subTest(workflow=workflow):
                matrix = self.matrix(workflow, '2026-11')
                self.assertEqual({(x['distro'], x['version']) for x in matrix}, targets)
                self.assertTrue(all(re.fullmatch(r'[^@]+@sha256:[0-9a-f]{64}',
                                                x['container']) for x in matrix))

    def test_retired_targets_removed_from_declaration(self):
        rows = json.loads((ROOT / 'supported-distributions.json').read_text())['distributions']
        self.assertTrue(all(x['eol_date'] >= '2026-11' for x in rows))
        keys = [(x['family'], x['version']) for x in rows]
        self.assertEqual(len(keys), len(set(keys)))

    def test_corrected_expiry_months(self):
        deb = self.matrix('build-debian-packages.yml', '2027-06')
        self.assertNotIn(('ubuntu', '22.04'), {(x['distro'], x['version']) for x in deb})
        rpm = self.matrix('build-rpm-packages.yml', '2027-01')
        self.assertNotIn(('fedora', '43'), {(x['distro'], x['version']) for x in rpm})


if __name__ == '__main__':
    unittest.main()
