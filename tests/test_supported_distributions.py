"""Exercise the package workflows against the November 2026 lifecycle audit."""

import configparser
import contextlib
import io
import json
import os
from pathlib import Path
import re
import subprocess
import sys
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
        december = self.matrix('build-rpm-packages.yml', '2026-12')
        self.assertIn(('fedora', '43'), {(x['distro'], x['version']) for x in december})
        rpm = self.matrix('build-rpm-packages.yml', '2027-01')
        self.assertNotIn(('fedora', '43'), {(x['distro'], x['version']) for x in rpm})


    def test_ubuntu_2610_dependency_names(self):
        def packages(version):
            result = subprocess.run(
                [sys.executable, str(ROOT / 'scripts/get-dependencies.py'), '--distro', 'ubuntu',
                 '--version', version, '--profile', 'pythonscad-qt6', '--list'],
                check=True, capture_output=True, text=True)
            return set(result.stdout.splitlines()[3:])
        current = packages('26.10')
        self.assertIn('qt6-5compat-dev', current)
        self.assertNotIn('libboost-system-dev', current)
        self.assertNotIn('libqt6core5compat6-dev', current)
        self.assertIn('libqt6core5compat6-dev', packages('24.04'))
        self.assertIn('catch2', packages('24.04'))
        self.assertIn('libqt6core5compat6-dev | qt6-5compat-dev',
                      (ROOT / 'debian/control').read_text())


    @unittest.skipIf(os.name == "nt", "YUM repository generator requires Unix Bash")
    def test_el_client_repository_matches_release(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, repo, tools = (root / name for name in ('rpms', 'repo', 'bin'))
            for path in (source, repo, tools):
                path.mkdir()
            for version in ('9', '10'):
                (source / f'pythonscad-1.1.2-1.el{version}.x86_64.rpm').touch()
            stub = tools / 'createrepo_c'
            stub.write_text('#!/bin/sh\nexit 0\n')
            stub.chmod(0o755)
            subprocess.run(['bash', str(ROOT / 'scripts/update-yum-repo.sh'), str(source)],
                           env={**os.environ, 'REPO_DIR': str(repo), 'GPG_KEY': '',
                                'PATH': str(tools) + os.pathsep + os.environ['PATH']},
                           check=True, capture_output=True, text=True)
            config = configparser.ConfigParser(interpolation=None)
            config.read(repo / 'pythonscad.repo')
            for version in ('9', '9.4', '10', '10.2'):
                major = version.split('.')[0]
                baseurl = config[f'pythonscad-el{major}']['baseurl']
                self.assertNotIn('$releasever', baseurl)
                url = baseurl.replace('$basearch', 'x86_64')
                self.assertTrue((repo / url.split('/yum/', 1)[1]).is_dir())
            instructions = (repo / 'index.html').read_text()
            self.assertIn('config-manager --set-disabled pythonscad', instructions)
            self.assertIn('config-manager --set-enabled pythonscad-el9\n', instructions)
            self.assertIn('<code>pythonscad-el10</code>', instructions)


if __name__ == '__main__':
    unittest.main()
