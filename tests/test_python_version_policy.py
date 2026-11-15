"""Check the supported Python floor without building the C++ application."""
from pathlib import Path
import tomllib
import unittest

ROOT = Path(__file__).resolve().parents[1]


class PythonVersionPolicyTests(unittest.TestCase):
    def test_wheel_metadata_requires_supported_python(self):
        metadata = tomllib.loads((ROOT / "pyproject.toml").read_text())
        self.assertEqual(metadata["project"]["requires-python"], ">=3.11")
        self.assertNotIn("Programming Language :: Python :: 3.10",
                         metadata["project"]["classifiers"])
        targets = metadata["tool"]["cibuildwheel"]["build"].split()
        self.assertNotIn("cp310-*", targets)
        self.assertIn("cp311-*", targets)
        self.assertEqual(metadata["project"]["dependencies"], ["ipython>=9,<10"])


if __name__ == "__main__":
    unittest.main()
