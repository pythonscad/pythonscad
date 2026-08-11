"""Resolve the pinned PyQt6 sdist for a given Qt major.minor version.

Usage: resolve_pyqt6_version.py <major.minor> [--version-only]
Prints: "<version> <sdist_url>"
"""
import json
import sys
import urllib.request
from pathlib import Path


def main() -> None:
    if len(sys.argv) not in (2, 3) or (
        len(sys.argv) == 3 and sys.argv[2] != "--version-only"
    ):
        print(
            "Usage: resolve_pyqt6_version.py <major.minor> [--version-only]",
            file=sys.stderr,
        )
        sys.exit(1)

    major_minor = sys.argv[1]
    versions_path = Path(__file__).with_name("pyqt6-versions.json")
    versions = json.loads(versions_path.read_text())
    try:
        pinned_version = versions[major_minor]
    except KeyError:
        sys.exit(
            f"No pinned PyQt6 version for Qt {major_minor}; "
            f"update {versions_path.name}"
        )

    if len(sys.argv) == 3:
        print(pinned_version)
        return

    with urllib.request.urlopen("https://pypi.org/pypi/PyQt6/json") as r:
        data = json.load(r)

    files = data["releases"].get(pinned_version, [])
    try:
        sdist = next(f for f in files if f["packagetype"] == "sdist")
    except StopIteration:
        sys.exit(f"PyQt6 {pinned_version} has no sdist on PyPI")
    print(pinned_version, sdist["url"])


if __name__ == "__main__":
    main()
