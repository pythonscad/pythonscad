"""Resolve the newest PyQt6 sdist release matching a given Qt
major.minor version, without triggering pip's build-isolation
metadata hook (which would otherwise kick off a full Qt6 module
build just to read package metadata).

Usage: resolve_pyqt6_version.py <major.minor>
Prints: "<version> <sdist_url>"
"""
import json
import sys
import urllib.request

from packaging.version import Version


def main() -> None:
    if len(sys.argv) != 2:
        print("Usage: resolve_pyqt6_version.py <major.minor>", file=sys.stderr)
        sys.exit(1)

    major_minor = sys.argv[1]
    with urllib.request.urlopen("https://pypi.org/pypi/PyQt6/json") as r:
        data = json.load(r)

    candidates = [Version(v) for v, files in data["releases"].items() if files]
    matching = [v for v in candidates if f"{v.major}.{v.minor}" == major_minor]
    if not matching:
        sys.exit(f"Keine PyQt6-Version fuer Qt {major_minor} gefunden")

    best = max(matching)
    sdist = next(f for f in data["releases"][str(best)] if f["packagetype"] == "sdist")
    print(best, sdist["url"])


if __name__ == "__main__":
    main()
