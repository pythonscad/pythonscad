#!/usr/bin/env python3
"""Regression test for the release flow (#1281).

New releases must not become public / "latest" the moment the release PR is
merged. Release Please therefore creates them as *drafts* (a "prerelease"
setting is ignored for 1.x versions), release-please.yml publishes the draft as
a prerelease that is not "latest", and every workflow that builds or deploys
release artifacts reacts to the `published` event (the `created` event is not
sent for drafts).
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

RELEASE_PLEASE_CONFIGS = [
    ".release-please-config.json",
    ".release-please-config.release.json",
]

# Workflows that must run once a release has been published (as a prerelease).
RELEASE_EVENT_WORKFLOWS = [
    "build-appimage.yml",
    "build-debian-packages.yml",
    "build-macos-release.yml",
    "build-rpm-packages.yml",
    "build-wasm.yml",
    "build-windows-native.yml",
    "deploy-website.yml",
    "notify-signing.yml",
]


def fail(msg: str) -> int:
    print(f"FAIL: {msg}", file=sys.stderr)
    return 1


def main() -> int:
    for name in RELEASE_PLEASE_CONFIGS:
        cfg = json.loads((ROOT / name).read_text(encoding="utf-8"))
        pkg = cfg["packages"]["."]
        if pkg.get("draft") is not True:
            return fail(f'{name}: "draft" must be true so releases are not public at creation')
        if pkg.get("force-tag-creation") is not True:
            return fail(
                f'{name}: "force-tag-creation" must be true, otherwise Release Please '
                "cannot find the previous release while the tag does not exist yet"
            )

    for name in RELEASE_EVENT_WORKFLOWS:
        text = (ROOT / ".github" / "workflows" / name).read_text(encoding="utf-8")
        m = re.search(r"^  release:\n    types: \[([^\]]*)\]", text, re.MULTILINE)
        if not m:
            return fail(f"{name}: no `release:` trigger with explicit types found")
        types = [t.strip() for t in m.group(1).split(",")]
        if types != ["published"]:
            return fail(
                f"{name}: release trigger types must be [published] (the `created` event "
                f"is not sent for draft releases), got {types}"
            )

    rp = (ROOT / ".github" / "workflows" / "release-please.yml").read_text(encoding="utf-8")
    for needle in ("draft: false", "prerelease: true", "make_latest: 'false'"):
        if needle not in rp:
            return fail(f"release-please.yml must publish the draft release with `{needle}`")

    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
