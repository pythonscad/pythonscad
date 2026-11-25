#!/usr/bin/env python3
"""Find the "Publish to PyPI" run that already built the wheels/sdist for a tag.

publish-to-pypi.yml builds the wheels and sdist when a release is *published*
(as a pre-release) and uploads them to PyPI when the release is later
*released* (promoted after the Windows installers were signed). The `released`
run uses this script to look up the build run of the same tag so it can upload
exactly those artifacts instead of rebuilding.

Prints the workflow run id, or nothing if no usable run exists (in which case
the caller falls back to rebuilding).

Usage: find-pypi-build-run.py --repo OWNER/REPO --tag vX.Y.Z --sha <tag commit>
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from typing import Any, Callable

WORKFLOW_FILE = "publish-to-pypi.yml"


def select_run(
    runs: list[dict[str, Any]],
    artifacts_for: Callable[[int], list[dict[str, Any]]],
    tag: str,
    sha: str,
) -> int | None:
    """Return the id of the newest usable build run for `tag`, or None.

    A run is usable when it is a completed, successful `release` run for the
    same tag *and* the same commit, and it still has a non-expired sdist
    artifact and at least one non-expired wheel artifact.
    """
    candidates = [
        r
        for r in runs
        if r.get("event") == "release"
        and r.get("head_branch") == tag
        and r.get("head_sha") == sha
        and r.get("status") == "completed"
        and r.get("conclusion") == "success"
    ]
    # Newest first (ISO 8601 timestamps sort lexicographically).
    candidates.sort(key=lambda r: r.get("created_at", ""), reverse=True)

    for run in candidates:
        live = [a["name"] for a in artifacts_for(run["id"]) if not a.get("expired")]
        if "pypi-sdist" in live and any(n.startswith("pypi-wheels-") for n in live):
            return int(run["id"])
    return None


def _gh_api(path: str) -> Any:
    out = subprocess.run(
        ["gh", "api", path],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    return json.loads(out)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--sha", required=True)
    args = parser.parse_args()

    try:
        runs = _gh_api(
            f"repos/{args.repo}/actions/workflows/{WORKFLOW_FILE}/runs"
            "?event=release&status=success&per_page=100"
        ).get("workflow_runs", [])

        def artifacts_for(run_id: int) -> list[dict[str, Any]]:
            return _gh_api(f"repos/{args.repo}/actions/runs/{run_id}/artifacts?per_page=100").get(
                "artifacts", []
            )

        run_id = select_run(runs, artifacts_for, args.tag, args.sha)
    except (subprocess.CalledProcessError, json.JSONDecodeError, OSError) as exc:
        # Never fail the release: no run id means "rebuild".
        print(f"warning: could not look up an earlier build run: {exc}", file=sys.stderr)
        run_id = None

    if run_id is not None:
        print(run_id)
    return 0


if __name__ == "__main__":
    sys.exit(main())
