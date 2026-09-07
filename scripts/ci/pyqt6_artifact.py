#!/usr/bin/env python3
"""Create and validate identities for immutable PyQt6 binding artifacts."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path


SCHEMA_VERSION = 1
IDENTITY_INPUTS = (
    ".github/workflows/build-pyqt6-bindings.yml",
    "scripts/ci/patch_pyqt6_free_operators.py",
    "scripts/ci/patch_pyqtbuild_win32_pythonlib.py",
    "scripts/ci/pyqt6-build-requirements.txt",
    "scripts/ci/pyqt6-versions.json",
    "scripts/ci/pyqt6_artifact.py",
    "scripts/ci/pyqt6_smoke_test.py",
    "scripts/ci/resolve_pyqt6_version.py",
)


def builder_revision(repo_root: Path) -> str:
    digest = hashlib.sha256()
    for relative_path in IDENTITY_INPUTS:
        path = repo_root / relative_path
        digest.update(relative_path.encode())
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def identity(args: argparse.Namespace) -> dict[str, str | int]:
    return {
        "schema_version": SCHEMA_VERSION,
        "platform": args.platform,
        "profile": args.profile,
        "qt_version": args.qt_version,
        "pyqt_version": args.pyqt_version,
        "python_version": args.python_version,
        "soabi": args.soabi,
        "builder_revision": builder_revision(args.repo_root),
    }


def artifact_id(metadata: dict[str, str | int]) -> str:
    canonical = json.dumps(metadata, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(canonical).hexdigest()[:16]


def descriptor(args: argparse.Namespace) -> dict[str, object]:
    metadata = identity(args)
    identifier = artifact_id(metadata)
    platform = str(metadata["platform"])
    return {
        "identity": metadata,
        "artifact_id": identifier,
        "release_tag": f"pyqt6-bindings-v1-{platform}-{identifier}",
        "asset_name": f"pyqt6-bindings-{platform}-{identifier}.zip",
    }


def write_manifest(args: argparse.Namespace) -> None:
    manifest = descriptor(args)
    manifest["pyqt_version"] = args.pyqt_version
    manifest["built_at"] = datetime.now(timezone.utc).isoformat()
    manifest["built_from_commit"] = os.environ.get("GITHUB_SHA") or git_revision(args.repo_root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def git_revision(repo_root: Path) -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo_root,
        check=False,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip() if result.returncode == 0 else "unknown"


def read_manifest(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise SystemExit(f"{path}: manifest root must be an object")
    return value


def validate(args: argparse.Namespace) -> None:
    actual = read_manifest(args.manifest)
    expected = descriptor(args)
    for key in ("identity", "artifact_id", "release_tag", "asset_name"):
        if actual.get(key) != expected[key]:
            raise SystemExit(
                f"{args.manifest}: {key} mismatch\n"
                f"expected: {json.dumps(expected[key], sort_keys=True)}\n"
                f"actual:   {json.dumps(actual.get(key), sort_keys=True)}"
            )
    print(f"Validated PyQt6 artifact {expected['artifact_id']} for {args.platform}")


def print_field(args: argparse.Namespace) -> None:
    value = descriptor(args)[args.field]
    if not isinstance(value, str):
        raise SystemExit(f"{args.field} is not a scalar field")
    print(value)


def add_identity_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--platform", required=True)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--qt-version", required=True)
    parser.add_argument("--pyqt-version", required=True)
    parser.add_argument("--python-version", required=True)
    parser.add_argument("--soabi", required=True)
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    create_parser = subparsers.add_parser("create")
    add_identity_arguments(create_parser)
    create_parser.add_argument("--output", type=Path, required=True)
    create_parser.set_defaults(handler=write_manifest)

    field_parser = subparsers.add_parser("field")
    add_identity_arguments(field_parser)
    field_parser.add_argument(
        "--field",
        choices=("artifact_id", "release_tag", "asset_name"),
        required=True,
    )
    field_parser.set_defaults(handler=print_field)

    validate_parser = subparsers.add_parser("validate")
    add_identity_arguments(validate_parser)
    validate_parser.add_argument("--manifest", type=Path, required=True)
    validate_parser.set_defaults(handler=validate)

    args = parser.parse_args()
    args.repo_root = args.repo_root.resolve()
    args.handler(args)


if __name__ == "__main__":
    main()
