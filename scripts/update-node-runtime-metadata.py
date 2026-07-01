#!/usr/bin/env python3

import argparse
import json
import pathlib
import urllib.request


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Refresh pinned bundled Node.js runtime URLs and SHA256 checksums."
    )
    parser.add_argument(
        "--metadata",
        default="packaging/node-runtime.json",
        help="Path to node runtime metadata JSON.",
    )
    parser.add_argument(
        "--version",
        help="Node.js version to pin, without leading 'v'. Defaults to the JSON nodeVersion.",
    )
    args = parser.parse_args()

    metadata_path = pathlib.Path(args.metadata)
    data = json.loads(metadata_path.read_text())
    version = args.version or data["nodeVersion"]
    version = version.removeprefix("v")

    shasums_url = f"https://nodejs.org/dist/v{version}/SHASUMS256.txt"
    with urllib.request.urlopen(shasums_url, timeout=30) as response:
        shasums_text = response.read().decode("utf-8")

    shasums = {}
    for line in shasums_text.splitlines():
        checksum, filename = line.split(maxsplit=1)
        shasums[filename] = checksum

    data["nodeVersion"] = version
    for platform, entry in data["platforms"].items():
        if platform == "win-x64":
            filename = "win-x64/node.exe"
            url = f"https://nodejs.org/dist/v{version}/win-x64/node.exe"
        else:
            extension = "".join(pathlib.PurePosixPath(entry["url"]).suffixes[-2:])
            filename = f"node-v{version}-{platform}{extension}"
            url = f"https://nodejs.org/dist/v{version}/{filename}"

        if filename not in shasums:
            raise SystemExit(f"missing checksum for {filename} in {shasums_url}")
        entry["url"] = url
        entry["sha256"] = shasums[filename]

    metadata_path.write_text(json.dumps(data, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
