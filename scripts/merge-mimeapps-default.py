#!/usr/bin/env python3
"""Merge or remove a single Default Applications entry in mimeapps.list.

Preserves unrelated associations. Used at cmake install time and from RPM
%post / %postun so PythonSCAD never overwrites a shared mimeapps.list.

Usage:
  merge-mimeapps-default.py merge  <mimeapps.list> <desktop-id> <mime-type>
  merge-mimeapps-default.py remove <mimeapps.list> <desktop-id> <mime-type>
"""

from __future__ import annotations

import pathlib
import sys


def _parse(text: str) -> list[tuple[str | None, list[str]]]:
    """Return [(section_name|None for preamble, lines_without_header), ...]."""
    sections: list[tuple[str | None, list[str]]] = [(None, [])]
    for raw in text.splitlines():
        line = raw.rstrip("\n")
        if line.startswith("[") and line.endswith("]") and len(line) > 2:
            sections.append((line[1:-1], []))
        else:
            sections[-1][1].append(line)
    return sections


def _serialize(sections: list[tuple[str | None, list[str]]]) -> str:
    out: list[str] = []
    for name, lines in sections:
        if name is not None:
            if out and out[-1] != "":
                out.append("")
            out.append(f"[{name}]")
        out.extend(lines)
    text = "\n".join(out)
    if text and not text.endswith("\n"):
        text += "\n"
    return text


def _default_apps(sections: list[tuple[str | None, list[str]]]) -> list[str]:
    for name, lines in sections:
        if name == "Default Applications":
            return lines
    sections.append(("Default Applications", []))
    return sections[-1][1]


def merge(path: pathlib.Path, desktop_id: str, mime_type: str) -> None:
    text = path.read_text(encoding="utf-8") if path.is_file() else ""
    sections = _parse(text)
    lines = _default_apps(sections)
    entry = f"{mime_type}={desktop_id};"
    key = f"{mime_type}="
    replaced = False
    for i, line in enumerate(lines):
        if line.startswith(key) and not line.lstrip().startswith("#"):
            lines[i] = entry
            replaced = True
            break
    if not replaced:
        lines.append(entry)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_serialize(sections), encoding="utf-8")


def remove(path: pathlib.Path, desktop_id: str, mime_type: str) -> None:
    if not path.is_file():
        return
    sections = _parse(path.read_text(encoding="utf-8"))
    for name, lines in sections:
        if name != "Default Applications":
            continue
        keep: list[str] = []
        for line in lines:
            if line.startswith(f"{mime_type}=") and desktop_id in line.split("=", 1)[1]:
                continue
            keep.append(line)
        lines[:] = keep
    path.write_text(_serialize(sections), encoding="utf-8")


def main(argv: list[str]) -> int:
    if len(argv) != 5 or argv[1] not in ("merge", "remove"):
        print(__doc__.strip(), file=sys.stderr)
        return 2
    action, path_s, desktop_id, mime_type = argv[1], argv[2], argv[3], argv[4]
    path = pathlib.Path(path_s)
    if action == "merge":
        merge(path, desktop_id, mime_type)
    else:
        remove(path, desktop_id, mime_type)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
