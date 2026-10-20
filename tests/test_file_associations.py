#!/usr/bin/env python3
"""Validate desktop/MIME/appdata/plist/MSIX file-association metadata.

Source-tree check (no build required). Ensures:
  - .scad is advertised as application/x-openscad (default-eligible)
  - .py is Open With only via text/x-python (never a custom *.py MIME claim)
"""

from __future__ import annotations

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]


def _read(rel: str) -> str:
    path = ROOT / rel
    if not path.is_file():
        raise AssertionError(f"missing file: {rel}")
    return path.read_text(encoding="utf-8")


def check_desktop() -> None:
    text = _read("resources/icons/pythonscad.desktop.in")
    match = re.search(r"^MimeType=(.+)$", text, re.MULTILINE)
    assert match, "desktop.in missing MimeType="
    mimes = {m.strip() for m in match.group(1).split(";") if m.strip()}
    assert "application/x-openscad" in mimes, mimes
    assert "text/x-python" in mimes, mimes


def check_mime_xml() -> None:
    text = _read("resources/icons/pythonscad.xml")
    assert 'type="application/x-openscad"' in text
    assert 'pattern="*.scad"' in text
    assert 'pattern="*.py"' not in text
    assert "application/x-pythonscad" not in text


def check_mimeapps() -> None:
    text = _read("resources/icons/pythonscad-mimeapps.list.in")
    assert "[Default Applications]" in text
    assert "application/x-openscad=org.pythonscad.pythonscad${SUFFIX_WITH_DASH}.desktop;" in text
    assert "text/x-python" not in text
    cmake = _read("CMakeLists.txt")
    # Must merge into shared mimeapps.list, never install(FILES ... RENAME mimeapps.list).
    assert "merge-mimeapps-default.py" in cmake
    assert "RENAME mimeapps.list" not in cmake


def check_merge_script() -> None:
    import tempfile

    script = ROOT / "scripts" / "merge-mimeapps-default.py"
    assert script.is_file()
    import importlib.util

    spec = importlib.util.spec_from_file_location("merge_mimeapps_default", script)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    with tempfile.TemporaryDirectory() as tmp:
        path = pathlib.Path(tmp) / "mimeapps.list"
        path.write_text(
            "[Default Applications]\n"
            "text/x-python=code.desktop;\n"
            "application/x-openscad=openscad.desktop;\n",
            encoding="utf-8",
        )
        mod.merge(path, "org.pythonscad.pythonscad.desktop", "application/x-openscad")
        text = path.read_text(encoding="utf-8")
        assert "text/x-python=code.desktop;" in text
        assert "application/x-openscad=org.pythonscad.pythonscad.desktop;" in text
        assert "openscad.desktop" not in text

        # Multi-handler value: remove only our exact desktop ID.
        path.write_text(
            "[Default Applications]\n"
            "application/x-openscad=openscad.desktop;org.pythonscad.pythonscad.desktop;\n"
            "text/x-python=code.desktop;\n",
            encoding="utf-8",
        )
        mod.remove(path, "org.pythonscad.pythonscad.desktop", "application/x-openscad")
        text = path.read_text(encoding="utf-8")
        assert "application/x-openscad=openscad.desktop;" in text
        assert "org.pythonscad.pythonscad.desktop" not in text
        assert "text/x-python=code.desktop;" in text
        # Must not match a different desktop ID that merely contains ours as a substring.
        path.write_text(
            "[Default Applications]\n"
            "application/x-openscad=custom-org.pythonscad.pythonscad.desktop;\n",
            encoding="utf-8",
        )
        mod.remove(path, "org.pythonscad.pythonscad.desktop", "application/x-openscad")
        text = path.read_text(encoding="utf-8")
        assert "application/x-openscad=custom-org.pythonscad.pythonscad.desktop;" in text


def check_appdata() -> None:
    text = _read("pythonscad.appdata.xml.in")
    assert "<mediatype>application/x-openscad</mediatype>" in text
    assert "<mediatype>text/x-python</mediatype>" in text
    assert "application/x-pythonscad" not in text


def check_info_plist() -> None:
    text = _read("Info.plist.in")
    assert "<string>scad</string>" in text
    assert "<string>py</string>" in text
    # .py must be Alternate (Open With), never Apple default-for-type.
    py_block = text.split("<string>py</string>", 1)[1].split("</dict>", 1)[0]
    assert "<string>Alternate</string>" in py_block
    assert "LSIsAppleDefaultForType" not in py_block
    scad_block = text.split("<string>scad</string>", 1)[1].split("</dict>", 1)[0]
    assert "LSIsAppleDefaultForType" in scad_block


def check_msix_action() -> None:
    text = _read(".github/actions/build-msix/action.yml")
    assert "windows.fileTypeAssociation" in text
    assert 'FileType = ".scad"' in text or '".scad"' in text
    assert 'FileType = ".py"' in text or '".py"' in text
    # Both associations must be constructed in the manifest builder.
    assert "PythonSCAD Design" in text
    assert "PythonSCAD Python File" in text
    # UAP FileTypeAssociation.Name must be 3–100 characters (not "py").
    assert 'Name = "python"' in text
    assert 'Name = "py"' not in text


def check_destdir_merge_guard() -> None:
    cmake = _read("CMakeLists.txt")
    assert "PYTHONSCAD_MERGE_MIMEAPPS" in cmake
    appimage = _read("scripts/create_appimage.sh")
    assert "PYTHONSCAD_MERGE_MIMEAPPS=1" in appimage


def main() -> int:
    checks = (
        check_desktop,
        check_mime_xml,
        check_mimeapps,
        check_merge_script,
        check_appdata,
        check_info_plist,
        check_msix_action,
        check_destdir_merge_guard,
    )
    failures = 0
    for check in checks:
        try:
            check()
            print(f"OK  {check.__name__}")
        except AssertionError as exc:
            failures += 1
            print(f"FAIL {check.__name__}: {exc}", file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
