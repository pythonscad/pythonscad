#!/usr/bin/env python3
"""Install IPython 9.x into a --target directory on MSYS2 UCRT64.

IPython 9.x requires psutil>=7, but PyPI has no prebuilt psutil wheel
for the mingw_x86_64_ucrt_gnu-cpython ABI and source builds fail against
MSYS2's winternl.h. MSYS2 ships a pacman-built psutil; this script
installs IPython's other dependencies via pip, copies psutil from the
system site-packages into the bundle, then installs IPython itself with
--no-deps from the downloaded wheel.
"""

from __future__ import annotations

import argparse
import email
import pathlib
import shutil
import site
import subprocess
import sys
import tempfile
import zipfile

from packaging.markers import default_environment
from packaging.requirements import Requirement


def _run(cmd: list[str]) -> None:
    print(f"[bundle-ipython-msys2] running: {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, check=True)


def _download_ipython_wheel(
    python: str, requirements_file: pathlib.Path, dest_dir: pathlib.Path
) -> pathlib.Path:
    _run(
        [
            python,
            "-m",
            "pip",
            "download",
            "--no-deps",
            "-r",
            str(requirements_file),
            "-d",
            str(dest_dir),
        ]
    )
    wheels = sorted(dest_dir.glob("ipython-*.whl"), key=lambda p: p.name, reverse=True)
    if not wheels:
        raise SystemExit(f"no ipython wheel downloaded to {dest_dir}")
    return wheels[0]


def _read_requires_dist(wheel_path: pathlib.Path) -> list[Requirement]:
    with zipfile.ZipFile(wheel_path) as zf:
        meta_name = next(n for n in zf.namelist() if n.endswith(".dist-info/METADATA"))
        meta_text = zf.read(meta_name).decode("utf-8")
    msg = email.message_from_string(meta_text)
    return [Requirement(val) for key, val in msg.items() if key.lower() == "requires-dist"]


def _filter_requirements(
    requires_dist: list[Requirement], exclude_names: set[str]
) -> list[str]:
    env = default_environment()
    exclude = {name.lower() for name in exclude_names}
    selected: list[str] = []
    for req in requires_dist:
        if req.name.lower() in exclude:
            continue
        if req.marker is not None and not req.marker.evaluate(env):
            continue
        selected.append(str(req).split(";")[0])
    return selected


def _install_deps(python: str, target: pathlib.Path, req_strings: list[str]) -> None:
    if not req_strings:
        return
    _run(
        [
            python,
            "-m",
            "pip",
            "install",
            "--target",
            str(target),
            "--upgrade",
            "--no-compile",
            *req_strings,
        ]
    )


def _vend_psutil_from_system(target: pathlib.Path) -> None:
    import importlib.metadata

    import psutil

    src_pkg = pathlib.Path(psutil.__file__).parent
    dest_pkg = target / src_pkg.name
    if dest_pkg.exists():
        shutil.rmtree(dest_pkg)
    shutil.copytree(src_pkg, dest_pkg)

    for sp in [*site.getsitepackages(), site.getusersitepackages()]:
        if not sp:
            continue
        sp_path = pathlib.Path(sp)
        if not sp_path.is_dir():
            continue
        for meta in sp_path.glob("psutil-*.dist-info"):
            dest_meta = target / meta.name
            if dest_meta.exists():
                shutil.rmtree(dest_meta)
            shutil.copytree(meta, dest_meta)
            return

    version = importlib.metadata.version("psutil")
    raise SystemExit(f"psutil-{version} dist-info not found in site-packages")


def _install_ipython_wheel(python: str, target: pathlib.Path, wheel_path: pathlib.Path) -> None:
    _run(
        [
            python,
            "-m",
            "pip",
            "install",
            "--target",
            str(target),
            "--no-deps",
            "--no-compile",
            str(wheel_path),
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python", required=True, help="Python interpreter to drive pip")
    parser.add_argument("--target", required=True, type=pathlib.Path, help="Bundle staging dir")
    parser.add_argument(
        "--requirements", required=True, type=pathlib.Path, help="requirements/runtime.txt"
    )
    args = parser.parse_args()

    try:
        subprocess.run([args.python, "-c", "import psutil"], check=True, capture_output=True)
    except subprocess.CalledProcessError:
        raise SystemExit(
            "psutil is not importable in the build Python. "
            "On MSYS2 UCRT64, run: pacboy -S python-psutil:p"
        ) from None

    args.target.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="ipython-wheel.") as tmp:
        wheel_dir = pathlib.Path(tmp)
        wheel = _download_ipython_wheel(args.python, args.requirements, wheel_dir)
        deps = _filter_requirements(_read_requires_dist(wheel), exclude_names={"psutil"})
        _install_deps(args.python, args.target, deps)
        _vend_psutil_from_system(args.target)
        _install_ipython_wheel(args.python, args.target, wheel)


if __name__ == "__main__":
    main()
