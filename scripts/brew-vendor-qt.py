#!/usr/bin/env python3
"""Vendor a consistent set of Qt formulae from a fixed homebrew-core commit
into a local Homebrew tap, so that macOS CI stays on one Qt minor version.

Why this exists
---------------
Homebrew splits Qt 6 into many formulae (`qtbase`, `qtsvg`, `qtmultimedia`,
...). The `qt` formula is only a meta package whose dependencies are plain,
unqualified names such as `depends_on "qtbase"`. Copying just `qt.rb` from an
old commit into a local tap therefore does *not* pin Qt: Homebrew resolves
those dependencies against homebrew/core and pours whatever version core ships
today (e.g. 6.12 when we want 6.11).

This script instead downloads every requested formula *and* its Qt runtime
dependencies (transitively) from the same homebrew-core commit, rewrites the
dependency names to the local tap, and points the bottle blocks back at the
homebrew/core bottle storage so the old bottles are poured instead of building
Qt from source.

Usage
-----
    brew-vendor-qt.py --commit <sha> qtbase qtsvg qscintilla2 ...

Afterwards install the formulae with their tap-qualified names, e.g.
`brew install local/pinned-qt/qtbase`.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Set

RAW_URL = "https://raw.githubusercontent.com/Homebrew/homebrew-core/{commit}/Formula/{shard}/{name}.rb"
CORE_BOTTLE_ROOT = "https://ghcr.io/v2/homebrew/core"

# Formulae that belong to the pinned Qt stack. Any runtime dependency whose
# name matches is vendored from the same commit; everything else (icu4c,
# freetype, ...) is left to homebrew/core.
QT_STACK_RE = re.compile(r"^(qt[a-z0-9]*|qscintilla2)$")

DEPENDS_RE = re.compile(r'^\s*depends_on\s+"([^"]+)"(.*)$')


def shard_for(name: str) -> str:
    return name[0].lower()


def fetch(commit: str, name: str) -> str:
    url = RAW_URL.format(commit=commit, shard=shard_for(name), name=name)
    try:
        return subprocess.check_output(["curl", "-fsSL", "--retry", "3", url], text=True)
    except subprocess.CalledProcessError as exc:
        raise SystemExit(f"error: could not download {url} (curl exit {exc.returncode})")


def runtime_qt_deps(source: str) -> List[str]:
    """Qt-stack dependencies that are needed at runtime (no `=> :build/:test`)."""
    deps: List[str] = []
    for line in source.splitlines():
        m = DEPENDS_RE.match(line)
        if not m:
            continue
        name, rest = m.group(1), m.group(2)
        if "=>" in rest:
            continue
        if QT_STACK_RE.match(name) and name not in deps:
            deps.append(name)
    return deps


def rewrite(source: str, tap: str, vendored: Set[str]) -> str:
    def qualify(name: str) -> str:
        return f"{tap}/{name}" if name in vendored else name

    # depends_on "qtbase"  ->  depends_on "local/pinned-qt/qtbase"
    source = re.sub(
        r'(depends_on\s+")([^"/]+)(")',
        lambda m: m.group(1) + qualify(m.group(2)) + m.group(3),
        source,
    )
    # Formula["qtbase"]  ->  Formula["local/pinned-qt/qtbase"]
    # (unqualified lookups become ambiguous once the name exists in two taps)
    source = re.sub(
        r'(Formula\[")([^"/]+)("\])',
        lambda m: m.group(1) + qualify(m.group(2)) + m.group(3),
        source,
    )
    # Bottles of a non-core tap are looked up under that tap's own registry
    # by default; point them back at homebrew/core where they were uploaded.
    if "root_url" not in source:
        source = re.sub(
            r"^(\s*)bottle do\n",
            lambda m: f'{m.group(0)}{m.group(1)}  root_url "{CORE_BOTTLE_ROOT}"\n',
            source,
            count=1,
            flags=re.MULTILINE,
        )
    return source


def tap_formula_dir(tap: str) -> Path:
    repo = Path(subprocess.check_output(["brew", "--repository", tap], text=True).strip())
    if not repo.is_dir():
        subprocess.check_call(["brew", "tap-new", "--no-git", tap])
    formula_dir = repo / "Formula"
    formula_dir.mkdir(parents=True, exist_ok=True)
    return formula_dir


def main(argv: List[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--commit", required=True, help="homebrew-core commit to take the formulae from")
    ap.add_argument("--tap", default="local/pinned-qt", help="local tap to write into (default: %(default)s)")
    ap.add_argument(
        "--formula-dir",
        type=Path,
        help="write formulae here instead of the tap's Formula directory (testing; skips `brew`)",
    )
    ap.add_argument("formulae", nargs="+", help="root formulae, e.g. qtbase qtsvg qscintilla2")
    args = ap.parse_args(argv)

    # Resolve the Qt-stack closure first, so rewriting knows every vendored name.
    sources: Dict[str, str] = {}
    queue = list(args.formulae)
    while queue:
        name = queue.pop(0)
        if name in sources:
            continue
        sources[name] = fetch(args.commit, name)
        queue.extend(d for d in runtime_qt_deps(sources[name]) if d not in sources)

    out_dir = args.formula_dir or tap_formula_dir(args.tap)
    out_dir.mkdir(parents=True, exist_ok=True)
    vendored = set(sources)
    for name, src in sorted(sources.items()):
        (out_dir / f"{name}.rb").write_text(rewrite(src, args.tap, vendored))
        ver = re.search(r'^\s*version\s+"([^"]+)"', src, re.MULTILINE) or re.search(
            r"(\d+\.\d+\.\d+)", src
        )
        print(f"vendored {args.tap}/{name} {ver.group(1) if ver else '?'}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
