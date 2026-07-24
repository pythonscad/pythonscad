#!/usr/bin/env python3
"""Entfernt freie ('bool operator==/!=/</<=/>/>=(...)') Vergleichsoperator-
Deklarationen aus allen PyQt6-.sip-Dateien. Diese Deklarationen gehen davon
aus, dass Qt dafuer einen passenden freien Operator bereitstellt -- neuere
Qt6-Patch-Versionen brechen die ADL-Aufloesung fuer mehrere davon (Kollision
mit QCborTag/QCborKnownTags-Operatoren aus unabhaengigen Headern). Das
Entfernen kostet nur die Python-seitige Vergleichbarkeit (==, !=, <, <=, >,
>=) fuer den jeweiligen Qt-Typ -- alles andere bleibt unberuehrt.
"""
import re
import sys
from pathlib import Path

# Nur unindentierte (Spalte 0), freie Operator-Deklarationen, die direkt mit
# ");" enden -- Member-Operatoren ("    bool operator==(...) const;")
# haben fuehrende Leerzeichen und zusaetzlichen Text vor dem ";" und werden
# durch den Anker "^bool" bzw. das fehlende " const" vor ";" nicht getroffen.
PATTERN = re.compile(r'^bool operator(==|!=|<=|>=|<|>)\([^)]*\);\s*$')


def patch_file(path: Path) -> int:
    lines = path.read_text().splitlines(keepends=True)
    kept = []
    removed = 0
    for line in lines:
        if PATTERN.match(line.strip()):
            removed += 1
            continue
        kept.append(line)
    if removed:
        path.write_text("".join(kept))
    return removed


def main():
    if len(sys.argv) != 2:
        print("Usage: patch_pyqt6_free_operators.py <sip-root-dir>", file=sys.stderr)
        sys.exit(1)

    root = Path(sys.argv[1])
    total = 0
    for sip_file in root.rglob("*.sip"):
        n = patch_file(sip_file)
        if n:
            print(f"{sip_file}: {n} Operator-Deklaration(en) entfernt")
            total += n
    print(f"Insgesamt entfernt: {total}")


if __name__ == "__main__":
    main()
