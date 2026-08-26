#!/usr/bin/env bash
set -euo pipefail
echo "DEBUG: pin-qt6-msys2.sh gestartet" >&2
which pacman curl awk grep || echo "DEBUG: eines der Tools fehlt im PATH" >&2

QT_PIN_VERSION="6.11.0"
BASE_URL="https://repo.msys2.org/mingw/ucrt64"

echo "DEBUG: installierte qt6-* Pakete:" >&2
pacman -Q | grep '^mingw-w64-ucrt-x86_64-qt6-' >&2 || echo "DEBUG: keine qt6-* Pakete installiert!" >&2

pacman -Q | awk '{print $1, $2}' | grep '^mingw-w64-ucrt-x86_64-qt6-' | while read -r pkg version; do
  if [[ "$version" != "${QT_PIN_VERSION}"* ]]; then
    echo "Pinne $pkg von $version auf $QT_PIN_VERSION..." >&2
    FILE=$(curl -s "${BASE_URL}/" | grep -oE "${pkg}-${QT_PIN_VERSION}-[0-9]+-any\.pkg\.tar\.zst" | head -1)
    if [ -z "$FILE" ]; then
      echo "WARNUNG: keine ${QT_PIN_VERSION}-Version von $pkg im Repo gefunden" >&2
      continue
    fi
    pacman -U --noconfirm "${BASE_URL}/${FILE}"
  fi
done || true

echo "DEBUG: pin-qt6-msys2.sh fertig" >&2
