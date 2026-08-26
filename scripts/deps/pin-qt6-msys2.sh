#!/usr/bin/env bash
set -euo pipefail

QT_PIN_VERSION="6.11.0"
BASE_URL="https://repo.msys2.org/mingw/ucrt64"

# Alle installierten qt6-* Pakete ermitteln, die NICHT schon auf der Pin-Version sind
pacman -Q | awk '{print $1, $2}' | grep '^mingw-w64-ucrt-x86_64-qt6-' | while read -r pkg version; do
  if [[ "$version" != "${QT_PIN_VERSION}"* ]]; then
    echo "Pinne $pkg von $version auf $QT_PIN_VERSION..."
    # Build-Nummer variiert pro Paket, daher zuerst die verfügbare Datei im Repo suchen
    FILE=$(curl -s "${BASE_URL}/" | grep -oE "${pkg}-${QT_PIN_VERSION}-[0-9]+-any\.pkg\.tar\.zst" | head -1)
    if [ -z "$FILE" ]; then
      echo "WARNUNG: keine ${QT_PIN_VERSION}-Version von $pkg im Repo gefunden" >&2
      continue
    fi
    pacman -U --noconfirm "${BASE_URL}/${FILE}"
  fi
done
