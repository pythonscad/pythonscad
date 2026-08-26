#!/usr/bin/env bash
set -euo pipefail
echo "DEBUG: pin-qt6-msys2.sh gestartet" >&2

QT_PIN_VERSION="6.11.0"
BASE_URL="https://repo.msys2.org/mingw/ucrt64"

echo "DEBUG: Phase 1 - explizite Kernpakete installieren/pinnen" >&2
for pkg in qt6-base qt6-5compat qt6-multimedia qt6-svg ; do
  FULL_PKG="mingw-w64-ucrt-x86_64-${pkg}"
  FILE=$(curl -s "${BASE_URL}/" | grep -oE "${FULL_PKG}-${QT_PIN_VERSION}-[0-9]+-any\.pkg\.tar\.zst" | head -1) || true
  if [ -z "$FILE" ]; then
    echo "WARNUNG: keine ${QT_PIN_VERSION}-Version von $FULL_PKG im Repo gefunden" >&2
    continue
  fi
  echo "Installiere/pinne $FULL_PKG auf $QT_PIN_VERSION..." >&2
  pacman -U --noconfirm "${BASE_URL}/${FILE}" || {
    status=$?
    echo "WARNUNG: pacman -U für $FULL_PKG gab Exit-Code $status zurück (evtl. nur 'up to date -- reinstalling')" >&2
  }
done

echo "DEBUG: qscintilla-qt6 installieren (eigenes Versionsschema, kein Qt-Pin nötig)" >&2
pacman -S --needed --noconfirm mingw-w64-ucrt-x86_64-qscintilla-qt6

echo "DEBUG: Phase 2 - alle qt6-* Pakete nachpinnen (transitive Deps)" >&2
pacman -Q | grep '^mingw-w64-ucrt-x86_64-qt6-' >&2 || echo "DEBUG: keine qt6-* Pakete gefunden" >&2

pacman -Q | awk '{print $1, $2}' | grep '^mingw-w64-ucrt-x86_64-qt6-' | while read -r pkg version; do
  if [[ "$version" != "${QT_PIN_VERSION}"* ]]; then
    echo "Pinne $pkg von $version auf $QT_PIN_VERSION..." >&2
    FILE=$(curl -s "${BASE_URL}/" | grep -oE "${pkg}-${QT_PIN_VERSION}-[0-9]+-any\.pkg\.tar\.zst" | head -1) || true
    if [ -z "$FILE" ]; then
      echo "WARNUNG: keine ${QT_PIN_VERSION}-Version von $pkg im Repo gefunden" >&2
      continue
    fi
    pacman -U --noconfirm "${BASE_URL}/${FILE}" || {
    status=$?
    echo "WARNUNG: pacman -U für $FULL_PKG gab Exit-Code $status zurück (evtl. nur 'up to date -- reinstalling')" >&2
  }
  fi
done || true

echo "DEBUG: pin-qt6-msys2.sh fertig" >&2
