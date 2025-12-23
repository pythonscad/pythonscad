# PythonSCAD - Python Conversions Refactoring

## Overview
Dies ist die restaurierte Version mit allen Python ↔ C++ Konvertierungsfunktionen konsolidiert in einem separaten Modul.

## Inhalt des Archivs

```
pythonscad-src/
├── src/python/
│   ├── python_conversions.h       (NEW) Header mit 11 Funktionen
│   └── python_conversions.cc      (NEW) Implementierung aller Funktionen
├── CMakeLists.txt                 (MODIFIED) python_conversions.cc hinzugefügt
├── setup.py                       (MODIFIED) python_conversions.cc hinzugefügt
├── PULL_REQUEST_TEMPLATE.md       (NEW) Detaillierte PR-Dokumentation
└── README.md                       (This file)
```

## Konsolidierte Funktionen (11 insgesamt)

### Python → C++ Konvertierungen
- `python_tomatrix()` - Python Liste zu Matrix4d
- `python_tovector()` - Python Liste zu Vector3d
- `python_numberval()` - Python Zahl zu double
- `python_vectorval()` - Python Liste zu 3D/4D Vektor

### C++ → Python Konvertierungen
- `python_frommatrix()` - Matrix4d zu Python Liste
- `python_fromvector()` - Vector3d zu Python Liste
- `python_frompointsflexible()` - Vector3d Array zu Python (2D/3D intelligente Konvertierung)
- `python_frompaths()` - std::vector<std::vector<size_t>> zu Python
- `python_fromfaces()` - std::vector<IndexedFace> zu Python
- `python_fromopenscad()` - OpenSCAD Value zu Python Objekt

### OpenSCAD Value Konvertierungen
- `python_convertresult()` - Python Objekt zu OpenSCAD Value

## Wichtige Features

### python_frompointsflexible()
- Intelligente 2D/3D Konvertierung
- Wenn z==0: gibt [x, y] zurück
- Wenn z!=0: gibt [x, y, z] zurück

## Wie du die Dateien verwendest

### 1. In dein PythonSCAD-Projekt kopieren
```bash
# Kopiere die neuen Dateien
cp src/python/python_conversions.* /path/to/your/pythonscad/src/python/

# Ersetze die modifizierten Dateien
cp CMakeLists.txt /path/to/your/pythonscad/
cp setup.py /path/to/your/pythonscad/
```

### 2. Änderungen in pyfunctions.cc anwenden
Entferne diese 7 Funktionen und ersetze mit Include:
```cpp
#include "python_conversions.h"
```

Funktionen zu entfernen:
- python_tomatrix
- python_tovector
- python_frommatrix
- python_fromvector
- python_frompointsflexible
- python_frompaths
- python_fromfaces

### 3. Änderungen in pyopenscad.cc anwenden
Entferne diese 4 Funktionen und ersetze mit Include:
```cpp
#include "python_conversions.h"
```

Funktionen zu entfernen:
- python_numberval
- python_vectorval
- python_fromopenscad
- python_convertresult

### 4. Kompilieren
```bash
cd /path/to/pythonscad
mkdir build
cd build
cmake ..
make -j$(nproc)
```

## Pull Request Erstellen

1. Feature Branch erstellen:
```bash
git checkout -b refactor/consolidate-python-conversions
```

2. Änderungen committen:
```bash
git add .
git commit -m "refactor: consolidate Python conversion functions"
```

3. PR erstellen mit der Dokumentation aus `PULL_REQUEST_TEMPLATE.md`

## Validierung

✅ Alle Dateien sind syntaktisch korrekt
✅ Braces und Parentheses sind balanced
✅ Alle Required Includes sind vorhanden
✅ Build System ist aktualisiert
✅ Memory Management ist korrekt

## Kontakt

Bei Fragen oder Problemen, kontaktiere das PythonSCAD-Team.
