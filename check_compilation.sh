#!/bin/bash

echo "=== Syntax-Prüfung der neuen Dateien ==="
echo ""

# Prüfe python_conversions.h
echo "1. Checking python_conversions.h..."
gcc -fsyntax-only -std=c++17 src/python/python_conversions.h 2>&1 | grep -E "error|warning" && echo "ERRORS FOUND" || echo "✓ No syntax errors"

# Prüfe python_conversions.cc (ohne Python.h)
echo ""
echo "2. Checking python_conversions.cc structure..."
if grep -q "#include.*Python.h" src/python/python_conversions.cc; then
  echo "✓ Python.h included"
else
  echo "✗ Missing Python.h include"
fi

if grep -q "#include.*matrix.h" src/python/python_conversions.cc; then
  echo "✓ matrix.h included"
else
  echo "✗ Missing matrix.h include"
fi

# Prüfe python_getitem_hier.h
echo ""
echo "3. Checking python_getitem_hier.h..."
gcc -fsyntax-only -std=c++17 src/python/python_getitem_hier.h 2>&1 | grep -E "error|warning" && echo "ERRORS FOUND" || echo "✓ No syntax errors"

# Prüfe python_getitem_hier.cc
echo ""
echo "4. Checking python_getitem_hier.cc structure..."
if grep -q "#include.*python_getitem_hier.h" src/python/python_getitem_hier.cc; then
  echo "✓ python_getitem_hier.h included"
else
  echo "✗ Missing python_getitem_hier.h include"
fi

if grep -q "#include.*python_conversions.h" src/python/python_getitem_hier.cc; then
  echo "✓ python_conversions.h included"
else
  echo "✗ Missing python_conversions.h include"
fi

# Kontrolliere auf Braces-Matching
echo ""
echo "5. Checking brace balance..."

for file in src/python/python_conversions.{h,cc} src/python/python_getitem_hier.{h,cc}; do
  if [ -f "$file" ]; then
    braces=$(grep -o '{' "$file" | wc -l)
    closing=$(grep -o '}' "$file" | wc -l)
    if [ "$braces" -eq "$closing" ]; then
      echo "✓ $file: Braces balanced ($braces)"
    else
      echo "✗ $file: Brace mismatch (open=$braces, close=$closing)"
    fi
  fi
done

# Kontrolliere auf Parentheses-Matching
echo ""
echo "6. Checking parenthesis balance..."

for file in src/python/python_conversions.{h,cc} src/python/python_getitem_hier.{h,cc}; do
  if [ -f "$file" ]; then
    open=$(grep -o '(' "$file" | wc -l)
    close=$(grep -o ')' "$file" | wc -l)
    if [ "$open" -eq "$close" ]; then
      echo "✓ $file: Parentheses balanced ($open)"
    else
      echo "✗ $file: Parenthesis mismatch (open=$open, close=$close)"
    fi
  fi
done

echo ""
echo "=== Detaillierte Fehlerprüfung ==="
echo ""

# Suche nach häufigen C++-Fehlern
echo "7. Looking for common C++ issues..."

# Check for missing semicolons before closing braces
if grep -E '}[[:space:]]*[}]' src/python/python_*.{h,cc} 2>/dev/null; then
  echo "⚠ Potential missing semicolon before closing brace"
else
  echo "✓ No obvious missing semicolons detected"
fi

# Check for unmatched quotes
echo ""
echo "8. Checking string literals..."
for file in src/python/python_*.{h,cc}; do
  if [ -f "$file" ]; then
    # Count double quotes
    quotes=$(grep -o '"' "$file" | wc -l)
    if [ $((quotes % 2)) -eq 0 ]; then
      echo "✓ $file: String quotes balanced ($quotes)"
    else
      echo "⚠ $file: Odd number of quotes ($quotes)"
    fi
  fi
done
