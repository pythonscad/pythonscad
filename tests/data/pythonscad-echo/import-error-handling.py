"""Echo test for osimport() and osuse() error handling.

Verifies that osimport() and osuse() raise proper Python exceptions immediately
at the call site rather than silently returning a valid object or producing a
cryptic message during geometry evaluation.

Scenarios covered:
- osimport(): empty filename       -> ValueError
- osimport(): nonexistent file     -> FileNotFoundError
- osimport(): path is a directory  -> OSError
- osuse():    empty filename       -> ValueError
- osuse():    nonexistent file     -> FileNotFoundError
- osuse():    path is a directory  -> OSError

Scenarios NOT covered here (unavoidable deferred errors):
- Corrupt/invalid file contents: these are detected during geometry evaluation
  after osimport() returns, so they still surface via LOG() rather than as
  Python exceptions.
"""

import tempfile

from openscad import osimport, osuse


def expect(label, fn, exc):
    try:
        fn()
    except exc as e:
        print(f"{label}: {type(e).__name__}")
    except Exception as e:
        print(f"{label}: UNEXPECTED {type(e).__name__}: {e}")
    else:
        print(f"{label}: NO EXCEPTION (expected {exc.__name__})")


# --- osimport() ---
expect("osimport empty filename", lambda: osimport(""), ValueError)
expect(
    "osimport nonexistent file",
    lambda: osimport("__nonexistent_file_that_cannot_exist_xyz__.stl"),
    FileNotFoundError,
)
with tempfile.TemporaryDirectory() as tmp:
    expect("osimport directory path", lambda: osimport(tmp), OSError)

# --- osuse() ---
expect("osuse empty filename", lambda: osuse(""), ValueError)
expect(
    "osuse nonexistent file",
    lambda: osuse("__nonexistent_file_that_cannot_exist_xyz__.scad"),
    FileNotFoundError,
)
with tempfile.TemporaryDirectory() as tmp:
    expect("osuse directory path", lambda: osuse(tmp), OSError)
