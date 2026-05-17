"""Echo test for import() error handling.

Verifies that osimport() raises proper Python exceptions immediately at the
call site rather than silently returning a valid object and logging a cryptic
message during geometry evaluation.

Scenarios covered:
- Empty filename       -> ValueError
- Nonexistent file     -> FileNotFoundError
- Path is a directory  -> OSError

Scenarios NOT covered here (unavoidable deferred errors):
- Corrupt/invalid file contents: these are detected during geometry evaluation
  after osimport() returns, so they still surface via LOG() rather than as
  Python exceptions.
"""

import os
import tempfile

from openscad import cube, osimport


def expect(label, fn, exc):
    try:
        fn()
    except exc as e:
        print(f"{label}: {type(e).__name__}")
    except Exception as e:
        print(f"{label}: UNEXPECTED {type(e).__name__}: {e}")
    else:
        print(f"{label}: NO EXCEPTION (expected {exc.__name__})")


# Empty filename
expect("empty filename", lambda: osimport(""), ValueError)

# Nonexistent file
expect(
    "nonexistent file",
    lambda: osimport("__nonexistent_file_that_cannot_exist_xyz__.stl"),
    FileNotFoundError,
)

# Path is a directory (use a temp dir to keep it portable)
with tempfile.TemporaryDirectory() as tmp:
    expect("directory path", lambda: osimport(tmp), OSError)
