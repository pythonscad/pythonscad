"""Echo test: export() must not write files when preview is true.

CLI echo tests run with -o *.echo and without --render, so preview=True
(the same flag as GUI F5). export() is a no-op in that mode: it returns
without creating the target file.
"""

import os
import tempfile

from openscad import cube, export

with tempfile.TemporaryDirectory() as tmp:
    out = os.path.join(tmp, "should-not-exist.stl")
    result = export(cube(10), out)
    exists = os.path.exists(out)
    print(f"return: {result!r}")
    print(f"wrote file: {exists}")
