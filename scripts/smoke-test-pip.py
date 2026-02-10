#!/usr/bin/env python3
"""Minimal smoke test for the pip-installed openscad module (no CMake/ctest)."""
from openscad import *

c = cube(5)
# show() requires the full GUI runtime (pythonMainModule), which is not
# initialised when openscad is loaded as a plain pip module.  Skip it here;
# the GUI path is already covered by the regular ctest suite.
export(c, "/tmp/pip-smoke.3mf")
print("smoke test OK")
