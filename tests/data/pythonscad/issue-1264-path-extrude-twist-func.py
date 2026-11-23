"""Regression for #1264: path_extrude(twist=callable) must not segfault.

In 1.2.0, twist/profile Python callbacks ran from the geometry worker without
the GIL (lockDepth is thread_local since #1024), so the first callback into the
Python C API crashed. This uses the reporter's minimal repro.
"""
from openscad import *


def tw(h):
  return 0.0


# Extrude a square along a path with rounded corners
p = path_extrude(
  square(1),
  [[0, 0, 0], [0, 0, 10, 3], [10, 0, 10, 3], [10, 10, 10]],
  twist=tw,
)
p.show()
