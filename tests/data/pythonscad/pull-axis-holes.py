"""Regression: pull() must not drop faces when edges are parallel to the pull axis.

A cube pulled apart through its center along X previously lost most straddling
faces (open holes) because the cut helper became singular for X-aligned edges.
Y and Z pulls on the same cube already worked; include them so all three axes
stay covered.
"""
from openscad import *

c = cube(10)
cx = c.pull([5, 5, 5], [4, 0, 0])
cy = c.pull([5, 5, 5], [0, 4, 0]).right(20)
cz = c.pull([5, 5, 5], [0, 0, 4]).right(40)
show(cx + cy + cz)
