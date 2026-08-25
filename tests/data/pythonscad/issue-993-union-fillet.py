"""Regression for issue #993: union(..., r=) should fillet the full boolean seam.

A hex prism sitting on a cube used to skip seam edges whose endpoints were
original cylinder vertices that already lay on the cube face.
"""
from openscad import *

plate = cube([40, 40, 20])
tower = cylinder(50, 10, fn=6).rotx(30) + [20, 20, 0]
union(plate, tower, fn=32, r=2).show()
