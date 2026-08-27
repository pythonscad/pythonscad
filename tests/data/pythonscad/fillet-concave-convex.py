"""Regression: nearby convex and concave rails fit on one face at small r."""

from openscad import *

outer = cube([40, 40, 20])
opening = cube([34, 20, 30]) + [3, 10, -5]
frame = outer - opening

frame.fillet(1, fn=8).show()
