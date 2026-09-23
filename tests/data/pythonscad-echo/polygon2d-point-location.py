"""Regression for Polygon2d point_location / point_inside boundary semantics."""
from openscad import *

s = square(10, center=True)
print("interior", s.inside([0, 0]))
print("edge", s.inside([5, 0]))
print("vertex", s.inside([5, 5]))
print("outside", s.inside([10, 10]))
