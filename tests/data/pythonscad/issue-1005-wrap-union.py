"""Regression for issue #1005: wrap(r=) then union must keep letter wedges.

Wrapping extruded text around a radius and unioning with a sphere used to
drop wedges from letters such as S in the render even though preview looked
fine.
"""
from openscad import *

t = text("S", size=8).linear_extrude(height=1.2).rotx(90).wrap(r=20)
s = sphere(r=21, fn=24)
union(t, s).rotz(270).show()
