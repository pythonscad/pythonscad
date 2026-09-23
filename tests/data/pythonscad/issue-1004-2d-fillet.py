"""Regression for issue #1004: union/difference(..., r=) must fillet 2D seams."""
from pythonscad import *

s1 = union(square(10, center=True), square(10, center=True).rotz(45), r=1.4, fn=8)
s2 = difference(square(10, center=True), square(10, center=True).rotz(45).back(10), r=1.4, fn=8).right(13)
s3 = intersection(square(10, center=True), square(10, center=True).rotz(45), r=1.4, fn=8).right(26)
show([s1, s2, s3])
