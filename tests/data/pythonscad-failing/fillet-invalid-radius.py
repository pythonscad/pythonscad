"""CLI rendering must fail cleanly when the requested fillet cannot fit."""

from openscad import *

outer = cube([40, 40, 20])
opening = cube([34, 20, 30]) + [3, 10, -5]
(outer - opening).fillet(2, fn=8).show()
