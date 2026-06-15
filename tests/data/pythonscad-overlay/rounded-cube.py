"""Render test: rounded_cube() via the pythonscad overlay package."""
from pythonscad import *

rounded_cube(20, r=2).show()
rounded_cube([30, 20, 10], d=4).translate([22, 0, 0]).show()
