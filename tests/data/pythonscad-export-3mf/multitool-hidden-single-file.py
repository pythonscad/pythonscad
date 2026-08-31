"""Regression fixture for hidden cutters in ``export(single_file=...)``.

The final hidden entry subtracts from earlier parts but must not appear as a
named object in the output 3MF.
"""
from pythonscad import MultiToolExporter, cube

red_part = cube([6, 6, 4]).color("red")
slot = cube([4, 4, 4]).translate([1, 1, 0])
blue_part = cube([6, 6, 4]).translate([8, 0, 0]).color("blue")

MultiToolExporter("", ".stl", items=[
    ("red", red_part),
    ("blue", blue_part),
    ("slot", slot, False),
]).export(single_file="multitool-hidden.3mf")
