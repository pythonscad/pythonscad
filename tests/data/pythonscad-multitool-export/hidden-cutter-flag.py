"""MultiToolExporter regression fixture: hidden cutter between two colors.

A final hidden ``export=False`` cutter subtracts from both printable colors but
produces no output file, so only ``blue.stl`` and ``red.stl`` are written.
"""
from pythonscad import MultiToolExporter, cube

background = cube([20, 20, 4])
slot = cube([2, 20, 4]).translate([9, 0, 0])
star = cube([8, 8, 4]).translate([6, 6, 0])

MultiToolExporter("", ".stl", items=[
    ("blue", background),
    ("red", star),
    ("slot", slot, False),
]).export()
