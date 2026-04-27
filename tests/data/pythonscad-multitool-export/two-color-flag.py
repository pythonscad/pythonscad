"""MultiToolExporter regression fixture: two overlapping bricks split into
red and blue parts.

The test driver runs PythonSCAD with this script's directory set to a
scratch run directory, so the exporter's relative filenames ("red.stl",
"blue.stl") land where the driver expects them. The driver then byte-
compares each produced file against a checked-in golden after running
post_process_progname for header normalization.

Geometry is intentionally small and axis-aligned so the resulting ASCII
STL stays compact and stable under --enable=predictible-output.
"""
from pythonscad import MultiToolExporter, cube

background = cube([20, 20, 4])
star = cube([8, 8, 4]).translate([6, 6, 0])

MultiToolExporter("", ".stl", items=[
    (background, "blue"),
    (star, "red"),
]).export()
