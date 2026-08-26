"""Geometry-time fillet failures become ValueError in eager Python APIs."""

from pathlib import Path
from tempfile import TemporaryDirectory

from openscad import *


def expect_value_error(label, operation):
    try:
        operation()
    except ValueError as error:
        print(f"{label}: ValueError: {error}")
    except Exception as error:
        print(f"{label}: UNEXPECTED {type(error).__name__}: {error}")
    else:
        print(f"{label}: NO EXCEPTION")


outer = cube([40, 40, 20])
opening = cube([34, 20, 30]) + [3, 10, -5]
frame = outer - opening
invalid_fillet = frame.fillet(2, fn=8)
print("lazy construction: ok")
expect_value_error("face collision mesh", invalid_fillet.mesh)

with TemporaryDirectory() as directory:
    output = Path(directory) / "invalid.stl"
    expect_value_error("face collision export", lambda: invalid_fillet.export(str(output)))

expect_value_error("edge radius", lambda: cube(10).fillet(6, fn=8).mesh())

horizontal = cube([20, 10, 10])
vertical = cube([10, 20, 10])
expect_value_error(
    "union radius",
    lambda: union(horizontal, vertical, r=30, fn=8).mesh(),
)
