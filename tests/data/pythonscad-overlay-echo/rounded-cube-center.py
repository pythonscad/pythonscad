"""Check rounded_cube(center=True) bounding-box placement."""
from pythonscad import *


def _assert_close(actual, expected, label):
    assert abs(actual - expected) < 1e-6, f"{label}: got {actual}, expected {expected}"


default_box = rounded_cube([30, 20, 10], r=2, fn=100)
none_center_box = rounded_cube([30, 20, 10], r=2, center=None, fn=100)
centered_box = rounded_cube([30, 20, 10], r=2, center=True, fn=100)

default_position = list(default_box.position)
default_size = list(default_box.size)
none_center_position = list(none_center_box.position)
none_center_size = list(none_center_box.size)
centered_position = list(centered_box.position)
centered_size = list(centered_box.size)

for index, value in enumerate(default_position):
    assert value >= -1e-6, f"default position {index} must stay in positive octant"

for index, (actual, expected) in enumerate(zip(none_center_position, default_position)):
    _assert_close(actual, expected, f"center=None position {index}")

for index, (actual, expected) in enumerate(zip(none_center_size, default_size)):
    _assert_close(actual, expected, f"center=None size {index}")

for index, (actual, expected) in enumerate(zip(centered_size, default_size)):
    _assert_close(actual, expected, f"centered size {index}")

for index, (actual, size) in enumerate(zip(centered_position, centered_size)):
    _assert_close(actual, size * -0.5, f"centered position {index}")

print("OK")
