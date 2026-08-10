"""Behavior checks for the PythonSCAD-only Customizer container."""

import openscad
import pythonscad
from pythonscad import Customizer


assert pythonscad.add_parameter is openscad.add_parameter

params = Customizer()
base = params.group("Base parameters")
assert base is params.group("Base parameters")

returned_diameter = base.add_parameter("diameter", 10, description="Diameter")
returned_scale = params.add_parameter("scale", 1.0, description="Scale")
label = params.add_parameter("label", "plate")

assert type(returned_diameter) is int
assert type(returned_scale) is float
assert returned_diameter == 10
assert returned_scale == 1.0
assert label == "plate"
assert dict(params) == {"diameter": 10, "scale": 1.0, "label": "plate"}
assert dict(base) == {"diameter": 10}
assert len(params) == 3
assert len(base) == 1

assert "diameter" not in globals()

try:
    params.add_parameter("diameter", 20)
except ValueError as error:
    assert "already registered" in str(error)
else:
    raise AssertionError("duplicate parameter name was accepted")

try:
    Customizer().add_parameter("diameter", 20)
except ValueError as error:
    assert "already registered" in str(error)
else:
    raise AssertionError("duplicate parameter name across containers was accepted")

try:
    Customizer().add_parameter("unsupported", object())
except TypeError as error:
    assert "unsupported default type" in str(error)
else:
    raise AssertionError("unsupported default type was accepted")

print("OK")
