"""Behavior checks for the PythonSCAD-only Customizer container."""

import openscad
import pythonscad
from pythonscad import Customizer


assert pythonscad.add_parameter is openscad.add_parameter

params = Customizer()
base = params.group("Base parameters")
assert base is params.group("Base parameters")

returned_diameter = base.add_parameter("diameter", 10, description="Diameter")
label = params.add_parameter("label", "plate")

assert returned_diameter == 10
assert label == "plate"
assert dict(params) == {"diameter": 10.0, "label": "plate"}
assert dict(base) == {"diameter": 10.0}
assert len(params) == 2
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

print("OK")
