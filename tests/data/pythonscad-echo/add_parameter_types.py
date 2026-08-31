"""Regression: int defaults must stay int when finished values are reapplied.

The GUI dry-run / re-preview path stores Customizer numbers as OpenSCAD
doubles and feeds them back through customizer_parameters_finished. This
test seeds that finished-parameter path via -D (see CMake ARGS) and
checks that Python int/float/str types are preserved on return.
"""

from openscad import add_parameter

int_value = add_parameter("int_param", 256)
float_value = add_parameter("float_param", 256.0)
str_value = add_parameter("str_param", "256")

print(f"int: {type(int_value).__name__} {int_value!r}")
print(f"float: {type(float_value).__name__} {float_value!r}")
print(f"str: {type(str_value).__name__} {str_value!r}")
