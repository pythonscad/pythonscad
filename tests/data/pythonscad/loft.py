from pythonscad import *


try:
    loft(square(1), square(2), height=0)
except ValueError as error:
    assert str(error) == "loft height must be non-zero"
else:
    raise AssertionError("loft() accepted a zero height")

profile = loft(
    square(10, center=True),
    circle(r=3, fn=24),
    height=20,
    n=24,
    rot=45,
)
linear_extrude(profile, height=20, slices=20).show()
