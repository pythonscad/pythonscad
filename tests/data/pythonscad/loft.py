from pythonscad import *


profile = loft(
    square(10, center=True),
    circle(r=3, fn=24),
    height=20,
    n=24,
    rot=45,
)
linear_extrude(profile, height=20, slices=20).show()
