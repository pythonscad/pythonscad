"""wrap() around a 2D outline without an experimental feature flag."""
from openscad import *

flat = text("S", size=5).linear_extrude(height=1).rotx(90)
flat.wrap(circle(r=10, fn=24)).show()
