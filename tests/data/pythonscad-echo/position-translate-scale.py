from openscad import *
fn=32

# 2D: scalar multiply on .position must negate components (not list-repeat)
# and recenter the circle at the origin.
c2 = circle(r=10)
c2 = c2.translate(c2.position * -1)
print(c2.position)

# Same pattern for 3D
c3 = cylinder(r=10, h=5)
c3 = c3.translate(c3.position * -1)
print(c3.position)

# 2D size scalar multiply
s = square([10, 20])
print(s.size * 2)
