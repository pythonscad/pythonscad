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

# 2D size/position must work as scale/extrude args (sequence + omitted-Z).
# scale() defaults missing Z to 1 — must not flatten the cube.
print(cube(1).scale(square([2, 3]).size).size)
sq = square([10, 20])
print(sq.linear_extrude(5, origin=sq.position).size)

# Native vectors are sequences for polygon points and solid-minus-vector translate
print(polygon([[0, 0], sq.size, [0, 20]]).size)
print((sq - sq.position).position)

# Scalar extrude args must TypeError (not crash)
try:
    square(1).linear_extrude(5, scale=2)
    print("scalar-scale-ok")
except TypeError:
    print("scalar-scale-typeerror")

# Iterator protocol: nested iter(it) must not crash, and repeated list()
# must not leak refs to the vector (refcount returns to baseline).
v = cube(1).position
it = iter(v)
print(list(it) == list(v) and list(iter(it)) == [])
import sys
base = sys.getrefcount(v)
for _ in range(50):
    list(v)
print(sys.getrefcount(v) == base)
