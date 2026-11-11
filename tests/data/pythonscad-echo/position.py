from openscad import *
fn=32

a = cube(10)
b = cube(10).translate([5, 5, 5])
c = cylinder(h=20, d=6)

print(a.position)         # vector(0,0,0)
print(b.position)         # vector(5,5,5)
print((a | b).position)   # vector(0,0,0)
print((a - b).position)   # vector(0,0,0)
print((a & b).position)   # vector(5,5,5)
print(c.position)         # vector(-3,-3,0)

# 2D objects
s = square([10, 20])
t = square([10, 20]).translate([5, 5])
u = circle(10)

print(s.position)         # vector(0,0)
print(t.position)         # vector(5,5)
print((s | t).position)   # vector(0,0)
print(u.position)         # vector(-10,-10)
