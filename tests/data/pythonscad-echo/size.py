from openscad import *
fn=32

a = cube(10)
b = cube(10).translate([5, 5, 5])
c = cylinder(h=20, d=6)

print(a.size)         # vector(10,10,10)
print(b.size)         # vector(10,10,10)
print((a | b).size)   # vector(15,15,15)
print((a - b).size)   # vector(10,10,10)
print((a & b).size)   # vector(5,5,5)
print(c.size)         # vector(6,6,20)

# 2D objects
s = square([10, 20])
t = square([10, 20]).translate([5, 5])
u = circle(10)

print(s.size)         # vector(10,20)
print(t.size)         # vector(10,20)
print((s | t).size)   # vector(15,25)
print(u.size)         # vector(20,20)
