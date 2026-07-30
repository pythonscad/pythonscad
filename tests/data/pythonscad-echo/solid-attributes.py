from openscad import *

# specific attributes
a=cube(2).rotz(30)
f=a.faces()[0]
print(f.matrix)
for c in f:
    print(c.paths)
    print(c.points)

# Generic Attributes
c=cube(10)
c.configuration=[1,2]
print(c.hasattr("myname"))
c.setattr("Hello","World")
print(c.hasattr("Hello"))
print(c.getattr("configuration"))
print(c.Hello)

transformed=c.multmatrix([
    [1,0,0,1],
    [0,1,0,2],
    [0,0,1,3],
    [0,0,0,1],
])
print(transformed.configuration)
print(transformed.Hello)
