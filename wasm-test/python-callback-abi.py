from pythonscad import *

assert model() is None
assert modelpath().endswith("python-callback-abi.py")
x = add_parameter("x", 1, range=range(0, 101))
assert x == 1
cube(1 + x).show()
