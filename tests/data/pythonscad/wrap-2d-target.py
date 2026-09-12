from pythonscad import *

t = text("PythonScad", size=5).scale([1,1.4,1]).linear_extrude(height=1.2).rotx(90).wrap(square(10)).down(2.25).rotz(215)
show(t)
