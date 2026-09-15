
from pythonscad import *

t=text("e", size=15, halign="center", valign="center")
t.linear_extrude(3).fillet(0.5,fn=14).show()
