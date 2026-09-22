from pythonscad import *

# Anzahl Kreissegmente fuer glatte Rundungen
nfn = 64

r0 = circle(40, fn=nfn)              # Boden, z=0
r1 = circle(50, fn=nfn).up(100)       # oben
r2 = circle(40, fn=nfn).up(100)       # oben
r3 = circle(35, fn=nfn).up(5)

tap1 = circle(10, fn=32).roty(-90).translate([42.5, 0, 30])
tap2 = circle(10, fn=32).roty(90).translate([42.5, 0, 70])

wall1 = patch(r1, holes=[r0,tap1, tap2], grid_spacing_uv=3)
wall2 = patch(r2, holes=[r1], grid_spacing_uv=3)
wall3 = patch(r3, holes=[r2], grid_spacing_uv=3)

handle = patch(tap2, holes=[tap1], grid_spacing_uv=3, use_tangents=True)

# geschlossene Scheibe genau bei z=0.
boden1 = patch(r0)
boden2 = patch(r3)

topf = concat(wall1,  wall2, wall3,boden1, boden2, handle )
topf.show()
