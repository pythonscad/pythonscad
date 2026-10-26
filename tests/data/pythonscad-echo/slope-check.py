import math

from pythonscad import *

# Face angle checks: slope() / overhang() / draft()

def rep(r):
    w = r["worst"]
    return (r["count"], r["angle"], r["undercut"], None if w is None else round(w, 2))

def box(o):
    b = o.bbox
    return [round(v, 2) for v in b.position], [round(v, 2) for v in b.size]

fn = 64
c = cube(10)
print("cube overhang", rep(c.overhang(report=True)))
print("cube draft", rep(c.draft(report=True)))

taper = cylinder(h=20, r1=10, r2=10 - 20 * math.tan(math.radians(5)))
print("taper5 draft2", rep(taper.draft(2, report=True)))
print("taper5 draft10", rep(taper.draft(10, report=True)))
print("taper5 upside down draft2", rep(taper.roty(180).draft(2, report=True)))

cone30 = cylinder(h=10, r1=1, r2=1 + 10 * math.tan(math.radians(30)))
cone60 = cylinder(h=10, r1=1, r2=1 + 10 * math.tan(math.radians(60)))
print("cone30 overhang45", rep(cone30.overhang(45, report=True)))
print("cone60 overhang45", rep(cone60.overhang(45, report=True)))

mush = cylinder(h=10, r=2) | cylinder(h=3, r=8).translate([0, 0, 10])
print("mushroom overhang", rep(mush.overhang(report=True)))
print("mushroom draft", rep(mush.draft(-0.5, report=True)))

C = cube([20, 20, 3]) | cube([20, 20, 3]).translate([0, 0, 10]) | cube([3, 20, 13])
print("C draft free", rep(C.draft(-0.5, report=True)))
print("C draft no undercut", rep(C.draft(-0.5, undercut=False, report=True)))
print("C draft parting 6.5", rep(C.draft(-0.5, parting=6.5, report=True)))
print("C draft pull x", rep(C.draft(-0.5, dir=[1, 0, 0], report=True)))
print("C trapped volume", box(C.draft(-0.5, grow=0)))
print("cube draft skin", box(c.draft(grow=0)))

print("slope window", rep(c.slope(dir=[0, 0, 1], min=-10, max=10, report=True)))
