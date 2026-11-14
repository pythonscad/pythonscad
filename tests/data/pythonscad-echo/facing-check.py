from openscad import *

# Width / spacing checks between facing surfaces: internal() / external()

fn = 64


def rep(r):
    m = r["min"]
    return r["count"] > 0, None if m is None else round(m, 4)


def box(o):
    b = o.bbox
    return [round(v, 3) for v in b.position], [round(v, 3) for v in b.size]


# plate 0.5 thick
plate = cube([10, 10, 0.5])
print("plate internal 0.8", rep(plate.internal(0.8, report=True)))
print("plate internal 0.3", rep(plate.internal(0.3, report=True)))
print("plate external 0.8", rep(plate.external(0.8, report=True)))
print("plate error solid", box(plate.internal(0.8, grow=0)))
print("plate error solid grown", box(plate.internal(0.8, grow=0.05)))

# gaps: 0.3 inside a, 0.4 between a and b
a = cube(1) | cube(1).translate([1.3, 0, 0])
b = cube(1).translate([2.7, 0, 0])
print("gap inside a", rep(external(a, 0.5, report=True)))
print("gap a-b", rep(external(a, 0.5, other=b, report=True)))
print("gap a-b 0.35", rep(external(a, 0.35, other=b, report=True)))
print("gap a-b error solid", box(a.external(0.5, other=b, grow=0)))

# tube R5 r4.5: wall 0.5, hole ~9
tube = cylinder(r=5, h=10) - cylinder(r=4.5, h=12).translate([0, 0, -1])
print("tube wall", rep(tube.internal(0.8, report=True)))
print("tube hole projecting", rep(tube.external(20, alpha=0, report=True)))
print("tube hole oblique", rep(tube.external(20, report=True)))
print("tube hole angle 175", rep(tube.external(20, angle=175, report=True)))

# wedge: base 2, height 10 -> width 1.0 at y = 5
wedge = polygon([[-1, 0], [1, 0], [0, 10]]).linear_extrude(height=5)
print("wedge", rep(wedge.internal(1.0, report=True)))
print("wedge angle 170", rep(wedge.internal(1.0, angle=170, report=True)))
print("wedge error solid", box(wedge.internal(1.0, grow=0)))

# step: all walls >= 4.5, edge to edge through air 2.236
step = cube([5, 5, 10]) | cube([5, 5, 4.5]).translate([6, 0, 8])
print("step", rep(step.internal(4.4, report=True)))
print("step no occlusion", rep(step.internal(4.4, occlusion=False, report=True)))
print("step no occlusion projecting", rep(step.internal(4.4, occlusion=False, alpha=0, report=True)))

# offset plates: lateral 0.2, vertical 0.3 -> diagonal 0.3606
plates = cube([10, 10, 1]) | cube([9.8, 10, 1]).translate([10.2, 0, 1.3])
for alpha in (90, 45, 30, 0):
    print("offset plates alpha", alpha, rep(plates.external(0.5, alpha=alpha, report=True)))

# clean part
print("clean", rep(cube(5).internal(1, report=True)))
