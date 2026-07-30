"""Echo test for ``MultiToolExporter`` hidden cutters (``export=False``).

Exercises:
  * mixed 2- and 3-tuples preserving stored shape
  * strict ``bool`` validation for the export flag
  * ``parts()`` / ``show()`` omitting hidden entries
  * cumulative subtraction including hidden cutters
  * hidden duplicate names allowed; exportable duplicates still rejected
  * all-hidden ``export(single_file=...)`` no-op
"""

import pythonscad
from pythonscad import MultiToolExporter, cube

red = cube(10)
blue = cube(10).right(5)
green = cube(10).right(10)
cutter = cube(4).right(2)

exp = MultiToolExporter("p-", ".stl", items=[
    ("red", red),
    ("blue", blue),
    ("slot", cutter, False),
])
print("stored shapes:", [len(item) for item in exp])
print("parts len:", len(exp.parts()))
print("parts names:", [name for name, _g in exp.parts()])
print("part0 uses hidden cutter:", exp._part(0) is not exp[0][1])
print("part1 uses hidden cutter:", exp._part(1) is not exp[1][1])
print("part2 is bare object:", exp._part(2) is exp[2][1])

def expect(label, fn, exc):
    try:
        fn()
    except exc as e:
        print(f"{label}: {type(e).__name__}")
    else:
        print(f"{label}: NO EXCEPTION (expected {exc.__name__})")

expect("append export int", lambda: exp.append(("x", red, 1)), TypeError)
expect("append export zero", lambda: exp.append(("x", red, 0)), TypeError)
expect("append export str", lambda: exp.append(("x", red, "y")), TypeError)
try:
    exp.append(("hidden", green, False))
    print("append export false ok: appended")
    exp.pop()
except Exception as e:
    print(f"append export false ok: {type(e).__name__}")

hidden_dup = MultiToolExporter("p-", ".stl", items=[
    ("a", red),
    ("a", cutter, False),
    ("b", blue),
    ("a", green, False),
])
hidden_dup._check_unique_filenames()
hidden_dup._check_unique_part_names()
print("hidden duplicate names export: ok")
duplicate = MultiToolExporter("p-", ".stl", items=[("a", red), ("a", blue)])
duplicate._part = lambda _i: (_ for _ in ()).throw(RuntimeError("geometry built"))
expect(
    "exportable duplicate names",
    lambda: duplicate.export(single_file="duplicate.3mf"),
    ValueError,
)

all_hidden = MultiToolExporter("p-", ".stl", items=[
    ("only", red, False),
    ("also", blue, False),
])
print("all-hidden parts len:", len(all_hidden.parts()))
all_hidden.export(single_file="empty.3mf")
print("all-hidden single-file no-op: ok")

show_calls = []
real_show = pythonscad.show
pythonscad.show = lambda obj: show_calls.append(obj is not None)
try:
    exp.show()
finally:
    pythonscad.show = real_show
print("show call count:", len(show_calls))

export_calls = []
real_export = pythonscad.export
def recording_export(obj, filename):
    if isinstance(obj, dict):
        export_calls.append((list(obj.keys()), filename))
    else:
        export_calls.append((obj is not None, filename))
pythonscad.export = recording_export
try:
    exp.export()
    exp.export(single_file="assembly.3mf")
finally:
    pythonscad.export = real_export

for obj_info, filename in export_calls:
    if filename == "assembly.3mf":
        print("single-file keys:", obj_info)
    elif filename.startswith("p-"):
        print("per-file export:", obj_info, filename)
