"""Echo test for ``MultiToolExporter`` (pythonscad-only helper).

Exercises:
  * constructor (positional + ``items=`` seeding) and ``len()``
  * shape and name validation in ``append`` / ``extend`` / ``insert`` /
    ``__setitem__``
  * ``_part(last)`` returning the underlying object as-is (no degenerate
    one-child ``difference`` node)
  * duplicate-name detection at ``export()`` time
  * ``mkdir=True`` with a directory-less filename (must not raise)
  * end-to-end ``export()`` calling the underlying ``pythonscad.export``
    once per part with the expected filename

The underlying ``pythonscad.export`` is monkey-patched to a deterministic
recorder so the test does not depend on a render backend, file system
state, or output formatting of any export format.
"""

import pythonscad
from pythonscad import MultiToolExporter, cube

# --- 1. Constructor and basic list semantics ----------------------------
red = cube(10)
blue = cube(10).right(5)
green = cube(10).right(10)

exp = MultiToolExporter("p-", ".stl")
exp.append((red, "red"))
exp.append((blue, "blue"))
print("len:", len(exp))
print("filename[0]:", exp._filename(0))
print("filename[1]:", exp._filename(1))

# --- 2. _part(last) returns the bare object -----------------------------
print("part(last) is bare object:", exp._part(len(exp) - 1) is exp[-1][0])

# --- 2b. parts() shape: ordered list of (name, geometry) pairs ----------
ps = exp.parts()
print("parts len:", len(ps))
print("parts names:", [name for name, _g in ps])
print("parts last is bare:", ps[-1][1] is exp[-1][0])

# --- 3. Validation: every entry path -----------------------------------
def expect(label, fn, exc):
    """Run ``fn`` and print whether it raised ``exc``."""
    try:
        fn()
    except exc as e:
        print(f"{label}: {type(e).__name__}")
    else:
        print(f"{label}: NO EXCEPTION (expected {exc.__name__})")

expect("append non-tuple", lambda: exp.append(red), TypeError)
expect("append wrong arity", lambda: exp.append((red, "x", "y")), TypeError)
expect("append non-str name", lambda: exp.append((red, 42)), TypeError)
expect("append empty name", lambda: exp.append((red, "")), ValueError)
expect("extend bad item", lambda: exp.extend([(red, "ok"), red]), TypeError)
expect("insert bad item", lambda: exp.insert(0, red), TypeError)
expect("setitem bad item", lambda: exp.__setitem__(0, red), TypeError)
expect("setitem slice bad", lambda: exp.__setitem__(slice(0, 1), [red]), TypeError)
expect(
    "constructor bad item",
    lambda: MultiToolExporter("p-", ".stl", items=[(red, "ok"), (blue, "")]),
    ValueError,
)

# Make sure no half-applied state was left behind by the failed extend.
print("len after failed extend:", len(exp))

# --- 4. Duplicate-name detection at export time ------------------------
dup = MultiToolExporter("p-", ".stl", items=[(red, "x"), (blue, "x")])
expect("export duplicate names", dup.export, ValueError)

# --- 5. End-to-end export(), with monkey-patched underlying export -----
calls = []
real_export = pythonscad.export
def recording_export(obj, filename):
    calls.append((obj is not None, filename))
pythonscad.export = recording_export
try:
    # 5a. Two parts, no directory in prefix, mkdir=True -- must not crash
    e1 = MultiToolExporter(
        "nodir-",
        ".stl",
        mkdir=True,
        items=[(red, "red"), (blue, "blue"), (green, "green")],
    )
    e1.export()
    # 5b. With a real directory portion, mkdir=True -- still works
    import tempfile, os
    with tempfile.TemporaryDirectory() as tmp:
        out_prefix = os.path.join(tmp, "nested", "x-")
        e2 = MultiToolExporter(
            out_prefix, ".3mf", mkdir=True, items=[(red, "a"), (blue, "b")]
        )
        e2.export()
        print("created nested dir:", os.path.isdir(os.path.join(tmp, "nested")))
finally:
    pythonscad.export = real_export

# Print the recorded export calls (last filename only -- temp paths vary).
for has_obj, filename in calls:
    if filename.startswith("nodir-"):
        print("export call:", has_obj, filename)
    else:
        print("export call:", has_obj, "tmp/.../" + os.path.basename(filename))
