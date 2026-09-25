"""Echo test for ``MultiToolExporter`` (pythonscad-only helper).

Exercises:
  * no-arg constructor and ``items=`` seeding without prefix/suffix
  * DeprecationWarning when constructor is given prefix/suffix (incl. ``""``)
  * shape and name validation in ``append`` / ``extend`` / ``insert`` /
    ``__setitem__`` / ``__iadd__`` (``+=``)
  * ``_part(last)`` returning the underlying object as-is (no degenerate
    one-child ``difference`` node)
  * duplicate-name detection at ``export()`` time
  * ``mkdir=True`` with a directory-less filename (must not raise)
  * end-to-end ``export(prefix=..., suffix=...)`` and legacy constructor
    path layout
  * ``export(single_file=...)`` calling the underlying ``pythonscad.export``
    once with the expected dict of named parts

The underlying ``pythonscad.export`` is monkey-patched to a deterministic
recorder so the test does not depend on a render backend, file system
state, or output formatting of any export format.
"""

import os
import tempfile
import warnings

import pythonscad
from pythonscad import MultiToolExporter, cube

# --- 1. Constructor and basic list semantics ----------------------------
red = cube(10)
blue = cube(10).right(5)
green = cube(10).right(10)

# Preferred: no path layout at construction.
exp = MultiToolExporter()
exp.append(("red", red))
exp.append(("blue", blue))
print("len:", len(exp))
# Instance defaults are empty; path layout comes from export().
print("filename[0] default:", exp._filename(0))
print("filename[0] override:", exp._filename(0, prefix="p-", suffix=".stl"))
print("filename[1] override:", exp._filename(1, prefix="p-", suffix=".stl"))

# Explicit "" / "" at construction must warn (None sentinel distinguishes
# MultiToolExporter() from MultiToolExporter("", "")).
with warnings.catch_warnings(record=True) as w_empty:
    warnings.simplefilter("always", DeprecationWarning)
    MultiToolExporter("", "")
print("ctor empty strings warns:", any(
    issubclass(w.category, DeprecationWarning) for w in w_empty
))
with warnings.catch_warnings(record=True) as w_none:
    warnings.simplefilter("always", DeprecationWarning)
    MultiToolExporter()
print("ctor no-arg warns:", any(
    issubclass(w.category, DeprecationWarning) for w in w_none
))
with warnings.catch_warnings(record=True) as w_legacy:
    warnings.simplefilter("always", DeprecationWarning)
    legacy = MultiToolExporter("p-", ".stl", items=[("red", red), ("blue", blue)])
print("ctor legacy warns:", any(
    issubclass(w.category, DeprecationWarning) for w in w_legacy
))
print("legacy filename[0]:", legacy._filename(0))
print("legacy filename[1]:", legacy._filename(1))

# --- 2. _part(last) returns the bare object -----------------------------
print("part(last) is bare object:", exp._part(len(exp) - 1) is exp[-1][1])

# --- 2b. parts() shape: ordered list of (name, geometry) pairs ----------
ps = exp.parts()
print("parts len:", len(ps))
print("parts names:", [name for name, _g in ps])
print("parts last is bare:", ps[-1][1] is exp[-1][1])

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
expect("append wrong arity", lambda: exp.append(("x", red, "y")), TypeError)
expect("append non-str name", lambda: exp.append((42, red)), TypeError)
expect("append empty name", lambda: exp.append(("", red)), ValueError)
expect("extend bad item", lambda: exp.extend([("ok", red), red]), TypeError)
expect("insert bad item", lambda: exp.insert(0, red), TypeError)
expect("setitem bad item", lambda: exp.__setitem__(0, red), TypeError)
expect("setitem slice bad", lambda: exp.__setitem__(slice(0, 1), [red]), TypeError)
# `+=` goes through list.__iadd__ at the C level, which bypasses the
# extend() override unless __iadd__ is also overridden -- exercise it.
expect("iadd bad item", lambda: exp.__iadd__([red]), TypeError)
expect("iadd mid-bad item", lambda: exp.__iadd__([("ok", red), red]), TypeError)
expect(
    "constructor bad item",
    lambda: MultiToolExporter(items=[("ok", red), ("", blue)]),
    ValueError,
)

# Make sure no half-applied state was left behind by the failed extend
# or the failed +=.
print("len after failed extend:", len(exp))

# --- 4. Duplicate-name detection at export time ------------------------
dup = MultiToolExporter(items=[("x", red), ("x", blue)])
expect(
    "export duplicate names",
    lambda: dup.export(prefix="p-", suffix=".stl"),
    ValueError,
)
expect(
    "single-file duplicate names",
    lambda: dup.export(single_file="assembly.3mf"),
    ValueError,
)
expect(
    "single-file non-3mf",
    lambda: exp.export(single_file="assembly.stl"),
    ValueError,
)
empty = MultiToolExporter()
empty.export(single_file="empty.3mf")
print("single-file empty no-op: ok")

# --- 5. End-to-end export(), with monkey-patched underlying export -----
calls = []
real_export = pythonscad.export
def recording_export(obj, filename):
    if isinstance(obj, dict):
        calls.append((list(obj.keys()), filename))
    else:
        calls.append((obj is not None, filename))
pythonscad.export = recording_export
try:
    # 5a. Path layout on export(); mkdir with no directory -- must not crash
    e1 = MultiToolExporter(items=[("red", red), ("blue", blue), ("green", green)])
    e1.export(prefix="nodir-", suffix=".stl", mkdir=True)
    # 5b. With a real directory portion, mkdir=True -- still works
    with tempfile.TemporaryDirectory() as tmp:
        out_prefix = os.path.join(tmp, "nested", "x-")
        e2 = MultiToolExporter(items=[("a", red), ("b", blue)])
        e2.export(prefix=out_prefix, suffix=".3mf", mkdir=True)
        print("created nested dir:", os.path.isdir(os.path.join(tmp, "nested")))
    # 5c. Single-file 3MF export -- no prefix/suffix needed
    with tempfile.TemporaryDirectory() as tmp:
        out_file = os.path.join(tmp, "assembly", "parts.3mf")
        e3 = MultiToolExporter(items=[("r", red), ("b", blue)])
        e3.export(single_file=out_file, mkdir=True)
        print("created assembly dir:", os.path.isdir(os.path.join(tmp, "assembly")))
    # 5d. Legacy constructor path layout still works
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        e4 = MultiToolExporter(
            "legacy-", ".stl", mkdir=True, items=[("a", red), ("b", blue)]
        )
        e4.export()
finally:
    pythonscad.export = real_export

# Print the recorded export calls (last filename only -- temp paths vary).
for obj_info, filename in calls:
    if filename.startswith("nodir-") or filename.startswith("legacy-"):
        print("export call:", obj_info, filename)
    else:
        print("export call:", obj_info, "tmp/.../" + os.path.basename(filename))
