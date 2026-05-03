"""Echo test for issue #609: ``python_numberval()`` must report
"not coercible" when a ``str`` argument cannot be parsed as a double,
instead of silently returning success and leaving ``*result`` carrying
the caller's prior (possibly uninitialised) value while still flipping
the dragflags bit.

The bug is observable through any builtin that takes a vector of
numbers, forwards it to ``python_vectorval()``, and **passes a non-null
``flags`` argument** -- the str branch in ``python_numberval()`` is
gated on ``flags != nullptr`` so callers that omit ``flags`` already
reject any str entry in the vector with ``return 1`` regardless of
this fix. ``cube()`` is the canonical caller that *does* pass
``flags`` (`&node->dragflags`), so it's the cleanest user-visible
surface for the bug:

  Before the fix
    cube([10, "two", 30])
      -> python_numberval("two", &dim[1], &dragflags, 2)
         -> sscanf returns 0, *result unchanged, dragflags |= 2,
            return 0 (success)
      -> python_vectorval returns 0 (success)
      -> CubeNode is built with whatever dim[1] was default-init'd
         to (often 0 or stack garbage). The "must be positive"
         downstream guard catches the 0 case but loses the original
         "non-numeric str" provenance, and the dragflags bit is
         falsely set on a parameter the user never actually
         numerically supplied.

  After the fix
    cube([10, "two", 30])
      -> python_numberval("two", &dim[1], &dragflags, 2)
         -> sscanf returns 0, return 1 (not coercible),
            *result and *flags both untouched
      -> python_vectorval returns 1
      -> cube() raises TypeError("Invalid Cube dimensions")

Each call below is wrapped in :func:`expect` -- the same helper used by
``issue587-export-non-str-keys.py`` -- which prints the exception class
so the fixture's stdout becomes a stable golden.
"""

from openscad import color, cube, scale, sphere, translate


def expect(label, fn, exc):
    """Run ``fn`` and print whether it raised ``exc``."""
    try:
        fn()
    except exc as e:
        print(f"{label}: {type(e).__name__}")
    except Exception as e:  # pragma: no cover - guards regressions
        print(f"{label}: UNEXPECTED {type(e).__name__}: {e}")
    else:
        print(f"{label}: NO EXCEPTION (expected {exc.__name__})")


def expect_no_exception(label, fn):
    """Run ``fn`` and assert that *no* exception is raised."""
    try:
        fn()
    except Exception as e:  # pragma: no cover - guards regressions
        print(f"{label}: UNEXPECTED {type(e).__name__}: {e}")
    else:
        print(f"{label}: OK")


# --- 1. cube() with a non-numeric str inside the size vector --------
# This is the primary surface for the fix: cube() passes &dragflags
# to python_vectorval, which forwards it into python_numberval, so
# the str branch actually runs. Pre-fix these would have parsed as
# 0 / stack-garbage with dragflags falsely set; post-fix they raise.
expect("cube str in size", lambda: cube([10, "two", 30]), TypeError)
expect("cube empty str in size", lambda: cube([10, "", 30]), TypeError)
expect("cube whitespace str in size", lambda: cube([10, "   ", 30]), TypeError)

# --- 2. cube() with a leading number followed by garbage ------------
# sscanf("1.5abc", ...) returns 1 (the "1.5" parses cleanly) so this
# is *not* a "not coercible" case. Lock in that the partial-match
# path still succeeds, so the fix doesn't accidentally tighten the
# coercion contract beyond "literally nothing matched".
expect_no_exception("cube partial-match str", lambda: cube([10, "1.5abc", 30]))
expect_no_exception("cube scientific-notation str", lambda: cube([10, "2.5e1", 30]))

# --- 3. cube() with a numeric str -- the str branch's intended path.
# Before the fix this also worked (sscanf parsed it cleanly), the
# fix only tightens the *failure* path. Lock in that the success
# path is still permissive.
expect_no_exception("cube numeric str", lambda: cube([10, "20", 30]))

# --- 4. Defensive coverage of the flags=nullptr callers
# (color / scale / translate). These always reject any str entry in
# their vector argument, regardless of this fix, because
# python_color_core / python_scale_core / python_translate_core do
# not pass `flags` to python_vectorval and the str branch in
# python_numberval is gated on `flags != nullptr`. Locking in the
# "always TypeError" contract here so a future caller that starts
# passing `flags` gets a deliberate review of which str inputs it
# wants to accept.
c = cube(10)
expect("color str in vector", lambda: c.color([1.0, "green", 0.5]), TypeError)
expect("color numeric str in vector", lambda: c.color([1.0, "0.5", 0.25]), TypeError)
expect("translate str in vector", lambda: translate([1, "y", 3], cube(1)), TypeError)
expect("scale str in vector", lambda: scale([1, "y", 3], cube(1)), TypeError)

# Pure-str colour name still works (it's the PyUnicode_Check branch
# in python_color_core, not python_vectorval).
expect_no_exception("color named str", lambda: c.color("red"))

# --- 5. sphere(r=str) -- documents that scalar `r` strings now also
# surface as "not coercible". sphere() initialises its r/d locals to
# NAN, so before the fix the bug was masked into "silently use the
# default radius 1" (and dragflags |= 1 was falsely set). After the
# fix r stays NAN cleanly *and* dragflags is no longer flipped.
# The user-visible surface from a Python script is still "silent
# radius=1" because sphere() does not itself raise on missing r/d,
# but the internal state is now consistent.
expect_no_exception("sphere non-numeric r", lambda: sphere(r="not-a-number"))

print("done")
