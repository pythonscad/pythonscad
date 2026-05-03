"""Echo test for issue #594: ``python__setitem__()`` must surface a
clean ``TypeError`` (mp_ass_subscript: ``return -1`` with the
helper's exception set) when the subscript key is not a ``str`` --
not silently swallow the error and report a successful assignment
that never actually wrote to the object's dict.

Pre-fix the bug shape was:

    int python__setitem__(...)
    {
      PyObject *keyname = PyUnicode_AsEncodedString(key, "utf-8", "~");
      if (keyname == nullptr) return 0;          // <-- swallows TypeError
      std::string keystr = PyBytes_AS_STRING(keyname);
      ...
      return 0;
    }

CPython's ``tp_as_mapping->mp_ass_subscript`` slot is documented to
return 0 only when no exception is set; returning 0-with-exception
is undefined behaviour and surfaces as a confusing ``SystemError:
<built-in function ...> returned a result with an exception set``
from whatever C-API call happens to run next.

PR #608 partially fixed the contract violation (changed ``return 0``
to ``return -1`` so CPython propagates the encoder's TypeError), but
the call site still used the bogus ``"~"`` error handler -- not a
registered Python error handler at all -- and the raw
``PyBytes_AS_STRING`` idiom that #587 / #595 swept everywhere else.
This change finishes the migration to ``python_pyobject_to_utf8``,
which produces the standard helper messages
("obj[key]: expected str, got <type>" for non-str keys, "obj[key]:
str cannot be UTF-8 encoded" for str keys the strict utf-8 handler
refuses) and consistent surrogate handling (strict, raises
TypeError) instead of the accidental ``LookupError`` from the
unregistered handler.

Each call below is wrapped in :func:`expect` -- the same helper used
by ``issue587-export-non-str-keys.py`` -- which prints the exception
class so the fixture's stdout becomes a stable golden.
"""

from openscad import cube


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


# --- 1. The issue #594 reproducer: c[non_str] = v ---------------
# Before the fix this silently appeared to succeed but left a stale
# TypeError on the error indicator (UB per the C-API contract).
# After the fix every non-str key path raises TypeError immediately.
c = cube(10)
expect("setitem int key", lambda: c.__setitem__(1, "good"), TypeError)
expect("setitem None key", lambda: c.__setitem__(None, "x"), TypeError)
expect("setitem bool key", lambda: c.__setitem__(True, "x"), TypeError)
expect("setitem tuple key", lambda: c.__setitem__((1, 2), "x"), TypeError)
expect("setitem bytes key", lambda: c.__setitem__(b"k", "x"), TypeError)
# frozenset is hashable so it's a legal dict key but still not a str.
expect(
    "setitem frozenset key",
    lambda: c.__setitem__(frozenset([1, 2]), "x"),
    TypeError,
)

# --- 2. Lone surrogates -- valid str but the strict utf-8 handler
# refuses to encode them, so the helper raises TypeError. Locks in
# that the new strict-handler path matches the rest of the codebase
# (#587 / #595) instead of the accidental LookupError the bogus
# "~" handler produced.
expect("setitem surrogate key", lambda: c.__setitem__("\udcff", "x"), TypeError)

# --- 3. del c[non_str] also routes through mp_ass_subscript with
# v == NULL. The same TypeError contract applies.
expect("delitem int key", lambda: c.__delitem__(1), TypeError)
expect("delitem None key", lambda: c.__delitem__(None), TypeError)

# --- 4. The successful str-key path is unchanged: a roundtrip set
# + read still works, and del removes the key. Locks in that the
# fix only tightens the failure contract.
c["meta"] = "ok"
print("get meta:", c["meta"])
expect_no_exception("delitem str key", lambda: c.__delitem__("meta"))
# After del, the get returns Py_None (PythonSCAD's documented
# missing-key contract -- not KeyError) so the success path stays
# observably distinct from the TypeError path above.
print("get meta after del:", c["meta"])

# --- 5. After every TypeError above, the object's dict must be
# unchanged (the assignment must not have partially landed) and no
# stale exception must leak. The simplest sanity check is that a
# follow-up valid assignment + read still works -- pre-fix this
# would have surfaced the stashed TypeError as a SystemError from
# the next C-API call, which is exactly the pathology #594 calls
# out.
c["after"] = "still-ok"
print("get after:", c["after"])

# --- 6. Attribute-style subscript (`c.attr = value`) routes through
# `tp_setattro` -> python__setattro__ -> python__setitem__, so the
# attribute path inherits the same fix. Python's syntax already
# constrains attr names to str at compile time, so the only way to
# hit the non-str path here is through the C API
# (`PyObject_SetAttr` with a non-str key) -- but route the str
# case through the public surface to make sure we didn't regress
# normal attribute assignment.
c.colour = "red"
print("get colour:", c.colour)

print("done")
