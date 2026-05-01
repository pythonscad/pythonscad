"""Echo test for issue #587: dict-form ``export()`` must raise
``TypeError`` (not SIGSEGV) when any key is not a ``str``.

Each call below is wrapped in :func:`expect` -- the same pattern used by
``multitool-exporter.py`` -- which prints the exception class so the
fixture's stdout becomes a stable golden:

  - The ``out`` path is in a ``tempfile.TemporaryDirectory`` so a partial
    or stray ``.3mf`` cannot litter the build tree (and the directory
    listing is printed at the end as a sanity check that nothing was
    written -- a correctly raised ``TypeError`` on every entry must
    fail before opening the output stream).
  - The fixture covers the four user-reachable crash classes from
    `pyfunctions.cc`:
      1. ``python_export_core`` dict branch  (outer-dict non-``str`` key)
      2. ``Export3mfPartInfo::writeProps``    (props_3mf non-``str`` key)
      3. ``python_show_core`` selection-handle dict (str-only by
         construction today, defensive coverage only)
      4. ``python_export_obj_att`` attribute dict (str-only by
         construction today, defensive coverage only)
"""

import os
import tempfile

from openscad import cube, export, show


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
    """Run ``fn`` and assert that *no* exception is raised.

    Prints ``OK`` on success and ``UNEXPECTED <ExcName>: <msg>`` on failure,
    so the golden output for a "should succeed" case is unambiguous instead
    of leaning on ``expect(..., Exception)`` and reading "NO EXCEPTION
    (expected Exception)" as the success line.
    """
    try:
        fn()
    except Exception as e:  # pragma: no cover - guards regressions
        print(f"{label}: UNEXPECTED {type(e).__name__}: {e}")
    else:
        print(f"{label}: OK")


with tempfile.TemporaryDirectory() as tmp:
    out_3mf = os.path.join(tmp, "x.3mf")
    out_stl = os.path.join(tmp, "x.stl")

    # --- 1. Outer dict: non-str keys (the issue #587 reproducer) -----
    expect("export int key", lambda: export({1: cube(10)}, out_3mf), TypeError)
    expect("export tuple key", lambda: export({(1, 2): cube(10)}, out_3mf), TypeError)
    expect("export None key", lambda: export({None: cube(10)}, out_3mf), TypeError)
    expect("export bytes key", lambda: export({b"x": cube(10)}, out_3mf), TypeError)
    # frozenset is hashable so it's a legal dict key but still not a str
    expect(
        "export frozenset key",
        lambda: export({frozenset([1, 2]): cube(10)}, out_3mf),
        TypeError,
    )
    # Same path is reachable for non-3MF formats too -- the dict branch
    # rejects the key well before the format-specific exporter runs.
    expect("export non-str key STL", lambda: export({1: cube(10)}, out_stl), TypeError)

    # --- 2. props_3mf dict with non-str keys (writeProps path) -------
    c_int_key = cube(10)
    c_int_key["props_3mf"] = {1: 1.5}
    expect(
        "props_3mf int key",
        lambda: export({"part": c_int_key}, out_3mf),
        TypeError,
    )

    c_tuple_key = cube(10)
    c_tuple_key["props_3mf"] = {(1, 2): "v"}
    expect(
        "props_3mf tuple key",
        lambda: export({"part": c_tuple_key}, out_3mf),
        TypeError,
    )

    # --- 3. props_3mf is 3MF-only: exporting the SAME bad-props_3mf
    # object to a non-3MF format must *not* raise, since props_3mf is
    # never consumed by the STL/OBJ/etc. exporters. Locking this in
    # so future pre-encode reshuffles don't accidentally regress
    # non-3MF exports on user-supplied props metadata.
    out_stl_part = os.path.join(tmp, "stl_part.stl")
    expect_no_exception(
        "props_3mf int key STL ignored",
        lambda: export({"part": c_int_key}, out_stl_part),
    )

    # --- 3b. Lone surrogate in a str key now raises TypeError instead
    # of silently substituting U+FFFD (which collides keys and produces
    # corrupted part names in the 3MF). The helper uses the "strict"
    # encoder so this fails cleanly at validation time.
    bad_str = "\udcff"  # lone low surrogate, valid str but not encodable
    expect(
        "export surrogate key",
        lambda: export({bad_str: cube(10)}, out_3mf),
        TypeError,
    )

    # --- 4. Sanity: directory listing. Only the legitimate STL we
    # wrote in step 3 should be present; every TypeError above must
    # have aborted before opening the output stream.
    leftover = sorted(os.listdir(tmp))
    print("leftover:", leftover)

# --- 4. Defensive coverage: show() on an object with a normal str
# attribute dict still works (regression check for the show_core sweep).
c_ok = cube(10)
c_ok["meta"] = "ok"
show(c_ok)
print("show ok")
