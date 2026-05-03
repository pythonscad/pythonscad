"""
Regression test for issue #596: PyOpenSCADObjectToNodeMulti dict leak.

Before the fix, every call that took the LIST branch of
`PyOpenSCADObjectToNodeMulti` allocated a fresh `PyDict_New()` whose
strong reference the caller never `Py_DECREF`-ed. That dict (plus the
borrowed values it pinned) accumulated one-per-call, so a tight loop
of e.g. `show([a, b])` leaked one dict per iteration -- visible in
`gc.get_objects()` as a steadily growing object population.

After the fix, the contract for `*dict` is "always a strong reference"
on every input shape, and every audited call site releases it
(typically via a `PyObjectUniquePtr` RAII wrapper). The
`gc.get_objects()` count therefore stays effectively flat across the
loop instead of growing linearly with N.

The script self-asserts: it prints `dict_leak_ok` on the fixed binary
and `dict_leak_FAIL` on the unfixed one, so the standard echo-test
comparator catches a regression deterministically by string match --
the *expected file* contains only that constant discriminator, no
numeric threshold. The script itself does pick a threshold (`N // 4`,
documented inline below) to decide which constant to print.
"""
import gc

from pythonscad import cube, cylinder, show

a = cube(10)
b = cylinder(r=5, h=10)

# Run a couple of warm-up iterations so any one-shot interpreter
# allocations (lazy module imports, jit-style caches, etc.) settle
# before we sample the baseline.
for _ in range(50):
    show([a, b])

gc.collect()
before = len(gc.get_objects())

N = 1000
for _ in range(N):
    show([a, b])

gc.collect()
after = len(gc.get_objects())
delta = after - before

# Pre-fix: delta scales ~linearly with N (dozens of objects per leaked
# dict because the merged dict pins the per-child borrowed values).
# Post-fix: delta should be a small constant. Use N/4 as the cutoff;
# that comfortably catches the one-dict-per-call leak while leaving
# headroom for unrelated, non-deterministic interpreter behaviour.
threshold = N // 4
if delta < threshold:
    print("dict_leak_ok")
else:
    print(f"dict_leak_FAIL delta={delta} threshold={threshold}")
