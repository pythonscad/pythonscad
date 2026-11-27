#pragma once

#include <Python.h>

#include <string>
#include <vector>

/**
 * Decide which `sys.modules` entries must be dropped when the embedded
 * interpreter is re-used for the next script run (re-init path of
 * `initPython()`).
 *
 * Entries are judged per *top-level package*, never individually: a package
 * and all of its submodules are either all kept or all dropped. Dropping only
 * part of a package leaves `sys.modules` and the package attributes
 * inconsistent (e.g. `importlib` dropped while the frozen
 * `importlib.machinery` stays cached, so a fresh `importlib` never gets its
 * `machinery` attribute back -> `AttributeError` in `pkgutil`).
 *
 * A top-level package is kept when any of its modules is
 *   - part of the standard library (`sys.stdlib_module_names`),
 *   - frozen or built-in,
 *   - a compiled extension (`.pyd`, `.so`, `.dylib`) — these cannot be
 *     re-imported reliably (PyQt6/sip, numpy, ...),
 *   - installed in a `site-packages` directory or below `usr/lib`.
 * Everything else (the user's own scripts and libraries) is dropped so edits
 * are picked up on the next run.
 *
 * @param modules       the `sys.modules` dict (borrowed reference)
 * @param forceDrop     top-level names that are always dropped, together with
 *                      all their submodules (the PythonSCAD pure-Python
 *                      overlays `openscad` and `pythonscad`)
 * @return              the keys to delete from `sys.modules`
 *
 * Never leaves a Python exception set. Caller must hold the GIL.
 */
std::vector<std::string> python_collect_stale_modules(PyObject *modules,
                                                      const std::vector<std::string>& forceDrop);
