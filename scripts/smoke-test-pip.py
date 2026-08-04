#!/usr/bin/env python3
"""Minimal smoke test for the pip-installed PythonSCAD module (no CMake/ctest).

Exercises the three-module layout shipped by the wheel:

* ``_openscad``  - C extension (low level).
* ``openscad``   - pure-Python overlay re-exporting ``_openscad`` (kept
  for compatibility with upstream OpenSCAD).
* ``pythonscad`` - strict superset of ``openscad``, the recommended
  import for PythonSCAD-only code.

Asserts the drop-in property: switching between ``from openscad import *``
and ``from pythonscad import *`` does not change which callables a script
binds (they are the same objects).
"""
import ast
import importlib.metadata
import json
import os
from pathlib import Path
import tempfile

import _openscad
import openscad
import pythonscad

distribution = importlib.metadata.distribution("pythonscad")
distribution_files = {
    str(path): path
    for path in (distribution.files or ())
}
required_typing_files = {
    "_openscad-stubs/__init__.pyi",
    "_openscad-stubs/py.typed",
    "openscad/py.typed",
    "pythonscad/py.typed",
}
missing_typing_files = required_typing_files - distribution_files.keys()
direct_url_text = distribution.read_text("direct_url.json")
direct_url = json.loads(direct_url_text) if direct_url_text else {}
is_editable = direct_url.get("dir_info", {}).get("editable", False)
if is_editable:
    project_root = Path(__file__).resolve().parents[1]
    editable_typing_files = {
        "_openscad-stubs/__init__.pyi": (
            project_root
            / "libraries/python/stubs/_openscad-stubs/__init__.pyi"
        ),
        "_openscad-stubs/py.typed": (
            project_root / "libraries/python/stubs/_openscad-stubs/py.typed"
        ),
        "openscad/py.typed": (
            project_root / "libraries/python/openscad/py.typed"
        ),
        "pythonscad/py.typed": (
            project_root / "libraries/python/pythonscad/py.typed"
        ),
    }
    missing_editable_typing_files = {
        name for name, path in editable_typing_files.items() if not path.is_file()
    }
    assert not missing_editable_typing_files, (
        "editable install is missing source type information: "
        f"{sorted(missing_editable_typing_files)}"
    )
    stub_path = editable_typing_files["_openscad-stubs/__init__.pyi"]
else:
    assert not missing_typing_files, (
        f"wheel is missing packaged type information: {sorted(missing_typing_files)}"
    )
    stub_path = distribution_files["_openscad-stubs/__init__.pyi"].locate()

stub_tree = ast.parse(stub_path.read_text(encoding="utf-8"), filename=str(stub_path))
stub_exports = {
    node.name
    for node in stub_tree.body
    if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef))
}
stub_exports.update(
    target.id
    for node in stub_tree.body
    if isinstance(node, ast.Assign)
    for target in node.targets
    if isinstance(target, ast.Name)
)
stub_exports.update(
    node.target.id
    for node in stub_tree.body
    if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name)
)
runtime_exports = {name for name in dir(_openscad) if not name.startswith("_")}
missing_stub_exports = runtime_exports - stub_exports
assert not missing_stub_exports, (
    "public _openscad exports missing from the bundled stub: "
    f"{sorted(missing_stub_exports)}"
)

assert openscad.cube is _openscad.cube, (
    "openscad.cube should be the same object as _openscad.cube"
)
assert pythonscad.cube is openscad.cube, (
    "pythonscad.cube should be the same object as openscad.cube"
)

assert set(dir(openscad)) >= set(n for n in dir(_openscad) if not n.startswith("_")), (
    "openscad must re-export every public name from _openscad"
)
assert set(dir(pythonscad)) >= set(n for n in dir(openscad) if not n.startswith("_")), (
    "pythonscad must re-export every public name from openscad"
)

# Hold an explicit reference to the TemporaryDirectory so it stays alive
# until the script exits; the directory (and therefore the export
# artifacts) is cleaned up automatically.  Using mkdtemp-style isolation
# keeps the smoke test safe to run in parallel and avoids leaving stale
# files behind in /tmp.
_workdir = tempfile.TemporaryDirectory(prefix="pythonscad-pip-smoke-")
WORKDIR = _workdir.name

from openscad import *  # noqa: F401,F403,E402

c = cube(5)
c.show()
export(c, os.path.join(WORKDIR, "pip-smoke-openscad.3mf"))

from pythonscad import *  # noqa: F401,F403,E402

c2 = cube(5)
c2.show()
export(c2, os.path.join(WORKDIR, "pip-smoke-pythonscad.3mf"))

print("smoke test OK")
