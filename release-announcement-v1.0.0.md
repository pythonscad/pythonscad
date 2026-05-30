# PythonSCAD v1.0.0 released

Hi everyone,

We are thrilled to announce **PythonSCAD v1.0.0** — our biggest release ever. This version marks a milestone: PythonSCAD has grown from a promising OpenSCAD fork with Python support into a mature, first-class tool for script-based 3D modeling. We have landed **83 pull requests** since v0.20.0, covering everything from fundamental API architecture to a real IPython shell, native vector math, a completely refreshed Windows experience, full rebranding, and dozens of bug fixes.

## Clean Python module layout: `openscad` / `pythonscad` / `_openscad` (new, breaking)

**v1.0.0 restructures the Python surface into three distinct modules** ([PR #579](https://github.com/pythonscad/pythonscad/pull/579)):

- **`_openscad`** — the C extension (renamed from `openscad`). Rarely imported directly; it is the engine room.
- **`openscad`** — a pure-Python overlay that re-exports `_openscad` 1:1. The home for OpenSCAD-compatible additions and future per-symbol deprecation notices.
- **`pythonscad`** — a pure-Python overlay that re-exports `openscad` 1:1. The home for PythonSCAD-only features.

The drop-in guarantee: `openscad.cube is _openscad.cube` and `pythonscad.cube is openscad.cube`, so **switching between `from openscad import *` and `from pythonscad import *` requires no other code change**. All three modules are shipped in the pip wheel, each with a PEP 561 type-stub package (`-stubs/`) so IDEs can provide autocompletion for whichever import style you prefer.

New scripts and templates use `from pythonscad import *`. Existing designs using `from openscad import *` continue to work unmodified.

*Note: the C extension rename from `openscad` to `_openscad` is a breaking change for code that imported `openscad` directly as a C extension (not via the overlay). Any code using `from openscad import *` is unaffected.*

## Real IPython shell (`--ipython`) and basic REPL (`--repl`) (new)

`pythonscad --ipython` now launches a **genuine IPython interactive shell** — the same rich `In [N]:` prompt, tab completion, history, and magic commands you get from the `ipython` command itself ([PR #600](https://github.com/pythonscad/pythonscad/pull/600)). IPython is bundled inside the AppImage, macOS .app, and Windows installer so it is available out of the box without any extra `pip install`.

```
$ pythonscad --ipython
In [1]: from pythonscad import *
In [2]: cube(10).color("Tomato")
```

For a lighter, dependency-free prompt, the new `--repl` flag opens the basic embedded Python REPL. Distro packages (`.deb`, `.rpm`) declare `python3-ipython` as a recommended dependency rather than bundling it, keeping the package footprint small.

## Native vector operations (new)

Solid geometry objects now support **arithmetic operators directly** ([PR #588](https://github.com/pythonscad/pythonscad/pull/588)): addition (`+`), subtraction (`-`), scaling (`*`), dot product, and cross product. This makes geometric expressions in Python feel natural without needing to call helper functions for every operation.

## MultiToolExporter (new)

A new `MultiToolExporter` helper in the `pythonscad` module ([PR #585](https://github.com/pythonscad/pythonscad/pull/585)) makes it straightforward to split a single model into multiple per-tool or per-color output files — a common need for multi-material or laser-cut workflows:

```python
from pythonscad import *

background = cube([200, 100, 1]).color("blue")
star       = cylinder(r=20, h=2, fn=5).translate([100, 50, -0.5]).color("red")

exporter = MultiToolExporter("out/flag-", ".stl", mkdir=True)
exporter.append(("blue", background))  # blue part, minus the star area
exporter.append(("red",  star))        # red part (later entries win)
exporter.export()
# writes out/flag-blue.stl and out/flag-red.stl

# Or: one 3MF file with two named parts
export(dict(exporter.parts()), "flag.3mf")
```

Later entries "win" over earlier ones, so overlapping volume is assigned to exactly one output file.

## Inline Python trust bar (new)

The **blocking modal trust dialog** that appeared before you could even read a Python file has been replaced with a **non-intrusive inline bar** above the editor ([PR #667](https://github.com/pythonscad/pythonscad/pull/667)). You can now open, scroll, and review a file before deciding to trust it. The bar appears at the top of the editor (pushing the text down, never overlapping it), Preview and Render are disabled until trust is granted, and once you click "Trust Design" the customizer populates immediately.

External-editor workflows are smoother too: if you edit a trusted file in VS Code while PythonSCAD is open, it auto-trusts on reload without asking again.

## Dynamic oversampling (new)

A new dynamic oversampling mode ([PR #554](https://github.com/pythonscad/pythonscad/pull/554)) improves render quality for curved surfaces by adapting sample density to the local geometry.

## Polyline improvements (new)

Polylines now participate correctly in **`difference()` and `intersection()` operations** ([PR #572](https://github.com/pythonscad/pythonscad/pull/572)), and polyline preview has been improved ([PR #575](https://github.com/pythonscad/pythonscad/pull/575)).

## `solid.c` color property (new)

A new read-only **`.c` property** on solid objects ([PR #562](https://github.com/pythonscad/pythonscad/pull/562)) exposes the RGBA color of the root `color()` wrapper as a `(r, g, b, a)` tuple — useful for reading back colors from a colored solid in Python code.

## Full PythonSCAD rebrand (new)

With v1.0.0, **PythonSCAD is fully rebranded throughout** ([PRs #661](https://github.com/pythonscad/pythonscad/pull/661), [#669](https://github.com/pythonscad/pythonscad/pull/669), [#672](https://github.com/pythonscad/pythonscad/pull/672), [#676](https://github.com/pythonscad/pythonscad/pull/676), [#684](https://github.com/pythonscad/pythonscad/pull/684), [#686](https://github.com/pythonscad/pythonscad/pull/686)): About dialog, application metadata, desktop file IDs, localization strings, homepage URLs, and the macOS dock icon all show PythonSCAD consistently. No more leftover "OpenSCAD" references in unexpected places.

## Windows improvements

v1.0.0 brings a significantly better Windows experience:

- **Per-user installer with on-demand UAC** ([PR #696](https://github.com/pythonscad/pythonscad/pull/696)): The installer no longer requires administrator rights by default. A "Just for me / For all users" page lets you choose. UAC elevation only fires if you explicitly request an all-users install — a UAC shield icon on the Next button signals when elevation will be requested.
- **Code-signed Windows installer** ([PR #703](https://github.com/pythonscad/pythonscad/pull/703)): The Windows installer is now signed with a Certum code-signing certificate (`CN=Open Source Developer Michael Postmann`), so Windows SmartScreen will recognize it as coming from a known publisher.
- **MSIX package** ([PR #638](https://github.com/pythonscad/pythonscad/pull/638)): A Windows Store–style `.msix` package is now produced alongside the traditional NSIS installer and ZIP archive.

## New platform support

Pre-built packages are now available for **Fedora 44** and **Ubuntu 26.04 LTS (Noble)** ([PR #624](https://github.com/pythonscad/pythonscad/pull/624)).

## Title bar shows full path (new)

The main window title bar now displays the **absolute file path** of the open design ([PR #576](https://github.com/pythonscad/pythonscad/pull/576)), so it's always clear exactly which file you're editing when working across multiple directories.

## Reliability and bug fixes

v1.0.0 includes a large number of stability and correctness fixes across the board:

- **Resize** now behaves consistently across all cases ([#557](https://github.com/pythonscad/pythonscad/pull/557)).
- **Fillet** calculation corrected ([#625](https://github.com/pythonscad/pythonscad/pull/625)).
- **`solid + solid`** is now a union as expected; `minkowski` test fixed ([#658](https://github.com/pythonscad/pythonscad/pull/658)).
- **SVG import**: title handling and `osimport:stroke` parameter fixed ([#570](https://github.com/pythonscad/pythonscad/pull/570), [#577](https://github.com/pythonscad/pythonscad/pull/577)).
- **Session restore** now persists per-window geometry ([#627](https://github.com/pythonscad/pythonscad/pull/627)) and the macOS IPC socket name is shortened to fit system limits ([#640](https://github.com/pythonscad/pythonscad/pull/640)).
- **Console output** from compile and parse operations is now routed correctly to the console widget ([#643](https://github.com/pythonscad/pythonscad/pull/643)).
- **Python binding robustness**: multiple memory leaks, UB, and error-handling gaps resolved ([#603](https://github.com/pythonscad/pythonscad/pull/603), [#608](https://github.com/pythonscad/pythonscad/pull/608), [#610](https://github.com/pythonscad/pythonscad/pull/610), [#611](https://github.com/pythonscad/pythonscad/pull/611), [#619](https://github.com/pythonscad/pythonscad/pull/619)).
- **CSG tree**: missing `CONCAT` case in switch handled ([#653](https://github.com/pythonscad/pythonscad/pull/653)).
- **Startup crash and double file dialog** from the welcome screen fixed ([#695](https://github.com/pythonscad/pythonscad/pull/695)).
- **Fontconfig 2.18.0** app-font scoring regression worked around ([#691](https://github.com/pythonscad/pythonscad/pull/691)).
- **macOS**: deprecated `kUTTypePNG` API replaced, dock icon corrected, theme mismatch in preferences fixed ([#651](https://github.com/pythonscad/pythonscad/pull/651), [#675](https://github.com/pythonscad/pythonscad/pull/675)).
- **`--define` on the command line** now works during animation too ([#699](https://github.com/pythonscad/pythonscad/pull/699)).
- **Welcome dialog** New button restored when Python support is disabled ([#637](https://github.com/pythonscad/pythonscad/pull/637)).

## Upstream changes

We have merged all changes from OpenSCAD into our codebase to bring in the latest bug fixes and enhancements from the upstream project ([#559](https://github.com/pythonscad/pythonscad/pull/559), [#583](https://github.com/pythonscad/pythonscad/pull/583)).

## Downloads

Pre-built packages are available for Linux (AppImage, `.deb`, `.rpm`), macOS (DMG), and Windows (installer, ZIP, and MSIX):

**https://pythonscad.org/downloads/**

As always, we would love to hear from you — feedback and bug reports are very welcome!
