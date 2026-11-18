# PythonSCAD v1.2.0 "Quadratum" released

Hi everyone,

We are happy to announce **PythonSCAD v1.2.0** "Quadratum".

Four months and more than 750 commits after v1.1.2, this is our biggest release
so far. The highlights:

- **[`patch()`](https://www.pythonscad.org/reference/primitives3d/#patch)** – a
  brand-new surface/solid construction primitive that builds geometry from a
  boundary plus holes (tube, panel and ring-bridge modes).
- **[PyQt6 inside PythonSCAD](https://www.pythonscad.org/development/pyqt6-and-forms/)**
  – bundled PyQt6 bindings let your Python classes open their own interactive
  dialogs from the editor, and you can extend the GUI with your own forms.
- **Interactive dragging, reworked** – drag `rotate`, `scale` and `translate`
  along their real paths, and drag the sizes of spheres, cylinders, cubes and
  polyhedra straight from the 3D view, with the code updated for you.
- **[NumPy everywhere](https://www.pythonscad.org/reference/math/#numpy-support)**
  – NumPy arrays and tuples are accepted wherever lists are.
- **[Better multi-object 3MF export](https://www.pythonscad.org/reference/display/#multi-object-3mf-export)**,
  a new **Export / Export as** workflow,
  **[Customizer parameter containers](https://www.pythonscad.org/reference/customizer/#customizer-containers)**
  and many more improvements.

Read on for the full list.

## Release name

Following our tradition, release names are 3D-printing terms translated
(loosely or directly) to Latin. v1.0.0 was "Impletio" (infill), v1.1 was
"Prominentio" (overhang).

This release is called "Quadratum" which means "square".

## Python API

### New: `patch()`

[`patch()`](https://www.pythonscad.org/reference/primitives3d/#patch)
(formerly developed as a generalised `loft()`) builds surfaces and solids from
an outer boundary and optional holes. Besides classic lofting it supports
tangent-driven curvature and three modes – **tube**, **panel** and
**ring-bridge** – and accepts 2D shapes as input
([PR #1054](https://github.com/pythonscad/pythonscad/pull/1054),
[PR #950](https://github.com/pythonscad/pythonscad/pull/950)).
`patch(displacement=None)` is accepted as well
([PR #1146](https://github.com/pythonscad/pythonscad/pull/1146)).

### NumPy support across the API

NumPy arrays and tuples can now be passed wherever lists are accepted. New
[`Vector1`/`Vector2`/`Vector3` and `Matrix4x4`](https://www.pythonscad.org/reference/math/#vector1-vector2-vector3-matrix4x4)
helper types are NumPy-backed when NumPy is installed and fall back to plain
lists otherwise, so scripts work either way
([PR #891](https://github.com/pythonscad/pythonscad/pull/891),
[PR #916](https://github.com/pythonscad/pythonscad/pull/916)).
See the [NumPy support](https://www.pythonscad.org/reference/math/#numpy-support)
section in the math reference.

### Customizer

- **[Parameter containers](https://www.pythonscad.org/reference/customizer/#customizer-containers)**:
  explicit, groupable parameter state without global injection; the classic
  [`add_parameter()`](https://www.pythonscad.org/reference/customizer/#add_parameter)
  API keeps working
  ([PR #956](https://github.com/pythonscad/pythonscad/pull/956)).
- **[Custom PyQt6 widgets in the Customizer](https://www.pythonscad.org/development/custom-customizer-widgets/)**
  ([PR #1024](https://github.com/pythonscad/pythonscad/pull/1024)).
- The Customizer panel now keeps its scroll position, focus and expanded groups
  when the design re-renders
  ([PR #976](https://github.com/pythonscad/pythonscad/pull/976)).
- Companion Customizer JSON is reloaded when it changes on disk.
- Integer parameters stay integers instead of turning into floats
  ([PR #972](https://github.com/pythonscad/pythonscad/pull/972)).

### Import and export

- **[Multi-object 3MF from Python](https://www.pythonscad.org/reference/display/#multi-object-3mf-export)**:
  top-level objects can be kept separate and are exported as individually named
  3MF objects (using the Python variable name where available)
  ([PR #1032](https://github.com/pythonscad/pythonscad/pull/1032)).
- **[`osimport(..., split_by_color=True)`](https://www.pythonscad.org/reference/io/#osimport)**
  splits a multi-color SVG into a dict of `{color: 2D object}` – perfect for
  multi-material prints
  ([PR #887](https://github.com/pythonscad/pythonscad/pull/887)):

  ```python
  parts = osimport("logo.svg", split_by_color=True)
  solids = {name: obj.linear_extrude(2) for name, obj in parts.items()}
  export(solids, "logo.3mf")
  ```

- **[`MultiToolExporter`](https://www.pythonscad.org/reference/multitool/#multitoolexporter)**
  - hidden "cutter" entries take part in subtraction without being exported
    ([PR #947](https://github.com/pythonscad/pythonscad/pull/947))
  - `prefix`, `suffix` and `mkdir` move to `export()`; the constructor
    arguments are deprecated but still work
    ([PR #1077](https://github.com/pythonscad/pythonscad/pull/1077))
  - parts are now rendered for F5 preview
- STEP export writes proper units
  ([PR #1128](https://github.com/pythonscad/pythonscad/pull/1128)).
- [`export()`](https://www.pythonscad.org/reference/display/#export) no longer
  writes files during F5 preview
  ([PR #1010](https://github.com/pythonscad/pythonscad/pull/1010)).

### Other API improvements

- [`rounded_cube()`](https://www.pythonscad.org/reference/primitives3d/#rounded_cube)
  gained a `center` option
  ([PR #912](https://github.com/pythonscad/pythonscad/pull/912)).
- 2D [`.position` / `.size`](https://www.pythonscad.org/reference/properties/#size)
  now return proper vectors, just like in 3D, so
  e.g. `obj.translate(obj.position * -1)` works for 2D objects
  ([PR #1162](https://github.com/pythonscad/pythonscad/pull/1162)).
- [`osuse()`](https://www.pythonscad.org/reference/io/#osuse) falls back to the
  OpenSCAD library path after the local lookup
  ([PR #967](https://github.com/pythonscad/pythonscad/pull/967)).
- Calltips are now available for all functions
  ([PR #926](https://github.com/pythonscad/pythonscad/pull/926)).
- Type information (stubs) now ships with the wheels for better IDE support
  ([PR #958](https://github.com/pythonscad/pythonscad/pull/958)).
- An installable **IDE shim** package provides a placeholder `openscad` module
  with type information for editor completion when you are not running inside
  PythonSCAD itself.

### Deprecations

- **[`nimport()` is deprecated](https://www.pythonscad.org/reference/io/#nimport).**
  Please publish and install libraries as regular, versioned Python packages via
  pip/venv instead. See
  [Creating and publishing a library](https://www.pythonscad.org/libraries/creating/)
  ([PR #1083](https://github.com/pythonscad/pythonscad/pull/1083)).
- [`output()`](https://www.pythonscad.org/reference/display/#show) is deprecated
  in favour of [`show()`](https://www.pythonscad.org/reference/display/#show);
  all website examples were updated
  ([PR #1209](https://github.com/pythonscad/pythonscad/pull/1209)).
- MultiToolExporter constructor arguments `prefix`, `suffix`, `mkdir`
  (see above).

## PyQt6 integration and forms

PythonSCAD now ships
[PyQt6 bindings](https://www.pythonscad.org/development/pyqt6-and-forms/) built
against its own Qt6
([PR #904](https://github.com/pythonscad/pythonscad/pull/904),
[PR #915](https://github.com/pythonscad/pythonscad/pull/915)). This opens up a
whole new way of designing:

- A Python class can define an
  **[interactive dialog](https://www.pythonscad.org/development/extend-with-interactive-form/)**
  (e.g. a point-editing canvas for a polygon) which opens when you click its
  calltip in the editor – the edited values are written back into your code.
- The new **`pyscadforms`** library provides ready-made building blocks.
- Together with
  [`.pythonscadrc` and `add_menuitem()`](https://www.pythonscad.org/development/pythonscadrc-and-menuitems/)
  you can extend the PythonSCAD GUI with your own tools. `.pythonscadrc` is now
  loaded from the PythonSCAD config directory.

On Linux, packages use the distribution's PyQt6, AppImages bundle it.
See the new developer documentation on
[pythonscad.org](https://www.pythonscad.org/development/pyqt6-and-forms/).

## GUI

- **Interactive dragging reworked**
  ([PR #1168](https://github.com/pythonscad/pythonscad/pull/1168)):
  - drag `rotate`, `scale` and `translate` along their real path
  - drag sphere, cylinder, polyhedron and scalar cube sizes
  - support for [`marked()`](https://www.pythonscad.org/reference/gui/#marked)
    values; variables are changed at their definition
  - smoother point dragging with a clearly visible (azure) drag handle
- **Export… / Export as…** (GIMP-style): remembers path, format and options per
  tab, replacing the per-format export submenu
  ([PR #1042](https://github.com/pythonscad/pythonscad/pull/1042)).
- **Render-first prompts for Export (F7) and 3D Print (F8)**: if the design has
  not been F6-rendered since the last change (or never rendered), a dialog now
  offers **Render and Export** / **Render and Print** instead of failing. See
  [Slicer integration](https://www.pythonscad.org/slicer-integration/) for
  configuring a local slicer with F8
  ([PR #1042](https://github.com/pythonscad/pythonscad/pull/1042)).
- **Fullscreen mode** (F11) and unified window-geometry persistence
  ([PR #995](https://github.com/pythonscad/pythonscad/pull/995)).
- Creating a
  [Python virtual environment from the GUI](https://www.pythonscad.org/libraries/using/#gui-file-python-virtual-environment)
  is much more robust (runs in a subprocess, bootstraps pip when needed,
  Windows fixes)
  ([PR #1219](https://github.com/pythonscad/pythonscad/pull/1219)).
- **Nord** color schemes and more consistent theme / face-color handling in
  preview and render
  ([upstream #7033](https://github.com/openscad/openscad/pull/7033),
  [PR #1217](https://github.com/pythonscad/pythonscad/pull/1217)).
- Search state is preserved across tabs, rendering is cancelled before closing
  a tab, and the icon theme survives OS theme changes.
- `.scad` files are registered as default and `.py` files as "Open With" on
  all platforms ([PR #1154](https://github.com/pythonscad/pythonscad/pull/1154)).

## Geometry fixes

- [Fillets](https://www.pythonscad.org/reference/fillet/) with `r=` now also work
  on **2D** union/difference
  ([PR #1069](https://github.com/pythonscad/pythonscad/pull/1069)).
- Fixed fillet self-intersections at 3-way corners
  ([issue #994](https://github.com/pythonscad/pythonscad/issues/994),
  [PR #1049](https://github.com/pythonscad/pythonscad/pull/1049)).
- [`wrap()`](https://www.pythonscad.org/reference/wrap/) no longer drops slices
  of wrapped solids
  ([PR #1043](https://github.com/pythonscad/pythonscad/pull/1043)).
- [`pull()`](https://www.pythonscad.org/reference/pull/) no longer drops faces on
  axis-aligned cuts and keeps pulled cavities in the final render
  ([PR #1131](https://github.com/pythonscad/pythonscad/pull/1131),
  [PR #1201](https://github.com/pythonscad/pythonscad/pull/1201)).
- Correct circular segment count in fillets
  ([issue #985](https://github.com/pythonscad/pythonscad/issues/985)).
- Degenerate triangles no longer cause asserts on export and are skipped in 3MF
  ([PR #1114](https://github.com/pythonscad/pythonscad/pull/1114),
  [PR #1163](https://github.com/pythonscad/pythonscad/pull/1163)).
- Deterministic results: no more `rand()` in geometry code
  ([PR #1053](https://github.com/pythonscad/pythonscad/pull/1053)).
- SVG `fill:none` is treated as unpainted instead of invisible geometry
  ([PR #987](https://github.com/pythonscad/pythonscad/pull/987)).
- Fixes for `multmatrix`, invalid color names and a crash in `color()`.

## Web playground (WebAssembly)

- Complete models render without a CGAL fallback
  ([PR #935](https://github.com/pythonscad/pythonscad/pull/935)).
- The playground stays responsive during evaluation, the debug log is hidden by
  default, and navigation is clearer.
- A reusable WASM toolchain image is published on GHCR.

## Translations

All translations were completed, and **Latin** and **Esperanto** were added
([PR #951](https://github.com/pythonscad/pythonscad/pull/951),
[PR #1078](https://github.com/pythonscad/pythonscad/pull/1078)).

## Merged from upstream OpenSCAD

We kept merging upstream OpenSCAD regularly. Notable changes include:

- Unicode (non-ASCII) identifiers following UAX #31
- `import()` derives the file type from the extension and allows specifying it
  explicitly
- [`sphere(style = "octa")`](https://www.pythonscad.org/reference/primitives3d/#sphere)
  octahedral tessellation (default remains `"orig"`)
- `fragments(r|d)` returns how many sides a circle of that size would get under
  the current [`fn` / `fa` / `fs`](https://www.pythonscad.org/reference/variables/)
  settings (also available from Python)
- `object()` and JSON `import()` are no longer experimental
- AMF import/export removed; OFF import improvements and header fix
- AI assistant: tool calling, chat export/import, viewport tools and agent modes
- Nord color schemes
- Keyboard shortcuts for all menu items, move line up/down shortcut
- Custom CA certificates via PEM file, IPv6 URL fix
- SVG: entity handling fix and limited `<use>` expansion
- Warnings when parameter limits are adjusted, NFC-normalised parameter set
  names
- Open files with non-ASCII paths on Windows
- Updated German, French, Swedish and Galician translations

## Packaging and platforms

- **macOS**: releases are built for macOS 15 on Apple Silicon (Intel runners
  were retired); fixed an Apple Silicon startup crash.
- **Linux**: AppImage is now built on Ubuntu 24.04 and is **Qt6-only** (the Qt5
  AppImage is gone); Debian 11 (bullseye) is no longer supported; supported
  distribution targets were refreshed (including newer Ubuntu and Enterprise
  Linux releases)
  ([PR #1252](https://github.com/pythonscad/pythonscad/pull/1252)).
- **Windows**: improved Python import library generation and venv handling.
- **PyPI** wheels ship a real lib3mf on manylinux and include type information.
- Updated Nix build (thanks to @mvirts) and macOS dependency build fixes
  (thanks to @javax0).
- New release smoke tests run against the packaged artifacts on clean runners.

## Thank you

A big thank you to everyone who contributed to this release – especially (in
alphabetical order) dasbidyendu, gsohler, javax0, John Neo, mvirts, nomike,
pinkfish and Shadowjump – and to everyone who reported bugs and tested
pre-releases.

Download PythonSCAD v1.2.0 from
[pythonscad.org](https://www.pythonscad.org/downloads/). If you've installed
PythonSCAD in a Python virtual environment, simply run
`pip install --upgrade pythonscad`.

Happy designing!
