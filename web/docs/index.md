# PythonSCAD

**Design precise 3D models with Python.**

PythonSCAD is a script-based 3D modeling application with a full GUI. Write
parametric, engineering-oriented models in Python, preview them live, and export
to STL, 3MF, and other formats for 3D printing and manufacturing.

<div id="hero-download">Loading latest release…</div>

[Get started](get_started.md){ .md-button }
[Tutorial](tutorial/getting_started.md){ .md-button .md-button--primary }
[All downloads](downloads.md)

![Step-by-step: building a camera housing in the PythonSCAD GUI](pictures/box_anim.gif)

## Is PythonSCAD for you?

### A great fit if you…

- think in code and want precise, repeatable parametric models
- already know Python — or want to **learn Python or programming** through something tangible
- need models for 3D printing, CNC, or engineering workflows
- want the full Python ecosystem: pip packages, your IDE, Jupyter notebooks

### Probably not the right tool if you…

- need organic sculpting, animation, or VFX → [Blender](https://www.blender.org/)
- prefer click-to-design CAD → [FreeCAD](https://www.freecad.org/)

## Why PythonSCAD?

- **Python-native** — write `.py` scripts directly in the app; no external code generator
- **Instant visual feedback** — see your code become a 3D model as you work; export and **3D-print something you wrote**
- **Rewarding for learners** — a few lines of Python can produce a useful physical object, which makes programming feel concrete
- **Live preview** — interactive 3D view, customizer, and export in one desktop app
- **Powerful geometry** — booleans, extrusions, fillets, surfaces, and more
- **Cross-platform** — Windows, Linux, macOS, plus a browser WebAssembly build
- **Open and extensible** — open source; integrate libraries, scripts, and tooling you already use

## From code to model

=== "Python"

    ```python
    from pythonscad import *

    c = cube([10, 20, 30]).color("Tomato")
    show(c)
    ```

=== "OpenSCAD-compatible"

    ```python
    from openscad import *

    c = cube([10, 20, 30]).color("Tomato")
    show(c)
    ```

![Example model rendered in PythonSCAD](pictures/red-box-example.png)

Browse the [Examples gallery](examples.md) for QR codes, gyroids, GDS import, and more.

## Learn more

1. [Download](downloads.md) for your platform — or use the button above
2. [Installation](installation.md) — APT/YUM repos, AppImage, PyPI, and build options
3. [Getting Started tutorial](tutorial/getting_started.md) — hands-on walkthrough
4. [Cheat sheet](cheatsheet.md) and [API reference](reference/primitives3d.md)

## Community

- [r/OpenPythonSCAD](https://www.reddit.com/r/OpenPythonSCAD/)
- [Google Group](https://groups.google.com/g/pythonscad)
- [Community wiki](http://old.reddit.com/r/openpythonscad/wiki/index)
- [Python type stubs for IDE support](https://raw.githubusercontent.com/pythonscad/pythonscad/refs/heads/master/libraries/python/stubs/openscad/__init__.pyi)
- [Contact & support](contact.md)

---

PythonSCAD is open source on [GitHub](https://github.com/pythonscad/pythonscad). It builds on
the solid foundation of [OpenSCAD](https://openscad.org) and stays closely synced with upstream;
Python-specific features are developed here. See [Upstream sync](development/upstream-sync.md).
