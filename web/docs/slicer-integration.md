# Slicer integration

You can send a rendered design straight to your slicer by configuring the
slicer as a local print service. After that, press **F8** (or use
**File → 3D Print**) to export and open the model in the slicer.

If the design has not been rendered with **F6** since the last change,
PythonSCAD offers to render it first (**Render and Print**) before sending
it to the slicer. The same prompt applies to **Export…** (**F7**).

## Setup

1. Open **Edit → Preferences** (or **PythonSCAD → Preferences** on macOS)
   and go to the **3D Print** tab.
2. Enable remote print services and set the default print service to
   **Local Application**.
3. Under **Local Application**, set the executable and parameters for your
   slicer.

The examples below configure [PrusaSlicer](https://www.prusa3d.com/page/prusaslicer_215/).
Other slicers work the same way: set the executable (on macOS the `.app`
bundle; elsewhere the slicer binary) and pass the output file path as a
parameter (use the **File** button to insert the placeholder).

## Linux (Flatpak)

Executable:

```text
/usr/bin/flatpak
```

Parameters (one per line):

- `run`
- `com.prusa3d.PrusaSlicer`
- `--single-instance`
- full path to the output file (use the **File** button)

## Windows

Executable (default installer location):

```text
C:\Program Files\Prusa3D\PrusaSlicer\prusa-slicer.exe
```

Parameters (one per line):

- `--single-instance`
- full path to the output file (use the **File** button)

If you installed PrusaSlicer elsewhere, browse to `prusa-slicer.exe` from
the Preferences dialog instead of typing the path by hand.

## macOS

On macOS, PythonSCAD launches the configured application with
`open -a`, so the executable must be the **application bundle** (the
`.app`), not the Mach-O binary inside `Contents/MacOS/`.

Executable:

```text
/Applications/PrusaSlicer.app
```

Parameters (one per line):

- `--single-instance`
- full path to the output file (use the **File** button)

If PrusaSlicer lives elsewhere (for example
`/Applications/Original Prusa Drivers/PrusaSlicer.app`), point the
executable at that `.app` bundle instead.
