"""Regression fixture for the in-script ``export()`` function: AMF.

AMF is deprecated -- ``src/io/export_amf.cc`` logs
"AMF export is deprecated. Please use 3mf instead." -- but the format
identifier is still wired and the writer is still functional, so we
test it for parity with the upstream ``export-amf`` ctest. The driver
runs ``post_process_amf`` on the file (rewriting the
``<metadata type="producer">PythonSCAD <version></metadata>`` line to a
fixed placeholder) before the normalized text compare so the golden
survives version bumps.
"""
from pythonscad import cube, export

export(cube(5), "cube.amf")
