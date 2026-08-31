"""Consumer-side typing contract for the packaged PythonSCAD API."""

from openscad import Openscad, cube
from pythonscad import MultiToolExporter, osimport, rounded_cube


solid: Openscad = cube(10)
rounded: Openscad = rounded_cube([20, 20, 20], r=2)
imported: Openscad = osimport("watermelon.svg", dpi=96)
colored_parts: dict[str, Openscad] = osimport(
    "watermelon.svg",
    split_by_color=True,
    dpi=96,
)

exporter = MultiToolExporter("out/watermelon-", ".stl", mkdir=True)
exporter.extend(colored_parts.items())
