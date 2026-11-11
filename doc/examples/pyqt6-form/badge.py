"""Model script for the editor-form demo.

The interactive pieces live in .pythonscadrc, which only the GUI loads, at
startup. This script defines `badge` itself so that it also renders from the
command line, where no editor and no startup file exist.
"""

from pythonscad import cylinder, show, text


if "badge" not in globals():

    def badge(label="PythonSCAD", radius=30.0):
        disc = cylinder(r=radius, h=2, center=True)
        glyph = text(label, size=radius / 3, halign="center", valign="center")
        return disc + glyph.linear_extrude(3).translate([0, 0, -1])


show(badge("PythonSCAD", 30))
