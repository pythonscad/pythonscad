"""Regression test for osuse() lookup through OPENSCADPATH."""

from openscad import osuse, scad

library = osuse("pythonscad_osuse_library_path_unique.scad")
library.library_box()
scad(f'echo("library path answer: {library.library_answer()}");')

local = osuse("pythonscad_osuse_library_path_precedence.scad")
scad(f'echo("local path answer: {local.library_answer()}");')
