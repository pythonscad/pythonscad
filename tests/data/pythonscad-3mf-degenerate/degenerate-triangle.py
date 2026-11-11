"""3MF export must drop triangles whose vertex indices are not unique.

lib3mf rejects those faces ("The provided Index is invalid") and aborts
the whole file, which left per-part exports as a 0-byte 3MF. Two ways a
face gets there:

* a face that already repeats an index
* vertex welding in ``createSortedPolySet`` (always on for lib3mf 1, and
  on for lib3mf 2 when ``--enable=predictible-output`` is set), which
  collapses coincident points onto one index

This fixture is registered with that flag so both lib3mf backends weld
the coincident vertex before writing. The degenerate faces sit between
real faces of a colored shell so a skipped face cannot shift the
lib3mf triangle ids used for the remaining per-triangle colors.
"""
from pythonscad import export, polyhedron

points = [
    [0, 0, 0],
    [10, 0, 0],
    [10, 10, 0],
    [0, 10, 0],
    [0, 0, 10],
    [10, 0, 10],
    [10, 10, 10],
    [0, 10, 10],
    [0, 0, 0],  # coincident with vertex 0; welding repeats an index
]
faces = [
    [0, 1, 8],  # welds to a repeated index
    [0, 1, 2],
    [0, 2, 3],
    [0, 1, 0],  # explicit repeated index, between real faces
    [4, 6, 5],
    [4, 7, 6],
    [0, 4, 5],
    [0, 5, 1],
    [1, 5, 6],
    [1, 6, 2],
    [2, 6, 7],
    [2, 7, 3],
    [3, 7, 4],
    [3, 4, 0],
]
export({"shell": polyhedron(points, faces).color("red")}, "degenerate-triangle.3mf")
