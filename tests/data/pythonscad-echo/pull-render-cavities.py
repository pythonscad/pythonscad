"""Regression #1200: pulled surfaces must survive F6 backend conversion.

Use only primitives: a square tube with a through-hole, pulled along each axis.
Check the raw mesh's orientation as well as its volume after a boolean union.
The old implementation reversed some faces and filled the hole with a hull.
"""
from collections import Counter
from openscad import *


def check_mesh(obj, expected_volume):
    points, faces = obj.mesh(triangulate=True)
    edges = Counter()
    volume = 0.0
    for face in faces:
        a, b, c = [points[i] for i in face]
        volume += (
            a[0] * (b[1] * c[2] - b[2] * c[1])
            + a[1] * (b[2] * c[0] - b[0] * c[2])
            + a[2] * (b[0] * c[1] - b[1] * c[0])
        ) / 6
        for i, j in zip(face, face[1:] + face[:1]):
            edges[i, j] += 1
    assert all(count == 1 and edges[j, i] == 1 for (i, j), count in edges.items()), \
        "mesh must be closed with consistently oriented faces"
    # CGAL quantizes decimal coordinates during backend conversion.
    assert abs(volume - expected_volume) < 1e-3, (volume, expected_volume)


tube = cube([10, 10, 10]) - cube([6, 6, 12]).translate([2, 2, -1])
for direction, volume in [([4, 0, 0], 800), ([0, 4, 0], 800), ([0, 0, 4], 896)]:
    pulled = tube.pull([5, 5, 5], direction)
    check_mesh(pulled, volume)
    # A second object forces conversion of the pulled PolySet into the selected
    # geometry backend, reproducing the final union used by F6.
    check_mesh(pulled + cube(1).translate([20, 0, 0]), volume + 1)
    print("pull", direction, "closed, oriented, cavity preserved")

# Decimal coordinates make opposite traversals of a shared cut edge round
# differently unless the intersection uses a canonical endpoint order.
offset = [0.123, 0.456, 0.789]
shifted = tube.translate(offset)
for direction, volume in [([4, 0, 0], 800), ([0, 4, 0], 800), ([0, 0, 4], 896)]:
    pulled = shifted.pull([5 + value for value in offset], direction)
    check_mesh(pulled, volume)
    check_mesh(pulled + cube(1).translate([20, 0, 0]), volume + 1)
    print("fractional pull", direction, "closed, oriented, cavity preserved")
