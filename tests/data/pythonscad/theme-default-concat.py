"""Unpainted faces retain their default role through nested concatenation."""
from pythonscad import cube, concat
mixed = concat(cube(1).color('red'), cube(1).translate([2, 0, 0]))
combined = concat(mixed, cube(1).translate([4, 0, 0]))
points, faces, colors, indices = combined.mesh(triangulate=True, color=True)
assert len(faces) == 36, len(faces)
assert len(indices) == len(faces), (len(indices), len(faces))
assert all(color == [1, 0, 0, 1] for color in colors), colors
assert indices[:12] == [0] * 12, indices
assert indices[12:] == [-1] * 24, indices
combined.show()
