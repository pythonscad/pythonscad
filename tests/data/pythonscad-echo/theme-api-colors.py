"""Public mesh colors resolve symbolic defaults without painting cached geometry."""
from pythonscad import cube, square

material = [249 / 255, 215 / 255, 44 / 255, 1]
cutout = [157 / 255, 203 / 255, 81 / 255, 1]
def matches(actual, expected):
    return all(abs(a - b) < 1e-6 for a, b in zip(actual, expected))

assert matches(square(10).mesh(color=True)[0][0], material)
assert square(10).color("red").mesh(color=True)[0][0] == [1, 0, 0, 1]
shape = cube(2) - cube(1).translate([1, 1, 1])
for _ in range(2):
    _, _, colors, indices = shape.mesh(triangulate=True, color=True)
    assert any(matches(color, material) for color in colors), colors
    assert any(matches(color, cutout) for color in colors), colors
    assert all(0 <= channel <= 1 for color in colors for channel in color), colors
print("Mesh and outline colors resolved")
