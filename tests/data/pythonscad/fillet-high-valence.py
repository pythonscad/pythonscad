"""Regression: fillet corners with four and five incident rounded edges."""

from math import cos, pi, sin

from openscad import *


def pyramid(sides, radius, height):
    points = [
        [
            radius * cos(2 * pi * i / sides),
            radius * sin(2 * pi * i / sides),
            0,
        ]
        for i in range(sides)
    ]
    points.append([0, 0, height])
    apex = sides
    faces = [[i, (i + 1) % sides, apex] for i in range(sides)]
    faces.append(list(reversed(range(sides))))
    return polyhedron(points, faces)


def octahedron(radius):
    points = [
        [radius, 0, 0],
        [-radius, 0, 0],
        [0, radius, 0],
        [0, -radius, 0],
        [0, 0, radius],
        [0, 0, -radius],
    ]
    faces = [
        [4, 0, 2],
        [4, 2, 1],
        [4, 1, 3],
        [4, 3, 0],
        [5, 2, 0],
        [5, 1, 2],
        [5, 3, 1],
        [5, 0, 3],
    ]
    return polyhedron(points, faces)


square_pyramid = pyramid(4, 12, 18)
pentagonal_pyramid = pyramid(5, 12, 18)
octahedral_solid = octahedron(12)
partial_mask = cube([20, 2, 30]) + [-1, -1, -1]
noncontiguous_mask = (
    cube(2, center=True).translate([0, 0, 18])
    | cube(2, center=True).translate([12, 0, 0])
    | cube(2, center=True).translate(
        [12 * cos(4 * pi / 5), 12 * sin(4 * pi / 5), 0]
    )
)
two_high_valence_mask = (
    cube(2, center=True).translate([0, 0, 12])
    | cube(2, center=True).translate([12, 0, 0])
)

square_pyramid.fillet(1.5, fn=8).translate([-18, -18, 0]).show()
pentagonal_pyramid.fillet(1.5, fn=8).translate([18, -18, 0]).show()
square_pyramid.fillet(1.5, partial_mask, fn=8).translate([-36, 18, 0]).show()
pentagonal_pyramid.fillet(1.5, noncontiguous_mask, fn=8).translate([0, 18, 0]).show()
octahedral_solid.fillet(1.5, two_high_valence_mask, fn=8).translate(
    [36, 18, 12]
).show()
