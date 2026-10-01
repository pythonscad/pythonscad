# Design Rule Checks

`internal()` and `external()` check a solid for walls that are too thin and
gaps that are too narrow. They are the 3D counterpart of the `INTERNAL` and
`EXTERNAL` rules known from chip layout DRC, and are useful wherever a
process has a minimum wall thickness or clearance: 3D printing, die casting,
injection molding, CNC.

Both return an **error solid**: the material (`internal`) or the air
(`external`) that violates the rule. It is empty when the design is clean, so
it can be colored and shown on top of the part, subtracted, measured or
exported like any other object. A summary is written to the console.

**How it measures:** two faces "face" each other when their normals are at
least `angle` degrees apart and each lies behind the other one. Both faces
are clipped exactly to the part that is closer than `d`, and the exact
distance between them is computed (vertex to face and edge to edge). There is
no sampling or voxelization, so results are exact for the mesh.

---

## internal

Minimum wall thickness, measured through the material.

**Syntax:**

=== "Python"

    ```python
    internal(obj, d, angle=120, alpha=90, occlusion=True, grow=None, report=False)
    obj.internal(d, angle=120, alpha=90, occlusion=True, grow=None, report=False)
    ```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `obj` | solid | — | The object to check |
| `d` | float | — | Minimum allowed wall thickness |
| `angle` | float | `120` | Minimum angle between the two face normals, `90 < angle <= 180`. `180` only measures exactly opposite faces. With `120`, wedges sharper than 60° count as thin |
| `alpha` | float | `90` | `0` measures strictly perpendicular to the faces (like `PROJECTING`), larger values allow oblique measurement up to that angle, `>= 90` any direction |
| `occlusion` | bool | `True` | Ignore pairs whose connecting line runs through air, e.g. from edge to edge across a step |
| `grow` | float | `None` | Lifts the error solid off the part by about this distance (each face moves outwards, vertices are moved along their normals, not an exact offset), so both can be shown together without z-fighting. `None`: 0.1 % of the part size, `0`: exact error solid |
| `report` | bool | `False` | Evaluate immediately and return `{"count": n, "min": d}` instead of the error solid. `min` is `None` when clean |

**Examples:**

=== "Python"

    ```python
    from pythonscad import *

    part = cube([20, 20, 10]) - cube([18, 18, 10]).translate([1, 1, 1.2])
    thin = part.internal(1.5)

    show([part.color("lightgray", 0.3), thin.color("red")])
    print(part.internal(1.5, report=True))   # {'count': 20, 'min': 1.0}
    ```

---

## external

Minimum spacing, measured through the air.

**Syntax:**

=== "Python"

    ```python
    external(obj, d, other=None, angle=120, alpha=90, occlusion=True, grow=None, report=False)
    obj.external(d, other=None, angle=120, alpha=90, occlusion=True, grow=None, report=False)
    ```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `obj` | solid | — | The object to check. Gaps inside it are checked, also between separate bodies it contains |
| `d` | float | — | Minimum allowed gap |
| `other` | solid | `None` | If given, only gaps between `obj` and `other` are checked, gaps inside either of them are ignored |
| `angle`, `alpha`, `occlusion`, `grow`, `report` | | | As for `internal` |

**Examples:**

=== "Python"

    ```python
    from pythonscad import *

    a = cube(10)
    b = cube(10).translate([10.3, 0, 0])

    gap = a.external(0.5, other=b)
    show([a, b, gap.color("red")])
    print(a.external(0.5, other=b, report=True))   # min ≈ 0.3
    ```

---

## Choosing angle and alpha

Allowing oblique measurement is stricter than measuring perpendicular only:
across a round hole of diameter 9, faces on opposite sides are 9 apart, but
oblique chords between faces 120° apart are only about 7.8 long.

| Setting | Measures |
|---------|----------|
| `alpha=0` | Only perpendicular distances, like a lateral DRC check |
| `alpha=90, angle=175` | Any direction, but only between nearly opposite faces |
| `alpha=90, angle=120` (default) | Any direction, including wedges sharper than 60° |

## Performance

The check uses a bounding volume hierarchy and only looks at face pairs
closer than `d`, so clean designs are fast (about 2 s for 500k triangles).
Large, finely meshed areas that violate the rule produce many face pairs.
The error solid uses one hull per violating triangle (to its nearest
partner); a fully too thin hollow sphere with 29k triangles takes about 15 s.
Above 20000 violating triangles only the most severe ones are shown and a
warning is printed. `report=True` never builds the error solid and needs
little memory.
