# Design Rule Checks

Design rule checks find the places where a part cannot be manufactured.
They are the 3D counterpart of the rules known from chip layout DRC:

| Function | Checks | DRC rule |
|----------|--------|----------|
| `internal()` | minimum wall thickness | `INTERNAL` |
| `external()` | minimum gap | `EXTERNAL` |
| `slope()` | face angle against a direction | `ANGLE` |
| `overhang()` | overhangs for 3D printing | |
| `draft()` | draft angle and undercuts for molding and casting | |

All of them return an **error solid**: the region that violates the rule,
empty when the design is clean. It can be colored and shown on top of the
part, subtracted, measured or exported like any other object. A summary is
written to the console.

**How internal and external measure:** two faces "face" each other when their normals are at
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

## slope

Face angle check against a direction, the 3D counterpart of an `ANGLE` rule.
For every face the angle `beta = asin(n · dir)` is computed:

| beta | Face |
|------|------|
| `0` | parallel to `dir`, e.g. a vertical wall when `dir` points up |
| `90` | looks along `dir` (top face) |
| `-90` | looks against `dir` (bottom face) |

Faces outside `[min, max]` are violations. `overhang` and `draft` are the
same test with fixed windows.

**Syntax:**

=== "Python"

    ```python
    slope(obj, dir=[0,0,1], min=None, max=None, undercut=False, grow=None, report=False)
    obj.slope(dir=[0,0,1], min=None, max=None, undercut=False, grow=None, report=False)
    ```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `dir` | vector | `[0,0,1]` | Reference direction |
| `min`, `max` | float | `-90`, `90` | Allowed range of `beta` in degrees |
| `undercut` | bool | `False` | Also report faces looking along `dir` that are hidden behind other material in that direction |
| `grow`, `report` | | | As for `internal`. The report is `{"count", "angle", "undercut", "worst"}` |

---

## overhang

Overhang check for 3D printing: faces that look down more than `angle`
degrees away from the vertical. Faces on the build plate (the lowest level
along `dir`) are ignored.

=== "Python"

    ```python
    overhang(obj, angle=45, dir=[0,0,1], grow=None, report=False)
    obj.overhang(angle=45, dir=[0,0,1], grow=None, report=False)
    ```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `angle` | float | `45` | Largest allowed overhang, measured from the vertical |
| `dir` | vector | `[0,0,1]` | Build direction |

In the report, `worst` is the largest overhang angle found.

=== "Python"

    ```python
    from pythonscad import *

    mushroom = cylinder(h=10, r=2) | cylinder(h=3, r=8).translate([0, 0, 10])
    show([mushroom.color("lightgray", 0.3), mushroom.overhang().color("red")])
    print(mushroom.overhang(report=True))   # worst: 90.0, the underside of the cap
    ```

---

## draft

Draft and undercut check for molding and casting (sand casting, die
casting, injection molding). The part is pulled out of a two-part mold:
the upper half along `dir`, the lower half along `-dir`.

Two kinds of errors are reported:

- **Angle:** a face has less than `angle` degrees of draft towards its mold
  half. Faces that lean the wrong way (negative draft) are local undercuts.
- **Undercut:** a face has a correct angle but is hidden behind other
  material along its pull direction, e.g. the inner faces of a C-profile
  lying on its side. The error solid then shows the trapped sand between
  the face and the material in front of it.

=== "Python"

    ```python
    draft(obj, angle=2, dir=[0,0,1], parting=None, undercut=True, grow=None, report=False)
    obj.draft(angle=2, dir=[0,0,1], parting=None, undercut=True, grow=None, report=False)
    ```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `angle` | float | `2` | Minimum draft angle in degrees |
| `dir` | vector | `[0,0,1]` | Pull direction of the upper mold half |
| `parting` | float | `None` | Position of a flat parting plane along `dir`. Faces above it belong to the upper half, faces below to the lower half; faces crossing it are split. `None`: free parting line, every face belongs to the half it looks at |
| `undercut` | bool | `True` | Detect undercuts |

=== "Python"

    ```python
    from pythonscad import *

    C = cube([20, 20, 3]) | cube([20, 20, 3]).translate([0, 0, 10]) | cube([3, 20, 13])
    print(C.draft(-0.5, report=True))                # 4 undercut faces
    print(C.draft(-0.5, dir=[1, 0, 0], report=True)) # pulled sideways: clean
    show([C.color("lightgray", 0.3), C.draft(-0.5).color("red")])
    ```

`angle=-0.5` in the example allows the vertical walls, so only the undercut
remains.

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
