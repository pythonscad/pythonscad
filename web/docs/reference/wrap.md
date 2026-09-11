# wrap

Wrap a flat object around a cylinder or a 2D outline. This transforms a
planar shape so that it conforms to the target surface.

**Syntax:**

=== "Python"

    ```python
    wrap(obj, target=None, r=None, d=None, fn=0, fa=0, fs=0)
    obj.wrap(target=None, r=None, d=None, fn=0, fa=0, fs=0)
    ```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `obj` | solid | — | The flat 3d object to wrap |
| `target` | solid | `None` | The target 2D object to wrap around. Optional if `r` or `d` is given |
| `r` | float | `None` | Cylinder radius (alternative to providing a target) |
| `d` | float | `None` | Cylinder diameter |
| `fn`, `fa`, `fs` | float | global | Curve discretization; defaults to the global `fn`/`fa`/`fs` values |

**Examples:**

=== "Python"

    ```python
    from pythonscad import *

    flat_text = text("Hello", size=5).linear_extrude(height=1).rotx(90)

    # Around a 2D outline
    cyl = circle(r=10)
    flat_text.wrap(cyl).show()

    # Around a cylinder of the given radius
    flat_text.wrap(r=10).show()
    ```
