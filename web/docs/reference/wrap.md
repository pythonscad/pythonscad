# wrap

Wrap a flat object around a cylinder. This transforms a planar shape so that it conforms to a cylindrical surface.

**Syntax:**

=== "Python"

    ```python
    wrap(obj, target=None, r=None, d=None, fn=0, fa=0, fs=0)
    obj.wrap(target=None, r=None, d=None, fn=0, fa=0, fs=0)
    ```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `obj` | solid | — | The flat object to wrap |
| `target` | solid | `None` | The target 2D object to wrap around. Optional if `r` or `d` is given |
| `r` | float | `None` | Cylinder radius (alternative to providing a target) |
| `d` | float | `None` | Cylinder diameter |
| `fn`, `fa`, `fs` | float | global | Curve discretization; defaults to the global `fn`/`fa`/`fs` values |

**Examples:**

=== "Python"

    ```python
    from pythonscad import *

    flat_text = text("Hello", size=5).linear_extrude(height=1).rotx(90)
    cyl = circle(r=10)
    flat_text.wrap(cyl).show()
    ```

    Side note: Wrapping around an 2D object is only working when "wrap" experimental feature is turned on.
    We aim not to distinguish here in the future anymore.
