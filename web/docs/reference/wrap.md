# wrap

Wrap a flat object around a cylinder. This transforms a planar shape so that it conforms to a cylindrical surface.

**Syntax:**

```python
wrap(obj, target, r=None, d=None, fn=0, fa=0, fs=0)
obj.wrap(target, r=None, d=None, fn=0, fa=0, fs=0)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `obj` | solid | — | The flat object to wrap |
| `target` | solid | — | The target cylinder to wrap around |
| `r` | float | `None` | Cylinder radius (alternative to providing a target) |
| `d` | float | `None` | Cylinder diameter |
| `fn`, `fa`, `fs` | float | `0` | Curve discretization parameters |

**Examples:**

```python
from openscad import *

flat_text = text("Hello", size=5).linear_extrude(height=1)
cyl = cylinder(r=10, h=20)
flat_text.wrap(cyl).show()
```
