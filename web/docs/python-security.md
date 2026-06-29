# Python Security

PythonSCAD designs can be written in Python. Python is powerful, but that also
means a design can do much more than describe geometry.

For example, a Python design running with normal desktop access could read,
change, or delete files that your user account can access:

```python
from pythonscad import *
import shutil

cube(10).show()
shutil.rmtree("/path/to/important/files")
```

It could also open network connections, inspect environment variables, run other
programs, or overwrite files:

```python
with open("/home/alex/model.stl", "w") as f:
    f.write("replaced")
```

## Sandboxed Python

Sandboxed Python runs the design inside the WebAssembly build of PythonSCAD.
That environment has its own virtual filesystem and cannot directly access your
home directory or system files. This is the recommended mode for designs you
download, receive from someone else, or have not reviewed carefully.

Sandboxed mode is designed so ordinary modeling code can run without asking you
to trust the design first.

## Native Python

Native Python uses the Python interpreter embedded in the desktop application.
It is useful for compatibility with scripts that need local Python packages,
custom filesystem access, or other desktop integration.

Only use native mode for designs you trust. In native mode, a design has the
same filesystem permissions as PythonSCAD itself.

## Exported Files

In sandboxed mode, files produced by Python code are first written inside the
virtual filesystem. PythonSCAD can then show those files to you and copy only
the files you choose to a destination on your computer.

This extra step is intentional. It prevents a design from silently overwriting a
file such as an existing model, document, or configuration file.
