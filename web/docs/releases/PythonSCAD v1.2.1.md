# PythonSCAD v1.2.1 released

Hi everyone,

**PythonSCAD v1.2.1** is a maintenance release addressing an issue that caused
a crash when using Python callbacks with `path_extrude()`.

The release also hardens `twist` and `profile` callbacks, improves callback
error handling, and fixes closed-path indexing.

Thanks to [@ianjghent](https://github.com/ianjghent) for reporting the issue
in [#1264](https://github.com/pythonscad/pythonscad/issues/1264) and to
[@gsohler](https://github.com/gsohler) and [@nomike](https://github.com/nomike)
for fixing it.

Download PythonSCAD v1.2.1 from
[pythonscad.org](https://www.pythonscad.org/downloads/). If you've installed
PythonSCAD in a Python virtual environment, simply run
`pip install --upgrade pythonscad`.

Happy designing!
