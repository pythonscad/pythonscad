from openscad import *

if not isinstance(qapp_ptr(), int):
    raise AssertionError("qapp_ptr() must return an integer pointer")
if not isinstance(mainwindow_ptr(), int):
    raise AssertionError("mainwindow_ptr() must return an integer pointer")

for callback in (qapp_ptr, mainwindow_ptr):
    try:
        callback(1)
    except TypeError:
        pass
    else:
        raise AssertionError("positional arguments must be rejected")

    try:
        callback(unexpected=True)
    except TypeError:
        pass
    else:
        raise AssertionError("keyword arguments must be rejected")

print("GUI callback ABI: OK")
