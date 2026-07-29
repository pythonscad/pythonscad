from openscad import *

assert isinstance(qapp_ptr(), int)
assert isinstance(mainwindow_ptr(), int)

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
