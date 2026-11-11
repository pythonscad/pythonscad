"""Regression for issue #1107: patch(displacement=None) must be accepted.

Documented signature allows displacement=None, but explicit None was
rejected because it is non-null and not a PyFunction.
"""
from pythonscad import *

outer = [
    [0, 0, 0],
    [20, 0, 0],
    [20, 20, 0],
    [0, 20, 0],
]
patch(outer, displacement=None).show()
