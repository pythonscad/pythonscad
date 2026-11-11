"""Echo test for ``MultiToolExporter.show(render=...)``.

Exercises:
  * default → one plain CSG part per exportable item (no render node)
  * ``render=True`` → every shown part wrapped in ``render()``
  * ``export=False`` cutters omitted from ``show`` in both modes
  * ``parts()`` geometry unchanged by ``show(render=True)``
"""

import pythonscad
from pythonscad import MultiToolExporter, cube

red = cube(10)
blue = cube(10).right(5)
cutter = cube(4).right(2)

exp = MultiToolExporter(items=[
    ("red", red),
    ("blue", blue),
    ("slot", cutter, False),
])


def capture_show(fn):
    shown = []
    real_show = pythonscad.show
    pythonscad.show = lambda obj: shown.append(repr(obj).lstrip())
    try:
        fn()
    finally:
        pythonscad.show = real_show
    return shown


def is_rendered(text):
    return text.startswith("render(")


shown_default = capture_show(lambda: exp.show())
print("default count:", len(shown_default))
print("default rendered:", [is_rendered(s) for s in shown_default])

shown_render = capture_show(lambda: exp.show(render=True))
print("render count:", len(shown_render))
print("render rendered:", [is_rendered(s) for s in shown_render])

print("parts rendered:", [is_rendered(repr(g).lstrip()) for _n, g in exp.parts()])
