# .pythonscadrc and add_menuitem

`.pythonscadrc` is PythonSCAD's Python startup file — the same idea as a
shell's `.bashrc`. It is loaded once, automatically, when PythonSCAD starts,
and anything it defines stays available for the entire session: functions,
classes, imports, custom menu items.

This page covers what `.pythonscadrc` is for, how its contents persist across
renders, and `add_menuitem`, the function it typically uses to add custom
commands to PythonSCAD's menus.

## Where It Lives and When It Runs

PythonSCAD looks for `.pythonscadrc` at startup and, if found, runs its
contents in the same Python namespace that your scripts run in. This happens
once per session — editing the file while PythonSCAD is already running has
no effect until you restart.

A minimal `.pythonscadrc`:

```python
from pythonscad import *

def double(x):
    return x * 2
```

After this loads, `double(21)` is available in any script you open, for the
rest of the session.

## Persistence Across Renders

PythonSCAD resets certain Python state between renders, so that leftover
variables from one script don't leak into the next. Anything defined while
`.pythonscadrc` loads is exempted from this reset — it is treated as part of
the session's baseline, not as script-local state. That is what makes it the
right place for anything you want to rely on being there every time,
regardless of which document is open or how many times it has been rendered.

Definitions made only inside an individual script are *not* protected this
way. They are visible to other open documents too (every open document shares
one Python namespace), but only after that script has been rendered, and only
until a different document without the same definitions is rendered next.
For anything you want reliably available — helper functions, custom classes,
menu items — put it in `.pythonscadrc`, or in a library `.pythonscadrc`
imports.

> **Watch out for repeated imports.** If a script itself starts with
> `from pythonscad import *`, rendering that script re-imports the native
> primitives into the same shared namespace `.pythonscadrc` uses. If you have
> redefined any of those names (see
> [Extend your design with an Interactive Form](extend-with-interactive-form.md)),
> this silently undoes the redefinition for every open document, not just the
> one being rendered. Keep that particular import in `.pythonscadrc` only,
> not repeated in every script.

## setup() — Where add_menuitem Belongs

When `.pythonscadrc` loads, PythonSCAD calls a function named `setup()` in
it, once, at startup. This is the intended place to register menu items —
`add_menuitem` calls made there run exactly once per session, when the menu
they add to already exists.

```python
from pythonscad import *

def my_action():
    print("Running my_action")

def setup():
    add_menuitem("Tools", "Run my action", "my_action()")
```

Calling `add_menuitem` outside `setup()` — for example directly at the top
level of `.pythonscadrc` — is not the supported pattern; use `setup()` for
anything menu-related.

> **Note:** Older versions of PythonSCAD reloaded `.pythonscadrc` on every
> menu click, which made it possible to edit the file and see changes without
> restarting. That is no longer the case — `.pythonscadrc` loads once at
> startup, `setup()` runs once as part of that, and menu clicks only run
> their callback. Editing `.pythonscadrc` requires restarting PythonSCAD to
> take effect, the same as any other change to it.

## Summary

- `.pythonscadrc` loads once at startup, into the same namespace every script
  uses.
- Definitions made there persist across renders; script-local definitions do
  not, reliably.
- `add_menuitem` adds a custom command to a menu; its callback string is
  evaluated in the shared namespace.
- Call `add_menuitem` from inside `setup()`, the function PythonSCAD calls
  once when `.pythonscadrc` loads at startup.
