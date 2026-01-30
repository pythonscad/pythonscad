# Python Editor Code Completion and Parameter Hints – Feasibility and Implementation Plan

## Summary

**Feasibility: Yes.** The same QScintilla API-based completion and call-tip mechanism used for SCAD is available for Python. The Python lexer is already used for `.py` files but has no `QsciAbstractAPIs` attached, so completion and parameter hints are currently inactive for Python. Adding a Python-specific API class and attaching it to the Python lexer will provide equivalent behaviour.

**Effort: Small–medium** (on the order of 1–3 days for an experienced C++/Qt developer), assuming we mirror the existing SCAD behaviour and reuse `Builtins::keywordList` for data.

---

## Current Behaviour (SCAD)

1. **Code completion**
   - `ScadApi` extends `QsciAbstractAPIs` and is created with the SCAD lexer (and editor).
   - It is associated with the lexer (lexer as parent; for SCAD the API is created in `setLexer(ScadLexer*)` / `setLexer(ScadLexer2*)`).
   - `updateAutoCompletionList()` filters `Builtins::keywordList` by the current word prefix and optionally handles `use`/`include` paths.
   - QScintilla uses `AcsAPIs`, so the completion list comes from the lexer’s APIs.
   - Behaviour: type “c” → list with “cube”, “circle”, “color”, etc.; cursor/tab to select.

2. **Parameter hints (call tips)**
   - `ScadApi::callTips()` receives the context (words before cursor). For `cube(`, the function name is `context.at(context.size()-2)`.
   - It returns the list of call-tip strings for that function from the same data used for completion (e.g. `cube(size)`, `cube([width, depth, height])`, `cube([width, depth, height], center = true)`).
   - Call tips are already configured in `ScintillaEditor::setupAutoComplete()`: `setCallTipsVisible(10)`, `setCallTipsStyle(CallTipsContext)`.

3. **Data source**
   - `Builtins::keywordList` is a static map: name → vector of call-tip strings.
   - It is filled from all `Builtins::init(..., calltipList)` calls (primitives, modules, functions) and is shared by SCAD and can be reused for Python.

---

## Why Python Has No Completion Today

- For `.py` files, `EditorInterface::recomputeLanguageActive()` sets `LANG_PYTHON` and `ScintillaEditor::onLanguageChanged(LANG_PYTHON)` runs `qsci->setLexer(this->pythonLexer)`.
- `pythonLexer` is a plain `QsciLexerPython` created in `ScintillaEditor` and **has no `QsciAbstractAPIs`** attached.
- Completion and call tips are driven by the **current lexer’s APIs**. So when the Python lexer is active, there are no APIs → no completion, no parameter hints.

---

## Proposed Solution

### 1. Add a Python API class (`PythonApi`)

- **Role:** Same as `ScadApi` but for the Python lexer: provide completion list and call tips for the `openscad` API used in Python files.
- **Base:** `QsciAbstractAPIs`, constructed with `QsciLexerPython*` (and a pointer to `ScintillaEditor` for current line/column if needed).
- **Data:**
  - Reuse **`Builtins::keywordList`** for all names that exist in both SCAD and Python (e.g. `cube`, `circle`, `sphere`, `color`, `translate`, …). Same call-tip strings (e.g. `cube(size)`, `cube([width, depth, height])`, …).
  - Add a small set of **Python-only** names with call tips, e.g.:
    - `show` → e.g. `show(obj)` or `show(obj, ...)`
    - `add_parameter` → e.g. `add_parameter(name, default, ...)`
    - Optionally: `version`, `version_num`, `scad`, `align`, etc., if desired for the first version.
- **Methods:**
  - `updateAutoCompletionList(context, list)`: same idea as `ScadApi` – filter by `context.last()` (current word prefix). No need for `use`/`include`; optionally skip inside string literals (Python: `'...'`, `"..."`, `'''...'''`) if we want to avoid completion inside strings.
  - `callTips(context, commas, style, shifts)`: same as `ScadApi` – use `context.at(context.size()-2)` as the function name and return the corresponding call-tip list from the same name→signatures map.
- **Files:** New `src/gui/PythonApi.h` and `src/gui/PythonApi.cc` (or similar). Can reuse `ApiFunc` and the same pattern as `ScadApi` for name + list of call-tip strings.

### 2. Attach Python API to the Python lexer

- **Where:** In `ScintillaEditor`, once `pythonLexer` exists (e.g. in constructor or in `initLexer()`), create `pythonApi = new PythonApi(this, pythonLexer)` and call `pythonLexer->setAPIs(pythonApi)` (QScintilla’s `QsciLexer::setAPIs(QsciAbstractAPIs*)`).
- **Lifetime:** `pythonApi` is owned by the editor (or by the lexer if we set the lexer as parent in the API ctor). Ensure it is deleted when the editor is destroyed if not owned by the lexer.
- No change to the way the editor switches lexers: when the user is in a `.py` file, `qsci->setLexer(pythonLexer)` is already called; with APIs set on `pythonLexer`, QScintilla will use them for completion and call tips.

### 3. Reuse existing editor settings

- `setupAutoComplete()` already sets:
  - `setAutoCompletionSource(QsciScintilla::AcsAPIs)`
  - `setAutoCompletionFillups("(")`
  - `setCallTipsVisible(10)` and `setCallTipsStyle(QsciScintilla::CallTipsContext)`
- These apply to whichever lexer is current, so **no change** is required for Python: once the Python lexer has APIs, completion and call tips will work when the Python lexer is active.

### 4. Optional refinements (later)

- **Skip completion inside strings:** In `PythonApi::updateAutoCompletionList`, use current line/column and a simple scan (single/double/triple quotes) to avoid offering completion inside string literals (mirroring SCAD’s `isInString`).
- **Python-only call tips:** Add more entries from the Python `openscad` module (e.g. from `PyOpenSCADFunctions` / `doc/openscad.pyi`) as needed.
- **Method chaining:** Completion after `obj.` (e.g. `obj.translate`, `obj.rotate`) would require context-aware completion (detecting that we are after a dot) and a list of method names; this can be a follow-up.

---

## Implementation Steps

1. **Add `PythonApi` class**
   - Create `PythonApi.h` / `PythonApi.cc` (or equivalent).
   - In constructor, build an internal list of (name, call-tip list) from `Builtins::keywordList` plus Python-only entries (e.g. `show`, `add_parameter` with one or two call-tip strings each).
   - Implement `updateAutoCompletionList()`: filter by prefix from `context.last()`; optionally add string-in-context check.
   - Implement `callTips()`: resolve function name from context (e.g. `context.size() >= 2 ? context.at(context.size()-2) : QString()`), return the corresponding call-tip list.
   - Add to `src/gui/CMakeLists.txt` (or build script) so the new sources are compiled.

2. **Wire Python API into `ScintillaEditor`**
   - In `ScintillaEditor`, after `pythonLexer` is created, instantiate `PythonApi(this, pythonLexer)` and call `pythonLexer->setAPIs(pythonApi)`.
   - Ensure `pythonApi` is kept alive and not double-deleted (e.g. if the API’s parent is the lexer, Qt’s parent-child will delete it with the lexer; otherwise delete in `ScintillaEditor` destructor).
   - Guard with `#ifdef ENABLE_PYTHON` so the feature is only built when Python support is enabled.

3. **Testing**
   - Open a `.py` file, type `c` and confirm completion list (e.g. `cube`, `circle`, `color`, …).
   - Type `cube(` and confirm parameter-hint tooltip with multiple signatures.
   - Test that existing SCAD completion and call tips still work in `.scad` files.
   - Test with “Enable Autocompletion” off and on in preferences.

4. **Documentation**
   - Short note in user docs or cheatsheet that Python files support the same completion and parameter hints as SCAD (and that they show the `openscad` API).

---

## Risk and Limitations

- **Context format:** QScintilla’s `context` for Python may differ slightly (e.g. word boundaries). If `callTips` or completion behaves oddly, verify `context` contents (e.g. for `cube(`) and adjust index (e.g. `context.size()-2`) if needed.
- **No runtime introspection:** Completion is static (from `Builtins::keywordList` + hard-coded Python names). We do not run the script to infer types or methods; that would be a much larger project (e.g. language server).
- **Method completion:** Completion for `obj.method()` is out of scope for this plan; only top-level names (e.g. after `from openscad import *`) are targeted.

---

## Effort Estimate

| Task | Effort |
|------|--------|
| PythonApi class (data + updateAutoCompletionList + callTips) | 0.5–1 day |
| Integrate into ScintillaEditor (create API, setAPIs, ENABLE_PYTHON) | 0.25 day |
| Testing (manual + regression) | 0.25–0.5 day |
| Optional: skip completion in Python strings | 0.25 day |
| **Total** | **~1–2.5 days** |

---

## References (in tree)

- SCAD completion / call tips: `src/gui/ScadApi.cc`, `src/gui/ScadApi.h`
- Lexer/API wiring: `src/gui/ScintillaEditor.cc` (e.g. `setLexer`, `setupAutoComplete`), `ScintillaEditor::onLanguageChanged`
- Builtins and call-tip data: `src/core/Builtins.h`, `Builtins::keywordList`; registration in `src/core/primitives.cc`, `builtin_functions.cc`, etc.
- Python exports: `src/python/pyfunctions.cc` (`PyOpenSCADFunctions`), `doc/openscad.pyi`
- QScintilla: `QsciAbstractAPIs`, `QsciLexer::setAPIs`, `QsciScintilla::setAutoCompletionSource`, `setCallTipsVisible` / `setCallTipsStyle`
