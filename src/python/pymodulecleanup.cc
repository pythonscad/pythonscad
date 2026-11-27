#include "pymodulecleanup.h"

#include <algorithm>
#include <cctype>
#include <cstring>
#include <map>
#include <utility>

namespace {

bool to_utf8(PyObject *obj, std::string& out)
{
  if (obj == nullptr || !PyUnicode_Check(obj)) return false;
  Py_ssize_t len = 0;
  const char *s = PyUnicode_AsUTF8AndSize(obj, &len);
  if (s == nullptr) {
    PyErr_Clear();
    return false;
  }
  out.assign(s, static_cast<size_t>(len));
  return true;
}

std::string top_level_name(const std::string& name)
{
  const auto dot = name.find('.');
  return dot == std::string::npos ? name : name.substr(0, dot);
}

bool ends_with_ci(const std::string& s, const char *suffix)
{
  const size_t n = std::strlen(suffix);
  if (s.size() < n) return false;
  return std::equal(s.end() - static_cast<std::ptrdiff_t>(n), s.end(), suffix,
                    [](char a, char b) { return std::tolower(static_cast<unsigned char>(a)) == b; });
}

bool starts_with_ci(const std::string& s, const std::string& prefix)
{
  if (s.size() < prefix.size()) return false;
  return std::equal(prefix.begin(), prefix.end(), s.begin(), [](char a, char b) {
    return std::tolower(static_cast<unsigned char>(a)) == std::tolower(static_cast<unsigned char>(b));
  });
}

bool is_stdlib_name(PyObject *stdlibNames, const std::string& top)
{
  if (stdlibNames == nullptr) return false;  // Python < 3.10: fall back to the path heuristics
  PyObject *topObj = PyUnicode_FromStringAndSize(top.data(), static_cast<Py_ssize_t>(top.size()));
  if (topObj == nullptr) {
    PyErr_Clear();
    return false;
  }
  const int found = PySequence_Contains(stdlibNames, topObj);
  Py_DECREF(topObj);
  if (found < 0) {
    PyErr_Clear();
    return false;
  }
  return found == 1;
}

bool module_file_path(PyObject *mod, std::string& path)
{
  PyObject *file = PyObject_GetAttrString(mod, "__file__");
  if (file == nullptr) {
    PyErr_Clear();
    return false;
  }
  const bool ok = to_utf8(file, path);
  Py_DECREF(file);
  return ok;
}

bool path_under_sys_prefix(const std::string& path, const char *sysAttr)
{
  PyObject *prefixObj = PySys_GetObject(sysAttr);  // borrowed
  std::string prefix;
  if (prefixObj == nullptr || !to_utf8(prefixObj, prefix) || prefix.empty()) return false;
  // Trim a trailing separator so "C:\\Python" and "C:\\Python\\" both match.
  while (!prefix.empty() && (prefix.back() == '/' || prefix.back() == '\\')) prefix.pop_back();
  if (!starts_with_ci(path, prefix)) return false;
  if (path.size() == prefix.size()) return true;
  const char next = path[prefix.size()];
  return next == '/' || next == '\\';
}

// A name in sys.stdlib_module_names is not enough: a user file next to the
// model (mailbox.py, profile.py, ...) shares that name but must still be
// dropped so edits are picked up on re-run. Only protect the group when the
// loaded module actually lives under the interpreter's install prefix.
bool is_stdlib_module_by_origin(PyObject *mod, PyObject *stdlibNames, const std::string& top)
{
  if (!is_stdlib_name(stdlibNames, top)) return false;
  std::string path;
  if (!module_file_path(mod, path)) return false;
  return path_under_sys_prefix(path, "base_prefix") || path_under_sys_prefix(path, "prefix");
}

bool is_extension_module(PyObject *mod)
{
  std::string path;
  if (!module_file_path(mod, path)) return false;
  return ends_with_ci(path, ".pyd") || ends_with_ci(path, ".so") || ends_with_ci(path, ".dylib");
}

// The historic repr()-based filters. Kept so that behaviour on Python < 3.10
// and for third-party packages stays as before.
bool repr_marks_protected(PyObject *mod)
{
  PyObject *repr = PyObject_Repr(mod);
  if (repr == nullptr) {
    PyErr_Clear();
    return true;  // cannot judge it -> leave it alone
  }
  std::string r;
  const bool ok = to_utf8(repr, r);
  Py_DECREF(repr);
  if (!ok) return true;
  static const char *const markers[] = {"(frozen)", "(built-in)",    "/encodings/", "\\encodings\\",
                                        "_frozen_", "site-packages", "usr/lib"};
  for (const char *m : markers) {
    if (r.find(m) != std::string::npos) return true;
  }
  return false;
}

}  // namespace

std::vector<std::string> python_collect_stale_modules(PyObject *modules,
                                                      const std::vector<std::string>& forceDrop)
{
  std::vector<std::string> result;
  if (modules == nullptr || !PyDict_Check(modules)) return result;

  PyObject *stdlibNames = PySys_GetObject("stdlib_module_names");  // borrowed, may be NULL

  PyObject *keys = PyDict_Keys(modules);
  if (keys == nullptr) {
    PyErr_Clear();
    return result;
  }

  std::vector<std::pair<std::string, std::string>> entries;  // (name, top-level)
  std::map<std::string, bool> keepGroup;                     // top-level -> keep?

  const Py_ssize_t n = PyList_GET_SIZE(keys);
  for (Py_ssize_t i = 0; i < n; ++i) {
    PyObject *key = PyList_GET_ITEM(keys, i);
    std::string name;
    if (!to_utf8(key, name)) continue;
    const std::string top = top_level_name(name);
    entries.emplace_back(name, top);

    bool& keep = keepGroup.emplace(top, false).first->second;
    if (keep) continue;
    if (top == "__main__") {
      keep = true;
      continue;
    }
    PyObject *value = PyDict_GetItem(modules, key);  // borrowed
    // Non-module values (e.g. None placeholders) neither protect nor block
    // their group; they are dropped together with it.
    if (value == nullptr || !PyModule_Check(value)) continue;
    // Do not protect on the stdlib name alone: a user module that shadows a
    // stdlib name (mailbox.py, profile.py, ...) must still be reloadable.
    if (is_extension_module(value) || repr_marks_protected(value) ||
        is_stdlib_module_by_origin(value, stdlibNames, top)) {
      keep = true;
    }
  }
  Py_DECREF(keys);

  for (const auto& f : forceDrop) keepGroup[f] = false;

  for (const auto& [name, top] : entries) {
    if (!keepGroup[top]) result.push_back(name);
  }
  return result;
}
