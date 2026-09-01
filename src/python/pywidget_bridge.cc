// Python.h MUSS die allererste Zeile dieser Datei sein (CPython-Vorgabe).
// Diese Datei includet im Gegenzug NIE einen Qt-Header und wird NIE
// durch moc geschickt - so kollidieren Python.h und Q_OBJECT nie in
// derselben Übersetzungseinheit.
#include <Python.h>

#include "src/python/pywidget_bridge.h"

#include <map>
#include <string>

// Diese Map wird von add_parameter_widget() befüllt, siehe
// src/python/pywidget_registry.cc. extern deklariert, damit es nur
// eine Quelle der Wahrheit gibt.
extern std::map<std::string, PyObject *> customizer_widget_factories;

namespace {

PyObject *json_to_pyobject(const json& j)
{
  if (j.is_boolean()) return PyBool_FromLong(j.get<bool>());
  if (j.is_number_integer()) return PyLong_FromLongLong(j.get<long long>());
  if (j.is_number()) return PyFloat_FromDouble(j.get<double>());
  if (j.is_string()) return PyUnicode_FromString(j.get<std::string>().c_str());
  if (j.is_array()) {
    PyObject *list = PyList_New(j.size());
    Py_ssize_t i = 0;
    for (const auto& item : j) {
      PyList_SET_ITEM(list, i++, json_to_pyobject(item));
    }
    return list;
  }
  Py_RETURN_NONE;
}

json pyobject_to_json(PyObject *obj)
{
  if (!obj || obj == Py_None) return json(nullptr);
  if (PyBool_Check(obj)) return json(obj == Py_True);
  if (PyLong_Check(obj)) return json(PyLong_AsLongLong(obj));
  if (PyFloat_Check(obj)) return json(PyFloat_AsDouble(obj));
  if (PyUnicode_Check(obj)) return json(std::string(PyUnicode_AsUTF8(obj)));
  if (PyList_Check(obj) || PyTuple_Check(obj)) {
    json arr = json::array();
    Py_ssize_t n = PySequence_Size(obj);
    for (Py_ssize_t i = 0; i < n; i++) {
      PyObject *item = PySequence_GetItem(obj, i);  // neue Referenz
      arr.push_back(pyobject_to_json(item));
      Py_DECREF(item);
    }
    return arr;
  }
  // Fallback: str(obj) - besser eine grobe Darstellung als ein Absturz.
  PyObject *str = PyObject_Str(obj);
  json result = str ? json(std::string(PyUnicode_AsUTF8(str))) : json(nullptr);
  Py_XDECREF(str);
  return result;
}

}  // namespace

PyWidgetHandle pywidget_create(const std::string& typeName, const std::string& paramName,
                                const json& defaultValue)
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  PyWidgetHandle result{nullptr, nullptr};

  auto it = customizer_widget_factories.find(typeName);
  if (it != customizer_widget_factories.end()) {
    PyObject *pyname = PyUnicode_FromString(paramName.c_str());
    PyObject *pydefault = json_to_pyobject(defaultValue);

    PyObject *obj = PyObject_CallFunctionObjArgs(it->second, pyname, pydefault, NULL);
    Py_XDECREF(pyname);
    Py_XDECREF(pydefault);

    if (obj) {
      // Gegenstück zu dem sip.wrapinstance(), das mainwindow_ptr() schon
      // in die andere Richtung benutzt (C++-Objekt -> PyQt6-Objekt).
      PyObject *sipmod = PyImport_ImportModule("sip");
      if (sipmod) {
        PyObject *addr = PyObject_CallMethod(sipmod, "unwrapinstance", "O", obj);
        if (addr) {
          result.qwidget_ptr = PyLong_AsVoidPtr(addr);
          Py_DECREF(addr);
        }
        Py_DECREF(sipmod);
      }
      if (result.qwidget_ptr) {
        result.py_object_handle = obj;  // Ownership geht an den Aufrufer über
      } else {
        Py_DECREF(obj);  // kein gültiges Widget -> Referenz nicht behalten
      }
    } else {
      PyErr_Print();  // Factory hat eine Exception geworfen - sichtbar machen statt zu verschlucken
    }
  }

  PyGILState_Release(gstate);
  return result;
}

json pywidget_get_value(void *pyObjectHandle)
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  json result = json(nullptr);

  PyObject *obj = reinterpret_cast<PyObject *>(pyObjectHandle);
  PyObject *val = PyObject_CallMethod(obj, "get_value", nullptr);
  if (val) {
    result = pyobject_to_json(val);
    Py_DECREF(val);
  } else {
    PyErr_Print();
  }

  PyGILState_Release(gstate);
  return result;
}

void pywidget_release(void *pyObjectHandle)
{
  if (!pyObjectHandle) return;
  PyGILState_STATE gstate = PyGILState_Ensure();
  Py_DECREF(reinterpret_cast<PyObject *>(pyObjectHandle));
  PyGILState_Release(gstate);
}
