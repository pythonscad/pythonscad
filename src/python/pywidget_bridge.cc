// pywidget_bridge.cc
#include <Python.h>  // MUSS die allererste Zeile sein
#include "pywidget_bridge.h"

extern std::map<std::string, PyObject *> customizer_widget_factories;

PyWidgetHandle pywidget_create(const std::string& typeName, const std::string& paramName,
                               PyObject *defaultVal)
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  PyWidgetHandle result{nullptr, nullptr};

  auto it = customizer_widget_factories.find(typeName);
  if (it != customizer_widget_factories.end()) {
    PyObject *pyname = PyUnicode_FromString(paramName.c_str());
    PyObject *obj = PyObject_CallFunctionObjArgs(it->second, pyname, defaultVal, NULL);
    Py_XDECREF(pyname);
    if (obj) {
      PyObject *sipmod = PyImport_ImportModule("sip");
      PyObject *addr = PyObject_CallMethod(sipmod, "unwrapinstance", "O", obj);
      result.qwidget_ptr = PyLong_AsVoidPtr(addr);
      result.py_object = obj;  // Ownership geht an den Aufrufer über
      Py_XDECREF(addr);
      Py_XDECREF(sipmod);
    }
  }
  PyGILState_Release(gstate);
  return result;
}

PyObject *pywidget_get_value(PyObject *obj)
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  PyObject *val = PyObject_CallMethod(obj, "get_value", nullptr);
  PyGILState_Release(gstate);
  return val;  // neue Referenz, Aufrufer muss DECREF (mit GIL!) machen
}

void pywidget_release(PyObject *obj)
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  Py_XDECREF(obj);
  PyGILState_Release(gstate);
}
