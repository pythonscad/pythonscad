// pywidget_bridge.h
#pragma once
#include <string>

#ifndef Py_PYTHON_H
struct _object;
typedef _object PyObject;
#endif

// Liefert rohen QWidget*-Pointer (als void*) + die PyObject*-Referenz zurück,
// damit der Aufrufer (Qt-Seite) beides ohne Python.h weiterreichen kann.
struct PyWidgetHandle {
  void *qwidget_ptr;    // -> reinterpret_cast<QWidget*> auf Qt-Seite
  PyObject *py_object;  // starke Referenz, muss später via pywidget_release freigegeben werden
};

PyWidgetHandle pywidget_create(const std::string& typeName, const std::string& paramName,
                               PyObject *defaultVal);
PyObject *pywidget_get_value(PyObject *obj);  // ruft obj.get_value() auf
void pywidget_release(PyObject *obj);
