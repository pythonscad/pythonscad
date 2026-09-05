#pragma once

// Diese Datei darf NIE Qt-Header includen (kein QWidget, kein QObject).
// Sie ist die einzige erlaubte Schnittstelle zwischen der Qt-Welt
// (ParameterPyQtWidget.h/.cc, läuft durch moc) und der Python-Welt
// (pywidget_bridge.cc, includet <Python.h> als erste Zeile).
//
// Der QWidget*-Pointer wird bewusst nur als void* durchgereicht - erst
// ParameterPyQtWidget.cc (Qt-Seite, kein Python.h) macht daraus wieder
// ein QWidget*. Genauso wird das PyObject* nur als opaker void*-Handle
// durchgereicht - erst pywidget_bridge.cc (Python-Seite, kein Qt)
// interpretiert ihn wieder als PyObject*.

#include <string>
#include "json/json.hpp"
using json = nlohmann::json;

struct PyWidgetHandle {
  void *qwidget_ptr;      // -> reinterpret_cast<QWidget*> auf Qt-Seite
  void *py_object_handle; // -> reinterpret_cast<PyObject*> auf Python-Seite
};

// Ruft die per add_parameter_widget(typeName, factory) registrierte
// Python-Factory auf: factory(name, default) -> QWidget (PyQt6).
// Bei Fehler/unbekanntem Typ ist das Ergebnis {nullptr, nullptr}.
PyWidgetHandle pywidget_create(const std::string& typeName, const std::string& paramName,
                                const json& defaultValue);

// Ruft obj.get_value() auf dem von der Factory gelieferten Objekt auf.
json pywidget_get_value(void *pyObjectHandle);

void pywidget_set_value(void *pyObjectHandle, const json& value);

// Gibt die gehaltene Python-Referenz frei (Py_DECREF unter GIL).
void pywidget_release(void *pyObjectHandle);
