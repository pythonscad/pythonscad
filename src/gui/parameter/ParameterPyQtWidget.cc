#include "gui/parameter/ParameterPyQtWidget.h"
#include "core/customizer/ParameterObject.h"
#include <QVBoxLayout>
#include <QMetaMethod>

extern std::map<std::string, PyObject *> customizer_widget_factories;  // aus pyfunctions.cc

ParameterPyQtWidget::ParameterPyQtWidget(QWidget *parent, CustomParameter *parameter,
                                         DescriptionStyle descriptionStyle)
  : ParameterVirtualWidget(parent, parameter), parameter(parameter)
{
  PyGILState_STATE gstate = PyGILState_Ensure();

  auto it = customizer_widget_factories.find(parameter->customTypeName);
  if (it != customizer_widget_factories.end()) {
    PyObject *pyname = PyUnicode_FromString(parameter->name().c_str());
    PyObject *pydefault = /* passendes PyObject aus parameter->currentValue() bauen */;
    pyWidgetObj = PyObject_CallFunctionObjArgs(it->second, pyname, pydefault, NULL);
    Py_XDECREF(pyname);
    Py_XDECREF(pydefault);

    if (pyWidgetObj) {
      PyObject *sipmod = PyImport_ImportModule("sip");
      PyObject *addr = PyObject_CallMethod(sipmod, "unwrapinstance", "O", pyWidgetObj);
      pyWidget = reinterpret_cast<QWidget *>(PyLong_AsVoidPtr(addr));
      Py_XDECREF(addr);
      Py_XDECREF(sipmod);
    }
  }

  PyGILState_Release(gstate);

  auto *layout = new QVBoxLayout(this);
  layout->setContentsMargins(0, 0, 0, 0);
  if (pyWidget) {
    layout->addWidget(pyWidget);

    // Generische Anbindung: das PyQt6-Widget muss ein Signal namens
    // "valueChanged()" deklarieren (pyqtSignal() in Python) - analog
    // dazu, dass jede bestehende Parameter*-Klasse ihr eigenes
    // Qt-Signal auf changed(bool) mappt.
    int idx = pyWidget->metaObject()->indexOfSignal("valueChanged()");
    if (idx >= 0) {
      QMetaMethod signal = pyWidget->metaObject()->method(idx);
      QMetaMethod slot =
        this->metaObject()->method(this->metaObject()->indexOfSlot("onPyWidgetValueChanged()"));
      connect(pyWidget, signal, this, slot);
    }
  }
}

ParameterPyQtWidget::~ParameterPyQtWidget()
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  Py_XDECREF(pyWidgetObj);
  PyGILState_Release(gstate);
}

void ParameterPyQtWidget::onPyWidgetValueChanged()
{
  emit changed(true);  // wie ParameterCheckBox/ParameterComboBox: sofort
}

void ParameterPyQtWidget::setValue()
{
  if (!pyWidgetObj) return;
  PyGILState_STATE gstate = PyGILState_Ensure();
  PyObject *val = PyObject_CallMethod(pyWidgetObj, "get_value", nullptr);
  if (val) {
    parameter->importValueFromPython(val);  // Pendant zu importValue() bei anderen Parametertypen
  }
  Py_XDECREF(val);
  PyGILState_Release(gstate);
}
