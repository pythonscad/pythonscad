#include "gui/parameter/ParameterPyQtWidget.h"
#include "core/customizer/ParameterObject.h"
#include "pywidget_bridge.h"
#include <QVBoxLayout>
#include <QTimer>
#include <QMetaMethod>


ParameterPyQtWidget::ParameterPyQtWidget(QWidget *parent, CustomParameter *parameter,
                                         DescriptionStyle descriptionStyle)
  : ParameterVirtualWidget(parent, parameter), parameter(parameter)
{
  PyWidgetHandle h =
    pywidget_create(parameter->customTypeName, parameter->name(), parameter->value);
  pyWidget = reinterpret_cast<QWidget *>(h.qwidget_ptr);
  pyWidgetObj = (PyObject *) h.py_object_handle;

  auto *layout = new QHBoxLayout(this);
  layout->setContentsMargins(0, 0, 0, 0);
  if (pyWidget) {
    auto *label = new QLabel(QString::fromStdString(parameter->name()), this);
      label->setToolTip(QString::fromStdString(parameter->description()));
    layout->addWidget(label);	  
    layout->addWidget(pyWidget);
    int idx = pyWidget->metaObject()->indexOfSignal("valueChanged()");
    if (idx >= 0) {
      connect(pyWidget, pyWidget->metaObject()->method(idx), this,
              this->metaObject()->method(this->metaObject()->indexOfSlot("onPyWidgetValueChanged()")));
    }
  }
}

void ParameterPyQtWidget::onPyWidgetValueChanged()
{
  if (pyWidgetObj) {
    parameter->value = pywidget_get_value(pyWidgetObj);  // Widget -> Parameter, HIER richtig
  }

  if (pyWidgetObj) {
    parameter->value = pywidget_get_value(pyWidgetObj);
  }
  // Nicht synchron emittieren: wir laufen hier noch innerhalb des sip-
  // Signal-Dispatch-Callstacks (GIL über sip's eigenes PyGILState_Ensure
  // gehalten, nicht über unser tstate-System). Ein synchrones
  // changed(true) würde MainWindow::compile() -> python_lock() re-entrant
  // im selben Stack auslösen und mit der eigenen tstate-Buchhaltung
  // kollidieren. Auf die nächste Event-Loop-Iteration verschieben, damit
  // der sip-Kontext sauber beendet ist, bevor compile() läuft.
  QTimer::singleShot(0, this, [this]() { emit changed(true); });
}
void ParameterPyQtWidget::setValue()
{
  //if (!pyWidgetObj) return;
  //parameter->value = pywidget_get_value(pyWidgetObj);
}

ParameterPyQtWidget::~ParameterPyQtWidget()
{
  pywidget_release(pyWidgetObj);
}
