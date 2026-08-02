#include "gui/parameter/ParameterPyQtWidget.h"
#include "pywidget_bridge.h"
#include <QVBoxLayout>
#include <QMetaMethod>

ParameterPyQtWidget::ParameterPyQtWidget(QWidget *parent, CustomParameter *parameter,
                                         DescriptionStyle descriptionStyle)
  : ParameterVirtualWidget(parent, parameter), parameter(parameter)
{
  PyWidgetHandle h =
    pywidget_create(parameter->customTypeName, parameter->name(), /* default als PyObject* */);
  pyWidget = reinterpret_cast<QWidget *>(h.qwidget_ptr);
  pyWidgetObj = h.py_object;

  auto *layout = new QVBoxLayout(this);
  layout->setContentsMargins(0, 0, 0, 0);
  if (pyWidget) {
    layout->addWidget(pyWidget);
    int idx = pyWidget->metaObject()->indexOfSignal("valueChanged()");
    if (idx >= 0) {
      connect(pyWidget, pyWidget->metaObject()->method(idx), this,
              this->metaObject()->method(this->metaObject()->indexOfSlot("onPyWidgetValueChanged()")));
    }
  }
}
