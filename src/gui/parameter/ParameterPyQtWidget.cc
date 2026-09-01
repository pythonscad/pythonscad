#include "gui/parameter/ParameterPyQtWidget.h"
#include "core/customizer/ParameterObject.h"
#include "pywidget_bridge.h"
#include <QVBoxLayout>
#include <QMetaMethod>


ParameterPyQtWidget::ParameterPyQtWidget(QWidget *parent, CustomParameter *parameter,
                                         DescriptionStyle descriptionStyle)
  : ParameterVirtualWidget(parent, parameter), parameter(parameter)
{
  PyWidgetHandle h =
    pywidget_create(parameter->customTypeName, parameter->name(), parameter->value);
  pyWidget = reinterpret_cast<QWidget *>(h.qwidget_ptr);
  pyWidgetObj = (PyObject *) h.py_object_handle;

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

void ParameterPyQtWidget::onPyWidgetValueChanged()
{
  // "true" = immediate, analog zu ParameterCheckBox/ParameterComboBox,
  // die ebenfalls diskrete statt kontinuierliche Änderungen liefern.
  emit changed(true);
}
void ParameterPyQtWidget::setValue()
{
  if (!pyWidgetObj) return;
  parameter->value = pywidget_get_value(pyWidgetObj);
}

ParameterPyQtWidget::~ParameterPyQtWidget()
{
  pywidget_release(pyWidgetObj);
}
