#pragma once
#include "gui/parameter/ParameterVirtualWidget.h"

// Opaque forward declaration statt #include <Python.h> - moc-sicher
#ifndef Py_PYTHON_H
struct _object;
typedef _object PyObject;
#endif

class CustomParameter;

class ParameterPyQtWidget : public ParameterVirtualWidget
{
  Q_OBJECT
public:
  ParameterPyQtWidget(QWidget *parent, CustomParameter *parameter, DescriptionStyle descriptionStyle);
  ~ParameterPyQtWidget() override;
  void setValue() override;

private slots:
  void onPyWidgetValueChanged();

private:
  CustomParameter *parameter;
  QWidget *pyWidget = nullptr;
  PyObject *pyWidgetObj = nullptr;  // nur als Zeiger bekannt, nie dereferenziert in dieser TU
};
