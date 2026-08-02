#pragma once

#include "gui/parameter/ParameterVirtualWidget.h"
#include <Python.h>

class CustomParameter;  // Erweiterung von ParameterObject, siehe unten

class ParameterPyQtWidget : public ParameterVirtualWidget
{
  Q_OBJECT

public:
  ParameterPyQtWidget(QWidget *parent, CustomParameter *parameter, DescriptionStyle descriptionStyle);
  ~ParameterPyQtWidget() override;

  void setValue() override;
  // valueApplied() -> Default-Implementierung (no-op) reicht meistens

private slots:
  void onPyWidgetValueChanged();

private:
  CustomParameter *parameter;
  QWidget *pyWidget = nullptr;      // via sip.unwrapinstance() gewonnen
  PyObject *pyWidgetObj = nullptr;  // starke Referenz auf das PyQt6-Objekt (für get_value())
};
