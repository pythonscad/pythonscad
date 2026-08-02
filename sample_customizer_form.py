from pythonscad import *

from PyQt6 import QtWidgets, QtCore

class ColorPicker(QtWidgets.QPushButton):
    valueChanged = QtCore.pyqtSignal()

    def __init__(self, name, default):
        super().__init__(default)
        self._value = default
        self.clicked.connect(self._pick)

    def _pick(self):
        c = QtWidgets.QColorDialog.getColor()
        if c.isValid():
            self._value = c.name()
            self.setText(self._value)
            self.valueChanged.emit()

    def get_value(self):
        return self._value

def color_picker_factory(name, default):
    return ColorPicker(name, default)


def setup():
    add_parameter_widget("color_picker", color_picker_factory)
