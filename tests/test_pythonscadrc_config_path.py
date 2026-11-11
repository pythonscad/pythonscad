"""Host-side regression checks for the platform-specific .pythonscadrc path."""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


class PythonScadRcConfigPathTests(unittest.TestCase):
    def test_gui_uses_python_config_path_without_fallback(self):
        source = (ROOT / "src/gui/MainWindow.cc").read_text()
        lookup = source.split("std::string MainWindow::loadInitFile", 1)[1].split(
            "void MainWindow::customSetup", 1
        )[0]
        self.assertIn('lookup_file(".pythonscadrc", PlatformUtils::pythonUserConfigPath(), "")', lookup)
        self.assertNotIn("userConfigPath()", lookup)

    def test_each_platform_has_a_pythonscad_specific_directory(self):
        posix = (ROOT / "src/platform/PlatformUtils-posix.cc").read_text()
        macos = (ROOT / "src/platform/PlatformUtils-mac.mm").read_text()
        windows = (ROOT / "src/platform/PlatformUtils-win.cc").read_text()

        self.assertIn('getXdgConfigDir() / "PythonSCAD"', posix)
        self.assertIn('std::string("/PythonSCAD")', macos)
        self.assertIn("CSIDL_LOCAL_APPDATA", windows)
        self.assertIn('std::string("/PythonSCAD")', windows)

    def test_form_model_preserves_the_startup_badge_class(self):
        model = (ROOT / "doc/examples/pyqt6-form/badge.py").read_text()
        self.assertIn('if "badge" not in globals():', model)

    def test_form_badge_returns_geometry(self):
        startup = (ROOT / "doc/examples/pyqt6-form/pythonscadrc.py").read_text()
        constructor = startup.split("def __new__", 1)[1].split("\n\n", 1)[0]
        self.assertIn("return disc +", constructor)
        self.assertNotIn("show(", constructor)


if __name__ == "__main__":
    unittest.main()
