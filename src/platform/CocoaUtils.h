#pragma once

#include <string>

class CocoaUtils
{
public:
  static void endApplication();
  static void nslog(const std::string& str, void *userdata);
  // Force the macOS appearance (dark/light) for the whole app so that native
  // Cocoa widgets match PythonSCAD's chosen theme even when it differs from
  // the system-wide setting.
  static void setAppearance(bool dark);
};
