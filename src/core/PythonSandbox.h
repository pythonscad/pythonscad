#pragma once

#include <string>

#ifdef ENABLE_PYTHON

struct PythonSandboxResult {
  bool ok = false;
  std::string csg;
  std::string error;
};

PythonSandboxResult evaluatePythonSandboxToCsg(const std::string& code, const std::string& scriptpath);

#endif  // ENABLE_PYTHON
