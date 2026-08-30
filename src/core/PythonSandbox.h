#pragma once

#include <string>
#include <vector>

#ifdef ENABLE_PYTHON

struct PythonSandboxOutputFile {
  std::string relativePath;
  std::string hostPath;
  unsigned long long size = 0;
};

struct PythonSandboxResult {
  bool ok = false;
  std::string csg;
  std::string error;
  std::string tempDir;
  std::string outputRoot;
  std::vector<PythonSandboxOutputFile> outputFiles;
};

PythonSandboxResult evaluatePythonSandboxToCsg(const std::string& code, const std::string& scriptpath);
void cleanupPythonSandboxResult(const PythonSandboxResult& result);

// Windows reserved device names (CON, NUL, ...) after stripping trailing
// spaces and dots that the Win32 path normalizer would otherwise remove.
bool isReservedWindowsDeviceNameComponent(const std::string& component);

#endif  // ENABLE_PYTHON
