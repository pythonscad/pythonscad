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

#endif  // ENABLE_PYTHON
