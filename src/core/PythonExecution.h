#pragma once

#include <cstdint>
#include <string>

#include "core/Settings.h"

#ifdef ENABLE_PYTHON

enum class PythonExecutionMode : std::uint8_t { Sandboxed, Native };

inline constexpr const char *PYTHON_EXECUTION_SANDBOXED = "sandboxed";
inline constexpr const char *PYTHON_EXECUTION_NATIVE = "native";

inline PythonExecutionMode pythonExecutionModeFromString(const std::string& mode)
{
  return mode == PYTHON_EXECUTION_NATIVE ? PythonExecutionMode::Native : PythonExecutionMode::Sandboxed;
}

inline const char *pythonExecutionModeToString(PythonExecutionMode mode)
{
  return mode == PythonExecutionMode::Native ? PYTHON_EXECUTION_NATIVE : PYTHON_EXECUTION_SANDBOXED;
}

inline PythonExecutionMode defaultPythonExecutionMode()
{
#ifdef __EMSCRIPTEN__
  // The Emscripten build's CPython interpreter is already inside the WASM sandbox.
  return PythonExecutionMode::Native;
#else
  return pythonExecutionModeFromString(Settings::SettingsPython::pythonExecutionMode.value());
#endif
}

inline bool pythonExecutionModeIsNative(PythonExecutionMode mode)
{
  return mode == PythonExecutionMode::Native;
}

#endif  // ENABLE_PYTHON
