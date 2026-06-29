#include "core/PythonSandbox.h"

#ifdef ENABLE_PYTHON

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>

#include "platform/PlatformUtils.h"

namespace fs = std::filesystem;

namespace {

std::string getenvString(const char *name)
{
  const char *value = std::getenv(name);
  return value ? std::string(value) : std::string();
}

std::string commandQuote(const fs::path& path)
{
#ifdef _WIN32
  std::string value = path.string();
  std::string quoted = "\"";
  for (const char ch : value) {
    if (ch == '"') quoted += "\\\"";
    else quoted += ch;
  }
  quoted += "\"";
  return quoted;
#else
  std::string value = path.string();
  std::string quoted = "'";
  for (const char ch : value) {
    if (ch == '\'') quoted += "'\\''";
    else quoted += ch;
  }
  quoted += "'";
  return quoted;
#endif
}

bool hasWasmBundle(const fs::path& dir)
{
  return fs::is_regular_file(dir / "pythonscad.js") && fs::is_regular_file(dir / "pythonscad.wasm") &&
         fs::is_regular_file(dir / "pythonscad.data");
}

fs::path findWasmBundleDir()
{
  const auto fromEnv = getenvString("PYTHONSCAD_WASM_DIR");
  if (!fromEnv.empty() && hasWasmBundle(fromEnv)) return fs::path(fromEnv);

  const fs::path cwd = fs::current_path();
  if (hasWasmBundle(cwd / "build-wasm-web")) return cwd / "build-wasm-web";

  const fs::path appPath = PlatformUtils::applicationPath();
  if (hasWasmBundle(appPath / "wasm")) return appPath / "wasm";
  if (hasWasmBundle(appPath)) return appPath;

  return {};
}

fs::path findSandboxRunner()
{
  const auto fromEnv = getenvString("PYTHONSCAD_SANDBOX_RUNNER");
  if (!fromEnv.empty() && fs::is_regular_file(fromEnv)) return fs::path(fromEnv);

  const fs::path cwdRunner = fs::current_path() / "scripts" / "python-sandbox-runner.mjs";
  if (fs::is_regular_file(cwdRunner)) return cwdRunner;

  const fs::path appRunner =
    fs::path(PlatformUtils::applicationPath()) / "scripts" / "python-sandbox-runner.mjs";
  if (fs::is_regular_file(appRunner)) return appRunner;

  return {};
}

fs::path makeTempDir()
{
  const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
  fs::path dir = fs::temp_directory_path() / ("pythonscad-sandbox-" + std::to_string(stamp));
  fs::create_directories(dir);
  return dir;
}

bool readFile(const fs::path& path, std::string& output)
{
  std::ifstream stream(path, std::ios::binary);
  if (!stream.is_open()) return false;
  output.assign(std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>());
  return true;
}

bool writeFile(const fs::path& path, const std::string& content)
{
  std::ofstream stream(path, std::ios::binary);
  if (!stream.is_open()) return false;
  stream.write(content.data(), static_cast<std::streamsize>(content.size()));
  return stream.good();
}

}  // namespace

PythonSandboxResult evaluatePythonSandboxToCsg(const std::string& code, const std::string& scriptpath)
{
  PythonSandboxResult result;

  const fs::path wasmDir = findWasmBundleDir();
  if (wasmDir.empty()) {
    result.error =
      "Sandboxed Python requires a WebAssembly bundle. Set PYTHONSCAD_WASM_DIR to a "
      "directory containing pythonscad.js, pythonscad.wasm, and pythonscad.data.";
    return result;
  }

  const fs::path runner = findSandboxRunner();
  if (runner.empty()) {
    result.error =
      "Sandboxed Python runner was not found. Set PYTHONSCAD_SANDBOX_RUNNER to "
      "scripts/python-sandbox-runner.mjs.";
    return result;
  }

  const std::string node =
    getenvString("PYTHONSCAD_NODE").empty() ? "node" : getenvString("PYTHONSCAD_NODE");
  const fs::path tempDir = makeTempDir();
  const fs::path inputFile =
    tempDir /
    (fs::path(scriptpath).filename().empty() ? fs::path("input.py") : fs::path(scriptpath).filename());
  const fs::path outputFile = tempDir / "output.csg";

  if (!writeFile(inputFile, code)) {
    result.error = "Could not write sandbox input file.";
    fs::remove_all(tempDir);
    return result;
  }

  std::ostringstream command;
  command << commandQuote(node) << " " << commandQuote(runner) << " --wasm-dir " << commandQuote(wasmDir)
          << " --input " << commandQuote(inputFile) << " --output " << commandQuote(outputFile);

  const int exitCode = std::system(command.str().c_str());
  if (exitCode != 0) {
    result.error = "Sandboxed Python process failed.";
    fs::remove_all(tempDir);
    return result;
  }

  if (!readFile(outputFile, result.csg)) {
    result.error = "Sandboxed Python did not produce CSG output.";
    fs::remove_all(tempDir);
    return result;
  }

  fs::remove_all(tempDir);
  result.ok = true;
  return result;
}

#endif  // ENABLE_PYTHON
