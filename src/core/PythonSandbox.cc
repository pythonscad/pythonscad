#include "core/PythonSandbox.h"

#ifdef ENABLE_PYTHON

#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

#include "platform/PlatformUtils.h"

namespace fs = std::filesystem;

namespace {

std::string getenvString(const char *name)
{
  const char *value = std::getenv(name);
  return value ? std::string(value) : std::string();
}

std::string windowsCommandLineQuote(const std::string& value)
{
  std::string quoted = "\"";
  size_t backslashes = 0;
  for (const char ch : value) {
    if (ch == '\\') {
      ++backslashes;
    } else if (ch == '"') {
      quoted.append(backslashes * 2 + 1, '\\');
      quoted += ch;
      backslashes = 0;
    } else {
      quoted.append(backslashes, '\\');
      quoted += ch;
      backslashes = 0;
    }
  }
  quoted.append(backslashes * 2, '\\');
  quoted += "\"";
  return quoted;
}

bool hasWasmBundle(const fs::path& dir)
{
  return fs::is_regular_file(dir / "pythonscad.js") && fs::is_regular_file(dir / "pythonscad.wasm") &&
         fs::is_regular_file(dir / "pythonscad.data");
}

std::string processStatusMessage(int status)
{
#ifdef _WIN32
  return "exit code " + std::to_string(status);
#else
  if (WIFEXITED(status)) return "exit code " + std::to_string(WEXITSTATUS(status));
  if (WIFSIGNALED(status)) return "signal " + std::to_string(WTERMSIG(status));
  return "status " + std::to_string(status);
#endif
}

struct ProcessResult {
  bool ok = false;
  std::string error;
};

ProcessResult runProcess(const std::vector<std::string>& args)
{
  ProcessResult result;
  if (args.empty()) {
    result.error = "missing process arguments";
    return result;
  }

#ifdef _WIN32
  std::string commandLine;
  for (const auto& arg : args) {
    if (!commandLine.empty()) commandLine += " ";
    commandLine += windowsCommandLineQuote(arg);
  }

  STARTUPINFOA startupInfo{};
  startupInfo.cb = sizeof(startupInfo);
  PROCESS_INFORMATION processInfo{};
  std::vector<char> mutableCommandLine(commandLine.begin(), commandLine.end());
  mutableCommandLine.push_back('\0');
  if (!CreateProcessA(nullptr, mutableCommandLine.data(), nullptr, nullptr, FALSE, 0, nullptr, nullptr,
                      &startupInfo, &processInfo)) {
    result.error = "could not start process";
    return result;
  }

  WaitForSingleObject(processInfo.hProcess, INFINITE);
  DWORD exitCode = 1;
  GetExitCodeProcess(processInfo.hProcess, &exitCode);
  CloseHandle(processInfo.hThread);
  CloseHandle(processInfo.hProcess);
  if (exitCode != 0) {
    result.error = processStatusMessage(static_cast<int>(exitCode));
    return result;
  }
  result.ok = true;
  return result;
#else
  const pid_t pid = fork();
  if (pid < 0) {
    result.error = "could not start process";
    return result;
  }
  if (pid == 0) {
    std::vector<char *> argv;
    argv.reserve(args.size() + 1);
    for (const auto& arg : args) {
      argv.push_back(const_cast<char *>(arg.c_str()));
    }
    argv.push_back(nullptr);
    execvp(argv[0], argv.data());
    _exit(127);
  }

  int status = 0;
  while (waitpid(pid, &status, 0) == -1) {
    if (errno != EINTR) {
      result.error = "could not wait for process";
      return result;
    }
  }
  if (status != 0) {
    result.error = processStatusMessage(status);
    return result;
  }
  result.ok = true;
  return result;
#endif
}

bool isPathInside(const fs::path& child, const fs::path& parent)
{
  std::error_code ec;
  const fs::path canonicalChild = fs::weakly_canonical(child, ec);
  if (ec) return false;
  const fs::path canonicalParent = fs::weakly_canonical(parent, ec);
  if (ec) return false;
  const fs::path relative = canonicalChild.lexically_relative(canonicalParent);
  const std::string relativeString = relative.generic_string();
  return !relative.empty() && relativeString.rfind("..", 0) != 0 && !relative.is_absolute();
}

bool hasManifestControlCharacter(const std::string& value)
{
  return value.find('\t') != std::string::npos || value.find('\r') != std::string::npos ||
         value.find('\n') != std::string::npos || value.find('\0') != std::string::npos;
}

bool isReservedWindowsManifestPathComponent(const std::string& component)
{
  std::string stem = component.substr(0, component.find('.'));
  std::transform(stem.begin(), stem.end(), stem.begin(),
                 [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
  static constexpr std::array<const char *, 22> reserved = {
    "con",  "prn",  "aux",  "nul",  "com1", "com2", "com3", "com4", "com5", "com6", "com7",
    "com8", "com9", "lpt1", "lpt2", "lpt3", "lpt4", "lpt5", "lpt6", "lpt7", "lpt8", "lpt9",
  };
  return std::find(reserved.begin(), reserved.end(), stem) != reserved.end();
}

bool isSafeManifestRelativePath(const std::string& value)
{
  if (hasManifestControlCharacter(value) || value.empty() || value[0] == '/') return false;
  if (value.size() >= 2 && std::isalpha(static_cast<unsigned char>(value[0])) && value[1] == ':') {
    return false;
  }
  if (value.rfind("\\\\", 0) == 0) return false;

  const fs::path normalized = fs::path(value).lexically_normal();
  if (normalized.empty() || normalized == "." || normalized.is_absolute()) return false;
  for (const auto& part : normalized) {
    const std::string component = part.generic_string();
    if (component.empty() || component == "." || component == "..") return false;
    if (isReservedWindowsManifestPathComponent(component)) return false;
  }
  return true;
}

fs::path findWasmBundleDir()
{
  const auto fromEnv = getenvString("PYTHONSCAD_WASM_DIR");
  if (!fromEnv.empty() && hasWasmBundle(fromEnv)) return fs::path(fromEnv);

  const fs::path cwd = fs::current_path();
  if (hasWasmBundle(cwd / "build-wasm-web")) return cwd / "build-wasm-web";

  fs::path appPath = PlatformUtils::applicationPath();
  if (hasWasmBundle(appPath / "wasm")) return appPath / "wasm";
  if (hasWasmBundle(appPath)) return appPath;

  return {};
}

fs::path findSandboxRunner()
{
  const auto fromEnv = getenvString("PYTHONSCAD_SANDBOX_RUNNER");
  if (!fromEnv.empty() && fs::is_regular_file(fromEnv)) return fs::path(fromEnv);

  fs::path cwdRunner = fs::current_path() / "scripts" / "python-sandbox-runner.mjs";
  if (fs::is_regular_file(cwdRunner)) return cwdRunner;

  fs::path appRunner =
    fs::path(PlatformUtils::applicationPath()) / "scripts" / "python-sandbox-runner.mjs";
  if (fs::is_regular_file(appRunner)) return appRunner;

  return {};
}

fs::path makeTempDir()
{
  std::random_device random;
  const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
  for (int attempt = 0; attempt < 100; ++attempt) {
    fs::path dir = fs::temp_directory_path() /
                   ("pythonscad-sandbox-" + std::to_string(stamp) + "-" + std::to_string(random()));
    std::error_code ec;
    if (fs::create_directory(dir, ec)) return dir;
    if (ec && ec != std::make_error_code(std::errc::file_exists)) break;
  }
  throw std::runtime_error("Could not create sandbox temporary directory.");
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

bool readManifest(const fs::path& manifestFile, const fs::path& outputRoot,
                  std::vector<PythonSandboxOutputFile>& files, std::string& error)
{
  std::string manifest;
  if (!readFile(manifestFile, manifest)) return true;

  std::istringstream stream(manifest);
  std::string line;
  while (std::getline(stream, line)) {
    if (line.empty()) continue;
    const size_t firstTab = line.find('\t');
    const size_t secondTab =
      firstTab == std::string::npos ? std::string::npos : line.find('\t', firstTab + 1);
    if (firstTab == std::string::npos || secondTab == std::string::npos) continue;

    PythonSandboxOutputFile file;
    file.relativePath = line.substr(0, firstTab);
    file.hostPath = line.substr(firstTab + 1, secondTab - firstTab - 1);
    if (!isSafeManifestRelativePath(file.relativePath) || !isPathInside(file.hostPath, outputRoot)) {
      error = "Sandboxed Python produced an invalid output manifest.";
      return false;
    }
    try {
      file.size = std::stoull(line.substr(secondTab + 1));
    } catch (const std::exception&) {
      file.size = 0;
    }
    files.push_back(std::move(file));
  }
  return true;
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
  fs::path tempDir;
  try {
    tempDir = makeTempDir();
  } catch (const std::exception& e) {
    result.error = e.what();
    return result;
  }
  result.tempDir = tempDir.string();
  const fs::path inputFile =
    tempDir /
    (fs::path(scriptpath).filename().empty() ? fs::path("input.py") : fs::path(scriptpath).filename());
  const fs::path outputFile = tempDir / "output.csg";
  const fs::path outputRoot = tempDir / "sandbox-outputs";
  const fs::path manifestFile = tempDir / "outputs.tsv";
  std::error_code ec;
  if (!fs::create_directories(outputRoot, ec) || ec) {
    result.error = "Could not create sandbox output directory.";
    fs::remove_all(tempDir);
    return result;
  }
  result.outputRoot = outputRoot.string();

  if (!writeFile(inputFile, code)) {
    result.error = "Could not write sandbox input file.";
    fs::remove_all(tempDir);
    return result;
  }

  const ProcessResult processResult =
    runProcess({node, runner.string(), "--wasm-dir", wasmDir.string(), "--input", inputFile.string(),
                "--output", outputFile.string(), "--host-output-dir", outputRoot.string(), "--manifest",
                manifestFile.string()});
  if (!processResult.ok) {
    result.error = "Sandboxed Python process failed with " + processResult.error + ".";
    fs::remove_all(tempDir);
    return result;
  }

  if (!readFile(outputFile, result.csg)) {
    result.error = "Sandboxed Python did not produce CSG output.";
    fs::remove_all(tempDir);
    return result;
  }

  if (!readManifest(manifestFile, outputRoot, result.outputFiles, result.error)) {
    fs::remove_all(tempDir);
    return result;
  }
  result.ok = true;
  return result;
}

void cleanupPythonSandboxResult(const PythonSandboxResult& result)
{
  if (!result.tempDir.empty()) {
    std::error_code ec;
    fs::remove_all(result.tempDir, ec);
  }
}

#endif  // ENABLE_PYTHON
