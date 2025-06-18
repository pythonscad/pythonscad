/*
 *  OpenSCAD (www.openscad.org)
 *  Copyright (C) 2009-2025 Clifford Wolf <clifford@clifford.at> and
 *                          Marius Kintel <marius@kintel.net>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  As a special exception, you have permission to link this program
 *  with the CGAL library and distribute executables, as long as you
 *  follow the requirements of the GNU GPL in regard to all of the
 *  software in the executable aside from CGAL.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */
#include <Python.h>

#include <string>
#include <vector>
#include <cstdlib>
#include <filesystem>

#include "core/Settings.h"
#include "platform/PlatformUtils.h"

#include "pyopenscad.h"

namespace fs = std::filesystem;

using SP = Settings::SettingsPython;

std::string venvBinDirFromSettings()
{
    const auto& venv = fs::path(SP::pythonVirtualEnv.value()) / "bin";
    if (fs::is_directory(venv)) {
        return venv.generic_string();
    }
    return "";
}

int pythonRunArgs(int argc, char **argv)
{
  PyStatus status;

  PyConfig config;
  PyConfig_InitPythonConfig(&config);

  status = PyConfig_SetBytesArgv(&config, argc, argv);
  if (PyStatus_Exception(status)) {
    goto fail;
  }

  status = Py_InitializeFromConfig(&config);
  if (PyStatus_Exception(status)) {
    goto fail;
  }
  PyConfig_Clear(&config);

  return Py_RunMain();

fail:
  PyConfig_Clear(&config);
  if (PyStatus_IsExit(status)) {
    return status.exitcode;
  }
  Py_ExitStatusException(status);
}

int pythonCreateVenv(const std::string& path)
{
  int result = pythonRunModule("", "venv", { path });
  if (result != 0) {
    return result;
  }

  // The created VENV points to the temporary mount point of the
  // AppImage, e.g. /tmp/.mount_OpenSCCpPaio - that is obviously
  // no good for any later runs.
  // To fix that, we point the link to the magic /proc/self/exe
  // so it can always just call itself as the python interpreter.
  const char *appdirenv = getenv("APPDIR");
  if (getenv("APPIMAGE") != nullptr && appdirenv != nullptr) {
    // Assume we are running as AppImage
    const std::string appdir = appdirenv;
    const auto vbin = fs::path{path} / "bin" / "openscad-python";
    if (fs::exists(vbin) && fs::is_symlink(vbin)) {
      const auto lbin = fs::read_symlink(vbin).generic_string();
      if (lbin.rfind(appdir, 0) == 0) {
        std::error_code ec;
        fs::remove(vbin, ec);
        if (ec.value() > 0) {
          return ec.value();
        }
        fs::create_symlink("/proc/self/exe", vbin, ec);
        if (ec.value() > 0) {
          return ec.value();
        }
      }
    }
  }

  return 0;
}

int pythonRunModule(const std::string& appPath, const std::string& module,
                    const std::vector<std::string>& args)
{
  PyStatus status;
  const auto name = "openscad-python";
  const auto exe = PlatformUtils::applicationPath() + "/" + name;

  PyPreConfig preconfig;
  PyPreConfig_InitPythonConfig(&preconfig);

  status = Py_PreInitialize(&preconfig);
  if (PyStatus_Exception(status)) {
    Py_ExitStatusException(status);
  }

  PyConfig config;
  PyConfig_InitPythonConfig(&config);

  status = PyConfig_SetBytesString(&config, &config.program_name, name);
  if (PyStatus_Exception(status)) {
    goto done;
  }

  status = PyConfig_SetBytesString(&config, &config.executable, exe.c_str());
  if (PyStatus_Exception(status)) {
    goto done;
  }

  status = PyConfig_SetBytesString(&config, &config.run_module, module.c_str());
  if (PyStatus_Exception(status)) {
    goto done;
  }

  /* Read all configuration at once */
  status = PyConfig_Read(&config);
  if (PyStatus_Exception(status)) {
    goto done;
  }

  for (const auto& arg : args) {
    std::wstring warg(arg.size(), L' ');
    warg.resize(std::mbstowcs(&warg[0], arg.c_str(), arg.size()));
    status = PyWideStringList_Append(&config.argv, warg.c_str());
    if (PyStatus_Exception(status)) {
      goto done;
    }
  }

  /* Override executable computed by PyConfig_Read() */
  status = PyConfig_SetBytesString(&config, &config.executable, exe.c_str());
  if (PyStatus_Exception(status)) {
    goto done;
  }

  status = Py_InitializeFromConfig(&config);
  if (PyStatus_Exception(status)) {
    goto done;
  }

  return Py_RunMain();

done:
  PyConfig_Clear(&config);
  return status.exitcode;
}

int curl_download(std::string url, std::string path);

int pip_bootstrap(void) // will be called from separate program
{
  std::string tempfile = std::tmpnam(nullptr); // TODO dangerous
  char *tempfile_wr = strdup(tempfile.c_str());

  // curl -sSL https://bootstrap.pypa.io/get-pip.py -o get-pip.py
  if(curl_download("https://bootstrap.pypa.io/get-pip.py", tempfile)) {
    LOG("Error downloading pip");
    return 1;
  }

  PyConfig config;
  PyConfig_InitPythonConfig(&config);
  std::string targetpath=PlatformUtils::resourceBasePath()+"/libraries/python/site-packages";
  char *targetpath_wr = strdup(targetpath.c_str());
  
  // python get-pip.py --target /home/gsohler/git/pythonscad/libraries/python/site-packages
  char * const argv[4]= { "", tempfile_wr, "--target", strdup(targetpath_wr) };
  PyConfig_SetBytesArgv(&config, 4, argv);
  config.parse_argv=1;

  PyStatus status;
  status = Py_InitializeFromConfig(&config);
  PyConfig_Clear(&config);

  FILE *fp =fopen(tempfile.c_str(),"rb");
  if(fp == nullptr) {
    LOG("Cannot open tempfile");
    return 1;
  }
  PyCompilerFlags cf = _PyCompilerFlags_INIT;
  PyRun_AnyFileFlags(fp, "<stdin>", &cf);
  fclose(fp);

  unlink(tempfile_wr); // TODO fail
  Py_Finalize();

  free(targetpath_wr);
  free(tempfile_wr);
  exit(0);
}
int pip_ensure(void) {
  std::string pip_path=PlatformUtils::resourceBasePath()+"/libraries/python/site-packages/pip";
  if(std::filesystem::exists(pip_path)) return 0;
  int status = system("./pythonscad --pip-bootstrap"); // TODO improve
  printf("status is %d\n",status);						       
  if(status) {
    LOG("Cannot bootstrap pip");
    return -1;
  }	  
  if(!std::filesystem::exists(pip_path)){
    LOG("Pip did not install as expected");
    return -1;
  }
  return 0;
}
int pip_install(const std::string &modname){
  PyConfig config;
  PyConfig_InitPythonConfig(&config);
  std::string targetpath=PlatformUtils::resourceBasePath()+"/libraries/python/site-packages";

  char *targetpath_wr = strdup(targetpath.c_str());

  char *modname_wr = strdup(modname.c_str());
  
  //python -m pip install --target ../libraries/python/site-packages numpy
  char * const argv[7]= { "", "pip", "install", "--upgrade", "--target", targetpath_wr,  strdup(modname_wr) };
  PyConfig_SetBytesArgv(&config, 7, argv);
  config.parse_argv=1;

  PyStatus status;
  status = Py_InitializeFromConfig(&config);
  PyConfig_Clear(&config);

  PyObject *module, *runpy, *runmodule, *runargs, *result;

  runpy = PyImport_ImportModule("runpy");
  if (runpy == NULL) {
      fprintf(stderr, "Could not import runpy module\n");
      exit(1);
  }
  runmodule = PyObject_GetAttrString(runpy, "_run_module_as_main");
  if (runmodule == NULL) {
      fprintf(stderr, "Could not access runpy._run_module_as_main\n");
      Py_DECREF(runpy);
      exit(1);
  }
  module = PyUnicode_FromString("pip");
  if (module == NULL) {
      fprintf(stderr, "Could not convert module name to unicode\n");
      Py_DECREF(runpy);
      Py_DECREF(runmodule);
      exit(1);
  }
  runargs = PyTuple_Pack(2, module, 1 ? Py_True : Py_False);
  if (runargs == NULL) {
      fprintf(stderr,
          "Could not create arguments for runpy._run_module_as_main\n");
      Py_DECREF(runpy);
      Py_DECREF(runmodule);
      Py_DECREF(module);
      exit(1);
  }
  result = PyObject_Call(runmodule, runargs, NULL);
  free(targetpath_wr);
  free(modname_wr);
  Py_DECREF(runpy);
  Py_DECREF(runmodule);
  Py_DECREF(module);
  Py_DECREF(runargs);
  if (result == NULL) {
	  exit(1);
  }
  Py_DECREF(result);
  exit(0);
  return 0;
}
// TODO make sure that only own modules are used

int pip_install_call(const std::string &package){
  std::string command="./pythonscad --pip-install "+package;
  int status = system(command.c_str()); // TODO improve
  return 0;
}


