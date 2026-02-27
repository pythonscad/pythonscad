/*
 *  OpenSCAD (www.openscad.org)
 *  Copyright (C) 2009-2011 Clifford Wolf <clifford@clifford.at> and
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

#include "openscad_gui.h"

#include <QtCore/qstringliteral.h>

#include <QByteArray>
#include <QCryptographicHash>
#include <QDialog>
#include <QDir>
#include <QFileInfo>
#include <QFutureWatcher>
#include <QGuiApplication>
#include <QIcon>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLocalServer>
#include <QLocalSocket>
#include <QLockFile>
#include <QMessageBox>
#include <QObject>
#include <QPalette>
#include <QDesktopServices>
#include <QFile>
#include <QSessionManager>
#include <QSaveFile>
#include <QSocketNotifier>
#include <QStandardPaths>
#include <QStringList>
#include <QTimer>
#include <QStyleHints>
#include <QVector>
#include <Qt>
#include <QtConcurrentRun>
#include <QtGlobal>
#include <array>
#include <cstring>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>
#include <QFile>

#ifdef Q_OS_UNIX
#include <fcntl.h>
#include <signal.h>
#include <unistd.h>
#endif

#include "Feature.h"
#include "FontCache.h"
#include "core/Settings.h"
#include "core/parsersettings.h"
#include "geometry/Geometry.h"
#include "gui/AppleEvents.h"
#include "gui/input/InputDriverManager.h"
#include "version.h"
#ifdef ENABLE_HIDAPI
#include "gui/input/HidApiInputDriver.h"
#endif
#ifdef ENABLE_SPNAV
#include "gui/input/SpaceNavInputDriver.h"
#endif
#ifdef ENABLE_JOYSTICK
#include "gui/input/JoystickInputDriver.h"
#endif
#ifdef ENABLE_DBUS
#include "gui/input/DBusInputDriver.h"
#endif
#ifdef ENABLE_QGAMEPAD
#include "gui/input/QGamepadInputDriver.h"
#endif
#include "gui/LaunchingScreen.h"
#include "gui/MainWindow.h"
#include "gui/OpenSCADApp.h"
#include "gui/Preferences.h"
#include "gui/QSettingsCached.h"
#include "openscad.h"
#include "platform/CocoaUtils.h"
#include "utils/printutils.h"

#ifdef ENABLE_GUI_TESTS
#include "guitests/guitests.h"
#endif

Q_DECLARE_METATYPE(Message);
Q_DECLARE_METATYPE(std::shared_ptr<const Geometry>);

extern std::string arg_colorscheme;

namespace {

constexpr int kAutosaveIntervalMs = 60000;

bool shouldOfferAutosaveRestore(const QString& autosavePath, const QString& sessionPath)
{
  const QFileInfo autosaveInfo(autosavePath);
  if (!autosaveInfo.exists()) return false;

  const QFileInfo sessionInfo(sessionPath);
  if (!sessionInfo.exists()) return true;

  return autosaveInfo.lastModified() > sessionInfo.lastModified();
}

bool promptAutosaveRestore(const QString& autosavePath)
{
  while (true) {
    QMessageBox box;
    box.setIcon(QMessageBox::Warning);
    box.setWindowTitle(QStringLiteral("PythonSCAD"));
    box.setText(_("Recovered session data was found."));
    box.setInformativeText(_("Restore the autosaved session?"));
    auto *restoreButton = box.addButton(_("Restore"), QMessageBox::AcceptRole);
    auto *discardButton = box.addButton(_("Discard"), QMessageBox::RejectRole);
    auto *showButton = box.addButton(_("Show File"), QMessageBox::ActionRole);
    box.setDefaultButton(restoreButton);
    box.exec();

    if (box.clickedButton() == showButton) {
      const QString dirPath = QFileInfo(autosavePath).absolutePath();
      QDesktopServices::openUrl(QUrl::fromLocalFile(dirPath));
      continue;
    }

    if (box.clickedButton() == restoreButton) {
      return true;
    }

    return false;
  }
}

void setupAutosaveTimer(OpenSCADApp *app)
{
  auto *timer = new QTimer(app);
  timer->setInterval(kAutosaveIntervalMs);

  auto *state = new uint64_t(0);
  auto *performed = new bool(false);

  QObject::connect(timer, &QTimer::timeout, app, [state, performed, timer]() {
    const int intervalSeconds = Settings::Settings::autosaveSessionIntervalSeconds.value();
    const int intervalMs = std::max(10, intervalSeconds) * 1000;
    if (timer->interval() != intervalMs) {
      timer->setInterval(intervalMs);
    }

    const bool enabled = Settings::Settings::autosaveSessionEnabled.value();

    if (!enabled) {
      if (*performed) {
        QFile::remove(TabManager::getAutosaveFilePath());
        *performed = false;
        *state = 0;
      }
      return;
    }

    const bool dirty = TabManager::hasDirtyTabs();
    if (!dirty) {
      if (*performed) {
        QFile::remove(TabManager::getAutosaveFilePath());
        *performed = false;
        *state = 0;
      }
      return;
    }

    const uint64_t generation = TabManager::sessionDirtyGeneration();
    if (!*performed || generation != *state) {
      TabManager::saveGlobalSession(TabManager::getAutosaveFilePath());
      *state = generation;
      *performed = true;
    }
  });

  timer->start();
}

bool saveSessionForShutdown()
{
  const auto& windows = scadApp->windowManager.getWindows();
  if (windows.isEmpty()) {
    return false;
  }
  if (TabManager::shouldSkipSessionSave()) {
    QFile::remove(TabManager::getAutosaveFilePath());
    return false;
  }
  for (auto *win : windows) {
    win->markSessionQuitting();
  }
  TabManager::saveGlobalSession(TabManager::getSessionFilePath());
  QFile::remove(TabManager::getAutosaveFilePath());
  return true;
}

constexpr int kIpcTimeoutMs = 1500;

QString lockFilePath()
{
  const QString baseDir = TabManager::getSessionFilePath();
  return QFileInfo(baseDir).absolutePath() + QStringLiteral("/pythonscad.lock");
}

QString serverNameFromPath(const QString& path)
{
  const QByteArray hash = QCryptographicHash::hash(path.toUtf8(), QCryptographicHash::Sha1).toHex();
  return QStringLiteral("pythonscad.instance.") + QString::fromUtf8(hash);
}

QString serverName()
{
  return serverNameFromPath(lockFilePath());
}

QString resolveOpenMode(const std::string& overrideMode)
{
  if (overrideMode == "new-window" || overrideMode == "active-window") {
    return QString::fromStdString(overrideMode);
  }
  return QString::fromStdString(Settings::Settings::singleInstanceOpenMode.value());
}

void focusWindow(MainWindow *window)
{
  if (!window) return;
  window->show();
  window->raise();
  window->activateWindow();
}

void openFilesInWindow(MainWindow *window, const QStringList& files)
{
  if (!window) return;
  for (const auto& file : files) {
    if (file.isEmpty()) continue;
    window->tabManager->open(file);
  }
  focusWindow(window);
}

MainWindow *getOrCreateActiveWindow()
{
  auto *active = scadApp->windowManager.getLastActive();
  if (active) return active;
  if (!scadApp->windowManager.getWindows().isEmpty()) {
    return *scadApp->windowManager.getWindows().begin();
  }
  return nullptr;
}

QJsonObject buildIpcMessage(const QString& action, const QStringList& files, const QString& openMode,
                            const QString& cwd)
{
  QJsonObject obj;
  obj.insert(QStringLiteral("action"), action);
  obj.insert(QStringLiteral("openMode"), openMode);
  obj.insert(QStringLiteral("cwd"), cwd);
  QJsonArray fileArray;
  for (const auto& file : files) fileArray.append(file);
  obj.insert(QStringLiteral("files"), fileArray);
  return obj;
}

bool sendIpcMessage(const QJsonObject& message)
{
  QLocalSocket socket;
  socket.connectToServer(serverName());
  if (!socket.waitForConnected(kIpcTimeoutMs)) {
    return false;
  }

  const QByteArray payload = QJsonDocument(message).toJson(QJsonDocument::Compact);
  socket.write(payload);
  if (!socket.waitForBytesWritten(kIpcTimeoutMs)) {
    return false;
  }
  socket.flush();
  if (!socket.waitForReadyRead(kIpcTimeoutMs)) {
    return false;
  }
  const QByteArray ack = socket.readAll();
  socket.disconnectFromServer();
  return ack.startsWith("ok");
}

void startIpcServer(QLocalServer *server)
{
  QObject::connect(server, &QLocalServer::newConnection, [server]() {
    while (auto *socket = server->nextPendingConnection()) {
      QObject::connect(socket, &QLocalSocket::readyRead, [socket]() {
        const QByteArray data = socket->readAll();
        const auto doc = QJsonDocument::fromJson(data);
        bool ok = doc.isObject();
        if (!doc.isObject()) {
          socket->write("error");
          socket->flush();
          socket->disconnectFromServer();
          return;
        }
        const auto obj = doc.object();
        const auto action = obj.value(QStringLiteral("action")).toString();
        const auto openMode = obj.value(QStringLiteral("openMode")).toString();
        const auto filesValue = obj.value(QStringLiteral("files"));

        if (action == QStringLiteral("focus")) {
          focusWindow(getOrCreateActiveWindow());
          socket->write("ok");
          socket->flush();
          socket->disconnectFromServer();
          return;
        }

        if (action != QStringLiteral("open") || !filesValue.isArray()) {
          ok = false;
          socket->write("error");
          socket->flush();
          socket->disconnectFromServer();
          return;
        }

        QStringList files;
        for (const auto& entry : filesValue.toArray()) {
          const auto path = entry.toString();
          if (!path.isEmpty()) files.append(path);
        }

        if (files.isEmpty()) {
          focusWindow(getOrCreateActiveWindow());
          socket->write("ok");
          socket->flush();
          socket->disconnectFromServer();
          return;
        }

        if (openMode == QStringLiteral("active-window")) {
          auto *active = getOrCreateActiveWindow();
          if (active) {
            openFilesInWindow(active, files);
          } else {
            new MainWindow(files);
          }
          socket->write("ok");
          socket->flush();
          socket->disconnectFromServer();
          return;
        }

        for (const auto& file : files) {
          new MainWindow(QStringList(file));
        }
        socket->write(ok ? "ok" : "error");
        socket->flush();
        socket->disconnectFromServer();
      });
    }
  });

  if (!server->listen(serverName())) {
    if (server->serverError() == QAbstractSocket::AddressInUseError) {
      QLocalServer::removeServer(serverName());
      server->listen(serverName());
    }
  }
}

#ifdef Q_OS_UNIX
int shutdownSignalPipe[2] = {-1, -1};

void shutdownSignalHandler(int)
{
  const char signalByte = 1;
  if (shutdownSignalPipe[1] != -1) {
    (void)::write(shutdownSignalPipe[1], &signalByte, sizeof(signalByte));
  }
}

void setupUnixSignalHandlers(OpenSCADApp *app)
{
  if (::pipe(shutdownSignalPipe) != 0) return;
  const int flags = fcntl(shutdownSignalPipe[1], F_GETFL, 0);
  if (flags != -1) {
    fcntl(shutdownSignalPipe[1], F_SETFL, flags | O_NONBLOCK);
  }

  auto *notifier = new QSocketNotifier(shutdownSignalPipe[0], QSocketNotifier::Read, app);
  QObject::connect(notifier, &QSocketNotifier::activated, app, [notifier](int) {
    notifier->setEnabled(false);
    char buffer[32];
    while (::read(shutdownSignalPipe[0], buffer, sizeof(buffer)) > 0) {
    }
    saveSessionForShutdown();
    QCoreApplication::quit();
    notifier->setEnabled(true);
  });

  struct sigaction action;
  memset(&action, 0, sizeof(action));
  action.sa_handler = shutdownSignalHandler;
  sigaction(SIGTERM, &action, nullptr);
  sigaction(SIGINT, &action, nullptr);
  sigaction(SIGHUP, &action, nullptr);
}
#endif

// Only if "fileName" is not absolute, prepend the "absoluteBase".
QString assemblePath(const std::filesystem::path& absoluteBaseDir, const std::string& fileName)
{
  if (fileName.empty()) return "";
  auto qsDir = QString::fromStdString(absoluteBaseDir.generic_string());
  auto qsFile = QString::fromStdString(fileName);
  // if qsfile is absolute, dir is ignored. (see documentation of QFileInfo)
  const QFileInfo fileInfo(qsDir, qsFile);
  return fileInfo.absoluteFilePath();
}

void dialogThreadFunc(FontCacheInitializer *initializer)
{
  initializer->run();
}

void dialogInitHandler(FontCacheInitializer *initializer, void *)
{
  QFutureWatcher<void> futureWatcher;
  QObject::connect(&futureWatcher, &QFutureWatcher<void>::finished, scadApp,
                   &OpenSCADApp::hideFontCacheDialog);

  auto future = QtConcurrent::run([initializer] { return dialogThreadFunc(initializer); });
  futureWatcher.setFuture(future);

  // We don't always get the started() signal, so we start manually
  QMetaObject::invokeMethod(scadApp, "showFontCacheDialog");

  // Block, in case we're in a separate thread, or the dialog was closed by the user
  futureWatcher.waitForFinished();

  // We don't always receive the finished signal. We still need the signal to break
  // out of the exec() though.
  QMetaObject::invokeMethod(scadApp, "hideFontCacheDialog");
}

#ifdef Q_OS_WIN
void registerDefaultIcon(QString applicationFilePath)
{
  // Not using cached instance here, so this needs to be in a
  // separate scope to ensure the QSettings instance is released
  // directly after use.
  QSettings reg_setting(QLatin1String("HKEY_CURRENT_USER"), QSettings::NativeFormat);
  auto appPath = QDir::toNativeSeparators(applicationFilePath + QLatin1String(",1"));
  reg_setting.setValue(QLatin1String("Software/Classes/OpenSCAD_File/DefaultIcon/Default"),
                       QVariant(appPath));
}
#else
void registerDefaultIcon(const QString&)
{
}
#endif

}  // namespace

#ifdef OPENSCAD_SUFFIX
#define DESKTOP_FILENAME "pythonscad" OPENSCAD_SUFFIX
#else
#define DESKTOP_FILENAME "pythonscad"
#endif

int gui(std::vector<std::string>& inputFiles, const std::filesystem::path& original_path, int argc,
        char **argv, const std::string& gui_test, const bool reset_window_settings,
        const std::string& open_in_override)
{
  OpenSCADApp app(argc, argv);

  // set up groups for QSettings
  QCoreApplication::setOrganizationName("PythonSCAD");
  QCoreApplication::setOrganizationDomain("pythonscad.org");
  QCoreApplication::setApplicationName("PythonSCAD");
  QCoreApplication::setApplicationVersion(QString::fromStdString(std::string(openscad_versionnumber)));
  QGuiApplication::setApplicationDisplayName("PythonSCAD");
  QGuiApplication::setDesktopFileName(DESKTOP_FILENAME);
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
  QCoreApplication::setAttribute(Qt::AA_UseHighDpiPixmaps);
#endif

#ifdef Q_OS_MACOS
  app.setWindowIcon(QIcon(":/icon-macos.png"));
#else
  app.setWindowIcon(QIcon(":/logo.png"));
#endif

  // Other global settings
  qRegisterMetaType<Message>();
  qRegisterMetaType<std::shared_ptr<const Geometry>>();

  FontCache::registerProgressHandler(dialogInitHandler);

  parser_init();

  const QString openMode = resolveOpenMode(open_in_override);
  const QString cwd = QString::fromStdString(original_path.generic_string());

  QLockFile lock(lockFilePath());
  lock.setStaleLockTime(0);
  if (!lock.tryLock()) {
    lock.removeStaleLockFile();
    if (!lock.tryLock()) {
      const QStringList ipcFiles = [&]() {
        QStringList files;
        for (const auto& infile : inputFiles) {
          if (!infile.empty()) {
            files.append(assemblePath(original_path, infile));
          }
        }
        return files;
      }();

      const QString action = ipcFiles.isEmpty() ? QStringLiteral("focus") : QStringLiteral("open");
      const QJsonObject message = buildIpcMessage(action, ipcFiles, openMode, cwd);

      while (true) {
        if (sendIpcMessage(message)) {
          return 0;
        }

        QMessageBox box;
        box.setIcon(QMessageBox::Warning);
        box.setWindowTitle(QStringLiteral("PythonSCAD"));
        box.setText(_("PythonSCAD is already running but is not responding."));
        box.setInformativeText(_("Retry to send the request, or exit."));
        const auto retryButton = box.addButton(_("Retry"), QMessageBox::AcceptRole);
        const auto exitButton = box.addButton(_("Exit"), QMessageBox::RejectRole);
        box.setDefaultButton(retryButton);
        box.exec();
        if (box.clickedButton() == exitButton) {
          return 1;
        }
      }
    }
  }

  QLocalServer ipcServer;
  startIpcServer(&ipcServer);

  QSettingsCached settings;
  if (settings.value("advanced/localization", true).toBool()) {
    localization_init();
  }
  if (reset_window_settings) {
    const auto keys = std::array<std::string, 20>{
      "editor/fontfamily",
      "editor/fontsize",
      "advanced/applicationFontSize",
      "advanced/applicationFontFamily",
      "advanced/consoleFontFamily",
      "advanced/consoleFontSize",
      "advanced/customizerFontFamily",
      "advanced/customizerFontSize",
      "advanced/undockableWindows",
      "window/state",
      "window/geometry",
      "window/position",
      "window/size",
      "view/hideEditor",
      "view/hideConsole",
      "view/hideErrorLog",
      "view/hideAnimate",
      "view/hideCustomizer",
      "view/hideFontList",
      "view/hideViewportControl",
    };
    for (const auto& key : keys) {
      settings.remove(QString::fromStdString(key));
    }
  }

#ifdef Q_OS_MACOS
  installAppleEventHandlers();
#endif

  registerDefaultIcon(app.applicationFilePath());
  app.setGuiTheme(GlobalPreferences::inst()->getValue("advanced/guiTheme").toString());

#ifdef OPENSCAD_UPDATER
  AutoUpdater *updater = new SparkleAutoUpdater;
  AutoUpdater::setUpdater(updater);
  if (updater->automaticallyChecksForUpdates()) updater->checkForUpdates();
  updater->init();
#endif

  QObject::connect(GlobalPreferences::inst(), &Preferences::applicationFontChanged, &app,
                   &OpenSCADApp::setApplicationFont);
  QObject::connect(GlobalPreferences::inst(), &Preferences::guiThemeChanged, &app,
                   &OpenSCADApp::setGuiTheme);
  QObject::connect(GlobalPreferences::inst(), &Preferences::renderBackend3DChanged, &app,
                   &OpenSCADApp::setRenderBackend3D);

  set_render_color_scheme(arg_colorscheme, false);
  auto noInputFiles = false;

  if (!inputFiles.size()) {
    noInputFiles = true;
    inputFiles.emplace_back("");
  }

  QStringList inputFilesList;
  for (const auto& infile : inputFiles) {
    inputFilesList.append(assemblePath(original_path, infile));
  }

  QVector<QStringList> windowsToOpen;

  if (noInputFiles) {
    const QString sessionPath = TabManager::getSessionFilePath();
    const QString autosavePath = TabManager::getAutosaveFilePath();
    if (shouldOfferAutosaveRestore(autosavePath, sessionPath)) {
      if (promptAutosaveRestore(autosavePath)) {
        QFile in(autosavePath);
        QSaveFile out(sessionPath);
        if (in.open(QIODevice::ReadOnly) && out.open(QIODevice::WriteOnly)) {
          out.write(in.readAll());
          const bool committed = out.commit();
          if (committed) {
            QFile::remove(autosavePath);
          }
        } else {
          out.cancelWriting();
        }
      } else {
        QFile::remove(autosavePath);
      }
    }
  }

  // When no files given, restore session if it exists (skip launcher)
  if (inputFilesList.size() == 1 && inputFilesList[0].isEmpty()) {
    const QString sessionPath = TabManager::getSessionFilePath();
    if (QFileInfo(sessionPath).exists()) {
      const int windowCount = TabManager::sessionWindowCount(sessionPath);
      if (windowCount > 0) {
        for (int i = 0; i < windowCount; ++i) {
          windowsToOpen.append(QStringList(QStringLiteral(":session:%1:").arg(i)));
        }
      } else {
        windowsToOpen.append(QStringList(QStringLiteral(":session:")));
      }
    }
  }

  if (windowsToOpen.isEmpty()) {
    // Show launcher only when no files and no session to restore
    if (noInputFiles && inputFilesList.size() == 1 && inputFilesList[0].isEmpty()) {
      auto showOnStartup = settings.value("launcher/showOnStartup");
      if (showOnStartup.isNull() || showOnStartup.toBool()) {
        LaunchingScreen launcher;
        if (launcher.exec() == QDialog::Accepted) {
          if (launcher.isForceShowEditor()) {
            settings.setValue("view/hideEditor", false);
          }
          const QStringList files = launcher.selectedFiles();
          if (!files.empty()) {
            inputFilesList.clear();
            for (const auto& f : files) {
              inputFilesList.append(assemblePath(original_path, f.toStdString()));
            }
          }
        } else {
          return 0;
        }
      }
    }
    windowsToOpen.append(inputFilesList);
  }

  for (const auto& files : windowsToOpen) {
    new MainWindow(files);
  }

  setupAutosaveTimer(&app);

  QObject::connect(&app, &QCoreApplication::aboutToQuit, []() {
    saveSessionForShutdown();
    QSettingsCached{}.release();
#ifdef Q_OS_MACOS
    CocoaUtils::endApplication();
#endif
  });

  QObject::connect(&app, &QGuiApplication::commitDataRequest, &app,
                   [](QSessionManager&) { saveSessionForShutdown(); });
  QObject::connect(&app, &QGuiApplication::saveStateRequest, &app,
                   [](QSessionManager&) { saveSessionForShutdown(); });

#ifdef Q_OS_UNIX
  setupUnixSignalHandlers(&app);
#endif

#ifdef ENABLE_HIDAPI
  if (Settings::Settings::inputEnableDriverHIDAPI.value()) {
    auto hidApi = new HidApiInputDriver();
    InputDriverManager::instance()->registerDriver(hidApi);
  }
#endif
#ifdef ENABLE_SPNAV
  if (Settings::Settings::inputEnableDriverSPNAV.value()) {
    auto spaceNavDriver = new SpaceNavInputDriver();
    bool spaceNavDominantAxisOnly = Settings::Settings::inputEnableDriverHIDAPI.value();
    spaceNavDriver->setDominantAxisOnly(spaceNavDominantAxisOnly);
    InputDriverManager::instance()->registerDriver(spaceNavDriver);
  }
#endif
#ifdef ENABLE_JOYSTICK
  if (Settings::Settings::inputEnableDriverJOYSTICK.value()) {
    std::string nr = STR(Settings::Settings::joystickNr.value());
    auto joyDriver = new JoystickInputDriver();
    joyDriver->setJoystickNr(nr);
    InputDriverManager::instance()->registerDriver(joyDriver);
  }
#endif
#ifdef ENABLE_QGAMEPAD
  if (Settings::Settings::inputEnableDriverQGAMEPAD.value()) {
    auto qGamepadDriver = new QGamepadInputDriver();
    InputDriverManager::instance()->registerDriver(qGamepadDriver);
  }
#endif
#ifdef ENABLE_DBUS
  if (Feature::ExperimentalInputDriverDBus.is_enabled()) {
    if (Settings::Settings::inputEnableDriverDBUS.value()) {
      auto dBusDriver = new DBusInputDriver();
      InputDriverManager::instance()->registerDriver(dBusDriver);
    }
  }
#endif

#ifdef ENABLE_GUI_TESTS
  // Adds a singleshot timer that will be executed when the application will be started.
  // the timer validates that each mainwindow respects the expected UX behavior.
  if (gui_test != "none") {
    QTimer::singleShot(0, [&]() {
      int failureCount = 0;
      for (auto w : app.windowManager.getWindows()) {
        failureCount += runAllTest(w);
      }
      app.exit(failureCount);
    });
  }
#endif  // ENABLE_GUI_TESTS

  InputDriverManager::instance()->init();
  return app.exec();
}
