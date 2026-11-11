#include "TestTabManager.h"

#include <QString>
#include <QStringList>
#include <QTest>
#include <QPointer>
#include <QSignalSpy>
#include <QTemporaryDir>
#include <QJsonDocument>
#include <QJsonArray>
#include <QJsonObject>
#include <QFile>
#include <QComboBox>
#include <QPushButton>
#include "gui/QWordSearchField.h"
#include "gui/ScintillaEditor.h"

#include "platform/PlatformUtils.h"

void TestTabManager::initTestCase()
{
}

void TestTabManager::checkOpenClose()
{
  // The window has only one editor with file default.scad
  restoreWindowInitialState();

  QString filename = fixturePath("basic-ux/empty.scad");
  QString filename2 = fixturePath("basic-ux/empty2.scad");

  window->tabManager->open(filename);
  // The active editor must have a filepath equal to the loaded file
  QCOMPARE(window->activeEditor->filepath, filename);

  window->tabManager->open(filename2);
  // The active editor must have a filepath equal to the loaded file
  QCOMPARE(window->activeEditor->filepath, filename2);

  // Close empty2.scad
  window->tabManager->closeCurrentTab();

  // Close empty.scad
  window->tabManager->closeCurrentTab();

  // Only default.scad remain.
  QCOMPARE(window->tabManager->count(), 1);
}

void TestTabManager::checkReOpen()
{
  restoreWindowInitialState();

  QString filename = fixturePath("basic-ux/empty.scad");
  auto numPanel = window->tabManager->count();

  // When we open a new file,
  window->tabManager->open(filename);
  QCOMPARE(numPanel + 1, window->tabManager->count());

  // When we re-open a new file, nothing should happens as the file is already there
  window->tabManager->open(filename);
  QCOMPARE(numPanel + 1, window->tabManager->count());

  // After the tests we close the current tab.
  window->tabManager->closeCurrentTab();
}

void TestTabManager::checkCloseDuringCompilation()
{
  restoreWindowInitialState();
  QTRY_VERIFY(!GuiLocker::isLocked());
  window->designActionAutoReload->setChecked(false);
  window->tabManager->open(fixturePath("basic-ux/empty.scad"));
  QTRY_VERIFY(!GuiLocker::isLocked());

  QPointer<EditorInterface> editor = window->activeEditor;
  QPointer<QWidget> parameters = editor->parameterWidget;
  const auto count = window->tabManager->count();
  QSignalSpy closing(window->tabManager, &TabManager::editorAboutToClose);
  {
    // The same lock spans parsing, event processing, and asynchronous rendering.
    const GuiLocker compiling;
    window->tabManager->closeCurrentTab();
    window->tabManager->closeCurrentTab();  // Duplicate requests must close only once.
    QCOMPARE(window->tabManager->count(), count);
    QVERIFY(editor);
    QVERIFY(parameters);
    QCOMPARE(window->activeEditor, editor.data());
    QCOMPARE(closing.count(), 0);
  }

  QTRY_COMPARE(window->tabManager->count(), count - 1);
  QVERIFY(editor.isNull());
  QVERIFY(parameters.isNull());
  QCOMPARE(closing.count(), 1);
}

void TestTabManager::checkCancelRenderingTab_data()
{
  QTest::addColumn<bool>("fullRender");
  QTest::addColumn<bool>("closeSource");
  QTest::newRow("F5-source-tab") << false << true;
  QTest::newRow("F6-source-tab") << true << true;
  QTest::newRow("F5-other-tab") << false << false;
  QTest::newRow("F6-other-tab") << true << false;
}

void TestTabManager::checkCancelRenderingTab()
{
  QFETCH(bool, fullRender);
  QFETCH(bool, closeSource);
  restoreWindowInitialState();
  QTRY_VERIFY(!GuiLocker::isLocked());
  window->designActionAutoReload->setChecked(false);

  auto *otherEditor = window->activeEditor;

  QTemporaryDir directory;
  QVERIFY(directory.isValid());
  QFile file(directory.filePath("cancel-render.scad"));
  QVERIFY(file.open(QIODevice::WriteOnly));
  file.write(
    "render() difference() { sphere(r=10, $fn=48); "
    "translate([0,0,5]) sphere(r=9, $fn=48); }\n");
  file.close();
  window->tabManager->open(file.fileName());
  auto *sourceEditor = window->activeEditor;
  QPointer<EditorInterface> editor = closeSource ? sourceEditor : otherEditor;
  const auto count = window->tabManager->count();
  QSignalSpy finished(window, &MainWindow::compilationDone);
  bool discarded = false;
  QObject observer;
  connect(window, &MainWindow::compilationDone, &observer, [&]() {
    // Observe completion before the queued tab close clears the renderers.
    discarded = !window->rootGeom && !window->thrownTogetherRenderer;
#ifdef ENABLE_OPENCSG
    discarded = discarded && !window->previewRenderer;
#endif
  });

  if (fullRender) window->designActionRender->trigger();
  else window->actionRenderPreview();
  QVERIFY(GuiLocker::isLocked());
  const auto root = window->tree.root();
  QVERIFY(root);

  if (!closeSource) window->tabManager->switchToEditor(otherEditor);
  window->tabManager->closeCurrentTab();
  // This used to delete the editor and clear the worker's tree immediately.
  QVERIFY(editor);
  QCOMPARE(window->tabManager->count(), count);
  QCOMPARE(window->tree.root(), root);

  QTRY_VERIFY_WITH_TIMEOUT(!GuiLocker::isLocked(), 30000);
  QVERIFY(!finished.isEmpty());
  QCOMPARE(discarded, closeSource);
  QTRY_VERIFY(editor.isNull());
  QCOMPARE(window->tabManager->count(), count - 1);
  if (closeSource) {
    QVERIFY(!window->tree.root());
    QVERIFY(!window->rootGeom);
    QVERIFY(!window->thrownTogetherRenderer);
#ifdef ENABLE_OPENCSG
    QVERIFY(!window->previewRenderer);
#endif
  } else {
    QVERIFY(window->tree.root());
    window->tabManager->open(fixturePath("basic-ux/default.scad"));
    window->tabManager->switchToEditor(sourceEditor);
    window->tabManager->closeCurrentTab();
  }
}

void TestTabManager::checkCloseLastRenderingTab()
{
  QTRY_VERIFY(!GuiLocker::isLocked());
  QTemporaryDir directory;
  QVERIFY(directory.isValid());
  QFile file(directory.filePath("last-tab.scad"));
  QVERIFY(file.open(QIODevice::WriteOnly));
  file.write("render() sphere(r=10, $fn=48);\n");
  file.close();

  QPointer<MainWindow> renderWindow = new MainWindow(QStringList{file.fileName()});
  renderWindow->setAttribute(Qt::WA_DeleteOnClose);
  renderWindow->designActionAutoReload->setChecked(false);
  QTRY_VERIFY(!GuiLocker::isLocked());
  QCOMPARE(renderWindow->tabManager->count(), size_t(1));
  QPointer<EditorInterface> editor = renderWindow->activeEditor;

  renderWindow->actionRenderPreview();
  QVERIFY(GuiLocker::isLocked());
  const auto root = renderWindow->tree.root();
  QVERIFY(root);
  renderWindow->tabManager->closeCurrentTab();
  QVERIFY(renderWindow);
  QVERIFY(editor);
  QCOMPARE(renderWindow->tree.root(), root);

  QTRY_VERIFY_WITH_TIMEOUT(renderWindow.isNull(), 30000);
  QVERIFY(editor.isNull());
  QVERIFY(!GuiLocker::isLocked());
}

void TestTabManager::checkWindowWideSearch()
{
  restoreWindowInitialState();
  window->designActionAutoReload->setChecked(false);
  window->tabManager->open(fixturePath("basic-ux/empty.scad"));
  auto *a = window->activeEditor;
  window->tabManager->open(fixturePath("basic-ux/empty2.scad"));
  auto *b = window->activeEditor;
  a->setPlainText("rainbow rainbow");
  b->setPlainText("rainbow");
  window->showFind(true);
  window->findInputField->setText("rainbow");
  window->replaceInputField->setText("flag");
  // A tab switch must not seed the query from either editor's selection.
  static_cast<ScintillaEditor *>(a)->qsci->setSelection(0, 0, 0, 3);
  static_cast<ScintillaEditor *>(b)->qsci->setSelection(0, 0, 0, 4);
  for (auto *editor : {a, b, a}) {
    window->tabManager->switchToEditor(editor);
    QVERIFY(!window->find_panel->isHidden());
    QCOMPARE(window->findTypeComboBox->currentIndex(), 1);
    QCOMPARE(window->findInputField->text(), QString("rainbow"));
    QCOMPARE(window->replaceInputField->text(), QString("flag"));
    QCOMPARE(window->findInputField->findCount(), editor == a ? 2 : 1);
  }
  QTest::mouseClick(window->findDoneButton, Qt::LeftButton);
  for (auto *editor : {b, a, b}) {
    window->tabManager->switchToEditor(editor);
    QVERIFY(window->find_panel->isHidden());
    QCOMPARE(window->findInputField->findCount(), 0);
  }
  window->showFind(false);
  // Explicit Find still seeds from the current editor's selected text.
  QCOMPARE(window->findInputField->text(), QString("rain"));
  window->tabManager->switchToEditor(a);
  QVERIFY(!window->find_panel->isHidden());
  QCOMPARE(window->findTypeComboBox->currentIndex(), 0);
  QTest::keyClick(window->find_panel, Qt::Key_Escape);
  window->tabManager->switchToEditor(b);
  QVERIFY(window->find_panel->isHidden());
  a->setContentModified(false);
  b->setContentModified(false);
}

void TestTabManager::checkSearchSession_data()
{
  QTest::addColumn<int>("state");
  QTest::addColumn<bool>("legacy");
  for (int state = TabManager::FIND_HIDDEN; state <= TabManager::FIND_REPLACE_VISIBLE; ++state) {
    QTest::newRow(qPrintable(QString("state-%1").arg(state))) << state << false;
    QTest::newRow(qPrintable(QString("legacy-%1").arg(state))) << state << true;
  }
}

void TestTabManager::checkSearchSession()
{
  QFETCH(int, state);
  QFETCH(bool, legacy);
  restoreWindowInitialState();
  window->designActionAutoReload->setChecked(false);
  window->tabManager->open(fixturePath("basic-ux/empty.scad"));
  window->findInputField->setText("rainbow");
  window->replaceInputField->setText("flag");
  window->setFindState(state);
  QTemporaryDir dir;
  QVERIFY(dir.isValid());
  const QString path = dir.filePath("session.json");
  window->tabManager->saveSession(path);
  QFile file(path);
  QVERIFY(file.open(QIODevice::ReadOnly));
  auto root = QJsonDocument::fromJson(file.readAll()).object();
  file.close();
  auto windows = root["windows"].toArray();
  auto win = windows[0].toObject();
  auto panel = win["findPanel"].toObject();
  QCOMPARE(panel["state"].toInt(-1), state);
  if (legacy) {
    panel.remove("state");
    win["findPanel"] = panel;
    auto tabs = win["tabs"].toArray();
    const int index = win["currentIndex"].toInt();
    // An inactive tab's old state must not influence the migrated window state.
    auto inactive = tabs[0].toObject();
    inactive["findState"] = (state + 1) % 3;
    tabs[0] = inactive;
    auto tab = tabs[index].toObject();
    tab["findState"] = state;
    tabs[index] = tab;
    win["tabs"] = tabs;
    windows[0] = win;
    root["windows"] = windows;
    QVERIFY(file.open(QIODevice::WriteOnly | QIODevice::Truncate));
    file.write(QJsonDocument(root).toJson());
    file.close();
  }
  window->setFindState(TabManager::FIND_HIDDEN);
  window->findInputField->clear();
  window->replaceInputField->clear();
  window->tabManager->closeCurrentTab();
  QVERIFY(window->tabManager->restoreSession(path));
  QCOMPARE(window->findState, state);
  QCOMPARE(window->find_panel->isHidden(), state == TabManager::FIND_HIDDEN);
  QCOMPARE(window->findInputField->text(), QString("rainbow"));
  QCOMPARE(window->replaceInputField->text(), QString("flag"));
  window->tabManager->prevTab();
  QCOMPARE(window->findState, state);
  QCOMPARE(window->find_panel->isHidden(), state == TabManager::FIND_HIDDEN);
  window->setFindState(TabManager::FIND_HIDDEN);
}
