#include "gui/Editor.h"
#include <QWheelEvent>
#include <QWidget>
#include "gui/Preferences.h"
#include "gui/QSettingsCached.h"
#include "genlang/genlang.h"

void EditorInterface::wheelEvent(QWheelEvent *event)
{
  QSettingsCached settings;
  bool wheelzoom_enabled = GlobalPreferences::inst()->getValue("editor/ctrlmousewheelzoom").toBool();
  if ((event->modifiers() == Qt::ControlModifier) && wheelzoom_enabled) {
    if (event->angleDelta().y() > 0) zoomIn();
    else if (event->angleDelta().y() < 0) zoomOut();
  } else {
    QWidget::wheelEvent(event);
  }
}

void EditorInterface::recomputeLanguageActive()
{
  printf("RecomputeLang\n");
  auto fnameba = filepath.toLocal8Bit();
  const char *fname = filepath.isEmpty() ? "" : fnameba;

  int oldLanguage = language;
  language = LANG_SCAD;
  if (fname != NULL) {
#ifdef ENABLE_PYTHON
    if (boost::algorithm::ends_with(fname, ".py")) {
      std::string content = toPlainText().toStdString();
      /* if (trust_python_file(std::string(fname), content)) */ language = LANG_PYTHON;  // TODO activate
      // else LOG(message_group::Warning, Location::NONE, "", "Python is not enabled");
    }
#endif
  }

#ifdef ENABLE_PYTHON
  if (oldLanguage != language) {
    onLanguageActiveChanged(language);
  }
#endif
}
