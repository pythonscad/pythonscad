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
  if (languageManuallySet) return;  // Don't override manual selection

  auto fnameba = filepath.toLocal8Bit();
  const char *fname = filepath.isEmpty() ? "" : fnameba;

  int oldLanguage = language;
  language = LANG_SCAD;
  if (fname != NULL) {
#ifdef ENABLE_PYTHON
    if (boost::algorithm::ends_with(fname, ".py")) {
      language = LANG_PYTHON;
    }
#endif
  }

#ifdef ENABLE_PYTHON
  if (oldLanguage != language) {
    onLanguageChanged(language);
  }
#endif
}

void EditorInterface::setLanguageManually(int lang)
{
  languageManuallySet = true;
  if (language != lang) {
    language = lang;
    onLanguageChanged(lang);
  }
}

void EditorInterface::resetLanguageDetection()
{
  languageManuallySet = false;
  recomputeLanguageActive();
}
