#pragma once

#include <QList>
#include <QString>
#include <QStringList>

#include <Qsci/qsciapis.h>
#include <Qsci/qscilexerpython.h>

class ScintillaEditor;

struct PythonApiFunc {
  QString name;
  QStringList params;
};

class PythonApi : public QsciAbstractAPIs
{
  Q_OBJECT

private:
  ScintillaEditor *editor;
  QList<PythonApiFunc> funcs;

  void autoCompleteFunctions(const QStringList& context, QStringList& list);

public:
  PythonApi(ScintillaEditor *editor, QsciLexerPython *lexer);

  void updateAutoCompletionList(const QStringList& context, QStringList& list) override;
  void autoCompletionSelected(const QString& selection) override;
  QStringList callTips(const QStringList& context, int commas, QsciScintilla::CallTipsStyle style,
                       QList<int>& shifts) override;
};
