#pragma once

#include <Qsci/qsciapis.h>

#include <QList>
#include <QObject>
#include <QString>
#include <QStringList>
#include <unordered_set>
#include <utility>

#include "core/SourceFile.h"
class ApiFunc;
class ScintillaEditor;

class PythonApi : public QsciAbstractAPIs
{
  Q_OBJECT

private:
  ScintillaEditor *editor;
  QList<ApiFunc> funcs;

protected:
  void autoCompleteFolder(const QStringList& context, const QString& text, const int col,
                          QStringList& list);
  void autoCompleteFunctions(const QStringList& context, QStringList& list);

public:
  PythonApi(ScintillaEditor *editor, QsciLexer *lexer);

  void updateAutoCompletionList(const QStringList& context, QStringList& list) override;
  void autoCompletionSelected(const QString& selection) override;
  QStringList callTips(const QStringList& context, int commas, QsciScintilla::CallTipsStyle style,
                       QList<int>& shifts) override;

  void correctUserVarNamesForCompletionFromSourceFile(const SourceFile *sourceFile,
                                                      bool flagAutoCompleteIncludeVariables,
                                                      bool flagAutoCompleteIncludeModules,
                                                      bool flagAutoCompleteIncludeFunctions);

  void correctUserVarNamesForCompletionFromInputText(bool flagAutoCompleteIncludeVariables,
                                                     bool flagAutoCompleteIncludeModules,
                                                     bool flagAutoCompleteIncludeFunctions);

private:
  QStringList userVariableNames;
};
