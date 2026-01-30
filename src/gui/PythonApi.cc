#include "gui/PythonApi.h"

#include <QString>
#include <QStringList>
#include <string>

#include "core/Builtins.h"
#include "gui/ScintillaEditor.h"

namespace {

bool isInPythonString(const QString& text, int col)
{
  bool lastWasEscape = false;
  bool inSingle = false;
  bool inDouble = false;
  int dx = 0;
  int count = col;
  while (count-- > 0 && dx < text.length()) {
    QChar ch = text[dx++];
    if (ch == '\\') {
      lastWasEscape = true;
    } else if (lastWasEscape) {
      lastWasEscape = false;
    } else if (ch == '\'' && !inDouble) {
      inSingle = !inSingle;
    } else if (ch == '"' && !inSingle) {
      inDouble = !inDouble;
    }
  }
  return inSingle || inDouble;
}

}  // namespace

PythonApi::PythonApi(ScintillaEditor *editor, QsciLexerPython *lexer)
  : QsciAbstractAPIs(lexer), editor(editor)
{
  for (const auto& iter : Builtins::keywordList) {
    QStringList calltipList;
    for (const auto& it : iter.second) {
      calltipList.append(QString::fromStdString(it));
    }
    funcs.append({QString::fromStdString(iter.first), calltipList});
  }
  // Python-only openscad names with call tips
  funcs.append(PythonApiFunc{"show", {"show(obj)", "show(obj, ...)"}});
  funcs.append(PythonApiFunc{"add_parameter", {"add_parameter(name, default, ...)"}});
}

void PythonApi::updateAutoCompletionList(const QStringList& context, QStringList& list)
{
  int line, col;
  editor->qsci->getCursorPosition(&line, &col);
  const QString text = editor->qsci->text(line);

  if (isInPythonString(text, col)) {
    return;
  }
  autoCompleteFunctions(context, list);
}

void PythonApi::autoCompleteFunctions(const QStringList& context, QStringList& list)
{
  if (context.isEmpty()) {
    return;
  }
  const QString& c = context.last();
  if (c.isEmpty()) {
    return;
  }

  for (const auto& func : funcs) {
    if (func.name.startsWith(c) && !list.contains(func.name)) {
      list << func.name;
    }
  }
}

void PythonApi::autoCompletionSelected(const QString& /*selection*/)
{
}

QStringList PythonApi::callTips(const QStringList& context, int /*commas*/,
                                QsciScintilla::CallTipsStyle /*style*/, QList<int>& /*shifts*/)
{
  QStringList tips;
  if (context.size() < 2) {
    return tips;
  }
  const QString name = context.at(context.size() - 2);
  for (const auto& func : funcs) {
    if (func.name == name) {
      tips = func.params;
      break;
    }
  }
  return tips;
}
