#include "gui/PythonApi.h"

#include <QString>
#include <QStringList>
#include <string>

#include "core/Builtins.h"
#include "gui/ScintillaEditor.h"

namespace {

// Returns true if the cursor at col is inside a string (single, double, or triple-quoted).
bool isInPythonString(const QString& text, int col)
{
  bool lastWasEscape = false;
  enum State { Normal, Single, Double, TripleSingle, TripleDouble };
  State state = Normal;

  for (int dx = 0; dx < col && dx < text.length(); ++dx) {
    QChar ch = text[dx];

    if (state == Normal) {
      if (ch == '\'') {
        if (dx + 2 < text.length() && text[dx + 1] == '\'' && text[dx + 2] == '\'') {
          state = TripleSingle;
          dx += 2;
        } else {
          state = Single;
        }
      } else if (ch == '"') {
        if (dx + 2 < text.length() && text[dx + 1] == '"' && text[dx + 2] == '"') {
          state = TripleDouble;
          dx += 2;
        } else {
          state = Double;
        }
      }
      continue;
    }

    if (state == Single) {
      if (ch == '\\') {
        lastWasEscape = true;
      } else if (lastWasEscape) {
        lastWasEscape = false;
      } else if (ch == '\'') {
        state = Normal;
      }
      continue;
    }

    if (state == Double) {
      if (ch == '\\') {
        lastWasEscape = true;
      } else if (lastWasEscape) {
        lastWasEscape = false;
      } else if (ch == '"') {
        state = Normal;
      }
      continue;
    }

    if (state == TripleSingle) {
      if (ch == '\\') {
        lastWasEscape = true;
      } else if (lastWasEscape) {
        lastWasEscape = false;
      } else if (ch == '\'' && dx + 2 < text.length() && text[dx + 1] == '\'' && text[dx + 2] == '\'') {
        state = Normal;
        dx += 2;
      }
      continue;
    }

    if (state == TripleDouble) {
      if (ch == '\\') {
        lastWasEscape = true;
      } else if (lastWasEscape) {
        lastWasEscape = false;
      } else if (ch == '"' && dx + 2 < text.length() && text[dx + 1] == '"' && text[dx + 2] == '"') {
        state = Normal;
        dx += 2;
      }
      continue;
    }
  }

  return state != Normal;
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
