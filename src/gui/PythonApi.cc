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
  funcs.append(PythonApiFunc{"right", {"right(x)", "right(obj, x)"}});
  funcs.append(PythonApiFunc{"left", {"left(x)", "left(obj, x)"}});
  funcs.append(PythonApiFunc{"back", {"back(x)", "back(obj, x)"}});
  funcs.append(PythonApiFunc{"front", {"front(x)", "front(obj, x)"}});
  funcs.append(PythonApiFunc{"up", {"up(x)", "up(obj, x)"}});
  funcs.append(PythonApiFunc{"down", {"down(x)", "down(obj, x)"}});
  funcs.append(PythonApiFunc{"rotx", {"rotx(angle)", "rotx(obj, angle)"}});
  funcs.append(PythonApiFunc{"roty", {"roty(angle)", "roty(obj, angle)"}});
  funcs.append(PythonApiFunc{"rotz", {"rotz(angle)", "rotz(obj, angle)"}});
  funcs.append(PythonApiFunc{"separate", {"separate(obj)"}});
  funcs.append(PythonApiFunc{"export", {"export(obj)", "export(obj, file=...)"}});
  funcs.append(PythonApiFunc{"find_face", {"find_face(obj, normal)"}});
  funcs.append(PythonApiFunc{"sitonto", {"sitonto(obj, x, y, z)"}});
  funcs.append(PythonApiFunc{"path_extrude", {"path_extrude(shape, path, ...)"}});
  funcs.append(PythonApiFunc{"skin", {"skin(profiles, ...)"}});
  funcs.append(PythonApiFunc{"concat", {"concat(*objects)", "concat(*objects, r=..., fn=...)"}});
  funcs.append(PythonApiFunc{"highlight", {"highlight(obj)"}});
  funcs.append(PythonApiFunc{"background", {"background(obj)"}});
  funcs.append(PythonApiFunc{"only", {"only(obj)"}});
  funcs.append(PythonApiFunc{"sheet", {"sheet(file, ...)"}});
  funcs.append(PythonApiFunc{"inside", {"inside(obj, point)"}});
  funcs.append(PythonApiFunc{"bbox", {"bbox(obj)"}});
  funcs.append(PythonApiFunc{"size", {"size(obj)"}});
  funcs.append(PythonApiFunc{"position", {"position(obj)"}});
  funcs.append(PythonApiFunc{"faces", {"faces(obj)"}});
  funcs.append(PythonApiFunc{"children", {"children(obj)"}});
  funcs.append(PythonApiFunc{"edges", {"edges(obj)"}});
  funcs.append(PythonApiFunc{"explode", {"explode(obj, vector)"}});
  funcs.append(PythonApiFunc{"oversample", {"oversample(obj, n, round=...)"}});
  funcs.append(PythonApiFunc{"debug", {"debug(obj, faces=...)"}});
  funcs.append(PythonApiFunc{"repair", {"repair(obj)"}});
  funcs.append(PythonApiFunc{"fillet", {"fillet(obj, r, sel, fn=...)"}});
  funcs.append(PythonApiFunc{"group", {"group(*objects)"}});
  funcs.append(PythonApiFunc{"osimport", {"osimport(file, ...)"}});
  funcs.append(PythonApiFunc{"osuse", {"osuse(path)"}});
  funcs.append(PythonApiFunc{"osinclude", {"osinclude(path)"}});
  funcs.append(PythonApiFunc{"version", {"version()"}});
  funcs.append(PythonApiFunc{"version_num", {"version_num()"}});
  funcs.append(PythonApiFunc{"scad", {"scad(code)"}});
  funcs.append(PythonApiFunc{"add_menuitem", {"add_menuitem(menuname, itemname, callback)"}});
  funcs.append(PythonApiFunc{"nimport", {"nimport(url)"}});
  funcs.append(PythonApiFunc{"model", {"model()"}});
  funcs.append(PythonApiFunc{"modelpath", {"modelpath()"}});
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
