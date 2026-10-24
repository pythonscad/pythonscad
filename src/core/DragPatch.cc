#include "core/DragPatch.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <locale>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

using Span = std::pair<size_t, size_t>;  // [first, second)

/// One text replacement [start, end) -> text in the source code.
struct DragEdit {
  size_t start;
  size_t end;
  std::string text;
};

bool isIdentChar(char c)
{
  return std::isalnum(static_cast<unsigned char>(c)) || c == '_';
}
bool isSpace(char c)
{
  return std::isspace(static_cast<unsigned char>(c));
}

void trimSpan(const std::string& src, size_t& start, size_t& end)
{
  while (start < end && isSpace(src[start])) start++;
  while (end > start && isSpace(src[end - 1])) end--;
}

/// Split the top level elements of a bracket starting at src[open] ('(' '[' '{').
/// Strings and python comments are skipped. close receives the matching bracket.
bool splitBracket(const std::string& src, size_t open, std::vector<Span>& elems, size_t& close)
{
  std::vector<char> stack;
  size_t elemstart = open + 1;
  for (size_t i = open; i < src.size(); i++) {
    const char c = src[i];
    if (c == '"' || c == '\'') {
      const size_t q = src.find(c, i + 1);
      if (q == std::string::npos) return false;
      i = q;
      continue;
    }
    if (c == '#') {  // python comment until end of line
      const size_t nl = src.find('\n', i);
      if (nl == std::string::npos) return false;
      i = nl;
      continue;
    }
    if (c == '(' || c == '[' || c == '{') {
      stack.push_back(c);
    } else if (c == ')' || c == ']' || c == '}') {
      if (stack.empty()) return false;
      stack.pop_back();
      if (stack.empty()) {
        size_t s = elemstart, e = i;
        trimSpan(src, s, e);
        if (s < e) elems.emplace_back(elemstart, i);  // skip empty trailing element
        close = i;
        return true;
      }
    } else if (c == ',' && stack.size() == 1) {
      elems.emplace_back(elemstart, i);
      elemstart = i + 1;
    }
  }
  return false;
}

/// If [start,end) is a plain variable name, replace the span by the right hand
/// side of the last "name = ..." assignment before limit (at the start of a line).
/// This lets e.g. polyhedron(pts, faces) patch the list assigned to pts.
void resolveVariable(const std::string& src, size_t limit, size_t& start, size_t& end)
{
  for (int hop = 0; hop < 8; hop++) {  // follow a = b chains, but never loop forever
    trimSpan(src, start, end);
    if (start >= end || !(std::isalpha(static_cast<unsigned char>(src[start])) || src[start] == '_'))
      return;
    for (size_t k = start; k < end; k++) {
      if (!isIdentChar(src[k])) return;  // not a bare identifier
    }
    const std::string name = src.substr(start, end - start);

    bool found = false;
    size_t pos = std::min(limit, start);
    while (pos > 0 && !found) {
      const size_t hit = src.rfind(name, pos - 1);
      if (hit == std::string::npos) break;
      pos = hit;
      // must start a line (after indentation) and be followed by a single '='
      size_t ls = hit;
      while (ls > 0 && (src[ls - 1] == ' ' || src[ls - 1] == '\t')) ls--;
      if (ls > 0 && src[ls - 1] != '\n') continue;
      size_t k = hit + name.size();
      if (k < src.size() && isIdentChar(src[k])) continue;
      while (k < src.size() && (src[k] == ' ' || src[k] == '\t')) k++;
      if (k >= src.size() || src[k] != '=' || (k + 1 < src.size() && src[k + 1] == '=')) continue;
      k++;
      while (k < src.size() && (src[k] == ' ' || src[k] == '\t')) k++;
      if (k >= src.size()) return;
      size_t rhs_end;
      if (src[k] == '[' || src[k] == '(') {
        std::vector<Span> dummy;
        size_t close;
        if (!splitBracket(src, k, dummy, close)) return;
        rhs_end = close + 1;
      } else {
        rhs_end = src.find_first_of("#\n", k);
        if (rhs_end == std::string::npos) rhs_end = src.size();
      }
      start = k;
      end = rhs_end;
      limit = hit;
      found = true;
    }
    if (!found) return;
  }
}

/// Inside [start,end) (a list literal) descend to the element given by path.
/// Variables are followed to their assignment before limit.
bool findListElement(const std::string& src, size_t limit, size_t& start, size_t& end,
                     const std::vector<int>& path)
{
  for (int index : path) {
    resolveVariable(src, limit, start, end);
    trimSpan(src, start, end);
    if (start >= end || (src[start] != '[' && src[start] != '(')) return false;
    std::vector<Span> elems;
    size_t close;
    if (!splitBracket(src, start, elems, close) || close >= end) return false;
    if (index < 0 || index >= static_cast<int>(elems.size())) return false;
    start = elems[index].first;
    end = elems[index].second;
  }
  resolveVariable(src, limit, start, end);
  trimSpan(src, start, end);
  return start < end;
}

/// If [start,end) is "name = value" with the given name, return the value span.
bool namedArgValue(const std::string& src, size_t start, size_t end, const std::string& name,
                   size_t& vstart, size_t& vend)
{
  trimSpan(src, start, end);
  if (end - start <= name.size() || src.compare(start, name.size(), name) != 0) return false;
  size_t k = start + name.size();
  if (isIdentChar(src[k])) return false;  // "sizes = " is not "size"
  while (k < end && isSpace(src[k])) k++;
  if (k >= end || src[k] != '=' || (k + 1 < end && src[k + 1] == '=')) return false;
  vstart = k + 1;
  vend = end;
  trimSpan(src, vstart, vend);
  return vstart < vend;
}

/// A draggable value is a quoted number ("10", '10') or a marked() call
/// (marked(10), marked(value=10), pythonscad.marked(10)).
/// On success [start,end) is narrowed to the bare number which gets replaced.
bool findMarkedNumber(const std::string& src, size_t& start, size_t& end)
{
  trimSpan(src, start, end);
  if (end - start < 2) return false;

  const char q = src[start];
  if ((q == '"' || q == '\'') && src[end - 1] == q) {
    start++;
    end--;
    return true;
  }

  // marked( ... )
  const size_t paren = src.find('(', start);
  if (paren == std::string::npos || paren >= end || src[end - 1] != ')') return false;
  size_t ns = start, ne = paren;
  trimSpan(src, ns, ne);
  const std::string name = src.substr(ns, ne - ns);
  static const std::string marked = "marked";
  const bool is_marked =
    name == marked || (name.size() > marked.size() &&
                       name.compare(name.size() - marked.size(), marked.size(), marked) == 0 &&
                       name[name.size() - marked.size() - 1] == '.');
  if (!is_marked) return false;

  std::vector<Span> args;
  size_t close;
  if (!splitBracket(src, paren, args, close) || close != end - 1 || args.size() != 1) return false;
  size_t vs, ve;
  if (namedArgValue(src, args[0].first, args[0].second, "value", vs, ve)) {
    start = vs;
    end = ve;
  } else {
    start = args[0].first;
    end = args[0].second;
    trimSpan(src, start, end);
  }
  return start < end;
}

/// Compute all edits for one call whose '(' is at src[open]; all mods must match.
bool editsForCall(const std::string& src, size_t open, bool method, const DragResult& result,
                  std::vector<DragEdit>& edits)
{
  std::vector<Span> args;
  size_t close;
  if (!splitBracket(src, open, args, close)) return false;

  std::vector<DragEdit> local;
  for (const auto& mod : result.mods) {
    size_t vstart = std::string::npos, vend = std::string::npos;
    for (const auto& arg : args) {
      if (namedArgValue(src, arg.first, arg.second, mod.name, vstart, vend)) break;
    }
    if (vstart == std::string::npos) {  // positional
      // obj.translate(v) has no obj argument, translate(obj, v) has
      const int index = mod.index - ((result.objFirst && method) ? 1 : 0);
      if (index < 0 || index >= static_cast<int>(args.size())) return false;
      vstart = args[index].first;
      vend = args[index].second;
    }
    if (!findListElement(src, open, vstart, vend, mod.arrinfo)) return false;
    if (!findMarkedNumber(src, vstart, vend)) return false;
    // two mods may resolve to the same text (e.g. a shared variable): patch once
    bool duplicate = false;
    for (const auto& e : local) duplicate |= e.start == vstart;
    if (!duplicate) local.push_back({vstart, vend, formatDragValue(mod.value)});
  }
  edits = std::move(local);
  return !edits.empty();
}

}  // namespace

std::string formatDragValue(double value)
{
  value = std::round(value * 1e4) / 1e4;
  if (value == 0) value = 0;  // no "-0"
  std::ostringstream ss;
  ss.imbue(std::locale::classic());
  ss << std::setprecision(10) << value;
  return ss.str();
}

bool applyDragResult(const std::string& src, const DragResult& result, std::string& out)
{
  if (result.modname.empty() || result.mods.empty()) return false;
  size_t pos = 0;
  while ((pos = src.find(result.modname, pos)) != std::string::npos) {
    const size_t after = pos + result.modname.size();
    const bool word_start = pos == 0 || !isIdentChar(src[pos - 1]);
    const bool word_end = after >= src.size() || !isIdentChar(src[after]);
    size_t open = after;
    while (open < src.size() && isSpace(src[open])) open++;
    pos = after;
    if (!word_start || !word_end || open >= src.size() || src[open] != '(') continue;

    size_t before = pos - result.modname.size();  // pos already points behind the name
    while (before > 0 && isSpace(src[before - 1])) before--;
    const bool method = before > 0 && src[before - 1] == '.';

    std::vector<DragEdit> edits;
    // "a.translate(" is a method call, but "pythonscad.translate(" is a module
    // qualified function: if the method reading does not fit, try the other one
    if (!editsForCall(src, open, method, result, edits) &&
        !(method && editsForCall(src, open, false, result, edits)))
      continue;
    // apply back to front so earlier offsets stay valid
    std::sort(edits.begin(), edits.end(),
              [](const DragEdit& a, const DragEdit& b) { return a.start > b.start; });
    out = src;
    for (const auto& e : edits) out.replace(e.start, e.end - e.start, e.text);
    return true;
  }
  return false;
}
