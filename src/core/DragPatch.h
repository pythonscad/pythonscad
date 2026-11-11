#pragma once

#include <string>

#include "core/DragResult.h"

/*!
   Write the values of an interactive drag back into the script source.

   The call named result.modname (e.g. "cube") is searched in src. For every
   DragMod the addressed argument (named or positional, optionally an element
   of a nested list) must be marked as draggable, i.e. written as
     "10"   '10'   marked(10)   marked(value=10)   pythonscad.marked(10)
   Only the number inside the marker is replaced. All mods of one call are
   applied together. If a call does not match all mods, the next occurrence
   of the name is tried.

   Returns true and sets out to the patched source on success.
 */
bool applyDragResult(const std::string& src, const DragResult& result, std::string& out);

/*! Format a dragged value without floating point noise (12.300000000000001 -> "12.3"). */
std::string formatDragValue(double value);
