#pragma once

#include <cstdint>

/*
 * This collects the settings enums used in GUI code. None of this must
 * depend on any GUI or Qt classes.
 */

// Where "Export as ..." starts its file dialog.
enum class ExportLocationMode : std::uint8_t {
  lastUsed,       // the folder used for this file type last, forgotten on exit
  nextToDesign,   // always beside the .scad / .py being edited
  perProject,     // whatever this particular design exported to last, remembered
  fixedFolder,    // one configured folder for everything
};

// Folder name for the optional dated subfolder. isoYearFirst is the default
// because it is the only one of the three that sorts chronologically.
enum class ExportDateFormat : std::uint8_t {
  isoYearFirst,
  dayFirst,
  monthFirst,
};

enum class ColorListFilterType : std::uint8_t {
  fixed,
  wildcard,
  regexp,
};

enum class ColorListSortType : std::uint8_t {
  alphabetical,
  by_color,
  by_color_warmth,
  by_lightness,
};
