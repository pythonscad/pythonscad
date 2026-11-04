#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include "geometry/Geometry.h"
#include "geometry/PolySet.h"
#include "geometry/linalg.h"

// Select check: spatial relationship filtering between two solid sets.
// A (the parts to filter) is separated into individual pieces.
// Each piece is tested against B (the reference solid) for spatial relations:
//   - inside: piece is completely inside B
//   - outside: piece is completely outside B
//   - straddle: piece crosses B's boundary (neither completely inside nor outside)
//
// The result is a union of parts that match (or don't match) the specified relation.
namespace SelectCheck {

enum class Relation {
  Inside,      // piece completely inside B
  NotInside,   // piece NOT completely inside B
  Outside,     // piece completely outside B
  NotOutside,  // piece NOT completely outside B
  Straddle,    // piece crosses B's boundary
  NotStraddle  // piece does NOT cross B's boundary
};

struct Options {
  Relation relation = Relation::Inside;
  unsigned threads = 0;  // 0 = hardware concurrency
};

struct PartInfo {
  int part_index = -1;
  bool matches = false;  // true if part matches the relation criterion
};

struct Result {
  std::vector<PartInfo> parts;                  // individual part results
  std::shared_ptr<const Geometry> error_solid;  // union of parts that don't match
  size_t count = 0;                             // number of parts that don't match the relation
  bool clean() const { return count == 0; }
};

// Perform select check: separate A into parts, test each against B for the relation
// Returns Result with error_solid containing parts that violate the relation
Result check(const PolySet& A, const PolySet& B, const Options& opt);

}  // namespace SelectCheck
