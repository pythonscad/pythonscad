#pragma once

#include <string>
#include <vector>

#include "core/node.h"
#include "geometry/linalg.h"

class LoftNode : public LeafNode
{
public:
  LoftNode(std::shared_ptr<const ModuleInstantiation> mi) : LeafNode(std::move(mi)) {}
  LoftNode(const LoftNode& other);  // eigene Copy-Ctor: haelt proj_func/
                                    // displacement_func korrekt am Leben (Py_XINCREF)
  ~LoftNode() override;             // Py_XDECREF auf proj_func/displacement_func

  std::string toString() const override;
  std::string name() const override { return "loft"; }
  std::unique_ptr<const class Geometry> createGeometry() const override;

  std::vector<Vector3d> outer;
  std::vector<std::vector<Vector3d>> holes;
  double grid_spacing_uv = 1.0;

  // Python-Funktionsobjekte (PyObject*), opak gehalten, damit dieser Header
  // kein Python.h braucht. Werden in LoftNode.cc korrekt refcounted.
  void *proj_func = nullptr;
  void *displacement_func = nullptr;

  // Inhaltsbasierte Hashes der beiden Funktionen (siehe python_func_content_hash),
  // einmalig bei Node-Erzeugung berechnet - fuer toString()/Caching.
  std::string proj_func_hash;
  std::string displacement_func_hash;
};
