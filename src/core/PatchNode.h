#pragma once

#include <string>
#include <vector>

#include "core/node.h"
#include "geometry/linalg.h"

class PatchNode : public LeafNode
{
public:
  PatchNode(std::shared_ptr<const ModuleInstantiation> mi) : LeafNode(std::move(mi)) {}
  PatchNode(const PatchNode& other);  // eigene Copy-Ctor: haelt proj_func/
                                      // displacement_func korrekt am Leben (Py_XINCREF)
  ~PatchNode() override;              // Py_XDECREF auf proj_func/displacement_func

  std::string toString() const override;
  std::string name() const override { return "loft"; }
  std::unique_ptr<const class Geometry> createGeometry() const override;

  std::vector<Vector3d> outer;
  std::vector<std::vector<Vector3d>> holes;
  double grid_spacing_uv = 1.0;

  // Optionale Tangenten-/Verlassrichtung je Randpunkt (parallel zu 'outer'
  // bzw. jedem Eintrag von 'holes'). Wird von python_loft_ring_from_shape()
  // aus der Ebenennormale eines 2D-Shapes befuellt; bei einer reinen
  // Punktliste bleiben die Arrays leer. Leer = altes rein lineares
  // Verhalten, siehe geometry/loft.h.
  std::vector<Vector3d> outer_normal;
  std::vector<std::vector<Vector3d>> holes_normal;

  // Ob die obigen Tangenten ueberhaupt angefordert wurden (Python:
  // use_tangents=True). MUSS Teil von toString() sein: 'outer'/'holes'
  // (die reinen 3D-Positionen) koennen fuer zwei Aufrufe identisch sein,
  // waehrend outer_normal/holes_normal sich unterscheiden (einmal mit,
  // einmal ohne Kruemmung angefordert) - ohne dieses Flag im Cache-Key
  // wuerden solche zwei PatchNodes denselben Geometrie-Cache-Eintrag
  // teilen und der falsche (linear statt gekruemmt, oder umgekehrt)
  // koennte zurueckgegeben werden.
  bool use_tangents = false;

  // Python-Funktionsobjekte (PyObject*), opak gehalten, damit dieser Header
  // kein Python.h braucht. Werden in PatchNode.cc korrekt refcounted.
  void *proj_func = nullptr;
  void *displacement_func = nullptr;

  // Inhaltsbasierte Hashes der beiden Funktionen (siehe python_func_content_hash),
  // einmalig bei Node-Erzeugung berechnet - fuer toString()/Caching.
  std::string proj_func_hash;
  std::string displacement_func_hash;
};
