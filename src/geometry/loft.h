#pragma once

#include <functional>
#include <memory>
#include <vector>

#include "geometry/PolySet.h"
#include "geometry/linalg.h"

// Erzeugt eine Flaeche zwischen einem geschlossenen Aussenring ('outer') und
// 0..N geschlossenen Innenringen ('holes'), beliebige Punktzahl pro Ring.
//
// 'proj' bildet einen 3D-Punkt auf eine 2D-Koordinate (u,v) ab - nur zur
// Bestimmung der Topologie (welcher Bereich gehoert zur Flaeche, welcher zu
// einem Loch). Muss fuer die Gesamtheit aller uebergebenen Ringe lokal
// ueberschneidungsfrei sein (z.B. proj_xy fuer eine flache Kappe, oder eine
// zylindrische (theta,z)-Abwicklung fuer eine Wand).
//
// 'grid_spacing_uv' ist die Zielaufloesung des inneren Stuetzpunkt-Gitters
// (dieselbe Einheit wie proj()'s Ausgabe).
//
// 'displacement' bekommt die (unverschobene) 3D-Basisposition jedes
// erzeugten Punktes und liefert einen Versatz entlang der lokalen
// Flaechennormalen zurueck - fuer texturierte/verbeulte Oberflaechen. Fuer
// eine glatte Flaeche einfach eine Funktion uebergeben, die immer 0
// zurueckgibt.
std::unique_ptr<PolySet> loft(const std::vector<Vector3d>& outer,
                              const std::vector<std::vector<Vector3d>>& holes,
                              const std::function<Vector2d(const Vector3d&)>& proj,
                              double grid_spacing_uv,
                              const std::function<double(const Vector3d&)>& displacement);
