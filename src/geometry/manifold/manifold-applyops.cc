// this file is split into many separate cgalutils* files
// in order to workaround gcc 4.9.1 crashing on systems with only 2GB of RAM

#ifdef ENABLE_MANIFOLD

#include <memory>
#include "geometry/manifold/manifoldutils.h"
#include "geometry/Geometry.h"
#include "core/AST.h"
#include "geometry/manifold/ManifoldGeometry.h"
#include "core/node.h"
#include "core/progress.h"
#include "utils/printutils.h"

namespace ManifoldUtils {

Location getLocation(const std::shared_ptr<const AbstractNode>& node)
{
  return node && node->modinst ? node->modinst->location() : Location::NONE;
}

/*!
   Applies op to all children and returns the result.
   The child list should be guaranteed to contain non-NULL 3D or empty Geometry objects
 */

int enable_negspace=1;

#define NEGSPACE_CLONE \
      arg = chN->copy1();\
      argneg = std::move(arg->neg_space);\
      if(argneg == nullptr) argneg = std::make_shared<ManifoldGeometry>();\

std::shared_ptr<ManifoldGeometry> applyOperator3DManifoldNegSpace(const Geometry::Geometries& children,
                                                          OpenSCADOperator op)
{
  std::shared_ptr<ManifoldGeometry> resultpp = nullptr, resultpn = nullptr, resultnp = nullptr, resultnn = nullptr;
  std::shared_ptr<ManifoldGeometry> firstp = nullptr, firstn = nullptr, arg, argneg, result;

  for (const auto& item : children) {
    auto chN = item.second ? createManifoldFromGeometry(item.second) : nullptr;

    // Intersecting something with nothing results in nothing
    if (!chN || chN->isEmpty()) {
      if (op == OpenSCADOperator::INTERSECTION) {
        break;
      }
      if (op == OpenSCADOperator::DIFFERENCE && resultpp == nullptr) {
        break;
      }
      continue;
    }

    // Initialize geom with first expected geometric object

    NEGSPACE_CLONE

    if (resultpp == nullptr) {
      if(op == OpenSCADOperator::INTERSECTION) {
        resultpp = arg;
        resultpn = argneg;
        NEGSPACE_CLONE
        resultnp = arg;
        resultnn = argneg;
      } else {	      
        resultpp = arg;
        resultnp = argneg;
        NEGSPACE_CLONE
        firstp = arg;
        firstn = argneg;
      }	
      continue;
    }

    switch (op) {
    case OpenSCADOperator::UNION:        
      *resultpp = *resultpp - *argneg;
      *resultnp = *resultnp - *arg;
      NEGSPACE_CLONE
      if(resultpn != nullptr) { 
        *resultpn = *resultpn + *arg;
        *resultnn = *resultnn + *argneg;
      } else {
        resultpn = arg;
        resultnn = argneg;
      }
      break;
    case OpenSCADOperator::INTERSECTION:        
      *resultpp = *resultpp * *arg;
      *resultpn = *resultpn * *argneg;
      *resultnp = *resultnp + *argneg;
      *resultnn = *resultnn + *arg;
      break;
    case OpenSCADOperator::DIFFERENCE:   
      *resultpp = *resultpp - *arg;
      *resultnp = *resultnp - *argneg;
      NEGSPACE_CLONE
      if(resultpn != nullptr) { 
        *resultpn = *resultpn + *argneg;
        *resultnn = *resultnn + *arg;
      } else {
        resultpn = argneg;
        resultnn = arg;
      }

	    break;
//    case OpenSCADOperator::MINKOWSKI:    *geom = geom->minkowski(*chN); break;
    default:                             LOG(message_group::Error, "Unsupported CGAL operator: %1$d", static_cast<int>(op));
    }
    if (item.first) item.first->progress_report();
  }
  switch(op){
    case OpenSCADOperator::UNION:   
    case OpenSCADOperator::DIFFERENCE:   
      if(resultpn != nullptr) {
        *resultpn = *resultpn - *firstn; 
        *resultnn = *resultnn - *firstp;      
      }	
      break;	    
  }
  // TODO switch off
  // DOCU A: arg1pos, B: arg1neg C: arg2pos D: arg2neg
  // UNION: POS = C-B | A-D NEG = D-A  | B-C
  // DIFF   POS = A-C | D-B NEG = B-D  | C-A
  // INTER  POS =A&C  | B&D NEG = A&D  | B*C
  result = resultpp;
  if(resultpn != nullptr) *result = *result + *resultpn; 

  result->neg_space = resultnp;

  if(resultnn != nullptr) *(result->neg_space) = *(result->neg_space)  + *resultnn; 
  return result;
}


std::shared_ptr<ManifoldGeometry> applyOperator3DManifold(const Geometry::Geometries& children,
                                                          OpenSCADOperator op)
{
  if(enable_negspace)  return applyOperator3DManifoldNegSpace(children, op);
  std::shared_ptr<ManifoldGeometry> geom;

  bool foundFirst = false;

  for (const auto& item : children) {
    auto chN = item.second ? createManifoldFromGeometry(item.second) : nullptr;

    // Intersecting something with nothing results in nothing
    if (!chN || chN->isEmpty()) {
      if (op == OpenSCADOperator::INTERSECTION) {
        geom = nullptr;
        break;
      }
      if (op == OpenSCADOperator::DIFFERENCE && !foundFirst) {
        geom = nullptr;
        break;
      }
      continue;
    }

    // Initialize geom with first expected geometric object
    if (!foundFirst) {
      geom = std::make_shared<ManifoldGeometry>(*chN);
      foundFirst = true;
      continue;
    }

    switch (op) {
    case OpenSCADOperator::UNION:        *geom = *geom + *chN; break;
    case OpenSCADOperator::INTERSECTION: *geom = *geom * *chN; break;
    case OpenSCADOperator::DIFFERENCE:   *geom = *geom - *chN; break;
    case OpenSCADOperator::MINKOWSKI:    *geom = geom->minkowski(*chN); break;
    default:                             LOG(message_group::Error, "Unsupported CGAL operator: %1$d", static_cast<int>(op));
    }
    if (item.first) item.first->progress_report();
  }
  return geom;
}

};  // namespace ManifoldUtils

#endif  // ENABLE_MANIFOLD
