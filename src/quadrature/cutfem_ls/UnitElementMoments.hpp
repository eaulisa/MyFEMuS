#ifndef FEMUS_UNIT_ELEMENT_MOMENTS_HPP
#define FEMUS_UNIT_ELEMENT_MOMENTS_HPP

#include <cstdlib>
#include "Hex27LevelSetMoments.hpp"

namespace FELevelSetMoments {

inline unsigned UnitMomentIndex(unsigned dimension, unsigned i, unsigned j = 0, unsigned k = 0) {
  if(dimension == 1) return i;
  if(dimension == 2) return (i+j)*(i+j+1U)/2U+j;
  return HexMomentIndex(i,j,k);
}

}

#endif
