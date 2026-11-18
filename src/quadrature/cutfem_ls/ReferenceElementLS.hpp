#ifndef FEMUS_REFERENCE_ELEMENT_LS_HPP
#define FEMUS_REFERENCE_ELEMENT_LS_HPP

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace FELevelSetMoments {

inline long double ReferencePower(long double value, unsigned degree) {
  long double result = 1.0L;
  for(unsigned i = 0; i < degree; ++i) result *= value;
  return result;
}

inline long double ReferenceBinomial(unsigned degree, unsigned index) {
  if(index > degree) return 0.0L;
  index = std::min(index,degree-index);
  long double result = 1.0L;
  for(unsigned i = 1; i <= index; ++i) result *= static_cast<long double>(degree-index+i)/i;
  return result;
}

inline unsigned ReferenceMomentIndex2D(unsigned degree, unsigned first, unsigned second) {
  return first*(2U*degree+3U-first)/2U+second;
}

inline long double ReferenceTriangleMoment(unsigned first, unsigned second) {
  long double result = 1.0L/((first+second+1.0L)*(first+second+2.0L));
  for(unsigned i = 1; i <= second; ++i) result *= static_cast<long double>(i)/(first+i);
  return result;
}

}

#endif
