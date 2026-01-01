#pragma once

#include <vector>

#include "Traits.hpp"

// struct PointLocatorResult {
//   unsigned elem = UMAX;       // UMAX if not found
//   std::vector<double> xi;      // size dim
//   //double xi[3];
//   bool ok = false;
// };

class PointLocatorResult {
  public:
    PointLocatorResult() {
      elem = UMAX;
      xi.reserve(3);
      xi.clear();
      ok = false;
    }
    unsigned elem;       // UMAX if not found
    std::vector<double> xi;      // size dim
    bool ok;
};
