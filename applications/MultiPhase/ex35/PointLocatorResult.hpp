#pragma once

#include <vector>

#include "Traits.hpp"

struct PointLocatorResult {
  unsigned elem = UMAX;       // UMAX if not found
  std::vector<double> xi;      // size dim
  bool ok = false;
};

