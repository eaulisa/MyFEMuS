#ifndef FEMUS_HEX27_DIRECTION_SELECTION_HPP
#define FEMUS_HEX27_DIRECTION_SELECTION_HPP

#include "Hex27LevelSetMoments.hpp"

namespace FELevelSetMoments {

inline void Hex27SelectedAxes(unsigned order, unsigned& inner, unsigned& middle, unsigned& outer) {
  if(order == 0 || order > 6) throw std::invalid_argument("Invalid explicit HEX27 integration order.");
  const unsigned axes[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
  inner = axes[order-1U][0];
  middle = axes[order-1U][1];
  outer = axes[order-1U][2];
}

}

#endif
