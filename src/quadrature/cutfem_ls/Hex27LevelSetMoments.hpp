#ifndef FEMUS_HEX27_LEVEL_SET_MOMENTS_HPP
#define FEMUS_HEX27_LEVEL_SET_MOMENTS_HPP

#include <string>
#include "FELevelSetMoments.hpp"

namespace FELevelSetMoments {

inline unsigned Hex27Index(unsigned i, unsigned j, unsigned k) {
  return 9U * i + 3U * j + k;
}

inline unsigned HexMomentIndex(unsigned i, unsigned j, unsigned k) {
  unsigned q = i + j + k;
  return q * (q + 1U) * (q + 2U) / 6U + (q - i) * (q - i + 1U) / 2U + k;
}

inline void Hex27OtherAxes(unsigned axis, unsigned& first, unsigned& second) {
  first = axis == 0 ? 1 : 0;
  second = axis == 2 ? 1 : 2;
}

inline unsigned Hex27PopCount(unsigned value) {
  unsigned count = 0;
  while(value != 0) { count += value & 1U; value >>= 1U; }
  return count;
}

inline void Hex27TrimBivariate(std::vector<std::vector<long double>>& value) {
  for(unsigned i = 0; i < value.size(); ++i) Trim(value[i]);
  while(value.size() > 1 && PolynomialNorm(value.back()) == 0.0L) value.pop_back();
}

inline std::vector<long double> Hex27Resultant(std::vector<std::vector<long double>> first, std::vector<std::vector<long double>> second) {
  Hex27TrimBivariate(first);
  Hex27TrimBivariate(second);
  unsigned m = static_cast<unsigned>(first.size()) - 1U;
  unsigned n = static_cast<unsigned>(second.size()) - 1U;
  if(m == 0 || n == 0) { std::vector<long double> result(1,1); for(unsigned i = 0; i < (m == 0 ? n : m); ++i) result = PolynomialProduct(result,m == 0 ? first[0] : second[0]); return result; }
  unsigned size = m + n;
  std::vector<std::vector<std::vector<long double>>> matrix(size, std::vector<std::vector<long double>>(size, std::vector<long double>(1, 0.0L)));
  for(unsigned row = 0; row < n; ++row) {
    for(unsigned j = 0; j <= m; ++j) matrix[row][row+j] = first[m-j];
  }
  for(unsigned row = 0; row < m; ++row) {
    for(unsigned j = 0; j <= n; ++j) matrix[n+row][row+j] = second[n-j];
  }
  std::vector<std::vector<long double>> determinant(1U << size, std::vector<long double>(1,0.0L));
  determinant[0][0] = 1.0L;
  for(unsigned mask = 0; mask + 1U < (1U << size); ++mask) {
    unsigned row = Hex27PopCount(mask);
    for(unsigned column = 0; column < size; ++column) {
      if((mask & (1U << column)) == 0 && PolynomialNorm(matrix[row][column]) != 0.0L) {
        long double sign = (Hex27PopCount(mask >> (column+1U)) & 1U) == 0 ? 1.0L : -1.0L;
        unsigned next = mask | (1U << column);
        determinant[next] = PolynomialSum(determinant[next], PolynomialProduct(determinant[mask],matrix[row][column]), sign);
      }
    }
  }
  return determinant.back();
}

}

#endif
