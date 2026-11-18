#ifndef FEMUS_HEX27_LEVEL_SET_MOMENTS_HPP
#define FEMUS_HEX27_LEVEL_SET_MOMENTS_HPP

#include "FELevelSetMoments.hpp"
#include <string>

namespace FELevelSetMoments {

inline unsigned Hex27Index(unsigned i, unsigned j, unsigned k) {
  return 9U * i + 3U * j + k;
}

inline unsigned HexMomentIndex(unsigned i, unsigned j, unsigned k) {
  unsigned q = i + j + k;
  return q * (q + 1U) * (q + 2U) / 6U + (q - i) * (q - i + 1U) / 2U + k;
}

inline void Hex27ReferenceCoefficientsFromNodalValues(const std::vector<double>& nodalValues, std::vector<std::vector<std::vector<double>>>& coefficient) {
  if(nodalValues.size() != 27) throw std::invalid_argument("HEX27 requires 27 nodal values in FEMuS order.");
  const unsigned node[27][3] = {{0,0,0},{2,0,0},{2,2,0},{0,2,0},{0,0,2},{2,0,2},{2,2,2},{0,2,2},{1,0,0},{2,1,0},{1,2,0},{0,1,0},{1,0,2},{2,1,2},{1,2,2},{0,1,2},{0,0,1},{2,0,1},{2,2,1},{0,2,1},{1,0,1},{2,1,1},{1,2,1},{0,1,1},{1,1,0},{1,1,2},{1,1,1}};
  const long double interpolation[3][3] = {{0.0L,1.0L,0.0L},{-0.5L,0.0L,0.5L},{0.5L,-1.0L,0.5L}};
  std::vector<long double> values(27, 0.0L);
  for(unsigned n = 0; n < 27; ++n) {
    if(!std::isfinite(nodalValues[n])) throw std::invalid_argument("A HEX27 nodal value is not finite.");
    values[Hex27Index(node[n][0], node[n][1], node[n][2])] = nodalValues[n];
  }
  for(unsigned axis = 0; axis < 3; ++axis) {
    std::vector<long double> result(27, 0.0L);
    for(unsigned i = 0; i < 3; ++i) {
      for(unsigned j = 0; j < 3; ++j) {
        for(unsigned k = 0; k < 3; ++k) {
          unsigned index[3] = {i,j,k};
          unsigned power = index[axis];
          for(unsigned sample = 0; sample < 3; ++sample) {
            index[axis] = sample;
            result[Hex27Index(i,j,k)] += interpolation[power][sample] * values[Hex27Index(index[0],index[1],index[2])];
          }
        }
      }
    }
    values.swap(result);
  }
  std::vector<std::vector<std::vector<double>>> result(3, std::vector<std::vector<double>>(3, std::vector<double>(3, 0.)));
  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = 0; j < 3; ++j) {
      for(unsigned k = 0; k < 3; ++k) {
        result[i][j][k] = static_cast<double>(values[Hex27Index(i,j,k)]);
        if(!std::isfinite(result[i][j][k])) throw std::overflow_error("A HEX27 coefficient exceeds double precision range.");
      }
    }
  }
  coefficient.swap(result);
}

inline void Hex27FaceNodalValues(const std::vector<double>& nodalValues, unsigned axis, unsigned side, std::vector<double>& faceValues) {
  if(nodalValues.size() != 27 || axis > 2 || side > 1) throw std::invalid_argument("Invalid HEX27 face request.");
  const unsigned node[27] = {0,16,4,11,23,15,3,19,7,8,20,12,24,26,25,10,22,14,1,17,5,9,21,13,2,18,6};
  const unsigned face[9][2] = {{0,0},{2,0},{2,2},{0,2},{1,0},{2,1},{1,2},{0,1},{1,1}};
  unsigned first = axis == 0 ? 1 : 0;
  unsigned second = axis == 2 ? 1 : 2;
  std::vector<double> values(9);
  for(unsigned n = 0; n < 9; ++n) {
    unsigned index[3] = {0,0,0}; index[axis] = 2U*side; index[first] = face[n][0]; index[second] = face[n][1];
    values[n] = nodalValues[node[Hex27Index(index[0],index[1],index[2])]];
    if(!std::isfinite(values[n])) throw std::invalid_argument("A HEX27 face value is not finite.");
  }
  faceValues.swap(values);
}

inline void Hex27Normalize(const std::vector<std::vector<std::vector<double>>>& coefficient, std::vector<long double>& normalized, long double& scale) {
  if(coefficient.size() != 3) throw std::invalid_argument("Expected a 3 by 3 by 3 HEX27 coefficient tensor.");
  scale = 0.0L;
  normalized.assign(27, 0.0L);
  for(unsigned i = 0; i < 3; ++i) {
    if(coefficient[i].size() != 3) throw std::invalid_argument("Expected a 3 by 3 by 3 HEX27 coefficient tensor.");
    for(unsigned j = 0; j < 3; ++j) {
      if(coefficient[i][j].size() != 3) throw std::invalid_argument("Expected a 3 by 3 by 3 HEX27 coefficient tensor.");
      for(unsigned k = 0; k < 3; ++k) {
        if(!std::isfinite(coefficient[i][j][k])) throw std::invalid_argument("A HEX27 coefficient is not finite.");
        scale = std::max(scale, std::abs(static_cast<long double>(coefficient[i][j][k])));
      }
    }
  }
  if(scale == 0.0L) scale = 1.0L;
  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = 0; j < 3; ++j) {
      for(unsigned k = 0; k < 3; ++k) normalized[Hex27Index(i,j,k)] = coefficient[i][j][k] / scale;
    }
  }
}

inline void Hex27OtherAxes(unsigned axis, unsigned& first, unsigned& second) {
  first = axis == 0 ? 1 : 0;
  second = axis == 2 ? 1 : 2;
}

// Substitute each native coordinate by 2*u-1 only inside the integration chart.
inline void ReferenceCubePolynomialToChart(const std::vector<long double>& coefficient, std::vector<long double>& chart) {
  const unsigned size = SurfacePolynomialCubeSize(static_cast<unsigned>(coefficient.size()));
  std::vector<std::vector<long double>> transform(size,std::vector<long double>(size,0.0L));
  for(unsigned power = 0; power < size; ++power) {
    for(unsigned source = power; source < size; ++source) transform[power][source] = ((source-power)%2U == 0 ? 1.0L : -1.0L)*Binomial(source,power)*Power(2.0L,power);
  }
  chart = coefficient;
  for(unsigned axis = 0; axis < 3; ++axis) {
    std::vector<long double> result(chart.size(),0.0L);
    for(unsigned i = 0; i < size; ++i) {
      for(unsigned j = 0; j < size; ++j) {
        for(unsigned k = 0; k < size; ++k) {
          unsigned index[3] = {i,j,k};
          const unsigned power = index[axis];
          long double value = 0.0L;
          long double correction = 0.0L;
          for(unsigned source = power; source < size; ++source) {
            index[axis] = source;
            const long double term = transform[power][source]*chart[size*size*index[0]+size*index[1]+index[2]]-correction;
            const long double next = value+term;
            correction = (next-value)-term;
            value = next;
          }
          if(!std::isfinite(value)) throw std::overflow_error("A reference polynomial chart exceeds long-double range.");
          result[size*size*i+size*j+k] = value;
        }
      }
    }
    chart.swap(result);
  }
}

inline void Hex27ReferenceSurfaceGradient(const std::vector<long double>& coefficient, long double scale, std::vector<std::vector<long double>>& gradient, const std::vector<std::vector<double>>* inverseJacobian = nullptr, const std::vector<std::vector<double>>* geometryCoordinates = nullptr) {
  gradient.assign(3,std::vector<long double>(27,0.0L));
  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = 0; j < 3; ++j) {
      for(unsigned k = 0; k < 3; ++k) {
        const long double value = scale*coefficient[Hex27Index(i,j,k)];
        if(i != 0) gradient[0][Hex27Index(i-1U,j,k)] = i*value;
        if(j != 0) gradient[1][Hex27Index(i,j-1U,k)] = j*value;
        if(k != 0) gradient[2][Hex27Index(i,j,k-1U)] = k*value;
      }
    }
  }
  if(inverseJacobian != nullptr) SurfaceGradientPhysical(gradient,*inverseJacobian);
  if(geometryCoordinates != nullptr) {
    std::vector<std::vector<long double>> coordinates(3,std::vector<long double>(27,0.0L));
    for(unsigned component = 0; component < 3; ++component) {
      std::vector<double> nodal(27);
      for(unsigned node = 0; node < 27; ++node) nodal[node] = (*geometryCoordinates)[component][node]-(*geometryCoordinates)[component][0];
      std::vector<std::vector<std::vector<double>>> polynomial;
      Hex27ReferenceCoefficientsFromNodalValues(nodal,polynomial);
      for(unsigned i = 0; i < 3; ++i) {
        for(unsigned j = 0; j < 3; ++j) {
          for(unsigned k = 0; k < 3; ++k) coordinates[component][Hex27Index(i,j,k)] = polynomial[i][j][k];
        }
      }
    }
    SurfaceGradientCofactor3D(gradient,coordinates,1.0L);
  }
}

inline std::vector<long double> Hex27Bernstein(const std::vector<long double>& coefficient) {
  std::vector<long double> result(27, 0.0L);
  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = 0; j < 3; ++j) {
      for(unsigned k = 0; k < 3; ++k) {
        for(unsigned r = 0; r <= i; ++r) {
          for(unsigned s = 0; s <= j; ++s) {
            for(unsigned t = 0; t <= k; ++t) result[Hex27Index(i,j,k)] += coefficient[Hex27Index(r,s,t)] * Binomial(i,r) * Binomial(j,s) * Binomial(k,t) / (Binomial(2,r) * Binomial(2,s) * Binomial(2,t));
          }
        }
      }
    }
  }
  return result;
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

inline void Hex27AddEventRoots(const std::vector<long double>& polynomial, std::vector<long double>& events) {
  if(PolynomialNorm(polynomial) != 0.0L) PolynomialRoots(polynomial, 0.0L, 1.0L, events);
}











}

#endif

