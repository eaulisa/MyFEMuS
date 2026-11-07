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

inline void Hex27UnitCoefficientsFromNodalValues(const std::vector<double>& nodalValues, std::vector<std::vector<std::vector<double>>>& coefficient) {
  if(nodalValues.size() != 27) throw std::invalid_argument("HEX27 requires 27 nodal values in FEMuS order.");
  const unsigned node[27][3] = {{0,0,0},{2,0,0},{2,2,0},{0,2,0},{0,0,2},{2,0,2},{2,2,2},{0,2,2},{1,0,0},{2,1,0},{1,2,0},{0,1,0},{1,0,2},{2,1,2},{1,2,2},{0,1,2},{0,0,1},{2,0,1},{2,2,1},{0,2,1},{1,0,1},{2,1,1},{1,2,1},{0,1,1},{1,1,0},{1,1,2},{1,1,1}};
  const long double interpolation[3][3] = {{1.0L,0.0L,0.0L},{-3.0L,4.0L,-1.0L},{2.0L,-4.0L,2.0L}};
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

inline void Hex27ReferenceSurfaceGradient(const std::vector<long double>& coefficient, long double scale, std::vector<std::vector<long double>>& gradient, const std::vector<std::vector<double>>* inverseJacobian = nullptr, const std::vector<std::vector<double>>* geometryCoordinates = nullptr) {
  gradient.assign(3,std::vector<long double>(27,0.0L));
  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = 0; j < 3; ++j) {
      for(unsigned k = 0; k < 3; ++k) {
        const long double value = 0.5L*scale*coefficient[Hex27Index(i,j,k)];
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
      Hex27UnitCoefficientsFromNodalValues(nodal,polynomial);
      for(unsigned i = 0; i < 3; ++i) {
        for(unsigned j = 0; j < 3; ++j) {
          for(unsigned k = 0; k < 3; ++k) coordinates[component][Hex27Index(i,j,k)] = polynomial[i][j][k];
        }
      }
    }
    SurfaceGradientCofactor3D(gradient,coordinates,0.5L);
  }
}

inline long double Hex27ReferenceSurfaceGradientBound(const std::vector<std::vector<long double>>& gradient) {
  long double bound = 0.0L;
  for(unsigned component = 0; component < gradient.size(); ++component) {
    long double magnitude = 0.0L;
    for(unsigned term = 0; term < gradient[component].size(); ++term) magnitude += std::abs(gradient[component][term]);
    bound = std::hypot(bound,magnitude);
  }
  return bound*(1.0L+256.0L*std::numeric_limits<long double>::epsilon());
}

inline void Hex27SliceReferenceSurfaceGradient(const std::vector<std::vector<long double>>& gradient, unsigned axis, long double point, std::vector<std::vector<std::vector<long double>>>& slice) {
  const unsigned size = SurfacePolynomialCubeSize(gradient[0].size());
  unsigned first;
  unsigned second;
  Hex27OtherAxes(axis,first,second);
  slice.assign(gradient.size(),std::vector<std::vector<long double>>(size,std::vector<long double>(size,0.0L)));
  for(unsigned component = 0; component < gradient.size(); ++component) {
    for(unsigned i = 0; i < size; ++i) {
      for(unsigned j = 0; j < size; ++j) {
        unsigned powers[3] = {0,0,0};
        powers[first] = i;
        powers[second] = j;
        long double value = 0.0L;
        for(unsigned k = size; k > 0; --k) {
          powers[axis] = k-1U;
          value = std::fma(value,point,gradient[component][size*size*powers[0]+size*powers[1]+powers[2]]);
        }
        slice[component][i][j] = value;
      }
    }
  }
}

inline bool Hex27SurfaceGradientIndependent(const std::vector<std::vector<long double>>& gradient, unsigned axis) {
  const unsigned size = SurfacePolynomialCubeSize(gradient[0].size());
  for(unsigned component = 0; component < gradient.size(); ++component) {
    for(unsigned i = 0; i < size; ++i) {
      for(unsigned j = 0; j < size; ++j) {
        for(unsigned k = 0; k < size; ++k) {
          const unsigned powers[3] = {i,j,k};
          if(powers[axis] != 0 && gradient[component][size*size*i+size*j+k] != 0.0L) return false;
        }
      }
    }
  }
  return true;
}

inline void Hex27Slice(const std::vector<long double>& coefficient, unsigned axis, long double point, std::vector<std::vector<double>>& slice) {
  unsigned first, second;
  Hex27OtherAxes(axis, first, second);
  slice.assign(3, std::vector<double>(3, 0.));
  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = 0; j < 3; ++j) {
      unsigned powers[3] = {0,0,0};
      powers[first] = i;
      powers[second] = j;
      long double value = 0.0L;
      for(unsigned k = 3; k > 0; --k) {
        powers[axis] = k - 1U;
        value = std::fma(value, point, coefficient[Hex27Index(powers[0],powers[1],powers[2])]);
      }
      slice[i][j] = static_cast<double>(value);
    }
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

inline void Hex27SplitBernstein(const std::vector<long double>& input, unsigned axis, long double point, std::vector<long double>& left, std::vector<long double>& right) {
  unsigned first, second;
  Hex27OtherAxes(axis, first, second);
  left.resize(27);
  right.resize(27);
  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = 0; j < 3; ++j) {
      unsigned index[3] = {0,0,0};
      index[first] = i;
      index[second] = j;
      long double value[3];
      for(unsigned k = 0; k < 3; ++k) {
        index[axis] = k;
        value[k] = input[Hex27Index(index[0],index[1],index[2])];
      }
      long double a = (1.0L-point)*value[0] + point*value[1];
      long double b = (1.0L-point)*value[1] + point*value[2];
      long double middle = (1.0L-point)*a + point*b;
      const long double lv[3] = {value[0],a,middle};
      const long double rv[3] = {middle,b,value[2]};
      for(unsigned k = 0; k < 3; ++k) {
        index[axis] = k;
        left[Hex27Index(index[0],index[1],index[2])] = lv[k];
        right[Hex27Index(index[0],index[1],index[2])] = rv[k];
      }
    }
  }
}

inline int Hex27BernsteinPhase(const std::vector<long double>& coefficient, unsigned first, unsigned second, unsigned depth, unsigned maximumDepth, unsigned& budget, bool computeDelta) {
  long double low = coefficient[0];
  long double high = low;
  for(unsigned k = 1; k < 27; ++k) {
    low = std::min(low, coefficient[k]);
    high = std::max(high, coefficient[k]);
  }
  long double padding = 256.0L * std::numeric_limits<long double>::epsilon() * std::max(std::abs(low),std::abs(high));
  if(low > padding || (!computeDelta && low >= 0.0L && high > 0.0L)) return 1;
  if(high < -padding || (!computeDelta && high <= 0.0L && low < 0.0L)) return -1;
  bool positive = false;
  bool negative = false;
  for(unsigned i = 0; i < 3; i += 2) {
    for(unsigned j = 0; j < 3; j += 2) {
      for(unsigned k = 0; k < 3; k += 2) {
        positive = positive || coefficient[Hex27Index(i,j,k)] > padding;
        negative = negative || coefficient[Hex27Index(i,j,k)] < -padding;
      }
    }
  }
  if((positive && negative) || depth >= maximumDepth || budget == 0) return 0;
  --budget;
  std::vector<long double> left, right;
  Hex27SplitBernstein(coefficient, (depth & 1U) == 0 ? first : second, 0.5L, left, right);
  int phase = Hex27BernsteinPhase(left, first, second, depth+1U, maximumDepth, budget, computeDelta);
  if(phase == 0) return 0;
  return Hex27BernsteinPhase(right, first, second, depth+1U, maximumDepth, budget, computeDelta) == phase ? phase : 0;
}

inline int Hex27SlabPhase(const std::vector<long double>& bernstein, unsigned axis, long double left, long double right, unsigned depth, bool computeDelta) {
  std::vector<long double> discarded, slab, restricted;
  Hex27SplitBernstein(bernstein, axis, left, discarded, slab);
  Hex27SplitBernstein(slab, axis, (right-left)/(1.0L-left), restricted, discarded);
  unsigned first, second;
  Hex27OtherAxes(axis, first, second);
  unsigned budget = 2048;
  return Hex27BernsteinPhase(restricted, first, second, 0, std::min(48U,16U+depth), budget, computeDelta);
}

inline long double Hex27DeltaBound(const std::vector<std::vector<long double>>& fields, const std::vector<long double>& widths, unsigned first, unsigned second, unsigned depth, unsigned& budget, long double target) {
  long double low = fields[0][0], high = low;
  for(unsigned k = 1; k < 27; ++k) { low = std::min(low,fields[0][k]); high = std::max(high,fields[0][k]); }
  long double padding = 512.0L*std::numeric_limits<long double>::epsilon()*(1.0L+std::max(std::abs(low),std::abs(high)));
  if(low > padding || high < -padding) return 0.0L;
  long double bound = std::numeric_limits<long double>::infinity();
  for(unsigned axis = 0; axis < 3; ++axis) {
    long double minimum = fields[axis+1U][0], maximum = minimum;
    for(unsigned k = 1; k < 27; ++k) { minimum = std::min(minimum,fields[axis+1U][k]); maximum = std::max(maximum,fields[axis+1U][k]); }
    long double roundoff = 512.0L*std::numeric_limits<long double>::epsilon()*(1.0L+std::max(std::abs(minimum),std::abs(maximum)));
    long double derivative = minimum > roundoff ? minimum-roundoff : (maximum < -roundoff ? -maximum-roundoff : 0.0L);
    if(derivative > 0.0L) {
      long double candidate = 1.0L/derivative;
      for(unsigned d = 0; d < 3; ++d) if(d != axis) candidate *= widths[d];
      bound = std::min(bound,candidate);
    }
  }
  if(bound <= target || depth >= 64 || budget == 0) return bound;
  --budget;
  unsigned splitAxis = (depth & 1U) == 0 ? first : second;
  std::vector<std::vector<long double>> left(4), right(4);
  for(unsigned i = 0; i < 4; ++i) Hex27SplitBernstein(fields[i],splitAxis,0.5L,left[i],right[i]);
  std::vector<long double> halfWidths(widths);
  halfWidths[splitAxis] *= 0.5L;
  long double result = Hex27DeltaBound(left,halfWidths,first,second,depth+1U,budget,target);
  result += Hex27DeltaBound(right,halfWidths,first,second,depth+1U,budget,target);
  return std::min(bound,result);
}

inline long double Hex27SlabDeltaBound(const std::vector<long double>& bernstein, unsigned axis, long double left, long double right, long double target) {
  std::vector<std::vector<long double>> fields(4,std::vector<long double>(27,0.0L));
  fields[0] = bernstein;
  for(unsigned d = 0; d < 3; ++d) {
    unsigned first, second;
    Hex27OtherAxes(d,first,second);
    for(unsigned i = 0; i < 3; ++i) {
      for(unsigned j = 0; j < 3; ++j) {
        unsigned index[3] = {0,0,0}; index[first] = i; index[second] = j;
        index[d] = 0; long double a = bernstein[Hex27Index(index[0],index[1],index[2])];
        index[d] = 1; long double b = bernstein[Hex27Index(index[0],index[1],index[2])];
        index[d] = 2; long double c = bernstein[Hex27Index(index[0],index[1],index[2])];
        const long double derivative[3] = {2.0L*(b-a),c-a,2.0L*(c-b)};
        for(unsigned k = 0; k < 3; ++k) { index[d] = k; fields[d+1U][Hex27Index(index[0],index[1],index[2])] = derivative[k]; }
      }
    }
  }
  for(unsigned i = 0; i < 4; ++i) {
    std::vector<long double> discarded, slab;
    Hex27SplitBernstein(fields[i],axis,left,discarded,slab);
    Hex27SplitBernstein(slab,axis,(right-left)/(1.0L-left),fields[i],discarded);
  }
  std::vector<long double> widths(3,1.0L);
  widths[axis] = right-left;
  unsigned first, second, budget = 16384;
  Hex27OtherAxes(axis,first,second);
  return Hex27DeltaBound(fields,widths,first,second,0,budget,target/64.0L);
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
  if(m == 0 || n == 0) return std::vector<long double>(1, 0.0L);
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

inline long double Hex27Derivative(const std::vector<long double>& coefficient, const long double* point, const unsigned* derivative) {
  long double value = 0.0L;
  for(unsigned i = derivative[0]; i < 3; ++i) {
    for(unsigned j = derivative[1]; j < 3; ++j) {
      for(unsigned k = derivative[2]; k < 3; ++k) {
        const unsigned power[3] = {i,j,k};
        long double term = coefficient[Hex27Index(i,j,k)];
        for(unsigned d = 0; d < 3; ++d) {
          for(unsigned r = 0; r < derivative[d]; ++r) term *= power[d]-r;
          term *= Power(point[d],power[d]-derivative[d]);
        }
        value += term;
      }
    }
  }
  return value;
}

inline bool Hex27RefineEvent(const std::vector<long double>& coefficient, unsigned outer, long double candidate, long double seedFirst, long double seedSecond, unsigned fixedAxis, unsigned fixedSide, long double& result) {
  unsigned first, second;
  Hex27OtherAxes(outer,first,second);
  unsigned variables[3] = {first,second,outer};
  unsigned size = 3;
  long double point[3] = {0.0L,0.0L,0.0L}; point[outer] = candidate; point[first] = seedFirst; point[second] = seedSecond;
  if(fixedAxis < 3) {
    point[fixedAxis] = fixedSide;
    variables[0] = fixedAxis == first ? second : first;
    variables[1] = outer;
    size = 2;
  }
  for(unsigned iteration = 0; iteration < 35; ++iteration) {
    long double matrix[3][4] = {};
    long double residual = 0.0L;
    for(unsigned row = 0; row < size; ++row) {
      unsigned derivative[3] = {0,0,0};
      if(row > 0) derivative[variables[row-1U]] = 1;
      long double value = Hex27Derivative(coefficient,point,derivative);
      residual = std::max(residual,std::abs(value));
      matrix[row][size] = -value;
      for(unsigned column = 0; column < size; ++column) {
        ++derivative[variables[column]];
        matrix[row][column] = Hex27Derivative(coefficient,point,derivative);
        --derivative[variables[column]];
      }
    }
    if(residual <= 32.0L*std::numeric_limits<long double>::epsilon()) {
      bool inside = true;
      for(unsigned d = 0; d < 3; ++d) inside = inside && point[d] >= 0.0L && point[d] <= 1.0L;
      if(inside && std::abs(point[outer]-candidate) < 1.0e-3L) { result = point[outer]; return true; }
      return false;
    }
    for(unsigned column = 0; column < size; ++column) {
      unsigned pivot = column;
      for(unsigned row = column+1U; row < size; ++row) if(std::abs(matrix[row][column]) > std::abs(matrix[pivot][column])) pivot = row;
      if(std::abs(matrix[pivot][column]) < 64.0L*std::numeric_limits<long double>::epsilon()) return false;
      for(unsigned j = column; j <= size; ++j) std::swap(matrix[pivot][j],matrix[column][j]);
      long double scale = matrix[column][column];
      for(unsigned j = column; j <= size; ++j) matrix[column][j] /= scale;
      for(unsigned row = 0; row < size; ++row) {
        if(row != column) {
          long double factor = matrix[row][column];
          for(unsigned j = column; j <= size; ++j) matrix[row][j] -= factor*matrix[column][j];
        }
      }
    }
    for(unsigned column = 0; column < size; ++column) point[variables[column]] += matrix[column][size];
    for(unsigned d = 0; d < 3; ++d) if(point[d] < -0.5L || point[d] > 1.5L) return false;
  }
  return false;
}

inline std::vector<long double> Hex27Events(const std::vector<long double>& coefficient, unsigned outer) {
  unsigned first, second;
  Hex27OtherAxes(outer, first, second);
  std::vector<long double> events = {0.0L,1.0L};
  for(unsigned choice = 0; choice < 2; ++choice) {
    unsigned fiber = choice == 0 ? first : second;
    unsigned remaining = choice == 0 ? second : first;
    std::vector<std::vector<std::vector<long double>>> field(3, std::vector<std::vector<long double>>(3, std::vector<long double>(3,0.0L)));
    for(unsigned i = 0; i < 3; ++i) {
      for(unsigned j = 0; j < 3; ++j) {
        for(unsigned k = 0; k < 3; ++k) {
          unsigned powers[3] = {0,0,0};
          powers[fiber] = i; powers[remaining] = j; powers[outer] = k;
          field[i][j][k] = coefficient[Hex27Index(powers[0],powers[1],powers[2])];
        }
      }
    }
    for(unsigned side = 0; side < 2; ++side) {
      std::vector<std::vector<long double>> face(3, std::vector<long double>(3,0.0L));
      for(unsigned i = 0; i < 3; ++i) {
        for(unsigned j = 0; j < 3; ++j) {
          if(j == 0 || side == 1) face[i] = PolynomialSum(face[i],field[i][j]);
        }
        Hex27AddEventRoots(face[i], events);
      }
      Hex27AddEventRoots(PolynomialSum(PolynomialProduct(face[1],face[1]),PolynomialProduct(face[0],face[2]),-4.0L), events);
      Hex27AddEventRoots(PolynomialSum(PolynomialSum(face[0],face[1]),face[2]), events);
    }
    std::vector<std::vector<long double>> discriminant(5, std::vector<long double>(1,0.0L));
    bool quadratic = false;
    bool linear = false;
    for(unsigned i = 0; i < 3; ++i) { quadratic = quadratic || PolynomialNorm(field[2][i]) != 0.0L; linear = linear || PolynomialNorm(field[1][i]) != 0.0L; }
    if(quadratic) {
      for(unsigned i = 0; i < 3; ++i) {
        for(unsigned j = 0; j < 3; ++j) discriminant[i+j] = PolynomialSum(discriminant[i+j],PolynomialSum(PolynomialProduct(field[1][i],field[1][j]),PolynomialProduct(field[0][i],field[2][j]),-4.0L));
      }
    }
    else if(linear) {
      Hex27AddEventRoots(Hex27Resultant(field[0],field[1]), events);
      discriminant = field[1];
    }
    else discriminant = field[0];
    Hex27TrimBivariate(discriminant);
    for(unsigned i = 0; i < discriminant.size(); ++i) Hex27AddEventRoots(discriminant[i], events);
    if(discriminant.size() > 1) {
      std::vector<std::vector<long double>> derivative(discriminant.size()-1U);
      for(unsigned i = 1; i < discriminant.size(); ++i) {
        derivative[i-1U] = discriminant[i];
        for(unsigned j = 0; j < derivative[i-1U].size(); ++j) derivative[i-1U][j] *= i;
      }
      Hex27AddEventRoots(Hex27Resultant(discriminant,derivative), events);
    }
  }
  SortDistinct(events);
  std::vector<long double> candidates(events);
  for(unsigned n = 1; n+1U < candidates.size(); ++n) {
    for(unsigned i = 1; i <= 3; ++i) {
      for(unsigned j = 1; j <= 3; ++j) {
        long double refined;
        if(Hex27RefineEvent(coefficient,outer,candidates[n],0.25L*i,0.25L*j,3,0,refined)) events.push_back(refined);
      }
    }
    for(unsigned axis = 0; axis < 3; ++axis) {
      if(axis != outer) {
        for(unsigned side = 0; side < 2; ++side) {
          for(unsigned seed = 1; seed <= 3; ++seed) {
            long double refined;
            if(Hex27RefineEvent(coefficient,outer,candidates[n],0.25L*seed,0.25L*seed,axis,side,refined)) events.push_back(refined);
          }
        }
      }
    }
  }
  SortDistinct(events);
  return events;
}

inline bool Hex27FlatSlice(const std::vector<long double>& coefficient, unsigned axis) {
  unsigned first, second;
  Hex27OtherAxes(axis,first,second);
  std::vector<std::vector<long double>> polynomials;
  std::vector<long double> candidates;
  bool selected = false;
  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = 0; j < 3; ++j) {
      unsigned powers[3] = {0,0,0}; powers[first] = i; powers[second] = j;
      std::vector<long double> polynomial(3);
      for(unsigned k = 0; k < 3; ++k) { powers[axis] = k; polynomial[k] = coefficient[Hex27Index(powers[0],powers[1],powers[2])]; }
      if(!selected && PolynomialNorm(polynomial) != 0.0L) { Hex27AddEventRoots(polynomial,candidates); selected = true; }
      polynomials.push_back(polynomial);
    }
  }
  for(unsigned k = 0; k < candidates.size(); ++k) {
    bool common = true;
    for(unsigned i = 0; i < polynomials.size(); ++i) common = common && std::abs(PolynomialValue(polynomials[i],candidates[k])) <= 256.0L * std::numeric_limits<long double>::epsilon() * PolynomialBound(polynomials[i],candidates[k]);
    if(common) return true;
  }
  return false;
}

inline unsigned Hex27GraphAxes(const std::vector<long double>& bernstein) {
  long double norm = 1.0L;
  for(unsigned k = 0; k < bernstein.size(); ++k) norm = std::max(norm,std::abs(bernstein[k]));
  long double padding = 1024.0L*std::numeric_limits<long double>::epsilon()*norm;
  long double conditioning = 32.0L*std::sqrt(std::numeric_limits<long double>::epsilon())*norm;
  unsigned mask = 0;
  for(unsigned axis = 0; axis < 3; ++axis) {
    unsigned first,second;
    Hex27OtherAxes(axis,first,second);
    long double lowMinimum = std::numeric_limits<long double>::infinity();
    long double lowMaximum = -lowMinimum;
    long double highMinimum = lowMinimum;
    long double highMaximum = lowMaximum;
    long double derivativeMinimum = lowMinimum;
    long double derivativeMaximum = lowMaximum;
    for(unsigned i = 0; i < 3; ++i) {
      for(unsigned j = 0; j < 3; ++j) {
        unsigned index[3] = {0,0,0};
        index[first] = i;
        index[second] = j;
        index[axis] = 0;
        long double low = bernstein[Hex27Index(index[0],index[1],index[2])];
        index[axis] = 1;
        long double middle = bernstein[Hex27Index(index[0],index[1],index[2])];
        index[axis] = 2;
        long double high = bernstein[Hex27Index(index[0],index[1],index[2])];
        lowMinimum = std::min(lowMinimum,low);
        lowMaximum = std::max(lowMaximum,low);
        highMinimum = std::min(highMinimum,high);
        highMaximum = std::max(highMaximum,high);
        const long double derivative[3] = {2.0L*(middle-low),high-low,2.0L*(high-middle)};
        for(unsigned k = 0; k < 3; ++k) {
          derivativeMinimum = std::min(derivativeMinimum,derivative[k]);
          derivativeMaximum = std::max(derivativeMaximum,derivative[k]);
        }
      }
    }
    if((lowMaximum < -padding && highMinimum > padding && derivativeMinimum > conditioning) || (lowMinimum > padding && highMaximum < -padding && derivativeMaximum < -conditioning)) mask |= 1U<<axis;
  }
  return mask;
}

inline void Hex27PlanarMoments(const std::vector<long double>& coefficient, unsigned degree, bool computeDelta, double absoluteTolerance, double relativeTolerance, std::vector<long double>& values) {
  unsigned count = (degree+1U)*(degree+2U)*(degree+3U)/6U;
  values.assign(2U*count,0.0L);
  const long double normal[3] = {coefficient[Hex27Index(1,0,0)],coefficient[Hex27Index(0,1,0)],coefficient[Hex27Index(0,0,1)]};
  unsigned pivot = 0;
  for(unsigned d = 1; d < 3; ++d) if(std::abs(normal[d]) > std::abs(normal[pivot])) pivot = d;
  if(normal[pivot] == 0.0L) {
    if(coefficient[0] == 0.0L && computeDelta) throw std::domain_error("The delta of an identically zero HEX27 field is undefined.");
    for(unsigned i = 0; i <= degree; ++i) {
      for(unsigned j = 0; j <= degree-i; ++j) {
        for(unsigned k = 0; k <= degree-i-j; ++k) values[HexMomentIndex(i,j,k)] = Heaviside(coefficient[0])/((i+1.0L)*(j+1.0L)*(k+1.0L));
      }
    }
    return;
  }
  std::vector<std::vector<std::vector<std::vector<double>>>> faces(3,std::vector<std::vector<std::vector<double>>>(2));
  for(unsigned axis = 0; axis < 3; ++axis) {
    for(unsigned side = 0; side < 2; ++side) {
      std::vector<std::vector<double>> slice, unused;
      Hex27Slice(coefficient,axis,side,slice);
      IntegrateQuad9Unit(slice,degree,faces[axis][side],unused,absoluteTolerance,relativeTolerance,false);
    }
  }
  for(unsigned q = 0; q <= degree; ++q) {
    for(int ii = q; ii >= 0; --ii) {
      for(int jj = q-static_cast<unsigned>(ii); jj >= 0; --jj) {
        unsigned i = static_cast<unsigned>(ii), j = static_cast<unsigned>(jj), k = q-i-j;
        unsigned powers[3] = {i,j,k};
        unsigned first, second;
        Hex27OtherAxes(pivot,first,second);
        long double delta = faces[pivot][1][powers[first]][powers[second]];
        if(powers[pivot] == 0) delta -= faces[pivot][0][powers[first]][powers[second]];
        else {
          --powers[pivot];
          delta -= (powers[pivot]+1U)*values[HexMomentIndex(powers[0],powers[1],powers[2])];
          ++powers[pivot];
        }
        delta /= normal[pivot];
        long double boundary = faces[0][1][j][k] + faces[1][1][i][k] + faces[2][1][i][j];
        unsigned index = HexMomentIndex(i,j,k);
        values[index] = (boundary+coefficient[0]*delta)/(q+3U);
        if(computeDelta) values[count+index] = delta;
      }
    }
  }
}

}

#endif

