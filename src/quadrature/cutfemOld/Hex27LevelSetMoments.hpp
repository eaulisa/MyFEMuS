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

inline int ClassifyHex27NodalValues(const std::vector<double>& nodalValues) {
  std::vector<std::vector<std::vector<double>>> coefficient;
  Hex27UnitCoefficientsFromNodalValues(nodalValues,coefficient);
  std::vector<long double> normalized;
  long double scale;
  Hex27Normalize(coefficient,normalized,scale);
  std::vector<long double> bernstein = Hex27Bernstein(normalized);
  unsigned budget = 0;
  return Hex27BernsteinPhase(bernstein,0,1,0,0,budget,true);
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

inline void Hex27Sample(const std::vector<long double>& coefficient, unsigned degree, unsigned axis, long double point, long double jacobian, long double absoluteTolerance, long double deltaTolerance, long double relativeTolerance, bool computeDelta, std::vector<long double>& value, std::vector<long double>& error, bool& cut, unsigned& evaluations, unsigned maximumEvaluations) {
  if(evaluations >= maximumEvaluations) throw std::runtime_error("HEX27 slice evaluation budget exceeded.");
  ++evaluations;
  std::vector<std::vector<double>> slice;
  Hex27Slice(coefficient,axis,point,slice);
  std::vector<std::vector<double>> positive, delta;
  std::vector<long double> sliceError;
  double innerTolerance = static_cast<double>(0.02L * std::min(absoluteTolerance,computeDelta ? deltaTolerance : absoluteTolerance));
  if(!(innerTolerance > 0.)) innerTolerance = std::numeric_limits<double>::min();
  IntegrateCanonicalPolynomial(slice,degree,false,positive,delta,innerTolerance,static_cast<double>(0.05L*relativeTolerance),computeDelta,false,32,28,0,false,nullptr,&sliceError);
  unsigned count = (degree+1U)*(degree+2U)*(degree+3U)/6U;
  value.assign(2U*count,0.0L);
  error.assign(2U*count,0.0L);
  unsigned first, second;
  Hex27OtherAxes(axis,first,second);
  unsigned sliceCount = (degree+1U)*(degree+2U)/2U;
  for(unsigned i = 0; i <= degree; ++i) {
    for(unsigned j = 0; j <= degree-i; ++j) {
      unsigned source = MomentIndex(degree,i,j);
      long double positiveError = sliceError[source]+std::numeric_limits<double>::epsilon()*std::abs(positive[i][j]);
      long double deltaError = computeDelta ? sliceError[sliceCount+source]+std::numeric_limits<double>::epsilon()*std::abs(delta[i][j]) : 0.0L;
      long double power = 1.0L;
      for(unsigned k = 0; k <= degree-i-j; ++k) {
        unsigned index[3] = {0,0,0}; index[first] = i; index[second] = j; index[axis] = k;
        unsigned target = HexMomentIndex(index[0],index[1],index[2]);
        value[target] = jacobian * power * positive[i][j];
        error[target] = jacobian*power*positiveError;
        if(computeDelta) {
          value[count+target] = jacobian * power * delta[i][j];
          error[count+target] = jacobian*power*deltaError;
        }
        power *= point;
      }
    }
  }
  cut = (positive[0][0] > 0.0 && positive[0][0] < 1.0) || (computeDelta && delta[0][0] != 0.0);
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

inline void Hex27IntegrateInterval(const std::vector<long double>& coefficient, const std::vector<long double>& bernstein, unsigned degree, unsigned axis, long double left, long double right, long double thetaLeft, long double thetaRight, long double absoluteTolerance, long double deltaTolerance, long double relativeTolerance, bool computeDelta, unsigned depth, unsigned maximumDepth, unsigned& evaluations, unsigned maximumEvaluations, unsigned& accepted, std::vector<long double>& total, std::vector<long double>& totalError, bool regularGraph = false) {
  // Gauss-Kronrod 15/7 constants: https://jacobwilliams.github.io/quadpack/proc/dqk15.html.
  const long double nodes[8] = {0.991455371120812639206854697526329L,0.949107912342758524526189684047851L,0.864864423359769072789712788640926L,0.741531185599394439863864773280788L,0.586087235467691130294144838258730L,0.405845151377397166906606412076961L,0.207784955007898467600689403773245L,0.0L};
  const long double kronrod[8] = {0.022935322010529224963732008058970L,0.063092092629978553290700663189204L,0.104790010322250183839876322541518L,0.140653259715525918745189590510238L,0.169004726639267902826583426598550L,0.190350578064785409913256402421014L,0.204432940075298892414161999234649L,0.209482141084727828012999174891714L};
  const long double gauss[4] = {0.129484966168869693270611432679082L,0.279705391489276667901467771423780L,0.381830050505118944950369775488975L,0.417959183673469387755102040816327L};
  long double width = right-left;
  long double zLeft = left + width*(regularGraph ? thetaLeft : Power(std::sin(thetaLeft),2));
  long double zRight = left + width*(regularGraph ? thetaRight : Power(std::sin(thetaRight),2));
  if(!(zRight > zLeft)) throw std::runtime_error("HEX27 slice interval reached floating point resolution.");
  unsigned count = (degree+1U)*(degree+2U)*(degree+3U)/6U;
  int phase = regularGraph ? 0 : Hex27SlabPhase(bernstein,axis,zLeft,zRight,depth,computeDelta);
  if(phase != 0) {
    if(phase > 0) {
      for(unsigned i = 0; i <= degree; ++i) {
        for(unsigned j = 0; j <= degree-i; ++j) {
          for(unsigned k = 0; k <= degree-i-j; ++k) {
            unsigned powers[3] = {i,j,k};
            long double value = MonomialIntegral(zLeft,zRight,powers[axis]);
            for(unsigned d = 0; d < 3; ++d) if(d != axis) value /= powers[d]+1U;
            total[HexMomentIndex(i,j,k)] += value;
          }
        }
      }
    }
    ++accepted;
    return;
  }
  long double middle = 0.5L*(thetaLeft+thetaRight);
  long double halfWidth = 0.5L*(thetaRight-thetaLeft);
  std::vector<long double> high(2U*count,0.0L), low(2U*count,0.0L), innerError(2U*count,0.0L), absolute(2U*count,0.0L);
  std::vector<std::vector<long double>> samples;
  std::vector<long double> sampleWeights;
  bool anyCut = false;
  bool failed = false;
  long double split = middle;
  for(unsigned node = 0; node < 8 && !failed; ++node) {
    unsigned sides = node == 7 ? 1 : 2;
    for(unsigned side = 0; side < sides; ++side) {
      long double theta = middle + (side == 0 ? -1.0L : 1.0L)*halfWidth*nodes[node];
      long double sine = std::sin(theta);
      long double point = left + width*(regularGraph ? theta : sine*sine);
      long double jacobian = regularGraph ? width : 2.0L*width*sine*std::cos(theta);
      std::vector<long double> value, error;
      bool sampleCut = false;
      try {
        Hex27Sample(coefficient,degree,axis,point,jacobian,absoluteTolerance,deltaTolerance,relativeTolerance,computeDelta,value,error,sampleCut,evaluations,maximumEvaluations);
      }
      catch(const std::domain_error&) { failed = true; split = theta; break; }
      anyCut = anyCut || sampleCut;
      samples.push_back(value);
      sampleWeights.push_back(kronrod[node]);
      for(unsigned k = 0; k < high.size(); ++k) {
        high[k] += kronrod[node]*value[k];
        absolute[k] += kronrod[node]*std::abs(value[k]);
        innerError[k] += kronrod[node]*error[k];
        if((node & 1U) != 0) low[k] += gauss[node/2U]*value[k];
      }
    }
  }
  bool converged = !failed && anyCut;
  std::vector<long double> errors(2U*count,0.0L);
  long double intervalFraction = std::max(1.0e-2L,width*(thetaRight-thetaLeft)/(regularGraph ? 1.0L : 0.5L*std::acos(-1.0L)));
  if(failed && zRight-zLeft <= 0.25L*absoluteTolerance*intervalFraction) {
    long double bound = computeDelta ? Hex27SlabDeltaBound(bernstein,axis,zLeft,zRight,0.25L*deltaTolerance*intervalFraction) : 0.0L;
    if(bound <= 0.25L*deltaTolerance*intervalFraction) {
      for(unsigned k = 0; k < count; ++k) { totalError[k] += zRight-zLeft; totalError[count+k] += bound; }
      ++accepted;
      return;
    }
  }
  if(!failed) {
    long double fraction = intervalFraction;
    // Charge certified unresolved-slab bounds to the same per-moment budget checked below.
    if(!anyCut && zRight-zLeft <= absoluteTolerance*fraction) {
      long double bound = computeDelta ? Hex27SlabDeltaBound(bernstein,axis,zLeft,zRight,deltaTolerance*fraction) : 0.0L;
      if(bound <= deltaTolerance*fraction) {
        converged = true;
        for(unsigned k = 0; k < count; ++k) { innerError[k] += (zRight-zLeft)/halfWidth; innerError[count+k] += bound/halfWidth; }
      }
    }
    for(unsigned k = 0; k < high.size(); ++k) {
      long double variation = 0.0L;
      for(unsigned node = 0; node < samples.size(); ++node) variation += sampleWeights[node]*std::abs(samples[node][k]-0.5L*high[k]);
      long double error = std::abs(high[k]-low[k]);
      if(variation > 0.0L && error > 0.0L) error = variation*std::min(1.0L,std::pow(200.0L*error/variation,1.5L));
      error = std::max(error,50.0L*std::numeric_limits<long double>::epsilon()*absolute[k]);
      errors[k] = halfWidth*(error+innerError[k]);
      high[k] *= halfWidth;
      long double tolerance = (k < count ? absoluteTolerance : deltaTolerance)*fraction + relativeTolerance*std::abs(high[k]);
      if(errors[k] > tolerance) converged = false;
    }
  }
  if(converged) {
    for(unsigned k = 0; k < high.size(); ++k) { total[k] += high[k]; totalError[k] += errors[k]; }
    ++accepted;
    return;
  }
  if(depth >= maximumDepth || !(split > thetaLeft && split < thetaRight)) throw std::runtime_error("HEX27 integration did not resolve the requested tolerance or slice topology.");
  Hex27IntegrateInterval(coefficient,bernstein,degree,axis,left,right,thetaLeft,split,absoluteTolerance,deltaTolerance,relativeTolerance,computeDelta,depth+1U,maximumDepth,evaluations,maximumEvaluations,accepted,total,totalError,regularGraph);
  Hex27IntegrateInterval(coefficient,bernstein,degree,axis,left,right,split,thetaRight,absoluteTolerance,deltaTolerance,relativeTolerance,computeDelta,depth+1U,maximumDepth,evaluations,maximumEvaluations,accepted,total,totalError,regularGraph);
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

inline void IntegrateHex27Unit(const std::vector<std::vector<std::vector<double>>>& coefficient, unsigned degree, std::vector<std::vector<std::vector<double>>>& positive, std::vector<std::vector<std::vector<double>>>& delta, double absoluteTolerance = 1.0e-10, double relativeTolerance = 1.0e-10, bool computeDelta = true, unsigned maximumDepth = 24, unsigned integrationAxis = 0, unsigned maximumEvaluations = 200000, std::vector<double>* estimatedError = nullptr) {
  if(degree > 64 || !(absoluteTolerance > 0.) || !(relativeTolerance >= 0.) || !std::isfinite(absoluteTolerance) || !std::isfinite(relativeTolerance) || maximumDepth > 48 || integrationAxis > 3 || maximumEvaluations == 0) throw std::invalid_argument("Invalid HEX27 integration input.");
  std::vector<long double> normalized;
  long double scale;
  Hex27Normalize(coefficient,normalized,scale);
  unsigned count = (degree+1U)*(degree+2U)*(degree+3U)/6U;
  std::vector<long double> values(2U*count,0.0L), errors(2U*count,0.0L);
  bool planar = true;
  long double derivative[3] = {0.0L,0.0L,0.0L};
  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = 0; j < 3; ++j) {
      for(unsigned k = 0; k < 3; ++k) {
        long double value = std::abs(normalized[Hex27Index(i,j,k)]);
        if(i+j+k > 1 && value != 0.0L) planar = false;
        derivative[0] += i*value; derivative[1] += j*value; derivative[2] += k*value;
      }
    }
  }
  unsigned accepted = 0, evaluations = 0, selectedAxis = 0;
  std::vector<long double> bernstein = Hex27Bernstein(normalized);
  unsigned graphAxes = planar ? 0 : Hex27GraphAxes(bernstein);
  unsigned phaseBudget = 0;
  int uniformPhase = Hex27BernsteinPhase(bernstein,0,1,0,0,phaseBudget,computeDelta);
  if(uniformPhase != 0) {
    if(uniformPhase > 0) {
      for(unsigned i = 0; i <= degree; ++i) {
        for(unsigned j = 0; j <= degree-i; ++j) {
          for(unsigned k = 0; k <= degree-i-j; ++k) values[HexMomentIndex(i,j,k)] = 1.0L/((i+1.0L)*(j+1.0L)*(k+1.0L));
        }
      }
    }
  }
  else if(planar) Hex27PlanarMoments(normalized,degree,computeDelta,absoluteTolerance,relativeTolerance,values);
  else {
    unsigned axes[3] = {0,1,2};
    for(unsigned i = 0; i < 3; ++i) {
      for(unsigned j = i+1U; j < 3; ++j) if(derivative[axes[j]] < derivative[axes[i]]) std::swap(axes[i],axes[j]);
    }
    if(integrationAxis != 0) axes[0] = integrationAxis-1U;
    bool completed = false;
    std::string failure = "No regular HEX27 slicing direction was found.";
    unsigned attempts = integrationAxis == 0 ? 3 : 1;
    for(unsigned attempt = 0; attempt < attempts && !completed; ++attempt) {
      unsigned axis = axes[attempt];
      bool regularGraph = (graphAxes & ~(1U<<axis)) != 0;
      if(computeDelta && !regularGraph && Hex27FlatSlice(normalized,axis)) { failure = "A slicing direction contains an entire zero slice."; continue; }
      std::fill(values.begin(),values.end(),0.0L);
      std::fill(errors.begin(),errors.end(),0.0L);
      evaluations = 0; accepted = 0;
      try {
        long double positiveWorkingTolerance = absoluteTolerance;
        long double deltaWorkingTolerance = std::max(static_cast<long double>(absoluteTolerance)*scale,static_cast<long double>(relativeTolerance));
        std::vector<long double> events;
        if(derivative[axis] != 0.0L) events = regularGraph ? std::vector<long double>{0.0L,1.0L} : Hex27Events(normalized,axis);
        for(unsigned refinement = 0; refinement < 4 && !completed; ++refinement) {
          std::fill(values.begin(),values.end(),0.0L);
          std::fill(errors.begin(),errors.end(),0.0L);
          accepted = 0;
          if(derivative[axis] == 0.0L) {
            bool cut = false;
            Hex27Sample(normalized,degree,axis,0.5L,1.0L,positiveWorkingTolerance,deltaWorkingTolerance,relativeTolerance,computeDelta,values,errors,cut,evaluations,maximumEvaluations);
            for(unsigned i = 0; i <= degree; ++i) {
              for(unsigned j = 0; j <= degree-i; ++j) {
                for(unsigned k = 0; k <= degree-i-j; ++k) {
                  unsigned powers[3] = {i,j,k};
                  long double factor = 1.0L/(Power(0.5L,powers[axis])*(powers[axis]+1U));
                  unsigned index = HexMomentIndex(i,j,k);
                  values[index] *= factor;
                  values[count+index] *= factor;
                  errors[index] *= factor;
                  errors[count+index] *= factor;
                }
              }
            }
            accepted = 1;
          }
          else {
            for(unsigned event = 1; event < events.size(); ++event) {
              if(events[event] > events[event-1U]) Hex27IntegrateInterval(normalized,bernstein,degree,axis,events[event-1U],events[event],0.0L,regularGraph ? 1.0L : 0.5L*std::acos(-1.0L),positiveWorkingTolerance,deltaWorkingTolerance,relativeTolerance,computeDelta,0,maximumDepth,evaluations,maximumEvaluations,accepted,values,errors,regularGraph);
            }
          }
          bool positiveConverged = true;
          bool deltaConverged = true;
          long double nextPositiveTolerance = 0.125L*positiveWorkingTolerance;
          long double nextDeltaTolerance = 0.125L*deltaWorkingTolerance;
          for(unsigned k = 0; k < values.size(); ++k) {
            // delta(G/scale) = scale*delta(G); relative tolerances remain dimensionless.
            long double tolerance = (k < count ? absoluteTolerance : static_cast<long double>(absoluteTolerance)*scale)+relativeTolerance*std::abs(values[k]);
            if(k < count) {
              positiveConverged = positiveConverged && errors[k] <= tolerance;
              nextPositiveTolerance = std::min(nextPositiveTolerance,0.25L*tolerance);
            }
            else if(computeDelta) {
              deltaConverged = deltaConverged && errors[k] <= tolerance;
              nextDeltaTolerance = std::min(nextDeltaTolerance,0.25L*tolerance);
            }
          }
          completed = positiveConverged && deltaConverged;
          if(!positiveConverged) positiveWorkingTolerance = nextPositiveTolerance;
          if(!deltaConverged) deltaWorkingTolerance = nextDeltaTolerance;
        }
        if(!completed) throw std::runtime_error("HEX27 accumulated error estimate exceeds the requested tolerance.");
        selectedAxis = axis+1U;
      }
      catch(const std::exception& exception) { failure = exception.what(); }
    }
    if(!completed) throw std::runtime_error("HEX27 interior integration failed: " + failure);
  }
  positive.resize(degree+1U);
  delta.resize(degree+1U);
  long double positiveError = 0.0L, deltaError = 0.0L;
  for(unsigned i = 0; i <= degree; ++i) {
    positive[i].resize(degree-i+1U); delta[i].resize(degree-i+1U);
    for(unsigned j = 0; j <= degree-i; ++j) {
      positive[i][j].resize(degree-i-j+1U); delta[i][j].resize(degree-i-j+1U);
      for(unsigned k = 0; k <= degree-i-j; ++k) {
        unsigned index = HexMomentIndex(i,j,k);
        positive[i][j][k] = static_cast<double>(values[index]);
        delta[i][j][k] = computeDelta ? static_cast<double>(values[count+index]/scale) : 0.;
        if(!std::isfinite(positive[i][j][k]) || !std::isfinite(delta[i][j][k])) throw std::overflow_error("A HEX27 moment exceeds double precision range.");
        positiveError = std::max(positiveError,errors[index]);
        deltaError = std::max(deltaError,errors[count+index]/scale);
      }
    }
  }
  if(estimatedError != nullptr) *estimatedError = {static_cast<double>(positiveError),static_cast<double>(deltaError),static_cast<double>(accepted),static_cast<double>(evaluations),static_cast<double>(selectedAxis)};
}

}

#endif

