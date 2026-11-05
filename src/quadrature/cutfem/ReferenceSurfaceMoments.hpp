#ifndef FEMUS_REFERENCE_SURFACE_MOMENTS_HPP
#define FEMUS_REFERENCE_SURFACE_MOMENTS_HPP

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace FELevelSetMoments {

inline void ValidateSurfaceInverseJacobian(const std::vector<std::vector<double>>& inverseJacobian, unsigned dimension) {
  if(dimension < 1 || dimension > 3 || inverseJacobian.size() != dimension) throw std::invalid_argument("Invalid surface inverse Jacobian dimensions.");
  long double factor[3][3] = {};
  for(unsigned reference = 0; reference < dimension; ++reference) {
    long double scale = 0.0L;
    if(inverseJacobian[reference].size() != dimension) throw std::invalid_argument("Invalid surface inverse Jacobian dimensions.");
    for(unsigned physical = 0; physical < dimension; ++physical) {
      if(!std::isfinite(inverseJacobian[reference][physical])) throw std::invalid_argument("A surface inverse Jacobian entry is not finite.");
      scale = std::max(scale,std::abs(static_cast<long double>(inverseJacobian[reference][physical])));
    }
    if(scale == 0.0L) throw std::invalid_argument("The surface inverse Jacobian is singular.");
    for(unsigned physical = 0; physical < dimension; ++physical) factor[reference][physical] = static_cast<long double>(inverseJacobian[reference][physical])/scale;
  }
  for(unsigned column = 0; column < dimension; ++column) {
    unsigned pivot = column;
    for(unsigned row = column+1U; row < dimension; ++row) {
      if(std::abs(factor[row][column]) > std::abs(factor[pivot][column])) pivot = row;
    }
    if(std::abs(factor[pivot][column]) <= 64.0L*std::numeric_limits<long double>::epsilon()) throw std::invalid_argument("The surface inverse Jacobian is numerically singular.");
    for(unsigned entry = column; entry < dimension; ++entry) std::swap(factor[column][entry],factor[pivot][entry]);
    for(unsigned row = column+1U; row < dimension; ++row) {
      const long double multiplier = factor[row][column]/factor[column][column];
      for(unsigned entry = column+1U; entry < dimension; ++entry) factor[row][entry] = std::fma(-multiplier,factor[column][entry],factor[row][entry]);
    }
  }
}

inline void SurfaceGradientPhysical(std::vector<std::vector<std::vector<long double>>>& gradient, const std::vector<std::vector<double>>& inverseJacobian) {
  const unsigned dimension = static_cast<unsigned>(gradient.size());
  ValidateSurfaceInverseJacobian(inverseJacobian,dimension);
  unsigned rows = 0;
  unsigned columns = 0;
  for(unsigned component = 0; component < dimension; ++component) {
    rows = std::max(rows,static_cast<unsigned>(gradient[component].size()));
    for(unsigned row = 0; row < gradient[component].size(); ++row) columns = std::max(columns,static_cast<unsigned>(gradient[component][row].size()));
  }
  std::vector<std::vector<std::vector<long double>>> transformed(dimension,std::vector<std::vector<long double>>(rows,std::vector<long double>(columns,0.0L)));
  for(unsigned reference = 0; reference < dimension; ++reference) {
    for(unsigned physical = 0; physical < dimension; ++physical) {
      for(unsigned row = 0; row < gradient[reference].size(); ++row) {
        for(unsigned column = 0; column < gradient[reference][row].size(); ++column) transformed[physical][row][column] = std::fma(static_cast<long double>(inverseJacobian[reference][physical]),gradient[reference][row][column],transformed[physical][row][column]);
      }
    }
  }
  gradient.swap(transformed);
}

inline void SurfaceGradientPhysical(std::vector<std::vector<long double>>& gradient, const std::vector<std::vector<double>>& inverseJacobian) {
  const unsigned dimension = static_cast<unsigned>(gradient.size());
  ValidateSurfaceInverseJacobian(inverseJacobian,dimension);
  unsigned count = 0;
  for(unsigned component = 0; component < dimension; ++component) count = std::max(count,static_cast<unsigned>(gradient[component].size()));
  std::vector<std::vector<long double>> transformed(dimension,std::vector<long double>(count,0.0L));
  for(unsigned reference = 0; reference < dimension; ++reference) {
    for(unsigned physical = 0; physical < dimension; ++physical) {
      for(unsigned term = 0; term < gradient[reference].size(); ++term) transformed[physical][term] = std::fma(static_cast<long double>(inverseJacobian[reference][physical]),gradient[reference][term],transformed[physical][term]);
    }
  }
  gradient.swap(transformed);
}

inline long double SurfacePolynomialValue(const std::vector<std::vector<long double>>& polynomial, long double first, long double second) {
  long double value = 0.0L;
  for(unsigned i = static_cast<unsigned>(polynomial.size()); i > 0; --i) {
    long double row = 0.0L;
    for(unsigned j = static_cast<unsigned>(polynomial[i-1U].size()); j > 0; --j) row = std::fma(row,second,polynomial[i-1U][j-1U]);
    value = std::fma(value,first,row);
  }
  return value;
}

inline long double SurfaceGradientNorm(const std::vector<std::vector<std::vector<long double>>>& gradient, long double first, long double second) {
  long double norm = 0.0L;
  for(unsigned component = 0; component < gradient.size(); ++component) norm = std::hypot(norm,SurfacePolynomialValue(gradient[component],first,second));
  return norm;
}

inline long double SurfaceGradientBound(const std::vector<std::vector<std::vector<long double>>>& gradient) {
  long double norm = 0.0L;
  for(unsigned component = 0; component < gradient.size(); ++component) {
    long double bound = 0.0L;
    for(unsigned i = 0; i < gradient[component].size(); ++i) {
      for(unsigned j = 0; j < gradient[component][i].size(); ++j) bound += std::abs(gradient[component][i][j]);
    }
    norm = std::hypot(norm,bound);
  }
  return norm;
}

inline void SurfaceGradientTranspose(const std::vector<std::vector<std::vector<long double>>>& input, std::vector<std::vector<std::vector<long double>>>& output) {
  output.resize(input.size());
  for(unsigned component = 0; component < input.size(); ++component) {
    unsigned columns = 0;
    for(unsigned i = 0; i < input[component].size(); ++i) columns = std::max(columns,static_cast<unsigned>(input[component][i].size()));
    output[component].assign(columns,std::vector<long double>(input[component].size(),0.0L));
    for(unsigned i = 0; i < input[component].size(); ++i) {
      for(unsigned j = 0; j < input[component][i].size(); ++j) output[component][j][i] = input[component][i][j];
    }
  }
}

inline long double SurfaceBinomial(unsigned n, unsigned k) {
  long double value = 1.0L;
  for(unsigned i = 1; i <= k; ++i) value *= static_cast<long double>(n-k+i)/i;
  return value;
}

inline long double SurfacePower(long double value, unsigned exponent) {
  long double result = 1.0L;
  for(unsigned i = 0; i < exponent; ++i) result *= value;
  return result;
}

inline void SurfaceGradientAffine(const std::vector<std::vector<std::vector<long double>>>& input, long double firstConstant, long double firstFirst, long double firstSecond, long double secondConstant, long double secondFirst, long double secondSecond, std::vector<std::vector<std::vector<long double>>>& output) {
  output.resize(input.size());
  for(unsigned component = 0; component < input.size(); ++component) {
    unsigned size = 1;
    for(unsigned i = 0; i < input[component].size(); ++i) size = std::max(size,i+static_cast<unsigned>(input[component][i].size()));
    output[component].assign(size,std::vector<long double>(size,0.0L));
    for(unsigned i = 0; i < input[component].size(); ++i) {
      for(unsigned j = 0; j < input[component][i].size(); ++j) {
        if(input[component][i][j] == 0.0L) continue;
        for(unsigned a = 0; a <= i; ++a) {
          for(unsigned b = 0; b <= i-a; ++b) {
            const long double first = SurfaceBinomial(i,a)*SurfaceBinomial(i-a,b)*SurfacePower(firstConstant,i-a-b)*SurfacePower(firstFirst,a)*SurfacePower(firstSecond,b);
            if(first == 0.0L) continue;
            for(unsigned c = 0; c <= j; ++c) {
              for(unsigned d = 0; d <= j-c; ++d) output[component][a+c][b+d] += input[component][i][j]*first*SurfaceBinomial(j,c)*SurfaceBinomial(j-c,d)*SurfacePower(secondConstant,j-c-d)*SurfacePower(secondFirst,c)*SurfacePower(secondSecond,d);
            }
          }
        }
      }
    }
  }
}

inline void SurfaceGradientDiagonal(const std::vector<std::vector<std::vector<long double>>>& input, std::vector<std::vector<std::vector<long double>>>& output) {
  output.resize(input.size());
  for(unsigned component = 0; component < input.size(); ++component) {
    unsigned size = 1;
    for(unsigned i = 0; i < input[component].size(); ++i) size = std::max(size,i+static_cast<unsigned>(input[component][i].size()));
    output[component].assign(size,std::vector<long double>(size,0.0L));
    for(unsigned i = 0; i < input[component].size(); ++i) {
      for(unsigned j = 0; j < input[component][i].size(); ++j) {
        for(unsigned a = 0; a <= i; ++a) {
          for(unsigned b = 0; b <= j; ++b) output[component][j-b][a+b] += (((a+b)%2U != 0) ? -1.0L : 1.0L)*SurfaceBinomial(i,a)*SurfaceBinomial(j,b)*input[component][i][j];
        }
      }
    }
  }
}

inline void SurfaceGradientFromPolynomial(const std::vector<std::vector<double>>& coefficient, long double firstFactor, long double secondFactor, std::vector<std::vector<std::vector<long double>>>& gradient) {
  unsigned columns = 0;
  for(unsigned i = 0; i < coefficient.size(); ++i) columns = std::max(columns,static_cast<unsigned>(coefficient[i].size()));
  gradient.assign(2,std::vector<std::vector<long double>>(coefficient.size(),std::vector<long double>(columns,0.0L)));
  for(unsigned i = 0; i < coefficient.size(); ++i) {
    for(unsigned j = 0; j < coefficient[i].size(); ++j) {
      if(i > 0) gradient[0][i-1U][j] = firstFactor*i*static_cast<long double>(coefficient[i][j]);
      if(j > 0) gradient[1][i][j-1U] = secondFactor*j*static_cast<long double>(coefficient[i][j]);
    }
  }
}

}

#endif
