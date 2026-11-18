#ifndef FEMUS_REFERENCE_SURFACE_MOMENTS_HPP
#define FEMUS_REFERENCE_SURFACE_MOMENTS_HPP

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace FELevelSetMoments {

inline void ValidateSurfaceGeometryCoordinates(const std::vector<std::vector<double>>& coordinates, unsigned dimension, unsigned nodeCount) {
  if(dimension < 1 || dimension > 3 || coordinates.size() != dimension || nodeCount == 0) throw std::invalid_argument("Invalid surface geometry dimensions.");
  for(unsigned component = 0; component < dimension; ++component) {
    if(coordinates[component].size() != nodeCount) throw std::invalid_argument("Invalid surface geometry node count.");
    for(unsigned node = 0; node < nodeCount; ++node) {
      if(!std::isfinite(coordinates[component][node])) throw std::invalid_argument("A surface geometry coordinate is not finite.");
    }
  }
}

inline unsigned SurfacePolynomialCubeSize(unsigned count) {
  if(count == 0) throw std::invalid_argument("A surface polynomial cube is empty.");
  unsigned side = 1;
  while(static_cast<unsigned long long>(side)*side*side < count) ++side;
  if(static_cast<unsigned long long>(side)*side*side != count) throw std::invalid_argument("Invalid surface polynomial cube size.");
  return side;
}

inline void SurfacePolynomialTrim2D(std::vector<std::vector<long double>>& polynomial) {
  unsigned rows = 1;
  unsigned columns = 1;
  for(unsigned i = 0; i < polynomial.size(); ++i) {
    for(unsigned j = 0; j < polynomial[i].size(); ++j) {
      if(polynomial[i][j] != 0.0L) {
        rows = std::max(rows,i+1U);
        columns = std::max(columns,j+1U);
      }
    }
  }
  polynomial.resize(rows);
  for(unsigned i = 0; i < rows; ++i) polynomial[i].resize(columns,0.0L);
}

inline std::vector<std::vector<long double>> SurfacePolynomialDerivative2D(const std::vector<std::vector<long double>>& polynomial, unsigned axis, long double factor) {
  if(axis > 1) throw std::invalid_argument("Invalid surface polynomial derivative axis.");
  unsigned columns = 0;
  for(unsigned i = 0; i < polynomial.size(); ++i) columns = std::max(columns,static_cast<unsigned>(polynomial[i].size()));
  std::vector<std::vector<long double>> derivative(polynomial.size(),std::vector<long double>(columns,0.0L));
  for(unsigned i = 0; i < polynomial.size(); ++i) {
    for(unsigned j = 0; j < polynomial[i].size(); ++j) {
      const unsigned exponent = axis == 0 ? i : j;
      if(exponent != 0) derivative[i-(axis == 0 ? 1U : 0U)][j-(axis == 1 ? 1U : 0U)] = factor*exponent*polynomial[i][j];
    }
  }
  SurfacePolynomialTrim2D(derivative);
  return derivative;
}

inline void SurfacePolynomialAccumulateProduct2D(const std::vector<std::vector<long double>>& first, const std::vector<std::vector<long double>>& second, long double factor, std::vector<std::vector<long double>>& result) {
  if(first.empty() || second.empty() || factor == 0.0L) return;
  unsigned firstColumns = 0;
  unsigned secondColumns = 0;
  for(unsigned i = 0; i < first.size(); ++i) firstColumns = std::max(firstColumns,static_cast<unsigned>(first[i].size()));
  for(unsigned i = 0; i < second.size(); ++i) secondColumns = std::max(secondColumns,static_cast<unsigned>(second[i].size()));
  if(firstColumns == 0 || secondColumns == 0) return;
  const unsigned rows = std::max(static_cast<unsigned>(result.size()),static_cast<unsigned>(first.size()+second.size()-1U));
  const unsigned columns = firstColumns+secondColumns-1U;
  result.resize(rows);
  for(unsigned i = 0; i < rows; ++i) result[i].resize(std::max(static_cast<unsigned>(result[i].size()),columns),0.0L);
  for(unsigned i = 0; i < first.size(); ++i) {
    for(unsigned j = 0; j < first[i].size(); ++j) {
      if(first[i][j] == 0.0L) continue;
      const long double coefficient = factor*first[i][j];
      for(unsigned k = 0; k < second.size(); ++k) {
        for(unsigned l = 0; l < second[k].size(); ++l) {
          if(second[k][l] != 0.0L) result[i+k][j+l] = std::fma(coefficient,second[k][l],result[i+k][j+l]);
        }
      }
    }
  }
}

inline void SurfaceGradientCofactor2D(std::vector<std::vector<std::vector<long double>>>& gradient, const std::vector<std::vector<std::vector<long double>>>& coordinatePolynomials, long double firstFactor, long double secondFactor) {
  if(gradient.size() != 2 || coordinatePolynomials.size() != 2) throw std::invalid_argument("A planar surface cofactor requires two components.");
  const std::vector<std::vector<long double>> firstFirst = SurfacePolynomialDerivative2D(coordinatePolynomials[0],0,firstFactor);
  const std::vector<std::vector<long double>> firstSecond = SurfacePolynomialDerivative2D(coordinatePolynomials[0],1,secondFactor);
  const std::vector<std::vector<long double>> secondFirst = SurfacePolynomialDerivative2D(coordinatePolynomials[1],0,firstFactor);
  const std::vector<std::vector<long double>> secondSecond = SurfacePolynomialDerivative2D(coordinatePolynomials[1],1,secondFactor);
  std::vector<std::vector<std::vector<long double>>> transformed(2);
  SurfacePolynomialAccumulateProduct2D(secondSecond,gradient[0],1.0L,transformed[0]);
  SurfacePolynomialAccumulateProduct2D(secondFirst,gradient[1],-1.0L,transformed[0]);
  SurfacePolynomialAccumulateProduct2D(firstFirst,gradient[1],1.0L,transformed[1]);
  SurfacePolynomialAccumulateProduct2D(firstSecond,gradient[0],-1.0L,transformed[1]);
  unsigned side = 1;
  for(unsigned component = 0; component < 2; ++component) {
    SurfacePolynomialTrim2D(transformed[component]);
    side = std::max(side,static_cast<unsigned>(transformed[component].size()));
    for(unsigned i = 0; i < transformed[component].size(); ++i) side = std::max(side,static_cast<unsigned>(transformed[component][i].size()));
  }
  for(unsigned component = 0; component < 2; ++component) {
    transformed[component].resize(side);
    for(unsigned i = 0; i < side; ++i) transformed[component][i].resize(side,0.0L);
  }
  gradient.swap(transformed);
}

inline void SurfacePolynomialResizeCube(std::vector<long double>& polynomial, unsigned side) {
  if(side == 0) throw std::invalid_argument("A surface polynomial cube is empty.");
  const unsigned previousSide = polynomial.empty() ? 0U : SurfacePolynomialCubeSize(static_cast<unsigned>(polynomial.size()));
  if(previousSide == side) return;
  std::vector<long double> result(static_cast<std::size_t>(side)*side*side,0.0L);
  const unsigned retained = std::min(previousSide,side);
  for(unsigned i = 0; i < retained; ++i) {
    for(unsigned j = 0; j < retained; ++j) {
      for(unsigned k = 0; k < retained; ++k) result[side*side*i+side*j+k] = polynomial[previousSide*previousSide*i+previousSide*j+k];
    }
  }
  polynomial.swap(result);
}

inline void SurfacePolynomialTrim3D(std::vector<long double>& polynomial) {
  const unsigned side = SurfacePolynomialCubeSize(static_cast<unsigned>(polynomial.size()));
  unsigned retained = 1;
  for(unsigned i = 0; i < side; ++i) {
    for(unsigned j = 0; j < side; ++j) {
      for(unsigned k = 0; k < side; ++k) {
        if(polynomial[side*side*i+side*j+k] != 0.0L) retained = std::max(retained,std::max(i,std::max(j,k))+1U);
      }
    }
  }
  SurfacePolynomialResizeCube(polynomial,retained);
}

inline std::vector<long double> SurfacePolynomialDerivative3D(const std::vector<long double>& polynomial, unsigned axis, long double factor) {
  if(axis > 2) throw std::invalid_argument("Invalid surface polynomial derivative axis.");
  const unsigned side = SurfacePolynomialCubeSize(static_cast<unsigned>(polynomial.size()));
  std::vector<long double> derivative(polynomial.size(),0.0L);
  for(unsigned i = 0; i < side; ++i) {
    for(unsigned j = 0; j < side; ++j) {
      for(unsigned k = 0; k < side; ++k) {
        const unsigned exponent = axis == 0 ? i : (axis == 1 ? j : k);
        if(exponent != 0) derivative[side*side*(i-(axis == 0 ? 1U : 0U))+side*(j-(axis == 1 ? 1U : 0U))+k-(axis == 2 ? 1U : 0U)] = factor*exponent*polynomial[side*side*i+side*j+k];
      }
    }
  }
  SurfacePolynomialTrim3D(derivative);
  return derivative;
}

inline void SurfacePolynomialAccumulateProduct3D(const std::vector<long double>& first, const std::vector<long double>& second, long double factor, std::vector<long double>& result) {
  const unsigned firstSide = SurfacePolynomialCubeSize(static_cast<unsigned>(first.size()));
  const unsigned secondSide = SurfacePolynomialCubeSize(static_cast<unsigned>(second.size()));
  const unsigned previousSide = result.empty() ? 0U : SurfacePolynomialCubeSize(static_cast<unsigned>(result.size()));
  const unsigned side = std::max(previousSide,firstSide+secondSide-1U);
  SurfacePolynomialResizeCube(result,side);
  for(unsigned firstIndex = 0; firstIndex < first.size(); ++firstIndex) {
    if(first[firstIndex] == 0.0L || factor == 0.0L) continue;
    const unsigned firstI = firstIndex/(firstSide*firstSide);
    const unsigned firstJ = (firstIndex/firstSide)%firstSide;
    const unsigned firstK = firstIndex%firstSide;
    const long double coefficient = factor*first[firstIndex];
    for(unsigned secondIndex = 0; secondIndex < second.size(); ++secondIndex) {
      if(second[secondIndex] == 0.0L) continue;
      const unsigned secondI = secondIndex/(secondSide*secondSide);
      const unsigned secondJ = (secondIndex/secondSide)%secondSide;
      const unsigned secondK = secondIndex%secondSide;
      const unsigned target = side*side*(firstI+secondI)+side*(firstJ+secondJ)+firstK+secondK;
      result[target] = std::fma(coefficient,second[secondIndex],result[target]);
    }
  }
}

inline void SurfaceGradientCofactor3D(std::vector<std::vector<long double>>& gradient, const std::vector<std::vector<long double>>& coordinatePolynomials, long double derivativeFactor) {
  if(gradient.size() != 3 || coordinatePolynomials.size() != 3) throw std::invalid_argument("A spatial surface cofactor requires three components.");
  std::vector<std::vector<std::vector<long double>>> jacobian(3,std::vector<std::vector<long double>>(3));
  for(unsigned physical = 0; physical < 3; ++physical) {
    SurfacePolynomialCubeSize(static_cast<unsigned>(gradient[physical].size()));
    for(unsigned reference = 0; reference < 3; ++reference) jacobian[physical][reference] = SurfacePolynomialDerivative3D(coordinatePolynomials[physical],reference,derivativeFactor);
  }
  std::vector<std::vector<long double>> transformed(3);
  for(unsigned physical = 0; physical < 3; ++physical) {
    for(unsigned reference = 0; reference < 3; ++reference) {
      const unsigned firstPhysical = (physical+1U)%3U;
      const unsigned secondPhysical = (physical+2U)%3U;
      const unsigned firstReference = (reference+1U)%3U;
      const unsigned secondReference = (reference+2U)%3U;
      std::vector<long double> minor;
      SurfacePolynomialAccumulateProduct3D(jacobian[firstPhysical][firstReference],jacobian[secondPhysical][secondReference],1.0L,minor);
      SurfacePolynomialAccumulateProduct3D(jacobian[firstPhysical][secondReference],jacobian[secondPhysical][firstReference],-1.0L,minor);
      SurfacePolynomialTrim3D(minor);
      SurfacePolynomialAccumulateProduct3D(minor,gradient[reference],1.0L,transformed[physical]);
    }
  }
  unsigned side = 1;
  for(unsigned component = 0; component < 3; ++component) {
    SurfacePolynomialTrim3D(transformed[component]);
    side = std::max(side,SurfacePolynomialCubeSize(static_cast<unsigned>(transformed[component].size())));
  }
  for(unsigned component = 0; component < 3; ++component) SurfacePolynomialResizeCube(transformed[component],side);
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

inline long double SurfaceGradientFactor(const std::vector<std::vector<std::vector<long double>>>& gradient, long double first, long double second, bool inverseSurface, int normalComponent = -1) {
  const long double norm = SurfaceGradientNorm(gradient,first,second);
  if(!inverseSurface && normalComponent < 0) return norm;
  if(!(norm > 0) || !std::isfinite(norm)) throw std::domain_error("An inverse surface moment requires a regular interface and geometry.");
  if(normalComponent >= 0) {
    if(static_cast<unsigned>(normalComponent) >= gradient.size()) throw std::invalid_argument("Invalid normal-weight component.");
    return SurfacePolynomialValue(gradient[normalComponent],first,second)/norm;
  }
  return 1.0L/norm;
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

}

#endif
