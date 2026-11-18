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

inline long double ReferenceIntervalMoment(unsigned degree) {
  return degree%2U == 0 ? 2.0L/(degree+1U) : 0.0L;
}

inline long double ReferenceTriangleMoment(unsigned first, unsigned second) {
  long double result = 1.0L/((first+second+1.0L)*(first+second+2.0L));
  for(unsigned i = 1; i <= second; ++i) result *= static_cast<long double>(i)/(first+i);
  return result;
}

template <class Type> void ReferencePolynomialAffine2D(const std::vector<std::vector<Type>>& input, long double firstOrigin, long double firstScale, long double secondOrigin, long double secondScale, std::vector<std::vector<long double>>& output) {
  unsigned columns = 0;
  for(unsigned i = 0; i < input.size(); ++i) columns = std::max(columns,static_cast<unsigned>(input[i].size()));
  std::vector<std::vector<long double>> result(input.size(),std::vector<long double>(columns,0.0L));
  std::vector<std::vector<long double>> correction(input.size(),std::vector<long double>(columns,0.0L));
  for(unsigned i = 0; i < input.size(); ++i) {
    for(unsigned j = 0; j < input[i].size(); ++j) {
      for(unsigned a = 0; a <= i; ++a) {
        const long double first = ReferenceBinomial(i,a)*ReferencePower(firstOrigin,i-a)*ReferencePower(firstScale,a);
        for(unsigned b = 0; b <= j; ++b) {
          const long double term = static_cast<long double>(input[i][j])*first*ReferenceBinomial(j,b)*ReferencePower(secondOrigin,j-b)*ReferencePower(secondScale,b);
          const long double adjusted = term-correction[a][b];
          const long double next = result[a][b]+adjusted;
          correction[a][b] = (next-result[a][b])-adjusted;
          result[a][b] = next;
        }
      }
    }
  }
  output.swap(result);
}

inline void ReferenceAffineMoments2D(unsigned degree, const std::vector<long double>& chartValues, const std::vector<long double>& chartErrors, long double firstOrigin, long double firstScale, long double secondOrigin, long double secondScale, long double measure, std::vector<long double>& values, std::vector<long double>& errors) {
  const unsigned count = (degree+1U)*(degree+2U)/2U;
  if(chartValues.size() != 2U*count || chartErrors.size() != 2U*count) throw std::invalid_argument("Invalid reference moment table.");
  values.assign(2U*count,0.0L);
  errors.assign(2U*count,0.0L);
  for(unsigned i = 0; i <= degree; ++i) {
    for(unsigned j = 0; j <= degree-i; ++j) {
      const unsigned target = ReferenceMomentIndex2D(degree,i,j);
      long double correction[2] = {0.0L,0.0L};
      long double absoluteSum[2] = {0.0L,0.0L};
      for(unsigned a = 0; a <= i; ++a) {
        const long double first = measure*ReferenceBinomial(i,a)*ReferencePower(firstOrigin,i-a)*ReferencePower(firstScale,a);
        for(unsigned b = 0; b <= j; ++b) {
          const unsigned source = ReferenceMomentIndex2D(degree,a,b);
          const long double multiplier = first*ReferenceBinomial(j,b)*ReferencePower(secondOrigin,j-b)*ReferencePower(secondScale,b);
          for(unsigned pair = 0; pair < 2; ++pair) {
            const unsigned destination = pair*count+target;
            const long double term = multiplier*chartValues[pair*count+source];
            const long double adjusted = term-correction[pair];
            const long double next = values[destination]+adjusted;
            correction[pair] = (next-values[destination])-adjusted;
            values[destination] = next;
            errors[destination] += std::abs(multiplier)*chartErrors[pair*count+source];
            absoluteSum[pair] += std::abs(term);
          }
        }
      }
      for(unsigned pair = 0; pair < 2; ++pair) errors[pair*count+target] += (4.0L+2.0L*(i+1U)*(j+1U))*std::numeric_limits<long double>::epsilon()*absoluteSum[pair];
    }
  }
}

inline void ReferenceTensorMoments2D(unsigned degree, const std::vector<long double>& chartValues, const std::vector<long double>& chartErrors, std::vector<long double>& values, std::vector<long double>& errors, long double measure = 4.0L) {
  ReferenceAffineMoments2D(degree,chartValues,chartErrors,-1.0L,2.0L,-1.0L,2.0L,measure,values,errors);
}

inline void Quad9ReferenceCoefficientsFromNodalValues(const std::vector<double>& nodalValues, std::vector<std::vector<double>>& coefficient) {
  if(nodalValues.size() != 9) throw std::invalid_argument("QUAD9 requires nine nodal values in FEMuS order.");
  for(unsigned i = 0; i < nodalValues.size(); ++i) if(!std::isfinite(nodalValues[i])) throw std::invalid_argument("A nodal value is not finite.");
  const unsigned node[3][3] = {{0,3,7},{1,2,5},{4,6,8}};
  const long double interpolation[3][3] = {{0.0L,0.0L,1.0L},{-0.5L,0.5L,0.0L},{0.5L,0.5L,-1.0L}};
  long double intermediate[3][3] = {};
  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned y = 0; y < 3; ++y) {
      for(unsigned x = 0; x < 3; ++x) intermediate[i][y] += interpolation[i][x]*nodalValues[node[x][y]];
    }
  }
  coefficient.assign(3,std::vector<double>(3,0.0));
  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = 0; j < 3; ++j) {
      long double value = 0.0L;
      for(unsigned y = 0; y < 3; ++y) value += interpolation[j][y]*intermediate[i][y];
      coefficient[i][j] = static_cast<double>(value);
      if(!std::isfinite(coefficient[i][j])) throw std::overflow_error("A QUAD9 coefficient exceeds double range.");
    }
  }
}

inline void Tri7ReferenceCoefficientsFromNodalValues(const std::vector<double>& nodalValues, std::vector<std::vector<double>>& coefficient) {
  if(nodalValues.size() != 7) throw std::invalid_argument("TRI7 requires seven nodal values in FEMuS order.");
  for(unsigned i = 0; i < nodalValues.size(); ++i) if(!std::isfinite(nodalValues[i])) throw std::invalid_argument("A nodal value is not finite.");
  const long double bubble = 3.0L*(static_cast<long double>(nodalValues[0])+nodalValues[1]+nodalValues[2])-12.0L*(static_cast<long double>(nodalValues[3])+nodalValues[4]+nodalValues[5])+27.0L*nodalValues[6];
  coefficient.assign(3,std::vector<double>(3,0.0));
  coefficient[0][0] = nodalValues[0];
  coefficient[1][0] = static_cast<double>(-3.0L*nodalValues[0]-nodalValues[1]+4.0L*nodalValues[3]);
  coefficient[0][1] = static_cast<double>(-3.0L*nodalValues[0]-nodalValues[2]+4.0L*nodalValues[5]);
  coefficient[2][0] = static_cast<double>(2.0L*nodalValues[0]+2.0L*nodalValues[1]-4.0L*nodalValues[3]);
  coefficient[0][2] = static_cast<double>(2.0L*nodalValues[0]+2.0L*nodalValues[2]-4.0L*nodalValues[5]);
  coefficient[1][1] = static_cast<double>(4.0L*nodalValues[0]-4.0L*nodalValues[3]+4.0L*nodalValues[4]-4.0L*nodalValues[5]+bubble);
  coefficient[2][1] = static_cast<double>(-bubble);
  coefficient[1][2] = static_cast<double>(-bubble);
  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = 0; j < 3; ++j) if(!std::isfinite(coefficient[i][j])) throw std::overflow_error("A TRI7 coefficient exceeds double range.");
  }
}

}

#endif
