#ifndef FEMUS_ENDPOINT_SIMPLEX_CHART_HPP
#define FEMUS_ENDPOINT_SIMPLEX_CHART_HPP

#include "ReferenceElementLS.hpp"

namespace FELevelSetMoments {

class EndpointSimplexChart {
public:
  // This barycentric permutation substitutes x=1-u-v exactly and retains every polynomial term.
  static std::vector<std::vector<long double>> Triangle(const std::vector<std::vector<long double>>& coefficient) {
    unsigned degree = 0;
    for(unsigned i = 0; i < coefficient.size(); ++i) for(unsigned j = 0; j < coefficient[i].size(); ++j) if(coefficient[i][j] != 0) degree = std::max(degree,i+j);
    std::vector<std::vector<long double>> result(std::max(3U,degree+1U),std::vector<long double>(std::max(3U,degree+1U),0));
    for(unsigned i = 0; i < coefficient.size(); ++i) for(unsigned j = 0; j < coefficient[i].size(); ++j) if(coefficient[i][j] != 0) for(unsigned a = 0; a <= i; ++a) for(unsigned b = 0; b <= i-a; ++b) result[a][b+j] += coefficient[i][j]*ReferenceBinomial(i,a)*ReferenceBinomial(i-a,b)*((a+b)%2U == 0 ? 1 : -1);
    return result;
  }

  static void RestoreTriangle(unsigned degree, const std::vector<long double>& source, const std::vector<long double>& sourceErrors, std::vector<long double>& values, std::vector<long double>& errors) {
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    values.assign(2U*count,0);
    errors.assign(2U*count,0);
    for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) for(unsigned a = 0; a <= i; ++a) for(unsigned b = 0; b <= i-a; ++b) {
      const unsigned target = ReferenceMomentIndex2D(degree,i,j);
      const unsigned origin = ReferenceMomentIndex2D(degree,a,b+j);
      const long double factor = ReferenceBinomial(i,a)*ReferenceBinomial(i-a,b)*((a+b)%2U == 0 ? 1 : -1);
      for(unsigned pair = 0; pair < 2; ++pair) { values[pair*count+target] += factor*source[pair*count+origin]; errors[pair*count+target] += std::abs(factor)*(sourceErrors[pair*count+origin]+16*std::numeric_limits<long double>::epsilon()*std::abs(source[pair*count+origin])); }
    }
  }

  static unsigned Index(unsigned i, unsigned j, unsigned k) {
    const unsigned degree = i+j+k;
    return degree*(degree+1U)*(degree+2U)/6U+(degree-i)*(degree-i+1U)/2U+k;
  }

  static void RestoreTetrahedron(unsigned degree, const std::vector<long double>& source, const std::vector<long double>& sourceErrors, std::vector<long double>& values, std::vector<long double>& errors) {
    const unsigned count = (degree+1U)*(degree+2U)*(degree+3U)/6U;
    values.assign(2U*count,0);
    errors.assign(2U*count,0);
    for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) for(unsigned k = 0; k <= degree-i-j; ++k) for(unsigned a = 0; a <= i; ++a) for(unsigned b = 0; b <= i-a; ++b) for(unsigned c = 0; c <= i-a-b; ++c) {
      const unsigned target = Index(i,j,k);
      const unsigned origin = Index(a,b+j,c+k);
      const long double factor = ReferenceBinomial(i,a)*ReferenceBinomial(i-a,b)*ReferenceBinomial(i-a-b,c)*((a+b+c)%2U == 0 ? 1 : -1);
      for(unsigned pair = 0; pair < 2; ++pair) { values[pair*count+target] += factor*source[pair*count+origin]; errors[pair*count+target] += std::abs(factor)*(sourceErrors[pair*count+origin]+16*std::numeric_limits<long double>::epsilon()*std::abs(source[pair*count+origin])); }
    }
  }
};

}

#endif
