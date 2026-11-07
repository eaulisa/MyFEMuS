#include "Quad9DirectionSelection.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

long double ScaledPower(long double value, unsigned degree) {
  long double result = 1.0L;
  for(unsigned i = 0; i < degree; ++i) result *= value;
  return result;
}

long double ScaledBinomial(unsigned degree, unsigned index) {
  long double result = 1.0L;
  for(unsigned i = 1; i <= index; ++i) result *= static_cast<long double>(degree-index+i)/i;
  return result;
}

long double ScaledInterval(unsigned degree) {
  return degree%2U == 0 ? 2.0L/(degree+1U) : 0.0L;
}

int main() {
  unsigned checks = 0;
  long double maximumNormalizedError = 0.0L;
  const double scale = std::ldexp(1.0,-20);
  const std::vector<std::vector<double>> coordinates = {{-1,1,1,-1,0,1,0,-1,0},{-1,-1,1,1,-1,0,1,0,0}};
  for(unsigned degree = 2; degree <= 5; degree += 3) {
    for(unsigned field = 0; field < 2; ++field) {
      for(unsigned mode = 0; mode < (field == 0 ? 3U : 1U); ++mode) {
        std::vector<std::vector<double>> coefficient(3,std::vector<double>(3,0.0));
        coefficient[0][0] = -0.125*scale;
        coefficient[1][0] = scale;
        coefficient[0][2] = field == 0 ? 0.0 : -0.25*scale;
        std::vector<std::vector<double>> positive;
        std::vector<std::vector<double>> delta;
        FELevelSetMoments::IntegrateQuad9ReferenceSelected(coefficient,degree,positive,delta,1.0e-12,1.0e-12,true,false,28,0,200000,nullptr,0,nullptr,false,mode != 0,nullptr,nullptr,mode == 2 ? &coordinates : nullptr);
        for(unsigned i = 0; i <= degree; ++i) {
          for(unsigned j = 0; j <= degree-i; ++j) {
            long double graphH = 0.0L;
            long double graphD = 0.0L;
            for(unsigned k = 0; k <= i+1U; ++k) graphH += ScaledBinomial(i+1U,k)*ScaledPower(0.125L,i+1U-k)*ScaledPower(field == 0 ? 0.0L : 0.25L,k)*ScaledInterval(j+2U*k);
            for(unsigned k = 0; k <= i; ++k) graphD += ScaledBinomial(i,k)*ScaledPower(0.125L,i-k)*ScaledPower(field == 0 ? 0.0L : 0.25L,k)*ScaledInterval(j+2U*k);
            const long double expected[2] = {(ScaledInterval(j)-graphH)/(i+1U),graphD/(mode == 0 ? scale : 1.0L)};
            const double actual[2] = {positive[i][j],delta[i][j]};
            for(unsigned pair = 0; pair < 2; ++pair) {
              const long double normalizedError = std::abs(static_cast<long double>(actual[pair])-expected[pair])/(1.0e-12L+1.0e-12L*std::abs(expected[pair]));
              maximumNormalizedError = std::max(maximumNormalizedError,normalizedError);
              if(!std::isfinite(actual[pair]) || normalizedError > 1.0L) throw std::runtime_error("A scaled native moment exceeds requested tolerance.");
              ++checks;
            }
          }
        }
      }
    }
  }
  std::cout << checks << " scaled plane/graph checks passed; maximum error fraction of requested budget " << maximumNormalizedError << '\n';
}
