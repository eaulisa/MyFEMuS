#include "LineIntegration.hpp"
#include <iostream>
#include <stdexcept>

int main() {
  const std::vector<long double> polynomial = {3.66848065654445018278e-17L,-1.16299486905685379029e-15L,1.12631006249140947052e-15L,-2.56245806002129715809e-31L,7.42051998969827515356e-32L};
  std::vector<long double> candidates;
  LineSliceUtilities::PolynomialRoots(polynomial,0,1,candidates);
  if(candidates.empty()) throw std::runtime_error("The recorded near-boundary root fixture lost all candidates.");
  const std::vector<long double> narrow = FELevelSetMoments::QuadraticRoots(1,-1e-24L,0);
  if(narrow.size() != 2 || narrow[0] == narrow[1]) throw std::runtime_error("Distinct narrow LINE roots were merged.");
  const long double gap = std::ldexp(1.0L,-30);
  const std::vector<long double> close = FELevelSetMoments::QuadraticRoots(1,-1,0.25L-gap*gap);
  if(close.size() != 2 || !(close[0] < close[1])) throw std::runtime_error("A close-root LINE control was lost.");
  std::vector<long double> positive;
  std::vector<long double> delta;
  FELevelSetMoments::Integrate1DInternal(1,-1,0.25L-gap*gap,0,0,1,true,false,positive,delta);
  if(!std::isfinite(delta[0]) || std::abs(delta[0]-1/gap) > 1e-12L/gap) throw std::runtime_error("A retained close-root LINE control lost its delta contribution.");
  std::cout << "UNSUPPORTED legacy LINE coefficient case: a root 1e-24 from a zero endpoint lies inside the preserved endpoint-snap tolerance; no accuracy claim is made for that integral\n";
  std::cout << "Near-boundary candidate and distinct narrow-root regressions passed\n";
}
