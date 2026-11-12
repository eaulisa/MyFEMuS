#include "SurfaceEndpointIntegration.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

long double ResolutionIntegral(long double endpoint, long double delta, int mapping, unsigned power) {
  long double result = 0;
  for(unsigned n = 0; n <= power; ++n) result += delta*FELevelSetMoments::ReferenceBinomial(power,n)*FELevelSetMoments::ReferencePower(endpoint,power-n)*FELevelSetMoments::ReferencePower(mapping*delta,n)/(n+1U);
  return result;
}

unsigned ResolutionCase(unsigned exponent, unsigned ulps, int mapping, bool selective) {
  const unsigned degree = 2;
  const unsigned count = 10;
  const long double endpoint = 0.54802043156957994162L;
  const long double delta = ulps*std::abs(std::nextafter(endpoint,1.0L)-endpoint);
  if(FELevelSetMoments::OuterIntervalQuadrature::Resolved<15>({endpoint,endpoint+mapping*delta,endpoint,delta,mapping,exponent,true},0,1)) throw std::runtime_error("The coordinate-resolution regression no longer exercises a narrow strip.");
  FELevelSetMoments::SurfaceEndpointData data;
  data.coefficient.assign(27,0);
  data.coefficient[0] = 0.75L;
  data.coefficient[FELevelSetMoments::Hex27Index(0,1,0)] = -1;
  data.degree = degree;
  data.graphAxis = 1;
  data.unitCube = true;
  data.absoluteH = data.absoluteD = 1e-12L;
  data.relative = 0;
  data.powers.resize(count);
  for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) for(unsigned k = 0; k <= degree-i-j; ++k) data.powers[FELevelSetMoments::HexMomentIndex(i,j,k)] = {i,j,k};
  FELevelSetMoments::SurfaceEndpointDomain cell;
  cell.vertices = {{0,0,0},{1,1,1}};
  std::vector<unsigned char> requested(2U*count,1);
  if(selective) for(unsigned n = 0; n < requested.size(); ++n) requested[n] = n%3U == 0;
  std::vector<long double> values(2U*count,0);
  std::vector<long double> errors(2U*count,0);
  FELevelSetMoments::SurfaceEndpointIntegration::Adaptive(data,cell,0,0,1,mapping,endpoint,delta,exponent,1,0,requested,values,errors);
  for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) for(unsigned k = 0; k <= degree-i-j; ++k) {
    const unsigned n = FELevelSetMoments::HexMomentIndex(i,j,k);
    const long double outer = ResolutionIntegral(endpoint,delta,mapping,i)/(k+1U);
    const long double expected[2] = {outer*FELevelSetMoments::ReferencePower(0.75L,j+1U)/(j+1U),outer*FELevelSetMoments::ReferencePower(0.75L,j)};
    for(unsigned pair = 0; pair < 2; ++pair) {
      const unsigned index = pair*count+n;
      if(!requested[index]) { if(values[index] != 0 || errors[index] != 0) throw std::runtime_error("A narrow strip evaluated an inactive moment."); continue; }
      if(!std::isfinite(values[index]) || std::abs(values[index]-expected[pair]) > errors[index]+128*std::numeric_limits<long double>::epsilon()*std::abs(expected[pair]) || errors[index] > 1e-12L) throw std::runtime_error("A resolved endpoint rule failed its analytical integral or error estimate.");
    }
  }
  bool rejected = false;
  try { FELevelSetMoments::SurfaceEndpointIntegration::Adaptive(data,cell,0,0,1,mapping,endpoint,delta/ulps,exponent,1,0,requested,values,errors); }
  catch(const std::runtime_error&) { rejected = true; }
  if(!rejected) throw std::runtime_error("An unresolved endpoint panel was accepted.");
  return 2U*count+1U;
}

int main() {
  unsigned checks = 0;
  for(const int mapping : {-1,1}) for(const bool selective : {false,true}) {
    checks += ResolutionCase(2,64,mapping,selective);
    checks += ResolutionCase(2,16384,mapping,selective);
    checks += ResolutionCase(4,65536,mapping,selective);
    checks += ResolutionCase(6,1048576,mapping,selective);
  }
  std::cout << checks << " endpoint-coordinate, transformation-Jacobian, selective-mask and precision-rejection checks passed\n";
}
