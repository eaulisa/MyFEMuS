#include "UnitElementMoments.hpp"
#include "Basis.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

long double OctantMoment(unsigned i, unsigned j, unsigned k, unsigned corner, bool surface) {
  long double value = 0;
  for(unsigned a = (corner&1U) != 0 ? 0U : i; a <= i; ++a) for(unsigned b = (corner&2U) != 0 ? 0U : j; b <= j; ++b) for(unsigned c = (corner&4U) != 0 ? 0U : k; c <= k; ++c) {
    long double factor = FELevelSetMoments::ReferenceBinomial(i,a)*FELevelSetMoments::ReferenceBinomial(j,b)*FELevelSetMoments::ReferenceBinomial(k,c);
    if((((corner&1U) != 0 ? a : 0U)+((corner&2U) != 0 ? b : 0U)+((corner&4U) != 0 ? c : 0U))%2U != 0) factor = -factor;
    value += factor*std::tgamma((a+1.0L)/2)*std::tgamma((b+1.0L)/2)*std::tgamma((c+1.0L)/2)/((surface ? 4 : 8)*std::tgamma((a+b+c+(surface ? 3.0L : 5.0L))/2));
  }
  return value;
}

void OctantCheck(long double actual, long double expected, long double tolerance) {
  if(!std::isfinite(actual) || std::abs(actual-expected) > 8*tolerance) throw std::runtime_error("A boundary-tangent sphere moment failed its analytical reference.");
}

void OctantRejectionChecks() {
  const std::vector<std::vector<long double>> isolated = {{0,0,1},{0},{1}};
  const std::vector<std::vector<long double>> saddle = {{0,0,-1},{0},{1}};
  const std::vector<std::vector<long double>> normal = {{1}};
  const std::vector<std::vector<long double>> zeroNormal = {{0}};
  if(FELevelSetMoments::EndpointProjection::ContactExponent(isolated,0,false,false,&normal) != 2) throw std::runtime_error("A regular isolated face tangency was rejected.");
  bool singularRejected = false;
  bool saddleRejected = false;
  try { FELevelSetMoments::EndpointProjection::ContactExponent(isolated,0,false,false,&zeroNormal); }
  catch(const std::runtime_error&) { singularRejected = true; }
  try { FELevelSetMoments::EndpointProjection::ContactExponent(saddle,0,false,false,&normal); }
  catch(const std::runtime_error&) { saddleRejected = true; }
  if(!singularRejected || !saddleRejected) throw std::runtime_error("An unsupported non-transverse contact was accepted.");
}

int main() {
  OctantRejectionChecks();
  femus::HexBiquadratic basis;
  const unsigned degree = 10;
  const unsigned count = (degree+1U)*(degree+2U)*(degree+3U)/6U;
  const long double radius = 0.25L;
  const std::vector<double> failing = {-0.125,-0.0625,0,-0.0625,-0.0625,0,0.0625,0,-0.078125,-0.015625,-0.015625,-0.078125,-0.015625,0.046875,0.046875,-0.015625,-0.078125,-0.015625,0.046875,-0.015625,-0.03125,0.03125,0.03125,-0.03125,-0.03125,0.03125,0.015625};
  unsigned cases = 0;
  unsigned checks = 3;
  long double maximumVolumeError = 0;
  long double maximumSurfaceError = 0;
  for(unsigned corner = 0; corner < 8; ++corner) for(unsigned metric = 0; metric < 3; ++metric) for(const long double tolerance : {1e-10L,1e-12L}) {
    std::vector<double> nodal(27);
    std::vector<std::vector<double>> coordinates(3,std::vector<double>(27));
    for(unsigned node = 0; node < 27; ++node) {
      const double* point = basis.GetXcoarse(node);
      nodal[node] = radius*radius;
      for(unsigned axis = 0; axis < 3; ++axis) {
        coordinates[axis][node] = radius*((point[axis]+1)/2-((corner&(1U<<axis)) != 0 ? 1 : 0));
        nodal[node] -= coordinates[axis][node]*coordinates[axis][node];
      }
    }
    if(corner == 7 && nodal != failing) throw std::runtime_error("The ex44_test element-6 regression data changed.");
    std::vector<long double> values;
    std::vector<long double> errors;
    FELevelSetMoments::IntegrateUnitTensorNodalMoments(nodal,3,degree,values,errors,tolerance,tolerance,metric != 0,metric == 2 ? &coordinates : nullptr,0,0);
    for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) for(unsigned k = 0; k <= degree-i-j; ++k) {
      const unsigned n = FELevelSetMoments::HexMomentIndex(i,j,k);
      const long double volume = OctantMoment(i,j,k,corner,false);
      const long double surface = OctantMoment(i,j,k,corner,true)*(metric == 0 ? 1/(2*radius*radius) : (metric == 1 ? 1 : radius*radius));
      OctantCheck(values[n],volume,tolerance);
      OctantCheck(values[count+n],surface,tolerance);
      maximumVolumeError = std::max(maximumVolumeError,radius*radius*radius*std::abs(values[n]-volume));
      if(metric == 2) maximumSurfaceError = std::max(maximumSurfaceError,std::abs(values[count+n]-surface));
      checks += 2;
    }
    ++cases;
  }
  std::cout << std::setprecision(8) << cases << " sphere-octant cases and " << checks << " checks passed; maximum physical volume error " << maximumVolumeError << ", surface error " << maximumSurfaceError << '\n';
}
