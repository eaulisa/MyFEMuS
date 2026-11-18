#include "UnitElementMoments.hpp"
#include "Basis.hpp"
#include <boost/math/quadrature/gauss.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

template<unsigned Order> std::vector<long double> CurvedArcReference(unsigned degree) {
  const long double pi = std::acos(-1.0L);
  std::vector<long double> result((degree+1U)*(degree+2U)/2U,0);
  for(unsigned q = 0; q < boost::math::quadrature::gauss<long double,Order>::abscissa().size(); ++q) for(unsigned side = 0; side < 2; ++side) {
    const long double angle = pi/4+(side == 0 ? -pi/4 : pi/4)*boost::math::quadrature::gauss<long double,Order>::abscissa()[q];
    const long double sine = std::sin(angle);
    const long double cosine = std::cos(angle);
    const long double weight = pi/16*boost::math::quadrature::gauss<long double,Order>::weights()[q]*std::sqrt(1-0.4L*cosine+0.04L);
    for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) result[FELevelSetMoments::UnitMomentIndex(2,i,j)] += weight*FELevelSetMoments::ReferencePower(cosine,i)*FELevelSetMoments::ReferencePower(sine,j);
  }
  return result;
}

long double CurvedQuadrantMoment(unsigned i, unsigned j, unsigned corner, const std::vector<long double>& base) {
  long double result = 0;
  for(unsigned a = (corner&1U) != 0 ? 0U : i; a <= i; ++a) for(unsigned b = (corner&2U) != 0 ? 0U : j; b <= j; ++b) {
    long double factor = FELevelSetMoments::ReferenceBinomial(i,a)*FELevelSetMoments::ReferenceBinomial(j,b);
    if((((corner&1U) != 0 ? a : 0U)+((corner&2U) != 0 ? b : 0U))%2U != 0) factor = -factor;
    result += factor*base[FELevelSetMoments::UnitMomentIndex(2,a,b)];
  }
  return result;
}

int main() {
  const std::vector<double> sourceNodal = {-0.06250000000000003, 0.0, 0.0625, -2.7755575615628914e-17, -0.015625000000000028, 0.046875, 0.04687499999999999, -0.01562500000000004, 0.03125};
  const std::vector<std::vector<double>> sourceCoordinates = {{-0.25, -0.025000000000000005, 0.0, -0.22500000000000003, -0.14375, -0.006250000000000002, -0.11875000000000001, -0.23125, -0.125}, {-0.2, -0.25, 0.0, 0.0, -0.22500000000000003, -0.125, 0.0, -0.10000000000000002, -0.11249999999999999}};
  const unsigned degree = 10;
  const unsigned count = (degree+1U)*(degree+2U)/2U;
  const std::vector<long double> physical = CurvedArcReference<128>(degree);
  const std::vector<long double> previous = CurvedArcReference<64>(degree);
  std::vector<long double> volume(count);
  std::vector<long double> surface(count);
  for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) {
    const unsigned n = FELevelSetMoments::UnitMomentIndex(2,i,j);
    const long double numerator = std::tgamma((i+1.0L)/2)*std::tgamma((j+1.0L)/2);
    volume[n] = numerator/(4*std::tgamma((i+j+4.0L)/2));
    surface[n] = numerator/(2*std::tgamma((i+j+2.0L)/2));
    if(std::abs(physical[n]-previous[n]) > 2e-18L) throw std::runtime_error("The independent curved-arc reference has not converged.");
  }
  femus::QuadBiquadratic basis;
  unsigned cases = 0;
  unsigned checks = count;
  long double maximumVolumeError = 0;
  long double maximumSurfaceError = 0;
  for(unsigned flip = 0; flip < 4; ++flip) for(unsigned metric = 0; metric < 3; ++metric) for(const long double tolerance : {1e-10L,1e-12L}) {
    std::vector<double> nodal(9);
    std::vector<std::vector<double>> coordinates(2,std::vector<double>(9));
    for(unsigned node = 0; node < 9; ++node) {
      const double* point = basis.GetXcoarse(node);
      unsigned source = 9;
      for(unsigned candidate = 0; candidate < 9; ++candidate) {
        const double* original = basis.GetXcoarse(candidate);
        bool match = true;
        for(unsigned axis = 0; axis < 2; ++axis) match = match && original[axis] == ((flip&(1U<<axis)) != 0 ? -point[axis] : point[axis]);
        if(match) source = candidate;
      }
      if(source == 9) throw std::runtime_error("Invalid reflected QUAD9 node.");
      nodal[node] = sourceNodal[source];
      for(unsigned axis = 0; axis < 2; ++axis) coordinates[axis][node] = sourceCoordinates[axis][source];
    }
    std::vector<long double> values;
    std::vector<long double> errors;
    FELevelSetMoments::IntegrateUnitTensorNodalMoments(nodal,2,degree,values,errors,tolerance,tolerance,metric != 0,metric == 2 ? &coordinates : nullptr,0,0);
    for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) {
      const unsigned n = FELevelSetMoments::UnitMomentIndex(2,i,j);
      const long double expectedH = CurvedQuadrantMoment(i,j,3U^flip,volume);
      const long double expectedD = CurvedQuadrantMoment(i,j,3U^flip,metric == 2 ? physical : surface)*(metric == 0 ? 8 : 1);
      if(!std::isfinite(values[n]) || !std::isfinite(values[count+n]) || std::abs(values[n]-expectedH) > 8*tolerance || std::abs(values[count+n]-expectedD) > 8*tolerance) throw std::runtime_error("The curved QUAD9 regression failed its independent reference.");
      maximumVolumeError = std::max(maximumVolumeError,std::abs(values[n]-expectedH));
      if(metric == 2) maximumSurfaceError = std::max(maximumSurfaceError,std::abs(values[count+n]-expectedD));
      checks += 2;
    }
    ++cases;
  }
  std::cout << std::setprecision(8) << cases << " curved QUAD9 cases and " << checks << " checks passed; maximum reference area error " << maximumVolumeError << ", physical curve error " << maximumSurfaceError << '\n';
}
