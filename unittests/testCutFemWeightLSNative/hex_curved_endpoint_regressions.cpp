#include "UnitElementMoments.hpp"
#include "Basis.hpp"
#include <boost/math/quadrature/gauss.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

template<unsigned Order> std::vector<long double> CurvedSurfaceReference(unsigned degree) {
  const long double pi = std::acos(-1.0L);
  std::vector<std::vector<long double>> theta(degree+1U,std::vector<long double>(degree+1U,0));
  for(unsigned q = 0; q < boost::math::quadrature::gauss<long double,Order>::abscissa().size(); ++q) for(unsigned side = 0; side < 2; ++side) {
    const long double angle = pi/4+(side == 0 ? -pi/4 : pi/4)*boost::math::quadrature::gauss<long double,Order>::abscissa()[q];
    const long double sine = std::sin(angle);
    const long double cosine = std::cos(angle);
    const long double weight = pi/4*boost::math::quadrature::gauss<long double,Order>::weights()[q]*sine*(1-0.2L*cosine)*std::sqrt(1-0.4L*cosine+0.04L);
    for(unsigned m = 0; m <= degree; ++m) for(unsigned k = 0; k <= degree-m; ++k) theta[m][k] += weight*FELevelSetMoments::ReferencePower(sine,m)*FELevelSetMoments::ReferencePower(cosine,k);
  }
  std::vector<long double> result((degree+1U)*(degree+2U)*(degree+3U)/6U);
  for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) for(unsigned k = 0; k <= degree-i-j; ++k) result[FELevelSetMoments::HexMomentIndex(i,j,k)] = 0.0625L*FELevelSetMoments::ReferencePower(0.5L,i+j+k)*std::tgamma((i+1.0L)/2)*std::tgamma((j+1.0L)/2)/(2*std::tgamma((i+j+2.0L)/2))*theta[i+j][k];
  return result;
}

long double CurvedOctantMoment(unsigned i, unsigned j, unsigned k, unsigned corner, const std::vector<long double>& base) {
  long double result = 0;
  for(unsigned a = (corner&1U) != 0 ? 0U : i; a <= i; ++a) for(unsigned b = (corner&2U) != 0 ? 0U : j; b <= j; ++b) for(unsigned c = (corner&4U) != 0 ? 0U : k; c <= k; ++c) {
    long double factor = FELevelSetMoments::ReferenceBinomial(i,a)*FELevelSetMoments::ReferenceBinomial(j,b)*FELevelSetMoments::ReferenceBinomial(k,c);
    if((((corner&1U) != 0 ? a : 0U)+((corner&2U) != 0 ? b : 0U)+((corner&4U) != 0 ? c : 0U))%2U != 0) factor = -factor;
    result += factor*base[FELevelSetMoments::HexMomentIndex(a,b,c)];
  }
  return result;
}

int main() {
  const std::vector<double> sourceNodal = {-0.4375,-0.187500000000000111,-0.4375,-0.687499999999999778,-0.1875,0.0625,-0.1875,-0.4375,-0.250000000000000111,-0.250000000000000111,-0.5,-0.5,0,0,-0.25,-0.25,-0.25,-2.77555756156289135e-17,-0.25,-0.5,-0.0625,-0.0625,-0.3125,-0.3125,-0.3125,-0.0625000000000000555,-0.124999999999999986};
  const std::vector<std::vector<double>> sourceCoordinates = {{-0.3,0,0,-0.3,-0.5,0,0,-0.5,-0.15,0,-0.15,-0.3,-0.25,0,-0.25,-0.5,-0.4,0,0,-0.4,-0.2,0,-0.2,-0.4,-0.15,-0.25,-0.2},{0,0,0.3,0.3,0,0,0.5,0.5,0,0.15,0.3,0.15,0,0.25,0.5,0.25,0,0,0.4,0.4,0,0.2,0.4,0.2,0.15,0.25,0.2},{-0.5,-0.4,-0.5,-0.6,-0.1,0,-0.1,-0.2,-0.425000000000000044,-0.425000000000000044,-0.525,-0.525,-0.025,-0.025,-0.125,-0.125,-0.325,-0.225000000000000033,-0.325,-0.425,-0.25,-0.25,-0.350000000000000033,-0.350000000000000033,-0.45,-0.05,-0.274999999999999967}};
  const unsigned degree = 10;
  const unsigned count = (degree+1U)*(degree+2U)*(degree+3U)/6U;
  const std::vector<long double> physical = CurvedSurfaceReference<128>(degree);
  const std::vector<long double> previous = CurvedSurfaceReference<64>(degree);
  std::vector<long double> volume(count);
  std::vector<long double> surface(count);
  for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) for(unsigned k = 0; k <= degree-i-j; ++k) {
    const unsigned n = FELevelSetMoments::HexMomentIndex(i,j,k);
    const long double numerator = std::tgamma((i+1.0L)/2)*std::tgamma((j+1.0L)/2)*std::tgamma((k+1.0L)/2);
    volume[n] = FELevelSetMoments::ReferencePower(0.5L,i+j+k+3U)*numerator/(8*std::tgamma((i+j+k+5.0L)/2));
    surface[n] = FELevelSetMoments::ReferencePower(0.5L,i+j+k+2U)*numerator/(4*std::tgamma((i+j+k+3.0L)/2));
    if(std::abs(physical[n]-previous[n]) > 2e-18L) throw std::runtime_error("The independent curved-surface reference has not converged.");
  }
  femus::HexBiquadratic basis;
  unsigned cases = 0;
  unsigned checks = count;
  long double maximumVolumeError = 0;
  long double maximumSurfaceError = 0;
  for(unsigned flip = 0; flip < 8; ++flip) for(unsigned metric = 0; metric < 3; ++metric) for(const long double tolerance : {1e-10L,1e-12L}) {
    std::vector<double> nodal(27);
    std::vector<std::vector<double>> coordinates(3,std::vector<double>(27));
    for(unsigned node = 0; node < 27; ++node) {
      const double* point = basis.GetXcoarse(node);
      unsigned source = 27;
      for(unsigned candidate = 0; candidate < 27; ++candidate) {
        const double* original = basis.GetXcoarse(candidate);
        bool match = true;
        for(unsigned axis = 0; axis < 3; ++axis) match = match && original[axis] == ((flip&(1U<<axis)) != 0 ? -point[axis] : point[axis]);
        if(match) source = candidate;
      }
      if(source == 27) throw std::runtime_error("Invalid reflected HEX27 node.");
      nodal[node] = sourceNodal[source];
      for(unsigned axis = 0; axis < 3; ++axis) coordinates[axis][node] = sourceCoordinates[axis][source];
    }
    std::vector<long double> values;
    std::vector<long double> errors;
    FELevelSetMoments::IntegrateUnitTensorNodalMoments(nodal,3,degree,values,errors,tolerance,tolerance,metric != 0,metric == 2 ? &coordinates : nullptr,0,0);
    for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) for(unsigned k = 0; k <= degree-i-j; ++k) {
      const unsigned n = FELevelSetMoments::HexMomentIndex(i,j,k);
      const long double expectedH = CurvedOctantMoment(i,j,k,5U^flip,volume);
      const long double expectedD = CurvedOctantMoment(i,j,k,5U^flip,metric == 2 ? physical : surface)*(metric == 0 ? 4 : 1);
      if(!std::isfinite(values[n]) || !std::isfinite(values[count+n]) || std::abs(values[n]-expectedH) > 8*tolerance || std::abs(values[count+n]-expectedD) > 8*tolerance) throw std::runtime_error("The curved HEX27 regression failed its independent reference.");
      maximumVolumeError = std::max(maximumVolumeError,std::abs(values[n]-expectedH));
      if(metric == 2) maximumSurfaceError = std::max(maximumSurfaceError,std::abs(values[count+n]-expectedD));
      checks += 2;
    }
    ++cases;
  }
  std::cout << std::setprecision(8) << cases << " curved HEX27 cases and " << checks << " checks passed; maximum reference volume error " << maximumVolumeError << ", physical surface error " << maximumSurfaceError << '\n';
}
