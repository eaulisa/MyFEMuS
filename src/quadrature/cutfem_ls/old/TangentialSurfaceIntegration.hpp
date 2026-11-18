#ifndef FEMUS_TANGENTIAL_SURFACE_INTEGRATION_HPP
#define FEMUS_TANGENTIAL_SURFACE_INTEGRATION_HPP

#include <stdexcept>
#include <vector>

namespace FELevelSetMoments {

// Returns |cof(J) grad_ref G|^2 and polynomial force and stiffness factors in native reference coordinates.
inline double TangentialSurfaceFactors(const std::vector<std::vector<double>>& jacobian, const std::vector<double>& gradient, const std::vector<double>& shapeGradient, std::vector<double>& force, std::vector<double>& tangent) {
  const unsigned dimension = gradient.size();
  if((dimension != 2 && dimension != 3) || jacobian.size() != dimension || shapeGradient.size()%dimension != 0) throw std::invalid_argument("Invalid tangential surface factor dimensions.");
  for(unsigned d = 0; d < dimension; ++d) if(jacobian[d].size() != dimension) throw std::invalid_argument("Invalid tangential surface Jacobian dimensions.");
  const unsigned count = shapeGradient.size()/dimension;
  const unsigned tangentSize = dimension == 2 ? 1U : 3U;
  std::vector<double> normal(dimension,0.);
  force.assign(count*dimension,0.);
  tangent.assign(count*tangentSize,0.);
  if(dimension == 2) {
    normal[0] = jacobian[1][1]*gradient[0]-jacobian[1][0]*gradient[1];
    normal[1] = jacobian[0][0]*gradient[1]-jacobian[0][1]*gradient[0];
    const double first = -jacobian[0][0]*gradient[1]+jacobian[0][1]*gradient[0];
    const double second = -jacobian[1][0]*gradient[1]+jacobian[1][1]*gradient[0];
    for(unsigned i = 0; i < count; ++i) {
      tangent[i] = -gradient[1]*shapeGradient[2U*i]+gradient[0]*shapeGradient[2U*i+1U];
      force[2U*i] = first*tangent[i];
      force[2U*i+1U] = second*tangent[i];
    }
  }
  else {
    for(unsigned physical = 0; physical < 3; ++physical) for(unsigned reference = 0; reference < 3; ++reference) normal[physical] += (jacobian[(physical+1U)%3U][(reference+1U)%3U]*jacobian[(physical+2U)%3U][(reference+2U)%3U]-jacobian[(physical+1U)%3U][(reference+2U)%3U]*jacobian[(physical+2U)%3U][(reference+1U)%3U])*gradient[reference];
    for(unsigned i = 0; i < count; ++i) {
      double cross[3];
      for(unsigned d = 0; d < 3; ++d) cross[d] = shapeGradient[3U*i+(d+1U)%3U]*gradient[(d+2U)%3U]-shapeGradient[3U*i+(d+2U)%3U]*gradient[(d+1U)%3U];
      for(unsigned d = 0; d < 3; ++d) for(unsigned k = 0; k < 3; ++k) tangent[3U*i+d] += jacobian[d][k]*cross[k];
      for(unsigned d = 0; d < 3; ++d) force[3U*i+d] = normal[(d+1U)%3U]*tangent[3U*i+(d+2U)%3U]-normal[(d+2U)%3U]*tangent[3U*i+(d+1U)%3U];
    }
  }
  double squared = 0.;
  for(unsigned d = 0; d < dimension; ++d) squared += normal[d]*normal[d];
  return squared;
}

}

#endif
