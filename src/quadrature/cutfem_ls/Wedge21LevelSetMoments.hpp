#ifndef FEMUS_WEDGE21_LEVEL_SET_MOMENTS_HPP
#define FEMUS_WEDGE21_LEVEL_SET_MOMENTS_HPP

#include "Hex27DirectionSelection.hpp"

namespace FELevelSetMoments {

// Reconstructs FEMuS TRI7 x LINE3 on the unit prism with t = (z+1)/2.
inline void Wedge21UnitPolynomialFromNodalValues(const std::vector<double>& nodal, std::vector<long double>& coefficient) {
  if(nodal.size() != 21) throw std::invalid_argument("WEDGE21 requires twenty-one nodal values in FEMuS order.");
  const unsigned nodes[3][7] = {{0,1,2,6,7,8,18},{12,13,14,15,16,17,20},{3,4,5,9,10,11,19}};
  const long double interpolation[3][3] = {{1,0,0},{-3,4,-1},{2,-4,2}};
  coefficient.assign(27,0.0L);
  for(unsigned plane = 0; plane < 3; ++plane) {
    std::vector<double> values(7);
    for(unsigned node = 0; node < 7; ++node) values[node] = nodal[nodes[plane][node]];
    std::vector<std::vector<double>> triangle;
    Tri7ReferenceCoefficientsFromNodalValues(values,triangle);
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) for(unsigned k = 0; k < 3; ++k) coefficient[Hex27Index(i,j,k)] += interpolation[k][plane]*triangle[i][j];
  }
}

inline void IntegrateWedge21NodalMoments(const std::vector<double>& nodal, unsigned degree, std::vector<long double>& values, std::vector<long double>& errors, long double absolute = 1.0e-12L, long double relative = 1.0e-12L, bool referenceSurface = false, const std::vector<std::vector<double>>* geometry = nullptr, unsigned direction = 0, unsigned validation = 1, const OuterQuadratureRules& rules = OuterQuadratureRules(), bool inverseSurface = false) {
  rules.Validate();
  if(inverseSurface && geometry == nullptr) throw std::invalid_argument("Inverse surface moments require nodal geometry.");
  if(degree > 64 || !(absolute > 0) || !(relative >= 0) || !std::isfinite(absolute) || !std::isfinite(relative) || direction > 6 || validation > 2 || (geometry != nullptr && !referenceSurface)) throw std::invalid_argument("Invalid WEDGE21 moment request.");
  if(geometry != nullptr) ValidateSurfaceGeometryCoordinates(*geometry,3,21);
  std::vector<long double> original;
  Wedge21UnitPolynomialFromNodalValues(nodal,original);
  long double scale = 0;
  for(unsigned n = 0; n < original.size(); ++n) scale = std::max(scale,std::abs(original[n]));
  if(scale == 0) throw std::domain_error("The delta of an identically zero WEDGE21 field is undefined.");
  int exponent = 0;
  std::frexp(scale,&exponent);
  scale = std::ldexp(1.0L,exponent-1);
  SurfaceEndpointData data;
  data.wedge = true;
  data.quadratureRules = rules;
  data.inverseSurface = inverseSurface;
  data.coefficient = original;
  for(unsigned n = 0; n < original.size(); ++n) data.coefficient[n] /= scale;
  data.degree = degree;
  data.absoluteH = absolute;
  data.absoluteD = absolute*scale;
  data.relative = relative;
  const unsigned faceNodes[5][9] = {{0,2,5,3,8,14,11,12,17},{0,1,4,3,6,13,9,12,15},{1,2,5,4,7,14,10,13,16},{0,1,2,6,7,8,18,18,18},{3,4,5,9,10,11,19,19,19}};
  data.zeroFaces.assign(5,1);
  for(unsigned face = 0; face < 5; ++face) for(unsigned node = 0; node < 9; ++node) data.zeroFaces[face] = data.zeroFaces[face] && nodal[faceNodes[face][node]] == 0;
  if(direction != 0) { unsigned middle; Hex27SelectedAxes(direction,data.preferredInner,middle,data.preferredOuter); }
  if(referenceSurface) {
    data.gradient.resize(3);
    for(unsigned axis = 0; axis < 3; ++axis) data.gradient[axis] = SurfaceEndpointIntegration::Derivative(original,axis);
    if(geometry == nullptr) for(unsigned n = 0; n < data.gradient[2].size(); ++n) data.gradient[2][n] *= 0.5L;
    else {
      std::vector<std::vector<long double>> coordinates(3);
      for(unsigned component = 0; component < 3; ++component) {
        std::vector<double> centered(21);
        for(unsigned node = 0; node < 21; ++node) centered[node] = (*geometry)[component][node]-(*geometry)[component][0];
        Wedge21UnitPolynomialFromNodalValues(centered,coordinates[component]);
      }
      SurfaceGradientCofactor3D(data.gradient,coordinates,1.0L);
    }
  }
  SurfaceEndpointIntegration::Integrate(data,false,values,errors);
  if(validation == 2) SurfaceEndpointIntegration::Compare(data,false,values,errors);
  const unsigned count = values.size()/2U;
  for(unsigned n = count; n < values.size(); ++n) { values[n] /= scale; errors[n] /= scale; }
}

}

#endif
