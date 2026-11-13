#ifndef FEMUS_TET15_LEVEL_SET_MOMENTS_HPP
#define FEMUS_TET15_LEVEL_SET_MOMENTS_HPP

#include "Tet15LevelSet.hpp"
#include "EndpointSimplexChart.hpp"
#include "SurfaceEndpointIntegration.hpp"

namespace FELevelSetMoments {

inline unsigned Tet15MomentIndex(unsigned i, unsigned j, unsigned k) {
  unsigned total = i+j+k;
  return total*(total+1U)*(total+2U)/6U+(total-i)*(total-i+1U)/2U+k;
}

inline unsigned Tet15TriangleIndex(unsigned degree, unsigned i, unsigned j) {
  return i*(2U*degree+3U-i)/2U+j;
}

inline void Tet15ReferenceSurfaceGradient(const std::vector<double>& nodalValues, std::vector<std::vector<long double>>& gradient, const std::vector<std::vector<double>>* inverseJacobian = nullptr, const std::vector<std::vector<double>>* geometryCoordinates = nullptr) {
  std::vector<long double> coefficient;
  Tet15PolynomialFromNodalValues(nodalValues,false,coefficient);
  gradient.assign(3,std::vector<long double>(125,0.0L));
  for(unsigned i = 0; i < 5; ++i) {
    for(unsigned j = 0; i+j < 5; ++j) {
      for(unsigned k = 0; i+j+k < 5; ++k) {
        const long double value = coefficient[Tet15Index(i,j,k)];
        if(i != 0) gradient[0][Tet15Index(i-1U,j,k)] = i*value;
        if(j != 0) gradient[1][Tet15Index(i,j-1U,k)] = j*value;
        if(k != 0) gradient[2][Tet15Index(i,j,k-1U)] = k*value;
      }
    }
  }
  if(inverseJacobian != nullptr) SurfaceGradientPhysical(gradient,*inverseJacobian);
  if(geometryCoordinates != nullptr) {
    std::vector<std::vector<long double>> coordinates(3);
    for(unsigned component = 0; component < 3; ++component) {
      std::vector<double> nodal(15);
      for(unsigned node = 0; node < 15; ++node) nodal[node] = (*geometryCoordinates)[component][node]-(*geometryCoordinates)[component][0];
      Tet15PolynomialFromNodalValues(nodal,false,coordinates[component]);
    }
    SurfaceGradientCofactor3D(gradient,coordinates,1.0L);
  }
}

inline void IntegrateTet15NodalByTriangleSlices(const std::vector<double>& nodalValues, unsigned degree, std::vector<std::vector<std::vector<double>>>& positive, std::vector<std::vector<std::vector<double>>>& delta, double absoluteTolerance = 1.0e-10, double relativeTolerance = 1.0e-10, bool computeDelta = true, unsigned maximumDepth = 24, unsigned integrationAxis = 0, unsigned maximumEvaluations = 200000, std::vector<double>* estimatedError = nullptr, unsigned directionValidation = 1, std::vector<double>* directionInfo = nullptr, bool referenceSurface = false, const std::vector<std::vector<double>>* inverseJacobian = nullptr, const std::vector<std::vector<double>>* geometryCoordinates = nullptr, const OuterQuadratureRules& rules = OuterQuadratureRules(), bool inverseSurface = false) {
  QuadratureWarningScope warningScope("TET15 moments");
  rules.Validate();
  if(inverseSurface && geometryCoordinates == nullptr) throw std::invalid_argument("Inverse surface moments require nodal geometry.");
  if(nodalValues.size() != 15 || &positive == &delta || degree > 64 || !(absoluteTolerance > 0.0) || !(relativeTolerance >= 0.0) || !std::isfinite(absoluteTolerance) || !std::isfinite(relativeTolerance) || maximumDepth > 48 || integrationAxis > 4 || maximumEvaluations < 3 || directionValidation > 2 || (directionInfo != nullptr && directionInfo == estimatedError)) throw std::invalid_argument("Invalid adaptive TET15 integration input.");
  if(inverseJacobian != nullptr && !referenceSurface) throw std::invalid_argument("A physical TET15 metric requires surface moments.");
  if(inverseJacobian != nullptr) ValidateSurfaceInverseJacobian(*inverseJacobian,3);
  if(geometryCoordinates != nullptr && (!referenceSurface || inverseJacobian != nullptr)) throw std::invalid_argument("TET15 surface geometry requires surface moments without an inverse Jacobian.");
  if(geometryCoordinates != nullptr) ValidateSurfaceGeometryCoordinates(*geometryCoordinates,3,15);
  long double scale = 0.0L;
  for(unsigned node = 0; node < 15; ++node) {
    if(!std::isfinite(nodalValues[node])) throw std::invalid_argument("A TET15 nodal value is not finite.");
    scale = std::max(scale,std::abs(static_cast<long double>(nodalValues[node])));
  }
  if(scale == 0.0L && computeDelta) throw std::domain_error("The Dirac delta of an identically zero TET15 field is undefined.");
  int exponent = 0;
  std::frexp(scale,&exponent);
  scale = scale == 0.0L ? 1.0L : std::ldexp(1.0L,exponent-1);

  const unsigned permutation[15] = {1,0,2,3,4,6,5,8,7,9,10,11,13,12,14};
  std::vector<double> fieldNodal = nodalValues;
  std::vector<std::vector<double>> chartGeometry;
  if(integrationAxis == 1) {
    for(unsigned n = 0; n < 15; ++n) fieldNodal[n] = nodalValues[permutation[n]];
    if(geometryCoordinates != nullptr) {
      chartGeometry = *geometryCoordinates;
      for(unsigned component = 0; component < 3; ++component) for(unsigned n = 0; n < 15; ++n) chartGeometry[component][n] = (*geometryCoordinates)[component][permutation[n]];
      geometryCoordinates = &chartGeometry;
    }
  }
  SurfaceEndpointData data;
  data.quadratureRules = rules;
  data.inverseSurface = inverseSurface;
  const unsigned boundaryNodes[4][7] = {{1,2,3,5,9,8,12},{0,2,3,6,9,7,13},{0,1,3,4,8,7,11},{0,1,2,4,5,6,10}};
  data.zeroFaces.assign(4,1);
  for(unsigned face = 0; face < 4; ++face) for(unsigned node = 0; node < 7; ++node) data.zeroFaces[face] = data.zeroFaces[face] && fieldNodal[boundaryNodes[face][node]] == 0;
  std::vector<long double> original;
  Tet15PolynomialFromNodalValues(fieldNodal,false,original);
  data.coefficient.assign(27,0);
  for(unsigned i = 0; i < 5; ++i) for(unsigned j = 0; j < 5; ++j) for(unsigned k = 0; k < 5; ++k) {
    if((i > 2 || j > 2 || k > 2) && original[Tet15Index(i,j,k)] != 0) throw std::logic_error("The exact enriched TET polynomial cannot be represented by this slice interface.");
    if(i < 3 && j < 3 && k < 3) data.coefficient[Hex27Index(i,j,k)] = original[Tet15Index(i,j,k)]/scale;
  }
  data.degree = degree;
  data.maximumDepth = maximumDepth;
  data.maximumEvaluations = maximumEvaluations;
  data.computeDelta = computeDelta;
  const long double amplification = integrationAxis == 1 ? ReferencePower(4.0L,degree) : 1;
  data.absoluteH = absoluteTolerance/amplification;
  data.absoluteD = static_cast<long double>(absoluteTolerance)*scale/amplification;
  data.relative = relativeTolerance/amplification;
  if(integrationAxis != 0) { data.preferredOuter = integrationAxis == 1 ? 0U : integrationAxis-2U; }
  if(referenceSurface && computeDelta) {
    Tet15ReferenceSurfaceGradient(fieldNodal,data.gradient,integrationAxis == 1 && geometryCoordinates == nullptr ? nullptr : inverseJacobian,geometryCoordinates);
    if(integrationAxis == 1 && geometryCoordinates == nullptr) {
      for(unsigned n = 0; n < data.gradient[0].size(); ++n) { const long double first = data.gradient[0][n]; data.gradient[0][n] = -first; data.gradient[1][n] -= first; data.gradient[2][n] -= first; }
      if(inverseJacobian != nullptr) SurfaceGradientPhysical(data.gradient,*inverseJacobian);
    }
  }
  std::vector<long double> values;
  std::vector<long double> errors;
  SurfaceEndpointIntegration::Integrate(data,true,values,errors);
  if(directionValidation == 2) SurfaceEndpointIntegration::Compare(data,true,values,errors);
  if(integrationAxis == 1) {
    std::vector<long double> source;
    std::vector<long double> sourceErrors;
    source.swap(values);
    sourceErrors.swap(errors);
    EndpointSimplexChart::RestoreTetrahedron(degree,source,sourceErrors,values,errors);
    data.absoluteH = absoluteTolerance;
    data.absoluteD = absoluteTolerance*scale;
    data.relative = relativeTolerance;
  }
  SurfaceEndpointIntegration::Output(data,scale,values,errors,positive,delta,estimatedError);
  if(directionInfo != nullptr) *directionInfo = {1,1,0};
}

}

#endif
