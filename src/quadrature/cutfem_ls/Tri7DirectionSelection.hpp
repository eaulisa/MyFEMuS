#ifndef FEMUS_TRI7_DIRECTION_SELECTION_HPP
#define FEMUS_TRI7_DIRECTION_SELECTION_HPP

#include "TriangleIntegration.hpp"
#include "EndpointSimplexChart.hpp"

namespace FELevelSetMoments {

inline void Tri7SurfaceCoordinatePolynomials(const std::vector<std::vector<double>>& geometryCoordinates, std::vector<std::vector<std::vector<long double>>>& coordinatePolynomials) {
  ValidateSurfaceGeometryCoordinates(geometryCoordinates,2,7);
  coordinatePolynomials.assign(2,std::vector<std::vector<long double>>(3,std::vector<long double>(3,0.0L)));
  for(unsigned component = 0; component < 2; ++component) {
    std::vector<double> shifted(7);
    for(unsigned node = 0; node < 7; ++node) shifted[node] = geometryCoordinates[component][node]-geometryCoordinates[component][0];
    std::vector<std::vector<double>> coefficient;
    Tri7ReferenceCoefficientsFromNodalValues(shifted,coefficient);
    for(unsigned i = 0; i < 3; ++i) {
      for(unsigned j = 0; j < 3; ++j) coordinatePolynomials[component][i][j] = coefficient[i][j];
    }
  }
}

inline void IntegrateTri7ReferenceSelected(const std::vector<std::vector<double>>& coefficient, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool computeDelta = true, bool fullBoundaryDelta = false, unsigned maximumDepth = 28, unsigned integrationDirection = 0, unsigned maximumEvaluations = 200000, std::vector<double>* estimatedError = nullptr, unsigned validationMode = 1, const std::vector<unsigned>* zeroBoundary = nullptr,  bool referenceSurface = false, const std::vector<std::vector<double>>* inverseJacobian = nullptr, const std::vector<std::vector<double>>* geometryCoordinates = nullptr, const OuterQuadratureRules& rules = OuterQuadratureRules(), bool inverseSurface = false) {
  rules.Validate();
  if(inverseSurface && geometryCoordinates == nullptr) throw std::invalid_argument("Inverse surface moments require nodal geometry.");
  if(geometryCoordinates != nullptr) {
    if(!referenceSurface || inverseJacobian != nullptr) throw std::invalid_argument("Surface geometry requires surface moments without an inverse Jacobian.");
    ValidateSurfaceGeometryCoordinates(*geometryCoordinates,2,7);
  }
  if(inverseJacobian != nullptr) {
    if(!referenceSurface) throw std::invalid_argument("A surface inverse Jacobian requires surface moments.");
    ValidateSurfaceInverseJacobian(*inverseJacobian,2);
  }
  if(&positive == &delta || coefficient.size() != 3 || degree > 64 || !(absoluteTolerance > 0.0) || !(relativeTolerance >= 0.0) || !std::isfinite(absoluteTolerance) || !std::isfinite(relativeTolerance) || maximumDepth > 60 || integrationDirection > 3 || validationMode > 2 || maximumEvaluations < 3) throw std::invalid_argument("Invalid TRI7 selected-direction input.");
  for(unsigned i = 0; i < 3; ++i) if(coefficient[i].size() != 3) throw std::invalid_argument("The TRI7 coefficient matrix must be 3 by 3.");
  if(coefficient[2][2] != 0.0 || coefficient[2][1] != coefficient[1][2]) throw std::invalid_argument("The TRI7 coefficients must represent a quadratic polynomial with a barycentric bubble.");
  if(zeroBoundary != nullptr && zeroBoundary->size() != 3) throw std::invalid_argument("TRI7 requires three boundary flags.");

  std::vector<unsigned> canonicalBoundary;
  if(zeroBoundary != nullptr) { canonicalBoundary = {(*zeroBoundary)[1],(*zeroBoundary)[0],(*zeroBoundary)[2]}; zeroBoundary = &canonicalBoundary; }
  std::vector<std::vector<long double>> field;
  ReferencePolynomialAffine2D(coefficient,0,1,0,1,field);
  std::vector<std::vector<long double>> normalized;
  long double scale;
  LineSliceUtilities::NormalizeCanonicalPolynomial(field,normalized,scale);
  std::vector<std::vector<std::vector<long double>>> gradient;
  if(referenceSurface) {
    SurfaceGradientFromPolynomial(coefficient,1,1,gradient);
    if(inverseJacobian != nullptr) SurfaceGradientPhysical(gradient,*inverseJacobian);
    if(geometryCoordinates != nullptr) {
      std::vector<std::vector<std::vector<long double>>> coordinates;
      Tri7SurfaceCoordinatePolynomials(*geometryCoordinates,coordinates);
      SurfaceGradientCofactor2D(gradient,coordinates,1,1);
    }
  }
  std::vector<long double> values;
  std::vector<long double> errors;
  unsigned evaluations = 0;
  unsigned accepted = 0;
  unsigned selected = 0;
  const long double amplification = integrationDirection == 3 ? ReferencePower(3.0L,degree) : 1;
  std::vector<unsigned> transformedBoundary;
  if(integrationDirection == 3) {
    const std::vector<std::vector<long double>> transformed = EndpointSimplexChart::Triangle(normalized);
    for(unsigned i = 0; i < transformed.size(); ++i) for(unsigned j = 0; j < transformed[i].size(); ++j) {
      if((i > 2 || j > 2) && transformed[i][j] != 0) throw std::logic_error("A TRI chart requires a higher-degree exact level-set representation.");
      if(i < 3 && j < 3) normalized[i][j] = transformed[i][j];
    }
    for(unsigned component = 0; component < gradient.size(); ++component) gradient[component] = EndpointSimplexChart::Triangle(gradient[component]);
    if(zeroBoundary != nullptr) { transformedBoundary = {(*zeroBoundary)[2],(*zeroBoundary)[1],(*zeroBoundary)[0]}; zeroBoundary = &transformedBoundary; }
  }
  TriangleIntegration::Integrate(normalized,degree,computeDelta,absoluteTolerance/amplification,absoluteTolerance*scale/amplification,relativeTolerance/amplification,maximumDepth,maximumEvaluations,values,errors,evaluations,accepted,selected,referenceSurface ? &gradient : nullptr,nullptr,nullptr,integrationDirection == 3 ? 1U : integrationDirection,fullBoundaryDelta,zeroBoundary,rules,inverseSurface);
  unsigned compared = 0;
  if(validationMode == 2) compared = PlanarEndpointIntegration::Compare(normalized,degree,true,computeDelta,absoluteTolerance/amplification,absoluteTolerance*scale/amplification,relativeTolerance/amplification,maximumDepth,maximumEvaluations,values,errors,evaluations,accepted,selected,referenceSurface ? &gradient : nullptr,nullptr,fullBoundaryDelta,zeroBoundary,rules,inverseSurface);
  if(integrationDirection == 3) {
    std::vector<long double> source;
    std::vector<long double> sourceErrors;
    source.swap(values);
    sourceErrors.swap(errors);
    EndpointSimplexChart::RestoreTriangle(degree,source,sourceErrors,values,errors);
    selected = 3;
  }
  const unsigned count = values.size()/2U;
  std::vector<std::vector<double>> resultH(degree+1U);
  std::vector<std::vector<double>> resultD(degree+1U);
  long double errorH = 0;
  long double errorD = 0;
  for(unsigned i = 0; i <= degree; ++i) {
    resultH[i].resize(degree-i+1U);
    resultD[i].resize(degree-i+1U);
    for(unsigned j = 0; j <= degree-i; ++j) {
      const unsigned n = ReferenceMomentIndex2D(degree,i,j);
      resultH[i][j] = static_cast<double>(values[n]);
      resultD[i][j] = computeDelta ? static_cast<double>(values[count+n]/scale) : 0;
      const long double hError = errors[n]+std::abs(values[n]-resultH[i][j]);
      const long double dError = computeDelta ? errors[count+n]/scale+std::abs(values[count+n]/scale-resultD[i][j]) : 0;
      if(!std::isfinite(resultH[i][j]) || !std::isfinite(resultD[i][j]) || hError > absoluteTolerance+relativeTolerance*std::abs(values[n]) || dError > absoluteTolerance+relativeTolerance*std::abs(resultD[i][j])) throw std::runtime_error("Rounded TRI moments exceed tolerance.");
      errorH = std::max(errorH,hError);
      errorD = std::max(errorD,dError);
    }
  }
  positive.swap(resultH);
  delta.swap(resultD);
  if(estimatedError != nullptr) *estimatedError = {static_cast<double>(errorH),static_cast<double>(errorD),static_cast<double>(accepted),static_cast<double>(evaluations),static_cast<double>(selected),static_cast<double>(compared),1};
}

inline void IntegrateTri7NodalSelected(const std::vector<double>& nodalValues, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool computeDelta = true, bool fullBoundaryDelta = false, unsigned validationMode = 1, unsigned integrationDirection = 0, std::vector<double>* estimatedError = nullptr,  bool referenceSurface = false, const std::vector<std::vector<double>>* inverseJacobian = nullptr, const std::vector<std::vector<double>>* geometryCoordinates = nullptr, const OuterQuadratureRules& rules = OuterQuadratureRules(), bool inverseSurface = false) {
  std::vector<std::vector<double>> coefficient;
  Tri7ReferenceCoefficientsFromNodalValues(nodalValues,coefficient);
  const std::vector<unsigned> boundary = {static_cast<unsigned>(nodalValues[0] == 0.0 && nodalValues[1] == 0.0 && nodalValues[3] == 0.0),static_cast<unsigned>(nodalValues[0] == 0.0 && nodalValues[2] == 0.0 && nodalValues[5] == 0.0),static_cast<unsigned>(nodalValues[1] == 0.0 && nodalValues[2] == 0.0 && nodalValues[4] == 0.0)};
  IntegrateTri7ReferenceSelected(coefficient,degree,positive,delta,absoluteTolerance,relativeTolerance,computeDelta,fullBoundaryDelta,28U,integrationDirection,200000,estimatedError,validationMode,&boundary,referenceSurface,inverseJacobian,geometryCoordinates,rules,inverseSurface);
}

}

#endif
