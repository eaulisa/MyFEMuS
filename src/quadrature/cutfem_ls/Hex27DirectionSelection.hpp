#ifndef FEMUS_HEX27_DIRECTION_SELECTION_HPP
#define FEMUS_HEX27_DIRECTION_SELECTION_HPP

#include "SurfaceEndpointIntegration.hpp"

namespace FELevelSetMoments {

inline void Hex27SelectedAxes(unsigned order, unsigned& inner, unsigned& middle, unsigned& outer) {
  if(order == 0 || order > 6) throw std::invalid_argument("Invalid explicit HEX27 integration order.");
  const unsigned axes[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
  inner = axes[order-1U][0];
  middle = axes[order-1U][1];
  outer = axes[order-1U][2];
}

inline void IntegrateHex27ReferenceSelected(const std::vector<std::vector<std::vector<double>>>& coefficient, unsigned degree, std::vector<std::vector<std::vector<double>>>& positive, std::vector<std::vector<std::vector<double>>>& delta, double absoluteTolerance = 1.0e-10, double relativeTolerance = 1.0e-10, bool computeDelta = true, unsigned maximumDepth = 24, unsigned integrationOrder = 0, unsigned maximumEvaluations = 200000, std::vector<double>* estimatedError = nullptr, unsigned validationMode = 1,  const std::vector<unsigned>* zeroFaces = nullptr, std::vector<double>* directionInfo = nullptr, bool referenceSurface = false, const std::vector<std::vector<double>>* inverseJacobian = nullptr, const std::vector<std::vector<double>>* geometryCoordinates = nullptr) {
  if(&positive == &delta || degree > 64 || !(absoluteTolerance > 0.0) || !(relativeTolerance >= 0.0) || !std::isfinite(absoluteTolerance) || !std::isfinite(relativeTolerance) || maximumDepth > 48 || integrationOrder > 6 || maximumEvaluations == 0 || validationMode > 2 || (directionInfo != nullptr && directionInfo == estimatedError)) throw std::invalid_argument("Invalid selected-direction HEX27 integration input.");
  if(inverseJacobian != nullptr && !referenceSurface) throw std::invalid_argument("A physical HEX27 metric requires surface moments.");
  if(inverseJacobian != nullptr) ValidateSurfaceInverseJacobian(*inverseJacobian,3);
  if(geometryCoordinates != nullptr && (!referenceSurface || inverseJacobian != nullptr)) throw std::invalid_argument("HEX27 surface geometry requires surface moments without an inverse Jacobian.");
  if(geometryCoordinates != nullptr) ValidateSurfaceGeometryCoordinates(*geometryCoordinates,3,27);
  if(zeroFaces != nullptr) {
    if(zeroFaces->size() != 6) throw std::invalid_argument("HEX27 requires six zero-face flags.");
    for(unsigned face = 0; face < 6; ++face) {
      if((*zeroFaces)[face] > 1) throw std::invalid_argument("Invalid HEX27 zero-face flag.");
    }
  }
  std::vector<long double> reference;
  long double scale;
  Hex27Normalize(coefficient,reference,scale);

  SurfaceEndpointData data;
  data.coefficient = reference;
  if(zeroFaces != nullptr) data.zeroFaces = *zeroFaces;
  data.degree = degree;
  data.maximumDepth = maximumDepth;
  data.maximumEvaluations = maximumEvaluations;
  data.computeDelta = computeDelta;
  data.absoluteH = absoluteTolerance;
  data.absoluteD = static_cast<long double>(absoluteTolerance)*scale;
  data.relative = relativeTolerance;
  if(integrationOrder != 0) { unsigned middle; Hex27SelectedAxes(integrationOrder,data.preferredInner,middle,data.preferredOuter); }
  if(referenceSurface && computeDelta) Hex27ReferenceSurfaceGradient(reference,scale,data.gradient,inverseJacobian,geometryCoordinates);
  std::vector<long double> values;
  std::vector<long double> errors;
  SurfaceEndpointIntegration::Integrate(data,false,values,errors);
  if(validationMode == 2) SurfaceEndpointIntegration::Compare(data,false,values,errors);
  SurfaceEndpointIntegration::Output(data,scale,values,errors,positive,delta,estimatedError);
  if(directionInfo != nullptr) *directionInfo = {static_cast<double>(data.graphAxis+1U),0,1,1,static_cast<double>(data.innerCounts[0]),static_cast<double>(data.innerCounts[1]),static_cast<double>(data.innerCounts[2])};
}

inline void IntegrateHex27NodalSelected(const std::vector<double>& nodalValues, unsigned degree, std::vector<std::vector<std::vector<double>>>& positive, std::vector<std::vector<std::vector<double>>>& delta, double absoluteTolerance = 1.0e-10, double relativeTolerance = 1.0e-10, bool computeDelta = true, unsigned validationMode = 1, unsigned integrationOrder = 0, std::vector<double>* estimatedError = nullptr,  std::vector<double>* directionInfo = nullptr, bool referenceSurface = false, const std::vector<std::vector<double>>* inverseJacobian = nullptr, const std::vector<std::vector<double>>* geometryCoordinates = nullptr) {
  std::vector<std::vector<std::vector<double>>> coefficient;
  Hex27ReferenceCoefficientsFromNodalValues(nodalValues,coefficient);
  std::vector<unsigned> zeroFaces(6,0);
  for(unsigned axis = 0; axis < 3; ++axis) {
    for(unsigned side = 0; side < 2; ++side) {
      std::vector<double> face;
      Hex27FaceNodalValues(nodalValues,axis,side,face);
      bool zero = true;
      for(unsigned node = 0; node < 9; ++node) zero = zero && face[node] == 0.0;
      zeroFaces[2U*axis+side] = zero;
    }
  }
  IntegrateHex27ReferenceSelected(coefficient,degree,positive,delta,absoluteTolerance,relativeTolerance,computeDelta,24,integrationOrder,200000,estimatedError,validationMode,&zeroFaces,directionInfo,referenceSurface,inverseJacobian,geometryCoordinates);
}

}

#endif
