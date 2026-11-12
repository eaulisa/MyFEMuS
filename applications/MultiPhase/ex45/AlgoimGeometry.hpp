#ifndef FEMUS_EX45_ALGOIM_GEOMETRY_HPP
#define FEMUS_EX45_ALGOIM_GEOMETRY_HPP

#include "ElemType.hpp"
#include <cmath>
#include <stdexcept>
#include <vector>

inline void EvaluateAlgoimGeometry(const femus::elem_type* finiteElement, const std::vector<std::vector<double>>& coordinates, const std::vector<double>& point, double& determinant, std::vector<std::vector<double>>& jacobian, std::vector<std::vector<double>>& inverseJacobian) {
  finiteElement->GetJacobianMatrix(coordinates,point,determinant,jacobian,inverseJacobian);
  if(!std::isfinite(determinant) || determinant == 0.) throw std::domain_error("Invalid element geometry at an integration point.");
}

inline void EvaluateAlgoimShape(const femus::elem_type* finiteElement, const std::vector<double>& point, const std::vector<std::vector<double>>& inverseJacobian, std::vector<double>& shape, std::vector<double>& physicalGradient) {
  const unsigned dimension = finiteElement->GetDim();
  const unsigned nodes = static_cast<unsigned>(finiteElement->GetNDofs());
  finiteElement->GetPhi(shape,point);
  physicalGradient.assign(nodes*dimension,0.);
  std::vector<double> derivative;
  for(unsigned reference = 0; reference < dimension; ++reference) {
    if(reference == 0) finiteElement->GetDPhiDXi(derivative,point);
    else if(reference == 1) finiteElement->GetDPhiDEta(derivative,point);
    else finiteElement->GetDPhiDZeta(derivative,point);
    for(unsigned node = 0; node < nodes; ++node) {
      for(unsigned physical = 0; physical < dimension; ++physical) physicalGradient[node*dimension+physical] += derivative[node]*inverseJacobian[reference][physical];
    }
  }
}

inline double MapAlgoimSurfaceWeight(double referenceWeight, double determinant, const std::vector<std::vector<double>>& inverseJacobian, const std::vector<double>& referenceNormal, std::vector<double>& physicalNormal) {
  const unsigned dimension = static_cast<unsigned>(referenceNormal.size());
  physicalNormal.assign(dimension,0.);
  double normalScale = 0.;
  for(unsigned physical = 0; physical < dimension; ++physical) {
    for(unsigned reference = 0; reference < dimension; ++reference) physicalNormal[physical] += inverseJacobian[reference][physical]*referenceNormal[reference];
    normalScale = std::hypot(normalScale,physicalNormal[physical]);
  }
  if(!(normalScale > 0.) || !std::isfinite(normalScale)) throw std::domain_error("Invalid interface normal at an integration point.");
  for(unsigned physical = 0; physical < dimension; ++physical) physicalNormal[physical] /= normalScale;
  return referenceWeight*std::abs(determinant)*normalScale;
}

#endif
