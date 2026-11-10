#ifndef FEMUS_QUAD9_DIRECTION_SELECTION_HPP
#define FEMUS_QUAD9_DIRECTION_SELECTION_HPP

#include "PlanarEndpointIntegration.hpp"

namespace FELevelSetMoments {

inline void Quad9SelectedQuadraticRange(long double a, long double b, long double c, long double& minimum, long double& maximum) {
  const long double endpoint = (a+b)+c;
  minimum = std::min(minimum,std::min(c,endpoint));
  maximum = std::max(maximum,std::max(c,endpoint));
  if(a != 0.0L && ((b < 0.0L && a > 0.0L) || (b > 0.0L && a < 0.0L)) && std::abs(b) < 2.0L*std::abs(a)) {
    const long double point = -b/(2.0L*a);
    const long double value = std::fma(std::fma(a,point,b),point,c);
    minimum = std::min(minimum,value);
    maximum = std::max(maximum,value);
  }
}

inline void Quad9SelectedDerivativeBounds(const std::vector<std::vector<long double>>& coefficient, unsigned axis, long double& minimum, long double& maximum) {
  if(coefficient.size() != 3 || axis > 1) throw std::invalid_argument("Invalid QUAD9 derivative input.");
  long double norm = 0.0L;
  for(unsigned i = 0; i < 3; ++i) {
    if(coefficient[i].size() != 3) throw std::invalid_argument("QUAD9 coefficients must be 3 by 3.");
    for(unsigned j = 0; j < 3; ++j) {
      if(!std::isfinite(coefficient[i][j])) throw std::invalid_argument("A QUAD9 coefficient is not finite.");
      norm += (axis == 0 ? i : j)*std::abs(coefficient[i][j]);
    }
  }
  minimum = std::numeric_limits<long double>::infinity();
  maximum = -minimum;
  for(unsigned side = 0; side < 2; ++side) {
    long double q[3];
    for(unsigned j = 0; j < 3; ++j) q[j] = axis == 0 ? coefficient[1][j]+2.0L*side*coefficient[2][j] : coefficient[j][1]+2.0L*side*coefficient[j][2];
    Quad9SelectedQuadraticRange(q[2],q[1],q[0],minimum,maximum);
  }
  const long double padding = 512.0L*std::numeric_limits<long double>::epsilon()*norm;
  minimum -= padding;
  maximum += padding;
}

inline unsigned Quad9SelectedOrder(const std::vector<std::vector<long double>>& coefficient, long double lower[2], long double upper[2]) {
  for(unsigned axis = 0; axis < 2; ++axis) {
    long double minimum;
    long double maximum;
    Quad9SelectedDerivativeBounds(coefficient,axis,minimum,maximum);
    upper[axis] = std::max(std::abs(minimum),std::abs(maximum));
    lower[axis] = minimum > 0.0L ? minimum : (maximum < 0.0L ? -maximum : 0.0L);
    const long double conditioning = 32.0L*std::sqrt(std::numeric_limits<long double>::epsilon())*upper[axis];
    if(lower[axis] <= conditioning) lower[axis] = 0.0L;
  }
  if(lower[0] != lower[1]) return lower[0] > lower[1] ? 2U : 1U;
  return upper[0] > upper[1] ? 2U : 1U;
}

inline unsigned Quad9SelectedWeakGraphOrder(const std::vector<std::vector<long double>>& coefficient, unsigned requestedOrder = 0) {
  long double bernstein[3][3];
  CanonicalBernstein(coefficient,bernstein);
  long double norm = 0.0L;
  for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) norm += std::abs(coefficient[i][j]);
  const long double padding = 1024.0L*std::numeric_limits<long double>::epsilon()*std::max(1.0L,norm);
  unsigned selected = 0;
  long double selectedMagnitude = 0.0L;
  for(unsigned order = 1; order <= 2; ++order) {
    if(requestedOrder != 0 && requestedOrder != order) continue;
    long double minimum = std::numeric_limits<long double>::infinity();
    long double maximum = -minimum;
    for(unsigned i = 0; i < 3; ++i) {
      const long double low = order == 1 ? bernstein[i][0] : bernstein[0][i];
      const long double middle = order == 1 ? bernstein[i][1] : bernstein[1][i];
      const long double high = order == 1 ? bernstein[i][2] : bernstein[2][i];
      const long double derivative[3] = {2.0L*(middle-low),high-low,2.0L*(high-middle)};
      for(unsigned j = 0; j < 3; ++j) {
        minimum = std::min(minimum,derivative[j]);
        maximum = std::max(maximum,derivative[j]);
      }
    }
    const bool nonnegative = minimum >= -padding && maximum > padding;
    const bool nonpositive = maximum <= padding && minimum < -padding;
    if(nonnegative || nonpositive) {
      const long double magnitude = std::max(std::abs(minimum),std::abs(maximum));
      if(magnitude > selectedMagnitude) {
        selectedMagnitude = magnitude;
        selected = order;
      }
    }
  }
  return selected;
}

inline int Quad9SelectedReflectionParity(const std::vector<std::vector<long double>>& coefficient, unsigned axis) {
  unsigned columns = 0;
  for(unsigned row = 0; row < coefficient.size(); ++row) columns = std::max(columns,static_cast<unsigned>(coefficient[row].size()));
  const unsigned longitudinal = axis == 0 ? coefficient.size() : columns;
  const unsigned transverse = axis == 0 ? columns : coefficient.size();
  bool even = true;
  bool odd = true;
  for(unsigned k = 0; k < transverse; ++k) {
    long double polynomial[3] = {0.0L,0.0L,0.0L};
    for(unsigned power = 0; power < longitudinal; ++power) {
      const unsigned row = axis == 0 ? power : k;
      const unsigned column = axis == 0 ? k : power;
      const long double value = column < coefficient[row].size() ? coefficient[row][column] : 0.0L;
      if(power > 2 && value != 0.0L) return 0;
      if(power <= 2) polynomial[power] = value;
    }
    even = even && polynomial[1] == -polynomial[2];
    odd = odd && polynomial[2] == 0.0L && polynomial[1] == -2.0L*polynomial[0];
  }
  return even ? 1 : (odd ? -1 : 0);
}

inline void Quad9SelectedReferenceParity(const std::vector<std::vector<long double>>& coefficient, unsigned degree, std::vector<long double>& values, std::vector<long double>& errors, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient) {
  const unsigned count = (degree+1U)*(degree+2U)/2U;
  for(unsigned axis = 0; axis < 2; ++axis) {
    const int parity = Quad9SelectedReflectionParity(coefficient,axis);
    if(parity == 0) continue;
    bool deltaSymmetric = true;
    if(surfaceGradient != nullptr) {
      for(unsigned component = 0; component < surfaceGradient->size(); ++component) deltaSymmetric = deltaSymmetric && Quad9SelectedReflectionParity((*surfaceGradient)[component],axis) != 0;
    }
    for(unsigned i = 0; i <= degree; ++i) {
      for(unsigned j = 0; j <= degree-i; ++j) {
        if((axis == 0 ? i : j)%2U == 0) continue;
        const unsigned index = LineSliceUtilities::MomentIndex(degree,i,j);
        if(parity == 1) {
          values[index] = 0.0L;
          errors[index] = 0.0L;
        }
        if(deltaSymmetric) {
          values[count+index] = 0.0L;
          errors[count+index] = 0.0L;
        }
      }
    }
  }
}

inline void Quad9SurfaceCoordinatePolynomials(const std::vector<std::vector<double>>& geometryCoordinates, std::vector<std::vector<std::vector<long double>>>& coordinatePolynomials) {
  ValidateSurfaceGeometryCoordinates(geometryCoordinates,2,9);
  coordinatePolynomials.assign(2,std::vector<std::vector<long double>>(3,std::vector<long double>(3,0.0L)));
  for(unsigned component = 0; component < 2; ++component) {
    std::vector<double> shifted(9);
    for(unsigned node = 0; node < 9; ++node) shifted[node] = geometryCoordinates[component][node]-geometryCoordinates[component][0];
    std::vector<std::vector<double>> coefficient;
    Quad9ReferenceCoefficientsFromNodalValues(shifted,coefficient);
    for(unsigned i = 0; i < 3; ++i) {
      for(unsigned j = 0; j < 3; ++j) coordinatePolynomials[component][i][j] = coefficient[i][j];
    }
  }
}

inline void Quad9SelectedDirect(const std::vector<std::vector<long double>>& coefficient, unsigned degree, unsigned direction, bool computeDelta, long double absoluteH, long double absoluteD, long double relative, unsigned maximumDepth, unsigned maximumEvaluations, std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted, const std::vector<std::vector<std::vector<long double>>>* gradient, const std::array<long double,5>& coordinates, const std::vector<unsigned char>* requested, const std::vector<unsigned>* zeroBoundary = nullptr) {
  unsigned selected = 0;
  PlanarEndpointIntegration::Integrate(coefficient,degree,false,computeDelta,absoluteH,absoluteD,relative,maximumDepth,maximumEvaluations,values,errors,evaluations,accepted,selected,gradient,&coordinates,requested,direction,false,zeroBoundary);
}

inline void IntegrateQuad9ReferenceSelected(const std::vector<std::vector<double>>& coefficient, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool computeDelta = true, bool fullBoundaryDelta = false, unsigned maximumDepth = 28, unsigned integrationOrder = 0, unsigned maximumEvaluations = 200000, std::vector<double>* estimatedError = nullptr, unsigned validationMode = 1, const std::vector<unsigned>* zeroBoundary = nullptr,  bool referenceSurface = false, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient = nullptr, const std::vector<std::vector<double>>* inverseJacobian = nullptr, const std::vector<std::vector<double>>* geometryCoordinates = nullptr) {
  if(geometryCoordinates != nullptr) {
    if(!referenceSurface || inverseJacobian != nullptr) throw std::invalid_argument("Surface geometry requires surface moments without an inverse Jacobian.");
    ValidateSurfaceGeometryCoordinates(*geometryCoordinates,2,9);
  }
  if(inverseJacobian != nullptr) {
    if(!referenceSurface) throw std::invalid_argument("A surface inverse Jacobian requires surface moments.");
    ValidateSurfaceInverseJacobian(*inverseJacobian,2);
  }
  if(&positive == &delta || coefficient.size() != 3 || degree > 64 || !(absoluteTolerance > 0.0) || !(relativeTolerance >= 0.0) || !std::isfinite(absoluteTolerance) || !std::isfinite(relativeTolerance) || maximumDepth > 60 || integrationOrder > 2 || validationMode > 2 || maximumEvaluations < 15) throw std::invalid_argument("Invalid QUAD9 selected-direction input.");
  for(unsigned i = 0; i < 3; ++i) if(coefficient[i].size() != 3) throw std::invalid_argument("The QUAD9 coefficient matrix must be 3 by 3.");
  if(zeroBoundary != nullptr && zeroBoundary->size() != 4) throw std::invalid_argument("QUAD9 requires four boundary flags.");
  std::vector<std::vector<long double>> chartCoefficient;
  ReferencePolynomialAffine2D(coefficient,-1.0L,2.0L,-1.0L,2.0L,chartCoefficient);
  std::vector<std::vector<long double>> normalized;
  long double scale;
  LineSliceUtilities::NormalizeCanonicalPolynomial(chartCoefficient,normalized,scale);
  std::vector<std::vector<std::vector<long double>>> referenceGradient;
  if(referenceSurface && surfaceGradient == nullptr) {
    SurfaceGradientFromPolynomial(coefficient,1.0L,1.0L,referenceGradient);
    surfaceGradient = &referenceGradient;
  }
  if(inverseJacobian != nullptr) {
    if(surfaceGradient != &referenceGradient) referenceGradient = *surfaceGradient;
    SurfaceGradientPhysical(referenceGradient,*inverseJacobian);
    surfaceGradient = &referenceGradient;
  }
  if(geometryCoordinates != nullptr) {
    if(surfaceGradient != &referenceGradient) referenceGradient = *surfaceGradient;
    std::vector<std::vector<std::vector<long double>>> coordinatePolynomials;
    Quad9SurfaceCoordinatePolynomials(*geometryCoordinates,coordinatePolynomials);
    SurfaceGradientCofactor2D(referenceGradient,coordinatePolynomials,1.0L,1.0L);
    surfaceGradient = &referenceGradient;
  }
  std::vector<std::vector<std::vector<long double>>> chartGradient;
  if(surfaceGradient != nullptr) {
    chartGradient.resize(surfaceGradient->size());
    for(unsigned component = 0; component < surfaceGradient->size(); ++component) ReferencePolynomialAffine2D((*surfaceGradient)[component],-1.0L,2.0L,-1.0L,2.0L,chartGradient[component]);
    surfaceGradient = &chartGradient;
  }
  const std::array<long double,5> parent = {-1,2,-1,2,4};
  std::vector<long double> values;
  std::vector<long double> errors;
  unsigned evaluations = 0;
  unsigned accepted = 0;
  unsigned selected = 0;
  PlanarEndpointIntegration::Integrate(normalized,degree,false,computeDelta,absoluteTolerance,absoluteTolerance*scale,relativeTolerance,maximumDepth,maximumEvaluations,values,errors,evaluations,accepted,selected,surfaceGradient,&parent,nullptr,integrationOrder,fullBoundaryDelta,zeroBoundary);
  unsigned compared = 0;
  if(validationMode == 2) compared = PlanarEndpointIntegration::Compare(normalized,degree,false,computeDelta,absoluteTolerance,absoluteTolerance*scale,relativeTolerance,maximumDepth,maximumEvaluations,values,errors,evaluations,accepted,selected,surfaceGradient,&parent,fullBoundaryDelta,zeroBoundary);
  Quad9SelectedReferenceParity(normalized,degree,values,errors,surfaceGradient);
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
      if(!std::isfinite(resultH[i][j]) || !std::isfinite(resultD[i][j]) || hError > absoluteTolerance+relativeTolerance*std::abs(values[n]) || dError > absoluteTolerance+relativeTolerance*std::abs(resultD[i][j])) throw std::runtime_error("Rounded QUAD moments exceed tolerance.");
      errorH = std::max(errorH,hError);
      errorD = std::max(errorD,dError);
    }
  }
  positive.swap(resultH);
  delta.swap(resultD);
  if(estimatedError != nullptr) *estimatedError = {static_cast<double>(errorH),static_cast<double>(errorD),static_cast<double>(accepted),static_cast<double>(evaluations),static_cast<double>(selected),static_cast<double>(compared),1};
}

inline void IntegrateQuad9NodalSelected(const std::vector<double>& nodalValues, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool computeDelta = true, bool fullBoundaryDelta = false, unsigned validationMode = 1, unsigned integrationOrder = 0, std::vector<double>* estimatedError = nullptr,  bool referenceSurface = false, const std::vector<std::vector<double>>* inverseJacobian = nullptr, const std::vector<std::vector<double>>* geometryCoordinates = nullptr) {
  std::vector<std::vector<double>> coefficient;
  Quad9ReferenceCoefficientsFromNodalValues(nodalValues,coefficient);
  std::vector<unsigned> boundary = {static_cast<unsigned>(nodalValues[0] == 0.0 && nodalValues[3] == 0.0 && nodalValues[7] == 0.0),static_cast<unsigned>(nodalValues[1] == 0.0 && nodalValues[2] == 0.0 && nodalValues[5] == 0.0),static_cast<unsigned>(nodalValues[0] == 0.0 && nodalValues[1] == 0.0 && nodalValues[4] == 0.0),static_cast<unsigned>(nodalValues[3] == 0.0 && nodalValues[2] == 0.0 && nodalValues[6] == 0.0)};
  IntegrateQuad9ReferenceSelected(coefficient,degree,positive,delta,absoluteTolerance,relativeTolerance,computeDelta,fullBoundaryDelta,28U,integrationOrder,200000,estimatedError,validationMode,&boundary,referenceSurface,nullptr,inverseJacobian,geometryCoordinates);
}

}

#endif
