#ifndef FEMUS_UNIT_ELEMENT_MOMENTS_HPP
#define FEMUS_UNIT_ELEMENT_MOMENTS_HPP

#include "Hex27DirectionSelection.hpp"
#include "TangencyStripIntegration.hpp"

namespace FELevelSetMoments {

inline unsigned UnitMomentIndex(unsigned dimension, unsigned i, unsigned j = 0, unsigned k = 0) {
  if(dimension == 1) return i;
  if(dimension == 2) return (i+j)*(i+j+1U)/2U+j;
  return HexMomentIndex(i,j,k);
}

// Reconstructs the tensor Lagrange polynomial directly at unit-cell nodes 0, 1/2, and 1.
inline void UnitTensorPolynomialFromNodalValues(const std::vector<double>& nodal, unsigned dimension, std::vector<long double>& coefficient) {
  if(dimension < 1 || dimension > 3 || nodal.size() != (dimension == 1 ? 3U : (dimension == 2 ? 9U : 27U))) throw std::invalid_argument("Invalid unit-cell nodal field.");
  const unsigned line[3][3] = {{0,0,0},{2,0,0},{1,0,0}};
  const unsigned quad[9][3] = {{0,0,0},{2,0,0},{2,2,0},{0,2,0},{1,0,0},{2,1,0},{1,2,0},{0,1,0},{1,1,0}};
  const unsigned hex[27][3] = {{0,0,0},{2,0,0},{2,2,0},{0,2,0},{0,0,2},{2,0,2},{2,2,2},{0,2,2},{1,0,0},{2,1,0},{1,2,0},{0,1,0},{1,0,2},{2,1,2},{1,2,2},{0,1,2},{0,0,1},{2,0,1},{2,2,1},{0,2,1},{1,0,1},{2,1,1},{1,2,1},{0,1,1},{1,1,0},{1,1,2},{1,1,1}};
  const long double interpolation[3][3] = {{1,0,0},{-3,4,-1},{2,-4,2}};
  coefficient.assign(27,0.0L);
  for(unsigned n = 0; n < nodal.size(); ++n) {
    if(!std::isfinite(nodal[n])) throw std::invalid_argument("A unit-cell nodal value is not finite.");
    const unsigned* node = dimension == 1 ? line[n] : (dimension == 2 ? quad[n] : hex[n]);
    coefficient[Hex27Index(node[0],node[1],node[2])] = nodal[n];
  }
  for(unsigned axis = 0; axis < dimension; ++axis) {
    std::vector<long double> result(27,0.0L);
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < (dimension > 1 ? 3U : 1U); ++j) for(unsigned k = 0; k < (dimension > 2 ? 3U : 1U); ++k) {
      unsigned index[3] = {i,j,k};
      const unsigned power = index[axis];
      for(unsigned sample = 0; sample < 3; ++sample) {
        index[axis] = sample;
        result[Hex27Index(i,j,k)] += interpolation[power][sample]*coefficient[Hex27Index(index[0],index[1],index[2])];
      }
    }
    coefficient.swap(result);
  }
}

inline void UnitTensorGradient(const std::vector<long double>& coefficient, unsigned dimension, const std::vector<std::vector<double>>* geometry, std::vector<std::vector<long double>>& gradient) {
  gradient.resize(dimension);
  for(unsigned axis = 0; axis < dimension; ++axis) gradient[axis] = GraphCellIntegration::Derivative(coefficient,axis);
  if(geometry == nullptr || dimension == 1) return;
  ValidateSurfaceGeometryCoordinates(*geometry,dimension,dimension == 2 ? 9U : 27U);
  std::vector<std::vector<long double>> coordinates(dimension);
  for(unsigned component = 0; component < dimension; ++component) UnitTensorPolynomialFromNodalValues((*geometry)[component],dimension,coordinates[component]);
  if(dimension == 3) SurfaceGradientCofactor3D(gradient,coordinates,1.0L);
  else {
    std::vector<std::vector<std::vector<long double>>> planarGradient(2,std::vector<std::vector<long double>>(3,std::vector<long double>(3)));
    std::vector<std::vector<std::vector<long double>>> planarCoordinates = planarGradient;
    for(unsigned component = 0; component < 2; ++component) for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) {
      planarGradient[component][i][j] = gradient[component][Hex27Index(i,j,0)];
      planarCoordinates[component][i][j] = coordinates[component][Hex27Index(i,j,0)];
    }
    SurfaceGradientCofactor2D(planarGradient,planarCoordinates,1,1);
    const unsigned side = planarGradient[0].size();
    gradient.assign(2,std::vector<long double>(side*side*side,0.0L));
    for(unsigned component = 0; component < 2; ++component) for(unsigned i = 0; i < side; ++i) for(unsigned j = 0; j < side; ++j) gradient[component][side*side*i+side*j] = planarGradient[component][i][j];
  }
}

inline void UnitTensorCheck(unsigned count, const std::vector<long double>& values, const std::vector<long double>& errors, long double absoluteH, long double absoluteD, long double relative) {
  for(unsigned n = 0; n < 2U*count; ++n) {
    const long double tolerance = (n < count ? absoluteH : absoluteD)+relative*std::abs(values[n]);
    if(!std::isfinite(values[n]) || !std::isfinite(errors[n]) || errors[n] < 0.0L || errors[n] > tolerance) throw std::runtime_error("A unit-coordinate moment exceeds its requested tolerance.");
  }
}

inline void UnitQuadMoments(const std::vector<long double>& coefficient, long double scale, unsigned degree, const std::vector<std::vector<long double>>& metric, long double absolute, long double relative, unsigned direction, unsigned validation, std::vector<long double>& values, std::vector<long double>& errors, const std::vector<std::vector<double>>* geometry = nullptr) {
  std::vector<std::vector<long double>> field(3,std::vector<long double>(3));
  for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) field[i][j] = coefficient[Hex27Index(i,j,0)];
  std::vector<std::vector<std::vector<long double>>> gradient;
  for(unsigned component = 0; component < metric.size(); ++component) {
    const unsigned side = SurfacePolynomialCubeSize(metric[component].size());
    gradient.push_back(std::vector<std::vector<long double>>(side,std::vector<long double>(side)));
    for(unsigned i = 0; i < side; ++i) for(unsigned j = 0; j < side; ++j) gradient.back()[i][j] = metric[component][side*side*i+side*j];
  }
  long double lower[2];
  long double upper[2];
  const unsigned automatic = Quad9SelectedOrder(field,lower,upper);
  const unsigned selected = direction == 0 ? automatic : direction;
  const std::array<long double,5> coordinates = {0,1,0,1,1};
  unsigned evaluations = 0;
  unsigned accepted = 0;
  std::vector<long double> native;
  std::vector<long double> nativeErrors;
  bool tangencyIntegrated = false;
  if(direction == 0 && validation != 2 && lower[0] == 0 && lower[1] == 0 && (upper[0] > 0 || upper[1] > 0) && Quad9SelectedWeakGraphOrder(field) == 0) {
    int phase = 0;
    LineSliceUtilities::CanonicalGraphOrder(field,0,phase);
    unsigned chosen = 0;
    if(phase == 0) {
      std::vector<std::vector<std::vector<long double>>> coordinates;
      if(geometry != nullptr) for(unsigned d = 0; d < 2; ++d) {
        std::vector<long double> polynomial;
        UnitTensorPolynomialFromNodalValues((*geometry)[d],2,polynomial);
        coordinates.push_back(std::vector<std::vector<long double>>(3,std::vector<long double>(3)));
        for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) coordinates[d][i][j] = polynomial[Hex27Index(i,j,0)];
      }
      tangencyIntegrated = TangencyStripIntegration::Try(field,degree,false,true,absolute,absolute*scale,relative,24,200000,native,nativeErrors,evaluations,accepted,chosen,gradient.empty() ? nullptr : &gradient,scale,coordinates.empty() ? nullptr : &coordinates);
    }
  }
  for(unsigned attempt = 0; !tangencyIntegrated && attempt < (validation == 2 ? 2U : 1U); ++attempt) {
    const unsigned order = attempt == 0 ? selected : 3U-selected;
    std::vector<long double> trial;
    std::vector<long double> trialErrors;
    try { Quad9SelectedDirect(field,degree,order,true,absolute,absolute*scale,relative,24,200000,trial,trialErrors,evaluations,accepted,gradient.empty() ? nullptr : &gradient,coordinates,nullptr); }
    catch(const std::runtime_error&) {
      std::vector<std::vector<long double>> transposed = field;
      for(unsigned i = 0; i < 3; ++i) for(unsigned j = i+1U; j < 3; ++j) std::swap(transposed[i][j],transposed[j][i]);
      const std::vector<unsigned> boundary = {static_cast<unsigned>(LineIntegration::BoundaryIsZero(transposed,0,0)),static_cast<unsigned>(LineIntegration::BoundaryIsZero(transposed,1,0)),static_cast<unsigned>(LineIntegration::BoundaryIsZero(field,0,0)),static_cast<unsigned>(LineIntegration::BoundaryIsZero(field,1,0))};
      Quad9SelectedOneDirection(field,scale,degree,order,lower[2U-order],true,false,absolute,relative,24,200000,&boundary,trial,trialErrors,evaluations,accepted,true,gradient.empty() ? nullptr : &gradient);
    }
    UnitTensorCheck(trial.size()/2U,trial,trialErrors,absolute,absolute*scale,relative);
    if(attempt == 0) { native.swap(trial); nativeErrors.swap(trialErrors); }
    else for(unsigned n = 0; n < native.size(); ++n) nativeErrors[n] = std::max(nativeErrors[n],std::max(trialErrors[n],std::abs(native[n]-trial[n])));
  }
  const unsigned count = (degree+1U)*(degree+2U)/2U;
  values.assign(2U*count,0);
  errors.assign(2U*count,0);
  for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) {
    const unsigned source = ReferenceMomentIndex2D(degree,i,j);
    const unsigned target = UnitMomentIndex(2,i,j);
    values[target] = native[source];
    errors[target] = nativeErrors[source];
    values[count+target] = native[count+source]/scale;
    errors[count+target] = nativeErrors[count+source]/scale;
  }
}

inline void IntegrateUnitTensorNodalMoments(const std::vector<double>& nodal, unsigned dimension, unsigned degree, std::vector<long double>& values, std::vector<long double>& errors, long double absolute = 1.0e-12L, long double relative = 1.0e-12L, bool surface = false, const std::vector<std::vector<double>>* geometry = nullptr, unsigned direction = 0, unsigned validation = 1) {
  if(degree > 64 || !(absolute > 0) || !(relative >= 0) || !std::isfinite(absolute) || !std::isfinite(relative) || validation > 2 || direction > (dimension == 1 ? 1U : (dimension == 2 ? 2U : 6U)) || (geometry != nullptr && !surface)) throw std::invalid_argument("Invalid unit-coordinate moment request.");
  std::vector<long double> original;
  UnitTensorPolynomialFromNodalValues(nodal,dimension,original);
  long double scale = 0;
  for(unsigned n = 0; n < original.size(); ++n) scale = std::max(scale,std::abs(original[n]));
  if(scale == 0) throw std::domain_error("The delta of an identically zero unit-cell field is undefined.");
  int exponent = 0;
  std::frexp(scale,&exponent);
  scale = std::ldexp(1.0L,exponent-1);
  std::vector<long double> coefficient = original;
  for(unsigned n = 0; n < coefficient.size(); ++n) coefficient[n] /= scale;
  std::vector<std::vector<long double>> metric;
  if(surface) UnitTensorGradient(original,dimension,geometry,metric);
  unsigned count = degree+1U;
  if(dimension == 2) count = (degree+1U)*(degree+2U)/2U;
  if(dimension == 3) count = (degree+1U)*(degree+2U)*(degree+3U)/6U;
  if(dimension == 1) {
    std::vector<long double> positive;
    std::vector<long double> delta;
    const long double a = coefficient[Hex27Index(2,0,0)];
    const long double b = coefficient[Hex27Index(1,0,0)];
    const long double c = coefficient[0];
    LineIntegration::IntegratePolynomial(a,b,c,degree,0,1,false,false,positive,delta);
    std::vector<long double> roots = QuadraticRoots(a,b,c);
    if(nodal[0] == 0) { roots = {0}; if(a != 0) roots.push_back(-b/a); }
    else if(nodal[1] == 0) { roots = {1}; if(a != 0) roots.push_back(c/a); }
    values.assign(2U*count,0);
    errors.assign(2U*count,0);
    for(unsigned n = 0; n < count; ++n) values[n] = positive[n];
    for(unsigned n = 0; n < roots.size(); ++n) {
      const long double root = roots[n];
      if(root < 0 || root > 1) continue;
      const long double derivative = std::abs(std::fma(2*a,root,b));
      if(derivative == 0) throw std::domain_error("A repeated line root has no regular unit-coordinate delta moment.");
      const long double weight = (root == 0 || root == 1 ? 0.5L : 1.0L)*(surface ? 1.0L : 1.0L/(scale*derivative));
      long double power = 1;
      for(unsigned i = 0; i <= degree; ++i) { values[count+i] += weight*power; power *= root; }
    }
    for(unsigned n = 0; n < values.size(); ++n) errors[n] = 128*std::numeric_limits<long double>::epsilon()*std::abs(values[n]);
  }
  else if(dimension == 2) UnitQuadMoments(coefficient,scale,degree,metric,absolute,relative,direction,validation,values,errors,geometry);
  else {
    GraphCellData data;
    data.coefficient = coefficient;
    data.gradient = metric;
    data.degree = degree;
    data.unitCube = true;
    data.absoluteH = absolute;
    data.absoluteD = absolute*scale;
    data.relative = relative;
    if(direction != 0) { unsigned middle; Hex27SelectedAxes(direction,data.preferredInner,middle,data.preferredOuter); }
    GraphCellIntegration::Integrate(data,false,values,errors);
    if(validation == 2) {
      GraphCellData independent = data;
      independent.evaluations = independent.accepted = 0;
      unsigned outer = 0;
      for(unsigned axis = 1; axis < 3; ++axis) if(data.outerCounts[axis] > data.outerCounts[outer]) outer = axis;
      independent.preferredOuter = 3U-data.graphAxis-outer;
      std::vector<long double> second;
      std::vector<long double> secondErrors;
      GraphCellIntegration::Integrate(independent,false,second,secondErrors);
      for(unsigned n = 0; n < values.size(); ++n) errors[n] = std::max(errors[n],std::max(secondErrors[n],std::abs(values[n]-second[n])));
    }
    for(unsigned n = count; n < values.size(); ++n) { values[n] /= scale; errors[n] /= scale; }
  }
  UnitTensorCheck(count,values,errors,absolute,absolute,relative);
}

}

#endif
