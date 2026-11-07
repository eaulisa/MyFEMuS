#include "FemusAlgoim.hpp"
#include "FELevelSetMoments.hpp"
#include "Hex27LevelSetMoments.hpp"
#include "Hex27DirectionSelection.hpp"
#include "Tet15LevelSet.hpp"
#include "algoim/quadrature_multipoly.hpp"
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

namespace FemusAlgoim {
namespace {

long double Binomial(unsigned degree, unsigned index) {
  long double result = 1.0L;
  for(unsigned i = 1; i <= index; ++i) {
    result *= static_cast<long double>(degree+1U-i)/i;
  }
  return result;
}

template<int N>
void AddMoments(const algoim::uvector<double,N>& point, long double weight, unsigned simplexDimension, unsigned degree, std::vector<long double>& powers, std::vector<long double>& moments) {
  for(int axis = 0; axis < N; ++axis) {
    long double coordinate = point(axis);
    if(simplexDimension == 2 && axis == 1) {
      coordinate *= point(0);
    }
    if(simplexDimension == 3 && axis >= 1) {
      coordinate *= point(0);
      if(axis == 2) coordinate *= point(1);
    }
    unsigned offset = static_cast<unsigned>(axis)*(degree+1U);
    powers[offset] = 1.0L;
    for(unsigned i = 1; i <= degree; ++i) {
      powers[offset+i] = powers[offset+i-1U]*coordinate;
    }
  }
  unsigned index = 0;
  for(unsigned total = 0; total <= degree; ++total) {
    if constexpr(N == 2) {
      for(unsigned j = 0; j <= total; ++j) {
        moments[index] += weight*powers[total-j]*powers[degree+1U+j];
        ++index;
      }
    }
    else {
      for(int i = static_cast<int>(total); i >= 0; --i) {
        for(int j = static_cast<int>(total)-i; j >= 0; --j) {
          unsigned k = total-static_cast<unsigned>(i+j);
          moments[index] += weight*powers[i]*powers[degree+1U+j]*powers[2U*(degree+1U)+k];
          ++index;
        }
      }
    }
  }
}

template<int N>
void VolumePoint(const algoim::uvector<double,N>& local, double weight, const algoim::xarray<double,N>& polynomial, const algoim::uvector<double,N>& low, double width, unsigned simplexDimension, unsigned degree, std::vector<long double>& powers, std::vector<long double>& moments) {
  algoim::uvector<double,N> point = low+width*local;
  double value = algoim::bernstein::evalBernsteinPoly(polynomial,point);
  if(!std::isfinite(value) || !std::isfinite(weight)) {
    throw std::runtime_error("Algoim returned a non-finite volume sample.");
  }
  if(value <= 0.0) {
    return;
  }
  long double mappedWeight = weight;
  for(int axis = 0; axis < N; ++axis) {
    mappedWeight *= width;
  }
  if(simplexDimension == 2) {
    mappedWeight *= point(0);
  }
  else if(simplexDimension == 3) {
    mappedWeight *= static_cast<long double>(point(0))*point(0)*point(1);
  }
  AddMoments<N>(point,mappedWeight,simplexDimension,degree,powers,moments);
}

void Line3SurfaceMoments(double a, double b, double c, unsigned degree, std::vector<double>& delta, bool leftEndpoint, bool rightEndpoint) {
  const long double scale = std::max(std::abs(static_cast<long double>(a)),std::max(std::abs(static_cast<long double>(b)),std::abs(static_cast<long double>(c))));
  if(scale == 0.0L) throw std::domain_error("The Dirac delta of an identically zero level set is undefined.");
  const long double aa = a/scale;
  const long double bb = b/scale;
  const long double cc = c/scale;
  std::vector<long double> roots;
  if(leftEndpoint && rightEndpoint) roots = {0.0L,1.0L};
  else if(leftEndpoint) {
    roots.push_back(0.0L);
    if(aa != 0.0L) roots.push_back(-bb/aa);
  }
  else if(rightEndpoint) {
    roots.push_back(1.0L);
    if(aa != 0.0L) roots.push_back(cc/aa);
  }
  else roots = FELevelSetMoments::QuadraticRoots(aa,bb,cc);
  std::vector<long double> moments(degree+1U,0.0L);
  for(unsigned index = 0; index < roots.size(); ++index) {
    const long double root = roots[index];
    if(root < 0.0L || root > 1.0L) continue;
    const long double derivative = rightEndpoint ? aa-cc : leftEndpoint ? bb : 2.0L*aa*root+bb;
    if(derivative == 0.0L) throw std::domain_error("A repeated 1D zero has no regular surface moment.");
    long double value = root == 0.0L || root == 1.0L ? 0.25L : 0.5L;
    for(unsigned power = 0; power <= degree; ++power) {
      moments[power] += value;
      value *= root;
    }
  }
  delta.resize(degree+1U);
  for(unsigned power = 0; power <= degree; ++power) delta[power] = static_cast<double>(moments[power]);
}

void ValidateGradientTransform(const std::vector<std::vector<double>>* inverseJacobian, unsigned dimension, bool referenceSurface) {
  if(inverseJacobian == nullptr) return;
  if(!referenceSurface) throw std::invalid_argument("A surface inverse Jacobian requires gradient-weighted moments.");
  FELevelSetMoments::ValidateSurfaceInverseJacobian(*inverseJacobian,dimension);
}

void ValidateSurfaceGeometry(const std::vector<std::vector<double>>* geometryCoordinates, const std::vector<std::vector<double>>* inverseJacobian, unsigned dimension, unsigned nodes, bool referenceSurface) {
  if(geometryCoordinates == nullptr) return;
  if(!referenceSurface || inverseJacobian != nullptr) throw std::invalid_argument("Surface geometry requires gradient-weighted moments without an inverse Jacobian.");
  FELevelSetMoments::ValidateSurfaceGeometryCoordinates(*geometryCoordinates,dimension,nodes);
}

std::vector<std::vector<long double>> SurfacePolynomial2D(const std::vector<long double>& coefficient, unsigned nodes) {
  std::vector<std::vector<long double>> polynomial(3,std::vector<long double>(3,0.0L));
  if(nodes == 7) {
    polynomial[0][0] = coefficient[0];
    polynomial[1][0] = coefficient[1];
    polynomial[0][1] = coefficient[2];
    polynomial[2][0] = coefficient[3];
    polynomial[1][1] = coefficient[4];
    polynomial[0][2] = coefficient[5];
    polynomial[2][1] = coefficient[6];
    polynomial[1][2] = -coefficient[6];
  }
  else {
    for(unsigned i = 0; i < 3; ++i) {
      for(unsigned j = 0; j < 3; ++j) polynomial[i][j] = coefficient[3U*i+j];
    }
  }
  return polynomial;
}

void SurfaceCoordinatePolynomial(const std::vector<double>& coordinate, unsigned nodes, std::vector<long double>& coefficient) {
  std::vector<double> translated(nodes);
  for(unsigned node = 0; node < nodes; ++node) translated[node] = coordinate[node]-coordinate[0];
  if(nodes == 7) {
    std::vector<double> data;
    FELevelSetMoments::Tri7ReflectedCoefficientsFromNodalValues(translated,data);
    coefficient.assign(data.begin(),data.end());
  }
  else if(nodes == 9) {
    std::vector<std::vector<double>> data;
    FELevelSetMoments::Quad9UnitCoefficientsFromNodalValues(translated,data);
    coefficient.resize(9);
    for(unsigned i = 0; i < 3; ++i) {
      for(unsigned j = 0; j < 3; ++j) coefficient[3U*i+j] = data[i][j];
    }
  }
  else if(nodes == 27) {
    std::vector<std::vector<std::vector<double>>> data;
    FELevelSetMoments::Hex27UnitCoefficientsFromNodalValues(translated,data);
    coefficient.resize(27);
    for(unsigned i = 0; i < 3; ++i) {
      for(unsigned j = 0; j < 3; ++j) {
        for(unsigned k = 0; k < 3; ++k) coefficient[9U*i+3U*j+k] = data[i][j][k];
      }
    }
  }
  else FELevelSetMoments::Tet15PolynomialFromNodalValues(translated,false,coefficient);
}

void SurfaceCofactorGradient(const std::vector<long double>& coefficient, unsigned nodes, const std::vector<std::vector<double>>* geometryCoordinates, std::vector<std::vector<long double>>& gradient, unsigned& extent) {
  extent = 0;
  if(geometryCoordinates == nullptr) return;
  if(nodes == 7 || nodes == 9) {
    const long double firstFactor = nodes == 7 ? -1.0L : 0.5L;
    const long double secondFactor = nodes == 7 ? 1.0L : 0.5L;
    const std::vector<std::vector<long double>> polynomial = SurfacePolynomial2D(coefficient,nodes);
    std::vector<std::vector<std::vector<long double>>> field(2);
    std::vector<std::vector<std::vector<long double>>> coordinates(2);
    field[0] = FELevelSetMoments::SurfacePolynomialDerivative2D(polynomial,0,firstFactor);
    field[1] = FELevelSetMoments::SurfacePolynomialDerivative2D(polynomial,1,secondFactor);
    for(unsigned physical = 0; physical < 2; ++physical) {
      std::vector<long double> data;
      SurfaceCoordinatePolynomial((*geometryCoordinates)[physical],nodes,data);
      coordinates[physical] = SurfacePolynomial2D(data,nodes);
    }
    FELevelSetMoments::SurfaceGradientCofactor2D(field,coordinates,firstFactor,secondFactor);
    extent = static_cast<unsigned>(field[0].size());
    gradient.assign(2,std::vector<long double>(extent*extent,0.0L));
    for(unsigned physical = 0; physical < 2; ++physical) {
      for(unsigned i = 0; i < extent; ++i) {
        for(unsigned j = 0; j < extent; ++j) gradient[physical][extent*i+j] = field[physical][i][j];
      }
    }
  }
  else {
    const long double factor = nodes == 15 ? 1.0L : 0.5L;
    std::vector<std::vector<long double>> coordinates(3);
    gradient.resize(3);
    for(unsigned physical = 0; physical < 3; ++physical) {
      gradient[physical] = FELevelSetMoments::SurfacePolynomialDerivative3D(coefficient,physical,factor);
      SurfaceCoordinatePolynomial((*geometryCoordinates)[physical],nodes,coordinates[physical]);
    }
    FELevelSetMoments::SurfaceGradientCofactor3D(gradient,coordinates,factor);
    extent = FELevelSetMoments::SurfacePolynomialCubeSize(static_cast<unsigned>(gradient[0].size()));
  }
}

long double SurfaceCofactorNorm(const std::vector<std::vector<long double>>& gradient, unsigned extent, const long double* point, unsigned dimension) {
  long double norm = 0.0L;
  for(unsigned component = 0; component < dimension; ++component) {
    long double value = 0.0L;
    for(unsigned i = extent; i > 0; --i) {
      long double row = 0.0L;
      for(unsigned j = extent; j > 0; --j) {
        long double column = gradient[component][extent*(i-1U)+j-1U];
        if(dimension == 3) {
          column = 0.0L;
          for(unsigned k = extent; k > 0; --k) column = std::fma(column,point[2],gradient[component][extent*extent*(i-1U)+extent*(j-1U)+k-1U]);
        }
        row = std::fma(row,point[1],column);
      }
      value = std::fma(value,point[0],row);
    }
    norm = std::hypot(norm,value);
  }
  return norm;
}

template<int N>
long double SurfaceCofactorNorm(const std::vector<std::vector<long double>>& gradient, unsigned extent, const algoim::uvector<double,N>& point, unsigned simplexDimension) {
  long double referencePoint[3] = {};
  for(int axis = 0; axis < N; ++axis) referencePoint[axis] = point(axis);
  if(simplexDimension == 2) referencePoint[1] *= referencePoint[0];
  if(simplexDimension == 3) {
    referencePoint[0] = static_cast<long double>(point(0))*(1.0L-point(1));
    referencePoint[1] = static_cast<long double>(point(0))*point(1)*(1.0L-point(2));
    referencePoint[2] = static_cast<long double>(point(0))*point(1)*point(2);
  }
  return SurfaceCofactorNorm(gradient,extent,referencePoint,N);
}

long double SurfaceGradientNorm(const long double* gradient, unsigned dimension, const std::vector<std::vector<double>>* inverseJacobian) {
  long double norm = 0.0L;
  for(unsigned physical = 0; physical < dimension; ++physical) {
    long double value = gradient[physical];
    if(inverseJacobian != nullptr) {
      value = 0.0L;
      for(unsigned reference = 0; reference < dimension; ++reference) value = std::fma(static_cast<long double>((*inverseJacobian)[reference][physical]),gradient[reference],value);
    }
    norm = std::hypot(norm,value);
  }
  return norm;
}

long double TensorGradientNorm(const std::vector<long double>& coefficient, const long double* point, unsigned dimension, unsigned extent, const std::vector<std::vector<double>>* inverseJacobian) {
  long double powers[3][5] = {};
  long double gradient[3] = {};
  for(unsigned axis = 0; axis < dimension; ++axis) {
    powers[axis][0] = 1.0L;
    for(unsigned power = 1; power < extent; ++power) powers[axis][power] = powers[axis][power-1U]*point[axis];
  }
  for(unsigned index = 0; index < coefficient.size(); ++index) {
    if(coefficient[index] == 0.0L) continue;
    unsigned exponent[3] = {};
    unsigned remaining = index;
    for(unsigned axis = dimension; axis > 0; --axis) {
      exponent[axis-1U] = remaining%extent;
      remaining /= extent;
    }
    for(unsigned derivative = 0; derivative < dimension; ++derivative) {
      if(exponent[derivative] == 0) continue;
      long double value = coefficient[index]*exponent[derivative];
      for(unsigned axis = 0; axis < dimension; ++axis) value *= powers[axis][exponent[axis]-(axis == derivative ? 1U : 0U)];
      gradient[derivative] += value;
    }
  }
  return SurfaceGradientNorm(gradient,dimension,inverseJacobian);
}

template<int N>
long double ReferenceGradientNorm(const algoim::uvector<double,N>& point, const std::vector<long double>& coefficient, unsigned simplexDimension, const std::vector<std::vector<double>>* inverseJacobian) {
  if(simplexDimension == 2) {
    const long double first = point(0);
    const long double second = first*point(1);
    const long double firstDerivative = coefficient[1]+2.0L*coefficient[3]*first+coefficient[4]*second+coefficient[6]*second*(2.0L*first-second);
    const long double secondDerivative = coefficient[2]+coefficient[4]*first+2.0L*coefficient[5]*second+coefficient[6]*first*(first-2.0L*second);
    if(inverseJacobian == nullptr) return std::hypot(firstDerivative,secondDerivative);
    const long double gradient[2] = {-firstDerivative,secondDerivative};
    return SurfaceGradientNorm(gradient,2,inverseJacobian);
  }
  long double referencePoint[3] = {};
  for(int axis = 0; axis < N; ++axis) referencePoint[axis] = point(axis);
  if(simplexDimension == 3) {
    referencePoint[0] = static_cast<long double>(point(0))*(1.0L-point(1));
    referencePoint[1] = static_cast<long double>(point(0))*point(1)*(1.0L-point(2));
    referencePoint[2] = static_cast<long double>(point(0))*point(1)*point(2);
  }
  return (simplexDimension == 0 ? 0.5L : 1.0L)*TensorGradientNorm(coefficient,referencePoint,N,simplexDimension == 3 ? 5U : 3U,inverseJacobian);
}

template<int N>
void SurfacePoint(const algoim::uvector<double,N>& local, double weight, const algoim::uvector<double,N>& normalWeight, const algoim::xarray<double,N>& polynomial, const algoim::uvector<double,N>& low, double width, long double scale, unsigned simplexDimension, unsigned removedFirstPower, unsigned removedSecondPower, unsigned degree, const std::vector<unsigned>& zeroBoundary, bool referenceSurface, const std::vector<long double>& surfaceCoefficient, const std::vector<std::vector<double>>* inverseJacobian, const std::vector<std::vector<long double>>* cofactorGradient, unsigned cofactorExtent, std::vector<long double>& powers, std::vector<long double>& moments) {
  static_cast<void>(normalWeight);
  algoim::uvector<double,N> point = low+width*local;
  for(unsigned boundary = 0; boundary < zeroBoundary.size(); ++boundary) {
    if(zeroBoundary[boundary] != 0 && std::abs(point(boundary/2U)-static_cast<double>(boundary%2U)) <= 64.0*std::numeric_limits<double>::epsilon()) return;
  }
  algoim::uvector<double,N> gradient = algoim::bernstein::evalBernsteinPolyGradient(polynomial,point);
  long double norm = 0.0L;
  for(int axis = 0; axis < N; ++axis) {
    norm = std::hypot(norm,static_cast<long double>(gradient(axis)));
  }
  if(!std::isfinite(weight) || !std::isfinite(norm) || !(norm > 64.0L*std::numeric_limits<double>::epsilon())) {
    throw std::domain_error("The Algoim Dirac integral requires a numerically regular zero level set.");
  }
  long double mappedWeight = static_cast<long double>(weight)/norm/scale;
  for(int axis = 1; axis < N; ++axis) {
    mappedWeight *= width;
  }
  for(int axis = 0; axis < N; ++axis) {
    if(local(axis) == 0.0 || local(axis) == 1.0) {
      mappedWeight *= 0.5L;
    }
  }
  if(simplexDimension == 2 && removedFirstPower == 0) {
    mappedWeight *= point(0);
  }
  else if(simplexDimension == 3) {
    for(unsigned power = removedFirstPower; power < 2; ++power) mappedWeight *= point(0);
    if(removedSecondPower == 0) mappedWeight *= point(1);
  }
  if(referenceSurface) {
    const long double surfaceNorm = cofactorGradient == nullptr ? ReferenceGradientNorm<N>(point,surfaceCoefficient,simplexDimension,inverseJacobian) : SurfaceCofactorNorm<N>(*cofactorGradient,cofactorExtent,point,simplexDimension);
    if(!(surfaceNorm > 0.0L) || !std::isfinite(surfaceNorm)) throw std::domain_error("The Algoim surface integral requires a regular reference level set.");
    mappedWeight *= surfaceNorm;
  }
  AddMoments<N>(point,mappedWeight,simplexDimension,degree,powers,moments);
}

template<int N>
void AddFullCell(const algoim::uvector<double,N>& low, const algoim::uvector<double,N>& high, unsigned simplexDimension, unsigned degree, std::vector<long double>& moments) {
  unsigned stride = degree+(simplexDimension == 3 ? 3U : 2U);
  std::vector<long double> integrals(static_cast<unsigned>(N)*stride);
  for(int axis = 0; axis < N; ++axis) {
    for(unsigned power = 0; power < stride; ++power) {
      integrals[static_cast<unsigned>(axis)*stride+power] = FELevelSetMoments::MonomialIntegral(low(axis),high(axis),power);
    }
  }
  unsigned index = 0;
  for(unsigned total = 0; total <= degree; ++total) {
    if constexpr(N == 2) {
      for(unsigned j = 0; j <= total; ++j) {
        moments[index] += integrals[simplexDimension == 2 ? total+1U : total-j]*integrals[stride+j];
        ++index;
      }
    }
    else {
      for(int i = static_cast<int>(total); i >= 0; --i) {
        for(int j = static_cast<int>(total)-i; j >= 0; --j) {
          unsigned k = total-static_cast<unsigned>(i+j);
          moments[index] += simplexDimension == 3 ? integrals[total+2U]*integrals[stride+j+k+1U]*integrals[2U*stride+k] : integrals[i]*integrals[stride+j]*integrals[2U*stride+k];
          ++index;
        }
      }
    }
  }
}

template<int N>
bool HasMonotoneDirection(const algoim::xarray<double,N>& polynomial) {
  for(int axis = 0; axis < N; ++axis) {
    unsigned stride = 1;
    for(int other = axis+1; other < N; ++other) {
      stride *= static_cast<unsigned>(polynomial.ext(other));
    }
    bool positive = true;
    bool negative = true;
    bool examined = false;
    for(int index = 0; index < polynomial.size(); ++index) {
      unsigned coordinate = (static_cast<unsigned>(index)/stride)%static_cast<unsigned>(polynomial.ext(axis));
      if(coordinate+1U < static_cast<unsigned>(polynomial.ext(axis))) {
        double derivative = (polynomial[index+stride]-polynomial[index])*(polynomial.ext(axis)-1);
        positive = positive && derivative > 256.0*std::numeric_limits<double>::epsilon();
        negative = negative && derivative < -256.0*std::numeric_limits<double>::epsilon();
        examined = true;
      }
    }
    if(examined && (positive || negative)) {
      return true;
    }
  }
  return false;
}

template<int N>
void IntegrateTensor(const std::vector<long double>& coefficient, algoim::uvector<int,N> extent, unsigned degree, unsigned simplexDimension, unsigned removedFirstPower, unsigned removedSecondPower, std::vector<double>& positive, std::vector<double>& delta, unsigned quadratureOrder, unsigned subdivisions, bool referenceSurface, const std::vector<long double>& surfaceCoefficient, const std::vector<std::vector<double>>* inverseJacobian, const std::vector<std::vector<long double>>* cofactorGradient, unsigned cofactorExtent, const std::vector<unsigned>* prescribedZeroBoundary = nullptr, const std::vector<unsigned>* surfaceZeroBoundary = nullptr) {
  if(&positive == &delta || degree > 64 || quadratureOrder == 0 || quadratureOrder > 100 || quadratureOrder < (degree+(simplexDimension == 0 ? 2U : simplexDimension+1U))/2U || subdivisions > 32) {
    throw std::invalid_argument("Invalid Algoim moment degree, quadrature order, subdivisions, or output vectors.");
  }
  unsigned size = 1;
  for(int axis = 0; axis < N; ++axis) {
    size *= static_cast<unsigned>(extent(axis));
  }
  if(coefficient.size() != size) {
    throw std::invalid_argument("Invalid Algoim polynomial coefficient tensor.");
  }
  long double scale = 0.0L;
  for(unsigned i = 0; i < size; ++i) {
    if(!std::isfinite(coefficient[i])) {
      throw std::invalid_argument("An Algoim polynomial coefficient is not finite.");
    }
    scale = std::max(scale,std::abs(static_cast<long double>(coefficient[i])));
  }
  if(scale == 0.0L) {
    throw std::domain_error("The Dirac delta of an identically zero level set is undefined.");
  }
  std::vector<long double> converted(size);
  for(unsigned i = 0; i < size; ++i) {
    converted[i] = static_cast<long double>(coefficient[i])/scale;
  }
  unsigned stride = size;
  for(int axis = 0; axis < N; ++axis) {
    unsigned count = static_cast<unsigned>(extent(axis));
    stride /= count;
    std::vector<long double> next(size,0.0L);
    for(unsigned i = 0; i < size; ++i) {
      unsigned power = (i/stride)%count;
      unsigned offset = i-power*stride;
      for(unsigned k = 0; k <= power; ++k) {
        next[i] += converted[offset+k*stride]*Binomial(power,k)/Binomial(count-1U,k);
      }
    }
    converted.swap(next);
  }
  std::vector<double> bernstein(size);
  for(unsigned i = 0; i < size; ++i) {
    bernstein[i] = static_cast<double>(converted[i]);
  }
  std::vector<unsigned> zeroBoundary;
  if(prescribedZeroBoundary != nullptr) {
    zeroBoundary = *prescribedZeroBoundary;
    for(unsigned boundary = 0; boundary < zeroBoundary.size(); ++boundary) {
      if(zeroBoundary[boundary] == 0) continue;
      unsigned stride = 1;
      unsigned axis = boundary/2U;
      for(unsigned other = axis+1U; other < static_cast<unsigned>(N); ++other) stride *= extent(other);
      unsigned coordinate = boundary%2U == 0 ? 0U : static_cast<unsigned>(extent(axis)-1);
      for(unsigned index = 0; index < size; ++index) {
        if((index/stride)%static_cast<unsigned>(extent(axis)) == coordinate) bernstein[index] = 0.0;
      }
    }
  }
  std::vector<unsigned> skippedSurfaceBoundary(zeroBoundary);
  if(referenceSurface && surfaceZeroBoundary != nullptr) skippedSurfaceBoundary = *surfaceZeroBoundary;
  algoim::xarray<double,N> polynomial(bernstein.data(),extent);
  unsigned momentCount = (degree+1U)*(degree+2U)/2U;
  if constexpr(N == 3) {
    momentCount = momentCount*(degree+3U)/3U;
  }
  std::vector<long double> h(momentCount,0.0L);
  std::vector<long double> d(momentCount,0.0L);
  std::vector<long double> powers(static_cast<unsigned>(N)*(degree+1U));
  int globalPhase = algoim::bernstein::uniformSign(polynomial);
  if(globalPhase == 0 && !zeroBoundary.empty()) {
    bool nonnegative = true;
    bool nonpositive = true;
    for(unsigned index = 0; index < size; ++index) {
      nonnegative = nonnegative && bernstein[index] >= 0.0;
      nonpositive = nonpositive && bernstein[index] <= 0.0;
    }
    if(nonnegative != nonpositive) globalPhase = nonnegative ? 1 : -1;
  }
  if(subdivisions == 0) {
    subdivisions = globalPhase != 0 || HasMonotoneDirection<N>(polynomial) ? 1 : 4;
  }
  if(globalPhase > 0) {
    AddFullCell<N>(algoim::uvector<double,N>(0.0),algoim::uvector<double,N>(1.0),simplexDimension,degree,h);
  }
  double width = 1.0/static_cast<double>(subdivisions);
  for(algoim::MultiLoop<N> cell(0,globalPhase == 0 ? static_cast<int>(subdivisions) : 0); ~cell; ++cell) {
    algoim::uvector<double,N> low;
    algoim::uvector<double,N> high;
    for(int axis = 0; axis < N; ++axis) {
      low(axis) = static_cast<double>(cell(axis))/subdivisions;
      high(axis) = static_cast<double>(cell(axis)+1)/subdivisions;
    }
    std::vector<double> localData(size);
    algoim::xarray<double,N> localPolynomial(localData.data(),extent);
    algoim::bernstein::deCasteljau(polynomial,low,high,localPolynomial);
    int phase = algoim::bernstein::uniformSign(localPolynomial);
    if(phase > 0) {
      AddFullCell<N>(low,high,simplexDimension,degree,h);
      continue;
    }
    if(phase < 0) {
      continue;
    }
    algoim::ImplicitPolyQuadrature<N> quadrature(localPolynomial);
    quadrature.integrate(algoim::AutoMixed,static_cast<int>(quadratureOrder),std::bind(&VolumePoint<N>,std::placeholders::_1,std::placeholders::_2,std::cref(polynomial),std::cref(low),width,simplexDimension,degree,std::ref(powers),std::ref(h)));
    quadrature.integrate_surf(algoim::AutoMixed,static_cast<int>(quadratureOrder),std::bind(&SurfacePoint<N>,std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::cref(polynomial),std::cref(low),width,scale,simplexDimension,removedFirstPower,removedSecondPower,degree,std::cref(skippedSurfaceBoundary),referenceSurface,std::cref(surfaceCoefficient),inverseJacobian,cofactorGradient,cofactorExtent,std::ref(powers),std::ref(d)));
  }
  std::vector<double> resultH(momentCount);
  std::vector<double> resultD(momentCount);
  for(unsigned i = 0; i < momentCount; ++i) {
    resultH[i] = static_cast<double>(h[i]);
    resultD[i] = static_cast<double>(d[i]);
    if(!std::isfinite(resultH[i]) || !std::isfinite(resultD[i])) {
      throw std::overflow_error("An Algoim moment exceeds double precision range.");
    }
  }
  positive.swap(resultH);
  delta.swap(resultD);
}

}

void IntegrateQuad9Moments(const std::vector<std::vector<double>>& coefficient, unsigned degree, std::vector<double>& positive, std::vector<double>& delta, unsigned quadratureOrder, unsigned subdivisions, bool referenceSurface, const std::vector<std::vector<double>>* inverseJacobian, const std::vector<std::vector<double>>* geometryCoordinates) {
  ValidateGradientTransform(inverseJacobian,2,referenceSurface);
  ValidateSurfaceGeometry(geometryCoordinates,inverseJacobian,2,9,referenceSurface);
  if(coefficient.size() != 3) {
    throw std::invalid_argument("QUAD9 requires a 3 by 3 coefficient tensor.");
  }
  std::vector<long double> data(9);
  for(unsigned i = 0; i < 3; ++i) {
    if(coefficient[i].size() != 3) {
      throw std::invalid_argument("QUAD9 requires a 3 by 3 coefficient tensor.");
    }
    for(unsigned j = 0; j < 3; ++j) {
      data[3U*i+j] = coefficient[i][j];
    }
  }
  std::vector<std::vector<long double>> cofactorGradient;
  unsigned cofactorExtent = 0;
  SurfaceCofactorGradient(data,9,geometryCoordinates,cofactorGradient,cofactorExtent);
  IntegrateTensor<2>(data,algoim::uvector<int,2>(3,3),degree,0,0,0,positive,delta,quadratureOrder,subdivisions,referenceSurface,data,inverseJacobian,geometryCoordinates == nullptr ? nullptr : &cofactorGradient,cofactorExtent);
}

static bool Tri7ExactSumZero(long double first, long double second, long double third) {
  const long double sum = first+second;
  const long double recoveredSecond = sum-first;
  const long double error = (first-(sum-recoveredSecond))+(second-recoveredSecond);
  return sum == -third && error == 0.0L;
}

static void Tri7GeometryZeroEdgeMoments(const std::vector<double>& coefficient, const std::vector<unsigned>& zeroBoundary, const std::vector<std::vector<long double>>& cofactorGradient, unsigned cofactorExtent, unsigned degree, unsigned order, unsigned subdivisions, std::vector<double>& delta) {
  std::vector<long double> moments(delta.size(),0.0L);
  std::vector<long double> powers(2U*(degree+1U));
  const unsigned cells = subdivisions == 0 ? 4U : subdivisions;
  long double scale = 0.0L;
  for(unsigned i = 0; i < coefficient.size(); ++i) scale = std::max(scale,std::abs(static_cast<long double>(coefficient[i])));
  for(unsigned edge = 1; edge < 4; ++edge) {
    if(zeroBoundary[edge] == 0) continue;
    for(unsigned cell = 0; cell < cells; ++cell) {
      for(unsigned node = 0; node < order; ++node) {
        const long double parameter = (cell+algoim::GaussQuad::x(order,node))/cells;
        const long double point[2] = {edge == 1 ? 1.0L : parameter,edge == 2 ? 0.0L : parameter};
        const long double firstDerivative = coefficient[1]+2.0L*coefficient[3]*point[0]+coefficient[4]*point[1]+coefficient[6]*point[1]*(2.0L*point[0]-point[1]);
        const long double secondDerivative = coefficient[2]+coefficient[4]*point[0]+2.0L*coefficient[5]*point[1]+coefficient[6]*point[0]*(point[0]-2.0L*point[1]);
        const long double derivative = std::abs(edge == 1 ? firstDerivative : secondDerivative);
        const long double norm = SurfaceCofactorNorm(cofactorGradient,cofactorExtent,point,2);
        if(!(derivative > 64.0L*std::numeric_limits<double>::epsilon()*scale) || !(norm > 0.0L) || !std::isfinite(norm)) throw std::domain_error("A TRI7 zero edge has a nonregular surface density.");
        const long double weight = 0.5L*algoim::GaussQuad::w(order,node)*norm/(cells*derivative);
        AddMoments<2>(algoim::uvector<double,2>(static_cast<double>(point[0]),static_cast<double>(point[1]/point[0])),weight,2,degree,powers,moments);
      }
    }
  }
  for(unsigned index = 0; index < delta.size(); ++index) {
    delta[index] = static_cast<double>(delta[index]+moments[index]);
    if(!std::isfinite(delta[index])) throw std::overflow_error("A TRI7 edge surface moment exceeds double precision range.");
  }
}

static void IntegrateTri7MomentsInternal(const std::vector<double>& coefficient, unsigned degree, std::vector<double>& positive, std::vector<double>& delta, unsigned quadratureOrder, unsigned subdivisions, bool referenceSurface, const std::vector<unsigned>* surfaceZeroBoundary, const std::vector<std::vector<double>>* inverseJacobian, const std::vector<std::vector<double>>* geometryCoordinates) {
  if(coefficient.size() != 7) {
    throw std::invalid_argument("TRI7 requires coefficients of 1,x,y,x*x,x*y,y*y,x*y*(x-y).");
  }
  if(referenceSurface && surfaceZeroBoundary != nullptr) {
    long double scale = 0.0L;
    for(unsigned i = 0; i < coefficient.size(); ++i) scale = std::max(scale,std::abs(static_cast<long double>(coefficient[i])));
    const long double tolerance = 64.0L*std::numeric_limits<double>::epsilon()*scale;
    const long double firstDerivative = std::abs(static_cast<long double>(coefficient[1])+2.0L*coefficient[3])+std::abs(static_cast<long double>(coefficient[4])+2.0L*coefficient[6])+std::abs(static_cast<long double>(coefficient[6]));
    const long double secondDerivative = std::abs(static_cast<long double>(coefficient[2]))+std::abs(static_cast<long double>(coefficient[4]))+std::abs(static_cast<long double>(coefficient[6]));
    const long double diagonalDerivative = std::abs(static_cast<long double>(coefficient[2]))+std::abs(static_cast<long double>(coefficient[4])+2.0L*coefficient[5])+std::abs(static_cast<long double>(coefficient[6]));
    if(((*surfaceZeroBoundary)[1] != 0 && firstDerivative <= tolerance) || ((*surfaceZeroBoundary)[2] != 0 && secondDerivative <= tolerance) || ((*surfaceZeroBoundary)[3] != 0 && diagonalDerivative <= tolerance)) throw std::domain_error("A TRI7 zero edge has a numerically critical level set.");
  }
  std::vector<long double> data(12,0.0);
  data[0] = coefficient[0];
  data[3] = coefficient[1];
  data[4] = coefficient[2];
  data[6] = coefficient[3];
  data[7] = coefficient[4];
  data[8] = coefficient[5];
  data[10] = coefficient[6];
  data[11] = -coefficient[6];
  unsigned removed = 0;
  unsigned rows = 4;
  while(rows > 1 && data[0] == 0.0 && data[1] == 0.0 && data[2] == 0.0) {
    data.erase(data.begin(),data.begin()+3);
    --rows;
    ++removed;
  }
  if(removed > 1) {
    throw std::domain_error("TRI7 has a critical zero at the collapsed vertex.");
  }
  const std::vector<long double> surfaceCoefficient(referenceSurface ? coefficient.begin() : coefficient.end(),coefficient.end());
  std::vector<std::vector<long double>> cofactorGradient;
  unsigned cofactorExtent = 0;
  SurfaceCofactorGradient(surfaceCoefficient,7,geometryCoordinates,cofactorGradient,cofactorExtent);
  IntegrateTensor<2>(data,algoim::uvector<int,2>(static_cast<int>(rows),3),degree,2,removed,0,positive,delta,quadratureOrder,subdivisions,referenceSurface,surfaceCoefficient,inverseJacobian,geometryCoordinates == nullptr ? nullptr : &cofactorGradient,cofactorExtent,nullptr,surfaceZeroBoundary);
  if(geometryCoordinates != nullptr && surfaceZeroBoundary != nullptr) Tri7GeometryZeroEdgeMoments(coefficient,*surfaceZeroBoundary,cofactorGradient,cofactorExtent,degree,quadratureOrder,subdivisions,delta);
  if(referenceSurface && surfaceZeroBoundary != nullptr && geometryCoordinates == nullptr) {
    const long double firstNormal[2] = {1.0L,0.0L};
    const long double secondNormal[2] = {0.0L,1.0L};
    const long double diagonalNormal[2] = {1.0L,1.0L};
    const double firstWeight = inverseJacobian == nullptr ? 0.5 : static_cast<double>(0.5L*SurfaceGradientNorm(firstNormal,2,inverseJacobian));
    const double secondWeight = inverseJacobian == nullptr ? 0.5 : static_cast<double>(0.5L*SurfaceGradientNorm(secondNormal,2,inverseJacobian));
    const double diagonalWeight = inverseJacobian == nullptr ? std::sqrt(0.5) : static_cast<double>(0.5L*SurfaceGradientNorm(diagonalNormal,2,inverseJacobian));
    unsigned index = 0;
    for(unsigned total = 0; total <= degree; ++total) {
      for(unsigned j = 0; j <= total; ++j) {
        if((*surfaceZeroBoundary)[1] != 0) delta[index] += firstWeight/static_cast<double>(j+1U);
        if((*surfaceZeroBoundary)[2] != 0 && j == 0) delta[index] += secondWeight/static_cast<double>(total+1U);
        if((*surfaceZeroBoundary)[3] != 0) delta[index] += diagonalWeight/static_cast<double>(total+1U);
        ++index;
      }
    }
  }
}

void IntegrateTri7Moments(const std::vector<double>& coefficient, unsigned degree, std::vector<double>& positive, std::vector<double>& delta, unsigned quadratureOrder, unsigned subdivisions, bool referenceSurface, const std::vector<std::vector<double>>* inverseJacobian, const std::vector<std::vector<double>>* geometryCoordinates) {
  ValidateGradientTransform(inverseJacobian,2,referenceSurface);
  ValidateSurfaceGeometry(geometryCoordinates,inverseJacobian,2,7,referenceSurface);
  std::vector<unsigned> surfaceZeroBoundary(4,0);
  if(referenceSurface) {
    if(coefficient.size() != 7) throw std::invalid_argument("TRI7 requires coefficients of 1,x,y,x*x,x*y,y*y,x*y*(x-y).");
    surfaceZeroBoundary[1] = Tri7ExactSumZero(coefficient[0],coefficient[1],coefficient[3]) && Tri7ExactSumZero(coefficient[2],coefficient[4],coefficient[6]) && coefficient[5] == coefficient[6];
    surfaceZeroBoundary[2] = coefficient[0] == 0.0 && coefficient[1] == 0.0 && coefficient[3] == 0.0;
    surfaceZeroBoundary[3] = coefficient[0] == 0.0 && coefficient[1] == -coefficient[2] && Tri7ExactSumZero(coefficient[3],coefficient[4],coefficient[5]);
  }
  IntegrateTri7MomentsInternal(coefficient,degree,positive,delta,quadratureOrder,subdivisions,referenceSurface,referenceSurface ? &surfaceZeroBoundary : nullptr,inverseJacobian,geometryCoordinates);
}

static bool Hex27ExactZeroFaces(const std::vector<long double>& coefficient, std::vector<unsigned>& zeroFaces) {
  zeroFaces.assign(6,0);
  bool any = false;
  for(unsigned axis = 0; axis < 3; ++axis) {
    const unsigned stride = axis == 0 ? 9U : axis == 1 ? 3U : 1U;
    for(unsigned side = 0; side < 2; ++side) {
      bool zero = true;
      for(unsigned index = 0; index < 27; ++index) {
        if((index/stride)%3U != 0) continue;
        zero = zero && (side == 0 ? coefficient[index] == 0.0L : Tri7ExactSumZero(coefficient[index],coefficient[index+stride],coefficient[index+2U*stride]));
      }
      zeroFaces[2U*axis+side] = zero;
      any = any || zero;
    }
  }
  return any;
}

static void Hex27GeometryZeroFaceMoments(const std::vector<std::vector<std::vector<double>>>& coefficient, const std::vector<unsigned>& zeroFaces, unsigned degree, std::vector<double>& positive, std::vector<double>& delta, unsigned quadratureOrder, unsigned subdivisions, const std::vector<std::vector<double>>& geometryCoordinates) {
  if(&positive == &delta || degree > 64 || quadratureOrder == 0 || quadratureOrder > 100 || quadratureOrder < (degree+2U)/2U || subdivisions > 32) throw std::invalid_argument("Invalid Algoim moment degree, quadrature order, subdivisions, or output vectors.");
  std::vector<std::vector<std::vector<double>>> h;
  std::vector<std::vector<std::vector<double>>> d;
  // Native slices avoid Algoim's zero-score recursion on exact HEX boundary faces in the geometry mode.
  FELevelSetMoments::IntegrateHex27UnitSelected(coefficient,degree,h,d,1.0e-12,1.0e-12,true,24,0,200000,nullptr,1,true,&zeroFaces,nullptr,true,nullptr,&geometryCoordinates);
  std::vector<double> resultH;
  std::vector<double> resultD;
  for(unsigned total = 0; total <= degree; ++total) {
    for(int i = static_cast<int>(total); i >= 0; --i) {
      for(int j = static_cast<int>(total)-i; j >= 0; --j) {
        const unsigned k = total-static_cast<unsigned>(i+j);
        resultH.push_back(h[i][j][k]);
        resultD.push_back(d[i][j][k]);
      }
    }
  }
  positive.swap(resultH);
  delta.swap(resultD);
}

void IntegrateHex27Moments(const std::vector<std::vector<std::vector<double>>>& coefficient, unsigned degree, std::vector<double>& positive, std::vector<double>& delta, unsigned quadratureOrder, unsigned subdivisions, bool referenceSurface, const std::vector<std::vector<double>>* inverseJacobian, const std::vector<std::vector<double>>* geometryCoordinates) {
  ValidateGradientTransform(inverseJacobian,3,referenceSurface);
  ValidateSurfaceGeometry(geometryCoordinates,inverseJacobian,3,27,referenceSurface);
  if(coefficient.size() != 3) {
    throw std::invalid_argument("HEX27 requires a 3 by 3 by 3 coefficient tensor.");
  }
  std::vector<long double> data(27);
  for(unsigned i = 0; i < 3; ++i) {
    if(coefficient[i].size() != 3) {
      throw std::invalid_argument("HEX27 requires a 3 by 3 by 3 coefficient tensor.");
    }
    for(unsigned j = 0; j < 3; ++j) {
      if(coefficient[i][j].size() != 3) {
        throw std::invalid_argument("HEX27 requires a 3 by 3 by 3 coefficient tensor.");
      }
      for(unsigned k = 0; k < 3; ++k) {
        data[9U*i+3U*j+k] = coefficient[i][j][k];
      }
    }
  }
  if(geometryCoordinates != nullptr) {
    std::vector<unsigned> zeroFaces;
    if(Hex27ExactZeroFaces(data,zeroFaces)) {
      Hex27GeometryZeroFaceMoments(coefficient,zeroFaces,degree,positive,delta,quadratureOrder,subdivisions,*geometryCoordinates);
      return;
    }
  }
  std::vector<std::vector<long double>> cofactorGradient;
  unsigned cofactorExtent = 0;
  SurfaceCofactorGradient(data,27,geometryCoordinates,cofactorGradient,cofactorExtent);
  IntegrateTensor<3>(data,algoim::uvector<int,3>(3,3,3),degree,0,0,0,positive,delta,quadratureOrder,subdivisions,referenceSurface,data,inverseJacobian,geometryCoordinates == nullptr ? nullptr : &cofactorGradient,cofactorExtent);
}

void Tet15ZeroFaceMoments(const std::vector<double>& nodalValues, const std::vector<unsigned>& zeroFaces, unsigned degree, unsigned order, unsigned subdivisions, std::vector<double>& delta, bool referenceSurface, const std::vector<std::vector<double>>* inverseJacobian, const std::vector<std::vector<long double>>* cofactorGradient, unsigned cofactorExtent) {
  std::vector<long double> coefficient;
  FELevelSetMoments::Tet15PolynomialFromNodalValues(nodalValues,false,coefficient);
  std::vector<long double> moments(delta.size(),0.0L);
  std::vector<long double> powers(3U*(degree+1U));
  unsigned cells = subdivisions == 0 ? 4U : subdivisions;
  for(unsigned face = 0; face < 4; ++face) {
    if(zeroFaces[face] == 0) continue;
    unsigned derivative = face == 0 || face == 3 ? 2U : face-1U;
    for(unsigned cellA = 0; cellA < cells; ++cellA) {
      for(unsigned cellB = 0; cellB < cells; ++cellB) {
        for(unsigned a = 0; a < order; ++a) {
          long double first = (cellA+algoim::GaussQuad::x(order,a))/cells;
          for(unsigned b = 0; b < order; ++b) {
            long double second = (cellB+algoim::GaussQuad::x(order,b))/cells;
            long double point[3] = {first,(1.0L-first)*second,0.0L};
            if(face == 0) point[2] = (1.0L-first)*(1.0L-second);
            if(face == 1) {
              point[0] = 0.0L;
              point[1] = first;
              point[2] = (1.0L-first)*second;
            }
            if(face == 2) {
              point[1] = 0.0L;
              point[2] = (1.0L-first)*second;
            }
            long double basis[3][5];
            for(unsigned axis = 0; axis < 3; ++axis) {
              basis[axis][0] = 1.0L;
              for(unsigned power = 1; power < 5; ++power) basis[axis][power] = basis[axis][power-1U]*point[axis];
            }
            long double normalDerivative = 0.0L;
            for(unsigned i = 0; i < 5; ++i) {
              for(unsigned j = 0; j < 5; ++j) {
                for(unsigned k = 0; k < 5; ++k) {
                  unsigned exponent[3] = {i,j,k};
                  long double value = coefficient[FELevelSetMoments::Tet15Index(i,j,k)]*exponent[derivative];
                  if(exponent[derivative] == 0) continue;
                  --exponent[derivative];
                  normalDerivative += value*basis[0][exponent[0]]*basis[1][exponent[1]]*basis[2][exponent[2]];
                }
              }
            }
            if(normalDerivative == 0.0L || !std::isfinite(normalDerivative)) throw std::domain_error("A TET15 zero face has a critical level set.");
            long double weight = 0.5L*algoim::GaussQuad::w(order,a)*algoim::GaussQuad::w(order,b)*(1.0L-first)/(cells*cells*std::abs(normalDerivative));
            if(referenceSurface) weight *= cofactorGradient == nullptr ? TensorGradientNorm(coefficient,point,3,5,inverseJacobian) : SurfaceCofactorNorm(*cofactorGradient,cofactorExtent,point,3);
            long double u = point[0]+point[1]+point[2];
            long double v = (point[1]+point[2])/u;
            long double w = point[1]+point[2] > 0.0L ? point[2]/(point[1]+point[2]) : 0.0L;
            AddMoments<3>(algoim::uvector<double,3>(u,v,w),weight,3,degree,powers,moments);
          }
        }
      }
    }
  }
  for(unsigned index = 0; index < delta.size(); ++index) {
    delta[index] = static_cast<double>(delta[index]+moments[index]);
    if(!std::isfinite(delta[index])) throw std::overflow_error("A TET15 face Dirac moment exceeds double precision range.");
  }
}

void IntegrateTet15NodalMoments(const std::vector<double>& nodalValues, unsigned degree, std::vector<double>& positive, std::vector<double>& delta, unsigned quadratureOrder, unsigned subdivisions, bool referenceSurface, const std::vector<std::vector<double>>* inverseJacobian, const std::vector<std::vector<double>>* geometryCoordinates) {
  ValidateGradientTransform(inverseJacobian,3,referenceSurface);
  ValidateSurfaceGeometry(geometryCoordinates,inverseJacobian,3,15,referenceSurface);
  std::vector<long double> polynomial;
  FELevelSetMoments::Tet15PolynomialFromNodalValues(nodalValues,true,polynomial);
  std::vector<long double> surfaceCoefficient;
  if(referenceSurface) FELevelSetMoments::Tet15PolynomialFromNodalValues(nodalValues,false,surfaceCoefficient);
  const unsigned faceNodes[4][7] = {{1,2,3,5,9,8,12},{0,2,3,6,9,7,13},{0,1,3,4,8,7,11},{0,1,2,4,5,6,10}};
  const unsigned faceBoundary[4] = {1,3,5,4};
  std::vector<unsigned> zeroFaces(4,0);
  std::vector<unsigned> zeroBoundary;
  unsigned numberOfZeroFaces = 0;
  for(unsigned face = 0; face < 4; ++face) {
    bool zero = true;
    for(unsigned node = 0; node < 7; ++node) zero = zero && nodalValues[faceNodes[face][node]] == 0.0;
    if(zero) {
      zeroFaces[face] = 1;
      if(zeroBoundary.empty()) zeroBoundary.assign(6,0);
      zeroBoundary[faceBoundary[face]] = 1;
      ++numberOfZeroFaces;
    }
  }
  if(numberOfZeroFaces > 1) throw std::domain_error("TET15 has intersecting zero faces and a nonregular Dirac integral.");
  unsigned removed[2] = {0,0};
  for(unsigned axis = 0; axis < 2; ++axis) {
    unsigned stride = axis == 0 ? 25U : 5U;
    bool factor = true;
    while(factor && removed[axis] < 4) {
      for(unsigned index = 0; index < 125; ++index) {
        if((index/stride)%5U == 0 && polynomial[index] != 0.0L) factor = false;
      }
      if(factor) {
        std::vector<long double> next(125,0.0L);
        for(unsigned index = 0; index < 125; ++index) {
          if((index/stride)%5U != 0) next[index-stride] = polynomial[index];
        }
        polynomial.swap(next);
        ++removed[axis];
      }
    }
  }
  if(removed[0] > 1 || removed[1] > 1) {
    throw std::domain_error("TET15 has a critical zero on a collapsed simplex boundary.");
  }
  // Algoim's elevated derivative requires at least two coefficients along every axis.
  unsigned maximum[3] = {1,1,1};
  for(unsigned i = 0; i < 5; ++i) {
    for(unsigned j = 0; j < 5; ++j) {
      for(unsigned k = 0; k < 5; ++k) {
        if(polynomial[FELevelSetMoments::Tet15Index(i,j,k)] != 0.0L) {
          maximum[0] = std::max(maximum[0],i);
          maximum[1] = std::max(maximum[1],j);
          maximum[2] = std::max(maximum[2],k);
        }
      }
    }
  }
  std::vector<long double> data((maximum[0]+1U)*(maximum[1]+1U)*(maximum[2]+1U));
  unsigned index = 0;
  for(unsigned i = 0; i <= maximum[0]; ++i) {
    for(unsigned j = 0; j <= maximum[1]; ++j) {
      for(unsigned k = 0; k <= maximum[2]; ++k) {
        data[index++] = polynomial[FELevelSetMoments::Tet15Index(i,j,k)];
      }
    }
  }
  std::vector<std::vector<long double>> cofactorGradient;
  unsigned cofactorExtent = 0;
  SurfaceCofactorGradient(surfaceCoefficient,15,geometryCoordinates,cofactorGradient,cofactorExtent);
  IntegrateTensor<3>(data,algoim::uvector<int,3>(maximum[0]+1U,maximum[1]+1U,maximum[2]+1U),degree,3,removed[0],removed[1],positive,delta,quadratureOrder,subdivisions,referenceSurface,surfaceCoefficient,inverseJacobian,geometryCoordinates == nullptr ? nullptr : &cofactorGradient,cofactorExtent,zeroBoundary.empty() ? nullptr : &zeroBoundary);
  if(numberOfZeroFaces != 0) Tet15ZeroFaceMoments(nodalValues,zeroFaces,degree,quadratureOrder,subdivisions,delta,referenceSurface,inverseJacobian,geometryCoordinates == nullptr ? nullptr : &cofactorGradient,cofactorExtent);
}

void IntegrateNodalMoments(const std::vector<double>& nodalValues, unsigned degree, std::vector<double>& positive, std::vector<double>& delta, unsigned quadratureOrder, unsigned subdivisions, bool referenceSurface, const std::vector<std::vector<double>>* inverseJacobian, const std::vector<std::vector<double>>* geometryCoordinates) {
  if(&positive == &delta || degree > 64 || quadratureOrder == 0 || quadratureOrder > 100 || subdivisions > 32) {
    throw std::invalid_argument("Invalid Algoim nodal integration request.");
  }
  if(nodalValues.size() == 27) {
    std::vector<std::vector<std::vector<double>>> coefficient;
    FELevelSetMoments::Hex27UnitCoefficientsFromNodalValues(nodalValues,coefficient);
    if(geometryCoordinates != nullptr) {
      ValidateSurfaceGeometry(geometryCoordinates,inverseJacobian,3,27,referenceSurface);
      std::vector<unsigned> zeroFaces(6,0);
      bool any = false;
      for(unsigned axis = 0; axis < 3; ++axis) {
        for(unsigned side = 0; side < 2; ++side) {
          std::vector<double> face;
          FELevelSetMoments::Hex27FaceNodalValues(nodalValues,axis,side,face);
          bool zero = true;
          for(unsigned node = 0; node < 9; ++node) zero = zero && face[node] == 0.0;
          zeroFaces[2U*axis+side] = zero;
          any = any || zero;
        }
      }
      if(any) {
        Hex27GeometryZeroFaceMoments(coefficient,zeroFaces,degree,positive,delta,quadratureOrder,subdivisions,*geometryCoordinates);
        return;
      }
    }
    IntegrateHex27Moments(coefficient,degree,positive,delta,quadratureOrder,subdivisions,referenceSurface,inverseJacobian,geometryCoordinates);
  }
  else if(nodalValues.size() == 15) {
    IntegrateTet15NodalMoments(nodalValues,degree,positive,delta,quadratureOrder,subdivisions,referenceSurface,inverseJacobian,geometryCoordinates);
  }
  else if(nodalValues.size() == 9) {
    std::vector<std::vector<double>> coefficient;
    FELevelSetMoments::Quad9UnitCoefficientsFromNodalValues(nodalValues,coefficient);
    IntegrateQuad9Moments(coefficient,degree,positive,delta,quadratureOrder,subdivisions,referenceSurface,inverseJacobian,geometryCoordinates);
  }
  else if(nodalValues.size() == 7) {
    ValidateGradientTransform(inverseJacobian,2,referenceSurface);
    ValidateSurfaceGeometry(geometryCoordinates,inverseJacobian,2,7,referenceSurface);
    std::vector<double> coefficient;
    FELevelSetMoments::Tri7ReflectedCoefficientsFromNodalValues(nodalValues,coefficient);
    std::vector<unsigned> surfaceZeroBoundary(4,0);
    if(referenceSurface) {
      surfaceZeroBoundary[1] = nodalValues[0] == 0.0 && nodalValues[2] == 0.0 && nodalValues[5] == 0.0;
      surfaceZeroBoundary[2] = nodalValues[0] == 0.0 && nodalValues[1] == 0.0 && nodalValues[3] == 0.0;
      surfaceZeroBoundary[3] = nodalValues[1] == 0.0 && nodalValues[2] == 0.0 && nodalValues[4] == 0.0;
    }
    IntegrateTri7MomentsInternal(coefficient,degree,positive,delta,quadratureOrder,subdivisions,referenceSurface,referenceSurface ? &surfaceZeroBoundary : nullptr,inverseJacobian,geometryCoordinates);
  }
  else if(nodalValues.size() == 3) {
    ValidateGradientTransform(inverseJacobian,1,referenceSurface);
    ValidateSurfaceGeometry(geometryCoordinates,inverseJacobian,1,3,referenceSurface);
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;
    FELevelSetMoments::Line3CoefficientsFromNodalValues(nodalValues,a,b,c);
    FELevelSetMoments::Integrate1D(a,b,c,degree,positive,delta,0.0,1.0,!referenceSurface);
    if(referenceSurface) Line3SurfaceMoments(a,b,c,degree,delta,nodalValues[0] == 0.0,nodalValues[1] == 0.0);
    if(inverseJacobian != nullptr) {
      for(unsigned power = 0; power <= degree; ++power) delta[power] *= std::abs((*inverseJacobian)[0][0]);
    }
  }
  else {
    throw std::invalid_argument("Algoim nodal integration supports LINE3, QUAD9, TRI7, HEX27, and TET15.");
  }
}

}

