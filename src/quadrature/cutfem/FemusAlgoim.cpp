#include "FemusAlgoim.hpp"
#include "FELevelSetMoments.hpp"
#include "Hex27LevelSetMoments.hpp"
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

template<int N>
void SurfacePoint(const algoim::uvector<double,N>& local, double weight, const algoim::uvector<double,N>& normalWeight, const algoim::xarray<double,N>& polynomial, const algoim::uvector<double,N>& low, double width, long double scale, unsigned simplexDimension, unsigned removedFirstPower, unsigned removedSecondPower, unsigned degree, const std::vector<unsigned>& zeroBoundary, std::vector<long double>& powers, std::vector<long double>& moments) {
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
void IntegrateTensor(const std::vector<long double>& coefficient, algoim::uvector<int,N> extent, unsigned degree, unsigned simplexDimension, unsigned removedFirstPower, unsigned removedSecondPower, std::vector<double>& positive, std::vector<double>& delta, unsigned quadratureOrder, unsigned subdivisions, const std::vector<unsigned>* prescribedZeroBoundary = nullptr) {
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
    quadrature.integrate_surf(algoim::AutoMixed,static_cast<int>(quadratureOrder),std::bind(&SurfacePoint<N>,std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::cref(polynomial),std::cref(low),width,scale,simplexDimension,removedFirstPower,removedSecondPower,degree,std::cref(zeroBoundary),std::ref(powers),std::ref(d)));
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

void IntegrateQuad9Moments(const std::vector<std::vector<double>>& coefficient, unsigned degree, std::vector<double>& positive, std::vector<double>& delta, unsigned quadratureOrder, unsigned subdivisions) {
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
  IntegrateTensor<2>(data,algoim::uvector<int,2>(3,3),degree,0,0,0,positive,delta,quadratureOrder,subdivisions);
}

void IntegrateTri7Moments(const std::vector<double>& coefficient, unsigned degree, std::vector<double>& positive, std::vector<double>& delta, unsigned quadratureOrder, unsigned subdivisions) {
  if(coefficient.size() != 7) {
    throw std::invalid_argument("TRI7 requires coefficients of 1,x,y,x*x,x*y,y*y,x*y*(x-y).");
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
  IntegrateTensor<2>(data,algoim::uvector<int,2>(static_cast<int>(rows),3),degree,2,removed,0,positive,delta,quadratureOrder,subdivisions);
}

void IntegrateHex27Moments(const std::vector<std::vector<std::vector<double>>>& coefficient, unsigned degree, std::vector<double>& positive, std::vector<double>& delta, unsigned quadratureOrder, unsigned subdivisions) {
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
  IntegrateTensor<3>(data,algoim::uvector<int,3>(3,3,3),degree,0,0,0,positive,delta,quadratureOrder,subdivisions);
}

void Tet15ZeroFaceMoments(const std::vector<double>& nodalValues, const std::vector<unsigned>& zeroFaces, unsigned degree, unsigned order, unsigned subdivisions, std::vector<double>& delta) {
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

void IntegrateTet15NodalMoments(const std::vector<double>& nodalValues, unsigned degree, std::vector<double>& positive, std::vector<double>& delta, unsigned quadratureOrder, unsigned subdivisions) {
  std::vector<long double> polynomial;
  FELevelSetMoments::Tet15PolynomialFromNodalValues(nodalValues,true,polynomial);
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
  IntegrateTensor<3>(data,algoim::uvector<int,3>(maximum[0]+1U,maximum[1]+1U,maximum[2]+1U),degree,3,removed[0],removed[1],positive,delta,quadratureOrder,subdivisions,zeroBoundary.empty() ? nullptr : &zeroBoundary);
  if(numberOfZeroFaces != 0) Tet15ZeroFaceMoments(nodalValues,zeroFaces,degree,quadratureOrder,subdivisions,delta);
}

void IntegrateNodalMoments(const std::vector<double>& nodalValues, unsigned degree, std::vector<double>& positive, std::vector<double>& delta, unsigned quadratureOrder, unsigned subdivisions) {
  if(&positive == &delta || degree > 64 || quadratureOrder == 0 || quadratureOrder > 100 || subdivisions > 32) {
    throw std::invalid_argument("Invalid Algoim nodal integration request.");
  }
  if(nodalValues.size() == 27) {
    std::vector<std::vector<std::vector<double>>> coefficient;
    FELevelSetMoments::Hex27UnitCoefficientsFromNodalValues(nodalValues,coefficient);
    IntegrateHex27Moments(coefficient,degree,positive,delta,quadratureOrder,subdivisions);
  }
  else if(nodalValues.size() == 15) {
    IntegrateTet15NodalMoments(nodalValues,degree,positive,delta,quadratureOrder,subdivisions);
  }
  else if(nodalValues.size() == 9) {
    std::vector<std::vector<double>> coefficient;
    FELevelSetMoments::Quad9UnitCoefficientsFromNodalValues(nodalValues,coefficient);
    IntegrateQuad9Moments(coefficient,degree,positive,delta,quadratureOrder,subdivisions);
  }
  else if(nodalValues.size() == 7) {
    std::vector<double> coefficient;
    FELevelSetMoments::Tri7ReflectedCoefficientsFromNodalValues(nodalValues,coefficient);
    IntegrateTri7Moments(coefficient,degree,positive,delta,quadratureOrder,subdivisions);
  }
  else if(nodalValues.size() == 3) {
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;
    FELevelSetMoments::Line3CoefficientsFromNodalValues(nodalValues,a,b,c);
    FELevelSetMoments::Integrate1D(a,b,c,degree,positive,delta);
  }
  else {
    throw std::invalid_argument("Algoim nodal integration supports LINE3, QUAD9, TRI7, HEX27, and TET15.");
  }
}

}
