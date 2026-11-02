#include "FemusAlgoim.hpp"
#include "FELevelSetMoments.hpp"
#include "Hex27LevelSetMoments.hpp"
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
void AddMoments(const algoim::uvector<double,N>& point, long double weight, bool triangle, unsigned degree, std::vector<long double>& powers, std::vector<long double>& moments) {
  for(int axis = 0; axis < N; ++axis) {
    long double coordinate = point(axis);
    if(triangle && axis == 1) {
      coordinate *= point(0);
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
void VolumePoint(const algoim::uvector<double,N>& local, double weight, const algoim::xarray<double,N>& polynomial, const algoim::uvector<double,N>& low, double width, bool triangle, unsigned degree, std::vector<long double>& powers, std::vector<long double>& moments) {
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
  if(triangle) {
    mappedWeight *= point(0);
  }
  AddMoments<N>(point,mappedWeight,triangle,degree,powers,moments);
}

template<int N>
void SurfacePoint(const algoim::uvector<double,N>& local, double weight, const algoim::uvector<double,N>& normalWeight, const algoim::xarray<double,N>& polynomial, const algoim::uvector<double,N>& low, double width, long double scale, bool triangle, unsigned removedTrianglePower, unsigned degree, std::vector<long double>& powers, std::vector<long double>& moments) {
  static_cast<void>(normalWeight);
  algoim::uvector<double,N> point = low+width*local;
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
  if(triangle && removedTrianglePower == 0) {
    mappedWeight *= point(0);
  }
  AddMoments<N>(point,mappedWeight,triangle,degree,powers,moments);
}

template<int N>
void AddFullCell(const algoim::uvector<double,N>& low, const algoim::uvector<double,N>& high, bool triangle, unsigned degree, std::vector<long double>& moments) {
  unsigned stride = degree+2U;
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
        moments[index] += integrals[triangle ? total+1U : total-j]*integrals[stride+j];
        ++index;
      }
    }
    else {
      for(int i = static_cast<int>(total); i >= 0; --i) {
        for(int j = static_cast<int>(total)-i; j >= 0; --j) {
          unsigned k = total-static_cast<unsigned>(i+j);
          moments[index] += integrals[i]*integrals[stride+j]*integrals[2U*stride+k];
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
void IntegrateTensor(const std::vector<double>& coefficient, algoim::uvector<int,N> extent, unsigned degree, bool triangle, unsigned removedTrianglePower, std::vector<double>& positive, std::vector<double>& delta, unsigned quadratureOrder, unsigned subdivisions) {
  if(&positive == &delta || degree > 64 || quadratureOrder == 0 || quadratureOrder > 100 || quadratureOrder < (degree+(triangle ? 3U : 2U))/2U || subdivisions > 32) {
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
  algoim::xarray<double,N> polynomial(bernstein.data(),extent);
  unsigned momentCount = (degree+1U)*(degree+2U)/2U;
  if constexpr(N == 3) {
    momentCount = momentCount*(degree+3U)/3U;
  }
  std::vector<long double> h(momentCount,0.0L);
  std::vector<long double> d(momentCount,0.0L);
  std::vector<long double> powers(static_cast<unsigned>(N)*(degree+1U));
  int globalPhase = algoim::bernstein::uniformSign(polynomial);
  if(subdivisions == 0) {
    subdivisions = globalPhase != 0 || HasMonotoneDirection<N>(polynomial) ? 1 : 4;
  }
  if(globalPhase > 0) {
    AddFullCell<N>(algoim::uvector<double,N>(0.0),algoim::uvector<double,N>(1.0),triangle,degree,h);
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
      AddFullCell<N>(low,high,triangle,degree,h);
      continue;
    }
    if(phase < 0) {
      continue;
    }
    algoim::ImplicitPolyQuadrature<N> quadrature(localPolynomial);
    quadrature.integrate(algoim::AutoMixed,static_cast<int>(quadratureOrder),std::bind(&VolumePoint<N>,std::placeholders::_1,std::placeholders::_2,std::cref(polynomial),std::cref(low),width,triangle,degree,std::ref(powers),std::ref(h)));
    quadrature.integrate_surf(algoim::AutoMixed,static_cast<int>(quadratureOrder),std::bind(&SurfacePoint<N>,std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::cref(polynomial),std::cref(low),width,scale,triangle,removedTrianglePower,degree,std::ref(powers),std::ref(d)));
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
  std::vector<double> data(9);
  for(unsigned i = 0; i < 3; ++i) {
    if(coefficient[i].size() != 3) {
      throw std::invalid_argument("QUAD9 requires a 3 by 3 coefficient tensor.");
    }
    for(unsigned j = 0; j < 3; ++j) {
      data[3U*i+j] = coefficient[i][j];
    }
  }
  IntegrateTensor<2>(data,algoim::uvector<int,2>(3,3),degree,false,0,positive,delta,quadratureOrder,subdivisions);
}

void IntegrateTri7Moments(const std::vector<double>& coefficient, unsigned degree, std::vector<double>& positive, std::vector<double>& delta, unsigned quadratureOrder, unsigned subdivisions) {
  if(coefficient.size() != 7) {
    throw std::invalid_argument("TRI7 requires coefficients of 1,x,y,x*x,x*y,y*y,x*y*(x-y).");
  }
  std::vector<double> data(12,0.0);
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
  IntegrateTensor<2>(data,algoim::uvector<int,2>(static_cast<int>(rows),3),degree,true,removed,positive,delta,quadratureOrder,subdivisions);
}

void IntegrateHex27Moments(const std::vector<std::vector<std::vector<double>>>& coefficient, unsigned degree, std::vector<double>& positive, std::vector<double>& delta, unsigned quadratureOrder, unsigned subdivisions) {
  if(coefficient.size() != 3) {
    throw std::invalid_argument("HEX27 requires a 3 by 3 by 3 coefficient tensor.");
  }
  std::vector<double> data(27);
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
  IntegrateTensor<3>(data,algoim::uvector<int,3>(3,3,3),degree,false,0,positive,delta,quadratureOrder,subdivisions);
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
    throw std::invalid_argument("Algoim nodal integration supports LINE3, QUAD9, TRI7, and HEX27.");
  }
}

}
