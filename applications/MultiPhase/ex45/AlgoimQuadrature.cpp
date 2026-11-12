#include "AlgoimQuadrature.hpp"
#include "algoim/quadrature_multipoly.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

namespace DirectAlgoim {
namespace {

long double Binomial(unsigned n, unsigned k) {
  long double value = 1.0L;
  for(unsigned i = 1; i <= k; ++i) value *= static_cast<long double>(n+1U-i)/i;
  return value;
}

unsigned IntegerPower(unsigned base, unsigned exponent) {
  unsigned value = 1;
  for(unsigned i = 0; i < exponent; ++i) value *= base;
  return value;
}

void TensorPolynomial(const std::vector<long double>& nodalValues, unsigned dimension, bool unitCoordinates, std::vector<long double>& coefficient) {
  const unsigned line[3] = {0,2,1};
  const unsigned quad[9][2] = {{0,0},{2,0},{2,2},{0,2},{1,0},{2,1},{1,2},{0,1},{1,1}};
  const unsigned hex[27][3] = {{0,0,0},{2,0,0},{2,2,0},{0,2,0},{0,0,2},{2,0,2},{2,2,2},{0,2,2},{1,0,0},{2,1,0},{1,2,0},{0,1,0},{1,0,2},{2,1,2},{1,2,2},{0,1,2},{0,0,1},{2,0,1},{2,2,1},{0,2,1},{1,0,1},{2,1,1},{1,2,1},{0,1,1},{1,1,0},{1,1,2},{1,1,1}};
  const long double unit[3][3] = {{1,0,0},{-3,4,-1},{2,-4,2}};
  const long double standard[3][3] = {{0,1,0},{-0.5L,0,0.5L},{0.5L,-1,0.5L}};
  coefficient.assign(nodalValues.size(),0.0L);
  for(unsigned node = 0; node < nodalValues.size(); ++node) {
    unsigned index = dimension == 1 ? line[node] : dimension == 2 ? 3U*quad[node][0]+quad[node][1] : 9U*hex[node][0]+3U*hex[node][1]+hex[node][2];
    coefficient[index] = nodalValues[node];
  }
  for(unsigned axis = 0; axis < dimension; ++axis) {
    unsigned stride = IntegerPower(3,dimension-axis-1U);
    std::vector<long double> next(coefficient.size(),0.0L);
    for(unsigned index = 0; index < coefficient.size(); ++index) {
      unsigned power = (index/stride)%3U;
      unsigned offset = index-power*stride;
      for(unsigned k = 0; k < 3; ++k) next[index] += (unitCoordinates ? unit[power][k] : standard[power][k])*coefficient[offset+k*stride];
    }
    coefficient.swap(next);
  }
}

void TriPolynomial(const std::vector<long double>& values, std::vector<long double>& original, std::vector<long double>& chart) {
  const long double bubble = 27.0L*values[6]+3.0L*(values[0]+values[1]+values[2])-12.0L*(values[3]+values[4]+values[5]);
  original.assign(9,0.0L);
  original[0] = values[0];
  original[3] = -3.0L*values[0]-values[1]+4.0L*values[3];
  original[1] = -3.0L*values[0]-values[2]+4.0L*values[5];
  original[6] = 2.0L*values[0]+2.0L*values[1]-4.0L*values[3];
  original[2] = 2.0L*values[0]+2.0L*values[2]-4.0L*values[5];
  original[4] = 4.0L*values[0]-4.0L*values[3]+4.0L*values[4]-4.0L*values[5]+bubble;
  original[7] = -bubble;
  original[5] = -bubble;
  chart.assign(16,0.0L);
  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = 0; j < 3; ++j) {
      if(original[3U*i+j] == 0.0L) continue;
      for(unsigned k = 0; k <= i; ++k) chart[4U*(j+k)+j] += original[3U*i+j]*Binomial(i,k)*(k%2U == 0 ? 1.0L : -1.0L);
    }
  }
}

void AddTetProduct(const std::vector<std::vector<long double>>& barycentric, const std::vector<unsigned>& factors, long double multiplier, std::vector<long double>& result) {
  if(multiplier == 0.0L) return;
  std::vector<long double> product(125,0.0L);
  product[0] = 1.0L;
  for(unsigned factor = 0; factor < factors.size(); ++factor) {
    std::vector<long double> next(125,0.0L);
    for(unsigned i = 0; i < 5; ++i) {
      for(unsigned j = 0; j < 5; ++j) {
        for(unsigned k = 0; k < 5; ++k) {
          if(product[25U*i+5U*j+k] == 0.0L) continue;
          for(unsigned a = 0; a+i < 5; ++a) {
            for(unsigned b = 0; b+j < 5; ++b) {
              for(unsigned c = 0; c+k < 5; ++c) next[25U*(i+a)+5U*(j+b)+k+c] += product[25U*i+5U*j+k]*barycentric[factors[factor]][25U*a+5U*b+c];
            }
          }
        }
      }
    }
    product.swap(next);
  }
  for(unsigned index = 0; index < 125; ++index) result[index] += multiplier*product[index];
}

void TetPolynomial(const std::vector<long double>& values, std::vector<long double>& original, std::vector<long double>& chart) {
  std::vector<std::vector<long double>> barycentric(4,std::vector<long double>(125,0.0L));
  barycentric[0][0] = 1.0L;
  barycentric[0][25] = -1.0L;
  barycentric[0][5] = -1.0L;
  barycentric[0][1] = -1.0L;
  barycentric[1][25] = 1.0L;
  barycentric[2][5] = 1.0L;
  barycentric[3][1] = 1.0L;
  const unsigned edges[6][2] = {{0,1},{1,2},{0,2},{0,3},{1,3},{2,3}};
  const unsigned faces[4][3] = {{0,1,2},{0,1,3},{1,2,3},{0,2,3}};
  original.assign(125,0.0L);
  long double interior = 256.0L*values[14];
  for(unsigned vertex = 0; vertex < 4; ++vertex) {
    AddTetProduct(barycentric,{vertex,vertex},2.0L*values[vertex],original);
    AddTetProduct(barycentric,{vertex},-values[vertex],original);
    interior -= 4.0L*values[vertex];
  }
  for(unsigned edge = 0; edge < 6; ++edge) {
    AddTetProduct(barycentric,{edges[edge][0],edges[edge][1]},4.0L*values[4U+edge],original);
    interior += 32.0L*values[4U+edge];
  }
  for(unsigned face = 0; face < 4; ++face) {
    long double value = 27.0L*values[10U+face];
    for(unsigned vertex = 0; vertex < 3; ++vertex) value += 3.0L*values[faces[face][vertex]];
    for(unsigned edge = 0; edge < 6; ++edge) {
      unsigned contained = 0;
      for(unsigned vertex = 0; vertex < 3; ++vertex) contained += edges[edge][0] == faces[face][vertex] || edges[edge][1] == faces[face][vertex];
      if(contained == 2) value -= 12.0L*values[4U+edge];
    }
    AddTetProduct(barycentric,{faces[face][0],faces[face][1],faces[face][2]},value,original);
    interior -= 108.0L*values[10U+face];
  }
  AddTetProduct(barycentric,{0,1,2,3},interior,original);
  chart.assign(125,0.0L);
  for(unsigned i = 0; i < 5; ++i) {
    for(unsigned j = 0; i+j < 5; ++j) {
      for(unsigned k = 0; i+j+k < 5; ++k) {
        if(original[25U*i+5U*j+k] == 0.0L) continue;
        for(unsigned a = 0; a <= i; ++a) {
          for(unsigned b = 0; b <= j; ++b) chart[25U*(i+j+k)+5U*(j+k+a)+k+b] += original[25U*i+5U*j+k]*Binomial(i,a)*Binomial(j,b)*((a+b)%2U == 0 ? 1.0L : -1.0L);
        }
      }
    }
  }
}

void BoundaryFlags(const std::vector<long double>& nodalValues, unsigned dimension, unsigned simplexDimension, std::vector<unsigned>& chartZero, std::vector<unsigned>& physicalZero) {
  chartZero.assign(2U*dimension,0);
  if(simplexDimension == 2) {
    const unsigned edge[3][3] = {{0,2,5},{0,1,3},{1,2,4}};
    physicalZero.assign(3,1);
    for(unsigned face = 0; face < 3; ++face) {
      for(unsigned node = 0; node < 3; ++node) physicalZero[face] = physicalZero[face] && nodalValues[edge[face][node]] == 0.0L;
    }
    chartZero[0] = nodalValues[1] == 0.0L;
    chartZero[1] = physicalZero[0];
    chartZero[2] = physicalZero[1];
    chartZero[3] = physicalZero[2];
  }
  else if(simplexDimension == 3) {
    const unsigned face[4][7] = {{1,2,3,5,9,8,12},{0,2,3,6,9,7,13},{0,1,3,4,8,7,11},{0,1,2,4,5,6,10}};
    physicalZero.assign(4,1);
    for(unsigned boundary = 0; boundary < 4; ++boundary) {
      for(unsigned node = 0; node < 7; ++node) physicalZero[boundary] = physicalZero[boundary] && nodalValues[face[boundary][node]] == 0.0L;
    }
    chartZero[0] = nodalValues[0] == 0.0L;
    chartZero[1] = physicalZero[0];
    chartZero[2] = nodalValues[0] == 0.0L && nodalValues[1] == 0.0L && nodalValues[4] == 0.0L;
    chartZero[3] = physicalZero[1];
    chartZero[4] = physicalZero[3];
    chartZero[5] = physicalZero[2];
  }
  else {
    const unsigned line[3] = {0,2,1};
    const unsigned quad[9][2] = {{0,0},{2,0},{2,2},{0,2},{1,0},{2,1},{1,2},{0,1},{1,1}};
    const unsigned hex[27][3] = {{0,0,0},{2,0,0},{2,2,0},{0,2,0},{0,0,2},{2,0,2},{2,2,2},{0,2,2},{1,0,0},{2,1,0},{1,2,0},{0,1,0},{1,0,2},{2,1,2},{1,2,2},{0,1,2},{0,0,1},{2,0,1},{2,2,1},{0,2,1},{1,0,1},{2,1,1},{1,2,1},{0,1,1},{1,1,0},{1,1,2},{1,1,1}};
    physicalZero.assign(2U*dimension,1);
    for(unsigned boundary = 0; boundary < 2U*dimension; ++boundary) {
      for(unsigned node = 0; node < nodalValues.size(); ++node) {
        unsigned coordinate = dimension == 1 ? line[node] : dimension == 2 ? quad[node][boundary/2U] : hex[node][boundary/2U];
        if(coordinate == 2U*(boundary%2U)) physicalZero[boundary] = physicalZero[boundary] && nodalValues[node] == 0.0L;
      }
    }
    chartZero = physicalZero;
  }
}

bool ZeroPowerFace(const std::vector<long double>& coefficient, unsigned dimension, unsigned extent, unsigned boundary) {
  const unsigned stride = IntegerPower(extent,dimension-boundary/2U-1U);
  for(unsigned index = 0; index < coefficient.size(); ++index) {
    if((index/stride)%extent != 0) continue;
    long double sum = coefficient[index];
    if(boundary%2U != 0) {
      for(unsigned power = 1; power < extent; ++power) sum += coefficient[index+power*stride];
    }
    if(sum != 0.0L) return false;
  }
  return true;
}

void RemovePowerFace(std::vector<long double>& coefficient, unsigned dimension, unsigned extent, unsigned boundary) {
  const unsigned stride = IntegerPower(extent,dimension-boundary/2U-1U);
  for(unsigned index = 0; index < coefficient.size(); ++index) {
    if((index/stride)%extent != 0) continue;
    if(boundary%2U == 0) {
      for(unsigned power = 0; power+1U < extent; ++power) coefficient[index+power*stride] = coefficient[index+(power+1U)*stride];
    }
    else {
      for(unsigned power = 1; power+1U < extent; ++power) coefficient[index+power*stride] += coefficient[index+(power-1U)*stride];
    }
    coefficient[index+(extent-1U)*stride] = 0.0L;
  }
}

void RemoveBoundaryFactors(std::vector<long double>& coefficient, unsigned dimension, unsigned extent, std::vector<unsigned>& removed) {
  for(unsigned boundary = 0; boundary < 2U*dimension; ++boundary) {
    bool forced = removed[boundary] != 0;
    removed[boundary] = 0;
    for(unsigned count = 0; count+1U < extent && (forced || ZeroPowerFace(coefficient,dimension,extent,boundary)); ++count) {
      RemovePowerFace(coefficient,dimension,extent,boundary);
      ++removed[boundary];
      forced = false;
    }
  }
}

void SurfaceNormal(const std::vector<double>& point, const std::vector<long double>& coefficient, unsigned extent, std::vector<double>& normal) {
  const unsigned dimension = point.size();
  std::vector<std::vector<long double>> powers(dimension,std::vector<long double>(extent,1.0L));
  for(unsigned axis = 0; axis < dimension; ++axis) {
    for(unsigned degree = 1; degree < extent; ++degree) powers[axis][degree] = powers[axis][degree-1U]*point[axis];
  }
  std::vector<long double> gradient(dimension,0.0L);
  long double bound = 0.0L;
  for(unsigned index = 0; index < coefficient.size(); ++index) {
    if(coefficient[index] == 0.0L) continue;
    for(unsigned derivative = 0; derivative < dimension; ++derivative) {
      long double term = coefficient[index];
      for(unsigned axis = 0; axis < dimension; ++axis) {
        unsigned degree = (index/IntegerPower(extent,dimension-axis-1U))%extent;
        term *= derivative == axis ? (degree == 0 ? 0.0L : degree*powers[axis][degree-1U]) : powers[axis][degree];
      }
      gradient[derivative] += term;
      bound += std::abs(term);
    }
  }
  long double magnitude = 0.0L;
  for(unsigned axis = 0; axis < dimension; ++axis) magnitude = std::hypot(magnitude,gradient[axis]);
  if(!std::isfinite(magnitude) || magnitude <= 64.0L*std::numeric_limits<long double>::epsilon()*bound || magnitude == 0.0L) throw std::domain_error("Direct Algoim requires a regular level-set gradient on the interface.");
  normal.resize(dimension);
  for(unsigned axis = 0; axis < dimension; ++axis) normal[axis] = static_cast<double>(gradient[axis]/magnitude);
}

template<int N>
std::vector<double> ReferencePoint(const algoim::uvector<double,N>& chart, unsigned simplexDimension) {
  std::vector<double> point(N);
  for(int axis = 0; axis < N; ++axis) point[axis] = 2.0*chart(axis)-1.0;
  if constexpr(N == 2) {
    if(simplexDimension == 2) point = {1.0-chart(0),chart(0)*chart(1)};
  }
  if constexpr(N == 3) {
    if(simplexDimension == 3) point = {chart(0)*(1.0-chart(1)),chart(0)*chart(1)*(1.0-chart(2)),chart(0)*chart(1)*chart(2)};
  }
  return point;
}

template<int N>
void VolumePoint(const algoim::uvector<double,N>& local, double weight, const algoim::xarray<double,N>& polynomial, const algoim::uvector<double,N>& low, double width, unsigned simplexDimension, std::vector<std::vector<double>>& positivePoints, std::vector<double>& positiveWeights, std::vector<std::vector<double>>& negativePoints, std::vector<double>& negativeWeights) {
  const algoim::uvector<double,N> point = low+width*local;
  const double value = algoim::bernstein::evalBernsteinPoly(polynomial,point);
  double mappedWeight = weight*std::pow(width,N);
  if(simplexDimension == 0) mappedWeight *= static_cast<double>(1U << N);
  if constexpr(N == 2) {
    if(simplexDimension == 2) mappedWeight *= point(0);
  }
  if constexpr(N == 3) {
    if(simplexDimension == 3) mappedWeight *= point(0)*point(0)*point(1);
  }
  if(!std::isfinite(mappedWeight) || !std::isfinite(value) || mappedWeight < 0.0) throw std::runtime_error("Algoim returned an invalid volume point or weight.");
  if(mappedWeight == 0.0) return;
  // Tanh-sinh endpoint samples retain their native weights and use the positive phase when G is exactly zero.
  if(value >= 0.0) {
    positivePoints.push_back(ReferencePoint<N>(point,simplexDimension));
    positiveWeights.push_back(mappedWeight);
  }
  else {
    negativePoints.push_back(ReferencePoint<N>(point,simplexDimension));
    negativeWeights.push_back(mappedWeight);
  }
}

void AppendSurface(const std::vector<double>& point, double weight, const std::vector<long double>& original, unsigned originalExtent, std::vector<std::vector<double>>& surfacePoints, std::vector<double>& surfaceWeights, std::vector<std::vector<double>>& surfaceNormals) {
  if(!std::isfinite(weight) || weight < 0.0) throw std::runtime_error("Algoim returned an invalid surface weight.");
  if(weight == 0.0) return;
  std::vector<double> normal;
  SurfaceNormal(point,original,originalExtent,normal);
  surfacePoints.push_back(point);
  surfaceWeights.push_back(weight);
  surfaceNormals.push_back(normal);
}

template<int N>
void SurfacePoint(const algoim::uvector<double,N>& local, double weight, const algoim::uvector<double,N>&, const algoim::xarray<double,N>& polynomial, const algoim::uvector<double,N>& low, double width, unsigned simplexDimension, const std::vector<unsigned>& removed, const std::vector<long double>& original, unsigned originalExtent, std::vector<std::vector<double>>& surfacePoints, std::vector<double>& surfaceWeights, std::vector<std::vector<double>>& surfaceNormals) {
  const algoim::uvector<double,N> point = low+width*local;
  for(int axis = 0; axis < N; ++axis) {
    if((removed[2U*axis] != 0 && point(axis) == 0.0) || (removed[2U*axis+1U] != 0 && point(axis) == 1.0)) return;
    if(local(axis) == 0.0 || local(axis) == 1.0) weight *= 0.5;
  }
  double metric = static_cast<double>(1U << (N-1));
  if(simplexDimension != 0) {
    algoim::uvector<double,N> normal = algoim::bernstein::evalBernsteinPolyGradient(polynomial,point);
    double magnitude = 0.0;
    for(int axis = 0; axis < N; ++axis) magnitude = std::hypot(magnitude,normal(axis));
    if(!std::isfinite(magnitude) || magnitude == 0.0) throw std::domain_error("Algoim produced a singular chart interface.");
    normal /= magnitude;
    if constexpr(N == 2) {
      metric = std::hypot(point(0)*normal(0)-point(1)*normal(1),normal(1));
    }
    if constexpr(N == 3) {
      const double common = point(0)*point(0)*point(1)*normal(0);
      const double first = common-point(0)*point(1)*point(1)*normal(1);
      const double second = common+point(0)*point(1)*(1.0-point(1))*normal(1)-point(0)*point(2)*normal(2);
      const double third = common+point(0)*point(1)*(1.0-point(1))*normal(1)+point(0)*(1.0-point(2))*normal(2);
      metric = std::hypot(std::hypot(first,second),third);
    }
  }
  AppendSurface(ReferencePoint<N>(point,simplexDimension),weight*std::pow(width,N-1)*metric,original,originalExtent,surfacePoints,surfaceWeights,surfaceNormals);
}

void BoundaryQuadrature(unsigned dimension, unsigned simplexDimension, const std::vector<unsigned>& physicalZero, const std::vector<long double>& original, unsigned originalExtent, unsigned quadratureOrder, unsigned subdivisions, std::vector<std::vector<double>>& surfacePoints, std::vector<double>& surfaceWeights, std::vector<std::vector<double>>& surfaceNormals) {
  const unsigned cells = subdivisions == 0 ? 1 : subdivisions;
  const unsigned count = quadratureOrder*cells;
  for(unsigned boundary = 0; boundary < physicalZero.size(); ++boundary) {
    if(physicalZero[boundary] == 0) continue;
    if(dimension == 1) {
      AppendSurface({boundary == 0 ? -1.0 : 1.0},0.5,original,originalExtent,surfacePoints,surfaceWeights,surfaceNormals);
      continue;
    }
    for(unsigned i = 0; i < count; ++i) {
      const double a = (static_cast<double>(i/quadratureOrder)+algoim::GaussQuad::x(quadratureOrder,i%quadratureOrder))/cells;
      const double wa = algoim::GaussQuad::w(quadratureOrder,i%quadratureOrder)/cells;
      for(unsigned j = 0; j < (dimension == 3 ? count : 1U); ++j) {
        const double b = (static_cast<double>(j/quadratureOrder)+algoim::GaussQuad::x(quadratureOrder,j%quadratureOrder))/cells;
        const double wb = dimension == 3 ? algoim::GaussQuad::w(quadratureOrder,j%quadratureOrder)/cells : 1.0;
        std::vector<double> point(dimension,0.0);
        double metric = static_cast<double>(1U << (dimension-1U));
        if(simplexDimension == 2) {
          point = boundary == 0 ? std::vector<double>{0.0,a} : boundary == 1 ? std::vector<double>{a,0.0} : std::vector<double>{1.0-a,a};
          metric = boundary == 2 ? std::sqrt(2.0) : 1.0;
        }
        else if(simplexDimension == 3) {
          if(boundary == 0) point = {a,(1.0-a)*b,(1.0-a)*(1.0-b)};
          else if(boundary == 1) point = {0.0,a,(1.0-a)*b};
          else if(boundary == 2) point = {a,0.0,(1.0-a)*b};
          else point = {a,(1.0-a)*b,0.0};
          metric = (1.0-a)*(boundary == 0 ? std::sqrt(3.0) : 1.0);
        }
        else {
          unsigned tangential = 0;
          for(unsigned axis = 0; axis < dimension; ++axis) {
            if(axis == boundary/2U) point[axis] = boundary%2U == 0 ? -1.0 : 1.0;
            else point[axis] = 2.0*(tangential++ == 0 ? a : b)-1.0;
          }
        }
        AppendSurface(point,0.5*metric*wa*wb,original,originalExtent,surfacePoints,surfaceWeights,surfaceNormals);
      }
    }
  }
}

template<int N>
void BernsteinPolynomial(const std::vector<long double>& coefficient, unsigned powerExtent, algoim::uvector<int,N>& extent, std::vector<double>& data) {
  extent = algoim::uvector<int,N>(2);
  for(unsigned index = 0; index < coefficient.size(); ++index) {
    if(coefficient[index] == 0.0L) continue;
    for(int axis = 0; axis < N; ++axis) extent(axis) = std::max(extent(axis),static_cast<int>((index/IntegerPower(powerExtent,N-axis-1))%powerExtent)+1);
  }
  unsigned size = 1;
  for(int axis = 0; axis < N; ++axis) size *= extent(axis);
  std::vector<long double> converted(size,0.0L);
  for(algoim::MultiLoop<N> index(0,extent); ~index; ++index) {
    unsigned source = 0;
    unsigned target = 0;
    for(int axis = 0; axis < N; ++axis) {
      source = source*powerExtent+index(axis);
      target = target*extent(axis)+index(axis);
    }
    converted[target] = coefficient[source];
  }
  for(int axis = 0; axis < N; ++axis) {
    unsigned stride = 1;
    for(int other = axis+1; other < N; ++other) stride *= extent(other);
    std::vector<long double> next(size,0.0L);
    for(unsigned index = 0; index < size; ++index) {
      unsigned power = (index/stride)%extent(axis);
      unsigned offset = index-power*stride;
      for(unsigned k = 0; k <= power; ++k) next[index] += converted[offset+k*stride]*Binomial(power,k)/Binomial(extent(axis)-1,k);
    }
    converted.swap(next);
  }
  long double scale = 0.0L;
  for(unsigned index = 0; index < size; ++index) scale = std::max(scale,std::abs(converted[index]));
  if(scale == 0.0L || !std::isfinite(scale)) throw std::domain_error("The level-set polynomial is identically zero or invalid.");
  data.resize(size);
  for(unsigned index = 0; index < size; ++index) data[index] = static_cast<double>(converted[index]/scale);
}

template<int N>
bool MonotoneDirection(const algoim::xarray<double,N>& polynomial) {
  std::vector<double> data(polynomial.size());
  algoim::xarray<double,N> derivative(data.data(),polynomial.ext());
  for(int axis = 0; axis < N; ++axis) {
    algoim::bernstein::elevatedDerivative(polynomial,axis,derivative);
    if(algoim::bernstein::uniformSign(derivative) != 0) return true;
  }
  return false;
}

template<int N>
void TensorQuadrature(const std::vector<long double>& coefficient, unsigned powerExtent, const std::vector<long double>& original, unsigned originalExtent, unsigned simplexDimension, const std::vector<unsigned>& removed, const std::vector<unsigned>& physicalZero, unsigned quadratureOrder, unsigned subdivisions, std::vector<std::vector<double>>& positivePoints, std::vector<double>& positiveWeights, std::vector<std::vector<double>>& negativePoints, std::vector<double>& negativeWeights, std::vector<std::vector<double>>& surfacePoints, std::vector<double>& surfaceWeights, std::vector<std::vector<double>>& surfaceNormals) {
  algoim::uvector<int,N> extent;
  std::vector<double> data;
  BernsteinPolynomial<N>(coefficient,powerExtent,extent,data);
  algoim::xarray<double,N> polynomial(data.data(),extent);
  const bool monotone = MonotoneDirection<N>(polynomial);
  if(subdivisions == 0) subdivisions = algoim::bernstein::uniformSign(polynomial) != 0 || monotone ? 1 : 4;
  const bool fullChartSurface = N > 1 && monotone && subdivisions > 1;
  const double width = 1.0/subdivisions;
  std::vector<double> localData(data.size());
  algoim::xarray<double,N> localPolynomial(localData.data(),extent);
  for(algoim::MultiLoop<N> cell(0,subdivisions); ~cell; ++cell) {
    const algoim::uvector<double,N> low = cell()*width;
    const algoim::uvector<double,N> high = low+width;
    algoim::bernstein::deCasteljau(polynomial,low,high,localPolynomial);
    algoim::ImplicitPolyQuadrature<N> quadrature(localPolynomial);
    quadrature.integrate(algoim::AutoMixed,quadratureOrder,std::bind(&VolumePoint<N>,std::placeholders::_1,std::placeholders::_2,std::cref(polynomial),std::cref(low),width,simplexDimension,std::ref(positivePoints),std::ref(positiveWeights),std::ref(negativePoints),std::ref(negativeWeights)));
    if constexpr(N > 1) {
      if(!fullChartSurface) quadrature.integrate_surf(algoim::AutoMixed,quadratureOrder,std::bind(&SurfacePoint<N>,std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::cref(polynomial),std::cref(low),width,simplexDimension,std::cref(removed),std::cref(original),originalExtent,std::ref(surfacePoints),std::ref(surfaceWeights),std::ref(surfaceNormals)));
    }
    else {
      std::vector<double> roots(extent(0)-1);
      const int count = algoim::bernstein::bernsteinUnitIntervalRealRoots(localData.data(),extent(0),roots.data());
      for(int index = 0; index < count; ++index) SurfacePoint<1>(algoim::uvector<double,1>(roots[index]),1.0,algoim::uvector<double,1>(0.0),polynomial,low,width,simplexDimension,removed,original,originalExtent,surfacePoints,surfaceWeights,surfaceNormals);
    }
  }
  if constexpr(N > 1) {
    if(fullChartSurface) {
      // A monotone interface uses one native surface rule to avoid duplicated roots on artificial subdivision faces.
      const algoim::uvector<double,N> low(0.0);
      algoim::ImplicitPolyQuadrature<N> quadrature(polynomial);
      quadrature.integrate_surf(algoim::AutoMixed,quadratureOrder,std::bind(&SurfacePoint<N>,std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::cref(polynomial),std::cref(low),1.0,simplexDimension,std::cref(removed),std::cref(original),originalExtent,std::ref(surfacePoints),std::ref(surfaceWeights),std::ref(surfaceNormals)));
    }
  }
  BoundaryQuadrature(N,simplexDimension,physicalZero,original,originalExtent,quadratureOrder,subdivisions,surfacePoints,surfaceWeights,surfaceNormals);
}

}

void BuildQuadrature(const std::vector<double>& nodalValues, std::vector<std::vector<double>>& positivePoints, std::vector<double>& positiveWeights, std::vector<std::vector<double>>& negativePoints, std::vector<double>& negativeWeights, std::vector<std::vector<double>>& surfacePoints, std::vector<double>& surfaceWeights, std::vector<std::vector<double>>& surfaceNormals, unsigned quadratureOrder, unsigned subdivisions) {
  const unsigned nodes = nodalValues.size();
  if(nodes != 3 && nodes != 7 && nodes != 9 && nodes != 15 && nodes != 27) throw std::invalid_argument("Direct Algoim supports LINE3, TRI7, QUAD9, TET15, and HEX27.");
  if(quadratureOrder == 0 || quadratureOrder > 100 || subdivisions > 32) throw std::invalid_argument("Direct Algoim requires order 1 to 100 and subdivisions 0 to 32.");
  if(&positivePoints == &negativePoints || &positivePoints == &surfacePoints || &positivePoints == &surfaceNormals || &negativePoints == &surfacePoints || &negativePoints == &surfaceNormals || &surfacePoints == &surfaceNormals || &positiveWeights == &negativeWeights || &positiveWeights == &surfaceWeights || &negativeWeights == &surfaceWeights || &nodalValues == &positiveWeights || &nodalValues == &negativeWeights || &nodalValues == &surfaceWeights) throw std::invalid_argument("Direct Algoim input and output vectors must be distinct.");
  long double scale = 0.0L;
  for(unsigned node = 0; node < nodes; ++node) {
    if(!std::isfinite(nodalValues[node])) throw std::invalid_argument("A level-set nodal value is not finite.");
    scale = std::max(scale,std::abs(static_cast<long double>(nodalValues[node])));
  }
  if(scale == 0.0L) throw std::domain_error("An identically zero level set has no regular interface quadrature.");
  const unsigned dimension = nodes == 3 ? 1 : nodes == 7 || nodes == 9 ? 2 : 3;
  const unsigned simplexDimension = nodes == 7 ? 2 : nodes == 15 ? 3 : 0;
  const unsigned originalExtent = nodes == 15 ? 5 : 3;
  const unsigned powerExtent = nodes == 7 ? 4 : originalExtent;
  std::vector<long double> normalized(nodes);
  for(unsigned node = 0; node < nodes; ++node) normalized[node] = static_cast<long double>(nodalValues[node])/scale;
  std::vector<long double> original;
  std::vector<long double> coefficient;
  if(nodes == 7) TriPolynomial(normalized,original,coefficient);
  else if(nodes == 15) TetPolynomial(normalized,original,coefficient);
  else {
    TensorPolynomial(normalized,dimension,false,original);
    TensorPolynomial(normalized,dimension,true,coefficient);
  }
  std::vector<unsigned> removed;
  std::vector<unsigned> physicalZero;
  BoundaryFlags(normalized,dimension,simplexDimension,removed,physicalZero);
  RemoveBoundaryFactors(coefficient,dimension,powerExtent,removed);
  positivePoints.clear();
  positiveWeights.clear();
  negativePoints.clear();
  negativeWeights.clear();
  surfacePoints.clear();
  surfaceWeights.clear();
  surfaceNormals.clear();
  if(dimension == 1) TensorQuadrature<1>(coefficient,powerExtent,original,originalExtent,simplexDimension,removed,physicalZero,quadratureOrder,subdivisions,positivePoints,positiveWeights,negativePoints,negativeWeights,surfacePoints,surfaceWeights,surfaceNormals);
  else if(dimension == 2) TensorQuadrature<2>(coefficient,powerExtent,original,originalExtent,simplexDimension,removed,physicalZero,quadratureOrder,subdivisions,positivePoints,positiveWeights,negativePoints,negativeWeights,surfacePoints,surfaceWeights,surfaceNormals);
  else TensorQuadrature<3>(coefficient,powerExtent,original,originalExtent,simplexDimension,removed,physicalZero,quadratureOrder,subdivisions,positivePoints,positiveWeights,negativePoints,negativeWeights,surfacePoints,surfaceWeights,surfaceNormals);
}

}
