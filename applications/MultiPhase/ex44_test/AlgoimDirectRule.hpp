#pragma once

// Minimal direct-point adapter for Algoim on LINE3/QUAD9/TRI7/HEX27/TET15.
// The returned coordinates are FEMuS reference coordinates, NOT unit-cube coordinates.
// volume.weight integrates d(xi); surfaceDelta.weight integrates delta(Psi_h) d(xi).
// To obtain physical surface measure, multiply surfaceDelta.weight by
// |det(Dx/Dxi)| * |grad_x Psi_h| at each returned point.
// This header deliberately does not reconstruct FEMuS moment-fitting weights.

#include "FELevelSetMoments.hpp"
#include "Hex27LevelSetMoments.hpp"
#include "Tet15LevelSet.hpp"
#include "algoim/quadrature_multipoly.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

namespace FemusAlgoimDirect {

struct Point {
  std::vector<double> xi;
  double weight;
};

struct Rule {
  std::vector<Point> volume;
  std::vector<Point> surfaceDelta;
};

namespace detail {

inline long double Binomial(unsigned n, unsigned k) {
  long double value = 1.0L;
  for(unsigned i = 1; i <= k; ++i) value *= static_cast<long double>(n+1U-i)/i;
  return value;
}

inline std::vector<double> ToReference(const algoim::uvector<double,2>& canonical, unsigned simplexDimension) {
  if(simplexDimension == 2) return {1.0-canonical(0),canonical(0)*canonical(1)};
  if(simplexDimension != 0) throw std::invalid_argument("AlgoimDirect: invalid 2D simplex mapping");
  return {2.0*canonical(0)-1.0,2.0*canonical(1)-1.0};
}

inline std::vector<double> ToReference(const algoim::uvector<double,3>& canonical, unsigned simplexDimension) {
  if(simplexDimension == 3) return {canonical(0)*(1.0-canonical(1)),canonical(0)*canonical(1)*(1.0-canonical(2)),canonical(0)*canonical(1)*canonical(2)};
  if(simplexDimension != 0) throw std::invalid_argument("AlgoimDirect: invalid 3D simplex mapping");
  return {2.0*canonical(0)-1.0,2.0*canonical(1)-1.0,2.0*canonical(2)-1.0};
}

template<int N>
long double ReferenceJacobianFactor(const algoim::uvector<double,N>& canonical, unsigned simplexDimension) {
  if(simplexDimension == 2) return canonical(0);
  if(simplexDimension == 3) return static_cast<long double>(canonical(0))*canonical(0)*canonical(1);
  return static_cast<long double>(1U << N);
}

template<int N>
void AddFullBox(const algoim::uvector<double,N>& low, const algoim::uvector<double,N>& high, unsigned order, unsigned simplexDimension, Rule& rule) {
  for(algoim::MultiLoop<N> node(0,static_cast<int>(order)); ~node; ++node) {
    algoim::uvector<double,N> canonical;
    long double weight = 1.0L;
    for(int d = 0; d < N; ++d) {
      const double width = high(d)-low(d);
      canonical(d) = low(d)+width*algoim::GaussQuad::x(static_cast<int>(order),node(d));
      weight *= static_cast<long double>(width)*algoim::GaussQuad::w(static_cast<int>(order),node(d));
    }
    weight *= ReferenceJacobianFactor<N>(canonical,simplexDimension);
    rule.volume.push_back(Point{ToReference(canonical,simplexDimension),static_cast<double>(weight)});
  }
}

template<int N>
bool HasMonotoneDirection(const algoim::xarray<double,N>& polynomial) {
  for(int axis = 0; axis < N; ++axis) {
    unsigned stride = 1;
    for(int other = axis+1; other < N; ++other) stride *= static_cast<unsigned>(polynomial.ext(other));
    bool positive = true;
    bool negative = true;
    bool examined = false;
    for(int index = 0; index < polynomial.size(); ++index) {
      const unsigned coordinate = (static_cast<unsigned>(index)/stride)%static_cast<unsigned>(polynomial.ext(axis));
      if(coordinate+1U < static_cast<unsigned>(polynomial.ext(axis))) {
        const double derivative = (polynomial[index+stride]-polynomial[index])*(polynomial.ext(axis)-1);
        positive = positive && derivative > 256.0*std::numeric_limits<double>::epsilon();
        negative = negative && derivative < -256.0*std::numeric_limits<double>::epsilon();
        examined = true;
      }
    }
    if(examined && (positive || negative)) return true;
  }
  return false;
}

template<int N>
void BuildTensor(const std::vector<long double>& powerCoefficient, algoim::uvector<int,N> extent, unsigned simplexDimension, unsigned removedFirstPower, unsigned removedSecondPower, unsigned order, unsigned subdivisions, Rule& rule) {
  unsigned size = 1;
  for(int axis = 0; axis < N; ++axis) size *= static_cast<unsigned>(extent(axis));
  if(powerCoefficient.size() != size) throw std::invalid_argument("AlgoimDirect: inconsistent coefficient extent");
  long double scale = 0.0L;
  for(unsigned i = 0; i < size; ++i) {
    if(!std::isfinite(powerCoefficient[i])) throw std::invalid_argument("AlgoimDirect: nonfinite polynomial coefficient");
    scale = std::max(scale,std::abs(powerCoefficient[i]));
  }
  if(scale == 0.0L) throw std::domain_error("AlgoimDirect: identically zero polynomial");
  std::vector<long double> converted(size);
  for(unsigned i = 0; i < size; ++i) converted[i] = powerCoefficient[i]/scale;
  unsigned stride = size;
  for(int axis = 0; axis < N; ++axis) {
    const unsigned count = static_cast<unsigned>(extent(axis));
    stride /= count;
    std::vector<long double> next(size,0.0L);
    for(unsigned i = 0; i < size; ++i) {
      const unsigned power = (i/stride)%count;
      const unsigned offset = i-power*stride;
      for(unsigned k = 0; k <= power; ++k) next[i] += converted[offset+k*stride]*Binomial(power,k)/Binomial(count-1U,k);
    }
    converted.swap(next);
  }
  std::vector<double> bernstein(size);
  for(unsigned i = 0; i < size; ++i) bernstein[i] = static_cast<double>(converted[i]);
  algoim::xarray<double,N> polynomial(bernstein.data(),extent);
  const int globalPhase = algoim::bernstein::uniformSign(polynomial);
  if(subdivisions == 0) subdivisions = (globalPhase != 0 || HasMonotoneDirection<N>(polynomial)) ? 1U : 4U;
  const double width = 1.0/static_cast<double>(subdivisions);
  if(globalPhase > 0) AddFullBox<N>(algoim::uvector<double,N>(0.0),algoim::uvector<double,N>(1.0),order,simplexDimension,rule);
  if(globalPhase != 0) return;

  for(algoim::MultiLoop<N> cell(0,static_cast<int>(subdivisions)); ~cell; ++cell) {
    algoim::uvector<double,N> low;
    algoim::uvector<double,N> high;
    for(int axis = 0; axis < N; ++axis) {
      low(axis) = static_cast<double>(cell(axis))/subdivisions;
      high(axis) = static_cast<double>(cell(axis)+1)/subdivisions;
    }
    std::vector<double> localData(size);
    algoim::xarray<double,N> localPolynomial(localData.data(),extent);
    algoim::bernstein::deCasteljau(polynomial,low,high,localPolynomial);
    const int phase = algoim::bernstein::uniformSign(localPolynomial);
    if(phase > 0) {
      AddFullBox<N>(low,high,order,simplexDimension,rule);
      continue;
    }
    if(phase < 0) continue;

    algoim::ImplicitPolyQuadrature<N> quadrature(localPolynomial);
    quadrature.integrate(algoim::AutoMixed,static_cast<int>(order),[&](const algoim::uvector<double,N>& local, double weight) {
      const algoim::uvector<double,N> point = low+width*local;
      if(algoim::bernstein::evalBernsteinPoly(polynomial,point) <= 0.0) return;
      long double mappedWeight = weight;
      for(int axis = 0; axis < N; ++axis) mappedWeight *= width;
      mappedWeight *= ReferenceJacobianFactor<N>(point,simplexDimension);
      rule.volume.push_back(Point{ToReference(point,simplexDimension),static_cast<double>(mappedWeight)});
    });
    quadrature.integrate_surf(algoim::AutoMixed,static_cast<int>(order),[&](const algoim::uvector<double,N>& local, double weight, const algoim::uvector<double,N>& normalWeight) {
      static_cast<void>(normalWeight);
      const algoim::uvector<double,N> point = low+width*local;
      const algoim::uvector<double,N> grad = algoim::bernstein::evalBernsteinPolyGradient(polynomial,point);
      long double norm = 0.0L;
      for(int axis = 0; axis < N; ++axis) norm = std::hypot(norm,static_cast<long double>(grad(axis)));
      if(!std::isfinite(norm) || !(norm > 64.0L*std::numeric_limits<double>::epsilon())) throw std::domain_error("AlgoimDirect: nonregular level-set gradient");
      long double mappedWeight = static_cast<long double>(weight)/norm/scale;
      for(int axis = 1; axis < N; ++axis) mappedWeight *= width;
      for(int axis = 0; axis < N; ++axis) if(local(axis) == 0.0 || local(axis) == 1.0) mappedWeight *= 0.5L;
      if(simplexDimension == 2 && removedFirstPower == 0) mappedWeight *= point(0);
      if(simplexDimension == 3) {
        for(unsigned power = removedFirstPower; power < 2U; ++power) mappedWeight *= point(0);
        if(removedSecondPower == 0) mappedWeight *= point(1);
      }
      if(simplexDimension == 0) mappedWeight *= static_cast<long double>(1U << N);
      // This is the delta measure in FEMuS reference coordinates. The caller
      // multiplies it by |det J| * |J^{-T} grad_ref Psi|, point by point.
      rule.surfaceDelta.push_back(Point{ToReference(point,simplexDimension),static_cast<double>(mappedWeight)});
    });
  }
}

} // namespace detail

inline Rule GetNodalQuadrature(const std::vector<double>& nodalPsi, unsigned order = 20, unsigned subdivisions = 0) {
  if(order == 0 || order > 100 || subdivisions > 32) throw std::invalid_argument("AlgoimDirect: invalid order/subdivisions");
  Rule rule;
  if(nodalPsi.size() == 9) {
    std::vector<std::vector<double>> coefficient;
    FELevelSetMoments::Quad9ReferenceCoefficientsFromNodalValues(nodalPsi,coefficient);
    std::vector<std::vector<long double>> chart;
    FELevelSetMoments::ReferencePolynomialAffine2D(coefficient,-1.0L,2.0L,-1.0L,2.0L,chart);
    std::vector<long double> data(9);
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) data[3U*i+j] = chart[i][j];
    detail::BuildTensor<2>(data,algoim::uvector<int,2>(3,3),0,0,0,order,subdivisions,rule);
  }
  else if(nodalPsi.size() == 7) {
    // Exact zero edges require half-boundary contributions. Handle them explicitly
    // instead of silently producing incorrect surface integrals.
    const bool zero01 = nodalPsi[0] == 0.0 && nodalPsi[1] == 0.0 && nodalPsi[3] == 0.0;
    const bool zero02 = nodalPsi[0] == 0.0 && nodalPsi[2] == 0.0 && nodalPsi[5] == 0.0;
    const bool zero12 = nodalPsi[1] == 0.0 && nodalPsi[2] == 0.0 && nodalPsi[4] == 0.0;
    if(zero01 || zero02 || zero12) throw std::domain_error("AlgoimDirect: exact zero TRI7 edges need boundary quadrature");
    std::vector<std::vector<double>> coefficient;
    FELevelSetMoments::Tri7ReferenceCoefficientsFromNodalValues(nodalPsi,coefficient);
    std::vector<std::vector<long double>> chart;
    FELevelSetMoments::ReferencePolynomialAffine2D(coefficient,1.0L,-1.0L,0.0L,1.0L,chart);
    std::vector<long double> data(12,0.0L);
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) if(i+j < 4) data[3U*(i+j)+j] += chart[i][j];
    unsigned removed = 0;
    unsigned rows = 4;
    while(rows > 1 && data[0] == 0.0L && data[1] == 0.0L && data[2] == 0.0L) {
      data.erase(data.begin(),data.begin()+3);
      --rows;
      ++removed;
    }
    if(removed > 1) throw std::domain_error("AlgoimDirect: critical TRI7 collapsed vertex");
    detail::BuildTensor<2>(data,algoim::uvector<int,2>(static_cast<int>(rows),3),2,removed,0,order,subdivisions,rule);
  }
  else if(nodalPsi.size() == 27) {
    std::vector<std::vector<std::vector<double>>> coefficient;
    FELevelSetMoments::Hex27ReferenceCoefficientsFromNodalValues(nodalPsi,coefficient);
    std::vector<long double> polynomial(27), data;
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) for(unsigned k = 0; k < 3; ++k) polynomial[9U*i+3U*j+k] = coefficient[i][j][k];
    FELevelSetMoments::ReferenceCubePolynomialToChart(polynomial,data);
    detail::BuildTensor<3>(data,algoim::uvector<int,3>(3,3,3),0,0,0,order,subdivisions,rule);
  }
  else if(nodalPsi.size() == 15) {
    const unsigned faces[4][7] = {{1,2,3,5,9,8,12},{0,2,3,6,9,7,13},{0,1,3,4,8,7,11},{0,1,2,4,5,6,10}};
    for(unsigned face = 0; face < 4; ++face) {
      bool zero = true;
      for(unsigned n = 0; n < 7; ++n) zero = zero && nodalPsi[faces[face][n]] == 0.0;
      if(zero) throw std::domain_error("AlgoimDirect: exact zero TET15 faces need boundary quadrature");
    }
    std::vector<long double> polynomial;
    FELevelSetMoments::Tet15PolynomialFromNodalValues(nodalPsi,true,polynomial);
    unsigned removed[2] = {0,0};
    for(unsigned axis = 0; axis < 2; ++axis) {
      const unsigned stride = (axis == 0) ? 25U : 5U;
      bool factor = true;
      while(factor && removed[axis] < 4U) {
        for(unsigned index = 0; index < 125; ++index) if((index/stride)%5U == 0 && polynomial[index] != 0.0L) factor = false;
        if(factor) {
          std::vector<long double> next(125,0.0L);
          for(unsigned index = 0; index < 125; ++index) if((index/stride)%5U != 0) next[index-stride] = polynomial[index];
          polynomial.swap(next);
          ++removed[axis];
        }
      }
    }
    if(removed[0] > 1 || removed[1] > 1) throw std::domain_error("AlgoimDirect: critical TET15 collapsed vertex");
    unsigned maximum[3] = {1,1,1};
    for(unsigned i = 0; i < 5; ++i) for(unsigned j = 0; j < 5; ++j) for(unsigned k = 0; k < 5; ++k) if(polynomial[FELevelSetMoments::Tet15Index(i,j,k)] != 0.0L) {
      maximum[0] = std::max(maximum[0],i);
      maximum[1] = std::max(maximum[1],j);
      maximum[2] = std::max(maximum[2],k);
    }
    std::vector<long double> data((maximum[0]+1U)*(maximum[1]+1U)*(maximum[2]+1U));
    unsigned index = 0;
    for(unsigned i = 0; i <= maximum[0]; ++i) for(unsigned j = 0; j <= maximum[1]; ++j) for(unsigned k = 0; k <= maximum[2]; ++k) data[index++] = polynomial[FELevelSetMoments::Tet15Index(i,j,k)];
    detail::BuildTensor<3>(data,algoim::uvector<int,3>(static_cast<int>(maximum[0]+1U),static_cast<int>(maximum[1]+1U),static_cast<int>(maximum[2]+1U)),3,removed[0],removed[1],order,subdivisions,rule);
  }
  else throw std::invalid_argument("AlgoimDirect: requires QUAD9, TRI7, HEX27, or TET15 nodal Psi");
  return rule;
}

} // namespace FemusAlgoimDirect
