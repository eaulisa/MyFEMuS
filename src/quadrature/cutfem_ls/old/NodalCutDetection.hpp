#ifndef FEMUS_NODAL_CUT_DETECTION_HPP
#define FEMUS_NODAL_CUT_DETECTION_HPP

#include "GeomElTypeEnum.hpp"
#include "UnitElementMoments.hpp"
#include "Tet15LevelSet.hpp"
#include "Wedge21LevelSetMoments.hpp"

namespace FELevelSetMoments {

// Converts the represented nodal polynomial exactly to Bernstein form, with a roundoff allowance.
inline void NodalLevelSetBernsteinControls(GeomElType geometry, const std::vector<double>& nodal, std::vector<long double>& controls, long double& padding) {
  std::vector<long double> coefficient;
  long double rounding = std::numeric_limits<long double>::epsilon();
  if(geometry == LINE || geometry == QUAD || geometry == HEX) {
    UnitTensorPolynomialFromNodalValues(nodal,geometry == LINE ? 1U : (geometry == QUAD ? 2U : 3U),coefficient);
    controls = Hex27Bernstein(coefficient);
  }
  else {
    const unsigned degree = geometry == TET ? 4U : 3U;
    const unsigned lineDegree = geometry == WEDGE ? 2U : 0U;
    std::vector<long double> polynomial(125,0);
    if(geometry == TRI) {
      std::vector<std::vector<double>> planar;
      Tri7ReferenceCoefficientsFromNodalValues(nodal,planar);
      for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) polynomial[Tet15Index(i,j,0)] = planar[i][j];
      rounding = std::numeric_limits<double>::epsilon();
    }
    else if(geometry == TET) Tet15PolynomialFromNodalValues(nodal,false,polynomial);
    else if(geometry == WEDGE) {
      Wedge21UnitPolynomialFromNodalValues(nodal,coefficient);
      for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) for(unsigned k = 0; k < 3; ++k) polynomial[Tet15Index(i,j,k)] = coefficient[Hex27Index(i,j,k)];
      rounding = std::numeric_limits<double>::epsilon();
    }
    else throw std::invalid_argument("Unsupported nodal cut-detection geometry.");
    coefficient = polynomial;
    controls.clear();
    for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) for(unsigned k = 0; k <= (geometry == TET ? degree-i-j : lineDegree); ++k) {
      long double value = 0;
      for(unsigned a = 0; a <= i; ++a) for(unsigned b = 0; b <= j; ++b) for(unsigned c = 0; c <= k; ++c) {
        long double factor = ReferenceBinomial(i,a)*ReferenceBinomial(j,b);
        if(geometry == TET) factor *= ReferenceBinomial(k,c)/(ReferenceBinomial(degree,a+b+c)*ReferenceBinomial(a+b+c,a)*ReferenceBinomial(b+c,b));
        else factor *= ReferenceBinomial(k,c)/(ReferenceBinomial(degree,a+b)*ReferenceBinomial(a+b,a)*ReferenceBinomial(lineDegree,c));
        value += factor*polynomial[Tet15Index(a,b,c)];
      }
      controls.push_back(value);
    }
  }
  long double norm = 0;
  for(unsigned n = 0; n < coefficient.size(); ++n) norm += std::abs(coefficient[n]);
  padding = 8192*rounding*norm;
}

inline bool NodalCutPositiveTensor(const std::vector<long double>& controls, unsigned degree, unsigned dimension, long double padding, unsigned depth) {
  if(*std::min_element(controls.begin(),controls.end()) > padding) return true;
  if(depth == 0 || degree == 0) return false;
  const unsigned side = degree+1U;
  for(unsigned i = 0; i <= degree; i += degree) for(unsigned j = 0; j <= degree; j += degree) for(unsigned k = 0; k <= degree; k += degree) if(controls[side*side*i+side*j+k] <= padding) return false;
  unsigned axis = 0;
  long double variation = 0;
  for(unsigned d = 0; d < dimension; ++d) for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree; ++j) for(unsigned k = 0; k <= degree; ++k) {
    const unsigned index[3] = {i,j,k};
    const unsigned stride = d == 0 ? side*side : (d == 1 ? side : 1U);
    if(index[d] < degree && std::abs(controls[side*side*i+side*j+k+stride]-controls[side*side*i+side*j+k]) > variation) { variation = std::abs(controls[side*side*i+side*j+k+stride]-controls[side*side*i+side*j+k]); axis = d; }
  }
  if(variation == 0) return false;
  std::vector<long double> left(controls.size()), right(controls.size());
  const unsigned stride = axis == 0 ? side*side : (axis == 1 ? side : 1U);
  for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree; ++j) for(unsigned k = 0; k <= degree; ++k) {
    const unsigned index[3] = {i,j,k};
    if(index[axis] != 0) continue;
    const unsigned base = side*side*i+side*j+k;
    std::vector<long double> row(side);
    for(unsigned n = 0; n <= degree; ++n) row[n] = controls[base+n*stride];
    left[base] = row[0];
    right[base+degree*stride] = row[degree];
    for(unsigned level = 1; level <= degree; ++level) {
      for(unsigned n = 0; n <= degree-level; ++n) row[n] = 0.5L*row[n]+0.5L*row[n+1U];
      left[base+level*stride] = row[0];
      right[base+(degree-level)*stride] = row[degree-level];
    }
  }
  return NodalCutPositiveTensor(left,degree,dimension,padding,depth-1U) && NodalCutPositiveTensor(right,degree,dimension,padding,depth-1U);
}

inline void NodalCutTensorControls(GeomElType geometry, const std::vector<double>& nodal, std::vector<long double>& controls, unsigned& degree) {
  if(geometry == LINE || geometry == QUAD || geometry == HEX) {
    std::vector<long double> coefficient;
    UnitTensorPolynomialFromNodalValues(nodal,geometry == LINE ? 1U : (geometry == QUAD ? 2U : 3U),coefficient);
    controls = Hex27Bernstein(coefficient);
    degree = 2;
    return;
  }
  degree = geometry == TET ? 4U : 3U;
  const unsigned side = degree+1U;
  std::vector<long double> polynomial(side*side*side,0);
  if(geometry == TET) Tet15PolynomialFromNodalValues(nodal,true,polynomial);
  else {
    std::vector<long double> source(27,0);
    if(geometry == WEDGE) Wedge21UnitPolynomialFromNodalValues(nodal,source);
    else {
      std::vector<std::vector<double>> planar;
      Tri7ReferenceCoefficientsFromNodalValues(nodal,planar);
      for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) source[Hex27Index(i,j,0)] = planar[i][j];
    }
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) for(unsigned k = 0; k < 3; ++k) for(unsigned n = 0; n <= j; ++n) if(source[Hex27Index(i,j,k)] != 0) polynomial[side*side*(i+n)+side*j+k] += (n%2U == 0 ? 1 : -1)*ReferenceBinomial(j,n)*source[Hex27Index(i,j,k)];
  }
  controls.assign(polynomial.size(),0);
  for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree; ++j) for(unsigned k = 0; k <= degree; ++k) for(unsigned a = 0; a <= i; ++a) for(unsigned b = 0; b <= j; ++b) for(unsigned c = 0; c <= k; ++c) controls[side*side*i+side*j+k] += polynomial[side*side*a+side*b+c]*ReferenceBinomial(i,a)*ReferenceBinomial(j,b)*ReferenceBinomial(k,c)/(ReferenceBinomial(degree,a)*ReferenceBinomial(degree,b)*ReferenceBinomial(degree,c));
}

// Subdivision removes false alarms; unresolved cells and contacts remain possible cuts.
inline int NodalLevelSetPhase(GeomElType geometry, const std::vector<double>& nodal) {
  std::vector<long double> controls;
  long double padding;
  NodalLevelSetBernsteinControls(geometry,nodal,controls,padding);
  const long double minimum = *std::min_element(controls.begin(),controls.end());
  const long double maximum = *std::max_element(controls.begin(),controls.end());
  if(minimum > padding) return 1;
  if(maximum < -padding) return -1;
  const int phase = nodal[0] > 0 ? 1 : -1;
  for(unsigned n = 0; n < nodal.size(); ++n) if(nodal[n] == 0 || (nodal[n] > 0 ? 1 : -1) != phase) return 0;
  unsigned degree;
  NodalCutTensorControls(geometry,nodal,controls,degree);
  for(unsigned n = 0; n < controls.size(); ++n) controls[n] *= phase;
  const unsigned dimension = geometry == LINE ? 1U : (geometry == QUAD || geometry == TRI ? 2U : 3U);
  if(NodalCutPositiveTensor(controls,degree,dimension,64*padding,12)) return phase;
  return 0;
}

}

#endif
