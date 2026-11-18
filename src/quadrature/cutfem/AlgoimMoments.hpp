#ifndef FEMUS_ALGOIM_MOMENTS_HPP
#define FEMUS_ALGOIM_MOMENTS_HPP

// JCP 2022 Algoim moment backend for FEMUS QUAD9 and HEX27.
//
// Only moments are returned; FEMUS keeps its existing moment-to-weight map.
// No additional subdivisions, phase-rule construction, or point filtering.
// Algoim sees only the Bernstein polynomial and the 1D quadrature order.
// All quadrature points are in the unit box [0,1]^Dim.
//
// Requires the existing FEMUS polynomial/geometry routines, Algoim 2022,
// C++17, and LAPACKE (as provided by FemusAlgoim.cmake).

#include "CutFemPreparedMoments.hpp"
#include "algoim/quadrature_multipoly.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace FELevelSetMoments {

enum class AlgoimSurfaceMode {
  Delta,             // delta(G) du
  ReferenceSurface,  // dS_u, normalized for FEMUS reference Gauss weights
  PhysicalSurface,   // physical surface metric, normalized for FEMUS Gauss weights
  InverseSurface     // delta(G)/|cof(J_u) grad_u G| du, normalized for FEMUS
};

struct AlgoimMomentResult {
  std::vector<long double> positive;
  std::vector<long double> surface;

  // Number of callback contributions emitted by Algoim, including
  // volume contributions discarded because they belong to G <= 0.
  std::uint64_t volumePointsGenerated = 0;
  std::uint64_t volumePointsUsed = 0;
  std::uint64_t surfacePoints = 0;

  std::uint64_t TotalPoints() const {
    return volumePointsGenerated + surfacePoints;
  }
};

namespace AlgoimMomentDetail {

// Evaluate a tensor power polynomial, with last index fastest.
// For QUAD9, the third coordinate is zero, and 2D polynomials
// are stored as a cube of side length s with nonzero k=0 entries.
inline long double EvalCube(
    const std::vector<long double>& coefficients,
    const std::array<long double, 3>& x) {
  const unsigned s = SurfacePolynomialCubeSize(
      static_cast<unsigned>(coefficients.size()));
  long double value = 0.0L;
  for(int i = static_cast<int>(s) - 1; i >= 0; --i) {
    long double row = 0.0L;
    for(int j = static_cast<int>(s) - 1; j >= 0; --j) {
      long double column = 0.0L;
      for(int k = static_cast<int>(s) - 1; k >= 0; --k)
        column = column * x[2] +
                 coefficients[(static_cast<unsigned>(i) * s +
                               static_cast<unsigned>(j)) * s +
                              static_cast<unsigned>(k)];
      row = row * x[1] + column;
    }
    value = value * x[0] + row;
  }
  return value;
}

// All monomials of TOTAL degree <= degree; exact FEMUS ordering.
// 'powers' is workspace allocated once per call, not per quadrature point.
template<int Dim>
inline void AddMoments(
    const algoim::uvector<double, Dim>& u,
    long double w,
    unsigned degree,
    std::vector<long double>& moments,
    std::array<std::vector<long double>, Dim>& powers) {
  for(int d = 0; d < Dim; ++d) {
    powers[d][0] = 1.0L;
    for(unsigned n = 1; n <= degree; ++n)
      powers[d][n] = powers[d][n-1U] * static_cast<long double>(u(d));
  }

  for(unsigned i = 0; i <= degree; ++i)
    for(unsigned j = 0; j <= degree-i; ++j)
      for(unsigned k = 0; k <= (Dim == 2 ? 0U : degree-i-j); ++k) {
        long double value = w * powers[0][i] * powers[1][j];
        if constexpr(Dim == 3) value *= powers[2][k];
        moments[UnitMomentIndex(Dim, i, j, k)] += value;
      }
}

// Tensor power basis -> quadratic Bernstein basis, with positive
// normalization. The domain and zero level set are unchanged.
template<int Dim>
inline std::vector<double> BuildBernstein(
    const std::vector<long double>& original) {
  static_assert(Dim == 2 || Dim == 3, "Only QUAD9 and HEX27");
  const unsigned count = Dim == 2 ? 9U : 27U;
  std::vector<long double> a(count, 0.0L);

  for(unsigned i = 0; i < 3; ++i)
    for(unsigned j = 0; j < 3; ++j)
      for(unsigned k = 0; k < (Dim == 2 ? 1U : 3U); ++k) {
        const unsigned n = Dim == 2 ? 3U*i+j : 9U*i+3U*j+k;
        a[n] = original[9U*i+3U*j+k];
      }

  long double scale = 0.0L;
  for(const long double v : a) {
    if(!std::isfinite(v))
      throw std::invalid_argument("Algoim: nonfinite level-set coefficient");
    scale = std::max(scale, std::abs(v));
  }
  if(!(scale > 0.0L))
    throw std::domain_error("Algoim: identically zero level-set polynomial");
  for(long double& v : a) v /= scale;

  // Degree-two power -> Bernstein, independently along each axis:
  // b0=a0, b1=a0+a1/2, b2=a0+a1+a2.
  for(int axis = 0; axis < Dim; ++axis) {
    const unsigned stride = axis == 0 ? (Dim == 2 ? 3U : 9U) :
                            axis == 1 ? (Dim == 2 ? 1U : 3U) : 1U;
    std::vector<long double> b(count);
    for(unsigned n = 0; n < count; ++n) {
      const unsigned power = (n / stride) % 3U;
      const unsigned start = n - power * stride;
      b[n] = a[start];
      if(power >= 1U) b[n] += (power == 1U ? 0.5L : 1.0L)*a[start+stride];
      if(power == 2U) b[n] += a[start+2U*stride];
    }
    a.swap(b);
  }

  std::vector<double> bernstein(count);
  for(unsigned i = 0; i < count; ++i) {
    bernstein[i] = static_cast<double>(a[i]);
    if(!std::isfinite(bernstein[i]))
      throw std::overflow_error("Algoim: nonfinite Bernstein coefficient");
  }
  return bernstein;
}

inline std::vector<long double> PlanarCube(const std::vector<long double>& coefficients) {
  unsigned side = 1;
  while(side*side < coefficients.size()) ++side;
  if(side*side != coefficients.size()) throw std::invalid_argument("Invalid planar polynomial coefficient count.");
  std::vector<long double> result(side*side*side,0);
  for(unsigned i = 0; i < side; ++i) for(unsigned j = 0; j < side; ++j) result[(i*side+j)*side] = coefficients[i*side+j];
  return result;
}

template<int Dim>
inline AlgoimMomentResult Integrate(
    const std::vector<double>& nodalValues,
    unsigned degree,
    unsigned order,
    AlgoimSurfaceMode mode,
    const std::vector<std::vector<double>>* geometry) {
  const bool usesGeometry = mode == AlgoimSurfaceMode::PhysicalSurface ||
                            mode == AlgoimSurfaceMode::InverseSurface;
  if(usesGeometry != (geometry != nullptr))
    throw std::invalid_argument("Algoim: geometry must be supplied only for physical/inverse surface modes");
  if(geometry != nullptr)
    ValidateSurfaceGeometryCoordinates(*geometry, Dim, Dim == 2 ? 9U : 27U);

  const CutFemLevelSetData field(Dim == 2 ? QUAD : HEX,nodalValues);
  const std::vector<long double> original = Dim == 2 ? PlanarCube(field.GetPowerCoefficients()) : field.GetPowerCoefficients();
  std::vector<double> control = BuildBernstein<Dim>(original);
  const algoim::uvector<int,Dim> extent(3);
  algoim::xarray<double,Dim> polynomial(control.data(), extent);
  algoim::ImplicitPolyQuadrature<Dim> quadrature(polynomial);

  // Algoim's surface callback measures dS_u, independently of the
  // positive normalization used to construct its Bernstein polynomial.
  std::vector<std::vector<long double>> gradG;
  std::vector<std::vector<long double>> cofGradG;
  if(mode != AlgoimSurfaceMode::ReferenceSurface) for(unsigned axis = 0; axis < Dim; ++axis) gradG.push_back(Dim == 2 ? PlanarCube(field.GetDerivativeCoefficients(axis)) : field.GetDerivativeCoefficients(axis));
  if(usesGeometry) {
    PreparedMomentWorkspace workspace;
    PrepareMetric(field,*geometry,workspace);
    cofGradG = workspace.metric;
    if(Dim == 2) for(std::vector<long double>& component : cofGradG) component = PlanarCube(component);
  }

  const unsigned count = Dim == 2 ?
      (degree+1U)*(degree+2U)/2U :
      (degree+1U)*(degree+2U)*(degree+3U)/6U;
  AlgoimMomentResult result;
  result.positive.assign(count,0.0L);
  result.surface.assign(count,0.0L);

  std::array<std::vector<long double>,Dim> powers;
  for(auto& p : powers) p.resize(degree+1U);

  // integrate() emits points in BOTH phases. Keep only G_h > 0.
  quadrature.integrate(algoim::AutoMixed,static_cast<int>(order),
    [&](const algoim::uvector<double,Dim>& u, double w) {
      ++result.volumePointsGenerated;
      if(algoim::bernstein::evalBernsteinPoly(polynomial,u) <= 0.0)
        return;
      ++result.volumePointsUsed;
      AddMoments<Dim>(u,static_cast<long double>(w),degree,
                      result.positive,powers);
    });

  // Normalization convention already used by CutFemWeight_LS::_fPair.
  const long double factor =
      mode == AlgoimSurfaceMode::Delta ? 1.0L :
      mode == AlgoimSurfaceMode::ReferenceSurface ? 0.5L :
      mode == AlgoimSurfaceMode::PhysicalSurface ?
          1.0L/static_cast<long double>(1U<<Dim) :
          static_cast<long double>(1U<<Dim);

  quadrature.integrate_surf(algoim::AutoMixed,static_cast<int>(order),
    [&](const algoim::uvector<double,Dim>& u, double w,
        const algoim::uvector<double,Dim>& /*weightedNormal*/) {
      ++result.surfacePoints;
      long double weight = static_cast<long double>(w) * factor;
      if(mode != AlgoimSurfaceMode::ReferenceSurface) {
        std::array<long double,3> x = {{
            static_cast<long double>(u(0)),
            static_cast<long double>(u(1)), 0.0L}};
        if constexpr(Dim == 3) x[2] = static_cast<long double>(u(2));
        long double g = 0.0L;
        for(const auto& component : gradG)
          g = std::hypot(g,EvalCube(component,x));
        if(!(g > 0.0L) || !std::isfinite(g))
          throw std::domain_error("Algoim: singular level-set gradient on interface");

        if(usesGeometry) {
          long double c = 0.0L;
          for(const auto& component : cofGradG)
            c = std::hypot(c,EvalCube(component,x));
          if(!(c > 0.0L) || !std::isfinite(c))
            throw std::domain_error("Algoim: singular cofactor-weighted gradient on interface");
          weight *= mode == AlgoimSurfaceMode::PhysicalSurface ? c/g : 1.0L/(g*c);
        } else {
          weight /= g; // delta(G) du from Algoim's dS_u
        }
      }
      if(!std::isfinite(weight))
        throw std::overflow_error("Algoim: nonfinite surface moment weight");
      AddMoments<Dim>(u,weight,degree,result.surface,powers);
    });

  return result;
}

} // namespace AlgoimMomentDetail

// The only public entry point. Exactly nine or twenty-seven nodal values.
inline AlgoimMomentResult IntegrateAlgoimMoments(
    const std::vector<double>& nodalValues,
    unsigned momentDegree,
    unsigned algoimOrder,
    AlgoimSurfaceMode surfaceMode,
    const std::vector<std::vector<double>>* geometry = nullptr) {
  if(momentDegree > 64U)
    throw std::invalid_argument("Algoim: moment degree exceeds 64");
  if(algoimOrder < 1U ||
     algoimOrder > static_cast<unsigned>(algoim::GaussQuad::p_max) ||
     algoimOrder > static_cast<unsigned>(algoim::TanhSinhQuadrature::n_max))
    throw std::invalid_argument("Algoim: quadrature order must be between 1 and 100");

  if(nodalValues.size() == 9U)
    return AlgoimMomentDetail::Integrate<2>(
        nodalValues,momentDegree,algoimOrder,surfaceMode,geometry);
  if(nodalValues.size() == 27U)
    return AlgoimMomentDetail::Integrate<3>(
        nodalValues,momentDegree,algoimOrder,surfaceMode,geometry);
  throw std::invalid_argument("Algoim: expected QUAD9 or HEX27 nodal values");
}

} // namespace FELevelSetMoments

#endif // FEMUS_ALGOIM_MOMENTS_HPP
