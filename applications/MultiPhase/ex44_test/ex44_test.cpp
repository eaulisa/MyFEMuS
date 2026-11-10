#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "Fem.hpp"
#include "GeomElTypeEnum.hpp"
#include "VTKWriter.hpp"
#include "../includeLS/Psi.hpp"
#include "IntegrateCut.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstddef>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace femus;

template <class PsiFunction>
void FlagFinestMeshLevel(MultiLevelMesh& mlMsh, const PsiFunction& psi);

template <class PsiFunction>
void InitLevelSet(MultiLevelSolution& mlSol, const std::string& name, const PsiFunction& psi);

// HEX27 Lauricella reference: multiprecision evaluation, independent of CutFEM.
namespace detail_hex {

using Real = boost::multiprecision::number<
    boost::multiprecision::cpp_dec_float<100>>;
using Polynomial = std::vector<Real>;

struct SurfaceResult {
    Real area;
    Real truncation_bound;
    std::size_t panels;
    int maxorder;
    Real target_absolute_truncation_bound;
    Real max_majorant;
};

inline Real abs_real(const Real& x) {
    using std::abs;
    return abs(x);
}

inline Real sqrt_real(const Real& x) {
    using std::sqrt;
    return sqrt(x);
}

inline Polynomial multiply(const Polynomial& p, const Polynomial& q) {
    Polynomial result(p.size() + q.size() - 1, Real(0));
    for (std::size_t i = 0; i < p.size(); ++i)
        for (std::size_t j = 0; j < q.size(); ++j)
            result[i+j] += p[i] * q[j];
    return result;
}

inline Polynomial derivative(const Polynomial& p) {
    Polynomial result(p.size() - 1, Real(0));
    for (std::size_t i = 1; i < p.size(); ++i)
        result[i-1] = Real(i) * p[i];
    return result;
}

struct LocalPolynomials {
    Polynomial p;
    Polynomial q;
    Real f0;
    Real p_majorant;
    Real q_majorant;
};

inline LocalPolynomials local_polynomials(
        const Real& center, const Real& radius,
        const Real& chi, const Real& aspect) {
    // First form coefficients in u=t-center, preserving the positive
    // unexpanded expression for P(center). Then normalize v=u/radius.
    const Polynomial t2{center*center, Real(2)*center, Real(1)};
    const Polynomial q{Real(1)+chi*t2[0], chi*t2[1], chi};
    const Polynomial one_minus_t2{Real(1)-t2[0], -t2[1], Real(-1)};
    Polynomial p = multiply(one_minus_t2, multiply(q, multiply(q, q)));
    const Real multiplier = aspect*(Real(1)+chi)*(Real(1)+chi);
    for (std::size_t j = 0; j < 3; ++j)
        p[j] += multiplier*t2[j];
    const Real p0 = p[0];
    const Real q0 = q[0];
    if (p0 <= 0 || q0 <= 0)
        throw std::runtime_error("Nonpositive HEX local polynomial value.");

    LocalPolynomials out;
    out.p.resize(p.size());
    out.q.resize(q.size());
    Real power = 1;
    for (std::size_t j = 0; j < p.size(); ++j) {
        out.p[j] = p[j]*power/p0;
        power *= radius;
    }
    power = 1;
    for (std::size_t j = 0; j < q.size(); ++j) {
        out.q[j] = q[j]*power/q0;
        power *= radius;
    }
    out.p[0] = out.q[0] = 1;
    out.p_majorant = out.q_majorant = 0;
    for (std::size_t j = 1; j < out.p.size(); ++j)
        out.p_majorant += abs_real(out.p[j]);
    for (std::size_t j = 1; j < out.q.size(); ++j)
        out.q_majorant += abs_real(out.q[j]);
    out.f0 = sqrt_real(p0)/(q0*q0);
    return out;
}

struct PanelResult {
    Real integral;
    Real tail_bound;
    int order;
};

inline PanelResult panel_series(
        const LocalPolynomials& loc, const Real& width,
        const Real& tolerance_per_unit_width) {
    // Normalized real interval: v in [-1/2,1/2]. On |v|=1,
    // |g(v)| <= M, where g=sqrt(P/P0)/(Q/Q0)^2 and g(0)=1.
    const Real rho_squared = Real(1)/4;
    const Real one_minus_rho_squared = Real(3)/4;
    const Real denominator = Real(1)-loc.q_majorant;
    const Real major_bound = sqrt_real(Real(1)+loc.p_majorant)
                           /(denominator*denominator);
    int order = 0;
    Real tail_power = rho_squared;  // rho^(order+2), rho=1/2.
    while (loc.f0*major_bound*tail_power
               /(Real(order+3)*one_minus_rho_squared)
           > tolerance_per_unit_width) {
        order += 2;
        tail_power *= rho_squared;
        if (order > 20000)
            throw std::runtime_error("Excessive HEX Taylor order.");
    }

    // Polynomial ODE: A*g'=B*g,
    // A=2*P*Q and B=Q*P'-4*P*Q', for normalized P and Q.
    Polynomial aa = multiply(loc.p, loc.q);
    for (Real& coefficient : aa)
        coefficient *= 2;
    const Polynomial b_left = multiply(loc.q, derivative(loc.p));
    const Polynomial b_right = multiply(loc.p, derivative(loc.q));
    Polynomial bb(std::max(b_left.size(), b_right.size()), Real(0));
    for (std::size_t j = 0; j < bb.size(); ++j) {
        if (j < b_left.size())
            bb[j] += b_left[j];
        if (j < b_right.size())
            bb[j] -= 4*b_right[j];
    }

    Polynomial coefficients;
    coefficients.reserve(static_cast<std::size_t>(order)+1);
    coefficients.push_back(Real(1));
    for (std::size_t m = 0; m < static_cast<std::size_t>(order); ++m) {
        Real right = 0;
        Real left = 0;
        const std::size_t b_max = std::min(m, bb.size()-1);
        for (std::size_t k = 0; k <= b_max; ++k)
            right += bb[k]*coefficients[m-k];
        const std::size_t a_max = std::min(m, aa.size()-1);
        for (std::size_t k = 1; k <= a_max; ++k)
            left += aa[k]*Real(m-k+1)*coefficients[m-k+1];
        coefficients.push_back((right-left)/(aa[0]*Real(m+1)));
    }
    // Odd coefficients integrate to zero on the symmetric panel.
    Real sum = 0;
    Real power = 1;
    for (int j = 0; j <= order; j += 2) {
        sum += coefficients[static_cast<std::size_t>(j)]*power/Real(j+1);
        power *= rho_squared;
    }
    PanelResult out;
    out.integral = width*loc.f0*sum;
    out.tail_bound = width*loc.f0*major_bound*tail_power
                   /(Real(order+3)*one_minus_rho_squared);
    out.order = order;
    return out;
}

inline SurfaceResult evaluate_lauricella_fd_hex(
        const Real& chi, const Real& a, const Real& b, int digits = 40) {
    if (digits < 10 || digits > 60)
        throw std::invalid_argument("HEX surface digits must be in [10,60].");
    if (!boost::math::isfinite(chi) || !boost::math::isfinite(a)
        || !boost::math::isfinite(b) || chi < 0 || a <= 0 || b <= 0)
        throw std::invalid_argument("Require finite chi>=0 and a,b>0.");

    const Real pi = boost::math::constants::pi<Real>();
    const Real aspect = b/a;
    const Real prefactor = 4*pi/sqrt_real(a*b);
    const Real lower_bound = 2*pi/a;
    Real relative_target = 1;
    for (int i = 0; i < digits; ++i)
        relative_target /= 10;
    const Real target_surface = lower_bound*relative_target;
    const Real target_integral = target_surface/prefactor;

    SurfaceResult result;
    result.area = result.truncation_bound = result.max_majorant = 0;
    result.panels = 0;
    result.maxorder = 0;
    result.target_absolute_truncation_bound = target_surface;
    if (chi == 0 && a == b) {
        result.area = 4*pi/a;
        return result;
    }

    const Real majorant_limit("0.4");
    Real left = 0;
    Real width = 1;
    Real total_integral = 0;
    Real total_bound = 0;
    while (left < 1) {
        const Real remaining = 1-left;
        if (width > remaining)
            width = remaining;
        LocalPolynomials loc;
        Real right;
        for (;;) {
            right = left+width;
            if (right == left)
                throw std::runtime_error("Insufficient precision for HEX panel width.");
            const Real center = (left+right)/2;
            // R=width=2*(real panel half-width).
            loc = local_polynomials(center, width, chi, aspect);
            if (loc.p_majorant <= majorant_limit
                && loc.q_majorant <= majorant_limit)
                break;
            width /= 2;
        }

        const PanelResult panel = panel_series(loc, width, target_integral);
        total_integral += panel.integral;
        total_bound += panel.tail_bound;
        ++result.panels;
        result.maxorder = std::max(result.maxorder, panel.order);
        if (loc.p_majorant > result.max_majorant)
            result.max_majorant = loc.p_majorant;
        if (loc.q_majorant > result.max_majorant)
            result.max_majorant = loc.q_majorant;
        if (result.panels > 20000)
            throw std::runtime_error("Too many HEX continuation panels.");
        left = right;
        const Real doubled_width = 2*width;
        const Real new_remaining = 1-left;
        width = (doubled_width < new_remaining) ? doubled_width : new_remaining;
    }
    result.area = prefactor*total_integral;
    result.truncation_bound = prefactor*total_bound;
    return result;
}

} // namespace detail_hex

// ============================================================================
// PUBLIC, CLEAN API.
// These are the functions/parameters you actually use.
// ============================================================================

namespace cutfem_reference {

using Real = boost::multiprecision::cpp_dec_float_100;

inline Real pi() { return boost::math::constants::pi<Real>(); }

// ----------------------------------------------------------------------------
// HEX27
// ----------------------------------------------------------------------------
// Reference cell:
//     [-1,1]^3
//
// LEVEL SET:
//     r2 = x^2 + y^2
//     phi(x,y,z) = a*r2 + b*z^2 + c*r2*z^2 - 1
//
// SHAPE PARAMETERS:
//     a > 0 : controls equatorial radius    R = 1/sqrt(a)
//     b > 0 : controls vertical semi-height H = 1/sqrt(b)
//     c >=0 : controls the pinching/non-convexity
//
// Useful dimensionless parameter:
//     chi = c/(a*b)
//
// Convexity:
//     chi <= 3 : convex
//     chi >  3 : non-convex
//
// To keep the full shape inside [-1,1]^3, use a>=1 and b>=1.
//
// Limit c -> infinity at fixed a,b:
//     volume  -> 0
//     surface -> 2*pi/a
// and the shape collapses to
//     disk z=0, x^2+y^2 <= 1/a
// plus the axial segment x=y=0, |z|<=1/sqrt(b).
//
// This is Q2 in each coordinate, hence exactly representable by HEX27.
// ----------------------------------------------------------------------------

inline double hex_levelset(double x, double y, double z,
                           double a, double b, double c) {
    const double r2 = x*x + y*y;
    return a*r2 + b*z*z + c*r2*z*z - 1.0;
}

struct HexReference {
    Real volume;
    Real surface;
    Real chi;
    bool nonconvex;
};

// Analytic volume formula.
inline Real hex_volume(const Real& a, const Real& b, const Real& c) {
    using boost::multiprecision::abs;
    using boost::multiprecision::atan;
    using boost::multiprecision::sqrt;

    if (a <= 0 || b <= 0 || c < 0)
        throw std::invalid_argument("HEX requires a>0, b>0, c>=0.");

    const Real chi = c/(a*b);
    if (chi == 0)
        return 4*pi()/(3*a*sqrt(b));

    // Stable analytic series for very small chi.
    if (chi <= Real("0.1")) {
        Real term = Real(1)/3;
        Real sum = term;
        for (std::size_t n=0; n<10000; ++n) {
            const Real next = -chi*Real(2*n+1)/Real(2*n+5)*term;
            if (abs(next) <= Real("1e-95")*abs(sum))
                return 4*pi()*sum/(a*sqrt(b));
            sum += next;
            term = next;
        }
        throw std::runtime_error("HEX volume series did not converge.");
    }

    const Real q = sqrt(chi);
    return 2*pi()/(a*sqrt(b))
           *((Real(1)+Real(1)/chi)*atan(q)/q - Real(1)/chi);
}

// Analytic surface formula:
//
//   chi = c/(a*b),  k=b/a,
//
//   A1 = 3*chi - 1 + k*(1+chi)^2
//   A2 = 3*chi*(chi-1)
//   A3 = chi^2*(chi-3)
//   A4 = -chi^3
//
// Let r1,...,r4 be the roots of
//   r^4 + A1*r^3 + A2*r^2 + A3*r + A4 = 0.
//
// Then exactly
//
//   S = 4*pi/sqrt(a*b)
//       * F_D^(5)(1/2;
//                 -1/2,-1/2,-1/2,-1/2,2;
//                 3/2;
//                 r1,r2,r3,r4,-chi).
//
// F_D^(5) is the Lauricella analytic special function.
// detail_hex::evaluate_lauricella_fd_hex() below evaluates this closed form
// in multiprecision. The evaluator internals are deliberately hidden from
// the normal output.
inline Real hex_surface(const Real& a, const Real& b, const Real& c) {
    if (a <= 0 || b <= 0 || c < 0)
        throw std::invalid_argument("HEX requires a>0, b>0, c>=0.");

    const Real chi = c/(a*b);
    const auto result = detail_hex::evaluate_lauricella_fd_hex(
        chi, a, b, 50);
    return result.area;
}

inline HexReference hex_reference(const Real& a, const Real& b, const Real& c) {
    if (a < 1 || b < 1)
        throw std::invalid_argument(
            "For a shape fully contained in [-1,1]^3 use a>=1 and b>=1.");
    const Real chi = c/(a*b);
    return {hex_volume(a,b,c), hex_surface(a,b,c), chi, chi>3};
}

} // namespace cutfem_reference



namespace {

// Exact-analytic curved BALL: psi(x) = r^2 - |F^{-1}(x)-c|^2.
// The polynomial F is implemented in Mesh for every element type.
struct AnalyticBallPsi {
  GeomElType geometry;
  double curvature;
  std::array<double,3> center;
  double radius;

  double operator()(const std::vector<double>& physical) const {
    const unsigned dim=(geometry==TET || geometry==HEX) ? 3u : 2u;
    if(physical.size()!=dim)
      throw std::invalid_argument("AnalyticBallPsi: incorrect dimension");

    std::array<double,3> x={{0.,0.,0.}};
    for(unsigned d=0;d<dim;++d) x[d]=physical[d];
    const auto q=Mesh::InvertAnalyticBallMap(geometry,center,curvature,x);

    double psi=radius*radius;
    for(unsigned d=0;d<dim;++d) {
      const double z=q[d]-center[d];
      psi-=z*z;
    }
    return psi;
  }
};

// Straight QUAD9 (Stella Q2 oval) and HEX27 (Lauricella surface).
// HEX analytical references use [-1,1]^3; the FEMuS element is
// [-0.5,0.5]^3, so the level-set evaluator uses (X,Y,Z)=2*(x,y,z).
// Positive Psi is the interior phase.
struct NonconvexPsi {
  GeomElType element;
  double quadA, quadB, quadC, hexA, hexB, hexC;

  double operator()(const std::vector<double>& x) const {
    if(element == QUAD) {
      if(x.size() != 2u)
        throw std::invalid_argument("Nonconvex QUAD requires 2D");
      const double x2 = x[0]*x[0], y2 = x[1]*x[1];
      return 1. - quadA*x2 - quadB*y2 - quadC*x2*y2;
    }
    if(element == HEX) {
      if(x.size() != 3u)
        throw std::invalid_argument("Nonconvex HEX requires 3D");
      const double X=2.*x[0], Y=2.*x[1], Z=2.*x[2];
      const double radial=X*X+Y*Y;
      return 1. - hexA*radial - hexB*Z*Z - hexC*radial*Z*Z;
    }
    throw std::invalid_argument("Nonconvex requires QUAD9 or HEX27");
  }
};

// Gauss-Legendre integration of a smooth, independent parametrization of Gamma.
// This is only a reference calculation: the CutFEM integrator is not used.
struct NonconvexQuadGaussLegendre {
  std::vector<long double> nodes, weights;

  explicit NonconvexQuadGaussLegendre(unsigned n) : nodes(n), weights(n) {
    if(n < 2u) throw std::invalid_argument("Gauss-Legendre order must be >= 2");
    const long double pi = std::acos(-1.L);

    for(unsigned i = 0; i < (n + 1u) / 2u; ++i) {
      long double z = std::cos(pi * (static_cast<long double>(i) + 0.75L)
                                     / (static_cast<long double>(n) + 0.5L));
      long double derivative = 0.L;
      for(unsigned iter = 0; iter < 100u; ++iter) {
        long double p0 = 1.L, p1 = z;
        for(unsigned k = 2; k <= n; ++k) {
          const long double p2 = ((2.L * k - 1.L) * z * p1
                                 - (k - 1.L) * p0) / k;
          p0 = p1;
          p1 = p2;
        }
        derivative = n * (z * p1 - p0) / (z * z - 1.L);
        const long double next = z - p1 / derivative;
        if(std::abs(next - z) < 3.e-19L) {
          z = next;
          break;
        }
        z = next;
      }

      // Recompute P_n'(z) at the final root for accurate weights.
      long double p0 = 1.L, p1 = z;
      for(unsigned k = 2; k <= n; ++k) {
        const long double p2 = ((2.L * k - 1.L) * z * p1
                                - (k - 1.L) * p0) / k;
        p0 = p1;
        p1 = p2;
      }
      derivative = n * (z * p1 - p0) / (z * z - 1.L);
      const long double w = 2.L / ((1.L - z * z) * derivative * derivative);
      nodes[i] = -z;
      nodes[n - 1u - i] = z;
      weights[i] = weights[n - 1u - i] = w;
    }
  }
};

void ValidateNonconvexQuad(double a, double b, double c) {
  if(!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c)
     || a <= 0. || b <= 0. || c < 0.)
    throw std::invalid_argument("Nonconvex QUAD requires finite a>0, b>0, c>=0");

  // Our single straight QUAD9 occupies [-0.5,0.5]^2, not [-1,1]^2.
  // If this condition fails, the full-oval formula is NOT the cut-cell reference.
  if(a < 4. || b < 4.)
    throw std::invalid_argument(
        "Nonconvex QUAD full-oval references require a>=4 and b>=4 "
        "in the FEMuS element [-0.5,0.5]^2");
}

// Closed form A = pi/sqrt(ab) * 2F1(1/2,1/2;2;-q), q=c/(ab).
// For q>0 use the equivalent real-modulus elliptic-integral formula:
// A = 4*sqrt(1+q)/(q*sqrt(ab)) * [K(sqrt(m))-E(sqrt(m))], m=q/(1+q).
// For small q, use the convergent hypergeometric Taylor series instead
// to avoid subtractive cancellation in K-E.
double NonconvexQuadExactArea(double a, double b, double c) {
  const long double aa = a, bb = b, cc = c;
  const long double q = cc / (aa * bb);
  const long double pi = std::acos(-1.L);
  const long double scale = pi / std::sqrt(aa * bb);
  if(q == 0.L) return static_cast<double>(scale);

  if(q < 0.05L) {
    long double sum = 1.L, term = 1.L;
    for(unsigned n = 1u; n <= 128u; ++n) {
      const long double nn = static_cast<long double>(n);
      term *= -q * (nn - 0.5L) * (nn - 0.5L) / (nn * (nn + 1.L));
      sum += term;
      if(std::abs(term) < 1.e-20L * std::abs(sum)) break;
    }
    return static_cast<double>(scale * sum);
  }

  const long double modulus = std::sqrt(q / (1.L + q));
  const long double K = std::comp_ellint_1(modulus);
  const long double E = std::comp_ellint_2(modulus);
  return static_cast<double>(4.L * std::sqrt(1.L + q)
                             * (K - E) / (q * std::sqrt(aa * bb)));
}

// The perimeter has an exact regular one-dimensional integral but generally
// no simple closed special-function formula.  Compute its reference value
// by independently converged high-order quadrature on that parameterization.
long double NonconvexQuadPerimeterAtOrder(double a, double b, double c, unsigned order) {
  const long double aa = a, bb = b, cc = c;
  const long double q = cc / (aa * bb);
  const long double pi = std::acos(-1.L);
  const long double halfPi = pi / 2.L;
  std::vector<long double> intervals = {
      0.L, pi / 8.L, pi / 4.L, 3.L * pi / 8.L, halfPi};
  if(q > 0.L) {
    const long double transition = std::atan(1.L / std::sqrt(q));
    for(const long double factor : {0.25L, 1.L, 4.L}) {
      const long double candidate = factor * transition;
      if(candidate > 0.L && candidate < halfPi) intervals.push_back(candidate);
    }
  }
  std::sort(intervals.begin(), intervals.end());
  intervals.erase(std::unique(intervals.begin(), intervals.end(),
      [](long double x, long double y) { return std::abs(x-y) < 1.e-17L; }),
      intervals.end());

  NonconvexQuadGaussLegendre rule(order);
  const auto speed = [&](long double t) {
    const long double st = std::sin(t), ct = std::cos(t);
    const long double denom = 1.L + q * st * st;
    const long double dx = ct / std::sqrt(aa);
    const long double dy = -(1.L + q) * st
        / (std::sqrt(bb) * denom * std::sqrt(denom));
    return std::sqrt(dx*dx + dy*dy);
  };

  long double result = 0.L;
  for(unsigned segment = 1u; segment < intervals.size(); ++segment) {
    const long double left = intervals[segment-1u];
    const long double right = intervals[segment];
    const long double middle = (left + right) / 2.L;
    const long double halfWidth = (right - left) / 2.L;
    for(unsigned i = 0u; i < order; ++i) {
      result += halfWidth * rule.weights[i]
          * speed(middle + halfWidth * rule.nodes[i]);
    }
  }
  return 4.L * result;
}

double NonconvexQuadReferencePerimeter(double a, double b, double c) {
  const long double p128 = NonconvexQuadPerimeterAtOrder(a, b, c, 128);
  const long double p256 = NonconvexQuadPerimeterAtOrder(a, b, c, 256);
  if(std::abs(p256-p128) > 3.e-13L * std::abs(p256))
    throw std::runtime_error("Nonconvex QUAD perimeter reference has not converged");
  return static_cast<double>(p256);
}

// NONCONVEX is intentionally restricted to straight QUAD9 and HEX27.
void ValidateNonconvex(const std::vector<GeomElType>& elements,
                       double quadA, double quadB, double quadC,
                       double hexA, double hexB, double hexC) {
  for(const GeomElType e : elements) {
    if(e != QUAD && e != HEX)
      throw std::invalid_argument(
          "--psi nonconvex supports only QUAD9 and HEX27; "
          "use --elements quad,hex");
  }
  if(std::find(elements.begin(), elements.end(), QUAD) != elements.end())
    ValidateNonconvexQuad(quadA, quadB, quadC);
  if(std::find(elements.begin(), elements.end(), HEX) != elements.end()) {
    if(!(std::isfinite(hexA) && std::isfinite(hexB) && std::isfinite(hexC)
         && hexA >= 1. && hexB >= 1. && hexC >= 0.))
      throw std::invalid_argument(
          "Nonconvex HEX requires a,b>=1 and c>=0 for full-domain references");
    if(!(hexC/(hexA*hexB) > 3.))
      throw std::invalid_argument("Nonconvex HEX requires chi=c/(a*b)>3");
  }
}

std::vector<IntegrationTest> MakeNonconvexTests(
    GeomElType element,
    double quadA, double quadB, double quadC,
    double hexA, double hexB, double hexC) {
  using cutfem_reference::Real;
  double volume=0., surface=0.;
  if(element == QUAD) {
    ValidateNonconvexQuad(quadA, quadB, quadC);
    volume=NonconvexQuadExactArea(quadA, quadB, quadC);
    surface=NonconvexQuadReferencePerimeter(quadA, quadB, quadC);
  } else if(element == HEX) {
    const auto r = cutfem_reference::hex_reference(
        Real(hexA), Real(hexB), Real(hexC));
    // The source cube is [-1,1]^3, while FEMuS uses [-0.5,0.5]^3.
    // Physical volume and surface measures are divided by 8 and 4.
    volume=(r.volume/Real(8)).convert_to<double>();
    surface=(r.surface/Real(4)).convert_to<double>();
  } else {
    throw std::invalid_argument("Nonconvex supports only QUAD9 and HEX27");
  }
  if(!std::isfinite(volume) || !std::isfinite(surface)
     || !(volume>0.) || !(surface>0.))
    throw std::runtime_error("Nonconvex: invalid analytical reference");
  return {{"1", [](const std::vector<double>&) {return 1.;}, volume, surface}};
}

struct DriverOptions {
  std::vector<GeomElType> elements = {HEX, QUAD};
  std::vector<std::string> psis = {"ball"};
  bool curved = false;
  unsigned levels = 3; // Number of AMR refinements; also evaluate level 0.
  double curvature = 0.40; // Analytical quadratic deformation strength.
  double nonconvexQuadA = 4.0, nonconvexQuadB = 4.0, nonconvexQuadC = 160.0;
  double nonconvexHexA = 4.0, nonconvexHexB = 4.0, nonconvexHexC = 160.0;
  bool vtk = true;
  IntegrationOptions integration;
};

std::string Lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return s;
}

std::vector<std::string> Split(const std::string& s) {
  std::vector<std::string> out;
  std::istringstream in(s);
  std::string token;
  while(std::getline(in, token, ',')) {
    token.erase(std::remove_if(token.begin(), token.end(), [](unsigned char c) { return std::isspace(c); }), token.end());
    if(!token.empty()) out.push_back(Lower(token));
  }
  if(out.empty()) throw std::invalid_argument("Empty comma-separated list");
  return out;
}

unsigned Unsigned(const std::string& value, const std::string& option) {
  std::size_t used = 0;
  if(value.empty() || value[0] == '-') throw std::invalid_argument("Invalid " + option + ": " + value);
  const unsigned long x = std::stoul(value, &used);
  if(used != value.size() || x > std::numeric_limits<unsigned>::max()) throw std::invalid_argument("Invalid " + option + ": " + value);
  return static_cast<unsigned>(x);
}

double Real(const std::string& value, const std::string& option) {
  std::size_t used = 0;
  const double x = std::stod(value, &used);
  if(used != value.size() || !std::isfinite(x)) throw std::invalid_argument("Invalid " + option + ": " + value);
  return x;
}

bool Boolean(const std::string& value) {
  const std::string s = Lower(value);
  if(s == "1" || s == "true" || s == "yes" || s == "on") return true;
  if(s == "0" || s == "false" || s == "no" || s == "off") return false;
  throw std::invalid_argument("Expected true/false or 1/0, got: " + value);
}

GeomElType Element(const std::string& value) {
  if(value == "tri") return TRI;
  if(value == "quad") return QUAD;
  if(value == "tet") return TET;
  if(value == "hex") return HEX;
  throw std::invalid_argument("Unknown element: " + value + " (expected tri,quad,tet,hex)");
}

const char* ElementName(GeomElType value) {
  switch(value) {
    case TRI: return "tri";
    case QUAD: return "quad";
    case TET: return "tet";
    case HEX: return "hex";
    default: throw std::invalid_argument("Unsupported element type");
  }
}

std::string PsiKind(const std::string& value) {
  if(value == "ball" || value == "nonconvex") return value;
  throw std::invalid_argument("Unknown psi type: " + value + " (expected ball,nonconvex)");
}

const std::string& PsiName(const std::string& value) { return value; }

void Usage() {
  std::cout << "Usage: MultiPhase_ex44_test [options]\n"
            << "  --elements tri,quad,tet,hex  Element types (nonconvex: quad,hex only)\n"
            << "  --psi ball,nonconvex          Level-set types (nonconvex: QUAD/HEX only)\n"
            << "  --curved 0|1                 0=straight, 1=analytic curved BALL (default 0)\n"
            << "  --levels N                   Number of AMR refinements, plus level 0 (default 3)\n"
            << "  --curvature X                Analytical BALL deformation (default 0.40)\n"
            << "  --nonconvex-quad-a X         QUAD x^2 coefficient (default 4)\n"
            << "  --nonconvex-quad-b X         QUAD y^2 coefficient (default 4)\n"
            << "  --nonconvex-quad-c X         QUAD x^2*y^2 coefficient (default 160)\n"
            << "  --nonconvex-hex-a X          HEX a>=1 in source [-1,1]^3 (default 4)\n"
            << "  --nonconvex-hex-b X          HEX b>=1 in source [-1,1]^3 (default 4)\n"
            << "  --nonconvex-hex-c X          HEX c with chi=c/(ab)>3 (default 160)\n"
            << "  --mode moment|algoim|compare Integration method (default moment)\n"
            << "  --degree N                   Moment-fitting degree (default 7)\n"
            << "  --cutfem-abs-tol X           CutFEM absolute tolerance (default 1e-5)\n"
            << "  --cutfem-rel-tol X           CutFEM relative tolerance (default 1e-5)\n"
            << "  --algoim-order N             Algoim quadrature order (default 5)\n"
            << "  --algoim-subdivisions N      Algoim subdivisions (default 0)\n"
            << "  --vtk 0|1                    Write VTK at each level (default 1)\n"
            << "  --help                       Show this help\n";
}

DriverOptions Parse(int argc, char** argv) {
  DriverOptions options;
  for(int i = 1; i < argc; ++i) {
    const std::string key = argv[i];
    if(key == "--help" || key == "-h") { Usage(); std::exit(0); }
    if(i + 1 >= argc) throw std::invalid_argument("Missing value for " + key);
    const std::string value = argv[++i];
    if(key == "--elements") {
      options.elements.clear();
      for(const auto& token : Split(value)) {
        const GeomElType e = Element(token);
        if(std::find(options.elements.begin(), options.elements.end(), e) == options.elements.end()) options.elements.push_back(e);
      }
    }
    else if(key == "--psi") {
      options.psis.clear();
      for(const auto& token : Split(value)) {
        const std::string p = PsiKind(token);
        if(std::find(options.psis.begin(), options.psis.end(), p) == options.psis.end()) options.psis.push_back(p);
      }
    }
    else if(key == "--curved") options.curved = Boolean(value);
    else if(key == "--levels") options.levels = Unsigned(value, key);
    else if(key == "--curvature") options.curvature = Real(value, key);
    else if(key == "--nonconvex-quad-a") options.nonconvexQuadA = Real(value, key);
    else if(key == "--nonconvex-quad-b") options.nonconvexQuadB = Real(value, key);
    else if(key == "--nonconvex-quad-c") options.nonconvexQuadC = Real(value, key);
    else if(key == "--nonconvex-hex-a") options.nonconvexHexA = Real(value,key);
    else if(key == "--nonconvex-hex-b") options.nonconvexHexB = Real(value,key);
    else if(key == "--nonconvex-hex-c") options.nonconvexHexC = Real(value,key);
    else if(key == "--degree") options.integration.momentDegree = Unsigned(value, key);
    else if(key == "--cutfem-abs-tol") options.integration.absoluteTolerance = Real(value, key);
    else if(key == "--cutfem-rel-tol") options.integration.relativeTolerance = Real(value, key);
    else if(key == "--algoim-order") options.integration.algoimOrder = Unsigned(value, key);
    else if(key == "--algoim-subdivisions") options.integration.algoimSubdivisions = Unsigned(value, key);
    else if(key == "--vtk") options.vtk = Boolean(value);
    else if(key == "--mode") {
      const auto s = Lower(value);
      if(s == "moment") options.integration.mode = IntegrationMode::MomentFitting;
      else if(s == "algoim") options.integration.mode = IntegrationMode::Algoim;
      else if(s == "compare") options.integration.mode = IntegrationMode::Compare;
      else throw std::invalid_argument("Unknown integration mode: " + value);
    }
    else throw std::invalid_argument("Unknown option: " + key);
  }
  if(options.levels > 12) throw std::invalid_argument("--levels must be <= 12");
  if(options.integration.momentDegree == 0) throw std::invalid_argument("--degree must be positive");
  if(!(options.integration.absoluteTolerance > 0.))
    throw std::invalid_argument("--cutfem-abs-tol must be positive");
  if(options.integration.relativeTolerance < 0.)
    throw std::invalid_argument("--cutfem-rel-tol must be nonnegative");
  if(options.curvature < 0.) throw std::invalid_argument("--curvature must be nonnegative");
  if(options.curved && std::any_of(options.psis.begin(), options.psis.end(),
                                  [](const std::string& kind) { return kind != "ball"; }))
    throw std::invalid_argument("--curved 1 requires --psi ball (analytic geometry)");
  if(std::find(options.psis.begin(), options.psis.end(), "nonconvex") != options.psis.end()) {
    ValidateNonconvex(options.elements,
                      options.nonconvexQuadA, options.nonconvexQuadB, options.nonconvexQuadC,
                      options.nonconvexHexA, options.nonconvexHexB, options.nonconvexHexC);
  }
  return options;
}

// Straight BALL: exact references for the physical affine ball. All integrands
// are ordinary physical polynomials; both tests have closed-form integrals.
std::vector<IntegrationTest> MakeStraightBallTests(
    unsigned dim, const std::vector<double>& center, double radius) {
  if(dim != 2u && dim != 3u)
    throw std::invalid_argument("Straight BALL requires dimension 2 or 3");
  if(!(radius > 0.) || !std::isfinite(radius))
    throw std::invalid_argument("Straight BALL requires finite positive radius");

  const double pi = std::acos(-1.);
  const double r2 = radius*radius;
  const double r3 = r2*radius;
  const double r4 = r2*r2;

  return {
    {"1", [](const std::vector<double>&) { return 1.; },
      dim == 2u ? pi*r2 : 4.*pi*r3/3.,
      dim == 2u ? 2.*pi*radius : 4.*pi*r2},
    {"r^2", [center](const std::vector<double>& x) {
        double distance2 = 0.;
        for(unsigned d=0; d<x.size(); ++d) {
          const double v = x[d]-center[d];
          distance2 += v*v;
        }
        return distance2;
      },
      dim == 2u ? pi*r4/2. : 4.*pi*r4*radius/5.,
      dim == 2u ? 2.*pi*r3 : 4.*pi*r4}
  };
}

// CLOSED-FORM references for F(B_r(c)), with the polynomial geometry
//   2D: (p0,p1) -> (p0+k*(p0*p0-p1*p1), p1+2*k*p0*p1)
//   3D: (p0,p1,p2) -> ((1+2*k*p2)*p0, (1+2*k*p2)*p1,
//                       p2+k*(p2*p2-p0*p0-p1*p1)).
// Both tests use ordinary PHYSICAL polynomial integrands, no reciprocal metrics.
// E and K are complete elliptic integrals with modulus k (C++17).
std::vector<IntegrationTest> MakeAnalyticBallTests(
    unsigned dim, const std::vector<double>& center,
    double radius, double curvature) {
  if(dim!=2u && dim!=3u)
    throw std::invalid_argument("Analytic BALL: unsupported dimension");
  if(!(radius>0.) || !(curvature>=0.) ||
     !std::isfinite(radius) || !std::isfinite(curvature))
    throw std::invalid_argument("Analytic BALL: invalid radius/curvature");

  const double pi=std::acos(-1.);
  const double r=radius;
  const double lambda=2.*curvature*r;
  // In this range the mapped circle/sphere is regular and the simple
  // analytic expressions below apply without branch/sign changes.
  if(!(lambda<1.))
    throw std::invalid_argument("Analytic BALL requires 2*curvature*radius < 1");

  const double L2=lambda*lambda;
  const double L4=L2*L2;
  const double L6=L4*L2;
  const double r2=r*r;
  const double r3=r2*r;
  const double r4=r2*r2;

  double omega1=0.,gamma1=0.,omegaR2=0.,gammaR2=0.;

  if(dim==2u) {
    const double m=4.*lambda/((1.+lambda)*(1.+lambda));
    const double ellipticK=std::comp_ellint_1(std::sqrt(m));
    const double ellipticE=std::comp_ellint_2(std::sqrt(m));
    const double Q=(2.*(2.-m)*ellipticE-(1.-m)*ellipticK)/3.;

    omega1=pi*r2*(1.+0.5*L2);
    gamma1=4.*r*(1.+lambda)*ellipticE;

    omegaR2=0.5*pi*r4*(1.+1.5*L2+L4/8.);
    gammaR2=2.*r3*(1.+lambda)*(
        (1.-0.5*L2)*ellipticE+(1.+lambda)*(1.+lambda)*Q);
  } else {
    omega1=4.*pi*r3/3.*(1.+L2);
    gamma1=4.*pi*r2*(1.+2.*L2/3.-L4/15.);

    omegaR2=4.*pi*r3*r2/5.*(1.+25.*L2/12.+5.*L4/12.);
    gammaR2=4.*pi*r4*(1.+19.*L2/12.+L4/10.+L6/420.);
  }

  return {
    {"1", [](const std::vector<double>&) {return 1.;},
      omega1,gamma1},
    {"r^2", [center](const std::vector<double>& x) {
       double value=0.;
       for(unsigned d=0;d<x.size();++d) {
         const double z=x[d]-center[d];
         value+=z*z;
       }
       return value;
     }, omegaR2,gammaR2}
  };
}

void PrintHeader() {
  std::cout << std::string(90, '-') << '\n';
  std::cout << std::left << std::setw(8) << "Test" << std::setw(9) << "Method"
            << std::setw(9) << "Domain" << std::right << std::setw(16) << "Numerical"
            << std::setw(16) << "Exact" << std::setw(16) << "Abs. error"
            << std::setw(16) << "Rel. error" << '\n';
  std::cout << std::string(90, '-') << '\n';
}

void PrintRow(const std::string& test, const char* method, const char* domain,
              double numerical, double exact) {
  std::cout << std::left << std::setw(8) << test << std::setw(9) << method
            << std::setw(9) << domain << std::right << std::scientific << std::setprecision(8)
            << std::setw(16) << numerical;
  if(std::isfinite(exact)) {
    const double absoluteError = std::abs(numerical - exact);
    std::cout << std::setw(16) << exact << std::setw(16) << absoluteError;
    if(exact != 0.) std::cout << std::setw(16) << absoluteError / std::abs(exact);
    else std::cout << std::setw(16) << "N/A";
  } else {
    std::cout << std::setw(16) << "N/A" << std::setw(16) << "N/A"
              << std::setw(16) << "N/A";
  }
  std::cout << '\n';
}

void SaveResults(std::ofstream& out, unsigned level, const IntegrationTest& test,
                 const IntegrationResult& result) {
  const auto write = [&](const char* method, const char* domain, double numerical, double exact) {
    out << level << ',' << test.name << ',' << method << ',' << domain << ','
        << std::setprecision(17) << numerical << ',';
    if(std::isfinite(exact)) {
      const double absoluteError = std::abs(numerical - exact);
      out << exact << ',' << absoluteError << ',';
      if(exact != 0.) out << absoluteError / std::abs(exact);
      else out << "NA";
    } else {
      out << "NA,NA,NA";
    }
    out << ',' << result.cutElements << ',' << result.algoimVolumePoints
        << ',' << result.algoimSurfacePoints << '\n';
  };
  if(result.momentComputed) {
    write("moment", "omega", result.omega1, test.exactOmega1);
    write("moment", "gamma", result.gamma, test.exactGamma);
  }
  if(result.algoimComputed) {
    write("algoim", "omega", result.omega1Algoim, test.exactOmega1);
    write("algoim", "gamma", result.gammaAlgoim, test.exactGamma);
  }
}

void PrintResults(const IntegrationTest& test, const IntegrationResult& result) {
  if(result.momentComputed) {
    PrintRow(test.name, "Moment", "Omega+", result.omega1, test.exactOmega1);
    PrintRow(test.name, "Moment", "Gamma", result.gamma, test.exactGamma);
  }
  if(result.algoimComputed) {
    PrintRow(test.name, "Algoim", "Omega+", result.omega1Algoim, test.exactOmega1);
    PrintRow(test.name, "Algoim", "Gamma", result.gammaAlgoim, test.exactGamma);
  }
}

} // namespace

int main(int argc, char** argv) {
  FemusInit mpinit(argc, argv, MPI_COMM_WORLD);
  int iproc = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);

  DriverOptions options;
  try {
    options = Parse(argc, argv);
  } catch(const std::exception& error) {
    if(iproc == 0) std::cerr << "Argument error: " << error.what() << "\nUse --help for options.\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
  }

  const std::string psiName = "Psi";
  for(const GeomElType element : options.elements) {
    const unsigned dim = (element == TET || element == HEX) ? 3u : 2u;
    for(const std::string& kind : options.psis) {
      std::vector<double> center(dim, 0.);
      double radius = 0.25;
      if(element == TRI || element == TET) {
        const double rin = 1. / (dim + std::sqrt(static_cast<double>(dim)));
        center.assign(dim, -0.5 + rin);
        radius = 0.5 * rin;
      }

      PsiConfig cfg;
      cfg.type = PsiType::BALL;
      cfg.center = center;
      cfg.radius = radius;
      Psi psi(cfg);
      const bool analytic = options.curved;

      MultiLevelMesh mlMsh;
      const std::array<double, 3> min = {{-0.5, -0.5, -0.5}};
      const std::array<double, 3> max = {{ 0.5,  0.5,  0.5}};
      const std::array<double, 3> ballCenter = {{center[0], center[1], dim == 3 ? center[2] : 0.}};
      const AnalyticBallPsi analyticPsi = {element, options.curvature, ballCenter, radius};
      const bool nonconvex = kind == "nonconvex";
      const NonconvexPsi nonconvexPsi = {
          element, options.nonconvexQuadA, options.nonconvexQuadB, options.nonconvexQuadC,
          options.nonconvexHexA, options.nonconvexHexB, options.nonconvexHexC};

      const auto tests = nonconvex
          ? MakeNonconvexTests(element,
              options.nonconvexQuadA, options.nonconvexQuadB, options.nonconvexQuadC,
              options.nonconvexHexA, options.nonconvexHexB, options.nonconvexHexC)
          : analytic ? MakeAnalyticBallTests(dim, center, radius, options.curvature)
                     : MakeStraightBallTests(dim, center, radius);

      if(iproc == 0) {
        std::cout<<std::endl;
        std::cout << "\n" << std::string(90, '=') << "\n"
                  << "Element: " << ElementName(element) << " | Psi: " << PsiName(kind)
                  << " | Mesh: " << (analytic ? "analytic curved ball" : "straight");
        if(analytic) std::cout << " | kappa: " << options.curvature;
        std::cout << "\n" << std::string(90, '=') << '\n';
      }

      if(analytic) {
        mlMsh.GenerateSingleAnalyticBallElementMesh(
          element, min, max, ballCenter, options.curvature, "fifth");
      } else {
        mlMsh.GenerateSingleElementMesh(element, min, max, "fifth");
      }

      std::string meshTag = "straight";
      if(nonconvex) {
        const auto tagCoeff = [](double v) {
          std::ostringstream os;
          os << std::setprecision(8) << v;
          std::string t=os.str();
          std::replace(t.begin(),t.end(),'.','p');
          return t;
        };
        meshTag = "straight_";
        if(element==QUAD) meshTag += "a"+tagCoeff(options.nonconvexQuadA)
                                   +"_b"+tagCoeff(options.nonconvexQuadB)
                                   +"_c"+tagCoeff(options.nonconvexQuadC);
        if(element==HEX) meshTag += "a"+tagCoeff(options.nonconvexHexA)
                                  +"_b"+tagCoeff(options.nonconvexHexB)
                                  +"_c"+tagCoeff(options.nonconvexHexC);
      }
      if(analytic) {
        std::ostringstream tag;
        tag << "analytic_k" << std::fixed << std::setprecision(3) << options.curvature;
        meshTag = tag.str();
        std::replace(meshTag.begin(), meshTag.end(), '.', 'p');
      }

      std::ofstream csv;
      if(iproc == 0) {
        const std::string file = "results_" + std::string(ElementName(element)) + "_" + PsiName(kind)
                               + "_" + meshTag + ".csv";
        csv.open(file.c_str(), std::ios::trunc);
        if(!csv) throw std::runtime_error("Cannot create CSV: " + file);
        csv << "level,test,method,domain,numerical,exact,absolute_error,relative_error,cut_elements,algoim_volume_points,algoim_surface_points\n";
        if(nonconvex) {
          std::cout << "Nonconvex exact references (interior interface only): "
                    << "volume=" << std::setprecision(17) << tests[0].exactOmega1
                    << " gamma=" << tests[0].exactGamma << '\n';
          if(element==QUAD && (options.nonconvexQuadA == 4. || options.nonconvexQuadB == 4.))
            std::cout << "NOTE: QUAD interface touches the element boundary "
                         "at x=+/-0.5 or y=+/-0.5; "
                         "set quad-a,quad-b>4 for strict interior.\n";
        }
        if(analytic) {
          std::cout << "Physical references (closed form):"
                    << " V=" << std::setprecision(17) << tests[0].exactOmega1
                    << " S=" << tests[0].exactGamma
                    << " Vr2=" << tests[1].exactOmega1
                    << " Sr2=" << tests[1].exactGamma << '\n';
        }
      }

      // Initialize a new solution on every level: Integrate() always reads the finest level.
      for(unsigned level = 0; level <= options.levels; ++level) {
        if(iproc == 0) {
          std::cout << std::string(90, '-') << '\n';
          std::cout << "Level " << level << " (" << mlMsh.GetNumberOfLevels()
                    << " mesh level(s))\n\n";
        }
        if(level != 0) {
          // AMR marking must use the same level set used for initialization.
          if(analytic) FlagFinestMeshLevel(mlMsh, analyticPsi);
          else if(nonconvex) FlagFinestMeshLevel(mlMsh, nonconvexPsi);
          else FlagFinestMeshLevel(mlMsh, psi);
          mlMsh.AddAMRMeshLevel(false);
        }
        {
          MultiLevelSolution mlSol(&mlMsh);
          mlSol.AddSolution(psiName.c_str(), LAGRANGE, SECOND, 0, false);
          mlSol.Initialize("All");
          // On curved elements, evaluate the quadratic reference BALL using
          // the inverse analytic map at each physical DoF, including AMR levels.
          if(analytic) InitLevelSet(mlSol, psiName, analyticPsi);
          else if(nonconvex) InitLevelSet(mlSol, psiName, nonconvexPsi);
          else InitLevelSet(mlSol, psiName, psi);

          if(options.vtk) {
            const std::string outputdir = "output_" + std::string(ElementName(element));
            VTKWriter writer(&mlSol);
            writer.SetDebugOutput(true);
            writer.Write(outputdir, "biquadratic", std::vector<std::string>{psiName}, 0);
          }

          if(iproc == 0) {
            std::cout<<std::endl;
            PrintHeader();
          }
          for(const auto& test : tests) {
            const IntegrationResult result = Integrate(mlSol, psiName, test.f, options.integration);
            if(iproc == 0) {
              PrintResults(test, result);
              SaveResults(csv, level, test, result);
            }
          }
        }
      }
    }
  }
  return 0;
}

template <class PsiType>
void FlagFinestMeshLevel(MultiLevelMesh& mlMsh, const PsiType& psi) {

  const unsigned level = mlMsh.GetNumberOfLevels() - 1;
  Mesh& msh = *mlMsh.GetLevel(level);

  const unsigned iproc = msh.processor_id();
  const unsigned dim = msh.GetDimension();
  const unsigned xType = 2;
  const unsigned amrIndex = msh.GetAmrIndex();

  auto& xv = msh._topology->_Sol;
  auto* amrFlag = msh._topology->_Sol[amrIndex];

  amrFlag->zero();

  const unsigned offset = msh._elementOffset[iproc];
  const unsigned offsetp1 = msh._elementOffset[iproc + 1];

  std::vector<double> x(dim);

  for(unsigned iel = offset; iel < offsetp1; ++iel) {

    const unsigned nDof = msh.GetElementDofNumber(iel, xType);

    const unsigned xDof0 = msh.GetSolutionDof(0, iel, xType);

    for(unsigned k = 0; k < dim; ++k)
      x[k] = (*xv[k])(xDof0);

    const double psi0 = psi(x);
    const int sign0 = (psi0 > 0.) ? 1 : -1;

    bool signChange = false;

    for(unsigned i = 1; i < nDof; ++i) {

      const unsigned xDof = msh.GetSolutionDof(i, iel, xType);

      for(unsigned k = 0; k < dim; ++k)
        x[k] = (*xv[k])(xDof);

      const double psii = psi(x);
      const int signi = (psii > 0.) ? 1 : -1;

      if(signi != sign0) {
        signChange = true;
        break;
      }
    }

    if(signChange)
      amrFlag->set(iel, 2);
  }

  amrFlag->close();

  for(unsigned iel = offset; iel < offsetp1; ++iel) {

    if((*amrFlag)(iel) == 2) {

      for(unsigned j = 0; j < msh.el->GetElementNearElementSize(iel, 1); ++j) {

        const unsigned jel = msh.el->GetElementNearElement(iel, j);

        if(offset <= jel && jel < offsetp1) {

          if((*amrFlag)(jel) < 0.5)
            amrFlag->set(jel, 1);
        }
        else {
          amrFlag->add(jel, 1);
        }
      }
    }
  }

  amrFlag->close();

  unsigned localRefined = 0;

  for(unsigned iel = offset; iel < offsetp1; ++iel) {

    if((*amrFlag)(iel) > 0.5) {

      amrFlag->set(iel, 1);

      ++localRefined;
    }
  }

  amrFlag->close();

  unsigned globalRefined = 0;

  MPI_Allreduce(
    &localRefined,
    &globalRefined,
    1,
    MPI_UNSIGNED,
    MPI_SUM,
    MPI_COMM_WORLD);

  msh.el->SetRefinedElementNumber(globalRefined);
}

template <class PsiType>
void InitLevelSet(MultiLevelSolution & mlSol, const std::string & name,
                  const PsiType & psi2D) {

  MultiLevelMesh &mlMsh = *mlSol.GetMultilevelMesh();
  const unsigned level = mlMsh.GetNumberOfLevels() - 1u;
  Mesh &msh = *mlMsh.GetLevel(level);
  Solution &sol = *mlSol.GetLevel(level);
  const unsigned iproc = msh.processor_id();
  const unsigned dim = msh.GetDimension();

  unsigned solIndex = mlSol.GetIndex(name.c_str());
  unsigned solType = mlSol.GetSolutionType(name.c_str());

  const unsigned xType = 2u; // coordinate field type

  auto &xv = msh._topology->_Sol;
  std::vector<double> x(dim);

  const unsigned offset = msh._elementOffset[iproc];
  const unsigned offsetp1 = msh._elementOffset[iproc + 1];

  auto &solVec = sol._Sol[solIndex];

  solVec->zero();

  // Loop over local elements and interpolate psi2D at solution DoFs
  for (unsigned iel = offset; iel < offsetp1; ++iel) {

    const unsigned nDof = msh.GetElementDofNumber(iel, solType);

    for (unsigned i = 0; i < nDof; ++i) {

      // Get physical coordinates of the current DoF
      const unsigned xDof = msh.GetSolutionDof(i, iel, xType);
      for (unsigned k = 0; k < dim; ++k) {
        x[k] = (*xv[k])(xDof);
      }

      // Evaluate and assign field value
      const unsigned solDof = msh.GetSolutionDof(i, iel, solType);
      solVec->set(solDof, psi2D(x));
    }
  }

  solVec->close();
}