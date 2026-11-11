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

// Straight QUAD9-only Q2 oval.  Positive phase is the interior.
// Psi = 1 - a*x^2 - b*y^2 - c*x^2*y^2.
// The x^2*y^2 monomial belongs to Q2([-1,1]^2), but not P2(TRI).
struct StellaPsi {
  double a, b, c;

  double operator()(const std::vector<double>& x) const {
    if(x.size() != 2u)
      throw std::invalid_argument("StellaPsi requires a 2D point");
    const double x2 = x[0] * x[0];
    const double y2 = x[1] * x[1];
    return 1. - a * x2 - b * y2 - c * x2 * y2;
  }
};

// Gauss-Legendre integration of a smooth, independent parametrization of Gamma.
// This is only a reference calculation: the CutFEM integrator is not used.
struct StellaGaussLegendre {
  std::vector<long double> nodes, weights;

  explicit StellaGaussLegendre(unsigned n) : nodes(n), weights(n) {
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

void ValidateStella(double a, double b, double c) {
  if(!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c)
     || a <= 0. || b <= 0. || c < 0.)
    throw std::invalid_argument("Stella requires finite a>0, b>0, c>=0");

  // Our single straight QUAD9 occupies [-0.5,0.5]^2, not [-1,1]^2.
  // If this condition fails, the full-oval formula is NOT the cut-cell reference.
  if(a < 4. || b < 4.)
    throw std::invalid_argument(
        "Stella full-oval references require a>=4 and b>=4 "
        "in the FEMuS element [-0.5,0.5]^2");
}

// Closed form A = pi/sqrt(ab) * 2F1(1/2,1/2;2;-q), q=c/(ab).
// For q>0 use the equivalent real-modulus elliptic-integral formula:
// A = 4*sqrt(1+q)/(q*sqrt(ab)) * [K(sqrt(m))-E(sqrt(m))], m=q/(1+q).
// For small q, use the convergent hypergeometric Taylor series instead
// to avoid subtractive cancellation in K-E.
double StellaExactArea(double a, double b, double c) {
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
long double StellaPerimeterAtOrder(double a, double b, double c, unsigned order) {
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

  StellaGaussLegendre rule(order);
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

double StellaReferencePerimeter(double a, double b, double c) {
  const long double p128 = StellaPerimeterAtOrder(a, b, c, 128);
  const long double p256 = StellaPerimeterAtOrder(a, b, c, 256);
  if(std::abs(p256-p128) > 3.e-13L * std::abs(p256))
    throw std::runtime_error("Stella perimeter reference has not converged");
  return static_cast<double>(p256);
}

std::vector<IntegrationTest> MakeStellaTests(double a, double b, double c) {
  ValidateStella(a, b, c);
  const double area = StellaExactArea(a, b, c);
  const double perimeter = StellaReferencePerimeter(a, b, c);
  return {{"1", [](const std::vector<double>&) { return 1.; }, area, perimeter}};
}

struct DriverOptions {
  std::vector<GeomElType> elements = {TET};
  std::vector<std::string> psis = {"ball"};
  bool curved = false;
  unsigned levels = 3; // Number of AMR refinements; also evaluate level 0.
  double curvature = 0.40; // Analytical quadratic deformation strength.
  double stellaA = 4.0, stellaB = 4.0, stellaC = 160.0;
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
  if(value == "ball" || value == "saye" || value == "stella") return value;
  throw std::invalid_argument("Unknown psi type: " + value + " (expected ball,saye,stella)");
}

const std::string& PsiName(const std::string& value) { return value; }

void Usage() {
  std::cout << FELevelSetMoments::GaussKronrodRules::Help();
  std::cout << "Usage: MultiPhase_ex44_test [options]\n"
            << "  --elements tri,quad,tet,hex  Element types (comma-separated)\n"
            << "  --psi ball,saye,stella       Level-set types (comma-separated)\n"
            << "  --curved 0|1                 0=straight, 1=analytic curved BALL (default 0)\n"
            << "  --levels N                   Number of AMR refinements, plus level 0 (default 3)\n"
            << "  --curvature X                Analytical BALL deformation (default 0.40)\n"
            << "  --stella-a X                 Stella x^2 coefficient (default 4)\n"
            << "  --stella-b X                 Stella y^2 coefficient (default 4)\n"
            << "  --stella-c X                 Stella x^2*y^2 coefficient (default 160)\n"
            << "  --mode moment|algoim|compare Integration method (default moment)\n"
            << "  --degree N                   Moment-fitting degree (default 7)\n"
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
    else if(FELevelSetMoments::GaussKronrodRules::IsOption(key)) options.integration.outerQuadratureRules.SetOption(key,value);
    else if(key == "--curved") options.curved = Boolean(value);
    else if(key == "--levels") options.levels = Unsigned(value, key);
    else if(key == "--curvature") options.curvature = Real(value, key);
    else if(key == "--stella-a") options.stellaA = Real(value, key);
    else if(key == "--stella-b") options.stellaB = Real(value, key);
    else if(key == "--stella-c") options.stellaC = Real(value, key);
    else if(key == "--degree") options.integration.momentDegree = Unsigned(value, key);
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
  if(options.curvature < 0.) throw std::invalid_argument("--curvature must be nonnegative");
  if(options.curved && std::any_of(options.psis.begin(), options.psis.end(),
                                  [](const std::string& kind) { return kind != "ball"; }))
    throw std::invalid_argument("--curved 1 requires --psi ball (analytic geometry)");
  if(std::find(options.psis.begin(), options.psis.end(), "stella") != options.psis.end()) {
    ValidateStella(options.stellaA, options.stellaB, options.stellaC);
    if(std::find(options.elements.begin(), options.elements.end(), QUAD) == options.elements.end())
      throw std::invalid_argument("--psi stella requires QUAD in --elements");
  }
  return options;
}

std::vector<IntegrationTest> MakeTests(const PsiType kind, const unsigned dim,
                                       const std::vector<double>& center, const double radius) {
  const double pi = std::acos(-1.);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const bool ball = kind == PsiType::BALL;
  return {
    {"1", [](const std::vector<double>&) { return 1.; },
      ball ? (dim == 2 ? pi * radius * radius : 4. * pi * std::pow(radius, 3) / 3.) : nan,
      ball ? (dim == 2 ? 2. * pi * radius : 4. * pi * radius * radius) : nan},
    {"r^4", [center](const std::vector<double>& x) {
        double r2 = 0.;
        for(unsigned d = 0; d < x.size(); ++d) r2 += (x[d] - center[d]) * (x[d] - center[d]);
        return r2 * r2;
      },
      ball ? (dim == 2 ? pi * std::pow(radius, 6) / 3. : 4. * pi * std::pow(radius, 7) / 7.) : nan,
      ball ? (dim == 2 ? 2. * pi * std::pow(radius, 5) : 4. * pi * std::pow(radius, 6)) : nan}
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
      if(kind == "stella" && element != QUAD) continue;
      std::vector<double> center(dim, 0.);
      double radius = 0.25;
      if(element == TRI || element == TET) {
        const double rin = 1. / (dim + std::sqrt(static_cast<double>(dim)));
        center.assign(dim, -0.5 + rin);
        radius = 0.5 * rin;
      }

      PsiConfig cfg;
      cfg.type = kind == "saye" ? PsiType::SAYE : PsiType::BALL;
      cfg.center = center;
      cfg.radius = radius;
      cfg.amplitude = 0.1;
      cfg.waveNumber = 2. * std::acos(-1.);
      Psi psi(cfg);
      const bool analytic = options.curved;

      MultiLevelMesh mlMsh;
      const std::array<double, 3> min = {{-0.5, -0.5, -0.5}};
      const std::array<double, 3> max = {{ 0.5,  0.5,  0.5}};
      const std::array<double, 3> ballCenter = {{center[0], center[1], dim == 3 ? center[2] : 0.}};
      const AnalyticBallPsi analyticPsi = {element, options.curvature, ballCenter, radius};
      const StellaPsi stellaPsi = {options.stellaA, options.stellaB, options.stellaC};
      const bool stella = kind == "stella";

      const auto tests = stella
          ? MakeStellaTests(options.stellaA, options.stellaB, options.stellaC)
          : analytic ? MakeAnalyticBallTests(dim, center, radius, options.curvature)
                     : MakeTests(cfg.type, dim, center, radius);

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
      if(stella) {
        const auto tagCoeff = [](double coefficient) {
          std::ostringstream out;
          out << std::setprecision(8) << coefficient;
          std::string tag = out.str();
          std::replace(tag.begin(), tag.end(), '.', 'p');
          return tag;
        };
        meshTag = "straight_a" + tagCoeff(options.stellaA)
                + "_b" + tagCoeff(options.stellaB)
                + "_c" + tagCoeff(options.stellaC);
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
        if(stella) {
          std::cout << "Stella Q2 references: a=" << options.stellaA
                    << " b=" << options.stellaB << " c=" << options.stellaC
                    << " area=" << std::setprecision(17) << tests[0].exactOmega1
                    << " perimeter=" << tests[0].exactGamma << '\n';
          if(options.stellaA == 4. || options.stellaB == 4.)
            std::cout << "NOTE: Stella interface is tangent to the QUAD boundary "
                         "at x=+/-0.5 or y=+/-0.5; "
                         "use a,b>4 for strictly interior interface.\n";
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
          else if(stella) FlagFinestMeshLevel(mlMsh, stellaPsi);
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
          else if(stella) InitLevelSet(mlSol, psiName, stellaPsi);
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