#pragma once

// Integration of the positive level-set phase and its zero interface.
// Default: in-house moment fitting only.
// Optional: direct Algoim quadrature only, or a side-by-side comparison.
// Based on MyFEMuSTTU/rifat b6932651b2bf42e09e3f36d96e0676485d57fab7.
// The direct Algoim adapter requires the ex44_test target to compile as C++17.

#include "CutFemWeightLS.hpp"
#include "AlgoimDirectRule.hpp"
#include "Fem.hpp"
#include "MultiLevelSolution.hpp"
#include "../includeLS/ElementTopology.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// These match the precision choices in applications/MultiPhase/includeLS/Utils.hpp.
// AlgoimDirectRule.hpp currently returns double-valued nodes and weights.
typedef double TypeIO;
// typedef long double TypeA;
typedef cpp_bin_float_oct TypeA;

enum class IntegrationMode {
  MomentFitting,
  Algoim,
  Compare
};

struct IntegrationOptions {
  IntegrationMode mode = IntegrationMode::MomentFitting;
  unsigned momentDegree = 10;
  double absoluteTolerance = 1.e-6;
  double relativeTolerance = 1.e-6;
  unsigned algoimOrder = 5;
  unsigned algoimSubdivisions = 0;
};

struct IntegrationResult {
  double omega1 = 0.;
  double gamma = 0.;
  double omega1Algoim = 0.;
  double gammaAlgoim = 0.;
  unsigned cutElements = 0;
  unsigned algoimVolumePoints = 0;
  unsigned algoimSurfacePoints = 0;
  bool momentComputed = false;
  bool algoimComputed = false;
};

struct IntegrationTest {
  std::string name;
  std::function<double(const std::vector<double>&)> f;
  double exactOmega1;
  double exactGamma;
};

inline bool EdgeIsCut(const double p0, const double pm, const double p1) {
  const double scale = std::max({std::abs(p0), std::abs(pm), std::abs(p1)});
  if(scale == 0.) return true;

  const double a = 2. * (p0 + p1 - 2. * pm) / scale;
  const double b = (4. * pm - 3. * p0 - p1) / scale;
  const double c = p0 / scale;
  const double tolerance = 1.e-12;
  const auto inside = [&](const double t) { return t > tolerance && t < 1. - tolerance; };

  if(std::abs(a) <= 1.e-14) return std::abs(b) > 1.e-14 && inside(-c / b);
  const double discriminant = b * b - 4. * a * c;
  if(discriminant < -1.e-13) return false;
  const double delta = std::sqrt(std::max(0., discriminant));
  return inside((-b - delta) / (2. * a)) || inside((-b + delta) / (2. * a));
}

inline bool IsCutElementQuadratic(const std::vector<double>& psi, const unsigned geom) {
  const auto& topology = ex40_topo::topologyFromTypeId(geom);
  if(psi.size() != topology.nNodes) throw std::runtime_error("IsCutElementQuadratic: incorrect nodal count");

  double maximum = 0.;
  for(const double value : psi) {
    if(!std::isfinite(value)) throw std::runtime_error("IsCutElementQuadratic: non-finite Psi");
    maximum = std::max(maximum, std::abs(value));
  }
  if(maximum == 0.) throw std::runtime_error("IsCutElementQuadratic: identically zero Psi");

  const double tolerance = 1.e-12 * maximum;
  bool positive = false;
  bool negative = false;
  for(const double value : psi) {
    if(std::abs(value) <= tolerance) return true;
    if(value > 0.) positive = true;
    else negative = true;
  }
  if(positive && negative) return true;

  for(unsigned e = 0; e < topology.nEdges; ++e) {
    const auto& edge = topology.edges[e];
    if(edge.nodes.n != 3 || edge.vertices.n != 2) throw std::runtime_error("IsCutElementQuadratic: expected a quadratic edge");
    if(EdgeIsCut(psi[edge.vertices[0]], psi[edge.nodes[2]], psi[edge.vertices[1]])) return true;
  }
  return false;
}

namespace ex44_integration_detail {

// Cache each geometry/degree rule to reuse its Gram-Schmidt matrices.
inline CutFemWeightLS<TypeIO, TypeA>& GetMomentRule(const unsigned geom, const unsigned degree) {
  static std::map<std::pair<unsigned, unsigned>, std::unique_ptr<CutFemWeightLS<TypeIO, TypeA>>> cache;
  const std::pair<unsigned, unsigned> key(geom, degree);
  std::unique_ptr<CutFemWeightLS<TypeIO, TypeA>>& rule = cache[key];
  if(!rule) {
    rule.reset(new CutFemWeightLS<TypeIO, TypeA>(static_cast<GeomElType>(geom), degree, "legendre"));
    rule->SetHex27IntegrationDirection(0);
    rule->SetTet15IntegrationDirection(0);
    rule->SetQuad9IntegrationDirection(0);
    rule->SetTri7IntegrationDirection(0);
    rule->SetHex27DirectionValidation(0);
    rule->SetTet15DirectionValidation(0);
    rule->SetQuad9DirectionValidation(0);
    rule->SetTri7DirectionValidation(0);
  }
  return *rule;
}

inline bool RunMoment(const IntegrationMode mode) {
  return mode == IntegrationMode::MomentFitting || mode == IntegrationMode::Compare;
}

inline bool RunAlgoim(const IntegrationMode mode) {
  return mode == IntegrationMode::Algoim || mode == IntegrationMode::Compare;
}

inline void PrintError(const char* domain, const double numerical, const double exact) {
  const double absoluteError = std::abs(numerical - exact);
  std::cout << "  " << domain << std::endl;
  std::cout << "    Numerical      = " << numerical << std::endl;
  std::cout << "    Exact          = " << exact << std::endl;
  std::cout << "    Absolute error = " << absoluteError << std::endl;
  if(exact == 0.) std::cout << "    Relative error = N/A (exact integral is zero)" << std::endl;
  else std::cout << "    Relative error = " << absoluteError / std::abs(exact) << std::endl;
}

inline void EvaluatePhysicalCoordinates(const femus::elem_type* fe, const std::vector<std::vector<double>>& coordX, const std::vector<double>& xi, std::vector<double>& phi, std::vector<double>& x) {
  fe->GetPhi(phi, xi);
  std::fill(x.begin(), x.end(), 0.);
  for(unsigned d = 0; d < x.size(); ++d) for(unsigned i = 0; i < phi.size(); ++i) x[d] += coordX[d][i] * phi[i];
}

inline double PhysicalLevelSetGradientNorm(const femus::elem_type* fe, const std::vector<double>& psi, const std::vector<double>& xi, const std::vector<std::vector<double>>& inverseJacobian, std::vector<double>& dphi, std::vector<double>& gradRef) {
  std::fill(gradRef.begin(), gradRef.end(), 0.);
  fe->GetDPhiDXi(dphi, xi);
  for(unsigned i = 0; i < psi.size(); ++i) gradRef[0] += psi[i] * dphi[i];
  if(gradRef.size() >= 2) {
    fe->GetDPhiDEta(dphi, xi);
    for(unsigned i = 0; i < psi.size(); ++i) gradRef[1] += psi[i] * dphi[i];
  }
  if(gradRef.size() == 3) {
    fe->GetDPhiDZeta(dphi, xi);
    for(unsigned i = 0; i < psi.size(); ++i) gradRef[2] += psi[i] * dphi[i];
  }
  double norm = 0.;
  for(unsigned d = 0; d < gradRef.size(); ++d) {
    double component = 0.;
    for(unsigned j = 0; j < gradRef.size(); ++j) component += inverseJacobian[j][d] * gradRef[j];
    norm = std::hypot(norm, component);
  }
  return norm;
}

} // namespace ex44_integration_detail

template <typename Function>
IntegrationResult Integrate(femus::MultiLevelSolution& mlSol, const std::string& psiName, const Function& f, const IntegrationOptions& options = IntegrationOptions()) {
  const bool runMoment = ex44_integration_detail::RunMoment(options.mode);
  const bool runAlgoim = ex44_integration_detail::RunAlgoim(options.mode);
  if(!runMoment && !runAlgoim) throw std::invalid_argument("Integrate: invalid integration mode");
  if(options.momentDegree == 0) throw std::invalid_argument("Integrate: momentDegree must be positive");
  if(runMoment && (!(options.absoluteTolerance > 0.) || !(options.relativeTolerance >= 0.) || !std::isfinite(options.absoluteTolerance) || !std::isfinite(options.relativeTolerance))) throw std::invalid_argument("Integrate: invalid moment-fitting tolerances");
  if(runAlgoim && (options.algoimOrder == 0 || options.algoimOrder > 100 || options.algoimSubdivisions > 32)) throw std::invalid_argument("Integrate: invalid Algoim parameters");

  femus::MultiLevelMesh& mlMsh = *mlSol.GetMultilevelMesh();
  const unsigned level = mlMsh.GetNumberOfLevels() - 1;
  femus::Mesh* msh = mlMsh.GetLevel(level);
  femus::Solution* sol = mlSol.GetLevel(level);
  const unsigned dim = msh->GetDimension();
  const unsigned iproc = msh->processor_id();
  const unsigned psiIndex = mlSol.GetIndex(psiName.c_str());
  const unsigned psiType = mlSol.GetSolutionType(psiIndex);
  const unsigned xType = 2;
  if(psiType != xType) throw std::runtime_error("Integrate: Psi must use LAGRANGE SECOND");
  if(dim < 2 || dim > 3) throw std::invalid_argument("Integrate: only 2D and 3D are supported");

  IntegrationResult result;
  result.momentComputed = runMoment;
  result.algoimComputed = runAlgoim;

  // No eager global quadrature initialization. FEMuS FE tables are local;
  // moment-fitting transforms are lazily cached across integration tests.
  femus::Fem fem(2U * options.momentDegree, dim);

  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; ++iel) {
    const unsigned geom = msh->GetElementType(iel);
    const unsigned nDofs = msh->GetElementDofNumber(iel, xType);
    const femus::elem_type* fe = fem.GetFiniteElement(geom, xType);
    const unsigned nGauss = fe->GetGaussPointNumber();

    std::vector<double> psi(nDofs);
    std::vector<std::vector<double>> coordX(dim, std::vector<double>(nDofs));
    for(unsigned i = 0; i < nDofs; ++i) {
      const unsigned dof = msh->GetSolutionDof(i, iel, xType);
      psi[i] = (*sol->_Sol[psiIndex])(dof);
      for(unsigned d = 0; d < dim; ++d) coordX[d][i] = (*msh->_topology->_Sol[d])(dof);
    }

    const bool isCut = IsCutElementQuadratic(psi, geom);
    if(isCut) ++result.cutElements;
    if(!isCut && psi[0] <= 0.) continue;

    std::vector<double> weightOmega;
    std::vector<double> weightGamma;

    if(isCut && runMoment) {
      CutFemWeightLS<TypeIO, TypeA>& rule = ex44_integration_detail::GetMomentRule(geom, options.momentDegree);
      if(nGauss != rule.GetGaussQuadraturePointNumber()) throw std::runtime_error("Integrate: moment-fitting Gauss rules do not match FEMuS");
      // The geometry-corrected surface weights already include
      // |cof(J) * grad_ref(Psi)|. The Gauss reference weight is applied below.
      rule.GetPhysicalSurfaceWeightsFromNodalGeometry(psi, coordX, weightOmega, weightGamma, options.absoluteTolerance, options.relativeTolerance);
    }

    // Direct Algoim nodes are evaluated independently; no FEMuS moment fitting.
    if(isCut && runAlgoim) {
      const FemusAlgoimDirect::Rule rule = FemusAlgoimDirect::GetNodalQuadrature(psi, options.algoimOrder, options.algoimSubdivisions);
      result.algoimVolumePoints += static_cast<unsigned>(rule.volume.size());
      result.algoimSurfacePoints += static_cast<unsigned>(rule.surfaceDelta.size());

      std::vector<std::vector<double>> J, Ji;
      std::vector<double> phi;
      std::vector<double> dphi;
      std::vector<double> x(dim, 0.);
      std::vector<double> gradRef(dim, 0.);

      for(const auto& qp : rule.volume) {
        double det = 0.;
        fe->GetJacobianMatrix(coordX, qp.xi, det, J, Ji);
        ex44_integration_detail::EvaluatePhysicalCoordinates(fe, coordX, qp.xi, phi, x);
        result.omega1Algoim += f(x) * std::abs(det) * qp.weight;
      }

      for(const auto& qp : rule.surfaceDelta) {
        double det = 0.;
        fe->GetJacobianMatrix(coordX, qp.xi, det, J, Ji);
        ex44_integration_detail::EvaluatePhysicalCoordinates(fe, coordX, qp.xi, phi, x);
        const double gradNorm = ex44_integration_detail::PhysicalLevelSetGradientNorm(fe, psi, qp.xi, Ji, dphi, gradRef);
        if(!(gradNorm > 0.)) throw std::domain_error("Integrate: zero physical gradient at an Algoim surface point");
        result.gammaAlgoim += f(x) * std::abs(det) * gradNorm * qp.weight;
      }
    }

    // FEMuS Gauss loop is needed for moment fitting, or for a completely
    // positive element (where the volume integral is the same for both modes).
    if(runMoment || !isCut) {
      std::vector<double> phi, phi_x;
      std::vector<double> x(dim, 0.);
      for(unsigned ig = 0; ig < nGauss; ++ig) {
        double physicalWeight = 0.;
        fe->Jacobian(coordX, ig, physicalWeight, phi, phi_x);
        std::fill(x.begin(), x.end(), 0.);
        for(unsigned d = 0; d < dim; ++d) for(unsigned i = 0; i < nDofs; ++i) x[d] += coordX[d][i] * phi[i];
        const double value = f(x);

        if(!isCut) {
          if(runMoment) result.omega1 += value * physicalWeight;
          if(runAlgoim) result.omega1Algoim += value * physicalWeight;
        }
        else {
          result.omega1 += value * physicalWeight * weightOmega[ig];
          if(weightGamma[ig] != 0.) result.gamma += value * fe->GetGaussWeight(ig) * weightGamma[ig];
        }
      }
    }
  }

  // Collective reductions: all ranks run the selected paths, but only rank 0
  // prints. Inactive method slots remain exactly zero.
  double integrals[4] = {result.omega1, result.gamma, result.omega1Algoim, result.gammaAlgoim};
  MPI_Allreduce(MPI_IN_PLACE, integrals, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  result.omega1 = integrals[0];
  result.gamma = integrals[1];
  result.omega1Algoim = integrals[2];
  result.gammaAlgoim = integrals[3];

  unsigned counts[3] = {result.cutElements, result.algoimVolumePoints, result.algoimSurfacePoints};
  MPI_Allreduce(MPI_IN_PLACE, counts, 3, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  result.cutElements = counts[0];
  result.algoimVolumePoints = counts[1];
  result.algoimSurfacePoints = counts[2];
  return result;
}

// Convenience overload: true means moment fitting + direct Algoim.
template <typename Function>
IntegrationResult Integrate(femus::MultiLevelSolution& mlSol, const std::string& psiName, const Function& f, const bool compareWithAlgoim) {
  IntegrationOptions options;
  options.mode = compareWithAlgoim ? IntegrationMode::Compare : IntegrationMode::MomentFitting;
  return Integrate(mlSol, psiName, f, options);
}

inline void RunIntegrationTest(femus::MultiLevelSolution& mlSol, const std::string& psiName, const IntegrationTest& test, const IntegrationOptions& options = IntegrationOptions()) {
  const IntegrationResult result = Integrate(mlSol, psiName, test.f, options);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank != 0) return;

  const std::streamsize oldPrecision = std::cout.precision();
  const std::ios::fmtflags oldFlags = std::cout.flags();
  std::cout << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10);
  std::cout << "\n========== " << test.name << " ==========" << std::endl;
  std::cout << "Cut elements = " << result.cutElements << std::endl;

  if(result.momentComputed) {
    std::cout << "  === Moment fitting ===" << std::endl;
    ex44_integration_detail::PrintError("Omega1", result.omega1, test.exactOmega1);
    ex44_integration_detail::PrintError("Gamma", result.gamma, test.exactGamma);
  }
  if(result.algoimComputed) {
    std::cout << "  === Direct Algoim ===" << std::endl;
    ex44_integration_detail::PrintError("Omega1", result.omega1Algoim, test.exactOmega1);
    ex44_integration_detail::PrintError("Gamma", result.gammaAlgoim, test.exactGamma);
    std::cout << "  Algoim volume points  = " << result.algoimVolumePoints << std::endl;
    std::cout << "  Algoim surface points = " << result.algoimSurfacePoints << std::endl;
  }
  if(result.momentComputed && result.algoimComputed) {
    std::cout << "  === Algoim - moment fitting ===" << std::endl;
    std::cout << "  Omega1 difference = " << result.omega1Algoim - result.omega1 << std::endl;
    std::cout << "  Gamma difference  = " << result.gammaAlgoim - result.gamma << std::endl;
  }
  std::cout.flags(oldFlags);
  std::cout.precision(oldPrecision);
}

// Convenience overload: RunIntegrationTest(..., true) runs both methods.
inline void RunIntegrationTest(femus::MultiLevelSolution& mlSol, const std::string& psiName, const IntegrationTest& test, const bool compareWithAlgoim) {
  IntegrationOptions options;
  options.mode = compareWithAlgoim ? IntegrationMode::Compare : IntegrationMode::MomentFitting;
  RunIntegrationTest(mlSol, psiName, test, options);
}

// Explicit direct Algoim order/subdivisions: compare both methods.
// Keeps the calling convention used by the previous comparison header.
inline void RunIntegrationTest(femus::MultiLevelSolution& mlSol, const std::string& psiName, const IntegrationTest& test, const unsigned algoimOrder, const unsigned algoimSubdivisions) {
  IntegrationOptions options;
  options.mode = IntegrationMode::Compare;
  options.algoimOrder = algoimOrder;
  options.algoimSubdivisions = algoimSubdivisions;
  RunIntegrationTest(mlSol, psiName, test, options);
}

inline void RunIntegrationTest(femus::MultiLevelSolution& mlSol, const std::string& psiName, const IntegrationTest& test, const bool compareWithAlgoim, const unsigned algoimOrder, const unsigned algoimSubdivisions) {
  IntegrationOptions options;
  options.mode = compareWithAlgoim ? IntegrationMode::Compare : IntegrationMode::MomentFitting;
  options.algoimOrder = algoimOrder;
  options.algoimSubdivisions = algoimSubdivisions;
  RunIntegrationTest(mlSol, psiName, test, options);
}
