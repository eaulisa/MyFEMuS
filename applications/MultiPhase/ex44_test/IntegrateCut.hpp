#pragma once

#include "CutFemWeightNew.hpp"
#include "Fem.hpp"
#include "MultiLevelSolution.hpp"
#include "../includeLS/ElementTopology.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

typedef double TypeIO;
// typedef cpp_bin_float_oct TypeA;
typedef long double TypeA;
typedef cpp_bin_float_oct oct;

CutFemWeight <TypeIO, TypeA> hex  = CutFemWeight<TypeIO, TypeA >(HEX, 5, "legendre");
CutFemWeight <TypeIO, TypeA> tet  = CutFemWeight<TypeIO, TypeA >(TET, 5, "legendre");
CutFemWeight <TypeIO, TypeA> wedge = CutFemWeight<TypeIO, TypeA >(WEDGE, 5, "legendre");
CutFemWeight <TypeIO, TypeA> quad = CutFemWeight<TypeIO, TypeA >(QUAD, 5, "legendre");
CutFemWeight <TypeIO, TypeA> tri = CutFemWeight<TypeIO, TypeA >(TRI, 5, "legendre");
CutFemWeight <TypeIO, TypeA> line = CutFemWeight<TypeIO, TypeA >(LINE, 5, "legendre");

const std::vector< CutFemWeight <TypeIO, TypeA> *> cfw = {&hex, &tet, &wedge, &quad, &tri, &line};

Fem fem = Fem(hex.GetGaussQuadratureOrder(), hex.GetDimension());

struct IntegrationResult {
  double omega1 = 0.;
  double gamma = 0.;
  unsigned cutElements = 0;
};

struct IntegrationTest {
  std::string name;
  std::function<double(const std::vector<double>&)> f;
  double exactOmega1;
  double exactGamma;
};

template <typename Function>
IntegrationResult Integrate(femus::MultiLevelSolution& mlSol, const std::string& psiName, const Function& f);



void RunIntegrationTest(MultiLevelSolution& mlSol, const std::string& psiName, const IntegrationTest& test) {

  const IntegrationResult result = Integrate(mlSol, psiName, test.f);

  auto PrintError = [](const std::string& domain, const double numerical, const double exact) {

    const double absError = std::abs(numerical - exact);
    const double relError = (std::abs(exact) > 0.) ? absError / std::abs(exact) : 0.;

    std::cout << "  " << domain << std::endl;
    std::cout << "    Numerical      = " << std::scientific << std::setprecision(12) << numerical << std::endl;
    std::cout << "    Exact          = " << exact << std::endl;
    std::cout << "    Absolute error = " << absError << std::endl;
    std::cout << "    Relative error = " << relError << std::endl;
  };

  std::cout << "\n========== " << test.name << " ==========" << std::endl;
  std::cout << "Cut elements = " << result.cutElements << std::endl;

  PrintError("Omega1", result.omega1, test.exactOmega1);
  PrintError("Gamma", result.gamma, test.exactGamma);
}

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

  if(psi.size() != topology.nNodes) throw std::runtime_error("IsCutElementQuadratic: wrong nodal count");
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
    if(edge.nodes.n != 3 || edge.vertices.n != 2) throw std::runtime_error("IsCutElementQuadratic: expected 3-node quadratic edges");
    if(EdgeIsCut(psi[edge.vertices[0]], psi[edge.nodes[2]], psi[edge.vertices[1]])) return true;
  }

  return false;
}


inline bool GetAffineSurfaceInverseJacobianTest(const femus::elem_type* fe, const std::vector<std::vector<double>>& coordinates, std::vector<std::vector<double>>& inverseJacobian) {
  const unsigned dim = fe->GetDim();
  const unsigned nodes = static_cast<unsigned>(fe->GetNDofs());
  inverseJacobian.clear();
  if(coordinates.size() != dim || nodes == 0) throw std::invalid_argument("GetAffineSurfaceInverseJacobianTest: invalid coordinates");
  for(unsigned d = 0; d < dim; ++d) if(coordinates[d].size() != nodes) throw std::invalid_argument("GetAffineSurfaceInverseJacobianTest: invalid nodes");
  const double* origin = fe->GetBasis()->GetXcoarse(0);
  if(origin == nullptr) return false;
  const std::vector<double> referenceOrigin(origin, origin + dim);
  std::vector<std::vector<double>> relative = coordinates;
  long double scale = 0.;
  for(unsigned d = 0; d < dim; ++d) {
    for(unsigned i = 0; i < nodes; ++i) {
      relative[d][i] -= coordinates[d][0];
      scale = std::max(scale, std::abs(static_cast<long double>(relative[d][i])));
    }
  }
  double determinant = 0.;
  std::vector<std::vector<double>> jacobian;
  fe->GetJacobianMatrix(relative, referenceOrigin, determinant, jacobian, inverseJacobian);
  if(!std::isfinite(determinant) || determinant == 0.) throw std::domain_error("GetAffineSurfaceInverseJacobianTest: singular geometry");
  const long double tolerance = 128.L * std::numeric_limits<double>::epsilon() * scale;
  for(unsigned i = 0; i < nodes; ++i) {
    const double* point = fe->GetBasis()->GetXcoarse(static_cast<int>(i));
    if(point == nullptr) { inverseJacobian.clear(); return false; }
    for(unsigned d = 0; d < dim; ++d) {
      long double displacement = 0.;
      for(unsigned j = 0; j < dim; ++j) displacement += static_cast<long double>(jacobian[d][j]) * (static_cast<long double>(point[j]) - origin[j]);
      if(std::abs(displacement - static_cast<long double>(relative[d][i])) > tolerance) { inverseJacobian.clear(); return false; }
    }
  }
  return true;
}

template <typename Function>
IntegrationResult Integrate(femus::MultiLevelSolution& mlSol, const std::string& psiName, const Function& f) {

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

  IntegrationResult result;

  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; ++iel) {

    const unsigned geom = msh->GetElementType(iel);
    const unsigned nDofs = msh->GetElementDofNumber(iel, xType);

    // CutFemWeight<double, cpp_bin_float_oct> cut(static_cast<GeomElType>(geom), 5, "legendre");
    // femus::Fem fem(cfw[geom]->GetGaussQuadratureOrder(), dim);
    const femus::elem_type* fe = fem.GetFiniteElement(geom, xType);
    const unsigned nGauss = fe->GetGaussPointNumber();

    if(nGauss != cfw[geom]->GetGaussQuadraturePointNumber()) throw std::runtime_error("Integrate: quadrature mismatch");

    std::vector<double> psi(nDofs);
    std::vector<std::vector<double>> coordX(dim, std::vector<double>(nDofs));
    for(unsigned i = 0; i < nDofs; ++i) {
      const unsigned dof = msh->GetSolutionDof(i, iel, xType);
      psi[i] = (*sol->_Sol[psiIndex])(dof);
      for(unsigned d = 0; d < dim; ++d) coordX[d][i] = (*msh->_topology->_Sol[d])(dof);
    }

    const bool isCut = IsCutElementQuadratic(psi, geom);
    std::vector<double> weightOmega(nGauss, 0.);
    std::vector<double> weightGamma(nGauss, 0.);
    bool physicalWeights = false;

    if(isCut) {

      ++result.cutElements;

      cfw[geom]->SetHex27IntegrationDirection(0);
      cfw[geom]->SetTet15IntegrationDirection(0);
      cfw[geom]->SetQuad9IntegrationDirection(0);
      cfw[geom]->SetTri7IntegrationDirection(0);
      cfw[geom]->SetHex27DirectionValidation(0);
      cfw[geom]->SetTet15DirectionValidation(0);
      cfw[geom]->SetQuad9DirectionValidation(0);
      cfw[geom]->SetTri7DirectionValidation(0);

      std::vector<std::vector<double>> surfaceInverseJacobian;
      physicalWeights = GetAffineSurfaceInverseJacobianTest(fe, coordX, surfaceInverseJacobian);

      if(physicalWeights) cfw[geom]->GetPhysicalSurfaceWeightsFromNodalValues(psi, surfaceInverseJacobian, weightOmega, weightGamma, 1.e-5, 1.e-5);
      else cfw[geom]->GetWeightsFromNodalValues(psi, weightOmega, weightGamma, 1.e-5, 1.e-5, true);

    }
    else if(psi[0] > 0.) {

      std::fill(weightOmega.begin(), weightOmega.end(), 1.);

    }

    std::vector<double> phi, phi_x, dphi;
    std::vector<double> x(dim, 0.), xi(dim, 0.), gradRef(dim, 0.);
    std::vector<std::vector<double>> J, Ji;

    for(unsigned ig = 0; ig < nGauss; ++ig) {

      double weight = 0.;
      fe->Jacobian(coordX, ig, weight, phi, phi_x);
      std::fill(x.begin(), x.end(), 0.);

      for(unsigned i = 0; i < nDofs; ++i) {
        for(unsigned d = 0; d < dim; ++d) x[d] += coordX[d][i] * phi[i];
      }

      const double value = f(x);
      result.omega1 += value * weight * weightOmega[ig];

      if(!isCut || weightGamma[ig] == 0.) continue;

      double dsN = 1.;
      if(!physicalWeights) {

        for(unsigned d = 0; d < dim; ++d) xi[d] = fe->GetGaussRule().GetGaussCoordinatePointer(d)[ig];

        std::fill(gradRef.begin(), gradRef.end(), 0.);
        fe->GetDPhiDXi(dphi, xi);
        for(unsigned i = 0; i < nDofs; ++i) gradRef[0] += psi[i] * dphi[i];
        fe->GetDPhiDEta(dphi, xi);
        for(unsigned i = 0; i < nDofs; ++i) gradRef[1] += psi[i] * dphi[i];
        if(dim == 3) {
          fe->GetDPhiDZeta(dphi, xi);
          for(unsigned i = 0; i < nDofs; ++i) gradRef[2] += psi[i] * dphi[i];
        }

        double normRef = 0.;
        for(unsigned d = 0; d < dim; ++d) normRef = std::hypot(normRef, gradRef[d]);
        if(normRef <= 1.e-14) throw std::runtime_error("Integrate: zero reference gradient on a nonaffine mapped element");

        double weightFake = 0.;
        fe->GetJacobianMatrix(coordX, ig, weightFake, J, Ji);
        double normPhys = 0.;
        for(unsigned d = 0; d < dim; ++d) {
          double g = 0.;
          for(unsigned j = 0; j < dim; ++j) g += Ji[j][d] * gradRef[j];
          normPhys = std::hypot(normPhys, g);
        }
        dsN = normPhys / normRef;

      }

      result.gamma += value * weight * weightGamma[ig] * dsN;

    }
  }

  double integralTotals[2] = {result.omega1, result.gamma};
  unsigned totalCut = result.cutElements;
  MPI_Allreduce(MPI_IN_PLACE, integralTotals, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &totalCut, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  
  result.omega1 = integralTotals[0];
  result.gamma = integralTotals[1];
  result.cutElements = totalCut;

  return result;

}
