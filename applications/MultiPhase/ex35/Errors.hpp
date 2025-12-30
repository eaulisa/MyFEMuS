#pragma once
#include <array>
#include <memory>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <cstddef>

#include "Field.hpp"
#include "FemProjection.hpp"

// phase 1 = Psi > 0
static inline double Hpos(const double psi) { return (psi > 0.0) ? 1.0 : 0.0; }

/**
 * Mass error:
 *   Em = | sum_e A_e C_e(t1) - A_e C_e(t0) | / | sum_e A_e C_e(t0) |
 *
 * Geometric error:
 *   Eg = sum_e A_e | C_e(t1) - C_e(t0) |
 *
 * Quadrature (Option B):
 *   A_e      = sum_g wJ_g
 *   Psi_g(t) = sum_a Psi_a(t) * phi_a(g)
 *   C_e(t)   = [sum_g wJ_g * H(Psi_g(t))] / A_e
 *
 * Notes:
 * - single mesh (from Field::mesh())
 * - PsiS and PsiE must be NODAL fields stored in Field with names "PsiS" and "PsiE"
 * - 'weight' returned by fem().Jacobian(...) is assumed to be (gauss_weight * detJ) for that qp.
 */
inline std::tuple<double, double, double>
computeMassAndGeometricError(
    const std::string Psi0,
    const std::string Psi1,
    const Field& F,
    const std::array<std::unique_ptr<FemProjection>, 6>& elProj)
{
  const Mesh& m = F.mesh();

  const std::size_t dim = m.dim();
  const std::size_t nEl = m.numElements();
  const std::size_t nN  = m.numNodes();

  const auto& elTplgy = m.elTplgy();
  const auto& elType  = m.elType();
  const auto& Xall    = m.X();

  if (elTplgy.size() != nEl) throw std::runtime_error("mass/geom: elTplgy size mismatch");
  if (elType.size()  != nEl) throw std::runtime_error("mass/geom: elType size mismatch");
  if (Xall.size()    != dim) throw std::runtime_error("mass/geom: X dimension mismatch");

  for (std::size_t d = 0; d < dim; ++d)
    if (Xall[d].size() != nN) throw std::runtime_error("mass/geom: mesh.X()[d] size mismatch");

  // Get PsiS / PsiE nodal vectors (Field API from your Field.hpp)
  const unsigned idS = F.id(Psi0);
  const unsigned idE = F.id(Psi1);

  const Field::Vec& PsiS = F.getNodalById(idS);
  const Field::Vec& PsiE = F.getNodalById(idE);

  if (PsiS.size() != nN) throw std::runtime_error("mass/geom: PsiS size mismatch vs mesh.numNodes()");
  if (PsiE.size() != nN) throw std::runtime_error("mass/geom: PsiE size mismatch vs mesh.numNodes()");

  double M0 = 0.0;
  double M1 = 0.0;
  double Eg = 0.0;

  // Workspaces reused to avoid allocations in hot loops
  std::vector<std::vector<double>> X;     // X[dim][nen]
  std::vector<double> phi;               // size nen
  std::vector<double> phi_x;             // size (dim*nen) or whatever your FEM fills

  X.resize(dim);

  for (std::size_t e = 0; e < nEl; ++e) {
    const unsigned et = elType[e];
    if (et >= elProj.size()) throw std::runtime_error("mass/geom: element type out of [0,5]");
    if (!elProj[et])         throw std::runtime_error("mass/geom: elProj[et] is null");

    const auto& conn = elTplgy[e];
    const std::size_t nen = conn.size();
    if (nen == 0) throw std::runtime_error("mass/geom: empty connectivity");

    // Build element coordinates X(dim, nen)
    for (std::size_t d = 0; d < dim; ++d) {
      X[d].assign(nen, 0.0);
      for (std::size_t a = 0; a < nen; ++a) {
        const unsigned node = conn[a];
        if (node >= nN) throw std::runtime_error("mass/geom: node id out of range");
        X[d][a] = Xall[d][node];
      }
    }

    // FEM handle (your Field uses elProj[et]->fem().GetPhi(...), so fem() is a ref)
    auto& fem = elProj[et]->fem();

    const unsigned nGp = fem.GetGaussPointNumber();
    if (nGp == 0u) throw std::runtime_error("mass/geom: fem.GetGaussPointNumber()==0");

    double A   = 0.0;  // element measure (area in 2D, volume in 3D) from quadrature
    double Ap0 = 0.0;  // positive-phase measure at t0
    double Ap1 = 0.0;  // positive-phase measure at t1

    for (unsigned ig = 0; ig < nGp; ++ig) {
      double weight = 0.0;

      phi.clear();
      phi_x.clear();
      fem.Jacobian(X, ig, weight, phi, phi_x);

      // weight is gauss_weight * detJ at this qp
      A += weight;

      if (phi.size() != nen) {
        throw std::runtime_error("mass/geom: phi.size() != nen (connectivity size)");
      }

      double psiS_q = 0.0;
      double psiE_q = 0.0;

      // Psi at quadrature point via shape functions
      for (std::size_t a = 0; a < nen; ++a) {
        const unsigned node = conn[a];
        // node range already checked above, but cheap to keep safe:
        if (node >= nN) throw std::runtime_error("mass/geom: node id out of range (qp loop)");

        const double ph = phi[a];
        psiS_q += PsiS[node] * ph;
        psiE_q += PsiE[node] * ph;
      }

      Ap0 += weight * Hpos(psiS_q);
      Ap1 += weight * Hpos(psiE_q);
    }

    if (!(A > 0.0)) throw std::runtime_error("mass/geom: non-positive element measure A");

    const double C0 = Ap0 / A;
    const double C1 = Ap1 / A;

    // M(t) = sum_e A*C (== sum_e Ap)
    M0 += A * C0;  // equals Ap0
    M1 += A * C1;  // equals Ap1

    // geometric error
    Eg += A * std::abs(C1 - C0);
  }

  const double denom = std::abs(M0);
  const double Em = (denom > 0.0) ? (std::abs(M1 - M0) / denom) : 0.0;
  const double Egs = (denom > 0.0) ? (Eg / denom) : 0.0;

  std::cout << "mass = " << M1 << " " << "Em = " << Em << " " << " Eg = " << Eg << " Egs = " << Egs << std::endl;

  return {Em, Eg, Egs};
}
