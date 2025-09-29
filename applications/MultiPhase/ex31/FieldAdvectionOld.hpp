#pragma once
#include "QuadTree2DNew.hpp"
#include <array>
#include <vector>
#include <functional>
#include <cassert>

namespace fem {

//----------------------------------------
// 1. Evaluate velocity in parent coords
//    Returns (dξ/dt, dη/dt)
//----------------------------------------
  inline std::array<double, 2> eval_velocity_parent(
    const QuadTree2D& qt,
    Basis velBasis,
    const std::vector<std::array<double, 2>>& vOld,
    const std::vector<std::array<double, 2>>& vNew,
    double xi, double eta,
    double t, double dt) {
    // Shape functions
    std::vector<double> N;
    if (velBasis == Basis::Q1_Quad4) {
      N.resize(4);
      Shapes::Q1(xi, eta, N.data());
    }
    else if (velBasis == Basis::Serendipity8) {
      N.resize(8);
      Shapes::Serendipity8(xi, eta, N.data());
    }
    else { // Q2_Quad9
      N.resize(9);
      Quad9Shape::N(xi, eta, N.data());
    }

    // Interpolate velocity in physical space
    double theta = (dt > 0.0) ? (t / dt) : (1.0 - t / dt);
    std::array<double, 2> v{0.0, 0.0};
    for (size_t a = 0; a < N.size(); ++a) {
      v[0] += N[a] * ((1.0 - theta) * vOld[a][0] + theta * vNew[a][0]);
      v[1] += N[a] * ((1.0 - theta) * vOld[a][1] + theta * vNew[a][1]);
    }

    // Build Jacobian at (xi,eta)
    double dNdxi[9], dNdeta[9];
    Quad9Shape::dN(xi, eta, dNdxi, dNdeta);
    double dXdxi = 0, dXdeta = 0, dYdxi = 0, dYdeta = 0;
    const auto& X = qt.get_X();
    const auto& Y = qt.get_Y();

    for (int a = 0; a < 9; ++a) {
      dXdxi  += dNdxi[a] * X[a];
      dXdeta += dNdeta[a] * X[a];
      dYdxi  += dNdxi[a] * Y[a];
      dYdeta += dNdeta[a] * Y[a];
    }

    double detJ = dXdxi * dYdeta - dXdeta * dYdxi;
    assert(std::abs(detJ) > 1e-14);

    // Map velocity to parent coords: (ξ̇,η̇) = J⁻¹ v
    std::array<double, 2> vParent;
    vParent[0] = (dYdeta * v[0] - dXdeta * v[1]) / detJ;
    vParent[1] = (-dYdxi * v[0] + dXdxi * v[1]) / detJ;

    return vParent;
  }

//----------------------------------------
// 2. RK4 in parent space
//----------------------------------------
  inline std::array<double, 2> rk4_advect_parent(
    const std::array<double, 2>& xi0eta0,
    double dt,
    const std::function<std::array<double, 2>(const std::array<double, 2>&, double)>& vel_eval) {
    auto add = [](const std::array<double, 2>& a, const std::array<double, 2>& b, double s) {
      return std::array<double, 2> { a[0] + s*b[0], a[1] + s*b[1] };
    };

    std::array<double, 2> k1 = vel_eval(xi0eta0, 0.0);
    std::array<double, 2> k2 = vel_eval(add(xi0eta0, k1, 0.5 * dt), 0.5 * dt);
    std::array<double, 2> k3 = vel_eval(add(xi0eta0, k2, 0.5 * dt), 0.5 * dt);
    std::array<double, 2> k4 = vel_eval(add(xi0eta0, k3, dt), dt);

    std::array<double, 2> xiEta = xi0eta0;
    xiEta[0] += dt / 6.0 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
    xiEta[1] += dt / 6.0 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);

    return xiEta;
  }

//----------------------------------------
// 3. Forward advection of markers (parent RK)
//----------------------------------------
  inline void advect_markers_forward(
    const QuadTree2D& qt0,
    Basis velBasis,
    const std::vector<std::array<double, 2>>& vOld,
    const std::vector<std::array<double, 2>>& vNew,
    double dt,
    std::vector<std::array<double, 2>>& coordsLeftOld,
    std::vector<std::array<double, 2>>& coordsStayedNew) {
    coordsLeftOld.clear();
    coordsStayedNew.clear();

    auto coords0 = qt0.extract_node_coords_in_level_range(qt0.max_depth(), qt0.max_depth(), Basis::Q2_Quad9);

    for (auto& p0 : coords0) {
      // 1. One inverse map at start
      double xi0, eta0;
      if (!qt0.inverse_map_quad9(p0[0], p0[1], xi0, eta0)) {
        coordsLeftOld.push_back(p0);
        continue;
      }

      auto vel_eval = [&](const std::array<double, 2>& xiEta, double t) {
        return eval_velocity_parent(qt0, velBasis, vOld, vNew, xiEta[0], xiEta[1], t, dt);
      };

      // 2. RK4 in parent space
      auto xiEta1 = rk4_advect_parent({xi0, eta0}, dt, vel_eval);

      // 3. Map back to physical
      double N9[9]; Quad9Shape::N(xiEta1[0], xiEta1[1], N9);
      double x1 = 0, y1 = 0;
      const auto& X = qt0.get_X();
      const auto& Y = qt0.get_Y();
      for (int a = 0; a < 9; ++a) {
        x1 += N9[a] * X[a];
        y1 += N9[a] * Y[a];
      }
      u32 leaf = qt0.locate_physical(x1, y1);

      if (leaf == npos32) coordsLeftOld.push_back(p0);
      else coordsStayedNew.push_back({x1, y1});
    }
  }

//----------------------------------------
// 4. Backward advection & field transport (parent RK)
//----------------------------------------
  inline void advect_nodes_backward_and_transport_field(
    const QuadTree2D& qt0,
    u32 fid0,
    Basis velBasis,
    const std::vector<std::array<double, 2>>& vOld,
    const std::vector<std::array<double, 2>>& vNew,
    double dt,
    QuadTree2D& qt1,
    u32 fid1) {
    auto coords1 = qt1.extract_node_coords_in_level_range(0, qt1.max_depth(), Basis::Q2_Quad9);
    Field& fld1 = qt1.field(fid1);
    fld1.coeffs.resize(coords1.size());

    for (size_t i = 0; i < coords1.size(); ++i) {
      auto p1 = coords1[i];

      // 1. Map target node to parent coords of qt0
      double xi1, eta1;
      if (!qt0.inverse_map_quad9(p1[0], p1[1], xi1, eta1)) {
        fld1.coeffs[i] = -1.0; continue;
      }

      auto vel_eval = [&](const std::array<double, 2>& xiEta, double t) {
        return eval_velocity_parent(qt0, velBasis, vOld, vNew, xiEta[0], xiEta[1], t, dt);
      };

      // 2. Backward RK in parent space
      auto xiEta0 = rk4_advect_parent({xi1, eta1}, -dt, vel_eval);

      // 3. Map back to physical
      double N9[9]; Quad9Shape::N(xiEta0[0], xiEta0[1], N9);
      double x0 = 0, y0 = 0;

      const auto& X = qt0.get_X();
      const auto& Y = qt0.get_Y();

      for (int a = 0; a < 9; ++a) {
        x0 += N9[a] * X[a];
        y0 += N9[a] * Y[a];
      }

      // 4. Evaluate source field at (x0,y0)
      double val;
      bool ok = qt0.evaluate_physical(fid0, x0, y0, val);
      fld1.coeffs[i] = ok ? val : -1.0;
    }
  }

} // namespace fem
