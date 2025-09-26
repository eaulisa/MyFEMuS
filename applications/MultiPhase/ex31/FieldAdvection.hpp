#pragma once
#include "QuadTree2DNew.hpp"
#include <array>
#include <vector>
#include <functional>
#include <cassert>

namespace fem {

  inline std::array<double, 2> eval_velocity_coarse(
    const QuadTree2D& qt,
    Basis velBasis,
    const std::vector<std::array<double, 2>>& vOld,
    const std::vector<std::array<double, 2>>& vNew,
    double x, double y, double t, double dt) {
    double xi, eta;
    if (!qt.inverse_map_quad9(x, y, xi, eta)) {
      return {0.0, 0.0}; // outside coarse element
    }

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

    double theta = (dt > 0.0) ? t / dt : 1. - t / dt;

    std::array<double, 2> v{0.0, 0.0};
    for (size_t a = 0; a < N.size(); ++a) {
      v[0] += N[a] * ((1.0 - theta) * vOld[a][0] + theta * vNew[a][0]);
      v[1] += N[a] * ((1.0 - theta) * vOld[a][1] + theta * vNew[a][1]);
    }

    return v;
  }


//----------------------------------------
// 2. Generic RK4 advector
//----------------------------------------
  inline std::array<double, 2> rk4_advect(
    const std::array<double, 2>& x0,
    double dt,
    const std::function<std::array<double, 2>(const std::array<double, 2>&, double)>& vel_eval) {
    auto add = [](const std::array<double, 2>& a, const std::array<double, 2>& b, double s) {
      return std::array<double, 2> { a[0] + s*b[0], a[1] + s*b[1] };
    };

    std::array<double, 2> k1 = vel_eval(x0, 0.0);
    std::array<double, 2> k2 = vel_eval(add(x0, k1, 0.5 * dt), 0.5 * dt);
    std::array<double, 2> k3 = vel_eval(add(x0, k2, 0.5 * dt), 0.5 * dt);
    std::array<double, 2> k4 = vel_eval(add(x0, k3, dt), dt);

    std::array<double, 2> x = x0;
    x[0] += dt / 6.0 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
    x[1] += dt / 6.0 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);

    return x;
  }

//----------------------------------------
// 3. Forward advection of finest-level markers
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
      auto vel_eval = [&](const std::array<double, 2>& x, double t) {
        return eval_velocity_coarse(qt0, velBasis, vOld, vNew, x[0], x[1], t, dt);
      };

      auto p1 = rk4_advect(p0, dt, vel_eval);

      u32 leaf = qt0.locate_physical(p1[0], p1[1]);
      if (leaf == npos32) {
        coordsLeftOld.push_back(p0);   // store old pos
      }
      else {
        coordsStayedNew.push_back(p1); // store advected pos
      }
    }
  }

//----------------------------------------
// 4. Backward advection & field transport
//----------------------------------------
  inline void advect_nodes_backward_and_transport_field(
    const QuadTree2D& qt0,
    u32 fid0,                  // index of source field in qt0
    Basis velBasis,
    const std::vector<std::array<double, 2>>& vOld,
    const std::vector<std::array<double, 2>>& vNew,
    double dt,
    QuadTree2D& qt1,
    u32 fid1) {                // index of destination field in qt1
    auto coords1 = qt1.extract_node_coords_in_level_range(0, qt1.max_depth(), Basis::Q2_Quad9);

    Field& fld1 = qt1.field(fid1);
    fld1.coeffs.resize(coords1.size());

    for (size_t i = 0; i < coords1.size(); ++i) {
      auto p1 = coords1[i];

      auto vel_eval = [&](const std::array<double, 2>& x, double t) {
        return eval_velocity_coarse(qt0, velBasis, vOld, vNew, x[0], x[1], t, dt);
      };

      // Backward advection with -dt
      auto p0 = rk4_advect(p1, -dt, vel_eval);

      // Check if inside coarse element
      u32 leaf = qt0.locate_physical(p0[0], p0[1]);
      if (leaf == npos32) {
        fld1.coeffs[i] = -1.0;
      }
      else {
        double val;
        bool ok = qt0.evaluate_physical(fid0, p0[0], p0[1], val);
        fld1.coeffs[i] = ok ? val : -1.0;
      }
    }
  }

} // namespace fem

