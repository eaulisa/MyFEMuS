#pragma once
#include "QuadTree2D.hpp"
#include <array>
#include <vector>
#include <cassert>
#include <cmath>

namespace fem {

//----------------------------------------
// 0) Utility: make a parent-space velocity evaluator
//    that converts RK local time offset τ to ABSOLUTE time.
//    start_time is absolute time at RK stage τ = 0.
//----------------------------------------
  template<class AnalyticVel>
  inline auto make_parent_vel_eval(
    const QuadTree2D& qt,
    double start_time,          // absolute time where RK local τ=0
    AnalyticVel&& vfun) {       // vfun(x,y,t_abs) -> {vx,vy}
    return [&](const std::array<double, 2>& xiEta, double tau) {
      const double t_abs = start_time + tau;
      // Map (xi,eta) -> (x,y), evaluate v(x,y,t_abs), and map to parent coords
      // using J^{-1}.
      // (Implemented below by eval_velocity_parent_analytic.)
      return eval_velocity_parent_analytic(qt, xiEta[0], xiEta[1], t_abs, vfun);
    };
  }

//----------------------------------------
// 1) Velocity evaluator from analytic function.
//    Evaluates v(x,y,t_abs) in physical space, then returns (ξ̇,η̇).
//----------------------------------------
  template<class AnalyticVel>
  inline std::array<double, 2> eval_velocity_parent_analytic(
    const QuadTree2D& qt,
    double xi, double eta,
    double t_abs,
    AnalyticVel&& vfun) {
    // -- Map to physical with Q2_Quad9 geometry
    double N9[9];
    Quad9Shape::N(xi, eta, N9);

    double x = 0.0, y = 0.0;
    const auto& X = qt.get_X();
    const auto& Y = qt.get_Y();
    for (int a = 0; a < 9; ++a) {
      x += N9[a] * X[a];
      y += N9[a] * Y[a];
    }

    // -- Physical velocity at absolute time t_abs
    std::array<double, 2> vPhys = vfun(x, y, t_abs); // {vx, vy}

    // -- Jacobian and its inverse
    double dNdxi[9], dNdeta[9];
    Quad9Shape::dN(xi, eta, dNdxi, dNdeta);

    double dXdxi = 0, dXdeta = 0, dYdxi = 0, dYdeta = 0;
    for (int a = 0; a < 9; ++a) {
      dXdxi  += dNdxi[a]  * X[a];
      dXdeta += dNdeta[a] * X[a];
      dYdxi  += dNdxi[a]  * Y[a];
      dYdeta += dNdeta[a] * Y[a];
    }
    const double detJ = dXdxi * dYdeta - dXdeta * dYdxi;
    assert(std::abs(detJ) > 1e-14);

    // -- Map to parent coords: (ξ̇,η̇) = J^{-1} v
    std::array<double, 2> vParent;
    vParent[0] = ( dYdeta * vPhys[0] - dXdeta * vPhys[1]) / detJ;
    vParent[1] = (-dYdxi  * vPhys[0] + dXdxi  * vPhys[1]) / detJ;
    return vParent;
  }

//----------------------------------------
// 3) Forward advection of markers (analytic velocity).
//    Integrates from absolute (time - dt) to (time).
//----------------------------------------
  template<class AnalyticVel>
  inline void advect_markers_forward_analytic(
    const QuadTree2D& qt0,
    double time,                 // absolute END time t^{n+1}
    double dt,                   // dt > 0
    AnalyticVel&& vfun,
    std::vector<std::array<double, 2>>& coordsLeftOld,
    std::vector<std::array<double, 2>>& coordsStayedNew) {
    assert(dt > 0.0 && "Forward advection expects dt > 0");
    coordsLeftOld.clear();
    coordsStayedNew.clear();

    // We advect the nodes at the finest level, like your original code.
    auto coords0 = qt0.extract_node_parent_coords_in_level_range(
                     qt0.max_depth(), qt0.max_depth(), Basis::Q2_Quad9);

    // RK start absolute time
    const double t0_abs = time - dt;

    for (const auto& xiEta0 : coords0) {
      // Build evaluator with absolute start time t0_abs.
      auto vel_eval = make_parent_vel_eval(qt0, t0_abs, vfun);

      // RK4 with +dt: integrates from t0_abs to time.
      const auto xiEta1 = rk4_advect_parent(xiEta0, dt, vel_eval);

      // Map back to physical, or record as "left"
      u32 leaf = qt0.locate_leaf_on_parent(xiEta1[0], xiEta1[1]);
      if (leaf == npos32) {
        coordsLeftOld.push_back(qt0.parent_to_physical(xiEta0[0], xiEta0[1]));
      }
      else {
        double N9[9];
        Quad9Shape::N(xiEta1[0], xiEta1[1], N9);
        double x1 = 0.0, y1 = 0.0;
        const auto& X = qt0.get_X();
        const auto& Y = qt0.get_Y();
        for (int a = 0; a < 9; ++a) {
          x1 += N9[a] * X[a];
          y1 += N9[a] * Y[a];
        }
        coordsStayedNew.push_back({x1, y1});
      }
    }
  }

//----------------------------------------
// 4) Backward advection & field transport (analytic).
//    Traces back from absolute time `time` to `time - dt`.
//    We integrate with -dt, and the RK local τ is added to START=time.
//----------------------------------------
  template<class AnalyticVel>
  inline void advect_nodes_backward_and_transport_field_analytic(
    const QuadTree2D& qt0,   // geometry + field at old time (time - dt)
    u32 fid0,                // field id in qt0 to sample from
    double time,             // absolute END time t^{n+1}
    double dt,               // dt > 0
    AnalyticVel&& vfun,      // v(x,y,t_abs)
    QuadTree2D& qt1,         // geometry at END time (time)
    u32 fid1) {              // output field id in qt1
    assert(dt > 0.0 && "Backward advection here expects dt > 0");
    auto coords1 = qt1.extract_node_parent_coords_in_level_range(
                     0, qt1.max_depth(), Basis::Q2_Quad9);

    Field& fld1 = qt1.field(fid1);
    fld1.coeffs.resize(coords1.size());

    // RK starts at absolute time `time` and steps with τ in {0, -dt/2, -dt/2, -dt}
    const double start_abs = time;

    for (size_t i = 0; i < coords1.size(); ++i) {
      const auto xiEta1 = coords1[i];

      // Build evaluator with absolute start time = time (end of step).
      auto vel_eval = make_parent_vel_eval(qt0, start_abs, vfun);

      // RK4 with -dt: integrates from time to time - dt
      const auto xiEta0 = rk4_advect_parent(xiEta1, -dt, vel_eval);

      // Sample field directly in parent coords of qt0
      double val;
      bool ok = qt0.evaluate_field_on_parent(fid0, xiEta0[0], xiEta0[1], val);
      fld1.coeffs[i] = ok ? val : -1.0;
    }
  }

} // namespace fem
