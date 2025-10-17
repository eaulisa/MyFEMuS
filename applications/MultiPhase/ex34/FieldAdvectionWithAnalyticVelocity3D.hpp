#pragma once
#include "HexTree3D.hpp"
#include <array>
#include <vector>
#include <cassert>
#include <cmath>
#include <limits>

namespace fem {

//----------------------------------------
// 0) Utility: parent-space velocity evaluator with absolute time
//    vfun(x,y,z,t_abs) -> {vx,vy,vz}
//----------------------------------------
  template<class AnalyticVel3D>
  inline auto make_parent_vel_eval(
    const HexTree3D& hex,
    double start_time,                 // absolute time where RK local τ=0
    AnalyticVel3D&& vfun) {
    return [&](const Point3 & s, double tau) {
      const double t_abs = start_time + tau;
      return eval_velocity_parent_analytic(hex, s, t_abs, vfun);
    };
  }

//----------------------------------------
// 1) Velocity evaluator from analytic function (3D).
//    Map (ξ,η,ζ)->(x,y,z), eval v(x,y,z,t_abs), then return (ξ̇,η̇,ζ̇).
//----------------------------------------
  template<class AnalyticVel3D>
  inline Point3 eval_velocity_parent_analytic(
    const HexTree3D& hex,
    const Point3& s,
    double t_abs,
    AnalyticVel3D&& vfun) {
    // -- Physical coordinates via H27
    double N27[27];
    Shapes3::H27(s, N27);

    const auto& X = hex.get_X();
    const auto& Y = hex.get_Y();
    const auto& Z = hex.get_Z();

    double x = 0.0, y = 0.0, z = 0.0;
    for (int a = 0; a < 27; ++a) {
      x += N27[a] * X[a];
      y += N27[a] * Y[a];
      z += N27[a] * Z[a];
    }

    // -- Physical velocity
    const Point3 v = vfun(x, y, z, t_abs); // {vx,vy,vz}

    // -- Jacobian J = d(x,y,z)/d(ξ,η,ζ) and inverse
    double dNdxi[27], dNdeta[27], dNdz[27];
    Shapes3::H27_dN(s, dNdxi, dNdeta, dNdz);

    double A11 = 0, A12 = 0, A13 = 0,
           A21 = 0, A22 = 0, A23 = 0,
           A31 = 0, A32 = 0, A33 = 0;

    for (int a = 0; a < 27; ++a) {
      A11 += dNdxi[a] * X[a];  A12 += dNdeta[a] * X[a];  A13 += dNdz[a] * X[a];
      A21 += dNdxi[a] * Y[a];  A22 += dNdeta[a] * Y[a];  A23 += dNdz[a] * Y[a];
      A31 += dNdxi[a] * Z[a];  A32 += dNdeta[a] * Z[a];  A33 += dNdz[a] * Z[a];
    }

    const double det =
      A11 * (A22 * A33 - A23 * A32)
      - A12 * (A21 * A33 - A23 * A31)
      + A13 * (A21 * A32 - A22 * A31);
    assert(std::abs(det) > 1e-20 && "singular hex mapping Jacobian");

    const double inv = 1.0 / det;

    // ṡ = J^{-1} v (adjugate * v / det)
    Point3 sdot;
    sdot[0] = ((A22 * A33 - A23 * A32) * v[0]
               - (A12 * A33 - A13 * A32) * v[1]
               + (A12 * A23 - A13 * A22) * v[2]) * inv;

    sdot[1] = (-(A21 * A33 - A23 * A31) * v[0]
               + (A11 * A33 - A13 * A31) * v[1]
               - (A11 * A23 - A13 * A21) * v[2]) * inv;

    sdot[2] = ((A21 * A32 - A22 * A31) * v[0]
               - (A11 * A32 - A12 * A31) * v[1]
               + (A11 * A22 - A12 * A21) * v[2]) * inv;

    return sdot;
  }

//----------------------------------------
// 2) Forward advection of markers (analytic velocity, 3D).
//    Integrates from (time - dt) to (time).
//    Uses existing rk4_advect_parent(Point3, dt, vel_eval).
//----------------------------------------
  template<class AnalyticVel3D>
  inline void advect_markers_forward_analytic(
    const HexTree3D& hex0,
    double time,                 // absolute END time t^{n+1}
    double dt,                   // dt > 0
    AnalyticVel3D&& vfun,
    std::vector<Point3>& coordsLeftOld,
    std::vector<Point3>& coordsStayedNew) {
    assert(dt > 0.0 && "Forward advection expects dt > 0");
    coordsLeftOld.clear();
    coordsStayedNew.clear();

    // sample nodes at finest level (H27)
    const auto s0_all =
      hex0.extract_node_parent_coords_in_level_range(hex0.max_depth(), hex0.max_depth(), Basis::H27);

    const double t0_abs = time - dt;

    for (const auto& s0 : s0_all) {
      auto vel_eval = make_parent_vel_eval(hex0, t0_abs, vfun);

      const Point3 s1 = rk4_advect_parent(s0, dt, vel_eval); // <- provided elsewhere

      const u32 leaf = hex0.locate_leaf_on_parent(s1);
      if (leaf == npos32) {
        coordsLeftOld.push_back(hex0.parent_to_physical(s0));
      }
      else {
        // physical at s1 via H27
        double N27[27];
        Shapes3::H27(s1, N27);
        const auto& X = hex0.get_X();
        const auto& Y = hex0.get_Y();
        const auto& Z = hex0.get_Z();
        double x = 0, y = 0, z = 0;
        for (int a = 0; a < 27; ++a) {
          x += N27[a] * X[a];
          y += N27[a] * Y[a];
          z += N27[a] * Z[a];
        }
        coordsStayedNew.push_back({x, y, z});
      }
    }
  }

//----------------------------------------
// 3) Backward advection & field transport (analytic, 3D).
//    Traces back from absolute time `time` to `time - dt`.
//    Uses existing rk4_advect_parent(Point3, -dt, vel_eval).
//----------------------------------------
  template<class AnalyticVel3D>
  inline void advect_nodes_backward_and_transport_field_analytic(
    const HexTree3D& hex0,   // geometry + field at old time (time - dt)
    u32 fid0,                // field id in hex0 to sample from
    double time,             // absolute END time t^{n+1}
    double dt,               // dt > 0
    AnalyticVel3D&& vfun,    // v(x,y,z,t_abs)
    HexTree3D& hex1,         // geometry at END time (time)
    u32 fid1) {              // output field id in hex1
    assert(dt > 0.0 && "Backward advection expects dt > 0");

    Field& dst = hex1.field(fid1);
    const Basis b = dst.basis;
    const auto& nodes = hex1.basis_nodes(b);
    dst.nodal.assign(nodes.size(), 0.0);

    const double start_abs = time;
    auto vel_eval = make_parent_vel_eval(hex0, start_abs, vfun);

    for (size_t gid = 0; gid < nodes.size(); ++gid) {
      const Point3 s1 = nodes[gid].parent;                 // parent coords in hex1
      const Point3 s0 = rk4_advect_parent(s1, -dt, vel_eval); // <- provided elsewhere

      double val = std::numeric_limits<double>::quiet_NaN();
      (void) hex0.evaluate_field_on_parent(fid0, s0, val);
      dst.nodal[gid] = val;
    }
  }

} // namespace fem

