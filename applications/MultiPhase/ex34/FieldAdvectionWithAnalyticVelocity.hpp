#pragma once
#include "OctTree.hpp"
#include <array>
#include <vector>
#include <cassert>
#include <cmath>
#include <limits>

namespace fem {

// ============================================================
// Analytic parent-velocity: 2D overload
// eval_parent_sdot(tree, s, t_abs, vfun)
// vfun: v(x,y,t) -> Point<2>
// ============================================================
  template<class AnalyticVel>
  inline Point<2> eval_parent_sdot(
    const OctTree<2>& tree,
    const Point<2>& s,
    double t_abs,
    const AnalyticVel& vfun) {
    // Map to physical with Q9
    double N9[9]; Shapes2D::Q9(s, N9);
    const auto& X = tree.get_X();
    const auto& Y = tree.get_Y();

    double x = 0.0, y = 0.0;
    for (int a = 0; a < 9; ++a) {
      x += N9[a] * X[a];
      y += N9[a] * Y[a];
    }

    // Physical velocity
    const Point<2> v = vfun(x, y, t_abs);

    // J = d(x,y)/d(ξ,η) via Q9_dN
    double dNdxi[9], dNdeta[9];
    Shapes2D::Q9_dN(s, dNdxi, dNdeta);

    double A11 = 0, A12 = 0, A21 = 0, A22 = 0;
    for (int a = 0; a < 9; ++a) {
      A11 += dNdxi[a] * X[a];  A12 += dNdeta[a] * X[a];
      A21 += dNdxi[a] * Y[a];  A22 += dNdeta[a] * Y[a];
    }

    const double det = A11 * A22 - A12 * A21;
    assert(std::fabs(det) > 1e-20 && "Singular 2D mapping Jacobian");
    const double inv = 1.0 / det;

    Point<2> sdot{};
    sdot[0] = (A22 * v[0] - A12 * v[1]) * inv;
    sdot[1] = (-A21 * v[0] + A11 * v[1]) * inv;
    return sdot;
  }

// ============================================================
// Analytic physical-velocity: 2D overload
// eval_physical(s, t_abs, vfun)
// vfun: v(x,y,t) -> Point<2>
// ============================================================
  template<class AnalyticVel>
  inline Point<2> eval_physical(
    const Point<2>& s,
    double t_abs,
    const AnalyticVel& vfun) {
    // Map to physical with H27

    // Physical velocity
    const Point<2> v = vfun(s[0], s[1], t_abs);

    return v;
  }

// ============================================================
// Analytic parent-velocity: 3D overload
// eval_parent_sdot(tree, s, t_abs, vfun)
// vfun: v(x,y,z,t) -> Point<3>
// ============================================================
  template<class AnalyticVel>
  inline Point<3> eval_parent_sdot(
    const OctTree<3>& tree,
    const Point<3>& s,
    double t_abs,
    const AnalyticVel& vfun) {
    // Map to physical with H27
    double N27[27]; Shapes3D::H27(s, N27);
    const auto& X = tree.get_X();
    const auto& Y = tree.get_Y();
    const auto& Z = tree.get_Z();

    double x = 0, y = 0, z = 0;
    for (int a = 0; a < 27; ++a) {
      x += N27[a] * X[a];
      y += N27[a] * Y[a];
      z += N27[a] * Z[a];
    }

    // Physical velocity
    const Point<3> v = vfun(x, y, z, t_abs);

    // J = d(x,y,z)/d(ξ,η,ζ) via H27_dN
    double dNdxi[27], dNdeta[27], dNdz[27];
    Shapes3D::H27_dN(s, dNdxi, dNdeta, dNdz);

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
    assert(std::fabs(det) > 1e-20 && "Singular 3D mapping Jacobian");
    const double inv = 1.0 / det;

    Point<3> sdot{};
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

// ============================================================
// Analytic physical-velocity: 3D overload
// eval_physical(s, t_abs, vfun)
// vfun: v(x,y,z,t) -> Point<3>
// ============================================================
  template<class AnalyticVel>
  inline Point<3> eval_physical(
    const Point<3>& s,
    double t_abs,
    const AnalyticVel& vfun) {
    // Map to physical with H27

    // Physical velocity
    const Point<3> v = vfun(s[0], s[1], s[2], t_abs);

    return v;
  }

// ============================================================
// Make parent-space velocity evaluator with absolute time
// Returns: lambda (s, tau) -> sdot
// Works for DIM=2 or 3 by picking the right overload above
// ============================================================
  template<class AnalyticVel, std::size_t DIM>
  inline auto make_parent_vel_eval(
    const OctTree<DIM>& tree,
    double start_time,
    AnalyticVel&& vfun) {
    // capture by value to avoid ref-binding issues with function pointers
    auto   vfun_captured          = std::forward<AnalyticVel>(vfun);
    const double start_time_local = start_time;

    return [&tree, start_time_local, vfun_captured](const Point<DIM>& s, double tau) {
      const double t_abs = start_time_local + tau;
      return eval_parent_sdot(tree, s, t_abs, vfun_captured); // overload picks 2D/3D
    };
  }

// ============================================================
// Make physical-space velocity evaluator with absolute time
// Returns: lambda (s, tau) -> sdot
// Works for DIM=2 or 3 by picking the right overload above
// ============================================================
  template<class AnalyticVel, std::size_t DIM>
  inline auto make_physical_vel_eval(
    double start_time,
    AnalyticVel&& vfun) {
    // capture by value to avoid ref-binding issues with function pointers
    auto   vfun_captured          = std::forward<AnalyticVel>(vfun);
    const double start_time_local = start_time;

    return [start_time_local, vfun_captured](const Point<DIM>& s, double tau) {
      const double t_abs = start_time_local + tau;
      return eval_physical(s, t_abs, vfun_captured); // overload picks 2D/3D
    };
  }

// ======================================================================
// 2) Forward advection of markers (analytic velocity, DIM=2/3)
//    Integrates from (time - dt) to (time).
// ======================================================================
  template<class AnalyticVel, std::size_t DIM>
  inline void advect_interface_markers_forward_analytic(
    const OctTree<DIM>& tree0,
    const u32 &fid0,
    double time,                 // absolute END time t^{n+1}
    double dt,                   // dt > 0
    AnalyticVel&& vfun,
    std::vector<Point<DIM>>& coordsLeftOld,
    std::vector<Point<DIM>>& coordsStayedNew) {
    assert(dt > 0.0 && "Forward advection expects dt > 0");
    coordsLeftOld.clear();
    coordsStayedNew.clear();

    //sample nodes at finest level using highest-order geometry basis (id==0 ⇒ Q4/H8)
    std::vector<Point<DIM>> s0_all;
    //tree0.extract_node_parent_coords_in_level_range(tree0.max_depth(), tree0.max_depth(), static_cast<BasisT<DIM>>(0), s0_all);

    tree0.extract_interface_leaf_parent_coords_in_level_range(tree0.max_depth(), tree0.max_depth(),
                                                              static_cast<BasisT<DIM>>(2), fid0,
                                                              s0_all);


    const double t0_abs = time - dt;

    for (const auto& s0 : s0_all) {
      auto vel_eval = make_parent_vel_eval(tree0, t0_abs, vfun);
      const Point<DIM> s1 = rk4_advect_parent<DIM>(s0, dt, vel_eval);

      const u32 leaf = tree0.locate_leaf_on_parent(s1);
      if (leaf == npos32) {
        coordsLeftOld.push_back(tree0.parent_to_physical(s0));
      }
      else {
        coordsStayedNew.push_back(tree0.parent_to_physical(s1));
      }
    }
  }

// ======================================================================
// 2) Forward advection of markers (analytic velocity, DIM=2/3)
//    Integrates from (time - dt) to (time).
// ======================================================================
  template<class AnalyticVel, std::size_t DIM>
  inline void advect_physical_markers_forward_analytic(
    const OctTree<DIM>& tree0,
    double time,                 // absolute END time t^{n+1}
    double dt,                   // dt > 0
    AnalyticVel&& vfun,
    const std::vector<Point<DIM>>& markers,
    std::vector<Point<DIM>>& coordsLeftOld,
    std::vector<Point<DIM>>& coordsStayedNew) {
    assert(dt > 0.0 && "Forward advection expects dt > 0");
    coordsLeftOld.clear();
    coordsStayedNew.clear();

    //sample nodes at finest level using highest-order geometry basis (id==0 ⇒ Q4/H8)
    const double t0_abs = time - dt;

    for (const auto& s0 : markers) {
      auto vel_eval = make_physical_vel_eval<decltype(vfun), DIM>(t0_abs, vfun);

      const Point<DIM> s1 = rk4_advect_parent<DIM>(s0, dt, vel_eval);

      const u32 leaf = tree0.locate_leaf_on_parent(s1);
      if (leaf == npos32) {
        coordsLeftOld.push_back(s0);
      }
      else {
        coordsStayedNew.push_back(s1);
      }
    }
  }

// ======================================================================
// 3) Backward advection & field transport (analytic, DIM=2/3).
//    Traces back from absolute time `time` to `time - dt`.
// ======================================================================
  template<class AnalyticVel, std::size_t DIM>
  inline void advect_nodes_backward_and_transport_field_analytic(
    const OctTree<DIM>& tree0,  // geometry + source field at old time (time - dt)
    u32 fid0,                   // field id in tree0 to sample from
    double time,                // absolute END time t^{n+1}
    double dt,                  // dt > 0
    AnalyticVel&& vfun,         // v(x,y[,z], t_abs) -> Point<DIM>
    OctTree<DIM>& tree1,        // geometry at END time (time)
    u32 fid1) {                 // output field id in tree1
    assert(dt > 0.0 && "Backward advection expects dt > 0");

    tree1.rebuild_connectivity_active_bases();

    Field& dst = tree1.field(fid1);
    const auto b = to_basis<DIM>(dst.basis_id);
    const auto& nodes = tree1.basis_nodes(b);
    dst.nodal.assign(nodes.size(), 0.0);

    // local time origin at t^{n+1}; integrate with negative dt
    auto vel_eval = make_parent_vel_eval(tree0, time, vfun);

    for (size_t gid = 0; gid < nodes.size(); ++gid) {
      const Point<DIM> s1 = nodes[gid].parent;                 // parent coords in tree1
      const Point<DIM> s0 = rk4_advect_parent<DIM>(s1, -dt, vel_eval);

      double val = std::numeric_limits<double>::quiet_NaN();
      (void)tree0.evaluate_field_on_parent(fid0, s0, val);
      dst.nodal[gid] = val;
    }
  }

} // namespace fem
