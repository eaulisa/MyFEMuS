#pragma once
#include "OctTree.hpp"
#include <array>
#include <vector>
#include <cassert>
#include <limits>

namespace fem {

//----------------------------------------
// 1. Evaluate velocity in parent coords (ξ,η,ζ)
//    vOld/vNew are nodal 3D velocities at the chosen velocity basis
//----------------------------------------
  template<std::size_t DIM>
  inline Point<DIM> eval_velocity_parent(
    const OctTree<DIM, HexTree>& hex,
    BasisT<DIM> velBasis,
    const std::vector<Point<DIM>>& vOld,
    const std::vector<Point<DIM>>& vNew,
    const Point<DIM>& s,          // parent coords (ξ,η,ζ)
    double t, double dt) {    // t ∈ {0, dt/2, dt/2, dt} during RK
    // ---- shapes at s ----
    double N8 [8];
    double N20[20];
    double N27[27];

    const double* N = nullptr;
    int nBasis = 0;
    switch (velBasis) {
    case BasisT<DIM>::H8:  Shapes3::H8(s, N8); N = N8;  nBasis = 8;  break;
    case BasisT<DIM>::H20: Shapes3::H20(s, N20); N = N20; nBasis = 20; break;
    case BasisT<DIM>::H27: Shapes3::H27(s, N27); N = N27; nBasis = 27; break;
    }

    // ---- time interpolation weight ----
    // same convention as your 2D file
    const double theta = (dt > 0.0) ? (t / dt) : (1.0 - t / dt);

    // ---- interpolate physical velocity v(x(s),t) ----
    double vx = 0.0, vy = 0.0, vz = 0.0;
    for (int a = 0; a < nBasis; ++a) {
      const double vxa = (1.0 - theta) * vOld[a][0] + theta * vNew[a][0];
      const double vya = (1.0 - theta) * vOld[a][1] + theta * vNew[a][1];
      const double vza = (1.0 - theta) * vOld[a][2] + theta * vNew[a][2];
      vx += N[a] * vxa;
      vy += N[a] * vya;
      vz += N[a] * vza;
    }

    // ---- geometry Jacobian J(s) from coarse HEX27 map ----
    // (your HexTree geometry is set with set_physical_hex27, so use H27_dN)
    double dNdxi[27], dNdeta[27], dNdz[27];
    Shapes3::H27_dN(s, dNdxi, dNdeta, dNdz);

    const auto& X = hex.get_X();
    const auto& Y = hex.get_Y();
    const auto& Z = hex.get_Z();

    double dXdxi = 0, dXdeta = 0, dXdz = 0,
           dYdxi = 0, dYdeta = 0, dYdz = 0,
           dZdxi = 0, dZdeta = 0, dZdz = 0;

    for (int a = 0; a < 27; ++a) {
      dXdxi  += dNdxi [a] * X[a];  dXdeta += dNdeta[a] * X[a];  dXdz += dNdz[a] * X[a];
      dYdxi  += dNdxi [a] * Y[a];  dYdeta += dNdeta[a] * Y[a];  dYdz += dNdz[a] * Y[a];
      dZdxi  += dNdxi [a] * Z[a];  dZdeta += dNdeta[a] * Z[a];  dZdz += dNdz[a] * Z[a];
    }

    // ---- invert J and return ṡ = J^{-1} v ----
    const double A11 = dXdxi, A12 = dXdeta, A13 = dXdz;
    const double A21 = dYdxi, A22 = dYdeta, A23 = dYdz;
    const double A31 = dZdxi, A32 = dZdeta, A33 = dZdz;

    const double det =
      A11 * (A22 * A33 - A23 * A32)
      - A12 * (A21 * A33 - A23 * A31)
      + A13 * (A21 * A32 - A22 * A31);
    assert(std::abs(det) > 1e-20 && "singular hex mapping Jacobian");

    const double inv = 1.0 / det;

    Point<DIM> sdot;
    sdot[0] = ((A22 * A33 - A23 * A32) * vx - (A12 * A33 - A13 * A32) * vy + (A12 * A23 - A13 * A22) * vz) * inv;
    sdot[1] = (-(A21 * A33 - A23 * A31) * vx + (A11 * A33 - A13 * A31) * vy - (A11 * A23 - A13 * A21) * vz) * inv;
    sdot[2] = ((A21 * A32 - A22 * A31) * vx - (A11 * A32 - A12 * A31) * vy + (A11 * A22 - A12 * A21) * vz) * inv;

    return sdot;
  }

//----------------------------------------
// 2. RK4 in parent space (3D)
//----------------------------------------
  template<class VelEval3D, std::size_t DIM>
  inline Point<DIM> rk4_advect_parent(
    const Point<DIM>& s0,
    double dt,
    VelEval3D&& vel_eval) {
    auto add = [](const Point<DIM> & a, const Point<DIM> & b, double s) {
      return Point<DIM> { a[0] + s*b[0], a[1] + s*b[1], a[2] + s*b[2] };
    };

    const Point<DIM> k1 = vel_eval(s0,              0.0);
    const Point<DIM> k2 = vel_eval(add(s0, k1, 0.5 * dt), 0.5 * dt);
    const Point<DIM> k3 = vel_eval(add(s0, k2, 0.5 * dt), 0.5 * dt);
    const Point<DIM> k4 = vel_eval(add(s0, k3,     dt),     dt);

    Point<DIM> s = s0;
    s[0] += dt / 6.0 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
    s[1] += dt / 6.0 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
    s[2] += dt / 6.0 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]);
    return s;
  }

//----------------------------------------
// 3. Forward advection of markers (parent RK, 3D)
//    Emits physical coords that left the domain and those that stayed
//----------------------------------------
  template<std::size_t DIM>
  inline void advect_markers_forward(
    const HexTree& hex0,
    BasisT<DIM> velBasis,
    const std::vector<Point<DIM>>& vOld,    // size 8/20/27 depending on velBasis
    const std::vector<Point<DIM>>& vNew,    // same size as vOld
    double dt,
    std::vector<Point<DIM>>& coordsLeftOld, // physical positions at t^n that left
    std::vector<Point<DIM>>& coordsStayedNew) { // physical positions at t^{n+1} that stayed
    coordsLeftOld.clear();
    coordsStayedNew.clear();

    // seed: parent nodes on the finest current level (you can change range if needed)
    const auto s0_all = hex0.extract_node_parent_coords_in_level_range(hex0.max_depth(), hex0.max_depth(), BasisT<DIM>::H27);

    for (const auto& s0 : s0_all) {
      // velocity evaluator in parent space (uses hex0 coarse geometry)
      auto vel_eval = [&](const Point<DIM> & s, double t) {
        return eval_velocity_parent(hex0, velBasis, vOld, vNew, s, t, dt);
      };

      // RK4 step in parent
      const Point<DIM> s1 = rk4_advect_parent(s0, dt, vel_eval);

      // classify + output physical coordinates
      const u32 leaf = hex0.locate_leaf_on_parent(s1);
      if (leaf == npos32) {
        coordsLeftOld.push_back(hex0.parent_to_physical(s0));
      }
      else {
        // physical at s1 via coarse H27 map
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
// 4. Backward advection & field transport (parent RK, 3D)
//    Samples source field in parent space on hex0
//----------------------------------------
  template<std::size_t DIM>
  inline void advect_nodes_backward_and_transport_field(
    const HexTree& hex0,    // geometry + source field at old time
    u32 fid0,                 // source field id in hex0
    BasisT<DIM> velBasis,
    const std::vector<Point<DIM>>& vOld,
    const std::vector<Point<DIM>>& vNew,
    double dt,                // > 0
    HexTree& hex1,          // geometry at end time
    u32 fid1) {               // destination field id in hex1
    // make sure destination connectivity/registries are up to date
    hex1.rebuild_connectivity_active_bases();

    // destination field (nodal)
    Field& fld1 = hex1.field(fid1);
    const BasisT<DIM> b = to_basis<DIM>(fld1.basis_id);

    // unique destination nodes (global)
    const auto& nodes = hex1.basis_nodes(b);
    fld1.nodal.resize(nodes.size());

    for (size_t gid = 0; gid < nodes.size(); ++gid) {
      const Point<DIM> s1 = nodes[gid].parent; // (ξ,η,ζ) at t^{n+1} in parent

      // velocity evaluator (uses hex0 geometry/Jacobian)
      auto vel_eval = [&](const Point<DIM> & s, double t) {
        return eval_velocity_parent(hex0, velBasis, vOld, vNew, s, t, dt);
      };

      // backward RK4 in parent: s1 -> s0
      const Point<DIM> s0 = rk4_advect_parent(s1, -dt, vel_eval);

      // sample source field on hex0 directly in parent space
      double val = std::numeric_limits<double>::quiet_NaN();
      (void)hex0.evaluate_field_on_parent(fid0, s0, val);

      fld1.nodal[gid] = val;
    }
  }

} // namespace fem
