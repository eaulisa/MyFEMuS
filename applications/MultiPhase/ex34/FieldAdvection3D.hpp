#pragma once
#include "OctTree.hpp"
#include <array>
#include <vector>
#include <cassert>
#include <limits>
#include <algorithm>
#include <cmath>

namespace fem {

// -------------------------------
// Helpers: shape dispatch by DIM
// -------------------------------
  template<std::size_t DIM>
  struct ShapeDispatch;

  template<>
  struct ShapeDispatch<2> {
    static inline int eval_N(Basis2D b, const Point2& s, double* N) {
      switch (b) {
      case Basis2D::Q4: Shapes2D::Q4(s, N); return 4;
      case Basis2D::Q8: Shapes2D::Q8(s, N); return 8;
      case Basis2D::Q9: Shapes2D::Q9(s, N); return 9;
      default: return 0;
      }
    }
    static inline int eval_dN_geom(const Point2& s, double* dNdxi, double* dNdeta) {
      Shapes2D::Q9_dN(s, dNdxi, dNdeta); return 9; // geometry via Q9
    }
  };

  template<>
  struct ShapeDispatch<3> {
    static inline int eval_N(Basis3D b, const Point3& s, double* N) {
      switch (b) {
      case Basis3D::H8:  Shapes3D::H8(s,  N); return 8;
      case Basis3D::H20: Shapes3D::H20(s, N); return 20;
      case Basis3D::H27: Shapes3D::H27(s, N); return 27;
      default: return 0;
      }
    }
    static inline int eval_dN_geom(const Point3& s, double* dNdxi, double* dNdeta, double* dNdz) {
      Shapes3D::H27_dN(s, dNdxi, dNdeta, dNdz); return 27; // geometry via H27
    }
  };

// ----------------------------------------
// 1) Evaluate velocity in parent coords s
// vOld/vNew: nodal velocities at 'velBasis' nodes (size = #basis nodes)
// t in [0,dt] for forward, and in [0,dt] while dt<0 for backward (see theta)
// ----------------------------------------
  template<std::size_t DIM>
  inline Point<DIM> eval_velocity_parent(
    const OctTree<DIM>& tree,
    BasisT<DIM> velBasis,
    const std::vector<Point<DIM>>& vOld,
    const std::vector<Point<DIM>>& vNew,
    const Point<DIM>& s_parent,  // (ξ,η[,ζ])
    double t, double dt) {
    // ---- shapes at s for the selected velocity basis ----
    double Nbuf[27]; // enough for Q9/H27
    int nBasis = 0;
    if (DIM == 2) {
      const auto b = static_cast<Basis2D>(velBasis);
      Point2 s{ s_parent[0], s_parent[1] };
      nBasis = ShapeDispatch<2>::eval_N(b, s, Nbuf);
    }
    else {
      const auto b = static_cast<Basis3D>(velBasis);
      Point3 s{ s_parent[0], s_parent[1], s_parent[2] };
      nBasis = ShapeDispatch<3>::eval_N(b, s, Nbuf);
    }
    assert(nBasis > 0 && "Unsupported velocity basis/DIM combo");

    // ---- time interpolation weight θ ----
    // Forward (dt>0): t goes 0→dt ⇒ θ = t/dt in [0,1]
    // Backward (dt<0): t goes 0→dt (negative) ⇒ θ = 1 + t/dt in [1,0]
    double theta;
    if (dt > 0.0)      theta = t / dt;
    else if (dt < 0.0) theta = 1.0 + t / dt;
    else               theta = 0.0;
    theta = std::max(0.0, std::min(1.0, theta));

    // ---- interpolate physical velocity v(x(s), t) ----
    Point<DIM> v{}; // zero
    for (int a = 0; a < nBasis; ++a) {
      Point<DIM> va{};
      for (std::size_t k = 0; k < DIM; ++k)
        va[k] = (1.0 - theta) * vOld[a][k] + theta * vNew[a][k];
      const double Na = Nbuf[a];
      for (std::size_t k = 0; k < DIM; ++k) v[k] += Na * va[k];
    }

    // ---- geometry Jacobian J(s) from Q9/H27 mapping ----
    double J[DIM][DIM];
    for (std::size_t i = 0; i < DIM; ++i)
      for (std::size_t j = 0; j < DIM; ++j)
        J[i][j] = 0.0;

    if (DIM == 2) {
      double dNdxi[9], dNdeta[9];
      Point2 s{ s_parent[0], s_parent[1] };
      const int ngeom = ShapeDispatch<2>::eval_dN_geom(s, dNdxi, dNdeta);
      (void)ngeom; // 9
      const auto& X = tree.get_X();
      const auto& Y = tree.get_Y();
      for (int a = 0; a < 9; ++a) {
        J[0][0] += dNdxi[a]  * X[a];  J[0][1] += dNdeta[a] * X[a];
        J[1][0] += dNdxi[a]  * Y[a];  J[1][1] += dNdeta[a] * Y[a];
      }
    }
    else {
      double dNdxi[27], dNdeta[27], dNdz[27];
      Point3 s{ s_parent[0], s_parent[1], s_parent[2] };
      const int ngeom = ShapeDispatch<3>::eval_dN_geom(s, dNdxi, dNdeta, dNdz);
      (void)ngeom; // 27
      const auto& X = tree.get_X();
      const auto& Y = tree.get_Y();
      const auto& Z = tree.get_Z();
      for (int a = 0; a < 27; ++a) {
        J[0][0] += dNdxi[a]  * X[a];  J[0][1] += dNdeta[a] * X[a];  J[0][2] += dNdz[a] * X[a];
        J[1][0] += dNdxi[a]  * Y[a];  J[1][1] += dNdeta[a] * Y[a];  J[1][2] += dNdz[a] * Y[a];
        J[2][0] += dNdxi[a]  * Z[a];  J[2][1] += dNdeta[a] * Z[a];  J[2][2] += dNdz[a] * Z[a];
      }
    }

    // ---- ṡ = J^{-1} v ----
    Point<DIM> sdot{};
    if (DIM == 2) {
      const double A11 = J[0][0], A12 = J[0][1];
      const double A21 = J[1][0], A22 = J[1][1];
      const double det = A11 * A22 - A12 * A21;
      assert(std::fabs(det) > 1e-20 && "Singular 2D mapping Jacobian");
      const double inv = 1.0 / det;
      sdot[0] = (A22 * v[0] - A12 * v[1]) * inv;
      sdot[1] = (-A21 * v[0] + A11 * v[1]) * inv;
    }
    else {
      const double A11 = J[0][0], A12 = J[0][1], A13 = J[0][2];
      const double A21 = J[1][0], A22 = J[1][1], A23 = J[1][2];
      const double A31 = J[2][0], A32 = J[2][1], A33 = J[2][2];
      const double det =
        A11 * (A22 * A33 - A23 * A32)
        - A12 * (A21 * A33 - A23 * A31)
        + A13 * (A21 * A32 - A22 * A31);
      assert(std::fabs(det) > 1e-20 && "Singular 3D mapping Jacobian");
      const double inv = 1.0 / det;
      sdot[0] = ((A22 * A33 - A23 * A32) * v[0]
                 - (A12 * A33 - A13 * A32) * v[1]
                 + (A12 * A23 - A13 * A22) * v[2]) * inv;
      sdot[1] = (-(A21 * A33 - A23 * A31) * v[0]
                 + (A11 * A33 - A13 * A31) * v[1]
                 - (A11 * A23 - A13 * A21) * v[2]) * inv;
      sdot[2] = ((A21 * A32 - A22 * A31) * v[0]
                 - (A11 * A32 - A12 * A31) * v[1]
                 + (A11 * A22 - A12 * A21) * v[2]) * inv;
    }
    return sdot;
  }

// ----------------------------------------
// 2) RK4 in parent space (DIM = 2 or 3)
// vel_eval(s, tlocal) returns ṡ
// ----------------------------------------
  template<std::size_t DIM, class VelEval>
  inline Point<DIM> rk4_advect_parent(const Point<DIM>& s0, double dt, VelEval&& vel_eval) {
    auto add = [](const Point<DIM>& a, const Point<DIM>& b, double s) {
      Point<DIM> r = a; for (std::size_t k = 0; k < DIM; ++k) r[k] += s * b[k]; return r;
    };
    const Point<DIM> k1 = vel_eval(s0,                 0.0);
    const Point<DIM> k2 = vel_eval(add(s0, k1, 0.5 * dt), 0.5 * dt);
    const Point<DIM> k3 = vel_eval(add(s0, k2, 0.5 * dt), 0.5 * dt);
    const Point<DIM> k4 = vel_eval(add(s0, k3,     dt),     dt);

    Point<DIM> s = s0;
    for (std::size_t k = 0; k < DIM; ++k)
      s[k] += (dt / 6.0) * (k1[k] + 2.0 * k2[k] + 2.0 * k3[k] + k4[k]);
    return s;
  }

// ----------------------------------------
// 3) Forward advection of markers (parent RK, DIM=2/3)
// Emits physical coords that left (at t^n) and those that stayed (at t^{n+1})
// ----------------------------------------
  template<std::size_t DIM>
  inline void advect_markers_forward(
    const OctTree<DIM>& tree0,
    BasisT<DIM> velBasis,
    const std::vector<Point<DIM>>& vOld,   // sized to velBasis nodes
    const std::vector<Point<DIM>>& vNew,   // sized to velBasis nodes
    double dt,
    std::vector<Point<DIM>>& coordsLeftOld,
    std::vector<Point<DIM>>& coordsStayedNew) {
    coordsLeftOld.clear();
    coordsStayedNew.clear();

    // seed: parent nodes on finest current level, using highest-order geometry basis (id=2)
    const auto s0_all =
      tree0.extract_node_parent_coords_in_level_range(
        tree0.max_depth(), tree0.max_depth(),
        static_cast<BasisT<DIM>>(2));

    for (const auto& s0 : s0_all) {
      auto vel_eval = [&](const Point<DIM>& s, double tl) {
        return eval_velocity_parent(tree0, velBasis, vOld, vNew, s, tl, dt);
      };
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

// ----------------------------------------
// 4) Backward advection & field transport (parent RK, DIM=2/3)
// Samples source field in parent space on tree0 into tree1
// ----------------------------------------
  template<std::size_t DIM>
  inline void advect_nodes_backward_and_transport_field(
    const OctTree<DIM>& tree0,   // geometry + source at old time
    u32 fid0,                    // source field id in tree0
    BasisT<DIM> velBasis,
    const std::vector<Point<DIM>>& vOld,
    const std::vector<Point<DIM>>& vNew,
    double dt,                   // > 0 (we will integrate with -dt)
    OctTree<DIM>& tree1,         // geometry at end time
    u32 fid1) {                  // destination field id in tree1
    tree1.rebuild_connectivity_active_bases();

    Field& dst = tree1.field(fid1);
    const auto b = to_basis<DIM>(dst.basis_id);
    const auto& nodes = tree1.basis_nodes(b);
    dst.nodal.resize(nodes.size());

    for (size_t gid = 0; gid < nodes.size(); ++gid) {
      const Point<DIM> s1 = nodes[gid].parent; // (ξ,η[,ζ]) at t^{n+1}

      auto vel_eval = [&](const Point<DIM>& s, double tl) {
        // negative step handled inside RK via tl ∈ [0,-dt]; eval uses theta logic
        return eval_velocity_parent(tree0, velBasis, vOld, vNew, s, tl, -dt);
      };
      const Point<DIM> s0 = rk4_advect_parent<DIM>(s1, -dt, vel_eval);

      double val = std::numeric_limits<double>::quiet_NaN();
      (void)tree0.evaluate_field_on_parent(fid0, s0, val);
      dst.nodal[gid] = val;
    }
  }

} // namespace fem
