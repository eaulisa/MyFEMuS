 
// HexTree.hpp
#pragma once
#include "OctTree.hpp"  // your file above

namespace fem {

class HexTree : public OctTree<3, HexTree> {
public:
  using Base   = OctTree<3, HexTree>;
  using Point3 = Point<3>;

  // Inherit all OctTree<3> constructors
  using Base::Base;

  // Convenience aliases that read nicer at call sites
  // using Base::get_X;
  // using Base::get_Y;
  // using Base::get_Z;
  using Base::parent_to_physical;
  using Base::inverse_map_hex27;       // if you want to expose it with the base name
  using Base::extract_leaf_parent_coords;
  using Base::leaf_physical_nodes;
  using Base::write_binary_vtu;

  // inline uint64_t interleave(u32Point<3> x) noexcept {
  //   uint64_t xx = expand3_21(x[0]);
  //   uint64_t yy = expand3_21(x[1]) << 1;
  //   uint64_t zz = expand3_21(x[2]) << 2;
  //   return (zz | yy | xx);
  // }

  // Optional: friendlier, HexTree-named wrappers
  // void set_geometry_hex27(const double X27[27], const double Y27[27], const double Z27[27]) {
  //   this->set_physical_hex27(X27, Y27, Z27);
  // }
  // bool inverse_map(const Point3& x, Point3& s, int it=30, double tol=1e-12) const {
  //   return this->inverse_map_hex27(x, s, it, tol);
  // }

  // // ---- 3D-specific helper you asked about earlier: s_dot = J^{-1}(s) * v(x(s), t) ----
  // template<class AnalyticVel3D>
  // inline Point3 eval_velocity_parent_analytic(const Point3& s,
  //                                             double t_abs,
  //                                             AnalyticVel3D&& vfun) const {
  //   // -- Physical coordinates via H27
  //   double N27[27];
  //   Shapes3::H27(s, N27);
  //
  //   const auto& X = this->get_X();
  //   const auto& Y = this->get_Y();
  //   const auto& Z = this->get_Z();
  //
  //   double x = 0.0, y = 0.0, z = 0.0;
  //   for (int a = 0; a < 27; ++a) {
  //     x += N27[a] * X[a];
  //     y += N27[a] * Y[a];
  //     z += N27[a] * Z[a];
  //   }
  //
  //   // -- Physical velocity
  //   const Point3 v = vfun(x, y, z, t_abs); // {vx, vy, vz}
  //
  //   // -- Jacobian J and inverse via adjugate
  //   double dNdxi[27], dNdeta[27], dNdz[27];
  //   Shapes3::H27_dN(s, dNdxi, dNdeta, dNdz);
  //
  //   double A11=0, A12=0, A13=0,
  //          A21=0, A22=0, A23=0,
  //          A31=0, A32=0, A33=0;
  //
  //   for (int a = 0; a < 27; ++a) {
  //     const double Xa = X[a], Ya = Y[a], Za = Z[a];
  //     A11 += dNdxi[a]*Xa;  A12 += dNdeta[a]*Xa;  A13 += dNdz[a]*Xa;
  //     A21 += dNdxi[a]*Ya;  A22 += dNdeta[a]*Ya;  A23 += dNdz[a]*Ya;
  //     A31 += dNdxi[a]*Za;  A32 += dNdeta[a]*Za;  A33 += dNdz[a]*Za;
  //   }
  //
  //   const double det =
  //     A11 * (A22*A33 - A23*A32) -
  //     A12 * (A21*A33 - A23*A31) +
  //     A13 * (A21*A32 - A22*A31);
  //   assert(std::abs(det) > 1e-20 && "singular hex mapping Jacobian");
  //   const double inv = 1.0 / det;
  //
  //   Point3 sdot;
  //   sdot[0] = ((A22*A33 - A23*A32)*v[0] - (A12*A33 - A13*A32)*v[1] + (A12*A23 - A13*A22)*v[2]) * inv;
  //   sdot[1] = (-(A21*A33 - A23*A31)*v[0] + (A11*A33 - A13*A31)*v[1] - (A11*A23 - A13*A21)*v[2]) * inv;
  //   sdot[2] = ((A21*A32 - A22*A31)*v[0] - (A11*A32 - A12*A31)*v[1] + (A11*A22 - A12*A21)*v[2]) * inv;
  //   return sdot;
  // }
};

} // namespace fem
