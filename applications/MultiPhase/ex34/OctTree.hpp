#pragma once
#include <vector>
#include <array>
#include <unordered_map>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <unordered_set>
#include <cstdint>
#include <functional>
#include <utility>
#include <deque>
#include <type_traits>
#include <variant>
#include "Encoder.hpp"
#include "Mollifier.hpp"

#include "Reinitializer.hpp"

namespace fem {

  using u32 = uint32_t;
  using u64 = uint64_t;
  constexpr u32 npos32 = std::numeric_limits<u32>::max();

  template<std::size_t DIM>
  using Point = std::array<double, DIM>;

  using Point2 = std::array<double, 2>;
  using Point3 = std::array<double, 3>;

  // Build a Point<DIM> from (x,y,z) consistently in 2D/3D.
  template<std::size_t DIM> struct PtMake;

  template<> struct PtMake<2> {
    static inline Point<2> mk(double x, double y, double /*z*/) {
      return Point<2> {x, y};
    }
  };

  template<> struct PtMake<3> {
    static inline Point<3> mk(double x, double y, double z) {
      return Point<3> {x, y, z};
    }
  };

  template<std::size_t DIM>
  using u32Point = std::array<u32, DIM>;

  inline constexpr std::size_t NDOFS[4][3] = {
    {0, 0, 0},     // DIM = 0 (unused)
    {0, 0, 0},     // DIM = 1 (unused)
    {4, 8, 9},     // DIM = 2 -> Q1, Serendipity8, Q2
    {8, 20, 27}    // DIM = 3 -> H8, H20, H27
  };

  // ==== Basis enums (outside any class) ====
  enum class Basis3D : uint8_t { H8 = 0, H20 = 1, H27 = 2 };
  enum class Basis2D : uint8_t { Q4 = 0, Q8  = 1, Q9  = 2 };

// ==== Dimension-selected basis alias ====
  template <std::size_t DIM>
  using BasisT = typename std::conditional<(DIM == 3), Basis3D, Basis2D>::type;

// ==== Free helpers (tree-agnostic) ====
  template <std::size_t DIM>
  inline BasisT<DIM> to_basis(uint8_t id) noexcept {
    static_assert(DIM == 2 || DIM == 3, "to_basis only supports DIM = 2 or 3");
    assert(id < 3);
    return static_cast<BasisT<DIM>>(id);
  }

  template <std::size_t DIM>
  inline uint8_t to_basis_id(BasisT<DIM> b) noexcept {
    static_assert(DIM == 2 || DIM == 3, "to_basis_id only supports DIM = 2 or 3");
    return static_cast<uint8_t>(b);
  }

// Optional: names for printing/debug
  inline const char* basis_name(Basis2D b) noexcept {
    switch (b) {
    case Basis2D::Q4: return "Q4";
    case Basis2D::Q8: return "Q8";
    case Basis2D::Q9: return "Q9";
    default: return "Unknown2D";
    }
  }
  inline const char* basis_name(Basis3D b) noexcept {
    switch (b) {
    case Basis3D::H8: return "H8";
    case Basis3D::H20: return "H20";
    case Basis3D::H27: return "H27";
    default: return "Unknown3D";
    }
  }
  template <std::size_t DIM>
  inline const char* basis_name(BasisT<DIM> b) noexcept {
    return (DIM == 3) ? basis_name(static_cast<Basis3D>(b))
           : basis_name(static_cast<Basis2D>(b));
  }

// ---------- Morton helpers (3-way) ----------
  static inline uint64_t expand3_21(uint64_t v) {
    v &= 0x1FFFFF;
    v = (v | (v << 32)) & 0x1F00000000FFFFULL;
    v = (v | (v << 16)) & 0x1F0000FF0000FFULL;
    v = (v | (v << 8))  & 0x100F00F00F00F00FULL;
    v = (v | (v << 4))  & 0x10C30C30C30C30C3ULL;
    v = (v | (v << 2))  & 0x1249249249249249ULL;
    return v;
  }

// 2D helper: expand 32 LSBits to 0b a0 0 a1 0 a2 ...
  static inline uint64_t expand2_32(uint64_t v) noexcept {
    v &= 0x00000000FFFFFFFFull;
    v = (v | (v << 16)) & 0x0000FFFF0000FFFFull;
    v = (v | (v << 8)) & 0x00FF00FF00FF00FFull;
    v = (v | (v << 4)) & 0x0F0F0F0F0F0F0F0Full;
    v = (v | (v << 2)) & 0x3333333333333333ull;
    v = (v | (v << 1)) & 0x5555555555555555ull;  // 0101...
    return v;
  }

// Tag-dispatch impls (C++14-friendly)
  template<std::size_t D>
  inline uint64_t interleave_impl(const u32Point<D>& x, std::integral_constant<std::size_t, 2>) noexcept {
    // Up to 32 bits per axis (2*32 = 64)
    const uint64_t xx = expand2_32(x[0]);
    const uint64_t yy = expand2_32(x[1]) << 1;
    return xx | yy;
  }

  template<std::size_t D>
  inline uint64_t interleave_impl(const u32Point<D>& x, std::integral_constant<std::size_t, 3>) noexcept {
    // Uses your existing expand3_21 (3 * 21 = 63 bits)
    const uint64_t xx = expand3_21(x[0]);
    const uint64_t yy = expand3_21(x[1]) << 1;
    const uint64_t zz = expand3_21(x[2]) << 2;
    return xx | yy | zz;
  }

  struct Shapes2D {
    // ---------- 1D Q2 Lagrange ----------
    static inline void q2_1d_vals(double s, double& L0, double& L1, double& L2) noexcept {
      const double s2 = s * s;
      L0 = 0.5 * (s2 - s);   // 0.5*s*(s-1)  -> node at s=-1
      L1 = 1.0 - s2;         //              -> node at s= 0
      L2 = 0.5 * (s2 + s);   // 0.5*s*(s+1)  -> node at s=+1
    }
    static inline void q2_1d_derivs(double s, double& dL0, double& dL1, double& dL2) noexcept {
      dL0 = s - 0.5;
      dL1 = -2.0 * s;
      dL2 = s + 0.5;
    }

    // ============================================================
    // Q4: bilinear quad (4 nodes: corners only, indices 0..3)
    // Corner order: (ξ,η) = (-1,-1),(+1,-1),(+1,+1),(-1,+1)
    // ============================================================
    static inline void Q4(const Point2& s, double* __restrict__ N) noexcept {
      const double xi = s[0], eta = s[1];
      const double a = (1.0 - xi), b = 1.0 + xi;
      const double c = 1.0 - eta,  d = 1.0 + eta;

      N[0] = 0.25 * a * c;            // 0.25*(1-ξ)*(1-η)
      N[1] = 0.25 * b * c;            // 0.25*(1+ξ)*(1-η)
      N[2] = 0.25 * b * d;            // 0.25*(1+ξ)*(1+η)
      N[3] = 0.25 * a * d;            // 0.25*(1-ξ)*(1+η)
    }
    static inline void Q4_dN(const Point2& s,
                             double* __restrict__ dNdxi,
                             double* __restrict__ dNdeta) noexcept {
      const double xi = s[0], eta = s[1];
      // dN/dξ
      dNdxi[0] = -0.25 * (1.0 - eta);
      dNdxi[1] = +0.25 * (1.0 - eta);
      dNdxi[2] = +0.25 * (1.0 + eta);
      dNdxi[3] = -0.25 * (1.0 + eta);
      // dN/dη
      dNdeta[0] = -0.25 * (1.0 - xi);
      dNdeta[1] = -0.25 * (1.0 + xi);
      dNdeta[2] = +0.25 * (1.0 + xi);
      dNdeta[3] = +0.25 * (1.0 - xi);

#ifndef NDEBUG
      {
        constexpr double tol = 1e-12;
        double sx = 0.0, se = 0.0;
        for (int i = 0; i < 4; ++i) {
          sx += dNdxi[i];
          se += dNdeta[i];
        }
        assert(std::fabs(sx) < tol && "Q4: sum(dN/dxi) must be 0");
        assert(std::fabs(se) < tol && "Q4: sum(dN/deta) must be 0");
      }
#endif

    }

    // ============================================================
    // Q9: biquadratic Lagrange (9 nodes)
    // Indices:
    //   0..3  corners  (-1,-1),(+1,-1),(+1,+1),(-1,+1)
    //   4..7  edge mids: bottom(ξ-mid,η=-1), right(ξ=+1,η-mid),
    //                    top(ξ-mid,η=+1),    left(ξ=-1,η-mid)
    //   8     center (ξ=0,η=0)
    // ============================================================
    static inline void Q9(const Point2& s, double* __restrict__ N) noexcept {
      double Lx0, Lx1, Lx2, Ly0, Ly1, Ly2;
      q2_1d_vals(s[0], Lx0, Lx1, Lx2);
      q2_1d_vals(s[1], Ly0, Ly1, Ly2);

      auto S = [&](double Lx, double Ly) {
        return Lx * Ly;
      };

      // corners 0..3
      N[0] = S(Lx0, Ly0);
      N[1] = S(Lx2, Ly0);
      N[2] = S(Lx2, Ly2);
      N[3] = S(Lx0, Ly2);

      // edges 4..7 (bottom, right, top, left)
      N[4] = S(Lx1, Ly0); // bottom
      N[5] = S(Lx2, Ly1); // right
      N[6] = S(Lx1, Ly2); // top
      N[7] = S(Lx0, Ly1); // left

      // center 8
      N[8] = S(Lx1, Ly1);
    }
    static inline void Q9_dN(const Point2& s,
                             double* __restrict__ dNdxi,
                             double* __restrict__ dNdeta) noexcept {
      double Lx0, Lx1, Lx2, Ly0, Ly1, Ly2;
      double dLx0, dLx1, dLx2, dLy0, dLy1, dLy2;
      q2_1d_vals(s[0], Lx0, Lx1, Lx2); q2_1d_derivs(s[0], dLx0, dLx1, dLx2);
      q2_1d_vals(s[1], Ly0, Ly1, Ly2); q2_1d_derivs(s[1], dLy0, dLy1, dLy2);

      auto fill = [&](int i, double Lx, double Ly, double dLx, double dLy) {
        dNdxi[i]  = dLx * Ly;
        dNdeta[i] = Lx  * dLy;
      };

      // corners
      fill(0, Lx0, Ly0, dLx0, dLy0);
      fill(1, Lx2, Ly0, dLx2, dLy0);
      fill(2, Lx2, Ly2, dLx2, dLy2);
      fill(3, Lx0, Ly2, dLx0, dLy2);

      // edges (bottom, right, top, left)
      fill(4, Lx1, Ly0, dLx1, dLy0);
      fill(5, Lx2, Ly1, dLx2, dLy1);
      fill(6, Lx1, Ly2, dLx1, dLy2);
      fill(7, Lx0, Ly1, dLx0, dLy1);

      // center
      fill(8, Lx1, Ly1, dLx1, dLy1);

#ifndef NDEBUG
      {
        constexpr double tol = 1e-12;
        double sx = 0.0, se = 0.0;
        for (int i = 0; i < 9; ++i) {
          sx += dNdxi[i];
          se += dNdeta[i];
        }
        assert(std::fabs(sx) < tol && "Q9: sum(dN/dxi) must be 0");
        assert(std::fabs(se) < tol && "Q9: sum(dN/deta) must be 0");
      }
#endif
    }

    // ============================================================
    // Q8: 8-node serendipity (no center node)
    // Indices:
    //   0..3 corners as above
    //   4..7 edge mids: bottom, right, top, left
    // ============================================================
    static inline void Q8(const Point2& s, double* __restrict__ N) noexcept {
      const double xi = s[0], eta = s[1];

      // corners
      N[0] = 0.25 * (1.0 - xi) * (1.0 - eta) * (-xi - eta - 1.0);
      N[1] = 0.25 * (1.0 + xi) * (1.0 - eta) * (xi - eta - 1.0);
      N[2] = 0.25 * (1.0 + xi) * (1.0 + eta) * (xi + eta - 1.0);
      N[3] = 0.25 * (1.0 - xi) * (1.0 + eta) * (-xi + eta - 1.0);

      // edge mids (bottom, right, top, left)
      N[4] = 0.5 * (1.0 - xi * xi) * (1.0 - eta);   // bottom (η=-1)
      N[5] = 0.5 * (1.0 + xi)     * (1.0 - eta * eta); // right  (ξ=+1)
      N[6] = 0.5 * (1.0 - xi * xi) * (1.0 + eta);   // top    (η=+1)
      N[7] = 0.5 * (1.0 - xi)     * (1.0 - eta * eta); // left   (ξ=-1)
    }

    // Q8 serendipity derivatives (node order: 0..3 corners, 4..7 edges)
    static inline void Q8_dN(const Point2& s,
                             double* __restrict__ dNdxi,
                             double* __restrict__ dNdeta) noexcept {
      const double xi  = s[0];
      const double eta = s[1];

      // Corners: N = 1/4 (1+sx*xi)(1+sy*eta)(sx*xi + sy*eta - 1)
      // dN/dxi  = 1/4 * sx * (1+sy*eta) * [ (sx*xi + sy*eta - 1) + (1+sx*xi) ]
      // dN/deta = 1/4 * sy * (1+sx*xi)  * [ (sx*xi + sy*eta - 1) + (1+sy*eta) ]
      auto dN_corner = [&](int sx, int sy, int i) {
        const double a = 1.0 + sx * xi;
        const double b = 1.0 + sy * eta;
        const double c = sx * xi + sy * eta - 1.0;
        dNdxi[i]  = 0.25 * sx * b * (c + a);
        dNdeta[i] = 0.25 * sy * a * (c + b);
      };

      // corners 0..3: (-,-), (+,-), (+,+), (-,+)
      dN_corner(-1, -1, 0);
      dN_corner(+1, -1, 1);
      dN_corner(+1, +1, 2);
      dN_corner(-1, +1, 3);

      // Edges (serendipity)
      // 4 bottom:  1/2 (1 - xi^2)(1 - eta)
      dNdxi[4]  = -xi * (1.0 - eta);
      dNdeta[4] = -0.5 * (1.0 - xi * xi);

      // 5 right:   1/2 (1 + xi)(1 - eta^2)
      dNdxi[5]  =  0.5 * (1.0 - eta * eta);
      dNdeta[5] = -(1.0 + xi) * eta;

      // 6 top:     1/2 (1 - xi^2)(1 + eta)
      dNdxi[6]  = -xi * (1.0 + eta);
      dNdeta[6] =  0.5 * (1.0 - xi * xi);

      // 7 left:    1/2 (1 - xi)(1 - eta^2)
      dNdxi[7]  = -0.5 * (1.0 - eta * eta);
      dNdeta[7] = -(1.0 - xi) * eta;

#ifndef NDEBUG
      // Partition of unity ⇒ ∑ dN/dξ = 0, ∑ dN/dη = 0 (use a loose tolerance)
      double sxsum = 0.0, setasum = 0.0;
      for (int i = 0; i < 8; ++i) {
        sxsum += dNdxi[i];
        setasum += dNdeta[i];
      }
      assert(std::fabs(sxsum)   < 1e-12 && "Q8: sum(dN/dxi) must be 0");
      assert(std::fabs(setasum) < 1e-12 && "Q8: sum(dN/deta) must be 0");
#endif
    }
  };

  struct Shapes3D {
    // 1D Q2 Lagrange
    static inline void q2_1d_vals(double s, double& L0, double& L1, double& L2) noexcept {
      const double s2 = s * s;
      L0 = 0.5 * (s2 - s); // 0.5*s*(s-1)
      L1 = 1.0 - s2;
      L2 = 0.5 * (s2 + s); // 0.5*s*(s+1)
    }
    static inline void q2_1d_derivs(double s, double& dL0, double& dL1, double& dL2) noexcept {
      dL0 = s - 0.5;
      dL1 = -2.0 * s;
      dL2 = s + 0.5;
    }

    // ---------- H8 (trilinear) ----------
    static inline void H8(const Point3 &s, double* __restrict__ N) noexcept {
      const double xi = s[0], eta = s[1], z = s[2];
      const double a0 = 0.125 * (1.0 - xi) * (1.0 - eta);
      const double a1 = 0.125 * (1.0 + xi) * (1.0 - eta);
      const double a2 = 0.125 * (1.0 + xi) * (1.0 + eta);
      const double a3 = 0.125 * (1.0 - xi) * (1.0 + eta);

      N[0] = a0 * (1.0 - z);
      N[1] = a1 * (1.0 - z);
      N[2] = a2 * (1.0 - z);
      N[3] = a3 * (1.0 - z);
      N[4] = a0 * (1.0 + z);
      N[5] = a1 * (1.0 + z);
      N[6] = a2 * (1.0 + z);
      N[7] = a3 * (1.0 + z);
    }

    static inline void H8_dN(const Point3& s,
                             double* __restrict__ dNdxi,
                             double* __restrict__ dNdeta,
                             double* __restrict__ dNdz) noexcept {
      const double xi = s[0], eta = s[1], z = s[2];

      const double a0 = 0.125 * (1 - eta) * (1 - z);
      const double a1 = 0.125 * (1 + eta) * (1 - z);
      const double a2 = 0.125 * (1 - eta) * (1 + z);
      const double a3 = 0.125 * (1 + eta) * (1 + z);
      dNdxi[0] = -a0; dNdxi[1] =  a0; dNdxi[2] =  a1; dNdxi[3] = -a1;
      dNdxi[4] = -a2; dNdxi[5] =  a2; dNdxi[6] =  a3; dNdxi[7] = -a3;

      const double b0 = 0.125 * (1 - xi) * (1 - z);
      const double b1 = 0.125 * (1 + xi) * (1 - z);
      const double b2 = 0.125 * (1 - xi) * (1 + z);
      const double b3 = 0.125 * (1 + xi) * (1 + z);
      dNdeta[0] = -b0; dNdeta[1] = -b1; dNdeta[2] =  b1; dNdeta[3] =  b0;
      dNdeta[4] = -b2; dNdeta[5] = -b3; dNdeta[6] =  b3; dNdeta[7] =  b2;

      const double c0 = 0.125 * (1 - xi) * (1 - eta);
      const double c1 = 0.125 * (1 + xi) * (1 - eta);
      const double c2 = 0.125 * (1 + xi) * (1 + eta);
      const double c3 = 0.125 * (1 - xi) * (1 + eta);
      dNdz[0] = -c0; dNdz[1] = -c1; dNdz[2] = -c2; dNdz[3] = -c3;
      dNdz[4] =  c0; dNdz[5] =  c1; dNdz[6] =  c2; dNdz[7] =  c3;

#ifndef NDEBUG
      {
        constexpr double tol = 1e-12;
        double sx = 0.0, se = 0.0, sz = 0.0;
        for (int i = 0; i < 8; ++i) {
          sx += dNdxi[i];
          se += dNdeta[i];
          sz += dNdz[i];
        }
        assert(std::fabs(sx) < tol && "H8: sum(dN/dxi) must be 0");
        assert(std::fabs(se) < tol && "H8: sum(dN/deta) must be 0");
        assert(std::fabs(sz) < tol && "H8: sum(dN/dz)   must be 0");
      }
#endif

    }

    // ---------- H27 (triquadratic) ----------
    static inline void H27(const Point3& s, double* __restrict__ N) noexcept {
      double Lx0, Lx1, Lx2, Ly0, Ly1, Ly2, Lz0, Lz1, Lz2;
      q2_1d_vals(s[0], Lx0, Lx1, Lx2);
      q2_1d_vals(s[1], Ly0, Ly1, Ly2);
      q2_1d_vals(s[2], Lz0, Lz1, Lz2);

      auto S = [&](double Lx, double Ly, double Lz) {
        return Lx * Ly * Lz;
      };

      // corners (0..7)
      N[0] = S(Lx0, Ly0, Lz0); N[1] = S(Lx2, Ly0, Lz0); N[2] = S(Lx2, Ly2, Lz0); N[3] = S(Lx0, Ly2, Lz0);
      N[4] = S(Lx0, Ly0, Lz2); N[5] = S(Lx2, Ly0, Lz2); N[6] = S(Lx2, Ly2, Lz2); N[7] = S(Lx0, Ly2, Lz2);

      // z=-1 edges (8..11)
      N[ 8] = S(Lx1, Ly0, Lz0); N[ 9] = S(Lx2, Ly1, Lz0);
      N[10] = S(Lx1, Ly2, Lz0); N[11] = S(Lx0, Ly1, Lz0);

      // z=+1 edges (12..15)
      N[12] = S(Lx1, Ly0, Lz2); N[13] = S(Lx2, Ly1, Lz2);
      N[14] = S(Lx1, Ly2, Lz2); N[15] = S(Lx0, Ly1, Lz2);

      // vertical edges (16..19)
      N[16] = S(Lx0, Ly0, Lz1); N[17] = S(Lx2, Ly0, Lz1);
      N[18] = S(Lx2, Ly2, Lz1); N[19] = S(Lx0, Ly2, Lz1);

      // face centers (20..25)
      N[20] = S(Lx1, Ly0, Lz1); N[21] = S(Lx2, Ly1, Lz1);
      N[22] = S(Lx1, Ly2, Lz1); N[23] = S(Lx0, Ly1, Lz1);
      N[24] = S(Lx1, Ly1, Lz0); N[25] = S(Lx1, Ly1, Lz2);

      // cell center (26)
      N[26] = S(Lx1, Ly1, Lz1);
    }
    static inline void H27_dN(const Point3& s,
                              double* __restrict__ dNdxi,
                              double* __restrict__ dNdeta,
                              double* __restrict__ dNdz) noexcept {
      double Lx0, Lx1, Lx2, Ly0, Ly1, Ly2, Lz0, Lz1, Lz2;
      double dLx0, dLx1, dLx2, dLy0, dLy1, dLy2, dLz0, dLz1, dLz2;
      q2_1d_vals(s[0], Lx0, Lx1, Lx2); q2_1d_derivs(s[0], dLx0, dLx1, dLx2);
      q2_1d_vals(s[1], Ly0, Ly1, Ly2); q2_1d_derivs(s[1], dLy0, dLy1, dLy2);
      q2_1d_vals(s[2], Lz0, Lz1, Lz2); q2_1d_derivs(s[2], dLz0, dLz1, dLz2);

      auto fill = [&](int i, double Lx, double Ly, double Lz,
      double dLx, double dLy, double dLz) {
        dNdxi[i]  = dLx * Ly * Lz;
        dNdeta[i] = Lx * dLy * Lz;
        dNdz[i]   = Lx * Ly * dLz;
      };

      // corners (0..7)
      fill(0, Lx0, Ly0, Lz0, dLx0, dLy0, dLz0);
      fill(1, Lx2, Ly0, Lz0, dLx2, dLy0, dLz0);
      fill(2, Lx2, Ly2, Lz0, dLx2, dLy2, dLz0);
      fill(3, Lx0, Ly2, Lz0, dLx0, dLy2, dLz0);
      fill(4, Lx0, Ly0, Lz2, dLx0, dLy0, dLz2);
      fill(5, Lx2, Ly0, Lz2, dLx2, dLy0, dLz2);
      fill(6, Lx2, Ly2, Lz2, dLx2, dLy2, dLz2);
      fill(7, Lx0, Ly2, Lz2, dLx0, dLy2, dLz2);

      // z=-1 edges (8..11)
      fill(8, Lx1, Ly0, Lz0, dLx1, dLy0, dLz0);
      fill(9, Lx2, Ly1, Lz0, dLx2, dLy1, dLz0);
      fill(10, Lx1, Ly2, Lz0, dLx1, dLy2, dLz0);
      fill(11, Lx0, Ly1, Lz0, dLx0, dLy1, dLz0);

      // z=+1 edges (12..15)
      fill(12, Lx1, Ly0, Lz2, dLx1, dLy0, dLz2);
      fill(13, Lx2, Ly1, Lz2, dLx2, dLy1, dLz2);
      fill(14, Lx1, Ly2, Lz2, dLx1, dLy2, dLz2);
      fill(15, Lx0, Ly1, Lz2, dLx0, dLy1, dLz2);

      // vertical edges (16..19)
      fill(16, Lx0, Ly0, Lz1, dLx0, dLy0, dLz1);
      fill(17, Lx2, Ly0, Lz1, dLx2, dLy0, dLz1);
      fill(18, Lx2, Ly2, Lz1, dLx2, dLy2, dLz1);
      fill(19, Lx0, Ly2, Lz1, dLx0, dLy2, dLz1);

      // face centers (20..25)
      fill(20, Lx1, Ly0, Lz1, dLx1, dLy0, dLz1);
      fill(21, Lx2, Ly1, Lz1, dLx2, dLy1, dLz1);
      fill(22, Lx1, Ly2, Lz1, dLx1, dLy2, dLz1);
      fill(23, Lx0, Ly1, Lz1, dLx0, dLy1, dLz1);
      fill(24, Lx1, Ly1, Lz0, dLx1, dLy1, dLz0);
      fill(25, Lx1, Ly1, Lz2, dLx1, dLy1, dLz2);

      // cell center (26)
      fill(26, Lx1, Ly1, Lz1, dLx1, dLy1, dLz1);

#ifndef NDEBUG
      {
        constexpr double tol = 1e-12;
        double sx = 0.0, se = 0.0, sz = 0.0;
        for (int i = 0; i < 27; ++i) {
          sx += dNdxi[i];
          se += dNdeta[i];
          sz += dNdz[i];
        }
        assert(std::fabs(sx) < tol && "H27: sum(dN/dxi) must be 0");
        assert(std::fabs(se) < tol && "H27: sum(dN/deta) must be 0");
        assert(std::fabs(sz) < tol && "H27: sum(dN/dz)   must be 0");
      }
#endif

    }

    // ---------- H20 (true serendipity) ----------
// Node order (matches your earlier code):
// 0..7  : corners in (+/-) order
// 8..11 : z = -1 edges   (ξ-mid on bottom face)   : (η,ζ) = (±1, -1)
// 12..15: z = +1 edges   (ξ-mid on top face)      : (η,ζ) = (±1, +1)
// 16..19: vertical edges (ζ-mid on vertical edges): (ξ,η) = corners (±1,±1)
    static inline void H20(const Point3& s, double* __restrict__ N) noexcept {
      const double xi = s[0], eta = s[1], z = s[2];

      // corners (use compact serendipity form)
      auto Nc = [&](int sx, int sy, int sz) -> double {
        // sx,sy,sz ∈ { -1, +1 } indicate the corner (signs of (ξ,η,ζ))
        // Ni = 1/8 (1 + sx*ξ)(1 + sy*η)(1 + sz*ζ) * (sx*ξ + sy*η + sz*ζ - 2)
        const double a = 0.125 * (1.0 + sx * xi) * (1.0 + sy * eta) * (1.0 + sz * z);
        const double b = double(sx) * xi + double(sy) * eta + double(sz) * z - 2.0;
        return a * b;
      };

      // edges (true serendipity; no face/center bubbles)
      auto Ne_x = [&](int sy, int sz) -> double {
        // mid-edge along ξ, with fixed (η=sy, ζ=sz)
        return 0.25 * (1.0 - xi * xi) * (1.0 + sy * eta) * (1.0 + sz * z);
      };
      auto Ne_y = [&](int sx, int sz) -> double {
        // mid-edge along η, with fixed (ξ=sx, ζ=sz)
        return 0.25 * (1.0 - eta * eta) * (1.0 + sx * xi) * (1.0 + sz * z);
      };
      auto Ne_z = [&](int sx, int sy) -> double {
        // mid-edge along ζ, with fixed (ξ=sx, η=sy)
        return 0.25 * (1.0 - z * z) * (1.0 + sx * xi) * (1.0 + sy * eta);
      };

      // corners 0..7
      N[0] = Nc(-1, -1, -1);
      N[1] = Nc(+1, -1, -1);
      N[2] = Nc(+1, +1, -1);
      N[3] = Nc(-1, +1, -1);
      N[4] = Nc(-1, -1, +1);
      N[5] = Nc(+1, -1, +1);
      N[6] = Nc(+1, +1, +1);
      N[7] = Nc(-1, +1, +1);

      // z = -1 edges 8..11 (η = -1,+1 around bottom face; ξ-mid)
      N[ 8] = Ne_x(-1, -1); // (-, -, z=-1)
      N[ 9] = Ne_y(+1, -1); // (ξ=+1, -, z=-1)   (along η)
      N[10] = Ne_x(+1, -1); // ( - , +, z=-1)
      N[11] = Ne_y(-1, -1); // (ξ=-1, +, z=-1)

      // z = +1 edges 12..15 (same pattern, z=+1)
      N[12] = Ne_x(-1, +1);
      N[13] = Ne_y(+1, +1);
      N[14] = Ne_x(+1, +1);
      N[15] = Ne_y(-1, +1);

      // vertical edges 16..19 (ζ-mid, fixed (ξ,η) at the four corners)
      N[16] = Ne_z(-1, -1);
      N[17] = Ne_z(+1, -1);
      N[18] = Ne_z(+1, +1);
      N[19] = Ne_z(-1, +1);
    }

    static inline void H20_dN(const Point3& s,
                              double* __restrict__ dNdxi,
                              double* __restrict__ dNdeta,
                              double* __restrict__ dNdz) noexcept {
      const double xi = s[0], eta = s[1], z = s[2];

      // Corner derivatives
      auto dNc = [&](int sx, int sy, int sz,
      double & dx, double & dy, double & dz) {
        // a = 1/8 (1+sx*ξ)(1+sy*η)(1+sz*ζ), b = sx*ξ + sy*η + sz*ζ - 2
        // dN/dξ = da/dξ*b + a*db/dξ,   da/dξ = 1/8*sx*(1+sy*η)(1+sz*ζ),   db/dξ = sx
        const double a = 0.125 * (1.0 + sx * xi) * (1.0 + sy * eta) * (1.0 + sz * z);
        const double b = double(sx) * xi + double(sy) * eta + double(sz) * z - 2.0;

        const double dax = 0.125 * double(sx) * (1.0 + sy * eta) * (1.0 + sz * z);
        const double day = 0.125 * double(sy) * (1.0 + sx * xi) * (1.0 + sz * z);
        const double daz = 0.125 * double(sz) * (1.0 + sx * xi) * (1.0 + sy * eta);

        dx = dax * b + a * double(sx);
        dy = day * b + a * double(sy);
        dz = daz * b + a * double(sz);
      };

      // Edge derivatives
      auto dNe_x = [&](int sy, int sz,
      double & dx, double & dy, double & dz) {
        // N = 1/4 (1-ξ^2)(1+sy*η)(1+sz*ζ)
        dx = 0.25 * (-2.0 * xi) * (1.0 + sy * eta) * (1.0 + sz * z); // -0.5*xi*(1+sy*η)(1+sz*ζ)
        dy = 0.25 * (1.0 - xi * xi) * double(sy) * (1.0 + sz * z);
        dz = 0.25 * (1.0 - xi * xi) * (1.0 + sy * eta) * double(sz);
      };
      auto dNe_y = [&](int sx, int sz,
      double & dx, double & dy, double & dz) {
        // N = 1/4 (1-η^2)(1+sx*ξ)(1+sz*ζ)
        dx = 0.25 * (1.0 - eta * eta) * double(sx) * (1.0 + sz * z);
        dy = 0.25 * (-2.0 * eta) * (1.0 + sx * xi) * (1.0 + sz * z); // -0.5*eta*(1+sx*ξ)(1+sz*ζ)
        dz = 0.25 * (1.0 - eta * eta) * (1.0 + sx * xi) * double(sz);
      };
      auto dNe_z = [&](int sx, int sy,
      double & dx, double & dy, double & dz) {
        // N = 1/4 (1-ζ^2)(1+sx*ξ)(1+sy*η)
        dx = 0.25 * (1.0 - z * z) * double(sx) * (1.0 + sy * eta);
        dy = 0.25 * (1.0 - z * z) * (1.0 + sx * xi) * double(sy);
        dz = 0.25 * (-2.0 * z) * (1.0 + sx * xi) * (1.0 + sy * eta); // -0.5*z*(1+sx*ξ)(1+sy*η)
      };

      // corners 0..7
      dNc(-1, -1, -1, dNdxi[0], dNdeta[0], dNdz[0]);
      dNc(+1, -1, -1, dNdxi[1], dNdeta[1], dNdz[1]);
      dNc(+1, +1, -1, dNdxi[2], dNdeta[2], dNdz[2]);
      dNc(-1, +1, -1, dNdxi[3], dNdeta[3], dNdz[3]);
      dNc(-1, -1, +1, dNdxi[4], dNdeta[4], dNdz[4]);
      dNc(+1, -1, +1, dNdxi[5], dNdeta[5], dNdz[5]);
      dNc(+1, +1, +1, dNdxi[6], dNdeta[6], dNdz[6]);
      dNc(-1, +1, +1, dNdxi[7], dNdeta[7], dNdz[7]);

      // z = -1 edges 8..11
      dNe_x(-1, -1, dNdxi[ 8], dNdeta[ 8], dNdz[ 8]);
      dNe_y(+1, -1, dNdxi[ 9], dNdeta[ 9], dNdz[ 9]);
      dNe_x(+1, -1, dNdxi[10], dNdeta[10], dNdz[10]);
      dNe_y(-1, -1, dNdxi[11], dNdeta[11], dNdz[11]);

      // z = +1 edges 12..15
      dNe_x(-1, +1, dNdxi[12], dNdeta[12], dNdz[12]);
      dNe_y(+1, +1, dNdxi[13], dNdeta[13], dNdz[13]);
      dNe_x(+1, +1, dNdxi[14], dNdeta[14], dNdz[14]);
      dNe_y(-1, +1, dNdxi[15], dNdeta[15], dNdz[15]);

      // vertical edges 16..19
      dNe_z(-1, -1, dNdxi[16], dNdeta[16], dNdz[16]);
      dNe_z(+1, -1, dNdxi[17], dNdeta[17], dNdz[17]);
      dNe_z(+1, +1, dNdxi[18], dNdeta[18], dNdz[18]);
      dNe_z(-1, +1, dNdxi[19], dNdeta[19], dNdz[19]);


#ifndef NDEBUG
      {
        constexpr double tol = 1e-12;
        double sx = 0.0, se = 0.0, sz = 0.0;
        for (int i = 0; i < 20; ++i) {
          sx += dNdxi[i];
          se += dNdeta[i];
          sz += dNdz[i];
        }
        assert(std::fabs(sx) < tol && "H20: sum(dN/dxi) must be 0");
        assert(std::fabs(se) < tol && "H20: sum(dN/deta) must be 0");
        assert(std::fabs(sz) < tol && "H20: sum(dN/dz)   must be 0");
      }
#endif
    }
  };

// ---- Tree node & helpers (DIM-independent)
  template<std::size_t DIM>
  struct TreeNode {
    TreeNode() {
      static_assert(DIM == 2 || DIM == 3, "TreeNode supports DIM=2 or 3");
      child.fill(npos32);
    }

    u64 morton{0};
    std::array<u32, DIM> ix{};                 // value-initialized: all zeros for any DIM
    u32 level{0};
    static constexpr std::size_t NCHILD = (1u << DIM);
    std::array<u32, NCHILD> child{};           // filled to npos32 in ctor

    bool is_leaf{true};
    u32 node2leafIdx{0};
    u32 parent{npos32};
  };

// ---- Field lives outside any class; stores an opaque basis_id (0..2) (DIM-independent) ----
  struct Field {
    uint8_t basis_id{0};           // 0..2
    std::vector<double> nodal;
    std::string name;
  };

// ---- Finite-element node with parent-space and physical coordinates (DIM-independent) ----
  template<std::size_t DIM>
  struct FEMNode {
    int         gid{ -1 };
    Point<DIM>  parent{};   // in [-1,1]^DIM
    Point<DIM>  physical{}; // mapped coordinates
  };

// ---- Per-basis registry of global nodes and element connectivity (DIM-independent) ----
  template<std::size_t DIM>
  struct BasisRegistry {
    std::unordered_map<u64, int> nodeMap;   // Morton key -> global id (gid)
    std::vector<FEMNode<DIM>> nodes;    // gid -> FEM node (parent + physical)
    std::vector<std::vector<int>> elem2glob; // per-element connectivity (local -> gid)

    // ---- Clears all registry data while keeping capacities (DIM-independent) ----
    inline void clear() noexcept {
      nodeMap.clear();
      nodes.clear();
      elem2glob.clear();
    }
  };

  template<std::size_t DIM> class Reinitializer;


// ---------- OctTree (minimum workable stub) ----------
  template<std::size_t DIM/*, class Derived*/>
  class OctTree {
    public:

      friend class Reinitializer<DIM>;

      using Basis = BasisT<DIM>;
      using PM = PtMake<DIM>;
      // ---- Constructs an empty tree with a single root covering [-1,1]^DIM; clamps minDepth; validates depth vs Morton capacity.
      OctTree(u32 maxDepth, u32 minDepth = 0)
        : _maxDepth(maxDepth),
          _minDepth(std::min(minDepth, maxDepth)) {
        static_assert(DIM == 2 || DIM == 3, "OctTree supports DIM = 2 or 3");
        recompute_quant_params();

        // 1) Build root TreeNode covering [-1,1]^DIM
        TreeNode<DIM> root{};

        // 2) Reset topology containers and seed with root
        _tree_nodes.clear();       _tree_nodes.push_back(root);
        _root = 0;
        _leaves.clear();           _leaves.push_back(_root);
        _leaf_ids.clear();


        // 3) Basis registries empty (fixed length = 3)
        for (int i = 0; i < 3; ++i) _basisReg[i].clear();

        // 4) Geometry not set yet
        _geom_ready = false;
      }

// ---- Virtual destructor (DIM-independent) ----
      virtual ~OctTree() = default;

// ---- Morton interleave wrapper (DIM-independent) ----
      inline uint64_t interleave(const u32Point<DIM>& x) const noexcept {
        static_assert(DIM == 2 || DIM == 3, "interleave supports DIM = 2 or 3");
#ifndef NDEBUG
        // 2D: up to 32 levels; 3D: up to 21 levels.
#endif
        return interleave_impl(x, std::integral_constant<std::size_t, DIM> {});
      }

// // ---- CRTP helpers to access the derived type (DIM-independent) ----
//       inline const Derived& derived() const noexcept {
//         return static_cast<const Derived&>(*this);
//       }
//       inline       Derived& derived()       noexcept {
//         return static_cast<      Derived&>(*this);
//       }

// ---- Resets the tree to a single root; optionally keeps geometry/fields; rebuilds registries (fixed _basisReg length = 3) (DIM-independent) ----
      inline void reset(bool keep_geometry = true, bool keep_fields = true) {
        // 1) optionally keep geometry
        if (!keep_geometry) {
          _X = {};
          _geom_ready = false;
        }

        // 2) optionally keep fields (values + active bases), otherwise clear them
        if (!keep_fields) {
          _fields.clear();
          _activeBases.clear();
          for (int i = 0; i < 3; ++i) _basisReg[i].clear();
        }
        else {
          // keep fields/active bases; registries will be rebuilt
          for (int i = 0; i < 3; ++i) _basisReg[i].clear();
        }

        // 3) reset topology to just the root (reuse capacity; no shrink)
        _leaves.clear();
        _tree_nodes.clear();
        TreeNode<DIM> root{};
        _tree_nodes.push_back(root);
        _root = 0;
        _leaves.push_back(_root);

        // 4) clear/refresh lookups (reuse capacity)
        _leaf_ids.clear();

        // 5) basis registries / field sizes depend on active bases
        post_topology_update(); // very cheap (there is only one element)
      }


      // ---- Defaulted special members (move-only; copy disabled) (DIM-independent) ----

      OctTree(const OctTree&)            = default;
      OctTree& operator=(const OctTree&) = default;
      OctTree(OctTree&&) noexcept        = default;
      OctTree& operator=(OctTree&&) noexcept = default;


      OctTree(const OctTree& oct1, u32 fid1,
              const OctTree& oct2, u32 fid2)
        : OctTree( /* choose a conservative cap & min */
            std::max(oct1.max_depth(), oct2.max_depth()),
            std::min(oct1.min_depth(), oct2.min_depth())) {

        static_assert(DIM == 2 || DIM == 3, "OctTree supports DIM = 2 or 3");
        recompute_quant_params();

        _X = oct1._X;
        _geom_ready = true;
        // Allow safe coarsening below min if either source allowed it
        this->set_allow_coarsen_below_min(oct1.allow_coarsen_below_min() || oct2.allow_coarsen_below_min());

        // Helper to pick the highest geometry basis (Q9/H27) for dense node sampling per level.
        const BasisT<DIM> geomHi =
          (DIM == 2) ? static_cast<BasisT<DIM>>(static_cast<int>(Basis2D::Q9))
          : static_cast<BasisT<DIM>>(static_cast<int>(Basis3D::H27));

        const u32 Lmax = std::max(oct1.max_depth(), oct2.max_depth());

        std::vector<Point<DIM>> s_all{};
        // Refine to match 'oct1'
        for (u32 L = 0; L <= Lmax; ++L) {
          // sample parent nodes of level L on oct1 using the highest geometry basis
          s_all.resize(0);
          oct1.extract_leaf_centers_in_level_range(L, L, s_all);
          if (s_all.empty()) continue;

          std::vector<Point<DIM>> pts_phys;
          pts_phys.reserve(s_all.size());
          for (const auto& s : s_all) pts_phys.push_back(oct1.parent_to_physical(s));

          // Refine this overlay so each point is contained in a leaf at level >= L
          this->refine_to_contain_points(pts_phys, L);
        }

        // Refine to match 'oct2'
        for (u32 L = 0; L <= Lmax; ++L) {
          s_all.resize(0);
          oct2.extract_leaf_centers_in_level_range(L, L, s_all);
          if (s_all.empty()) continue;

          std::vector<Point<DIM>> pts_phys;
          pts_phys.reserve(s_all.size());
          for (const auto& s : s_all) pts_phys.push_back(oct2.parent_to_physical(s));

          this->refine_to_contain_points(pts_phys, L);
        }

        // 2) Add two fields on the overlay using the SAME bases as the sources (keep original bases)
        const Field& src1 = oct1.field(fid1);
        const Field& src2 = oct2.field(fid2);
        const auto b1 = to_basis<DIM>(src1.basis_id);
        const auto b2 = to_basis<DIM>(src2.basis_id);

        const u32 dst_fid1 = this->add_field(b1, oct1.field(fid1).name + "1"); // overlay field for src1
        const u32 dst_fid2 = this->add_field(b2, oct1.field(fid2).name + "2"); // overlay field for src2

        // Build connectivity for both active bases and size nodal arrays
        this->rebuild_connectivity_active_bases();
        this->resize_fields_to_nodes();

        // 3) Project (sample) src fields into overlay at overlay nodes of each basis
        {
          Field& D = this->field(dst_fid1);
          const auto bb = to_basis<DIM>(D.basis_id);
          const auto& nodes = this->basis_nodes(bb);
          D.nodal.assign(nodes.size(), 0.0);
          for (size_t gid = 0; gid < nodes.size(); ++gid) {
            double val = 0.0;
            // evaluate *on oct1* at overlay's parent coordinate
            (void) oct1.evaluate_field_on_parent(fid1, nodes[gid].parent, val);
            D.nodal[gid] = val;
          }
        }
        {
          Field& D = this->field(dst_fid2);
          const auto bb = to_basis<DIM>(D.basis_id);
          const auto& nodes = this->basis_nodes(bb);
          D.nodal.assign(nodes.size(), 0.0);
          for (size_t gid = 0; gid < nodes.size(); ++gid) {
            double val = 0.0;
            // evaluate *on oct2* at overlay's parent coordinate
            (void) oct2.evaluate_field_on_parent(fid2, nodes[gid].parent, val);
            D.nodal[gid] = val;
          }
        }
      }

// --- swap Octrees (DIM-independent) ----
      friend void swap(OctTree& a, OctTree& b) noexcept {
        using std::swap;
        swap(a._maxDepth, b._maxDepth);
        a.recompute_quant_params();
        b.recompute_quant_params();
        swap(a._minDepth, b._minDepth);
        swap(a._allowCoarsenBelowMinDepth, b._allowCoarsenBelowMinDepth);

        swap(a._X, b._X);
        swap(a._geom_ready, b._geom_ready);

        swap(a._leaves, b._leaves);
        swap(a._leaf_ids, b._leaf_ids);
        swap(a._fields, b._fields);

        swap(a._basisReg[0], b._basisReg[0]);
        swap(a._basisReg[1], b._basisReg[1]);
        swap(a._basisReg[2], b._basisReg[2]);

        swap(a._tree_nodes, b._tree_nodes);
        swap(a._root, b._root);
        //swap(a._node2leafpos, b._node2leafpos);
        swap(a._activeBases, b._activeBases);
        swap(a._level_offset, b._level_offset);
      }

// ---- Sets the minimum refinement depth (clamped to _maxDepth; DIM-independent).
      inline void set_min_depth(u32 d) noexcept {
        _minDepth = (d <= _maxDepth) ? d : _maxDepth;
      }

// ---- Returns the current minimum refinement depth (DIM-independent) ----
      inline u32 min_depth() const noexcept {
        return _minDepth;
      }

// ---- Enables/disables coarsening below the minimum depth (DIM-independent) ----
      inline void set_allow_coarsen_below_min(bool v) noexcept {
        _allowCoarsenBelowMinDepth = v;
      }

// ---- Reports whether coarsening below the minimum depth is allowed (DIM-independent) ----
      inline bool allow_coarsen_below_min() const noexcept {
        return _allowCoarsenBelowMinDepth;
      }

// --- geometry setters (DIM-independent) ----
      void set_physical_coordinates(const std::array<std::array<double, NDOFS[DIM][2]>, DIM> X) {
        _X = X;
        _geom_ready = true;
      }

// ---- Ensures geometry has been initialized before any geometric operation (DIM-independent) ----
      inline void require_geometry() const {
        assert(_geom_ready && "Geometry not initialized: call the appropriate set_*_geometry(...) before using geometric ops.");
      }

// ---- Order _leaves by level, and then morton. (DIM-independnet) ---
      void sort_leaves_by_level_then_morton() {
        std::sort(_leaves.begin(), _leaves.end(),
        [this](u32 a, u32 b) {
          const auto& na = _tree_nodes[a];
          const auto& nb = _tree_nodes[b];
          if (na.level != nb.level) return na.level < nb.level;   // primary: level
          return na.morton < nb.morton;                           // secondary: morton
        });
      }

// ---- Order _leaves by morton. (DIM-independnet) ---
      void sort_leaves_by_morton() {
        std::sort(_leaves.begin(), _leaves.end(),
        [this](u32 a, u32 b) {
          return _tree_nodes[a].morton < _tree_nodes[b].morton;
        });
      }


// ---- Build offsets so that leaves of level L are in the half-open range [offsets[L], offsets[L+1]) of _leaves (DIM -independent) ----
      void compute_level_offsets() {
        _level_offset.clear();
        if (_leaves.empty()) return;

        // Find the maximum level present among current leaves
        u32 maxL = 0;
        for (u32 node_id : _leaves) {
          maxL = std::max(maxL, _tree_nodes[node_id].level);
        }

        // Count leaves per level
        _level_offset.assign(maxL + 2, 0u);               // +1 for sentinel
        for (u32 node_id : _leaves) {
          ++_level_offset[_tree_nodes[node_id].level];
        }

        // Prefix sum → offsets[L] = start index, offsets[L+1] = end
        u32 run = 0;
        for (u32 l = 0; l < _level_offset.size(); ++l) {
          u32 c = _level_offset[l];
          _level_offset[l] = run;
          run += c;
        }
        // Now: offsets[maxL+1] == _leaves.size() (the sentinel)
      }


// ---- Re-sorts leaves by level/Morton order, refreshes lookups/registries, and resizes field storage (DIM-independent) ----
      inline void post_topology_update() {
        sort_leaves_by_level_then_morton();
        compute_level_offsets();
        rebuild_leafpos_lookup();            // node index -> position in _leaves
        rebuild_connectivity_active_bases(); // per-basis registries/connectivity
        resize_fields_to_nodes();            // field storage size matches registries
      }

// ---- Rebuilds the node-index → leaf-position lookup table (DIM-independent) ----
      inline void rebuild_leafpos_lookup() {

        for (u32 i = 0; i < _tree_nodes.size(); i++) {
          _tree_nodes[i].node2leafIdx = npos32;
        }

        //_node2leafpos.assign(_tree_nodes.size(), npos32);
        for (u32 i = 0; i < static_cast<u32>(_leaves.size()); ++i) {
          const u32 node_idx = _leaves[i];
          assert(node_idx < _tree_nodes.size());
          _tree_nodes[node_idx].node2leafIdx = i;
        }
      }

// ---- Rebuilds connectivity/registries for all active bases if geometry is ready (DIM-independent) ----
      inline void rebuild_connectivity_active_bases() {
        if (!_geom_ready) return;
        for (Basis b : _activeBases) {
          rebuild_connectivity(b);
        }
      }

// ---- Rebuilds per-basis unique nodes and element-to-global connectivity (DIM-independent) ----
      inline void rebuild_connectivity(Basis b) {

        BasisRegistry<DIM>& R = _basisReg[static_cast<int>(b)];
        R.clear();

        // --- Reserve to avoid rehash/realloc churn ---
        const u32 nleaves = leaf_count();

        // Nodes per element by basis id (0..2) for each DIM
        const int per = leaf_dof_number(b);

        // A decent estimate for unique nodes is ~(per * leaves) * 0.6..0.8 due to sharing.
        // Keep it conservative (1.5x leaves*per) to minimize rehashes:
        const size_t est_nodes =
          std::max<size_t>(32, size_t(nleaves) * size_t(per) * size_t(3) / size_t(2)); // 1.5x

        R.nodeMap.reserve(est_nodes);
        R.nodes.reserve(est_nodes);
        R.elem2glob.reserve(nleaves);

        std::vector<Point<DIM>> s; // parent-space nodes per leaf

        for (u32 e = 0; e < nleaves; ++e) {
          // Fill parent-space element nodes for leaf e, for the requested basis
          extract_leaf_parent_coords(b, e, s);

          std::vector<int> conn;
          conn.reserve(s.size());

          // Insert/lookup each node; get_or_insert_gid is DIM-generic (quantizes by _maxDepth+1)
          for (const auto& p : s) {
            const int gid = get_or_insert_gid(R, p);
            conn.push_back(gid);
          }

          R.elem2glob.push_back(std::move(conn));
        }
      }

// ---- Resizes every field’s nodal vector to match its basis registry (DIM-independent); preserves existing values and zero-fills any new entries.
      inline void resize_fields_to_nodes() noexcept {
        for (auto& f : _fields) {
          const auto& R = _basisReg[static_cast<int>(f.basis_id)];
          const std::size_t need = R.nodes.size();
          if (f.nodal.size() != need) {
            f.nodal.resize(need, 0.0); // keeps existing entries; appends zeros if growing
          }
        }
      }

// ---- midpoint helper for any axis (DIM-independet)
      static inline double mid(double a, double b) {
        return 0.5 * (a + b);
      }

// ---- populate out_pts with leaf interpolation nodes for the given basis (DIM-independet)
      void extract_leaf_parent_coords(Basis basis, u32 leaf_idx,
                                      std::vector<Point<DIM>>& out_pts) const {
        const TreeNode<DIM>& leaf = _tree_nodes[_leaves[leaf_idx]];

        Point<DIM> X0, X1;
        leaf_bounds(leaf, X0, X1);

        out_pts.clear();

        if (DIM == 2) {
          const double x0 = X0[0], y0 = X0[1];
          const double x1 = X1[0], y1 = X1[1];
          const double xm = mid(x0, x1), ym = mid(y0, y1);

          // Explicit static_cast to the 2D enum so case labels match the type
          switch (static_cast<Basis2D>(basis)) {
          case Basis2D::Q4: {
            out_pts.reserve(4);
            out_pts = {
              Point<DIM>{x0, y0}, Point<DIM>{x1, y0},
              Point<DIM>{x1, y1}, Point<DIM>{x0, y1}
            };
          } break;

          case Basis2D::Q8: {
            out_pts.reserve(8);
            // corners 0..3
            out_pts = {
              Point<DIM>{x0, y0}, Point<DIM>{x1, y0},
              Point<DIM>{x1, y1}, Point<DIM>{x0, y1}
            };
            // edges 4..7: bottom, right, top, left
            out_pts.push_back(Point<DIM> {xm, y0});
            out_pts.push_back(Point<DIM> {x1, ym});
            out_pts.push_back(Point<DIM> {xm, y1});
            out_pts.push_back(Point<DIM> {x0, ym});
          } break;

          case Basis2D::Q9: {
            out_pts.reserve(9);
            // corners 0..3
            out_pts = {
              Point<DIM>{x0, y0}, Point<DIM>{x1, y0},
              Point<DIM>{x1, y1}, Point<DIM>{x0, y1}
            };
            // edges 4..7
            out_pts.push_back(Point<DIM> {xm, y0});
            out_pts.push_back(Point<DIM> {x1, ym});
            out_pts.push_back(Point<DIM> {xm, y1});
            out_pts.push_back(Point<DIM> {x0, ym});
            // center 8
            out_pts.push_back(Point<DIM> {xm, ym});
          } break;

          default:
            out_pts.clear(); // unsupported 2D basis
          }

        }
        else if (DIM == 3) {
          const double x0 = X0[0], y0 = X0[1], z0 = X0[2];
          const double x1 = X1[0], y1 = X1[1], z1 = X1[2];
          const double xm = mid(x0, x1), ym = mid(y0, y1), zm = mid(z0, z1);

          // Explicit static_cast to the 3D enum so case labels match the type
          switch (static_cast<Basis3D>(basis)) {
          case Basis3D::H8: {
            out_pts.reserve(8);
            out_pts = {
              PM::mk(x0, y0, z0), PM::mk(x1, y0, z0),
              PM::mk(x1, y1, z0), PM::mk(x0, y1, z0),
              PM::mk(x0, y0, z1), PM::mk(x1, y0, z1),
              PM::mk(x1, y1, z1), PM::mk(x0, y1, z1)
            };
          } break;

          case Basis3D::H20: {
            out_pts.reserve(20);
            // corners 0..7
            out_pts = {
              PM::mk(x0, y0, z0), PM::mk(x1, y0, z0),
              PM::mk(x1, y1, z0), PM::mk(x0, y1, z0),
              PM::mk(x0, y0, z1), PM::mk(x1, y0, z1),
              PM::mk(x1, y1, z1), PM::mk(x0, y1, z1)
            };
            // bottom z=-1 edges 8..11
            out_pts.push_back(PM::mk(xm, y0, z0));
            out_pts.push_back(PM::mk(x1, ym, z0));
            out_pts.push_back(PM::mk(xm, y1, z0));
            out_pts.push_back(PM::mk(x0, ym, z0));
            // top z=+1 edges 12..15
            out_pts.push_back(PM::mk(xm, y0, z1));
            out_pts.push_back(PM::mk(x1, ym, z1));
            out_pts.push_back(PM::mk(xm, y1, z1));
            out_pts.push_back(PM::mk(x0, ym, z1));
            // vertical edges 16..19
            out_pts.push_back(PM::mk(x0, y0, zm));
            out_pts.push_back(PM::mk(x1, y0, zm));
            out_pts.push_back(PM::mk(x1, y1, zm));
            out_pts.push_back(PM::mk(x0, y1, zm));
          } break;

          case Basis3D::H27: {
            out_pts.reserve(27);
            // corners 0..7
            out_pts = {
              PM::mk(x0, y0, z0), PM::mk(x1, y0, z0),
              PM::mk(x1, y1, z0), PM::mk(x0, y1, z0),
              PM::mk(x0, y0, z1), PM::mk(x1, y0, z1),
              PM::mk(x1, y1, z1), PM::mk(x0, y1, z1)
            };
            // bottom z=-1 edges 8..11
            out_pts.push_back(PM::mk(xm, y0, z0));
            out_pts.push_back(PM::mk(x1, ym, z0));
            out_pts.push_back(PM::mk(xm, y1, z0));
            out_pts.push_back(PM::mk(x0, ym, z0));
            // top z=+1 edges 12..15
            out_pts.push_back(PM::mk(xm, y0, z1));
            out_pts.push_back(PM::mk(x1, ym, z1));
            out_pts.push_back(PM::mk(xm, y1, z1));
            out_pts.push_back(PM::mk(x0, ym, z1));
            // vertical edges 16..19
            out_pts.push_back(PM::mk(x0, y0, zm));
            out_pts.push_back(PM::mk(x1, y0, zm));
            out_pts.push_back(PM::mk(x1, y1, zm));
            out_pts.push_back(PM::mk(x0, y1, zm));
            // face centers 20..25
            out_pts.push_back(PM::mk(xm, y0, zm)); // y=y0
            out_pts.push_back(PM::mk(x1, ym, zm)); // x=x1
            out_pts.push_back(PM::mk(xm, y1, zm)); // y=y1
            out_pts.push_back(PM::mk(x0, ym, zm)); // x=x0
            out_pts.push_back(PM::mk(xm, ym, z0)); // bottom face center
            out_pts.push_back(PM::mk(xm, ym, z1)); // top face center
            // cell center 26
            out_pts.push_back(PM::mk(xm, ym, zm));
          } break;

          default:
            out_pts.clear(); // unsupported 3D basis
          }

        }
        else {
          out_pts.clear(); // unsupported DIM at compile time
        }
      }

// ---- Parent-space center of a leaf (DIM-independent) ----
      void extract_leaf_parent_center_coord(u32 leaf_idx, Point<DIM>& out_pt) const {
        const TreeNode<DIM>& leaf = _tree_nodes[_leaves[leaf_idx]];

        Point<DIM> X0, X1;
        leaf_bounds(leaf, X0, X1);           // parent-space AABB of the leaf

        for (std::size_t d = 0; d < DIM; ++d)
          out_pt[d] = 0.5 * (X0[d] + X1[d]); // midpoint along each dimension
      }

      // ---- map parent → physical with DIM-appropriate shape (DIM-independet)
      inline Point<DIM> parent_to_physical(const Point<DIM>& s) const noexcept {
        require_geometry();

        double N[NDOFS[DIM][2]];

        if (DIM == 2) {
          // Shapes2 expects Point2
          Point2 s2{ s[0], s[1] };
          Shapes2D::Q9(s2, N);
        }
        else if (DIM == 3) {
          // Shapes3 expects Point3
          Point3 s3{ s[0], s[1], s[2] };
          Shapes3D::H27(s3, N);
        }
        else {
          static_assert(DIM == 2 || DIM == 3, "Unsupported DIM");
        }

        Point<DIM> x{}; // zero-initialize
        constexpr std::size_t Nnodes = NDOFS[DIM][2];
        for (std::size_t a = 0; a < Nnodes; ++a) {
          const double Na = N[a];
          for (std::size_t k = 0; k < DIM; ++k) x[k] += Na * _X[k][a]; // _X[axis][node]
        }
        return x;
      }


// ---- Returns the current number of leaves (DIM-independent) ----
      inline u32 leaf_count() const noexcept {
        return static_cast<u32>(_leaves.size());
      }

// // ---- Returns a cached vector of leaf positions [0..leaf_count-1] (DIM-independent) ----
//       inline const std::vector<u32>& leaf_indices() const {
//         const u32 n = static_cast<u32>(_leaves.size());
//         _leaf_ids.resize(n);
//         for (u32 i = 0; i < n; ++i) _leaf_ids[i] = i;
//         return _leaf_ids;
//       }

// ---- Adds a new field for basis `b`; sizes nodal storage to current topology (DIM-independent) ----
      inline u32 add_field(const Basis &b, const std::string &name) {
        // Insert returns {iterator, bool}. bool==true iff b was newly inserted.
        const bool first_time_for_basis = _activeBases.insert(b).second;
        if (first_time_for_basis) {
          // Build connectivity if not already there
          rebuild_connectivity(b);
        }
        // Now size the new field to whatever registry currently has
        const auto& R = _basisReg[static_cast<int>(b)];

        Field f;
        f.basis_id = to_basis_id<DIM>(b);
        f.nodal.assign(R.nodes.size(), 0.0);
        f.name = name;

        _fields.emplace_back(std::move(f));
        return static_cast<u32>(_fields.size() - 1);
      }

// ---- Returns a mutable reference to the field with id `fid` (bounds-checked in debug) (DIM-independent) ----
      inline Field& field(u32 fid) {
        assert(fid < _fields.size());
        return _fields[fid];
      }

// ---- Returns a const reference to the field with id `fid` (bounds-checked in debug) (DIM-independent) ----
      inline const Field& field(u32 fid) const {
        assert(fid < _fields.size());
        return _fields[fid];
      }

// ---- Newton inverse map using Q4/H8 warm-up then Q9/H27; generic J-solve (DIM-independent) ----
      inline bool inverse_map(const Point<DIM>& x,
                              Point<DIM>& s,
                              int maxIts = 30, double tol = 1e-12) const {
        require_geometry(); // ensure _X is set: _X[axis][node]

        // init at center
        s = {};
        const double tol2 = tol * tol;

        // Node counts per DIM
        const std::size_t Nnodes0 = NDOFS[DIM][0]; // Q4 (4) or H8 (9)
        const std::size_t Nnodes1 = NDOFS[DIM][2]; // Q9 (9) or H27 (27)

        // Buffers (z-components unused in 2D)
        // warm-up shapes
        std::vector<double> N0(Nnodes0), dN0dx(Nnodes0), dN0dy(Nnodes0), dN0dz(Nnodes0, 0.0);
        // full shapes
        std::vector<double> N1(Nnodes1), dN1dx(Nnodes1), dN1dy(Nnodes1), dN1dz(Nnodes1, 0.0);

        // --- Warm-up with Q4 (2D) or H8 (3D) at current s ---
        if (DIM == 2) {
          Point2 s2{ s[0], s[1] };
          Shapes2D::Q4(s2, N0.data());
        }
        else if (DIM == 3) {
          Point3 s3{ s[0], s[1], s[2] };
          Shapes3D::H8(s3, N0.data());
        }
        else {
          return false; // unsupported DIM
        }

        // Xp = sum_a N0[a] * X[:,a]
        Point<DIM> Xp{}; // zero-init
        for (std::size_t a = 0; a < Nnodes0; ++a) {
          const double Na = N0[a];
          for (std::size_t k = 0; k < DIM; ++k) Xp[k] += Na * _X[k][a];
        }

        // J at warm-up shape
        if (DIM == 2) {
          Point2 s2{ s[0], s[1] };
          Shapes2D::Q4_dN(s2, dN0dx.data(), dN0dy.data());
        }
        else { /* DIM==3 */
          Point3 s3{ s[0], s[1], s[2] };
          Shapes3D::H8_dN(s3, dN0dx.data(), dN0dy.data(), dN0dz.data());
        }

        double J[DIM][DIM] = {{0.0}};
        for (std::size_t a = 0; a < Nnodes0; ++a) {
          const double dNa[3] = { dN0dx[a], dN0dy[a], dN0dz[a] }; // only first DIM used
          for (std::size_t i = 0; i < DIM; ++i) {
            const double Xi_a = _X[i][a];
            for (std::size_t j = 0; j < DIM; ++j) J[i][j] += dNa[j] * Xi_a;
          }
        }

        // Residual r = Xp - x
        double r[DIM];
        for (std::size_t k = 0; k < DIM; ++k) r[k] = Xp[k] - x[k];

        for (int it = 0; it < maxIts; ++it) {
          // ---- Solve J * ds = r (supports DIM=2 and DIM=3)
          double ds[DIM];

          if (DIM == 2) {
            const double A11 = J[0][0], A12 = J[0][1];
            const double A21 = J[1][0], A22 = J[1][1];
            const double det = A11 * A22 - A12 * A21;
            if (std::fabs(det) < 1e-20) break;
            const double inv = 1.0 / det;
            ds[0] = (A22 * r[0] - A12 * r[1]) * inv;
            ds[1] = (-A21 * r[0] + A11 * r[1]) * inv;
          }
          else { /* DIM==3 */
            const double A11 = J[0][0], A12 = J[0][1], A13 = J[0][2];
            const double A21 = J[1][0], A22 = J[1][1], A23 = J[1][2];
            const double A31 = J[2][0], A32 = J[2][1], A33 = J[2][2];
            const double det =
              A11 * (A22 * A33 - A23 * A32) -
              A12 * (A21 * A33 - A23 * A31) +
              A13 * (A21 * A32 - A22 * A31);
            if (std::fabs(det) < 1e-20) break;
            const double inv = 1.0 / det;
            ds[0] = ((A22 * A33 - A23 * A32) * r[0]
                     - (A12 * A33 - A13 * A32) * r[1]
                     + (A12 * A23 - A13 * A22) * r[2]) * inv;
            ds[1] = (-(A21 * A33 - A23 * A31) * r[0]
                     + (A11 * A33 - A13 * A31) * r[1]
                     - (A11 * A23 - A13 * A21) * r[2]) * inv;
            ds[2] = ((A21 * A32 - A22 * A31) * r[0]
                     - (A11 * A32 - A12 * A31) * r[1]
                     + (A11 * A22 - A12 * A21) * r[2]) * inv;
          }

          // Update s := s - ds with clamp to [-1,1]
          for (std::size_t k = 0; k < DIM; ++k)
            s[k] = std::max(-1.0, std::min(1.0, s[k] - ds[k]));

          // Recompute residual with Q9/H27 at updated s
          if (DIM == 2) {
            Point2 s2{ s[0], s[1] };
            Shapes2D::Q9(s2, N1.data());
          }
          else { /* DIM==3 */
            Point3 s3{ s[0], s[1], s[2] };
            Shapes3D::H27(s3, N1.data());
          }

          for (std::size_t k = 0; k < DIM; ++k) Xp[k] = 0.0;
          for (std::size_t a = 0; a < Nnodes1; ++a) {
            const double Na = N1[a];
            for (std::size_t k = 0; k < DIM; ++k) Xp[k] += Na * _X[k][a];
          }

          double nrm2 = 0.0;
          for (std::size_t k = 0; k < DIM; ++k) {
            const double rk = (Xp[k] - x[k]);
            r[k] = rk;
            nrm2 += rk * rk;
          }
          if (nrm2 < tol2) return true;

          // Next Jacobian with Q9/H27 derivatives
          if (DIM == 2) {
            Point2 s2{ s[0], s[1] };
            Shapes2D::Q9_dN(s2, dN1dx.data(), dN1dy.data());
            for (std::size_t i = 0; i < 2; ++i)
              for (std::size_t j = 0; j < 2; ++j) J[i][j] = 0.0;
            for (std::size_t a = 0; a < Nnodes1; ++a) {
              const double dNa0 = dN1dx[a], dNa1 = dN1dy[a];
              for (std::size_t i = 0; i < 2; ++i) {
                const double Xi_a = _X[i][a];
                J[i][0] += dNa0 * Xi_a;
                J[i][1] += dNa1 * Xi_a;
              }
            }
          }
          else { /* DIM==3 */
            Point3 s3{ s[0], s[1], s[2] };
            Shapes3D::H27_dN(s3, dN1dx.data(), dN1dy.data(), dN1dz.data());
            for (std::size_t i = 0; i < 3; ++i)
              for (std::size_t j = 0; j < 3; ++j) J[i][j] = 0.0;
            for (std::size_t a = 0; a < Nnodes1; ++a) {
              const double dNa0 = dN1dx[a], dNa1 = dN1dy[a], dNa2 = dN1dz[a];
              for (std::size_t i = 0; i < 3; ++i) {
                const double Xi_a = _X[i][a];
                J[i][0] += dNa0 * Xi_a;
                J[i][1] += dNa1 * Xi_a;
                J[i][2] += dNa2 * Xi_a;
              }
            }
          }
        }

        return false; // not converged
      }

// ---- Locate the leaf containing parent-space point s in [-1,1]^DIM (DIM-independent) ----
      inline u32 locate_leaf_on_parent(const Point<DIM>& s) const noexcept {
        // Bounds check in [-1,1]^DIM
        for (std::size_t k = 0; k < DIM; ++k)
          if (s[k] < -1.0 || s[k] > 1.0) return npos32;

        // Degenerate tree: no refinement levels
        if (_maxDepth == 0) return _root;

        // Map s to integer grid indices on a 2^_maxDepth grid per axis.
        // scale = 2^{_maxDepth-1}; clamp to max index so +1.0 maps to the last cell.
        const double scale = std::ldexp(1.0, int(_maxDepth) - 1);
        const u32 max_idx  = (u32(1) << _maxDepth) - 1;

        u32Point<DIM> gix; // integer coordinates per axis
        for (std::size_t k = 0; k < DIM; ++k) {
          const double t  = (s[k] + 1.0) * scale;      // [-1,1] -> [0, 2^_maxDepth]
          const u32     i = (t <= 0.0) ? 0u : (u32)t;  // fast floor for non-negative t
          gix[k] = (i > max_idx) ? max_idx : i;        // clamp
        }

        // Descend from root using bit k of gix[k] at each level.
        u32 node = _root;
        for (u32 level = 0; level < _maxDepth; ++level) {
          const TreeNode<DIM>& n = _tree_nodes[node];
          if (n.is_leaf) return node;

          const int shift = int(_maxDepth - level - 1);
          // Assemble child code: bit k = ((gix[k] >> shift) & 1)
          u32 code = 0;
          for (std::size_t k = 0; k < DIM; ++k)
            code |= (((gix[k] >> shift) & 1u) << k);

          const u32 c = n.child[code];
          if (c == npos32) return node; // stop at last existing ancestor
          node = c;
        }
        return node;
      }

// ---- physical nodes for a leaf (DIM-independent) ----
      void leaf_physical_nodes(Basis basis,
                               u32 leaf_idx,
                               std::vector<Point<DIM>>& out_xyz) const {
        require_geometry();

        // Get parent-space nodes
        std::vector<Point<DIM>> s;
        s.reserve(32); // small heuristic; avoids a couple reallocations
        extract_leaf_parent_coords(basis, leaf_idx, s);

        // Map to physical space
        out_xyz.resize(s.size());
        std::transform(s.begin(), s.end(), out_xyz.begin(),
        [this](const Point<DIM>& p) {
          return parent_to_physical(p);
        });
      }

// ---- Runs one adapt cycle (coarsen then refine passes) until stable or max_passes (DIM-independent) ----
      template<class CoarsenPred, class RefinePred>
      inline std::size_t adapt_cycle(CoarsenPred&& should_coarsen,
                                     RefinePred&&  should_refine,
                                     u32   max_passes  = 10) {
        ensure_min_depth(); // enforce refinement floor only
        std::size_t total = 0;

        for (u32 pass = 0; pass < max_passes; ++pass) {
          std::size_t changed = 0;
          changed += coarsen_pass(should_coarsen);
          changed += refine_pass(should_refine);
          if (changed == 0) break;
          total += changed;
        }

        enforce_balance();
        post_topology_update();
        return total;
      }

// ---- Perform one refinement pass: split leaves if predicate says so (DIM-independent) ----
      template<class RefinePred>
      std::size_t refine_pass(RefinePred&& should_refine) {
        std::vector<u32> newLeaves;
        newLeaves.reserve(_leaves.size() * 2u); // heuristic (not all leaves refine)
        std::size_t refined = 0;

        const Basis basisHi = (DIM == 2)
                              ? static_cast<Basis>(Basis2D::Q9)
                              : static_cast<Basis>(Basis3D::H27);

        const u32 n0 = static_cast<u32>(leaf_count()); // snapshot leaf count

        for (u32 idx = 0; idx < n0; ++idx) {
          const u32 leaf_node_idx = _leaves[idx];                 // index into _tree_nodes
          const TreeNode<DIM> leaf_copy = _tree_nodes[leaf_node_idx]; // snapshot (safe vs reallocation)

          // Probe geometry for refinement (keeps current API semantics: probe by leaf position 'idx')
          std::vector<Point<DIM>> pts_s, pts_xyz;
          extract_leaf_parent_coords(basisHi, idx, pts_s);
          leaf_physical_nodes(basisHi, idx, pts_xyz);
          std::vector<std::array<double, NDOFS[DIM][2]>> dummy; // placeholder if predicate needs extra data

          if (leaf_copy.is_leaf &&
              leaf_copy.level < _maxDepth &&
              should_refine(idx, pts_s, pts_xyz, dummy)) {
            // Turn parent internal and clear its child slots
            _tree_nodes[leaf_node_idx].is_leaf = false;
            for (std::size_t c = 0; c < TreeNode<DIM>::NCHILD; ++c)
              _tree_nodes[leaf_node_idx].child[c] = npos32;

            // Precompute base indices shifted by 1 (child = 2*parent + {0,1} per axis)
            u32Point<DIM> base2 = leaf_copy.ix;
            for (std::size_t k = 0; k < DIM; ++k) base2[k] <<= 1;

            // Create all 2^DIM children; 'mask' is the child code with bit-k selecting hi/lo along axis k
            const u32 child_count = static_cast<u32>(TreeNode<DIM>::NCHILD);
            for (u32 mask = 0; mask < child_count; ++mask) {
              u32Point<DIM> cix;
              for (std::size_t k = 0; k < DIM; ++k)
                cix[k] = base2[k] + ((mask >> k) & 1u);

              TreeNode<DIM> child{};
              child.ix      = cix;
              child.level   = leaf_copy.level + 1;
              child.is_leaf = true;
              child.parent  = leaf_node_idx;
              for (std::size_t c = 0; c < TreeNode<DIM>::NCHILD; ++c)
                child.child[c] = npos32;

              child.morton = interleave(cix); // or derived().interleave(cix) if CRTP

              const u32 child_idx = static_cast<u32>(_tree_nodes.size());
              _tree_nodes.push_back(child);

              _tree_nodes[leaf_node_idx].child[mask] = child_idx; // assign into parent
              newLeaves.push_back(child_idx);
            }

            ++refined;
          }
          else {
            // Keep existing leaf
            newLeaves.push_back(leaf_node_idx);
          }
        }

        _leaves.swap(newLeaves);
        return refined;
      }

// ---- coarsen pass with DIM-appropriate grouping and probe points (DIM-independet) ----
      template<class CoarsenPred>
      std::size_t coarsen_pass(CoarsenPred&& should_coarsen) {
        struct Group {
          u32      pos[(1u << DIM)]; // positions in _leaves of the children (4 in 2D, 8 in 3D)
          uint8_t  cnt   = 0;
          uint32_t tag   = 0;        // generation tag
        };

        constexpr u32 NCH = (1u << DIM);         // number of children per parent
        const u32 L = (u32)_leaves.size();
        const u32 N = (u32)_tree_nodes.size();

        // --- thread-local scratch to avoid alloc churn ---
        static thread_local std::vector<u32> touched;
        static thread_local std::vector<u32> to_add;
        static thread_local std::vector<u32> newLeaves;
        static thread_local std::vector<Group> per;
        static thread_local uint32_t per_epoch = 1;

        if (per.size() < N) per.resize(N);        // grow once, never shrink
        ++per_epoch;                              // bump generation

        touched.clear(); touched.reserve(L / NCH + NCH);
        to_add.clear();  to_add.reserve(L / NCH + NCH);
        newLeaves.clear(); newLeaves.reserve(L);

        // removal marks (epoch-tagged, no clears)
        static thread_local std::vector<uint32_t> mark_tag;
        static thread_local uint32_t mark_epoch = 1;
        if (mark_tag.size() < L) mark_tag.resize(L, 0);
        ++mark_epoch;

        auto mark_pos = [&](u32 i) {
          mark_tag[i] = mark_epoch;
        };
        auto is_marked = [&](u32 i) -> bool {
          return (i < mark_tag.size()) && (mark_tag[i] == mark_epoch);
        };

        // tiny fixed buffers for parent-space samples
        std::array<Point<DIM>, NDOFS[DIM][2]> s_buf, xyz_buf;

        // Reuse vector wrappers for the predicate (avoid allocs)
        static thread_local std::vector<Point<DIM>> pts_s_vec, pts_xyz_vec;

        auto fill_probe_nodes = [&](const Point<DIM>& X0, const Point<DIM>& X1,
        int nprobe, std::array<Point<DIM>, NDOFS[DIM][2]>& sbuf) -> int {
          if (DIM == 2) {
            const double x0 = X0[0], y0 = X0[1];
            const double x1 = X1[0], y1 = X1[1];
            const double xm = 0.5 * (x0 + x1);
            const double ym = 0.5 * (y0 + y1);
            int n = 0;
            // corners 0..3
            sbuf[n++] = Point<DIM> {x0, y0};
            sbuf[n++] = Point<DIM> {x1, y0};
            sbuf[n++] = Point<DIM> {x1, y1};
            sbuf[n++] = Point<DIM> {x0, y1};
            if (nprobe >= 8) {
              // edges 4..7: bottom, right, top, left
              sbuf[n++] = Point<DIM> {xm, y0};
              sbuf[n++] = Point<DIM> {x1, ym};
              sbuf[n++] = Point<DIM> {xm, y1};
              sbuf[n++] = Point<DIM> {x0, ym};
            }
            if (nprobe == 9) {
              // center 8
              sbuf[n++] = Point<DIM> {xm, ym};
            }
            return n;
          }
          else {
            const double x0 = X0[0], y0 = X0[1], z0 = X0[2];
            const double x1 = X1[0], y1 = X1[1], z1 = X1[2];
            const double xm = 0.5 * (x0 + x1);
            const double ym = 0.5 * (y0 + y1);
            const double zm = 0.5 * (z0 + z1);
            int n = 0;
            // corners 0..7
            sbuf[n++] = PM::mk(x0, y0, z0); sbuf[n++] = PM::mk(x1, y0, z0);
            sbuf[n++] = PM::mk(x1, y1, z0); sbuf[n++] = PM::mk(x0, y1, z0);
            sbuf[n++] = PM::mk(x0, y0, z1); sbuf[n++] = PM::mk(x1, y0, z1);
            sbuf[n++] = PM::mk(x1, y1, z1); sbuf[n++] = PM::mk(x0, y1, z1);
            if (nprobe >= 20) {
              // z=-1 edges (8..11)
              sbuf[n++] = PM::mk(xm, y0, z0);
              sbuf[n++] = PM::mk(x1, ym, z0);
              sbuf[n++] = PM::mk(xm, y1, z0);
              sbuf[n++] = PM::mk(x0, ym, z0);
              // z=+1 edges (12..15)
              sbuf[n++] = PM::mk(xm, y0, z1);
              sbuf[n++] = PM::mk(x1, ym, z1);
              sbuf[n++] = PM::mk(xm, y1, z1);
              sbuf[n++] = PM::mk(x0, ym, z1);
              // vertical edges (16..19)
              sbuf[n++] = PM::mk(x0, y0, zm);
              sbuf[n++] = PM::mk(x1, y0, zm);
              sbuf[n++] = PM::mk(x1, y1, zm);
              sbuf[n++] = PM::mk(x0, y1, zm);
            }
            if (nprobe == 27) {
              // face centers (20..25)
              sbuf[n++] = PM::mk(xm, y0, zm);
              sbuf[n++] = PM::mk(x1, ym, zm);
              sbuf[n++] = PM::mk(xm, y1, zm);
              sbuf[n++] = PM::mk(x0, ym, zm);
              sbuf[n++] = PM::mk(xm, ym, z0);
              sbuf[n++] = PM::mk(xm, ym, z1);
              // cell center (26)
              sbuf[n++] = PM::mk(xm, ym, zm);
            }
            return n;
          }
        };

        // 1) group leaves by parent (store LEAF POSITIONS, tag parent lazily)
        for (u32 pos = 0; pos < L; ++pos) {
          const u32 c = _leaves[pos];
          const u32 p = _tree_nodes[c].parent;
          if (p == npos32) continue;

          Group& g = per[p];
          if (g.tag != per_epoch) {        // first time we see this parent in this pass
            g.tag = per_epoch;
            g.cnt = 0;
            touched.push_back(p);
          }
          if (g.cnt < NCH) g.pos[g.cnt++] = pos;
        }

        // 2) coarsen candidates
        std::size_t coarsened = 0;

        for (u32 pidx : touched) {
          Group& g = per[pidx];
          if (g.cnt != NCH) continue;

          // level/minDepth check
          const u32 kid0 = _leaves[g.pos[0]];
          const u32 lev  = _tree_nodes[kid0].level;
          if (lev <= _minDepth && !_allowCoarsenBelowMinDepth) continue;

          // probe parent cell bounds
          TreeNode<DIM>& P = _tree_nodes[pidx];
          Point<DIM> X0, X1;
          leaf_bounds(P, X0, X1);

          // precompute midpoints once
          Point<DIM> M;
          for (u32 d = 0; d < DIM; ++d) M[d] = 0.5 * (X0[d] + X1[d]);

          auto child_bounds_from_parent = [&](u32 child_ord, Point<DIM>& C0, Point<DIM>& C1) {
            if constexpr(DIM == 2) {
              const bool rx = child_ord & 1u; // bit 0 = x high
              const bool ry = child_ord & 2u; // bit 1 = y high
              C0 = { rx ? M[0] : X0[0], ry ? M[1] : X0[1] };
              C1 = { rx ? X1[0] : M[0], ry ? X1[1] : M[1] };
            }
            else {
              const bool rx = child_ord & 1u;        // x high?
              const bool ry = (child_ord >> 1) & 1u; // y high?
              const bool rz = (child_ord >> 2) & 1u; // z high?
              C0 = { rx ? M[0] : X0[0], ry ? M[1] : X0[1], rz ? M[2] : X0[2] };
              C1 = { rx ? X1[0] : M[0], ry ? X1[1] : M[1], rz ? X1[2] : M[2] };
            }
          };

          // ---- CHILD-ONLY CHECK: all must pass ----
          bool all_children_ok = true;

          for (u32 k = 0; k < NCH; ++k) {
            const u32 leaf_pos = g.pos[k];
            const u32 child_id = _leaves[leaf_pos];
            const TreeNode<DIM>& C = _tree_nodes[child_id];

            // // min-depth guard per child (explicit)
            // if (C.level <= _minDepth && !_allowCoarsenBelowMinDepth) {
            //   all_children_ok = false;
            //   break;
            // }

            // derive this child's bbox from the parent bbox
            // (Assumes canonical child ordering; if not, compute ordinal from C.morton ^ P.morton.)
            Point<DIM> C0, C1;
            child_bounds_from_parent(/*child ordinal*/ k, C0, C1);

            // sample inside the child (PHYSICAL sampling to match your current fill_probe_nodes)
            const int nC = fill_probe_nodes(C0, C1, NDOFS[DIM][2], s_buf);

            // map to physical
            for (int i = 0; i < nC; ++i) xyz_buf[i] = parent_to_physical(s_buf[i]);

            pts_s_vec.resize(nC);
            pts_xyz_vec.resize(nC);
            for (int i = 0; i < nC; ++i) {
              pts_s_vec[i] = s_buf[i];   // parent points
              pts_xyz_vec[i] = xyz_buf[i];   // physical-points)
            }

            // reuse the same predicate for children; pass child's morton and level
            if (!should_coarsen(C.morton, C.level, pts_s_vec, pts_xyz_vec, /*dummy*/{})) {
              all_children_ok = false;
              break; // short-circuit: one fails => all fail
            }
          }

          if (!all_children_ok) continue; // do not coarsen this parent

          // ---- reach here: ALL children passed -> mark and promote parent ----
          for (u32 k = 0; k < NCH; ++k) {
            const u32 pos = g.pos[k];
            _tree_nodes[_leaves[pos]].is_leaf = false;
            mark_pos(pos);
          }
          P.is_leaf = true;
          for (u32 c = 0; c < NCH; ++c) P.child[c] = npos32;
          to_add.push_back(pidx);
          ++coarsened;
        }

        // 3) rebuild leaf set (single pass)
        if (coarsened) {
          for (u32 i = 0; i < L; ++i) if (!is_marked(i)) newLeaves.push_back(_leaves[i]);
          newLeaves.insert(newLeaves.end(), to_add.begin(), to_add.end());
          _leaves.swap(newLeaves);
        }

        return coarsened;
      }


// ---- Refines the given leaf into 2^DIM children and updates topology; returns true if a split occurred (DIM-independent) ----
      inline bool refine_leaf_once(u32 leaf_node_idx) noexcept {
        if (leaf_node_idx == npos32) return false;

        const TreeNode<DIM> parent_copy = _tree_nodes[leaf_node_idx];
        if (!parent_copy.is_leaf || parent_copy.level >= _maxDepth) return false;

        // Remove old leaf from list (swap-erase)
        {
          auto it = std::find(_leaves.begin(), _leaves.end(), leaf_node_idx);
          if (it != _leaves.end()) {
            *it = _leaves.back();
            _leaves.pop_back();
          }
        }

        // Turn parent into an internal node and clear its child slots
        _tree_nodes[leaf_node_idx].is_leaf = false;
        for (std::size_t c = 0; c < TreeNode<DIM>::NCHILD; ++c)
          _tree_nodes[leaf_node_idx].child[c] = npos32;

        // Precompute base indices shifted by 1
        u32Point<DIM> base2 = parent_copy.ix;
        for (std::size_t k = 0; k < DIM; ++k) base2[k] <<= 1;

        // Create all 2^DIM children
        const u32 child_count = static_cast<u32>(TreeNode<DIM>::NCHILD);
        for (u32 mask = 0; mask < child_count; ++mask) {
          // Child logical indices: add 0/1 per axis according to mask bit k
          u32Point<DIM> cix;
          for (std::size_t k = 0; k < DIM; ++k)
            cix[k] = base2[k] + ((mask >> k) & 1u);

          TreeNode<DIM> child{};
          child.ix      = cix;
          child.level   = parent_copy.level + 1;
          child.is_leaf = true;
          child.parent  = leaf_node_idx;
          // init child slots
          for (std::size_t c = 0; c < TreeNode<DIM>::NCHILD; ++c)
            child.child[c] = npos32;

          // Morton key (use your DIM-generic interleave)
          child.morton = interleave(cix);

          const u32 cidx = static_cast<u32>(_tree_nodes.size());
          _tree_nodes.push_back(child);
          _tree_nodes[leaf_node_idx].child[mask] = cidx;
          _leaves.push_back(cidx);
        }

        return true;
      }

// ---- neighbor by face for quad/octant trees (DIM-independent) ----
      u32 neighbor_leaf_by_face_any(u32 leaf_node_idx, int dir) const {
        // Faces are encoded as: 0..(2*DIM-1), where axis = dir>>1 and side = dir&1 (0: -, 1: +).
        if (dir < 0 || dir >= int(2 * DIM)) return npos32;

        const auto& me   = _tree_nodes[leaf_node_idx];
        const u32   axis = u32(dir >> 1);
        const bool  plus = (dir & 1) != 0;

        // Child code within parent (Morton-style): bit a = parity of ix[a].
        auto child_code = [&](u32 node_idx) -> u32 {
          const auto& nd = _tree_nodes[node_idx];
          u32 code = 0;
          for (u32 a = 0; a < DIM; ++a) code |= ((nd.ix[a] & 1u) << a);
          return code;
        };

        // Try stepping to sibling across the requested face at the same parent.
        auto try_local = [&](u32 code, u32 & sib_code) -> bool {
          const u32 bit = (1u << axis);
          if (plus) {
            if ((code & bit) == 0u) {
              sib_code = code | bit;
              return true;
            }
          }
          else {
            if ((code & bit) != 0u) {
              sib_code = code & ~bit;
              return true;
            }
          }
          return false;
        };

        // Descend from a subtree root to the leaf that touches the given face midpoint of `me`.
        auto descend_to_face_mid = [this, &me, dir](u32 subtree_root) -> u32 {
          assert(dir >= 0 && dir < int(2 * DIM) && "dir must be in [0, 2*DIM)");
          Point<DIM> X0, X1; leaf_bounds(me, X0, X1);

          // Build face midpoint q on `me` for the requested face.
          const u32  ax = u32(dir >> 1);
          const bool hi = (dir & 1) != 0;
          Point<DIM> q;
          for (u32 k = 0; k < DIM; ++k) q[k] = 0.5 * (X0[k] + X1[k]);
          q[ax] = hi ? X1[ax] : X0[ax];

          // Walk down choosing the child that contains q; stop at leaf or last existing ancestor.
          u32 node = subtree_root;
          Point<DIM> N0, N1;
          for (;;) {
            const auto& n = _tree_nodes[node];
            if (n.is_leaf) return node;

            leaf_bounds(n, N0, N1);
            u32 code = 0;
            for (u32 k = 0; k < DIM; ++k) {
              const double ck = 0.5 * (N0[k] + N1[k]);
              code |= ((q[k] >= ck) ? 1u : 0u) << k;
            }

            const u32 ch = n.child[code];
            if (ch == npos32) return node; // stop at last existing ancestor
            node = ch;
          }
        };

        // 1) Same-parent sibling?
        if (me.parent != npos32) {
          const u32 code = child_code(leaf_node_idx);
          u32 sib_code = code;
          if (try_local(code, sib_code)) {
            const u32 cand = _tree_nodes[me.parent].child[sib_code];
            if (cand != npos32) return descend_to_face_mid(cand);
          }
        }

        // 2) Climb until a parent where a local step across the requested face is possible; then descend.
        u32 cur = leaf_node_idx;
        for (;;) {
          const auto& c = _tree_nodes[cur];
          if (c.parent == npos32) return npos32;

          const u32 pidx = c.parent;
          const u32 code = child_code(cur);
          u32 sib_code = code;

          if (try_local(code, sib_code)) {
            const u32 sib = _tree_nodes[pidx].child[sib_code];
            if (sib == npos32) return npos32;
            return descend_to_face_mid(sib);
          }

          cur = pidx; // keep climbing
        }
      }

// ---- Returns the vector of current leaf node indices (DIM-independent) ----
      inline const std::vector<u32>& leaves() const noexcept {
        return _leaves;
      }

// ---- Returns a const reference to the tree node at index idx (DIM-independent) ----
      inline const TreeNode<DIM>& node(u32 idx) const noexcept {
        return _tree_nodes[idx];
      }

// ---- evaluate scalar field at parent coords for quad/octant bases (DIM-independent) ----
      bool evaluate_field_on_parent(u32 fid, const Point<DIM>& s, double& value) const {
        if (fid >= _fields.size()) return false;

        u32        leaf_node_idx = npos32;
        Point<DIM> shat; // local reference coords in [-1,1]^DIM

        // 1) locate leaf and local reference coords in [-1,1]^DIM
        if (!locate_leaf_on_parent_and_ref(s, leaf_node_idx, shat)) {
          value = -1.;// if the node is outside the bounding box [-1,+1]^DIM
          return false;
        }


        // 2) leaf node index -> position in coefficient storage
        const u32 leaf_pos = (leaf_node_idx < _tree_nodes.size()) ? _tree_nodes[leaf_node_idx].node2leafIdx : npos32;
        if (leaf_pos == npos32) return false;

        // 3) access field + connectivity
        const Field&              f = _fields[fid];
        const BasisRegistry<DIM>& R = _basisReg[(int)f.basis_id];
        const auto&               conn = R.elem2glob[leaf_pos];

        value = 0.0;

        // Runtime DIM branch (C++14): construct both s2/s3 safely, use enum casts per DIM
        if (DIM == 2) {
          Point2 s2{ shat[0], shat[1] };

          switch (static_cast<Basis2D>(to_basis<DIM>(f.basis_id))) {
          case Basis2D::Q4: {
            double N[4];  Shapes2D::Q4(s2, N);
            for (int a = 0; a < 4; ++a) value += N[a] * f.nodal[conn[a]];
          } break;
          case Basis2D::Q8: {
            double N[8];  Shapes2D::Q8(s2, N);
            for (int a = 0; a < 8; ++a) value += N[a] * f.nodal[conn[a]];
          } break;
          case Basis2D::Q9: {
            double N[9];  Shapes2D::Q9(s2, N);
            for (int a = 0; a < 9; ++a) value += N[a] * f.nodal[conn[a]];
          } break;
          default: return false;
          }

        }
        else if (DIM == 3) {
          // Only access shat[2] in the DIM==3 path
          Point3 s3{ shat[0], shat[1], shat[2] };

          switch (static_cast<Basis3D>(to_basis<DIM>(f.basis_id))) {
          case Basis3D::H8: {
            double N[8];   Shapes3D::H8(s3, N);
            for (int a = 0; a < 8;  ++a) value += N[a] * f.nodal[conn[a]];
          } break;
          case Basis3D::H20: {
            double N[20];  Shapes3D::H20(s3, N);
            for (int a = 0; a < 20; ++a) value += N[a] * f.nodal[conn[a]];
          } break;
          case Basis3D::H27: {
            double N[27];  Shapes3D::H27(s3, N);
            for (int a = 0; a < 27; ++a) value += N[a] * f.nodal[conn[a]];
          } break;
          default: return false;
          }

        }
        else {
          return false; // unsupported DIM
        }

        return true;
      }

// ---- evaluate gradient at parent coords for quad/octant bases (DIM-independent) ----
      bool evaluate_gradient_on_parent(u32 fid, const Point<DIM>& s,
                                       std::array<double, DIM>& gradient) const {
        if (fid >= _fields.size()) return false;

        u32        leaf_node_idx = npos32;
        Point<DIM> shat; // local reference coords in [-1,1]^DIM

        // 1) locate leaf and local reference coords in [-1,1]^DIM
        if (!locate_leaf_on_parent_and_ref(s, leaf_node_idx, shat)) return false;

        // 2) leaf node index -> position in coefficient storage

        const u32 leaf_pos = (leaf_node_idx < _tree_nodes.size()) ? _tree_nodes[leaf_node_idx].node2leafIdx : npos32;
        if (leaf_pos == npos32) return false;

        // 3) access field + connectivity
        const Field&              f = _fields[fid];
        const BasisRegistry<DIM>& R = _basisReg[(int)f.basis_id];
        const auto&               conn = R.elem2glob[leaf_pos];

        for (int idim = 0; idim < DIM; idim++)
          gradient[idim] = 0;

        if (DIM == 2) {
          Point2 s2{ shat[0], shat[1] };

          switch (static_cast<Basis2D>(to_basis<DIM>(f.basis_id))) {
          case Basis2D::Q4: {
            std::array<std::array<double, 4>, 2> dN{};
            Shapes2D::Q4_dN(s2, dN[0].data(), dN[1].data());
            for (int a = 0; a < 4; ++a)
              for (int idim = 0; idim < 2; idim++)
                gradient[idim] += dN[idim][a] * f.nodal[conn[a]];
          } break;
          case Basis2D::Q8: {
            std::array<std::array<double, 8>, 2> dN{};
            Shapes2D::Q8_dN(s2, dN[0].data(), dN[1].data());
            for (int a = 0; a < 8; ++a)
              for (int idim = 0; idim < 2; idim++)
                gradient[idim] += dN[idim][a] * f.nodal[conn[a]];
          } break;
          case Basis2D::Q9: {
            std::array<std::array<double, 9>, 2> dN{};
            Shapes2D::Q9_dN(s2, dN[0].data(), dN[1].data());
            for (int a = 0; a < 9; ++a)
              for (int idim = 0; idim < 2; idim++)
                gradient[idim] += dN[idim][a] * f.nodal[conn[a]];
          } break;
          default: return false;
          }

        }
        else if (DIM == 3) {
          // Only access shat[2] in the DIM==3 path
          Point3 s3{ shat[0], shat[1], shat[2] };

          switch (static_cast<Basis3D>(to_basis<DIM>(f.basis_id))) {
          case Basis3D::H8: {
            std::array<std::array<double, 8>, 3> dN{};
            Shapes3D::H8_dN(s3, dN[0].data(), dN[1].data(), dN[2].data());
            for (int a = 0; a < 8; ++a)
              for (int idim = 0; idim < 3; idim++)
                gradient[idim] += dN[idim][a] * f.nodal[conn[a]];
          } break;
          case Basis3D::H20: {
            std::array<std::array<double, 20>, 3> dN{};
            Shapes3D::H20_dN(s3, dN[0].data(), dN[1].data(), dN[2].data());
            for (int a = 0; a < 20; ++a)
              for (int idim = 0; idim < 3; idim++)
                gradient[idim] += dN[idim][a] * f.nodal[conn[a]];
          } break;
          case Basis3D::H27: {
            std::array<std::array<double, 27>, 3> dN{};
            Shapes3D::H27_dN(s3, dN[0].data(), dN[1].data(), dN[2].data());
            for (int a = 0; a < 27; ++a)
              for (int idim = 0; idim < 3; idim++)
                gradient[idim] += dN[idim][a] * f.nodal[conn[a]];
          } break;
          default: return false;
          }

        }
        else {
          return false; // unsupported DIM
        }

        const u32   L     = _tree_nodes[leaf_pos].level;
        const double scale = std::ldexp(1.0, int(L)); // 2^L
        for (std::size_t d = 0; d < DIM; ++d)
          gradient[d] = scale * gradient[d] ;

        return true;
      }



// ---- Locates the leaf for parent-space point s and maps it to that leaf’s local [-1,1]^DIM (DIM-independent) ----
      inline bool locate_leaf_on_parent_and_ref(const Point<DIM>& s,
                                                u32& leaf_node_idx,
                                                Point<DIM>& shat) const noexcept {
        // 1) Find the leaf containing s in parent space.
        leaf_node_idx = locate_leaf_on_parent(s);
        if (leaf_node_idx == npos32) return false;

        // 2) Convert to local [-1,1]^DIM using level + indices (no divisions).
        local_ref_fast(_tree_nodes[leaf_node_idx], s, shat);
        return true;
      }

// ---- refine so each point lies in a leaf by at least maxDepthTarget — caches only inverse-map per point (DIM-independent) ----
      void refine_to_contain_points(const std::vector<Point<DIM>>& pts, u32 maxDepthTarget) {
        if (pts.empty()) return;

        // ---- Cache inverse map once for all points ----
        struct XiEtaZeta {
          Point<DIM> xi;
          bool       valid;
        };
        std::vector<XiEtaZeta> pm(pts.size());

        for (size_t i = 0; i < pts.size(); ++i) {
          Point<DIM> s;
          bool ok = false;

          ok = inverse_map(pts[i], s);

          if (!ok) {
            pm[i] = { Point<DIM>{}, false };
          }
          else {
            pm[i] = { s, true };
          }
        }

        // ---- Refinement loop ----
        for (;;) {
          const u32 nleaves = (u32)leaf_count();
          if (nleaves == 0) break;

          // Keep node -> leaf_pos lookup synced with current _leaves.
          rebuild_leafpos_lookup(); // fills _node2leafpos

          // Mark which leaf positions need refinement this pass
          std::vector<uint8_t> mark_by_leaf(nleaves, 0);

          for (size_t i = 0; i < pm.size(); ++i) {
            if (!pm[i].valid) continue;

            // Reuse cached (ξ,η[,ζ]); do only the tree lookup per pass
            const u32 node = locate_leaf_on_parent(pm[i].xi);
            if (node == npos32) continue;

            const u32 lvl = _tree_nodes[node].level;
            if (lvl >= maxDepthTarget) continue;

            const u32 pos = _tree_nodes[node].node2leafIdx; // fast, no per-pass map build
            if (pos != npos32) mark_by_leaf[pos] = 1;
          }

          // Nothing to do? we're done.
          bool any = false;
          for (uint8_t b : mark_by_leaf) {
            if (b) {
              any = true;
              break;
            }
          }
          if (!any) break;

          // Minimal predicate: only checks the per-pass mark (by leaf position)
          auto should_refine_pos = [&](u32 leaf_pos) -> bool {
            return mark_by_leaf[leaf_pos] != 0;
          };

          // Count how many leaves are marked (we’ll refine at most that many)
          u32 to_split = 0;
          for (u32 pos = 0; pos < nleaves; ++pos) if (mark_by_leaf[pos]) ++to_split;

          // Net new nodes/leaves per split: (1<<DIM) - 1  (children minus the parent-as-leaf)
          const u32 add_per_split = (1u << DIM) - 1u;
          if (to_split) {
            _tree_nodes.reserve(_tree_nodes.size() + to_split * add_per_split);
            _leaves.reserve(_leaves.size() + to_split * add_per_split);
          }

          // Refine marked leaves (by position)
          const std::size_t nsplit = refine_pass_min(should_refine_pos);
          if (nsplit == 0) break; // reached max depth locally or topology unchanged

          // _leaves changed inside refine_pass_min; _node2leafpos will be rebuilt
          // at the top of the next loop iteration.
        }

        // Keep mesh balanced and refresh bookkeeping

        enforce_balance();
        post_topology_update();
      }

// ---- refine pass — splits each selected leaf into 2^DIM children — (DIM-independent) ----
      template<class RefinePred>
      std::size_t refine_pass_min(RefinePred&& should_refine) {
        std::vector<u32> newLeaves;
        newLeaves.reserve(_leaves.size() * 2u); // heuristic
        std::size_t refined = 0;

        const u32 n0 = static_cast<u32>(leaf_count()); // snapshot leaf count

        for (u32 idx = 0; idx < n0; ++idx) {
          const u32 leaf_node_idx = _leaves[idx];                      // index into _tree_nodes
          const TreeNode<DIM> leaf_copy = _tree_nodes[leaf_node_idx];  // snapshot (safe vs reallocation)

          if (leaf_copy.is_leaf && leaf_copy.level < _maxDepth && should_refine(idx)) {
            // Parent becomes internal; reset its child slots
            _tree_nodes[leaf_node_idx].is_leaf = false;
            for (std::size_t c = 0; c < TreeNode<DIM>::NCHILD; ++c)
              _tree_nodes[leaf_node_idx].child[c] = npos32;

            // Precompute base indices shifted by 1
            u32Point<DIM> base2 = leaf_copy.ix;
            for (std::size_t k = 0; k < DIM; ++k) base2[k] <<= 1;

            const u32 child_count = static_cast<u32>(TreeNode<DIM>::NCHILD); // = 1<<DIM
            for (u32 mask = 0; mask < child_count; ++mask) {
              // Child Morton index components
              u32Point<DIM> cix;
              for (std::size_t k = 0; k < DIM; ++k)
                cix[k] = base2[k] + ((mask >> k) & 1u);

              // Create child
              TreeNode<DIM> child{};
              child.ix      = cix;
              child.level   = leaf_copy.level + 1;
              child.morton  = interleave(cix);           // or derived().interleave(cix) in CRTP
              child.is_leaf = true;
              child.parent  = leaf_node_idx;
              for (std::size_t c = 0; c < TreeNode<DIM>::NCHILD; ++c)
                child.child[c] = npos32;

              const u32 child_idx = static_cast<u32>(_tree_nodes.size());
              _tree_nodes.push_back(child);

              // 'mask' is exactly the child code: bit k = side along axis k
              _tree_nodes[leaf_node_idx].child[mask] = child_idx;
              newLeaves.push_back(child_idx);
            }
            ++refined;
          }
          else {
            // Keep existing leaf
            newLeaves.push_back(leaf_node_idx);
          }
        }
        _leaves.swap(newLeaves);
        return refined;
      }

// ---- Simple cache-friendly FIFO ring buffer with advancing head (DIM-independent) ----
      struct RingQ {
        std::vector<u32> buf;
        std::size_t head = 0; // index of current front

        inline void clear() noexcept {
          buf.clear();
          head = 0;
        }
        inline bool empty() const noexcept {
          return head >= buf.size();
        }
        inline void push_back(u32 v) {
          buf.push_back(v);
        }
        inline u32 front() const noexcept {
          return buf[head];
        }
        inline void pop_front() noexcept {
          ++head;
        }
      };

// ---- Enforces 2:1 balance (1-irregularity): adjacent leaves may differ by at most one level (DIM-independent) ----
      inline void enforce_balance() {
        sort_leaves_by_level_then_morton(); // this is not needed, but it speeds up the leaf neighboor search
        compute_level_offsets();
        // local enqueue that avoids duplicates
        auto enqueue = [](RingQ & q, std::vector<char>& inq, u32 n) {
          if (n == npos32) return;
          if (n >= inq.size()) return;
          if (!inq[n]) {
            inq[n] = 1;
            q.push_back(n);
          }
        };

        // Seed queue with current leaves (bulk copy, no O(n) shifts later)
        RingQ q;
        q.buf.clear();
        q.head = 0;
        q.buf.reserve(_leaves.size() * (2 * DIM));   // heuristic proportional to face count
        q.buf.assign(_leaves.begin(), _leaves.end());

        // Track which nodes are in the queue (local to this call)
        std::vector<char> inq(_tree_nodes.size(), 0);
        for (u32 l : _leaves) if (l < inq.size()) inq[l] = 1;

        while (!q.empty()) {
          u32 leaf = q.front(); q.pop_front();
          if (leaf == npos32 || leaf >= _tree_nodes.size()) continue;
          inq[leaf] = 0;

          // Might have been refined since it was enqueued
          const auto& node = _tree_nodes[leaf];
          if (!node.is_leaf) continue;

          const u32 lev = node.level;

          // directions: 0..2*DIM-1  (pairs: -axis,+axis)
          for (int dir = 0; dir < int(2 * DIM); ++dir) {
            u32 nb = neighbor_leaf_by_face_any(leaf, dir);
            if (nb == npos32 || nb >= _tree_nodes.size()) continue;

            const auto& nbn = _tree_nodes[nb];
            if (!nbn.is_leaf) continue;

            const u32 lnb = nbn.level;

            // 2:1 rule: |lev - lnb| <= 1
            if (lev > lnb + 1) {
              // neighbor too coarse -> refine neighbor
              if (refine_leaf_once(nb)) {
                // containers may have grown; keep inq sized
                if (inq.size() < _tree_nodes.size()) inq.resize(_tree_nodes.size(), 0);

                const auto& refined_node = _tree_nodes[nb];  // `nb` is the index of the refined node
                for (int i = 0; i < (1 << DIM); ++i) {
                  u32 child_idx = refined_node.child[i];
                  if (child_idx == npos32) continue; // if uninitialized or invalid
                  enqueue(q, inq, child_idx);
                }
              }
            }
            else if (lnb > lev + 1) {
              // this leaf too coarse -> refine this leaf
              if (refine_leaf_once(leaf)) {
                if (inq.size() < _tree_nodes.size()) inq.resize(_tree_nodes.size(), 0);

                const auto& refined_node = _tree_nodes[leaf];  // `leaf` is the index of the refined node
                for (int i = 0; i < (1 << DIM); ++i) {
                  u32 child_idx = refined_node.child[i];
                  if (child_idx == npos32) continue; // if uninitialized or invalid
                  enqueue(q, inq, child_idx);
                }
                // no point checking other dirs for the old leaf; it’s not a leaf now
                break;
              }
            }
          }
        }
      }


// ---- Fast copy that reuses capacity without recomputing registries (DIM-independent)
      void assign_from(const OctTree& rhs) {
        if (this == &rhs) return;

        // ---- local helpers to copy while reusing capacity ----
        auto copy_vec = [](auto & dst, const auto & src) {
          if (dst.capacity() < src.size()) dst.reserve(src.size());
          dst.resize(src.size());
          std::copy(src.begin(), src.end(), dst.begin());
        };
        auto copy_vecvec_int = [](std::vector<std::vector<int>>& dst,
        const std::vector<std::vector<int>>& src) {
          if (dst.capacity() < src.size()) dst.reserve(src.size());
          dst.resize(src.size());
          for (size_t i = 0; i < src.size(); ++i) {
            if (dst[i].capacity() < src[i].size()) dst[i].reserve(src[i].size());
            dst[i].resize(src[i].size());
            std::copy(src[i].begin(), src[i].end(), dst[i].begin());
          }
        };
        auto copy_umap = [](auto & dst, const auto & src) {
          dst.clear(); // keep existing buckets
          if (dst.bucket_count() < src.size()) dst.rehash(src.size());
          for (const auto& kv : src) dst.insert(kv);
        };
        auto copy_uset = [](auto & dst, const auto & src) {
          dst.clear(); // keep existing buckets
          if (dst.bucket_count() < src.size()) dst.rehash(src.size());
          for (const auto& v : src) dst.insert(v);
        };

        // ---- scalars / config ----
        _maxDepth  = rhs._maxDepth;
        recompute_quant_params();

        _minDepth  = rhs._minDepth;
        _allowCoarsenBelowMinDepth = rhs._allowCoarsenBelowMinDepth;

        _root       = rhs._root;
        _geom_ready = rhs._geom_ready;

        // ---- geometry ----
        _X = rhs._X;

        // ---- active bases and fields ----
        copy_uset(_activeBases, rhs._activeBases);
        copy_vec(_fields,       rhs._fields);

        // ---- topology ----
        copy_vec(_tree_nodes,    rhs._tree_nodes);    // vector<TreeNode<DIM>>
        copy_vec(_leaves,        rhs._leaves);        // indices of leaf nodes
        copy_vec(_leaf_ids,      rhs._leaf_ids);      // cached leaf ids [0..n-1]
        //copy_vec(_node2leafpos,  rhs._node2leafpos);  // node index -> leaf position
        copy_vec(_level_offset,  rhs._level_offset);
        // ---- per-basis registries ----
        for (int b = 0; b < 3; ++b) {
          copy_umap(_basisReg[b].nodeMap, rhs._basisReg[b].nodeMap);            // u64 -> int
          copy_vec(_basisReg[b].nodes,    rhs._basisReg[b].nodes);              // vector<FEMNode<DIM>>
          copy_vecvec_int(_basisReg[b].elem2glob, rhs._basisReg[b].elem2glob);  // connectivity
        }
      }

// ---- conservative coarsen cycle using snapshot + parent coords; rebuild all fields — (DIM-independent) ----
      std::size_t coarsen_only_cycle_safe(u32 fid,
                                          std::vector <double> tolerance,
                                          OctTree& snapshot,
                                          u32 max_passes = 10) { // highest-order for DIM (Q9/H27)
        // Freeze current state for conservative evaluation/transfer
        snapshot.assign_from(*this);

        // Coarsen predicate (evaluated on the snapshot, in parent coords)
        auto pred = [&](u64 parent_morton, u32 level,
                        const std::vector<Point<DIM>>& pts_s,     // parent (ξ,η[,ζ])
                        const std::vector<Point<DIM>>& pts_xyz,
        const std::vector<std::array<double, NDOFS[DIM][2]>>& /*Nvals*/) -> bool {
          if (level <= min_depth() + 1u || level == max_depth()) return false;
          if (pts_s.empty()) return false;

          double v0;
          if (!snapshot.evaluate_field_on_parent(fid, pts_s[0], v0)) return false;

          double mn = v0, mx = v0;
          for (size_t i = 1; i < pts_s.size(); ++i) {
            double val;
            if (snapshot.evaluate_field_on_parent(fid, pts_s[i], val)) {
              mn = std::min(mn, val);
              mx = std::max(mx, val);
            }
          }
          return (mn > +tolerance[0]) || (mx < -tolerance[0]);
        };



        auto pred1 = [&](u64 /*parent_morton*/, u32 level,
                         const std::vector<Point<DIM>>& pts_s,     // Q9 (2D): 0..3 verts, 4..7 mids, 8 center
                         const std::vector<Point<DIM>>& /*pts_xyz*/,
        const std::vector<std::array<double, NDOFS[DIM][2]>>& /*Nvals*/) noexcept -> bool {
          if (level <= min_depth() + 1u) return false;
          if (pts_s.empty()) return false;

          // ---- small helpers (hot, inlinable) ----
          auto sq   = [](double x) noexcept {
            return x * x;
          };
          auto avg2 = [](double a, double b) noexcept {
            return 0.5 * (a + b);
          };

          double mn = std::numeric_limits<double>::infinity();
          double mx = -mn;

          // sample one point; updates mn/mx; returns false on failure
          auto sample = [&](size_t i, double & out) noexcept -> bool {
            double v;
            if (!snapshot.evaluate_field_on_parent(fid, pts_s[i], v)) return false;
            out = v;
            mn = (v < mn ? v : mn);
            mx = (v > mx ? v : mx);
            return true;
          };




          if constexpr(DIM == 2) {
            if (pts_s.size() < 9) return false;

            // sample all 9 (abort on first failure)
            std::array<double, 9> vQ9;
            for (size_t i = 0; i < 9; ++i) {
              if (!sample(i, vQ9[i])) return false;
            }

            // if sign change and not at deepest allowable level, do not accept
            if (mn * mx < 0.0 && level + 2u < max_depth()) return false;
            if (((mn > +tolerance[0]) || (mx < -tolerance[0])) && level + 2u < max_depth()) return true;

            // vertices
            const double v0 = vQ9[0], v1 = vQ9[1], v2 = vQ9[2], v3 = vQ9[3];

            // Q4 “reference” from vertex averages (edges 4..7), center 8
            double vQ4_edge[4];
            vQ4_edge[0] = avg2(v0, v1); // edge(0-1) -> node 4
            vQ4_edge[1] = avg2(v1, v2); // edge(1-2) -> node 5
            vQ4_edge[2] = avg2(v2, v3); // edge(2-3) -> node 6
            vQ4_edge[3] = avg2(v3, v0); // edge(3-0) -> node 7

            // center: 0.25*(0+1+2+3) == 0.5*(4+6) == 0.5*(5+7)
            const double vQ4_center = avg2(vQ4_edge[0], vQ4_edge[2]);

            // error: mids (weight 1) + center (weight 2)
            double e = 0.0;
            e += sq(vQ9[4] - vQ4_edge[0]);
            e += sq(vQ9[5] - vQ4_edge[1]);
            e += sq(vQ9[6] - vQ4_edge[2]);
            e += sq(vQ9[7] - vQ4_edge[3]);
            e += 2.0 * sq(vQ9[8] - vQ4_center);

            // element measure (constant, as in your code)
            const double meas = 6.0 / sq(mx - mn);
            e *= meas;
            //std::cout<<e<<" ";
            return e < tolerance[1];
          }
          else {
            if (pts_s.size() < 27) return false;

            // sample all 27 (abort on first failure)
            std::array<double, 27> vH27;
            for (size_t i = 0; i < 27; ++i) {
              if (!sample(i, vH27[i])) return false;
            }

            // if sign change and not at deepest allowable level, do not accept
            if (mn * mx < 0.0 && level + 2u < max_depth()) return false;
            if (((mn > +tolerance[0]) || (mx < -tolerance[0])) && level + 2u < max_depth()) return true;

            // vertices (0..7)
            const double v0 = vH27[0], v1 = vH27[1], v2 = vH27[2], v3 = vH27[3];
            const double v4 = vH27[4], v5 = vH27[5], v6 = vH27[6], v7 = vH27[7];

            // H8 edge mids from vertex avgs (reference), H27 edge DOFs are 8..19
            double vH8_edge[12];
            vH8_edge[0]  = avg2(v0, v1);
            vH8_edge[1]  = avg2(v1, v2);
            vH8_edge[2]  = avg2(v2, v3);
            vH8_edge[3]  = avg2(v3, v0);
            vH8_edge[4]  = avg2(v4, v5);
            vH8_edge[5]  = avg2(v5, v6);
            vH8_edge[6]  = avg2(v6, v7);
            vH8_edge[7]  = avg2(v7, v4);
            vH8_edge[8]  = avg2(v0, v4);
            vH8_edge[9]  = avg2(v1, v5);
            vH8_edge[10] = avg2(v2, v6);
            vH8_edge[11] = avg2(v3, v7);

            // H8 face centers from edge mids (reference), H27 face DOFs are 20..25
            double vH8_face[6];
            vH8_face[0] = avg2(vH8_edge[8],  vH8_edge[9]);   // (0,1,5,4)
            vH8_face[1] = avg2(vH8_edge[9],  vH8_edge[10]);  // (1,2,6,5)
            vH8_face[2] = avg2(vH8_edge[10], vH8_edge[11]);  // (2,3,7,6)
            vH8_face[3] = avg2(vH8_edge[11], vH8_edge[8]);   // (3,0,4,7)
            vH8_face[4] = avg2(vH8_edge[0],  vH8_edge[2]);   // bottom (0,1,2,3)
            vH8_face[5] = avg2(vH8_edge[4],  vH8_edge[6]);   // top    (4,5,6,7)

            // H8 center
            const double vH8_center = avg2(vH8_face[4], vH8_face[5]);

            // error accumulators (weights as in your code)
            double e = 0.0;
            auto add_edge_err = [&](unsigned i) noexcept {
              const double d = vH27[i + 8] - vH8_edge[i];
              e += d * d;
            };
            auto add_face_err = [&](unsigned i) noexcept {
              const double d = vH27[i + 20] - vH8_face[i];
              e += 2.0 * d * d;
            };

            for (unsigned i = 0; i < 12; ++i) add_edge_err(i);
            for (unsigned i = 0; i < 6;  ++i) add_face_err(i);
            {
              const double d = vH27[26] - vH8_center;
              e += 8.0 * d * d;
            }

            const double meas = 32. / sq(mx - mn);
            e *= meas;

            return e < tolerance[1];
          }
        };

        std::size_t total = 0;

        if (tolerance.size() == 1) {
          // // Coarsen passes; (balance+topology update happen after the loop)
          for (u32 pass = 0; pass < max_passes; ++pass) {
            std::size_t c = coarsen_pass(pred);
            if (c == 0) break;
            total += c;
          }
        }
        else if (tolerance.size() > 1) {
          // Coarsen passes; (balance+topology update happen after the loop)
          for (u32 pass = 0; pass < max_passes; ++pass) {
            std::size_t c = coarsen_pass(pred1);
            if (c == 0) break;
            total += c;
          }
        }



        // Enforce 1-irregularity and refresh bookkeeping

        enforce_balance();

        post_topology_update();

        // Rebuild all fields from the snapshot (conservative transfer) on nodal sets
        for (u32 f = 0; f < _fields.size(); ++f) {
          rebuild_field_from(snapshot, f);   // ensure rebuild_field_from is DIM-aware (2D/3D)
        }

        return total;
      }

// ---- rebuild field from another tree by sampling at parent nodes — (DIM-independent) ----
      void rebuild_field_from(const OctTree& src, u32 fid) {
        static_assert(DIM == 2 || DIM == 3, "rebuild_field_from supports DIM = 2 or 3");
        assert(fid < _fields.size() && fid < src._fields.size());

        Field&       dst = _fields[fid];
        const Field& s   = src._fields[fid];
        assert(dst.basis_id == s.basis_id);

        const BasisRegistry<DIM>& Rdst = _basisReg[(int)dst.basis_id];
        const BasisRegistry<DIM>& Rsrc = src._basisReg[(int)s.basis_id];

        dst.nodal.resize(Rdst.nodes.size());

        for (size_t gid = 0; gid < Rdst.nodes.size(); ++gid) {
          const auto& pr = Rdst.nodes[gid].parent;             // parent coord (ξ,η[,ζ])
          const u64 key = morton_node_key_from_parent(pr);

          // Fast path: direct copy if the node exists in src's registry
          auto it = Rsrc.nodeMap.find(key);
          if (it != Rsrc.nodeMap.end()) {
            dst.nodal[gid] = s.nodal[it->second];
            continue;
          }

          // Fallback: sample from src at the exact parent coord (works for DIM=2 or 3)
          double val = std::numeric_limits<double>::quiet_NaN();
          Point<DIM> p{}; for (std::size_t k = 0; k < DIM; ++k) p[k] = pr[k];
          const bool ok = src.evaluate_field_on_parent(fid, p, val);
          dst.nodal[gid] = ok ? val : std::numeric_limits<double>::quiet_NaN();
        }
      }

// ---- Returns the element-to-global connectivity for the given basis (DIM-independent) ----
      inline const std::vector<std::vector<int>>& basis_connectivity(Basis b) const {
        return _basisReg[static_cast<int>(b)].elem2glob;
      }

// ---- Returns the unique global FEM nodes (with parent & physical coords) for the given basis (DIM-independent) ----
      inline const std::vector<FEMNode<DIM>>& basis_nodes(Basis b) const {
        return _basisReg[static_cast<int>(b)].nodes;
      }

// ---- VTU writer: automatic 2D→3D vector padding for ParaView + cell "level" ----
      bool write_binary_vtu(const std::string &filename, Basis b,
                            const std::vector<std::variant<u32, std::vector<u32>>> &field_groups,
                            bool cell_centered = false) const {

        require_geometry();

        const BasisRegistry<DIM> &R = _basisReg[(int)b];
        const size_t numCells = R.elem2glob.size();
        if (numCells == 0) {
          std::cerr << "[VTU] Warning: No elements registered for basis " << (int)b << "\n";
          return false;
        }

        // -------- Points (global) --------
        std::vector<double> flatPoints;
        flatPoints.reserve(R.nodes.size() * 3);
        for (const auto &n : R.nodes) {
          flatPoints.push_back(n.physical[0]);
          flatPoints.push_back(n.physical[1]);
          if constexpr(DIM == 3)
            flatPoints.push_back(n.physical[2]);
          else
            flatPoints.push_back(0.0); // 2D → embed in z=0 plane
        }

        // -------- Connectivity --------
        std::vector<int> connectivity;
        std::vector<int> offsets;
        std::vector<unsigned char> types;

        auto vtk_cell_type = [&](Basis bb) -> unsigned char {
          if constexpr(DIM == 3) {
            switch (static_cast<Basis3D>(bb)) {
            case Basis3D::H8:  return 12;
            case Basis3D::H20: return 25;
            case Basis3D::H27: return 29;
            }
          }
          else {
            switch (static_cast<Basis2D>(bb)) {
            case Basis2D::Q4: return 9;
            case Basis2D::Q8: return 23;
            case Basis2D::Q9: return 28;
            }
          }
          return 0;
        };

        const int perCell = leaf_dof_number(b);
        connectivity.reserve(numCells * perCell);
        offsets.reserve(numCells);
        types.reserve(numCells);

        const std::array<int, 7> tail_map{{3, 1, 0, 2, 4, 5, 6}};
        int offset = 0;

        for (size_t e = 0; e < numCells; ++e) {
          const auto &conn = R.elem2glob[e];
          if (DIM == 3 && static_cast<Basis3D>(b) == Basis3D::H27) {
            connectivity.insert(connectivity.end(), conn.begin(), conn.begin() + 20);
            for (int i = 0; i < 7; ++i) connectivity.push_back(conn[20 + tail_map[i]]);
            offset += 27;
          }
          else {
            connectivity.insert(connectivity.end(), conn.begin(), conn.end());
            offset += static_cast<int>(conn.size());
          }
          offsets.push_back(offset);
          types.push_back(vtk_cell_type(b));
        }

        // -------- Prepare per-cell "level" from _leaves --------
        std::vector<int> cellLevel;
        cellLevel.reserve(numCells);
        const size_t nmap = std::min(numCells, _leaves.size());
        for (size_t e = 0; e < nmap; ++e) {
          const u32 leaf_node = _leaves[e];
          cellLevel.push_back(static_cast<int>(_tree_nodes[leaf_node].level));
        }
        if (nmap < numCells) {
          std::cerr << "[VTU] Warning: _leaves.size() (" << _leaves.size()
                    << ") < NumberOfCells (" << numCells
                    << "); padding remaining levels with 0.\n";
          cellLevel.insert(cellLevel.end(), numCells - nmap, 0);
        }

        // -------- XML output --------
        std::ofstream os(filename);
        if (!os) {
          std::cerr << "[VTU] Error: cannot open file " << filename << "\n";
          return false;
        }

        os << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        os << "  <UnstructuredGrid>\n";
        os << "    <Piece NumberOfPoints=\"" << R.nodes.size()
           << "\" NumberOfCells=\"" << numCells << "\">\n";

        os << "      <Points>\n";
        write_binary_array(os, "Float64", "", 3, flatPoints);
        os << "      </Points>\n";

        os << "      <Cells>\n";
        write_binary_array(os, "Int32",  "connectivity", 1, connectivity);
        write_binary_array(os, "Int32",  "offsets",      1, offsets);
        write_binary_array(os, "UInt8",  "types",        1, types);
        os << "      </Cells>\n";

        // -------- Attribute names for ParaView (per the section we will write) --------
        auto detect_first_names = [&](bool for_cell_section) {
          std::string first_scalar_name, first_vector_name;
          for (const auto &entry : field_groups) {
            if (std::holds_alternative<u32>(entry)) {
              const u32 fid = std::get<u32>(entry);
              if (fid < _fields.size()) {
                const Field &fld = _fields[fid];
                if (to_basis<DIM>(fld.basis_id) == b && first_scalar_name.empty()) {
                  first_scalar_name = fld.name; // same label even if averaged to cells
                }
              }
            }
            else {
              const auto &group = std::get<std::vector<u32>>(entry);
              if (group.size() == DIM && first_vector_name.empty()) {
                const Field &f0 = _fields[group[0]];
                std::string nm = f0.name;
                for (size_t i = 1; i < group.size(); ++i) nm += "-" + _fields[group[i]].name;
                first_vector_name = nm; // same label; layout differs by section
              }
            }
          }
          return std::pair<std::string, std::string> {first_scalar_name, first_vector_name};
        };

        // -------- Write Data Sections --------
        if (!cell_centered) {
          // --- PointData: nodal fields ---
          auto [p_scalar, p_vector] = detect_first_names(false);
          os << "      <PointData";
          if (!p_scalar.empty()) os << " Scalars=\"" << p_scalar << "\"";
          if (!p_vector.empty()) os << " Vectors=\"" << p_vector << "\"";
          os << ">\n";

          // Scalars/vectors at nodes
          for (const auto &entry : field_groups) {
            if (std::holds_alternative<u32>(entry)) {
              const u32 fid = std::get<u32>(entry);
              if (fid >= _fields.size()) continue;
              const Field &fld = _fields[fid];
              if (to_basis<DIM>(fld.basis_id) != b) continue;
              write_binary_array(os, "Float64", fld.name, 1, fld.nodal);
            }
            else {
              const auto &group = std::get<std::vector<u32>>(entry);
              if (group.size() != DIM) continue;

              std::vector<const Field*> vecFld;
              vecFld.reserve(DIM);
              bool valid = true;
              for (u32 fid : group) {
                if (fid >= _fields.size()) {
                  valid = false;
                  break;
                }
                const Field &f = _fields[fid];
                if (to_basis<DIM>(f.basis_id) != b || f.nodal.size() != R.nodes.size()) {
                  valid = false; break;
                }
                vecFld.push_back(&f);
              }
              if (!valid) continue;

              std::string vname = vecFld.front()->name;
              for (size_t i = 1; i < vecFld.size(); ++i) vname += "-" + vecFld[i]->name;

              std::vector<double> combined;
              combined.reserve(R.nodes.size() * 3);
              for (size_t i = 0; i < R.nodes.size(); ++i) {
                combined.push_back(vecFld[0]->nodal[i]);
                if constexpr(DIM > 1) combined.push_back(vecFld[1]->nodal[i]);
                else                   combined.push_back(0.0);
                if constexpr(DIM == 3) combined.push_back(vecFld[2]->nodal[i]);
                else                    combined.push_back(0.0);
              }
              write_binary_array(os, "Float64", vname, 3, combined);
            }
          }
          os << "      </PointData>\n";

          // --- CellData: just the per-cell "level" ---
          os << "      <CellData Scalars=\"level\">\n";
          write_binary_array(os, "Int32", "level", 1, cellLevel);
          os << "      </CellData>\n";
        }
        else {
          // --- CellData: averaged cell fields + per-cell "level" ---
          auto [c_scalar, c_vector] = detect_first_names(true);
          os << "      <CellData";
          // We always have 'level'; choose user scalar name if present, else "level"
          if (!c_scalar.empty()) os << " Scalars=\"" << c_scalar << "\"";
          else                   os << " Scalars=\"level\"";
          if (!c_vector.empty()) os << " Vectors=\"" << c_vector << "\"";
          os << ">\n";

          // Scalars: average nodal to cell
          for (const auto &entry : field_groups) {
            if (std::holds_alternative<u32>(entry)) {
              const u32 fid = std::get<u32>(entry);
              if (fid >= _fields.size()) continue;
              const Field &fld = _fields[fid];
              if (to_basis<DIM>(fld.basis_id) != b) continue;

              std::vector<double> cellData;
              cellData.reserve(numCells);
              for (const auto &conn : R.elem2glob) {
                double v = 0.0;
                for (int gid : conn) v += fld.nodal[gid];
                cellData.push_back(v / static_cast<double>(conn.size()));
              }
              write_binary_array(os, "Float64", fld.name, 1, cellData);
            }
            else {
              const auto &group = std::get<std::vector<u32>>(entry);
              if (group.size() != DIM) continue;

              std::vector<const Field*> vecFld;
              vecFld.reserve(DIM);
              bool valid = true;
              for (u32 fid : group) {
                if (fid >= _fields.size()) {
                  valid = false;
                  break;
                }
                const Field &f = _fields[fid];
                if (to_basis<DIM>(f.basis_id) != b || f.nodal.size() != R.nodes.size()) {
                  valid = false; break;
                }
                vecFld.push_back(&f);
              }
              if (!valid) continue;

              std::string vname = vecFld.front()->name;
              for (size_t i = 1; i < vecFld.size(); ++i) vname += "-" + vecFld[i]->name;

              // Average components to cell centers, then 3-pad
              std::vector<double> cx, cy, cz;
              cx.reserve(numCells);
              if constexpr(DIM > 1) cy.reserve(numCells);
              if constexpr(DIM == 3) cz.reserve(numCells);

              for (const auto &conn : R.elem2glob) {
                auto avg_comp = [&](const Field * fptr) {
                  double v = 0.0; for (int gid : conn) v += fptr->nodal[gid];
                  return v / static_cast<double>(conn.size());
                };
                cx.push_back(avg_comp(vecFld[0]));
                if constexpr(DIM > 1) cy.push_back(avg_comp(vecFld[1]));
                if constexpr(DIM == 3) cz.push_back(avg_comp(vecFld[2]));
              }

              std::vector<double> combined3;
              combined3.reserve(numCells * 3);
              for (size_t i = 0; i < numCells; ++i) {
                combined3.push_back(cx[i]);
                if constexpr(DIM > 1) combined3.push_back(cy[i]); else combined3.push_back(0.0);
                if constexpr(DIM == 3) combined3.push_back(cz[i]); else combined3.push_back(0.0);
              }
              write_binary_array(os, "Float64", vname, 3, combined3);
            }
          }

          // Append the per-cell "level"
          write_binary_array(os, "Int32", "level", 1, cellLevel);
          os << "      </CellData>\n";
        }

        os << "    </Piece>\n";
        os << "  </UnstructuredGrid>\n";
        os << "</VTKFile>\n";

        std::cout << "[VTU] Wrote " << filename << " with "
                  << field_groups.size() << " field groups and cell levels.\n";
        return true;
      }

// ---- write VTU for all fields with the given Basis (simple scalar listing) ----
      bool write_binary_vtu(const std::string &filename, Basis b,
                            bool cell_centered = false) const {
        using FG = std::variant<u32, std::vector<u32>>;

        std::vector<FG> field_groups;
        field_groups.reserve(_fields.size());

        for (u32 i = 0; i < _fields.size(); ++i) {
          if (to_basis<DIM>(_fields[i].basis_id) == b) {
            field_groups.emplace_back(i); // add as scalar
          }
        }

        // Delegate to the main overload (which also writes per-cell "level")
        return write_binary_vtu(filename, b, field_groups, cell_centered);
      }

// ---- global maximum depth allowed (DIM-independent) ----
      u32 max_depth() const {
        return _maxDepth;
      }

// ---- Returns the X-coordinates array of the coarse geometry (DIM-independent) ----
      inline const std::array<double, NDOFS[DIM][2]>& get_X() const noexcept {
        return _X[0];
      }

// ---- Returns the Y-coordinates array of the coarse geometry (DIM-independent) ----
      inline const std::array<double, NDOFS[DIM][2]>& get_Y() const noexcept {
        return _X[1];
      }

// ---- Returns the Z-coordinates array of the coarse geometry; asserts DIM == 3 (DIM-independent) ----
// ---- Only participates in overload resolution when DIM==3
      template <std::size_t D = DIM, typename std::enable_if<D == 3, int>::type = 0>
      inline const std::array<double, NDOFS[DIM][2]>& get_Z() const noexcept {
        return _X[2];
      }

// ---- Collect parent-space node coordinates from leaves with level in [lev_min, lev_max] (DIM-independent)
      inline void extract_node_parent_coords_in_level_range(u32 lev_min,
                                                            u32 lev_max,
                                                            Basis basis,
                                                            std::vector<Point<DIM>>& coords) const {
        coords.clear();

        // Need level offsets prepared
        if (_level_offset.empty()) return;

        // Levels present are [0 .. max_level_present]; index max_level_present+1 is the sentinel
        const u32 max_level_present =
          static_cast<u32>(_level_offset.size() >= 2 ? _level_offset.size() - 2 : 0);

        // Clamp and normalize input range
        u32 Lmin = (lev_min > max_level_present) ? max_level_present : lev_min;
        u32 Lmax = (lev_max > max_level_present) ? max_level_present : lev_max;
        if (Lmin > Lmax) std::swap(Lmin, Lmax);

        // Compute half-open span [lo, hi)
        const u32 lo = _level_offset[Lmin];
        const u32 hi = _level_offset[Lmax + 1];

        // Empty span? nothing to do
        if (lo >= hi) return;

        const u32 ndofs_per_leaf = leaf_dof_number(basis);

        // Conservative reserve
        coords.reserve(static_cast<size_t>(hi - lo) * ndofs_per_leaf);

        // Scratch buffer for per-leaf parent nodes; avoid realloc by reserving once
        std::vector<Point<DIM>> s;
        s.reserve(ndofs_per_leaf);

        for (u32 leaf_pos = lo; leaf_pos < hi; ++leaf_pos) {
          s.clear();
          extract_leaf_parent_coords(basis, leaf_pos, s);
          if (!s.empty()) {
            coords.insert(coords.end(), s.begin(), s.end());
          }
        }
      }

// ---- For each leaf (level in [lev_min, lev_max]), return parent-space coords of the Q9/H27 nodes where field fid reaches min and max on that leaf (DIM-independent) ----
      inline void extract_parent_coords_by_field_extremes_in_level_range(u32 lev_min,
                                                                         u32 lev_max,
                                                                         u32 fid,
                                                                         std::vector<Point<DIM>>& out) const {
        out.clear();

        // Must have offsets computed: offsets[L] = start of level L, offsets[maxL+1] = sentinel
        if (_level_offset.empty()) return;

        // Valid level range present in this tree:
        // If size is K, levels are [0 .. K-2], and index K-1 is the sentinel (== _leaves.size()).
        const u32 max_level_present = static_cast<u32>(_level_offset.size() >= 2 ? _level_offset.size() - 2 : 0);

        // Clamp and normalize input range
        u32 Lmin = (lev_min > max_level_present) ? max_level_present : lev_min;
        u32 Lmax = (lev_max > max_level_present) ? max_level_present : lev_max;
        if (Lmin > Lmax) std::swap(Lmin, Lmax);

        // Nothing to do?
        if (_level_offset[Lmin] == _level_offset[Lmax + 1]) return;

        // Reserve: at most 2 points per leaf (min+max), for all leaves in [Lmin, Lmax]
        const u32 lo = _level_offset[Lmin];
        const u32 hi = _level_offset[Lmax + 1];
        const u32 nleafs = hi - lo;
        out.reserve(static_cast<size_t>(nleafs) * 2u);

        // Always sample extremes on the highest geometry basis (Q9/H27)
        const Basis basisHi = (DIM == 2)
                              ? static_cast<Basis>(Basis2D::Q9)
                              : static_cast<Basis>(Basis3D::H27);

        for (u32 leaf_pos = lo; leaf_pos < hi; ++leaf_pos) {
          std::vector<Point<DIM>> s;
          extract_leaf_parent_coords(basisHi, leaf_pos, s);
          if (s.empty()) continue;

          double vmin =  std::numeric_limits<double>::infinity();
          double vmax = -std::numeric_limits<double>::infinity();
          u32 imin = 0, imax = 0;

          for (u32 i = 0; i < static_cast<u32>(s.size()); ++i) {
            double val = std::numeric_limits<double>::quiet_NaN();
            (void)evaluate_field_on_parent(fid, s[i], val);
            if (val < vmin) {
              vmin = val;
              imin = i;
            }
            if (val > vmax) {
              vmax = val;
              imax = i;
            }
          }

          out.push_back(s[imin]);
          if (imax != imin) out.push_back(s[imax]);
        }
      }

// ---- Collect parent-space *centers* from leaves with level in [lev_min, lev_max] (DIM-independent)
      inline void extract_leaf_centers_in_level_range(u32 lev_min,
                                                      u32 lev_max,
                                                      std::vector<Point<DIM>>& centers) const {
        centers.clear();

        // Need level offsets available
        if (_level_offset.empty()) return;

        // Levels present are [0 .. max_level_present]; index max_level_present+1 is the sentinel
        const u32 max_level_present =
          static_cast<u32>(_level_offset.size() >= 2 ? _level_offset.size() - 2 : 0);

        // Clamp and normalize input range
        u32 Lmin = (lev_min > max_level_present) ? max_level_present : lev_min;
        u32 Lmax = (lev_max > max_level_present) ? max_level_present : lev_max;
        if (Lmin > Lmax) std::swap(Lmin, Lmax);

        // No leaves in the range? (start == end)
        if (_level_offset[Lmin] == _level_offset[Lmax + 1]) return;

        const u32 lo = _level_offset[Lmin];
        const u32 hi = _level_offset[Lmax + 1];

        centers.reserve(static_cast<size_t>(hi - lo));

        for (u32 leaf_pos = lo; leaf_pos < hi; ++leaf_pos) {
          Point<DIM> s_center{};
          extract_leaf_parent_center_coord(leaf_pos, s_center); // Q9/H27 center internally
          centers.push_back(s_center);
        }
      }

      int leaf_dof_number(Basis b) const noexcept {
        return NDOFS[DIM][static_cast<uint8_t>(b)];
      }

    private:
// ---- Computes the axis-aligned bounds [x0, x1] of leaf node n in parent space [-1,1]^DIM (DIM-independent) ----
      inline void leaf_bounds(const TreeNode<DIM>& n,
                              Point<DIM>& x0,
                              Point<DIM>& x1) const noexcept {
        // cell size d = 2 / 2^level = 2^(1 - level)
        const double d = std::ldexp(1.0, 1 - int(n.level));
        for (std::size_t k = 0; k < DIM; ++k) {
          const double a = std::fma(static_cast<double>(n.ix[k]), d, -1.0); // -1 + ix[k]*d
          x0[k] = a;
          x1[k] = a + d;
        }
      }

// ---- Maps parent-space point s to the local coordinates shat in the leaf n’s reference cube [-1,1]^DIM (DIM-independent) ----
      inline void local_ref_fast(const TreeNode<DIM>& n,
                                 const Point<DIM>& s,
                                 Point<DIM>& shat) const noexcept {
        const int    L = static_cast<int>(n.level); // refinement level
        const double N = std::ldexp(1.0, L);        // N = 2^L
        for (std::size_t k = 0; k < DIM; ++k) {
          const double bias = N - double((n.ix[k] << 1) + 1); // N - (2*ix[k]+1)
          shat[k] = std::ldexp(s[k], L) + bias;               // N*s[k] + bias
        }
      }

// ---- Ensures all leaves are refined down to at least _minDepth (DIM-independent) ----
      inline void ensure_min_depth() {
        if (_minDepth == 0) return;

        bool changed = true;
        while (changed) {
          changed = false;

          std::vector<u32> newLeaves;
          newLeaves.reserve(_leaves.size() * 2u); // heuristic

          for (u32 leaf_idx : _leaves) {
            // Copy out (vector may reallocate during refinement)
            const TreeNode<DIM> node_copy = _tree_nodes[leaf_idx];

            if (node_copy.level < _minDepth) {
              // Mark parent as internal and clear child slots
              _tree_nodes[leaf_idx].is_leaf = false;
              for (std::size_t c = 0; c < TreeNode<DIM>::NCHILD; ++c)
                _tree_nodes[leaf_idx].child[c] = npos32;

              // Precompute base indices shifted by 1
              u32Point<DIM> base2 = node_copy.ix;
              for (std::size_t k = 0; k < DIM; ++k) base2[k] <<= 1;

              // Create all 2^DIM children
              const u32 child_count = static_cast<u32>(TreeNode<DIM>::NCHILD);
              for (u32 mask = 0; mask < child_count; ++mask) {
                // Child logical index
                u32Point<DIM> cix;
                for (std::size_t k = 0; k < DIM; ++k)
                  cix[k] = base2[k] + ((mask >> k) & 1u);

                // New child node
                TreeNode<DIM> child{};
                child.ix      = cix;
                child.level   = node_copy.level + 1;
                child.is_leaf = true;
                child.parent  = leaf_idx;
                for (std::size_t c = 0; c < TreeNode<DIM>::NCHILD; ++c)
                  child.child[c] = npos32;

                // Morton key (CRTP: Derived may provide interleave; or use this->interleave)
                child.morton = interleave(cix);

                const u32 cindex = static_cast<u32>(_tree_nodes.size());
                _tree_nodes.push_back(child);
                _tree_nodes[leaf_idx].child[mask] = cindex;
                newLeaves.push_back(cindex);
              }

              changed = true;
            }
            else {
              // Already deep enough
              newLeaves.push_back(leaf_idx);
            }
          }

          if (changed) _leaves.swap(newLeaves);
        }
      }

// ---- Returns the refinement level of a given leaf position (DIM-independent) ----
      inline u32 level_of(u32 leaf_pos) const noexcept {
        assert(leaf_pos < _leaves.size());
        const u32 node_idx = _leaves[leaf_pos];
        return _tree_nodes[node_idx].level;
      }

// ---- Inserts (or retrieves) a grid node for parent point s in [-1,1]^DIM; returns its global id (DIM-independent) ----
      inline int get_or_insert_gid(BasisRegistry<DIM>& R, const Point<DIM>& s) {
        const u64 key = morton_node_key_from_parent(s);

        auto it = R.nodeMap.find(key);
        if (it != R.nodeMap.end()) return it->second;

        const int gid = (int)R.nodes.size();
        R.nodes.push_back(FEMNode<DIM> {gid, s, parent_to_physical(s)});
        R.nodeMap.emplace(key, gid);
        return gid;
      }

// ---- Hash functor for dimension-selected Basis enum (DIM-independent) ----
      struct BasisHasher {
        template <typename E>
        std::size_t operator()(E b) const noexcept {
          static_assert(std::is_enum<E>::value, "BasisHasher expects an enum type");
          using U = typename std::underlying_type<E>::type;
          return std::hash<U> {}(static_cast<U>(b));
        }
      };

// ---- Recompute node-quantization parameters from _maxDepth (DIM-independent) ---
      inline void recompute_quant_params() noexcept {
        constexpr u32 kMaxBitsPerAxis = (DIM == 3) ? 21u : 32u; // bits per axis supported by 64-bit Morton
        assert(_maxDepth <= kMaxBitsPerAxis - 1u && "Depth exceeds Morton capacity for this DIM");
        _nodeBitsQ  = std::max<u32>(2u, std::min<u32>(_maxDepth + 1u, kMaxBitsPerAxis));
        _nodesNQ    = (1u << _nodeBitsQ);
        _halfNodesN = double(_nodesNQ) * 0.5;
        _maxIdxQ    = _nodesNQ - 1u;
      }

// ---- Quantize parent coord s∈[-1,1] to integer node index on 2^{_maxDepth+1} grid (DIM-independent) ---
      inline u32 quantize_parent_coord(double s) const noexcept {
        if (s <= -1.0) return 0u;
        if (s >=  1.0) return _maxIdxQ;

        // Map [-1,1] → [0, _nodesNQ] and round-to-nearest.
        // We avoid llround; t is non-negative so floor(t+0.5) is fine.
        double t = (s + 1.0) * _halfNodesN;

        // Guard against tiny FP drift at the ends.
        if (t <= 0.0) return 0u;
        if (t >= double(_nodesNQ)) return _maxIdxQ;

        u64 iu = (u64)(t + 0.5);
        if (iu > _maxIdxQ) iu = _maxIdxQ; // handle t==_nodesNQ after rounding
        return (u32)iu;
      }

// ---- Build Morton key for a parent-space point using nodal quantization (DIM-independent) ---
      inline u64 morton_node_key_from_parent(const Point<DIM>& s) const noexcept {
        u32Point<DIM> ix{};
        for (std::size_t k = 0; k < DIM; ++k) ix[k] = quantize_parent_coord(s[k]);
        return interleave(ix);
      }

// ---- Determinante di J in s (usa _X e le tue shape derivate Q9/H27)
      double detJ_at(const Point<DIM>& s) const noexcept {
        require_geometry();

        // buffer derivati
        constexpr std::size_t Nnodes = NDOFS[DIM][2]; // Q9/H27
        std::vector<double> dNx(Nnodes), dNy(Nnodes), dNz((DIM == 3) ? Nnodes : 0);

        if constexpr(DIM == 2) {
          Point2 s2{ s[0], s[1] };
          Shapes2D::Q9_dN(s2, dNx.data(), dNy.data());
          double J[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
          for (std::size_t a = 0; a < Nnodes; ++a) {
            const double dNa0 = dNx[a], dNa1 = dNy[a];
            const double X0 = _X[0][a], X1 = _X[1][a];
            J[0][0] += dNa0 * X0;  J[0][1] += dNa1 * X0;
            J[1][0] += dNa0 * X1;  J[1][1] += dNa1 * X1;
          }
          return J[0][0] * J[1][1] - J[0][1] * J[1][0];
        }
        else {   // DIM == 3
          Point3 s3{ s[0], s[1], s[2] };
          Shapes3D::H27_dN(s3, dNx.data(), dNy.data(), dNz.data());
          double J[3][3] = {{0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0}
          };
          for (std::size_t a = 0; a < Nnodes; ++a) {
            const double dNa0 = dNx[a], dNa1 = dNy[a], dNa2 = dNz[a];
            const double X0 = _X[0][a], X1 = _X[1][a], X2 = _X[2][a];
            // J(i,j) = sum_a dN_j(a) * X_i(a)
            J[0][0] += dNa0 * X0;  J[0][1] += dNa1 * X0;  J[0][2] += dNa2 * X0;
            J[1][0] += dNa0 * X1;  J[1][1] += dNa1 * X1;  J[1][2] += dNa2 * X1;
            J[2][0] += dNa0 * X2;  J[2][1] += dNa1 * X2;  J[2][2] += dNa2 * X2;
          }
          // det 3x3
          const double A11 = J[0][0], A12 = J[0][1], A13 = J[0][2];
          const double A21 = J[1][0], A22 = J[1][1], A23 = J[1][2];
          const double A31 = J[2][0], A32 = J[2][1], A33 = J[2][2];
          return A11 * (A22 * A33 - A23 * A32)
                 - A12 * (A21 * A33 - A23 * A31)
                 + A13 * (A21 * A32 - A22 * A31);
        }
      }

    public:
// ---- Determinante di J in un punto parent s (wrapper pubblico, const) ----
      double detJ_parent(const Point<DIM>& s) const noexcept {
        return detJ_at(s);  // chiama la privata
      }

    private: //data
// config
      u32  _maxDepth;
      u32  _minDepth{0};
      bool _allowCoarsenBelowMinDepth{true};

// Store coordinates as _G[axis][a], where axis = 0:x,1:y,(2:z if DIM==3) and a = node index.
      std::array<std::array<double, NDOFS[DIM][2]>, DIM> _X{};
      bool _geom_ready{false};

// mesh leaves
      std::vector<u32> _leaves;
      std::vector<u32> _level_offset;
      mutable std::vector<u32> _leaf_ids;

// fields (per-leaf coefficients for API compatibility)
      std::vector<Field> _fields;

// per-basis global node registries + connectivity
      std::vector<TreeNode<DIM>> _tree_nodes;        // full tree hierarchy
      u32 _root{0};

      //mutable std::vector<u32> _node2leafpos; // size == _tree_nodes.size(), npos32 default
      std::unordered_set<Basis, BasisHasher> _activeBases;

      BasisRegistry<DIM> _basisReg[3]; // 0:Q4/H8, 1:Q8/H20, 2:Q9/H27

// Quantization params (depend only on DIM and _maxDepth)
      u32 _nodeBitsQ{0};
      u32 _nodesNQ{0};
      double _halfNodesN{0.0};   // = _nodesNQ * 0.5
      u32 _maxIdxQ{0};           // = _nodesNQ - 1
  };

// ============================================================
// Compute integral of sign(phi) over the domain using public API
// ============================================================
// template<std::size_t DIM>
// double compute_sign_integral(const fem::OctTree<DIM>& tree, fem::u32 field_id) {
//   using namespace fem;
//   double integral = 0.0;
//
//   // 1. Ottieni riferimento al campo
//   const Field& fld = tree.field(field_id);
//   const auto basis = to_basis<DIM>(fld.basis_id);
//
//   // 2. Ottieni riferimento al registro del basis
//   const auto& registry = tree.basis_nodes(basis);
//   const auto& elem2glob = tree.basis_connectivity(basis);
//
//   // 3. Cicla sugli elementi (leaves)
//   for (const auto& conn : elem2glob) {
//     // Ottieni le coordinate fisiche dei nodi
//     std::vector<Point<DIM>> x(conn.size());
//     std::vector<double> phi(conn.size());
//     for (size_t i = 0; i < conn.size(); ++i) {
//       x[i] = registry[conn[i]].physical;
//       phi[i] = fld.nodal[conn[i]];
//     }
//
//     // 4. Calcola un'approssimazione semplice dell'integrale sul leaf
//     //    (es. media del segno * area/volume approssimato)
//     double avg_sign = 0.0;
//     for (double p : phi) avg_sign += (p >= 0.0 ? 1.0 : -1.0);
//     avg_sign /= static_cast<double>(phi.size());
//
//     // Stima dimensionale dell'elemento
//     double measure = 1.0;
//     if constexpr (DIM == 2) {
//       // area ≈ (Δx * Δy)
//       double dx = std::abs(x[1][0] - x[0][0]);
//       double dy = std::abs(x[3][1] - x[0][1]);
//       measure = dx * dy;
//     } else {
//       // volume ≈ (Δx * Δy * Δz)
//       double dx = std::abs(x[1][0] - x[0][0]);
//       double dy = std::abs(x[3][1] - x[0][1]);
//       double dz = std::abs(x[4][2] - x[0][2]);
//       measure = dx * dy * dz;
//     }
//
//     integral += avg_sign * measure;
//   }
//
//   return integral;
// }
// ------- mass and geometric error -------
  template<std::size_t DIM>
  std::pair<double, double>
  compute_mass_and_geom_errors_overlay(
    const OctTree<DIM>& overlay,
    const OctTree<DIM>& tree_t0, u32 fid_t0,
    const OctTree<DIM>& tree_t1, u32 fid_t1
  ) {

    const u32 Nel = overlay.leaf_count();
    double M0 = 0.0, M1 = 0.0;
    double Eg = 0.0;
    double Vtot = 0.0;

    std::vector<fem::Point<DIM>> gauss3_pts;
    std::vector<double> w3;
    std::vector<fem::Point<DIM>> gauss7_pts;
    std::vector<double> w7;

    if constexpr(DIM == 2) {
      gauss3_pts.reserve(9);  w3.reserve(9);
      gauss7_pts.reserve(49); w7.reserve(49);
    }
    else if constexpr(DIM == 3) {
      gauss3_pts.reserve(27);  w3.reserve(27);
      gauss7_pts.reserve(343); w7.reserve(343);
    }

    // ===============================
    // COSTRUZIONE QUADRATURE
    // ===============================
    const double a = std::sqrt(3.0 / 5.0);
    if constexpr(DIM == 2) {
      {
        // ---- 3x3 ----
        const double pts[3] = {-a, 0.0, +a};
        const double w[3] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j) {
            gauss3_pts.push_back({pts[i], pts[j]});
            w3.push_back(w[i]*w[j]);
          }
      }
      {
        // ---- 7x7 ----
        const double x[7] = {-0.9491079123, -0.7415311856, -0.4058451514, 0.0, 0.4058451514,  0.7415311856,  0.9491079123};
        const double w[7] = {0.1294849662, 0.2797053915, 0.3818300505, 0.4179591837, 0.3818300505, 0.2797053915, 0.1294849662};
        for (int i = 0; i < 7; ++i)
          for (int j = 0; j < 7; ++j) {
            gauss7_pts.push_back({x[i], x[j]});
            w7.push_back(w[i]*w[j]);
          }
      }
    }
    else if constexpr(DIM == 3) {
      {
        // ---- 3x3x3 ----
        const double pts[3] = {-a, 0.0, +a};
        const double w[3] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k) {
              gauss3_pts.push_back({pts[i], pts[j], pts[k]});
              w3.push_back(w[i]*w[j]*w[k]);
            }
      }
      {
        // ---- 7x7x7 ----
        const double x[7] = {-0.9491079123, -0.7415311856, -0.4058451514, 0.0, 0.4058451514,  0.7415311856,  0.9491079123};
        const double w[7] = {0.1294849662, 0.2797053915, 0.3818300505, 0.4179591837, 0.3818300505, 0.2797053915, 0.1294849662};
        for (int i = 0; i < 7; ++i)
          for (int j = 0; j < 7; ++j)
            for (int k = 0; k < 7; ++k) {
              gauss7_pts.push_back({x[i], x[j], x[k]});
              w7.push_back(w[i]*w[j]*w[k]);
            }
      }
    }

    // ===============================
    // LOOP SULLE CELLE
    // ===============================
    for (u32 e = 0; e < Nel; ++e) {
      Point<DIM> X0, X1;
      std::vector<Point<DIM>> pts;
      if constexpr(DIM == 2)
        overlay.extract_leaf_parent_coords(Basis2D::Q4, e, pts);
      else if constexpr(DIM == 3)
        overlay.extract_leaf_parent_coords(Basis3D::H8, e, pts);

      // Bounds in parent space (min/max sui nodi Q4/H8 del leaf)
      for (int d = 0; d < DIM; ++d) {
        X0[d] = std::numeric_limits<double>::max();
        X1[d] = -std::numeric_limits<double>::max();
      }
      for (const auto& p : pts) {
        for (int d = 0; d < DIM; ++d) {
          X0[d] = std::min(X0[d], p[d]);
          X1[d] = std::max(X1[d], p[d]);
        }
      }

      // Guardie anti-degenerazione (facoltative ma utili)
      if ((X1[0] <= X0[0]) || (X1[1] <= X0[1]) || (DIM == 3 && (X1[2] <= X0[2]))) {
        continue; // salta cella degenerata
      }


      // === Determina se la cella è tagliata ===
      double phi_min0 = 1e30, phi_max0 = -1e30;
      double phi_min1 = 1e30, phi_max1 = -1e30;
      for (const auto& p : pts) {
        double phi0, phi1;
        tree_t0.evaluate_field_on_parent(fid_t0, p, phi0);
        tree_t1.evaluate_field_on_parent(fid_t1, p, phi1);
        phi_min0 = std::min(phi_min0, phi0);
        phi_max0 = std::max(phi_max0, phi0);
        phi_min1 = std::min(phi_min1, phi1);
        phi_max1 = std::max(phi_max1, phi1);
      }
      bool cut_cell = ((phi_min0 * phi_max0 < 0.0) || (phi_min1 * phi_max1 < 0.0));

      // === Seleziona quadratura ===
      const auto& gauss_pts = cut_cell ? gauss7_pts : gauss3_pts;
      const auto& weights   = cut_cell ? w7        : w3;

      // === Loop sui punti di Gauss ===
      double integ_sign0 = 0.0, integ_sign1 = 0.0;
      double Jacc = 0.0;
      const double Jparent = 0.5 * (X1[0] - X0[0])
                             * 0.5 * (X1[1] - X0[1])
                             * ((DIM == 3) ? 0.5 * (X1[2] - X0[2]) : 1.0);


      for (size_t q = 0; q < gauss_pts.size(); ++q) {
        const auto& s = gauss_pts[q];
        const double wq = weights[q];
        Point<DIM> s_parent;
        for (int d = 0; d < DIM; ++d) s_parent[d] = 0.5 * ((1.0 - s[d]) * X0[d] + (1.0 + s[d]) * X1[d]);

        const double detJ_phys = std::fabs(overlay.detJ_parent(s_parent));
        const double w_phys = wq * detJ_phys * Jparent;

        double phi0 = 0.0, phi1 = 0.0;
        tree_t0.evaluate_field_on_parent(fid_t0, s_parent, phi0);
        tree_t1.evaluate_field_on_parent(fid_t1, s_parent, phi1);

        integ_sign0 += w_phys * ((phi0 >= 0.0) ? 1.0 : -1.0);
        integ_sign1 += w_phys * ((phi1 >= 0.0) ? 1.0 : -1.0);
        Jacc        += w_phys;
      }

      // Media del segno → color function
      const double C0 = 0.5 * (1.0 + integ_sign0 / Jacc);
      const double C1 = 0.5 * (1.0 + integ_sign1 / Jacc);
      const double Ai = Jacc;
      M0 += Ai * C0;
      M1 += Ai * C1;
      Eg += Ai * std::fabs(C1 - C0);
      Vtot += Ai;
    }

    const double Em = (std::fabs(M0) > 1e-14) ? std::fabs(M1 - M0) / std::fabs(M0) : 0.0;
    const double Eg_rel = (std::fabs(M0) > 1e-14) ? Eg / std::fabs(M0) : 0.0;
    return {Em, Eg_rel};
  }

} // namespace fem























