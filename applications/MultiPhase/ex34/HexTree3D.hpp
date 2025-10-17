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


#include "Encoder.hpp"

namespace fem {

  using u32 = uint32_t;
  using u64 = uint64_t;
  constexpr u32 npos32 = std::numeric_limits<u32>::max();

  using Point3 = std::array<double, 3>;

  enum class Basis : uint8_t { H8 = 0, H20 = 1, H27 = 2 };

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
  inline uint64_t interleave3(u32 x, u32 y, u32 z) {
    uint64_t xx = expand3_21(x);
    uint64_t yy = expand3_21(y) << 1;
    uint64_t zz = expand3_21(z) << 2;
    return (zz | yy | xx);
  }


  struct Shapes3 {
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
    static inline void H8(const Point3& s, double* __restrict__ N) noexcept {
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
    }
  };

// ---------- Tree node & helpers ----------
  struct TreeNode3D {
    u64 morton{0};
    u32 ix{0}, iy{0}, iz{0};
    u32 level{0};
    u32 child[8] {npos32, npos32, npos32, npos32, npos32, npos32, npos32, npos32};
    bool is_leaf{true};
    u32 parent{npos32};
  };

  struct Field {
    Basis basis{Basis::H8};
    // Nodal values, one per unique global node in the basis registry.
    std::vector<double> nodal;  // size == _basisReg[(int)basis].nodes.size()

  };

  struct FEMNode3D {
    int gid;
    Point3 parent;
    Point3 physical;
  };

// ---------- HexTree3D (minimum workable stub) ----------
  class HexTree3D {
    public:
      HexTree3D(u32 maxDepth, u32 minDepth = 0)
        : _maxDepth(maxDepth),
          _minDepth(std::min(minDepth, maxDepth)) {

        // We interleave 21 bits per axis: limit depth accordingly.
        assert(_maxDepth <= 20 && "HexTree3D supports up to 20 levels per axis (expand3_21).");

        // --- 1) Build root TreeNode covering [-1,1]^3 ---
        TreeNode3D root{};
        root.ix = 0;
        root.iy = 0;
        root.iz = 0;
        root.level  = 0;
        root.morton = interleave3(0, 0, 0);
        root.is_leaf = true;
        root.parent  = npos32;
        for (int c = 0; c < 8; ++c) root.child[c] = npos32;

        // --- 2) Reset topology containers and seed with root ---
        _tree_nodes.clear();
        _tree_nodes.push_back(root);
        _root = 0;

        _leaves.clear();
        _leaves.push_back(_root);

        _leaf_ids.clear();
        _node2leafpos.clear();
        _leafpos_valid = false;

        // --- 3) Basis registries empty ---
        for (int b = 0; b < 3; ++b) {
          _basisReg[b].clear();
        }

        // --- 4) Geometry not set yet ---
        _geom_ready = false;
      }

      // --- reset (3D) ---
      void reset(bool keep_geometry = true, bool keep_fields = true) {
        // 1) optionally keep geometry
        if (!keep_geometry) {
          _X = {}; _Y = {}; _Z = {};
          _geom_ready = false;
        }

        // 2) optionally keep fields (values + active bases), otherwise nuke them
        if (!keep_fields) {
          _fields.clear();
          _activeBases.clear();
          for (int b = 0; b < 3; ++b) _basisReg[b].clear();
        }
        else {
          // we keep _activeBases and _fields, but their registries will be rebuilt
          for (int b = 0; b < 3; ++b) _basisReg[b].clear();
        }

        // 3) reset topology to just the root (reuse capacity; no shrink)
        _leaves.clear();
        _tree_nodes.clear();

        TreeNode3D root{};
        root.ix = 0; root.iy = 0; root.iz = 0;
        root.morton  = interleave3(0, 0, 0);
        root.level   = 0;
        root.is_leaf = true;
        root.parent  = npos32;
        for (int c = 0; c < 8; ++c) root.child[c] = npos32;

        _tree_nodes.push_back(root);
        _root = 0;
        _leaves.push_back(_root);

        // 4) clear/refresh lookups (reuse capacity)
        _leaf_ids.clear();
        _node2leafpos.clear();
        rebuild_leafpos_lookup();

        // 5) basis registries / field sizes depend on active bases
        post_topology_update(); // rebuild connectivity for active bases, resize fields if kept
      }

// --- special members (3D) ---
      HexTree3D(HexTree3D&&) noexcept = default;
      HexTree3D& operator=(HexTree3D&&) noexcept = default;
      HexTree3D(const HexTree3D&) = delete;
      HexTree3D& operator=(const HexTree3D&) = delete;

// --- swap (3D) ---
      friend void swap(HexTree3D& a, HexTree3D& b) noexcept {
        using std::swap;
        swap(a._maxDepth, b._maxDepth);
        swap(a._minDepth, b._minDepth);
        swap(a._allowCoarsenBelowMinDepth, b._allowCoarsenBelowMinDepth);
        swap(a._leafpos_valid, b._leafpos_valid);

        swap(a._X, b._X);
        swap(a._Y, b._Y);
        swap(a._Z, b._Z);
        swap(a._geom_ready, b._geom_ready);

        swap(a._leaves, b._leaves);
        swap(a._leaf_ids, b._leaf_ids);
        swap(a._fields, b._fields);

        swap(a._basisReg[0], b._basisReg[0]);
        swap(a._basisReg[1], b._basisReg[1]);
        swap(a._basisReg[2], b._basisReg[2]);

        swap(a._tree_nodes, b._tree_nodes);
        swap(a._root, b._root);
        swap(a._node2leafpos, b._node2leafpos);
        swap(a._activeBases, b._activeBases);
      }

// ---- config (3D) ----
      void set_min_depth(u32 d) {
        _minDepth = std::min(d, _maxDepth);
      }
      u32 min_depth() const {
        return _minDepth;
      }
      void set_allow_coarsen_below_min(bool v) {
        _allowCoarsenBelowMinDepth = v;
      }
      bool allow_coarsen_below_min() const {
        return _allowCoarsenBelowMinDepth;
      }

// --- geometry setters / guards (3D) ---
      void set_physical_hex27(const double X27[27],
                              const double Y27[27],
                              const double Z27[27]) {
        for (int i = 0; i < 27; ++i) {
          _X[i] = X27[i];
          _Y[i] = Y27[i];
          _Z[i] = Z27[i];
        }
        _geom_ready = true;
      }

      inline void require_geometry() const {
        assert(_geom_ready && "Call set_physical_hex27(...) before geometric operations.");
      }

      inline void post_topology_update() {
        // Keep leaves ordered by 3D Morton (Z-order)
        std::sort(_leaves.begin(), _leaves.end(),
        [&](u32 a, u32 b) {
          return _tree_nodes[a].morton < _tree_nodes[b].morton;
        });

        rebuild_leafpos_lookup();            // refresh node -> leaf position map
        rebuild_connectivity_active_bases(); // rebuild per-basis registries/connectivity
        resize_fields_to_nodes();            // sync field storage to current node registries
      }

      inline void rebuild_leafpos_lookup() {
        _node2leafpos.assign(_tree_nodes.size(), npos32);
        for (u32 i = 0; i < static_cast<u32>(_leaves.size()); ++i) {
          const u32 node_idx = _leaves[i];
          if (node_idx < _node2leafpos.size()) _node2leafpos[node_idx] = i;
        }
        _leafpos_valid = true;
      }

      void rebuild_connectivity_active_bases() {
        if (!_geom_ready) return;
        for (Basis b : _activeBases) {
          rebuild_connectivity(b);
        }
      }

      void rebuild_connectivity(Basis b) {
        BasisRegistry &R = _basisReg[(int)b];
        R.clear();

        // --- Reserve to avoid rehash/realloc churn ---
        const u32 nleaves = leaf_count();
        const int per = (b == Basis::H8  ? 8 :
                         b == Basis::H20 ? 20 : 27); // H27

        // A decent estimate for unique nodes is ~ (per * leaves) * 0.6..0.8 due to sharing.
        // Keep it conservative (1.5x leaves*per) to minimize rehashes:
        const size_t est_nodes =
          std::max<size_t>(32, size_t(nleaves) * size_t(per) * size_t(3) / size_t(2)); // 1.5x
        R.nodeMap.reserve(est_nodes);
        R.nodes.reserve(est_nodes);
        R.elem2glob.reserve(nleaves);

        std::vector<Point3> s; // parent-space (xi,eta,zeta) nodes per leaf

        for (u32 e = 0; e < nleaves; ++e) {
          extract_leaf_parent_coords(b, e, s); // fills s with per-node parent coords for leaf e

          std::vector<int> conn;
          conn.reserve(s.size());

          for (const auto& p : s) {
            // 3D registry insert/lookup (quantized by _maxDepth+1)
            // Assumes you have: int get_or_insert_gid(BasisRegistry&, double xi, double eta, double zeta)

            int gid = get_or_insert_gid(R, p[0], p[1], p[2]);
            conn.push_back(gid);
          }

          R.elem2glob.push_back(std::move(conn));
        }
      }



      // Resize all fields so their nodal arrays match current basis registries
      void resize_fields_to_nodes() {
        for (auto& f : _fields) {
          const auto& R = _basisReg[(int)f.basis];
          f.nodal.resize(R.nodes.size(), 0.0);
        }
      }


      // ---- node generators (parent & physical) ----
// Fill vector with parent-space coordinates of interpolation nodes (by basis) for a leaf
      void extract_leaf_parent_coords(Basis basis, u32 leaf_idx,
                                      std::vector<Point3>& out_pts) const {
        // Map leaf_idx → actual TreeNode3D
        const TreeNode3D& leaf = _tree_nodes[_leaves[leaf_idx]];

        // Leaf bounds in parent space
        double x0, y0, z0, x1, y1, z1;
        leaf_bounds_3d(leaf, x0, y0, z0, x1, y1, z1);

        const double xm = 0.5 * (x0 + x1);
        const double ym = 0.5 * (y0 + y1);
        const double zm = 0.5 * (z0 + z1);

        out_pts.clear();

        switch (basis) {
        case Basis::H8: {
          out_pts.reserve(8);
          // corners 0..7
          out_pts = {
            Point3{x0, y0, z0}, Point3{x1, y0, z0}, Point3{x1, y1, z0}, Point3{x0, y1, z0},
            Point3{x0, y0, z1}, Point3{x1, y0, z1}, Point3{x1, y1, z1}, Point3{x0, y1, z1}
          };
        } break;

        case Basis::H20: {
          out_pts.reserve(20);
          // corners 0..7
          out_pts = {
            Point3{x0, y0, z0}, Point3{x1, y0, z0}, Point3{x1, y1, z0}, Point3{x0, y1, z0},
            Point3{x0, y0, z1}, Point3{x1, y0, z1}, Point3{x1, y1, z1}, Point3{x0, y1, z1}
          };
          // z = -1 edges (8..11)
          out_pts.push_back(Point3{xm, y0, z0});
          out_pts.push_back(Point3{x1, ym, z0});
          out_pts.push_back(Point3{xm, y1, z0});
          out_pts.push_back(Point3{x0, ym, z0});
          // z = +1 edges (12..15)
          out_pts.push_back(Point3{xm, y0, z1});
          out_pts.push_back(Point3{x1, ym, z1});
          out_pts.push_back(Point3{xm, y1, z1});
          out_pts.push_back(Point3{x0, ym, z1});
          // vertical edges (16..19)
          out_pts.push_back(Point3{x0, y0, zm});
          out_pts.push_back(Point3{x1, y0, zm});
          out_pts.push_back(Point3{x1, y1, zm});
          out_pts.push_back(Point3{x0, y1, zm});
        } break;

        case Basis::H27: {
          out_pts.reserve(27);
          // corners 0..7
          out_pts = {
            Point3{x0, y0, z0}, Point3{x1, y0, z0}, Point3{x1, y1, z0}, Point3{x0, y1, z0},
            Point3{x0, y0, z1}, Point3{x1, y0, z1}, Point3{x1, y1, z1}, Point3{x0, y1, z1}
          };
          // z = -1 edges (8..11)
          out_pts.push_back(Point3{xm, y0, z0});
          out_pts.push_back(Point3{x1, ym, z0});
          out_pts.push_back(Point3{xm, y1, z0});
          out_pts.push_back(Point3{x0, ym, z0});
          // z = +1 edges (12..15)
          out_pts.push_back(Point3{xm, y0, z1});
          out_pts.push_back(Point3{x1, ym, z1});
          out_pts.push_back(Point3{xm, y1, z1});
          out_pts.push_back(Point3{x0, ym, z1});
          // vertical edges (16..19)
          out_pts.push_back(Point3{x0, y0, zm});
          out_pts.push_back(Point3{x1, y0, zm});
          out_pts.push_back(Point3{x1, y1, zm});
          out_pts.push_back(Point3{x0, y1, zm});
          // face centers (20..25)
          out_pts.push_back(Point3{xm, y0, zm}); // x-mid on y=− face, z-mid
          out_pts.push_back(Point3{x1, ym, zm}); // y-mid on x=+ face, z-mid
          out_pts.push_back(Point3{xm, y1, zm}); // x-mid on y=+ face, z-mid
          out_pts.push_back(Point3{x0, ym, zm}); // y-mid on x=− face, z-mid
          out_pts.push_back(Point3{xm, ym, z0}); // center of bottom face
          out_pts.push_back(Point3{xm, ym, z1}); // center of top face
          // cell center (26)
          out_pts.push_back(Point3{xm, ym, zm});
        } break;
        }
      }

      // map parent -> physical using H27 geometry
      Point3 parent_to_physical(const Point3& s) const {
        require_geometry();
        double N[27];
        Shapes3::H27(s, N);
        double X = 0, Y = 0, Z = 0;
        for (int a = 0; a < 27; ++a) {
          X += N[a] * _X[a];
          Y += N[a] * _Y[a];
          Z += N[a] * _Z[a];
        }
        return {X, Y, Z};
      }

      // ---- leaves / indices ----
      u32 leaf_count() const {
        return static_cast<u32>(_leaves.size());
      }

      const std::vector<u32>& leaf_indices() const {
        _leaf_ids.resize(_leaves.size());
        for (u32 i = 0; i < static_cast<u32>(_leaves.size()); ++i) _leaf_ids[i] = i;
        return _leaf_ids;
      }

// ---- fields (per-leaf coefficients) ----
      u32 add_field(Basis b) {
        Field f;
        f.basis = b;

        // ensure registry exists so we can size nodal storage
        _activeBases.insert(b);
        post_topology_update();

        const auto& R = _basisReg[static_cast<int>(b)];
        f.nodal.assign(R.nodes.size(), 0.0);

        _fields.push_back(std::move(f));
        return static_cast<u32>(_fields.size() - 1);
      }

      Field& field(u32 fid) {
        return _fields[fid];
      }
      const Field& field(u32 fid) const {
        return _fields[fid];
      }







      // inverse map (Newton, H27)
      // Newton inverse for H27 with H8 warm-up (analogous to quad9 w/ quad4 warm-up)
      bool inverse_map_hex27(const std::array<double, 3>& x,
                             std::array<double, 3>& s,
                             int maxIts = 30, double tol = 1e-12) const {
        require_geometry(); // make sure _X,_Y,_Z are set

        s = {0.0, 0.0, 0.0};
        const double tol2 = tol * tol;

        // Buffers
        double N8[8], dN8dx[8], dN8dy[8], dN8dz[8];
        double N27[27], dN27dx[27], dN27dy[27], dN27dz[27];

        // --- Warm-up with H8 at current s (start center) ---
        Shapes3::H8(s, N8);

        double Xp = 0.0, Yp = 0.0, Zp = 0.0;
        for (int a = 0; a < 8; ++a) {
          Xp += N8[a] * _X[a];
          Yp += N8[a] * _Y[a];
          Zp += N8[a] * _Z[a];
        }

        double rx = Xp - x[0], ry = Yp - x[1], rz = Zp - x[2];

        Shapes3::H8_dN(s, dN8dx, dN8dy, dN8dz);
        double dXdx = 0, dXdy = 0, dXdz = 0, dYdx = 0, dYdy = 0, dYdz = 0, dZdx = 0, dZdy = 0, dZdz = 0;
        for (int a = 0; a < 8; ++a) {
          dXdx += dN8dx[a] * _X[a]; dXdy += dN8dy[a] * _X[a]; dXdz += dN8dz[a] * _X[a];
          dYdx += dN8dx[a] * _Y[a]; dYdy += dN8dy[a] * _Y[a]; dYdz += dN8dz[a] * _Y[a];
          dZdx += dN8dx[a] * _Z[a]; dZdy += dN8dy[a] * _Z[a]; dZdz += dN8dz[a] * _Z[a];
        }

        for (int it = 0; it < maxIts; ++it) {
          // Solve J * ds = r   (same layout as your 2D: we step by subtracting ds)
          // J = [dXdx dXdy dXdz; dYdx dYdy dYdz; dZdx dZdy dZdz]
          const double A11 = dXdx, A12 = dXdy, A13 = dXdz;
          const double A21 = dYdx, A22 = dYdy, A23 = dYdz;
          const double A31 = dZdx, A32 = dZdy, A33 = dZdz;

          const double det =
            A11 * (A22 * A33 - A23 * A32) -
            A12 * (A21 * A33 - A23 * A31) +
            A13 * (A21 * A32 - A22 * A31);

          if (std::fabs(det) < 1e-20) break;
          const double inv = 1.0 / det;

          // adjugate * r
          const double dxs =
            ((A22 * A33 - A23 * A32) * rx - (A12 * A33 - A13 * A32) * ry + (A12 * A23 - A13 * A22) * rz) * inv;
          const double dys =
            (-(A21 * A33 - A23 * A31) * rx + (A11 * A33 - A13 * A31) * ry - (A11 * A23 - A13 * A21) * rz) * inv;
          const double dzs =
            ((A21 * A32 - A22 * A31) * rx - (A11 * A32 - A12 * A31) * ry + (A11 * A22 - A12 * A21) * rz) * inv;

          // update iterate
          s[0] = std::max(-1.0, std::min(1.0, s[0] - dxs));
          s[1] = std::max(-1.0, std::min(1.0, s[1] - dys));
          s[2] = std::max(-1.0, std::min(1.0, s[2] - dzs));

          // Recompute residual with H27 at updated s
          Shapes3::H27(s, N27);
          Xp = Yp = Zp = 0.0;
          for (int a = 0; a < 27; ++a) {
            Xp += N27[a] * _X[a];
            Yp += N27[a] * _Y[a];
            Zp += N27[a] * _Z[a];
          }
          rx = Xp - x[0]; ry = Yp - x[1]; rz = Zp - x[2];
          const double nrm2 = rx * rx + ry * ry + rz * rz;
          if (nrm2 < tol2) return true;

          // Build NEXT Jacobian with H27 derivatives
          Shapes3::H27_dN(s, dN27dx, dN27dy, dN27dz);
          dXdx = dXdy = dXdz = dYdx = dYdy = dYdz = dZdx = dZdy = dZdz = 0.0;
          for (int a = 0; a < 27; ++a) {
            dXdx += dN27dx[a] * _X[a]; dXdy += dN27dy[a] * _X[a]; dXdz += dN27dz[a] * _X[a];
            dYdx += dN27dx[a] * _Y[a]; dYdy += dN27dy[a] * _Y[a]; dYdz += dN27dz[a] * _Y[a];
            dZdx += dN27dx[a] * _Z[a]; dZdy += dN27dy[a] * _Z[a]; dZdz += dN27dz[a] * _Z[a];
          }
        }
        return false; // not converged
      }


      // locate leaf containing parent point s ([-1,1]^3)
      u32 locate_leaf_on_parent(const Point3& s) const {
        if (s[0] < -1.0 || s[0] > 1.0 || s[1] < -1.0 || s[1] > 1.0 || s[2] < -1.0 || s[2] > 1.0) return npos32;
        const double scale = double(1u << _maxDepth) / 2.0;
        u32 ix = std::min<u32>((u32)((s[0] + 1.0) * scale), (1u << _maxDepth) - 1);
        u32 iy = std::min<u32>((u32)((s[1] + 1.0) * scale), (1u << _maxDepth) - 1);
        u32 iz = std::min<u32>((u32)((s[2] + 1.0) * scale), (1u << _maxDepth) - 1);
        u32 node = _root;
        for (u32 level = 0; level < _maxDepth; ++level) {
          const TreeNode3D& n = _tree_nodes[node];
          if (n.is_leaf) return node;
          const int shift = int(_maxDepth - level - 1);
          const int qx = (ix >> shift) & 1u, qy = (iy >> shift) & 1u, qz = (iz >> shift) & 1u;
          const int q = (qz << 2) | (qy << 1) | qx;
          u32 c = n.child[q];
          if (c == npos32) return node;
          node = c;
        }
        return node;
      }


// ---- physical nodes for a leaf ----
      void leaf_physical_nodes(Basis basis, u32 leaf_idx,
                               std::vector<Point3>& out_xyz) const {
        require_geometry();
        std::vector<Point3> s;
        extract_leaf_parent_coords(basis, leaf_idx, s);
        out_xyz.resize(s.size());
        for (size_t i = 0; i < s.size(); ++i) out_xyz[i] = parent_to_physical(s[i]);
      }

// ---- adaptation ---------------------------------------------------
      template<class CoarsenPred, class RefinePred>
      std::size_t adapt_cycle(CoarsenPred&& should_coarsen,
                              RefinePred&&  should_refine,
                              Basis probe_basis = Basis::H27,
                              u32 max_passes   = 10) {
        ensure_min_depth(); // refine floor only
        std::size_t total = 0;
        for (u32 pass = 0; pass < max_passes; ++pass) {
          std::size_t changed = 0;
          changed += coarsen_pass(should_coarsen, probe_basis);
          changed += refine_pass(should_refine, probe_basis);
          if (changed == 0) break;
          total += changed;
        }
        post_topology_update();
        return total;
      }

// Perform one refinement pass: split leaves if predicate says so
      template<class RefinePred>
      std::size_t refine_pass(RefinePred&& should_refine, Basis probe_basis) {
        std::vector<u32> newLeaves;
        newLeaves.reserve(_leaves.size() * 2u); // heuristic (not all leaves refine)
        std::size_t refined = 0;

        for (u32 idx = 0; idx < leaf_count(); ++idx) {
          u32 leaf_node_idx = _leaves[idx];                // index into _tree_nodes
          const TreeNode3D leaf_copy = _tree_nodes[leaf_node_idx]; // snapshot

          // probe geometry for refinement
          std::vector<Point3> pts_s, pts_xyz;
          extract_leaf_parent_coords(probe_basis, idx, pts_s);
          leaf_physical_nodes(probe_basis, idx, pts_xyz);
          std::vector<std::array<double, 27>> dummy; // placeholder, if predicate wants it

          if (should_refine(idx, pts_s, pts_xyz, dummy) && leaf_copy.level < _maxDepth) {
            // split into 8 children
            const u32 ix = leaf_copy.ix;
            const u32 iy = leaf_copy.iy;
            const u32 iz = leaf_copy.iz;

            for (int dz = 0; dz < 2; ++dz) {
              for (int dy = 0; dy < 2; ++dy) {
                for (int dx = 0; dx < 2; ++dx) {
                  const int q = (dz << 2) | (dy << 1) | dx;

                  TreeNode3D child{};
                  child.ix = (ix << 1) + (u32)dx;
                  child.iy = (iy << 1) + (u32)dy;
                  child.iz = (iz << 1) + (u32)dz;
                  child.morton  = interleave3(child.ix, child.iy, child.iz);
                  child.level   = leaf_copy.level + 1;
                  child.is_leaf = true;
                  child.parent  = leaf_node_idx;
                  for (int c = 0; c < 8; ++c) child.child[c] = npos32;

                  const u32 child_idx = (u32)_tree_nodes.size();
                  _tree_nodes.push_back(child);

                  // re-fetch parent AFTER push_back
                  _tree_nodes[leaf_node_idx].child[q] = child_idx;
                  newLeaves.push_back(child_idx);
                }
              }
            }

            _tree_nodes[leaf_node_idx].is_leaf = false;  // mark parent
            ++refined;
          }
          else {
            newLeaves.push_back(leaf_node_idx);
          }
        }

        _leaves.swap(newLeaves);
        return refined;
      }

      template<class CoarsenPred>
      std::size_t coarsen_pass(CoarsenPred&& should_coarsen, Basis probe_basis) {
        struct Group {
          u32      pos[8];      // positions in _leaves of the 8 children
          uint8_t  cnt   = 0;
          uint32_t tag   = 0;   // generation tag
        };

        const u32 L = (u32)_leaves.size();
        const u32 N = (u32)_tree_nodes.size();

        // --- thread-local scratch to avoid alloc churn ---
        static thread_local std::vector<u32> touched;
        static thread_local std::vector<u32> to_add;
        static thread_local std::vector<u32> newLeaves;
        static thread_local std::vector<Group> per;
        static thread_local uint32_t per_epoch = 1;

        if (per.size() < N) per.resize(N);        // grow once, never shrink
        ++per_epoch;                               // bump generation

        touched.clear(); touched.reserve(L / 8 + 8);
        to_add.clear();  to_add.reserve(L / 8 + 8);
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
        std::array<Point3, 27> s_buf, xyz_buf;

        // Reuse vector wrappers for the predicate (avoid allocs)
        static thread_local std::vector<Point3> pts_s_vec, pts_xyz_vec;

        // how many probe points?
        const int nprobe = (probe_basis == Basis::H8  ? 8 :
                            (probe_basis == Basis::H20 ? 20 : 27));

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
          if (g.cnt < 8) g.pos[g.cnt++] = pos;
        }

        // 2) coarsen candidates
        std::size_t coarsened = 0;

        for (u32 pidx : touched) {
          Group& g = per[pidx];
          if (g.cnt != 8) continue;

          // level/minDepth check
          const u32 kid0 = _leaves[g.pos[0]];
          const u32 lev  = _tree_nodes[kid0].level;
          if (lev <= _minDepth && !_allowCoarsenBelowMinDepth) continue;

          // probe parent cell
          TreeNode3D& P = _tree_nodes[pidx];

          double x0, y0, z0, x1, y1, z1;
          leaf_bounds_3d(P, x0, y0, z0, x1, y1, z1);
          const double xm = 0.5 * (x0 + x1);
          const double ym = 0.5 * (y0 + y1);
          const double zm = 0.5 * (z0 + z1);

          // fill parent-space points into stack buffer
          int n = 0;
          // corners (0..7)
          s_buf[n++] = {x0, y0, z0}; s_buf[n++] = {x1, y0, z0};
          s_buf[n++] = {x1, y1, z0}; s_buf[n++] = {x0, y1, z0};
          s_buf[n++] = {x0, y0, z1}; s_buf[n++] = {x1, y0, z1};
          s_buf[n++] = {x1, y1, z1}; s_buf[n++] = {x0, y1, z1};
          if (nprobe >= 20) {
            // z = -1 edges (8..11)
            s_buf[n++] = {xm, y0, z0};
            s_buf[n++] = {x1, ym, z0};
            s_buf[n++] = {xm, y1, z0};
            s_buf[n++] = {x0, ym, z0};
            // z = +1 edges (12..15)
            s_buf[n++] = {xm, y0, z1};
            s_buf[n++] = {x1, ym, z1};
            s_buf[n++] = {xm, y1, z1};
            s_buf[n++] = {x0, ym, z1};
            // vertical edges (16..19)
            s_buf[n++] = {x0, y0, zm};
            s_buf[n++] = {x1, y0, zm};
            s_buf[n++] = {x1, y1, zm};
            s_buf[n++] = {x0, y1, zm};
          }
          if (nprobe == 27) {
            // face centers (20..25)
            s_buf[n++] = {xm, y0, zm};
            s_buf[n++] = {x1, ym, zm};
            s_buf[n++] = {xm, y1, zm};
            s_buf[n++] = {x0, ym, zm};
            s_buf[n++] = {xm, ym, z0};
            s_buf[n++] = {xm, ym, z1};
            // cell center (26)
            s_buf[n++] = {xm, ym, zm};
          }

          // map to physical
          for (int i = 0; i < n; ++i) xyz_buf[i] = parent_to_physical(s_buf[i]);

          // wrap with pre-sized vectors (no alloc)
          pts_s_vec.resize(n);
          pts_xyz_vec.resize(n);
          for (int i = 0; i < n; ++i) {
            pts_s_vec[i]   = s_buf[i];
            pts_xyz_vec[i] = xyz_buf[i];
          }

          if (!should_coarsen(P.morton, lev, pts_s_vec, pts_xyz_vec, /*dummy*/{}))
            continue;

          // mark 8 kids by leaf position; don’t mutate _leaves yet
          for (int k = 0; k < 8; ++k) {
            const u32 pos = g.pos[k];
            _tree_nodes[_leaves[pos]].is_leaf = false;
            mark_pos(pos);
          }

          // parent becomes leaf
          P.is_leaf = true;
          for (int c = 0; c < 8; ++c) P.child[c] = npos32;
          to_add.push_back(pidx);
          ++coarsened;
        }

        // 3) rebuild leaf set (single pass)
        if (coarsened) {
          for (u32 i = 0; i < L; ++i) if (!is_marked(i)) newLeaves.push_back(_leaves[i]);
          newLeaves.insert(newLeaves.end(), to_add.begin(), to_add.end());
          _leaves.swap(newLeaves);
          // (optional) defer _node2leafpos rebuild; post_topology_update() will handle it.
        }

        return coarsened;
      }

















      // refine a leaf into 8 children (returns true on actual split)
      bool refine_leaf_once(u32 leaf_node_idx) {
        if (leaf_node_idx == npos32) return false;
        const TreeNode3D parent_copy = _tree_nodes[leaf_node_idx];
        if (!parent_copy.is_leaf || parent_copy.level >= _maxDepth) return false;

        // remove old leaf from list
        auto it = std::find(_leaves.begin(), _leaves.end(), leaf_node_idx);
        if (it != _leaves.end()) {
          *it = _leaves.back();
          _leaves.pop_back();
        }

        _tree_nodes[leaf_node_idx].is_leaf = false;

        for (int dz = 0; dz < 2; ++dz)
          for (int dy = 0; dy < 2; ++dy)
            for (int dx = 0; dx < 2; ++dx) {
              const int q = (dz << 2) | (dy << 1) | dx;
              TreeNode3D child{};
              child.ix = (parent_copy.ix << 1) + (u32)dx;
              child.iy = (parent_copy.iy << 1) + (u32)dy;
              child.iz = (parent_copy.iz << 1) + (u32)dz;
              child.level  = parent_copy.level + 1;
              child.morton = interleave3(child.ix, child.iy, child.iz);
              child.parent = leaf_node_idx;
              const u32 cidx = (u32)_tree_nodes.size();
              _tree_nodes.push_back(child);
              _tree_nodes[leaf_node_idx].child[q] = cidx;
              _leaves.push_back(cidx);
            }
        return true;
      }

      // neighbor by face (0=-x,1=+x,2=-y,3=+y,4=-z,5=+z), returns node index or npos32
      u32 neighbor_leaf_by_face_any(u32 leaf_node_idx, int dir) const {
        const auto& me = _tree_nodes[leaf_node_idx];
        const u32 bx = me.ix & 1u, by = me.iy & 1u, bz = me.iz & 1u;
        auto try_local = [&](u32 & sib_code)->bool {
          const u32 code = (bz << 2) | (by << 1) | bx;
          sib_code = code;
          switch (dir) {
          case 0: if (bx == 1) {
              sib_code = (code & ~1u);
              return true;
            } break;
          case 1: if (bx == 0) {
              sib_code = (code & ~1u) | 1u;
              return true;
            } break;
          case 2: if (by == 1) {
              sib_code = (code & ~2u);
              return true;
            } break;
          case 3: if (by == 0) {
              sib_code = (code & ~2u) | 2u;
              return true;
            } break;
          case 4: if (bz == 1) {
              sib_code = (code & ~4u);
              return true;
            } break;
          case 5: if (bz == 0) {
              sib_code = (code & ~4u) | 4u;
              return true;
            } break;
          }
          return false;
        };

        auto descend_to_face_mid = [this, &me, dir](u32 subtree_root)->u32 {
          // compute the shared-face midpoint in parent space
          double x0, y0, z0, x1, y1, z1;
          leaf_bounds_3d(me, x0, y0, z0, x1, y1, z1);

          Point3 q;
          switch (dir) {
          case 0: q = {x0, 0.5 * (y0 + y1), 0.5 * (z0 + z1)}; break; // -x
          case 1: q = {x1, 0.5 * (y0 + y1), 0.5 * (z0 + z1)}; break; // +x
          case 2: q = {0.5 * (x0 + x1), y0, 0.5 * (z0 + z1)}; break; // -y
          case 3: q = {0.5 * (x0 + x1), y1, 0.5 * (z0 + z1)}; break; // +y
          case 4: q = {0.5 * (x0 + x1), 0.5 * (y0 + y1), z0}; break; // -z
          default: q = {0.5 * (x0 + x1), 0.5 * (y0 + y1), z1}; break; // +z
          }

          u32 node = subtree_root;
          while (true) {
            const auto& n = _tree_nodes[node];
            if (n.is_leaf) return node;

            double bx0, by0, bz0, bx1, by1, bz1;
            leaf_bounds_3d(n, bx0, by0, bz0, bx1, by1, bz1);
            const double cx = 0.5 * (bx0 + bx1), cy = 0.5 * (by0 + by1), cz = 0.5 * (bz0 + bz1);

            const int qx = (q[0] >= cx) ? 1 : 0;
            const int qy = (q[1] >= cy) ? 1 : 0;
            const int qz = (q[2] >= cz) ? 1 : 0;
            const int code = (qz << 2) | (qy << 1) | qx;

            u32 ch = n.child[code];
            if (ch == npos32) return node;
            node = ch;
          }
        };


        // same-parent sibling?
        if (me.parent != npos32) {
          u32 sib_code;
          if (try_local(sib_code)) {
            u32 cand = _tree_nodes[me.parent].child[sib_code];
            if (cand != npos32) return descend_to_face_mid(cand);
          }
        }

        // climb
        u32 cur = leaf_node_idx;
        while (true) {
          const auto& c = _tree_nodes[cur];
          if (c.parent == npos32) return npos32;
          const auto& p = _tree_nodes[c.parent];
          const u32 cx = _tree_nodes[cur].ix & 1u, cy = _tree_nodes[cur].iy & 1u, cz = _tree_nodes[cur].iz & 1u;
          u32 code = (cz << 2) | (cy << 1) | cx, sib_code = code;
          bool can = false;
          switch (dir) {
          case 0: if (cx == 1) {
              sib_code = (code & ~1u);
              can = true;
            } break;
          case 1: if (cx == 0) {
              sib_code = (code & ~1u) | 1u;
              can = true;
            } break;
          case 2: if (cy == 1) {
              sib_code = (code & ~2u);
              can = true;
            } break;
          case 3: if (cy == 0) {
              sib_code = (code & ~2u) | 2u;
              can = true;
            } break;
          case 4: if (cz == 1) {
              sib_code = (code & ~4u);
              can = true;
            } break;
          case 5: if (cz == 0) {
              sib_code = (code & ~4u) | 4u;
              can = true;
            } break;
          }
          if (can) {
            u32 sib = p.child[sib_code];
            if (sib == npos32) return npos32;
            return descend_to_face_mid(sib);
          }
          cur = c.parent;
        }
      }

      const std::vector<u32>& leaves() const {
        return _leaves;
      }
      const TreeNode3D& node(u32 idx) const {
        return _tree_nodes[idx];
      }


      // Evaluate field directly in parent coordinates (xi, eta, zeta).
// Locates the leaf, maps to local [-1,1]^3, and interpolates.
      bool evaluate_field_on_parent(u32 fid, const Point3& s, double& value) const {
        if (fid >= _fields.size()) return false;

        u32    leaf_node_idx = npos32;
        Point3 shat;

        // 1) locate leaf and local reference coords in [-1,1]^3
        if (!locate_leaf_on_parent_and_ref(s, leaf_node_idx, shat)) {
          return false;
        }

        // 2) leaf node index -> position in coefficient storage
        const u32 leaf_pos = (leaf_node_idx < _node2leafpos.size()) ? _node2leafpos[leaf_node_idx] : npos32;
        if (leaf_pos == npos32) return false;

        // 3) access field + connectivity
        const Field&          f = _fields[fid];
        const BasisRegistry&  R = _basisReg[(int)f.basis];
        const auto&           conn = R.elem2glob[leaf_pos];

        // 4) interpolate by basis
        switch (f.basis) {
        case Basis::H8: {
          double N[8];
          Shapes3::H8(shat, N);
          value = 0.0;
          for (int a = 0; a < 8; ++a) value += N[a] * f.nodal[conn[a]];
        } break;

        case Basis::H20: {
          double N[20];
          Shapes3::H20(shat, N);
          value = 0.0;
          for (int a = 0; a < 20; ++a) value += N[a] * f.nodal[conn[a]];
        } break;

        case Basis::H27: {
          double N[27];
          Shapes3::H27(shat, N);
          value = 0.0;
          for (int a = 0; a < 27; ++a) value += N[a] * f.nodal[conn[a]];
        } break;
        }

        return true;
      }


      inline bool locate_leaf_on_parent_and_ref(const Point3& s,
                                                u32& leaf_node_idx,
                                                Point3& shat) const {
        // 1) Find the leaf containing (xi,eta,zeta) in parent space.
        leaf_node_idx = locate_leaf_on_parent(s);
        if (leaf_node_idx == npos32) return false;

        // 2) Convert to local [-1,1]^3 using level + indices (no divisions).
        local_ref_fast_3d(_tree_nodes[leaf_node_idx], s, shat);

        return true;
      }



// Refine tree so that each physical point lies in a leaf of at least 'maxDepthTarget'.
// Caches ONLY the inverse map (xi,eta,zeta) per point; uses _node2leafpos instead of per-pass node_to_pos.
      void refine_to_contain_points(const std::vector<Point3>& pts,
                                    u32 maxDepthTarget) {
        if (pts.empty()) return;

        // ---- Cache inverse map (xi,eta,zeta) once for all points ----
        struct XiEtaZeta {
          Point3 xi;
          bool   valid;
        };
        std::vector<XiEtaZeta> pm(pts.size());
        for (size_t i = 0; i < pts.size(); ++i) {
          Point3 s;
          if (!inverse_map_hex27(pts[i], s)) {
            pm[i] = { {0.0, 0.0, 0.0}, false };
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

            // Reuse cached (xi,eta,zeta); do only the tree lookup per pass
            const u32 node = locate_leaf_on_parent(pm[i].xi);
            if (node == npos32) continue;

            const u32 lvl = _tree_nodes[node].level;
            if (lvl >= maxDepthTarget) continue;

            const u32 pos = _node2leafpos[node]; // fast, no per-pass map build
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

          // Net new nodes per split in 3D octree: +7 tree nodes (8 kids minus parent as leaf)
          if (to_split) {
            _tree_nodes.reserve(_tree_nodes.size() + to_split * 7u);
            _leaves.reserve(_leaves.size() + to_split * 7u); // parent replaced by 8 children → +7 leaves
          }

          // Assumes you have a 3D version that refines leaves by *position* mark:
          // - It should read marks via should_refine_pos(leaf_pos)
          // - It should split a leaf into 8 children
          // - It should update _leaves and _tree_nodes accordingly
          const std::size_t nsplit = refine_pass_min(should_refine_pos);
          if (nsplit == 0) break; // reached max depth locally or topology unchanged

          // _leaves changed inside refine_pass_min; _node2leafpos will be rebuilt
          // at the top of the next loop iteration.
        }

        // Keep mesh 1-irregular and refresh bookkeeping
        enforce_balance();
        post_topology_update();
      }

// Perform one refinement pass with a minimal predicate: bool(u32 leaf_pos)
// No geometry probing; fastest path. 3D (octree) version: splits into 8 children.
      template<class RefinePred>
      std::size_t refine_pass_min(RefinePred&& should_refine) {
        std::vector<u32> newLeaves;
        newLeaves.reserve(_leaves.size() * 2);
        std::size_t refined = 0;

        const u32 n0 = (u32)leaf_count(); // snapshot leaf count

        for (u32 idx = 0; idx < n0; ++idx) {
          const u32 leaf_node_idx = _leaves[idx];               // index into _tree_nodes
          const TreeNode3D leaf_copy = _tree_nodes[leaf_node_idx]; // snapshot

          if (leaf_copy.is_leaf && leaf_copy.level < _maxDepth && should_refine(idx)) {
            // split into 8 children
            const u32 ix = leaf_copy.ix;
            const u32 iy = leaf_copy.iy;
            const u32 iz = leaf_copy.iz;

            for (int dz = 0; dz < 2; ++dz)
              for (int dy = 0; dy < 2; ++dy)
                for (int dx = 0; dx < 2; ++dx) {
                  const int q = (dz << 2) | (dy << 1) | dx;

                  TreeNode3D child{};
                  child.ix = (ix << 1) + (u32)dx;
                  child.iy = (iy << 1) + (u32)dy;
                  child.iz = (iz << 1) + (u32)dz;
                  child.level  = leaf_copy.level + 1;
                  child.morton = interleave3(child.ix, child.iy, child.iz);
                  child.is_leaf = true;
                  child.parent  = leaf_node_idx;
                  for (int c = 0; c < 8; ++c) child.child[c] = npos32;

                  const u32 child_idx = (u32)_tree_nodes.size();
                  _tree_nodes.push_back(child);
                  _tree_nodes[leaf_node_idx].child[q] = child_idx;
                  newLeaves.push_back(child_idx);
                }

            _tree_nodes[leaf_node_idx].is_leaf = false;
            refined++;
          }
          else {
            newLeaves.push_back(leaf_node_idx);
          }
        }

        _leaves.swap(newLeaves);
        return refined;
      }




// Simple ring (really: advancing head) FIFO: contiguous + cache-friendly.
      struct RingQ {
        std::vector<u32> buf;
        size_t head = 0; // index of current front
        inline void clear() {
          buf.clear();
          head = 0;
        }
        inline bool empty() const {
          return head >= buf.size();
        }
        inline void push_back(u32 v) {
          buf.push_back(v);
        }
        inline u32 front() const {
          return buf[head];
        }
        inline void pop_front() {
          ++head;
        }
      };

// Enforce 1-irregularity in 3D: adjacent leaves may differ by at most 1 level.
      void enforce_balance() {
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
        q.buf.reserve(_leaves.size() * 6);           // heuristic
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

          for (int dir = 0; dir < 6; ++dir) { // 0=−x,1=+x,2=−y,3=+y,4=−z,5=+z
            u32 nb = neighbor_leaf_by_face_any(leaf, dir);
            if (nb == npos32 || nb >= _tree_nodes.size()) continue;

            const auto& nbn = _tree_nodes[nb];
            if (!nbn.is_leaf) continue;

            const u32 lnb = nbn.level;

            // 2:1 rule: |lev - lnb| <= 1
            if (lev > lnb + 1) {
              // neighbor too coarse -> refine neighbor
              if (refine_leaf_once(nb)) {
                // _tree_nodes/_leaves may have grown; keep inq in range
                if (inq.size() < _tree_nodes.size()) inq.resize(_tree_nodes.size(), 0);
                // enqueue neighbor's neighbors (only ones affected)
                for (int d = 0; d < 6; ++d) {
                  u32 nb2 = neighbor_leaf_by_face_any(nb, d);
                  if (nb2 != npos32) enqueue(q, inq, nb2);
                }
                // also re-check the current leaf again (it may still violate with new children)
                enqueue(q, inq, leaf);
              }
            }
            else if (lnb > lev + 1) {
              // this leaf too coarse -> refine this leaf
              if (refine_leaf_once(leaf)) {
                if (inq.size() < _tree_nodes.size()) inq.resize(_tree_nodes.size(), 0);
                // enqueue this (former) leaf's neighbors
                for (int d = 0; d < 6; ++d) {
                  u32 nb2 = neighbor_leaf_by_face_any(leaf, d);
                  if (nb2 != npos32) enqueue(q, inq, nb2);
                }
                // no point checking other dirs for the old leaf; it’s not a leaf now
                break;
              }
            }
          }
        }
      }


      // Copia "veloce": riusa la memoria già allocata, cresce solo se serve.
// NON chiama post_topology_update() (stai copiando strutture già consistenti).
      void assign_from(const HexTree3D& rhs) {
        if (this == &rhs) return;

        // ---- helper locali per copiare riusando capacità ----
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
          dst.clear(); // mantiene i bucket esistenti
          if (dst.bucket_count() < src.size()) dst.rehash(src.size());
          for (const auto& kv : src) dst.insert(kv);
        };
        auto copy_uset = [](auto & dst, const auto & src) {
          dst.clear();
          if (dst.bucket_count() < src.size()) dst.rehash(src.size());
          for (const auto& v : src) dst.insert(v);
        };

        // ---- scalari / config ----
        _maxDepth  = rhs._maxDepth;
        _minDepth  = rhs._minDepth;
        _allowCoarsenBelowMinDepth = rhs._allowCoarsenBelowMinDepth;
        _leafpos_valid = rhs._leafpos_valid;

        _root       = rhs._root;
        _geom_ready = rhs._geom_ready;

        // ---- geometria (H27) ----
        _X = rhs._X;
        _Y = rhs._Y;
        _Z = rhs._Z;

        // ---- basi attive e campi ----
        copy_uset(_activeBases, rhs._activeBases);
        copy_vec(_fields,       rhs._fields);

        // ---- topologia ----
        copy_vec(_tree_nodes,    rhs._tree_nodes);   // std::vector<TreeNode3D>
        copy_vec(_leaves,        rhs._leaves);       // vettore di indici-nodo foglia
        copy_vec(_leaf_ids,      rhs._leaf_ids);     // cache degli id foglia [0..n-1]
        copy_vec(_node2leafpos,  rhs._node2leafpos); // mappa nodo -> pos foglia

        // ---- registri per-basis ----
        for (int b = 0; b < 3; ++b) {
          copy_umap(_basisReg[b].nodeMap, rhs._basisReg[b].nodeMap);     // u64 -> int
          copy_vec(_basisReg[b].nodes,    rhs._basisReg[b].nodes);       // std::vector<FEMNode3D>
          copy_vecvec_int(_basisReg[b].elem2glob, rhs._basisReg[b].elem2glob); // connettività
        }
      }


      // Conservative coarsen cycle using snapshot + parent coords; rebuild all fields (3D)
      std::size_t coarsen_only_cycle_safe(u32 fid,
                                          double tau_coarse,
                                          HexTree3D& snapshot,
                                          u32 max_passes = 10,
                                          Basis probe_basis = Basis::H27) {
        // Freeze current state for conservative evaluation/transfer
        snapshot.assign_from(*this);

        // Coarsen predicate (evaluated on the snapshot, in parent coords)
        auto pred = [&](u64 /*parent_morton*/, u32 level,
                        const std::vector<Point3>& pts_s,     // parent (xi,eta,zeta)
                        const std::vector<Point3>& /*pts_xyz*/,
        const std::vector<std::array<double, 27>>& /*Nvals*/) -> bool {
          if (level <= min_depth()) return false;
          if (pts_s.empty()) return false;

          double v0;
          if (!snapshot.evaluate_field_on_parent(fid, pts_s[0], v0))
            return false;

          double mn = v0, mx = v0;
          for (size_t i = 1; i < pts_s.size(); ++i) {
            double val;
            if (snapshot.evaluate_field_on_parent(fid, pts_s[i], val)) {
              mn = std::min(mn, val);
              mx = std::max(mx, val);
            }
          }
          return (mn > +tau_coarse) || (mx < -tau_coarse);
        };

        std::size_t total = 0;

        // Coarsen passes; (balance+topology update happen after the loop)
        for (u32 pass = 0; pass < max_passes; ++pass) {
          std::size_t c = coarsen_pass(pred, probe_basis);
          if (c == 0) break;
          total += c;
        }

        // Enforce 1-irregularity and refresh bookkeeping
        enforce_balance();
        post_topology_update();

        // Rebuild all fields from the snapshot (conservative transfer) on nodal sets
        for (u32 f = 0; f < _fields.size(); ++f) {
          rebuild_field_from(snapshot, f);
        }

        return total;
      }


      // Rebuild field fid on *this* from source tree 'src' by sampling at parent nodes (global nodal storage) — 3D
      void rebuild_field_from(const HexTree3D& src, u32 fid) {
        assert(fid < _fields.size() && fid < src._fields.size());

        Field&       dst = _fields[fid];
        const Field& s   = src._fields[fid];
        assert(dst.basis == s.basis);

        const BasisRegistry& Rdst = _basisReg[(int)dst.basis];
        const BasisRegistry& Rsrc = src._basisReg[(int)s.basis];

        dst.nodal.resize(Rdst.nodes.size());

        // Quantization helper that matches src.get_or_insert_gid(...)
        const u32 nodeBits_src = src._maxDepth + 1;                  // include midpoints at last level
        assert(nodeBits_src <= 21 && "node index packing requires (_maxDepth + 1) <= 21");
        const u32 nodesN_src   = (1u << nodeBits_src);

        auto to_idx_src = [nodesN_src](double ss) -> u32 {
          if (ss <= -1.0) return 0u;
          if (ss >=  1.0) return nodesN_src - 1u;                    // keep within range, match src
          double t = (ss + 1.0) * double(nodesN_src) * 0.5;          // [-1,1] -> [0, nodesN)
          long long li = llround(t);
          if (li < 0) li = 0;
          if (li > (long long)(nodesN_src - 1)) li = (long long)(nodesN_src - 1);
          return (u32)li;
        };

        for (size_t gid = 0; gid < Rdst.nodes.size(); ++gid) {
          const auto& pr = Rdst.nodes[gid].parent;                   // (xi,eta,zeta) at destination
          const u32 ix = to_idx_src(pr[0]);
          const u32 iy = to_idx_src(pr[1]);
          const u32 iz = to_idx_src(pr[2]);
          const u64 key = interleave3(ix, iy, iz);                   // same 3-way Morton packing as src

          auto it = Rsrc.nodeMap.find(key);
          if (it != Rsrc.nodeMap.end()) {
            // Fast path: direct copy from src nodal storage
            dst.nodal[gid] = s.nodal[it->second];
          }
          else {
            // Fallback: sample from src at parent coords (rare if depths/registries differ)
            double val = std::numeric_limits<double>::quiet_NaN();
            const bool ok = src.evaluate_field_on_parent(fid, Point3{pr[0], pr[1], pr[2]}, val);
            dst.nodal[gid] = ok ? val : std::numeric_limits<double>::quiet_NaN();
          }
        }
      }


// Returns the element-to-global connectivity (per-basis). (3D)
      const std::vector<std::vector<int>>& basis_connectivity(Basis b) const {
        return _basisReg[(int)b].elem2glob;
      }

// Returns the unique global FEM nodes (per-basis), with parent & physical coords.
      const std::vector<FEMNode3D>& basis_nodes(Basis b) const {
        return _basisReg[(int)b].nodes;
      }

// Print field on a vtu binary mesh to visualize on Paraview
      bool write_binary_vtu(const std::string &filename, u32 fid, const std::string &name,
                            bool cell_centered = false) const {
        require_geometry();
        if (fid >= _fields.size()) return false;

        const Field &fld = _fields[fid];
        Basis b = fld.basis;
        const BasisRegistry &R = _basisReg[(int)b];

        const size_t numCells = R.elem2glob.size();
        if (numCells == 0) return false;

        // -------- Points (global) --------
        std::vector<double> flatPoints;
        flatPoints.reserve(R.nodes.size() * 3);
        for (const auto &n : R.nodes) {
          flatPoints.push_back(n.physical[0]);
          flatPoints.push_back(n.physical[1]);
          flatPoints.push_back(n.physical[2]);
        }

        // -------- Connectivity --------
        std::vector<int> connectivity;
        std::vector<int> offsets;
        std::vector<unsigned char> types;

        const int perCell = (b == Basis::H8 ? 8 : (b == Basis::H20 ? 20 : 27));
        connectivity.reserve(numCells * perCell);
        offsets.reserve(numCells);
        types.reserve(numCells);

        auto vtk_cell_type = [](Basis bb) -> unsigned char {
          // VTK: HEX=12, QUADRATIC_HEX=25, TRIQUADRATIC_HEX=29
          switch (bb) {
          case Basis::H8:  return 12;
          case Basis::H20: return 25;
          case Basis::H27: return 29;
          }
          return 12;
        };

        // H27 tail remap (indices 20..26 in our internal ordering → VTK ordering)
        // Confirmed correct: {20+3, 20+1, 20+0, 20+2, 20+4, 20+5, 20+6}
        const std::array<int, 7> tail_map{{3, 1, 0, 2, 4, 5, 6}};

        int offset = 0;
        for (size_t e = 0; e < numCells; ++e) {
          const auto &conn = R.elem2glob[e];

          if (b == Basis::H27) {
            // Expect 27 nodes. First 20 match VTK; remap the last 7 via tail_map.
            assert(conn.size() == 27 && "H27 element must have 27 nodes");
            connectivity.insert(connectivity.end(), conn.begin(), conn.begin() + 20);
            for (int i = 0; i < 7; ++i) {
              connectivity.push_back(conn[20 + tail_map[i]]);
            }
            offset += 27;
          }
          else {
            connectivity.insert(connectivity.end(), conn.begin(), conn.end());
            offset += (int)conn.size();
          }

          offsets.push_back(offset);
          types.push_back(vtk_cell_type(b));
        }

        // -------- Data --------
        std::vector<double> pointData, cellData;
        if (cell_centered) {
          cellData.reserve(numCells);
          for (size_t e = 0; e < numCells; ++e) {
            const auto &conn = R.elem2glob[e];
            double v = 0.0; for (int gid : conn) v += fld.nodal[gid];
            cellData.push_back(v / (double)conn.size());
          }
        }
        else {
          pointData = fld.nodal; // unchanged: points array order == registry order
        }

        // -------- XML --------
        std::ofstream os(filename);
        if (!os) return false;

        os << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        os << "  <UnstructuredGrid>\n";
        os << "    <Piece NumberOfPoints=\"" << R.nodes.size()
           << "\" NumberOfCells=\"" << numCells << "\">\n";

        os << "      <Points>\n";
        write_binary_array(os, "Float64", "", 3, flatPoints);
        os << "      </Points>\n";

        os << "      <Cells>\n";
        write_binary_array(os, "Int32", "connectivity", 1, connectivity);
        write_binary_array(os, "Int32", "offsets", 1, offsets);
        write_binary_array(os, "UInt8", "types", 1, types);
        os << "      </Cells>\n";

        if (cell_centered) {
          os << "      <CellData Scalars=\"" << name << "\">\n";
          write_binary_array(os, "Float64", name, 1, cellData);
          os << "      </CellData>\n";
        }
        else {
          os << "      <PointData Scalars=\"" << name << "\">\n";
          write_binary_array(os, "Float64", name, 1, pointData);
          os << "      </PointData>\n";
        }

        os << "    </Piece>\n";
        os << "  </UnstructuredGrid>\n";
        os << "</VTKFile>\n";
        return true;
      }

// ---- global maximum depth allowed (unchanged) ----
      u32 max_depth() const {
        return _maxDepth;
      }

      // ---- coarse geometry accessors (3D) ----
      inline const std::array<double, 27>& get_X() const {
        return _X;
      }
      inline const std::array<double, 27>& get_Y() const {
        return _Y;
      }
      inline const std::array<double, 27>& get_Z() const {
        return _Z;
      }

      // ---- collect reference (parent) coordinates across a level range (3D) ----
      std::vector<Point3>
      extract_node_parent_coords_in_level_range(u32 lev_min, u32 lev_max, Basis basis) const {
        std::vector<Point3> coords;

        const auto& ids = leaf_indices();          // assumes you have the 3D version returning [0..n-1]
        for (u32 leaf_pos : ids) {
          const u32 lev = level_of(leaf_pos);      // leaf position -> level
          if (lev < lev_min || lev > lev_max) continue;

          std::vector<Point3> s;
          extract_leaf_parent_coords(basis, leaf_pos, s);  // fills parent-space nodes for this leaf
          coords.insert(coords.end(), s.begin(), s.end());
        }
        return coords;
      }

    private:

      // bounds of a leaf in parent space [-1,1]^3
      inline void leaf_bounds_3d(const TreeNode3D& n,
                                 double& x0, double& y0, double& z0,
                                 double& x1, double& y1, double& z1) const {
        const double N = double(1u << n.level);
        const double d = 2.0 / N;
        x0 = -1.0 + n.ix * d; x1 = x0 + d;
        y0 = -1.0 + n.iy * d; y1 = y0 + d;
        z0 = -1.0 + n.iz * d; z1 = z0 + d;
      }

      // fast local reference map to [-1,1]^3 for a given leaf node
      inline void local_ref_fast_3d(const TreeNode3D& n,
                                    const Point3& s,
                                    Point3& shat) const {
        const double N = double(1u << n.level);
        shat[0] = N * s[0] + N - double((n.ix << 1) + 1);
        shat[1] = N * s[1] + N - double((n.iy << 1) + 1);
        shat[2] = N * s[2] + N - double((n.iz << 1) + 1);
      }

      struct BasisRegistry {
        std::unordered_map<u64, int> nodeMap;            // key -> gid
        std::vector<FEMNode3D> nodes;                      // gid -> FEM node (parent + physical)
        std::vector<std::vector<int>> elem2glob;         // per-element connectivity
        void clear() {
          nodeMap.clear();
          nodes.clear();
          elem2glob.clear();
        }
      };


      // Ensure all leaves are refined down to at least _minDepth (3D, 8-way splits)
      void ensure_min_depth() {
        if (_minDepth == 0) return;

        bool changed = true;
        while (changed) {
          changed = false;
          std::vector<u32> newLeaves;
          newLeaves.reserve(_leaves.size() * 2u); // heuristic

          for (u32 leaf_idx : _leaves) {
            const TreeNode3D node_copy = _tree_nodes[leaf_idx];

            if (node_copy.level < _minDepth) {
              // refine this leaf into 8 children
              const u32 ix = node_copy.ix;
              const u32 iy = node_copy.iy;
              const u32 iz = node_copy.iz;

              _tree_nodes[leaf_idx].is_leaf = false;

              for (int dz = 0; dz < 2; ++dz) {
                for (int dy = 0; dy < 2; ++dy) {
                  for (int dx = 0; dx < 2; ++dx) {
                    const int q = (dz << 2) | (dy << 1) | dx;
                    const u32 cindex = (u32)_tree_nodes.size();

                    TreeNode3D child{};
                    child.ix = (ix << 1) + (u32)dx;
                    child.iy = (iy << 1) + (u32)dy;
                    child.iz = (iz << 1) + (u32)dz;
                    child.morton  = interleave3(child.ix, child.iy, child.iz);
                    child.level   = node_copy.level + 1;
                    child.is_leaf = true;
                    child.parent  = leaf_idx;
                    for (int c = 0; c < 8; ++c) child.child[c] = npos32;

                    _tree_nodes.push_back(child);
                    _tree_nodes[leaf_idx].child[q] = cindex;
                    newLeaves.push_back(cindex);
                  }
                }
              }

              changed = true;
            }
            else {
              newLeaves.push_back(leaf_idx);
            }
          }

          if (changed) {
            _leaves.swap(newLeaves);
          }
        }

        post_topology_update();
      }


// ---- return refinement level of a given leaf *position* (3D) ----
      u32 level_of(u32 leaf_pos) const {
        const u32 node_idx = _leaves[leaf_pos];    // leaf position -> node index
        return _tree_nodes[node_idx].level;
      }


      int get_or_insert_gid(BasisRegistry &R, double xi, double eta, double zeta) {
        // Need +1 bit so level=_maxDepth mid/center nodes (half-steps) map to integers.
        const u32 nodeBits = std::max<u32>(2, std::min<u32>(_maxDepth + 1, 21));


        // Our 3-way Morton interleaver packs 21 bits per axis.
        assert(nodeBits <= 21 && "node index packing requires (_maxDepth + 1) <= 21");
        const u32 nodesN = (1u << nodeBits); // indices in [0..nodesN-1]

        auto to_idx = [nodesN](double s)->u32 {
          if (s <= -1.0) return 0u;
          if (s >=  1.0) return nodesN - 1u;        // clamp top to fit 21-bit interleaver
          double t = (s + 1.0) * double(nodesN) * 0.5; // [-1,1] -> [0,nodesN]
          long long li = llround(t);
          if (li < 0) li = 0;
          if (li > (long long)(nodesN - 1)) li = (long long)(nodesN - 1);
          return (u32)li;
        };



        const u32 ix = to_idx(xi);
        const u32 iy = to_idx(eta);
        const u32 iz = to_idx(zeta);

        // 3-way Morton key (21-bit per axis)
        const u64 key = interleave3(ix, iy, iz);

        auto it = R.nodeMap.find(key);
        if (it != R.nodeMap.end()) return it->second;

        const int gid = (int)R.nodes.size();
        const Point3 parent = {xi, eta, zeta};
        R.nodes.push_back(FEMNode3D{gid, parent, parent_to_physical(parent)});
        R.nodeMap.emplace(key, gid);
        return gid;
      }

      struct BasisHasher {
        size_t operator()(fem::Basis b) const noexcept {
          return std::hash<uint8_t>()(static_cast<uint8_t>(b));
        }
      };

    private: //data
// config
      u32  _maxDepth;
      u32  _minDepth{0};
      bool _allowCoarsenBelowMinDepth{true};
      bool _leafpos_valid{false};

// geometry (H27)
      std::array<double, 27> _X {}, _Y {}, _Z{};
      bool _geom_ready{false};

// mesh leaves
      std::vector<u32> _leaves;
      mutable std::vector<u32> _leaf_ids;

// fields (per-leaf coefficients for API compatibility)
      std::vector<Field> _fields;

// per-basis global node registries + connectivity
      BasisRegistry _basisReg[3]; // 0:H8, 1:H20, 2:H27

      std::vector<TreeNode3D> _tree_nodes;        // full tree hierarchy
      u32 _root{0};

      mutable std::vector<u32> _node2leafpos; // size == _tree_nodes.size(), npos32 default
      std::unordered_set<Basis, BasisHasher> _activeBases;

  };

} // namespace fem




