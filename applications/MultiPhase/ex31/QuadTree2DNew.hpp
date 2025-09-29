#pragma once
#include <vector>
#include <array>
#include <cstdint>
#include <limits>
#include <cmath>
#include <cassert>
#include <functional>
#include <fstream>
#include <iomanip>
#include <algorithm>

/*
  QuadTree2D.hpp (ASCII VTU)
  - 2D quadtree living on parent domain [-1,1]^2
  - Quad9 (Q2) geometry map for the coarse/root element
  - Supports Q1, Serendipity8, and Q2 fields per leaf
  - Refinement/coarsening, field evaluation, and snapshot-safe coarsening
  - Parent-space helpers to avoid repeated inverse maps
*/

namespace fem {

  using u32 = uint32_t;
  using u64 = uint64_t;
  constexpr u32 npos32 = std::numeric_limits<u32>::max();

  // Supported element bases
  enum class Basis : uint8_t { Q1_Quad4 = 0, Serendipity8 = 1, Q2_Quad9 = 2 };

  //----------------------------------------
  // Quad9 shapes (parent space [-1,1]^2)
  //----------------------------------------
  struct Quad9Shape {
    // 1D Q2 Lagrange at s in {-1,0,1}: returns basis and derivatives
    static inline void q2_1d(double s, double& L0, double& L1, double& L2,
                             double& dL0, double& dL1, double& dL2) {
      L0 = 0.5 * s * (s - 1.0);   // 0.5*(s^2 - s)
      L1 = 1.0 - s * s;           // 1 - s^2
      L2 = 0.5 * s * (s + 1.0);   // 0.5*(s^2 + s)
      dL0 = s - 0.5;
      dL1 = -2.0 * s;
      dL2 = s + 0.5;
    }

    // Compute Q2 shape functions at (xi,eta) in parent space
    static inline void N(double xi, double eta, double N9[9]) {
      double Lx[3], Ly[3], dx[3], dy[3];
      q2_1d(xi,  Lx[0], Lx[1], Lx[2], dx[0], dx[1], dx[2]);
      q2_1d(eta, Ly[0], Ly[1], Ly[2], dy[0], dy[1], dy[2]);
      (void)dx;
      (void)dy;

      // corners
      N9[0] = Lx[0] * Ly[0]; // (-1,-1)
      N9[1] = Lx[2] * Ly[0]; // ( 1,-1)
      N9[2] = Lx[2] * Ly[2]; // ( 1, 1)
      N9[3] = Lx[0] * Ly[2]; // (-1, 1)
      // mids
      N9[4] = Lx[1] * Ly[0]; // (0,-1)
      N9[5] = Lx[2] * Ly[1]; // (1,0)
      N9[6] = Lx[1] * Ly[2]; // (0,1)
      N9[7] = Lx[0] * Ly[1]; // (-1,0)
      // center
      N9[8] = Lx[1] * Ly[1]; // (0,0)
    }

    // Compute derivative wrt xi and eta for Q2 at (xi,eta)
    static inline void dN(double xi, double eta, double dN_dxi[9], double dN_deta[9]) {
      double Lx[3], Ly[3], dx[3], dy[3];
      q2_1d(xi,  Lx[0], Lx[1], Lx[2], dx[0], dx[1], dx[2]);
      q2_1d(eta, Ly[0], Ly[1], Ly[2], dy[0], dy[1], dy[2]);

      dN_dxi[0]  = dx[0] * Ly[0];
      dN_deta[0] = Lx[0] * dy[0];
      dN_dxi[1]  = dx[2] * Ly[0];
      dN_deta[1] = Lx[2] * dy[0];
      dN_dxi[2]  = dx[2] * Ly[2];
      dN_deta[2] = Lx[2] * dy[2];
      dN_dxi[3]  = dx[0] * Ly[2];
      dN_deta[3] = Lx[0] * dy[2];
      dN_dxi[4]  = dx[1] * Ly[0];
      dN_deta[4] = Lx[1] * dy[0];
      dN_dxi[5]  = dx[2] * Ly[1];
      dN_deta[5] = Lx[2] * dy[1];
      dN_dxi[6]  = dx[1] * Ly[2];
      dN_deta[6] = Lx[1] * dy[2];
      dN_dxi[7]  = dx[0] * Ly[1];
      dN_deta[7] = Lx[0] * dy[1];
      dN_dxi[8]  = dx[1] * Ly[1];
      dN_deta[8] = Lx[1] * dy[1];
    }
  };

  //----------------------------------------
  // Q1 (Quad4) and Serendipity8 shapes on [-1,1]^2
  //----------------------------------------
  struct Shapes {
    // Compute Q1 bilinear shape functions at (xi,eta)
    static inline void Q1(double xi, double eta, double N4[4]) {
      const double a = 0.25 * (1 - xi), b = 0.25 * (1 + xi);
      N4[0] = a * (1 - eta); // (-1,-1)
      N4[1] = b * (1 - eta); // ( 1,-1)
      N4[2] = b * (1 + eta); // ( 1, 1)
      N4[3] = a * (1 + eta); // (-1, 1)
    }

    // Compute serendipity-8 shape functions at (xi,eta)
    static inline void Serendipity8(double xi, double eta, double N8[8]) {
      const double xm = xi, em = eta;
      N8[0] = 0.25 * (1 - xm) * (1 - em) * (-xm - em - 1);
      N8[1] = 0.25 * (1 + xm) * (1 - em) * (xm - em - 1);
      N8[2] = 0.25 * (1 + xm) * (1 + em) * (xm + em - 1);
      N8[3] = 0.25 * (1 - xm) * (1 + em) * (-xm + em - 1);
      N8[4] = 0.5 * (1 - xm * xm) * (1 - em);
      N8[5] = 0.5 * (1 + xm) * (1 - em * em);
      N8[6] = 0.5 * (1 - xm * xm) * (1 + em);
      N8[7] = 0.5 * (1 - xm) * (1 - em * em);
    }
  };

  //----------------------------------------
  // Morton helpers (2D)
  //----------------------------------------
  // Interleave bits of (x,y) into a 64-bit Morton code
  static inline u64 interleave2(u32 x, u32 y) {
    u64 xx = x, yy = y;
    xx = (xx | (xx << 16)) & 0x0000FFFF0000FFFFULL;
    xx = (xx | (xx << 8)) & 0x00FF00FF00FF00FFULL;
    xx = (xx | (xx << 4)) & 0x0F0F0F0F0F0F0F0FULL;
    xx = (xx | (xx << 2)) & 0x3333333333333333ULL;
    xx = (xx | (xx << 1)) & 0x5555555555555555ULL;

    yy = (yy | (yy << 16)) & 0x0000FFFF0000FFFFULL;
    yy = (yy | (yy << 8)) & 0x00FF00FF00FF00FFULL;
    yy = (yy | (yy << 4)) & 0x0F0F0F0F0F0F0F0FULL;
    yy = (yy | (yy << 2)) & 0x3333333333333333ULL;
    yy = (yy | (yy << 1)) & 0x5555555555555555ULL;

    return (yy << 1) | xx;
  }

  // Reverse Morton code to (x,y)
  static inline void deinterleave2(u64 m, u32& x, u32& y) {
    u64 xx = m, yy = m >> 1;
    xx &= 0x5555555555555555ULL;
    yy &= 0x5555555555555555ULL;

    xx = (xx | (xx >> 1)) & 0x3333333333333333ULL;
    xx = (xx | (xx >> 2)) & 0x0F0F0F0F0F0F0F0FULL;
    xx = (xx | (xx >> 4)) & 0x00FF00FF00FF00FFULL;
    xx = (xx | (xx >> 8)) & 0x0000FFFF0000FFFFULL;
    xx = (xx | (xx >> 16)) & 0x00000000FFFFFFFFULL;

    yy = (yy | (yy >> 1)) & 0x3333333333333333ULL;
    yy = (yy | (yy >> 2)) & 0x0F0F0F0F0F0F0F0FULL;
    yy = (yy | (yy >> 4)) & 0x00FF00FF00FF00FFULL;
    yy = (yy | (yy >> 8)) & 0x0000FFFF0000FFFFULL;
    yy = (yy | (yy >> 16)) & 0x00000000FFFFFFFFULL;

    x = (u32)xx;
    y = (u32)yy;
  }

  //----------------------------------------
  // Quadtree node
  //----------------------------------------
  struct QuadNode {
    u32 parent{npos32};
    std::array<u32, 4> child{npos32, npos32, npos32, npos32};
    u32 level{0};
    u64 morton{0};
    uint8_t flags{0};

    // True if node has no children
    inline bool is_leaf() const {
      return child[0] == npos32;
    }
  };

  //----------------------------------------
  // Field storage (per-leaf coefficients)
  //----------------------------------------
  struct Field {
    Basis basis{Basis::Q1_Quad4};
    u32   dofs_per_cell{4};
    std::vector<double> coeffs;
  };

  //----------------------------------------
  // QuadTree2D: quadtree over [-1,1]^2 with Quad9 geometry
  //----------------------------------------
  class QuadTree2D {
    public:
      // Construct with maxDepth; tree starts as a single root
      QuadTree2D(u32 maxDepth = 28)
        : _maxDepth(maxDepth) {
        _nodes.reserve(1024);
        _alive.reserve(1024);
        _free.reserve(256);
        _root = alloc_node();
        _nodes[_root].level = 0;
        _nodes[_root].morton = 0;
      }

      // Construct with maxDepth and initial minDepth (soft floor)
      QuadTree2D(u32 maxDepth, u32 minDepth)
        : _maxDepth(maxDepth), _minDepth(std::min(minDepth, maxDepth)) {
        _nodes.reserve(1024);
        _alive.reserve(1024);
        _free.reserve(256);
        _root = alloc_node();
        _nodes[_root].level = 0;
        _nodes[_root].morton = 0;
      }

      // Set minimum allowed depth for automatic floor
      void set_min_depth(u32 d) {
        _minDepth = std::min(d, _maxDepth);
      }

      // Get current minimum depth floor
      u32  min_depth() const {
        return _minDepth;
      }

      // Allow/disallow coarsening below min depth
      void set_allow_coarsen_below_min(bool v) {
        _allowCoarsenBelowMinDepth = v;
      }

      // Query flag for coarsening below min depth
      bool allow_coarsen_below_min() const {
        return _allowCoarsenBelowMinDepth;
      }

      // Set physical coordinates (Q2 geometry) for the coarse/root element
      void set_physical_quad9(const double x9[9], const double y9[9]) {
        for (int i = 0; i < 9; ++i) {
          _X[i] = x9[i];
          _Y[i] = y9[i];
        }
        _geom_ready = true;
      }

      // Assert geometry is defined
      inline void require_geometry() const {
        assert(_geom_ready && "Call set_physical_quad9(...) before geometric operations.");
      }

      // Refine a leaf into 4 children (Morton ordering), returns false if not allowed
      bool refine(u32 i) {
        if (i == npos32 || !_alive[i]) return false;
        if (!_nodes[i].is_leaf() || _nodes[i].level >= _maxDepth) return false;

        const u32 parent_level  = _nodes[i].level;
        const u64 parent_morton = _nodes[i].morton;

        for (int k = 0; k < 4; ++k) {
          const u32 c = alloc_node();
          _nodes[i].child[k] = c;
          QuadNode& cn = _nodes[c];
          cn.parent = i;
          cn.level  = parent_level + 1;
          cn.morton = (parent_morton << 2)
                      | (u64)((k & 2) ? 2 : 0)
                      | (u64)((k & 1) ? 1 : 0);
        }
        _leaf_dirty = true;
        return true;
      }

      // Coarsen a parent whose 4 children are leaves; returns false if not possible
      bool coarsen(u32 i) {
        if (i == npos32 || !_alive[i]) return false;
        if (_nodes[i].is_leaf()) return false;
        if (_nodes[i].level <= _minDepth && !_allowCoarsenBelowMinDepth) return false;

        std::array<u32, 4> ch = _nodes[i].child;

        for (int k = 0; k < 4; ++k) {
          const u32 c = ch[k];
          if (c == npos32) return false;
          if (!_alive[c])  return false;
          if (!_nodes[c].is_leaf()) return false;
          if (_nodes[c].parent != i) return false;
        }

        for (int k = 0; k < 4; ++k) free_node(ch[k]);
        _nodes[i].child = {npos32, npos32, npos32, npos32};
        _leaf_dirty = true;
        return true;
      }

      // Build/refresh the compact list of leaves (indices into _nodes)
      const std::vector<u32>& leaves() const {
        if (_leaf_dirty) {
          _leaves.clear();
          _leaves.reserve(_nodes.size());
          for (u32 i = 0; i < (u32)_nodes.size(); ++i)
            if (_alive[i] && _nodes[i].is_leaf())
              _leaves.push_back(i);

          _leaf_pos.assign(_nodes.size(), npos32);
          for (u32 pos = 0; pos < (u32)_leaves.size(); ++pos) _leaf_pos[_leaves[pos]] = pos;
          _leaf_dirty = false;
        }
        return _leaves;
      }

      // Perform one refinement pass using a predicate on current leaves
      template<class Pred>
      std::size_t refine_pass(Pred&& should_refine, Basis probe_basis = Basis::Q2_Quad9) {
        const auto& L = leaves();
        std::vector<u32> to_refine;
        to_refine.reserve(L.size());

        for (u32 k = 0; k < (u32)L.size(); ++k) {
          const u32 leaf = L[k];

          std::vector<std::array<double, 2>> pts_xi, pts_xy;
          std::vector<std::array<double, 9>> Nvals;
          leaf_parent_nodes(probe_basis, leaf, pts_xi);
          leaf_physical_nodes(probe_basis, leaf, pts_xy);
          quad9_shapes_at(pts_xi, Nvals);

          if (should_refine(leaf, pts_xi, pts_xy, Nvals)) to_refine.push_back(leaf);
        }

        std::size_t refined = 0;
        for (u32 leaf : to_refine) if (refine(leaf)) ++refined;
        (void)leaves();
        return refined;
      }

      // Return axis-aligned bounds of a node in parent coordinates
      inline void leaf_bounds(u32 i, double& xi0, double& eta0, double& xi1, double& eta1) const {
        const QuadNode& n = _nodes[i];
        u32 ix, iy;
        deinterleave2(n.morton, ix, iy);
        const double N = double(1u << n.level);
        const double dx = 2.0 / N, dy = 2.0 / N;
        xi0  = -1.0 + ix * dx;
        eta0 = -1.0 + iy * dy;
        xi1  = xi0 + dx;
        eta1 = eta0 + dy;
      }

      // Locate leaf containing (xi,eta) in the parent space [-1,1]^2 using Morton codes
      u32 locate_leaf_on_parent(double xi, double eta) const {
        if (xi < -1.0 || xi > 1.0 || eta < -1.0 || eta > 1.0) return npos32;

        // Scale to integer grid at max depth
        const double scale = double(1u << _maxDepth) / 2.0; // maps [-1,1] -> [0,2^maxDepth)
        u32 ix = std::min<u32>((u32)((xi + 1.0) * scale), (1u << _maxDepth) - 1);
        u32 iy = std::min<u32>((u32)((eta + 1.0) * scale), (1u << _maxDepth) - 1);

        u32 node = _root;
        for (u32 level = 0; level < _maxDepth; ++level) {
          const auto& n = _nodes[node];
          if (n.is_leaf()) return node;

          // pick bit for this level (from MSB downwards)
          const int shift = _maxDepth - level - 1;
          const int qx = (ix >> shift) & 1u;
          const int qy = (iy >> shift) & 1u;
          const int q  = (qy << 1) | qx; // 0=LL,1=LR,2=UL,3=UR

          const u32 c = n.child[q];
          if (c == npos32) return node;
          node = c;
        }
        return node;
      }


      // Locate leaf containing (xi,eta) in parent space [-1,1]^2
      // and compute leaf-local coordinates (shat,that) in [-1,1]^2.
      // Returns false if (xi,eta) outside domain.
      bool locate_leaf_on_parent_and_ref(double xi, double eta,
                                         u32& leaf,
                                         double& shat, double& that) const {
        // Reject outside parent element
        if (xi < -1.0 || xi > 1.0 || eta < -1.0 || eta > 1.0) {
          leaf = npos32;
          return false;
        }

        // Step 1. Find the leaf index (fast Morton-based descent)
        leaf = locate_leaf_on_parent(xi, eta);
        if (leaf == npos32) return false;

        // Step 2. Compute leaf bounds once
        double xi0, eta0, xi1, eta1;
        leaf_bounds(leaf, xi0, eta0, xi1, eta1);

        // Step 3. Map to leaf reference coordinates ([-1,1]²)
        to_leaf_ref(xi, eta, xi0, eta0, xi1, eta1, shat, that);

        return true;
      }



      // Locate leaf containing physical point (x,y)
      u32 locate_leaf_on_physical(double x, double y) const {
        double xi, eta;
        if (!inverse_map_quad9(x, y, xi, eta)) return npos32;
        return locate_leaf_on_parent(xi, eta);
      }

      // Inverse map physical (x,y) -> (xi,eta) using Newton on Q2 geometry
      bool inverse_map_quad9(double x, double y, double& xi, double& eta,
                             double tol = 1e-12, int maxit = 25) const {
        require_geometry();
        xi  = 0.0;
        eta = 0.0;
        for (int it = 0; it < maxit; ++it) {
          double N9[9], dNdxi[9], dNdeta[9];
          Quad9Shape::N(xi, eta, N9);
          Quad9Shape::dN(xi, eta, dNdxi, dNdeta);

          double Xh = 0.0, Yh = 0.0, dXdxi = 0.0, dXdeta = 0.0, dYdxi = 0.0, dYdeta = 0.0;
          for (int a = 0; a < 9; ++a) {
            Xh += N9[a] * _X[a];
            Yh += N9[a] * _Y[a];
            dXdxi  += dNdxi[a] * _X[a];
            dXdeta += dNdeta[a] * _X[a];
            dYdxi  += dNdxi[a] * _Y[a];
            dYdeta += dNdeta[a] * _Y[a];
          }

          const double rx = Xh - x, ry = Yh - y;
          const double norm = std::sqrt(rx * rx + ry * ry);
          if (norm < tol) return true;

          const double J00 = dXdxi, J01 = dXdeta;
          const double J10 = dYdxi, J11 = dYdeta;
          const double detJ = J00 * J11 - J01 * J10;
          if (std::abs(detJ) <= 1e-30) return false;

          const double dxi  = (J11 * rx - J01 * ry) / detJ;
          const double deta = (-J10 * rx + J00 * ry) / detJ;
          xi  -= dxi;
          eta -= deta;
          if (std::abs(dxi) + std::abs(deta) < tol * 1e-6) return true;
        }
        return false;
      }

      // Map parent (xi,eta) to leaf-local (shat,that) in [-1,1]^2 given leaf bounds
      static inline void to_leaf_ref(double xi, double eta,
                                     double xi0, double eta0, double xi1, double eta1,
                                     double& shat, double& that) {
        const double cx = 0.5 * (xi0 + xi1), hx = 0.5 * (xi1 - xi0);
        const double cy = 0.5 * (eta0 + eta1), hy = 0.5 * (eta1 - eta0);
        shat = (xi  - cx) / hx;
        that = (eta - cy) / hy;
      }

      // Add a field and return its id; caller later resizes coeffs appropriately
      u32 add_field(Basis b) {
        Field f;
        f.basis = b;
        f.dofs_per_cell = (b == Basis::Q1_Quad4) ? 4 : (b == Basis::Serendipity8 ? 8 : 9);
        _fields.push_back(std::move(f));
        return (u32)(_fields.size() - 1);
      }

      // Mutable access to field by id
      Field& field(u32 fid) {
        return _fields[fid];
      }

      // Const access to field by id
      const Field& field(u32 fid) const {
        return _fields[fid];
      }

      // Get pointer into leaf's coefficient block for field fid (mutable)
      double* leaf_coeff_ptr(u32 fid, u32 leaf_pos) {
        Field& f = _fields[fid];
        return f.coeffs.data() + size_t(leaf_pos) * f.dofs_per_cell;
      }

      // Get pointer into leaf's coefficient block for field fid (const)
      const double* leaf_coeff_ptr(u32 fid, u32 leaf_pos) const {
        const Field& f = _fields[fid];
        return f.coeffs.data() + size_t(leaf_pos) * f.dofs_per_cell;
      }

      // Evaluate field fid at physical (x,y); returns false if outside or inverse fails
      bool evaluate_field_on_physical(u32 fid, double x, double y, double& value) const {
        require_geometry();
        double xi, eta;
        if (!inverse_map_quad9(x, y, xi, eta)) return false;
        const u32 leaf = locate_leaf_on_parent(xi, eta);
        if (leaf == npos32) return false;

        double xi0, eta0, xi1, eta1;
        leaf_bounds(leaf, xi0, eta0, xi1, eta1);
        double shat, that;
        to_leaf_ref(xi, eta, xi0, eta0, xi1, eta1, shat, that);

        const u32 leaf_pos = leaf_position(leaf);
        if (leaf_pos == npos32) return false;

        const Field& f = _fields[fid];
        const double* c = leaf_coeff_ptr(fid, leaf_pos);

        switch (f.basis) {
          case Basis::Q1_Quad4: {
            double N4[4];
            Shapes::Q1(shat, that, N4);
            value = N4[0] * c[0] + N4[1] * c[1] + N4[2] * c[2] + N4[3] * c[3];
          }
          break;
          case Basis::Serendipity8: {
            double N8[8];
            Shapes::Serendipity8(shat, that, N8);
            double v = 0.0;
            for (int i = 0; i < 8; ++i) v += N8[i] * c[i];
            value = v;
          }
          break;
          case Basis::Q2_Quad9: {
            double N9[9];
            Quad9Shape::N(shat, that, N9);
            double v = 0.0;
            for (int i = 0; i < 9; ++i) v += N9[i] * c[i];
            value = v;
          }
          break;
        }
        return true;
      }

      // Return current leaf count
      u32 leaf_count() {
        return (u32)leaves().size();
      }

      // Resize all field coefficient vectors to match current leaves (keeps values where possible)
      void resize_fields_to_leaves() {
        const auto nL = leaf_count();
        for (auto& f : _fields) f.coeffs.resize(size_t(nL) * f.dofs_per_cell);
      }

      // Expose leaf indices (const reference to compact leaf list)
      const std::vector<u32>& leaf_indices() const {
        return leaves();
      }

      // Fill vector with parent coordinates of interpolation nodes for given leaf+basis
      void leaf_parent_nodes(Basis basis, u32 leaf,
                             std::vector<std::array<double, 2>>& out_pts) const {
        double xi0, eta0, xi1, eta1;
        leaf_bounds(leaf, xi0, eta0, xi1, eta1);
        const double cx = 0.5 * (xi0 + xi1), hx = 0.5 * (xi1 - xi0);
        const double cy = 0.5 * (eta0 + eta1), hy = 0.5 * (eta1 - eta0);

        auto toParent = [&](double s, double t) {
          return std::array<double, 2> { cx + s * hx, cy + t * hy };
        };

        out_pts.clear();
        switch (basis) {
          case Basis::Q1_Quad4: {
            static const double S[4][2] = {{-1, -1}, {+1, -1}, {+1, +1}, {-1, +1}};
            out_pts.reserve(4);
            for (int i = 0; i < 4; ++i) out_pts.push_back(toParent(S[i][0], S[i][1]));
          }
          break;
          case Basis::Serendipity8: {
            static const double S[8][2] = {
              {-1, -1}, {+1, -1}, {+1, +1}, {-1, +1},
              {0, -1}, {+1, 0}, {0, +1}, {-1, 0}
            };
            out_pts.reserve(8);
            for (int i = 0; i < 8; ++i) out_pts.push_back(toParent(S[i][0], S[i][1]));
          }
          break;
          case Basis::Q2_Quad9: {
            static const double S[9][2] = {
              {-1, -1}, {+1, -1}, {+1, +1}, {-1, +1},
              {0, -1}, {+1, 0}, {0, +1}, {-1, 0}, {0, 0}
            };
            out_pts.reserve(9);
            for (int i = 0; i < 9; ++i) out_pts.push_back(toParent(S[i][0], S[i][1]));
          }
          break;
        }
      }

      // Map given parent nodes to physical coordinates via Q2 geometry
      void leaf_physical_nodes(Basis basis, u32 leaf,
                               std::vector<std::array<double, 2>>& out_xy) const {
        require_geometry();
        std::vector<std::array<double, 2>> pts;
        leaf_parent_nodes(basis, leaf, pts);
        out_xy.resize(pts.size());
        for (size_t i = 0; i < pts.size(); ++i) {
          const double xi = pts[i][0], eta = pts[i][1];
          double N9[9];
          Quad9Shape::N(xi, eta, N9);
          double X = 0.0, Y = 0.0;
          for (int a = 0; a < 9; ++a) {
            X += N9[a] * _X[a];
            Y += N9[a] * _Y[a];
          }
          out_xy[i] = {X, Y};
        }
      }

      // Evaluate Q2 shape arrays for a list of parent points
      void quad9_shapes_at(const std::vector<std::array<double, 2>>& points,
                           std::vector<std::array<double, 9>>& out_N) const {
        out_N.resize(points.size());
        for (size_t i = 0; i < points.size(); ++i) {
          double N9[9];
          Quad9Shape::N(points[i][0], points[i][1], N9);
          for (int a = 0; a < 9; ++a) out_N[i][a] = N9[a];
        }
      }

      // Refine all leaves where predicate is true
      template<class Pred>
      void refine_where(Pred&& should_refine, Basis probe_basis = Basis::Q2_Quad9) {
        const auto& L = leaves();
        std::vector<u32> to_refine;
        to_refine.reserve(L.size());
        for (u32 k = 0; k < (u32)L.size(); ++k) {
          u32 leaf = L[k];
          std::vector<std::array<double, 2>> pts_xi, pts_xy;
          std::vector<std::array<double, 9>> Nvals;
          leaf_parent_nodes(probe_basis, leaf, pts_xi);
          leaf_physical_nodes(probe_basis, leaf, pts_xy);
          quad9_shapes_at(pts_xi, Nvals);
          if (should_refine(leaf, pts_xi, pts_xy, Nvals)) to_refine.push_back(leaf);
        }
        for (u32 leaf : to_refine) refine(leaf);
        (void)leaves();
      }

      // Coarsen-pass over parents with 4 leaf children using provided predicate
      template<class Pred>
      std::size_t coarsen_pass(Pred&& should_coarsen,
                               Basis probe_basis = Basis::Q2_Quad9) {
        std::vector<u32> to_coarsen;
        to_coarsen.reserve(_nodes.size() / 4 + 1);

        for (u32 i = 0; i < (u32)_nodes.size(); ++i) {
          if (!_alive[i]) continue;
          const QuadNode& n = _nodes[i];
          if (n.is_leaf()) continue;

          bool ok = true;
          for (int k = 0; k < 4; ++k) {
            u32 c = n.child[k];
            if (c == npos32 || !_alive[c] || !_nodes[c].is_leaf() || _nodes[c].parent != i) {
              ok = false;
              break;
            }
          }
          if (!ok) continue;

          std::vector<std::array<double, 2>> pts_xi, pts_xy;
          std::vector<std::array<double, 9>> Nvals;
          leaf_parent_nodes(probe_basis, i, pts_xi);
          leaf_physical_nodes(probe_basis, i, pts_xy);
          quad9_shapes_at(pts_xi, Nvals);

          bool ok_min = !(n.level <= _minDepth && !_allowCoarsenBelowMinDepth);
          if (ok_min && should_coarsen(i, n.level, pts_xi, pts_xy, Nvals)) {
            to_coarsen.push_back(i);
          }
        }

        std::size_t done = 0;
        for (u32 p : to_coarsen) if (coarsen(p)) ++done;
        (void)leaves();
        return done;
      }

      // Multilevel refinement until convergence or pass limit
      template<class Pred>
      std::size_t adapt_refine_until(Pred&& should_refine,
                                     Basis probe_basis = Basis::Q2_Quad9,
                                     u32 max_passes = 10) {
        ensure_min_depth();
        std::size_t total_refined = 0;
        for (u32 pass = 0; pass < max_passes; ++pass) {
          const std::size_t r = refine_pass(should_refine, probe_basis);
          total_refined += r;
          if (r == 0) break;
        }
        return total_refined;
      }

      // Full adaptivity cycle: coarsen then refine for up to max_passes
      template<class PredCoarsen, class PredRefine>
      std::size_t adapt_cycle(PredCoarsen&& should_coarsen,
                              PredRefine&& should_refine,
                              Basis probe_basis = Basis::Q2_Quad9,
                              u32 max_passes = 10) {
        ensure_min_depth();
        std::size_t total = 0;
        for (u32 pass = 0; pass < max_passes; ++pass) {
          std::size_t c = coarsen_pass(should_coarsen, probe_basis);
          std::size_t r = refine_pass(should_refine, probe_basis);
          total += c + r;
          if (c == 0 && r == 0) break;
        }
        return total;
      }




      // Fill vector with parent-space coordinates of interpolation nodes (by basis) for leaf
      void leaf_reference_nodes(Basis basis, u32 leaf,
                                std::vector<std::array<double, 2>>& xi) const {
        xi.clear();
        double xi0, eta0, xi1, eta1;
        leaf_bounds(leaf, xi0, eta0, xi1, eta1);
        switch (basis) {
          case Basis::Q1_Quad4: {
            xi.resize(4);
            xi[0] = {xi0, eta0};
            xi[1] = {xi1, eta0};
            xi[2] = {xi1, eta1};
            xi[3] = {xi0, eta1};
          }
          break;
          case Basis::Serendipity8: {
            xi.resize(8);
            double xm = 0.5 * (xi0 + xi1), ym = 0.5 * (eta0 + eta1);
            xi[0] = {xi0, eta0};
            xi[1] = {xi1, eta0};
            xi[2] = {xi1, eta1};
            xi[3] = {xi0, eta1};
            xi[4] = {xm,  eta0};
            xi[5] = {xi1, ym };
            xi[6] = {xm,  eta1};
            xi[7] = {xi0, ym };
          }
          break;
          case Basis::Q2_Quad9: {
            xi.resize(9);
            double xm = 0.5 * (xi0 + xi1), ym = 0.5 * (eta0 + eta1);
            xi[0] = {xi0, eta0};
            xi[1] = {xi1, eta0};
            xi[2] = {xi1, eta1};
            xi[3] = {xi0, eta1};
            xi[4] = {xm,  eta0};
            xi[5] = {xi1, ym };
            xi[6] = {xm,  eta1};
            xi[7] = {xi0, ym };
            xi[8] = {xm,  ym };
          }
          break;
        }
      }

      // Evaluate a field on a known leaf using (shat,that) in leaf-local [-1,1]^2
      bool evaluate_field_on_leaf(u32 fid, u32 leaf, double shat, double that, double& value) const {
        const u32 leaf_pos = leaf_position(leaf);
        if (leaf_pos == npos32) return false;
        const Field& f = _fields[fid];
        const double* c = _fields[fid].coeffs.data() + size_t(leaf_pos) * f.dofs_per_cell;
        switch (f.basis) {
          case Basis::Q1_Quad4: {
            double N4[4];
            Shapes::Q1(shat, that, N4);
            value = N4[0] * c[0] + N4[1] * c[1] + N4[2] * c[2] + N4[3] * c[3];
          }
          break;
          case Basis::Serendipity8: {
            double N8[8];
            Shapes::Serendipity8(shat, that, N8);
            double v = 0.0;
            for (int i = 0; i < 8; ++i) v += N8[i] * c[i];
            value = v;
          }
          break;
          case Basis::Q2_Quad9: {
            double N9[9];
            Quad9Shape::N(shat, that, N9);
            double v = 0.0;
            for (int i = 0; i < 9; ++i) v += N9[i] * c[i];
            value = v;
          }
          break;
        }
        return true;
      }

      // Evaluate field directly in parent coordinates (xi,eta)
      // Uses fused locate+ref mapping for efficiency
      bool evaluate_field_on_parent(u32 fid, double xi, double eta, double& value) const {
        u32 leaf;
        double shat, that;

        // Step 1. Locate leaf and compute local ref coords in one shot
        if (!locate_leaf_on_parent_and_ref(xi, eta, leaf, shat, that)) {
          return false;
        }

        // Step 2. Map leaf index to position in coefficient storage
        const u32 leaf_pos = leaf_position(leaf);
        if (leaf_pos == npos32) return false;

        // Step 3. Access field coefficients
        const Field& f = _fields[fid];
        const double* c = leaf_coeff_ptr(fid, leaf_pos);

        // Step 4. Evaluate basis interpolation
        switch (f.basis) {
          case Basis::Q1_Quad4: {
            double N4[4];
            Shapes::Q1(shat, that, N4);
            value = N4[0] * c[0] + N4[1] * c[1] + N4[2] * c[2] + N4[3] * c[3];
          }
          break;

          case Basis::Serendipity8: {
            double N8[8];
            Shapes::Serendipity8(shat, that, N8);
            double v = 0.0;
            for (int i = 0; i < 8; ++i) v += N8[i] * c[i];
            value = v;
          }
          break;

          case Basis::Q2_Quad9: {
            double N9[9];
            Quad9Shape::N(shat, that, N9);
            double v = 0.0;
            for (int i = 0; i < 9; ++i) v += N9[i] * c[i];
            value = v;
          }
          break;
        }
        return true;
      }



      // Rebuild field fid on *this* from source tree 'src' by sampling at parent nodes
      void rebuild_field_from(const QuadTree2D& src, u32 fid) {
        resize_fields_to_leaves();
        Field& f = field(fid);

        const auto& L = leaf_indices();
        for (u32 k = 0; k < L.size(); ++k) {
          u32 leaf = L[k];
          std::vector<std::array<double, 2>> xi;
          leaf_reference_nodes(f.basis, leaf, xi);

          double* coeffs = leaf_coeff_ptr(fid, k);
          for (size_t j = 0; j < xi.size(); ++j) {
            double val;
            if (!src.evaluate_field_on_parent(fid, xi[j][0], xi[j][1], val)) {
              val = 0.0;
              std::cout << "error!";
            }
            coeffs[j] = val;
          }
        }
      }


      // Conservative coarsen cycle using snapshot + parent coords; rebuild all fields
      std::size_t coarsen_only_cycle_safe(u32 fid,
                                          double tau_coarse,
                                          u32 max_passes = 10,
                                          Basis probe_basis = Basis::Q2_Quad9) {


        QuadTree2D snapshot = *this;



        // Predicate: coarsening criterion
        auto pred = [&](u32 /*parent*/, u32 level,
                        const std::vector<std::array<double, 2>>& pts_xi,
                        const std::vector<std::array<double, 2>>& /*pts_xy*/,
        const std::vector<std::array<double, 9>>& /*Nvals*/) -> bool {
          if (level <= min_depth()) return false;
          if (pts_xi.empty()) return false;

          double v0;
          if (!snapshot.evaluate_field_on_parent(fid, pts_xi[0][0], pts_xi[0][1], v0))
            return false;

          double mn = v0, mx = v0;
          for (size_t i = 1; i < pts_xi.size(); ++i) {
            double val;
            if (snapshot.evaluate_field_on_parent(fid, pts_xi[i][0], pts_xi[i][1], val)) {
              mn = std::min(mn, val);
              mx = std::max(mx, val);
            }
          }
          return (mn > +tau_coarse) || (mx < -tau_coarse);
        };

        std::size_t total = 0;

        // Perform coarsening passes
        for (u32 pass = 0; pass < max_passes; ++pass) {
          std::size_t c = coarsen_pass(pred, probe_basis);
          if (c == 0) break;

          total += c;

          // 🔧 Enforce 1-irregularity after each coarsening pass
          enforce_balance();
        }

        // Rebuild all fields from snapshot (conservative transfer)
        for (u32 f = 0; f < _fields.size(); ++f) {
          rebuild_field_from(snapshot, f);
        }

        //enforce_hanging_constraints();

        return total;
      }


      u32 find_neighbor(u32 leaf, int edge) const {
        double xi0, eta0, xi1, eta1;
        leaf_bounds(leaf, xi0, eta0, xi1, eta1);

        double eps = 1.0E-10;

        // local element sizes in parent space
        double dx = xi1 - xi0;
        double dy = eta1 - eta0;

        // relative perturbation factor (safe wrt machine precision)


        // test point starts at element center
        double xm = 0.5 * (xi0 + xi1);
        double ym = 0.5 * (eta0 + eta1);

        // shift across the requested edge
        switch (edge) {
          case 0:
            ym = eta0 - eps * dy;
            break; // bottom
          case 1:
            xm = xi1 + eps * dx;
            break; // right
          case 2:
            ym = eta1 + eps * dy;
            break; // top
          case 3:
            xm = xi0 - eps * dx;
            break; // left
          default:
            break;
        }

        return locate_leaf_on_parent(xm, ym);



      }


      // Collect edge nodes of a leaf and return both coordinates (xi,eta) and local indices.
      // The coordinates are shifted slightly outside the fine element so evaluation
      // comes from the coarse neighbor.
      void edge_reference_nodes_with_mapping(Basis basis, u32 leaf, int edge,
                                             std::vector<std::array<double, 2>>& out_coords,
                                             std::vector<int>& out_indices) const {
        out_coords.clear();
        out_indices.clear();

        std::vector<std::array<double, 2>> xi;
        leaf_reference_nodes(basis, leaf, xi);

        double eps = 1e-10;  // small outward shift
        double xi0, eta0, xi1, eta1;
        leaf_bounds(leaf, xi0, eta0, xi1, eta1);
        const double dx = xi1 - xi0, dy = eta1 - eta0;

        for (int i = 0; i < (int)xi.size(); ++i) {
          auto p = xi[i];
          bool isOnEdge = false;

          if (edge == 0) { // bottom edge
            p[1] -= eps * dy;
            switch (i) {
              case 0: // bottom-left corner
                p[0] += eps * dx;
                isOnEdge = true;
                break;
              case 1: // bottom-right corner
                p[0] -= eps * dx;
                isOnEdge = true;
                break;
              case 4: // bottom midpoint
                isOnEdge = true;
                break;
              default:
                break;
            }
          }
          else if (edge == 1) { // right edge
            p[0] += eps * dx;
            switch (i) {
              case 1: // bottom-right corner
                p[1] += eps * dy;
                isOnEdge = true;
                break;
              case 2: // top-right corner
                p[1] -= eps * dy;
                isOnEdge = true;
                break;
              case 5: // right midpoint
                isOnEdge = true;
                break;
              default:
                break;
            }
          }
          else if (edge == 2) { // top edge
            p[1] += eps * dy;
            switch (i) {
              case 2: // top-right corner
                p[0] -= eps * dx;
                isOnEdge = true;
                break;
              case 3: // top-left corner
                p[0] += eps * dx;
                isOnEdge = true;
                break;
              case 6: // top midpoint
                isOnEdge = true;
                break;
              default:
                break;
            }
          }
          else if (edge == 3) { // left edge
            p[0] -= eps * dx;
            switch (i) {
              case 0: // bottom-left corner
                p[1] += eps * dy;
                isOnEdge = true;
                break;
              case 3: // top-left corner
                p[1] -= eps * dy;
                isOnEdge = true;
                break;
              case 7: // left midpoint
                isOnEdge = true;
                break;
              default:
                break;
            }
          }

          if (isOnEdge) {
            out_coords.push_back(p);   // store nudged point
            out_indices.push_back(i); // store local index of node
          }
        }
      }

      // Enforce hanging-node constraints for all fields
      void enforce_hanging_constraints() {
        for (u32 leaf : leaf_indices()) {
          for (int e = 0; e < 4; ++e) {
            u32 neigh = find_neighbor(leaf, e);
            if (neigh == npos32) continue;

            if (level_of(leaf) > level_of(neigh)) {
              enforce_edge_constraints(neigh, leaf, e);
            }
          }
        }
      }


      void enforce_edge_constraints(u32 coarse, u32 fine, int edge) {
        for (u32 fid = 0; fid < _fields.size(); ++fid) {
          const Field& f = _fields[fid];
          double* coeffF = leaf_coeff_ptr(fid, leaf_position(fine));

          std::vector<std::array<double, 2>> xiFine;
          std::vector<int> idxFine;
          edge_reference_nodes_with_mapping(f.basis, fine, edge, xiFine, idxFine);

          for (size_t j = 0; j < xiFine.size(); ++j) {
            double val;
            if (!evaluate_field_on_parent(fid, xiFine[j][0], xiFine[j][1], val)) {
              val = 0.0;
              std::cout << "error!";
            }




            coeffF[idxFine[j]] = val;  // write into correct slot
          }
        }
      }





      // // Conservative coarsen cycle using snapshot + parent coords; rebuild all fields
      // std::size_t coarsen_only_cycle_safe(u32 fid,
      //                                     double tau_coarse,
      //                                     u32 max_passes = 10,
      //                                     Basis probe_basis = Basis::Q2_Quad9) {
      //   QuadTree2D snapshot = *this;
      //
      //   auto pred = [&](u32 /*parent*/, u32 level,
      //                   const std::vector<std::array<double, 2>>& pts_xi,
      //                   const std::vector<std::array<double, 2>>& /*pts_xy*/,
      //   const std::vector<std::array<double, 9>>& /*Nvals*/) -> bool {
      //     if (level <= min_depth()) return false;
      //     if (pts_xi.empty()) return false;
      //
      //     double v0;
      //     if (!snapshot.evaluate_field_on_parent(fid, pts_xi[0][0], pts_xi[0][1], v0))
      //       return false;
      //     double mn = v0, mx = v0;
      //     for (size_t i = 1; i < pts_xi.size(); ++i) {
      //       double val;
      //       if (snapshot.evaluate_field_on_parent(fid, pts_xi[i][0], pts_xi[i][1], val)) {
      //         mn = std::min(mn, val);
      //         mx = std::max(mx, val);
      //       }
      //     }
      //     return (mn > +tau_coarse) || (mx < -tau_coarse);
      //   };
      //
      //   std::size_t total = 0;
      //   for (u32 pass = 0; pass < max_passes; ++pass) {
      //     std::size_t c = coarsen_pass(pred, probe_basis);
      //     total += c;
      //     if (c == 0) break;
      //   }
      //
      //   for (u32 f = 0; f < _fields.size(); ++f) rebuild_field_from(snapshot, f);
      //   return total;
      // }
      //
      u32 level_of(u32 leaf) const {
        return _nodes[leaf].level;
      }


      // Map a parent coordinate (xi,eta) in [-1,1]^2 to physical (x,y)
// using the Quad9 isoparametric geometry map.
      std::array<double, 2> parent_to_physical(double xi, double eta) const {
        require_geometry();

        double N9[9];
        Quad9Shape::N(xi, eta, N9);

        double X = 0.0, Y = 0.0;
        for (int a = 0; a < 9; ++a) {
          X += N9[a] * _X[a];
          Y += N9[a] * _Y[a];
        }
        return {X, Y};
      }



      // Collect reference coordinates (xi,eta) of nodes for all leaves
// in a level range. Much like extract_node_coords_in_level_range,
// but stays in parent reference space (avoids inverse_map).
      std::vector<std::array<double, 2>>
      extract_node_parent_coords_in_level_range(u32 lev_min, u32 lev_max, Basis basis) const {
        std::vector<std::array<double, 2>> coords;

        for (u32 leaf : leaf_indices()) {
          u32 lev = level_of(leaf);          // use accessor, not _nodes
          if (lev < lev_min || lev > lev_max) continue;

          std::vector<std::array<double, 2>> xi;
          leaf_reference_nodes(basis, leaf, xi);

          coords.insert(coords.end(), xi.begin(), xi.end());
        }

        return coords;
      }



      // ---------------------------------------------------------
// Neighbor lookup (axis-aligned): dir = 0:left, 1:right, 2:down, 3:up
// Returns the leaf covering the neighbor cell, or npos32 if outside.
// ---------------------------------------------------------
      u32 neighbor_leaf(u32 leaf, int dir) const {
        const QuadNode& n = _nodes[leaf];

        // Decode Morton index -> (ix, iy)
        u32 ix, iy;
        deinterleave2(n.morton, ix, iy);

        // Step one cell at this level
        const u32 N = 1u << n.level;
        if (dir == 0) {
          if (ix == 0) return npos32;  // left
          else ix -= 1;
        }
        if (dir == 1) {
          if (ix + 1 >= N) return npos32;  // right
          else ix += 1;
        }
        if (dir == 2) {
          if (iy == 0) return npos32;  // down
          else iy -= 1;
        }
        if (dir == 3) {
          if (iy + 1 >= N) return npos32;  // up
          else iy += 1;
        }

        // Map to parent-space coordinate at neighbor cell center
        const double dx = 2.0 / double(N);
        const double dy = 2.0 / double(N);
        const double xi  = -1.0 + (ix + 0.5) * dx;
        const double eta = -1.0 + (iy + 0.5) * dy;

        // Locate leaf covering that center point
        return locate_leaf_on_parent(xi, eta);
      }

// ---------------------------------------------------------
// Enforce 1-irregularity: no adjacent leaves differ by >1 level
// ---------------------------------------------------------
      void enforce_balance() {
        bool changed = true;
        while (changed) {
          changed = false;
          const auto& L = leaves();

          for (u32 leaf : L) {
            u32 lev = _nodes[leaf].level;

            for (int dir = 0; dir < 4; ++dir) {
              u32 nb = neighbor_leaf(leaf, dir);
              if (nb == npos32) continue;

              u32 lev_nb = _nodes[nb].level;
              if (lev > lev_nb + 1) {
                if (refine(nb)) changed = true;
              }
            }
          }
          (void)leaves(); // refresh compact list
        }
      }



      // // Write current mesh + field to VTK UnstructuredGrid ASCII (.vtu)
      // bool write_vtu(const std::string & filename, u32 fid, const std::string & name,
      //                bool cell_centered = false) const {
      //   require_geometry();
      //   const auto& L = leaves();
      //   const size_t numCells = L.size();
      //   if (numCells == 0) return false;
      //
      //   const Field& fld = _fields[fid];
      //   const size_t need = numCells * (size_t)fld.dofs_per_cell;
      //   if (fld.coeffs.size() < need) {
      //     assert(false && "Field coefficients not sized to current leaves. Call resize_fields_to_leaves().");
      //     return false;
      //   }
      //
      //   struct Pt {
      //     double x, y, z;
      //   };
      //   std::vector<Pt> points;
      //   points.reserve(numCells * 4);
      //
      //   std::vector<int> connectivity;
      //   connectivity.reserve(numCells * 4);
      //   std::vector<int> offsets;
      //   offsets.reserve(numCells);
      //   std::vector<unsigned char> types;
      //   types.reserve(numCells);
      //
      //   std::vector<double> pointData;
      //   pointData.reserve(numCells * 4);
      //   std::vector<double> cellData;
      //   cellData.reserve(numCells);
      //
      //   std::vector<int> levels;
      //   levels.reserve(numCells);
      //
      //   static const double ST[4][2] = { {-1, -1}, {+1, -1}, {+1, +1}, {-1, +1} };
      //
      //   const int ndof = (int)fld.dofs_per_cell;
      //   const double* coeff_base = fld.coeffs.data();
      //
      //   for (size_t k = 0; k < numCells; ++k) {
      //     const u32 leaf = L[k];
      //     levels.push_back((int)_nodes[leaf].level);
      //
      //     const double* c = coeff_base + size_t(k) * ndof;
      //
      //     double xi0, eta0, xi1, eta1;
      //     leaf_bounds(leaf, xi0, eta0, xi1, eta1);
      //     const double cx = 0.5 * (xi0 + xi1), hx = 0.5 * (xi1 - xi0);
      //     const double cy = 0.5 * (eta0 + eta1), hy = 0.5 * (eta1 - eta0);
      //
      //     for (int q = 0; q < 4; ++q) {
      //       const double sh = ST[q][0], th = ST[q][1];
      //       const double xi  = cx + sh * hx;
      //       const double eta = cy + th * hy;
      //
      //       double N9[9];
      //       Quad9Shape::N(xi, eta, N9);
      //       double X = 0.0, Y = 0.0;
      //       for (int a = 0; a < 9; ++a) {
      //         X += N9[a] * _X[a];
      //         Y += N9[a] * _Y[a];
      //       }
      //       points.push_back({X, Y, 0.0});
      //
      //       if (!cell_centered) {
      //         double v = 0.0;
      //         switch (fld.basis) {
      //           case Basis::Q1_Quad4: {
      //             double N4[4];
      //             Shapes::Q1(sh, th, N4);
      //             v = N4[0] * c[0] + N4[1] * c[1] + N4[2] * c[2] + N4[3] * c[3];
      //           }
      //           break;
      //           case Basis::Serendipity8: {
      //             double N8[8];
      //             Shapes::Serendipity8(sh, th, N8);
      //             for (int a = 0; a < 8; ++a) v += N8[a] * c[a];
      //           }
      //           break;
      //           case Basis::Q2_Quad9: {
      //             double Nq[9];
      //             Quad9Shape::N(sh, th, Nq);
      //             for (int a = 0; a < 9; ++a) v += Nq[a] * c[a];
      //           }
      //           break;
      //         }
      //         pointData.push_back(v);
      //       }
      //     }
      //
      //     const int base = (int)(k * 4);
      //     connectivity.push_back(base + 0);
      //     connectivity.push_back(base + 1);
      //     connectivity.push_back(base + 2);
      //     connectivity.push_back(base + 3);
      //     offsets.push_back(base + 4);
      //     types.push_back(9); // VTK_QUAD
      //
      //     if (cell_centered) {
      //       double v = 0.0;
      //       switch (fld.basis) {
      //         case Basis::Q1_Quad4: {
      //           double N4[4];
      //           Shapes::Q1(0.0, 0.0, N4);
      //           v = N4[0] * c[0] + N4[1] * c[1] + N4[2] * c[2] + N4[3] * c[3];
      //         }
      //         break;
      //         case Basis::Serendipity8: {
      //           double N8[8];
      //           Shapes::Serendipity8(0.0, 0.0, N8);
      //           for (int a = 0; a < 8; ++a) v += N8[a] * c[a];
      //         }
      //         break;
      //         case Basis::Q2_Quad9: {
      //           double Nq[9];
      //           Quad9Shape::N(0.0, 0.0, Nq);
      //           for (int a = 0; a < 9; ++a) v += Nq[a] * c[a];
      //         }
      //         break;
      //       }
      //       cellData.push_back(v);
      //     }
      //   }
      //
      //   std::ofstream os(filename);
      //   if (!os) return false;
      //   os << std::setprecision(16);
      //   os << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
      //   os << "  <UnstructuredGrid>\n";
      //   os << "    <Piece NumberOfPoints=\"" << points.size()
      //      << "\" NumberOfCells=\"" << numCells << "\">\n";
      //
      //   os << "      <Points>\n";
      //   os << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
      //   for (const auto& p : points) os << "          " << p.x << " " << p.y << " " << p.z << "\n";
      //   os << "        </DataArray>\n";
      //   os << "      </Points>\n";
      //
      //   os << "      <Cells>\n";
      //   os << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
      //   for (size_t i = 0; i < connectivity.size(); i += 4)
      //     os << "          " << connectivity[i] << " " << connectivity[i + 1]
      //        << " " << connectivity[i + 2] << " " << connectivity[i + 3] << "\n";
      //   os << "        </DataArray>\n";
      //   os << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
      //   for (int off : offsets) os << "          " << off << "\n";
      //   os << "        </DataArray>\n";
      //   os << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
      //   for (unsigned char t : types) os << "          " << (int)t << "\n";
      //   os << "        </DataArray>\n";
      //   os << "      </Cells>\n";
      //
      //   os << "      <CellData Scalars=\"" << name << "\">\n";
      //   os << "        <DataArray type=\"Int32\" Name=\"level\" format=\"ascii\">\n";
      //   for (int lv : levels) os << "          " << lv << "\n";
      //   os << "        </DataArray>\n";
      //
      //   if (cell_centered) {
      //     os << "        <DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
      //     for (double v : cellData) os << "          " << v << "\n";
      //     os << "        </DataArray>\n";
      //     os << "      </CellData>\n";
      //   }
      //   else {
      //     os << "      </CellData>\n";
      //     os << "      <PointData Scalars=\"" << name << "\">\n";
      //     os << "        <DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
      //     for (double v : pointData) os << "          " << v << "\n";
      //     os << "        </DataArray>\n";
      //     os << "      </PointData>\n";
      //   }
      //
      //   os << "    </Piece>\n";
      //   os << "  </UnstructuredGrid>\n";
      //   os << "</VTKFile>\n";
      //   return true;
      // }

      // Write current mesh + field to VTK UnstructuredGrid ASCII (.vtu)
      bool write_vtu(const std::string & filename, u32 fid, const std::string & name,
                     bool cell_centered = false) const {
        require_geometry();
        const auto& L = leaves();
        const size_t numCells = L.size();
        if (numCells == 0) return false;

        const Field& fld = _fields[fid];
        const size_t need = numCells * (size_t)fld.dofs_per_cell;
        if (fld.coeffs.size() < need) {
          assert(false && "Field coefficients not sized to current leaves. Call resize_fields_to_leaves().");
          return false;
        }

        struct Pt {
          double x, y, z;
        };
        std::vector<Pt> points;
        std::vector<int> connectivity;
        std::vector<int> offsets;
        std::vector<unsigned char> types;
        std::vector<double> pointData;
        std::vector<double> cellData;
        std::vector<int> levels;

        const int ndof = (int)fld.dofs_per_cell;
        const double* coeff_base = fld.coeffs.data();
        int pointCounter = 0;

        for (size_t k = 0; k < numCells; ++k) {
          const u32 leaf = L[k];
          levels.push_back((int)_nodes[leaf].level);
          const double* c = coeff_base + size_t(k) * ndof;

          // Get geometry nodes depending on basis
          std::vector<std::array<double, 2>> xy;
          leaf_physical_nodes(fld.basis, leaf, xy);

          // Append nodes to point list
          for (auto &p : xy) {
            points.push_back({p[0], p[1], 0.0});
          }

          // Connectivity: sequential nodes
          for (int i = 0; i < ndof; ++i)
            connectivity.push_back(pointCounter + i);
          pointCounter += ndof;
          offsets.push_back((int)connectivity.size());

          // Cell type
          switch(fld.basis) {
            case Basis::Q1_Quad4:
              types.push_back(9);
              break;  // VTK_QUAD
            case Basis::Serendipity8:
              types.push_back(23);
              break;  // VTK_QUADRATIC_QUAD
            case Basis::Q2_Quad9:
              types.push_back(28);
              break;  // VTK_BIQUADRATIC_QUAD
          }

          if (cell_centered) {
            // Just evaluate at element center (0,0)
            double v = 0.0;
            switch(fld.basis) {
              case Basis::Q1_Quad4: {
                double N4[4];
                Shapes::Q1(0, 0, N4);
                for (int i = 0; i < 4; ++i) v += N4[i] * c[i];
              }
              break;
              case Basis::Serendipity8: {
                double N8[8];
                Shapes::Serendipity8(0, 0, N8);
                for (int i = 0; i < 8; ++i) v += N8[i] * c[i];
              }
              break;
              case Basis::Q2_Quad9: {
                double N9[9];
                Quad9Shape::N(0, 0, N9);
                for (int i = 0; i < 9; ++i) v += N9[i] * c[i];
              }
              break;
            }
            cellData.push_back(v);
          }
          else {
            // Push values at nodes
            for (int i = 0; i < ndof; ++i) pointData.push_back(c[i]);
          }
        }

        // Write ASCII VTK
        std::ofstream os(filename);
        if (!os) return false;
        os << std::setprecision(16);
        os << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        os << "  <UnstructuredGrid>\n";
        os << "    <Piece NumberOfPoints=\"" << points.size()
           << "\" NumberOfCells=\"" << numCells << "\">\n";

        os << "      <Points>\n";
        os << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (auto &p : points) os << "          " << p.x << " " << p.y << " " << p.z << "\n";
        os << "        </DataArray>\n";
        os << "      </Points>\n";

        os << "      <Cells>\n";
        os << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
        for (size_t i = 0; i < connectivity.size(); ++i) {
          os << connectivity[i] << ((i + 1) % ndof == 0 ? "\n" : " ");
        }
        os << "        </DataArray>\n";
        os << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
        for (int off : offsets) os << off << "\n";
        os << "        </DataArray>\n";
        os << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
        for (auto t : types) os << (int)t << "\n";
        os << "        </DataArray>\n";
        os << "      </Cells>\n";

        os << "      <CellData Scalars=\"" << name << "\">\n";
        os << "        <DataArray type=\"Int32\" Name=\"level\" format=\"ascii\">\n";
        for (int lv : levels) os << lv << "\n";
        os << "        </DataArray>\n";

        if (cell_centered) {
          os << "        <DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
          for (double v : cellData) os << v << "\n";
          os << "        </DataArray>\n";
          os << "      </CellData>\n";
        }
        else {
          os << "      </CellData>\n";
          os << "      <PointData Scalars=\"" << name << "\">\n";
          os << "        <DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
          for (double v : pointData) os << v << "\n";
          os << "        </DataArray>\n";
          os << "      </PointData>\n";
        }

        os << "    </Piece>\n";
        os << "  </UnstructuredGrid>\n";
        os << "</VTKFile>\n";
        return true;
      }


      // === Data extraction utilities ===

      // Test if two points are the same within tolerance (used for global gather)
      static inline bool same_point(const std::array<double, 2>& a,
                                    const std::array<double, 2>& b,
                                    double tol = 1e-12) {
        return (std::fabs(a[0] - b[0]) < tol && std::fabs(a[1] - b[1]) < tol);
      }

      // Return number of nodes for a given basis
      static int basis_nodes(Basis b) {
        switch (b) {
          case Basis::Q1_Quad4:
            return 4;
          case Basis::Serendipity8:
            return 8;
          case Basis::Q2_Quad9:
            return 9;
        }
        return 0;
      }

      // Gather unique coords/values across leaves into global arrays
      void gather_field(u32 fid) {
        require_geometry();
        struct Flat {
          std::array<double, 2> xy;
          double val;
          int elem, local;
        };
        std::vector<Flat> flat;

        const auto& L = leaves();
        const Field& f = _fields[fid];
        int nNodes = f.dofs_per_cell;

        for (int e = 0; e < (int)L.size(); ++e) {
          u32 leaf = L[e];
          std::vector<std::array<double, 2>> pts_xy;
          leaf_physical_nodes(f.basis, leaf, pts_xy);
          const double* coeffs = leaf_coeff_ptr(fid, e);
          for (int i = 0; i < nNodes; ++i) flat.push_back({pts_xy[i], coeffs[i], e, i});
        }

        std::sort(flat.begin(), flat.end(),
        [](auto & a, auto & b) {
          return (a.xy[0] < b.xy[0]) ||
                 (a.xy[0] == b.xy[0] && a.xy[1] < b.xy[1]);
        });

        global_coords.clear();
        global_field.clear();
        elem2glob.assign(L.size(), {});

        int gid = -1;
        for (size_t k = 0; k < flat.size(); ++k) {
          if (k == 0 || !same_point(flat[k].xy, flat[k - 1].xy)) {
            ++gid;
            global_coords.push_back(flat[k].xy);
            global_field.push_back(flat[k].val);
          }
          if (elem2glob[flat[k].elem].empty()) elem2glob[flat[k].elem].resize(nNodes);
          elem2glob[flat[k].elem][flat[k].local] = gid;
        }
      }

      // Scatter modified global_field back to leaf coefficient storage
      void scatter_field(u32 fid) {
        auto& f = _fields[fid];
        int nNodes = f.dofs_per_cell;
        const auto& L = leaves();
        for (int e = 0; e < (int)L.size(); ++e) {
          double* coeffs = leaf_coeff_ptr(fid, e);
          for (int i = 0; i < nNodes; ++i) {
            int gid = elem2glob[e][i];
            coeffs[i] = global_field[gid];
          }
        }
      }

      // Refine tree so that each point lies in a leaf of at least 'maxDepthTarget'
      void refine_to_contain_points(const std::vector<std::array<double, 2>>& pts,
                                    u32 maxDepthTarget) {
        for (const auto& p : pts) {
          double xi, eta;
          if (!inverse_map_quad9(p[0], p[1], xi, eta)) continue;

          u32 leaf = locate_leaf_on_parent(xi, eta);
          while (leaf != npos32 && _nodes[leaf].level < maxDepthTarget) {
            if (!refine(leaf)) break;
            leaf = locate_leaf_on_parent(xi, eta);
          }
        }
        (void)leaves();
        enforce_balance();
      }

      // Extract node coordinates (physical) for all leaves in [minLevel,maxLevel]
      std::vector<std::vector<std::array<double, 2>>>
      extract_leaf_nodes_in_level_range(u32 minLevel, u32 maxLevel, Basis basis) const {
        require_geometry();
        std::vector<std::vector<std::array<double, 2>>> result;
        const auto& L = leaves();
        for (u32 leaf : L) {
          u32 lev = _nodes[leaf].level;
          if (lev >= minLevel && lev <= maxLevel) {
            std::vector<std::array<double, 2>> xy;
            leaf_physical_nodes(basis, leaf, xy);
            result.push_back(std::move(xy));
          }
        }
        return result;
      }

      // Extract all node coordinates (physical) for leaves in [minLevel,maxLevel] (flat list)
      std::vector<std::array<double, 2>>
      extract_node_coords_in_level_range(u32 minLevel, u32 maxLevel, Basis basis) const {
        require_geometry();
        std::vector<std::array<double, 2>> coords;
        const auto& L = leaves();
        for (u32 leaf : L) {
          u32 lev = _nodes[leaf].level;
          if (lev >= minLevel && lev <= maxLevel) {
            std::vector<std::array<double, 2>> xy;
            leaf_physical_nodes(basis, leaf, xy);
            coords.insert(coords.end(), xy.begin(), xy.end());
          }
        }
        return coords;
      }

      // Reset to single root cell; keep existing fields but size to 1 leaf
      void reset() {
        _nodes.clear();
        _alive.clear();
        _free.clear();
        _leaves.clear();
        _leaf_dirty = true;

        _root = alloc_node();
        _nodes[_root].level = 0;
        _nodes[_root].morton = 0;

        for (auto& f : _fields) {
          f.coeffs.clear();
          f.coeffs.resize(f.dofs_per_cell);
        }
      }

      // Return global maximum depth allowed
      u32 max_depth() const {
        return _maxDepth;
      }

      // Accessors for coarse geometry coordinates (const-ref)
      inline const std::array<double, 9>& get_X() const {
        return _X;
      }
      inline const std::array<double, 9>& get_Y() const {
        return _Y;
      }

      // === Public state used by gather/scatter utilities ===
      std::vector<std::array<double, 2>> global_coords;     // gathered unique coords
      std::vector<double>                global_field;      // gathered corresponding values
      std::vector<std::vector<int>>      elem2glob;         // map elem-local -> global ids

    private:
      // Coarse geometry (Q2 nodes)
      std::array<double, 9> _X{};
      std::array<double, 9> _Y{};

      // Ensure all leaves satisfy min depth floor (by refining only)
      void ensure_min_depth() {
        if (_minDepth == 0) return;
        bool changed = true;
        while (changed) {
          changed = false;
          const auto& L = leaves();
          for (u32 leaf : L) {
            if (_nodes[leaf].level < _minDepth) {
              if (refine(leaf)) changed = true;
            }
          }
        }
      }

      // Allocate a new node id (reuse from freelist if available)
      u32 alloc_node() {
        u32 i;
        if (!_free.empty()) {
          i = _free.back();
          _free.pop_back();
        }
        else {
          i = (u32)_nodes.size();
          _nodes.emplace_back();
          _alive.push_back(0);
        }
        _alive[i] = 1;
        _nodes[i] = QuadNode();
        return i;
      }

      // Free a node id (push to freelist), mark leaves dirty
      void free_node(u32 i) {
        _alive[i] = 0;
        _free.push_back(i);
        _leaf_dirty = true;
      }

      // Bounds helper for any node
      inline void bounds_of(const QuadNode & n, double & xi0, double & eta0, double & xi1, double & eta1) const {
        u32 ix, iy;
        deinterleave2(n.morton, ix, iy);
        const double N = double(1u << n.level);
        const double dx = 2.0 / N, dy = 2.0 / N;
        xi0  = -1.0 + ix * dx;
        eta0 = -1.0 + iy * dy;
        xi1  = xi0 + dx;
        eta1 = eta0 + dy;
      }

      // Return compact-leaf position for node id, or npos32 if not a leaf
      inline u32 leaf_position(u32 leaf) const {
        (void)leaves();
        return (leaf < _leaf_pos.size()) ? _leaf_pos[leaf] : npos32;
      }

      // --- Internal data ---
      u32  _root{npos32};
      u32  _maxDepth;
      u32  _minDepth{0};
      bool _allowCoarsenBelowMinDepth{true};
      std::vector<QuadNode> _nodes;
      std::vector<uint8_t>  _alive;
      std::vector<u32>      _free;
      mutable std::vector<u32> _leaves;
      mutable bool _leaf_dirty{true};
      std::vector<Field> _fields;
      mutable std::vector<u32> _leaf_pos;
      bool _geom_ready{false};
  };

} // namespace fem


