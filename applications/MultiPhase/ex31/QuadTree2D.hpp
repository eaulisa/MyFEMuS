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


namespace fem {

  using u32 = uint32_t;
  using u64 = uint64_t;
  constexpr u32 npos32 = std::numeric_limits<u32>::max();

  enum class Basis : uint8_t { Q1_Quad4 = 0, Serendipity8 = 1, Q2_Quad9 = 2 };

//----------------------------------------
// Quad9 shapes (parent space [-1,1]^2)
//----------------------------------------
  struct Quad9Shape {
    // 1D Q2 Lagrange at xi: nodes {-1,0,1} -> [L0,L1,L2] and derivs
    static inline void q2_1d(double s, double& L0, double& L1, double& L2,
                             double& dL0, double& dL1, double& dL2) {
      // Nodes at s = {-1, 0, 1}
      L0 = 0.5 * s * (s - 1.0);   // 0.5*(s^2 - s)
      L1 = 1.0 - s * s;           // 1 - s^2
      L2 = 0.5 * s * (s + 1.0);   // 0.5*(s^2 + s)
      dL0 = s - 0.5;
      dL1 = -2.0 * s;
      dL2 = s + 0.5;
    }


    // N[0..8] in the standard 2D tensor product ordering:
    // corners: (-1,-1)=0, (1,-1)=1, (1,1)=2, (-1,1)=3
    // mid-edges: (0,-1)=4, (1,0)=5, (0,1)=6, (-1,0)=7
    // center: (0,0)=8
    static inline void N(double xi, double eta, double N9[9]) {
      double Lx[3], Ly[3], dx[3], dy[3];
      q2_1d(xi,  Lx[0], Lx[1], Lx[2], dx[0], dx[1], dx[2]);
      q2_1d(eta, Ly[0], Ly[1], Ly[2], dy[0], dy[1], dy[2]);
      // (dx,dy) unused here; silence -Wall if desired:
      (void)dx;
      (void)dy;

      // corners
      N9[0] = Lx[0] * Ly[0];
      N9[1] = Lx[2] * Ly[0];
      N9[2] = Lx[2] * Ly[2];
      N9[3] = Lx[0] * Ly[2];
      // mids
      N9[4] = Lx[1] * Ly[0];
      N9[5] = Lx[2] * Ly[1];
      N9[6] = Lx[1] * Ly[2];
      N9[7] = Lx[0] * Ly[1];
      // center
      N9[8] = Lx[1] * Ly[1];
    }


    static inline void dN(double xi, double eta, double dN_dxi[9], double dN_deta[9]) {
      double Lx[3], Ly[3], dx[3], dy[3];
      q2_1d(xi,  Lx[0], Lx[1], Lx[2], dx[0], dx[1], dx[2]);
      q2_1d(eta, Ly[0], Ly[1], Ly[2], dy[0], dy[1], dy[2]);
      // corners
      dN_dxi[0]  = dx[0] * Ly[0];
      dN_deta[0] = Lx[0] * dy[0];
      dN_dxi[1]  = dx[2] * Ly[0];
      dN_deta[1] = Lx[2] * dy[0];
      dN_dxi[2]  = dx[2] * Ly[2];
      dN_deta[2] = Lx[2] * dy[2];
      dN_dxi[3]  = dx[0] * Ly[2];
      dN_deta[3] = Lx[0] * dy[2];
      // mids
      dN_dxi[4]  = dx[1] * Ly[0];
      dN_deta[4] = Lx[1] * dy[0];
      dN_dxi[5]  = dx[2] * Ly[1];
      dN_deta[5] = Lx[2] * dy[1];
      dN_dxi[6]  = dx[1] * Ly[2];
      dN_deta[6] = Lx[1] * dy[2];
      dN_dxi[7]  = dx[0] * Ly[1];
      dN_deta[7] = Lx[0] * dy[1];
      // center
      dN_dxi[8]  = dx[1] * Ly[1];
      dN_deta[8] = Lx[1] * dy[1];
    }
  };

//----------------------------------------
// Q1 (Quad4) and Serendipity8 shapes on [-1,1]^2
//----------------------------------------
  struct Shapes {
    static inline void Q1(double xi, double eta, double N4[4]) {
      const double a = 0.25 * (1 - xi), b = 0.25 * (1 + xi);
      N4[0] = a * (1 - eta); // (-1,-1)
      N4[1] = b * (1 - eta); // ( 1,-1)
      N4[2] = b * (1 + eta); // ( 1, 1)
      N4[3] = a * (1 + eta); // (-1, 1)
    }

    static inline void Serendipity8(double xi, double eta, double N8[8]) {
      const double xm = xi, em = eta; // to match common notation
      // Standard serendipity 8-node shape functions
      N8[0] = 0.25 * (1 - xm) * (1 - em) * (-xm - em - 1);
      N8[1] = 0.25 * (1 + xm) * (1 - em) * (xm - em - 1);
      N8[2] = 0.25 * (1 + xm) * (1 + em) * (xm + em - 1);
      N8[3] = 0.25 * (1 - xm) * (1 + em) * (-xm + em - 1);
      N8[4] = 0.5 * (1 - xm * xm) * (1 - em); // mid-bottom
      N8[5] = 0.5 * (1 + xm) * (1 - em * em); // mid-right
      N8[6] = 0.5 * (1 - xm * xm) * (1 + em); // mid-top
      N8[7] = 0.5 * (1 - xm) * (1 - em * em); // mid-left
    }
  };

//----------------------------------------
// Morton helpers (2D)
//----------------------------------------
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
    u64 morton{0};   // cell origin index in grid at this level
    uint8_t flags{0};
    inline bool is_leaf() const {
      return child[0] == npos32;
    }
  };

//----------------------------------------
// Field storage (per-leaf coefficients)
//----------------------------------------
  struct Field {
    Basis basis{Basis::Q1_Quad4};
    u32   dofs_per_cell{4}; // 4, 8, or 9
    // Flat coefficient array laid out as blocks of dofs_per_cell for each leaf in leaf order.
    // You manage the values; class provides indexing helpers.
    std::vector<double> coeffs;
  };

//----------------------------------------
// Quadtree living on parent [-1,1]^2, with Quad9 geometry map
//----------------------------------------
  class QuadTree2D {
    public:
      QuadTree2D(u32 maxDepth = 28)
        : _maxDepth(maxDepth) {
        _nodes.reserve(1024);
        _alive.reserve(1024);
        _free.reserve(256);
        _root = alloc_node();
        _nodes[_root].level = 0;
        _nodes[_root].morton = 0;
      }


      // Overload: set minDepth (soft floor)
      QuadTree2D(u32 maxDepth, u32 minDepth)
        : _maxDepth(maxDepth), _minDepth(std::min(minDepth, maxDepth)) {
        _nodes.reserve(1024);
        _alive.reserve(1024);
        _free.reserve(256);
        _root = alloc_node();
        _nodes[_root].level = 0;
        _nodes[_root].morton = 0;
      }

      // Min-depth controls
      void set_min_depth(u32 d) {
        _minDepth = std::min(d, _maxDepth);
      }
      u32  min_depth() const {
        return _minDepth;
      }

      void set_allow_coarsen_below_min(bool v) {
        _allowCoarsenBelowMinDepth = v;
      }
      bool allow_coarsen_below_min() const {
        return _allowCoarsenBelowMinDepth;
      }
// ---- geometry setup (Quad9 physical coords) ----
      // Provide x[9], y[9] in the standard Quad9 node order (see Quad9Shape::N comment).
      void set_physical_quad9(const double x9[9], const double y9[9]) {
        for (int i = 0; i < 9; ++i) {
          _X[i] = x9[i];
          _Y[i] = y9[i];
        }
      }

      // ---- refinement ----
      bool refine(u32 i) {
        // Basic guards
        if (i == npos32 || !_alive[i]) return false;
        if (!_nodes[i].is_leaf() || _nodes[i].level >= _maxDepth) return false;

        // Snapshot parent info BEFORE allocations (avoid invalidated refs)
        const u32 parent_level  = _nodes[i].level;
        const u64 parent_morton = _nodes[i].morton;

        // Allocate and initialize 4 children
        for (int k = 0; k < 4; ++k) {
          const u32 c = alloc_node();              // may reallocate _nodes
          _nodes[i].child[k] = c;                  // use index, not a stale reference
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


      bool coarsen(u32 i) {
        // Basic guards
        if (i == npos32 || !_alive[i]) return false;
        if (_nodes[i].is_leaf()) return false;

        // If parent is at/below floor and dropping below is disabled, block
        if (_nodes[i].level <= _minDepth && !_allowCoarsenBelowMinDepth) return false;

        // Snapshot children (avoid using a reference across mutations)
        std::array<u32, 4> ch = _nodes[i].child;

        // Validate: all 4 children must exist, be alive, and be leaves
        for (int k = 0; k < 4; ++k) {
          const u32 c = ch[k];
          if (c == npos32) return false;
          if (!_alive[c])  return false;
          if (!_nodes[c].is_leaf()) return false;
          // Optional: ensure parent matches
          if (_nodes[c].parent != i) return false;
        }

        // Free children first
        for (int k = 0; k < 4; ++k) {
          free_node(ch[k]);   // marks _alive[c]=0 and pushes to freelist; no reallocation
        }

        // Now clear the parent's child pointers
        _nodes[i].child = {npos32, npos32, npos32, npos32};

        _leaf_dirty = true;
        return true;
      }


      // Build/refresh the compact leaf list (indices into _nodes)

      const std::vector<u32>& leaves() const {
        if (_leaf_dirty) {
          _leaves.clear();
          _leaves.reserve(_nodes.size());
          for (u32 i = 0; i < (u32)_nodes.size(); ++i)
            if (_alive[i] && _nodes[i].is_leaf())
              _leaves.push_back(i);
          _leaf_dirty = false;
        }
        return _leaves;
      }


// One refinement "pass": evaluate the predicate on the CURRENT leaves only,
// collect those that need refinement, refine them, refresh leaf list.
// Returns the number of leaves refined in this pass.
      template<class Pred>
      std::size_t refine_pass(Pred&& should_refine, Basis probe_basis = Basis::Q2_Quad9) {
        const auto& L = leaves();               // snapshot of current leaves
        std::vector<u32> to_refine;
        to_refine.reserve(L.size());

        for (u32 k = 0; k < (u32)L.size(); ++k) {
          const u32 leaf = L[k];

          std::vector<std::array<double, 2>> pts_xi, pts_xy;
          std::vector<std::array<double, 9>> Nvals;
          leaf_parent_nodes(probe_basis, leaf, pts_xi);
          leaf_physical_nodes(probe_basis, leaf, pts_xy);
          quad9_shapes_at(pts_xi, Nvals);

          if (should_refine(leaf, pts_xi, pts_xy, Nvals)) {
            // respect maxDepth inside refine(); it returns false if at max
            to_refine.push_back(leaf);
          }
        }

        std::size_t refined = 0;
        for (u32 leaf : to_refine) if (refine(leaf)) ++refined;

        (void)leaves(); // refresh compact leaf list
        return refined;
      }




      // Get leaf bounding box in parent coords (xi,eta)
      inline void leaf_bounds(u32 i, double& xi0, double& eta0, double& xi1, double& eta1) const {
        const QuadNode& n = _nodes[i];
        u32 ix, iy;
        deinterleave2(n.morton, ix, iy);
        const double N = double(1u << n.level);
        const double dx = 2.0 / N; // width in [-1,1]
        const double dy = 2.0 / N;
        xi0  = -1.0 + ix * dx;
        eta0 = -1.0 + iy * dy;
        xi1  = xi0 + dx;
        eta1 = eta0 + dy;
      }

      // Locate leaf containing (xi,eta) in parent space (assumes inside [-1,1]^2)
      u32 locate_parent(double xi, double eta) const {
        if (xi < -1.0 || xi > 1.0 || eta < -1.0 || eta > 1.0) return npos32;
        u32 i = _root;
        while (true) {
          const auto& n = _nodes[i];
          if (n.is_leaf()) return i;
          // compute child quadrant by midpoint test in parent coords
          double xi0, eta0, xi1, eta1;
          bounds_of(n, xi0, eta0, xi1, eta1);
          const double cx = 0.5 * (xi0 + xi1);
          const double cy = 0.5 * (eta0 + eta1);
          const int qx = (xi >= cx) ? 1 : 0;
          const int qy = (eta >= cy) ? 1 : 0;
          const int q  = (qy << 1) | qx; // 0:LL,1:LR,2:UL,3:UR
          const u32 c = n.child[q];
          if (c == npos32) return i;
          i = c;
        }
      }

      // Inverse map physical (x,y) -> (xi,eta) using Quad9 isoparametric map (Newton)
      // Returns false if fails to converge or Jacobian singular. On success sets (xi,eta).
      bool inverse_map_quad9(double x, double y, double& xi, double& eta,
                             double tol = 1e-12, int maxit = 25) const {
        // Start from center or a quick affine guess:
        xi  = 0.0;
        eta = 0.0;

        for (int it = 0; it < maxit; ++it) {
          double N9[9], dNdxi[9], dNdeta[9];
          Quad9Shape::N(xi, eta, N9);
          Quad9Shape::dN(xi, eta, dNdxi, dNdeta);

          // xhat(xi,eta), yhat(xi,eta)
          double Xh = 0.0, Yh = 0.0, dXdxi = 0.0, dXdeta = 0.0, dYdxi = 0.0, dYdeta = 0.0;
          for (int a = 0; a < 9; ++a) {
            Xh += N9[a] * _X[a];
            Yh += N9[a] * _Y[a];
            dXdxi  += dNdxi[a] * _X[a];
            dXdeta += dNdeta[a] * _X[a];
            dYdxi  += dNdxi[a] * _Y[a];
            dYdeta += dNdeta[a] * _Y[a];
          }

          // residual
          const double rx = Xh - x;
          const double ry = Yh - y;
          const double norm = std::sqrt(rx * rx + ry * ry);
          if (norm < tol) return true;

          // Jacobian and its inverse
          const double J00 = dXdxi, J01 = dXdeta;
          const double J10 = dYdxi, J11 = dYdeta;
          const double detJ = J00 * J11 - J01 * J10;
          if (std::abs(detJ) <= 1e-30) return false;

          // Newton step: [d_xi; d_eta] = J^{-1} * [rx; ry]
          const double dxi  = (J11 * rx - J01 * ry) / detJ;
          const double deta = (-J10 * rx + J00 * ry) / detJ;

          xi  -= dxi;
          eta -= deta;

          // Optional damping for robustness on highly distorted quads
          if (std::abs(dxi) + std::abs(deta) < tol * 1e-6) return true;
        }
        return false;
      }

      // Map global parent (xi,eta) into a leaf-local reference \hat{s},\hat{t} in [-1,1]^2
      static inline void to_leaf_ref(double xi, double eta,
                                     double xi0, double eta0, double xi1, double eta1,
                                     double& shat, double& that) {
        const double cx = 0.5 * (xi0 + xi1), hx = 0.5 * (xi1 - xi0);
        const double cy = 0.5 * (eta0 + eta1), hy = 0.5 * (eta1 - eta0);
        shat = (xi  - cx) / hx; // in [-1,1]
        that = (eta - cy) / hy; // in [-1,1]
      }

      // -------- Field API --------
      // Add a field (returns field id). You must later size coeffs:
      //   coeffs.resize(num_leaves * dofs_per_cell)
      u32 add_field(Basis b) {
        Field f;
        f.basis = b;
        f.dofs_per_cell = (b == Basis::Q1_Quad4) ? 4 : (b == Basis::Serendipity8 ? 8 : 9);
        _fields.push_back(std::move(f));
        return (u32)(_fields.size() - 1);
      }

      Field& field(u32 fid) {
        return _fields[fid];
      }
      const Field& field(u32 fid) const {
        return _fields[fid];
      }

      // Index helper: returns pointer to leaf's coefficient block for field fid
      double* leaf_coeff_ptr(u32 fid, u32 leaf_pos) {
        Field& f = _fields[fid];
        return f.coeffs.data() + size_t(leaf_pos) * f.dofs_per_cell;
      }
      const double* leaf_coeff_ptr(u32 fid, u32 leaf_pos) const {
        const Field& f = _fields[fid];
        return f.coeffs.data() + size_t(leaf_pos) * f.dofs_per_cell;
      }

      // Evaluate field fid at physical (x,y). Returns false if inverse fails or outside.
      bool evaluate_physical(u32 fid, double x, double y, double& value) {
        double xi, eta;
        if (!inverse_map_quad9(x, y, xi, eta)) return false;
        const u32 leaf = locate_parent(xi, eta);
        if (leaf == npos32) return false;

        // Get leaf bounds and map to leaf reference
        double xi0, eta0, xi1, eta1;
        leaf_bounds(leaf, xi0, eta0, xi1, eta1);
        double shat, that;
        to_leaf_ref(xi, eta, xi0, eta0, xi1, eta1, shat, that);

        // Locate leaf position in compact list
        const auto& L = leaves();
        // For speed, cache a map; here we do linear search (ok for prototype).
        u32 leaf_pos = npos32;
        for (u32 k = 0; k < (u32)L.size(); ++k) if (L[k] == leaf) {
            leaf_pos = k;
            break;
          }
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

      // Utility: number of leaves
      u32 leaf_count() {
        return (u32)leaves().size();
      }

      // Rebuild field storage sizes to match current leaves (does not touch values)
      void resize_fields_to_leaves() {
        const auto nL = leaf_count();
        for (auto& f : _fields) {
          f.coeffs.resize(size_t(nL) * f.dofs_per_cell);
        }
      }

      // Expose nodes/leaves if you want to drive refinement externally
      const std::vector<u32>& leaf_indices() {
        return leaves();
      }


// Return the basis node positions (xi,eta) for a given leaf subcell
// out_pts.size() will be 4, 8, or 9 depending on 'basis'.
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
        if (basis == Basis::Q1_Quad4) {
          // corners in [-1,1]^2 of the leaf
          static const double S[4][2] = {{-1, -1}, {+1, -1}, {+1, +1}, {-1, +1}};
          out_pts.reserve(4);
          for (int i = 0; i < 4; ++i) out_pts.push_back(toParent(S[i][0], S[i][1]));
        }
        else if (basis == Basis::Serendipity8) {
          static const double S[8][2] = {{-1, -1}, {+1, -1}, {+1, +1}, {-1, +1}, {0, -1}, {+1, 0}, {0, +1}, {-1, 0}};
          out_pts.reserve(8);
          for (int i = 0; i < 8; ++i) out_pts.push_back(toParent(S[i][0], S[i][1]));
        }
        else {   // Q2_Quad9
          static const double S[9][2] = {{-1, -1}, {+1, -1}, {+1, +1}, {-1, +1}, {0, -1}, {+1, 0}, {0, +1}, {-1, 0}, {0, 0}};
          out_pts.reserve(9);
          for (int i = 0; i < 9; ++i) out_pts.push_back(toParent(S[i][0], S[i][1]));
        }
      }

// Map given leaf basis nodes to physical (x,y) via Quad9 isoparametric map
      void leaf_physical_nodes(Basis basis, u32 leaf,
                               std::vector<std::array<double, 2>>& out_xy) const {
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

// Evaluate Quad9 shape values at arbitrary parent points
// out_N: length = points.size(), each entry is array<double,9>
      void quad9_shapes_at(const std::vector<std::array<double, 2>>& points,
                           std::vector<std::array<double, 9>>& out_N) const {
        out_N.resize(points.size());
        for (size_t i = 0; i < points.size(); ++i) {
          double N9[9];
          Quad9Shape::N(points[i][0], points[i][1], N9);
          for (int a = 0; a < 9; ++a) out_N[i][a] = N9[a];
        }
      }

// Predicate-based refinement pass over current leaves.
// The predicate receives: leaf index, parent-basis node coords (xi,eta),
// physical node coords (x,y), and Quad9 shape values at those nodes.
// Return true to refine the leaf.
      template<class Pred>
      void refine_where(Pred&& should_refine, Basis probe_basis = Basis::Q2_Quad9) {
        const auto& L = leaves(); // compact leaf list
        std::vector<u32> to_refine;
        to_refine.reserve(L.size());
        for (u32 k = 0; k < (u32)L.size(); ++k) {
          u32 leaf = L[k];
          std::vector<std::array<double, 2>> pts_xi, pts_xy;
          std::vector<std::array<double, 9>> Nvals;
          leaf_parent_nodes(probe_basis, leaf, pts_xi);
          leaf_physical_nodes(probe_basis, leaf, pts_xy);
          quad9_shapes_at(pts_xi, Nvals);
          if (should_refine(leaf, pts_xi, pts_xy, Nvals)) {
            to_refine.push_back(leaf);
          }
        }
        // apply
        for (u32 leaf : to_refine) refine(leaf);
        (void)leaves(); // refresh leaf list
      }



      // Coarsen-pass: visit parents whose 4 children are leaves.
      // Predicate signature:
      //   bool(u32 parent, u32 level,
      //        const std::vector<std::array<double,2>>& parent_pts_xi,
      //        const std::vector<std::array<double,2>>& parent_pts_xy,
      //        const std::vector<std::array<double,9>>& parent_Nvals)
      template<class Pred>
      std::size_t coarsen_pass(Pred&& should_coarsen,
                               Basis probe_basis = Basis::Q2_Quad9) {
        std::vector<u32> to_coarsen;
        to_coarsen.reserve(_nodes.size() / 4 + 1);

        for (u32 i = 0; i < (u32)_nodes.size(); ++i) {
          if (!_alive[i]) continue;
          const QuadNode& n = _nodes[i];
          if (n.is_leaf()) continue;

          // must have 4 alive leaf-children
          bool ok = true;
          for (int k = 0; k < 4; ++k) {
            u32 c = n.child[k];
            if (c == npos32 || !_alive[c] || !_nodes[c].is_leaf() || _nodes[c].parent != i) {
              ok = false;
              break;
            }
          }
          if (!ok) continue;

          // Build probe data on the PARENT cell
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
        (void)leaves(); // refresh cache
        return done;
      }

      // Multilevel adaptivity: run up to max_passes passes until no new refinements.
      template<class Pred>
      std::size_t adapt_refine_until(Pred&& should_refine,
                                     Basis probe_basis = Basis::Q2_Quad9,
                                     u32 max_passes = 10) {


        ensure_min_depth();
        std::size_t total_refined = 0;
        for (u32 pass = 0; pass < max_passes; ++pass) {
          const std::size_t r = refine_pass(should_refine, probe_basis);
          total_refined += r;
          if (r == 0) break; // converged
        }
        return total_refined;
      }


      // Full adaptivity cycle: coarsen (should_coarsen) then refine (should_refine), iterate.
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

      // Evaluate a field on a known leaf using leaf-local (shat,that) in [-1,1]^2
      bool evaluate_on_leaf(u32 fid, u32 leaf, double shat, double that, double& value) const {
        // Find leaf position in the compact leaf list
        // (For speed, you can maintain a leaf->position map. This is simple & safe.)
        const auto& L = leaves();
        u32 leaf_pos = npos32;
        for (u32 k = 0; k < (u32)L.size(); ++k) if (L[k] == leaf) {
            leaf_pos = k;
            break;
          }
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

      bool write_vtu(const std::string& filename, u32 fid, const std::string& name,
                     bool cell_centered = false) const {
        const auto& L = leaves();
        const size_t numCells = L.size();
        if (numCells == 0) return false;

        struct Pt {
          double x, y, z;
        };
        std::vector<Pt> points;
        points.reserve(numCells * 4);

        std::vector<int> connectivity;
        connectivity.reserve(numCells * 4);
        std::vector<int> offsets;
        offsets.reserve(numCells);
        std::vector<unsigned char> types;
        types.reserve(numCells);

        std::vector<double> pointData;
        pointData.reserve(numCells * 4);
        std::vector<double> cellData;
        cellData.reserve(numCells);

        // NEW: per-cell refinement levels
        std::vector<int> levels;
        levels.reserve(numCells);

        static const double ST[4][2] = { {-1, -1}, {+1, -1}, {+1, +1}, {-1, +1} };

        const Field& fld = _fields[fid];
        const int ndof = (int)fld.dofs_per_cell;
        const double* coeff_base = fld.coeffs.data();

        for (size_t k = 0; k < numCells; ++k) {
          const u32 leaf = L[k];
          levels.push_back((int)_nodes[leaf].level); // <<< record level

          const double* c = coeff_base + size_t(k) * ndof;

          double xi0, eta0, xi1, eta1;
          leaf_bounds(leaf, xi0, eta0, xi1, eta1);
          const double cx = 0.5 * (xi0 + xi1), hx = 0.5 * (xi1 - xi0);
          const double cy = 0.5 * (eta0 + eta1), hy = 0.5 * (eta1 - eta0);

          for (int q = 0; q < 4; ++q) {
            const double sh = ST[q][0], th = ST[q][1];
            const double xi  = cx + sh * hx;
            const double eta = cy + th * hy;

            double N9[9];
            Quad9Shape::N(xi, eta, N9);
            double X = 0.0, Y = 0.0;
            for (int a = 0; a < 9; ++a) {
              X += N9[a] * _X[a];
              Y += N9[a] * _Y[a];
            }
            points.push_back({X, Y, 0.0});

            if (!cell_centered) {
              double v = 0.0;
              switch (fld.basis) {
                case Basis::Q1_Quad4: {
                  double N4[4];
                  Shapes::Q1(sh, th, N4);
                  v = N4[0] * c[0] + N4[1] * c[1] + N4[2] * c[2] + N4[3] * c[3];
                }
                break;
                case Basis::Serendipity8: {
                  double N8[8];
                  Shapes::Serendipity8(sh, th, N8);
                  for (int a = 0; a < 8; ++a) v += N8[a] * c[a];
                }
                break;
                case Basis::Q2_Quad9: {
                  double Nq[9];
                  Quad9Shape::N(sh, th, Nq);
                  for (int a = 0; a < 9; ++a) v += Nq[a] * c[a];
                }
                break;
              }
              pointData.push_back(v);
            }
          }

          const int base = (int)(k * 4);
          connectivity.push_back(base + 0);
          connectivity.push_back(base + 1);
          connectivity.push_back(base + 2);
          connectivity.push_back(base + 3);
          offsets.push_back(base + 4);
          types.push_back(9);

          if (cell_centered) {
            double v = 0.0;
            switch (fld.basis) {
              case Basis::Q1_Quad4: {
                double N4[4];
                Shapes::Q1(0.0, 0.0, N4);
                v = N4[0] * c[0] + N4[1] * c[1] + N4[2] * c[2] + N4[3] * c[3];
              }
              break;
              case Basis::Serendipity8: {
                double N8[8];
                Shapes::Serendipity8(0.0, 0.0, N8);
                for (int a = 0; a < 8; ++a) v += N8[a] * c[a];
              }
              break;
              case Basis::Q2_Quad9: {
                double Nq[9];
                Quad9Shape::N(0.0, 0.0, Nq);
                for (int a = 0; a < 9; ++a) v += Nq[a] * c[a];
              }
              break;
            }
            cellData.push_back(v);
          }
        }

        // --- ASCII VTU ---
        std::ofstream os(filename);
        if (!os) return false;
        os << std::setprecision(16);
        os << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        os << "  <UnstructuredGrid>\n";
        os << "    <Piece NumberOfPoints=\"" << points.size()
           << "\" NumberOfCells=\"" << numCells << "\">\n";

        // Points
        os << "      <Points>\n";
        os << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (const auto& p : points) os << "          " << p.x << " " << p.y << " " << p.z << "\n";
        os << "        </DataArray>\n";
        os << "      </Points>\n";

        // Cells
        os << "      <Cells>\n";
        os << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
        for (size_t i = 0; i < connectivity.size(); i += 4)
          os << "          " << connectivity[i] << " " << connectivity[i + 1]
             << " " << connectivity[i + 2] << " " << connectivity[i + 3] << "\n";
        os << "        </DataArray>\n";
        os << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
        for (int off : offsets) os << "          " << off << "\n";
        os << "        </DataArray>\n";
        os << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
        for (unsigned char t : types) os << "          " << (int)t << "\n";
        os << "        </DataArray>\n";
        os << "      </Cells>\n";

        // Data
        // Always write CellData "level" (Int32)
        os << "      <CellData Scalars=\"" << name << "\">\n";
        os << "        <DataArray type=\"Int32\" Name=\"level\" format=\"ascii\">\n";
        for (int lv : levels) os << "          " << lv << "\n";
        os << "        </DataArray>\n";

        if (cell_centered) {
          // cell-centered scalar in the same CellData section
          os << "        <DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
          for (double v : cellData) os << "          " << v << "\n";
          os << "        </DataArray>\n";
          os << "      </CellData>\n";
          // no PointData
        }
        else {
          os << "      </CellData>\n";
          // point-centered scalar
          os << "      <PointData Scalars=\"" << name << "\">\n";
          os << "        <DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
          for (double v : pointData) os << "          " << v << "\n";
          os << "        </DataArray>\n";
          os << "      </PointData>\n";
        }

        os << "    </Piece>\n";
        os << "  </UnstructuredGrid>\n";
        os << "</VTKFile>\n";
        return true;
      }



    private:
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

      // Physical geometry (Quad9 nodes)
      double _X[9] {}, _Y[9] {};

      // Ensure all leaves have level >= _minDepth by refining only
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
      void free_node(u32 i) {
        _alive[i] = 0;
        _free.push_back(i);
      }

      inline void bounds_of(const QuadNode& n, double& xi0, double& eta0, double& xi1, double& eta1) const {
        u32 ix, iy;
        deinterleave2(n.morton, ix, iy);
        const double N = double(1u << n.level);
        const double dx = 2.0 / N, dy = 2.0 / N;
        xi0  = -1.0 + ix * dx;
        eta0 = -1.0 + iy * dy;
        xi1  = xi0 + dx;
        eta1 = eta0 + dy;
      }


    public:


      // === Unified node gathering ===
      std::vector<std::array<double, 2>> global_coords;
      std::vector<double> global_field;
      std::vector<std::vector<int>> elem2glob;

      static inline bool same_point(const std::array<double, 2>& a,
                                    const std::array<double, 2>& b,
                                    double tol = 1e-12) {
        return (std::fabs(a[0] - b[0]) < tol && std::fabs(a[1] - b[1]) < tol);
      }

      static int basis_nodes(Basis b) {
        switch(b) {
          case Basis::Q1_Quad4:
            return 4;
          case Basis::Serendipity8:
            return 8;
          case Basis::Q2_Quad9:
            return 9;
        }
        return 0;
      }

// Gather unique coordinates and associated field values
      void gather_field(u32 fid) {
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
          for (int i = 0; i < nNodes; ++i) {
            flat.push_back({pts_xy[i], coeffs[i], e, i});
          }
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
          if (elem2glob[flat[k].elem].empty())
            elem2glob[flat[k].elem].resize(nNodes);
          elem2glob[flat[k].elem][flat[k].local] = gid;
        }
      }

// Scatter modified global_field back into the field coeffs
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


      // Refine tree so that all given physical points are inside leaves
// at least as deep as maxDepthTarget (or until global _maxDepth).
      void refine_to_contain_points(const std::vector<std::array<double, 2>>& pts,
                                    u32 maxDepthTarget) {
        for (auto& p : pts) {
          double xi, eta;
          if (!inverse_map_quad9(p[0], p[1], xi, eta)) continue; // skip if outside domain

          u32 leaf = locate_parent(xi, eta);
          while (leaf != npos32 && _nodes[leaf].level < maxDepthTarget) {
            if (!refine(leaf)) break;   // refine returns false if not refinable
            leaf = locate_parent(xi, eta); // find the new child containing p
          }
        }
        (void)leaves(); // refresh cached leaf list
      }

// Extract node coordinates of all leaves whose refinement level is in [minLevel, maxLevel].
// Each entry in the result corresponds to one leaf and contains its node coordinates in physical space.
      std::vector<std::vector<std::array<double, 2>>>
      extract_leaf_nodes_in_level_range(u32 minLevel, u32 maxLevel, Basis basis) const {
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

      // Extract all node coordinates from leaves whose level is in [minLevel, maxLevel].
// Returns a flat vector of coordinates (duplicates possible).
      std::vector<std::array<double, 2>>
      extract_node_coords_in_level_range(u32 minLevel, u32 maxLevel, Basis basis) const {
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


      // Reset the tree to a single root cell (remove all refinements)
void reset() {
    _nodes.clear();
    _alive.clear();
    _free.clear();
    _leaves.clear();
    _leaf_dirty = true;

    _root = alloc_node();
    _nodes[_root].level = 0;
    _nodes[_root].morton = 0;

    // Keep existing fields but resize coefficients to 1 cell
    for (auto& f : _fields) {
        f.coeffs.clear();
        f.coeffs.resize(f.dofs_per_cell);
    }
}





  };

} // namespace fem


















