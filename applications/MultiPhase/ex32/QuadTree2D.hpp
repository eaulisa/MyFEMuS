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

  enum class Basis : uint8_t { Q1_Quad4 = 0, Serendipity8 = 1, Q2_Quad9 = 2 };


//--------------------------------------------
// VTK cell-type helper
//--------------------------------------------
  inline unsigned char vtk_cell_type(Basis b) {
    switch (b) {
    case Basis::Q1_Quad4:
      return 9;   // VTK_QUAD
    case Basis::Serendipity8:
      return 23;  // VTK_QUADRATIC_QUAD
    case Basis::Q2_Quad9:
      return 28;  // VTK_BIQUADRATIC_QUAD
    }
    return 9; // sensible default
  }



//--------------------------------------------
// Morton helpers
//--------------------------------------------
  inline u64 interleave2(u32 x, u32 y) {
    u64 xx = x, yy = y;
    xx = (xx | (xx << 16)) & 0x0000FFFF0000FFFFULL;
    xx = (xx | (xx << 8))  & 0x00FF00FF00FF00FFULL;
    xx = (xx | (xx << 4))  & 0x0F0F0F0F0F0F0F0FULL;
    xx = (xx | (xx << 2))  & 0x3333333333333333ULL;
    xx = (xx | (xx << 1))  & 0x5555555555555555ULL;

    yy = (yy | (yy << 16)) & 0x0000FFFF0000FFFFULL;
    yy = (yy | (yy << 8))  & 0x00FF00FF00FF00FFULL;
    yy = (yy | (yy << 4))  & 0x0F0F0F0F0F0F0F0FULL;
    yy = (yy | (yy << 2))  & 0x3333333333333333ULL;
    yy = (yy | (yy << 1))  & 0x5555555555555555ULL;
    return (yy << 1) | xx;
  }

  inline void deinterleave2(u64 m, u32 &x, u32 &y) {
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

//--------------------------------------------
// Shapes / geometry
//--------------------------------------------
  struct Quad9Shape {
    // 1D quadratic Lagrange basis (values only)
    static inline void q2_1d_vals(double s, double& L0, double& L1, double& L2) {
      L0 = 0.5 * s * (s - 1.0);
      L1 = 1.0 - s * s;
      L2 = 0.5 * s * (s + 1.0);
    }

    // 1D quadratic Lagrange basis (derivatives only)
    static inline void q2_1d_derivs(double s, double& dL0, double& dL1, double& dL2) {
      dL0 = s - 0.5;
      dL1 = -2.0 * s;
      dL2 = s + 0.5;
    }

    // 2D Q2 shape values (no derivative work here)
    static inline void N(double xi, double eta, double N9[9]) {
      double Lx0, Lx1, Lx2;
      double Ly0, Ly1, Ly2;
      q2_1d_vals(xi,  Lx0, Lx1, Lx2);
      q2_1d_vals(eta, Ly0, Ly1, Ly2);

      N9[0] = Lx0 * Ly0;  // (-1,-1)
      N9[1] = Lx2 * Ly0;  // (+1,-1)
      N9[2] = Lx2 * Ly2;  // (+1,+1)
      N9[3] = Lx0 * Ly2;  // (-1,+1)
      N9[4] = Lx1 * Ly0;  // ( 0,-1)
      N9[5] = Lx2 * Ly1;  // (+1, 0)
      N9[6] = Lx1 * Ly2;  // ( 0,+1)
      N9[7] = Lx0 * Ly1;  // (-1, 0)
      N9[8] = Lx1 * Ly1;  // ( 0, 0)
    }

    // 2D Q2 shape derivatives
    static inline void dN(double xi, double eta, double dN_dxi[9], double dN_deta[9]) {
      double Lx0, Lx1, Lx2, dx0, dx1, dx2;
      double Ly0, Ly1, Ly2, dy0, dy1, dy2;

      q2_1d_vals(xi,  Lx0, Lx1, Lx2);
      q2_1d_derivs(xi,  dx0, dx1, dx2);
      q2_1d_vals(eta, Ly0, Ly1, Ly2);
      q2_1d_derivs(eta, dy0, dy1, dy2);

      dN_dxi[0]  = dx0 * Ly0;
      dN_deta[0] = Lx0 * dy0;
      dN_dxi[1]  = dx2 * Ly0;
      dN_deta[1] = Lx2 * dy0;
      dN_dxi[2]  = dx2 * Ly2;
      dN_deta[2] = Lx2 * dy2;
      dN_dxi[3]  = dx0 * Ly2;
      dN_deta[3] = Lx0 * dy2;
      dN_dxi[4]  = dx1 * Ly0;
      dN_deta[4] = Lx1 * dy0;
      dN_dxi[5]  = dx2 * Ly1;
      dN_deta[5] = Lx2 * dy1;
      dN_dxi[6]  = dx1 * Ly2;
      dN_deta[6] = Lx1 * dy2;
      dN_dxi[7]  = dx0 * Ly1;
      dN_deta[7] = Lx0 * dy1;
      dN_dxi[8]  = dx1 * Ly1;
      dN_deta[8] = Lx1 * dy1;
    }

    // Optional: compute N and dN together without duplicating 1D work
    static inline void N_and_dN(double xi, double eta,
                                double N9[9], double dN_dxi[9], double dN_deta[9]) {
      double Lx0, Lx1, Lx2, dx0, dx1, dx2;
      double Ly0, Ly1, Ly2, dy0, dy1, dy2;

      q2_1d_vals(xi,  Lx0, Lx1, Lx2);
      q2_1d_derivs(xi,  dx0, dx1, dx2);
      q2_1d_vals(eta, Ly0, Ly1, Ly2);
      q2_1d_derivs(eta, dy0, dy1, dy2);

      // N
      N9[0] = Lx0 * Ly0;
      N9[1] = Lx2 * Ly0;
      N9[2] = Lx2 * Ly2;
      N9[3] = Lx0 * Ly2;
      N9[4] = Lx1 * Ly0;
      N9[5] = Lx2 * Ly1;
      N9[6] = Lx1 * Ly2;
      N9[7] = Lx0 * Ly1;
      N9[8] = Lx1 * Ly1;

      // dN
      dN_dxi[0]  = dx0 * Ly0;
      dN_deta[0] = Lx0 * dy0;
      dN_dxi[1]  = dx2 * Ly0;
      dN_deta[1] = Lx2 * dy0;
      dN_dxi[2]  = dx2 * Ly2;
      dN_deta[2] = Lx2 * dy2;
      dN_dxi[3]  = dx0 * Ly2;
      dN_deta[3] = Lx0 * dy2;
      dN_dxi[4]  = dx1 * Ly0;
      dN_deta[4] = Lx1 * dy0;
      dN_dxi[5]  = dx2 * Ly1;
      dN_deta[5] = Lx2 * dy1;
      dN_dxi[6]  = dx1 * Ly2;
      dN_deta[6] = Lx1 * dy2;
      dN_dxi[7]  = dx0 * Ly1;
      dN_deta[7] = Lx0 * dy1;
      dN_dxi[8]  = dx1 * Ly1;
      dN_deta[8] = Lx1 * dy1;
    }
  };

  struct Shapes {
    // ---------- Values (yours) ----------
    static inline void Q1(double xi, double eta, double N4[4]) {
      const double a = 0.25 * (1 - xi), b = 0.25 * (1 + xi);
      N4[0] = a * (1 - eta);
      N4[1] = b * (1 - eta);
      N4[2] = b * (1 + eta);
      N4[3] = a * (1 + eta);
    }
    static inline void Q8(double xi, double eta, double N8[8]) {
      const double xm = xi, em = eta;
      N8[0] = 0.25 * (1 - xm) * (1 - em) * (-xm - em - 1);
      N8[1] = 0.25 * (1 + xm) * (1 - em) * (xm - em - 1);
      N8[2] = 0.25 * (1 + xm) * (1 + em) * (xm + em - 1);
      N8[3] = 0.25 * (1 - xm) * (1 + em) * (-xm + em - 1);
      N8[4] = 0.5  * (1 - xm * xm) * (1 - em);
      N8[5] = 0.5  * (1 + xm) * (1 - em * em);
      N8[6] = 0.5  * (1 - xm * xm) * (1 + em);
      N8[7] = 0.5  * (1 - xm) * (1 - em * em);
    }

    // ---------- Derivatives ----------
    // Derivatives of Q1 (Quad4) shape functions in the standard node order:
// 0:(-1,-1), 1:(1,-1), 2:(1,1), 3:(-1,1)
    static inline void Q1_dN(double xi, double eta,
                             double* __restrict dNdxi,
                             double* __restrict dNdeta) {
      // Precompute (1 ± xi), (1 ± eta)
      const double sm = 1.0 - xi;
      const double sp = 1.0 + xi;
      const double tm = 1.0 - eta;
      const double tp = 1.0 + eta;

      // Common scale
      const double c = 0.25;

      // ∂N/∂xi
      dNdxi[0] = -c * tm;   // -0.25 * (1 - eta)
      dNdxi[1] =  c * tm;   //  0.25 * (1 - eta)
      dNdxi[2] =  c * tp;   //  0.25 * (1 + eta)
      dNdxi[3] = -c * tp;   // -0.25 * (1 + eta)

      // ∂N/∂eta
      dNdeta[0] = -c * sm;  // -0.25 * (1 - xi)
      dNdeta[1] = -c * sp;  // -0.25 * (1 + xi)
      dNdeta[2] =  c * sp;  //  0.25 * (1 + xi)
      dNdeta[3] =  c * sm;  //  0.25 * (1 - xi)
    }


    // Q8 (serendipity) derivatives in the same node order as Q8():
    // corners 0..3, then edge mids 4..7 (bottom,right,top,left)
    static inline void Q8_dN(double xi, double eta,
                             double dNdxi[8], double dNdeta[8]) {
      const double x = xi, e = eta;

      // corners
      dNdxi[0]  = 0.25 * (1 - e) * (2 * x + e);  // ∂N0/∂x
      dNdeta[0] = 0.25 * (1 - x) * (x + 2 * e);  // ∂N0/∂e

      dNdxi[1]  = 0.25 * (1 - e) * (2 * x - e);
      dNdeta[1] = 0.25 * (1 + x) * (-x + 2 * e);

      dNdxi[2]  = 0.25 * (1 + e) * (2 * x + e);
      dNdeta[2] = 0.25 * (1 + x) * (x + 2 * e);

      dNdxi[3]  = 0.25 * (1 + e) * (2 * x - e);
      dNdeta[3] = 0.25 * (1 - x) * (-x + 2 * e);

      // edges
      dNdxi[4]  =      -x * (1 - e);             // ∂/∂x of 0.5*(1-x^2)*(1-e)
      dNdeta[4] = -0.5 * (1 - x * x);

      dNdxi[5]  =  0.5 * (1 - e * e);
      dNdeta[5] =    -(1 + x) * e;

      dNdxi[6]  =      -x * (1 + e);
      dNdeta[6] =  0.5 * (1 - x * x);

      dNdxi[7]  = -0.5 * (1 - e * e);
      dNdeta[7] =    -(1 - x) * e;
    }
  };


//--------------------------------------------
// Core data
//--------------------------------------------

// Quadtree node (internal or leaf)
  struct TreeNode {
    u64 morton{0};            // Morton code of this node
    u32 level{0};             // Refinement level
    u32 child[4] {npos32, npos32, npos32, npos32}; // Children indices (npos32 = none)
    bool is_leaf{true};       // True if this node is a leaf
    u32 parent{npos32};       // Parent node index (npos32 if root)
  };

  struct Leaf {
    u64 morton{0};
    u32 level{0};
  };

  struct Field {
    Basis basis{Basis::Q1_Quad4};
    // Nodal values, one per unique global node in the basis registry.
    std::vector<double> nodal;  // size == _basisReg[(int)basis].nodes.size()

  };

  struct FEMNode {
    int gid;
    std::array<double, 2> parent;
    std::array<double, 2> physical;
  };

//--------------------------------------------
// QuadTree2D
//--------------------------------------------
  class QuadTree2D {
    public:
      QuadTree2D(u32 maxDepth, u32 minDepth = 0)
        : _maxDepth(maxDepth),
          _minDepth(std::min(minDepth, maxDepth)) {
        // --- 1. Build root TreeNode covering [-1,1]^2 ---
        TreeNode root;
        root.morton = interleave2(0, 0);
        root.level  = 0;
        root.is_leaf = true;
        root.parent  = npos32;
        root.child[0] = root.child[1] = root.child[2] = root.child[3] = npos32;

        _tree_nodes.clear();
        _tree_nodes.push_back(root);
        _root = 0;

        // --- 2. Root is the only active leaf at start ---
        _leaves.clear();
        _leaves.push_back(_root);

        // --- 3. Basis registries empty ---
        for (int b = 0; b < 3; ++b) {
          _basisReg[b].clear();
        }

        // --- 4. Geometry not set yet ---
        _geom_ready = false;
      }

      // ---- config ----
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

      void set_physical_quad9(const double x9[9], const double y9[9]) {
        for (int i = 0; i < 9; i++) {
          _X[i] = x9[i];
          _Y[i] = y9[i];
        }
        _geom_ready = true;
        // post_topology_update();

      }
      inline void require_geometry() const {
        assert(_geom_ready && "Call set_physical_quad9(...) before geometric operations.");
      }

// ---- leaves / indices ----
      u32 leaf_count() const {
        return (u32)_leaves.size();
      }

      const std::vector<u32>& leaf_indices() const {
        // expose [0..n-1] as leaf "ids" (stable across a pass)
        _leaf_ids.resize(_leaves.size());
        for (u32 i = 0; i < (u32)_leaves.size(); ++i) _leaf_ids[i] = i;
        return _leaf_ids;
      }

// ---- fields (per-leaf coefficients, compatible with your main) ----
      u32 add_field(Basis b) {
        Field f;
        f.basis = b;

        // Ensure registry exists so we can size nodal storage
        _activeBases.insert(b);

        post_topology_update();
        const auto& R = _basisReg[(int)b];
        f.nodal.assign(R.nodes.size(), 0.0);

        _fields.push_back(std::move(f));
        return (u32)_fields.size() - 1;
      }

      Field& field(u32 fid) {
        return _fields[fid];
      }
      const Field& field(u32 fid) const {
        return _fields[fid];
      }

// Resize all fields so their nodal arrays match current basis registries
      void resize_fields_to_nodes() {
        for (auto& f : _fields) {
          const auto& R = _basisReg[(int)f.basis];
          f.nodal.resize(R.nodes.size(), 0.0);
        }
      }

// ---- node generators (parent & physical) ----
// Fill vector with parent-space coordinates of interpolation nodes (by basis) for leaf
      void extract_node_parent_coords_in_level_range(Basis basis, u32 leaf_idx,
                                                     std::vector<std::array<double, 2>>& out_pts) const {

// Map leaf_idx → actual TreeNode
        const TreeNode& leaf = _tree_nodes[_leaves[leaf_idx]];

        double xi0, eta0, xi1, eta1;
        leaf_bounds(leaf, xi0, eta0, xi1, eta1);

        out_pts.clear();
        out_pts.reserve(9);
        switch (basis) {
        case Basis::Q1_Quad4: {
          out_pts = {
            {xi0, eta0}, {xi1, eta0}, {xi1, eta1}, {xi0, eta1}
          };
        }
        break;

        case Basis::Serendipity8: {
          const double xm = 0.5 * (xi0 + xi1);
          const double ym = 0.5 * (eta0 + eta1);
          out_pts = {
            {xi0, eta0}, {xi1, eta0}, {xi1, eta1}, {xi0, eta1},
            {xm, eta0}, {xi1, ym}, {xm, eta1}, {xi0, ym}
          };
        }
        break;

        case Basis::Q2_Quad9: {
          const double xm = 0.5 * (xi0 + xi1);
          const double ym = 0.5 * (eta0 + eta1);
          out_pts = {
            {xi0, eta0}, {xi1, eta0}, {xi1, eta1}, {xi0, eta1},
            {xm, eta0}, {xi1, ym}, {xm, eta1}, {xi0, ym}, {xm, ym}
          };
        }
        break;
        }
      }


      void leaf_physical_nodes(Basis basis, u32 leaf_idx,
                               std::vector<std::array<double, 2>>& out_xy) const {
        require_geometry();
        std::vector<std::array<double, 2>> xi;
        extract_node_parent_coords_in_level_range(basis, leaf_idx, xi);
        out_xy.resize(xi.size());
        for (size_t i = 0; i < xi.size(); ++i) out_xy[i] = parent_to_physical(xi[i][0], xi[i][1]);
      }

      // ---- adaptation ---------------------------------------------------
      template<class CoarsenPred, class RefinePred>
      std::size_t adapt_cycle(CoarsenPred&& should_coarsen,
                              RefinePred&&  should_refine,
                              Basis probe_basis = Basis::Q2_Quad9,
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
        newLeaves.reserve(_leaves.size() * 2);
        std::size_t refined = 0;

        for (u32 idx = 0; idx < leaf_count(); ++idx) {
          u32 leaf_node_idx = _leaves[idx];          // index into _tree_nodes
          const TreeNode leaf_copy = _tree_nodes[leaf_node_idx]; // snapshot (safe!)

          // probe geometry for refinement
          std::vector<std::array<double, 2>> pts_xi, pts_xy;
          extract_node_parent_coords_in_level_range(probe_basis, idx, pts_xi);
          leaf_physical_nodes(probe_basis, idx, pts_xy);
          std::vector<std::array<double, 9>> dummy;

          if (should_refine(idx, pts_xi, pts_xy, dummy) && leaf_copy.level < _maxDepth) {
            // split into 4 children
            u32 ix, iy;
            deinterleave2(leaf_copy.morton, ix, iy);

            for (int dy = 0; dy < 2; ++dy) {
              for (int dx = 0; dx < 2; ++dx) {
                u64 cmorton = interleave2(2 * ix + dx, 2 * iy + dy);

                TreeNode child;
                child.morton  = cmorton;
                child.level   = leaf_copy.level + 1;
                child.is_leaf = true;
                child.parent  = leaf_node_idx;
                child.child[0] = child.child[1] = child.child[2] = child.child[3] = npos32;


                u32 child_idx = (u32)_tree_nodes.size();
                _tree_nodes.push_back(child);

                // ✅ re-fetch parent AFTER push_back
                _tree_nodes[leaf_node_idx].child[(dy << 1) | dx] = child_idx;
                newLeaves.push_back(child_idx);
              }
            }

            _tree_nodes[leaf_node_idx].is_leaf = false;  // mark parent
            refined++;
          }
          else {
            newLeaves.push_back(leaf_node_idx);
          }
        }

        _leaves.swap(newLeaves);
        //post_topology_update();

        return refined;
      }

// Perform one coarsening pass: merge 4 siblings into parent if predicate allows
      template<class CoarsenPred>
      std::size_t coarsen_pass(CoarsenPred&& should_coarsen, Basis probe_basis) {
        struct Kid {
          u32 leaf_idx;
          u32 node_idx;
        };
        std::unordered_map<u64, std::vector<Kid>> groups;

        // group leaves by parent Morton
        for (u32 idx = 0; idx < leaf_count(); ++idx) {
          u32 node_idx = _leaves[idx];
          const TreeNode& lf = _tree_nodes[node_idx];
          if (lf.level == 0) continue;

          u32 ix, iy;
          deinterleave2(lf.morton, ix, iy);
          u64 pm = interleave2(ix / 2, iy / 2);
          groups[pm].push_back({idx, node_idx});
        }

        std::vector<char> mark_remove(leaf_count(), 0);
        std::vector<u32> to_add;
        to_add.reserve(groups.size());
        std::size_t coarsened = 0;

        for (auto& kv : groups) {
          const u64 pm = kv.first;
          auto& kids = kv.second;
          if (kids.size() != 4) continue;

          // Ensure all 4 kids share the same existing parent node
          const u32 pidx0 = _tree_nodes[kids[0].node_idx].parent;
          if (pidx0 == npos32) continue; // safety
          bool same_parent = true;
          for (int i = 1; i < 4; ++i)
            if (_tree_nodes[kids[i].node_idx].parent != pidx0) {
              same_parent = false;
              break;
            }
          if (!same_parent) continue;

          // Level checks
          u32 lev = _tree_nodes[kids[0].node_idx].level; // kids' level
          if (lev <= _minDepth && !_allowCoarsenBelowMinDepth) continue;

          // Probe coarse parent cell (use the *existing* parent)
          TreeNode& parent = _tree_nodes[pidx0];

          // parent.morton should equal pm; keep it as is
          double xi0, eta0, xi1, eta1;
          leaf_bounds(parent, xi0, eta0, xi1, eta1);

          std::vector<std::array<double, 2>> pts_xi;
          pts_xi.reserve(9);
          const double xm = 0.5 * (xi0 + xi1), ym = 0.5 * (eta0 + eta1);
          switch (probe_basis) {
          case Basis::Q1_Quad4:
            pts_xi = {{ {xi0, eta0}, {xi1, eta0}, {xi1, eta1}, {xi0, eta1} }};
            break;
          case Basis::Serendipity8:
            pts_xi = {{ {xi0, eta0}, {xi1, eta0}, {xi1, eta1}, {xi0, eta1},
                {xm, eta0}, {xi1, ym}, {xm, eta1}, {xi0, ym}
              }
            };
            break;
          case Basis::Q2_Quad9:
            pts_xi = {{ {xi0, eta0}, {xi1, eta0}, {xi1, eta1}, {xi0, eta1},
                {xm, eta0}, {xi1, ym}, {xm, eta1}, {xi0, ym}, {xm, ym}
              }
            };
            break;
          }

          std::vector<std::array<double, 2>> pts_xy;
          pts_xy.resize(pts_xi.size());
          for (size_t i = 0; i < pts_xi.size(); ++i)
            pts_xy[i] = parent_to_physical(pts_xi[i][0], pts_xi[i][1]);

          std::vector<std::array<double, 9>> dummy;

          if (should_coarsen((u32)pm, lev, pts_xi, pts_xy, dummy)) {
            // (1) remove kids from leaves
            for (auto& k : kids) {
              mark_remove[k.leaf_idx] = 1;
              _tree_nodes[k.node_idx].is_leaf = false; // they are no longer active leaves
            }

            // (2) make the existing parent a true leaf
            parent.is_leaf = true;
            parent.child[0] = parent.child[1] = parent.child[2] = parent.child[3] = npos32;

            // (3) add parent to new leaf set
            to_add.push_back(pidx0);
            coarsened++;
          }
        }

        // rebuild leaf set
        if (coarsened > 0) {
          std::vector<u32> newLeaves;
          newLeaves.reserve(_leaves.size());
          for (u32 i = 0; i < leaf_count(); ++i)
            if (!mark_remove[i]) newLeaves.push_back(_leaves[i]);
          newLeaves.insert(newLeaves.end(), to_add.begin(), to_add.end());
          _leaves.swap(newLeaves);

          //post_topology_update();
        }
        return coarsened;
      }

// ---- geometry map ----
      std::array<double, 2> parent_to_physical(double xi, double eta) const {
        require_geometry();
        double N9[9];
        Quad9Shape::N(xi, eta, N9);
        double X = 0, Y = 0;
        for (int a = 0; a < 9; a++) {
          X += N9[a] * _X[a];
          Y += N9[a] * _Y[a];
        }
        return {X, Y};
      }

// Add to QuadTree2D: call this after any leaf rebuild
      inline void post_topology_update() {
        rebuild_connectivity_active_bases();
        rebuild_leafpos_lookup();
        resize_fields_to_nodes();
      }


      bool write_binary_vtu(const std::string &filename, u32 fid, const std::string &name,
                            bool cell_centered = false) const {
        require_geometry();
        if (fid >= _fields.size()) return false;

        const Field &fld = _fields[fid];
        Basis b = fld.basis;
        const BasisRegistry &R = _basisReg[(int)b];

        const size_t numCells = R.elem2glob.size();
        if (numCells == 0) return false;

        // -------------------------
        // Points (global nodes)
        // -------------------------
        std::vector<double> flatPoints;
        flatPoints.reserve(R.nodes.size() * 3);
        for (auto &n : R.nodes) {
          flatPoints.push_back(n.physical[0]);
          flatPoints.push_back(n.physical[1]);
          flatPoints.push_back(0.0); // 2D
        }

        // -------------------------
        // Connectivity
        // -------------------------
        std::vector<int> connectivity;
        std::vector<int> offsets;
        std::vector<unsigned char> types;
        //connectivity.reserve(numCells * fld.dofs_per_cell);

        // Reserve approximately (safe upper bound if uniform)
        connectivity.reserve(numCells * (b == Basis::Q1_Quad4 ? 4 : (b == Basis::Serendipity8 ? 8 : 9)));

        offsets.reserve(numCells);
        types.reserve(numCells);

        int offset = 0;
        for (size_t e = 0; e < numCells; ++e) {
          const auto &conn = R.elem2glob[e];
          for (int gid : conn) connectivity.push_back(gid);
          offset += (int)conn.size();
          offsets.push_back(offset);

          types.push_back(vtk_cell_type(fld.basis));
        }

        // -------------------------
        // Field values
        // -------------------------
        std::vector<double> pointData;
        std::vector<double> cellData;
        if (cell_centered) {
          cellData.reserve(numCells);
          for (size_t e = 0; e < numCells; ++e) {
            const auto &conn = R.elem2glob[e];
            double v = 0.0;
            for (int gid : conn) v += fld.nodal[gid];
            cellData.push_back(v / (double)conn.size());
          }
        }
        else {
          // Point data is directly the nodal storage
          pointData = fld.nodal;
        }


        // -------------------------
        // Write XML
        // -------------------------
        std::ofstream os(filename);
        if (!os) return false;

        os << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        os << "  <UnstructuredGrid>\n";
        os << "    <Piece NumberOfPoints=\"" << R.nodes.size()
           << "\" NumberOfCells=\"" << numCells << "\">\n";

        // Points
        os << "      <Points>\n";
        write_binary_array(os, "Float64", "", 3, flatPoints);
        os << "      </Points>\n";

        // Cells
        os << "      <Cells>\n";
        write_binary_array(os, "Int32", "connectivity", 1, connectivity);
        write_binary_array(os, "Int32", "offsets", 1, offsets);
        write_binary_array(os, "UInt8", "types", 1, types);
        os << "      </Cells>\n";

        // Data
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



// Accessors for coarse geometry coordinates (const-ref)
      inline const std::array<double, 9>& get_X() const {
        return _X;
      }
      inline const std::array<double, 9>& get_Y() const {
        return _Y;
      }

// Collect reference coordinates (xi,eta) of nodes for all leaves in a level range.
      std::vector<std::array<double, 2>>
      extract_node_parent_coords_in_level_range(u32 lev_min, u32 lev_max, Basis basis) const {
        std::vector<std::array<double, 2>> coords;

        for (u32 leaf : leaf_indices()) {
          u32 lev = level_of(leaf);          // use accessor, not _nodes
          if (lev < lev_min || lev > lev_max) continue;

          std::vector<std::array<double, 2>> xi;
          extract_node_parent_coords_in_level_range(basis, leaf, xi);

          coords.insert(coords.end(), xi.begin(), xi.end());
        }
        return coords;
      }

// Return refinement level of a given leaf index
      u32 level_of(u32 leaf_idx) const {
        u32 node_idx = _leaves[leaf_idx];          // leaf index → node index
        return _tree_nodes[node_idx].level;        // get level from TreeNode
      }

//Return global maximum depth allowed
      u32 max_depth() const {
        return _maxDepth;
      }

// Locate leaf containing (xi,eta) in parent space [-1,1]^2 using Morton traversal
      u32 locate_leaf_on_parent(double xi, double eta) const {
        if (xi < -1.0 || xi > 1.0 || eta < -1.0 || eta > 1.0) return npos32;

        // Scale to integer grid at max depth
        const double scale = double(1u << _maxDepth) / 2.0;  // maps [-1,1] -> [0, 2^maxDepth)
        u32 ix = std::min<u32>((u32)((xi + 1.0) * scale), (1u << _maxDepth) - 1);
        u32 iy = std::min<u32>((u32)((eta + 1.0) * scale), (1u << _maxDepth) - 1);

        u32 node = _root;
        for (u32 level = 0; level < _maxDepth; ++level) {
          const TreeNode& n = _tree_nodes[node];
          if (n.is_leaf) return node;  // found a leaf

          // Extract bit for this level (from MSB downwards)
          const int shift = _maxDepth - level - 1;
          const int qx = (ix >> shift) & 1u;
          const int qy = (iy >> shift) & 1u;
          const int q  = (qy << 1) | qx; // quadrant index 0=LL,1=LR,2=UL,3=UR

          u32 c = n.child[q];
          if (c == npos32) return node; // child not refined -> current is leaf
          node = c;
        }
        return node;
      }


// Evaluate field directly in parent coordinates (xi, eta).
// Locates the leaf, maps to local [-1,1]^2, and interpolates.
      bool evaluate_field_on_parent(u32 fid, double xi, double eta, double& value) const {
        if (fid >= _fields.size()) return false;

        u32 leaf_node_idx;
        double shat, that;

        // 1) locate leaf and local reference coords
        if (!locate_leaf_on_parent_and_ref(xi, eta, leaf_node_idx, shat, that)) {
          return false;
        }


        // 2) leaf node index -> position in coefficient storage
        const u32 leaf_pos = _node2leafpos[leaf_node_idx]; //BIG WINNER!!!
        if (leaf_pos == npos32) return false;

        // 3) access field + coefficients
        const Field& f = _fields[fid];
        const BasisRegistry& R = _basisReg[(int)f.basis];

        // Connectivity for this element (leaf_pos is aligned with rebuild_connectivity)
        const auto& conn = R.elem2glob[leaf_pos];

        // 4) interpolate
        switch (f.basis) {
        case Basis::Q1_Quad4: {
          double N4[4];
          Shapes::Q1(shat, that, N4);
          value = 0.0;
          for (int a = 0; a < 4; ++a) value += N4[a] * f.nodal[conn[a]];
        }
        break;
        case Basis::Serendipity8: {
          double N8[8];
          Shapes::Q8(shat, that, N8);
          value = 0.0;
          for (int a = 0; a < 8; ++a) value += N8[a] * f.nodal[conn[a]];
        }
        break;
        case Basis::Q2_Quad9: {
          double N9[9];
          Quad9Shape::N(shat, that, N9);
          value = 0.0;
          for (int a = 0; a < 9; ++a) value += N9[a] * f.nodal[conn[a]];
        }
        break;
        }
        return true;
      }


      inline void local_ref_fast(u32 leaf_node_idx, double xi, double eta,
                                 double& shat, double& that) const {
        const TreeNode& n = _tree_nodes[leaf_node_idx];
        u32 ix, iy; deinterleave2(n.morton, ix, iy);
        const double N = double(1u << n.level);

        shat = N * xi  + N - double((ix << 1) + 1);
        that = N * eta + N - double((iy << 1) + 1);
      }

      inline bool locate_leaf_on_parent_and_ref(double xi, double eta,
                                                u32& leaf_node_idx,
                                                double& shat, double& that) const {
        // 1) Find the leaf containing (xi,eta) in parent space.
        leaf_node_idx = locate_leaf_on_parent(xi, eta);
        if (leaf_node_idx == npos32) return false;

        // 2) Convert to local [-1,1]^2 using level+Morton (no divisions).
        local_ref_fast(leaf_node_idx, xi, eta, shat, that);

        return true;
      }


// Locate leaf containing (xi,eta) and compute local reference coords (shat,that) in [-1,1]^2.
      bool locate_leaf_on_parent_and_ref2(double xi, double eta,
                                          u32& leaf_node_idx,
                                          double& shat, double& that) const {
        leaf_node_idx = locate_leaf_on_parent(xi, eta);
        if (leaf_node_idx == npos32) return false;

        double xi0, eta0, xi1, eta1;
        leaf_bounds(_tree_nodes[leaf_node_idx], xi0, eta0, xi1, eta1);
        const double dx = xi1 - xi0;
        const double dy = eta1 - eta0;
        if (dx <= 0.0 || dy <= 0.0) return false;

        const double xm = 0.5 * (xi0 + xi1);
        const double ym = 0.5 * (eta0 + eta1);
        shat = 2.0 * (xi  - xm) / dx;
        that = 2.0 * (eta - ym) / dy;
        return true;
      }





// Refine tree so that each physical point lies in a leaf of at least 'maxDepthTarget'
      void refine_to_contain_points(const std::vector<std::array<double, 2>>& pts,
                                    u32 maxDepthTarget) {
        if (pts.empty()) return;

        // Iterate until all points are in leaves with level >= target (or no more splits possible)
        for (;;) {
          // Mark the *node indices* (not positions) that need refinement this pass
          std::unordered_set<u32> nodes_to_refine;
          nodes_to_refine.reserve(pts.size());

          for (const auto& p : pts) {
            double xi, eta;
            if (!inverse_map_quad9(p[0], p[1], xi, eta)) continue;

            const u32 node = locate_leaf_on_parent(xi, eta);
            if (node == npos32) continue;

            if (_tree_nodes[node].level < maxDepthTarget)
              nodes_to_refine.insert(node);
          }

          if (nodes_to_refine.empty()) break;

          // Refine exactly those leaves; use the cheapest probe basis since we ignore pts_xi/pts_xy
          auto should_refine = [&](u32 leaf_pos,
                                   const std::vector<std::array<double, 2>>&,
                                   const std::vector<std::array<double, 2>>&,
          const std::vector<std::array<double, 9>>&) -> bool {
            const u32 node = _leaves[leaf_pos];
            return nodes_to_refine.count(node) &&
            (_tree_nodes[node].level < maxDepthTarget);
          };

          const std::size_t nsplit = refine_pass(should_refine, Basis::Q1_Quad4);
          if (nsplit == 0) break; // nothing changed (e.g., already at/above target)
        }

        // Keep mesh 1-irregular
        enforce_balance();
        post_topology_update();
      }

// Enforce 1-irregularity: adjacent leaves may differ by at most 1 level.
      // void enforce_balance() {
      //   bool changed = true;
      //   while (changed) {
      //     changed = false;
      //     // work on a snapshot to avoid invalidating indices while refining
      //     const std::vector<u32> leaves_snapshot = _leaves;
      //
      //     for (u32 leaf_node_idx : leaves_snapshot) {
      //       if (_tree_nodes[leaf_node_idx].is_leaf == false) continue; // may have been refined already
      //       const u32 lev = _tree_nodes[leaf_node_idx].level;
      //
      //       for (int dir = 0; dir < 4; ++dir) {
      //         u32 nb = neighbor_leaf_by_dir(leaf_node_idx, dir);
      //         if (nb == npos32) continue;
      //         if (_tree_nodes[nb].is_leaf == false) continue; // very unlikely, but safe
      //
      //         const u32 lev_nb = _tree_nodes[nb].level;
      //
      //         // If neighbor is too coarse relative to this leaf, refine neighbor
      //         if (lev > lev_nb + 1) {
      //           changed |= refine_leaf_once(nb);
      //         }
      //         // Symmetric check (in case we iterate this leaf before the finer neighbor does)
      //         else if (lev_nb > lev + 1) {
      //           changed |= refine_leaf_once(leaf_node_idx);
      //         }
      //       }
      //     }
      //   }
      //   // Rebuild once at the end
      //   post_topology_update();
      // }

      void enforce_balance() {
        // local enqueue that avoids duplicates
        auto enqueue = [](std::deque<u32>& q, std::vector<char>& inq, u32 n) {
          if (n == npos32) return;
          if (n >= inq.size()) return;
          if (!inq[n]) {
            inq[n] = 1;
            q.push_back(n);
          }
        };

        // Seed queue with current leaves
        std::deque<u32> q;
        q.insert(q.end(), _leaves.begin(), _leaves.end());

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

          for (int dir = 0; dir < 4; ++dir) {
            u32 nb = neighbor_leaf_by_dir(leaf, dir);
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
                // enqueue neighbor's neighbors (they are the only ones that can be affected)
                for (int d = 0; d < 4; ++d) {
                  u32 nb2 = neighbor_leaf_by_dir(nb, d);
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
                for (int d = 0; d < 4; ++d) {
                  u32 nb2 = neighbor_leaf_by_dir(leaf, d);
                  if (nb2 != npos32) enqueue(q, inq, nb2);
                }
                // no point checking other dirs for the old leaf; it’s not a leaf now
                break;
              }
            }
          }
        }

        // Rebuild derived data once
      }




// Conservative coarsen cycle using snapshot + parent coords; rebuild all fields
      std::size_t coarsen_only_cycle_safe(u32 fid,
                                          double tau_coarse,
                                          u32 max_passes = 10,
                                          Basis probe_basis = Basis::Q2_Quad9) {
        // Freeze current state for conservative evaluation/transfer
        QuadTree2D snapshot = *this;

        // Coarsen predicate (evaluated on the snapshot, in parent coords)
        auto pred = [&](u32 /*parent_morton*/, u32 level,
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

        // Coarsen passes + enforce 1-irregularity after each pass
        for (u32 pass = 0; pass < max_passes; ++pass) {
          std::size_t c = coarsen_pass(pred, probe_basis);
          if (c == 0) break;
          total += c;

        }
        // enforce 1-level irregularity
        enforce_balance();
        post_topology_update();

        // Rebuild all fields from the snapshot (conservative transfer) on nodal sets
        for (u32 f = 0; f < _fields.size(); ++f) {
          rebuild_field_from(snapshot, f);
        }

        return total;
      }

      void rebuild_connectivity_active_bases() {
        if (!_geom_ready) return;
        for (Basis b : _activeBases) {
          rebuild_connectivity(b);
        }
      }

// Returns the unique global FEM nodes (per-basis), with parent & physical coords.
      const std::vector<FEMNode>& basis_nodes(Basis b) const {
        return _basisReg[(int)b].nodes;
      }

// Returns the element-to-global connectivity (per-basis).
      const std::vector<std::vector<int>>& basis_connectivity(Basis b) const {
        return _basisReg[(int)b].elem2glob;
      }



    private: //methods

// Rebuild field fid on *this* from source tree 'src' by sampling at parent nodes (global nodal storage)
      void rebuild_field_from(const QuadTree2D& src, u32 fid) {

        assert(fid < _fields.size() && fid < src._fields.size());
        Field& dst = _fields[fid];
        const Field& s = src._fields[fid];
        assert(dst.basis == s.basis && "Source/dest field bases must match");

        const BasisRegistry& Rdst = _basisReg[(int)dst.basis];
        dst.nodal.resize(Rdst.nodes.size());

        // Sample source field at each destination nodal parent coordinate
        for (size_t gid = 0; gid < Rdst.nodes.size(); ++gid) {
          const auto& pr = Rdst.nodes[gid].parent;
          double val;
          const bool ok = src.evaluate_field_on_parent(fid, pr[0], pr[1], val);
          dst.nodal[gid] = ok ? val : std::numeric_limits<double>::quiet_NaN();
        }
      }


// Compute a neighbor leaf by probing an edge midpoint in parent space.
// dir: 0=left (−x), 1=right (+x), 2=down (−y), 3=up (+y).
      u32 neighbor_leaf_by_dir(u32 leaf_node_idx, int dir) const {
        const TreeNode& n = _tree_nodes[leaf_node_idx];

        double xi0, eta0, xi1, eta1;
        leaf_bounds(n, xi0, eta0, xi1, eta1);

        const double xm = 0.5 * (xi0 + xi1);
        const double ym = 0.5 * (eta0 + eta1);

        // step slightly outside the edge to guarantee we land in the neighbor
        const double dx = (xi1 - xi0);
        const double dy = (eta1 - eta0);
        const double epsx = 1e-6 * dx;
        const double epsy = 1e-6 * dy;

        double qx, qy;
        switch (dir) {
        case 0:
          qx = xi0 - epsx;
          qy = ym;
          break; // left
        case 1:
          qx = xi1 + epsx;
          qy = ym;
          break; // right
        case 2:
          qx = xm;
          qy = eta0 - epsy;
          break; // down
        default:/*3*/
          qx = xm;
          qy = eta1 + epsy;
          break; // up
        }

        // outside global domain? no neighbor
        if (qx < -1.0 || qx > 1.0 || qy < -1.0 || qy > 1.0) return npos32;

        return locate_leaf_on_parent(qx, qy);
      }

// Refine a single leaf node into 4 children, update _leaves in place.
// Returns true if a refinement actually occurred.
      bool refine_leaf_once(u32 leaf_node_idx) {
        TreeNode& parent = _tree_nodes[leaf_node_idx];
        if (!parent.is_leaf) return false;
        if (parent.level >= _maxDepth) return false;

        // mark parent as internal
        parent.is_leaf = false;

        u32 ix, iy;
        deinterleave2(parent.morton, ix, iy);

        // Remove parent from _leaves (swap-erase)
        for (u32 i = 0; i < _leaves.size(); ++i) {
          if (_leaves[i] == leaf_node_idx) {
            _leaves[i] = _leaves.back();
            _leaves.pop_back();
            break;
          }
        }

        // create 2x2 children
        for (int dy = 0; dy < 2; ++dy) {
          for (int dx = 0; dx < 2; ++dx) {
            const u64 cmorton = interleave2(2 * ix + dx, 2 * iy + dy);

            TreeNode child;
            child.morton  = cmorton;
            child.level   = parent.level + 1;
            child.is_leaf = true;
            child.parent  = leaf_node_idx;
            child.child[0] = child.child[1] = child.child[2] = child.child[3] = npos32;

            const u32 child_idx = (u32)_tree_nodes.size();
            _tree_nodes.push_back(child);

            const int q = (dy << 1) | dx;
            _tree_nodes[leaf_node_idx].child[q] = child_idx;

            // add to leaves
            _leaves.push_back(child_idx);
          }
        }
        return true;
      }

// Ensure all leaves are refined down to at least _minDepth
      void ensure_min_depth() {
        if (_minDepth == 0) return;

        bool changed = true;
        while (changed) {
          changed = false;
          std::vector<u32> newLeaves;

          for (u32 leaf_idx : _leaves) {
            const TreeNode node_copy = _tree_nodes[leaf_idx];

            if (node_copy.level < _minDepth) {
              // refine this leaf into 4 children
              u32 ix, iy;
              deinterleave2(node_copy.morton, ix, iy);

              _tree_nodes[leaf_idx].is_leaf = false;

              for (int dy = 0; dy < 2; ++dy) {
                for (int dx = 0; dx < 2; ++dx) {
                  u64 cmorton = interleave2(2 * ix + dx, 2 * iy + dy);
                  u32 cindex = (u32)_tree_nodes.size();

                  TreeNode child;
                  child.morton  = cmorton;
                  child.level   = node_copy.level + 1;
                  child.is_leaf = true;
                  child.parent  = leaf_idx;
                  child.child[0] = child.child[1] = child.child[2] = child.child[3] = npos32;


                  _tree_nodes.push_back(child);
                  _tree_nodes[leaf_idx].child[(dy << 1) | dx] = cindex;
                  newLeaves.push_back(cindex);
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


// Compute physical bounds in parent space [-1,1]^2 for a given TreeNode
      inline void leaf_bounds(const TreeNode& node,
                              double& xi0, double& eta0,
                              double& xi1, double& eta1) const {
        u32 ix, iy;
        deinterleave2(node.morton, ix, iy);

        const double N  = double(1u << node.level); // grid resolution at this level
        const double dx = 2.0 / N;
        const double dy = 2.0 / N;

        xi0  = -1.0 + ix * dx;
        eta0 = -1.0 + iy * dy;
        xi1  = xi0 + dx;
        eta1 = eta0 + dy;
      }

      struct BasisRegistry {
        std::unordered_map<u64, int> nodeMap;            // key -> gid
        std::vector<FEMNode> nodes;                      // gid -> FEM node (parent + physical)
        std::vector<std::vector<int>> elem2glob;         // per-element connectivity
        void clear() {
          nodeMap.clear();
          nodes.clear();
          elem2glob.clear();
        }
      };

      int get_or_insert_gid(BasisRegistry &R, double xi, double eta) {
        // Need +1 bit so level=_maxDepth mid/center nodes (half-steps) map to integers.
        const u32 nodeBits = _maxDepth + 1;
        const u32 nodesN   = (1u << nodeBits); // indices 0..nodesN

        // Optional: clamp xi,eta early to avoid NaNs/inf sneaking in
        if (!std::isfinite(xi) || !std::isfinite(eta)) {
          // handle or assert as you prefer
          assert(false && "Non-finite xi/eta");
        }

        auto to_idx = [nodesN](double s)->u32 {
          // Clamp to kill 1.0 ± eps and -1.0 ± eps cases
          if (s <= -1.0) return 0u;
          if (s >=  1.0) return nodesN;
          double t = (s + 1.0) * double(nodesN) * 0.5; // [-1,1] -> [0,nodesN]
          long long li = llround(t);
          if (li < 0) li = 0;
          if (li > (long long)nodesN) li = (long long)nodesN;
          return (u32)li;
        };

        const u32 ix = to_idx(xi);
        const u32 iy = to_idx(eta);

        const u64 key = (u64(ix) << 32) | u64(iy); // or Morton(bits=nodeBits)

        auto it = R.nodeMap.find(key);
        if (it != R.nodeMap.end()) return it->second;

        const int gid = (int)R.nodes.size();
        R.nodes.push_back(FEMNode{gid, {xi, eta}, parent_to_physical(xi, eta)});
        R.nodeMap.emplace(key, gid);
        return gid;
      }

      void rebuild_connectivity(Basis b) {
        BasisRegistry &R = _basisReg[(int)b];
        R.clear();

        std::vector<std::array<double, 2>> xi;
        for (u32 e = 0; e < leaf_count(); ++e) {
          extract_node_parent_coords_in_level_range(b, e, xi);
          std::vector<int> conn;
          conn.reserve(xi.size());

          for (auto &p : xi) {
            int gid = get_or_insert_gid(R, p[0], p[1]); // now uses registry
            conn.push_back(gid);
          }
          R.elem2glob.push_back(std::move(conn));
        }
      }

// ======= VTU (ASCII) =======
      bool write_ascii_vtu(const std::string& filename, u32 fid, const std::string& name,
                           bool cell_centered) const {
        if (!_geom_ready) return false;
        if (fid >= _fields.size()) return false;

        const Field& fld = _fields[fid];
        const Basis b = fld.basis;
        const BasisRegistry& R = _basisReg[(int)b];

        // point data is directly the nodal values
        const std::vector<double>& pointData = fld.nodal;

        std::ofstream os(filename);
        if (!os) return false;
        os << std::setprecision(16);
        os << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        os << "  <UnstructuredGrid>\n";
        os << "    <Piece NumberOfPoints=\"" << R.nodes.size()
           << "\" NumberOfCells=\"" << R.elem2glob.size() << "\">\n";

        // Points
        os << "      <Points>\n";
        os << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (auto &n : R.nodes)
          os << "          " << n.physical[0] << " " << n.physical[1] << " 0\n";
        os << "        </DataArray>\n";
        os << "      </Points>\n";

        // Cells
        os << "      <Cells>\n";
        os << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
        for (auto &conn : R.elem2glob) {
          for (size_t i = 0; i < conn.size(); ++i)
            os << conn[i] << (i + 1 == conn.size() ? '\n' : ' ');
        }
        os << "        </DataArray>\n";
        os << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
        int off = 0;
        for (auto &conn : R.elem2glob) {
          off += (int)conn.size();
          os << "          " << off << "\n";
        }
        os << "        </DataArray>\n";
        os << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
        const int vtkType = vtk_cell_type(b);
        for (size_t k = 0; k < R.elem2glob.size(); ++k)
          os << "          " << vtkType << "\n";
        os << "        </DataArray>\n";
        os << "      </Cells>\n";

        // Data
        if (!cell_centered) {
          os << "      <PointData Scalars=\"" << name << "\">\n";
          os << "        <DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
          for (double v : pointData) os << "          " << v << "\n";
          os << "        </DataArray>\n";
          os << "      </PointData>\n";
          os << "      <CellData>\n";
          os << "      </CellData>\n";
        }
        else {
          // simple center eval: take average of nodal values
          std::vector<double> cellData(R.elem2glob.size(), 0.0);
          for (size_t e = 0; e < R.elem2glob.size(); ++e) {
            double s = 0;
            for (int gid : R.elem2glob[e]) s += pointData[gid];
            cellData[e] = s / R.elem2glob[e].size();
          }
          os << "      <CellData Scalars=\"" << name << "\">\n";
          os << "        <DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
          for (double v : cellData) os << "          " << v << "\n";
          os << "        </DataArray>\n";
          os << "      </CellData>\n";
          os << "      <PointData>\n";
          os << "      </PointData>\n";
        }

        os << "    </Piece>\n";
        os << "  </UnstructuredGrid>\n";
        os << "</VTKFile>\n";
        return true;
      }


// Insert a unique FEM node (global interpolation node) based on parent coords (xi, eta)
// into the registry for the given basis
      int get_or_insert_node(BasisRegistry& R, double xi, double eta, u32 maxDepth) {
        // snap to dyadic grid at maxDepth to build a unique Morton key
        u32 ix = (u32)std::llround((xi + 1.0) * (1u << maxDepth) / 2.0);
        u32 iy = (u32)std::llround((eta + 1.0) * (1u << maxDepth) / 2.0);
        u64 key = interleave2(ix, iy);

        // check if already exists
        auto it = R.nodeMap.find(key);
        if (it != R.nodeMap.end())
          return it->second;

        // assign new global id
        int gid = static_cast<int>(R.nodes.size());
        R.nodeMap[key] = gid;

        // create FEM node with parent and physical coordinates
        FEMNode node;
        node.gid      = gid;
        node.parent   = {xi, eta};
        node.physical = parent_to_physical(xi, eta);


        R.nodes.push_back(std::move(node));
        return gid;
      }


// Q9 inverse with Q1 warm-up; same control flow & loops (no unrolling)
      bool inverse_map_quad9(double x, double y, double& xi, double& eta,
                             int maxIts = 25, double tol = 1e-12) const {
        require_geometry();

        // start at center
        xi = 0.0; eta = 0.0;
        const double tol2 = tol * tol;

        // buffers
        double N[9], dNdxi[9], dNdeta[9];

        // --- Warm-up with Q1: residual and Jacobian from corners only ---
        Shapes::Q1(xi, eta, N);           // N[0..3] valid
        double Xp = 0.0, Yp = 0.0;
        for (int a = 0; a < 4; ++a) {
          Xp += N[a] * _X[a];
          Yp += N[a] * _Y[a];
        }
        double rx = Xp - x, ry = Yp - y;

        Shapes::Q1_dN(xi, eta, dNdxi, dNdeta); // fills [0..3]
        double dXdxi = 0.0, dXdeta = 0.0, dYdxi = 0.0, dYdeta = 0.0;
        for (int a = 0; a < 4; ++a) {
          dXdxi  += dNdxi[a] * _X[a];  dXdeta += dNdeta[a] * _X[a];
          dYdxi  += dNdxi[a] * _Y[a];  dYdeta += dNdeta[a] * _Y[a];
        }

        // --- Newton iterations ---
        for (int it = 0; it < maxIts; ++it) {
          // step using CURRENT residual (rx,ry) and CURRENT J
          const double detJ = dXdxi * dYdeta - dXdeta * dYdxi;
          if (std::fabs(detJ) < 1e-20) break;

          const double inv = 1.0 / detJ;
          const double dXi  = (dYdeta * rx - dXdeta * ry) * inv;
          const double dEta = (-dYdxi * rx + dXdxi * ry) * inv;

          // update iterate (light clamp)
          xi  = std::max(-1.0, std::min(1.0, xi  - dXi));
          eta = std::max(-1.0, std::min(1.0, eta - dEta));

          // recompute N at updated (xi,eta), update residual for convergence check
          Quad9Shape::N(xi, eta, N);
          Xp = 0.0; Yp = 0.0;
          for (int a = 0; a < 9; ++a) {
            Xp += N[a] * _X[a];
            Yp += N[a] * _Y[a];
          }
          rx = Xp - x;  ry = Yp - y;
          if ((rx * rx + ry * ry) < tol2) return true;

          // build NEXT Jacobian at updated (xi,eta) using Q9 derivatives
          Quad9Shape::dN(xi, eta, dNdxi, dNdeta);
          dXdxi = dXdeta = dYdxi = dYdeta = 0.0;   // reset accumulators
          for (int a = 0; a < 9; ++a) {
            dXdxi  += dNdxi[a] * _X[a];  dXdeta += dNdeta[a] * _X[a];
            dYdxi  += dNdxi[a] * _Y[a];  dYdeta += dNdeta[a] * _Y[a];
          }
        }

        return false; // not converged
      }




      // --- private: Q1 (Quad4) inverse map (Newton in parent space), using Shapes::Q1 ---
      bool inverse_map_quad4(double x, double y, double& xi, double& eta,
                             int maxIts = 12, double tol = 1e-12) const {
        require_geometry();

        // Start at center
        xi = 0.0; eta = 0.0;
        const double tol2    = tol * tol;
        const double tol_dx2 = 1e-20; // ~1e-10 step

        // // Cache vertices (assumed first 4 entries, in standard Quad4 order)
        // const double X0 = _X[0], Y0 = _Y[0];
        // const double X1 = _X[1], Y1 = _Y[1];
        // const double X2 = _X[2], Y2 = _Y[2];
        // const double X3 = _X[3], Y3 = _Y[3];

        for (int it = 0; it < maxIts; ++it) {

          // Values from your Shapes::Q1
          double N4[4], dNdxi[4], dNdeta[4];
          Shapes::Q1(xi, eta, N4);
          Shapes::Q1_dN(xi, eta, dNdxi, dNdeta);/*

          // Derivatives (standard bilinear basis, node order:
          // 0:(-1,-1), 1:(1,-1), 2:(1,1), 3:(-1,1))
          dNdxi[0] = -0.25 * (1.0 - eta);
          dNdeta[0] = -0.25 * (1.0 - xi);
          dNdxi[1] =  0.25 * (1.0 - eta);
          dNdeta[1] = -0.25 * (1.0 + xi);
          dNdxi[2] =  0.25 * (1.0 + eta);
          dNdeta[2] =  0.25 * (1.0 + xi);
          dNdxi[3] = -0.25 * (1.0 + eta);
          dNdeta[3] =  0.25 * (1.0 - xi);*/

          const auto& X = _X;
          const auto& Y = _Y;

          // current mapped point and Jacobian
          double Xp = 0.0, Yp = 0.0, dXdxi = 0.0, dXdeta = 0.0, dYdxi = 0.0, dYdeta = 0.0;
          for (int a = 0; a < 4; ++a) {
            Xp     += N4[a]    * X[a];
            Yp     += N4[a]    * Y[a];
            dXdxi  += dNdxi[a] * X[a];
            dXdeta += dNdeta[a] * X[a];
            dYdxi  += dNdxi[a] * Y[a];
            dYdeta += dNdeta[a] * Y[a];
          }

          double rx = Xp - x, ry = Yp - y;
          double nrm2 = (rx * rx + ry * ry);
          if (nrm2 < tol2) return true;

          double detJ = dXdxi * dYdeta - dXdeta * dYdxi;
          if (std::abs(detJ) < 1e-20) break; // singular

          // Newton step: [dXi, dEta] = J^{-1} * r
          double dXi  = (dYdeta * rx - dXdeta * ry) / detJ;
          double dEta = (-dYdxi * rx + dXdxi * ry) / detJ;

          xi  -= dXi;
          eta -= dEta;

          // keep iterates reasonable (light clamp)
          xi  = std::max(-1.5, std::min(1.5, xi));
          eta = std::max(-1.5, std::min(1.5, eta));
        }
        return false; // not converged
      }



// --- private: split one leaf and update _leaves ---
      bool refine_leaf_and_update_leaves(u32 leaf_node_idx) {
        if (leaf_node_idx == npos32) return false;
        TreeNode& parent = _tree_nodes[leaf_node_idx];
        if (!parent.is_leaf || parent.level >= _maxDepth) return false;

        u32 ix, iy;
        deinterleave2(parent.morton, ix, iy);

        for (int dy = 0; dy < 2; ++dy) {
          for (int dx = 0; dx < 2; ++dx) {
            const u64 cmorton = interleave2(2 * ix + dx, 2 * iy + dy);
            TreeNode child;
            child.morton  = cmorton;
            child.level   = parent.level + 1;
            child.is_leaf = true;
            child.parent  = leaf_node_idx;
            child.child[0] = child.child[1] = child.child[2] = child.child[3] = npos32;

            const u32 cidx = (u32)_tree_nodes.size();
            _tree_nodes.push_back(child);
            parent.child[(dy << 1) | dx] = cidx;
          }
        }
        parent.is_leaf = false;

        // replace parent in _leaves with its 4 children
        auto it = std::find(_leaves.begin(), _leaves.end(), leaf_node_idx);
        if (it != _leaves.end()) {
          // erase parent
          _leaves.erase(it);
          // append children
          for (int q = 0; q < 4; ++q) _leaves.push_back(parent.child[q]);
        }
        else {
          // if for some reason parent wasn't in _leaves, just append children
          for (int q = 0; q < 4; ++q) _leaves.push_back(parent.child[q]);
        }
        return true;
      }

      mutable std::vector<u32> _node2leafpos; // size == _tree_nodes.size(), npos32 default

      inline void rebuild_leafpos_lookup() {
        _node2leafpos.assign(_tree_nodes.size(), npos32);
        for (u32 i = 0; i < _leaves.size(); ++i) _node2leafpos[_leaves[i]] = i;
      }

    private: //data
// config
      u32  _maxDepth;
      u32  _minDepth{0};
      bool _allowCoarsenBelowMinDepth{true};

// geometry (Q2)
      std::array<double, 9> _X {}, _Y {};
      bool _geom_ready{false};

// mesh leaves
      std::vector<u32> _leaves;
      mutable std::vector<u32> _leaf_ids;

// fields (per-leaf coefficients for API compatibility)
      std::vector<Field> _fields;

// per-basis global node registries + connectivity
      BasisRegistry _basisReg[3]; // 0:Q1, 1:S8, 2:Q2

      std::vector<TreeNode> _tree_nodes;        // full tree hierarchy
      u32 _root{0};

      struct BasisHasher {
        size_t operator()(fem::Basis b) const noexcept {
          return std::hash<uint8_t>()(static_cast<uint8_t>(b));
        }
      };
      std::unordered_set<Basis, BasisHasher> _activeBases;
  };

} // namespace fem

