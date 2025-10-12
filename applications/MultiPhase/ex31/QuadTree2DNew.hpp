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
    static inline void Q1(double xi, double eta, double N4[4]) {
      const double a = 0.25 * (1 - xi), b = 0.25 * (1 + xi);
      N4[0] = a * (1 - eta);
      N4[1] = b * (1 - eta);
      N4[2] = b * (1 + eta);
      N4[3] = a * (1 + eta);
    }
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
    // Basis basis{Basis::Q1_Quad4};
    // u32   dofs_per_cell{4};
    // std::vector<double> coeffs; // per-leaf coefficients (kept for API compatibility)

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
        //rebuild_connectivity_active_bases(); // geometry known -> can build physical coords

        post_topology_update();

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
        // f.dofs_per_cell = (b == Basis::Q1_Quad4 ? 4 : (b == Basis::Serendipity8 ? 8 : 9));
        // f.coeffs.resize(leaf_count() * f.dofs_per_cell, 0.0);


        // Ensure registry exists so we can size nodal storage
        _activeBases.insert(b);
        //rebuild_connectivity_active_bases();

        post_topology_update();
        const auto& R = _basisReg[(int)b];
        f.nodal.assign(R.nodes.size(), 0.0);

        // _fields.push_back(std::move(f));
        //_activeBases.insert(b);   // mark basis as active

        _fields.push_back(std::move(f));
        return (u32)_fields.size() - 1;
      }

      Field& field(u32 fid) {
        return _fields[fid];
      }
      const Field& field(u32 fid) const {
        return _fields[fid];
      }

      // void resize_fields_to_leaves() {
      //   const u32 nL = leaf_count();
      //   for (auto& f : _fields) f.coeffs.resize(size_t(nL)*f.dofs_per_cell);
      // }
      //
      // double* leaf_coeff_ptr(u32 fid, u32 leaf_pos) {
      //   Field& f = _fields[fid];
      //   return f.coeffs.data() + size_t(leaf_pos) * f.dofs_per_cell;
      // }
      // const double* leaf_coeff_ptr(u32 fid, u32 leaf_pos) const {
      //   const Field& f = _fields[fid];
      //   return f.coeffs.data() + size_t(leaf_pos) * f.dofs_per_cell;
      // }

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
        //rebuild_connectivity_active_bases();
        //resize_fields_to_leaves();

        //rebuild_connectivity_active_bases();
        //resize_fields_to_nodes();

        post_topology_update();

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

          //rebuild_connectivity_active_bases();
          //resize_fields_to_nodes();
          post_topology_update();
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
        // std::vector<double> pointData;
        // std::vector<double> cellData;
        // if (cell_centered) {
        //   cellData.reserve(numCells);
        //   for (size_t e = 0; e < numCells; ++e) {
        //     const double *c = fld.coeffs.data() + e * fld.dofs_per_cell;
        //     double v = 0.0;
        //     for (int i = 0; i < fld.dofs_per_cell; ++i) v += c[i];
        //     cellData.push_back(v / fld.dofs_per_cell);
        //   }
        // }
        // else {
        //   pointData.resize(R.nodes.size(), 0.0);
        //   std::vector<int> counts(R.nodes.size(), 0);
        //   for (size_t e = 0; e < numCells; ++e) {
        //     const auto &conn = R.elem2glob[e];
        //     const double *c = fld.coeffs.data() + e * fld.dofs_per_cell;
        //     for (size_t i = 0; i < conn.size(); ++i) {
        //       int gid = conn[i];
        //       pointData[gid] += c[i];
        //       counts[gid]++;
        //     }
        //   }
        //   for (size_t i = 0; i < pointData.size(); ++i)
        //     if (counts[i] > 0) pointData[i] /= counts[i];
        // }

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


        // const u32 leaf_pos = leaf_position(leaf_node_idx);
        // if (leaf_pos == npos32) return false;

        // 3) access field + coefficients
        const Field& f = _fields[fid];
        //const double* c = leaf_coeff_ptr(fid, leaf_pos);
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
          Shapes::Serendipity8(shat, that, N8);
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

        // shat = N*xi + N - (2*ix + 1)
        // that = N*eta + N - (2*iy + 1)
        shat = N * xi  + N - double((ix << 1) + 1);
        that = N * eta + N - double((iy << 1) + 1);
      }



// Map a leaf's node index (into _tree_nodes) to its position in _leaves.
      // u32 leaf_position(u32 leaf_node_idx) const {
      //   for (u32 i = 0; i < (u32)_leaves.size(); ++i)
      //     if (_leaves[i] == leaf_node_idx) return i;
      //   return npos32;
      // }

      inline bool locate_leaf_on_parent_and_ref(double xi, double eta,
                                                u32& leaf_node_idx,
                                                double& shat, double& that) const {
        // 1) Find the leaf containing (xi,eta) in parent space.
        leaf_node_idx = locate_leaf_on_parent(xi, eta);
        if (leaf_node_idx == npos32) return false;

        // 2) Convert to local [-1,1]^2 using level+Morton (no divisions).
        local_ref_fast(leaf_node_idx, xi, eta, shat, that);

        // Optional: tiny clamp to counter FP jitter at boundaries.
        // const double eps = 1e-12;
        // shat = std::max(-1.0 - eps, std::min(1.0 + eps, shat));
        // that = std::max(-1.0 - eps, std::min(1.0 + eps, that));

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
      }


// // --- public: refine so that each physical point is inside a leaf at least 'maxDepthTarget' ---
//       void refine_to_contain_points(const std::vector<std::array<double, 2>>& pts,
//                                     u32 maxDepthTarget) {
//         if (pts.empty()) return;
//         require_geometry();
//
//         maxDepthTarget = std::min(maxDepthTarget, _maxDepth);
//
//         for (const auto& p : pts) {
//           double xi, eta;
//           if (!inverse_map_quad9(p[0], p[1], xi, eta)) continue;
//
//           u32 node_idx = locate_leaf_on_parent(xi, eta);
//           // Refine until we reach the requested depth (or maxDepth)
//           while (node_idx != npos32 &&
//                  _tree_nodes[node_idx].is_leaf &&
//                  _tree_nodes[node_idx].level < maxDepthTarget) {
//             if (!refine_leaf_and_update_leaves(node_idx)) break;
//             // after refinement, locate again (the point is in one of the new children)
//             node_idx = locate_leaf_on_parent(xi, eta);
//           }
//         }
//
//         // refresh connectivity/fields after topology changes
//         rebuild_connectivity_active_bases();
//         resize_fields_to_leaves();
//       }


      // ---- Public: enforce 1-irregularity ---------------------------------
    public:
      // Enforce 1-irregularity: adjacent leaves may differ by at most 1 level.
      // This refines only the coarser side until the constraint is satisfied.
      void enforce_balance() {
        bool changed = true;
        while (changed) {
          changed = false;

          // work on a snapshot to avoid invalidating indices while refining
          const std::vector<u32> leaves_snapshot = _leaves;

          for (u32 leaf_node_idx : leaves_snapshot) {
            if (_tree_nodes[leaf_node_idx].is_leaf == false) continue; // may have been refined already
            const u32 lev = _tree_nodes[leaf_node_idx].level;

            for (int dir = 0; dir < 4; ++dir) {
              u32 nb = neighbor_leaf_by_dir(leaf_node_idx, dir);
              if (nb == npos32) continue;
              if (_tree_nodes[nb].is_leaf == false) continue; // very unlikely, but safe

              const u32 lev_nb = _tree_nodes[nb].level;

              // If neighbor is too coarse relative to this leaf, refine neighbor
              if (lev > lev_nb + 1) {
                changed |= refine_leaf_once(nb);
              }
              // Symmetric check (in case we iterate this leaf before the finer neighbor does)
              else if (lev_nb > lev + 1) {
                changed |= refine_leaf_once(leaf_node_idx);
              }
            }
          }
        }

        // Rebuild once at the end

        post_topology_update();

        // rebuild_connectivity_active_bases();
        // resize_fields_to_nodes();
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
          // enforce 1-level irregularity
          enforce_balance();
        }


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

      // ------------------------------
      // Public read-only registry API
      // ------------------------------
      // Returns the unique global FEM nodes (per-basis), with parent & physical coords.
      const std::vector<FEMNode>& basis_nodes(Basis b) const {
        return _basisReg[(int)b].nodes;
      }
      // Returns the element-to-global connectivity (per-basis).
      const std::vector<std::vector<int>>& basis_connectivity(Basis b) const {
        return _basisReg[(int)b].elem2glob;
      }


      // ---- Private helpers -------------------------------------------------
    private:

      // // Rebuild field fid on *this* from source tree 'src' by sampling at parent nodes
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


// ======= internal helpers =======
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

        // rebuild_connectivity_active_bases();
        // resize_fields_to_nodes();
      }

      // // In class QuadTree2D (public:)
      // const std::vector<FEMNode>& basis_nodes(Basis b) const {
      //   return _basisReg[(int)b].nodes;
      // }
      //
      // const std::vector<std::vector<int>>& basis_connectivity(Basis b) const {
      //   return _basisReg[(int)b].elem2glob;
      // }


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
        const u32 gridN = (1u << _maxDepth);
        long long ix_l = llround((xi + 1.0) * gridN / 2.0);
        long long iy_l = llround((eta + 1.0) * gridN / 2.0);
        if (ix_l < 0) ix_l = 0;
        if (iy_l < 0) iy_l = 0;
        if (ix_l > (long long)gridN) ix_l = (long long)gridN;
        if (iy_l > (long long)gridN) iy_l = (long long)gridN;
        u32 ix = (u32)ix_l;
        u32 iy = (u32)iy_l;

        u64 key = interleave2(ix, iy);

        auto it = R.nodeMap.find(key);
        if (it != R.nodeMap.end()) return it->second;

        int gid = (int)R.nodes.size();
        FEMNode node{gid, {xi, eta}, parent_to_physical(xi, eta)};
        R.nodes.push_back(node);
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

        // // gather point data from per-leaf coeffs onto unique global nodes
        // std::vector<double> pointData(R.nodes.size(), 0.0);
        // std::vector<char>   filled(R.nodes.size(), 0);
        //
        // const int ndof = (int)fld.dofs_per_cell;
        // for (u32 e = 0; e < leaf_count(); ++e) {
        //   const double* c = fld.coeffs.data() + size_t(e) * ndof;
        //   const auto& conn = R.elem2glob[e];
        //   for (int a = 0; a < ndof; ++a) {
        //     int gid = conn[a];
        //     pointData[gid] = c[a];
        //     filled[gid] = 1;
        //   }
        // }

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

// --- private: Q2 inverse map (Newton in parent space) ---
      bool inverse_map_quad9(double x, double y, double& xi, double& eta,
                             int maxIts = 25, double tol = 1e-12) const {
        require_geometry();
        // start at center
        xi = 0.0;
        eta = 0.0;

        for (int it = 0; it < maxIts; ++it) {
          double N9[9], dNdxi[9], dNdeta[9];
          Quad9Shape::N(xi, eta, N9);
          Quad9Shape::dN(xi, eta, dNdxi, dNdeta);

          const auto& X = _X;
          const auto& Y = _Y;

          // current mapped point and Jacobian
          double Xp = 0.0, Yp = 0.0, dXdxi = 0.0, dXdeta = 0.0, dYdxi = 0.0, dYdeta = 0.0;
          for (int a = 0; a < 9; ++a) {
            Xp     += N9[a]    * X[a];
            Yp     += N9[a]    * Y[a];
            dXdxi  += dNdxi[a] * X[a];
            dXdeta += dNdeta[a] * X[a];
            dYdxi  += dNdxi[a] * Y[a];
            dYdeta += dNdeta[a] * Y[a];
          }

          double rx = Xp - x, ry = Yp - y;
          double nrm = std::sqrt(rx * rx + ry * ry);
          if (nrm < tol) return true;

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

    private:
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


      // in QuadTree2D
      mutable std::vector<u32> _node2leafpos; // size == _tree_nodes.size(), npos32 default

      inline void rebuild_leafpos_lookup() {
        _node2leafpos.assign(_tree_nodes.size(), npos32);
        for (u32 i = 0; i < _leaves.size(); ++i) _node2leafpos[_leaves[i]] = i;
      }



  };

} // namespace fem

























// #pragma once
// #include <vector>
// #include <array>
// #include <cstdint>
// #include <limits>
// #include <cmath>
// #include <cassert>
// #include <functional>
// #include <fstream>
// #include <iomanip>
// #include <algorithm>
//
// #include "Encoder.hpp"
//
// /*
//   QuadTree2D.hpp (ASCII VTU)
//   - 2D quadtree living on parent domain [-1,1]^2
//   - Quad9 (Q2) geometry map for the coarse/root element
//   - Supports Q1, Serendipity8, and Q2 fields per leaf
//   - Refinement/coarsening, field evaluation, and snapshot-safe coarsening
//   - Parent-space helpers to avoid repeated inverse maps
// */
//
// namespace fem {
//
//   using u32 = uint32_t;
//   using u64 = uint64_t;
//   constexpr u32 npos32 = std::numeric_limits<u32>::max();
//
//   // Supported element bases
//   enum class Basis : uint8_t { Q1_Quad4 = 0, Serendipity8 = 1, Q2_Quad9 = 2 };
//
//   //----------------------------------------
//   // Quad9 shapes (parent space [-1,1]^2)
//   //----------------------------------------
//   struct Quad9Shape {
//     // 1D Q2 Lagrange at s in {-1,0,1}: returns basis and derivatives
//     static inline void q2_1d(double s, double& L0, double& L1, double& L2,
//                              double& dL0, double& dL1, double& dL2) {
//       L0 = 0.5 * s * (s - 1.0);   // 0.5*(s^2 - s)
//       L1 = 1.0 - s * s;           // 1 - s^2
//       L2 = 0.5 * s * (s + 1.0);   // 0.5*(s^2 + s)
//       dL0 = s - 0.5;
//       dL1 = -2.0 * s;
//       dL2 = s + 0.5;
//     }
//
//     // Compute Q2 shape functions at (xi,eta) in parent space
//     static inline void N(double xi, double eta, double N9[9]) {
//       double Lx[3], Ly[3], dx[3], dy[3];
//       q2_1d(xi,  Lx[0], Lx[1], Lx[2], dx[0], dx[1], dx[2]);
//       q2_1d(eta, Ly[0], Ly[1], Ly[2], dy[0], dy[1], dy[2]);
//       (void)dx;
//       (void)dy;
//
//       // corners
//       N9[0] = Lx[0] * Ly[0]; // (-1,-1)
//       N9[1] = Lx[2] * Ly[0]; // ( 1,-1)
//       N9[2] = Lx[2] * Ly[2]; // ( 1, 1)
//       N9[3] = Lx[0] * Ly[2]; // (-1, 1)
//       // mids
//       N9[4] = Lx[1] * Ly[0]; // (0,-1)
//       N9[5] = Lx[2] * Ly[1]; // (1,0)
//       N9[6] = Lx[1] * Ly[2]; // (0,1)
//       N9[7] = Lx[0] * Ly[1]; // (-1,0)
//       // center
//       N9[8] = Lx[1] * Ly[1]; // (0,0)
//     }
//
//     // Compute derivative wrt xi and eta for Q2 at (xi,eta)
//     static inline void dN(double xi, double eta, double dN_dxi[9], double dN_deta[9]) {
//       double Lx[3], Ly[3], dx[3], dy[3];
//       q2_1d(xi,  Lx[0], Lx[1], Lx[2], dx[0], dx[1], dx[2]);
//       q2_1d(eta, Ly[0], Ly[1], Ly[2], dy[0], dy[1], dy[2]);
//
//       dN_dxi[0]  = dx[0] * Ly[0];
//       dN_deta[0] = Lx[0] * dy[0];
//       dN_dxi[1]  = dx[2] * Ly[0];
//       dN_deta[1] = Lx[2] * dy[0];
//       dN_dxi[2]  = dx[2] * Ly[2];
//       dN_deta[2] = Lx[2] * dy[2];
//       dN_dxi[3]  = dx[0] * Ly[2];
//       dN_deta[3] = Lx[0] * dy[2];
//       dN_dxi[4]  = dx[1] * Ly[0];
//       dN_deta[4] = Lx[1] * dy[0];
//       dN_dxi[5]  = dx[2] * Ly[1];
//       dN_deta[5] = Lx[2] * dy[1];
//       dN_dxi[6]  = dx[1] * Ly[2];
//       dN_deta[6] = Lx[1] * dy[2];
//       dN_dxi[7]  = dx[0] * Ly[1];
//       dN_deta[7] = Lx[0] * dy[1];
//       dN_dxi[8]  = dx[1] * Ly[1];
//       dN_deta[8] = Lx[1] * dy[1];
//     }
//   };
//
//   //----------------------------------------
//   // Q1 (Quad4) and Serendipity8 shapes on [-1,1]^2
//   //----------------------------------------
//   struct Shapes {
//     // Compute Q1 bilinear shape functions at (xi,eta)
//     static inline void Q1(double xi, double eta, double N4[4]) {
//       const double a = 0.25 * (1 - xi), b = 0.25 * (1 + xi);
//       N4[0] = a * (1 - eta); // (-1,-1)
//       N4[1] = b * (1 - eta); // ( 1,-1)
//       N4[2] = b * (1 + eta); // ( 1, 1)
//       N4[3] = a * (1 + eta); // (-1, 1)
//     }
//
//     // Compute serendipity-8 shape functions at (xi,eta)
//     static inline void Serendipity8(double xi, double eta, double N8[8]) {
//       const double xm = xi, em = eta;
//       N8[0] = 0.25 * (1 - xm) * (1 - em) * (-xm - em - 1);
//       N8[1] = 0.25 * (1 + xm) * (1 - em) * (xm - em - 1);
//       N8[2] = 0.25 * (1 + xm) * (1 + em) * (xm + em - 1);
//       N8[3] = 0.25 * (1 - xm) * (1 + em) * (-xm + em - 1);
//       N8[4] = 0.5 * (1 - xm * xm) * (1 - em);
//       N8[5] = 0.5 * (1 + xm) * (1 - em * em);
//       N8[6] = 0.5 * (1 - xm * xm) * (1 + em);
//       N8[7] = 0.5 * (1 - xm) * (1 - em * em);
//     }
//   };
//
//   //----------------------------------------
//   // Morton helpers (2D)
//   //----------------------------------------
//   // Interleave bits of (x,y) into a 64-bit Morton code
//   static inline u64 interleave2(u32 x, u32 y) {
//     u64 xx = x, yy = y;
//     xx = (xx | (xx << 16)) & 0x0000FFFF0000FFFFULL;
//     xx = (xx | (xx << 8)) & 0x00FF00FF00FF00FFULL;
//     xx = (xx | (xx << 4)) & 0x0F0F0F0F0F0F0F0FULL;
//     xx = (xx | (xx << 2)) & 0x3333333333333333ULL;
//     xx = (xx | (xx << 1)) & 0x5555555555555555ULL;
//
//     yy = (yy | (yy << 16)) & 0x0000FFFF0000FFFFULL;
//     yy = (yy | (yy << 8)) & 0x00FF00FF00FF00FFULL;
//     yy = (yy | (yy << 4)) & 0x0F0F0F0F0F0F0F0FULL;
//     yy = (yy | (yy << 2)) & 0x3333333333333333ULL;
//     yy = (yy | (yy << 1)) & 0x5555555555555555ULL;
//
//     return (yy << 1) | xx;
//   }
//
//   // Reverse Morton code to (x,y)
//   static inline void deinterleave2(u64 m, u32& x, u32& y) {
//     u64 xx = m, yy = m >> 1;
//     xx &= 0x5555555555555555ULL;
//     yy &= 0x5555555555555555ULL;
//
//     xx = (xx | (xx >> 1)) & 0x3333333333333333ULL;
//     xx = (xx | (xx >> 2)) & 0x0F0F0F0F0F0F0F0FULL;
//     xx = (xx | (xx >> 4)) & 0x00FF00FF00FF00FFULL;
//     xx = (xx | (xx >> 8)) & 0x0000FFFF0000FFFFULL;
//     xx = (xx | (xx >> 16)) & 0x00000000FFFFFFFFULL;
//
//     yy = (yy | (yy >> 1)) & 0x3333333333333333ULL;
//     yy = (yy | (yy >> 2)) & 0x0F0F0F0F0F0F0F0FULL;
//     yy = (yy | (yy >> 4)) & 0x00FF00FF00FF00FFULL;
//     yy = (yy | (yy >> 8)) & 0x0000FFFF0000FFFFULL;
//     yy = (yy | (yy >> 16)) & 0x00000000FFFFFFFFULL;
//
//     x = (u32)xx;
//     y = (u32)yy;
//   }
//
//   //----------------------------------------
//   // Quadtree node
//   //----------------------------------------
//   struct QuadNode {
//     u32 parent{npos32};
//     std::array<u32, 4> child{npos32, npos32, npos32, npos32};
//     u32 level{0};
//     u64 morton{0};
//     uint8_t flags{0};
//
//     // True if node has no children
//     inline bool is_leaf() const {
//       return child[0] == npos32;
//     }
//   };
//
//   //----------------------------------------
//   // Field storage (per-leaf coefficients)
//   //----------------------------------------
//   struct Field {
//     Basis basis{Basis::Q1_Quad4};
//     u32   dofs_per_cell{4};
//     std::vector<double> coeffs;
//   };
//
//   //----------------------------------------
//   // QuadTree2D: quadtree over [-1,1]^2 with Quad9 geometry
//   //----------------------------------------
//   class QuadTree2D {
//     public:
//       // Construct with maxDepth; tree starts as a single root
//       QuadTree2D(u32 maxDepth = 28)
//         : _maxDepth(maxDepth) {
//         _nodes.reserve(1024);
//         _alive.reserve(1024);
//         _free.reserve(256);
//         _root = alloc_node();
//         _nodes[_root].level = 0;
//         _nodes[_root].morton = 0;
//       }
//
//       // Construct with maxDepth and initial minDepth (soft floor)
//       QuadTree2D(u32 maxDepth, u32 minDepth)
//         : _maxDepth(maxDepth), _minDepth(std::min(minDepth, maxDepth)) {
//         _nodes.reserve(1024);
//         _alive.reserve(1024);
//         _free.reserve(256);
//         _root = alloc_node();
//         _nodes[_root].level = 0;
//         _nodes[_root].morton = 0;
//       }
//
//       // Set minimum allowed depth for automatic floor
//       void set_min_depth(u32 d) {
//         _minDepth = std::min(d, _maxDepth);
//       }
//
//       // Get current minimum depth floor
//       u32  min_depth() const {
//         return _minDepth;
//       }
//
//       // Allow/disallow coarsening below min depth
//       void set_allow_coarsen_below_min(bool v) {
//         _allowCoarsenBelowMinDepth = v;
//       }
//
//       // Query flag for coarsening below min depth
//       bool allow_coarsen_below_min() const {
//         return _allowCoarsenBelowMinDepth;
//       }
//
//       // Set physical coordinates (Q2 geometry) for the coarse/root element
//       void set_physical_quad9(const double x9[9], const double y9[9]) {
//         for (int i = 0; i < 9; ++i) {
//           _X[i] = x9[i];
//           _Y[i] = y9[i];
//         }
//         _geom_ready = true;
//       }
//
//       // Assert geometry is defined
//       inline void require_geometry() const {
//         assert(_geom_ready && "Call set_physical_quad9(...) before geometric operations.");
//       }
//
//       // Refine a leaf into 4 children (Morton ordering), returns false if not allowed
//       bool refine(u32 i) {
//         if (i == npos32 || !_alive[i]) return false;
//         if (!_nodes[i].is_leaf() || _nodes[i].level >= _maxDepth) return false;
//
//         const u32 parent_level  = _nodes[i].level;
//         const u64 parent_morton = _nodes[i].morton;
//
//         for (int k = 0; k < 4; ++k) {
//           const u32 c = alloc_node();
//           _nodes[i].child[k] = c;
//           QuadNode& cn = _nodes[c];
//           cn.parent = i;
//           cn.level  = parent_level + 1;
//           cn.morton = (parent_morton << 2)
//                       | (u64)((k & 2) ? 2 : 0)
//                       | (u64)((k & 1) ? 1 : 0);
//         }
//         _leaf_dirty = true;
//         return true;
//       }
//
//       // Coarsen a parent whose 4 children are leaves; returns false if not possible
//       bool coarsen(u32 i) {
//         if (i == npos32 || !_alive[i]) return false;
//         if (_nodes[i].is_leaf()) return false;
//         if (_nodes[i].level <= _minDepth && !_allowCoarsenBelowMinDepth) return false;
//
//         std::array<u32, 4> ch = _nodes[i].child;
//
//         for (int k = 0; k < 4; ++k) {
//           const u32 c = ch[k];
//           if (c == npos32) return false;
//           if (!_alive[c])  return false;
//           if (!_nodes[c].is_leaf()) return false;
//           if (_nodes[c].parent != i) return false;
//         }
//
//         for (int k = 0; k < 4; ++k) free_node(ch[k]);
//         _nodes[i].child = {npos32, npos32, npos32, npos32};
//         _leaf_dirty = true;
//         return true;
//       }
//
//       // Build/refresh the compact list of leaves (indices into _nodes)
//       const std::vector<u32>& leaves() const {
//         if (_leaf_dirty) {
//           _leaves.clear();
//           _leaves.reserve(_nodes.size());
//           for (u32 i = 0; i < (u32)_nodes.size(); ++i)
//             if (_alive[i] && _nodes[i].is_leaf())
//               _leaves.push_back(i);
//
//           _leaf_pos.assign(_nodes.size(), npos32);
//           for (u32 pos = 0; pos < (u32)_leaves.size(); ++pos) _leaf_pos[_leaves[pos]] = pos;
//           _leaf_dirty = false;
//         }
//         return _leaves;
//       }
//
//       // Perform one refinement pass using a predicate on current leaves
//       template<class Pred>
//       std::size_t refine_pass(Pred&& should_refine, Basis probe_basis = Basis::Q2_Quad9) {
//         const auto& L = leaves();
//         std::vector<u32> to_refine;
//         to_refine.reserve(L.size());
//
//         for (u32 k = 0; k < (u32)L.size(); ++k) {
//           const u32 leaf = L[k];
//
//           std::vector<std::array<double, 2>> pts_xi, pts_xy;
//           std::vector<std::array<double, 9>> Nvals;
//           leaf_parent_nodes(probe_basis, leaf, pts_xi);
//           leaf_physical_nodes(probe_basis, leaf, pts_xy);
//           quad9_shapes_at(pts_xi, Nvals);
//
//           if (should_refine(leaf, pts_xi, pts_xy, Nvals)) to_refine.push_back(leaf);
//         }
//
//         std::size_t refined = 0;
//         for (u32 leaf : to_refine) if (refine(leaf)) ++refined;
//         (void)leaves();
//         return refined;
//       }
//
//       // Return axis-aligned bounds of a node in parent coordinates
//       inline void leaf_bounds(u32 i, double& xi0, double& eta0, double& xi1, double& eta1) const {
//         const QuadNode& n = _nodes[i];
//         u32 ix, iy;
//         deinterleave2(n.morton, ix, iy);
//         const double N = double(1u << n.level);
//         const double dx = 2.0 / N, dy = 2.0 / N;
//         xi0  = -1.0 + ix * dx;
//         eta0 = -1.0 + iy * dy;
//         xi1  = xi0 + dx;
//         eta1 = eta0 + dy;
//       }
//
//       // Locate leaf containing (xi,eta) in the parent space [-1,1]^2 using Morton codes
//       u32 locate_leaf_on_parent(double xi, double eta) const {
//         if (xi < -1.0 || xi > 1.0 || eta < -1.0 || eta > 1.0) return npos32;
//
//         // Scale to integer grid at max depth
//         const double scale = double(1u << _maxDepth) / 2.0; // maps [-1,1] -> [0,2^maxDepth)
//         u32 ix = std::min<u32>((u32)((xi + 1.0) * scale), (1u << _maxDepth) - 1);
//         u32 iy = std::min<u32>((u32)((eta + 1.0) * scale), (1u << _maxDepth) - 1);
//
//         u32 node = _root;
//         for (u32 level = 0; level < _maxDepth; ++level) {
//           const auto& n = _nodes[node];
//           if (n.is_leaf()) return node;
//
//           // pick bit for this level (from MSB downwards)
//           const int shift = _maxDepth - level - 1;
//           const int qx = (ix >> shift) & 1u;
//           const int qy = (iy >> shift) & 1u;
//           const int q  = (qy << 1) | qx; // 0=LL,1=LR,2=UL,3=UR
//
//           const u32 c = n.child[q];
//           if (c == npos32) return node;
//           node = c;
//         }
//         return node;
//       }
//
//
//       // Locate leaf containing (xi,eta) in parent space [-1,1]^2
//       // and compute leaf-local coordinates (shat,that) in [-1,1]^2.
//       // Returns false if (xi,eta) outside domain.
//       bool locate_leaf_on_parent_and_ref(double xi, double eta,
//                                          u32& leaf,
//                                          double& shat, double& that) const {
//         // Reject outside parent element
//         if (xi < -1.0 || xi > 1.0 || eta < -1.0 || eta > 1.0) {
//           leaf = npos32;
//           return false;
//         }
//
//         // Step 1. Find the leaf index (fast Morton-based descent)
//         leaf = locate_leaf_on_parent(xi, eta);
//         if (leaf == npos32) return false;
//
//         // Step 2. Compute leaf bounds once
//         double xi0, eta0, xi1, eta1;
//         leaf_bounds(leaf, xi0, eta0, xi1, eta1);
//
//         // Step 3. Map to leaf reference coordinates ([-1,1]²)
//         to_leaf_ref(xi, eta, xi0, eta0, xi1, eta1, shat, that);
//
//         return true;
//       }
//
//
//
//       // Locate leaf containing physical point (x,y)
//       u32 locate_leaf_on_physical(double x, double y) const {
//         double xi, eta;
//         if (!inverse_map_quad9(x, y, xi, eta)) return npos32;
//         return locate_leaf_on_parent(xi, eta);
//       }
//
//       // Inverse map physical (x,y) -> (xi,eta) using Newton on Q2 geometry
//       bool inverse_map_quad9(double x, double y, double& xi, double& eta,
//                              double tol = 1e-12, int maxit = 25) const {
//         require_geometry();
//         xi  = 0.0;
//         eta = 0.0;
//         for (int it = 0; it < maxit; ++it) {
//           double N9[9], dNdxi[9], dNdeta[9];
//           Quad9Shape::N(xi, eta, N9);
//           Quad9Shape::dN(xi, eta, dNdxi, dNdeta);
//
//           double Xh = 0.0, Yh = 0.0, dXdxi = 0.0, dXdeta = 0.0, dYdxi = 0.0, dYdeta = 0.0;
//           for (int a = 0; a < 9; ++a) {
//             Xh += N9[a] * _X[a];
//             Yh += N9[a] * _Y[a];
//             dXdxi  += dNdxi[a] * _X[a];
//             dXdeta += dNdeta[a] * _X[a];
//             dYdxi  += dNdxi[a] * _Y[a];
//             dYdeta += dNdeta[a] * _Y[a];
//           }
//
//           const double rx = Xh - x, ry = Yh - y;
//           const double norm = std::sqrt(rx * rx + ry * ry);
//           if (norm < tol) return true;
//
//           const double J00 = dXdxi, J01 = dXdeta;
//           const double J10 = dYdxi, J11 = dYdeta;
//           const double detJ = J00 * J11 - J01 * J10;
//           if (std::abs(detJ) <= 1e-30) return false;
//
//           const double dxi  = (J11 * rx - J01 * ry) / detJ;
//           const double deta = (-J10 * rx + J00 * ry) / detJ;
//           xi  -= dxi;
//           eta -= deta;
//           if (std::abs(dxi) + std::abs(deta) < tol * 1e-6) return true;
//         }
//         return false;
//       }
//
//       // Map parent (xi,eta) to leaf-local (shat,that) in [-1,1]^2 given leaf bounds
//       static inline void to_leaf_ref(double xi, double eta,
//                                      double xi0, double eta0, double xi1, double eta1,
//                                      double& shat, double& that) {
//         const double cx = 0.5 * (xi0 + xi1), hx = 0.5 * (xi1 - xi0);
//         const double cy = 0.5 * (eta0 + eta1), hy = 0.5 * (eta1 - eta0);
//         shat = (xi  - cx) / hx;
//         that = (eta - cy) / hy;
//       }
//
//       // Add a field and return its id; caller later resizes coeffs appropriately
//       u32 add_field(Basis b) {
//         Field f;
//         f.basis = b;
//         f.dofs_per_cell = (b == Basis::Q1_Quad4) ? 4 : (b == Basis::Serendipity8 ? 8 : 9);
//         _fields.push_back(std::move(f));
//         return (u32)(_fields.size() - 1);
//       }
//
//       // Mutable access to field by id
//       Field& field(u32 fid) {
//         return _fields[fid];
//       }
//
//       // Const access to field by id
//       const Field& field(u32 fid) const {
//         return _fields[fid];
//       }
//
//       // Get pointer into leaf's coefficient block for field fid (mutable)
//       double* leaf_coeff_ptr(u32 fid, u32 leaf_pos) {
//         Field& f = _fields[fid];
//         return f.coeffs.data() + size_t(leaf_pos) * f.dofs_per_cell;
//       }
//
//       // Get pointer into leaf's coefficient block for field fid (const)
//       const double* leaf_coeff_ptr(u32 fid, u32 leaf_pos) const {
//         const Field& f = _fields[fid];
//         return f.coeffs.data() + size_t(leaf_pos) * f.dofs_per_cell;
//       }
//
//       // Evaluate field fid at physical (x,y); returns false if outside or inverse fails
//       bool evaluate_field_on_physical(u32 fid, double x, double y, double& value) const {
//         require_geometry();
//         double xi, eta;
//         if (!inverse_map_quad9(x, y, xi, eta)) return false;
//         const u32 leaf = locate_leaf_on_parent(xi, eta);
//         if (leaf == npos32) return false;
//
//         double xi0, eta0, xi1, eta1;
//         leaf_bounds(leaf, xi0, eta0, xi1, eta1);
//         double shat, that;
//         to_leaf_ref(xi, eta, xi0, eta0, xi1, eta1, shat, that);
//
//         const u32 leaf_pos = leaf_position(leaf);
//         if (leaf_pos == npos32) return false;
//
//         const Field& f = _fields[fid];
//         const double* c = leaf_coeff_ptr(fid, leaf_pos);
//
//         switch (f.basis) {
//         case Basis::Q1_Quad4: {
//           double N4[4];
//           Shapes::Q1(shat, that, N4);
//           value = N4[0] * c[0] + N4[1] * c[1] + N4[2] * c[2] + N4[3] * c[3];
//         }
//         break;
//         case Basis::Serendipity8: {
//           double N8[8];
//           Shapes::Serendipity8(shat, that, N8);
//           double v = 0.0;
//           for (int i = 0; i < 8; ++i) v += N8[i] * c[i];
//           value = v;
//         }
//         break;
//         case Basis::Q2_Quad9: {
//           double N9[9];
//           Quad9Shape::N(shat, that, N9);
//           double v = 0.0;
//           for (int i = 0; i < 9; ++i) v += N9[i] * c[i];
//           value = v;
//         }
//         break;
//         }
//         return true;
//       }
//
//       // Return current leaf count
//       u32 leaf_count() {
//         return (u32)leaves().size();
//       }
//
//       // Resize all field coefficient vectors to match current leaves (keeps values where possible)
//       void resize_fields_to_leaves() {
//         const auto nL = leaf_count();
//         for (auto& f : _fields) f.coeffs.resize(size_t(nL) * f.dofs_per_cell);
//       }
//
//       // Expose leaf indices (const reference to compact leaf list)
//       const std::vector<u32>& leaf_indices() const {
//         return leaves();
//       }
//
//       // Fill vector with parent coordinates of interpolation nodes for given leaf+basis
//       void leaf_parent_nodes(Basis basis, u32 leaf,
//                              std::vector<std::array<double, 2>>& out_pts) const {
//         double xi0, eta0, xi1, eta1;
//         leaf_bounds(leaf, xi0, eta0, xi1, eta1);
//         const double cx = 0.5 * (xi0 + xi1), hx = 0.5 * (xi1 - xi0);
//         const double cy = 0.5 * (eta0 + eta1), hy = 0.5 * (eta1 - eta0);
//
//         auto toParent = [&](double s, double t) {
//           return std::array<double, 2> { cx + s * hx, cy + t * hy };
//         };
//
//         out_pts.clear();
//         switch (basis) {
//         case Basis::Q1_Quad4: {
//           static const double S[4][2] = {{-1, -1}, {+1, -1}, {+1, +1}, {-1, +1}};
//           out_pts.reserve(4);
//           for (int i = 0; i < 4; ++i) out_pts.push_back(toParent(S[i][0], S[i][1]));
//         }
//         break;
//         case Basis::Serendipity8: {
//           static const double S[8][2] = {
//             {-1, -1}, {+1, -1}, {+1, +1}, {-1, +1},
//             {0, -1}, {+1, 0}, {0, +1}, {-1, 0}
//           };
//           out_pts.reserve(8);
//           for (int i = 0; i < 8; ++i) out_pts.push_back(toParent(S[i][0], S[i][1]));
//         }
//         break;
//         case Basis::Q2_Quad9: {
//           static const double S[9][2] = {
//             {-1, -1}, {+1, -1}, {+1, +1}, {-1, +1},
//             {0, -1}, {+1, 0}, {0, +1}, {-1, 0}, {0, 0}
//           };
//           out_pts.reserve(9);
//           for (int i = 0; i < 9; ++i) out_pts.push_back(toParent(S[i][0], S[i][1]));
//         }
//         break;
//         }
//       }
//
//       // Map given parent nodes to physical coordinates via Q2 geometry
//       void leaf_physical_nodes(Basis basis, u32 leaf,
//                                std::vector<std::array<double, 2>>& out_xy) const {
//         require_geometry();
//         std::vector<std::array<double, 2>> pts;
//         leaf_parent_nodes(basis, leaf, pts);
//         out_xy.resize(pts.size());
//         for (size_t i = 0; i < pts.size(); ++i) {
//           const double xi = pts[i][0], eta = pts[i][1];
//           double N9[9];
//           Quad9Shape::N(xi, eta, N9);
//           double X = 0.0, Y = 0.0;
//           for (int a = 0; a < 9; ++a) {
//             X += N9[a] * _X[a];
//             Y += N9[a] * _Y[a];
//           }
//           out_xy[i] = {X, Y};
//         }
//       }
//
//       // Evaluate Q2 shape arrays for a list of parent points
//       void quad9_shapes_at(const std::vector<std::array<double, 2>>& points,
//                            std::vector<std::array<double, 9>>& out_N) const {
//         out_N.resize(points.size());
//         for (size_t i = 0; i < points.size(); ++i) {
//           double N9[9];
//           Quad9Shape::N(points[i][0], points[i][1], N9);
//           for (int a = 0; a < 9; ++a) out_N[i][a] = N9[a];
//         }
//       }
//
//       // Refine all leaves where predicate is true
//       template<class Pred>
//       void refine_where(Pred&& should_refine, Basis probe_basis = Basis::Q2_Quad9) {
//         const auto& L = leaves();
//         std::vector<u32> to_refine;
//         to_refine.reserve(L.size());
//         for (u32 k = 0; k < (u32)L.size(); ++k) {
//           u32 leaf = L[k];
//           std::vector<std::array<double, 2>> pts_xi, pts_xy;
//           std::vector<std::array<double, 9>> Nvals;
//           leaf_parent_nodes(probe_basis, leaf, pts_xi);
//           leaf_physical_nodes(probe_basis, leaf, pts_xy);
//           quad9_shapes_at(pts_xi, Nvals);
//           if (should_refine(leaf, pts_xi, pts_xy, Nvals)) to_refine.push_back(leaf);
//         }
//         for (u32 leaf : to_refine) refine(leaf);
//         (void)leaves();
//       }
//
//       // Coarsen-pass over parents with 4 leaf children using provided predicate
//       template<class Pred>
//       std::size_t coarsen_pass(Pred&& should_coarsen,
//                                Basis probe_basis = Basis::Q2_Quad9) {
//         std::vector<u32> to_coarsen;
//         to_coarsen.reserve(_nodes.size() / 4 + 1);
//
//         for (u32 i = 0; i < (u32)_nodes.size(); ++i) {
//           if (!_alive[i]) continue;
//           const QuadNode& n = _nodes[i];
//           if (n.is_leaf()) continue;
//
//           bool ok = true;
//           for (int k = 0; k < 4; ++k) {
//             u32 c = n.child[k];
//             if (c == npos32 || !_alive[c] || !_nodes[c].is_leaf() || _nodes[c].parent != i) {
//               ok = false;
//               break;
//             }
//           }
//           if (!ok) continue;
//
//           std::vector<std::array<double, 2>> pts_xi, pts_xy;
//           std::vector<std::array<double, 9>> Nvals;
//           leaf_parent_nodes(probe_basis, i, pts_xi);
//           leaf_physical_nodes(probe_basis, i, pts_xy);
//           quad9_shapes_at(pts_xi, Nvals);
//
//           bool ok_min = !(n.level <= _minDepth && !_allowCoarsenBelowMinDepth);
//           if (ok_min && should_coarsen(i, n.level, pts_xi, pts_xy, Nvals)) {
//             to_coarsen.push_back(i);
//           }
//         }
//
//         std::size_t done = 0;
//         for (u32 p : to_coarsen) if (coarsen(p)) ++done;
//         (void)leaves();
//         return done;
//       }
//
//       // Multilevel refinement until convergence or pass limit
//       template<class Pred>
//       std::size_t adapt_refine_until(Pred&& should_refine,
//                                      Basis probe_basis = Basis::Q2_Quad9,
//                                      u32 max_passes = 10) {
//         ensure_min_depth();
//         std::size_t total_refined = 0;
//         for (u32 pass = 0; pass < max_passes; ++pass) {
//           const std::size_t r = refine_pass(should_refine, probe_basis);
//           total_refined += r;
//           if (r == 0) break;
//         }
//         return total_refined;
//       }
//
//       // Full adaptivity cycle: coarsen then refine for up to max_passes
//       template<class PredCoarsen, class PredRefine>
//       std::size_t adapt_cycle(PredCoarsen&& should_coarsen,
//                               PredRefine&& should_refine,
//                               Basis probe_basis = Basis::Q2_Quad9,
//                               u32 max_passes = 10) {
//         ensure_min_depth();
//         std::size_t total = 0;
//         for (u32 pass = 0; pass < max_passes; ++pass) {
//           std::size_t c = coarsen_pass(should_coarsen, probe_basis);
//           std::size_t r = refine_pass(should_refine, probe_basis);
//           total += c + r;
//           if (c == 0 && r == 0) break;
//         }
//         return total;
//       }
//
//
//
//
//       // Fill vector with parent-space coordinates of interpolation nodes (by basis) for leaf
//       void leaf_reference_nodes(Basis basis, u32 leaf,
//                                 std::vector<std::array<double, 2>>& xi) const {
//         xi.clear();
//         double xi0, eta0, xi1, eta1;
//         leaf_bounds(leaf, xi0, eta0, xi1, eta1);
//         switch (basis) {
//         case Basis::Q1_Quad4: {
//           xi.resize(4);
//           xi[0] = {xi0, eta0};
//           xi[1] = {xi1, eta0};
//           xi[2] = {xi1, eta1};
//           xi[3] = {xi0, eta1};
//         }
//         break;
//         case Basis::Serendipity8: {
//           xi.resize(8);
//           double xm = 0.5 * (xi0 + xi1), ym = 0.5 * (eta0 + eta1);
//           xi[0] = {xi0, eta0};
//           xi[1] = {xi1, eta0};
//           xi[2] = {xi1, eta1};
//           xi[3] = {xi0, eta1};
//           xi[4] = {xm,  eta0};
//           xi[5] = {xi1, ym };
//           xi[6] = {xm,  eta1};
//           xi[7] = {xi0, ym };
//         }
//         break;
//         case Basis::Q2_Quad9: {
//           xi.resize(9);
//           double xm = 0.5 * (xi0 + xi1), ym = 0.5 * (eta0 + eta1);
//           xi[0] = {xi0, eta0};
//           xi[1] = {xi1, eta0};
//           xi[2] = {xi1, eta1};
//           xi[3] = {xi0, eta1};
//           xi[4] = {xm,  eta0};
//           xi[5] = {xi1, ym };
//           xi[6] = {xm,  eta1};
//           xi[7] = {xi0, ym };
//           xi[8] = {xm,  ym };
//         }
//         break;
//         }
//       }
//
//       // Evaluate a field on a known leaf using (shat,that) in leaf-local [-1,1]^2
//       bool evaluate_field_on_leaf(u32 fid, u32 leaf, double shat, double that, double& value) const {
//         const u32 leaf_pos = leaf_position(leaf);
//         if (leaf_pos == npos32) return false;
//         const Field& f = _fields[fid];
//         const double* c = _fields[fid].coeffs.data() + size_t(leaf_pos) * f.dofs_per_cell;
//         switch (f.basis) {
//         case Basis::Q1_Quad4: {
//           double N4[4];
//           Shapes::Q1(shat, that, N4);
//           value = N4[0] * c[0] + N4[1] * c[1] + N4[2] * c[2] + N4[3] * c[3];
//         }
//         break;
//         case Basis::Serendipity8: {
//           double N8[8];
//           Shapes::Serendipity8(shat, that, N8);
//           double v = 0.0;
//           for (int i = 0; i < 8; ++i) v += N8[i] * c[i];
//           value = v;
//         }
//         break;
//         case Basis::Q2_Quad9: {
//           double N9[9];
//           Quad9Shape::N(shat, that, N9);
//           double v = 0.0;
//           for (int i = 0; i < 9; ++i) v += N9[i] * c[i];
//           value = v;
//         }
//         break;
//         }
//         return true;
//       }
//
//       // Evaluate field directly in parent coordinates (xi,eta)
//       // Uses fused locate+ref mapping for efficiency
//       bool evaluate_field_on_parent(u32 fid, double xi, double eta, double& value) const {
//         u32 leaf;
//         double shat, that;
//
//         // Step 1. Locate leaf and compute local ref coords in one shot
//         if (!locate_leaf_on_parent_and_ref(xi, eta, leaf, shat, that)) {
//           return false;
//         }
//
//         // Step 2. Map leaf index to position in coefficient storage
//         const u32 leaf_pos = leaf_position(leaf);
//         if (leaf_pos == npos32) return false;
//
//         // Step 3. Access field coefficients
//         const Field& f = _fields[fid];
//         const double* c = leaf_coeff_ptr(fid, leaf_pos);
//
//         // Step 4. Evaluate basis interpolation
//         switch (f.basis) {
//         case Basis::Q1_Quad4: {
//           double N4[4];
//           Shapes::Q1(shat, that, N4);
//           value = N4[0] * c[0] + N4[1] * c[1] + N4[2] * c[2] + N4[3] * c[3];
//         }
//         break;
//
//         case Basis::Serendipity8: {
//           double N8[8];
//           Shapes::Serendipity8(shat, that, N8);
//           double v = 0.0;
//           for (int i = 0; i < 8; ++i) v += N8[i] * c[i];
//           value = v;
//         }
//         break;
//
//         case Basis::Q2_Quad9: {
//           double N9[9];
//           Quad9Shape::N(shat, that, N9);
//           double v = 0.0;
//           for (int i = 0; i < 9; ++i) v += N9[i] * c[i];
//           value = v;
//         }
//         break;
//         }
//         return true;
//       }
//
//
//
//       // Rebuild field fid on *this* from source tree 'src' by sampling at parent nodes
//       void rebuild_field_from(const QuadTree2D& src, u32 fid) {
//         resize_fields_to_leaves();
//         Field& f = field(fid);
//
//         const auto& L = leaf_indices();
//         for (u32 k = 0; k < L.size(); ++k) {
//           u32 leaf = L[k];
//           std::vector<std::array<double, 2>> xi;
//           leaf_reference_nodes(f.basis, leaf, xi);
//
//           double* coeffs = leaf_coeff_ptr(fid, k);
//           for (size_t j = 0; j < xi.size(); ++j) {
//             double val;
//             if (!src.evaluate_field_on_parent(fid, xi[j][0], xi[j][1], val)) {
//               val = 0.0;
//               std::cout << "error!";
//             }
//             coeffs[j] = val;
//           }
//         }
//       }
//
//
//       // Conservative coarsen cycle using snapshot + parent coords; rebuild all fields
//       std::size_t coarsen_only_cycle_safe(u32 fid,
//                                           double tau_coarse,
//                                           u32 max_passes = 10,
//                                           Basis probe_basis = Basis::Q2_Quad9) {
//
//
//         QuadTree2D snapshot = *this;
//
//
//
//         // Predicate: coarsening criterion
//         auto pred = [&](u32 /*parent*/, u32 level,
//                         const std::vector<std::array<double, 2>>& pts_xi,
//                         const std::vector<std::array<double, 2>>& /*pts_xy*/,
//         const std::vector<std::array<double, 9>>& /*Nvals*/) -> bool {
//           if (level <= min_depth()) return false;
//           if (pts_xi.empty()) return false;
//
//           double v0;
//           if (!snapshot.evaluate_field_on_parent(fid, pts_xi[0][0], pts_xi[0][1], v0))
//             return false;
//
//           double mn = v0, mx = v0;
//           for (size_t i = 1; i < pts_xi.size(); ++i) {
//             double val;
//             if (snapshot.evaluate_field_on_parent(fid, pts_xi[i][0], pts_xi[i][1], val)) {
//               mn = std::min(mn, val);
//               mx = std::max(mx, val);
//             }
//           }
//           return (mn > +tau_coarse) || (mx < -tau_coarse);
//         };
//
//         std::size_t total = 0;
//
//         // Perform coarsening passes
//         for (u32 pass = 0; pass < max_passes; ++pass) {
//           std::size_t c = coarsen_pass(pred, probe_basis);
//           if (c == 0) break;
//
//           total += c;
//
//           // 🔧 Enforce 1-irregularity after each coarsening pass
//           enforce_balance();
//         }
//
//         // Rebuild all fields from snapshot (conservative transfer)
//         for (u32 f = 0; f < _fields.size(); ++f) {
//           rebuild_field_from(snapshot, f);
//         }
//
//         //enforce_hanging_constraints();
//
//         return total;
//       }
//
//
//       u32 find_neighbor(u32 leaf, int edge) const {
//         double xi0, eta0, xi1, eta1;
//         leaf_bounds(leaf, xi0, eta0, xi1, eta1);
//
//         double eps = 1.0E-10;
//
//         // local element sizes in parent space
//         double dx = xi1 - xi0;
//         double dy = eta1 - eta0;
//
//         // relative perturbation factor (safe wrt machine precision)
//
//
//         // test point starts at element center
//         double xm = 0.5 * (xi0 + xi1);
//         double ym = 0.5 * (eta0 + eta1);
//
//         // shift across the requested edge
//         switch (edge) {
//         case 0:
//           ym = eta0 - eps * dy;
//           break; // bottom
//         case 1:
//           xm = xi1 + eps * dx;
//           break; // right
//         case 2:
//           ym = eta1 + eps * dy;
//           break; // top
//         case 3:
//           xm = xi0 - eps * dx;
//           break; // left
//         default:
//           break;
//         }
//
//         return locate_leaf_on_parent(xm, ym);
//
//
//
//       }
//
//
//       // Collect edge nodes of a leaf and return both coordinates (xi,eta) and local indices.
//       // The coordinates are shifted slightly outside the fine element so evaluation
//       // comes from the coarse neighbor.
//       void edge_reference_nodes_with_mapping(Basis basis, u32 leaf, int edge,
//                                              std::vector<std::array<double, 2>>& out_coords,
//                                              std::vector<int>& out_indices) const {
//         out_coords.clear();
//         out_indices.clear();
//
//         std::vector<std::array<double, 2>> xi;
//         leaf_reference_nodes(basis, leaf, xi);
//
//         double eps = 1e-10;  // small outward shift
//         double xi0, eta0, xi1, eta1;
//         leaf_bounds(leaf, xi0, eta0, xi1, eta1);
//         const double dx = xi1 - xi0, dy = eta1 - eta0;
//
//         for (int i = 0; i < (int)xi.size(); ++i) {
//           auto p = xi[i];
//           bool isOnEdge = false;
//
//           if (edge == 0) { // bottom edge
//             p[1] -= eps * dy;
//             switch (i) {
//             case 0: // bottom-left corner
//               p[0] += eps * dx;
//               isOnEdge = true;
//               break;
//             case 1: // bottom-right corner
//               p[0] -= eps * dx;
//               isOnEdge = true;
//               break;
//             case 4: // bottom midpoint
//               isOnEdge = true;
//               break;
//             default:
//               break;
//             }
//           }
//           else if (edge == 1) { // right edge
//             p[0] += eps * dx;
//             switch (i) {
//             case 1: // bottom-right corner
//               p[1] += eps * dy;
//               isOnEdge = true;
//               break;
//             case 2: // top-right corner
//               p[1] -= eps * dy;
//               isOnEdge = true;
//               break;
//             case 5: // right midpoint
//               isOnEdge = true;
//               break;
//             default:
//               break;
//             }
//           }
//           else if (edge == 2) { // top edge
//             p[1] += eps * dy;
//             switch (i) {
//             case 2: // top-right corner
//               p[0] -= eps * dx;
//               isOnEdge = true;
//               break;
//             case 3: // top-left corner
//               p[0] += eps * dx;
//               isOnEdge = true;
//               break;
//             case 6: // top midpoint
//               isOnEdge = true;
//               break;
//             default:
//               break;
//             }
//           }
//           else if (edge == 3) { // left edge
//             p[0] -= eps * dx;
//             switch (i) {
//             case 0: // bottom-left corner
//               p[1] += eps * dy;
//               isOnEdge = true;
//               break;
//             case 3: // top-left corner
//               p[1] -= eps * dy;
//               isOnEdge = true;
//               break;
//             case 7: // left midpoint
//               isOnEdge = true;
//               break;
//             default:
//               break;
//             }
//           }
//
//           if (isOnEdge) {
//             out_coords.push_back(p);   // store nudged point
//             out_indices.push_back(i); // store local index of node
//           }
//         }
//       }
//
//       // Enforce hanging-node constraints for all fields
//       void enforce_hanging_constraints() {
//         for (u32 leaf : leaf_indices()) {
//           for (int e = 0; e < 4; ++e) {
//             u32 neigh = find_neighbor(leaf, e);
//             if (neigh == npos32) continue;
//
//             if (level_of(leaf) > level_of(neigh)) {
//               enforce_edge_constraints(neigh, leaf, e);
//             }
//           }
//         }
//       }
//
//
//       void enforce_edge_constraints(u32 coarse, u32 fine, int edge) {
//         for (u32 fid = 0; fid < _fields.size(); ++fid) {
//           const Field& f = _fields[fid];
//           double* coeffF = leaf_coeff_ptr(fid, leaf_position(fine));
//
//           std::vector<std::array<double, 2>> xiFine;
//           std::vector<int> idxFine;
//           edge_reference_nodes_with_mapping(f.basis, fine, edge, xiFine, idxFine);
//
//           for (size_t j = 0; j < xiFine.size(); ++j) {
//             double val;
//             if (!evaluate_field_on_parent(fid, xiFine[j][0], xiFine[j][1], val)) {
//               val = 0.0;
//               std::cout << "error!";
//             }
//             coeffF[idxFine[j]] = val;  // write into correct slot
//           }
//         }
//       }
//
//
//       u32 level_of(u32 leaf) const {
//         return _nodes[leaf].level;
//       }
//
//
//       // Map a parent coordinate (xi,eta) in [-1,1]^2 to physical (x,y)
//       // using the Quad9 isoparametric geometry map.
//       std::array<double, 2> parent_to_physical(double xi, double eta) const {
//         require_geometry();
//
//         double N9[9];
//         Quad9Shape::N(xi, eta, N9);
//
//         double X = 0.0, Y = 0.0;
//         for (int a = 0; a < 9; ++a) {
//           X += N9[a] * _X[a];
//           Y += N9[a] * _Y[a];
//         }
//         return {X, Y};
//       }
//
//
//
//       // Collect reference coordinates (xi,eta) of nodes for all leaves
//       // in a level range. Much like extract_node_coords_in_level_range,
//       // but stays in parent reference space (avoids inverse_map).
//       std::vector<std::array<double, 2>>
//       extract_node_parent_coords_in_level_range(u32 lev_min, u32 lev_max, Basis basis) const {
//         std::vector<std::array<double, 2>> coords;
//
//         for (u32 leaf : leaf_indices()) {
//           u32 lev = level_of(leaf);          // use accessor, not _nodes
//           if (lev < lev_min || lev > lev_max) continue;
//
//           std::vector<std::array<double, 2>> xi;
//           leaf_reference_nodes(basis, leaf, xi);
//
//           coords.insert(coords.end(), xi.begin(), xi.end());
//         }
//
//         return coords;
//       }
//
//
//
//       // ---------------------------------------------------------
//       // Neighbor lookup (axis-aligned): dir = 0:left, 1:right, 2:down, 3:up
//       // Returns the leaf covering the neighbor cell, or npos32 if outside.
//       // ---------------------------------------------------------
//       u32 neighbor_leaf(u32 leaf, int dir) const {
//         const QuadNode& n = _nodes[leaf];
//
//         // Decode Morton index -> (ix, iy)
//         u32 ix, iy;
//         deinterleave2(n.morton, ix, iy);
//
//         // Step one cell at this level
//         const u32 N = 1u << n.level;
//         if (dir == 0) {
//           if (ix == 0) return npos32;  // left
//           else ix -= 1;
//         }
//         if (dir == 1) {
//           if (ix + 1 >= N) return npos32;  // right
//           else ix += 1;
//         }
//         if (dir == 2) {
//           if (iy == 0) return npos32;  // down
//           else iy -= 1;
//         }
//         if (dir == 3) {
//           if (iy + 1 >= N) return npos32;  // up
//           else iy += 1;
//         }
//
//         // Map to parent-space coordinate at neighbor cell center
//         const double dx = 2.0 / double(N);
//         const double dy = 2.0 / double(N);
//         const double xi  = -1.0 + (ix + 0.5) * dx;
//         const double eta = -1.0 + (iy + 0.5) * dy;
//
//         // Locate leaf covering that center point
//         return locate_leaf_on_parent(xi, eta);
//       }
//
//       // ---------------------------------------------------------
//       // Enforce 1-irregularity: no adjacent leaves differ by >1 level
//       // ---------------------------------------------------------
//       void enforce_balance() {
//         bool changed = true;
//         while (changed) {
//           changed = false;
//           const auto& L = leaves();
//
//           for (u32 leaf : L) {
//             u32 lev = _nodes[leaf].level;
//
//             for (int dir = 0; dir < 4; ++dir) {
//               u32 nb = neighbor_leaf(leaf, dir);
//               if (nb == npos32) continue;
//
//               u32 lev_nb = _nodes[nb].level;
//               if (lev > lev_nb + 1) {
//                 if (refine(nb)) changed = true;
//               }
//             }
//           }
//           (void)leaves(); // refresh compact list
//         }
//       }
//
//       // Write current mesh + field to VTK UnstructuredGrid ASCII (.vtu)
//       bool write_vtu(const std::string & filename, u32 fid, const std::string & name,
//                      bool cell_centered = false) const {
//         require_geometry();
//         const auto& L = leaves();
//         const size_t numCells = L.size();
//         if (numCells == 0) return false;
//
//         const Field& fld = _fields[fid];
//         const size_t need = numCells * (size_t)fld.dofs_per_cell;
//         if (fld.coeffs.size() < need) {
//           assert(false && "Field coefficients not sized to current leaves. Call resize_fields_to_leaves().");
//           return false;
//         }
//
//         struct Pt {
//           double x, y, z;
//         };
//         std::vector<Pt> points;
//         std::vector<int> connectivity;
//         std::vector<int> offsets;
//         std::vector<unsigned char> types;
//         std::vector<double> pointData;
//         std::vector<double> cellData;
//         std::vector<int> levels;
//
//         const int ndof = (int)fld.dofs_per_cell;
//         const double* coeff_base = fld.coeffs.data();
//         int pointCounter = 0;
//
//         for (size_t k = 0; k < numCells; ++k) {
//           const u32 leaf = L[k];
//           levels.push_back((int)_nodes[leaf].level);
//           const double* c = coeff_base + size_t(k) * ndof;
//
//           // Get geometry nodes depending on basis
//           std::vector<std::array<double, 2>> xy;
//           leaf_physical_nodes(fld.basis, leaf, xy);
//
//           // Append nodes to point list
//           for (auto &p : xy) {
//             points.push_back({p[0], p[1], 0.0});
//           }
//
//           // Connectivity: sequential nodes
//           for (int i = 0; i < ndof; ++i)
//             connectivity.push_back(pointCounter + i);
//           pointCounter += ndof;
//           offsets.push_back((int)connectivity.size());
//
//           // Cell type
//           switch (fld.basis) {
//           case Basis::Q1_Quad4:
//             types.push_back(9);
//             break;  // VTK_QUAD
//           case Basis::Serendipity8:
//             types.push_back(23);
//             break;  // VTK_QUADRATIC_QUAD
//           case Basis::Q2_Quad9:
//             types.push_back(28);
//             break;  // VTK_BIQUADRATIC_QUAD
//           }
//
//           if (cell_centered) {
//             // Just evaluate at element center (0,0)
//             double v = 0.0;
//             switch (fld.basis) {
//             case Basis::Q1_Quad4: {
//               double N4[4];
//               Shapes::Q1(0, 0, N4);
//               for (int i = 0; i < 4; ++i) v += N4[i] * c[i];
//             }
//             break;
//             case Basis::Serendipity8: {
//               double N8[8];
//               Shapes::Serendipity8(0, 0, N8);
//               for (int i = 0; i < 8; ++i) v += N8[i] * c[i];
//             }
//             break;
//             case Basis::Q2_Quad9: {
//               double N9[9];
//               Quad9Shape::N(0, 0, N9);
//               for (int i = 0; i < 9; ++i) v += N9[i] * c[i];
//             }
//             break;
//             }
//             cellData.push_back(v);
//           }
//           else {
//             // Push values at nodes
//             for (int i = 0; i < ndof; ++i) pointData.push_back(c[i]);
//           }
//         }
//
//         // Write ASCII VTK
//         std::ofstream os(filename);
//         if (!os) return false;
//         os << std::setprecision(16);
//         os << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
//         os << "  <UnstructuredGrid>\n";
//         os << "    <Piece NumberOfPoints=\"" << points.size()
//            << "\" NumberOfCells=\"" << numCells << "\">\n";
//
//         os << "      <Points>\n";
//         os << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
//         for (auto &p : points) os << "          " << p.x << " " << p.y << " " << p.z << "\n";
//         os << "        </DataArray>\n";
//         os << "      </Points>\n";
//
//         os << "      <Cells>\n";
//         os << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
//         for (size_t i = 0; i < connectivity.size(); ++i) {
//           os << connectivity[i] << ((i + 1) % ndof == 0 ? "\n" : " ");
//         }
//         os << "        </DataArray>\n";
//         os << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
//         for (int off : offsets) os << off << "\n";
//         os << "        </DataArray>\n";
//         os << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
//         for (auto t : types) os << (int)t << "\n";
//         os << "        </DataArray>\n";
//         os << "      </Cells>\n";
//
//         os << "      <CellData Scalars=\"" << name << "\">\n";
//         os << "        <DataArray type=\"Int32\" Name=\"level\" format=\"ascii\">\n";
//         for (int lv : levels) os << lv << "\n";
//         os << "        </DataArray>\n";
//
//         if (cell_centered) {
//           os << "        <DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
//           for (double v : cellData) os << v << "\n";
//           os << "        </DataArray>\n";
//           os << "      </CellData>\n";
//         }
//         else {
//           os << "      </CellData>\n";
//           os << "      <PointData Scalars=\"" << name << "\">\n";
//           os << "        <DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
//           for (double v : pointData) os << v << "\n";
//           os << "        </DataArray>\n";
//           os << "      </PointData>\n";
//         }
//
//         os << "    </Piece>\n";
//         os << "  </UnstructuredGrid>\n";
//         os << "</VTKFile>\n";
//         return true;
//       }
//
//
//       // -----------------------------------------------------------------------------
// // Main VTU writer (binary)
// // -----------------------------------------------------------------------------
//
//
//       bool write_binary_vtu(const std::string &filename, u32 fid, const std::string &name,
//                             bool cell_centered = false) const {
//         require_geometry();
//         const auto& L = leaves();
//         const size_t numCells = L.size();
//         if (numCells == 0) return false;
//
//         const Field& fld = _fields[fid];
//         const size_t need = numCells * (size_t)fld.dofs_per_cell;
//         if (fld.coeffs.size() < need) {
//           assert(false && "Field coefficients not sized to current leaves. Call resize_fields_to_leaves().");
//           return false;
//         }
//
//         struct Pt {
//           double x, y, z;
//         };
//         std::vector<Pt> points;
//         std::vector<int> connectivity;
//         std::vector<int> offsets;
//         std::vector<unsigned char> types;
//         std::vector<double> pointData;
//         std::vector<double> cellData;
//         std::vector<int> levels;
//
//         const int ndof = (int)fld.dofs_per_cell;   // 9 for Q2
//         const double* coeff_base = fld.coeffs.data();
//
//         // ---- Reserve memory up front ----
//         points.reserve(numCells * ndof);
//         connectivity.reserve(numCells * ndof);
//         offsets.reserve(numCells);
//         types.reserve(numCells);
//         levels.reserve(numCells);
//         if (cell_centered) {
//           cellData.reserve(numCells);
//         }
//         else {
//           pointData.reserve(numCells * ndof);
//         }
//
//         int pointCounter = 0;
//
//         for (size_t k = 0; k < numCells; ++k) {
//           const u32 leaf = L[k];
//           levels.push_back((int)_nodes[leaf].level);
//           const double* c = coeff_base + size_t(k) * ndof;
//
//           std::vector<std::array<double, 2>> xy;
//           leaf_physical_nodes(Basis::Q2_Quad9, leaf, xy);
//
//           for (auto &p : xy) {
//             points.push_back({p[0], p[1], 0.0});
//           }
//
//           for (int i = 0; i < ndof; ++i)
//             connectivity.push_back(pointCounter + i);
//           pointCounter += ndof;
//           offsets.push_back((int)connectivity.size());
//           types.push_back(28); // VTK_BIQUADRATIC_QUAD
//
//           if (cell_centered) {
//             // Evaluate at element center (0,0)
//             double v = 0.0;
//             double N9[9];
//             Quad9Shape::N(0, 0, N9);
//             for (int i = 0; i < 9; ++i) v += N9[i] * c[i];
//             cellData.push_back(v);
//           }
//           else {
//             for (int i = 0; i < ndof; ++i)
//               pointData.push_back(c[i]);
//           }
//         }
//
//         // Flatten points
//         std::vector<double> flatPoints;
//         flatPoints.reserve(points.size() * 3);
//         for (auto &p : points) {
//           flatPoints.push_back(p.x);
//           flatPoints.push_back(p.y);
//           flatPoints.push_back(p.z);
//         }
//
//         // ---- Write file ----
//         std::ofstream os(filename);
//         if (!os) return false;
//
//         os << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
//         os << "  <UnstructuredGrid>\n";
//         os << "    <Piece NumberOfPoints=\"" << points.size()
//            << "\" NumberOfCells=\"" << numCells << "\">\n";
//
//         // Points
//         os << "      <Points>\n";
//         write_binary_array(os, "Float64", "", 3, flatPoints);
//         os << "      </Points>\n";
//
//         // Cells
//         os << "      <Cells>\n";
//         write_binary_array(os, "Int32", "connectivity", 1, connectivity);
//         write_binary_array(os, "Int32", "offsets", 1, offsets);
//         write_binary_array(os, "UInt8", "types", 1, types);
//         os << "      </Cells>\n";
//
//         // CellData: refinement level always written
//         os << "      <CellData Scalars=\"" << name << "\">\n";
//         write_binary_array(os, "Int32", "level", 1, levels);
//
//         if (cell_centered) {
//           write_binary_array(os, "Float64", name, 1, cellData);
//           os << "      </CellData>\n";
//         }
//         else {
//           os << "      </CellData>\n";
//           os << "      <PointData Scalars=\"" << name << "\">\n";
//           write_binary_array(os, "Float64", name, 1, pointData);
//           os << "      </PointData>\n";
//         }
//
//         os << "    </Piece>\n";
//         os << "  </UnstructuredGrid>\n";
//         os << "</VTKFile>\n";
//
//         return true;
//       }
//
//
//
//
//
//
//       bool write_binary_vtu_mesh(const std::string &filename) const {
//         require_geometry();
//         const auto& L = leaves();
//         if (L.empty()) return false;
//
//         struct Pt {
//           double x, y, z;
//         };
//         std::vector<Pt> points;
//         std::vector<int> connectivity;
//         std::vector<int> offsets;
//         std::vector<unsigned char> types;
//
//         const int vtk_type = 28;  // VTK_BIQUADRATIC_QUAD
//         const int ndof = 9;       // Q2 element has 9 nodes
//
//         int pointCounter = 0;
//
//         for (size_t k = 0; k < L.size(); ++k) {
//           const u32 leaf = L[k];
//           std::vector<std::array<double, 2>> xy;
//           leaf_physical_nodes(Basis::Q2_Quad9, leaf, xy);
//
//           // add nodes
//           for (auto &p : xy) {
//             points.push_back({p[0], p[1], 0.0});
//           }
//
//           // connectivity
//           for (int i = 0; i < ndof; ++i) {
//             connectivity.push_back(pointCounter + i);
//           }
//           pointCounter += ndof;
//
//           offsets.push_back((int)connectivity.size());
//           types.push_back((unsigned char)vtk_type);
//         }
//
//         // Flatten points into xyz array
//         std::vector<double> flatPoints;
//         flatPoints.reserve(points.size() * 3);
//         for (auto &p : points) {
//           flatPoints.push_back(p.x);
//           flatPoints.push_back(p.y);
//           flatPoints.push_back(p.z);
//         }
//
//         // ---- Write file ----
//         std::ofstream os(filename);
//         if (!os) return false;
//
//         os << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
//         os << "  <UnstructuredGrid>\n";
//         os << "    <Piece NumberOfPoints=\"" << points.size()
//            << "\" NumberOfCells=\"" << L.size() << "\">\n";
//
//         // Points
//         os << "      <Points>\n";
//         write_binary_array(os, "Float64", "", 3, flatPoints);
//         os << "      </Points>\n";
//
//         // Cells
//         os << "      <Cells>\n";
//         write_binary_array(os, "Int32", "connectivity", 1, connectivity);
//         write_binary_array(os, "Int32", "offsets", 1, offsets);
//         write_binary_array(os, "UInt8", "types", 1, types);
//         os << "      </Cells>\n";
//
//         os << "    </Piece>\n";
//         os << "  </UnstructuredGrid>\n";
//         os << "</VTKFile>\n";
//
//         return true;
//       }
//
//
//
//
//       bool write_binary_vtu_points(const std::string &filename) const {
//         require_geometry();
//         const auto& L = leaves();
//         if (L.empty()) return false;
//
//         struct Pt {
//           double x, y, z;
//         };
//         std::vector<Pt> points;
//
//         // Collect Q1 quad corner nodes just for testing
//         for (size_t k = 0; k < L.size(); ++k) {
//           const u32 leaf = L[k];
//           std::vector<std::array<double, 2>> xy;
//           leaf_physical_nodes(Basis::Q1_Quad4, leaf, xy);
//           for (auto &p : xy)
//             points.push_back({p[0], p[1], 0.0});
//         }
//
//         // Flatten points
//         std::vector<double> flatPoints;
//         flatPoints.reserve(points.size() * 3);
//         for (auto &p : points) {
//           flatPoints.push_back(p.x);
//           flatPoints.push_back(p.y);
//           flatPoints.push_back(p.z);
//         }
//
//         // Connectivity: one vertex per point
//         std::vector<int> connectivity;
//         std::vector<int> offsets;
//         std::vector<unsigned char> types;
//         for (size_t i = 0; i < points.size(); ++i) {
//           connectivity.push_back((int)i);
//           offsets.push_back((int)(i + 1));
//           types.push_back(1); // VTK_VERTEX
//         }
//
//         // ---- Write file ----
//         std::ofstream os(filename);
//         if (!os) return false;
//
//         os << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
//         os << "  <UnstructuredGrid>\n";
//         os << "    <Piece NumberOfPoints=\"" << points.size()
//            << "\" NumberOfCells=\"" << points.size() << "\">\n";
//
//         // Points
//         os << "      <Points>\n";
//         write_binary_array(os, "Float64", "", 3, flatPoints);
//         os << "      </Points>\n";
//
//         // Cells
//         os << "      <Cells>\n";
//         write_binary_array(os, "Int32", "connectivity", 1, connectivity);
//         write_binary_array(os, "Int32", "offsets", 1, offsets);
//         write_binary_array(os, "UInt8", "types", 1, types);
//         os << "      </Cells>\n";
//
//         os << "    </Piece>\n";
//         os << "  </UnstructuredGrid>\n";
//         os << "</VTKFile>\n";
//
//         return true;
//       }
//
//
//
//
//
//
//
//       // === Data extraction utilities ===
//
//       // Test if two points are the same within tolerance (used for global gather)
//       static inline bool same_point(const std::array<double, 2>& a,
//                                     const std::array<double, 2>& b,
//                                     double tol = 1e-12) {
//         return (std::fabs(a[0] - b[0]) < tol && std::fabs(a[1] - b[1]) < tol);
//       }
//
//       // Return number of nodes for a given basis
//       static int basis_nodes(Basis b) {
//         switch (b) {
//         case Basis::Q1_Quad4:
//           return 4;
//         case Basis::Serendipity8:
//           return 8;
//         case Basis::Q2_Quad9:
//           return 9;
//         }
//         return 0;
//       }
//
//       // Gather unique coords/values across leaves into global arrays
//       void gather_field(u32 fid) {
//         require_geometry();
//         struct Flat {
//           std::array<double, 2> xy;
//           double val;
//           int elem, local;
//         };
//         std::vector<Flat> flat;
//
//         const auto& L = leaves();
//         const Field& f = _fields[fid];
//         int nNodes = f.dofs_per_cell;
//
//         for (int e = 0; e < (int)L.size(); ++e) {
//           u32 leaf = L[e];
//           std::vector<std::array<double, 2>> pts_xy;
//           leaf_physical_nodes(f.basis, leaf, pts_xy);
//           const double* coeffs = leaf_coeff_ptr(fid, e);
//           for (int i = 0; i < nNodes; ++i) flat.push_back({pts_xy[i], coeffs[i], e, i});
//         }
//
//         std::sort(flat.begin(), flat.end(),
//         [](auto & a, auto & b) {
//           return (a.xy[0] < b.xy[0]) ||
//                  (a.xy[0] == b.xy[0] && a.xy[1] < b.xy[1]);
//         });
//
//         global_coords.clear();
//         global_field.clear();
//         elem2glob.assign(L.size(), {});
//
//         int gid = -1;
//         for (size_t k = 0; k < flat.size(); ++k) {
//           if (k == 0 || !same_point(flat[k].xy, flat[k - 1].xy)) {
//             ++gid;
//             global_coords.push_back(flat[k].xy);
//             global_field.push_back(flat[k].val);
//           }
//           if (elem2glob[flat[k].elem].empty()) elem2glob[flat[k].elem].resize(nNodes);
//           elem2glob[flat[k].elem][flat[k].local] = gid;
//         }
//       }
//
//       // Scatter modified global_field back to leaf coefficient storage
//       void scatter_field(u32 fid) {
//         auto& f = _fields[fid];
//         int nNodes = f.dofs_per_cell;
//         const auto& L = leaves();
//         for (int e = 0; e < (int)L.size(); ++e) {
//           double* coeffs = leaf_coeff_ptr(fid, e);
//           for (int i = 0; i < nNodes; ++i) {
//             int gid = elem2glob[e][i];
//             coeffs[i] = global_field[gid];
//           }
//         }
//       }
//
//       // Refine tree so that each point lies in a leaf of at least 'maxDepthTarget'
//       void refine_to_contain_points(const std::vector<std::array<double, 2>>& pts,
//                                     u32 maxDepthTarget) {
//         for (const auto& p : pts) {
//           double xi, eta;
//           if (!inverse_map_quad9(p[0], p[1], xi, eta)) continue;
//
//           u32 leaf = locate_leaf_on_parent(xi, eta);
//           while (leaf != npos32 && _nodes[leaf].level < maxDepthTarget) {
//             if (!refine(leaf)) break;
//             leaf = locate_leaf_on_parent(xi, eta);
//           }
//         }
//         (void)leaves();
//         enforce_balance();
//       }
//
//       // Extract node coordinates (physical) for all leaves in [minLevel,maxLevel]
//       std::vector<std::vector<std::array<double, 2>>>
//       extract_leaf_nodes_in_level_range(u32 minLevel, u32 maxLevel, Basis basis) const {
//         require_geometry();
//         std::vector<std::vector<std::array<double, 2>>> result;
//         const auto& L = leaves();
//         for (u32 leaf : L) {
//           u32 lev = _nodes[leaf].level;
//           if (lev >= minLevel && lev <= maxLevel) {
//             std::vector<std::array<double, 2>> xy;
//             leaf_physical_nodes(basis, leaf, xy);
//             result.push_back(std::move(xy));
//           }
//         }
//         return result;
//       }
//
//       // Extract all node coordinates (physical) for leaves in [minLevel,maxLevel] (flat list)
//       std::vector<std::array<double, 2>>
//       extract_node_coords_in_level_range(u32 minLevel, u32 maxLevel, Basis basis) const {
//         require_geometry();
//         std::vector<std::array<double, 2>> coords;
//         const auto& L = leaves();
//         for (u32 leaf : L) {
//           u32 lev = _nodes[leaf].level;
//           if (lev >= minLevel && lev <= maxLevel) {
//             std::vector<std::array<double, 2>> xy;
//             leaf_physical_nodes(basis, leaf, xy);
//             coords.insert(coords.end(), xy.begin(), xy.end());
//           }
//         }
//         return coords;
//       }
//
//       // Reset to single root cell; keep existing fields but size to 1 leaf
//       void reset() {
//         _nodes.clear();
//         _alive.clear();
//         _free.clear();
//         _leaves.clear();
//         _leaf_dirty = true;
//
//         _root = alloc_node();
//         _nodes[_root].level = 0;
//         _nodes[_root].morton = 0;
//
//         for (auto& f : _fields) {
//           f.coeffs.clear();
//           f.coeffs.resize(f.dofs_per_cell);
//         }
//       }
//
//       // Return global maximum depth allowed
//       u32 max_depth() const {
//         return _maxDepth;
//       }
//
//       // Accessors for coarse geometry coordinates (const-ref)
//       inline const std::array<double, 9>& get_X() const {
//         return _X;
//       }
//       inline const std::array<double, 9>& get_Y() const {
//         return _Y;
//       }
//
//       // === Public state used by gather/scatter utilities ===
//       std::vector<std::array<double, 2>> global_coords;     // gathered unique coords
//       std::vector<double>                global_field;      // gathered corresponding values
//       std::vector<std::vector<int>>      elem2glob;         // map elem-local -> global ids
//
//     private:
//       // Coarse geometry (Q2 nodes)
//       std::array<double, 9> _X{};
//       std::array<double, 9> _Y{};
//
//       // Ensure all leaves satisfy min depth floor (by refining only)
//       void ensure_min_depth() {
//         if (_minDepth == 0) return;
//         bool changed = true;
//         while (changed) {
//           changed = false;
//           const auto& L = leaves();
//           for (u32 leaf : L) {
//             if (_nodes[leaf].level < _minDepth) {
//               if (refine(leaf)) changed = true;
//             }
//           }
//         }
//       }
//
//       // Allocate a new node id (reuse from freelist if available)
//       u32 alloc_node() {
//         u32 i;
//         if (!_free.empty()) {
//           i = _free.back();
//           _free.pop_back();
//         }
//         else {
//           i = (u32)_nodes.size();
//           _nodes.emplace_back();
//           _alive.push_back(0);
//         }
//         _alive[i] = 1;
//         _nodes[i] = QuadNode();
//         return i;
//       }
//
//       // Free a node id (push to freelist), mark leaves dirty
//       void free_node(u32 i) {
//         _alive[i] = 0;
//         _free.push_back(i);
//         _leaf_dirty = true;
//       }
//
//       // Bounds helper for any node
//       inline void bounds_of(const QuadNode & n, double & xi0, double & eta0, double & xi1, double & eta1) const {
//         u32 ix, iy;
//         deinterleave2(n.morton, ix, iy);
//         const double N = double(1u << n.level);
//         const double dx = 2.0 / N, dy = 2.0 / N;
//         xi0  = -1.0 + ix * dx;
//         eta0 = -1.0 + iy * dy;
//         xi1  = xi0 + dx;
//         eta1 = eta0 + dy;
//       }
//
//       // Return compact-leaf position for node id, or npos32 if not a leaf
//       inline u32 leaf_position(u32 leaf) const {
//         (void)leaves();
//         return (leaf < _leaf_pos.size()) ? _leaf_pos[leaf] : npos32;
//       }
//
//       // --- Internal data ---
//       u32  _root{npos32};
//       u32  _maxDepth;
//       u32  _minDepth{0};
//       bool _allowCoarsenBelowMinDepth{true};
//       std::vector<QuadNode> _nodes;
//       std::vector<uint8_t>  _alive;
//       std::vector<u32>      _free;
//       mutable std::vector<u32> _leaves;
//       mutable bool _leaf_dirty{true};
//       std::vector<Field> _fields;
//       mutable std::vector<u32> _leaf_pos;
//       bool _geom_ready{false};
//   };
//
// } // namespace fem
//
//


