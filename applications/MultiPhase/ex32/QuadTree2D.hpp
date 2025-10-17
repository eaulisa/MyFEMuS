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

  using Point2 = std::array<double, 2>;

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

  // inline void deinterleave2(u64 m, u32 &x, u32 &y) {
  //   u64 xx = m, yy = m >> 1;
  //   xx &= 0x5555555555555555ULL;
  //   yy &= 0x5555555555555555ULL;
  //   xx = (xx | (xx >> 1)) & 0x3333333333333333ULL;
  //   xx = (xx | (xx >> 2)) & 0x0F0F0F0F0F0F0F0FULL;
  //   xx = (xx | (xx >> 4)) & 0x00FF00FF00FF00FFULL;
  //   xx = (xx | (xx >> 8)) & 0x0000FFFF0000FFFFULL;
  //   xx = (xx | (xx >> 16)) & 0x00000000FFFFFFFFULL;
  //   yy = (yy | (yy >> 1)) & 0x3333333333333333ULL;
  //   yy = (yy | (yy >> 2)) & 0x0F0F0F0F0F0F0F0FULL;
  //   yy = (yy | (yy >> 4)) & 0x00FF00FF00FF00FFULL;
  //   yy = (yy | (yy >> 8)) & 0x0000FFFF0000FFFFULL;
  //   yy = (yy | (yy >> 16)) & 0x00000000FFFFFFFFULL;
  //   x = (u32)xx;
  //   y = (u32)yy;
  // }

//--------------------------------------------
// Shapes / geometry
//--------------------------------------------
  struct Shapes {
    public:
      // ---------- Q1 impl ----------
      static inline void Q1(const Point2& se, double* __restrict__ N) noexcept {
        const double xi  = se[0];
        const double eta = se[1];
        const double a = 0.25 * (1.0 - xi);
        const double b = 0.25 * (1.0 + xi);
        const double tm = (1.0 - eta);
        const double tp = (1.0 + eta);
        N[0] = a * tm;  // (-1,-1)
        N[1] = b * tm;  // ( 1,-1)
        N[2] = b * tp;  // ( 1, 1)
        N[3] = a * tp;  // (-1, 1)
      }

      static inline void Q1_dN(const Point2& se,
                               double* __restrict__ dNdxi,
                               double* __restrict__ dNdeta) noexcept {
        const double xi  = se[0];
        const double eta = se[1];
        const double sm = 1.0 - xi, sp = 1.0 + xi;
        const double tm = 1.0 - eta, tp = 1.0 + eta;
        const double c  = 0.25;

        dNdxi[0] = -c * tm;  dNdeta[0] = -c * sm;
        dNdxi[1] =  c * tm;  dNdeta[1] = -c * sp;
        dNdxi[2] =  c * tp;  dNdeta[2] =  c * sp;
        dNdxi[3] = -c * tp;  dNdeta[3] =  c * sm;
      }

      // ---------- Q8 impl ----------
      static inline void Q8(const Point2& se, double* __restrict__ N) noexcept {
        const double xi  = se[0];
        const double eta = se[1];

        const double xm = xi, em = eta;
        const double one_m_x = 1.0 - xm;
        const double one_p_x = 1.0 + xm;
        const double one_m_e = 1.0 - em;
        const double one_p_e = 1.0 + em;
        const double x2 = xm * xm;
        const double e2 = em * em;

        N[0] = 0.25 * one_m_x * one_m_e * (-xm - em - 1.0);
        N[1] = 0.25 * one_p_x * one_m_e * (xm - em - 1.0);
        N[2] = 0.25 * one_p_x * one_p_e * (xm + em - 1.0);
        N[3] = 0.25 * one_m_x * one_p_e * (-xm + em - 1.0);
        N[4] = 0.5  * (1.0 - x2) * one_m_e;
        N[5] = 0.5  * one_p_x * (1.0 - e2);
        N[6] = 0.5  * (1.0 - x2) * one_p_e;
        N[7] = 0.5  * one_m_x * (1.0 - e2);
      }

      static inline void Q8_dN(const Point2& se,
                               double* __restrict__ dNdxi,
                               double* __restrict__ dNdeta) noexcept {
        const double x = se[0];
        const double e = se[1];

        const double one_m_x = 1.0 - x;
        const double one_p_x = 1.0 + x;
        const double one_m_e = 1.0 - e;
        const double one_p_e = 1.0 + e;

        // corners
        dNdxi[0]  = 0.25 * one_m_e * (2.0 * x + e);
        dNdeta[0] = 0.25 * one_m_x * (x + 2.0 * e);

        dNdxi[1]  = 0.25 * one_m_e * (2.0 * x - e);
        dNdeta[1] = 0.25 * one_p_x * (-x + 2.0 * e);

        dNdxi[2]  = 0.25 * one_p_e * (2.0 * x + e);
        dNdeta[2] = 0.25 * one_p_x * (x + 2.0 * e);

        dNdxi[3]  = 0.25 * one_p_e * (2.0 * x - e);
        dNdeta[3] = 0.25 * one_m_x * (-x + 2.0 * e);

        // edges
        dNdxi[4]  =    -x * one_m_e;            // d/dx of 0.5*(1-x^2)*(1-e)
        dNdeta[4] = -0.5 * (1.0 - x * x);

        dNdxi[5]  =  0.5 * (1.0 - e * e);
        dNdeta[5] =    -(1.0 + x) * e;

        dNdxi[6]  =    -x * one_p_e;
        dNdeta[6] =  0.5 * (1.0 - x * x);

        dNdxi[7]  = -0.5 * (1.0 - e * e);
        dNdeta[7] =    -(1.0 - x) * e;
      }

      static inline void Q9(const Point2& se, double* __restrict__ N) noexcept {
        const double xi  = se[0];
        const double eta = se[1];

        double Lx0, Lx1, Lx2;
        double Ly0, Ly1, Ly2;
        q2_1d_vals(xi,  Lx0, Lx1, Lx2);
        q2_1d_vals(eta, Ly0, Ly1, Ly2);

        // corners 0..3, edges 4..7 (bottom,right,top,left), center 8
        N[0] = Lx0 * Ly0;  // (-1,-1)
        N[1] = Lx2 * Ly0;  // ( 1,-1)
        N[2] = Lx2 * Ly2;  // ( 1, 1)
        N[3] = Lx0 * Ly2;  // (-1, 1)
        N[4] = Lx1 * Ly0;  // ( 0,-1)
        N[5] = Lx2 * Ly1;  // ( 1, 0)
        N[6] = Lx1 * Ly2;  // ( 0, 1)
        N[7] = Lx0 * Ly1;  // (-1, 0)
        N[8] = Lx1 * Ly1;  // ( 0, 0)
      }

      static inline void Q9_dN(const Point2& se,
                               double* __restrict__ dNdxi,
                               double* __restrict__ dNdeta) noexcept {
        const double xi  = se[0];
        const double eta = se[1];

        double Lx0, Lx1, Lx2, dx0, dx1, dx2;
        double Ly0, Ly1, Ly2, dy0, dy1, dy2;

        q2_1d_vals(xi,  Lx0, Lx1, Lx2);
        q2_1d_derivs(xi,  dx0, dx1, dx2);
        q2_1d_vals(eta, Ly0, Ly1, Ly2);
        q2_1d_derivs(eta, dy0, dy1, dy2);

        dNdxi[0] = dx0 * Ly0;  dNdeta[0] = Lx0 * dy0;
        dNdxi[1] = dx2 * Ly0;  dNdeta[1] = Lx2 * dy0;
        dNdxi[2] = dx2 * Ly2;  dNdeta[2] = Lx2 * dy2;
        dNdxi[3] = dx0 * Ly2;  dNdeta[3] = Lx0 * dy2;

        dNdxi[4] = dx1 * Ly0;  dNdeta[4] = Lx1 * dy0;
        dNdxi[5] = dx2 * Ly1;  dNdeta[5] = Lx2 * dy1;
        dNdxi[6] = dx1 * Ly2;  dNdeta[6] = Lx1 * dy2;
        dNdxi[7] = dx0 * Ly1;  dNdeta[7] = Lx0 * dy1;

        dNdxi[8] = dx1 * Ly1;  dNdeta[8] = Lx1 * dy1;
      }

    private:
      // ---------- Q9 impl ----------
      static inline void q2_1d_vals(double s, double& L0, double& L1, double& L2) noexcept {
        const double s2 = s * s;
        L0 = 0.5 * (s2 - s);  // 0.5*s*(s-1)
        L1 = 1.0 - s2;
        L2 = 0.5 * (s2 + s);  // 0.5*s*(s+1)
      }
      static inline void q2_1d_derivs(double s, double& dL0, double& dL1, double& dL2) noexcept {
        dL0 = s - 0.5;
        dL1 = -2.0 * s;
        dL2 = s + 0.5;
      }
  };


//--------------------------------------------
// Core data
//--------------------------------------------

// Quadtree node (internal or leaf)
  struct TreeNode {
    u64 morton{0};            // Morton code of this node
    u32 ix{0}, iy{0};
    u32 level{0};             // Refinement level
    u32 child[4] {npos32, npos32, npos32, npos32}; // Children indices (npos32 = none)
    bool is_leaf{true};       // True if this node is a leaf
    u32 parent{npos32};       // Parent node index (npos32 if root)
  };

  struct Field {
    Basis basis{Basis::Q1_Quad4};
    // Nodal values, one per unique global node in the basis registry.
    std::vector<double> nodal;  // size == _basisReg[(int)basis].nodes.size()

  };

  struct FEMNode {
    int gid;
    Point2 parent;
    Point2 physical;
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
        root.ix = 0;
        root.iy = 0;
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

      void reset(bool keep_geometry = true, bool keep_fields = true) {
        // 1) optionally keep geometry
        if (!keep_geometry) {
          _X = {}; _Y = {};
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

        TreeNode root;
        root.ix = 0;
        root.iy = 0;
        root.morton  = interleave2(0, 0);
        root.level   = 0;
        root.is_leaf = true;
        root.parent  = npos32;
        root.child[0] = root.child[1] = root.child[2] = root.child[3] = npos32;

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

      QuadTree2D(QuadTree2D&&) noexcept = default;
      QuadTree2D& operator=(QuadTree2D&&) noexcept = default;
      QuadTree2D(const QuadTree2D&) = delete;
      QuadTree2D& operator=(const QuadTree2D&) = delete;

      friend void swap(QuadTree2D& a, QuadTree2D& b) noexcept {
        using std::swap;
        swap(a._maxDepth, b._maxDepth);
        swap(a._minDepth, b._minDepth);
        swap(a._allowCoarsenBelowMinDepth, b._allowCoarsenBelowMinDepth);
        swap(a._leafpos_valid, b._leafpos_valid);
        swap(a._X, b._X);
        swap(a._Y, b._Y);
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
      void extract_leaf_parent_coords(Basis basis, u32 leaf_idx,
                                      std::vector<Point2>& out_pts) const {

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
                               std::vector<Point2>& out_xy) const {
        require_geometry();
        std::vector<Point2> xi;
        extract_leaf_parent_coords(basis, leaf_idx, xi);
        out_xy.resize(xi.size());
        for (size_t i = 0; i < xi.size(); ++i) out_xy[i] = parent_to_physical(xi[i]);
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
          std::vector<Point2> pts_xi, pts_xy;
          extract_leaf_parent_coords(probe_basis, idx, pts_xi);
          leaf_physical_nodes(probe_basis, idx, pts_xy);
          std::vector<std::array<double, 9>> dummy;

          if (should_refine(idx, pts_xi, pts_xy, dummy) && leaf_copy.level < _maxDepth) {
            // split into 4 children

            u32 ix = leaf_copy.ix;
            u32 iy = leaf_copy.iy;
            //u32 ix, iy; deinterleave2(leaf_copy.morton, ix, iy);

            for (int dy = 0; dy < 2; ++dy) {
              for (int dx = 0; dx < 2; ++dx) {
                u64 cmorton = interleave2(2 * ix + dx, 2 * iy + dy);

                TreeNode child;
                child.ix = 2 * ix + dx;
                child.iy = 2 * iy + dy;
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

        return refined;
      }


      template<class CoarsenPred>
      std::size_t coarsen_pass(CoarsenPred&& should_coarsen, Basis probe_basis) {
        struct Group {
          u32      pos[4];
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
        ++per_epoch;                               // bump generation (wrap is astronomically far)

        touched.clear(); touched.reserve(L / 4 + 8);
        to_add.clear();  to_add.reserve(L / 4 + 8);
        newLeaves.clear(); newLeaves.reserve(L);

        // removal marks (epoch-tagged, no clears)
        static thread_local std::vector<uint32_t> mark_tag;
        static thread_local uint32_t mark_epoch = 1;
        if (mark_tag.size() < L) mark_tag.resize(L, 0);
        ++mark_epoch;

        auto mark_pos = [&](u32 i)       {
          mark_tag[i] = mark_epoch;
        };

        auto is_marked = [&](u32 i) -> bool {
          return (i < mark_tag.size()) && (mark_tag[i] == mark_epoch);
        };


        // tiny fixed buffers for parent-space samples
        std::array<Point2, 9> xi_buf, xy_buf;

        // Reuse vector wrappers for the predicate (avoid allocs)
        static thread_local std::vector<Point2> pts_xi_vec, pts_xy_vec;

        // how many probe points?
        const int nprobe = (probe_basis == Basis::Q1_Quad4 ? 4 :
                            (probe_basis == Basis::Serendipity8 ? 8 : 9));

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
          if (g.cnt < 4) g.pos[g.cnt++] = pos;
        }

        // 2) coarsen candidates
        std::size_t coarsened = 0;

        for (u32 pidx : touched) {
          Group& g = per[pidx];
          if (g.cnt != 4) continue;

          // level/minDepth check
          const u32 kid0 = _leaves[g.pos[0]];
          const u32 lev  = _tree_nodes[kid0].level;
          if (lev <= _minDepth && !_allowCoarsenBelowMinDepth) continue;

          // probe parent cell
          TreeNode& P = _tree_nodes[pidx];
          double xi0, eta0, xi1, eta1; leaf_bounds(P, xi0, eta0, xi1, eta1);
          const double xm = 0.5 * (xi0 + xi1), ym = 0.5 * (eta0 + eta1);

          // fill parent-space points into stack buffer
          int n = 0;
          xi_buf[n++] = {xi0, eta0}; xi_buf[n++] = {xi1, eta0};
          xi_buf[n++] = {xi1, eta1}; xi_buf[n++] = {xi0, eta1};
          if (nprobe >= 8) {
            xi_buf[n++] = {xm,  eta0}; xi_buf[n++] = {xi1, ym};
            xi_buf[n++] = {xm,  eta1}; xi_buf[n++] = {xi0, ym};
          }
          if (nprobe == 9) {
            xi_buf[n++] = {xm,  ym};
          }

          // map to physical (cheap Q9 mat-vec)
          for (int i = 0; i < n; ++i) {
            xy_buf[i] = parent_to_physical(xi_buf[i]);
          }

          // wrap with pre-sized vectors (no alloc)
          pts_xi_vec.resize(n);
          pts_xy_vec.resize(n);
          // memcpy would also be fine — this is tiny
          for (int i = 0; i < n; ++i) {
            pts_xi_vec[i] = xi_buf[i];
            pts_xy_vec[i] = xy_buf[i];
          }

          if (!should_coarsen(P.morton, lev, pts_xi_vec, pts_xy_vec, /*dummy*/{}))
            continue;

          // mark 4 kids by leaf position; don’t mutate _leaves yet
          for (int k = 0; k < 4; ++k) {
            const u32 pos = g.pos[k];
            _tree_nodes[_leaves[pos]].is_leaf = false;
            mark_pos(pos);
          }

          // parent becomes leaf
          P.is_leaf = true;
          P.child[0] = P.child[1] = P.child[2] = P.child[3] = npos32;
          to_add.push_back(pidx);
          ++coarsened;
        }

        // 3) rebuild leaf set (single pass)
        if (coarsened) {
          for (u32 i = 0; i < L; ++i) if (!is_marked(i)) newLeaves.push_back(_leaves[i]);
          newLeaves.insert(newLeaves.end(), to_add.begin(), to_add.end());
          _leaves.swap(newLeaves);
          // defer _node2leafpos rebuild to your pass boundary if you like
        }

        return coarsened;
      }

// ---- geometry map ----
      Point2 parent_to_physical(Point2 xi) const {
        require_geometry();
        double N9[9];
        Shapes::Q9(xi, N9);
        double X = 0, Y = 0;
        for (int a = 0; a < 9; a++) {
          X += N9[a] * _X[a];
          Y += N9[a] * _Y[a];
        }
        return {X, Y};
      }

      inline void post_topology_update() {
        // Keep leaves in Z-order
        std::sort(_leaves.begin(), _leaves.end(),
        [&](u32 a, u32 b) {
          return _tree_nodes[a].morton < _tree_nodes[b].morton;
        });

        rebuild_leafpos_lookup();
        rebuild_connectivity_active_bases();
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
      std::vector<Point2>
      extract_node_parent_coords_in_level_range(u32 lev_min, u32 lev_max, Basis basis) const {
        std::vector<Point2> coords;

        for (u32 leaf : leaf_indices()) {
          u32 lev = level_of(leaf);          // use accessor, not _nodes
          if (lev < lev_min || lev > lev_max) continue;

          std::vector<Point2> xi;
          extract_leaf_parent_coords(basis, leaf, xi);

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
      u32 locate_leaf_on_parent(Point2 xi) const {
        if (xi[0] < -1.0 || xi[0] > 1.0 || xi[1] < -1.0 || xi[1] > 1.0) return npos32;

        // Scale to integer grid at max depth
        const double scale = double(1u << _maxDepth) / 2.0;  // maps [-1,1] -> [0, 2^maxDepth)
        u32 ix = std::min<u32>((u32)((xi[0] + 1.0) * scale), (1u << _maxDepth) - 1);
        u32 iy = std::min<u32>((u32)((xi[1] + 1.0) * scale), (1u << _maxDepth) - 1);

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
      bool evaluate_field_on_parent(u32 fid, Point2 xi, double& value) const {
        if (fid >= _fields.size()) return false;

        u32 leaf_node_idx;
        Point2 shat;

        // 1) locate leaf and local reference coords
        if (!locate_leaf_on_parent_and_ref(xi, leaf_node_idx, shat)) {
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
          Shapes::Q1(shat, N4);
          value = 0.0;
          for (int a = 0; a < 4; ++a) value += N4[a] * f.nodal[conn[a]];
        }
        break;
        case Basis::Serendipity8: {
          double N8[8];
          Shapes::Q8(shat, N8);
          value = 0.0;
          for (int a = 0; a < 8; ++a) value += N8[a] * f.nodal[conn[a]];
        }
        break;
        case Basis::Q2_Quad9: {
          double N9[9];
          Shapes::Q9(shat, N9);
          value = 0.0;
          for (int a = 0; a < 9; ++a) value += N9[a] * f.nodal[conn[a]];
        }
        break;
        }
        return true;
      }


      inline void local_ref_fast(u32 leaf_node_idx,
                                 Point2 xi,
                                 Point2 &shat) const {
        const TreeNode& n = _tree_nodes[leaf_node_idx];
        u32 ix = n.ix;
        u32 iy = n.iy;
        //u32 ix, iy; deinterleave2(n.morton, ix, iy);
        const double N = double(1u << n.level);

        shat[0] = N * xi[0] + N - double((ix << 1) + 1);
        shat[1] = N * xi[1] + N - double((iy << 1) + 1);
      }

      inline bool locate_leaf_on_parent_and_ref(Point2 xi,
                                                u32& leaf_node_idx,
                                                Point2& shat) const {
        // 1) Find the leaf containing (xi,eta) in parent space.
        leaf_node_idx = locate_leaf_on_parent(xi);
        if (leaf_node_idx == npos32) return false;

        // 2) Convert to local [-1,1]^2 using level + Morton (no divisions).
        local_ref_fast(leaf_node_idx, xi, shat);

        return true;
      }


      // Refine tree so that each physical point lies in a leaf of at least 'maxDepthTarget'.
// Caches ONLY the inverse map (xi, eta) per point; uses _node2leafpos instead of per-pass node_to_pos.
      void refine_to_contain_points(const std::vector<Point2>& pts,
                                    u32 maxDepthTarget) {
        if (pts.empty()) return;

        // ---- Cache inverse map (xi, eta) once for all points ----
        struct XiEta {
          Point2 xi;
          bool valid;
        };
        std::vector<XiEta> pm(pts.size());
        for (size_t i = 0; i < pts.size(); ++i) {
          Point2 xi;
          if (!inverse_map_quad9(pts[i], xi)) {
            pm[i] = {{0.0, 0.0}, false};
          }
          else {
            pm[i] = {xi, true};
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

            // Reuse cached (xi, eta); do only the tree lookup per pass
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

          // Minimal predicate: only checks the per-pass mark
          auto should_refine_idx = [&](u32 leaf_pos) -> bool {
            return mark_by_leaf[leaf_pos] != 0;
          };


          // Count how many leaves are marked (we’ll refine at most that many)
          u32 to_split = 0;
          for (u32 pos = 0; pos < nleaves; ++pos) if (mark_by_leaf[pos]) ++to_split;

// Net new nodes per split: +3 tree nodes (4 kids minus parent as leaf)
          if (to_split) {
            _tree_nodes.reserve(_tree_nodes.size() + to_split * 3u);
            _leaves.reserve(_leaves.size() + to_split * 3u); // parent replaced by 4 children → +3 leaves
          }



          const std::size_t nsplit = refine_pass_min(should_refine_idx);
          if (nsplit == 0) break; // reached max depth locally or topology unchanged

          // _leaves changed inside refine_pass_min; _node2leafpos will be rebuilt
          // at the top of the next loop iteration.
        }

        // Keep mesh 1-irregular and refresh bookkeeping
        enforce_balance();
        post_topology_update();
      }



// Perform one refinement pass with a minimal predicate: bool(u32 leaf_pos)
// No geometry probing; fastest path.
      template<class RefinePred>
      std::size_t refine_pass_min(RefinePred&& should_refine) {
        std::vector<u32> newLeaves;
        newLeaves.reserve(_leaves.size() * 2);
        std::size_t refined = 0;

        const u32 n0 = (u32)leaf_count(); // snapshot leaf count

        for (u32 idx = 0; idx < n0; ++idx) {
          const u32 leaf_node_idx = _leaves[idx];             // index into _tree_nodes
          const TreeNode leaf_copy = _tree_nodes[leaf_node_idx]; // snapshot

          if (leaf_copy.is_leaf && leaf_copy.level < _maxDepth && should_refine(idx)) {
            // split into 4 children

            u32 ix = leaf_copy.ix;
            u32 iy = leaf_copy.iy;
            //u32 ix, iy; deinterleave2(leaf_copy.morton, ix, iy);

            for (int dy = 0; dy < 2; ++dy) {
              for (int dx = 0; dx < 2; ++dx) {
                const u64 cmorton = interleave2(2 * ix + dx, 2 * iy + dy);

                TreeNode child;
                child.ix = 2 * ix + dx;
                child.iy = 2 * iy + dy;
                child.morton  = cmorton;
                child.level   = leaf_copy.level + 1;
                child.is_leaf = true;
                child.parent  = leaf_node_idx;
                child.child[0] = child.child[1] = child.child[2] = child.child[3] = npos32;

                const u32 child_idx = (u32)_tree_nodes.size();
                _tree_nodes.push_back(child);
                _tree_nodes[leaf_node_idx].child[(dy << 1) | dx] = child_idx;
                newLeaves.push_back(child_idx);
              }
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

// Enforce 1-irregularity: adjacent leaves may differ by at most 1 level.
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
        q.buf.reserve(_leaves.size() * 5);   // or exactly _leaves.size()
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

          for (int dir = 0; dir < 4; ++dir) {
            u32 nb = neighbor_leaf_by_dir_any(leaf, dir);
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
                for (int d = 0; d < 4; ++d) {
                  u32 nb2 = neighbor_leaf_by_dir_any(nb, d);
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
                  u32 nb2 = neighbor_leaf_by_dir_any(leaf, d);
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
      void assign_from(const QuadTree2D& rhs) {
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

        _root       = rhs._root;
        _geom_ready = rhs._geom_ready;

        // ---- geometria ----
        _X = rhs._X;
        _Y = rhs._Y;

        // ---- basi attive e campi ----
        copy_uset(_activeBases, rhs._activeBases);
        copy_vec(_fields,      rhs._fields);

        // ---- topologia ----
        copy_vec(_tree_nodes,   rhs._tree_nodes);
        copy_vec(_leaves,       rhs._leaves);
        copy_vec(_leaf_ids,     rhs._leaf_ids);
        copy_vec(_node2leafpos, rhs._node2leafpos);

        // ---- registri per-basis ----
        for (int b = 0; b < 3; ++b) {
          copy_umap(_basisReg[b].nodeMap,   rhs._basisReg[b].nodeMap);
          copy_vec(_basisReg[b].nodes,     rhs._basisReg[b].nodes);
          copy_vecvec_int(_basisReg[b].elem2glob, rhs._basisReg[b].elem2glob);
        }
      }



// Conservative coarsen cycle using snapshot + parent coords; rebuild all fields
      std::size_t coarsen_only_cycle_safe(u32 fid,
                                          double tau_coarse,
                                          QuadTree2D& snapshot,
                                          u32 max_passes = 10,
                                          Basis probe_basis = Basis::Q2_Quad9) {
        // Freeze current state for conservative evaluation/transfer
        snapshot.assign_from(*this);

        // Coarsen predicate (evaluated on the snapshot, in parent coords)
        auto pred = [&](u32 /*parent_morton*/, u32 level,
                        const std::vector<Point2>& pts_xi,
                        const std::vector<Point2>& /*pts_xy*/,
        const std::vector<std::array<double, 9>>& /*Nvals*/) -> bool {
          if (level <= min_depth()) return false;
          if (pts_xi.empty()) return false;

          double v0;
          if (!snapshot.evaluate_field_on_parent(fid, pts_xi[0], v0))
            return false;

          double mn = v0, mx = v0;
          for (size_t i = 1; i < pts_xi.size(); ++i) {
            double val;
            if (snapshot.evaluate_field_on_parent(fid, pts_xi[i], val)) {
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

      inline u64 node_key_from_parent(Point2 xi) const {
        const u32 nodeBits = _maxDepth + 1;
        const u32 nodesN   = (1u << nodeBits);

        auto to_idx = [nodesN](double s)->u32 {
          if (s <= -1.0) return 0u;
          if (s >=  1.0) return nodesN;
          double t = (s + 1.0) * double(nodesN) * 0.5;
          long long li = llround(t);
          if (li < 0) li = 0;
          if (li > (long long)nodesN) li = (long long)nodesN;
          return (u32)li;
        };

        const u32 ix = to_idx(xi[0]);
        const u32 iy = to_idx(xi[1]);
        return (u64(ix) << 32) | u64(iy);
      }

// Rebuild field fid on *this* from source tree 'src' by sampling at parent nodes (global nodal storage)
      void rebuild_field_from(const QuadTree2D& src, u32 fid) {
        assert(fid < _fields.size() && fid < src._fields.size());
        Field& dst = _fields[fid];
        const Field& s = src._fields[fid];
        assert(dst.basis == s.basis);

        const BasisRegistry& Rdst = _basisReg[(int)dst.basis];
        const BasisRegistry& Rsrc = src._basisReg[(int)s.basis];

        dst.nodal.resize(Rdst.nodes.size());

        for (size_t gid = 0; gid < Rdst.nodes.size(); ++gid) {
          const auto& pr = Rdst.nodes[gid].parent;               // (xi,eta) at destination
          const u64 key  = node_key_from_parent({pr[0], pr[1]});   // same quantization as src
          auto it = Rsrc.nodeMap.find(key);
          if (it != Rsrc.nodeMap.end()) {
            dst.nodal[gid] = s.nodal[it->second];                // direct copy (fast path)
          }
          else {
            double val;                                          // rare fallback
            const bool ok = src.evaluate_field_on_parent(fid, {pr[0], pr[1]}, val);
            dst.nodal[gid] = ok ? val : std::numeric_limits<double>::quiet_NaN();
          }
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

        Point2 q;
        switch (dir) {
        case 0:
          q[0] = xi0 - epsx;
          q[1] = ym;
          break; // left
        case 1:
          q[0] = xi1 + epsx;
          q[1] = ym;
          break; // right
        case 2:
          q[0] = xm;
          q[1] = eta0 - epsy;
          break; // down
        default:/*3*/
          q[0] = xm;
          q[1] = eta1 + epsy;
          break; // up
        }

        // outside global domain? no neighbor
        if (q[0] < -1.0 || q[0] > 1.0 || q[1] < -1.0 || q[1] > 1.0) return npos32;

        return locate_leaf_on_parent(q);
      }

// dir: 0=left (−x), 1=right (+x), 2=down (−y), 3=up (+y).
// Returns leaf node index or npos32 if no neighbor (domain boundary).
      u32 neighbor_leaf_by_dir_any(u32 leaf_node_idx, int dir) const {
        const auto side_dx = (dir == 0 ? -1 : dir == 1 ? +1 : 0);
        const auto side_dy = (dir == 2 ? -1 : dir == 3 ? +1 : 0);

        // 1) Fast same-parent sibling if not on that side boundary.
        u32 nidx = leaf_node_idx;
        const TreeNode& n = _tree_nodes[nidx];
        u32 ix = n.ix, iy = n.iy;
        u32 L  = n.level;

        if (L > 0) {
          const u32 dx_bit = ix & 1u;
          const u32 dy_bit = iy & 1u;
          const u32 child_code = (dy_bit << 1) | dx_bit;

          // Is there a sibling across 'dir' inside the same parent?
          bool has_local = false;
          u32 sib_code = child_code;

          if (dir == 0 && dx_bit == 1) {
            sib_code = (child_code & ~1u) | 0u;  // left -> dx 1->0
            has_local = true;
          }
          if (dir == 1 && dx_bit == 0) {
            sib_code = (child_code & ~1u) | 1u;  // right -> dx 0->1
            has_local = true;
          }
          if (dir == 2 && dy_bit == 1) {
            sib_code = (child_code & ~2u) | 0u;  // down -> dy 1->0
            has_local = true;
          }
          if (dir == 3 && dy_bit == 0) {
            sib_code = (child_code & ~2u) | 2u;  // up   -> dy 0->1
            has_local = true;
          }

          if (has_local) {
            u32 p = n.parent;
            if (p != npos32) {
              u32 cand = _tree_nodes[p].child[sib_code];
              if (cand != npos32) {
                // descend along the shared edge to the leaf that contains the edge midpoint
                return descend_to_edge_midpoint_leaf(cand, dir, leaf_node_idx);
              }
            }
          }
        }

        // 2) No same-parent neighbor: climb until we can cross the edge at that ancestor.
        u32 cur = leaf_node_idx;
        while (true) {
          const TreeNode& c = _tree_nodes[cur];
          if (c.parent == npos32) return npos32; // hit domain boundary

          const TreeNode& p = _tree_nodes[c.parent];

          // child position (within parent)
          const u32 dx_bit = _tree_nodes[cur].ix & 1u;
          const u32 dy_bit = _tree_nodes[cur].iy & 1u;
          const u32 child_code = (dy_bit << 1) | dx_bit;

          bool can_cross_here = false;
          u32 sib_code = child_code;

          if (dir == 0 && dx_bit == 1) {
            sib_code = (child_code & ~1u) | 0u;
            can_cross_here = true;
          }
          if (dir == 1 && dx_bit == 0) {
            sib_code = (child_code & ~1u) | 1u;
            can_cross_here = true;
          }
          if (dir == 2 && dy_bit == 1) {
            sib_code = (child_code & ~2u) | 0u;
            can_cross_here = true;
          }
          if (dir == 3 && dy_bit == 0) {
            sib_code = (child_code & ~2u) | 2u;
            can_cross_here = true;
          }

          if (can_cross_here) {
            u32 sib = p.child[sib_code];
            if (sib == npos32) return npos32; // logically shouldn’t happen if parent exists, but guard
            // 3) We have the neighbor subtree root at some level; now descend to the leaf
            return descend_to_edge_midpoint_leaf(sib, dir, leaf_node_idx);
          }

          // Keep climbing
          cur = c.parent;
        }
      }

// Descend inside 'subtree_root' to the leaf that contains the shared-edge midpoint
// relative to 'ref_leaf_idx'. Works for any level difference (no 1-irregularity needed).
      u32 descend_to_edge_midpoint_leaf(u32 subtree_root, int dir, u32 ref_leaf_idx) const {
        // Compute the shared-edge midpoint in parent coordinates (cheap; uses bounds).
        double xi0, eta0, xi1, eta1;
        leaf_bounds(_tree_nodes[ref_leaf_idx], xi0, eta0, xi1, eta1);
        Point2 q;
        switch (dir) {
        case 0: q = {xi0, 0.5 * (eta0 + eta1)}; break; // left  edge midpoint
        case 1: q = {xi1, 0.5 * (eta0 + eta1)}; break; // right edge midpoint
        case 2: q = {0.5 * (xi0 + xi1), eta0}; break; // down  edge midpoint
        default: q = {0.5 * (xi0 + xi1), eta1}; break; // up    edge midpoint
        }

        // Walk down inside the neighbor subtree toward q.
        u32 node = subtree_root;
        while (true) {
          const TreeNode& n = _tree_nodes[node];
          if (n.is_leaf) return node;

          // Choose child quadrant by comparing q to this node's center.
          double bx0, by0, bx1, by1;
          leaf_bounds(n, bx0, by0, bx1, by1);
          const double cx = 0.5 * (bx0 + bx1);
          const double cy = 0.5 * (by0 + by1);

          const int qx = (q[0] >= cx) ? 1 : 0;
          const int qy = (q[1] >= cy) ? 1 : 0;
          const int code = (qy << 1) | qx;

          u32 child = n.child[code];
          if (child == npos32) {
            // Neighbor subtree isn’t refined this deep; current node is the leaf w.r.t. mesh.
            return node;
          }
          node = child;
        }
      }
















// Refine a single leaf node into 4 children, update _leaves in place.
// Returns true if a refinement actually occurred.
      bool refine_leaf_once(u32 leaf_node_idx) {
        if (leaf_node_idx == npos32) return false;

        const TreeNode parent_copy = _tree_nodes[leaf_node_idx];
        if (!parent_copy.is_leaf || parent_copy.level >= _maxDepth) return false;

        // Try O(1) position if cache is valid; else fallback scan.
        u32 pos = npos32;
        if (_leafpos_valid &&
            leaf_node_idx < _node2leafpos.size() &&
            _node2leafpos[leaf_node_idx] < _leaves.size() &&
            _leaves[_node2leafpos[leaf_node_idx]] == leaf_node_idx) {
          pos = _node2leafpos[leaf_node_idx];
        }
        else {
          // Rare slow path; still O(n) but infrequent.
          for (u32 i = 0; i < (u32)_leaves.size(); ++i) {
            if (_leaves[i] == leaf_node_idx) {
              pos = i;
              break;
            }
          }
          if (pos == npos32) return false; // not a leaf anymore
        }

        // Parent becomes internal
        _tree_nodes[leaf_node_idx].is_leaf = false;

        // Remove parent from leaves (swap-erase). Do NOT touch the map here.
        const u32 last = _leaves.back();
        _leaves[pos] = last;
        _leaves.pop_back();

        // Children
        u32 ix = parent_copy.ix;
        u32 iy = parent_copy.iy;
        //u32 ix, iy; deinterleave2(parent_copy.morton, ix, iy);

        // Optional micro-opt: reserve a little to reduce realloc churn (not required)
        // if (_tree_nodes.capacity() < _tree_nodes.size() + 4) _tree_nodes.reserve(_tree_nodes.size() + 4);
        // if (_leaves.capacity()     < _leaves.size()     + 4) _leaves.reserve(_leaves.size()     + 4);

        for (int dy = 0; dy < 2; ++dy) {
          for (int dx = 0; dx < 2; ++dx) {
            const int q = (dy << 1) | dx;
            const u64 cm = interleave2(2 * ix + (u32)dx, 2 * iy + (u32)dy);

            TreeNode child;
            child.ix = 2 * ix + (u32)dx;
            child.iy = 2 * iy + (u32)dy;
            child.morton  = cm;
            child.level   = parent_copy.level + 1;
            child.is_leaf = true;
            child.parent  = leaf_node_idx;
            child.child[0] = child.child[1] = child.child[2] = child.child[3] = npos32;

            const u32 cidx = (u32)_tree_nodes.size();
            _tree_nodes.push_back(child);                  // parent index remains valid
            _tree_nodes[leaf_node_idx].child[q] = cidx;    // write by index, safe
            _leaves.push_back(cidx);                       // no map writes
          }
        }

        // Map is now stale; mark dirty and rebuild later when needed.
        _leafpos_valid = false;
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

              u32 ix = node_copy.ix;
              u32 iy = node_copy.iy;
              //u32 ix, iy; deinterleave2(node_copy.morton, ix, iy);

              _tree_nodes[leaf_idx].is_leaf = false;

              for (int dy = 0; dy < 2; ++dy) {
                for (int dx = 0; dx < 2; ++dx) {
                  u64 cmorton = interleave2(2 * ix + dx, 2 * iy + dy);
                  u32 cindex = (u32)_tree_nodes.size();

                  TreeNode child;
                  child.ix = 2 * ix + dx;
                  child.iy = 2 * iy + dy;
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
        u32 ix = node.ix;
        u32 iy = node.iy;
        //u32 ix, iy; deinterleave2(node.morton, ix, iy);

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
        R.nodes.push_back(FEMNode{gid, {xi, eta}, parent_to_physical({xi, eta})});
        R.nodeMap.emplace(key, gid);
        return gid;
      }

      void rebuild_connectivity(Basis b) {
        BasisRegistry &R = _basisReg[(int)b];
        R.clear();

        // --- NEW: reserve to avoid rehash ---
        const u32 nleaves = leaf_count();
        const int per = (b == Basis::Q1_Quad4 ? 4 : (b == Basis::Serendipity8 ? 8 : 9));

// A decent estimate for unique nodes is ~ (per * leaves)*0.6..0.8 due to sharing.
// Keep it conservative but not tiny to minimize rehashes:
        const size_t est_nodes = std::max<size_t>(32, size_t(nleaves) * size_t(per) * size_t(3) / size_t(2)); // 1.5x leaves*per
        R.nodeMap.reserve(est_nodes);             // avoids _M_rehash spikes
        R.nodes.reserve(est_nodes);               // avoids _M_realloc_insert spikes
        R.elem2glob.reserve(nleaves);             // avoid vector growth



        std::vector<Point2> xi;
        for (u32 e = 0; e < leaf_count(); ++e) {
          extract_leaf_parent_coords(b, e, xi);
          std::vector<int> conn;
          conn.reserve(xi.size());

          for (auto &p : xi) {
            int gid = get_or_insert_gid(R, p[0], p[1]); // now uses registry
            conn.push_back(gid);
          }
          R.elem2glob.push_back(std::move(conn));
        }
      }

// Q9 inverse with Q1 warm-up; same control flow & loops (no unrolling)
      bool inverse_map_quad9(Point2 x, Point2 &xi, int maxIts = 25, double tol = 1e-12) const {
        require_geometry();

        // start at center
        xi = {0.0, 0.0};
        const double tol2 = tol * tol;

        // buffers
        double N[9], dNdxi[9], dNdeta[9];

        // --- Warm-up with Q1: residual and Jacobian from corners only ---
        Shapes::Q1(xi, N);           // N[0..3] valid
        double Xp = 0.0, Yp = 0.0;
        for (int a = 0; a < 4; ++a) {
          Xp += N[a] * _X[a];
          Yp += N[a] * _Y[a];
        }
        double rx = Xp - x[0], ry = Yp - x[1];

        Shapes::Q1_dN(xi, dNdxi, dNdeta); // fills [0..3]
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
          xi[0]  = std::max(-1.0, std::min(1.0, xi[0]  - dXi));
          xi[1] = std::max(-1.0, std::min(1.0, xi[1] - dEta));

          // recompute N at updated (xi,eta), update residual for convergence check
          Shapes::Q9(xi, N);
          Xp = 0.0; Yp = 0.0;
          for (int a = 0; a < 9; ++a) {
            Xp += N[a] * _X[a];
            Yp += N[a] * _Y[a];
          }
          rx = Xp - x[0];  ry = Yp - x[1];
          if ((rx * rx + ry * ry) < tol2) return true;

          // build NEXT Jacobian at updated (xi,eta) using Q9 derivatives
          Shapes::Q9_dN(xi, dNdxi, dNdeta);
          dXdxi = dXdeta = dYdxi = dYdeta = 0.0;   // reset accumulators
          for (int a = 0; a < 9; ++a) {
            dXdxi  += dNdxi[a] * _X[a];  dXdeta += dNdeta[a] * _X[a];
            dYdxi  += dNdxi[a] * _Y[a];  dYdeta += dNdeta[a] * _Y[a];
          }
        }

        return false; // not converged
      }


      // --- private: Q1 (Quad4) inverse map (Newton in parent space), using Shapes::Q1 ---
      bool inverse_map_quad4(Point2 x, Point2 &xi,
                             int maxIts = 12, double tol = 1e-12) const {
        require_geometry();

        // Start at center
        xi = {0.0, 0.0};
        const double tol2    = tol * tol;
        const double tol_dx2 = 1e-20; // ~1e-10 step

        for (int it = 0; it < maxIts; ++it) {

          // Values from your Shapes::Q1
          double N4[4], dNdxi[4], dNdeta[4];
          Shapes::Q1(xi, N4);
          Shapes::Q1_dN(xi, dNdxi, dNdeta);

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

          double rx = Xp - x[0], ry = Yp - x[1];
          double nrm2 = (rx * rx + ry * ry);
          if (nrm2 < tol2) return true;

          double detJ = dXdxi * dYdeta - dXdeta * dYdxi;
          if (std::abs(detJ) < 1e-20) break; // singular

          // Newton step: [dXi, dEta] = J^{-1} * r
          double dXi  = (dYdeta * rx - dXdeta * ry) / detJ;
          double dEta = (-dYdxi * rx + dXdxi * ry) / detJ;

          xi[0] -= dXi;
          xi[1] -= dEta;

          // keep iterates reasonable (light clamp)
          xi[0] = std::max(-1.5, std::min(1.5, xi[0]));
          xi[1] = std::max(-1.5, std::min(1.5, xi[1]));
        }
        return false; // not converged
      }



// --- private: split one leaf and update _leaves ---
      bool refine_leaf_and_update_leaves(u32 leaf_node_idx) {
        if (leaf_node_idx == npos32) return false;
        TreeNode& parent = _tree_nodes[leaf_node_idx];
        if (!parent.is_leaf || parent.level >= _maxDepth) return false;

        u32 ix = parent.ix;
        u32 iy = parent.iy;
        //u32 ix, iy; deinterleave2(parent.morton, ix, iy);

        for (int dy = 0; dy < 2; ++dy) {
          for (int dx = 0; dx < 2; ++dx) {
            const u64 cmorton = interleave2(2 * ix + dx, 2 * iy + dy);
            TreeNode child;
            child.ix = 2 * ix + dx;
            child.iy = 2 * iy + dy;
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




      inline void rebuild_leafpos_lookup() {
        _node2leafpos.assign(_tree_nodes.size(), npos32);
        for (u32 i = 0; i < _leaves.size(); ++i) _node2leafpos[_leaves[i]] = i;
        _leafpos_valid = true;
      }


    private: //data
// config
      u32  _maxDepth;
      u32  _minDepth{0};
      bool _allowCoarsenBelowMinDepth{true};
      bool _leafpos_valid{false};


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

      mutable std::vector<u32> _node2leafpos; // size == _tree_nodes.size(), npos32 default
      //std::unordered_map<u64, u32> _leaf_by_morton; // morton -> node index


      struct BasisHasher {
        size_t operator()(fem::Basis b) const noexcept {
          return std::hash<uint8_t>()(static_cast<uint8_t>(b));
        }
      };
      std::unordered_set<Basis, BasisHasher> _activeBases;
  };

} // namespace fem






