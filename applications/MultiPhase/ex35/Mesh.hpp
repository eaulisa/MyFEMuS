// #pragma once
//
// #include <vector>
#pragma once

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <limits>
#include <string>

class Mesh {
  public:
    using ElemConnectivity = std::vector<std::vector<unsigned>>;
    using ElemType         = std::vector<unsigned>;
    using ElemLevel        = std::vector<unsigned>;
    using Coordinates      = std::vector<std::vector<double>>;
    using AMRVector        = std::vector<int>;                   // -1,0,1 or 0/1, etc.
    using FatherVector     = std::vector<unsigned>;              // index on coarser level
    using ChildrenVector   = std::vector<std::vector<unsigned>>; // indices on finer level
    using NodeElemAdj      = std::vector<std::vector<unsigned>>; // node -> elements
    using ElemNeighbors    = std::vector<std::vector<unsigned>>; // neighbors[e][f]

    enum class CheckMode {
      Basic,  // cheap global checks
      Full    // Basic + heavier hierarchy checks
    };

    static constexpr unsigned UMAX = std::numeric_limits<unsigned>::max();

    Mesh(ElemConnectivity &elTplgy,
         ElemType         &elType,
         ElemLevel        &elLevel,
         Coordinates      &X,
         AMRVector        &AMR,
         FatherVector     &father,
         ChildrenVector   &children,
         NodeElemAdj      &nodeElems,
         ElemNeighbors    &neighbors)
      : _elTplgy(elTplgy),
        _elType(elType),
        _elLevel(elLevel),
        _X(X),
        _AMR(AMR),
        _father(father),
        _children(children),
        _nodeElems(nodeElems),
        _neighbors(neighbors) {
      checkSanity("Mesh::Mesh", CheckMode::Full);
    }


    void clearAllData() {
      _elTplgy.clear();
      _elType.clear();
      _elLevel.clear();
      _X.clear();
      _AMR.clear();
      _father.clear();
      _children.clear();
      _nodeElems.clear();
      _neighbors.clear(); // if you added neighbors
    }


    // ==================================
    // Assignment: copy DATA (not binding)
    // ==================================
    Mesh& operator=(const Mesh &other) {
      if (this != &other) {
        _elTplgy    = other._elTplgy;    // topology
        _elType     = other._elType;     // element types
        _elLevel    = other._elLevel;    // refinement levels
        _X          = other._X;          // coordinates
        _AMR        = other._AMR;        // AMR flags
        _father     = other._father;     // father indices (coarser level)
        _children   = other._children;   // children indices (finer level)
        _nodeElems  = other._nodeElems;  // node->elements adjacency
        _neighbors  = other._neighbors;  // face neighbors

        checkSanity("Mesh::operator=", CheckMode::Full);
      }
      return *this;
    }

    // Optional: friend swap (swaps underlying containers)
    friend void swap(Mesh &a, Mesh &b) noexcept {
      using std::swap;
      swap(a._elTplgy,    b._elTplgy);
      swap(a._elType,     b._elType);
      swap(a._elLevel,    b._elLevel);
      swap(a._X,          b._X);
      swap(a._AMR,        b._AMR);
      swap(a._father,     b._father);
      swap(a._children,   b._children);
      swap(a._nodeElems,  b._nodeElems);
      swap(a._neighbors,  b._neighbors);
    }

    // ============================
    // Getters (non-const)
    // ============================
    ElemConnectivity & elTplgy()   {
      return _elTplgy;
    }
    ElemType         & elType()    {
      return _elType;
    }
    ElemLevel        & elLevel()   {
      return _elLevel;
    }
    Coordinates      & X()         {
      return _X;
    }
    AMRVector        & AMR()       {
      return _AMR;
    }
    FatherVector     & father()    {
      return _father;
    }
    ChildrenVector   & children()  {
      return _children;
    }
    NodeElemAdj      & nodeElems() {
      return _nodeElems;
    }
    ElemNeighbors    & neighbors() {
      return _neighbors;
    }

    // ============================
    // Getters (const)
    // ============================
    const ElemConnectivity & elTplgy()   const {
      return _elTplgy;
    }
    const ElemType         & elType()    const {
      return _elType;
    }
    const ElemLevel        & elLevel()   const {
      return _elLevel;
    }
    const Coordinates      & X()         const {
      return _X;
    }
    const AMRVector        & AMR()       const {
      return _AMR;
    }
    const FatherVector     & father()    const {
      return _father;
    }
    const ChildrenVector   & children()  const {
      return _children;
    }
    const NodeElemAdj      & nodeElems() const {
      return _nodeElems;
    }
    const ElemNeighbors    & neighbors() const {
      return _neighbors;
    }

    // ============================
    // Mesh info interface
    // ============================

    std::size_t dim() const {
      return _X.size();
    }

    std::size_t numNodes() const {
      return _X.empty() ? 0 : _X[0].size();
    }

    std::size_t numElements() const {
      return _elType.size();
    }

    // ============================
    // AMR interface
    // ============================

    /// Uniform refinement: AMR cleared => everyone refined
    void setUniformRefinement() {
      _AMR.clear();  // AMR.size()==0 => uniform refinement
    }

    /// Debug: refine only even elements, keep odd ones
    /// AMR[i] = 1 (refine) for even i
    /// AMR[i] = 0 (keep)   for odd i
    void setRefineEvenElements() {
      const std::size_t nEl = numElements();
      _AMR.assign(nEl, 0);
      for (std::size_t i = 0; i < nEl; ++i) {
        if ((i % 2u) == 0u) {
          _AMR[i] = 1;
        }
      }
      checkSanity("Mesh::setRefineEvenElements", CheckMode::Basic);
    }

    void setRefinementFromBallLevelSetCrossing_OneRing(const std::vector<double>& center,
                                                       double r,
                                                       unsigned neighMode,   // 0,1,2,3 (see below)
                                                       double eps = -1.0) {
      // neighMode:
      // 0 = all connectivity nodes
      // 1 = vertices only
      // 2 = faces only (3D face-centers / 2D edge-mids / 1D endpoints)
      // 3 = hybrid: Hex27 & Quad9 -> faces; Tet15, Wedge21, Tri7 -> vertices; Line3 -> vertices

      const std::size_t d    = dim();
      const std::size_t nEl  = numElements();
      const std::size_t nNod = numNodes();

      if (neighMode > 3u) {
        throw std::runtime_error("setRefinementFromBallLevelSetCrossing_OneRing: neighMode must be 0,1,2,3");
      }
      if (d == 0) {
        throw std::runtime_error("setRefinementFromBallLevelSetCrossing_OneRing: dim()==0");
      }
      if (center.size() != d) {
        std::cout << center.size() << std::endl;
        throw std::runtime_error("setRefinementFromBallLevelSetCrossing_OneRing: center.size() != dim()");
      }
      if (!(r > 0.0)) {
        throw std::runtime_error("setRefinementFromBallLevelSetCrossing_OneRing: r must be > 0");
      }

      if (eps < 0.0) {
        eps = 1e-12 * std::max(1.0, std::abs(r));
      }

      if (_nodeElems.size() != nNod) {
        throw std::runtime_error(
          "setRefinementFromBallLevelSetCrossing_OneRing: _nodeElems is not built or wrong size");
      }

      _AMR.assign(nEl, 0);

      auto phi_of_node = [&](unsigned node) -> double {
        if (node >= nNod) {
          throw std::runtime_error("setRefinementFromBallLevelSetCrossing_OneRing: node index out of range");
        }
        double s2 = 0.0;
        for (std::size_t k = 0; k < d; ++k) {
          const double dx = _X[k][node] - center[k];
          s2 += dx * dx;
        }
        return std::sqrt(s2) - r;
      };

      auto vertex_range = [&](std::size_t connSize) -> std::pair<std::size_t, std::size_t> {
        switch (connSize) {
        case 27: return std::pair<std::size_t, std::size_t>(0, 8); // Hex27 vertices 0..7
        case 15: return std::pair<std::size_t, std::size_t>(0, 4); // Tet15 vertices 0..3
        case 21: return std::pair<std::size_t, std::size_t>(0, 6); // Wedge21 vertices 0..5
        case 9:  return std::pair<std::size_t, std::size_t>(0, 4); // Quad9 vertices 0..3
        case 7:  return std::pair<std::size_t, std::size_t>(0, 3); // Tri7 vertices 0..2
        case 3:  return std::pair<std::size_t, std::size_t>(0, 2); // Line3 vertices 0..1
        default:
          throw std::runtime_error(
            "setRefinementFromBallLevelSetCrossing_OneRing: unsupported element connectivity size (vertices)");
        }
      };

      auto face_range = [&](std::size_t connSize) -> std::pair<std::size_t, std::size_t> {
        switch (connSize) {
        case 27: return std::pair<std::size_t, std::size_t>(20, 26); // Hex27 face centers 20..25
        case 15: return std::pair<std::size_t, std::size_t>(10, 14); // Tet15 face centers 10..13
        case 21: return std::pair<std::size_t, std::size_t>(15, 20); // Wedge21 face centers 15..19
        case 9:  return std::pair<std::size_t, std::size_t>(4,  8);  // Quad9 edge mids 4..7
        case 7:  return std::pair<std::size_t, std::size_t>(3,  6);  // Tri7 edge mids 3..5
        case 3:  return std::pair<std::size_t, std::size_t>(0,  2);  // Line3 endpoints 0..1 (treat as "faces")
        default:
          throw std::runtime_error(
            "setRefinementFromBallLevelSetCrossing_OneRing: unsupported element connectivity size (faces)");
        }
      };

      auto pick_range = [&](std::size_t connSize, unsigned mode) -> std::pair<std::size_t, std::size_t> {
        if (mode == 0u) {
          return std::pair<std::size_t, std::size_t>(0, connSize); // all connectivity nodes
        }
        if (mode == 1u) {
          return vertex_range(connSize); // vertices only
        }
        if (mode == 2u) {
          return face_range(connSize); // faces only (as defined above)
        }
        // mode == 3 hybrid:
        // Hex27 & Quad9 -> faces
        // Tet15, Wedge21, Tri7 -> vertices
        // Line3 -> vertices
        if (connSize == 27 || connSize == 9) {
          return face_range(connSize);
        }
        if (connSize == 15 || connSize == 21 || connSize == 7 || connSize == 3) {
          return vertex_range(connSize);
        }
        throw std::runtime_error(
          "setRefinementFromBallLevelSetCrossing_OneRing: unsupported element connectivity size (hybrid)");
      };

      // -------------------------
      // 1) First sweep: crossing elements => 2
      // -------------------------
      for (std::size_t e = 0; e < nEl; ++e) {
        const auto& conn = _elTplgy[e];
        if (conn.empty()) continue;

        bool hasNeg = false;
        bool hasPos = false;
        bool hasZeroBand = false;

        for (unsigned node : conn) {
          const double phi = phi_of_node(node);

          if (phi <= -eps) hasNeg = true;
          else if (phi >= +eps) hasPos = true;
          else hasZeroBand = true;
        }

        const bool crosses =
          (hasNeg && hasPos) ||
          (hasZeroBand && (hasNeg || hasPos));

        if (crosses) _AMR[e] = 2;
      }

      // -------------------------
      // 2) Second sweep: expand from 2-elements to neighbors via _nodeElems
      //    - only 2-elements seed the halo (no cascading from new 1s)
      //    - neighbor candidates defined by neighMode range
      // -------------------------
      for (std::size_t e = 0; e < nEl; ++e) {
        if (_AMR[e] != 2) continue;

        const auto& conn = _elTplgy[e];
        if (conn.empty()) continue;

        const std::pair<std::size_t, std::size_t> liRange = pick_range(conn.size(), neighMode);
        const std::size_t liBeg = liRange.first;
        const std::size_t liEnd = liRange.second;

        if (liEnd > conn.size()) {
          throw std::runtime_error(
            "setRefinementFromBallLevelSetCrossing_OneRing: chosen local range exceeds connectivity size");
        }

        for (std::size_t li = liBeg; li < liEnd; ++li) {
          const unsigned gnode = conn[li];
          if (gnode >= nNod) {
            throw std::runtime_error(
              "setRefinementFromBallLevelSetCrossing_OneRing: global node out of range");
          }

          const auto& incident = _nodeElems[gnode];
          for (unsigned ne_u : incident) {
            const std::size_t ne = static_cast<std::size_t>(ne_u);

            if (ne >= nEl) {
              throw std::runtime_error(
                "setRefinementFromBallLevelSetCrossing_OneRing: _nodeElems contains element out of range");
            }
            if (_AMR[ne] == 0) _AMR[ne] = 1;
          }
        }
      }

      // -------------------------
      // 3) Final sweep: everything nonzero -> 1
      // -------------------------
      for (std::size_t e = 0; e < nEl; ++e) {
        if (_AMR[e] != 0) _AMR[e] = 1;
      }

      checkSanity("Mesh::setRefinementFromBallLevelSetCrossing_OneRing", CheckMode::Basic);
    }

    /// Query: does element i need refinement?
    /// - if AMR.size() == 0  -> uniform refinement: return true
    /// - else                -> return (AMR[i] == 1)
    bool needsRefinement(std::size_t i) const {
      if (_AMR.empty()) {
        return true;  // uniform refinement
      }
      if (i >= _AMR.size()) {
        throw std::out_of_range("Mesh::needsRefinement: index out of range");
      }
      return (_AMR[i] == 1);
    }

    // ============================
    // Hierarchy interface (cross-level)
    // ============================

    /// Sentinel value meaning "no father" (for level 0)
    static unsigned noFather() {
      return std::numeric_limits<unsigned>::max();
    }

    /// Set children of element i (indices live on the *finer* mesh).
    /// Option 2 semantics:
    ///   - ch.size() == 1  -> unrefined (identity child)
    ///   - ch.size() == 4  -> refined in 2D
    ///   - ch.size() == 8  -> refined in 3D
    void setChildren(std::size_t i, const std::vector<unsigned> &ch) {
      if (i >= numElements()) {
        throw std::out_of_range("Mesh::setChildren: index out of range");
      }

      if (_children.size() != numElements()) {
        _children.resize(numElements());
      }
      _children[i] = ch;
    }

    /// Access children of element i (const)
    const std::vector<unsigned>& childrenOf(std::size_t i) const {
      if (i >= numElements()) {
        throw std::out_of_range("Mesh::childrenOf: index out of range");
      }
      return _children[i];
    }

    /// Access children of element i (non-const)
    std::vector<unsigned>& childrenOf(std::size_t i) {
      if (i >= numElements()) {
        throw std::out_of_range("Mesh::childrenOf: index out of range");
      }
      return _children[i];
    }

    /// Access father of element i (index on coarser level or noFather())
    unsigned fatherOf(std::size_t i) const {
      if (i >= numElements()) {
        throw std::out_of_range("Mesh::fatherOf: index out of range");
      }
      return _father[i];
    }

    /// Clear all children information for this mesh level.
    void clearAllChildren() {
      _children.clear();
    }

    /// Reset all fathers to `noFather()` for this mesh level.
    /// Use this only when you *know* this mesh has no coarser level
    /// (i.e., it's the global coarsest mesh).
    void resetAllFathersToNoFather() {
      const std::size_t nEl = numElements();
      const unsigned no_father = noFather();
      _father.assign(nEl, no_father);
    }

    // ============================
    // Overall mesh check (heavy)
    // ============================
    void checkAll(const char* where = "Mesh::checkAll") const {
      checkSanity(where, CheckMode::Full);
    }

    // ============================
    // Node -> elements adjacency
    // ============================

    /// Build node -> elements adjacency into the class data _nodeElems:
    /// _nodeElems[n] = list of element indices that contain node n.
    void buildNodeToElementAdjacency() {
      const std::size_t nElems  = numElements();
      const std::size_t nNodes  = numNodes();

      _nodeElems.clear();
      _nodeElems.resize(nNodes);

      // Pre-count to reserve exact sizes (more cache-friendly).
      std::vector<unsigned> counts(nNodes, 0u);

      for (std::size_t e = 0; e < nElems; ++e) {
        const auto &conn = _elTplgy[e];
        for (unsigned node : conn) {
          if (node >= nNodes) {
            throw std::runtime_error("Mesh::buildNodeToElementAdjacency: node index out of range (count)");
          }
          ++counts[node];
        }
      }

      for (std::size_t n = 0; n < nNodes; ++n) {
        _nodeElems[n].reserve(counts[n]);
      }

      // Fill adjacency
      for (std::size_t e = 0; e < nElems; ++e) {
        const auto &conn = _elTplgy[e];
        for (unsigned node : conn) {
          if (node >= nNodes) {
            throw std::runtime_error("Mesh::buildNodeToElementAdjacency: node index out of range (fill)");
          }
          _nodeElems[node].push_back(static_cast<unsigned>(e));
        }
      }
    }

    // ============================================================
    // Build face neighbors via face-center nodes
    // ============================================================
    //
    // Conventions (by connectivity size):
    //   Hex27   (27): 6 faces, centers conn[20..25]
    //   Tet15   (15): 4 faces, centers conn[10..13]
    //   Wedge21 (21): 5 faces, centers conn[15..19]
    //   Quad9    (9): 4 edges, centers conn[ 4.. 7]
    //   Tri7     (7): 3 edges, centers conn[ 3.. 5]
    //   Line3    (3): 2 "faces" endpoints conn[0],conn[1]
    //
    // For such a face-center node n:
    //   _nodeElems[n] has size 1 (boundary) or 2 (interior).
    //   one entry is e, the other (if present) is the neighbor.
    void buildFaceNeighborsFromNodeToElement() {
      const std::size_t nEl  = numElements();
      const std::size_t nNod = numNodes();

      if (_nodeElems.size() != nNod) {
        throw std::runtime_error(
          "Mesh::buildFaceNeighborsFromNodeToElement: _nodeElems is not built or wrong size");
      }

      _neighbors.clear();
      _neighbors.resize(nEl);

      for (std::size_t e = 0; e < nEl; ++e) {
        const auto& conn = _elTplgy[e];
        if (conn.empty()) {
          _neighbors[e].clear();
          continue;
        }

        const unsigned nf = numFacesFromConnSize(conn.size());
        const unsigned i0 = faceCenterStartFromConnSize(conn.size());

        if (static_cast<std::size_t>(i0 + nf) > conn.size()) {
          throw std::runtime_error("Mesh::buildFaceNeighborsFromNodeToElement: conn too small for face-center range");
        }

        _neighbors[e].assign(nf, UMAX);

        for (unsigned f = 0; f < nf; ++f) {
          const unsigned node = conn[i0 + f];
          if (node >= nNod) {
            throw std::runtime_error("Mesh::buildFaceNeighborsFromNodeToElement: face-center node out of range");
          }

          const auto& inc = _nodeElems[node];

          if (inc.size() == 1u) {
            // boundary face: should be itself
            if (inc[0] != static_cast<unsigned>(e)) {
              throw std::runtime_error("Mesh::buildFaceNeighborsFromNodeToElement: boundary face-center not owned by e");
            }
            _neighbors[e][f] = UMAX;
          }
          else if (inc.size() == 2u) {
            const unsigned a0 = inc[0], a1 = inc[1];
            if (a0 == static_cast<unsigned>(e)) _neighbors[e][f] = a1;
            else if (a1 == static_cast<unsigned>(e)) _neighbors[e][f] = a0;
            else {
              throw std::runtime_error("Mesh::buildFaceNeighborsFromNodeToElement: face-center does not reference e");
            }
          }
          else {
            throw std::runtime_error("Mesh::buildFaceNeighborsFromNodeToElement: face-center shared by >2 elements");
          }
        }
      }

      checkSanity("Mesh::buildFaceNeighborsFromNodeToElement", CheckMode::Basic);
    }

    /// Enforce 1-level discontinuity for the *next* mesh level by
    /// adjusting the AMR flags (assumed 0 or 1).
    ///
    /// After this call, for all elements ei, ej that share at least one node,
    /// we have: | (level[ei] + AMR[ei]) - (level[ej] + AMR[ej]) | <= 1.
    ///
    /// This may turn some AMR[e] from 0 to 1 (propagating refinement),
    /// but never from 1 to 0.
    void adjustAMRForOneLevelDiscontinuity() {
      const std::size_t nElems  = numElements();
      const std::size_t nNodes  = numNodes();

      // AMR empty => uniform refinement; if current mesh is 1-level balanced,
      // the next mesh is, too.
      if (_AMR.empty()) {
        return;
      }

      // Sanity: AMR size and values
      if (_AMR.size() != nElems) {
        throw std::runtime_error(
          "adjustAMRForOneLevelDiscontinuity: AMR.size() != numElements()");
      }
      for (std::size_t e = 0; e < nElems; ++e) {
        if (_AMR[e] != 0 && _AMR[e] != 1) {
          throw std::runtime_error(
            "adjustAMRForOneLevelDiscontinuity: AMR entries must be 0 or 1");
        }
      }

      auto &nodeElems = _nodeElems;

      bool changedGlobal;
      do {
        changedGlobal = false;

        // For each node (patch)
        for (std::size_t n = 0; n < nNodes; ++n) {
          auto &els = nodeElems[n];
          const std::size_t m = els.size();
          if (m < 2) continue; // no neighbors at this node

          // Compute Lmin and Lmax for this patch
          unsigned Lmin = std::numeric_limits<unsigned>::max();
          unsigned Lmax = 0u;

          for (std::size_t k = 0; k < m; ++k) {
            const unsigned e = els[k];
            const unsigned L = _elLevel[e] + static_cast<unsigned>(_AMR[e]);
            if (L < Lmin) Lmin = L;
            if (L > Lmax) Lmax = L;
          }

          if (Lmax <= Lmin + 1u) {
            // Patch already satisfies 1-level condition for next mesh
            continue;
          }
          else if (Lmax > Lmin + 2u) {
            throw std::runtime_error("adjustAMRForOneLevelDiscontinuity: "
                                     "Lmax-Lmin > 2 at a node (inconsistent input)");
          }

          // Here Lmax == Lmin + 2 is the only possible violation.
          // Refine all elements at Lmin (if AMR==0); single pass is enough.
          bool refinedSomething = false;
          for (std::size_t k = 0; k < m; ++k) {
            const unsigned e = els[k];
            const unsigned L = _elLevel[e] + static_cast<unsigned>(_AMR[e]);
            if (L == Lmin && _AMR[e] == 0) {
              _AMR[e] = 1;
              refinedSomething = true;
              changedGlobal    = true;
            }
          }

          if (!refinedSomething) {
            throw std::runtime_error(
              "adjustAMRForOneLevelDiscontinuity: cannot fix Lmax-Lmin>1 at a node");
          }
        }

      }
      while (changedGlobal);
    }

  private:
    ElemConnectivity &_elTplgy;
    ElemType         &_elType;
    ElemLevel        &_elLevel;
    Coordinates      &_X;
    AMRVector        &_AMR;
    FatherVector     &_father;     // index on level-1 (or noFather)
    ChildrenVector   &_children;   // indices on level+1
    NodeElemAdj      &_nodeElems;  // node -> elements adjacency
    ElemNeighbors    &_neighbors;  // NEW: element -> face neighbors

    static unsigned numFacesFromConnSize(std::size_t connSize) {
      switch (connSize) {
      case 27: return 6; // Hex27
      case 15: return 4; // Tet15
      case 21: return 5; // Wedge21
      case 9:  return 4; // Quad9
      case 7:  return 3; // Tri7
      case 3:  return 2; // Line3
      default:
        throw std::runtime_error("Mesh::numFacesFromConnSize: unsupported element connectivity size");
      }
    }

    static unsigned faceCenterStartFromConnSize(std::size_t connSize) {
      switch (connSize) {
      case 27: return 20; // Hex27 face centers 20..25
      case 15: return 10; // Tet15 face centers 10..13
      case 21: return 15; // Wedge21 face centers 15..19
      case 9:  return 4;  // Quad9 edge mids 4..7
      case 7:  return 3;  // Tri7 edge mids 3..5
      case 3:  return 0;  // Line3 endpoints 0..1
      default:
        throw std::runtime_error("Mesh::faceCenterStartFromConnSize: unsupported element connectivity size");
      }
    }

    // ============================
    // Sanity checks
    // ============================
    void checkSanity(const char *where,
                     CheckMode mode = CheckMode::Basic) const {
      const std::size_t nEl = _elType.size();

      // ========== BASIC CHECKS (always run) ==========

      // 1) element array sizes
      if (_elLevel.size() != nEl || _elTplgy.size() != nEl) {
        std::cerr << where << ": element arrays size mismatch: "
                  << "elType="    << _elType.size()
                  << ", elLevel=" << _elLevel.size()
                  << ", elTplgy=" << _elTplgy.size() << "\n";
        throw std::runtime_error("Mesh sanity check failed: element arrays size mismatch");
      }

      // 1b) AMR size: either 0 (uniform) or nEl
      if (!_AMR.empty() && _AMR.size() != nEl) {
        std::cerr << where << ": AMR size invalid: "
                  << "AMR.size() = " << _AMR.size()
                  << ", numElements = " << nEl << "\n";
        throw std::runtime_error("Mesh sanity check failed: AMR size must be 0 or numElements");
      }

      // 1c) neighbors size: either 0 (not built) or nEl
      if (!_neighbors.empty() && _neighbors.size() != nEl) {
        std::cerr << where << ": neighbors size invalid: "
                  << "neighbors.size() = " << _neighbors.size()
                  << ", numElements = " << nEl << "\n";
        throw std::runtime_error("Mesh sanity check failed: neighbors size must be 0 or numElements");
      }

      // 2) coordinate sizes
      const std::size_t d = dim();
      if (d > 0) {
        const std::size_t nNodes0 = _X[0].size();
        for (std::size_t k = 1; k < d; ++k) {
          if (_X[k].size() != nNodes0) {
            std::cerr << where << ": coordinate arrays size mismatch: "
                      << "X[" << k << "].size() = " << _X[k].size()
                      << ", X[0].size() = " << nNodes0 << "\n";
            throw std::runtime_error("Mesh sanity check failed: coordinate arrays size mismatch");
          }
        }
      }

      // 3) dim vs element type (only checks first type; extend if needed)
      if (!_elType.empty()) {
        const unsigned t  = _elType[0];
        const std::size_t d0 = dim();
        if (t < 3) { // 3D: Hex27, Tet15, Wedge21
          if (d0 != 3) {
            std::cerr << where << ": element type " << t
                      << " implies dim=3, but dim()=" << d0 << "\n";
            throw std::runtime_error("Mesh sanity check failed: 3D element with non-3D coordinates");
          }
        }
        else if (t < 5) {   // 2D: Quad9, Tri7
          if (d0 != 2) {
            std::cerr << where << ": element type " << t
                      << " implies dim=2, but dim()=" << d0 << "\n";
            throw std::runtime_error("Mesh sanity check failed: 2D element with non-2D coordinates");
          }
        }
      }

      if (mode == CheckMode::Basic) {
        return; // skip heavy hierarchy checks
      }

      // ========== FULL / HEAVY CHECKS (hierarchy etc.) ==========

      // 4) hierarchy container sizes: either 0 or nEl
      if (!_father.empty() && _father.size() != nEl) {
        std::cerr << where << ": father size mismatch: "
                  << "father.size() = " << _father.size()
                  << ", numElements = " << nEl << "\n";
        throw std::runtime_error("Mesh sanity check failed: father size must be 0 or numElements");
      }
      if (!_children.empty() && _children.size() != nEl) {
        std::cerr << where << ": children size mismatch: "
                  << "children.size() = " << _children.size()
                  << ", numElements = " << nEl << "\n";
        throw std::runtime_error("Mesh sanity check failed: children size must be 0 or numElements");
      }

      // 4b) neighbors per-element size check (if built)
      if (_neighbors.size() == nEl && !_neighbors.empty()) {
        for (std::size_t e = 0; e < nEl; ++e) {
          const auto& conn = _elTplgy[e];
          if (conn.empty()) continue;
          const unsigned nf = numFacesFromConnSize(conn.size());
          if (_neighbors[e].size() != nf) {
            std::cerr << where << ": neighbors[" << e << "].size()=" << _neighbors[e].size()
                      << " but expected " << nf << " for conn.size()=" << conn.size() << "\n";
            throw std::runtime_error("Mesh sanity check failed: wrong neighbors[e] size");
          }
        }
      }

      // 5) Hierarchy invariants (local, not cross-level range checks)
      if (_children.size() == nEl && !_children.empty()) {
        const std::size_t d0 = dim();
        const unsigned expectedRefChildren =
          (d0 == 2 ? 4u :
           d0 == 3 ? 8u : 0u);

        for (std::size_t e = 0; e < nEl; ++e) {
          const auto &ch = _children[e];
          if (!ch.empty() && expectedRefChildren != 0u) {
            const std::size_t k = ch.size();
            // Option 2: 1 child (unrefined) or 4/8 children (refined)
            if (k != 1u && k != expectedRefChildren) {
              std::cerr << where << ": element " << e
                        << " has " << k
                        << " children; expected 1 or " << expectedRefChildren << "\n";
              throw std::runtime_error("Mesh sanity check failed: wrong number of children");
            }
          }
        }
      }

      if (_father.size() == nEl && !_father.empty()) {
        for (std::size_t e = 0; e < nEl; ++e) {
          const unsigned lvl = _elLevel[e];
          if (lvl > 0 && _father[e] == noFather()) {
            std::cerr << where << ": element " << e
                      << " has level " << lvl
                      << " but father[e] == noFather()\n";
            throw std::runtime_error(
              "Mesh sanity check failed: element with level>0 has no father");
          }
        }
      }
    }
};





// #include <algorithm>
// #include <stdexcept>
// #include <iostream>
// #include <limits>
// #include <string>
//
// class Mesh {
//   public:
//     using ElemConnectivity = std::vector<std::vector<unsigned>>;
//     using ElemType         = std::vector<unsigned>;
//     using ElemLevel        = std::vector<unsigned>;
//     using Coordinates      = std::vector<std::vector<double>>;
//     using AMRVector        = std::vector<int>;                 // -1,0,1 or 0/1, etc.
//     using FatherVector     = std::vector<unsigned>;            // index on coarser level
//     using ChildrenVector   = std::vector<std::vector<unsigned>>; // indices on finer level
//     using NodeElemAdj      = std::vector<std::vector<unsigned>>; // node -> elements
//
//     enum class CheckMode {
//       Basic,  // cheap global checks
//       Full    // Basic + heavier hierarchy checks
//     };
//
//     Mesh(ElemConnectivity &elTplgy,
//          ElemType         &elType,
//          ElemLevel        &elLevel,
//          Coordinates      &X,
//          AMRVector        &AMR,
//          FatherVector     &father,
//          ChildrenVector   &children,
//          NodeElemAdj      &nodeElems)
//       : _elTplgy(elTplgy),
//         _elType(elType),
//         _elLevel(elLevel),
//         _X(X),
//         _AMR(AMR),
//         _father(father),
//         _children(children),
//         _nodeElems(nodeElems) {
//       checkSanity("Mesh::Mesh", CheckMode::Full);
//     }
//
//     // ==================================
//     // Assignment: copy DATA (not binding)
//     // ==================================
//     Mesh& operator=(const Mesh &other) {
//       if (this != &other) {
//         _elTplgy   = other._elTplgy;    // topology
//         _elType    = other._elType;     // element types
//         _elLevel   = other._elLevel;    // refinement levels
//         _X         = other._X;          // coordinates
//         _AMR       = other._AMR;        // AMR flags
//         _father    = other._father;     // father indices (coarser level)
//         _children  = other._children;   // children indices (finer level)
//         _nodeElems = other._nodeElems;  // node->elements adjacency
//
//         checkSanity("Mesh::operator=", CheckMode::Full);
//       }
//       return *this;
//     }
//
//     // Optional: friend swap (swaps underlying containers)
//     friend void swap(Mesh &a, Mesh &b) noexcept {
//       using std::swap;
//       swap(a._elTplgy,   b._elTplgy);
//       swap(a._elType,    b._elType);
//       swap(a._elLevel,   b._elLevel);
//       swap(a._X,         b._X);
//       swap(a._AMR,       b._AMR);
//       swap(a._father,    b._father);
//       swap(a._children,  b._children);
//       swap(a._nodeElems, b._nodeElems);
//     }
//
//     // ============================
//     // Getters (non-const)
//     // ============================
//     ElemConnectivity & elTplgy()   {
//       return _elTplgy;
//     }
//     ElemType         & elType()    {
//       return _elType;
//     }
//     ElemLevel        & elLevel()   {
//       return _elLevel;
//     }
//     Coordinates      & X()         {
//       return _X;
//     }
//     AMRVector        & AMR()       {
//       return _AMR;
//     }
//     FatherVector     & father()    {
//       return _father;
//     }
//     ChildrenVector   & children()  {
//       return _children;
//     }
//     NodeElemAdj      & nodeElems() {
//       return _nodeElems;
//     }
//
//     // ============================
//     // Getters (const)
//     // ============================
//     const ElemConnectivity & elTplgy()   const {
//       return _elTplgy;
//     }
//     const ElemType         & elType()    const {
//       return _elType;
//     }
//     const ElemLevel        & elLevel()   const {
//       return _elLevel;
//     }
//     const Coordinates      & X()         const {
//       return _X;
//     }
//     const AMRVector        & AMR()       const {
//       return _AMR;
//     }
//     const FatherVector     & father()    const {
//       return _father;
//     }
//     const ChildrenVector   & children()  const {
//       return _children;
//     }
//     const NodeElemAdj      & nodeElems() const {
//       return _nodeElems;
//     }
//
//     // ============================
//     // Mesh info interface
//     // ============================
//
//     std::size_t dim() const {
//       return _X.size();
//     }
//
//     std::size_t numNodes() const {
//       return _X.empty() ? 0 : _X[0].size();
//     }
//
//     std::size_t numElements() const {
//       return _elType.size();
//     }
//
//     // ============================
//     // AMR interface
//     // ============================
//
//     /// Uniform refinement: AMR cleared => everyone refined
//     void setUniformRefinement() {
//       _AMR.clear();  // AMR.size()==0 => uniform refinement
//     }
//
//     /// Debug: refine only even elements, keep odd ones
//     /// AMR[i] = 1 (refine) for even i
//     /// AMR[i] = 0 (keep)   for odd i
//     void setRefineEvenElements() {
//       const std::size_t nEl = numElements();
//       _AMR.assign(nEl, 0);
//       for (std::size_t i = 0; i < nEl; ++i) {
//         if ((i % 2u) == 0u) {
//           _AMR[i] = 1;
//         }
//       }
//       checkSanity("Mesh::setRefineEvenElements", CheckMode::Basic);
//     }
//
//
//     void setRefinementFromBallLevelSetCrossing_OneRing(const std::vector<double>& center,
//                                                        double r,
//                                                        unsigned neighMode,   // 0,1,2,3 (see below)
//                                                        double eps = -1.0) {
//       // neighMode:
//       // 0 = all connectivity nodes
//       // 1 = vertices only
//       // 2 = faces only (3D face-centers / 2D edge-mids / 1D endpoints)
//       // 3 = hybrid: Hex27 & Quad9 -> faces; Tet15, Wedge21, Tri7 -> vertices; Line3 -> vertices
//
//       const std::size_t d    = dim();
//       const std::size_t nEl  = numElements();
//       const std::size_t nNod = numNodes();
//
//       if (neighMode < 0 || neighMode > 3) {
//         throw std::runtime_error("setRefinementFromBallLevelSetCrossing_OneRing: neighMode must be 0,1,2,3");
//       }
//       if (d == 0) {
//         throw std::runtime_error("setRefinementFromBallLevelSetCrossing_OneRing: dim()==0");
//       }
//       if (center.size() != d) {
//         std::cout << center.size() << std::endl;
//         throw std::runtime_error("setRefinementFromBallLevelSetCrossing_OneRing: center.size() != dim()");
//       }
//       if (!(r > 0.0)) {
//         throw std::runtime_error("setRefinementFromBallLevelSetCrossing_OneRing: r must be > 0");
//       }
//
//       if (eps < 0.0) {
//         eps = 1e-12 * std::max(1.0, std::abs(r));
//       }
//
//       if (_nodeElems.size() != nNod) {
//         throw std::runtime_error(
//           "setRefinementFromBallLevelSetCrossing_OneRing: _nodeElems is not built or wrong size");
//       }
//
//       _AMR.assign(nEl, 0);
//
//       auto phi_of_node = [&](unsigned node) -> double {
//         if (node >= nNod) {
//           throw std::runtime_error("setRefinementFromBallLevelSetCrossing_OneRing: node index out of range");
//         }
//         double s2 = 0.0;
//         for (std::size_t k = 0; k < d; ++k) {
//           const double dx = _X[k][node] - center[k];
//           s2 += dx * dx;
//         }
//         return std::sqrt(s2) - r;
//       };
//
//       auto vertex_range = [&](std::size_t connSize) -> std::pair<std::size_t, std::size_t> {
//         switch (connSize) {
//         case 27: return std::pair<std::size_t, std::size_t>(0, 8); // Hex27 vertices 0..7
//         case 15: return std::pair<std::size_t, std::size_t>(0, 4); // Tet15 vertices 0..3
//         case 21: return std::pair<std::size_t, std::size_t>(0, 6); // Wedge21 vertices 0..5
//         case 9:  return std::pair<std::size_t, std::size_t>(0, 4); // Quad9 vertices 0..3
//         case 7:  return std::pair<std::size_t, std::size_t>(0, 3); // Tri7 vertices 0..2
//         case 3:  return std::pair<std::size_t, std::size_t>(0, 2); // Line3 vertices 0..1
//         default:
//           throw std::runtime_error(
//             "setRefinementFromBallLevelSetCrossing_OneRing: unsupported element connectivity size (vertices)");
//         }
//       };
//
//       auto face_range = [&](std::size_t connSize) -> std::pair<std::size_t, std::size_t> {
//         switch (connSize) {
//         case 27: return std::pair<std::size_t, std::size_t>(20, 26); // Hex27 face centers 20..25
//         case 15: return std::pair<std::size_t, std::size_t>(10, 14); // Tet15 face centers 10..13
//         case 21: return std::pair<std::size_t, std::size_t>(15, 20); // Wedge21 face centers 15..19
//         case 9:  return std::pair<std::size_t, std::size_t>(4,  8);  // Quad9 edge mids 4..7
//         case 7:  return std::pair<std::size_t, std::size_t>(3,  6);  // Tri7 edge mids 3..5
//         case 3:  return std::pair<std::size_t, std::size_t>(0,  2);  // Line3 endpoints 0..1 (treat as "faces")
//         default:
//           throw std::runtime_error(
//             "setRefinementFromBallLevelSetCrossing_OneRing: unsupported element connectivity size (faces)");
//         }
//       };
//
//       auto pick_range = [&](std::size_t connSize, int mode) -> std::pair<std::size_t, std::size_t> {
//         if (mode == 0) {
//           return std::pair<std::size_t, std::size_t>(0, 0); // no neighbor
//         }
//         if (mode == 1) {
//           return vertex_range(connSize); // vertices only
//         }
//         if (mode == 2) {
//           return face_range(connSize); // faces only (as defined above)
//         }
//         // mode == 3 hybrid:
//         // Hex27 & Quad9 -> faces
//         // Tet15, Wedge21, Tri7 -> vertices
//         // Line3 -> vertices (your call; endpoints also fine, but this matches "tri only on vertices")
//         if (connSize == 27 || connSize == 9) {
//           return face_range(connSize);
//         }
//         if (connSize == 15 || connSize == 21 || connSize == 7 || connSize == 3) {
//           return vertex_range(connSize);
//         }
//         throw std::runtime_error(
//           "setRefinementFromBallLevelSetCrossing_OneRing: unsupported element connectivity size (hybrid)");
//       };
//
//       // -------------------------
//       // 1) First sweep: crossing elements => 2
//       // -------------------------
//       for (std::size_t e = 0; e < nEl; ++e) {
//         const auto& conn = _elTplgy[e];
//         if (conn.empty()) continue;
//
//         bool hasNeg = false;
//         bool hasPos = false;
//         bool hasZeroBand = false;
//
//         for (unsigned node : conn) {
//           const double phi = phi_of_node(node);
//
//           if (phi <= -eps) hasNeg = true;
//           else if (phi >= +eps) hasPos = true;
//           else hasZeroBand = true;
//         }
//
//         const bool crosses =
//           (hasNeg && hasPos) ||
//           (hasZeroBand && (hasNeg || hasPos));
//
//         if (crosses) _AMR[e] = 2;
//       }
//
//       // -------------------------
//       // 2) Second sweep: expand from 2-elements to neighbors via _nodeElems
//       //    - only 2-elements seed the halo (no cascading from new 1s)
//       //    - neighbor candidates defined by neighMode range
//       // -------------------------
//       for (std::size_t e = 0; e < nEl; ++e) {
//         if (_AMR[e] != 2) continue;
//
//         const auto& conn = _elTplgy[e];
//         if (conn.empty()) continue;
//
//         const std::pair<std::size_t, std::size_t> liRange = pick_range(conn.size(), neighMode);
//         const std::size_t liBeg = liRange.first;
//         const std::size_t liEnd = liRange.second;
//
//         if (liEnd > conn.size()) {
//           throw std::runtime_error(
//             "setRefinementFromBallLevelSetCrossing_OneRing: chosen local range exceeds connectivity size");
//         }
//
//         for (std::size_t li = liBeg; li < liEnd; ++li) {
//           const unsigned gnode = conn[li];
//           if (gnode >= nNod) {
//             throw std::runtime_error(
//               "setRefinementFromBallLevelSetCrossing_OneRing: global node out of range");
//           }
//
//           const auto& incident = _nodeElems[gnode];
//           for (unsigned ne_u : incident) {
//             const std::size_t ne = static_cast<std::size_t>(ne_u);
//
//             if (ne >= nEl) {
//               throw std::runtime_error(
//                 "setRefinementFromBallLevelSetCrossing_OneRing: _nodeElems contains element out of range");
//             }
//             if (_AMR[ne] == 0) _AMR[ne] = 1;
//           }
//         }
//       }
//
//       // -------------------------
//       // 3) Final sweep: everything nonzero -> 1
//       // -------------------------
//       for (std::size_t e = 0; e < nEl; ++e) {
//         if (_AMR[e] != 0) _AMR[e] = 1;
//       }
//
//       checkSanity("Mesh::setRefinementFromBallLevelSetCrossing_OneRing", CheckMode::Basic);
//     }
//
//
//     /// Query: does element i need refinement?
//     /// - if AMR.size() == 0  -> uniform refinement: return true
//     /// - else                -> return (AMR[i] == 1)
//     bool needsRefinement(std::size_t i) const {
//       if (_AMR.empty()) {
//         return true;  // uniform refinement
//       }
//       if (i >= _AMR.size()) {
//         throw std::out_of_range("Mesh::needsRefinement: index out of range");
//       }
//       return (_AMR[i] == 1);
//     }
//
//     // ============================
//     // Hierarchy interface (cross-level)
//     // ============================
//
//     /// Sentinel value meaning "no father" (for level 0)
//     static unsigned noFather() {
//       return std::numeric_limits<unsigned>::max();
//     }
//
//     /// Set children of element i (indices live on the *finer* mesh).
//     /// Option 2 semantics:
//     ///   - ch.size() == 1  -> unrefined (identity child)
//     ///   - ch.size() == 4  -> refined in 2D
//     ///   - ch.size() == 8  -> refined in 3D
//     void setChildren(std::size_t i, const std::vector<unsigned> &ch) {
//       if (i >= numElements()) {
//         throw std::out_of_range("Mesh::setChildren: index out of range");
//       }
//
//       if (_children.size() != numElements()) {
//         _children.resize(numElements());
//       }
//       _children[i] = ch;
//     }
//
//     /// Access children of element i (const)
//     const std::vector<unsigned>& childrenOf(std::size_t i) const {
//       if (i >= numElements()) {
//         throw std::out_of_range("Mesh::childrenOf: index out of range");
//       }
//       return _children[i];
//     }
//
//     /// Access children of element i (non-const)
//     std::vector<unsigned>& childrenOf(std::size_t i) {
//       if (i >= numElements()) {
//         throw std::out_of_range("Mesh::childrenOf: index out of range");
//       }
//       return _children[i];
//     }
//
//     /// Access father of element i (index on coarser level or noFather())
//     unsigned fatherOf(std::size_t i) const {
//       if (i >= numElements()) {
//         throw std::out_of_range("Mesh::fatherOf: index out of range");
//       }
//       return _father[i];
//     }
//
//     /// Clear all children information for this mesh level.
//     void clearAllChildren() {
//       _children.clear();
//     }
//
//     /// Reset all fathers to `noFather()` for this mesh level.
//     /// Use this only when you *know* this mesh has no coarser level
//     /// (i.e., it's the global coarsest mesh).
//     void resetAllFathersToNoFather() {
//       const std::size_t nEl = numElements();
//       const unsigned no_father = noFather();
//       _father.assign(nEl, no_father);
//     }
//
//     // ============================
//     // Overall mesh check (heavy)
//     // ============================
//     void checkAll(const char* where = "Mesh::checkAll") const {
//       checkSanity(where, CheckMode::Full);
//     }
//
//     // ============================
//     // Node -> elements adjacency
//     // ============================
//
//     /// Build node -> elements adjacency into the class data _nodeElems:
//     /// _nodeElems[n] = list of element indices that contain node n.
//     void buildNodeToElementAdjacency() {
//       const std::size_t nElems  = numElements();
//       const std::size_t nNodes  = numNodes();
//
//       _nodeElems.clear();
//       _nodeElems.resize(nNodes);
//
//       // Pre-count to reserve exact sizes (more cache-friendly).
//       std::vector<unsigned> counts(nNodes, 0u);
//
//       for (std::size_t e = 0; e < nElems; ++e) {
//         const auto &conn = _elTplgy[e];
//         for (unsigned node : conn) {
//           if (node >= nNodes) {
//             throw std::runtime_error("Mesh::buildNodeToElementAdjacency: node index out of range (count)");
//           }
//           ++counts[node];
//         }
//       }
//
//       for (std::size_t n = 0; n < nNodes; ++n) {
//         _nodeElems[n].reserve(counts[n]);
//       }
//
//       // Fill adjacency
//       for (std::size_t e = 0; e < nElems; ++e) {
//         const auto &conn = _elTplgy[e];
//         for (unsigned node : conn) {
//           if (node >= nNodes) {
//             throw std::runtime_error("Mesh::buildNodeToElementAdjacency: node index out of range (fill)");
//           }
//           _nodeElems[node].push_back(static_cast<unsigned>(e));
//         }
//       }
//     }
//
//     /// Enforce 1-level discontinuity for the *next* mesh level by
//     /// adjusting the AMR flags (assumed 0 or 1).
//     ///
//     /// After this call, for all elements ei, ej that share at least one node,
//     /// we have: | (level[ei] + AMR[ei]) - (level[ej] + AMR[ej]) | <= 1.
//     ///
//     /// This may turn some AMR[e] from 0 to 1 (propagating refinement),
//     /// but never from 1 to 0.
//     void adjustAMRForOneLevelDiscontinuity() {
//       const std::size_t nElems  = numElements();
//       const std::size_t nNodes  = numNodes();
//
//       // AMR empty => uniform refinement; if current mesh is 1-level balanced,
//       // the next mesh is, too.
//       if (_AMR.empty()) {
//         return;
//       }
//
//       // Sanity: AMR size and values
//       if (_AMR.size() != nElems) {
//         throw std::runtime_error(
//           "adjustAMRForOneLevelDiscontinuity: AMR.size() != numElements()");
//       }
//       for (std::size_t e = 0; e < nElems; ++e) {
//         if (_AMR[e] != 0 && _AMR[e] != 1) {
//           throw std::runtime_error(
//             "adjustAMRForOneLevelDiscontinuity: AMR entries must be 0 or 1");
//         }
//       }
//
//       auto &nodeElems = _nodeElems;
//
//       bool changedGlobal;
//       do {
//         changedGlobal = false;
//
//         // For each node (patch)
//         for (std::size_t n = 0; n < nNodes; ++n) {
//           auto &els = nodeElems[n];
//           const std::size_t m = els.size();
//           if (m < 2) continue; // no neighbors at this node
//
//           // Compute Lmin and Lmax for this patch
//           unsigned Lmin = std::numeric_limits<unsigned>::max();
//           unsigned Lmax = 0u;
//
//           for (std::size_t k = 0; k < m; ++k) {
//             const unsigned e = els[k];
//             const unsigned L = _elLevel[e] + static_cast<unsigned>(_AMR[e]);
//             if (L < Lmin) Lmin = L;
//             if (L > Lmax) Lmax = L;
//           }
//
//           if (Lmax <= Lmin + 1u) {
//             // Patch already satisfies 1-level condition for next mesh
//             continue;
//           }
//           else if (Lmax > Lmin + 2u) {
//             throw std::runtime_error("adjustAMRForOneLevelDiscontinuity: "
//                                      "Lmax-Lmin > 2 at a node (inconsistent input)");
//           }
//
//           // Here Lmax == Lmin + 2 is the only possible violation.
//           // Refine all elements at Lmin (if AMR==0); single pass is enough.
//           bool refinedSomething = false;
//           for (std::size_t k = 0; k < m; ++k) {
//             const unsigned e = els[k];
//             const unsigned L = _elLevel[e] + static_cast<unsigned>(_AMR[e]);
//             if (L == Lmin && _AMR[e] == 0) {
//               _AMR[e] = 1;
//               refinedSomething = true;
//               changedGlobal    = true;
//             }
//           }
//
//           if (!refinedSomething) {
//             // In theory this shouldn’t happen under your assumptions
//             // (you should always be able to raise Lmin by refining some AMR=0).
//             throw std::runtime_error(
//               "adjustAMRForOneLevelDiscontinuity: cannot fix Lmax-Lmin>1 at a node");
//           }
//         }
//
//       }
//       while (changedGlobal);
//     }
//
//   private:
//     ElemConnectivity &_elTplgy;
//     ElemType         &_elType;
//     ElemLevel        &_elLevel;
//     Coordinates      &_X;
//     AMRVector        &_AMR;
//     FatherVector     &_father;     // index on level-1 (or noFather)
//     ChildrenVector   &_children;   // indices on level+1
//     NodeElemAdj      &_nodeElems;  // node -> elements adjacency
//
//     // ============================
//     // Sanity checks
//     // ============================
//     void checkSanity(const char *where,
//                      CheckMode mode = CheckMode::Basic) const {
//       const std::size_t nEl = _elType.size();
//
//       // ========== BASIC CHECKS (always run) ==========
//
//       // 1) element array sizes
//       if (_elLevel.size() != nEl || _elTplgy.size() != nEl) {
//         std::cerr << where << ": element arrays size mismatch: "
//                   << "elType="    << _elType.size()
//                   << ", elLevel=" << _elLevel.size()
//                   << ", elTplgy=" << _elTplgy.size() << "\n";
//         throw std::runtime_error("Mesh sanity check failed: element arrays size mismatch");
//       }
//
//       // 1b) AMR size: either 0 (uniform) or nEl
//       if (!_AMR.empty() && _AMR.size() != nEl) {
//         std::cerr << where << ": AMR size invalid: "
//                   << "AMR.size() = " << _AMR.size()
//                   << ", numElements = " << nEl << "\n";
//         throw std::runtime_error("Mesh sanity check failed: AMR size must be 0 or numElements");
//       }
//
//       // 2) coordinate sizes
//       const std::size_t d = dim();
//       if (d > 0) {
//         const std::size_t nNodes0 = _X[0].size();
//         for (std::size_t k = 1; k < d; ++k) {
//           if (_X[k].size() != nNodes0) {
//             std::cerr << where << ": coordinate arrays size mismatch: "
//                       << "X[" << k << "].size() = " << _X[k].size()
//                       << ", X[0].size() = " << nNodes0 << "\n";
//             throw std::runtime_error("Mesh sanity check failed: coordinate arrays size mismatch");
//           }
//         }
//       }
//
//       // 3) dim vs element type (only checks first type; extend if needed)
//       if (!_elType.empty()) {
//         const unsigned t  = _elType[0];
//         const std::size_t d0 = dim();
//         if (t < 3) { // 3D: Hex27, Tet15, Wedge21
//           if (d0 != 3) {
//             std::cerr << where << ": element type " << t
//                       << " implies dim=3, but dim()=" << d0 << "\n";
//             throw std::runtime_error("Mesh sanity check failed: 3D element with non-3D coordinates");
//           }
//         }
//         else if (t < 5) {   // 2D: Quad9, Tri7
//           if (d0 != 2) {
//             std::cerr << where << ": element type " << t
//                       << " implies dim=2, but dim()=" << d0 << "\n";
//             throw std::runtime_error("Mesh sanity check failed: 2D element with non-2D coordinates");
//           }
//         }
//       }
//
//       if (mode == CheckMode::Basic) {
//         return; // skip heavy hierarchy checks
//       }
//
//       // ========== FULL / HEAVY CHECKS (hierarchy etc.) ==========
//
//       // 4) hierarchy container sizes: either 0 or nEl
//       if (!_father.empty() && _father.size() != nEl) {
//         std::cerr << where << ": father size mismatch: "
//                   << "father.size() = " << _father.size()
//                   << ", numElements = " << nEl << "\n";
//         throw std::runtime_error("Mesh sanity check failed: father size must be 0 or numElements");
//       }
//       if (!_children.empty() && _children.size() != nEl) {
//         std::cerr << where << ": children size mismatch: "
//                   << "children.size() = " << _children.size()
//                   << ", numElements = " << nEl << "\n";
//         throw std::runtime_error("Mesh sanity check failed: children size must be 0 or numElements");
//       }
//
//       // 5) Hierarchy invariants (local, not cross-level range checks)
//       if (_children.size() == nEl && !_children.empty()) {
//         const std::size_t d0 = dim();
//         const unsigned expectedRefChildren =
//           (d0 == 2 ? 4u :
//            d0 == 3 ? 8u : 0u);
//
//         for (std::size_t e = 0; e < nEl; ++e) {
//           const auto &ch = _children[e];
//           if (!ch.empty() && expectedRefChildren != 0u) {
//             const std::size_t k = ch.size();
//             // Option 2: 1 child (unrefined) or 4/8 children (refined)
//             if (k != 1u && k != expectedRefChildren) {
//               std::cerr << where << ": element " << e
//                         << " has " << k
//                         << " children; expected 1 or " << expectedRefChildren << "\n";
//               throw std::runtime_error("Mesh sanity check failed: wrong number of children");
//             }
//           }
//         }
//       }
//
//       if (_father.size() == nEl && !_father.empty()) {
//         for (std::size_t e = 0; e < nEl; ++e) {
//           const unsigned lvl = _elLevel[e];
//           if (lvl > 0 && _father[e] == noFather()) {
//             std::cerr << where << ": element " << e
//                       << " has level " << lvl
//                       << " but father[e] == noFather()\n";
//             throw std::runtime_error(
//               "Mesh sanity check failed: element with level>0 has no father");
//           }
//         }
//       }
//     }
// };
//
//
//
//
//
