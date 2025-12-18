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
    using AMRVector        = std::vector<int>;                 // -1,0,1 or 0/1, etc.
    using FatherVector     = std::vector<unsigned>;            // index on coarser level
    using ChildrenVector   = std::vector<std::vector<unsigned>>; // indices on finer level
    using NodeElemAdj      = std::vector<std::vector<unsigned>>; // node -> elements

    enum class CheckMode {
      Basic,  // cheap global checks
      Full    // Basic + heavier hierarchy checks
    };

    Mesh(ElemConnectivity &elTplgy,
         ElemType         &elType,
         ElemLevel        &elLevel,
         Coordinates      &X,
         AMRVector        &AMR,
         FatherVector     &father,
         ChildrenVector   &children,
         NodeElemAdj      &nodeElems)
      : _elTplgy(elTplgy),
        _elType(elType),
        _elLevel(elLevel),
        _X(X),
        _AMR(AMR),
        _father(father),
        _children(children),
        _nodeElems(nodeElems)
    {
      checkSanity("Mesh::Mesh", CheckMode::Full);
    }

    // ==================================
    // Assignment: copy DATA (not binding)
    // ==================================
    Mesh& operator=(const Mesh &other) {
      if (this != &other) {
        _elTplgy   = other._elTplgy;    // topology
        _elType    = other._elType;     // element types
        _elLevel   = other._elLevel;    // refinement levels
        _X         = other._X;          // coordinates
        _AMR       = other._AMR;        // AMR flags
        _father    = other._father;     // father indices (coarser level)
        _children  = other._children;   // children indices (finer level)
        _nodeElems = other._nodeElems;  // node->elements adjacency

        checkSanity("Mesh::operator=", CheckMode::Full);
      }
      return *this;
    }

    // Optional: friend swap (swaps underlying containers)
    friend void swap(Mesh &a, Mesh &b) noexcept {
      using std::swap;
      swap(a._elTplgy,   b._elTplgy);
      swap(a._elType,    b._elType);
      swap(a._elLevel,   b._elLevel);
      swap(a._X,         b._X);
      swap(a._AMR,       b._AMR);
      swap(a._father,    b._father);
      swap(a._children,  b._children);
      swap(a._nodeElems, b._nodeElems);
    }

    // ============================
    // Getters (non-const)
    // ============================
    ElemConnectivity & elTplgy()   { return _elTplgy;   }
    ElemType         & elType()    { return _elType;    }
    ElemLevel        & elLevel()   { return _elLevel;   }
    Coordinates      & X()         { return _X;         }
    AMRVector        & AMR()       { return _AMR;       }
    FatherVector     & father()    { return _father;    }
    ChildrenVector   & children()  { return _children;  }
    NodeElemAdj      & nodeElems() { return _nodeElems; }

    // ============================
    // Getters (const)
    // ============================
    const ElemConnectivity & elTplgy()   const { return _elTplgy;   }
    const ElemType         & elType()    const { return _elType;    }
    const ElemLevel        & elLevel()   const { return _elLevel;   }
    const Coordinates      & X()         const { return _X;         }
    const AMRVector        & AMR()       const { return _AMR;       }
    const FatherVector     & father()    const { return _father;    }
    const ChildrenVector   & children()  const { return _children;  }
    const NodeElemAdj      & nodeElems() const { return _nodeElems; }

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
            // In theory this shouldn’t happen under your assumptions
            // (you should always be able to raise Lmin by refining some AMR=0).
            throw std::runtime_error(
              "adjustAMRForOneLevelDiscontinuity: cannot fix Lmax-Lmin>1 at a node");
          }
        }

      } while (changedGlobal);
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
