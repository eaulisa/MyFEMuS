#pragma once

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cstddef>
#include <iostream>
#include <limits>
#include <string>
#include "Traits.hpp"
#include "PointLocatorResult.hpp"

class Mesh {
  public:
    using ElemConnectivity = std::vector<std::vector<unsigned>>;
    using ElemType         = std::vector<unsigned>;              // MUST stay unsigned
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
      _neighbors.clear();
    }

    // ==================================
    // Assignment: copy DATA (not binding)
    // ==================================
    Mesh& operator=(const Mesh &other) {
      if (this != &other) {
        _elTplgy    = other._elTplgy;
        _elType     = other._elType;
        _elLevel    = other._elLevel;
        _X          = other._X;
        _AMR        = other._AMR;
        _father     = other._father;
        _children   = other._children;
        _nodeElems  = other._nodeElems;
        _neighbors  = other._neighbors;

        checkSanity("Mesh::operator=", CheckMode::Full);
      }
      return *this;
    }

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
    void setUniformRefinement() {
      _AMR.clear();
    }

    void setRefineEvenElements() {
      const std::size_t nEl = numElements();
      _AMR.assign(nEl, 0);
      for (std::size_t i = 0; i < nEl; ++i) {
        if ((i % 2u) == 0u) _AMR[i] = 1;
      }
      checkSanity("Mesh::setRefineEvenElements", CheckMode::Basic);
    }

    void setRefinementFromBallLevelSetCrossing_OneRing(const std::vector<double>& center,
                                                       double r,
                                                       unsigned neighMode,   // 0,1,2,3
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

      auto vertex_range = [&](unsigned elTypeU) -> std::pair<std::size_t, std::size_t> {
        const ElemTraits &tr = getTraits(elTypeU);
        const std::size_t b = tr.vert.start;
        const std::size_t e = tr.vert.start + tr.vert.count;
        return std::pair<std::size_t, std::size_t>(b, e);
      };

      auto face_range = [&](unsigned elTypeU) -> std::pair<std::size_t, std::size_t> {
        const ElemTraits &tr = getTraits(elTypeU);

        // 3D: face-centers range
        if (tr.dim == 3) {
          const std::size_t b = tr.faceCenters.start;
          const std::size_t e = tr.faceCenters.start + tr.faceCenters.count;
          return std::pair<std::size_t, std::size_t>(b, e);
        }

        // 2D: "faces" are edges -> use edge mids range
        if (tr.dim == 2) {
          const std::size_t b = tr.edge.start;
          const std::size_t e = tr.edge.start + tr.edge.count;
          return std::pair<std::size_t, std::size_t>(b, e);
        }

        // 1D: "faces" are endpoints -> use vertices range (0..2)
        if (tr.dim == 1) {
          const std::size_t b = tr.vert.start;
          const std::size_t e = tr.vert.start + tr.vert.count;
          return std::pair<std::size_t, std::size_t>(b, e);
        }

        throw std::runtime_error("setRefinementFromBallLevelSetCrossing_OneRing: unsupported dim in traits");
      };

      auto pick_range = [&](unsigned elTypeU, std::size_t connSize, unsigned mode)
      -> std::pair<std::size_t, std::size_t> {

        if (mode == 0u) {
          return std::pair<std::size_t, std::size_t>(0, connSize);
        }
        if (mode == 1u) {
          return vertex_range(elTypeU);
        }
        if (mode == 2u) {
          return face_range(elTypeU);
        }

        // mode == 3 hybrid:
        // Hex27 & Quad9 -> faces; others -> vertices
        if (elTypeU == static_cast<unsigned>(Hex27) ||
            elTypeU == static_cast<unsigned>(Quad9)) {
          return face_range(elTypeU);
        }
        if (elTypeU == static_cast<unsigned>(Tet15)   ||
            elTypeU == static_cast<unsigned>(Wedge21) ||
            elTypeU == static_cast<unsigned>(Tri7)    ||
            elTypeU == static_cast<unsigned>(Line3)) {
          return vertex_range(elTypeU);
        }

        throw std::runtime_error("setRefinementFromBallLevelSetCrossing_OneRing: unsupported elType (hybrid)");
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
      // -------------------------
      for (std::size_t e = 0; e < nEl; ++e) {
        if (_AMR[e] != 2) continue;

        const auto& conn = _elTplgy[e];
        if (conn.empty()) continue;

        const unsigned et = _elType[e]; // unsigned as requested
        const std::pair<std::size_t, std::size_t> liRange = pick_range(et, conn.size(), neighMode);
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


    void setRefinementFromLocatedPoints_OneRing(
      const std::vector<PointLocatorResult>& results,
      unsigned neighMode   // 0,1,2,3
    ) {
      const std::size_t d    = dim();
      const std::size_t nEl  = numElements();
      const std::size_t nNod = numNodes();

      if (neighMode > 3u) {
        throw std::runtime_error(
          "setRefinementFromLocatedPoints_OneRing: neighMode must be 0,1,2,3");
      }

      if (_nodeElems.size() != nNod) {
        throw std::runtime_error(
          "setRefinementFromLocatedPoints_OneRing: _nodeElems is not built or wrong size");
      }

      _AMR.assign(nEl, 0);

      // --------------------------------------------------
      // Helpers (same spirit as your existing method)
      // --------------------------------------------------
      auto vertex_range = [&](unsigned elTypeU)
      -> std::pair<std::size_t, std::size_t> {
        const ElemTraits& tr = getTraits(elTypeU);
        return {tr.vert.start, tr.vert.start + tr.vert.count};
      };

      auto face_range = [&](unsigned elTypeU)
      -> std::pair<std::size_t, std::size_t> {
        const ElemTraits& tr = getTraits(elTypeU);

        if (tr.dim == 3) {
          return {
            tr.faceCenters.start,
            tr.faceCenters.start + tr.faceCenters.count};
        }
        if (tr.dim == 2) {
          return {tr.edge.start,
                  tr.edge.start + tr.edge.count};
        }
        if (tr.dim == 1) {
          return {tr.vert.start,
                  tr.vert.start + tr.vert.count};
        }

        throw std::runtime_error(
          "setRefinementFromLocatedPoints_OneRing: unsupported dim in traits");
      };

      auto pick_range =
        [&](unsigned elTypeU,
            std::size_t connSize,
            unsigned mode)
      -> std::pair<std::size_t, std::size_t> {

        if (mode == 0u) return {0, connSize};
        if (mode == 1u) return vertex_range(elTypeU);
        if (mode == 2u) return face_range(elTypeU);

        // hybrid (mode == 3)
        if (elTypeU == static_cast<unsigned>(Hex27) ||
            elTypeU == static_cast<unsigned>(Quad9)) {
          return face_range(elTypeU);
        }

        return vertex_range(elTypeU);
      };

      // --------------------------------------------------
      // 1) Seed elements from locator results
      // --------------------------------------------------
      for (const auto& r : results) {
        if (!r.ok) continue;  // skip not-found
        const std::size_t e = static_cast<std::size_t>(r.elem);
        if (e >= nEl) {
          throw std::runtime_error(
            "setRefinementFromLocatedPoints_OneRing: element index out of range");
        }
        _AMR[e] = 2;
      }

      // --------------------------------------------------
      // 2) One-ring expansion via node adjacency
      // --------------------------------------------------
      for (std::size_t e = 0; e < nEl; ++e) {
        if (_AMR[e] != 2) continue;

        const auto& conn = _elTplgy[e];
        if (conn.empty()) continue;

        const unsigned et = _elType[e];
        auto liRange = pick_range(et, conn.size(), neighMode);

        if (liRange.second > conn.size()) {
          throw std::runtime_error(
            "setRefinementFromLocatedPoints_OneRing: local range exceeds connectivity");
        }

        for (std::size_t li = liRange.first; li < liRange.second; ++li) {
          const unsigned gnode = conn[li];
          if (gnode >= nNod) {
            throw std::runtime_error(
              "setRefinementFromLocatedPoints_OneRing: node index out of range");
          }

          for (unsigned ne_u : _nodeElems[gnode]) {
            const std::size_t ne = static_cast<std::size_t>(ne_u);
            if (ne >= nEl) {
              throw std::runtime_error(
                "setRefinementFromLocatedPoints_OneRing: _nodeElems out of range");
            }
            if (_AMR[ne] == 0) _AMR[ne] = 1;
          }
        }
      }

      // --------------------------------------------------
      // 3) Normalize: everything nonzero -> 1
      // --------------------------------------------------
      for (std::size_t e = 0; e < nEl; ++e) {
        if (_AMR[e] != 0) _AMR[e] = 1;
      }

      checkSanity(
        "Mesh::setRefinementFromLocatedPoints_OneRing",
        CheckMode::Basic);
    }



    bool needsRefinement(std::size_t i) const {
      if (_AMR.empty()) return true;
      if (i >= _AMR.size()) throw std::out_of_range("Mesh::needsRefinement: index out of range");
      return (_AMR[i] == 1);
    }

    // ============================
    // Hierarchy interface (cross-level)
    // ============================
    static unsigned noFather() {
      return std::numeric_limits<unsigned>::max();
    }

    void setChildren(std::size_t i, const std::vector<unsigned> &ch) {
      if (i >= numElements()) throw std::out_of_range("Mesh::setChildren: index out of range");
      if (_children.size() != numElements()) _children.resize(numElements());
      _children[i] = ch;
    }

    const std::vector<unsigned>& childrenOf(std::size_t i) const {
      if (i >= numElements()) throw std::out_of_range("Mesh::childrenOf: index out of range");
      return _children[i];
    }

    std::vector<unsigned>& childrenOf(std::size_t i) {
      if (i >= numElements()) throw std::out_of_range("Mesh::childrenOf: index out of range");
      return _children[i];
    }

    unsigned fatherOf(std::size_t i) const {
      if (i >= numElements()) throw std::out_of_range("Mesh::fatherOf: index out of range");
      return _father[i];
    }

    void clearAllChildren() {
      _children.clear();
    }

    void resetAllFathersToNoFather() {
      const std::size_t nEl = numElements();
      const unsigned no_father = noFather();
      _father.assign(nEl, no_father);
    }

    void checkAll(const char* where = "Mesh::checkAll") const {
      checkSanity(where, CheckMode::Full);
    }

    // ============================
    // Node -> elements adjacency
    // ============================
    void buildNodeToElementAdjacency() {
      const std::size_t nElems  = numElements();
      const std::size_t nNodes  = numNodes();


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


      _nodeElems.resize(nNodes);           // keep existing inner vectors & capacities
      for (auto& v : _nodeElems) v.clear(); // clear sizes but keep capacity

      for (std::size_t n = 0; n < nNodes; ++n) {
        const auto need = (std::size_t)counts[n];
        if (_nodeElems[n].capacity() < need) {
          _nodeElems[n].reserve(need);
        }
      }


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
    // Build face neighbors via face-center nodes (UNIFORM MESH ASSUMPTION)
    // ============================================================
    //
    // Uses traits:
    //   3D: faceCenters range
    //   2D: edge range (edge mids)
    //   1D: vert range (endpoints)
    //
    // For such a "side-center" node n:
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

        const unsigned et = _elType[e];
        const ElemTraits &tr = getTraits(et);

        const unsigned nf = tr.nSides;

        unsigned i0 = 0u;
        if (tr.dim == 3) {
          i0 = tr.faceCenters.start;
        }
        else if (tr.dim == 2) {
          i0 = tr.edge.start;
        }
        else if (tr.dim == 1) {
          i0 = tr.vert.start; // endpoints
        }
        else {
          throw std::runtime_error("Mesh::buildFaceNeighborsFromNodeToElement: unsupported dim in traits");
        }

        if (static_cast<std::size_t>(i0 + nf) > conn.size()) {
          throw std::runtime_error("Mesh::buildFaceNeighborsFromNodeToElement: conn too small for side-center range");
        }

        _neighbors[e].assign(nf, UMAX);

        for (unsigned f = 0; f < nf; ++f) {
          const unsigned node = conn[i0 + f];
          if (node >= nNod) {
            throw std::runtime_error("Mesh::buildFaceNeighborsFromNodeToElement: side-center node out of range");
          }

          const auto& inc = _nodeElems[node];

          if (inc.size() == 1u) {
            if (inc[0] != static_cast<unsigned>(e)) {
              throw std::runtime_error("Mesh::buildFaceNeighborsFromNodeToElement: boundary side-center not owned by e");
            }
            _neighbors[e][f] = UMAX;
          }
          else if (inc.size() == 2u) {
            const unsigned a0 = inc[0], a1 = inc[1];
            if (a0 == static_cast<unsigned>(e)) _neighbors[e][f] = a1;
            else if (a1 == static_cast<unsigned>(e)) _neighbors[e][f] = a0;
            else {
              throw std::runtime_error("Mesh::buildFaceNeighborsFromNodeToElement: side-center does not reference e");
            }
          }
          else {
            throw std::runtime_error("Mesh::buildFaceNeighborsFromNodeToElement: side-center shared by >2 elements");
          }
        }
      }

      checkSanity("Mesh::buildFaceNeighborsFromNodeToElement", CheckMode::Basic);
    }

    void adjustAMRForOneLevelDiscontinuity() {
      const std::size_t nElems  = numElements();
      const std::size_t nNodes  = numNodes();

      if (_AMR.empty()) {
        return;
      }

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

        for (std::size_t n = 0; n < nNodes; ++n) {
          auto &els = nodeElems[n];
          const std::size_t m = els.size();
          if (m < 2) continue;

          unsigned Lmin = std::numeric_limits<unsigned>::max();
          unsigned Lmax = 0u;

          for (std::size_t k = 0; k < m; ++k) {
            const unsigned e = els[k];
            const unsigned L = _elLevel[e] + static_cast<unsigned>(_AMR[e]);
            if (L < Lmin) Lmin = L;
            if (L > Lmax) Lmax = L;
          }

          if (Lmax <= Lmin + 1u) {
            continue;
          }
          else if (Lmax > Lmin + 2u) {
            throw std::runtime_error("adjustAMRForOneLevelDiscontinuity: "
                                     "Lmax-Lmin > 2 at a node (inconsistent input)");
          }

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



    // Push-back for reference triangle: r>=0, s>=0, r+s<=1
    static inline void pushback_tri(double &r, double &s, const double eps) {
      if (r < 0.0 && r > -8.0 * eps) r = 0.0;
      if (s < 0.0 && s > -8.0 * eps) s = 0.0;

      const double rs = r + s;
      if (rs > 1.0 && rs < 1.0 + 8.0 * eps) {
        // scale back to the boundary r+s=1
        r /= rs;
        s /= rs;
      }
    }

// Push-back for line z in [-1,1]
    static inline void pushback_line(double &z, const double eps) {
      if (z < -1.0 && z > -1.0 - 8.0 * eps) z = -1.0;
      if (z >  1.0 && z <  1.0 + 8.0 * eps) z =  1.0;
    }

// Push-back for reference tet: x>=0,y>=0,z>=0,x+y+z<=1
    static inline void pushback_tet(double &x, double &y, double &z, const double eps) {
      if (x < 0.0 && x > -8.0 * eps) x = 0.0;
      if (y < 0.0 && y > -8.0 * eps) y = 0.0;
      if (z < 0.0 && z > -8.0 * eps) z = 0.0;

      const double s = x + y + z;
      if (s > 1.0 && s < 1.0 + 8.0 * eps) {
        // scale back to the boundary x+y+z=1
        x /= s;
        y /= s;
        z /= s;
      }
    }



// ------------------------------------------------------------
// Project PointLocatorResult from this mesh to mesh1
// ------------------------------------------------------------
    void projectPointLocatorResultsToNextLevel(
      const Mesh& mesh1,
      const std::vector<PointLocatorResult>& in0,
      std::vector<PointLocatorResult>& out1,
      double eps = 1e-12) const {
      const std::size_t nEl0 = numElements();
      const std::size_t nEl1 = mesh1.numElements();

      auto inside_tet = [&](double u, double v, double w) -> bool {
        return (u >= -eps) && (v >= -eps) && (w >= -eps) && (u + v + w <= 1.0 + eps);
      };

      auto inside_tri = [&](double r, double s) -> bool {
        return (r >= -eps) && (s >= -eps) && (r + s <= 1.0 + eps);
      };

      auto inside_line = [&](double t) -> bool {
        return (t >= -1.0 - eps) && (t <= 1.0 + eps);
      };

      // small helpers: avoid xi = {..} (which tends to show up as vector::operator=)
      auto set1 = [&](std::vector<double>& xi, double a) {
        //xi.resize(1);
        xi[0] = a;
      };
      auto set2 = [&](std::vector<double>& xi, double a, double b) {
        //xi.resize(2);
        xi[0] = a; xi[1] = b;
      };
      auto set3 = [&](std::vector<double>& xi, double a, double b, double c) {
        //xi.resize(3);
        xi[0] = a; xi[1] = b; xi[2] = c;
      };

      auto pick_child_by_contains = [&](unsigned et,
                                        const std::vector<double>& xi0,
                                        unsigned & childIndex,
      std::vector<double>& xi1) {
        // Sets childIndex and xi1 (child reference coords)
        // Child ordering MUST match your verts tables.

        if (et == static_cast<unsigned>(Line3)) {
          if (xi0.size() != 1) throw std::runtime_error("project: Line3 expects xi size 1");
          const double x = xi0[0];

          // order: c0 = left (x<=0), c1 = right (x>0); tie goes to c0
          if (x <= 0.0) {
            childIndex = 0;
            set1(xi1, 2.0 * x + 1.0);
          }
          else          {
            childIndex = 1;
            set1(xi1, 2.0 * x - 1.0);
          }

          if (!inside_line(xi1[0])) throw std::runtime_error("project: Line3 mapped xi out of [-1,1]");
          return;
        }

        if (et == static_cast<unsigned>(Quad9)) {
          if (xi0.size() != 2) throw std::runtime_error("project: Quad9 expects xi size 2");
          const double x = xi0[0];
          const double y = xi0[1];

          const bool right = (x > 0.0);
          const bool top   = (y > 0.0);

          // c0 BL, c1 BR, c2 TR, c3 TL (top row swapped)
          if (!top) childIndex = right ? 1u : 0u;
          else      childIndex = right ? 2u : 3u;

          const double x1 = right ? (2.0 * x - 1.0) : (2.0 * x + 1.0);
          const double y1 = top   ? (2.0 * y - 1.0) : (2.0 * y + 1.0);
          set2(xi1, x1, y1);

          if (!inside_line(xi1[0]) || !inside_line(xi1[1])) {
            throw std::runtime_error("project: Quad9 mapped xi out of [-1,1]^2");
          }
          return;
        }

        if (et == static_cast<unsigned>(Hex27)) {
          if (xi0.size() != 3) throw std::runtime_error("project: Hex27 expects xi size 3");
          const double x = xi0[0];
          const double y = xi0[1];
          const double z = xi0[2];

          const bool right = (x > 0.0);
          const bool top   = (y > 0.0);
          const bool up    = (z > 0.0);

          unsigned c2d;
          // same swapped top row as Quad9
          if (!top) c2d = right ? 1u : 0u;
          else      c2d = right ? 2u : 3u;

          childIndex = c2d + (up ? 4u : 0u);

          const double x1 = right ? (2.0 * x - 1.0) : (2.0 * x + 1.0);
          const double y1 = top   ? (2.0 * y - 1.0) : (2.0 * y + 1.0);
          const double z1 = up    ? (2.0 * z - 1.0) : (2.0 * z + 1.0);
          set3(xi1, x1, y1, z1);

          if (!inside_line(xi1[0]) || !inside_line(xi1[1]) || !inside_line(xi1[2])) {
            throw std::runtime_error("project: Hex27 mapped xi out of [-1,1]^3");
          }
          return;
        }


        if (et == static_cast<unsigned>(Tri7)) {
          // Reference triangle: (r,s) with r>=0, s>=0, r+s<=1
          if (xi0.size() != 2) throw std::runtime_error("project: Tri7 expects xi size 2");
          const double r  = xi0[0];
          const double s  = xi0[1];
          const double rs = r + s;

          if (!inside_tri(r, s)) {
            throw std::runtime_error("project: Tri7 input xi not inside reference triangle");
          }

          // Precompute mapped coordinates (cheap, avoids repeated ops)
          const double r2  = 2.0 * r;
          const double s2  = 2.0 * s;
          const double r2m = r2 - 1.0;   // 2r - 1
          const double s2m = s2 - 1.0;   // 2s - 1
          const double r1  = 1.0 - r2;   // 1 - 2r
          const double s1  = 1.0 - s2;   // 1 - 2s

          // Robust, consistent classification near split lines r=0.5, s=0.5, r+s=0.5
          const double a = r  - 0.5;
          const double b = s  - 0.5;
          const double c = rs - 0.5;

          // Child order: c0 near V0, c1 near V1, c2 near V2, c3 central (inverted)
          if (a <= eps && b <= eps) {
            if (c <= eps) {
              childIndex = 0;
              set2(xi1, r2, s2);
            }
            else {
              childIndex = 3;
              set2(xi1, r1, s1);
            }
          }
          else if (a > eps) {
            childIndex = 1;
            set2(xi1, r2m, s2);
          }
          else {
            childIndex = 2;
            set2(xi1, r2, s2m);
          }

          // Optional: snap tiny negatives/overshoot due to roundoff (always on)
          if (xi1[0] < 0.0 && xi1[0] > -8.0 * eps) xi1[0] = 0.0;
          if (xi1[1] < 0.0 && xi1[1] > -8.0 * eps) xi1[1] = 0.0;
          const double sum = xi1[0] + xi1[1];
          if (sum > 1.0 && sum < 1.0 + 8.0 * eps) {
            xi1[0] /= sum;
            xi1[1] /= sum;
          }

#ifndef NDEBUG
          if (!inside_tri(xi1[0], xi1[1])) {
            std::cout << r << " " << s << std::endl;
            std::cout << xi1[0] << " " << xi1[1] << std::endl;
            throw std::runtime_error("project: Tri7 mapped xi not inside child reference triangle");
          }
#endif
          return;
        }


        if (et == static_cast<unsigned>(Wedge21)) {
          // Wedge: (r,s) triangle, z in [-1,1]
          if (xi0.size() != 3) throw std::runtime_error("project: Wedge21 expects xi size 3");
          const double r = xi0[0];
          const double s = xi0[1];
          const double z = xi0[2];

          if (!inside_tri(r, s) || !inside_line(z)) {
            throw std::runtime_error("project: Wedge21 input xi not inside reference wedge");
          }

          const bool up = (z > 0.0);
          const unsigned cz = up ? 1u : 0u;
          const double z1 = up ? (2.0 * z - 1.0) : (2.0 * z + 1.0);

          unsigned ctri = 0;
          double r1 = 0.0, s1 = 0.0;

          if (r <= 0.5 + eps && s <= 0.5 + eps) {
            if (r + s <= 0.5 + eps) {
              ctri = 0; r1 = 2.0 * r;       s1 = 2.0 * s;
            }
            else {
              ctri = 3; r1 = 1.0 - 2.0 * r; s1 = 1.0 - 2.0 * s;
            }
          }
          else if (r > 0.5) {
            ctri = 1; r1 = 2.0 * r - 1.0;   s1 = 2.0 * s;
          }
          else {
            ctri = 2; r1 = 2.0 * r;         s1 = 2.0 * s - 1.0;
          }

          childIndex = ctri + 4u * cz;


          set3(xi1, r1, s1, z1);

          // Push back to boundaries if tiny overshoot due to roundoff
          pushback_tri(xi1[0], xi1[1], eps);
          pushback_line(xi1[2], eps);

#ifndef NDEBUG
          if (!inside_tri(xi1[0], xi1[1]) || !inside_line(xi1[2])) {
            std::cout << "Wedge21 input  : " << r << " " << s << " " << z << "\n";
            std::cout << "Wedge21 mapped : " << xi1[0] << " " << xi1[1] << " " << xi1[2] << "\n";
            throw std::runtime_error("project: Wedge21 mapped xi not inside child reference wedge");
          }
#endif

          return;



          if (!inside_tri(xi1[0], xi1[1]) || !inside_line(xi1[2])) {
            throw std::runtime_error("project: Wedge21 mapped xi not inside child reference wedge");
          }
          return;
        }

        if (et == static_cast<unsigned>(Tet15)) {
          // Tet: (x,y,z) with x>=0,y>=0,z>=0,x+y+z<=1
          if (xi0.size() != 3) throw std::runtime_error("project: Tet15 expects xi size 3");
          const double x = xi0[0];
          const double y = xi0[1];
          const double z = xi0[2];

          if (!inside_tet(x, y, z)) {
            throw std::runtime_error("project: Tet15 input xi not inside reference tet");
          }

          auto map_child = [&](unsigned c, double & u, double & v, double & w) {
            switch (c) {
            case 0: u = 2.0 * x;          v = 2.0 * y;          w = 2.0 * z;          return;
            case 1: u = 2.0 * x - 1.0;    v = 2.0 * y;          w = 2.0 * z;          return;
            case 2: u = 2.0 * x;          v = 2.0 * y - 1.0;    w = 2.0 * z;          return;
            case 3: u = 2.0 * x;          v = 2.0 * y;          w = 2.0 * z - 1.0;    return;

            case 4: u = 1.0 - 2.0 * x - 2.0 * z;
              v = 1.0 - 2.0 * y - 2.0 * z;
              w = 2.0 * z;
              return;

            case 5: u = 1.0 - 2.0 * x;
              v = 2.0 * y;
              w = 1.0 - 2.0 * y - 2.0 * z;
              return;

            case 6: u = 2.0 * y + 2.0 * z - 1.0;
              v = 2.0 * x + 2.0 * z - 1.0;
              w = 1.0 - 2.0 * z;
              return;

            case 7: u = 2.0 * x;
              v = 1.0 - 2.0 * y;
              w = 1.0 - 2.0 * x - 2.0 * z;
              return;

            default: throw std::runtime_error("project: Tet15 child index out of range");
            }
          };


          for (unsigned c = 0; c < 8u; ++c) {
            double u, v, w;
            map_child(c, u, v, w);

            // Push back to boundaries if tiny overshoot due to roundoff
            pushback_tet(u, v, w, eps);

            if (inside_tet(u, v, w)) {
              childIndex = c;
              set3(xi1, u, v, w);
#ifndef NDEBUG
              if (!inside_tet(xi1[0], xi1[1], xi1[2])) {
                std::cout << "Tet15 input  : " << x << " " << y << " " << z << "\n";
                std::cout << "Tet15 mapped : " << xi1[0] << " " << xi1[1] << " " << xi1[2] << "\n";
                throw std::runtime_error("project: Tet15 mapped xi not inside child reference tet");
              }
#endif
              return;
            }
          }

          throw std::runtime_error("project: Tet15 could not select a child (no candidate contained)");
        }
        throw std::runtime_error("project: unsupported element type for projection");
      };

      // ------------------------------------------------------------
      // Main loop: map each PointLocatorResult
      // ------------------------------------------------------------
      out1.resize(in0.size());

      const unsigned dim = this->dim();

      for (std::size_t i = 0; i < in0.size(); ++i) {
        const PointLocatorResult& r0 = in0[i];
        PointLocatorResult&       r1 = out1[i];

        r1.ok   = r0.ok;
        r1.elem = r0.elem;

        if (!r0.ok || r0.elem == UMAX) {
          r1.xi.clear();          // or resize(0)
          continue;
        }

#ifndef NDEBUG
        if (r0.xi.size() < dim) {
          throw std::runtime_error("PointLocatorResult: r0.xi too small");
        }
#endif

        r1.xi.resize(dim);        // no alloc (reserve(3) already)
        std::copy_n(r0.xi.data(), dim, r1.xi.data());



        const std::size_t e0 = static_cast<std::size_t>(r0.elem);
        if (e0 >= nEl0) {
          throw std::runtime_error("projectPointLocator: input elem out of range on mesh0");
        }

        const auto& ch = _children[e0];
        if (ch.empty()) {
          throw std::runtime_error("projectPointLocator: _children[e0] is empty");
        }

        const unsigned et = _elType[e0];

        std::size_t nChild = 0;
        switch (static_cast<ElType>(et)) {
        case Line3:   nChild = 2; break;
        case Quad9:   nChild = 4; break;
        case Tri7:    nChild = 4; break;
        case Hex27:   nChild = 8; break;
        case Wedge21: nChild = 8; break;
        case Tet15:   nChild = 8; break;
        default:
          throw std::runtime_error("projectPointLocator: unsupported element type");
        }

        if (ch.size() == 1) {
          const std::size_t e1 = static_cast<std::size_t>(ch[0]);
          if (e1 >= nEl1) throw std::runtime_error("projectPointLocator: child elem out of range on mesh1");

          r1.elem = static_cast<unsigned>(e1);
          // r1.xi already equals r0.xi due to r1=r0 above (identity mapping)
          r1.ok   = true;
          continue;
        }

        if (ch.size() != nChild) {
          throw std::runtime_error("projectPointLocator: children size does not match element refinement arity");
        }

        // Refined: pick child and map xi (WRITE DIRECTLY INTO r1.xi)
        unsigned childIndex = 0;
        pick_child_by_contains(et, r0.xi, childIndex, r1.xi);

        const std::size_t e1 = static_cast<std::size_t>(ch[childIndex]);
        if (e1 >= nEl1) throw std::runtime_error("projectPointLocator: child elem out of range on mesh1");

        r1.elem = static_cast<unsigned>(e1);
        r1.ok   = true;
      }
    }


  private:
    ElemConnectivity &_elTplgy;
    ElemType         &_elType;
    ElemLevel        &_elLevel;
    Coordinates      &_X;
    AMRVector        &_AMR;
    FatherVector     &_father;
    ChildrenVector   &_children;
    NodeElemAdj      &_nodeElems;
    ElemNeighbors    &_neighbors;

    // ============================
    // Sanity checks
    // ============================
    void checkSanity(const char *where,
                     CheckMode mode = CheckMode::Basic) const {
      const std::size_t nEl = _elType.size();

      // BASIC CHECKS

      if (_elLevel.size() != nEl || _elTplgy.size() != nEl) {
        std::cerr << where << ": element arrays size mismatch: "
                  << "elType="    << _elType.size()
                  << ", elLevel=" << _elLevel.size()
                  << ", elTplgy=" << _elTplgy.size() << "\n";
        throw std::runtime_error("Mesh sanity check failed: element arrays size mismatch");
      }

      if (!_AMR.empty() && _AMR.size() != nEl) {
        std::cerr << where << ": AMR size invalid: "
                  << "AMR.size() = " << _AMR.size()
                  << ", numElements = " << nEl << "\n";
        throw std::runtime_error("Mesh sanity check failed: AMR size must be 0 or numElements");
      }

      if (!_neighbors.empty() && _neighbors.size() != nEl) {
        std::cerr << where << ": neighbors size invalid: "
                  << "neighbors.size() = " << _neighbors.size()
                  << ", numElements = " << nEl << "\n";
        throw std::runtime_error("Mesh sanity check failed: neighbors size must be 0 or numElements");
      }

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

      // dim vs element type (traits-driven)
      if (!_elType.empty()) {
        const unsigned t0 = _elType[0];
        const ElemTraits &tr0 = getTraits(t0);
        if (dim() != static_cast<std::size_t>(tr0.dim)) {
          std::cerr << where << ": element type " << t0
                    << " implies dim=" << tr0.dim
                    << ", but dim()=" << dim() << "\n";
          throw std::runtime_error("Mesh sanity check failed: dim mismatch with traits");
        }
      }

      if (mode == CheckMode::Basic) {
        return;
      }

      // FULL / HEAVY CHECKS

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

      // neighbors per-element size check (if built) -> traits-driven
      if (_neighbors.size() == nEl && !_neighbors.empty()) {
        for (std::size_t e = 0; e < nEl; ++e) {
          const auto& conn = _elTplgy[e];
          if (conn.empty()) continue;

          const unsigned et = _elType[e];
          const ElemTraits &tr = getTraits(et);
          const unsigned nf = tr.nSides;

          if (_neighbors[e].size() != nf) {
            std::cerr << where << ": neighbors[" << e << "].size()=" << _neighbors[e].size()
                      << " but expected " << nf << " for elType=" << et << "\n";
            throw std::runtime_error("Mesh sanity check failed: wrong neighbors[e] size");
          }
        }
      }

      // hierarchy checks (unchanged)
      if (_children.size() == nEl && !_children.empty()) {
        const std::size_t d0 = dim();
        const unsigned expectedRefChildren =
          (d0 == 2 ? 4u :
           d0 == 3 ? 8u : 0u);

        for (std::size_t e = 0; e < nEl; ++e) {
          const auto &ch = _children[e];
          if (!ch.empty() && expectedRefChildren != 0u) {
            const std::size_t k = ch.size();
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




