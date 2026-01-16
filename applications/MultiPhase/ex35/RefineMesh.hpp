#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
#include <unordered_map>

#include "Mesh.hpp"
#include "FemProjection.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif


double computeHSnapElementBased(const Mesh& mesh, double rel = 1e-2) {

  const std::vector<std::vector<double>>& X = mesh.X();
  const std::vector<unsigned>& elType = mesh.elType();
  const std::vector<std::vector<unsigned>>& elTplgy = mesh.elTplgy();
  const std::size_t dim = mesh.dim();
  const std::size_t numElem = mesh.numElements();

  static const unsigned cornerCount[6] = {8, 4, 6, 4, 3, 2};//Hex, Tet, Wedge, Quad, Tri, Line number of Corners

  std::vector<double> xCorner;
  xCorner.reserve(8);

  double hMin = std::numeric_limits<double>::max();

  for (std::size_t e = 0; e < numElem; ++e) {
    const unsigned type    = elType[e];
    const unsigned nCorner = cornerCount[type];
    const std::vector<unsigned>& conn = elTplgy[e];

    xCorner.resize(nCorner);

    // element characteristic size = min range across used dimensions
    double h_e = std::numeric_limits<double>::max();

    for (unsigned d = 0; d < dim; ++d) {
      for (unsigned j = 0; j < nCorner; ++j) {
        xCorner[j] = X[d][conn[j]];
      }

      std::pair<std::vector<double>::iterator,
          std::vector<double>::iterator> minmax_it =
            std::minmax_element(xCorner.begin(), xCorner.end());

      double range_d = *minmax_it.second - *minmax_it.first;
      if (range_d < h_e) h_e = range_d;
    }
    if (h_e < hMin) hMin = h_e;
  }

  if (hMin <= 0.0) { //degenerate mesh
    // fallback: absolute snap length
    hMin = 1e-12;
  }

  return hMin * rel;
}

// hash-based dedup using snapped coordinates as key
void dedupNodesHash(Mesh &mesh,
                    std::vector<unsigned>& candidateIndices
                   ) {

  std::vector<std::vector<double>>& X = mesh.X();
  std::vector<std::vector<unsigned>>& elTplgy = mesh.elTplgy();

  unsigned dim = mesh.dim();
  if (dim != 2 && dim != 3) {
    throw std::runtime_error("dedupNodesHash: dim must be 2 or 3");
  }

  const unsigned N = static_cast<unsigned>(X[0].size());
  if (N == 0) {
    return;
  }
  for (unsigned d = 1; d < dim; ++d) {
    if (X[d].size() != X[0].size()) {
      throw std::runtime_error("dedupNodesHash: XF[d] size mismatch among components");
    }
  }

  // 1) unique candidateIndices
  std::sort(candidateIndices.begin(), candidateIndices.end());
  candidateIndices.erase(
    std::unique(candidateIndices.begin(), candidateIndices.end()),
    candidateIndices.end()
  );
  const unsigned K = static_cast<unsigned>(candidateIndices.size());
  if (K == 0) {
    return;
  }

  // 2) snapping length from element geometry
  double h_snap = computeHSnapElementBased(mesh, 1.e-2);
  if (h_snap <= 0.0) {
    // degenerate case fallback: absolute tiny length
    h_snap = 1e-12;
  }

  // optional: comment out for fair timing
  //std::cout << "characteristic size (hash) = " << h_snap << std::endl;

  // snap[d][i] = snapped integer coord of node i in dimension d
  // std::vector<std::vector<long long>> snap(dim, std::vector<long long>(N));
  // for (unsigned d = 0; d < dim; ++d) {
  //   const double inv_h = 1.0 / h_snap;
  //   for (unsigned i = 0; i < N; ++i) {
  //     snap[d][i] = static_cast<long long>(std::llround(X[d][i] * inv_h));
  //   }
  // }


  const double inv_h = 1.0 / h_snap;
  // Key for hash map
  struct Key {
    long long s0, s1, s2;
  };

  struct KeyHash {
    std::size_t operator()(const Key& k) const noexcept {
      std::size_t h = 1469598103934665603ull; // FNV-ish mix
      auto mix = [&](long long v) {
        std::size_t x = static_cast<std::size_t>(v);
        h ^= x + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
      };
      mix(k.s0);
      mix(k.s1);
      mix(k.s2);
      return h;
    }
  };

  struct KeyEq {
    bool operator()(const Key& a, const Key& b) const noexcept {
      return (a.s0 == b.s0) && (a.s1 == b.s1) && (a.s2 == b.s2);
    }
  };

  // 3) Build representative mapping rep[i] with hash buckets
  std::vector<unsigned> rep(N);
  std::iota(rep.begin(), rep.end(), 0);

  std::unordered_map<Key, unsigned, KeyHash, KeyEq> cellRep;
  cellRep.reserve(K * 2); // a bit of slack

  for (unsigned idx : candidateIndices) {
    Key key;
    key.s0 = static_cast<long long>(std::llround(X[0][idx] * inv_h)); //snap[0][idx];
    key.s1 = (dim > 1 ? static_cast<long long>(std::llround(X[1][idx] * inv_h))/*snap[1][idx]*/ : 0);
    key.s2 = (dim > 2 ? static_cast<long long>(std::llround(X[2][idx] * inv_h))/*;snap[2][idx]*/ : 0);

    auto it = cellRep.find(key);
    if (it == cellRep.end()) {
      // first node in this snapped cell becomes representative
      cellRep.emplace(key, idx);
      // rep[idx] already == idx from iota
    }
    else {
      // merge into existing representative
      rep[idx] = it->second;
    }
  }

  // 4) Assign new indices 0..numUnique-1
  std::vector<unsigned> old2new(N);
  unsigned numUnique = 0;

  for (unsigned i = 0; i < N; ++i) {
    if (rep[i] == i) {
      old2new[i] = numUnique++;
    }
  }


  #pragma omp parallel for schedule(static)
  for (long long i = 0; i < (long long)N; ++i) {
    unsigned ui = (unsigned)i;
    if (rep[ui] != ui) {
      old2new[ui] = old2new[ rep[ui] ];
    }
  }

  // for (unsigned i = 0; i < N; ++i) {
  //   if (rep[i] != i) {
  //     old2new[i] = old2new[rep[i]];
  //   }
  // }

  if (numUnique == N) {
    return;
  }



  const unsigned* map = old2new.data();

  #pragma omp parallel for schedule(static)
  for (long long e = 0; e < (long long)elTplgy.size(); ++e) {
    auto& conn = elTplgy[(size_t)e];
    for (unsigned& v : conn) {
      v = map[v];
    }
  }




  // // 5) Update connectivity
  // for (auto& conn : elTplgy) {
  //   for (auto& v : conn) {
  //     v = old2new[v];
  //   }
  // }


  // 6) Compact X to size numUnique (faster + OpenMP)
  std::vector<unsigned> reps;
  reps.reserve(numUnique);
  for (unsigned i = 0; i < N; ++i) {
    if (rep[i] == i) reps.push_back(i);
  }

  for (unsigned d = 0; d < dim; ++d) {
    std::vector<double> Xnew(numUnique);
    const double* Xd = X[d].data();

    #pragma omp parallel for schedule(static)
    for (long long j = 0; j < (long long)numUnique; ++j) {
      unsigned i = reps[(size_t)j];
      Xnew[(size_t)j] = Xd[i];
    }
    X[d].swap(Xnew);
  }



  /*
    // 6) Compact XF to size numUnique
    for (unsigned d = 0; d < dim; ++d) {
      std::vector<double> Xnew(numUnique);
      for (unsigned i = 0; i < N; ++i) {
        if (rep[i] == i) {
          unsigned j = old2new[i];
          Xnew[j] = X[d][i];
        }
      }
      X[d].swap(Xnew);
    }*/
}


void refineAndProjectMesh(
  const std::array<std::unique_ptr<FemProjection>, 6>& elProj,
  Mesh &meshC,
  Mesh &meshF) {

  const auto& XC       = meshC.X();
  const auto& elTypeC  = meshC.elType();
  const auto& elTplgyC = meshC.elTplgy();
  const auto& elLevelC = meshC.elLevel();

  auto& XF       = meshF.X();
  auto& elTypeF  = meshF.elType();
  auto& elTplgyF = meshF.elTplgy();
  auto& elLevelF = meshF.elLevel();
  auto& fatherF  = meshF.father();

  const std::size_t numElemC       = elTplgyC.size();
  const unsigned    dim            = static_cast<unsigned>(meshC.dim());
  const std::size_t numCoarseNodes = meshC.numNodes();

  if (dim != 2 && dim != 3) {
    throw std::runtime_error("refineAndProjectMesh: dim must be 2 or 3");
  }

  if (numCoarseNodes == 0) {
    throw std::runtime_error("refineAndProjectMesh: XC[0] is empty");
  }

  for (unsigned d = 1; d < dim; ++d) {
    if (XC[d].size() != numCoarseNodes) {
      throw std::runtime_error("refineAndProjectMesh: XC[d] size mismatch among components");
    }
  }

  const unsigned maxCoarseNodesPerElem = (dim == 2) ? 9u  : 27u;   // Quad9 vs Hex27
  const unsigned maxNewNodesPerElem    = (dim == 2) ? 16u : 98u;
  const unsigned maxChildrenPerElem    = (dim == 2) ? 4u  : 8u;

  // reset hierarchy info
  meshC.clearAllChildren();
  meshF.clearAllChildren();

  // reset fine mesh containers
  elTypeF.clear();
  elTplgyF.clear();
  elLevelF.clear();
  fatherF.clear();

  elTypeF.reserve(maxChildrenPerElem * numElemC);
  elTplgyF.reserve(maxChildrenPerElem * numElemC);
  elLevelF.reserve(maxChildrenPerElem * numElemC);
  fatherF.reserve(maxChildrenPerElem * numElemC);

  XF.resize(dim);
  for (unsigned d = 0; d < dim; ++d) {
    XF[d].clear();
    XF[d].reserve(numCoarseNodes + maxNewNodesPerElem * numElemC);
    XF[d].insert(XF[d].end(), XC[d].begin(), XC[d].end());
  }

  // local coarse/fine buffers
  std::vector<std::vector<double>> U0(dim), U1(dim);
  for (unsigned d = 0; d < dim; ++d) {
    U0[d].reserve(maxCoarseNodesPerElem);
    U1[d].reserve(maxNewNodesPerElem);
  }

  // candidate boundary nodes for dedup
  std::vector<unsigned> dupLocIndicesF;
  dupLocIndicesF.reserve(numElemC * maxNewNodesPerElem);

  // temporary per-element refined connectivity
  std::vector<std::vector<unsigned>> elTplgy;
  elTplgy.reserve(maxChildrenPerElem);

  std::vector<unsigned> children;
  children.reserve(maxChildrenPerElem);

  // 'n' = next free node index, passed by reference to refineConnectivity
  unsigned n = static_cast<unsigned>(numCoarseNodes);

  for (std::size_t i = 0; i < numElemC; ++i) {

    const unsigned type   = elTypeC[i];
    const auto&    connC  = elTplgyC[i];
    const unsigned levelC = elLevelC[i];

    if (meshC.needsRefinement(i)) {
      // ================================
      // REFINED ELEMENT
      // ================================

      elTplgy.clear();
      elProj[type]->refineConnectivity(connC, n, elTplgy, dupLocIndicesF);

      std::size_t e = elTypeF.size();
      children.clear();
      for (auto &conn : elTplgy) {
        elTypeF.push_back(type);
        fatherF.push_back(static_cast<unsigned>(i));  // father = coarse index
        elTplgyF.emplace_back(std::move(conn));
        children.push_back(static_cast<unsigned>(e++));
      }
      meshC.setChildren(i, children);

      // set level of children: parent level + 1
      elLevelF.resize(elLevelF.size() + elProj[type]->childNumber(), levelC + 1u);

      const unsigned cs = elProj[type]->coarseSize();

      for (unsigned d = 0; d < dim; ++d) {
        U0[d].resize(cs);
        for (unsigned j = 0; j < cs; ++j) {
          U0[d][j] = XC[d][ connC[j] ];
        }

        elProj[type]->project(U0[d], U1[d]);

        const std::size_t oldSize = XF[d].size();
        XF[d].resize(oldSize + U1[d].size());
        std::copy(U1[d].begin(), U1[d].end(), XF[d].begin() + oldSize);
      }
    }
    else {
      // ================================
      // UNREFINED ELEMENT
      // ================================

      const std::size_t e = elTypeF.size();          // index of identity child

      elTypeF.push_back(type);
      elTplgyF.push_back(connC);
      elLevelF.push_back(levelC);                    // same level
      fatherF.push_back(static_cast<unsigned>(i));   // father = coarse index

      children.clear();
      children.push_back(static_cast<unsigned>(e));
      meshC.setChildren(i, children);

      // boundary coarse nodes as candidates for dedup
      const unsigned cs = static_cast<unsigned>(connC.size());
      if (cs > 0) {
        for (unsigned j = 0; j + 1 < cs; ++j) {      // skip last (center)
          dupLocIndicesF.push_back(connC[j]);
        }
      }
    }
  }

  const unsigned N_before_dedup = static_cast<unsigned>(XF[0].size());

  dedupNodesHash(meshF, dupLocIndicesF);

  const unsigned N_after = static_cast<unsigned>(XF[0].size());
  //std::cout << "Node Final Number = "        << N_after              << std::endl;
  //std::cout << "Node Deduplicated Number = " << (N_before_dedup - N_after) << std::endl;

  // On the new fine mesh, next step default = uniform refinement
  meshF.setUniformRefinement();
  meshF.buildNodeToElementAdjacency();

  // Full sanity checks (can be ifdef'ed out in release)
  meshC.checkAll("refineAndProjectMesh: coarse");
  meshF.checkAll("refineAndProjectMesh: fine");
}









