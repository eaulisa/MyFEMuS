void dedupNodes(
  unsigned dim,
  std::vector<std::vector<double>>& XF,
  std::vector<std::vector<unsigned>>& elTplgyF,
  std::vector<unsigned>& elTypeF,
  unsigned& n,
  std::vector<unsigned>& candidateIndices
) {
  // --- basic sanity on dim / XF ---
  if (dim != 2 && dim != 3) {
    throw std::runtime_error("dedupNodes: dim must be 2 or 3");
  }
  if (XF.size() < dim) {
    throw std::runtime_error("dedupNodes: XF has fewer components than dim");
  }

  const unsigned N = static_cast<unsigned>(XF[0].size());
  if (N == 0) {
    n = 0;
    return;
  }

  for (unsigned d = 1; d < dim; ++d) {
    if (XF[d].size() != XF[0].size()) {
      throw std::runtime_error("dedupNodes: XF[d] size mismatch among components");
    }
  }

  // -------- 1) Make candidateIndices unique by index --------
  std::sort(candidateIndices.begin(), candidateIndices.end());
  candidateIndices.erase(
    std::unique(candidateIndices.begin(), candidateIndices.end()),
    candidateIndices.end()
  );

  const unsigned K = static_cast<unsigned>(candidateIndices.size());
  if (K == 0) {
    n = N;
    return;
  }

  // -------- 2) Build snapped integer coordinates (lexicographic with tolerance) --------
  //
  // Any two nodes whose coordinates are within ~h_snap/2 in each component
  // will snap to the same integer in that component.
  //
  const double h_snap = computeHSnapElementBased(XF, dim, elTplgyF, elTypeF, 1.e-2);   // relative to smallest element

  std::cout << "characteristic size = " << h_snap << std::endl;

  // snap[d][i] = snapped integer coord of node i in dimension d
  std::vector<std::vector<long long>> snap(dim, std::vector<long long>(N));
  for (unsigned d = 0; d < dim; ++d) {
    for (unsigned i = 0; i < N; ++i) {
      snap[d][i] = static_cast<long long>(std::llround(XF[d][i] / h_snap));
    }
  }

  // perm is the list of candidate indices to be ordered geometrically (on snapped coords)
  std::vector<unsigned> perm = candidateIndices;

  auto coordLess = [&](unsigned a, unsigned b) {
    // lexicographic compare on snapped coordinates
    for (unsigned d = 0; d < dim; ++d) {
      if (snap[d][a] < snap[d][b]) return true;
      if (snap[d][a] > snap[d][b]) return false;
    }
    return a < b; // tie-break by original index
  };

  std::sort(perm.begin(), perm.end(), coordLess);

  auto coordEqual = [&](unsigned a, unsigned b) {
    // equality is exact equality of snapped ints → transitive, no tolerance chains
    for (unsigned d = 0; d < dim; ++d) {
      if (snap[d][a] != snap[d][b]) return false;
    }
    return true;
  };

// -------- 3) Build representative mapping rep[i] --------
// Initialise all nodes as their own representative
  std::vector<unsigned> rep(N);
  std::iota(rep.begin(), rep.end(), 0);

// Only candidates can be merged. Interior nodes keep rep[i] = i.
  unsigned k = 0;
  while (k < K) {
    // perm[k] is the first node in this snapped-equal group
    const unsigned base = perm[k];

    // Find run [k, k2) of candidates equal (in snapped space) to perm[k]
    unsigned k2 = k + 1;
    while (k2 < K && coordEqual(perm[k], perm[k2])) {
      ++k2;
    }

    // Assign all nodes in the group to this representative
    for (unsigned j = k; j < k2; ++j) {
      const unsigned idx = perm[j];
      rep[idx] = base;
    }

    k = k2;
  }


  // -------- 4) Assign new indices 0..numUnique-1 to each representative --------
  std::vector<unsigned> old2new(N);
  unsigned numUnique = 0;

  for (unsigned i = 0; i < N; ++i) {
    if (rep[i] == i) {
      // i is a representative (includes interior nodes and boundary masters)
      old2new[i] = numUnique++;
    }
  }
  for (unsigned i = 0; i < N; ++i) {
    if (rep[i] != i) {
      old2new[i] = old2new[rep[i]];
    }
  }

  // If no duplicates were actually found, numUnique == N, nothing to do
  if (numUnique == N) {
    n = N;
    return;
  }

  // -------- 5) Update connectivity --------
  for (auto& conn : elTplgyF) {
    for (auto& v : conn) {
      v = old2new[v];
    }
  }

  // -------- 6) Compact XF to size numUnique --------
  // Keep the representative's coordinates explicitly
  for (unsigned d = 0; d < dim; ++d) {
    std::vector<double> Xnew(numUnique);

    for (unsigned i = 0; i < N; ++i) {
      if (rep[i] == i) {           // i is representative
        unsigned j = old2new[i];
        Xnew[j] = XF[d][i];
      }
    }

    XF[d].swap(Xnew);
  }

  n = numUnique;
}
