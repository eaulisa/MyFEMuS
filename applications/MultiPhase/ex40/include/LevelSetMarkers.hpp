#pragma once

#include "MultiLevelSolution.hpp"
#include "MyMatrix.hpp"

#include <gperftools/profiler.h>
#include <stdexcept>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

using namespace femus;

#include "ElementTopology.hpp"
#include "RefineElement.hpp"

class LevelSetMarkers {
public:
  LevelSetMarkers(std::string name, uint dim, uint lmax1 = 5)
      : _name(name), _lmax(lmax1){
        if(dim == 3) {
          _refineElement[0][0] = new RefineElement(lmax1, "hex", "linear", "fifth", "fifth", "legendre");
          _refineElement[0][1] = new RefineElement(lmax1, "hex", "quadratic", "fifth", "fifth", "legendre");
          _refineElement[0][2] = new RefineElement(lmax1, "hex", "biquadratic", "fifth", "fifth", "legendre");

          _refineElement[1][0] = new RefineElement(lmax1, "tet", "linear", "third", "third", "legendre");
          _refineElement[1][1] = new RefineElement(lmax1, "tet", "quadratic", "third", "third", "legendre");
          _refineElement[1][2] = new RefineElement(lmax1, "tet", "biquadratic", "third", "third", "legendre");
        }
        else if (dim == 2) {
          _refineElement[3][0] = new RefineElement(lmax1, "quad", "linear", "fifth", "fifth", "legendre");
          _refineElement[3][1] = new RefineElement(lmax1, "quad", "quadratic", "fifth", "fifth", "legendre");
          _refineElement[3][2] = new RefineElement(lmax1, "quad", "biquadratic", "fifth", "fifth", "legendre");

          _refineElement[4][0] = new RefineElement(lmax1, "tri", "linear", "fifth", "fifth", "legendre");
          _refineElement[4][1] = new RefineElement(lmax1, "tri", "quadratic", "fifth", "fifth", "legendre");
          _refineElement[4][2] = new RefineElement(lmax1, "tri", "biquadratic", "fifth", "fifth", "legendre");

        }
      };

  inline void GetCutElementPoints(MultiLevelSolution &mlSol, std::vector<MyVector<double>> &X, MyVector<unsigned> &Xiel) {

    MultiLevelMesh &mlMsh = *mlSol.GetMultilevelMesh();
    const unsigned level = mlMsh.GetNumberOfLevels() - 1u;
    Mesh &msh = *mlMsh.GetLevel(level);
    Solution &sol = *mlSol.GetLevel(level);
    const unsigned iproc = msh.processor_id();
    const unsigned dim = msh.GetDimension();

    const unsigned solIndex = mlSol.GetIndex(_name.c_str());
    const unsigned solType = mlSol.GetSolutionType(_name.c_str());
    const unsigned xType = 2u;
    if (solType != 2u) {
      throw std::runtime_error("LevelSetMarkers::GetCutElementPoints: solType is not biquadratic ");
    }

    auto &xv = msh._topology->_Sol;
    auto &solVec = sol._Sol[solIndex];

    const unsigned offset = msh._elementOffset[iproc];
    const unsigned offsetp1 = msh._elementOffset[iproc + 1];

    std::vector<std::vector<double>> x(dim);

    std::vector<double> phi;
    std::vector<double> gradPhi(dim);

    std::vector<std::vector<double>> Y(dim);
    std::vector<unsigned> Yiel;

    // TODO Add reservation to std::vector

    // Loop over local elements and collect interior points associated with cut
    // elements
    for (unsigned iel = offset; iel < offsetp1; ++iel) {

      const unsigned nDof = msh.GetElementDofNumber(iel, solType);

      const unsigned solDof0 = msh.GetSolutionDof(0, iel, solType);
      const double val0 = (*solVec)(solDof0);
      const int sign0 = (val0 > 0.) - (val0 < 0.);

      for (unsigned i = 1; i < nDof; ++i) {
        const unsigned solDofi = msh.GetSolutionDof(i, iel, solType);
        const double vali = (*solVec)(solDofi);
        const int signi = (vali > 0.) - (vali < 0.);

        if (signi != sign0) {

          // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          // fill dof values and coordinates
          phi.resize(nDof + 1u);
          for (unsigned k = 0; k < dim; ++k) {
            x[k].resize(nDof + 1u);
          }
          for (unsigned j = 0; j < nDof; ++j) {
            const unsigned solDofj = msh.GetSolutionDof(j, iel, solType);
            phi[j] = (*solVec)(solDofj);
            const unsigned xDof = msh.GetSolutionDof(j, iel, xType);
            for (unsigned k = 0; k < dim; ++k) {
              x[k][j] = (*xv[k])(xDof);
            }
          }
          // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

          // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          // get edge intersections
          auto el_type = msh.el->GetElementType(iel);
          const ex40_topo::ElementTopology &topo =
              ex40_topo::topologyFromTypeId(el_type);

          std::vector<std::vector<double>> element_roots;
          element_roots.resize(dim);
          for (unsigned k = 0; k < topo.nEdges; ++k) {
            const ex40_topo::EdgeInfo &E = topo.edges[k];
            if (E.nodes.n != 3 || E.vertices.n != 2)
              throw std::runtime_error("LevelSetMarkers::GetCutElementPoints: expected quadratic "
                                       "edge (3 nodes, 2 vertices)");

            const auto edge_roots = computeEdgeIntersections(x, phi, E);
            for (int idim = 0; idim < dim; idim++)
              element_roots[idim].insert(element_roots[idim].end(),
                                         edge_roots[idim].begin(),
                                         edge_roots[idim].end());
          }
          // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

          // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          // compute markers
          const auto element_markers =
              computeElementMarkers(element_roots, x, phi, topo);

          for (unsigned j = 0; j < element_markers[0].size(); ++j) {
            for (unsigned k = 0; k < dim; ++k) {
              Y[k].push_back(element_markers[k][j]);
            }
            Yiel.push_back(iel);
          }
          // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          break;
        }
      }
    }

    X.resize(dim);
    for (unsigned k = 0; k < dim; ++k) {
      X[k].buildFromLocal(Y[k]);
    }
    Xiel.buildFromLocal(Yiel);
  }

  inline std::vector<std::vector<double>>
  computeEdgeIntersections(const std::vector<std::vector<double>> &x,
                           const std::vector<double> &phi,
                           const ex40_topo::EdgeInfo &E) {
    int dim = x.size();
    std::vector<std::vector<double>> edge_roots(dim);

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // initialize edge values
    const std::array<uint, 3> edge_id = {E.vertices[0], E.nodes.p[2],
                                         E.vertices[1]};

    std::vector<double> edge_values(3);
    std::vector<std::vector<double>> edge_coords(dim,
                                                 std::vector<double>(3, 0));

    for (int idof = 0; idof < 3; idof++) {
      for (int idim = 0; idim < dim; idim++) {
        edge_coords[idim][idof] = x[idim][edge_id[idof]];
      }
      edge_values[idof] = phi[edge_id[idof]];
    }
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // compute edge intersections
    const double a =
        2.0 * (edge_values[0] + edge_values[2] - 2.0 * edge_values[1]);
    const double b =
        4.0 * edge_values[1] - 3.0 * edge_values[0] - edge_values[2];
    const double c = edge_values[0];

    double roots[2];
    unsigned nRoots = 0;
    const double epsT = 1e-12;

    auto pushRootIfIn01 = [&](double t) {
      if (t >= -epsT && t <= 1.0 + epsT) {
        if (t < 0.0)
          t = 0.0;
        if (t > 1.0)
          t = 1.0;
        roots[nRoots++] = t;
      }
    };

    if (std::abs(a) < 1e-16) {
      if (std::abs(b) > 1e-16) {
        const double t = -c / b;
        pushRootIfIn01(t);
      }
    } else {
      const double disc = b * b - 4.0 * a * c;
      if (disc >= -1e-16) {
        const double sdisc = std::sqrt(std::max(0.0, disc));
        const double inv2a = 1.0 / (2.0 * a);
        pushRootIfIn01((-b - sdisc) * inv2a);
        pushRootIfIn01((-b + sdisc) * inv2a);
      }
    }

    if (nRoots == 0)
      return edge_roots;

    for (unsigned iroot = 0; iroot < nRoots; iroot++) {

      double t = roots[iroot];
      const double N0 = 2.0 * (t - 0.5) * (t - 1.0);
      const double Nm = 4.0 * t * (1.0 - t);
      const double N1 = 2.0 * t * (t - 0.5);

      const std::size_t ip = edge_roots[0].size();

      for (std::size_t idim = 0; idim < dim; ++idim) {
        const double x0 = x[idim][edge_id[0]];
        const double xm = x[idim][edge_id[1]];
        const double x1 = x[idim][edge_id[2]];
        edge_roots[idim].push_back(N0 * x0 + Nm * xm + N1 * x1);
      }
    }
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    return edge_roots;
  }

  inline std::vector<std::vector<double>>
  computeElementMarkers(std::vector<std::vector<double>> element_roots,
                        const std::vector<std::vector<double>> &x,
                        const std::vector<double> &phi,
                        const ex40_topo::ElementTopology &topo) {
    int dim = element_roots.size();
    std::vector<std::vector<double>> element_markers(dim);

    double longest_edge = 0.;
    for (unsigned k = 0; k < topo.nEdges; ++k) {
      double edge_length = 0.;
      const ex40_topo::EdgeInfo &E = topo.edges[k];
      const unsigned edge_id0 = E.vertices[0];
      const unsigned edge_id1 = E.vertices[1];
      for (unsigned idim = 0; idim < dim; ++idim)        
        edge_length += (x[idim][edge_id1] - x[idim][edge_id0])*(x[idim][edge_id1] - x[idim][edge_id0]);
      longest_edge = std::max(longest_edge, std::sqrt(edge_length));
    }

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // deduplicate roots
    std::vector<std::size_t> old_to_new;
    auto unique_element_roots =
        deduplicatePoints(element_roots, &old_to_new, 1e-12);
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    bool simple = (dim == 2 && unique_element_roots[0].size() <= 2) ||
                  (dim == 3 && (unique_element_roots[0].size() >= 3 &&
                                unique_element_roots[0].size() <= 5));

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // simple configuration
    if (simple) {
      placeSimpleCutMarkers(unique_element_roots, element_markers, longest_edge);
    }
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // complex configuration
    else {

      RefineElement *refine = _refineElement[topo.typeId][2];
      const auto &PMatrix   = refine->GetProlongationMatrix();
      unsigned nChildren    = refine->GetNumberOfChildren();

      struct RefinementNode {
        std::vector<std::vector<double>> x;
        std::vector<double>              phi;
        int                              level;
      };

      std::vector<RefinementNode> frontier;
      frontier.push_back({x, phi, 0});

      std::vector<RefinementNode> simple_leaves;

      while (!frontier.empty()) {
        RefinementNode parent = std::move(frontier.back());
        frontier.pop_back();

        auto x_children   = prolongCoords(parent.x, PMatrix);
        auto phi_children = prolongPhi(parent.phi, PMatrix);

        for (unsigned i = 0; i < nChildren; i++) {

          std::vector<std::vector<double>> child_roots;
          child_roots.resize(dim);
          for (unsigned k = 0; k < topo.nEdges; ++k) {
            const ex40_topo::EdgeInfo &E = topo.edges[k];
            if (E.nodes.n != 3 || E.vertices.n != 2)
              throw std::runtime_error("LevelSetMarkers::GetCutElementPoints: expected quadratic "
                                       "edge (3 nodes, 2 vertices)");
            const auto edge_roots = computeEdgeIntersections(x_children[i], phi_children[i], E);
            for (int idim = 0; idim < dim; idim++)
              child_roots[idim].insert(child_roots[idim].end(),
                                       edge_roots[idim].begin(),
                                       edge_roots[idim].end());
          }

          std::vector<std::size_t> old_to_new_child;
          auto unique_child_roots = deduplicatePoints(child_roots, &old_to_new_child, 1e-12);

          bool child_simple =
              (dim == 2 && unique_child_roots[0].size() <= 2) ||
              (dim == 3 && unique_child_roots[0].size() >= 3 && unique_child_roots[0].size() <= 5);

          if (child_simple) {
            simple_leaves.push_back({x_children[i], phi_children[i], parent.level + 1});
          } else if (parent.level + 1 < _lmax) {
            frontier.push_back({x_children[i], phi_children[i], parent.level + 1});
          } else {
            throw std::runtime_error("levelSetMarkers::computeElementMarkers: COMPLEX CUT ARRIVED AT MAX LEVEL");
          }
        }
      }

      for (unsigned leaf_idx = 0; leaf_idx < simple_leaves.size(); leaf_idx++) {
        const auto &leaf = simple_leaves[leaf_idx];

        std::vector<std::vector<double>> leaf_roots;
        leaf_roots.resize(dim);
        for (unsigned k = 0; k < topo.nEdges; ++k) {
          const ex40_topo::EdgeInfo &E = topo.edges[k];
          if (E.nodes.n != 3 || E.vertices.n != 2)
            throw std::runtime_error("LevelSetMarkers::computeElementMarkers: expected quadratic "
                                     "edge (3 nodes, 2 vertices)");
          const auto edge_roots = computeEdgeIntersections(leaf.x, leaf.phi, E);
          for (int idim = 0; idim < dim; idim++)
            leaf_roots[idim].insert(leaf_roots[idim].end(),
                                    edge_roots[idim].begin(),
                                    edge_roots[idim].end());
        }

        std::vector<std::size_t> old_to_new_leaf;
        auto unique_leaf_roots = deduplicatePoints(leaf_roots, &old_to_new_leaf, 1e-12);

        if (unique_leaf_roots[0].size() > 0) {
          placeSimpleCutMarkers(unique_leaf_roots, element_markers, longest_edge);
        }
      }
    }
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    return element_markers;
  }

  inline void placeSimpleCutMarkers( const std::vector<std::vector<double>> &roots, std::vector<std::vector<double>> &markers, double longest_edge) {

    std::vector<double> b(roots.size(), 0.);
    for (int idim = 0; idim < (int)roots.size(); idim++) {
      for (int iroot = 0; iroot < (int)roots[0].size(); iroot++) {
        b[idim] += roots[idim][iroot] / roots[idim].size();
      }

      markers[idim].push_back(b[idim]);
    }

    
    for (int iroot = 0; iroot < (int)roots[0].size(); iroot++) {
      double root_2_b_dist = 0.;
      for (int idim = 0; idim < (int)roots.size(); idim++) {
        root_2_b_dist += (roots[idim][iroot] - b[idim]) * (roots[idim][iroot] - b[idim]);
      }
      root_2_b_dist = std::sqrt(root_2_b_dist);

      const unsigned nmarkers = (root_2_b_dist > longest_edge/3.) ? 2 : (root_2_b_dist > longest_edge/6.) ? 1 : 0;
      const double offset0 = 1 / (nmarkers + 0.5);
      double offset = offset0;

      for (int imarker = 0; imarker < (int)nmarkers; imarker++) {
        for (int idim = 0; idim < (int)roots.size(); idim++) {
          markers[idim].push_back(
              b[idim] + offset * (roots[idim][iroot] - b[idim]));
        }
        offset += offset0;
      }
    }

  }

  // ================================================================
  // AUXILIARY FUNCTIONS
  // ================================================================

  std::vector<std::vector<double>>
  deduplicatePoints(const std::vector<std::vector<double>> &pts,
                    std::vector<std::size_t> *old_to_new = nullptr,
                    double tol = 1e-12) {
    if (pts.empty() || pts[0].empty())
      return pts;

    const std::size_t dim = pts.size();
    const std::size_t nPts = pts[0].size();

    std::vector<std::size_t> idx(nPts);
    std::iota(idx.begin(), idx.end(), 0);

    std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) {
      for (std::size_t d = 0; d < dim; ++d) {
        if (pts[d][a] < pts[d][b] - tol)
          return true;
        if (pts[d][a] > pts[d][b] + tol)
          return false;
      }
      return false;
    });

    std::vector<std::size_t> map(nPts);
    std::vector<std::size_t> uniqueIdx;
    uniqueIdx.reserve(nPts);

    auto samePoint = [&](std::size_t a, std::size_t b) {
      for (std::size_t d = 0; d < dim; ++d)
        if (std::abs(pts[d][a] - pts[d][b]) > tol)
          return false;
      return true;
    };

    for (std::size_t k = 0; k < nPts; ++k) {
      const std::size_t i = idx[k];
      if (k == 0 || !samePoint(i, idx[k - 1])) {
        uniqueIdx.push_back(i);
      }
      map[i] = uniqueIdx.size() - 1;
    }

    const std::size_t nUnique = uniqueIdx.size();
    std::vector<std::vector<double>> result(dim, std::vector<double>(nUnique));
    for (std::size_t d = 0; d < dim; ++d)
      for (std::size_t k = 0; k < nUnique; ++k)
        result[d][k] = pts[d][uniqueIdx[k]];

    if (old_to_new)
      *old_to_new = std::move(map);

    return result;
  }

  // ---------------------------------------------------------------
  // Apply PMatrix to prol phi values from parent to children.
  // phi_parent[j] = value on node j of the parent
  // Returns phi_children[i][j] = value on node j of child i
  // ---------------------------------------------------------------
  static std::vector<std::vector<double>> prolongPhi(
      const std::vector<double> &phi_parent,
      const std::vector<std::vector<std::vector<std::pair<unsigned, double>>>> &PMatrix) {

    unsigned nChildren = PMatrix.size();
    unsigned nNodes    = PMatrix[0].size();
    std::vector<std::vector<double>> phi_children(nChildren,
                                                  std::vector<double>(nNodes, 0.));
    for (unsigned i = 0; i < nChildren; i++)
      for (unsigned j = 0; j < nNodes; j++)
        for (const auto &p : PMatrix[i][j])
          phi_children[i][j] += p.second * phi_parent[p.first];

    return phi_children;
  }

  // ---------------------------------------------------------------
  // Apply PMatrix to prol coords values from parent to children.
  // x_parent[k][j] = coordinate k of node j of the parent
  // Returns x_children[i][k][j] = coordinate k of node j of child i
  // ---------------------------------------------------------------
  static std::vector<std::vector<std::vector<double>>> prolongCoords(
      const std::vector<std::vector<double>> &x_parent,
      const std::vector<std::vector<std::vector<std::pair<unsigned, double>>>> &PMatrix) {

    unsigned nChildren = PMatrix.size();
    unsigned nNodes    = PMatrix[0].size();
    unsigned dim       = x_parent.size();

    std::vector<std::vector<std::vector<double>>> x_children(
        nChildren, std::vector<std::vector<double>>(dim, std::vector<double>(nNodes, 0.)));

    for (unsigned i = 0; i < nChildren; i++)
      for (unsigned k = 0; k < dim; k++)
        for (unsigned j = 0; j < nNodes; j++)
          for (const auto &p : PMatrix[i][j])
            x_children[i][k][j] += p.second * x_parent[k][p.first];

    return x_children;
  }

private:
  const std::string _name;
  const unsigned _lmax;
  RefineElement *_refineElement[6][3];
  int it = 0;
};