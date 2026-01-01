#pragma once

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <stdexcept>
#include <memory>
#include <string>
#include <array>
#include <cmath>
#include <iostream>
#include "Fem.hpp"


// convenience aliases for Eigen
using SpMat   = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

// ====================== Base class ======================
//
// - Owns _P (projection matrix), _coarseSize, _fineSize.
// - Owns _pattern (refinement pattern) and _femName.
// - Derived classes fill _P and _pattern in their constructors.
// - User only calls coarseSize(), fineSize(), femName(),
//   project(...), refineConnectivity(...) on FemProjection*.
//
// ========================================================

class FemProjection {
  public:
    virtual ~FemProjection() {}

    unsigned coarseSize() const {
      return _coarseSize;
    }
    unsigned fineSize()   const {
      return _fineSize;
    }
    const std::string& femName() const {
      return _femName;
    }
    unsigned childNumber()   const {
      return _numChildren;
    }


    const SpMat& GetProjection() const {
      return _P;
    }

    // Project coarse -> fine (values).
    // Requires: coarse.size() == coarseSize().
    // On return: fine.size() == fineSize().
    void project(const std::vector<double>& coarse,
                 std::vector<double>& fine) const {
      if (coarse.size() != _coarseSize) {
        throw std::runtime_error(
          "FemProjection::project(" + _femName +
          "): coarse size must be " + std::to_string(_coarseSize) + "."
        );
      }

      Eigen::Map<const Eigen::VectorXd> uc(coarse.data(), _coarseSize);
      Eigen::VectorXd uf = _P * uc;  // _fineSize entries

      fine.resize(_fineSize);
      Eigen::Map<Eigen::VectorXd>(fine.data(), _fineSize) = uf;
    }

    // Build refined connectivity (topology).
    void refineConnectivity(const std::vector<unsigned>& coarse_nodes,
                            unsigned& n_nodes,
                            std::vector<std::vector<unsigned>>& children,
                            std::vector<unsigned>& dupGlobIndices) const {
      if (coarse_nodes.size() != _coarseSize) {
        throw std::runtime_error(
          "FemProjection::refineConnectivity(" + _femName +
          "): coarse_nodes size must be " +
          std::to_string(_coarseSize) + "."
        );
      }

      const unsigned base       = n_nodes;                      // first new global node
      const unsigned C          = _numChildren;
      const unsigned m          = _coarseSize;                  // nodes per child (Q2)
      const unsigned firstNew   = _coarseSize;

      children.resize(C);

      for (unsigned c = 0; c < C; ++c) {
        std::vector<unsigned>& child = children[c];
        child.resize(m);

        const unsigned offset = c * m;
        for (unsigned j = 0; j < m; ++j) {
          const unsigned v = _childLocalConn[offset + j];
          if (v < _coarseSize) {
            child[j] = coarse_nodes[v];
          }
          else {
            child[j] = base + (v - firstNew);
          }
        }
      }

      // boundary nodes
      const unsigned start = dupGlobIndices.size();
      dupGlobIndices.resize(start + _boundaryNodes.size());
      for (unsigned i = 0; i < _boundaryNodes.size(); ++i) {
        const unsigned loc = _boundaryNodes[i];
        unsigned g;
        if (loc < _coarseSize) {
          g = coarse_nodes[loc];
        }
        else {
          g = base + (loc - firstNew);
        }
        dupGlobIndices[start + i] = g;
      }

      n_nodes += _fineSize;
    }

    virtual const femus::elem_type& fem() = 0;

  protected:
    // Only derived classes can construct.

    FemProjection(unsigned coarse_size,
                  const std::string& femName)
      : _P()
      , _coarseSize(coarse_size)
      , _fineSize(0u)
      , _femName(femName)
      , _numChildren(0u)
    {}

    SpMat       _P;           // (_fineSize) x (_coarseSize)
    unsigned    _coarseSize;  // # coarse dofs (the nodes in the coarse element)
    unsigned    _fineSize;    // # fine dofs (only the new nodes in the fine patch)
    std::string _femName;

    unsigned              _numChildren;       // # fine elements in patch
    std::vector<unsigned> _childLocalConn;    // flattened [child][a]
    std::vector<unsigned> _boundaryNodes;     // local ids in [0, coarseSize+fineSize)

};

// =============================
// Base class: Q2RefinementProjection
// =============================

class Q2RefinementProjection : public FemProjection {
  public:
    const femus::elem_type& fem() {
      if (!_fem) throw std::runtime_error("Field::fem(): _fem is null");
      return *_fem;
    }
  protected:
    using Node     = std::array<double, 3>;
    using ChildConn = std::vector<std::vector<unsigned>>;

    const unsigned _dim;         // spatial dimension
    const unsigned _nChild;      // number of children
    const unsigned _nVert;       // number of vertex nodes
    const unsigned _nEdge;       // number of edge nodes
    const unsigned _nFace;       // number of face nodes
    const unsigned _nCenter;     // number of center nodes
    const unsigned _nCoarse;     // total coarse nodes
    const unsigned _vertOff;     // offset of vertices
    const unsigned _edgeOff;     // offset of edges
    const unsigned _faceOff;     // offset of faces
    const unsigned _centerOff;   // offset of center
    std::unique_ptr<femus::elem_type>_fem;

    ChildConn _childVerts;       // child vertex connectivity (nChild x nVert)

    Q2RefinementProjection(unsigned dim,
                           unsigned nChild,
                           unsigned nVert,
                           unsigned nEdge,
                           unsigned nFace,
                           unsigned nCenter,
                           const std::string &name)
      : FemProjection(nVert + nEdge + nFace + nCenter, name),
        _dim(dim),
        _nChild(nChild),
        _nVert(nVert),
        _nEdge(nEdge),
        _nFace(nFace),
        _nCenter(nCenter),
        _nCoarse(nVert + nEdge + nFace + nCenter),
        _vertOff(0),
        _edgeOff(_vertOff + nVert),
        _faceOff(_edgeOff + nEdge),
        _centerOff(_faceOff + nFace),
        _childVerts(nChild, std::vector<unsigned>(nVert, 0u)) {}

    virtual void LinearToQuadraticCoordinates(const Node &X0,
                                              const Node &X1,
                                              std::vector<Node> &X) const = 0;

    virtual Node getParentMin() const = 0;
    virtual Node getParentMax() const = 0;

    virtual void buildChildX(const std::vector<Node> &X_C,
                             std::vector<std::vector<Node>> &childX) const = 0;

    virtual bool isBoundaryNode(const std::vector<std::vector<double>> &X_F,
                                unsigned i) const = 0;

    void build() {

      std::vector<Node> X_C;
      LinearToQuadraticCoordinates(getParentMin(), getParentMax(), X_C);

      std::vector<std::vector<Node>> childX;
      buildChildX(X_C, childX);

      std::vector<std::vector<double>> X_F(_dim);
      for (unsigned d = 0; d < _dim; ++d) X_F[d].reserve(_nCoarse * _nChild);

      for (unsigned d = 0; d < _dim; ++d) {
        for (unsigned i = 0; i < _nCoarse; ++i) {
          X_F[d].push_back(X_C[i][d]);
        }
      }

      std::vector<std::vector<unsigned>> childNodes(_nChild, std::vector<unsigned>(_nCoarse, 0u));

      for (unsigned c = 0; c < _nChild; ++c) {
        for (unsigned v = 0; v < _nVert; ++v) {
          childNodes[c][_vertOff + v] = _childVerts[c][v];
        }
      }

      auto same_point = [this](const Node & a, const Node & b) {
        for (unsigned d = 0; d < _dim; ++d) {
          if (std::fabs(a[d] - b[d]) > 1.0e-8) return false;
        }
        return true;
      };

      unsigned nodes = _nCoarse;

      if (_nEdge > 0) {
        for (unsigned e = 0; e < _nEdge; ++e) {
          const Node &Xi = childX[0][_edgeOff + e];
          for (unsigned d = 0; d < _dim; ++d)  X_F[d].push_back(Xi[d]);
          childNodes[0][_edgeOff + e] = nodes++;
        }

        for (unsigned c = 1; c < _nChild; ++c) {
          for (unsigned edgei = 0; edgei < _nEdge; ++edgei) {
            const Node &Xi = childX[c][_edgeOff + edgei];
            bool found = false;

            for (unsigned cp = 0; cp < c && !found; ++cp) {
              for (unsigned edgej = 0; edgej < _nEdge; ++edgej) {
                const Node &Xj = childX[cp][_edgeOff + edgej];
                if (same_point(Xi, Xj)) {
                  childNodes[c][_edgeOff + edgei] = childNodes[cp][_edgeOff + edgej];
                  found = true;
                  break;
                }
              }
            }
            for (unsigned j = 0; j < X_C.size() && !found; ++j) {
              if (same_point(Xi, X_C[j])) {
                childNodes[c][_edgeOff + edgei] = j;
                found = true;
                break;
              }
            }

            if (!found) {
              for (unsigned d = 0; d < _dim; ++d) X_F[d].push_back(Xi[d]);
              childNodes[c][_edgeOff + edgei] = nodes++;
            }
          }
        }
      }

      if (_nFace > 0) {
        for (unsigned f = 0; f < _nFace; ++f) {
          const Node &Xi = childX[0][_faceOff + f];
          for (unsigned d = 0; d < _dim; ++d) X_F[d].push_back(Xi[d]);
          childNodes[0][_faceOff + f] = nodes++;
        }

        for (unsigned c = 1; c < _nChild; ++c) {
          for (unsigned facei = 0; facei < _nFace; ++facei) {
            const Node &Xi = childX[c][_faceOff + facei];
            bool found = false;

            for (unsigned cp = 0; cp < c && !found; ++cp) {
              for (unsigned facej = 0; facej < _nFace; ++facej) {
                const Node &Xj = childX[cp][_faceOff + facej];
                if (same_point(Xi, Xj)) {
                  childNodes[c][_faceOff + facei] = childNodes[cp][_faceOff + facej];
                  found = true;
                  break;
                }
              }
            }
            for (unsigned j = 0; j < X_C.size() && !found; ++j) {
              if (same_point(Xi, X_C[j])) {
                childNodes[c][_faceOff + facei] = j;
                found = true;
                break;
              }
            }

            if (!found) {
              for (unsigned d = 0; d < _dim; ++d) X_F[d].push_back(Xi[d]);
              childNodes[c][_faceOff + facei] = nodes++;
            }
          }
        }
      }

      if (_nCenter > 0) {
        for (unsigned c = 0; c < _nChild; ++c) {
          const Node &Xi = childX[c][_centerOff];
          bool found = false;
          for (unsigned j = 0; j < X_C.size() && !found; ++j) {
            if (same_point(Xi, X_C[j])) {
              childNodes[c][_centerOff] = j;
              found = true;
              break;
            }
          }
          if (!found) {
            for (unsigned d = 0; d < _dim; ++d) X_F[d].push_back(Xi[d]);
            childNodes[c][_centerOff] = nodes++;
          }
        }
      }

      //std::cout << "Number of refined nodes = " << nodes << std::endl;

      _boundaryNodes.clear();
      for (unsigned i = 0; i < nodes; ++i) {
        if (isBoundaryNode(X_F, i)) {
          _boundaryNodes.push_back(i);
        }
      }

      std::vector<double> phiC;
      phiC.reserve(_nCoarse);
      std::vector<double> xi(_dim);

      unsigned nNew = nodes - _nCoarse;

      _P.resize(nNew, _nCoarse);
      std::vector<Triplet> T;
      T.reserve(nNew * _nCoarse);

      for (unsigned i = _nCoarse; i < nodes; ++i) {
        for (unsigned d = 0; d < _dim; ++d) xi[d] = X_F[d][i];

        _fem->GetPhi(phiC, xi);

        const unsigned row = i - _nCoarse;
        for (unsigned j = 0; j < _nCoarse; ++j) {
          const double val = phiC[j];
          if (std::fabs(val) > 1.0e-10) {
            T.emplace_back(static_cast<int>(row),
                           static_cast<int>(j),
                           val);
          }
        }
      }

      _P.setFromTriplets(T.begin(), T.end());

      _coarseSize   = _nCoarse;
      _fineSize     = nNew;
      _numChildren  = _nChild;

      // Flatten childNodes into _childLocalConn:
      _childLocalConn.clear();
      _childLocalConn.reserve(_nChild * _nCoarse);
      for (unsigned c = 0; c < _nChild; ++c) {
        for (unsigned a = 0; a < _nCoarse; ++a) {
          _childLocalConn.push_back(childNodes[c][a]);
        }
      }
    }
};


// =============================
// Intermediate box base: BoxQ2RefinementProjection
// shared by Hex27 and Quad9
// =============================

class BoxQ2RefinementProjection : public Q2RefinementProjection {
  protected:
    using Node = Q2RefinementProjection::Node;

    const unsigned _oppVertIndex;  // only used for box-type elements

    BoxQ2RefinementProjection(unsigned dim,
                              unsigned nChild,
                              unsigned nVert,
                              unsigned nEdge,
                              unsigned nFace,
                              unsigned nCenter,
                              // unsigned nFine,
                              unsigned oppVertIndex,
                              const std::string &name)
      : Q2RefinementProjection(dim, nChild, nVert, nEdge, nFace, nCenter, name),
        _oppVertIndex(oppVertIndex) {}

    void buildChildX(const std::vector<Node> &X_C,
                     std::vector<std::vector<Node>> &childX) const override {
      childX.assign(_nChild, std::vector<Node>());
      for (unsigned c = 0; c < _nChild; ++c) {
        const unsigned v0  = _childVerts[c][0];
        const unsigned vOp = _childVerts[c][_oppVertIndex];
        LinearToQuadraticCoordinates(X_C[v0], X_C[vOp], childX[c]);
      }
    }

    bool isBoundaryNode(const std::vector<std::vector<double>> &X_F,
                        unsigned i) const override {
      for (unsigned d = 0; d < _dim; ++d) {
        if (X_F[d][i] < -0.999 || X_F[d][i] > 0.999) return true;
      }
      return false;
    }
};



// =============================
// Hex27Projection
// =============================

class Hex27Projection : public BoxQ2RefinementProjection {
  public:
    Hex27Projection()
      : BoxQ2RefinementProjection(3,          // dim
                                  8,          // nChild
                                  8,          // nVert
                                  12,         // nEdge
                                  6,          // nFace
                                  1,          // nCenter
                                  // 125,        // nFine
                                  6,          // oppVertIndex
                                  "Hex27") {

      static const unsigned verts[8][8] = {
        {0,  8, 24, 11, 16, 20, 26, 23},
        {8,  1,  9, 24, 20, 17, 21, 26},
        {24, 9,  2, 10, 26, 21, 18, 22},
        {11, 24, 10,  3, 23, 26, 22, 19},
        {16, 20, 26, 23,  4, 12, 25, 15},
        {20, 17, 21, 26, 12,  5, 13, 25},
        {26, 21, 18, 22, 25, 13,  6, 14},
        {23, 26, 22, 19, 15, 25, 14,  7}
      };

      for (unsigned c = 0; c < _nChild; ++c) {
        for (unsigned v = 0; v < _nVert; ++v) {
          _childVerts[c][v] = verts[c][v];
        }
      }
      _fem = std::make_unique<femus::elem_type_3D>("hex", "biquadratic", "eleventh");

      build();
    }

  private:
    Node getParentMin() const override {
      return Node{-1.0, -1.0, -1.0};
    }
    Node getParentMax() const override {
      return Node{ 1.0,  1.0,  1.0};
    }

    void LinearToQuadraticCoordinates(const Node &X0,
                                      const Node &X1,
                                      std::vector<Node> &X) const override {

      const double x0 = X0[0]; const double y0 = X0[1]; const double z0 = X0[2];
      const double x1 = X1[0]; const double y1 = X1[1]; const double z1 = X1[2];

      const double xm = 0.5 * (x0 + x1);
      const double ym = 0.5 * (y0 + y1);
      const double zm = 0.5 * (z0 + z1);

      X = {
        Node{x0, y0, z0}, Node{x1, y0, z0}, Node{x1, y1, z0}, Node{x0, y1, z0},
        Node{x0, y0, z1}, Node{x1, y0, z1}, Node{x1, y1, z1}, Node{x0, y1, z1},
        Node{xm, y0, z0}, Node{x1, ym, z0}, Node{xm, y1, z0}, Node{x0, ym, z0},
        Node{xm, y0, z1}, Node{x1, ym, z1}, Node{xm, y1, z1}, Node{x0, ym, z1},
        Node{x0, y0, zm}, Node{x1, y0, zm}, Node{x1, y1, zm}, Node{x0, y1, zm},
        Node{xm, y0, zm}, Node{x1, ym, zm}, Node{xm, y1, zm}, Node{x0, ym, zm},
        Node{xm, ym, z0}, Node{xm, ym, z1},
        Node{xm, ym, zm}
      };
    }
};



// =============================
// Quad9Projection
// =============================

class Quad9Projection : public BoxQ2RefinementProjection {
  public:
    Quad9Projection()
      : BoxQ2RefinementProjection(2,          // dim
                                  4,          // nChild
                                  4,          // nVert
                                  4,          // nEdge
                                  0,          // nFace
                                  1,          // nCenter
                                  //  25,         // nFine
                                  2,          // oppVertIndex
                                  "Quad9") {

      static const unsigned verts[4][4] = {
        {0, 4, 8, 7},
        {4, 1, 5, 8},
        {8, 5, 2, 6},
        {7, 8, 6, 3}
      };

      for (unsigned c = 0; c < _nChild; ++c) {
        for (unsigned v = 0; v < _nVert; ++v) {
          _childVerts[c][v] = verts[c][v];
        }
      }

      _fem = std::make_unique<femus::elem_type_2D>("quad", "biquadratic", "eleventh");

      build();
    }

  private:
    Node getParentMin() const override {
      return Node{-1.0, -1.0, 0.0};
    }
    Node getParentMax() const override {
      return Node{ 1.0,  1.0, 0.0};
    }

    void LinearToQuadraticCoordinates(const Node &X0,
                                      const Node &X1,
                                      std::vector<Node> &X) const override {

      const double x0 = X0[0];
      const double y0 = X0[1];
      const double x1 = X1[0];
      const double y1 = X1[1];

      const double xm = 0.5 * (x0 + x1);
      const double ym = 0.5 * (y0 + y1);
      const double z  = 0.0;

      X = {
        Node{x0, y0, z}, Node{x1, y0, z}, Node{x1, y1, z}, Node{x0, y1, z},
        Node{xm, y0, z}, Node{x1, ym, z}, Node{xm, y1, z}, Node{x0, ym, z},
        Node{xm, ym, z}
      };
    }
};

class Line3Projection : public BoxQ2RefinementProjection {
  public:
    Line3Projection()
      : BoxQ2RefinementProjection(1,          // dim
                                  2,          // nChild
                                  2,          // nVert
                                  1,          // nEdge
                                  0,          // nFace
                                  0,          // nCenter
                                  1,          // oppVertIndex
                                  "Line3") {

      static const unsigned verts[2][2] = {
        {0, 2},
        {2, 1}
      };

      for (unsigned c = 0; c < _nChild; ++c) {
        for (unsigned v = 0; v < _nVert; ++v) {
          _childVerts[c][v] = verts[c][v];
        }
      }

      _fem = std::make_unique<femus::elem_type_1D>("line", "biquadratic", "eleventh");

      build();
    }

  private:
    Node getParentMin() const override {
      return Node{-1.0, 0.0, 0.0};
    }
    Node getParentMax() const override {
      return Node{ 1.0, 0.0, 0.0};
    }

    void LinearToQuadraticCoordinates(const Node &X0,
                                      const Node &X1,
                                      std::vector<Node> &X) const override {

      const double x0 = X0[0];
      const double x1 = X1[0];

      const double xm = 0.5 * (x0 + x1);
      const double y  = 0.0;
      const double z  = 0.0;

      X = {
        Node{x0, y, z},  Node{x1, y, z}, Node{xm, y, z}
      };
    }
};

// =============================
// Tri7Projection
// =============================

class Tri7Projection : public Q2RefinementProjection {
  public:
    Tri7Projection()
      : Q2RefinementProjection(2,          // dim
                               4,          // nChild
                               3,          // nVert
                               3,          // nEdge
                               0,          // nFace
                               1,          // nCenter
                               //   19,         // nFine
                               "Tri7") {

      static const unsigned verts[4][3] = {
        {0, 3, 5},
        {3, 1, 4},
        {5, 4, 2},
        {4, 5, 3}
      };

      for (unsigned c = 0; c < _nChild; ++c) {
        for (unsigned v = 0; v < _nVert; ++v) {
          _childVerts[c][v] = verts[c][v];
        }
      }

      _fem = std::make_unique<femus::elem_type_2D>("tri", "biquadratic", "eleventh");

      build();
    }

  private:
    using Node = Q2RefinementProjection::Node;

    Node getParentMin() const override {
      return Node{0.0, 0.0, 0.0};
    }
    Node getParentMax() const override {
      return Node{1.0, 1.0, 0.0};
    }

    void buildChildX(const std::vector<Node> &X_C,
                     std::vector<std::vector<Node>> &childX) const override {
      childX.assign(_nChild, std::vector<Node>());
      for (unsigned c = 0; c < _nChild; ++c) {
        const Node v0 = X_C[_childVerts[c][0]];
        const Node v1 = X_C[_childVerts[c][1]];
        const Node v2 = X_C[_childVerts[c][2]];

        LinearToQuadraticCoordinatesUsingVertices(v0, v1, v2, childX[c]);

      }
    }

    bool isBoundaryNode(const std::vector<std::vector<double>> &X_F,
                        unsigned i) const override {
      const double x    = X_F[0][i];
      const double y    = X_F[1][i];
      const double zeta = 1.0 - x - y;
      return (x < 1.0e-4 || y < 1.0e-4 || zeta < 1.0e-4);
    }

    void LinearToQuadraticCoordinates(const Node &X0,
                                      const Node &X1,
                                      std::vector<Node> &X) const override {

      const Node v0 = X0;
      const Node v1 = {X1[0], X0[1], 0.};
      const Node v2 = {X0[0], X1[1], 0.};

      LinearToQuadraticCoordinatesUsingVertices(v0, v1, v2, X);
    }

    void LinearToQuadraticCoordinatesUsingVertices(const Node &v0, const Node &v1,
                                                   const Node &v2,
                                                   std::vector<Node> &X) const {

      auto Barycenter = [this](std::initializer_list<Node> pts) {
        Node x = {0., 0., 0.};
        const double n = static_cast<double>(pts.size());
        for (unsigned d = 0; d < _dim; ++d) {
          double sum = 0.0;
          for (const auto &p : pts) {
            sum += p[d];          // sum coordinate d over all points
          }
          x[d] = sum / n;         // average in dimension d
        }
        return x;
      };

      const Node e0  = Barycenter({v0, v1});
      const Node e1  = Barycenter({v1, v2});
      const Node e2  = Barycenter({v2, v0});
      const Node c0 = Barycenter({v0, v1, v2});

      X = { v0, v1, v2,
            e0, e1, e2, c0
          };
    }
};


class Tet15Projection : public Q2RefinementProjection {
  public:
    Tet15Projection()
      : Q2RefinementProjection(3,          // dim
                               8,          // nChild
                               4,          // nVert
                               6,          // nEdge
                               4,          // nFace
                               1,          // nCenter
                               "Tet15") {

      static const unsigned verts[8][4] = { //nuova numerazione
        {0, 4, 6, 7},
        {4, 1, 5, 8},
        {6, 5, 2, 9},
        {7, 8, 9, 3},
        {5, 6, 4, 7},
        {8, 7, 5, 4},
        {7, 9, 8, 5},
        {9, 5, 7, 6}
      };
      for (unsigned c = 0; c < _nChild; ++c) {
        for (unsigned v = 0; v < _nVert; ++v) {
          _childVerts[c][v] = verts[c][v];
        }
      }

      _fem = std::make_unique<femus::elem_type_3D>("tet", "biquadratic", "eleventh");

      build();
    }

  private:
    using Node = Q2RefinementProjection::Node;

    Node getParentMin() const override {
      return Node{0.0, 0.0, 0.0};
    }
    Node getParentMax() const override {
      return Node{1.0, 1.0, 1.0};
    }

    bool isBoundaryNode(const std::vector<std::vector<double>> &X_F,
                        unsigned i) const override {
      const double x    = X_F[0][i];
      const double y    = X_F[1][i];
      const double z    = X_F[2][i];
      const double zeta = 1.0 - x - y - z;
      return (x < 1.0e-4 || y < 1.0e-4 || z < 1.0e-4 || zeta < 1.0e-4);
    }

    void buildChildX(const std::vector<Node> &X_C,
                     std::vector<std::vector<Node>> &childX) const override {
      childX.assign(_nChild, std::vector<Node>());
      for (unsigned c = 0; c < _nChild; ++c) {
        const Node v0 = X_C[_childVerts[c][0]];
        const Node v1 = X_C[_childVerts[c][1]];
        const Node v2 = X_C[_childVerts[c][2]];
        const Node v3 = X_C[_childVerts[c][3]];

        LinearToQuadraticCoordinatesUsingVertices(v0, v1, v2, v3, childX[c]);

      }
    }

    void LinearToQuadraticCoordinates(const Node &X0,
                                      const Node &X1,
                                      std::vector<Node> &X) const override {

      const Node v0 = X0;
      const Node v1 = {X1[0], X0[1], X0[2]};
      const Node v2 = {X0[0], X1[1], X0[2]};
      const Node v3 = {X0[0], X0[1], X1[2]};

      LinearToQuadraticCoordinatesUsingVertices(v0, v1, v2, v3, X);
    }

    void LinearToQuadraticCoordinatesUsingVertices(const Node &v0, const Node &v1,
                                                   const Node &v2, const Node &v3,
                                                   std::vector<Node> &X) const {

      auto Barycenter = [this](std::initializer_list<Node> pts) {
        Node x = {0., 0., 0.};
        const double n = static_cast<double>(pts.size());
        for (unsigned d = 0; d < _dim; ++d) {
          double sum = 0.0;
          for (const auto &p : pts) {
            sum += p[d];          // sum coordinate d over all points
          }
          x[d] = sum / n;         // average in dimension d
        }
        return x;
      };

      const Node e0  = Barycenter({v0, v1});
      const Node e1  = Barycenter({v1, v2});
      const Node e2  = Barycenter({v2, v0});
      const Node e3  = Barycenter({v0, v3});
      const Node e4  = Barycenter({v1, v3});
      const Node e5  = Barycenter({v2, v3});
      const Node f0 = Barycenter({v0, v1, v2});
      const Node f1 = Barycenter({v0, v1, v3});
      const Node f2 = Barycenter({v1, v2, v3});
      const Node f3 = Barycenter({v2, v0, v3});
      const Node c0 = Barycenter({v0, v1, v2, v3});

      X = { v0, v1, v2, v3,
            e0, e1, e2, e3, e4, e5,
            f0, f1, f2, f3, c0
          };
    }
};

class Wedge21Projection : public Q2RefinementProjection {
  public:
    Wedge21Projection()
      : Q2RefinementProjection(3,          // dim
                               8,          // nChild
                               6,          // nVert
                               9,          // nEdge
                               5,          // nFace
                               1,          // nCenter
                               "Wedge21") {

      static const unsigned verts[8][6] = {
        {0,  6,  8, 12, 15, 17},
        {6,  1,  7, 15, 13, 16},
        {8,  7,  2, 17, 16, 14},
        {7,  8,  6, 16, 17, 15},
        {12, 15, 17, 3,  9,  11},
        {15, 13, 16, 9,  4,  10},
        {17, 16, 14, 11, 10, 5},
        {16, 17, 15, 10, 11, 9}
      };

      for (unsigned c = 0; c < _nChild; ++c) {
        for (unsigned v = 0; v < _nVert; ++v) {
          _childVerts[c][v] = verts[c][v];
        }
      }

      _fem = std::make_unique<femus::elem_type_3D>("wedge", "biquadratic", "eleventh");

      build();
    }

  private:
    using Node = Q2RefinementProjection::Node;

    Node getParentMin() const override {
      return Node{0.0, 0.0, -1.0};
    }

    Node getParentMax() const override {
      return Node{1.0, 1.0, 1.0};
    }

    bool isBoundaryNode(const std::vector<std::vector<double>> &X_F,
                        unsigned i) const override {
      const double x    = X_F[0][i];
      const double y    = X_F[1][i];
      const double z    = X_F[2][i];
      const double zeta = 1.0 - x - y;
      return (x < 1.0e-4 || y < 1.0e-4 || zeta < 1.0e-4
              || z < -1.0 + 1.0e-4 || z > 1.0 - 1.0e-4);
    }

    void buildChildX(const std::vector<Node> &X_C,
                     std::vector<std::vector<Node>> &childX) const override {
      childX.assign(_nChild, std::vector<Node>());
      for (unsigned c = 0; c < _nChild; ++c) {
        const Node v0 = X_C[_childVerts[c][0]];
        const Node v1 = X_C[_childVerts[c][1]];
        const Node v2 = X_C[_childVerts[c][2]];
        const Node v3 = X_C[_childVerts[c][3]];
        const Node v4 = X_C[_childVerts[c][4]];
        const Node v5 = X_C[_childVerts[c][5]];

        LinearToQuadraticCoordinatesUsingVertices(v0, v1, v2, v3, v4, v5, childX[c]);
      }
    }

    void LinearToQuadraticCoordinates(const Node &X0,
                                      const Node &X1,
                                      std::vector<Node> &X) const override {

      const Node v0 = X0;
      const Node v1 = {X1[0], X0[1], X0[2]};
      const Node v2 = {X0[0], X1[1], X0[2]};

      const Node v3 = {X0[0], X0[1], X1[2]};
      const Node v4 = {X1[0], X0[1], X1[2]};
      const Node v5 = {X0[0], X1[1], X1[2]};

      LinearToQuadraticCoordinatesUsingVertices(v0, v1, v2, v3, v4, v5, X);
    }

    void LinearToQuadraticCoordinatesUsingVertices(const Node &v0, const Node &v1,
                                                   const Node &v2, const Node &v3,
                                                   const Node &v4, const Node &v5,
                                                   std::vector<Node> &X) const {

      auto Barycenter = [this](std::initializer_list<Node> pts) {
        Node x = {0., 0., 0.};
        const double n = static_cast<double>(pts.size());
        for (unsigned d = 0; d < _dim; ++d) {
          double sum = 0.0;
          for (const auto &p : pts) sum += p[d];
          x[d] = sum / n;
        }
        return x;
      };

      const Node e0 = Barycenter({v0, v1});
      const Node e1 = Barycenter({v1, v2});
      const Node e2 = Barycenter({v2, v0});

      const Node e3 = Barycenter({v3, v4});
      const Node e4 = Barycenter({v4, v5});
      const Node e5 = Barycenter({v5, v3});

      const Node e6 = Barycenter({v0, v3});
      const Node e7 = Barycenter({v1, v4});
      const Node e8 = Barycenter({v2, v5});

      const Node f0 = Barycenter({v0, v1, v3, v4});
      const Node f1 = Barycenter({v1, v2, v4, v5});
      const Node f2 = Barycenter({v0, v2, v3, v5});
      const Node f3 = Barycenter({v0, v1, v2});
      const Node f4 = Barycenter({v3, v4, v5});

      const Node c0 = Barycenter({v0, v1, v2, v3, v4, v5});

      X = { v0, v1, v2, v3, v4, v5,
            e0, e1, e2, e3, e4, e5, e6, e7, e8,
            f0, f1, f2, f3, f4, c0
          };
    }
};
