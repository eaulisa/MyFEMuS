#pragma once

#include <vector>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <algorithm>

#include "PolynomialBases.hpp"
#include "Traits.hpp"
#include "PointLocatorResult.hpp"

#include "Mesh.hpp"

class PointLocator {
  public:

    explicit PointLocator(const Mesh& meshLevel0, double binScale = 0.25)
      : _mesh(meshLevel0), _binScale(binScale) {
      const unsigned d = static_cast<unsigned>(_mesh.dim());
      if (d < 1 || d > 3) throw std::runtime_error("PointLocator: dim must be 1..3");
      if (_binScale <= 0.0) throw std::runtime_error("PointLocator: binScale must be > 0");
    }

    void setBinScale(double binScale) {
      if (binScale <= 0.0) throw std::runtime_error("PointLocator::setBinScale: binScale must be > 0");
      _binScale = binScale;
      _isInitialized = false; // force rebuild of bins/seeds
    }

    // Precompute:
    //  - element geometry caches (_Xall, _aPall, bboxes)
    //  - spatial bins
    //  - per-(bin,element) xiSeed at bin center
    void initializeElementData() {
      const unsigned nEl = static_cast<unsigned>(_mesh.numElements());
      const unsigned d   = static_cast<unsigned>(_mesh.dim());

      // -------------------------
      // global mesh scale + bbox tol
      // -------------------------
      _meshScale = 1.0;
      const std::size_t nNodes = _mesh.numNodes();
      if (nNodes > 0) {
        for (unsigned k = 0; k < d; ++k) {
          for (std::size_t gn = 0; gn < nNodes; ++gn) {
            _meshScale = std::max(_meshScale, std::abs(_mesh.X()[k][gn]));
          }
        }
      }
      _bboxTol = 1e-12 * std::max(1.0, _meshScale) + 1e-14;

      // -------------------------
      // element caches
      // -------------------------
      _Xall.clear();
      _aPall.clear();
      _bboxMin.clear();
      _bboxMax.clear();

      _Xall.resize(nEl);
      _aPall.resize(nEl);
      _bboxMin.resize(nEl, std::vector<double>(d, 0.0));
      _bboxMax.resize(nEl, std::vector<double>(d, 0.0));

      // For bin sizing (avg element bbox edge length)
      std::vector<double> sumEdge(d, 0.0);

      // For global domain bbox of element bboxes
      _gridMin.assign(d, 0.0);
      _gridMax.assign(d, 0.0);
      bool first = true;

      for (unsigned e = 0; e < nEl; ++e) {
        const unsigned et = _mesh.elType()[e];
        const ElemTraits& tr = getTraits(et);
        const unsigned nloc  = tr.connSize;

        const std::vector<unsigned>& conn = _mesh.elTplgy()[e];
        if (conn.size() != nloc) {
          throw std::runtime_error("PointLocator::initializeElementData: conn.size() != traits.connSize");
        }

        // build Xall[e] as dim x nloc
        _Xall[e].assign(d, std::vector<double>(nloc, 0.0));
        for (unsigned li = 0; li < nloc; ++li) {
          const unsigned gn = conn[li];
          for (unsigned k = 0; k < d; ++k) {
            _Xall[e][k][li] = _mesh.X()[k][gn];
          }
        }

        // bbox from all local nodes
        for (unsigned k = 0; k < d; ++k) {
          double mn = _Xall[e][k][0];
          double mx = _Xall[e][k][0];
          for (unsigned li = 1; li < nloc; ++li) {
            mn = std::min(mn, _Xall[e][k][li]);
            mx = std::max(mx, _Xall[e][k][li]);
          }
          _bboxMin[e][k] = mn;
          _bboxMax[e][k] = mx;

          sumEdge[k] += (mx - mn);

          if (first) {
            _gridMin[k] = mn;
            _gridMax[k] = mx;
          }
          else {
            _gridMin[k] = std::min(_gridMin[k], mn);
            _gridMax[k] = std::max(_gridMax[k], mx);
          }
        }
        first = false;

        // allocate 3 orders for this element
        _aPall[e].resize(3);

        // build all orders 0,1,2 from the SAME Q2 geometry nodes
        for (unsigned ord = 0; ord <= 2; ++ord) {
          _aPall[e][ord].clear();
          femus::ProjectNodalToPolynomialCoefficients(_aPall[e][ord], _Xall[e], et, ord);
        }
      }

      // If no elements, done
      if (nEl == 0) {
        _binN.clear();
        _binH.clear();
        _binElems.clear();
        _binXi.clear();
        _isInitialized = true;
        return;
      }

      // -------------------------
      // build bins (small vs element size)
      // -------------------------
      _binN.assign(d, 1u);
      _binH.assign(d, 1.0);

      for (unsigned k = 0; k < d; ++k) {
        const double avgEdge = sumEdge[k] / std::max(1u, nEl);
        const double L = _gridMax[k] - _gridMin[k];

        double hk = _binScale * avgEdge;
        if (hk <= 0.0) hk = (L > 0.0 ? L : 1.0);

        // If L is 0 in some dimension (degenerate), keep N=1 and h=1
        if (L <= 0.0) {
          _binN[k] = 1u;
          _binH[k] = 1.0;
        }
        else {
          unsigned Nk = static_cast<unsigned>(std::ceil(L / hk));
          if (Nk < 1u) Nk = 1u;
          _binN[k] = Nk;
          _binH[k] = L / static_cast<double>(Nk);
        }
      }

      const unsigned nBins = numBins();
      _binElems.assign(nBins, std::vector<unsigned> {});
      _binXi.assign(nBins, std::vector<std::vector<double>> {});

      // -------------------------
      // fill binElems by element bbox overlap
      // -------------------------
      for (unsigned e = 0; e < nEl; ++e) {
        std::vector<unsigned> i0(d, 0u), i1(d, 0u);
        bboxToBinRange(e, i0, i1);

        if (d == 1) {
          for (unsigned ix = i0[0]; ix <= i1[0]; ++ix) {
            _binElems[binId(ix, 0u, 0u)].push_back(e);
          }
        }
        else if (d == 2) {
          for (unsigned iy = i0[1]; iy <= i1[1]; ++iy) {
            for (unsigned ix = i0[0]; ix <= i1[0]; ++ix) {
              _binElems[binId(ix, iy, 0u)].push_back(e);
            }
          }
        }
        else {   // d==3
          for (unsigned iz = i0[2]; iz <= i1[2]; ++iz) {
            for (unsigned iy = i0[1]; iy <= i1[1]; ++iy) {
              for (unsigned ix = i0[0]; ix <= i1[0]; ++ix) {
                _binElems[binId(ix, iy, iz)].push_back(e);
              }
            }
          }
        }
      }

      // -------------------------
      // precompute xi seed for each (bin, element) using bin center
      // and reorder elements so that the element containing the bin center comes first
      // -------------------------
      for (unsigned id = 0; id < nBins; ++id) {
        if (_binElems[id].empty()) continue;

        std::vector<unsigned> ijk(3, 0u);
        decodeBinId(id, ijk);

        std::vector<double> xc(d, 0.0);
        for (unsigned k = 0; k < d; ++k) {
          xc[k] = _gridMin[k] + (static_cast<double>(ijk[k]) + 0.5) * _binH[k];
        }

        // We'll build a reordered list:
        // first those elements for which xc is inside (reference inside),
        // then the rest. We keep ALL elements in the bin.
        std::vector<unsigned> elems_in;
        std::vector<std::vector<double>> xi_in;

        std::vector<unsigned> elems_out;
        std::vector<std::vector<double>> xi_out;

        elems_in.reserve(_binElems[id].size());
        xi_in.reserve(_binElems[id].size());
        elems_out.reserve(_binElems[id].size());
        xi_out.reserve(_binElems[id].size());

        for (unsigned j = 0; j < _binElems[id].size(); ++j) {
          const unsigned e = _binElems[id][j];

          const unsigned et_u = _mesh.elType()[e];
          unsigned short et_su = static_cast<unsigned short>(et_u);

          // compute xi for the bin center in this element
          std::vector<double> xi(d, 0.0);
          femus::GetClosestPointInReferenceElement(_Xall[e], xc, et_su, xi);


          const bool ok = femus::GetInverseMapping(2u, et_su, _aPall[e], xc, xi, 100u);

          // If inverse fails, we still keep the element in the bin list,
          // but mark "no precomputed seed" by storing an EMPTY xi.
          if (!ok) {
            xi.clear();              // <-- empty means invalid seed
            elems_out.push_back(e);
            xi_out.push_back(xi);
            continue;
          }

          // If xc maps inside this element, prioritize it
          if (isInsideReference(et_u, xi)) {
            elems_in.push_back(e);
            xi_in.push_back(xi);
          }
          else {
            elems_out.push_back(e);
            xi_out.push_back(xi);
          }
        }

        // Merge: inside-first
        std::vector<unsigned> newElems;
        std::vector<std::vector<double>> newXi;

        newElems.reserve(elems_in.size() + elems_out.size());
        newXi.reserve(xi_in.size() + xi_out.size());

        for (unsigned i = 0; i < elems_in.size(); ++i) {
          newElems.push_back(elems_in[i]);
          newXi.push_back(xi_in[i]);
        }
        for (unsigned i = 0; i < elems_out.size(); ++i) {
          newElems.push_back(elems_out[i]);
          newXi.push_back(xi_out[i]);
        }

        _binElems[id].swap(newElems);
        _binXi[id].swap(newXi);
      }

      _isInitialized = true;
    }

    void locateAll(std::vector<PointLocatorResult>& out,
                   const std::vector<std::vector<double>>& points) {
      if (!_isInitialized) {
        initializeElementData();
      }

      const unsigned d = static_cast<unsigned>(_mesh.dim());
      if (points.size() != d) {
        throw std::runtime_error("PointLocator::locateAll: point has wrong dimension");
      }

      const std::size_t nPts = points[0].size();
      for (unsigned k = 1; k < d; ++k) {
        if (points[k].size() != nPts) {
          throw std::runtime_error(
            "PointLocator::locateAll: inconsistent point array sizes"
          );
        }
      }

      out.resize(nPts);
      std::vector<double> Xp(d);

      for (std::size_t i = 0; i < nPts; ++i) {
        for (unsigned k = 0; k < d; ++k) {
          Xp[k] = points[k][i];
        }
        out[i] = locateOne(Xp);
      }
    }


    double meshScale() const {
      return _meshScale;
    }
    double bboxTol()   const {
      return _bboxTol;
    }

  private:
    const Mesh& _mesh;

    // Per-element cached data (coarse mesh only)
    std::vector<std::vector<std::vector<double>>> _Xall; // [e][dim][nloc]
    std::vector<std::vector<std::vector<std::vector<double>>>> _aPall; // [e][3][...][...][...]
    std::vector<std::vector<double>> _bboxMin; // [e][dim]
    std::vector<std::vector<double>> _bboxMax; // [e][dim]

    double _meshScale = 1.0;
    double _bboxTol   = 1e-12;

    bool _isInitialized = false;

    // Bin data
    double _binScale = 0.25;
    std::vector<double> _gridMin;           // [dim]
    std::vector<double> _gridMax;           // [dim]
    std::vector<unsigned> _binN;            // [dim] (for dim<3, unused entries absent)
    std::vector<double> _binH;              // [dim]
    std::vector<std::vector<unsigned>> _binElems;                 // [binId] -> list of element ids
    std::vector<std::vector<std::vector<double>>> _binXi;         // [binId] -> list of xi vectors parallel to _binElems

    // ------------------------------
    // bbox utilities
    // ------------------------------
    bool pointInBBox(unsigned e, const std::vector<double>& xp) const {
      const unsigned d = static_cast<unsigned>(xp.size());
      for (unsigned k = 0; k < d; ++k) {
        if (xp[k] < _bboxMin[e][k] - _bboxTol) return false;
        if (xp[k] > _bboxMax[e][k] + _bboxTol) return false;
      }
      return true;
    }

    // ------------------------------
    // bin helpers
    // ------------------------------
    unsigned numBins() const {
      const unsigned d = static_cast<unsigned>(_mesh.dim());
      if (d == 1) return _binN[0];
      if (d == 2) return _binN[0] * _binN[1];
      return _binN[0] * _binN[1] * _binN[2];
    }

    unsigned clampIndex(int i, unsigned N) const {
      if (N == 0u) return 0u;
      if (i < 0) return 0u;
      if (static_cast<unsigned>(i) >= N) return N - 1u;
      return static_cast<unsigned>(i);
    }

    unsigned binId(unsigned ix, unsigned iy, unsigned iz) const {
      const unsigned d = static_cast<unsigned>(_mesh.dim());
      if (d == 1) return ix;
      if (d == 2) return ix + _binN[0] * iy;
      return ix + _binN[0] * (iy + _binN[1] * iz);
    }

    void decodeBinId(unsigned id, std::vector<unsigned>& ijk) const {
      const unsigned d = static_cast<unsigned>(_mesh.dim());
      ijk.assign(3, 0u);

      if (d == 1) {
        ijk[0] = id;
        return;
      }
      if (d == 2) {
        ijk[1] = id / _binN[0];
        ijk[0] = id % _binN[0];
        return;
      }
      // d==3
      const unsigned Nx = _binN[0];
      const unsigned Ny = _binN[1];
      const unsigned Nxy = Nx * Ny;

      ijk[2] = id / Nxy;
      const unsigned rem = id % Nxy;
      ijk[1] = rem / Nx;
      ijk[0] = rem % Nx;
    }

    void bboxToBinRange(unsigned e, std::vector<unsigned>& i0, std::vector<unsigned>& i1) const {
      const unsigned d = static_cast<unsigned>(_mesh.dim());
      i0.assign(d, 0u);
      i1.assign(d, 0u);

      for (unsigned k = 0; k < d; ++k) {
        const double L = _gridMax[k] - _gridMin[k];
        if (L <= 0.0) {
          i0[k] = 0u;
          i1[k] = 0u;
          continue;
        }
        const double invH = 1.0 / _binH[k];

        const double a = (_bboxMin[e][k] - _gridMin[k]) * invH;
        const double b = (_bboxMax[e][k] - _gridMin[k]) * invH;

        const int ia = static_cast<int>(std::floor(a));
        const int ib = static_cast<int>(std::floor(b));

        i0[k] = clampIndex(ia, _binN[k]);
        i1[k] = clampIndex(ib, _binN[k]);

        if (i1[k] < i0[k]) std::swap(i0[k], i1[k]);
      }
    }

    bool pointToBin(const std::vector<double>& xp, unsigned& outId) const {
      const unsigned d = static_cast<unsigned>(_mesh.dim());
      std::vector<unsigned> idx(3, 0u);

      for (unsigned k = 0; k < d; ++k) {
        const double L = _gridMax[k] - _gridMin[k];
        if (L <= 0.0) {
          idx[k] = 0u;
          continue;
        }
        // outside domain bbox => outside
        if (xp[k] < _gridMin[k] - _bboxTol) return false;
        if (xp[k] > _gridMax[k] + _bboxTol) return false;

        const double t = (xp[k] - _gridMin[k]) / _binH[k];
        int it = static_cast<int>(std::floor(t));
        idx[k] = clampIndex(it, _binN[k]);
      }

      outId = binId(idx[0], idx[1], idx[2]);
      return true;
    }

    // ------------------------------
    // per-point locate using (bin -> elements) and (bin,element)->xiSeed
    // ------------------------------
    PointLocatorResult locateOne(const std::vector<double>& xp) const {
      PointLocatorResult out;

      const unsigned nEl = static_cast<unsigned>(_mesh.numElements());
      if (nEl == 0) return out;

      unsigned id = 0u;
      if (!pointToBin(xp, id)) {
        return out; // outside domain bbox
      }

      if (id >= _binElems.size()) return out;
      if (_binElems[id].empty()) return out;

      const unsigned d = static_cast<unsigned>(_mesh.dim());

      // Try each candidate element in this bin, starting from precomputed xiSeed
      for (unsigned j = 0; j < _binElems[id].size(); ++j) {
        const unsigned e = _binElems[id][j];

        // safety: bbox prune (cheap)
        if (!pointInBBox(e, xp)) continue;

        const unsigned et_u = _mesh.elType()[e];
        unsigned short et_su = static_cast<unsigned short>(et_u);


        std::vector<double> xi = _binXi[id][j]; // seed at bin center; may be EMPTY

        // If no precomputed seed, build a seed for this point now.
        if (xi.empty()) {
          xi.assign(d, 0.0);
          femus::GetClosestPointInReferenceElement(_Xall[e], xp, et_su, xi);
        }

        const bool ok = femus::GetInverseMapping(2u, et_su, _aPall[e], xp, xi, 50u);
        if (!ok) continue;
        if (!isInsideReference(et_u, xi)) continue;

        out.elem = e;
        out.xi   = xi;
        out.ok   = true;
        return out;
      }

      return out;
    }

    // ------------------------------
    // inside tests in reference space (as you specified)
    // ------------------------------
    bool isInsideReference(unsigned elType, const std::vector<double>& xi) const {
      const unsigned d = static_cast<unsigned>(xi.size());

      // Box elements: [-1,1]^dim
      if (elType == static_cast<unsigned>(Hex27) ||
          elType == static_cast<unsigned>(Quad9) ||
          elType == static_cast<unsigned>(Line3)) {
        for (unsigned k = 0; k < d; ++k) {
          if (xi[k] < -1.0 - 1e-12) return false;
          if (xi[k] >  1.0 + 1e-12) return false;
        }
        return true;
      }

      // Tri: (0,0),(1,0),(0,1)
      if (elType == static_cast<unsigned>(Tri7)) {
        if (d != 2) return false;
        const double r = xi[0];
        const double s = xi[1];
        if (r < -1e-12) return false;
        if (s < -1e-12) return false;
        if (r + s > 1.0 + 1e-12) return false;
        return true;
      }

      // Tet: (0,0,0),(1,0,0),(0,1,0),(0,0,1)
      if (elType == static_cast<unsigned>(Tet15)) {
        if (d != 3) return false;
        const double r = xi[0];
        const double s = xi[1];
        const double t = xi[2];
        if (r < -1e-12) return false;
        if (s < -1e-12) return false;
        if (t < -1e-12) return false;
        if (r + s + t > 1.0 + 1e-12) return false;
        return true;
      }

      // Wedge: (r,s) tri with r>=0,s>=0,r+s<=1 and z in [-1,1]
      if (elType == static_cast<unsigned>(Wedge21)) {
        if (d != 3) return false;
        const double r = xi[0];
        const double s = xi[1];
        const double z = xi[2];
        if (r < -1e-12) return false;
        if (s < -1e-12) return false;
        if (r + s > 1.0 + 1e-12) return false;
        if (z < -1.0 - 1e-12) return false;
        if (z >  1.0 + 1e-12) return false;
        return true;
      }

      return false;
    }
};
