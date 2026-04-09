#pragma once

#include <vector>
#include <cmath>
#include <cfloat>
#include <stdexcept>
#include <iostream>

#include "PolynomialBases.hpp"

class BBoxToIel {
  public:

    BBoxToIel(const MultiLevelMesh& mlMsh,
              const unsigned nPartition,
              const unsigned level): _mlMsh(mlMsh), _msh(*mlMsh.GetLevel(level)), _dim(_msh.GetDimension()), _level(level) {

      const unsigned iproc = _msh.processor_id();
      _aPall.clear();
      const unsigned offset   = _msh._elementOffset[iproc];
      const unsigned offsetp1 = _msh._elementOffset[iproc + 1];
      _aPall.resize(offsetp1 - offset);
      _xAll.resize(offsetp1 - offset);
      auto& xv = _msh._topology->_Sol;

      static constexpr unsigned XType = 2;

      for(unsigned iel = offset; iel < offsetp1; ++iel) {
        unsigned ielType = _msh.GetElementType(iel);

        _xAll[iel - offset].resize(_dim);

        const unsigned nDof = _msh.GetElementDofNumber(iel, XType);
        for(unsigned k = 0; k < _dim; ++k) {
          _xAll[iel - offset][k].resize(nDof);
        }
        for(unsigned i = 0; i < nDof; i++) {
          const unsigned xDofi = _msh.GetSolutionDof(i, iel, XType);
          for(unsigned k = 0; k < _dim; ++k) {
            _xAll[iel - offset][k][i] =  (*xv[k])(xDofi);
          }
        }

        static constexpr unsigned NSolType = 3u;
        _aPall[iel - offset].resize(NSolType);
        for (unsigned solType = 0; solType < NSolType; ++solType) {
          _aPall[iel - offset][solType].clear();
          ProjectNodalToPolynomialCoefficients(_aPall[iel - offset][solType], _xAll[iel - offset], ielType, solType);
        }
      }
      Build(nPartition);
    }

    void Build(const unsigned nPartition) {
      if(nPartition == 0u) {
        throw std::runtime_error("BBoxToIel::Build: nPartition must be > 0");
      }

      const unsigned iproc = _msh.processor_id();

      if(_dim < 1u || _dim > 3u) {
        throw std::runtime_error("BBoxToIel::Build: dim must be 1, 2, or 3");
      }

      const unsigned offset   = _msh._elementOffset[iproc];
      const unsigned offsetp1 = _msh._elementOffset[iproc + 1];

      _xMin.assign(_dim,  DBL_MAX);
      _xMax.assign(_dim, -DBL_MAX);
      _hBox.assign(_dim, 0.0);
      _N.assign(_dim, 0u);
      _bboxToIel.clear();

      if(offset == offsetp1) {
        return;
      }

      std::vector<double> hCellMin(_dim, DBL_MAX);
      std::vector<double> xMinIel(_dim, 0.0);
      std::vector<double> xMaxIel(_dim, 0.0);

      // First pass: local partition bounding box and minimum local size
      for(unsigned iel = offset; iel < offsetp1; ++iel) {
        GetElementBoundingBox(_xAll[iel - offset], xMinIel, xMaxIel);

        for(unsigned k = 0; k < _dim; ++k) {
          if(_xMin[k] > xMinIel[k]) _xMin[k] = xMinIel[k];
          if(_xMax[k] < xMaxIel[k]) _xMax[k] = xMaxIel[k];

          const double hIel = xMaxIel[k] - xMinIel[k];
          if(hCellMin[k] > hIel) hCellMin[k] = hIel;
        }
      }

      static constexpr double PaddingFactor = 0.1;
      static constexpr double EpsFactor = 1.e-12;

      for(unsigned k = 0; k < _dim; ++k) {
        if(hCellMin[k] <= 0.0) {
          throw std::runtime_error("BBoxToIel::Build: nonpositive local element size");
        }

        _hBox[k] = hCellMin[k] / static_cast<double>(nPartition);

        _xMin[k] -= _hBox[k] * PaddingFactor;
        _xMax[k] += _hBox[k] * PaddingFactor;

        _N[k] = static_cast<unsigned>(std::ceil((_xMax[k] - _xMin[k]) / _hBox[k]));
        if(_N[k] == 0u) _N[k] = 1u;

        _hBox[k] = (_xMax[k] - _xMin[k]) / static_cast<double>(_N[k]);
      }

      const unsigned nBoxElem = GetNumberOfBoxElements();
      _bboxToIel.clear();
      _bboxToIel.resize(nBoxElem);

      for(unsigned ibox = 0; ibox < nBoxElem; ++ibox) {
        _bboxToIel[ibox].reserve(10);
      }

      std::vector<double> eps(_dim, 0.0);
      for(unsigned k = 0; k < _dim; ++k) {
        eps[k] = EpsFactor * _hBox[k];
      }

      std::vector<unsigned> iMin(_dim, 0u);
      std::vector<unsigned> iMax(_dim, 0u);

      // Second pass: element -> candidate BBox cells
      for(unsigned iel = offset; iel < offsetp1; ++iel) {
        GetElementBoundingBox(_xAll[iel - offset], xMinIel, xMaxIel);

        for(unsigned k = 0; k < _dim; ++k) {
          int ikMin = static_cast<int>(std::floor((xMinIel[k] - _xMin[k] - eps[k]) / _hBox[k]));
          int ikMax = static_cast<int>(std::floor((xMaxIel[k] - _xMin[k] + eps[k]) / _hBox[k]));

          if(ikMin < 0) ikMin = 0;
          if(ikMax < 0) ikMax = 0;
          if(ikMin >= static_cast<int>(_N[k])) ikMin = static_cast<int>(_N[k]) - 1;
          if(ikMax >= static_cast<int>(_N[k])) ikMax = static_cast<int>(_N[k]) - 1;

          if(ikMax < ikMin) std::swap(ikMin, ikMax);

          iMin[k] = static_cast<unsigned>(ikMin);
          iMax[k] = static_cast<unsigned>(ikMax);
        }

        if(_dim == 1u) {
          for(unsigned i = iMin[0]; i <= iMax[0]; ++i) {
            const unsigned ibox = Flatten1D(i);
            _bboxToIel[ibox].push_back(iel);
          }
        }
        else if(_dim == 2u) {
          for(unsigned j = iMin[1]; j <= iMax[1]; ++j) {
            for(unsigned i = iMin[0]; i <= iMax[0]; ++i) {
              const unsigned ibox = Flatten2D(i, j);
              _bboxToIel[ibox].push_back(iel);
            }
          }
        }
        else {
          for(unsigned k = iMin[2]; k <= iMax[2]; ++k) {
            for(unsigned j = iMin[1]; j <= iMax[1]; ++j) {
              for(unsigned i = iMin[0]; i <= iMax[0]; ++i) {
                const unsigned ibox = Flatten3D(i, j, k);
                _bboxToIel[ibox].push_back(iel);
              }
            }
          }
        }
      }



      _xi.resize(_bboxToIel.size());
      std::vector<double> xg(_dim);
      std::vector<unsigned> idx(_dim);
      for(unsigned ibox = 0; ibox < _bboxToIel.size(); ++ibox) {
        if (_bboxToIel[ibox].empty()) {
          _xi[ibox].clear();
          continue;
        }
        Unflatten(ibox, idx);
        for(unsigned k = 0; k < _dim; k++) {
          xg[k] = _xMin[k] + (idx[k] + 0.5) * _hBox[k];
        }

        _xi[ibox].resize(_dim);
        bool ok = false;
        for(unsigned j = 0; j < _bboxToIel[ibox].size(); ++j) {
          const unsigned jel = _bboxToIel[ibox][j];
          short unsigned elType = _msh.GetElementType(jel);
          GetClosestPointInReferenceElement(_xAll[jel - offset], xg, elType, _xi[ibox]);
          ok = GetInverseMapping_fast(2u, elType, _aPall[jel - offset], xg, _xi[ibox], 100u, _phi, _gradPhi);
          if(ok) {
            ok = isInsideReference(elType, _xi[ibox]);
            if(ok) {
              if(j != 0) {
                std::swap( _bboxToIel[ibox][0], _bboxToIel[ibox][j]);
              }
              break;
            }
          }
        }
        if(ok == false) _xi[ibox].clear();

      }

      for(unsigned ibox = 0; ibox < _bboxToIel.size(); ++ibox) {

        Unflatten(ibox, idx);
        for(unsigned k = 0; k < _dim; k++) {
          xg[k] = _xMin[k] + (idx[k] + 0.5) * _hBox[k];
        }

        std::cout << "iproc = " << iproc
                  << ", ibox = " << ibox
                  << ", nCandidates = " << _bboxToIel[ibox].size();
        if(_xi[ibox].size() == _dim) std::cout << ", mesh el = "<<_bboxToIel[ibox][0] << ", barycenter "<< _xi[ibox][0] << " "<<_xi[ibox][1] <<"     "<< xg[0] << " "<<xg[1];;
        std::cout<<std::endl;

      }


    }

    unsigned GetDim() const {
      return _dim;
    }

    const std::vector<double>& GetXMin() const {
      return _xMin;
    }

    const std::vector<double>& GetXMax() const {
      return _xMax;
    }

    const std::vector<double>& GetHBox() const {
      return _hBox;
    }

    const std::vector<unsigned>& GetN() const {
      return _N;
    }

    unsigned GetNumberOfBoxElements() const {
      if(_dim == 0u) return 0u;

      unsigned nBoxElem = 1u;
      for(unsigned k = 0; k < _dim; ++k) {
        nBoxElem *= _N[k];
      }
      return nBoxElem;
    }

    const std::vector<unsigned>& GetCandidates(const unsigned ibox) const {
      if(ibox >= _bboxToIel.size()) {
        throw std::out_of_range("BBoxToIel::GetCandidates: ibox out of range");
      }
      return _bboxToIel[ibox];
    }

    const std::vector<std::vector<unsigned>>& GetAllCandidates() const {
      return _bboxToIel;
    }

  private:
    unsigned Flatten1D(const unsigned i) const {
      return i;
    }

    unsigned Flatten2D(const unsigned i,
                       const unsigned j) const {
      return i + _N[0] * j;
    }

    unsigned Flatten3D(const unsigned i,
                       const unsigned j,
                       const unsigned k) const {
      return i + _N[0] * (j + _N[1] * k);
    }

    void Unflatten(const unsigned ibox, std::vector<unsigned>& idx) const {
      if (_dim == 0u) {
        throw std::runtime_error("BBoxToIel::Unflatten: _dim == 0");
      }

      if (_N.size() != _dim) {
        throw std::runtime_error("BBoxToIel::Unflatten: inconsistent _N size");
      }

      if (ibox >= GetNumberOfBoxElements()) {
        throw std::out_of_range("BBoxToIel::Unflatten: ibox out of range");
      }

      if (idx.size() != _dim) {
        throw std::runtime_error("BBoxToIel::Unflatten: idx has wrong size");
      }

      unsigned q = ibox;
      for (unsigned d = 0; d < _dim; ++d) {
        if (_N[d] == 0u) {
          throw std::runtime_error("BBoxToIel::Unflatten: _N[d] == 0");
        }
        idx[d] = q % _N[d];
        q /= _N[d];
      }
    }

    void GetElementBoundingBox(const std::vector<std::vector<double>>& xv,
                               std::vector<double>& xMinIel,
                               std::vector<double>& xMaxIel) const {
      if (xv.size() != _dim) {
        throw std::runtime_error("BBoxToIel::GetElementBoundingBox: wrong dimension");
      }
      for (unsigned k = 0; k < _dim; ++k) {
        if (xv[k].empty()) {
          throw std::runtime_error("BBoxToIel::GetElementBoundingBox: empty coordinate list");
        }
        xMinIel[k] = xMaxIel[k] = xv[k][0];
        for (unsigned i = 1; i < xv[k].size(); ++i) {
          const double xval = xv[k][i];
          if (xMinIel[k] > xval) xMinIel[k] = xval;
          if (xMaxIel[k] < xval) xMaxIel[k] = xval;
        }
      }
    }


    // ------------------------------
    // inside tests in reference space (as you specified)
    // ------------------------------
    bool isInsideReference(unsigned elType, const std::vector<double>& xi) const {
      const unsigned d = static_cast<unsigned>(xi.size());

      // Box elements: [-1,1]^dim
      if (elType == static_cast<unsigned>(HEX) ||
          elType == static_cast<unsigned>(QUAD) ||
          elType == static_cast<unsigned>(LINE)) {
        for (unsigned k = 0; k < d; ++k) {
          if (xi[k] < -1.0 - 1e-12) return false;
          if (xi[k] >  1.0 + 1e-12) return false;
        }
        return true;
      }

      // Tri: (0,0),(1,0),(0,1)
      if (elType == static_cast<unsigned>(TRI)) {
        if (d != 2) return false;
        const double r = xi[0];
        const double s = xi[1];
        if (r < -1e-12) return false;
        if (s < -1e-12) return false;
        if (r + s > 1.0 + 1e-12) return false;
        return true;
      }

      // Tet: (0,0,0),(1,0,0),(0,1,0),(0,0,1)
      if (elType == static_cast<unsigned>(TET)) {
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
      if (elType == static_cast<unsigned>(WEDGE)) {
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


  private:
    const MultiLevelMesh &_mlMsh;
    const Mesh &_msh;
    const unsigned _dim;
    const unsigned _level;

    std::vector<std::vector<std::vector<std::vector<double>>>> _aPall; // [e][3][...][...][...]
    std::vector<std::vector<std::vector<double>>> _xAll;

    std::vector<double> _xMin;
    std::vector<double> _xMax;
    std::vector<double> _hBox;
    std::vector<unsigned> _N;

    std::vector<std::vector<unsigned>> _bboxToIel;
    std::vector<std::vector<double>> _xi;

    std::vector<double> _phi;
    std::vector<std::vector<double>>_gradPhi;


};
