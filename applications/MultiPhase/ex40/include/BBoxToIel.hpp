#pragma once

#include <vector>
#include <cmath>
#include <cfloat>
#include <stdexcept>
#include <iostream>

class BBoxToIel {
  public:
    BBoxToIel() = default;

    BBoxToIel(const MultiLevelMesh& mlMsh,
              const unsigned nPartition,
              const unsigned level) {
      Build(mlMsh, nPartition, level);
    }

    void Build(const MultiLevelMesh& mlMsh,
               const unsigned nPartition,
               const unsigned level) {
      if(nPartition == 0u) {
        throw std::runtime_error("BBoxToIel::Build: nPartition must be > 0");
      }

      const Mesh& msh      = *mlMsh.GetLevel(level);
      const unsigned iproc = msh.processor_id();
      _dim                 = msh.GetDimension();
      const unsigned xType = 2u;

      if(_dim < 1u || _dim > 3u) {
        throw std::runtime_error("BBoxToIel::Build: dim must be 1, 2, or 3");
      }

      auto& xv = msh._topology->_Sol;

      const unsigned offset   = msh._elementOffset[iproc];
      const unsigned offsetp1 = msh._elementOffset[iproc + 1];

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
        GetElementBoundingBox(msh, xv, iel, xType, xMinIel, xMaxIel);

        for(unsigned k = 0; k < _dim; ++k) {
          if(_xMin[k] > xMinIel[k]) _xMin[k] = xMinIel[k];
          if(_xMax[k] < xMaxIel[k]) _xMax[k] = xMaxIel[k];

          const double hIel = xMaxIel[k] - xMinIel[k];
          if(hCellMin[k] > hIel) hCellMin[k] = hIel;
        }
      }

      for(unsigned k = 0; k < _dim; ++k) {
        if(hCellMin[k] <= 0.0) {
          throw std::runtime_error("BBoxToIel::Build: nonpositive local element size");
        }

        _hBox[k] = hCellMin[k] / static_cast<double>(nPartition);
        _xMin[k] -= _hBox[k] / 10.0;
        _xMax[k] += _hBox[k] / 10.0;

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
        eps[k] = 1.e-12 * _hBox[k];
      }

      std::vector<unsigned> iMin(_dim, 0u);
      std::vector<unsigned> iMax(_dim, 0u);

      // Second pass: element -> candidate BBox cells
      for(unsigned iel = offset; iel < offsetp1; ++iel) {
        GetElementBoundingBox(msh, xv, iel, xType, xMinIel, xMaxIel);

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

      for(unsigned ibox = 0; ibox < _bboxToIel.size(); ++ibox) {
        std::cout << "iproc = " << iproc
                  << ", ibox = " << ibox
                  << ", nCandidates = " << _bboxToIel[ibox].size() << std::endl;
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

  private:
    static void GetElementBoundingBox(const Mesh& msh,
                                      const std::vector<NumericVector*>& xv,
                                      const unsigned iel,
                                      const unsigned xType,
                                      std::vector<double>& xMinIel,
                                      std::vector<double>& xMaxIel) {
      const unsigned dim  = msh.GetDimension();
      const unsigned nDof = msh.GetElementDofNumber(iel, xType);

      const unsigned xDof0 = msh.GetSolutionDof(0, iel, xType);
      for(unsigned k = 0; k < dim; ++k) {
        const double x0 = (*xv[k])(xDof0);
        xMinIel[k] = x0;
        xMaxIel[k] = x0;
      }

      for(unsigned i = 1; i < nDof; ++i) {
        const unsigned xDof = msh.GetSolutionDof(i, iel, xType);
        for(unsigned k = 0; k < dim; ++k) {
          const double xval = (*xv[k])(xDof);
          if(xMinIel[k] > xval) xMinIel[k] = xval;
          if(xMaxIel[k] < xval) xMaxIel[k] = xval;
        }
      }
    }

  private:
    unsigned _dim = 0u;

    std::vector<double> _xMin;
    std::vector<double> _xMax;
    std::vector<double> _hBox;
    std::vector<unsigned> _N;

    std::vector<std::vector<unsigned>> _bboxToIel;
};
