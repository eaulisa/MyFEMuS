#pragma once

#include <vector>
#include <cmath>
#include <cfloat>
#include <stdexcept>
#include <iostream>

#include "PolynomialBases.hpp"

class BBoxToIel {
  public:

    BBoxToIel(MultiLevelMesh& mlMsh,
              const unsigned level,
              const unsigned nPartition
             ): _msh(mlMsh.GetLevel(level)), _dim(_msh->GetDimension()), _level(level), _iproc(_msh->processor_id()), _nprocs(_msh->n_processors()) {


      _aPall.clear();
      const unsigned offset   = _msh->_elementOffset[_iproc];
      const unsigned offsetp1 = _msh->_elementOffset[_iproc + 1];
      _aPall.resize(offsetp1 - offset);
      _xAll.resize(offsetp1 - offset);
      auto& xv = _msh->_topology->_Sol;

      static constexpr unsigned XType = 2;

      for(unsigned iel = offset; iel < offsetp1; ++iel) {
        unsigned ielType = _msh->GetElementType(iel);

        _xAll[iel - offset].resize(_dim);

        const unsigned nDof = _msh->GetElementDofNumber(iel, XType);
        for(unsigned k = 0; k < _dim; ++k) {
          _xAll[iel - offset][k].resize(nDof);
        }
        for(unsigned i = 0; i < nDof; i++) {
          const unsigned xDofi = _msh->GetSolutionDof(i, iel, XType);
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

      if(_dim < 1u || _dim > 3u) {
        throw std::runtime_error("BBoxToIel::Build: dim must be 1, 2, or 3");
      }

      const unsigned offset   = _msh->_elementOffset[_iproc];
      const unsigned offsetp1 = _msh->_elementOffset[_iproc + 1];

      _xMinMemory.resize(_nprocs * _dim);
      _xMaxMemory.resize(_nprocs * _dim);
      _xMin.resize(_nprocs);
      _xMax.resize(_nprocs);
      for(unsigned kproc = 0; kproc < _nprocs; ++kproc) {
        _xMin[kproc] = _xMinMemory.data() + kproc * _dim;
        _xMax[kproc] = _xMaxMemory.data() + kproc * _dim;
      }

      std::fill(_xMin[_iproc], _xMin[_iproc] + _dim, DBL_MAX);
      std::fill(_xMax[_iproc], _xMax[_iproc] + _dim, -DBL_MAX);
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
          if(_xMin[_iproc][k] > xMinIel[k]) _xMin[_iproc][k] = xMinIel[k];
          if(_xMax[_iproc][k] < xMaxIel[k]) _xMax[_iproc][k] = xMaxIel[k];

          const double hIel = xMaxIel[k] - xMinIel[k];
          if(hCellMin[k] > hIel) hCellMin[k] = hIel;
        }
      }

      MPI_Allgather(_xMin[_iproc], _dim, MPI_DOUBLE,
                    _xMinMemory.data(), _dim, MPI_DOUBLE,
                    MPI_COMM_WORLD);

      MPI_Allgather(_xMax[_iproc], _dim, MPI_DOUBLE,
                    _xMaxMemory.data(), _dim, MPI_DOUBLE,
                    MPI_COMM_WORLD);

      static constexpr double PaddingFactor = 0.1;
      static constexpr double EpsFactor = 1.e-12;

      for(unsigned k = 0; k < _dim; ++k) {
        if(hCellMin[k] <= 0.0) {
          throw std::runtime_error("BBoxToIel::Build: nonpositive local element size");
        }

        _hBox[k] = hCellMin[k] / static_cast<double>(nPartition);

        _xMin[_iproc][k] -= _hBox[k] * PaddingFactor;
        _xMax[_iproc][k] += _hBox[k] * PaddingFactor;

        _N[k] = static_cast<unsigned>(std::ceil((_xMax[_iproc][k] - _xMin[_iproc][k]) / _hBox[k]));
        if(_N[k] == 0u) _N[k] = 1u;

        _hBox[k] = (_xMax[_iproc][k] - _xMin[_iproc][k]) / static_cast<double>(_N[k]);
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
          int ikMin = static_cast<int>(std::floor((xMinIel[k] - _xMin[_iproc][k] - eps[k]) / _hBox[k]));
          int ikMax = static_cast<int>(std::floor((xMaxIel[k] - _xMin[_iproc][k] + eps[k]) / _hBox[k]));

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
          xg[k] = _xMin[_iproc][k] + (idx[k] + 0.5) * _hBox[k];
        }

        _xi[ibox].resize(_dim);
        bool ok = false;
        for(unsigned j = 0; j < _bboxToIel[ibox].size(); ++j) {
          const unsigned jel = _bboxToIel[ibox][j];
          short unsigned elType = _msh->GetElementType(jel);

          assert(jel >= offset);
          assert(jel - offset < _xAll.size());
          assert(jel - offset < _aPall.size());

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
          xg[k] = _xMin[_iproc][k] + (idx[k] + 0.5) * _hBox[k];
        }
      }
    }

    unsigned GetDim() const {
      return _dim;
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

    void GetInverseMappingOnCoarseLevel(std::vector<MyVector<double>>&X, LevelMarkers &lX, LevelMarkers &lY);
    void Project(const MultiLevelMesh &_mlMsh, LevelMarkers &lX, LevelMarkers &lY);

    bool inside_tet(double u, double v, double w, double eps) const {
      return (u >= -eps) && (v >= -eps) && (w >= -eps) && (u + v + w <= 1.0 + eps);
    };

    bool inside_tri(double r, double s, double eps)const {
      return (r >= -eps) && (s >= -eps) && (r + s <= 1.0 + eps);
    };

    bool inside_line(double t, double eps)const {
      return (t >= -1.0 - eps) && (t <= 1.0 + eps);
    };

    // small helpers: avoid xi = {..} (which tends to show up as vector::operator=)
    void set1(std::vector<double>& xi, double a) const {
      //xi.resize(1);
      xi[0] = a;
    };
    void set2(std::vector<double>& xi, double a, double b) const {
      //xi.resize(2);
      xi[0] = a;
      xi[1] = b;
    };
    void set3(std::vector<double>& xi, double a, double b, double c) const {
      //xi.resize(3);
      xi[0] = a;
      xi[1] = b;
      xi[2] = c;
    };

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

    void GetChildIndexAndLocalXi(const unsigned et, const std::vector<double>& xi0, unsigned & childIndex, std::vector<double>& xi1, const double eps = 1.0e-12) const;

    const unsigned& GetLevel() const {
      return _level;
    }

    void SetMesh(Mesh *msh) {
      _msh = msh;
    }


  private:
    unsigned Flatten(const std::vector<unsigned> &idx) const {
      if(_dim == 1) return Flatten1D(idx[0]);
      if(_dim == 2) return Flatten2D(idx[0], idx[1]);
      if(_dim == 3) return Flatten3D(idx[0], idx[1], idx[2]);
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
    Mesh *_msh;
    const unsigned _dim;
    const unsigned _level;
    const unsigned _iproc;
    const unsigned _nprocs;

    std::vector<std::vector<std::vector<std::vector<double>>>> _aPall; // [e][3][...][...][...]
    std::vector<std::vector<std::vector<double>>> _xAll;


    std::vector<double> _xMinMemory;
    std::vector<double*> _xMin;

    std::vector<double> _xMaxMemory;
    std::vector<double*> _xMax;

    std::vector<double> _hBox;
    std::vector<unsigned> _N;

    std::vector<std::vector<unsigned>> _bboxToIel;
    std::vector<std::vector<double>> _xi;

    std::vector<double> _phi;
    std::vector<std::vector<double>>_gradPhi;


};

void BBoxToIel::GetInverseMappingOnCoarseLevel(std::vector<MyVector<double>>&X, LevelMarkers &lX, LevelMarkers &lY) {


  if(X.size() != _dim) {
    throw std::runtime_error("GetInverseMappingOnCoarseLevel: X has wrong dimension");
  }

  const unsigned nPoints = X[0].size();
  lX.SetFieldLocalSize(nPoints);

  for(unsigned k = 1; k < _dim; ++k) {
    if(X[k].size() != nPoints || X[k].begin() != X[0].begin() || X[k].end() != X[0].end()) {
      throw std::runtime_error("GetInverseMappingOnCoarseLevel: X[k] indexing is inconsistent");
    }
  }

  std::vector<std::vector<double>> Z_s(_nprocs);
  std::vector<std::vector<unsigned>> &map_s = lX.GetMap_s();
  map_s.resize(_nprocs);


  std::vector<std::vector<unsigned>> ZisInside_s(_nprocs);
  std::vector<std::vector<double>> Z_r(_nprocs);
  std::vector<std::vector<unsigned>> ZisInside_r;// = lY.GetPointIsInside_r();
  ZisInside_r.resize(_nprocs);
  std::vector<std::vector<unsigned>> &map_r = lY.GetMap_r();
  map_r.resize(_nprocs);

  std::vector<std::vector<double>> field_r(_nprocs);

  std::vector<unsigned> size_s(_nprocs, 0);
  std::vector<unsigned> size_r(_nprocs, 0);

// Build packed send buffers:
// Z_s[j] = [x0,y0,(z0), x1,y1,(z1), ...] for points sent to process j
  for(unsigned j = 0; j < _nprocs; ++j) {
    Z_s[j].reserve(2 * _dim * nPoints / _nprocs);
  }

  for(unsigned j = 0; j < _nprocs; ++j) {
    map_s[j].reserve(2 * nPoints / _nprocs);
  }

  for(unsigned i = X[0].begin(); i < X[0].end(); i++) {
    for(unsigned j = 0; j < _nprocs; j++) {
      bool insideBB = true;
      for(unsigned k = 0; k < _dim; k++) {
        if(X[k][i] < _xMin[j][k] || X[k][i] > _xMax[j][k]) {
          insideBB = false;
          break;
        }
      }
      if(insideBB) {
        for(unsigned k = 0; k < _dim; k++) {
          Z_s[j].push_back(X[k][i]);
        }
        map_s[j].push_back(i - X[0].begin());
      }
    }
  }
  for(unsigned j = 0; j < _nprocs; j++) ZisInside_s[j].resize(map_s[j].size());

  for(unsigned j = 0; j < _nprocs; ++j) {
    assert(Z_s[j].size() % _dim == 0);
  }

// Number of doubles sent to each process
  for(unsigned j = 0; j < _nprocs; ++j) {
    size_s[j] = Z_s[j].size() / _dim;
  }

// Exchange sizes
  MPI_Alltoall(size_s.data(), 1, MPI_UNSIGNED,
               size_r.data(), 1, MPI_UNSIGNED,
               MPI_COMM_WORLD);

  const unsigned uint_max = std::numeric_limits<unsigned>::max();

// Resize receive buffers
  for(unsigned j = 0; j < _nprocs; ++j) {
    Z_r[j].resize(size_r[j] * _dim);
    ZisInside_r[j].assign(size_r[j], 0u);
    map_r[j].assign(size_r[j], uint_max);
  }

// Nonblocking receives first
  std::vector<MPI_Request> reqs;
  reqs.reserve(4 * _nprocs);

  for(int j = 0; j < static_cast<int>(_nprocs); ++j) {
    if(j == static_cast<int>(_iproc)) {
      Z_r[_iproc] = Z_s[_iproc];
    }
    else {
      if(Z_r[j].size() > 0) {
        MPI_Request req;
        MPI_Irecv(Z_r[j].data(), static_cast<int>(Z_r[j].size()), MPI_DOUBLE,
                  j, 100, MPI_COMM_WORLD, &req);
        reqs.push_back(req);
      }
    }
  }

// Then nonblocking sends
  for(int j = 0; j < static_cast<int>(_nprocs); ++j) {
    if(j != _iproc) {
      if(Z_s[j].size() > 0) {
        MPI_Request req;
        MPI_Isend(Z_s[j].data(), static_cast<int>(Z_s[j].size()), MPI_DOUBLE,
                  j, 100, MPI_COMM_WORLD, &req);
        reqs.push_back(req);
      }
    }
  }

// Complete communication
  if(!reqs.empty()) {
    MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
  }

  unsigned YLocalSize = 0;
  for(unsigned j = 0; j < _nprocs; ++j) {
    YLocalSize += size_r[j];
  }


  //std::vector < std::vector<double >> YFieldLocal(nFields);
  std::vector < std::vector<double >> YiLocal(_dim);
  std::vector<unsigned > YIelLocal;
  for( unsigned k = 0; k < _dim; k++) {
    YiLocal[k].reserve(YLocalSize);
  }

  YIelLocal.reserve(YLocalSize);

  const unsigned offset   = _msh->_elementOffset[_iproc];

  std::vector<unsigned> idx(_dim);
  std::vector<double> x(_dim);
  std::vector<double> xi(_dim);

  for(unsigned jproc = 0; jproc < _nprocs; ++jproc) {
    if(Z_r[jproc].size() % _dim != 0) {
      throw std::runtime_error("GetInverseMappingOnCoarseLevel: invalid packed receive size");
    }
    unsigned Npoints = Z_r[jproc].size() / _dim;
    for(unsigned i = 0; i < Npoints; i++) {
      for( unsigned k = 0; k < _dim; k++) {
        x[k] = Z_r[jproc][i * _dim + k];

        int ik = static_cast<int>(std::floor((x[k] - _xMin[_iproc][k]) / _hBox[k]));
        if(ik < 0) ik = 0;
        if(ik >= static_cast<int>(_N[k])) ik = static_cast<int>(_N[k]) - 1;
        idx[k] = static_cast<unsigned>(ik);

      }
      unsigned ibox = Flatten(idx);
      for(unsigned j = 0; j < _bboxToIel[ibox].size(); ++j) {
        const unsigned jel = _bboxToIel[ibox][j];
        short unsigned elType = _msh->GetElementType(jel);

        assert(jel >= offset);
        assert(jel - offset < _xAll.size());
        assert(jel - offset < _aPall.size());

        if( j != 0 || _xi[ibox].size() == 0) {
          GetClosestPointInReferenceElement(_xAll[jel - offset], x, elType, xi);
        }
        else xi = _xi[ibox];
        bool ok = GetInverseMapping_fast(2u, elType, _aPall[jel - offset], x, xi, 100u, _phi, _gradPhi);
        if(ok) {
          ok = isInsideReference(elType, xi);
          if(ok) {
            for( unsigned k = 0; k < _dim; k++) {
              YiLocal[k].push_back(xi[k]);
            }
            YIelLocal.push_back(jel);
            ZisInside_r[jproc][i] = 1u;
            map_r[jproc][i] = YIelLocal.size() - 1u;
            break;
          }
        }
      }
    }
  }

  reqs.clear();
  for(int j = 0; j < static_cast<int>(_nprocs); ++j) {
    if(j == static_cast<int>(_iproc)) {
      ZisInside_s[_iproc] = ZisInside_r[_iproc];
    }
    else {
      if(ZisInside_s[j].size() > 0) {
        MPI_Request req;
        MPI_Irecv(ZisInside_s[j].data(), static_cast<int>(ZisInside_s[j].size()), MPI_UNSIGNED,
                  j, 100, MPI_COMM_WORLD, &req);
        reqs.push_back(req);
      }
    }
  }

// Then nonblocking sends
  for(int j = 0; j < static_cast<int>(_nprocs); ++j) {
    if(j != _iproc && ZisInside_r[j].size() > 0) {
      MPI_Request req;
      MPI_Isend(ZisInside_r[j].data(), static_cast<int>(ZisInside_r[j].size()), MPI_UNSIGNED,
                j, 100, MPI_COMM_WORLD, &req);
      reqs.push_back(req);
    }
  }
// Complete communication
  if(!reqs.empty()) {
    MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
  }

  std::vector <bool> &ZisInside = lX.GetPointInsideDomain();
  ZisInside.assign(nPoints, false);
  for(unsigned kproc = 0; kproc < _nprocs; ++kproc) {
    for(unsigned i = 0; i < ZisInside_s[kproc].size(); ++i) {
      unsigned idx = map_s[kproc][i];
      ZisInside[idx] = ZisInside[idx] || static_cast<bool>(ZisInside_s[kproc][i]);
    }
  }



  for(unsigned kproc = 0; kproc < _nprocs; ++kproc) {
    for(unsigned i = 0; i < ZisInside_s[kproc].size(); ++i) {
      if(ZisInside_s[kproc][i] == 0) map_s[kproc][i] = uint_max;
    }
  }

  std::vector<MyVector<double>>&YField = lY.GetFields();
  std::vector<MyVector<double>>&Yi = lY.GetLocalCoordinates();
  MyVector<unsigned>&YIel = lY.GetElements();
  lY.SetLevel(_level);
  lY.SetFieldLocalSize(YIelLocal.size());

  Yi.resize(_dim);
  for( unsigned k = 0; k < _dim; k++) {
    Yi[k].buildFromLocal(YiLocal[k]);
  }
  YIel.buildFromLocal(YIelLocal);
}





void BBoxToIel::Project(const MultiLevelMesh &_mlMsh, LevelMarkers &lX, LevelMarkers &lY) {

  std::vector<MyVector<double>>&Xi = lX.GetLocalCoordinates();
  MyVector<unsigned>&XIel = lX.GetElements();

  const unsigned levelX = lX.GetLevel();
  const unsigned nPoints = XIel.size();


  if (Xi.size() != _dim) {
    throw std::runtime_error("Project: Xi has wrong dimension");
  }

  for (unsigned k = 0; k < _dim; ++k) {
    if (Xi[k].size() != nPoints) {
      throw std::runtime_error("Project: inconsistent point array sizes");
    }
  }

  if (XIel.size() != nPoints) {
    throw std::runtime_error("Project: inconsistent element index array size");
  }

  const Mesh* mshX = _mlMsh.GetLevel(levelX);

  if(levelX + 1 >= _mlMsh.GetNumberOfLevels()) {
    throw std::runtime_error("MultiLevelMesh has not enough levels to projec the points");
  };

  const Mesh* mshY = _mlMsh.GetLevel(levelX + 1);

  std::vector<std::vector<double>> Yi_s(_nprocs);
  std::vector<std::vector<unsigned>> YIel_s(_nprocs);

  std::vector<std::vector<unsigned>> &map_s = lX.GetMap_s();
  map_s.resize(_nprocs);

  for(unsigned j = 0; j < _nprocs; ++j) {

    map_s[j].clear();
    map_s[j].reserve(2 * nPoints / _nprocs);

    Yi_s[j].clear();
    Yi_s[j].reserve(2 * _dim * nPoints / _nprocs);

    YIel_s[j].clear();
    YIel_s[j].reserve(2 * nPoints / _nprocs);

  }

  const unsigned PWCtype = 3u; // PieceWiseConstant fem type enumerator

  std::vector<double> xi(_dim);
  std::vector<double> yi(_dim);
  unsigned childIndex;
  for(unsigned i = Xi[0].begin(); i < Xi[0].end(); i++) {
    for(unsigned k = 0; k < _dim; k++) {
      xi[k] = Xi[k][i];
    }
    unsigned xiel = XIel[i];
    const unsigned ielType = mshX->GetElementType(xiel);
    unsigned yiel = 0;
    if(mshX->GetRefinedElementIndex(xiel)) {
      GetChildIndexAndLocalXi(ielType, xi, childIndex, yi);
      yiel = mshX->el->GetChildElement(xiel, childIndex);
    }
    else {
      yi = xi;
      yiel = mshX->el->GetChildElement(xiel, 0);
    }
    unsigned kproc = mshY->IsdomBisectionSearch(yiel, PWCtype);

    for(unsigned k = 0; k < _dim; k++) {
      Yi_s[kproc].push_back(yi[k]);
    }

    YIel_s[kproc].push_back(yiel);
    map_s[kproc].push_back(i - XIel.begin());
  }

  std::vector<unsigned> size_s(_nprocs);
  std::vector<unsigned> size_r(_nprocs);
  for(unsigned kproc = 0; kproc < _nprocs; kproc++) {
    size_s[kproc] = YIel_s[kproc].size();
  }

  MPI_Alltoall(size_s.data(), 1, MPI_UNSIGNED,
               size_r.data(), 1, MPI_UNSIGNED,
               MPI_COMM_WORLD);

  std::vector<std::vector<double>> Yi_r(_nprocs);
  std::vector<std::vector<unsigned>> YIel_r(_nprocs);
  for(unsigned kproc = 0; kproc < _nprocs; kproc++) {
    Yi_r[kproc].resize(size_r[kproc] * _dim);
    YIel_r[kproc].resize(size_r[kproc]);
  }

  const int tagYi   = 101;
  const int tagYIel = 102;

  std::vector<MPI_Request> requests;
  requests.reserve(4 * _nprocs);

// Self-copy first, no MPI needed
  if (size_s[_iproc] > 0) {
    Yi_r[_iproc]  = Yi_s[_iproc];
    YIel_r[_iproc] = YIel_s[_iproc];
  }

// Post receives
  for (unsigned kproc = 0; kproc < _nprocs; ++kproc) {
    if (kproc == _iproc) continue;

    if (size_r[kproc] > 0) {
      MPI_Request req1;
      MPI_Irecv(Yi_r[kproc].data(),
                static_cast<int>(Yi_r[kproc].size()),
                MPI_DOUBLE,
                static_cast<int>(kproc),
                tagYi,
                MPI_COMM_WORLD,
                &req1);
      requests.push_back(req1);

      MPI_Request req2;
      MPI_Irecv(YIel_r[kproc].data(),
                static_cast<int>(YIel_r[kproc].size()),
                MPI_UNSIGNED,
                static_cast<int>(kproc),
                tagYIel,
                MPI_COMM_WORLD,
                &req2);
      requests.push_back(req2);
    }
  }

// Post sends
  for (unsigned kproc = 0; kproc < _nprocs; ++kproc) {
    if (kproc == _iproc) continue;

    if(size_s[kproc] > 0) {
      MPI_Request req1;
      MPI_Isend(Yi_s[kproc].data(),
                static_cast<int>(Yi_s[kproc].size()),
                MPI_DOUBLE,
                static_cast<int>(kproc),
                tagYi,
                MPI_COMM_WORLD,
                &req1);
      requests.push_back(req1);

      MPI_Request req2;
      MPI_Isend(YIel_s[kproc].data(),
                static_cast<int>(YIel_s[kproc].size()),
                MPI_UNSIGNED,
                static_cast<int>(kproc),
                tagYIel,
                MPI_COMM_WORLD,
                &req2);
      requests.push_back(req2);
    }
  }

  if (!requests.empty()) {
    MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
  }

  std::vector<MyVector<double>>&Yi = lY.GetLocalCoordinates();
  MyVector<unsigned>&YIel = lY.GetElements();
  lY.SetLevel(levelX + 1);

  std::vector<std::vector<unsigned>> &map_r = lY.GetMap_r();
  map_r.resize(_nprocs);
  for(unsigned j = 0; j < _nprocs; ++j) {
    map_r[j].resize(size_r[j]);
  }

  Yi.resize(_dim);

  unsigned Zsize = 0;
  for(unsigned kproc = 0; kproc < _nprocs; kproc++) {
    Zsize += YIel_r[kproc].size();
  }

  std::vector<std::vector<double>> Zi(_dim);
  std::vector<unsigned> ZIel;

  for(unsigned k = 0; k < _dim; k++) {
    Zi[k].reserve(Zsize);
  }

  ZIel.reserve(Zsize);

  for(unsigned kproc = 0; kproc < _nprocs; kproc++) {
    for(unsigned i = 0; i < YIel_r[kproc].size(); i++) {
      for(unsigned k = 0; k < _dim; k++) {
        Zi[k].push_back(Yi_r[kproc][i * _dim + k]);
      }
      ZIel.push_back(YIel_r[kproc][i]);
      map_r[kproc][i] = ZIel.size() - 1u;
    }
  }

  lY.SetFieldLocalSize(ZIel.size());

  for(unsigned k = 0; k < _dim; k++) {
    Yi[k].buildFromLocal(Zi[k]);
  }
  YIel.buildFromLocal(ZIel);


  // //BEGIN TEST
  //
  // std::vector<MyVector<double>>&Xfield = lX.GetFields();
  // unsigned nFields = Xfield.size();
  //
  // std::vector<std::vector<double>> Wfield_r, Wfield_s;
  //
  // LevelMarkers lXc = lX;
  // bool backward = true;
  //
  // lXc.RebuildLocalFromField(Wfield_s, nFields, !backward);
  // lY.SendLocalField(Wfield_s, Wfield_r);
  // lY.RebuildFieldFromLocal(Wfield_r, nFields, !backward);
  //
  // lY.RebuildLocalFromField(Wfield_r, nFields, backward);
  // lY.SendLocalField(Wfield_r, Wfield_s);
  // lXc.RebuildFieldFromLocal(Wfield_s, nFields, backward);
  //
  // std::vector<MyVector<double>>& Xcfield = lXc.GetFields();
  // for(unsigned k = 0; k < nFields; k++) {
  //   for(unsigned i = Xfield[k].begin(); i < Xfield[k].end(); i++) {
  //     if(fabs(Xfield[k][i] - Xcfield[k][i]) > 1.0e-12) std::cerr << "error ";
  //   }
  // }
  //
  // //END TEST

}

void BBoxToIel::GetChildIndexAndLocalXi(const unsigned et,
                                        const std::vector<double>& xi0,
                                        unsigned & childIndex,
                                        std::vector<double>& xi1, const double eps) const {
  // Sets childIndex and xi1 (child reference coords)
  // Child ordering MUST match your verts tables
  if (et == static_cast<unsigned>(LINE)) {
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
#ifndef NDEBUG
    if (!inside_line(xi1[0], eps)) {
      std::cout << xi1[0] << std::endl;
      throw std::runtime_error("project: Line3 mapped xi out of [-1,1]");
    }
#endif
    return;
  }

  if (et == static_cast<unsigned>(QUAD)) {
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
#ifndef NDEBUG
    if (!inside_line(xi1[0], eps) || !inside_line(xi1[1], eps)) {
      std::cout << xi1[0] << " " << xi1[1] << std::endl;
      throw std::runtime_error("project: Quad9 mapped xi out of [-1,1]^2");
    }
#endif
    return;
  }

  if (et == static_cast<unsigned>(HEX)) {
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

#ifndef NDEBUG
    if (!inside_line(xi1[0], eps) || !inside_line(xi1[1], eps) || !inside_line(xi1[2], eps)) {
      std::cout << xi1[0] << " " << xi1[1] << " " << xi1[2] << std::endl;
      throw std::runtime_error("project: Hex27 mapped xi out of [-1,1]^3");
    }
#endif
    return;
  }


  if (et == static_cast<unsigned>(TRI)) {
    // Reference triangle: (r,s) with r>=0, s>=0, r+s<=1
    if (xi0.size() != 2) throw std::runtime_error("project: Tri7 expects xi size 2");
    const double r  = xi0[0];
    const double s  = xi0[1];
    const double rs = r + s;

    if (!inside_tri(r, s, eps)) {
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
    if (!inside_tri(xi1[0], xi1[1], eps)) {
      std::cout << r << " " << s << std::endl;
      std::cout << xi1[0] << " " << xi1[1] << std::endl;
      throw std::runtime_error("project: Tri7 mapped xi not inside child reference triangle");
    }
#endif
    return;
  }


  if (et == static_cast<unsigned>(WEDGE)) {
    // Wedge: (r,s) triangle, z in [-1,1]
    if (xi0.size() != 3) throw std::runtime_error("project: Wedge21 expects xi size 3");
    const double r = xi0[0];
    const double s = xi0[1];
    const double z = xi0[2];

    if (!inside_tri(r, s, eps) || !inside_line(z, eps)) {
      throw std::runtime_error("project: Wedge21 input xi not inside reference wedge");
    }

    const bool up = (z > 0.0);
    const unsigned cz = up ? 1u : 0u;
    const double z1 = up ? (2.0 * z - 1.0) : (2.0 * z + 1.0);

    unsigned ctri = 0;
    double r1 = 0.0, s1 = 0.0;

    if (r <= 0.5 + eps && s <= 0.5 + eps) {
      if (r + s <= 0.5 + eps) {
        ctri = 0;
        r1 = 2.0 * r;
        s1 = 2.0 * s;
      }
      else {
        ctri = 3;
        r1 = 1.0 - 2.0 * r;
        s1 = 1.0 - 2.0 * s;
      }
    }
    else if (r > 0.5) {
      ctri = 1;
      r1 = 2.0 * r - 1.0;
      s1 = 2.0 * s;
    }
    else {
      ctri = 2;
      r1 = 2.0 * r;
      s1 = 2.0 * s - 1.0;
    }

    childIndex = ctri + 4u * cz;


    set3(xi1, r1, s1, z1);

    // Push back to boundaries if tiny overshoot due to roundoff
    pushback_tri(xi1[0], xi1[1], eps);
    pushback_line(xi1[2], eps);

#ifndef NDEBUG
    if (!inside_tri(xi1[0], xi1[1], eps) || !inside_line(xi1[2], eps)) {
      std::cout << "Wedge21 input  : " << r << " " << s << " " << z << "\n";
      std::cout << "Wedge21 mapped : " << xi1[0] << " " << xi1[1] << " " << xi1[2] << "\n";
      throw std::runtime_error("project: Wedge21 mapped xi not inside child reference wedge");
    }
#endif
    return;
  }

  if (et == static_cast<unsigned>(TET)) {
    // Tet: (x,y,z) with x>=0,y>=0,z>=0,x+y+z<=1
    if (xi0.size() != 3) throw std::runtime_error("project: Tet15 expects xi size 3");
    const double x = xi0[0];
    const double y = xi0[1];
    const double z = xi0[2];

    if (!inside_tet(x, y, z, eps)) {
      throw std::runtime_error("project: Tet15 input xi not inside reference tet");
    }

    auto map_child = [&](unsigned c, double & u, double & v, double & w) {
      switch (c) {
        case 0:
          u = 2.0 * x;
          v = 2.0 * y;
          w = 2.0 * z;
          return;
        case 1:
          u = 2.0 * x - 1.0;
          v = 2.0 * y;
          w = 2.0 * z;
          return;
        case 2:
          u = 2.0 * x;
          v = 2.0 * y - 1.0;
          w = 2.0 * z;
          return;
        case 3:
          u = 2.0 * x;
          v = 2.0 * y;
          w = 2.0 * z - 1.0;
          return;

        case 4:
          u = 1.0 - 2.0 * x - 2.0 * z;
          v = 1.0 - 2.0 * y - 2.0 * z;
          w = 2.0 * z;
          return;

        case 5:
          u = 1.0 - 2.0 * x;
          v = 2.0 * y;
          w = 1.0 - 2.0 * y - 2.0 * z;
          return;

        case 6:
          u = 2.0 * y + 2.0 * z - 1.0;
          v = 2.0 * x + 2.0 * z - 1.0;
          w = 1.0 - 2.0 * z;
          return;

        case 7:
          u = 2.0 * x;
          v = 1.0 - 2.0 * y;
          w = 1.0 - 2.0 * x - 2.0 * z;
          return;

        default:
          throw std::runtime_error("project: Tet15 child index out of range");
      }
    };

    // Child priority rule for boundary points:
    // if multiple children contain the mapped point within eps,
    // choose the first child in index order.
    for (unsigned c = 0; c < 8u; ++c) {
      double u, v, w;
      map_child(c, u, v, w);

      // Push back to boundaries if tiny overshoot due to roundoff
      pushback_tet(u, v, w, eps);

      if (inside_tet(u, v, w, eps)) {
        childIndex = c;
        set3(xi1, u, v, w);
#ifndef NDEBUG
        if (!inside_tet(xi1[0], xi1[1], xi1[2], eps)) {
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




