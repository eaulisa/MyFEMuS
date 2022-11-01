#ifndef __femus_Cloud_hpp__
#define __femus_Cloud_hpp__

#include <fstream>
#include <iostream>     // std::cout, std::ios
#include <sstream>      // std::ostringstream

#include "MultiLevelSolution.hpp"
#include "./MyMarker/MyMarker.hpp"
#include "MyEigenFunctions.hpp"


namespace femus {

  template <class Type>
  int mysign(Type x) {
    if(x > 0) return 1;
    if(x < 0) return -1;
    return 0;
  }


  class Cloud {
  public:
    Cloud(Solution* sol) {
      _mrk = MyMarker();
      _sol = sol;
    };
    ~Cloud() {};
    void SetNumberOfMarker(const unsigned &nMax);
    void AddQuadric(const std::vector<double> &A, const unsigned &npt);
    void AddInteriorQuadric(const std::vector<double> &A);

    void AddCircle(const std::vector<double> &xc, const double &R, const unsigned &npt);
    void AddEllipse(const std::vector<double> &xc, const std::vector<double> &a, const unsigned &npt);
    void AddEllipses(const std::vector<std::vector<double>> &xc, const std::vector<std::vector<double>> &a, const std::vector<unsigned> &npt);

    void AddInteriorEllipse(const std::vector<double> &xc, const std::vector<double> &a);
    void AddInteriorEllipses(const std::vector<std::vector<double>> &xc, const std::vector<std::vector<double>> &a);

    void PrintNoOrder(const unsigned &t);
    void PrintWithOrder(const unsigned &t);
    void PrintCSV(const std::string &filename, const unsigned &t);

    void ComputeQuadraticBestFit();

    double GetCost(const std::vector<std::vector<double>>&x,  const std::vector<double> &dotProd, const std::vector<double>&w, const unsigned &iel, const unsigned &nPoints);

    std::pair<std::vector<std::vector<double>>, std::vector<double>> GetCellPointsFromQuadric(const std::vector<std::vector<double>> &xv, const unsigned &iel, unsigned npt, unsigned &nInt, unsigned level = 0);

    void RebuildMarkers(const unsigned &nMin, const unsigned &nMax, const unsigned &npt);

    void RebuildInteriorMarkers(Cloud &intCloud, const std::string &C, const std::string &Cn);

    const std::map<unsigned, std::vector<double>> GetQuadraticBestFitCoefficients() {
      return _A;
    }

    void SetQuadraticBestFitCoefficients(const unsigned &iel, const std::vector<double> &A) {
      _A[iel] = A;
    }


    const std::vector<double> GetQuadraticBestFitCoefficients(const unsigned &iel) {
      if(_A.find(iel) != _A.end()) {
        return _A.at(iel);
      }
      else {
        return {};
      }
    }

    unsigned GetNumberOfMarker(const unsigned &iel) {
      if(_elMrkIdx.find(iel) != _elMrkIdx.end()) {
        return _elMrkIdx[iel][1] - _elMrkIdx[iel][0];
      }
      else {
        return 0;
      }
    }

    double GetValue(const unsigned &iel, const std::vector<double> &xp) {
      return  _A[iel][0] * xp[0] * xp[0] + _A[iel][1] * xp[0] * xp[1] + _A[iel][2] * xp[1] * xp[1] + _A[iel][3] * xp[0] + _A[iel][4] * xp[1] + _A[iel][5];
    }

    double GetValue(const std::vector<double> &A, const std::vector<double> &xp) {
      return  A[0] * xp[0] * xp[0] + A[1] * xp[0] * xp[1] + A[2] * xp[1] * xp[1] + A[3] * xp[0] + A[4] * xp[1] + A[5];
    }


    double GetCurvature(const unsigned &iel, const std::vector<double> &xp) {
      return (8 * _A[iel][0] * _A[iel][2] * _A[iel][2] * xp[1] * xp[1] + 2 * _A[iel][2] * ((_A[iel][3] + 2 * _A[iel][0] * xp[0]) * (_A[iel][3] + 2 * _A[iel][0] * xp[0]) + 4 * _A[iel][0] * (_A[iel][4] + _A[iel][1] * xp[0]) * xp[1] - _A[iel][1] * _A[iel][1] * xp[1] * xp[1]) - 2 * (_A[iel][4] + _A[iel][1] * xp[0]) * (-_A[iel][0] * _A[iel][4] + _A[iel][1] * (_A[iel][3] + _A[iel][0] * xp[0] + _A[iel][1] * xp[1]))) / pow(((_A[iel][4] + _A[iel][1] * xp[0] + 2 * _A[iel][2] * xp[1]) * (_A[iel][4] + _A[iel][1] * xp[0] + 2 * _A[iel][2] * xp[1]) + (_A[iel][3] + 2 * _A[iel][0] * xp[0] + _A[iel][1] * xp[1]) * (_A[iel][3] + 2 * _A[iel][0] * xp[0] + _A[iel][1] * xp[1])), 3. / 2.);
    }

    std::vector<double> GetNormal(const unsigned &iel, const std::vector<double> &xp) {
      std::vector<double> N(xp.size());

      N[0] = 2 * _A[iel][0] * xp[0] + _A[iel][1] * xp[1] + _A[iel][3];
      N[1] = 2 * _A[iel][2] * xp[1] + _A[iel][1] * xp[0] + _A[iel][4];

      double norm2 = 0.;
      for(unsigned i = 0; i < N.size(); i++) norm2 += N[i] * N[i];
      for(unsigned i = 0; i < N.size(); i++) N[i] /= sqrt(norm2);

      return N;
    }

    std::vector<double> GetNormal(const std::vector<double> &A, const std::vector<double> &xp) {
      std::vector<double> N(xp.size());

      N[0] = 2 * A[0] * xp[0] + A[1] * xp[1] + A[3];
      N[1] = 2 * A[2] * xp[1] + A[1] * xp[0] + A[4];

      double norm2 = 0.;
      for(unsigned i = 0; i < N.size(); i++) norm2 += N[i] * N[i];
      for(unsigned i = 0; i < N.size(); i++) N[i] /= sqrt(norm2);

      return N;
    }


    std::vector<double> GetCloudBaricenterInParentElement(const unsigned &iel) {
      unsigned dim = _sol->GetMesh()->GetDimension();
      std::vector <double> yg(dim, 0.);
      if(_elMrkIdx.find(iel) != _elMrkIdx.end()) {
        for(unsigned i = _elMrkIdx[iel][0]; i < _elMrkIdx[iel][1]; i++) {
          for(unsigned k = 0; k < dim; k++)  {
            yg[k] += _yi[_map[i]][k];
          }
        }
        for(unsigned k = 0; k < dim; k++)  yg[k] /= (_elMrkIdx[iel][1] - _elMrkIdx[iel][0]);
      }
      else {
        std::cerr << "In function Cloud::GetCloudBaricenterInParentElement, this element has no marker!!!!!!\n";
        abort();
      }
      return yg;
    }

    double GetAverageCurvature(const unsigned &iel) {
      unsigned i1 = _elMrkIdx[iel][1];
      unsigned i0 = _elMrkIdx[iel][0];
      double avgK = 0.;
      if(i1 - i0 > 0) {
        for(unsigned i = i0; i < i1; i++) {
          avgK += _kappa[_map[i]];
        }
        avgK /= (i1 - i0);
      }
      else {
        std::cerr << "No marker found in function getAverageCurvature \n";
        abort();
      }
      return avgK;
    }

    void RKAdvection(const unsigned & stages, const std::vector<std::string> &U, const double & dt);
    void GetLinearFit(const unsigned & iel, const std::vector<std::vector<double>> &Jac, std::vector < double > &a, double & d);
    void BuildColorFunction(const char C);

  private:

    void CreateMap();
    bool ParallelElementSearch(const std::vector<double> &xp, const unsigned previousElem);

    Solution *_sol;
    unsigned _nMrk;
    std::ofstream _fout;
    std::vector<std::vector<double>> _yp;
    std::vector<std::vector<double>> _N;
    std::vector<double> _kappa;
    std::vector<double> _ds;
    std::vector<std::vector<double>> _yi;

    std::vector<std::vector<double>> _ypNew;
    std::vector<std::vector<double>> _yiNew;
    std::vector<std::vector<double>> _NNew;
    std::vector<double> _kappaNew;
    std::vector<double> _dsNew;


    std::vector<unsigned> _elem;
    std::vector<unsigned> _elemNew;
    std::vector<unsigned> _map;
    MyMarker _mrk;
    std::map<unsigned, std::vector<double>> _A;
    std::map<unsigned, unsigned [2] > _elMrkIdx;
    std::map<unsigned, unsigned [2] >::iterator _itElMrkIdx;

    const std::vector<std::vector<double>> _yig = {{0., 0., 0.}, {1. / 3., 1. / 3., 1. / 3.}, {1. / 3., 1. / 3., 0.}, {0., 0.}, {1. / 3., 1. / 3.}, {0.}};


  };

  void Cloud::SetNumberOfMarker(const unsigned &nMax) {
    _nMrk = nMax;
  }

  void Cloud::AddQuadric(const std::vector<double> &A, const unsigned &npt) {
    Mesh *msh = _sol->GetMesh();
    unsigned dim = _sol->GetMesh()->GetDimension();

    unsigned coordXType = 2;
    std::vector< std::vector < double > > xv;
    std::pair<std::vector<std::vector<double>>, std::vector<double>> sol;


    unsigned solType = 2;
    _mrk.ClearElement();

    unsigned cnt = _elem.size();
    unsigned cnt0 = _elem.size();
    _ypNew.swap(_yp);
    _yiNew.swap(_yi);
    _NNew.swap(_N);
    _kappaNew.swap(_kappa);
    _dsNew.swap(_ds);


    const unsigned &nel = A.size();
    _ypNew.resize(cnt + 2 * nel * npt, std::vector<double>(dim));
    _yiNew.resize(cnt + 2 * nel * npt, std::vector<double>(dim));
    _elem.resize(cnt + 2 * nel * npt);
    _NNew.resize(cnt + 2 * nel * npt, std::vector<double> (dim));
    _kappaNew.resize(cnt + 2 * nel * npt);
    _dsNew.resize(cnt + 2 * nel * npt);

    xv.resize(dim);
    unsigned elCnt = 0;

    unsigned iproc  = msh->processor_id();

    unsigned offset = msh->_elementOffset[iproc];
    unsigned offsetp1 = msh->_elementOffset[iproc + 1];

    for(unsigned iel = offset; iel < offsetp1; iel++) {
//     for(_itElMrkIdx = _elMrkIdx.begin(); _itElMrkIdx != _elMrkIdx.end(); _itElMrkIdx++, elCnt++) {
      bool ielIsDefined = true;
      if(_A.find(iel) == _A.end()) ielIsDefined = false;

      for(unsigned i = 0; i < A.size(); i++) {
        _A[iel].resize(A.size());
        _A[iel][i] = A[i];
      }
      unsigned nDof = msh->GetElementDofNumber(iel, 0);

      for(unsigned k = 0; k < dim; k++) {
        xv[k].resize(nDof);
      }
      for(unsigned k = 0; k < dim; k++) {
        for(unsigned i = 0; i < nDof; i++) {
          unsigned xDof  = msh->GetSolutionDof(i, iel, coordXType);
          xv[k][i] = (*msh->_topology->_Sol[k])(xDof);
        }
      }

      unsigned nInt = 0;

      sol = GetCellPointsFromQuadric(xv, iel, npt, nInt);

      if(nInt == 2) {
        if(ielIsDefined) {
          std::cerr << "Same cell initialized two times in AddQuadric!\n";
          abort();
        }
        elCnt++;
        unsigned newSize = cnt0 + npt * elCnt ; //cnt + sol.first.size() + 2 * (nel - elCnt) * nMax;
        _ypNew.resize(newSize, std::vector<double>(dim));
        _yiNew.resize(newSize, std::vector<double>(dim));
        _elem.resize(newSize);
        _NNew.resize(newSize, std::vector<double>(dim));
        _kappaNew.resize(newSize);
        _dsNew.resize(newSize);

        for(unsigned i = 0; i < sol.first.size(); i++) {
          for(unsigned k = 0; k < dim; k++) {
            _ypNew[cnt][k] = sol.first[i][k];
          }
          _yiNew[cnt]  = _mrk.GetIprocLocalCoordinates(_ypNew[cnt], _sol, solType,  iel);
          _elem[cnt] = iel;
          _NNew[cnt] = GetNormal(iel, _ypNew[cnt]);
          _kappaNew[cnt] = GetCurvature(iel, _ypNew[cnt]);
          _dsNew[cnt] = sol.second[i];
          cnt++;
        }
      }
      else if(!ielIsDefined) _A.erase(iel);

    }

    _ypNew.resize(cnt);
    _yiNew.resize(cnt);
    _elem.resize(cnt);
    _NNew.resize(cnt);
    _kappaNew.resize(cnt);
    _dsNew.resize(cnt);

    _yp.swap(_ypNew);
    _yi.swap(_yiNew);
    _kappa.swap(_kappaNew);
    _ds.swap(_dsNew);
    _N.swap(_NNew);

    CreateMap();
  }


  void Cloud::AddInteriorQuadric(const std::vector<double> &A) {
    Mesh *msh = _sol->GetMesh();

    unsigned dim = msh->GetDimension();

    unsigned iproc  = msh->processor_id();
    unsigned nprocs  = msh->n_processors();

    unsigned offset = msh->_elementOffset[iproc];
    unsigned offsetp1 = msh->_elementOffset[iproc + 1];

    unsigned nMax = offsetp1 - offset;

    unsigned cnt = _elem.size();

    _yp.resize(nMax, std::vector<double>(dim));
    _yi.resize(nMax);
    _N.resize(nMax);
    _kappa.resize(nMax);
    _ds.resize(nMax);
    _elem.resize(nMax);

    for(unsigned iel = offset; iel < offsetp1; iel++) {
      short unsigned ielType = msh->GetElementType(iel);
      unsigned nDof = msh->GetElementDofNumber(iel, 0);  // number of coordinate linear element dofs

      bool ielIsDefined = true;
      if(_A.find(iel) == _A.end()) ielIsDefined = false;

      std::vector<double> x(dim);
      std::vector<double> xm(dim, 0.);

      unsigned cntNode = 0;
      for(unsigned i = 0; i < nDof; i++) {
        unsigned xDof  = msh->GetSolutionDof(i, iel, 2);
        for(unsigned k = 0; k < dim; k++) {
          x[k] = (*msh->_topology->_Sol[k])(xDof);
          xm[k] += x[k];
        }

        if(A[0] * x[0] * x[0] + A[1] * x[0] * x[1] + A[2] * x[1] * x[1] + A[3] * x[0] + A[4] * x[1] + A[5] < 0) cntNode++;

      }

      if(cntNode == nDof) {
        if(ielIsDefined) {
          std::cerr << "Same cell initialized two times in AddInteriorQuadric!\n";
          abort();
        }
        for(unsigned  k = 0; k < dim; k++) {
          _yp[cnt][k] = xm[k] / nDof;
        }
        _yi[cnt]  = _yig[ielType];
        _elem[cnt] = iel;

        _N[cnt] = std::vector<double>(dim, 0.);
        _kappa[cnt] = 0.; //TODO
        _ds[cnt] = 0.;

        cnt++;
      }
    }

    _yp.resize(cnt);
    _yi.resize(cnt);
    _elem.resize(cnt);
    _N.resize(cnt);
    _kappa.resize(cnt);
    _ds.resize(cnt);
    CreateMap();
  }

  void Cloud::AddEllipse(const std::vector<double> &xc, const std::vector<double> &a, const unsigned &npt) {
    std::vector<double> A(6, 0.);
    A[0] = a[1] * a[1];
    A[1] = 0.;
    A[2] = a[0] * a[0];
    A[3] = - 2 * xc[0] * a[1] * a[1];
    A[4] = - 2 * xc[1] * a[0] * a[0];
    A[5] = a[1] * a[1] * xc[0] * xc[0] + a[0] * a[0] * xc[1] * xc[1] - a[0] * a[0] * a[1] * a[1];
    AddQuadric(A, npt);
  }

  void Cloud::AddEllipses(const std::vector<std::vector<double>> &xc, const std::vector<std::vector<double>> &a, const std::vector<unsigned> &npt) {
    if(xc.size() != a.size() || npt.size() != a.size()) {
      std::cerr << "Non-matching vectors in ellipses initialization \n";
      abort();
    }
    else {
      for(unsigned it = 0; it < npt.size(); it++) {
        AddEllipse(xc[it], a[it], npt[it]);
      }
    }
  }

  void Cloud::AddInteriorEllipse(const std::vector<double> &xc, const std::vector<double> &a) {
    std::vector<double> A(6, 0.);
    A[0] = a[1] * a[1];
    A[1] = 0.;
    A[2] = a[0] * a[0];
    A[3] = - 2 * xc[0] * a[1] * a[1];
    A[4] = - 2 * xc[1] * a[0] * a[0];
    A[5] = a[1] * a[1] * xc[0] * xc[0] + a[0] * a[0] * xc[1] * xc[1] - a[0] * a[0] * a[1] * a[1];
    AddInteriorQuadric(A);
  }

  void Cloud::AddInteriorEllipses(const std::vector<std::vector<double>> &xc, const std::vector<std::vector<double>> &a) {
    if(xc.size() != a.size()) {
      std::cerr << "InitMultipleInteriorEllipses arguments not matching!\n";
      abort();
    }
    else {
      for(unsigned j = 0; j < xc.size(); j++) {
        AddInteriorEllipse(xc[j], a[j]);
      }
    }
  }

  void Cloud::AddCircle(const std::vector<double> &xc, const double &R, const unsigned &npt)  {
    AddEllipse(xc, {R, R}, npt);
  }

  void Cloud::CreateMap() {
    std::vector<unsigned*> vec(_elem.size());
    for(unsigned i = 0; i < _elem.size(); i++) {
      vec[i] = &_elem[i];
    }
    std::sort(vec.begin(), vec.end(), [](const unsigned * a, const unsigned * b) {
      return *a < *b;
    });
    _map.resize(_elem.size());
    for(unsigned i = 0; i < _map.size(); i++) {
      _map[i] =  static_cast<unsigned>(vec[i] - &_elem[0]);
    }

    _elMrkIdx.clear();
    unsigned i = 0;
    while(i < _elem.size()) {
      unsigned iel = _elem[_map[i]];
      _elMrkIdx[iel][0] = i;
      while(i < _elem.size() && _elem[_map[i]] == iel) {
        i++;
      }
      _elMrkIdx[iel][1] = i;
    }
  }

  bool Cloud::ParallelElementSearch(const std::vector<double> &xp, const unsigned previousElem = UINT_MAX) {
    bool elemSearch;
    if(previousElem == UINT_MAX) elemSearch = _mrk.ParallelElementSearchWithInverseMapping(xp, _sol, 2);
    else elemSearch = _mrk.ParallelElementSearchWithInverseMapping(xp, _sol, 2, previousElem);

    return elemSearch;
  }

  void Cloud::PrintNoOrder(const unsigned &t) {
    unsigned iproc = _sol->processor_id();
    unsigned nprocs = _sol->n_processors();
    unsigned dim = _sol->GetMesh()->GetDimension();
    for(unsigned kp = 0; kp < nprocs; kp++) {
      if(kp == iproc) {
        if(kp == 0) _fout.open("markerno.dat", std::fstream::out);
        else _fout.open("markerno.dat", std::fstream::app);
        for(unsigned i = 0; i < _yp.size(); i++) {
          for(unsigned k = 0; k < dim; k++) {
            _fout << _yp[i][k] << " ";
          }
          for(unsigned k = 0; k < dim; k++) {
            _fout << _yi[i][k] << " ";
          }
          for(unsigned k = 0; k < dim; k++) {
            _fout << _N[i][k] << " ";
          }
          _fout << _kappa[i] << " " << _ds[i] << " " << iproc << " " << _elem[i] << std::endl;
        }
        _fout.close();
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  void Cloud::PrintWithOrder(const unsigned &t) {
    unsigned iproc = _sol->processor_id();
    unsigned nprocs = _sol->n_processors();
    unsigned dim = _sol->GetMesh()->GetDimension();

    for(unsigned kp = 0; kp < nprocs; kp++) {
      if(kp == iproc) {
        if(kp == 0) _fout.open("markerTest.dat", std::fstream::out);
        else _fout.open("markerTest.dat", std::fstream::app);
        for(unsigned i = 0; i < _yp.size(); i++) {
          for(unsigned k = 0; k < dim; k++) {
            _fout << _yp[_map[i]][k] << " ";
          }
          for(unsigned k = 0; k < dim; k++) {
            _fout << _yi[_map[i]][k] << " ";
          }
          for(unsigned k = 0; k < dim; k++) {
            _fout << _N[_map[i]][k] << " ";
          }
          _fout << _kappa[_map[i]] << " "  << _ds[_map[i]] << " " << iproc << " " << _elem[_map[i]] << std::endl;
        }
        _fout.close();
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  void Cloud::PrintCSV(const std::string &filename, const unsigned &t) {
    unsigned iproc = _sol->processor_id();
    unsigned nprocs = _sol->n_processors();
    unsigned dim = _sol->GetMesh()->GetDimension();

    for(unsigned kp = 0; kp < nprocs; kp++) {
      if(kp == iproc) {
        std::ostringstream foo(std::ostringstream::ate);
        foo.str("./output/");
        foo << filename << t;
        foo << ".csv";
        if(kp == 0) _fout.open(foo.str(), std::fstream::out);
        else _fout.open(foo.str(), std::fstream::app);

        if(kp == 0) {
          _fout << "\"X\",\"Y\",\"Z\",\"xi\",\"eta\",\"zeta\",\"Nx\",\"Ny\",\"Nz\",\"kappa\",\"ds\",\"ipoc\",\"elem\"" << std::endl;
        }
        for(unsigned i = 0; i < _yp.size(); i++) {
          for(unsigned k = 0; k < dim; k++) {
            _fout << _yp[_map[i]][k] << ",";
          }
          _fout << "0.,";
          for(unsigned k = 0; k < dim; k++) {
            _fout << _yi[_map[i]][k] << ",";
          }
          _fout << "0.,";
          for(unsigned k = 0; k < dim; k++) {
            _fout << _N[_map[i]][k] << ",";
          }
          _fout << "0.,";
          _fout << _kappa[_map[i]] << "," << _ds[_map[i]] << "," << iproc << "," << _elem[_map[i]] << std::endl;
        }
        _fout.close();
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  void Cloud::ComputeQuadraticBestFit() {
    _A.clear();

    map<unsigned, bool> pSearch;

    Mesh *msh = _sol->GetMesh();

    unsigned dim = _sol->GetMesh()->GetDimension();
    std::vector<std::vector<double>> coord;
    std::vector<double> weight;
    coord.reserve(_elem.size());
    std::vector<double> norm;
    std::vector<double> xn;

    unsigned iproc = _sol->processor_id();
    unsigned nprocs = _sol->n_processors();

    for(_itElMrkIdx = _elMrkIdx.begin(); _itElMrkIdx != _elMrkIdx.end(); _itElMrkIdx++) {
      unsigned iel = _itElMrkIdx->first;
      unsigned i0 = _itElMrkIdx->second[0];
      unsigned i1 = _itElMrkIdx->second[1];

      unsigned cnt0 = i1 - i0;
      coord.resize(i1 - i0, std::vector<double> (dim));
      weight.resize(i1 - i0);
      norm.assign(dim, 0);
      unsigned cnt = 0;
      double weightSum = 0;

      for(unsigned i = i0; i < i1; i++, cnt++) {
        for(unsigned k = 0; k < dim; k++) {
          coord[cnt][k] = _yp[_map[i]][k];
          norm[k] += _N[_map[i]][k];
        }
        weight[cnt] = _ds[_map[i]];
        weightSum += weight[cnt];
      }

      double det = 0;
      for(unsigned k = 0; k < dim; k++) {
        det += norm[k] * norm[k];
      }
      det = sqrt(det);
      for(unsigned k = 0; k < dim; k++) {
        norm[k] /= det;
      }

      xn.assign(dim, 0);

      for(unsigned i = 0; i < cnt; i++) {
        for(unsigned k = 0; k < dim; k++) {
          xn[k] += weight[i] * coord[i][k];
        }
      }
      for(unsigned k = 0; k < dim; k++) {
        xn[k] /= weightSum;
      }

      unsigned nFaces = msh->GetElementFaceNumber(iel);
      std::vector<unsigned> jelFace(nFaces);
      unsigned cntF = 0;
      for(unsigned iface = 0; iface < nFaces; iface++) {
        int jel = msh->el->GetFaceElementIndex(iel, iface) - 1;
        if(jel >= 0) {
          jelFace[cntF] = jel;
          cntF++;
        }
      }
      jelFace.resize(cntF);

      bool parallelSearch = false;

      for(unsigned i = 1; i < msh->el->GetElementNearElementSize(iel, 1); i++) {
        int jel = msh->el->GetElementNearElement(iel, i);
        unsigned jproc = msh->IsdomBisectionSearch(jel, 3);
        if(jproc != iproc) {
          pSearch[iel] = true;
          parallelSearch = true;
          break;
        }

        bool isFaceElement = false;

        for(unsigned iFace = 0; iFace < jelFace.size(); iFace++) {
          if(jel == jelFace[iFace]) {
            isFaceElement = true;
            break;
          }
        }

        if(_elMrkIdx.find(jel) != _elMrkIdx.end()) { //jel is a cut fem inside iproc
          unsigned j0 = _elMrkIdx[jel][0];
          unsigned j1 = _elMrkIdx[jel][1];

          coord.resize(coord.size() + (j1 - j0), std::vector<double> (dim));
          weight.resize(coord.size() + (j1 - j0), 0.250 * !isFaceElement + isFaceElement);

          for(unsigned j = j0; j < j1; j++) {
            double dotProduct = 0.;
            for(unsigned k = 0; k < dim; k++) {
              dotProduct += norm[k] * _N[_map[j]][k];
            }
            if(acos(dotProduct) <= M_PI) {
              for(unsigned k = 0; k < dim; k++) {
                coord[cnt][k] = _yp[_map[j]][k];
              }
              weight[cnt] = _ds[_map[j]] * (0.250 * !isFaceElement + isFaceElement);
              cnt++;
            }
          }

          coord.resize(cnt);
          weight.resize(cnt);
        }
      }
      if(parallelSearch == false) {

        double sigma2 = 0.;
        double sigma = 0.;
        if(cnt > 1) {
          for(unsigned i = 0; i < cnt; i++) {
            for(unsigned k = 0; k < dim; k++) {
              sigma2 += (coord[i][k] - xn[k]) * (coord[i][k] - xn[k]);
            }
          }
          sigma2 /= cnt;
          sigma2 /= 2;
          sigma = sqrt(sigma2);
          for(unsigned i = 0; i < cnt; i++) {
            double a = 0;
            for(unsigned k = 0; k < dim; k++) {
              a += -0.5 / sigma2 * (coord[i][k] - xn[k]) * (coord[i][k] - xn[k]);
            }
            weight[i] *= exp(a) / (sigma * sqrt(2. * M_PI));
          }
        }
        else {
          std::cerr << "Abbiamo solo un marker!!!!!!!!!!!!!!!!!!!!!!!!!\n" ;
          abort();
        }

        FindQuadraticBestFit(coord, weight, norm, _A[iel]); //conica

        std::vector<double> dotProduct(i1 - i0, 0.);
        double n1Dotn = 0;
        for(unsigned i = i0; i < i1; i++) {
          std::vector <double> n1 = GetNormal(iel,  _yp[_map[i]]);
          for(unsigned k = 0; k < dim; k++) {
            dotProduct[i - i0] += _N[_map[i]][k] * n1[k];
          }
          n1Dotn += dotProduct[i - i0] * weight[i - i0];
        }
        if(n1Dotn < 0) {
          for(unsigned  i = 0; i < _A[iel].size(); i++) {
            _A[iel][i] *= -1.;
          }
          for(unsigned i = 0; i < cnt0; i++) {
            dotProduct[i] *= -1.;
          }
        }

        if(iel == 5588) {
          for(unsigned i = 0; i < cnt0; i++) {
            std::cout << dotProduct[i] << " ";
          }
          std::cout << std::endl;

        }


        double cost1 = GetCost(coord, dotProduct, weight, iel, cnt0);
        std::vector<double> Acon = _A[iel];
        if(cnt0 <= 2) {
          unsigned nEP = 3;
          nEP = (nEP > coord.size() - cnt0) ? (coord.size() - cnt0) : nEP;
          if(coord.size() > 6) {
            std::vector<unsigned> j(nEP);

            std::vector<double*> weightP(nEP);
            for(unsigned i = 0; i < nEP; i++) {
              weightP[i] = &weight[cnt0 + i];
            }
            std::sort(weightP.begin(), weightP.end(), [](const double * a, const double * b) {
              return *a < *b;
            });
            for(unsigned i = 0; i < nEP; i++) {
              j[i] =  cnt0 + static_cast<unsigned>(weightP[i] - &weight[cnt0]);
            }

            for(unsigned k = cnt0 + nEP; k < coord.size(); k++) {
              bool keepWhile = true;
              unsigned kn = nEP - 1;
              while(keepWhile) {
                keepWhile = false;
                if(weight[k] > weight[j[kn]]) {
                  std::swap(k, j[kn]);
                  if(kn > 0) {
                    keepWhile = true;
                    kn--;
                  }
                }
              }
            }
            for(unsigned k = 0; k < nEP; k++) {
              std::swap(weight[j[k]], weight[cnt0 + k]);
              coord[j[k]].swap(coord[cnt0 + k]);
            }


//             unsigned j1 = cnt0 + 1, j2 = cnt0 + 2, j3 = cnt0 + 3;
//             if(weight[j1] < weight[j2]) std::swap(j1, j2);
//             if(weight[j1] < weight[j3]) std::swap(j1, j3);
//             if(weight[j2] < weight[j3]) std::swap(j2, j3);
//             for(unsigned j = cnt0 + 4; j < coord.size(); j++) {
//               if(weight[j] > weight[j3]) {
//                 std::swap(j, j3);
//                 if(weight[j2] < weight[j3]) {
//                   std::swap(j2, j3);
//                   if(weight[j1] < weight[j2]) {
//                     std::swap(j1, j2);
//                   }
//                 }
//               }
//             }
//             std::vector<double> weight1 = weight;
//             std::swap(weight1[j1], weight1[cnt0 + 1]);
//             std::swap(weight1[j2], weight1[cnt0 + 2]);
//             std::swap(weight1[j3], weight1[cnt0 + 3]);
//
//             coord[j1].swap(coord[cnt0 + 1]);
//             coord[j2].swap(coord[cnt0 + 2]);
//             coord[j3].swap(coord[cnt0 + 3]);

            femus::GetQuadricBestFit(coord, weight, norm, _A[iel], cnt0 + 3); //parabola

          }
          else {
            femus::GetQuadricBestFit(coord, weight, norm, _A[iel], coord.size()); //parabola
          }
        }
        else {
          femus::GetQuadricBestFit(coord, weight, norm, _A[iel], cnt0); //parabola
        }

        dotProduct.assign(i1 - i0, 0);
        n1Dotn = 0;
        for(unsigned i = i0; i < i1; i++) {
          std::vector <double> n1 = GetNormal(iel,  _yp[_map[i]]);
          for(unsigned k = 0; k < dim; k++) {
            dotProduct[i - i0] += _N[_map[i]][k] * n1[k];
          }
          n1Dotn += dotProduct[i - i0] * weight[i - i0]; //TODO is this ok if we swap the weights above?
        }
        if(iel == 7736) {
          for(unsigned i = 0; i < cnt0; i++) {
            std::cout << dotProduct[i] << " ";
          }
          std::cout << std::endl << n1Dotn << std::endl;
        }

        if(n1Dotn < 0) {
          for(unsigned  i = 0; i < _A[iel].size(); i++) {
            _A[iel][i] *= -1.;
          }
          for(unsigned i = 0; i < cnt0; i++) {
            dotProduct[i] *= -1.;
          }
        }
        if(iel == 7736) {
          for(unsigned i = 0; i < cnt0; i++) {
            std::cout << dotProduct[i] << " ";
          }
          std::cout << std::endl;

        }

        double cost2 = GetCost(coord, dotProduct, weight, iel, cnt0);
        std::vector<double> Apar = _A[iel];

        double cost3 = 1.e10;
//         if(cnt0 > 2) {
//           std::vector<double> a(dim, 0.);
//           double d = 0.;
//           femus::FindBestFit(coord, weight, norm, a, d);
//           _A[iel] = {0., 0., 0., a[0], a[1], d};
//
//           dotProduct.assign(i1 - i0, 0);
//           n1Dotn = 0;
//           for(unsigned i = i0; i < i1; i++) {
//             std::vector <double> n1 = GetNormal(iel,  _yp[_map[i]]);
//             for(unsigned k = 0; k < dim; k++) {
//               dotProduct[i - i0] += _N[_map[i]][k] * n1[k];
//             }
//             n1Dotn += dotProduct[i - i0] * weight[i - i0];
//           }
//           if(n1Dotn < 0) {
//             for(unsigned  i = 0; i < _A[iel].size(); i++) {
//               _A[iel][i] *= -1.;
//             }
//             for(unsigned i = 0; i < cnt0; i++) {
//               dotProduct[i] *= -1.;
//             }
//           }
//           if(iel == 7243) {
//             for(unsigned i = 0; i < cnt0; i++) {
//               std::cout << dotProduct[i] << " ";
//             }
//             std::cout << std::endl;
//
//           }
//           cost3 = GetCost(coord, dotProduct, weight, iel, cnt0);
//         }

        if((cost2 < 0.0001 && cost1 / cost2 > 0.00001) && cost2 / cost1 > 0.00001) {

          double counter = 0.;
          std::vector <double> xp(dim);
          unsigned nDofs =  msh->GetElementDofNumber(iel, 2);
          for(unsigned i = 0; i < nDofs; i++) {
            unsigned xDof  = msh->GetSolutionDof(i, iel, 2);
            for(unsigned k = 0; k < dim; k++) {
              xp[k] = (*msh->_topology->_Sol[k])(xDof);
            }
            std::vector<double> n1 = GetNormal(Apar, xp);
            std::vector<double> n2 = GetNormal(Acon, xp);
            double dotProd = 0;
            for(unsigned k = 0; k < dim; k++) {
              dotProd += n1[k] * n2[k];
            }

            if(iel == 7736) {
              std::cout << i << " " <<  GetValue(Apar, xp) << " " << GetValue(Acon, xp) << " " << counter << std::endl;
            }

            if(dotProd < 0 && GetValue(Apar, xp) * GetValue(Acon, xp) < 0) {
              double distMin = 1.0E10;
              unsigned jMin = 0;
              for(unsigned j = 1; j < cnt0; j++) {
                double distj = 0.;
                for(unsigned k = 0; k < dim; k++) {
                  distj += (coord[j][k] - xp[k]) * (coord[j][k] - xp[k]);
                }
                if(distj > distMin) {
                  jMin = j;
                  distMin = distj;
                }
              }

              double dotProd1 = 0;
              for(unsigned k = 0; k < dim; k++) {
                dotProd1 += n1[k] * _N[_map[i0 + jMin]][k];
              }

              double dotProd2 = 0;
              for(unsigned k = 0; k < dim; k++) {
                dotProd2 += n2[k] * _N[_map[i0 + jMin]][k];
              }

              counter += dotProd1 - dotProd2;

              if(iel == 7736) {
                std::cout << i << " " <<  dotProd1 << " " << dotProd2 << " " << counter << std::endl;
              }

//               if( dotProd1 > dotProd2 ) {
//                 counter += (dotProd1 - dotProd2);
//               }
//               else if {
//                 counter -= (dotProd2 - dotProd1);
//               }
            }
          }
          if(iel == 7736) {


            std::cout << " counter = " << counter << std::endl;
            std::cout << "a = " << Acon[0] << "; b=" << Acon[1] << "; c=" << Acon[2] << "; d=" << Acon[3] << "; e=" << Acon[4] << "; f= " << Acon[5] << ";\n";
            std::cout << "a = " << Apar[0] << "; b=" << Apar[1] << "; c=" << Apar[2] << "; d=" << Apar[3] << "; e=" << Apar[4] << "; f= " << Apar[5] << ";\n";

          }
          if(counter > 0) {
            Acon = Apar;
          }
          else if(counter < 0) {
            Apar = Acon;
          }
        }


//         unsigned xDof0  = msh->GetSolutionDof(0, iel, 2);
//         unsigned xDof2  = msh->GetSolutionDof(2, iel, 2);
//         double treshold = 0;
//         for(unsigned k = 0; k < dim; k++) {
//           double h = (*msh->_topology->_Sol[k])(xDof2) - (*msh->_topology->_Sol[k])(xDof0);
//           treshold += 0.0001 * h * h;
//         }
//         treshold *= cnt0;
//         if(cost2 < treshold) {
//           _A[iel] = Apar;
//         }
//         else {
        _A[iel] = (cost1 < cost2) ? ((cost1 < cost3) ? Acon : _A[iel]) : ((cost2 < cost3) ? Apar : _A[iel]);
//         }

//         if(cost3 < cost1 && cost3 < cost2) {
//           std::cout << "CCCCC " << iel << std::endl;
//         }



        if(iel == 7736) {
          double xMin[2] = { 1.0e10, 1.0e10};
          double xMax[2] = { -1.0e10, -1.0e10};

          unsigned nDofs =  msh->GetElementDofNumber(iel, 2);
          for(unsigned i = 0; i < nDofs; i++) {
            unsigned xDof  = msh->GetSolutionDof(i, iel, 2);
            for(unsigned k = 0; k < dim; k++) {
              double x = (*msh->_topology->_Sol[k])(xDof);
              if(x < xMin[k]) xMin[k] = x;
              if(x > xMax[k]) xMax[k] = x;
            }
          }

          std::cout << "BBBBBBB " ;//<< treshold <<  " ";
          std::cout << cost1 << " " << cost2 << " " << cost3 << std::endl;
          std::cout << "a = " << _A[iel][0] << "; b=" << _A[iel][1] << "; c=" << _A[iel][2] << "; d=" << _A[iel][3] << "; e=" << _A[iel][4] << "; f= " << _A[iel][5] << ";\n";
          std::cout << " {x," << xMin[0] << ", " << xMax[0] << "},{" << "{y," << xMin[1] << ", " << xMax[1] << "}" << std::endl;
        }
        /*
                n1Dotn = 0;

                for(unsigned i = i0; i < i1; i++) {
                  std::vector <double> n1 = GetNormal(iel,  _yp[_map[i]]);
                  for(unsigned k = 0; k < dim; k++) {
                    n1Dotn += _N[_map[i]][k] * n1[k];
                  }
                }
                if(n1Dotn < 0) {
                  for(unsigned  i = 0; i < _A[iel].size(); i++) {
                    _A[iel][i] *= -1.;
                  }
                }*/
      }
    }

    map<unsigned, bool>::iterator it;


    if(nprocs > 1) {

      for(unsigned kp = 0; kp < nprocs; kp++) {

        unsigned nel;
        if(iproc == kp) {
          nel = pSearch.size();
        }
        MPI_Bcast(&nel, 1, MPI_UNSIGNED, kp, MPI_COMM_WORLD);

        if(nel > 0) {
          if(iproc == kp) {
            it =  pSearch.begin();
          }
          for(unsigned cntEl = 0; cntEl < nel; cntEl++) {
            unsigned kel;
            unsigned nNgbElms;
            unsigned i0;
            unsigned i1;

            std::vector<unsigned> jelFace;
            unsigned cnt0 = 0;

            if(iproc == kp) {
              kel = it->first;
              i0 = _elMrkIdx[kel][0];
              i1 = _elMrkIdx[kel][1];
              coord.resize(i1 - i0, std::vector<double> (dim));
              weight.resize(i1 - i0);
              cnt0 = i1 - i0;
              norm.assign(dim, 0);
              unsigned cnt = 0;
              double weightSum = 0;
              for(unsigned i = i0; i < i1; i++, cnt++) {
                for(unsigned k = 0; k < dim; k++) {
                  coord[cnt][k] = _yp[_map[i]][k];
                  norm[k] += _N[_map[i]][k];
                }
                weight[cnt] = _ds[_map[i]];
                weightSum += weight[cnt];
              }

              xn.assign(dim, 0);

              for(unsigned i = 0; i < cnt; i++) {
                for(unsigned k = 0; k < dim; k++) {
                  xn[k] += weight[i] * coord[i][k];
                }
              }
              for(unsigned k = 0; k < dim; k++) {
                xn[k] /= weightSum;
              }

              unsigned nFaces = msh->GetElementFaceNumber(kel);
              jelFace.resize(nFaces);
              unsigned cntF = 0;
              for(unsigned iface = 0; iface < nFaces; iface++) {
                int jel = msh->el->GetFaceElementIndex(kel, iface) - 1;
                if(jel >= 0) {
                  jelFace[cntF] = jel;
                  cntF++;
                }
              }
              jelFace.resize(cntF);
              nNgbElms = msh->el->GetElementNearElementSize(kel, 1);
            }
            MPI_Bcast(&nNgbElms, 1, MPI_UNSIGNED, kp, PETSC_COMM_WORLD);
            jelFace.resize(nNgbElms);
            MPI_Bcast(jelFace.data(), jelFace.size(), MPI_UNSIGNED, kp, PETSC_COMM_WORLD);
//             MPI_Bcast(&i1, 1, MPI_UNSIGNED, kp, PETSC_COMM_WORLD);
//             MPI_Bcast(&i0, 1, MPI_UNSIGNED, kp, PETSC_COMM_WORLD);

            for(unsigned i = 1; i < nNgbElms; i++) {

              int jel;
              if(iproc == kp) {
                jel = msh->el->GetElementNearElement(kel, i);
              }
              MPI_Bcast(&jel, 1, MPI_INT, kp, PETSC_COMM_WORLD);

              unsigned jp = msh->IsdomBisectionSearch(jel, 3);  // return  jproc for piece-wise constant discontinuous type (3)

              bool isFaceElement = false;
              for(unsigned iFace = 0; iFace < jelFace.size(); iFace++) {
                if(jel == jelFace[iFace]) {
                  isFaceElement = true;
                  break;
                }
              }

              std::vector<std::vector<double>> coordJel;
              std::vector<double> weightJel;
              unsigned cntJel = 0;
              if(iproc == jp) {
                if(_elMrkIdx.find(jel) != _elMrkIdx.end()) {   // if jel is cut cell
                  unsigned j0 = _elMrkIdx[jel][0];
                  unsigned j1 = _elMrkIdx[jel][1];
                  coordJel.resize(dim, std::vector<double> (j1 - j0));
                  weightJel.resize(j1 - j0);
                  for(unsigned j = j0; j < j1; j++, cntJel++) {
                    for(unsigned k = 0; k < dim; k++) {
                      coordJel[k][cntJel] = _yp[_map[j]][k];
                    }
                    weightJel[cntJel] = _ds[_map[j]] * (0.250 * !isFaceElement + isFaceElement);
                  }
                }

                if(jp != kp) {



                  MPI_Send(&cntJel, 1, MPI_UNSIGNED, kp, 0, MPI_COMM_WORLD);
                  if(cntJel != 0) {
                    for(unsigned k = 0; k < dim; k++) {
                      MPI_Send(coordJel[k].data(), coordJel[k].size(), MPI_DOUBLE, kp, 1 + k, MPI_COMM_WORLD);
                    }
                    MPI_Send(weightJel.data(), weightJel.size(), MPI_DOUBLE, kp, dim + 1, MPI_COMM_WORLD);
                  }
                }
              }

              if(iproc == kp) {
                if(kp != jp) {

                  MPI_Recv(&cntJel, 1, MPI_UNSIGNED, jp, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                  if(cntJel != 0) {
                    coordJel.resize(dim, std::vector<double> (cntJel));
                    weightJel.resize(cntJel);
                    for(unsigned k = 0; k < dim; k++) {
                      MPI_Recv(coordJel[k].data(), coordJel[k].size(), MPI_DOUBLE, jp, 1 + k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                    MPI_Recv(weightJel.data(), weightJel.size(), MPI_DOUBLE, jp, dim + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                  }
                }
                if(cntJel != 0) {

                  unsigned size0 = coord.size();
                  coord.resize(size0 + cntJel, std::vector<double> (dim));
                  weight.resize(size0 + cntJel);
                  for(unsigned j = 0; j < cntJel; j++) {
                    for(unsigned k = 0; k < dim; k++) {
                      coord[size0 + j][k] = coordJel[k][j];
                    }
                    weight[size0 + j] = weightJel[j];
                  }
                }
              }
            }//face loop

            if(iproc == kp) {

              unsigned cnt = coord.size();
              double sigma2 = 0.;
              double sigma = 0.;
              if(cnt > 1) {
                for(unsigned i = 0; i < cnt; i++) {
                  for(unsigned k = 0; k < dim; k++) {
                    sigma2 += (coord[i][k] - xn[k]) * (coord[i][k] - xn[k]);
                  }
                }
                sigma2 /= cnt;
                sigma2 /= 2;
                sigma = sqrt(sigma2);
                for(unsigned i = 0; i < cnt; i++) {
                  double a = 0;
                  for(unsigned k = 0; k < dim; k++) {
                    a += -0.5 / sigma2 * (coord[i][k] - xn[k]) * (coord[i][k] - xn[k]);
                  }
                  weight[i] *= exp(a) / (sigma * sqrt(2. * M_PI));
                }
              }
              else {
                std::cerr << "Abbiamo solo un marker!!!!!!!!!!!!!!!!!!!!!!!!!\n" ;
                abort();
              }

              FindQuadraticBestFit(coord, weight, norm, _A[kel]); //conica

              std::vector<double> dotProduct(i1 - i0, 0.);
              double n1Dotn = 0;
              for(unsigned i = i0; i < i1; i++) {
                std::vector <double> n1 = GetNormal(kel,  _yp[_map[i]]);
                for(unsigned k = 0; k < dim; k++) {
                  dotProduct[i - i0] += _N[_map[i]][k] * n1[k];
                }
                n1Dotn += dotProduct[i - i0] * weight[i - i0];
              }
              if(n1Dotn < 0) {
                for(unsigned  i = 0; i < _A[kel].size(); i++) {
                  _A[kel][i] *= -1.;
                }
                for(unsigned i = 0; i < cnt0; i++) {
                  dotProduct[i] *= -1.;
                }
              }

              double cost1 = GetCost(coord, dotProduct, weight, kel, cnt0);
              std::vector<double> Acon = _A[kel];
              if(cnt0 <= 2) {
                unsigned nEP = 3;
                nEP = (nEP > coord.size() - cnt0) ? (coord.size() - cnt0) : nEP;
                if(coord.size() > 6) {
                  std::vector<unsigned> j(nEP);

                  std::vector<double*> weightP(nEP);
                  for(unsigned i = 0; i < nEP; i++) {
                    weightP[i] = &weight[cnt0 + i];
                  }
                  std::sort(weightP.begin(), weightP.end(), [](const double * a, const double * b) {
                    return *a < *b;
                  });
                  for(unsigned i = 0; i < nEP; i++) {
                    j[i] =  cnt0 + static_cast<unsigned>(weightP[i] - &weight[cnt0]);
                  }

                  for(unsigned k = cnt0 + nEP; k < coord.size(); k++) {
                    bool keepWhile = true;
                    unsigned kn = nEP - 1;
                    while(keepWhile) {
                      keepWhile = false;
                      if(weight[k] > weight[j[kn]]) {
                        std::swap(k, j[kn]);
                        if(kn > 0) {
                          keepWhile = true;
                          kn--;
                        }
                      }
                    }
                  }
                  for(unsigned k = 0; k < nEP; k++) {
                    std::swap(weight[j[k]], weight[cnt0 + k]);
                    coord[j[k]].swap(coord[cnt0 + k]);
                  }

                  femus::GetQuadricBestFit(coord, weight, norm, _A[kel], cnt0 + 3); //parabola

                }
                else {
                  femus::GetQuadricBestFit(coord, weight, norm, _A[kel], coord.size()); //parabola
                }
              }
              else {
                femus::GetQuadricBestFit(coord, weight, norm, _A[kel], cnt0); //parabola
              }

              dotProduct.assign(i1 - i0, 0);
              n1Dotn = 0;
              for(unsigned i = i0; i < i1; i++) {
                std::vector <double> n1 = GetNormal(kel,  _yp[_map[i]]);
                for(unsigned k = 0; k < dim; k++) {
                  dotProduct[i - i0] += _N[_map[i]][k] * n1[k];
                }
                n1Dotn += dotProduct[i - i0] * weight[i - i0]; //TODO is this ok if we swap the weights above?
              }

              if(n1Dotn < 0) {
                for(unsigned  i = 0; i < _A[kel].size(); i++) {
                  _A[kel][i] *= -1.;
                }
                for(unsigned i = 0; i < cnt0; i++) {
                  dotProduct[i] *= -1.;
                }
              }

              double cost2 = GetCost(coord, dotProduct, weight, kel, cnt0);
              std::vector<double> Apar = _A[kel];

              double cost3 = 1.e10;

              if((cost2 < 0.0001 && cost1 / cost2 > 0.00001) && cost2 / cost1 > 0.00001) {

                double counter = 0.;
                std::vector <double> xp(dim);
                unsigned nDofs =  msh->GetElementDofNumber(kel, 2);
                for(unsigned i = 0; i < nDofs; i++) {
                  unsigned xDof  = msh->GetSolutionDof(i, kel, 2);
                  for(unsigned k = 0; k < dim; k++) {
                    xp[k] = (*msh->_topology->_Sol[k])(xDof);
                  }
                  std::vector<double> n1 = GetNormal(Apar, xp);
                  std::vector<double> n2 = GetNormal(Acon, xp);
                  double dotProd = 0;
                  for(unsigned k = 0; k < dim; k++) {
                    dotProd += n1[k] * n2[k];
                  }

                  if(dotProd < 0 && GetValue(Apar, xp) * GetValue(Acon, xp) < 0) {
                    double distMin = 1.0E10;
                    unsigned jMin = 0;
                    for(unsigned j = 1; j < cnt0; j++) {
                      double distj = 0.;
                      for(unsigned k = 0; k < dim; k++) {
                        distj += (coord[j][k] - xp[k]) * (coord[j][k] - xp[k]);
                      }
                      if(distj > distMin) {
                        jMin = j;
                        distMin = distj;
                      }
                    }

                    double dotProd1 = 0;
                    for(unsigned k = 0; k < dim; k++) {
                      dotProd1 += n1[k] * _N[_map[i0 + jMin]][k];
                    }

                    double dotProd2 = 0;
                    for(unsigned k = 0; k < dim; k++) {
                      dotProd2 += n2[k] * _N[_map[i0 + jMin]][k];
                    }

                    counter += dotProd1 - dotProd2;

                  }
                }
                if(counter > 0) {
                  Acon = Apar;
                }
                else if(counter < 0) {
                  Apar = Acon;
                }
              }

              _A[kel] = (cost1 < cost2) ? ((cost1 < cost3) ? Acon : _A[kel]) : ((cost2 < cost3) ? Apar : _A[kel]);

              it++;
            }

          }//element loop
        }
      }


    }


  }


  const double PJ[4][4][4] = {
    {{1.}, {0.5, 0.5}, {0.25, 0.25, 0.25, 0.25}, {0.5, 0., 0., 0.5}},
    {{0.5, 0.5}, {0., 1.}, {0., 0.5, 0.5}, {0.25, 0.25, 0.25, 0.25}},
    {{0.25, 0.25, 0.25, 0.25}, {0., 0.5, 0.5}, {0., 0., 1.}, {0., 0., 0.5, 0.5}},
    {{0.5, 0., 0., 0.5}, {0.25, 0.25, 0.25, 0.25}, {0., 0., 0.5, 0.5}, {0., 0., 0., 1.}}
  };
#include <boost/math/special_functions/math_fwd.hpp>
  std::pair<std::vector<std::vector<double>>, std::vector<double>> Cloud::GetCellPointsFromQuadric(const std::vector<std::vector<double>> &xv, const unsigned & iel, unsigned npt, unsigned & nInt, unsigned level) {

    unsigned cnt = 0;
    const unsigned dim = xv.size();
    std::vector < std::vector <double> > xe(((8 < npt) ? npt : 8), std::vector<double>(dim));
    std::vector <double> ds(npt);

    if(_A.find(iel) != _A.end()) {

      const unsigned nve = xv[0].size();
      const std::vector<double> &Cf = _A[iel];
      std::vector<double> v(dim, 0.);

      for(unsigned i = 0; i < nve; i++) {
        unsigned ip1 = (i + 1) % nve;
        for(unsigned k = 0; k < dim; k++) v[k] = xv[k][ip1] - xv[k][i];

        const double &x0 = xv[0][i];
        const double &y0 = xv[1][i];

        oct a = Cf[0] * v[0] * v[0] + Cf[1] * v[0] * v[1] + Cf[2] * v[1] * v[1];
        oct b = 2 * Cf[0] * v[0] * x0 + Cf[1] * v[1] * x0 + Cf[1] * v[0] * y0 + 2 * Cf[2] * v[1] * y0 + Cf[3] * v[0] + Cf[4] * v[1];
        oct c = Cf[0] * x0 * x0 + Cf[1] * x0 * y0 + Cf[2] * y0 * y0 + Cf[3] * x0 + Cf[4] * y0 + Cf[5];

        oct norm = sqrt(a * a + b * b + c * c);
        a /= norm;
        b /= norm;
        c /= norm;

        if(fabs(a) > 1.e-8) {
          oct delta = b * b - 4 * a * c;
          if(delta > 0) {
            for(unsigned j = 0; j < 2; j++) {
              double t = static_cast<double>((- b + pow(-1, j) * sqrt(delta)) / (2 * a));
              if(t >= 0 && t <= 1) {
                for(unsigned  k = 0; k < dim; k++) {
                  xe[cnt][k] = xv[k][i]  + t * v[k];
                }
                cnt++;
              }
            }
          }
        }
        else if(b != 0) {
          double t = static_cast<double>(-c / b);
          if(t >= 0 && t <= 1) {
            for(unsigned  k = 0; k < dim; k++) {
              xe[cnt][k] = xv[k][i]  + t * v[k];
            }
            cnt++;
          }
        }
      }

      nInt = cnt;

      if(cnt == 2) {

        cnt = 1;
        xe.resize(npt, std::vector<double>(dim));
        xe[npt - 1][0] = xe[1][0];
        xe[npt - 1][1] = xe[1][1];


        double kappa, dt0;

        for(unsigned i = 0; i < npt - 2; i++) {

          double x0 = xe[cnt - 1][0];
          double y0 = xe[cnt - 1][1];
          double x1 = xe[npt - 1][0];
          double y1 = xe[npt - 1][1];

          std::vector<double> xm = {0.5 * (x0 + x1), 0.5 * (y0 + y1)};
          std::vector<double> N = GetNormal(iel, xm);
          kappa = GetCurvature(iel, xm);
          if(fabs(kappa) < 1.e-5) kappa = 1.e-5; //TODO
          std::vector<double> xc = {xm[0] - N[0] / kappa, xm[1] - N[1] / kappa};

          double theta0 = atan2(y0 - xc[1], x0 - xc[0]);
          double theta1 = atan2(y1 - xc[1], x1 - xc[0]);

          dt0 = (theta1 > theta0) ? theta1 - theta0 : 2 * M_PI + theta1 - theta0;
          double dt1 = (theta0 > theta1) ? theta0 - theta1 : 2 * M_PI + theta0 - theta1;

          //std::cout << i << " " << theta0 << " " << theta1 << " " << x0 << " " << x1 << " " << y0 << " " << y1 << std::endl;

          if(dt0 > dt1) {
            dt0 = dt1;
            theta0 = theta1;
            xe[npt - 1] = {x0, y0};
            xe[cnt - 1] = {x1, y1};
            x0 = xe[cnt - 1][0];
            y0 = xe[cnt - 1][1];
            x1 = xe[npt - 1][0];
            y1 = xe[npt - 1][1];
          }
          double R = sqrt((x0 - xc[0]) * (x0 - xc[0]) + (y0 - xc[1]) * (y0 - xc[1]));


//           v[0] = R * cos(theta0 + (i + 1) * dt0 / (npt - 1));
//           v[1] = R * sin(theta0 + (i + 1) * dt0 / (npt - 1));

          v[0] = R * cos(theta0 + dt0 / (npt - i - 1));
          v[1] = R * sin(theta0 + dt0 / (npt - i - 1));

          oct a = Cf[0] * v[0] * v[0] + Cf[1] * v[0] * v[1] + Cf[2] * v[1] * v[1];
          oct b = 2 * Cf[0] * v[0] * xc[0] + Cf[1] * v[1] * xc[0] + Cf[1] * v[0] * xc[1] + 2 * Cf[2] * v[1] * xc[1] + Cf[3] * v[0] + Cf[4] * v[1];
          oct c = Cf[0] * xc[0] * xc[0] + Cf[1] * xc[0] * xc[1] + Cf[2] * xc[1] * xc[1] + Cf[3] * xc[0] + Cf[4] * xc[1] + Cf[5];

          oct norm = sqrt(a * a + b * b + c * c);
          a /= norm;
          b /= norm;
          c /= norm;

          //unsigned cnti = cnt;

          if(fabs(a) > 1.e-8) {
            oct delta = b * b - 4 * a * c;
            if(delta >= 0.) {
              double t[2];
              for(unsigned j = 0; j < 2; j++) {
                t[j] = static_cast<double>((- b + pow(-1, j) * sqrt(delta)) / (2 * a));
              }
              double ti = (fabs(t[0] - 1.) < fabs(t[1] - 1.)) ? t[0] : t[1];
              for(unsigned  k = 0; k < dim; k++) {
                xe[cnt][k] = xc[k]  + ti * v[k];
              }
              cnt++;
            }
          }
          else if(b != 0) {
            double t = static_cast<double>(-c / b);
            for(unsigned  k = 0; k < dim; k++) {
              xe[cnt][k] = xc[k]  + t * v[k];
            }
            cnt++;
          }
        }
        double dsi = 2. * dt0 / (fabs(kappa) * cnt);
        ds.assign(cnt + 1, dsi);
        ds[0] = ds[cnt] = 0.5 * dsi;

        if(cnt < npt - 1) {
          xe[cnt] = xe[npt - 1];
          xe.resize(cnt + 1);
        }









        cnt++;
        npt = cnt;
      }
//       else {
//
// //         std::cerr << iel << " " << level << " " << cnt << std::endl;
//
//         xe.resize(0);
//         if(cnt > 2) {
//           xe.reserve(4 * npt);
//           std::vector<std::vector<double> > xvj(dim, std::vector<double>(nve));
//           for(unsigned j = 0; j < 4; j++) {
//             xvj.assign(dim, std::vector<double>(nve, 0.));
//             for(unsigned k = 0; k < dim; k++) {
//               for(unsigned I = 0; I < nve; I++) {
//                 for(unsigned J = 0 ; J < nve; J++) {
//                   xvj[k][I] += PJ[j][I][J] * xv[k][J];
//                 }
//               }
//             }
//             unsigned nInt = 0;
//             std::vector <std::vector<double>> xej = GetCellPointsFromQuadric(xvj, iel, npt, nInt, level + 1);
//             xe.insert(xe.end(), xej.begin(), xej.end());
//           }
//         }
//         return xe;
//       }
    }

    if(cnt < npt) {
      xe.resize(cnt);
      npt = cnt;
    }

    return std::pair<std::vector<std::vector<double>>, std::vector<double>>(xe, ds);
  }

  void Cloud::RebuildMarkers(const unsigned & nMin, const unsigned & nMax, const unsigned & npt) {
    Mesh *msh = _sol->GetMesh();
    unsigned dim = _sol->GetMesh()->GetDimension();
    unsigned coordXType = 2;
    std::vector< std::vector < double > > xv;
    std::pair<std::vector<std::vector<double>>, std::vector<double>> sol;

    unsigned cnt = 0;

    const unsigned &nel = _A.size();
    _ypNew.resize(2 * nel * nMax, std::vector<double>(dim));
    _elem.resize(2 * nel * nMax);
    _NNew.resize(2 * nel * nMax, std::vector<double> (dim));
    _kappaNew.resize(2 * _A.size() * nMax);
    _dsNew.resize(2 * _A.size() * nMax);

    xv.resize(dim);
    unsigned elCnt = 0;
    for(_itElMrkIdx = _elMrkIdx.begin(); _itElMrkIdx != _elMrkIdx.end(); _itElMrkIdx++, elCnt++) {
      unsigned iel = _itElMrkIdx->first;
      unsigned i0 = _itElMrkIdx->second[0];
      unsigned i1 = _itElMrkIdx->second[1];

      unsigned nDof = msh->GetElementDofNumber(iel, 0);

      double t;
      if(fabs(_A[iel][1]) < 1.e-4) t = 0.;
      else if(fabs(_A[iel][0] - _A[iel][2]) < 1.e-4) t = M_PI / 4.;
      else t = 0.5 * atan(_A[iel][1] / (_A[iel][0] - _A[iel][2]));

      double ap = _A[iel][0] * cos(t) * cos(t) + _A[iel][1] * cos(t) * sin(t) + _A[iel][2] * sin(t) * sin(t);
      double cp = _A[iel][2] * cos(t) * cos(t) - _A[iel][1] * cos(t) * sin(t) + _A[iel][0] * sin(t) * sin(t);

      bool keepMrk = (fabs(ap / cp) > 100 || fabs(cp / ap) > 100) ? true : false;
      keepMrk = false;
      if(i1 - i0 < 2) keepMrk = true;

      keepMrk = false;


      for(unsigned k = 0; k < dim; k++) {
        xv[k].resize(nDof);
      }
      for(unsigned k = 0; k < dim; k++) {
        for(unsigned i = 0; i < nDof; i++) {
          unsigned xDof  = msh->GetSolutionDof(i, iel, coordXType);
          xv[k][i] = (*msh->_topology->_Sol[k])(xDof);
        }
      }

      unsigned nInt = 0;
      if(((i1 - i0) < nMin || (i1 - i0) > nMax) && !keepMrk) {
        sol = GetCellPointsFromQuadric(xv, iel, npt, nInt);
        if(nInt == 2) {
          if(cnt + sol.first.size() > _ypNew.size()) {
            unsigned newSize = cnt + sol.first.size() + 2 * (nel - elCnt) * nMax;
            _ypNew.resize(newSize, std::vector<double>(dim));
            _elem.resize(newSize);
            _NNew.resize(newSize, std::vector<double>(dim));
            _kappaNew.resize(newSize);
            _dsNew.resize(newSize);
          }

          for(unsigned i = 0; i < sol.first.size(); i++) {
            for(unsigned k = 0; k < dim; k++) {
              _ypNew[cnt][k] = sol.first[i][k];
            }
            _elem[cnt] = iel;
            _NNew[cnt] = GetNormal(iel, _ypNew[cnt]);
            _kappaNew[cnt] = GetCurvature(iel, _ypNew[cnt]);
            _dsNew[cnt] = sol.second[i];
            cnt++;
          }
        }
        else {
          keepMrk = true;
        }
      }
      if(((i1 - i0) >= nMin && (i1 - i0) <= nMax) || keepMrk) {

        if(cnt + (i1 - i0) > _ypNew.size()) {
          unsigned newSize = cnt + (i1 - i0) + 2 * (nel - elCnt) * nMax;
          _ypNew.resize(newSize, std::vector<double>(dim));
          _elem.resize(newSize);
          _NNew.resize(newSize, std::vector<double>(dim));
          _kappaNew.resize(newSize);
          _dsNew.resize(newSize);
        }

        for(unsigned i = i0; i < i1; i++) {
          _ypNew[cnt] = _yp[_map[i]];
          if(keepMrk) {
            _NNew[cnt] = _N[_map[i]];
            _kappaNew[cnt] = _kappa[_map[i]];
          }
          else {
            _NNew[cnt] = GetNormal(iel, _ypNew[cnt]);
            _kappaNew[cnt] = GetCurvature(iel, _ypNew[cnt]);
          }
          _dsNew[cnt] = _ds[_map[i]];
          _elem[cnt] = iel;
          cnt++;
        }
      }
    }

    _ypNew.resize(cnt);
    _elem.resize(cnt);
    _NNew.resize(cnt);
    _kappaNew.resize(cnt);
    _dsNew.resize(cnt);

    _yp.swap(_ypNew);
    _kappa.swap(_kappaNew);
    _ds.swap(_dsNew);
    _N.swap(_NNew);

    CreateMap();
    _yi.resize(cnt);

    unsigned solType = 2;
    _mrk.ClearElement();
    for(_itElMrkIdx = _elMrkIdx.begin(); _itElMrkIdx != _elMrkIdx.end(); _itElMrkIdx++) {
      for(unsigned i = _itElMrkIdx->second[0]; i < _itElMrkIdx->second[1]; i++) {
        _yi[_map[i]]  = _mrk.GetIprocLocalCoordinates(_yp[_map[i]], _sol, solType,  _itElMrkIdx->first);
      }
    }

  }

  void Cloud::RebuildInteriorMarkers(Cloud & intCloud, const std::string & C, const std::string & Cn) {

    std::vector<std::vector<double>> xn;

    std::vector<double> phi;
    std::vector<double> phi_x;

    Mesh *msh = _sol->GetMesh();
    unsigned SolCIndex = _sol->GetIndex(C.c_str());
    unsigned SolCnIndex = _sol->GetIndex(Cn.c_str());

    unsigned solType = _sol->GetSolutionType(SolCnIndex);

    unsigned iproc  = msh->processor_id();
    unsigned nprocs  = msh->n_processors();

    unsigned offset = msh->_elementOffset[iproc];
    unsigned offsetp1 = msh->_elementOffset[iproc + 1];

    unsigned dim = msh->GetDimension();

    _sol->_Sol[SolCIndex]->zero();
    _sol->_Sol[SolCnIndex]->zero();

    for(_itElMrkIdx = intCloud._elMrkIdx.begin(); _itElMrkIdx != intCloud._elMrkIdx.end(); _itElMrkIdx++) {
      unsigned iel = _itElMrkIdx->first;
      unsigned ielType = msh->GetElementType(iel);

      unsigned solTypeL = 0;
      unsigned nDofsL = msh->GetElementDofNumber(iel, solTypeL);  // number of coordinate linear element dofs

      unsigned nDofs = msh->GetElementDofNumber(iel, solType);

      xn.resize(dim);
      for(unsigned k = 0; k < dim; k++) {
        xn[k].resize(nDofs);
      }

      for(unsigned i = 0; i < nDofs; i++) {
        unsigned xDof  = msh->GetSolutionDof(i, iel, 2);
        for(unsigned k = 0; k < dim; k++) {
          xn[k][i] = (*msh->_topology->_Sol[k])(xDof);
        }
      }


      const elem_type *femL = fem.GetFiniteElement(ielType, solTypeL);

      std::vector<std::vector<double>> Jacob, JacI;
      double weight;
      femL->GetJacobianMatrix(xn, intCloud.GetCloudBaricenterInParentElement(iel), weight, Jacob, JacI);
      std::vector<double> a;
      double d;
      intCloud.GetLinearFit(iel, Jacob, a, d);

      d = -d;
      for(unsigned k = 0; k < dim; k++) a[k] = -a[k];

      std::vector <TypeIO> weightCF;
      if(ielType == 3) {
        quad.GetWeightWithMap(0, a, d, weightCF);
      }
      else if(ielType == 4) {
        tri.GetWeightWithMap(0, a, d, weightCF);
      }
      else {
        abort();
      }


      double area = 0.;
      double areaC = 0.;
      for(unsigned ig = 0; ig < femL->GetGaussPointNumber(); ig++) {
        femL->Jacobian(xn, ig, weight, phi, phi_x);
        area += weight;
        areaC += weight * weightCF[ig];
      }

      _sol->_Sol[SolCIndex]->set(iel, areaC / area);

      for(unsigned i = 0; i < nDofs; i++) {

        std::vector<double> xp(dim);
        for(unsigned k = 0; k < dim; k++) {
          xp[k] = xn[k][i];
        }
        double value = intCloud.GetValue(iel, xp);

        if(iel == 4964) {
          std::cout << "AAAAA ";
          std::cout << xp[0] << " " << xp[1] << " " << value << std::endl;
        }

        unsigned inode = msh->GetSolutionDof(i, iel, solType);
        _sol->_Sol[SolCnIndex]->set(inode, (value < 0) ? 0.75 : 0.25);
      }
    }

    _sol->_Sol[SolCnIndex]->close();



    const unsigned &nel = offsetp1 - offset;
    _ypNew.resize(nel, std::vector<double>(dim));
    _yiNew.resize(nel, std::vector<double>(dim));
    _NNew.resize(nel, std::vector<double> (dim));
    _kappaNew.resize(nel);
    _dsNew.resize(nel);
    _elem.resize(nel);

    unsigned cnt = 0;

    for(_itElMrkIdx = _elMrkIdx.begin(); _itElMrkIdx != _elMrkIdx.end(); _itElMrkIdx++) {
      unsigned iel = _itElMrkIdx->first;

      if((*_sol->_Sol[SolCIndex])(iel) == 0 || (*_sol->_Sol[SolCIndex])(iel) == 1.) {

        unsigned ielType = msh->GetElementType(iel);
        unsigned nDofsL = msh->GetElementDofNumber(iel, 0);  // number of coordinate linear element dofs

        std::vector<double> x(dim);
        std::vector<double> xm(dim, 0.);

        for(unsigned i = 0; i < nDofsL; i++) {
          unsigned xDof  = msh->GetSolutionDof(i, iel, 2);
          for(unsigned k = 0; k < dim; k++) {
            x[k] = (*msh->_topology->_Sol[k])(xDof);
            xm[k] += x[k];
          }
        }

        for(unsigned  k = 0; k < dim; k++) {
          _ypNew[cnt][k] = xm[k] / nDofsL;
        }
        _yiNew[cnt]  = _yig[ielType]; //TODO
        _elem[cnt] = iel;

        _NNew[cnt] = std::vector<double>(dim, 0.); //TODO
        _kappaNew[cnt] = 0.; //TODO
        _dsNew[cnt] = 0.; //TODO

        cnt++;
        _sol->_Sol[SolCIndex]->set(iel, 1.);

        unsigned nDofs = msh->GetElementDofNumber(iel, solType);
        for(unsigned i = 0; i < nDofs; i++) {
          unsigned inode = msh->GetSolutionDof(i, iel, solType);
          if((*_sol->_Sol[SolCnIndex])(inode) == 0.) {
            _sol->_Sol[SolCnIndex]->set(inode, 1.);
          }
        }

        for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
          int jel = msh->el->GetFaceElementIndex(iel, iface) - 1; // porcata ma fallo cosi' se negativo e' un boundary
          if(jel >= 0) {
            unsigned jproc = msh->IsdomBisectionSearch(jel, 3);
            if(jproc == iproc && ((*_sol->_Sol[SolCIndex])(jel) == 0. && GetNumberOfMarker(jel) == 0)) {
              unsigned jelType = msh->GetElementType(jel);
              unsigned nDofsL = msh->GetElementDofNumber(jel, 0);
              std::vector<double> xj(dim);
              std::vector<double> xmj(dim, 0.);
              for(unsigned j = 0; j < nDofsL; j++) {
                unsigned xDofj  = msh->GetSolutionDof(j, jel, 2);
                for(unsigned k = 0; k < dim; k++) {
                  xj[k] = (*msh->_topology->_Sol[k])(xDofj);
                  xmj[k] += xj[k];
                }
              }

              for(unsigned  k = 0; k < dim; k++) {
                _ypNew[cnt][k] = xmj[k] / nDofsL;
              }
              _yiNew[cnt]  = _yig[jelType]; //TODO
              _elem[cnt] = jel;

              _NNew[cnt] = std::vector<double>(dim, 0.); //TODO
              _kappaNew[cnt] = 0.;   //TODO
              _dsNew[cnt] = 0.;   //TODO

              cnt++;
              _sol->_Sol[SolCIndex]->set(jel, 1.);
              unsigned nDofs = msh->GetElementDofNumber(jel, solType);
              for(unsigned j = 0; j < nDofs; j++) {
                unsigned jnode = msh->GetSolutionDof(j, jel, solType);
                if((*_sol->_Sol[SolCnIndex])(jnode) == 0.) {
                  _sol->_Sol[SolCnIndex]->set(jnode, 1.);
                }
              }
            }
          }
        }
      }
    }

    _sol->_Sol[SolCnIndex]->close();

    unsigned newElemNumber = 1;
    while(newElemNumber != 0) {
      unsigned newElemNumberLocal = 0;
      for(unsigned iel = offset; iel < offsetp1; iel++) {
        if((*_sol->_Sol[SolCIndex])(iel) == 0) {

          bool atLeastOneOne = false; // C is zero, but at least one of its nodes is 1
          double allSurrounded = 1.;  // C is zero, but all of its nodes are either 0 or 1

          unsigned nDofs = msh->GetElementDofNumber(iel, solType);
          for(unsigned i = 0; i < nDofs - 1; i++) {
            unsigned inode = msh->GetSolutionDof(i, iel, solType);
            if((*_sol->_Sol[SolCnIndex])(inode) == 1.) {
              atLeastOneOne = true;
              break;
            }
            allSurrounded *= ((*_sol->_Sol[SolCnIndex])(inode) > 0.5) ? 1. : 0.;
          }

          if(atLeastOneOne || allSurrounded > 1.e-10) {


            unsigned ielType = msh->GetElementType(iel);
            newElemNumberLocal++;
            unsigned nDofsL = msh->GetElementDofNumber(iel, 0);  // number of coordinate linear element dofs

            std::vector<double> x(dim);
            std::vector<double> xm(dim, 0.);

            for(unsigned i = 0; i < nDofsL; i++) {
              unsigned xDof  = msh->GetSolutionDof(i, iel, 2);
              for(unsigned k = 0; k < dim; k++) {
                x[k] = (*msh->_topology->_Sol[k])(xDof);
                xm[k] += x[k];
              }
            }

            for(unsigned  k = 0; k < dim; k++) {
              _ypNew[cnt][k] = xm[k] / nDofsL;
            }
            _yiNew[cnt]  = _yig[ielType]; //TODO
            _elem[cnt] = iel;

            _NNew[cnt] = std::vector<double>(dim, 0.); //TODO
            _kappaNew[cnt] = 0.; //TODO
            _dsNew[cnt] = 0.; //TODO
            cnt++;

            _sol->_Sol[SolCIndex]->set(iel, 1.);
            unsigned nDofs = msh->GetElementDofNumber(iel, solType);
            for(unsigned i = 0; i < nDofs; i++) {
              unsigned inode = msh->GetSolutionDof(i, iel, solType);
              if((*_sol->_Sol[SolCnIndex])(inode) == 0.) {
                _sol->_Sol[SolCnIndex]->set(inode, 1.);
              }
            }
          }
        }
      }
      MPI_Allreduce(&newElemNumberLocal, &newElemNumber, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
      _sol->_Sol[SolCnIndex]->close();

    }

    _sol->_Sol[SolCIndex]->close();

    _ypNew.resize(cnt, std::vector<double>(dim));
    _elem.resize(cnt);
    _NNew.resize(cnt, std::vector<double>(dim));
    _kappaNew.resize(cnt);
    _dsNew.resize(cnt);

    _yp.swap(_ypNew);
    _yi.swap(_yiNew);
    _kappa.swap(_kappaNew);
    _ds.swap(_dsNew);
    _N.swap(_NNew);
    CreateMap();

  }


  void Cloud::RKAdvection(const unsigned & stages, const std::vector<std::string> &U, const double & dt) {

    Mesh *msh = _sol->GetMesh();
    unsigned dim = msh->GetDimension();
    elem* el = msh->el;

    _ypNew.resize(_yp.size());
    _elemNew.resize(_yp.size());
    _NNew.resize(_yp.size());
    _kappaNew.resize(_yp.size());
    _dsNew.resize(_yp.size());


    map<unsigned, bool> pSearch;

    unsigned cnt = 0;

    std::vector < unsigned > solUIndex(dim);
    for(unsigned k = 0; k < dim; k++) {
      solUIndex[k] = _sol->GetIndex(U[k].c_str());
    }


    unsigned solType = _sol->GetSolutionType(solUIndex[0]);
    std::vector<std::vector<std::vector<double>>> solU(stages);
    std::vector<std::vector<std::vector<double>>> solUOld(stages);
    std::vector<std::vector <double>> F(stages);
    std::vector<std::vector <double>> X(stages);
    std::vector<unsigned> iel(stages);
    std::vector<unsigned> ielType(stages);
    std::vector<unsigned> nDofs(stages);

    for(unsigned j  = 0; j < stages; j++) {
      solU[j].resize(dim);
      solUOld[j].resize(dim);
      F[j].resize(dim);
      X[j].resize(dim);
    }
    std::vector<std::vector<double>> c = {{}, {}, {}, {0., 0.5, 0.5, 1.}};
    std::vector<std::vector<std::vector<double>>> a = {
      {},
      {},
      {},
      {{}, {0.5}, {0., 0.5}, {0., 0., 1.}}
    };
    std::vector<std::vector<double>> b = {{}, {}, {}, {1. / 6., 1. / 3., 1. / 3., 1. / 6.}};

    std::vector<double> phi;

    for(_itElMrkIdx = _elMrkIdx.begin(); _itElMrkIdx != _elMrkIdx.end(); _itElMrkIdx++) {
      unsigned j0 = _itElMrkIdx->second[0];
      unsigned j1 = _itElMrkIdx->second[1];

      iel[0] = _itElMrkIdx->first;
      ielType[0] = msh->GetElementType(iel[0]);
      nDofs[0] = msh->GetElementDofNumber(iel[0], solType);

      for(unsigned k = 0; k < dim; k++) {
        solU[0][k].resize(nDofs[0]);
        solUOld[0][k].resize(nDofs[0]);
      }

      for(unsigned i = 0; i < nDofs[0]; i++) {
        unsigned iDof = msh->GetSolutionDof(i, iel[0], solType);
        for(unsigned k = 0; k < dim; k++) {
          solU[0][k][i] = (*_sol->_Sol[solUIndex[k]])(iDof);
          solUOld[0][k][i] = (*_sol->_SolOld[solUIndex[k]])(iDof);
        }
      }
      for(unsigned j = j0; j < j1; j++) {
        msh->_finiteElement[ielType[0]][solType]->GetPhi(phi, _yi[_map[j]]);
        X[0] = _yp[_map[j]];
        F[0].assign(dim, 0);
        for(unsigned i = 0; i < nDofs[0]; i++) {
          for(unsigned k = 0; k < dim; k++) {
            F[0][k] += ((1. - c[stages - 1][0]) * solUOld[0][k][i] + c[stages - 1][0] * solU[0][k][i]) * phi[i];
          }
        }

        bool insideLocalDomain = true;
        for(unsigned rk = 1; rk < stages; rk++) {
          X[rk] = X[0];
          for(unsigned jk = 0; jk < rk; jk++) {
            for(unsigned k = 0; k < dim; k++) {
              X[rk][k] += a[stages - 1][rk][jk] * F[jk][k] * dt;
            }
          }
          insideLocalDomain = _mrk.SerialElementSearchWithInverseMapping(X[rk], _sol, solType, iel[rk - 1]);
          if(!insideLocalDomain) {
            pSearch[j]  = true;
            break;
          }
          iel[rk] = _mrk.GetElement();
          bool sameElement = false;
          for(unsigned jk = 0; jk < rk; jk++) {
            if(iel[rk] == iel[jk]) {
              sameElement = true;
              ielType[rk] = ielType[jk];
              nDofs[rk] =  nDofs[jk];
              solU[rk] = solU[jk];
              solUOld[rk] = solUOld[jk];
              break;
            }
          }
          if(sameElement == false) {
            ielType[rk] = msh->GetElementType(iel[rk]);
            nDofs[rk] = msh->GetElementDofNumber(iel[rk], solType);
            for(unsigned k = 0; k < dim; k++) {
              solU[rk][k].resize(nDofs[rk]);
              solUOld[rk][k].resize(nDofs[rk]);
            }
            for(unsigned i = 0; i < nDofs[rk]; i++) {
              unsigned iDof = msh->GetSolutionDof(i, iel[rk], solType);
              for(unsigned k = 0; k < dim; k++) {
                solU[rk][k][i] = (*_sol->_Sol[solUIndex[k]])(iDof);
                solUOld[rk][k][i] = (*_sol->_SolOld[solUIndex[k]])(iDof);
              }
            }
          }
          msh->_finiteElement[ielType[rk]][solType]->GetPhi(phi, _mrk.GetIprocLocalCoordinates());
          F[rk].assign(dim, 0);
          for(unsigned i = 0; i < nDofs[rk]; i++) {
            for(unsigned k = 0; k < dim; k++) {
              F[rk][k] += ((1. - c[stages - 1][rk]) * solUOld[rk][k][i] + c[stages - 1][rk] * solU[rk][k][i]) * phi[i];
            }
          }
        }
        if(insideLocalDomain) {
          _ypNew[cnt] = _yp[_map[j]];
          for(unsigned k = 0; k < dim; k++) {
            for(unsigned rk = 0; rk < stages; rk++) {
              _ypNew[cnt][k] += b[stages - 1][rk] * F[rk][k] * dt;
            }
          }
          insideLocalDomain = _mrk.SerialElementSearch(_ypNew[cnt], _sol, solType, iel[0]);
          if(insideLocalDomain) {
            _elemNew[cnt] = _mrk.GetElement();
            _NNew[cnt] = _N[_map[j]];

            _kappaNew[cnt] = _kappa[_map[j]];
            _dsNew[cnt] = _ds[_map[j]];
            cnt++;
          }
          else {
            pSearch[j]  = true;
          }
        }
      }
    }
    _ypNew.resize(cnt);
    _elemNew.resize(cnt);
    _NNew.resize(cnt);
    _kappaNew.resize(cnt);
    _dsNew.resize(cnt);

    map<unsigned, bool>::iterator it;

    unsigned iproc = _sol->processor_id();
    unsigned nprocs = _sol->n_processors();

    if(nprocs > 1) {
      for(unsigned kp = 0; kp < nprocs; kp++) {

        unsigned np;
        if(iproc == kp) {
          np = pSearch.size();
        }
        MPI_Bcast(&np, 1, MPI_UNSIGNED, kp, MPI_COMM_WORLD);

        if(np > 0) {
          if(iproc == kp) {
            it =  pSearch.begin();
          }

          for(unsigned jcnt = 0; jcnt < np; jcnt++) {
            unsigned j;
            if(iproc == kp) {
              j = it->first;

              iel[0] = _elem[_map[j]];
              ielType[0] = msh->GetElementType(iel[0]);
              nDofs[0] = msh->GetElementDofNumber(iel[0], solType);

              for(unsigned k = 0; k < dim; k++) {
                solU[0][k].resize(nDofs[0]);
                solUOld[0][k].resize(nDofs[0]);
              }

              for(unsigned i = 0; i < nDofs[0]; i++) {
                unsigned iDof = msh->GetSolutionDof(i, iel[0], solType);
                for(unsigned k = 0; k < dim; k++) {
                  solU[0][k][i] = (*_sol->_Sol[solUIndex[k]])(iDof);
                  solUOld[0][k][i] = (*_sol->_SolOld[solUIndex[k]])(iDof);
                }
              }
              msh->_finiteElement[ielType[0]][solType]->GetPhi(phi, _yi[_map[j]]);
              X[0] = _yp[_map[j]];
              F[0].assign(dim, 0);
              for(unsigned i = 0; i < nDofs[0]; i++) {
                for(unsigned k = 0; k < dim; k++) {
                  F[0][k] += ((1. - c[stages - 1][0]) * solUOld[0][k][i] + c[stages - 1][0] * solU[0][k][i]) * phi[i];
                }
              }
            }
            MPI_Bcast(&iel[0], 1, MPI_UNSIGNED, kp, MPI_COMM_WORLD);
            MPI_Bcast(X[0].data(), X[0].size(), MPI_DOUBLE, kp, MPI_COMM_WORLD);
            MPI_Bcast(F[0].data(), F[0].size(), MPI_DOUBLE, kp, MPI_COMM_WORLD);

            bool insideLocalDomain = true;
            for(unsigned rk = 1; rk < stages; rk++) {
              X[rk] = X[0];
              for(unsigned jk = 0; jk < rk; jk++) {
                for(unsigned k = 0; k < dim; k++) {
                  X[rk][k] += a[stages - 1][rk][jk] * F[jk][k] * dt;
                }
              }
              insideLocalDomain = _mrk.ParallelElementSearchWithInverseMapping(X[rk], _sol, solType, iel[rk - 1]);
              if(!insideLocalDomain) {
                break;
              }

              iel[rk] = _mrk.GetElement();
              if(_mrk.GetProc() == iproc) {
                ielType[rk] = msh->GetElementType(iel[rk]);
                nDofs[rk] = msh->GetElementDofNumber(iel[rk], solType);
                for(unsigned k = 0; k < dim; k++) {
                  solU[rk][k].resize(nDofs[rk]);
                  solUOld[rk][k].resize(nDofs[rk]);
                }
                for(unsigned i = 0; i < nDofs[rk]; i++) {
                  unsigned iDof = msh->GetSolutionDof(i, iel[rk], solType);
                  for(unsigned k = 0; k < dim; k++) {
                    solU[rk][k][i] = (*_sol->_Sol[solUIndex[k]])(iDof);
                    solUOld[rk][k][i] = (*_sol->_SolOld[solUIndex[k]])(iDof);
                  }
                }

                msh->_finiteElement[ielType[rk]][solType]->GetPhi(phi, _mrk.GetIprocLocalCoordinates());
                F[rk].assign(dim, 0);
                for(unsigned i = 0; i < nDofs[rk]; i++) {
                  for(unsigned k = 0; k < dim; k++) {
                    F[rk][k] += ((1. - c[stages - 1][rk]) * solUOld[rk][k][i] + c[stages - 1][rk] * solU[rk][k][i]) * phi[i];
                  }
                }
              }
              MPI_Bcast(F[rk].data(), F[rk].size(), MPI_DOUBLE, _mrk.GetProc(), MPI_COMM_WORLD);
            }

            if(insideLocalDomain) {
              std::vector<double> ypNew = X[0];
              for(unsigned k = 0; k < dim; k++) {
                for(unsigned rk = 0; rk < stages; rk++) {
                  ypNew[k] += b[stages - 1][rk] * F[rk][k] * dt;
                }
              }
              insideLocalDomain = _mrk.ParallelElementSearch(ypNew, _sol, solType, iel[0]);
              if(insideLocalDomain) {
                if(kp == iproc) {
                  MPI_Send(_N[_map[j]].data(), _N[_map[j]].size(), MPI_DOUBLE, _mrk.GetProc(), 1, MPI_COMM_WORLD);
                  MPI_Send(&_kappa[_map[j]], 1, MPI_DOUBLE, _mrk.GetProc(), 2, MPI_COMM_WORLD);
                  MPI_Send(&_ds[_map[j]], 1, MPI_DOUBLE, _mrk.GetProc(), 3, MPI_COMM_WORLD);
                }

                if(_mrk.GetProc() == iproc) {
                  _ypNew.resize(cnt + 1);
                  _ypNew[cnt] = ypNew;
                  _elemNew.resize(cnt + 1);
                  _elemNew[cnt] = _mrk.GetElement();
                  _NNew.resize(cnt + 1, std::vector<double> (dim));
                  MPI_Recv(_NNew[cnt].data(), _NNew[cnt].size(), MPI_DOUBLE, kp, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                  _kappaNew.resize(cnt + 1);
                  MPI_Recv(&_kappaNew[cnt], 1, MPI_DOUBLE, kp, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                  _dsNew.resize(cnt + 1);
                  MPI_Recv(&_dsNew[cnt], 1, MPI_DOUBLE, kp, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                  cnt++;
                }
              }
            }
            if(iproc == kp) {
              it++;
            }
          }
        }
      }
    }




    _yp.swap(_ypNew);
    _elem.swap(_elemNew);
    _kappa.swap(_kappaNew);
    _ds.swap(_dsNew);
    _N.swap(_NNew);

    _yi.assign(cnt, std::vector<double>(dim, 0));

    CreateMap();

  }

  void Cloud::GetLinearFit(const unsigned & iel, const std::vector<std::vector<double>> &Jac, std::vector < double > &a, double & d) {
    unsigned dim = _sol->GetMesh()->GetDimension();
    a.resize(dim);

    if(_elMrkIdx.find(iel) != _elMrkIdx.end()) {
      unsigned i0 = _elMrkIdx[iel][0];
      unsigned i1 = _elMrkIdx[iel][1];
      std::vector<double> N(dim, 0.);
      std::vector< std::vector<double>> xl(i1 - i0, std::vector<double>(dim));
      unsigned cnt;
      for(unsigned i = i0, cnt = 0; i < i1; i++, cnt++) {
        xl[cnt] = _yi[_map[i]];
        for(unsigned k = 0; k < dim; k++) {
          for(unsigned j = 0; j < dim; j++) {
            N[k] += Jac[j][k] * _N[_map[i]][j];
          }
        }
      }
      double det = 0.;
      for(unsigned k = 0; k < dim; k++) det += N[k] * N[k];
      for(unsigned k = 0; k < dim; k++) N[k] /= det;

      if(xl.size() > 1) {
        FindBestFit(xl, boost::none, N, a, d);
      }
      else if(xl.size() == 1) {
        a = N;
        d = - a[0] * xl[0][0] - a[1] * xl[0][1];
      }
    }
    else {
      std::cerr << "In function Cloud::GetGetLinearFit, this element has no marker!!!!!!\n";
      abort();

    }
  }

  void Cloud::BuildColorFunction(const char C) {
    Mesh *msh = _sol->GetMesh();
    unsigned SolCIndex = _sol->GetIndex(&C);
    unsigned solType = _sol->GetSolutionType(SolCIndex);

    unsigned iproc  = msh->processor_id();
    unsigned nprocs  = msh->n_processors();

    unsigned offset = msh->_elementOffset[iproc];
    unsigned offsetp1 = msh->_elementOffset[iproc + 1];

    for(unsigned iel = offset; iel < offsetp1 - offset; iel++) {
      unsigned nDofs = msh->GetElementDofNumber(iel, solType);
      for(unsigned i = 0; i < nDofs; i++) {
        unsigned solDof = msh->GetSolutionDof(i, iel, solType);
//         (*_sol->_Sol[SolCIndex])(solDof);
      }
    }
  }



  double Cloud::GetCost(const std::vector<std::vector<double>>&x,  const std::vector<double> &y, const std::vector<double>&w, const unsigned & iel, const unsigned & nPoints) {
    double cost = 0;
    double h2 = 0.;
    Mesh* msh = _sol->GetMesh();
    unsigned dim = _sol->GetMesh()->GetDimension();

    unsigned xDof0  = msh->GetSolutionDof(0, iel, 2);
    unsigned xDof2  = msh->GetSolutionDof(2, iel, 2);
    for(unsigned k = 0; k < dim; k++) {
      double h = (*msh->_topology->_Sol[k])(xDof2) - (*msh->_topology->_Sol[k])(xDof0);
      h2 += h * h;
    }

    const std::vector<double> &Cf = _A[iel];
    std::vector<std::vector<double>> xe(2, std::vector<double>(dim));

    for(unsigned i = 0; i < nPoints; i++) {

      const double &x0 = x[i][0];
      const double &y0 = x[i][1];

      std::vector<double> v = GetNormal(iel, x[i]);
      bool test = false;
    vct:

      oct a = Cf[0] * v[0] * v[0] + Cf[1] * v[0] * v[1] + Cf[2] * v[1] * v[1];
      oct b = 2 * Cf[0] * v[0] * x0 + Cf[1] * v[1] * x0 + Cf[1] * v[0] * y0 + 2 * Cf[2] * v[1] * y0 + Cf[3] * v[0] + Cf[4] * v[1];
      oct c = Cf[0] * x0 * x0 + Cf[1] * x0 * y0 + Cf[2] * y0 * y0 + Cf[3] * x0 + Cf[4] * y0 + Cf[5];

      oct norm = sqrt(a * a + b * b + c * c);
      a /= norm;
      b /= norm;
      c /= norm;

      unsigned cnt = 0;
      if(fabs(a) > 1.e-8) {
        oct delta = b * b - 4 * a * c;
        if(delta > 0) {
          for(unsigned j = 0; j < 2; j++) {
            double t = static_cast<double>((- b + pow(-1, j) * sqrt(delta)) / (2 * a));
            //if(t >= 0 && t <= 1) {
            for(unsigned  k = 0; k < dim; k++) {
              xe[cnt][k] = x[i][k]  + t * v[k];
            }
            cnt++;
            //}
          }
        }
      }
      else if(b != 0) {
        double t = static_cast<double>(-c / b);
        //if(t >= 0 && t <= 1) {
        for(unsigned  k = 0; k < dim; k++) {
          xe[cnt][k] = x[i][k]  + t * v[k];
        }
        cnt++;
        //}
      }

      if(cnt == 0) {
        cost += h2;
      }
      else if(cnt == 1) {
        double cost1 = 0.;
        for(unsigned k = 0; k < dim; k++) {
          cost1 += (xe[0][k] - x[i][k]) * (xe[0][k] - x[i][k]);
        }
        cost += /*(45.5 * y[i] * y[i] - 49.5 * y[i] + 5) * w[i] **/ cost1;
        if(y[i] < -0.25) cost += h2;
      }
      else {
        double cost1 = 0.;
        double cost2 = 0.;

        for(unsigned k = 0; k < dim; k++) {
          cost1 += (xe[0][k] - x[i][k]) * (xe[0][k] - x[i][k]);
          cost2 += (xe[1][k] - x[i][k]) * (xe[1][k] - x[i][k]);
        }
        cost += /*(45.5 * y[i] * y[i] - 49.5 * y[i] + 5)*/ /** w[i] * */ ((cost1 < cost2) ? cost1 : cost2);
        if(y[i] < -0.25) cost += h2;
      }
    }
    cost /= (nPoints * h2);


    return cost;
  }




} // end namespace femus


#endif













