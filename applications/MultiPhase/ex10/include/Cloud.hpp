#ifndef __femus_Cloud_hpp__
#define __femus_Cloud_hpp__

#include <fstream>
#include <iostream>     // std::cout, std::ios
#include <sstream>      // std::ostringstream

#include "MultiLevelSolution.hpp"
#include "./MyMarker/MyMarker.hpp"
#include "MyEigenFunctions.hpp"

namespace femus {

  class Cloud {
    public:
      Cloud() {
        _mrk = MyMarker();
      };
      ~Cloud() {};
      void SetNumberOfMarker( const unsigned &nMax );
      void InitCircle(const std::vector<double> &xc, const double &R, const unsigned &nMax, Solution* sol);
      void InitEllipse(const std::vector<double> &xc, const std::vector<double> &a, const unsigned &nMax, Solution* sol);

      void PrintNoOrder(const unsigned &t);
      void PrintWithOrder(const unsigned &t);
      void PrintCSV(const unsigned &t);

      void ComputeQuadraticBestFit();

      const std::map<unsigned, std::vector<double>> GetQuadraticBestFitCoefficients() {
        return _A;
      }

      const std::vector<double> GetQuadraticBestFitCoefficients(const unsigned &iel) {
        if(_A.find(iel) != _A.end()) {
          return _A.at(iel);
        }
        else {
          return {};
        }
      }

    private:

      void CreateMap();
      bool ParallelElementSearch(const std::vector<double> &xp, const unsigned previousElem);

      Solution *_sol;
      unsigned _nMrk;
      std::ofstream _fout;
      std::vector<std::vector<double>> _yp;
      std::vector<std::vector<double>> _N;
      std::vector<double> _kappa;
      std::vector<std::vector<double>> _yi;
      std::vector<unsigned> _elem;
      std::vector<unsigned> _map;
      MyMarker _mrk;
      std::map<unsigned, std::vector<double>> _A;

  };

  void Cloud::SetNumberOfMarker( const unsigned &nMax ) {
    _nMrk = nMax;
  }

  void Cloud::InitEllipse(const std::vector<double> &xc, const std::vector<double> &a, const unsigned &nMax, Solution* sol) {

    _sol = sol;
    SetNumberOfMarker( nMax );
    double dt = 2. * M_PI / _nMrk;
    std::vector<double> xp(2);

    bool elemSearch;
    unsigned previousElem = UINT_MAX;
    unsigned iel;
    unsigned iproc = sol->processor_id();
    unsigned cnt = 0;

    _yp.resize(_nMrk);
    _yi.resize(_nMrk);
    _elem.resize(_nMrk);
    _N.resize(_nMrk);
    _kappa.resize(_nMrk);
    for(unsigned i = 0; i < _nMrk; i++) {
      double t = i * dt;
      xp[0] = xc[0] + a[0] * cos(t);
      xp[1] = xc[1] + a[1] * sin(t);

      elemSearch = ParallelElementSearch(xp, previousElem);
      if(elemSearch) {
        iel = _mrk.GetElement();
        if(_mrk.GetProc() == iproc) {
          double NNorm = sqrt( a[0] * a[0] * cos(t) * cos(t) + a[1] * a[1] * sin(t) * sin(t));
          _yp[cnt] = xp;
          _yi[cnt]  = _mrk.GetIprocLocalCoordinates();
          _N[cnt] = {a[0] * cos(t) / NNorm, a[1] * sin(t) / NNorm};
          _kappa[cnt] = a[0] * a[1] / (pow(sqrt(a[0] * a[0] * sin(t) * sin(t) + a[1] * a[1] * cos(t) * cos(t)), 3));
          _elem[cnt] = iel;
          cnt++;
        }
        previousElem = iel;
      }
      else {
        previousElem = UINT_MAX;
      }
    }
    _yp.resize(cnt);
    _yi.resize(cnt);
    _elem.resize(cnt);
    _N.resize(cnt);
    _kappa.resize(cnt);

    CreateMap();
  }

  void Cloud::InitCircle(const std::vector<double> &xc, const double &R, const unsigned &nMax, Solution* sol)  {
    InitEllipse(xc, {R, R}, nMax, sol);
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
  }

  bool Cloud::ParallelElementSearch(const std::vector<double> &xp, const unsigned previousElem = UINT_MAX ) {
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
          _fout << _kappa[i] << " " << iproc << " " << _elem[i] << std::endl;
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
          _fout << _kappa[_map[i]] << " " << iproc << " " << _elem[_map[i]] << std::endl;
        }
        _fout.close();
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  void Cloud::PrintCSV(const unsigned &t) {
    unsigned iproc = _sol->processor_id();
    unsigned nprocs = _sol->n_processors();
    unsigned dim = _sol->GetMesh()->GetDimension();

    for(unsigned kp = 0; kp < nprocs; kp++) {
      if(kp == iproc) {
        std::ostringstream foo (std::ostringstream::ate);
        foo.str("./output/marker");
        foo << t;
        foo << ".csv";
        if(kp == 0) _fout.open(foo.str(), std::fstream::out);
        else _fout.open(foo.str(), std::fstream::app);

        if(kp == 0) {
          _fout << "\"X\",\"Y\",\"Z\",\"xi\",\"eta\",\"zeta\",\"Nx\",\"Ny\",\"Nz\",\"kappa\",\"ipoc\",\"elem\"" << std::endl;
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
          _fout << _kappa[_map[i]] << "," << iproc << "," << _elem[_map[i]] << std::endl;
        }
        _fout.close();
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  void Cloud::ComputeQuadraticBestFit() {

    _A.clear();

    unsigned dim = _sol->GetMesh()->GetDimension();
    std::vector<std::vector<double>> coord;
    coord.reserve(_elem.size());
    std::vector<double> norm;

    coord.resize(0);
    unsigned i = 0;
    while(i < _elem.size()) {
      unsigned iel = _elem[_map[i]];
      unsigned cnt = 0;
      while (i < _elem.size() && _elem[_map[i]] == iel) {
        coord.resize(cnt + 1, std::vector<double> (dim));
        norm.assign(dim, 0);
        for( unsigned k = 0; k < dim; k++) {
          coord[cnt][k] = _yp[_map[i]][k];
          norm[k] += _N[_map[i]][k];
        }
        cnt++;
        i++;
      }
      femus::FindQuadraticBestFit(coord, boost::none, norm, _A[iel]);
    }
  }

}

#endif
