#ifndef __femus_Cloud_hpp__
#define __femus_Cloud_hpp__

#include <fstream>
#include <iostream>     // std::cout, std::ios
#include <sstream>      // std::ostringstream

#include "MultiLevelSolution.hpp"
#include "../MyMarker/MyMarker.hpp"

namespace femus {

class Cloud {
  public: 
    Cloud() {
      _mrk = MyMarker();
    };  
    ~Cloud(){};
    void SetNumberOfMarker( const unsigned &nMax );
    void InitCircle(std::vector<double> &xc, const double &R, const unsigned &nMax, Solution* sol);
    void InitEllipse(std::vector<double> &xc, const double &a, const double &b, const unsigned &nMax, Solution* sol);
    void CreateMap();
    bool ElementSearchWithGuess(Solution* sol, const std::vector<double> &xp, const unsigned previousElem);
    void PrintNoOrder(Solution* sol, const unsigned &dim);
    void PrintWithOrder(Solution* sol, const unsigned &dim);
    void PrintCSV(Solution* sol, const unsigned &dim, const unsigned &t);
  
  private:
    unsigned _nMrk;
    double _a;
    double _b;
    std::ofstream _fout;
    std::vector<std::vector<double>> _yp;
    std::vector<std::vector<double>> _N;
    std::vector<double> _kappa;
    std::vector<std::vector<double>> _yi;
    std::vector<unsigned> _elem;
    std::vector<unsigned> _map;
    MyMarker _mrk; 
    
};

void Cloud::SetNumberOfMarker( const unsigned &nMax ){
    _nMrk = nMax;
}

void Cloud::InitEllipse(std::vector<double> &xc, const double &a, const double &b, const unsigned &nMax, Solution* sol){
  SetNumberOfMarker( nMax );
  _a = a;
  _b = b;
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
    xp[0] = xc[0] + _a * cos(t);
    xp[1] = xc[1] + _b * sin(t);
    
    elemSearch = ElementSearchWithGuess(sol, xp, previousElem);
    if(elemSearch) {
      iel = _mrk.GetElement(); 
      if(_mrk.GetProc() == iproc) {
        double NNorm = sqrt(_a * _a * cos(t) * cos(t) + _b * _b * sin(t) * sin(t));
        _yp[cnt] = xp;
        _yi[cnt]  = _mrk.GetIprocLocalCoordinates();
        _N[cnt] = {_a * cos(t) / NNorm, _b * sin(t) / NNorm};
        _kappa[cnt] = _a * _b / (pow(sqrt(_a * _a * sin(t) * sin(t) + _b * _b * cos(t) * cos(t)), 3));
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
}

void Cloud::InitCircle(std::vector<double> &xc, const double &R, const unsigned &nMax, Solution* sol)  {
  InitEllipse(xc, R, R, nMax, sol);
}

void Cloud::CreateMap(){
std::vector<unsigned*> vec(_elem.size());
  for(unsigned i = 0; i < _elem.size(); i++) {
    vec[i] = &_elem[i];
  }
  std::sort(vec.begin(), vec.end(), [](const unsigned * a, const unsigned * b) { return *a < *b;});
  _map.resize(_elem.size());
  for(unsigned i = 0; i < _map.size(); i++) {
    _map[i] =  static_cast<unsigned>(vec[i] - &_elem[0]);
  }
}

bool Cloud::ElementSearchWithGuess(Solution* sol, const std::vector<double> &xp, const unsigned previousElem = UINT_MAX ){
  bool elemSearch;
  if(previousElem == UINT_MAX) elemSearch = _mrk.ParallelElementSearchWithInverseMapping(xp, sol, 2);
  else elemSearch = _mrk.ParallelElementSearchWithInverseMapping(xp, sol, 2, previousElem);    
  
  return elemSearch;
}

void Cloud::PrintNoOrder(Solution* sol, const unsigned &dim){
  unsigned iproc = sol->processor_id();
  unsigned nprocs = sol->n_processors();
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

void Cloud::PrintWithOrder(Solution* sol, const unsigned &dim){
  unsigned iproc = sol->processor_id();
  unsigned nprocs = sol->n_processors();
  CreateMap();
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

void Cloud::PrintCSV(Solution* sol, const unsigned &dim, const unsigned &t){
  unsigned iproc = sol->processor_id();
  unsigned nprocs = sol->n_processors();
  CreateMap();
  for(unsigned kp = 0; kp < nprocs; kp++) {
    if(kp == iproc) {
      std::ostringstream foo (std::ostringstream::ate);  
      foo.str("./output/marker");
      foo << t;
      foo << ".csv"; 
      if(kp == 0) _fout.open(foo.str(), std::fstream::out);
      else _fout.open(foo.str(), std::fstream::app);

      if(kp == 0){
        _fout<<"\"X\",\"Y\",\"Z\",\"xi\",\"eta\",\"zeta\",\"Nx\",\"Ny\",\"Nz\",\"kappa\",\"ipoc\",\"elem\""<<std::endl;
      }
      for(unsigned i = 0; i < _yp.size(); i++) {
        for(unsigned k = 0; k < dim; k++) {
          _fout << _yp[_map[i]][k] << ",";
        }
        _fout<<"0.,";
        for(unsigned k = 0; k < dim; k++) {
          _fout << _yi[_map[i]][k] << ",";
        }
        _fout<<"0.,";
        for(unsigned k = 0; k < dim; k++) {
          _fout << _N[_map[i]][k] << ",";
        }
        _fout<<"0.,";
        _fout << _kappa[_map[i]] << "," << iproc << "," << _elem[_map[i]] << std::endl;
      }
      _fout.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}


}

#endif
