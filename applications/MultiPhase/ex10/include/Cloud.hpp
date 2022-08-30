#ifndef __femus_Cloud_hpp__
#define __femus_Cloud_hpp__

#include <fstream>


class Cloud {
  public: 
    Cloud() {};  
    void InitCircle(std::vector<double> &xc, double &R, unsigned &nMax);
    void SetNumberOfMarker( unsigned &nMax );
    void CreateMap();
    void PrintNoOrder(unsigned &nprocs, unsigned &iproc, unsigned &dim);
    void PrintWithOrder(unsigned &nprocs, unsigned &iproc, unsigned &dim);
  
  private:
    unsigned _nMrk;
    std::ofstream _fout;
    std::vector<std::vector<double>> _xp;
    std::vector<std::vector<double>> _N;
    std::vector<double> _kappa;
    std::vector<std::vector<double>> _yi;
    std::vector<unsigned> _elem;
    std::vector<unsigned> _map;
  
    
};

void Cloud::InitCircle(std::vector<double> &xc, double &R, unsigned &nMax)  {
    SetNumberOfMarker( nMax );
    double dt = 2. * M_PI / _nMrk;
    _xp.resize(_nMrk);
    for(unsigned i = 0; i < _nMrk; i++) {
      _xp[i].resize(2);
      _xp[i][0] = xc[0] + R * cos(i * dt);
      _xp[i][1] = xc[1] + R * sin(i * dt);
    }
}

void Cloud::SetNumberOfMarker( unsigned &nMax ){
    _nMrk = nMax;
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
  

void Cloud::PrintNoOrder(unsigned &nprocs, unsigned &iproc, unsigned &dim){
  for(unsigned kp = 0; kp < nprocs; kp++) {
    if(kp == iproc) {
      if(kp == 0) _fout.open("markerno.dat", std::fstream::out);
      else _fout.open("markerno.dat", std::fstream::app);
      for(unsigned i = 0; i < _xp.size(); i++) {
        for(unsigned k = 0; k < dim; k++) {
          _fout << _xp[i][k] << " ";
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

void Cloud::PrintWithOrder(unsigned &nprocs, unsigned &iproc, unsigned &dim){
   for(unsigned kp = 0; kp < nprocs; kp++) {
    if(kp == iproc) {
      if(kp == 0) _fout.open("marker.dat", std::fstream::out);
      else _fout.open("marker.dat", std::fstream::app);
      for(unsigned i = 0; i < _xp.size(); i++) {
        for(unsigned k = 0; k < dim; k++) {
          _fout << _xp[_map[i]][k] << " ";
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

#endif
