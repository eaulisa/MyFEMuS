
#ifndef __femus_CDWeights_hpp__
#define __femus_CDWeights_hpp__

#include "CutFemWeight.hpp"

template <class TypeA>
class CDWeighsQUAD {
  public:

    CDWeighsQUAD(const unsigned &qM, const double &dx, const double &dt) {

      CutFemWeight <double, TypeA> quad  = CutFemWeight<double, TypeA >(QUAD, qM, "legendre");

      _dx = dx;
      _dt = dt;

      _nx = floor(2. / _dx);
      _nt = floor(90. / _dt);
      _dx = 2. / _nx;
      _dt = 90. / _nt;
      xi0 = - 1. - dx;
      t0 = -dt;
      _nx += 3;
      _nt += 3;

      for(unsigned k = 0; k < 4; k++) {
        _weight[k].resize(_nx);
        for(unsigned i = 0; i < _nx; ++i) {
          _weight[k][i].resize(_nt);
        }
      }
      for(unsigned k = 0; k < 4; k++) {
        std::vector<double> xi(2);
        for(unsigned i = 0; i < _nx; ++i) {
          xi[0] = xi0 + i * _dx;
          xi[1] = (k % 2 == 0) ? xi[0] : -xi[0];
          for(unsigned t = 0; t < _nt; ++t) {
            std::vector<double> a = {cos((k * 90 + t0 + t * _dt) * M_PI / 180), sin((k * 90 + t0 + t * _dt) * M_PI / 180)};
            double d = -a[0] * xi[0] - a[1] * xi[1];
            quad.clear();
            quad(0, a, d, _weight[k][i][t]);
          }
        }
      }
    }

    void GetWeight(const std::vector<double> &a, const double &d, std::vector<double> &weight) {

      unsigned k;
      double xi;
      double t;
      if(a[0] >= 0) {
        if(a[1] >= 0) {
          k = 0;
          t = atan2(a[1], a[0]) / M_PI * 180;
          xi = -d / (a[0] + a[1]);
        }
        else {
          k = 3;
          t = atan2(a[1], a[0]) / M_PI * 180 + 90;
          xi = -d / (a[0] - a[1]);
        }
      }
      else {
        if(a[1] <= 0) {
          k = 2;
          t = atan2(a[1], a[0]) / M_PI * 180 + 180;
          xi = -d / (a[0] + a[1]);
        }
        else {
          k = 1;
          t = atan2(a[1], a[0]) / M_PI * 180 - 90;
          xi = -d / (a[0] - a[1]);
        }
      }

//       unsigned i1 = floor((xi - xi0) / _dx);
//       double s1 = (xi - (xi0 + i1 * _dx)) / _dx;
//       unsigned i2 = floor((t - t0) / _dt);
//       double s2 = (t - (t0 + i2 * _dt)) / _dt;

//       weight.resize(_weight[k][i1][i2].size());
//       for(unsigned j = 0; j < weight.size(); j++) {
//         weight[j] = _weight[k][i1][i2][j] * (1 - s1) * (1 - s2) +
//                     _weight[k][i1][i2 + 1][j] * (1 - s1) * s2 +
//                     _weight[k][i1 + 1][i2][j] * s1 * (1 - s2) +
//                     _weight[k][i1 + 1][i2 + 1][j] * s1 * s2 ;
//       }
//       return;

      int i1 = static_cast <int>(floor((xi - xi0) / _dx));
      double s1 = (xi - (xi0 + i1 * _dx)) / _dx;
      if(i1 < 1) {
        i1 = 1;
        s1 = 0;
      }
      else if(i1 > _nx - 3) {
        i1 = _nx - 3;
        s1 = 1;
      }


      int i2 = static_cast <int>(floor((t - t0) / _dt));
      double s2 = (t - (t0 + i2 * _dt)) / _dt;
      if(i2 < 1) {
        i2 = 1;
        s2 = 0;
      }
      else if(i2 > _nt - 3) {
        i2 = _nt - 3;
        s2 = 1;
      }



      _phi1.resize(4);
      _phi2.resize(4);

      _phi1[0] =             s1 * (s1 - 1.) * (s1 - 2.) / (-6.);
      _phi1[1] = (s1 + 1.)      * (s1 - 1.) * (s1 - 2.) / (2.);
      _phi1[2] = (s1 + 1.) * s1             * (s1 - 2.) / (-2.);
      _phi1[3] = (s1 + 1.) * s1 * (s1 - 1.)             / (6.);

      _phi2[0] =             s2 * (s2 - 1.) * (s2 - 2.) / (-6.);
      _phi2[1] = (s2 + 1.)      * (s2 - 1.) * (s2 - 2.) / (2.);
      _phi2[2] = (s2 + 1.) * s2             * (s2 - 2.) / (-2.);
      _phi2[3] = (s2 + 1.) * s2 * (s2 - 1.) / (6.);

      weight.assign(_weight[k][i1][i2].size(), 0.);

      for(unsigned i = 0; i < weight.size(); i++) {
        for(int j1 = 0; j1 < 4; j1++) {
          for(int j2 = 0; j2 < 4; j2++) {
            weight[i] += _phi1[j1] * _phi2[j2] * _weight[k][i1 + j1 - 1][i2 + j2 - 1][i];
          }
        }
      }
    }
  private:
    std::vector < std::vector<std::vector <double> > > _weight[4];
    double _dx;
    double _dt;
    unsigned _nx;
    unsigned _nt;
    double xi0;
    double t0;

    std::vector<double> _phi1;
    std::vector<double> _phi2;
};


template <class TypeA>
class CDWeighsTRI {
  public:

    CDWeighsTRI(const unsigned &qM, const double &dx, const double &dt) {

      CutFemWeight <double, TypeA> tri  = CutFemWeight<double, TypeA >(TRI, qM, "legendre");

      _dx = dx;
      _dt = dt;

      _nx = floor(1. / _dx);
      _nt = floor(90. / _dt);
      _dx = 1. / _nx;
      _dt = 90. / _nt;
      xi0 = - dx;
      t0 = -dt;
      _nx += 3;
      _nt += 3;

      for(unsigned k = 0; k < 4; k++) {
        _weight[k].resize(_nx);
        for(unsigned i = 0; i < _nx; ++i) {
          _weight[k][i].resize(_nt);
        }
      }
      for(unsigned k = 0; k < 4; k++) {
        std::vector<double> xi(2);
        for(unsigned i = 0; i < _nx; ++i) {
          xi[0] = xi0 + i * _dx;
          xi[1] = (k % 2 == 0) ? xi[0] : 1 - xi[0];
          for(unsigned t = 0; t < _nt; ++t) {
            std::vector<double> a = {cos((k * 90 + t0 + t * _dt) * M_PI / 180), sin((k * 90 + t0 + t * _dt) * M_PI / 180)};
            double d = -a[0] * xi[0] - a[1] * xi[1];
            tri.clear();
            tri(0, a, d, _weight[k][i][t]);
          }
        }
      }
    }

    void GetWeight(const std::vector<double> &a, const double &d, std::vector<double> &weight) {

      unsigned k;
      double xi;
      double t;
      if(a[0] >= 0) {
        if(a[1] >= 0) {
          k = 0;
          t = atan2(a[1], a[0]) / M_PI * 180;
          xi = -d / (a[0] + a[1]);
        }
        else {
          k = 3;
          t = atan2(a[1], a[0]) / M_PI * 180 + 90;
          xi = (-d - a[1]) / (a[0] - a[1]);
        }
      }
      else {
        if(a[1] <= 0) {
          k = 2;
          t = atan2(a[1], a[0]) / M_PI * 180 + 180;
          xi = -d / (a[0] + a[1]);
        }
        else {
          k = 1;
          t = atan2(a[1], a[0]) / M_PI * 180 - 90;
          xi = (-d - a[1]) / (a[0] - a[1]);
        }
      }

//       unsigned i1 = floor((xi - xi0) / _dx);
//       double s1 = (xi - (xi0 + i1 * _dx)) / _dx;
//       unsigned i2 = floor((t - t0) / _dt);
//       double s2 = (t - (t0 + i2 * _dt)) / _dt;

//       weight.resize(_weight[k][i1][i2].size());
//       for(unsigned j = 0; j < weight.size(); j++) {
//         weight[j] = _weight[k][i1][i2][j] * (1 - s1) * (1 - s2) +
//                     _weight[k][i1][i2 + 1][j] * (1 - s1) * s2 +
//                     _weight[k][i1 + 1][i2][j] * s1 * (1 - s2) +
//                     _weight[k][i1 + 1][i2 + 1][j] * s1 * s2 ;
//       }
//       return;

      int i1 = static_cast <int>(floor((xi - xi0) / _dx));
      double s1 = (xi - (xi0 + i1 * _dx)) / _dx;
      if(i1 < 1) {
        i1 = 1;
        s1 = 0;
      }
      else if(i1 > _nx - 3) {
        i1 = _nx - 3;
        s1 = 1;
      }


      int i2 = static_cast <int>(floor((t - t0) / _dt));
      double s2 = (t - (t0 + i2 * _dt)) / _dt;
      if(i2 < 1) {
        i2 = 1;
        s2 = 0;
      }
      else if(i2 > _nt - 3) {
        i2 = _nt - 3;
        s2 = 1;
      }
      
      _phi1.resize(4);
      _phi2.resize(4);

      _phi1[0] =             s1 * (s1 - 1.) * (s1 - 2.) / (-6.);
      _phi1[1] = (s1 + 1.)      * (s1 - 1.) * (s1 - 2.) / (2.);
      _phi1[2] = (s1 + 1.) * s1             * (s1 - 2.) / (-2.);
      _phi1[3] = (s1 + 1.) * s1 * (s1 - 1.)             / (6.);

      _phi2[0] =             s2 * (s2 - 1.) * (s2 - 2.) / (-6.);
      _phi2[1] = (s2 + 1.)      * (s2 - 1.) * (s2 - 2.) / (2.);
      _phi2[2] = (s2 + 1.) * s2             * (s2 - 2.) / (-2.);
      _phi2[3] = (s2 + 1.) * s2 * (s2 - 1.) / (6.);

      weight.assign(_weight[k][i1][i2].size(), 0.);

      for(unsigned i = 0; i < weight.size(); i++) {
        for(int j1 = 0; j1 < 4; j1++) {
          for(int j2 = 0; j2 < 4; j2++) {
            weight[i] += _phi1[j1] * _phi2[j2] * _weight[k][i1 + j1 - 1][i2 + j2 - 1][i];
          }
        }
      }
    }
  private:
    std::vector < std::vector<std::vector <double> > > _weight[4];
    double _dx;
    double _dt;
    unsigned _nx;
    unsigned _nt;
    double xi0;
    double t0;

    std::vector<double> _phi1;
    std::vector<double> _phi2;
};


#endif

