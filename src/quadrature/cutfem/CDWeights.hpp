
#ifndef __femus_CDWeights_hpp__
#define __femus_CDWeights_hpp__

#include "CutFemWeight.hpp"



template <class TypeA>
class CDWeight {
  public:
    virtual ~CDWeight() {};
    virtual void GetWeight(const std::vector<double> &a, const double & d, std::vector<double> &weight) = 0;
};

template <class TypeA>
class CDWeightQUAD :
  public CDWeight <TypeA> {
  public:

    CDWeightQUAD(const unsigned & qM, const double & dx, const double & dt) {

      CutFemWeight <double, TypeA> quad  = CutFemWeight<double, TypeA >(QUAD, qM, "legendre");

      _dx = dx;
      _dt = dt;

      _nx = floor(2. / _dx);
      _nt = floor(90. / _dt);
      _dx = 2. / _nx;
      _dt = 90. / _nt;
      xi0 = - 1. - _dx;
      t0 = -_dt;
      _nx += 3;
      _nt += 3;
      _ng = quad.GetGaussQuadraturePointNumber();

      _weight.resize(4 * _nx * _nt * _ng);
      
      std::ostringstream fileName;
      fileName << "./save/gaussQUAD_" << _nx << "_" << _nt << "_" << _ng << ".dat";
      FILE *fp;

      fp = fopen(fileName.str().c_str(), "r");
      if(fp != NULL) {
        fread(_weight.data(), sizeof(double), _weight.size(), fp);
        fclose(fp);
      }
      else {
      std::vector<double> weight;
      unsigned cnt = 0;
      for(unsigned k = 0; k < 4; k++) {
        std::vector<double> xi(2);
        for(unsigned i = 0; i < _nx; ++i) {
          xi[0] = xi0 + i * _dx;
          xi[1] = (k % 2 == 0) ? xi[0] : -xi[0];
          for(unsigned t = 0; t < _nt; ++t) {
            std::vector<double> a = {cos((k * 90 + t0 + t * _dt) * M_PI / 180), sin((k * 90 + t0 + t * _dt) * M_PI / 180)};
            double d = -a[0] * xi[0] - a[1] * xi[1];
            quad.clear();
            quad(0, a, d, weight);
            for(unsigned ig = 0; ig < _ng; ig++, cnt++) {
              _weight[cnt] = weight[ig];
            }
          }
        }
      }
      fp = fopen(fileName.str().c_str(), "w");
      fwrite(_weight.data(), sizeof(double), _weight.size(), fp);
      fclose(fp);
      }
    }

    void GetWeight(const std::vector<double> &a, const double & d, std::vector<double> &weight) {

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

      weight.assign(_ng, 0.);

      for(unsigned i = 0; i < weight.size(); i++) {
        for(int j1 = 0; j1 < 4; j1++) {
          for(int j2 = 0; j2 < 4; j2++) {
            weight[i] += _phi1[j1] * _phi2[j2] * _weight[k * (_nx * _nt * _ng) +
                                                        (i1 + j1 - 1) * (_nt * _ng) +
                                                        (i2 + j2 - 1) * _ng +
                                                        i];
          }
        }
      }
    }
  private:
    std::vector<double> _weight;
    double _dx;
    double _dt;
    unsigned _nx;
    unsigned _nt;
    unsigned _ng;
    double xi0;
    double t0;

    std::vector<double> _phi1;
    std::vector<double> _phi2;
};

template <class TypeA> class CDWeightTRI :
  public CDWeight <TypeA> {
  public:

    CDWeightTRI(const unsigned & qM, const double & dx, const double & dt) {

      CutFemWeight <double, TypeA> tri  = CutFemWeight<double, TypeA >(TRI, qM, "legendre");

      _dx = dx;
      _dt = dt;

      _nx = floor(1. / _dx);
      _nt = floor(90. / _dt);
      _dx = 1. / _nx;
      _dt = 90. / _nt;
      xi0 = - _dx;
      t0 = -_dt;
      _nx += 3;
      _nt += 3;
      _ng = tri.GetGaussQuadraturePointNumber();

      _weight.resize(4 * _nx * _nt * _ng);

      std::ostringstream fileName;
      fileName << "./save/gaussTRI_" << _nx << "_" << _nt << "_" << _ng << ".dat";
      FILE *fp;
      
      fp = fopen(fileName.str().c_str(), "r");
      if(fp != NULL) {
        fread(_weight.data(), sizeof(double), _weight.size(), fp);
        fclose(fp);
      }
      else {
      unsigned cnt = 0;
      std::vector<double> weight;
      
      for(unsigned k = 0; k < 4; k++) {
        std::vector<double> xi(2);
        for(unsigned i = 0; i < _nx; ++i) {
          xi[0] = xi0 + i * _dx;
          xi[1] = (k % 2 == 0) ? xi[0] : 1 - xi[0];
          for(unsigned t = 0; t < _nt; ++t) {
            std::vector<double> a = {cos((k * 90 + t0 + t * _dt) * M_PI / 180), sin((k * 90 + t0 + t * _dt) * M_PI / 180)};
            double d = -a[0] * xi[0] - a[1] * xi[1];
            tri.clear();
            tri(0, a, d, weight);
            for(unsigned ig = 0; ig < _ng; ig++, cnt++) {
              _weight[cnt] = weight[ig];
            }
          }
        }
      }
      fp = fopen(fileName.str().c_str(), "w");
      fwrite(_weight.data(), sizeof(double), _weight.size(), fp);
      fclose(fp);
      }
    }

    void GetWeight(const std::vector<double> &a, const double & d, std::vector<double> &weight) {

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

      weight.assign(_ng, 0.);

      for(unsigned i = 0; i < weight.size(); i++) {
        for(int j1 = 0; j1 < 4; j1++) {
          for(int j2 = 0; j2 < 4; j2++) {
            weight[i] += _phi1[j1] * _phi2[j2] * _weight[k * (_nx * _nt * _ng) +
                                                        (i1 + j1 - 1) * (_nt * _ng) +
                                                        (i2 + j2 - 1) * _ng +
                                                        i];
          }
        }
      }
    }
  private:
    std::vector<double> _weight;  
    double _dx;
    double _dt;
    unsigned _nx;
    unsigned _nt;
    unsigned _ng;
    double xi0;
    double t0;

    std::vector<double> _phi1;
    std::vector<double> _phi2;
};


template <class TypeA>
class CDWeightTET :
  public CDWeight <TypeA> {
  public:

    CDWeightTET(const unsigned & qM, const double & dx, const double & dt) {

      CutFemWeight <double, TypeA> tet  = CutFemWeight<double, TypeA >(TET, qM, "legendre");

      //tet.SetTetBaseType(0);

      _dx = dx;
      _dt = dt;
      _df = dt;

      _nx = floor(1. / _dx);
      _nt = floor(90. / _dt);
      _nf = floor(90. / _df);

      _dx = 1. / _nx;
      _dt = 90. / _nt;
      _df = 90. / _nf;

      xi0 = 0;
      t0 = 0;
      f0 = 0;
      _nx += 1;
      _nt += 1;
      _nf += 1;
      _ng = tet.GetGaussQuadraturePointNumber();

      _weight.resize(8 * _nx * _nt * _nf * _ng);


      std::ostringstream fileName;
      fileName << "./save/gaussTET_" << _nx << "_" << _nt << "_" << _nf << "_" << _ng << ".dat";
      FILE *fp;

      fp = fopen(fileName.str().c_str(), "r");
      if(fp != NULL) {
        fread(_weight.data(), sizeof(double), _weight.size(), fp);
        fclose(fp);
      }
      else {
        unsigned cnt = 0;
        std::vector<double> weight;
        for(unsigned k = 0; k < 8; k++) {
          std::vector<double> xi(3);
          for(unsigned i = 0; i < _nx; ++i) {
            xi[0] = xi0 + i * _dx;
            switch(k) {
              case 0:
              case 6:
                xi[1] = xi[0];
                xi[2] = xi[0];
                break;
              case 1:
              case 7:
                xi[1] = 1. - xi[0];
                xi[2] = 1. - xi[0];
                break;
              case 2:
              case 4:
                xi[1] = xi[0];
                xi[2] = 1. - xi[0];
                break;
              case 3:
              case 5:
                xi[1] = 1. - xi[0];
                xi[2] = xi[0];
                break;
            }

            std::cout << k << " " << i << " " << std::flush;
            for(unsigned t = 0; t < _nt; ++t) {
              for(unsigned f = 0; f < _nf; ++f) {

                std::vector<double> a = {sin((k / 4 * 90 + f0 + f * _df) * M_PI / 180.) * cos((k % 4 * 90 + t0 + t * _dt) * M_PI / 180.),
                                         sin((k / 4 * 90 + f0 + f * _df) * M_PI / 180.) * sin((k % 4 * 90 + t0 + t * _dt) * M_PI / 180.),
                                         cos((k / 4 * 90 + f0 + f * _df) * M_PI / 180.)
                                        };
                double d = -a[0] * xi[0] - a[1] * xi[1] - a[2] * xi[2];
                tet.clear();
                tet(0, a, d, weight);

                for(unsigned ig = 0; ig < _ng; ig++, cnt++) {
                  _weight[cnt] = weight[ig];
                }
              }
            }
          }
        }

        fp = fopen(fileName.str().c_str(), "w");
        fwrite(_weight.data(), sizeof(double), _weight.size(), fp);
        fclose(fp);
      }

    }

    void GetWeight(const std::vector<double> &a, const double &d, std::vector<double> &weight) {

      unsigned k;
      double xi;
      double t;
      double f = acos(a[2] / sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])) / M_PI * 180;

      if(a[0] >= 0) {
        if(a[1] >= 0) {
          t = atan2(a[1], a[0]) / M_PI * 180;
          if(a[2] >= 0) {
            k = 0;
            xi = -d / (a[0] + a[1] + a[2]);
          }
          else {
            k = 4;
            xi = (-d - a[2]) / (a[0] + a[1] - a[2]);
          }
        }
        else {
          t = atan2(a[1], a[0]) / M_PI * 180 + 90;
          if(a[2] >= 0) {
            k = 3;
            xi = (-d - a[1]) / (a[0] - a[1] + a[2]);
          }
          else {
            k = 7;
            xi = (-d - a[1] - a[2]) / (a[0] - a[1] - a[2]);
          }
        }
      }
      else {
        if(a[1] <= 0) {
          t = atan2(a[1], a[0]) / M_PI * 180 + 180;
          if(a[2] >= 0) {
            k = 2;
            xi = (-d - a[2]) / (a[0] + a[1] - a[2]);
          }
          else {
            k = 6;
            xi = -d / (a[0] + a[1] + a[2]);
          }
        }
        else {
          t = atan2(a[1], a[0]) / M_PI * 180 - 90;
          if(a[2] >= 0) {
            k = 1;
            xi = (-d - a[1] - a[2]) / (a[0] - a[1] - a[2]);
          }
          else {
            k = 5;
            xi = (-d - a[1]) / (a[0] - a[1] + a[2]);
          }
        }
      }

      if(k >= 4) f -= 90;

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

//       std::cout << k << " " << a[0] << " " << a[1] << " " << a[2] << " " << d << std::endl;
//       std::cout << xi << std::endl;
//       std::cout << t << std::endl;
//       std::cout << f << std::endl << std::endl << std::flush;

      int i1 = static_cast <int>(floor((xi - xi0) / _dx));
      double s1 = (xi - (xi0 + i1 * _dx)) / _dx;
      if(i1 < 0) {
        i1 = 1;
        s1 = - 1;
      }
      else if(i1 == 0) {
        i1 = 1;
        s1 = s1 - 1;
      }
      else if(i1 == _nx - 2) {
        i1 = _nx - 3;
        s1 = s1 + 1;
      }
      else if(i1 >= _nx - 1) {
        i1 = _nx - 3;
        s1 = 2;
      }


      int i2 = static_cast <int>(floor((t - t0) / _dt));
      double s2 = (t - (t0 + i2 * _dt)) / _dt;
      if(i2 < 0) {
        i2 = 1;
        s2 = - 1;
      }
      else if(i2 == 0) {
        i2 = 1;
        s2 = s2 - 1;
      }
      else if(i2 == _nt - 2) {
        i2 = _nt - 3;
        s2 = s2 + 1;
      }
      else if(i2 >= _nt - 1) {
        i2 = _nt - 3;
        s2 = 2;
      }

      int i3 = static_cast <int>(floor((f - f0) / _df));
      double s3 = (f - (f0 + i3 * _df)) / _df;
      if(i3 < 0) {
        i3 = 1;
        s3 = - 1;
      }
      else if(i3 == 0) {
        i3 = 1;
        s3 = s3 - 1;
      }
      else if(i3 == _nf - 2) {
        i3 = _nf - 3;
        s3 = s3 + 1;
      }
      else if(i3 >= _nf - 1) {
        i3 = _nf - 3;
        s3 = 2;
      }

//       std::cout << k << std::endl;
//       std::cout << i1 << " " << s1 << std::endl;
//       std::cout << i2 << " " << s2 << std::endl;
//       std::cout << i3 << " " << s3 << std::endl << std::endl << std::flush;

      _phi1.resize(4);
      _phi2.resize(4);
      _phi3.resize(4);

      _phi1[0] =             s1 * (s1 - 1.) * (s1 - 2.) / (-6.);
      _phi1[1] = (s1 + 1.)      * (s1 - 1.) * (s1 - 2.) / (2.);
      _phi1[2] = (s1 + 1.) * s1             * (s1 - 2.) / (-2.);
      _phi1[3] = (s1 + 1.) * s1 * (s1 - 1.)             / (6.);

      _phi2[0] =             s2 * (s2 - 1.) * (s2 - 2.) / (-6.);
      _phi2[1] = (s2 + 1.)      * (s2 - 1.) * (s2 - 2.) / (2.);
      _phi2[2] = (s2 + 1.) * s2             * (s2 - 2.) / (-2.);
      _phi2[3] = (s2 + 1.) * s2 * (s2 - 1.) / (6.);

      _phi3[0] =             s3 * (s3 - 1.) * (s3 - 2.) / (-6.);
      _phi3[1] = (s3 + 1.)      * (s3 - 1.) * (s3 - 2.) / (2.);
      _phi3[2] = (s3 + 1.) * s3             * (s3 - 2.) / (-2.);
      _phi3[3] = (s3 + 1.) * s3 * (s3 - 1.) / (6.);

      weight.assign(_ng, 0.);

      for(unsigned i = 0; i < weight.size(); i++) {
        for(int j1 = 0; j1 < 4; j1++) {
          for(int j2 = 0; j2 < 4; j2++) {
            for(int j3 = 0; j3 < 4; j3++) {
              weight[i] += _phi1[j1] * _phi2[j2] * _phi3[j3] *
                           _weight[k * (_nx * _nt * _nf * _ng) +
                                     (i1 + j1 - 1) * (_nt * _nf * _ng) +
                                     (i2 + j2 - 1) * (_nf * _ng) +
                                     (i3 + j3 - 1) * _ng +
                                     i];
            }
          }
        }
      }
    }

  private:
    std::vector < double > _weight;
    double _dx;
    double _dt;
    double _df;
    unsigned _nx;
    unsigned _nt;
    unsigned _nf;
    unsigned _ng;
    double xi0;
    double t0;
    double f0;

    std::vector<double> _phi1;
    std::vector<double> _phi2;
    std::vector<double> _phi3;
};

template <class TypeA>
class CDWeightHEX :
  public CDWeight <TypeA> {
  public:

    CDWeightHEX(const unsigned & qM, const double & dx, const double & dt) {

      CutFemWeight <double, TypeA> hex  = CutFemWeight<double, TypeA >(HEX, qM, "legendre");

      //tet.SetTetBaseType(0);

      _dx = dx;
      _dt = dt;
      _df = dt;

      _nx = floor(2. / _dx);
      _nt = floor(90. / _dt);
      _nf = floor(90. / _df);

      _dx = 2. / _nx;
      _dt = 90. / _nt;
      _df = 90. / _nf;

      xi0 = -1.;
      t0 = 0;
      f0 = 0;
      _nx += 1;
      _nt += 1;
      _nf += 1;
      _ng = hex.GetGaussQuadraturePointNumber();

      _weight.resize(8 * _nx * _nt * _nf * _ng);


      std::ostringstream fileName;
      fileName << "./save/gaussHEX_" << _nx << "_" << _nt << "_" << _nf << "_" << _ng << ".dat";
      FILE *fp;

      fp = fopen(fileName.str().c_str(), "r");
      if(fp != NULL) {
        fread(_weight.data(), sizeof(double), _weight.size(), fp);
        fclose(fp);
      }
      else {
        unsigned cnt = 0;
        std::vector<double> weight;
        for(unsigned k = 0; k < 8; k++) {
          std::vector<double> xi(3);
          for(unsigned i = 0; i < _nx; ++i) {
            xi[0] = xi0 + i * _dx;
            switch(k) {
              case 0:
              case 6:
                xi[1] = xi[0];
                xi[2] = xi[0];
                break;
              case 1:
              case 7:
                xi[1] = - xi[0];
                xi[2] = - xi[0];
                break;
              case 2:
              case 4:
                xi[1] = xi[0];
                xi[2] = - xi[0];
                break;
              case 3:
              case 5:
                xi[1] = - xi[0];
                xi[2] = xi[0];
                break;
            }

            std::cout << k << "-" << i << " " << std::flush;
            for(unsigned t = 0; t < _nt; ++t) {
              for(unsigned f = 0; f < _nf; ++f) {

                std::vector<double> a = {sin((k / 4 * 90 + f0 + f * _df) * M_PI / 180.) * cos((k % 4 * 90 + t0 + t * _dt) * M_PI / 180.),
                                         sin((k / 4 * 90 + f0 + f * _df) * M_PI / 180.) * sin((k % 4 * 90 + t0 + t * _dt) * M_PI / 180.),
                                         cos((k / 4 * 90 + f0 + f * _df) * M_PI / 180.)
                                        };
                double d = -a[0] * xi[0] - a[1] * xi[1] - a[2] * xi[2];
                hex.clear();
                hex(0, a, d, weight);

                for(unsigned ig = 0; ig < _ng; ig++, cnt++) {
                  _weight[cnt] = weight[ig];
                }
              }
            }
          }
        }

        fp = fopen(fileName.str().c_str(), "w");
        fwrite(_weight.data(), sizeof(double), _weight.size(), fp);
        fclose(fp);
      }

    }

    void GetWeight(const std::vector<double> &a, const double &d, std::vector<double> &weight) {

      unsigned k;
      double xi;
      double t;
      double f = acos(a[2] / sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])) / M_PI * 180;

      if(a[0] >= 0) {
        if(a[1] >= 0) {
          t = atan2(a[1], a[0]) / M_PI * 180;
          if(a[2] >= 0) {
            k = 0;
            xi = -d / (a[0] + a[1] + a[2]);
          }
          else {
            k = 4;
            xi = -d / (a[0] + a[1] - a[2]);
          }
        }
        else {
          t = atan2(a[1], a[0]) / M_PI * 180 + 90;
          if(a[2] >= 0) {
            k = 3;
            xi = -d / (a[0] - a[1] + a[2]);
          }
          else {
            k = 7;
            xi = -d / (a[0] - a[1] - a[2]);
          }
        }
      }
      else {
        if(a[1] <= 0) {
          t = atan2(a[1], a[0]) / M_PI * 180 + 180;
          if(a[2] >= 0) {
            k = 2;
            xi = -d / (a[0] + a[1] - a[2]);
          }
          else {
            k = 6;
            xi = -d / (a[0] + a[1] + a[2]);
          }
        }
        else {
          t = atan2(a[1], a[0]) / M_PI * 180 - 90;
          if(a[2] >= 0) {
            k = 1;
            xi = -d / (a[0] - a[1] - a[2]);
          }
          else {
            k = 5;
            xi = -d / (a[0] - a[1] + a[2]);
          }
        }
      }

      if(k >= 4) f -= 90;

      std::cout<<" AAA  " <<k<<" "<<xi<<" "<<t<<" "<<f<<" "<<std::endl;

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

//       std::cout << k << " " << a[0] << " " << a[1] << " " << a[2] << " " << d << std::endl;
//       std::cout << xi << std::endl;
//       std::cout << t << std::endl;
//       std::cout << f << std::endl << std::endl << std::flush;

      int i1 = static_cast <int>(floor((xi - xi0) / _dx));
      double s1 = (xi - (xi0 + i1 * _dx)) / _dx;
      if(i1 < 0) {
        i1 = 1;
        s1 = - 1;
      }
      else if(i1 == 0) {
        i1 = 1;
        s1 = s1 - 1;
      }
      else if(i1 == _nx - 2) {
        i1 = _nx - 3;
        s1 = s1 + 1;
      }
      else if(i1 >= _nx - 1) {
        i1 = _nx - 3;
        s1 = 2;
      }


      int i2 = static_cast <int>(floor((t - t0) / _dt));
      double s2 = (t - (t0 + i2 * _dt)) / _dt;
      if(i2 < 0) {
        i2 = 1;
        s2 = - 1;
      }
      else if(i2 == 0) {
        i2 = 1;
        s2 = s2 - 1;
      }
      else if(i2 == _nt - 2) {
        i2 = _nt - 3;
        s2 = s2 + 1;
      }
      else if(i2 >= _nt - 1) {
        i2 = _nt - 3;
        s2 = 2;
      }

      int i3 = static_cast <int>(floor((f - f0) / _df));
      double s3 = (f - (f0 + i3 * _df)) / _df;
      if(i3 < 0) {
        i3 = 1;
        s3 = - 1;
      }
      else if(i3 == 0) {
        i3 = 1;
        s3 = s3 - 1;
      }
      else if(i3 == _nf - 2) {
        i3 = _nf - 3;
        s3 = s3 + 1;
      }
      else if(i3 >= _nf - 1) {
        i3 = _nf - 3;
        s3 = 2;
      }

      std::cout << k << std::endl;
      std::cout <<_nx <<" "<<_dx<<" "<< i1 << " " << s1 << std::endl;
      std::cout <<_nt <<" "<<_dt<<" "<< i2 << " " << s2 << std::endl;
      std::cout <<_nf <<" "<<_df<<" "<< i3 << " " << s3 << std::endl << std::endl << std::flush;

      _phi1.resize(4);
      _phi2.resize(4);
      _phi3.resize(4);

      _phi1[0] =             s1 * (s1 - 1.) * (s1 - 2.) / (-6.);
      _phi1[1] = (s1 + 1.)      * (s1 - 1.) * (s1 - 2.) / (2.);
      _phi1[2] = (s1 + 1.) * s1             * (s1 - 2.) / (-2.);
      _phi1[3] = (s1 + 1.) * s1 * (s1 - 1.)             / (6.);

      _phi2[0] =             s2 * (s2 - 1.) * (s2 - 2.) / (-6.);
      _phi2[1] = (s2 + 1.)      * (s2 - 1.) * (s2 - 2.) / (2.);
      _phi2[2] = (s2 + 1.) * s2             * (s2 - 2.) / (-2.);
      _phi2[3] = (s2 + 1.) * s2 * (s2 - 1.) / (6.);

      _phi3[0] =             s3 * (s3 - 1.) * (s3 - 2.) / (-6.);
      _phi3[1] = (s3 + 1.)      * (s3 - 1.) * (s3 - 2.) / (2.);
      _phi3[2] = (s3 + 1.) * s3             * (s3 - 2.) / (-2.);
      _phi3[3] = (s3 + 1.) * s3 * (s3 - 1.) / (6.);

      weight.assign(_ng, 0.);

      for(unsigned i = 0; i < weight.size(); i++) {
        for(int j1 = 0; j1 < 4; j1++) {
          for(int j2 = 0; j2 < 4; j2++) {
            for(int j3 = 0; j3 < 4; j3++) {
              weight[i] += _phi1[j1] * _phi2[j2] * _phi3[j3] *
                           _weight[k * (_nx * _nt * _nf * _ng) +
                                     (i1 + j1 - 1) * (_nt * _nf * _ng) +
                                     (i2 + j2 - 1) * (_nf * _ng) +
                                     (i3 + j3 - 1) * _ng +
                                     i];
            }
          }
        }
      }
    }

  private:
    std::vector < double > _weight;
    double _dx;
    double _dt;
    double _df;
    unsigned _nx;
    unsigned _nt;
    unsigned _nf;
    unsigned _ng;
    double xi0;
    double t0;
    double f0;

    std::vector<double> _phi1;
    std::vector<double> _phi2;
    std::vector<double> _phi3;
};

#endif


