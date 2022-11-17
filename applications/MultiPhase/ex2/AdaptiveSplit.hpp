#ifndef __femus_AdaptiveSplit_hpp__
#define __femus_AdaptiveSplit_hpp__

#include "CutFemWeight.hpp"
#include "CDWeights.hpp"
#include "Fem.hpp"

namespace femus {

  class AdaptiveSplit {
    public:

      AdaptiveSplit(const unsigned &qOrder) {
        _quad = new CutFemWeight<double, cpp_bin_float_oct>(QUAD, qOrder, "legendre");
        _tri = new CutFemWeight<double, cpp_bin_float_oct>(TRI, qOrder, "legendre");
        _femFine = new Fem(quad.GetGaussQuadratureOrder(), quad.GetDimension());
        _femCoarse = new Fem(qOrder, quad.GetDimension());
      }

      ~AdaptiveSplit() {
        delete _quad;
        delete _tri;
        delete _femCoarse;
        delete _femFine;
      }

      void Split(const std::vector<std::vector<double>> &xv, const unsigned ielType, const std::vector<std::vector<double>> &JacI, const unsigned &level, const unsigned &father, const unsigned &granFather = 0);
      const std::vector<double>& GetWeight1() {
        return _weight1;
      };
      const std::vector<double>& GetWeight2() {
        return _weight2;
      };
      const std::vector<double>& GetWeightI() {
        return _weightI;
      };

      std::vector<std::vector<double>> &GetXpFather() {
        return _xp0;
      }

      std::vector<std::vector<double>> &GetNpFather() {
        return _Np0;
      }

      std::vector<std::vector<double>> &GetNiFather() {
        if(_Ni.size() == 0) _Ni.resize(1);
        if(_Ni[0].size() == 0) _Ni[0].resize(1);
        return _Ni[0][0];
      }

      std::vector<std::vector<double>> &GetXiFather() {
        if(_xi.size() == 0) _xi.resize(1);
        if(_xi[0].size() == 0) _xi[0].resize(1);
        return _xi[0][0];
      }


    private:

      std::vector<std::vector<double>> _xp0;
      std::vector<std::vector<double>> _Np0;
      std::vector<std::vector<double>> _Ni0;

      std::vector<std::vector<std::vector<std::vector<double>>>> _xi;
      std::vector<std::vector<std::vector<std::vector<double>>>> _Ni;

      std::vector<std::vector<std::vector<unsigned>>> _map;

      CutFemWeight <double, cpp_bin_float_oct> *_quad;
      CutFemWeight <double, cpp_bin_float_oct> *_tri;
      Fem *_femFine;
      Fem *_femCoarse;

      std::vector <double> _weight1;
      std::vector <double> _weight2;
      std::vector <double> _weightI;

      std::vector <double> _weightCut1;
      std::vector <double> _weightCut2;
      std::vector <double> _weightCutI;

      std::vector <double> _phi;
      std::vector <double> _dphidx;
  };

  void AdaptiveSplit::Split(const std::vector<std::vector<double>> &xv, const unsigned ielType, const std::vector<std::vector<double>> &JacI, const unsigned &level, const unsigned &father, const unsigned &granFather) {

    const unsigned &dim = xv.size();
    if(level == 0) {

      _weight1.resize(0);
      _weight2.resize(0);
      _weightI.resize(0);

      if(_map.size() < 1) _map.resize(1);
      _map[0].assign(1, std::vector<unsigned> (_xp0.size()));
      for(unsigned j = 0; j < _xp0.size(); j++) _map[0][0][j] = j;
    }

    const unsigned &np = _map[level][father].size();
    if(np > 2) {

      //std::cout << level << " " << child << std::endl;

      unsigned nChilds = (dim == 2) ? 4 : 8;
      std::vector<double> j(nChilds, 0);

      if(_xi.size() < level + 2) _xi.resize(level + 2);
      _xi[level + 1].assign(nChilds, std::vector<std::vector<double>>(np, std::vector<double>(dim)));
      std::vector<std::vector<double>> &xi = _xi[level][father];
      std::vector<std::vector<double>>&xi0 = _xi[level + 1][0];
      std::vector<std::vector<double>> &xi1 = _xi[level + 1][1];
      std::vector<std::vector<double>> &xi2 = _xi[level + 1][2];
      std::vector<std::vector<double>> &xi3 = _xi[level + 1][3];

      if(_Ni.size() < level + 2) _Ni.resize(level + 2);
      _Ni[level + 1].assign(nChilds, std::vector<std::vector<double>>(np, std::vector<double>(dim)));
      std::vector<std::vector<double>> &Ni = _Ni[level][father];
      std::vector<std::vector<double>> &N0 = _Ni[level + 1][0];
      std::vector<std::vector<double>> &N1 = _Ni[level + 1][1];
      std::vector<std::vector<double>> &N2 = _Ni[level + 1][2];
      std::vector<std::vector<double>> &N3 = _Ni[level + 1][3];

      if(_map.size() < level + 2) _map.resize(level + 2);
      _map[level + 1].assign(nChilds, std::vector<unsigned> (_map[level][father].size()));

      std::vector<unsigned> &map = _map[level][father];
      std::vector<unsigned> &map0 = _map[level + 1][0];
      std::vector<unsigned> &map1 = _map[level + 1][1];
      std::vector<unsigned> &map2 = _map[level + 1][2];
      std::vector<unsigned> &map3 = _map[level + 1][3];


      //std::vector<std::vector<std::vector<double>>> xij(nChilds, std::vector<std::vector<double>>(xp.size(), std::vector<double>(dim)));

      if(ielType == 3) {
        for(unsigned i = 0; i < np; i++) {
          if(xi[i][0] < 0) {
            if(xi[i][1] < 0) {
              N0[j[0]] = Ni[i];
              xi0[j[0]][0] = -1 + 2 * (xi[i][0] + 1);
              xi0[j[0]][1] = -1 + 2 * (xi[i][1] + 1);
              map0[j[0]] = map[i];

              j[0]++;
            }
            else {
              N3[j[3]] = Ni[i];
              xi3[j[3]][0] = -1 + 2 * (xi[i][0] + 1);
              xi3[j[3]][1] = -1 + 2 * xi[i][1];
              map3[j[3]] = map[i];

              j[3]++;
            }
          }
          else {
            if(xi[i][1] < 0) {
              N1[j[1]] = Ni[i];
              xi1[j[1]][0] = -1 + 2 * xi[i][0];
              xi1[j[1]][1] = -1 + 2 * (xi[i][1] + 1);
              map1[j[1]] = map[i];

              j[1]++;
            }
            else {
              N2[j[2]] = Ni[i];
              xi2[j[2]][0] = -1 + 2 * xi[i][0];
              xi2[j[2]][1] = -1 + 2 * xi[i][1];
              map2[j[2]] = map[i];
              j[2]++;
            }

          }
        }
      }

      else if(ielType == 4) {
        for(unsigned i = 0; i < np; i++) {
          if(xi[i][0] > 0.5) {
            N1[j[1]] = Ni[i];
            xi1[j[1]][0] = 2 * (xi[i][0] - 0.5);
            xi1[j[1]][1] = 2 * xi[i][1];
            map1[j[1]] = map[i];
            j[1]++;
          }
          else if(xi[i][1] > 0.5) {
            N2[j[2]] = Ni[i];
            xi2[j[2]][0] = 2 * xi[i][0];
            xi2[j[2]][1] = 2 * (xi[i][1] - 0.5);
            map2[j[2]] = map[i];
            j[2]++;
          }
          else if(xi[i][0] + xi[i][1] < 0.5) {

            N0[j[0]] = Ni[i];
            xi0[j[0]][0] = 2 * xi[i][0];
            xi0[j[0]][1] = 2 * xi[i][1];
            map0[j[0]] = map[i];
            j[0]++;
          }
          else {
            N3[j[3]][0] = -Ni[i][0];
            N3[j[3]][1] = -Ni[i][1];
            xi3[j[3]][0] = 1 - 2 * xi[i][0];
            xi3[j[3]][1] = 1 - 2 * xi[i][1];
            map3[j[3]] = map[i];
            j[3]++;
          }
        }
      }

      for(unsigned l = 0; l < nChilds; l++) {

        _Ni[level + 1][l].resize(j[l]);
        _xi[level + 1][l].resize(j[l]);
        _map[level + 1][l].resize(j[l]);
        unsigned nve = (ielType == 3) ? 4 : 3;
        std::vector<std::vector<double> > xvj(dim, std::vector<double>(nve));


        xvj.assign(dim, std::vector<double>(nve, 0.));
        for(unsigned k = 0; k < dim; k++) {
          for(unsigned I = 0; I < nve; I++) {
            for(unsigned J = 0 ; J < nve; J++) {
              xvj[k][I] += PJ[ielType][l][I][J] * xv[k][J];
            }
          }
        }

        this->Split(xvj, ielType, JacI, level + 1, l, father);
      }
    }
    else if(np > 0) {
//       unsigned nve = (ielType == 3) ? 4 : 3;
//
//       double xmin = 100;
//       double xmax = -100;
//       double ymin = 100;
//       double ymax = -100;
//
//
//       for(unsigned i = 0; i < nve; i++) {
//         // std::cout << xv[0][i] << " " << xv[1][i] << std::endl;
//         if(xmin > xv[0][i]) xmin = xv[0][i];
//         if(xmax < xv[0][i]) xmax = xv[0][i];
//         if(ymin > xv[1][i]) ymin = xv[1][i];
//         if(ymax < xv[1][i]) ymax = xv[1][i];
//       }
//       //  std::cout << xv[0][0] << " " << xv[1][0] << std::endl;
//       //  std::cout << std::endl;
//
//
//       unsigned i0 = _map[level][father][0];
//       double x1 = _xp0[i0][0];
//       double y1 = _xp0[i0][1];
//
//       double N1 = _Np0[i0][0];
//       double N2 = _Np0[i0][1];
//
//       if(x1 < xmin || x1 > xmax) std::cout << "AAAAAAAAAAAAAA " << i0 << " " << father << " " << xmin <<  " " << x1 << " " << xmax << std::endl;
//       if(y1 < ymin || y1 > ymax) std::cout << "BBBBBBBBBBBBBB " << i0 << " " << father << " " << ymin <<  " " << y1 << " " << ymax << std::endl;


      std::vector <double> a(dim);
      std::vector <double> xi(dim);
      double deta = 0.;
      for(unsigned k = 0; k < dim; k++) {
        if(np == 1) {
          a[k] = _Ni[level][father][0][k];
          xi[k] = _xi[level][father][0][k];
        }
        else {
          a[k] = 0.5 * (_Ni[level][father][0][k] + _Ni[level][father][1][k]);
          xi[k] = 0.5 * (_xi[level][father][0][k] + _xi[level][father][1][k]);
          deta += a[k] * a[k];
        }
      }
      deta /= sqrt(deta);
      if(np == 2) for(unsigned k = 0; k < dim; k++) a[k] /= deta;

      double d = 0;
      for(unsigned k = 0; k < dim; k++) d -= a[k] * xi[k];

      double dsN = 0.;
      std::vector <double> Np(dim);
      for(unsigned k = 0; k < dim; k++) {
        for(unsigned j = 0; j < dim; j++) {
          Np[k] += JacI[j][k] * a[j];
        }
        dsN += Np[k] * Np[k];
      }
      dsN = sqrt(dsN);



      CutFemWeight <double, cpp_bin_float_oct> *cutElemeFem = (ielType == 3) ? _quad : _tri;
      cutElemeFem->GetWeightWithMap(0, a, d, _weightCut2);
      for(unsigned k = 0; k < dim; k++) a[k] = - a[k];
      d = -d;
      cutElemeFem->GetWeightWithMap(0, a, d, _weightCut1);
      cutElemeFem->GetWeightWithMap(-1, a, d, _weightCutI);

      const elem_type *elemFem = _femFine->GetFiniteElement(ielType, 0);
      unsigned ng = elemFem->GetGaussPointNumber();
      unsigned size0 = _weight1.size();
      _weight1.resize(size0 + ng, 0.);
      _weight2.resize(size0 + ng, 0.);
      _weightI.resize(size0 + ng, 0.);


      for(unsigned ig = 0; ig < ng; ig++) {
        double weight;
        elemFem->Jacobian(xv, ig, weight, _phi, _dphidx);

        _weight1[size0 + ig] = weight * _weightCut1[ig];
        _weight2[size0 + ig] = weight * _weightCut2[ig];
        _weightI[size0 + ig] = weight * _weightCutI[ig] * dsN;

      }
    }
    else {
      unsigned nve = (ielType == 3) ? 4 : 3;


      //std::cout << level << " " << child << " empty or full\n";

      double d2min = 1.e10;
      unsigned jmin = 0;
      unsigned jjmin = 0;
      for(unsigned i = 0; i < nve; i++) {
        for(unsigned jj = 0; jj < _map[level - 1][granFather].size(); jj++) {
          unsigned j = _map[level - 1][granFather][jj];
          double d2 = 0.;
          for(unsigned k = 0; k < dim; k++) {
            d2 += (xv[k][i] - _xp0[j][k]) * (xv[k][i] - _xp0[j][k]);
          }
          if(d2 < d2min) {
            jjmin = jj;
            jmin = j;
            d2min = d2;
          }
        }
      }
      std::vector<double> a = _Np0[jmin];
      double d = 0;
      for(unsigned k = 0; k < dim; k++) d -= a[k] * _xp0[jmin][k];

      std::vector<double> dist(nve, d);
      for(unsigned i = 0; i < nve; i++) {
        for(unsigned k = 0; k < dim; k++) {
          dist[i] += a[k] * xv[k][i];
        }
      }

      unsigned i0 = 0;
      while(dist[i0] == 0) {
        i0 = (i0 + 1) / nve;
      }

      bool sameSign = true;
      for(unsigned i = i0 + 1; i < nve; i++) {
        if(dist[i] != 0.) {
          if(dist[i0] * dist[i] < 0) sameSign = false;
        }
      }


      if(sameSign) {

        const elem_type *elemFem = _femCoarse->GetFiniteElement(ielType, 0);
        unsigned ng = elemFem->GetGaussPointNumber();
        unsigned size0 = _weight1.size();
        _weight1.resize(size0 + ng, 0.);
        _weight2.resize(size0 + ng, 0.);
        _weightI.resize(size0 + ng, 0.);



        if(dist[i0] < 0) {
          for(unsigned ig = 0; ig < ng; ig++) {
            elemFem->Jacobian(xv, ig, _weight1[size0 + ig], _phi, _dphidx);
          }
        }
        else {
          for(unsigned ig = 0; ig < ng; ig++) {
            elemFem->Jacobian(xv, ig, _weight2[size0 + ig], _phi, _dphidx);
          }
        }
      }
      else { //is a cut cell

        a = _Ni[level - 1][granFather][jjmin];
        double d = 0;
        for(unsigned k = 0; k < dim; k++) d -= a[k] * _xi[level - 1][granFather][jjmin][k];




        double dsN = 0.;
        std::vector <double> Np(dim);
        for(unsigned k = 0; k < dim; k++) {
          for(unsigned j = 0; j < dim; j++) {
            Np[k] += JacI[j][k] * a[j];
          }
          dsN += Np[k] * Np[k];
        }
        dsN = sqrt(dsN);



        CutFemWeight <double, cpp_bin_float_oct> *cutElemeFem = (ielType == 3) ? _quad : _tri;
        cutElemeFem->GetWeightWithMap(0, a, d, _weightCut2);
        for(unsigned k = 0; k < dim; k++) a[k] = - a[k];
        d = -d;
        cutElemeFem->GetWeightWithMap(0, a, d, _weightCut1);
        cutElemeFem->GetWeightWithMap(-1, a, d, _weightCutI);

        const elem_type *elemFem = _femFine->GetFiniteElement(ielType, 0);
        unsigned ng = elemFem->GetGaussPointNumber();
        unsigned size0 = _weight1.size();
        _weight1.resize(size0 + ng, 0.);
        _weight2.resize(size0 + ng, 0.);
        _weightI.resize(size0 + ng, 0.);


        for(unsigned ig = 0; ig < ng; ig++) {
          double weight;
          elemFem->Jacobian(xv, ig, weight, _phi, _dphidx);

          _weight1[size0 + ig] = weight * _weightCut1[ig];
          _weight2[size0 + ig] = weight * _weightCut2[ig];
          _weightI[size0 + ig] = weight * _weightCutI[ig] * dsN;

        }
      }
    }
  }



}
#endif
