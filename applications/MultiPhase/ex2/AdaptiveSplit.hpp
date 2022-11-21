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

      void Split(const std::vector<std::vector<double>> &xv, const unsigned ielType, const unsigned &level, const unsigned &father);
      const std::vector<double>& GetWeight1() {
        return _weight1;
      };
      const std::vector<double>& GetWeight2() {
        return _weight2;
      };
      const std::vector<double>& GetWeightI() {
        return _weightI;
      };
      const std::vector<std::vector<double>>& GetXv() {
        return _xvig;
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

      void RemapTriangles(const std::vector<double> &xi, const std::vector<double> &Ni,
                          const unsigned &i0, const unsigned &i1, const unsigned &i2, const unsigned &i3,
                          std::vector<std::vector<std::vector<double>>> &yi, std::vector<std::vector<std::vector<double>>> &Mi) {

        Mi[0][i0] = Ni;
        yi[0][i0][0] = 2 * xi[0];
        yi[0][i0][1] = 2 * xi[1];

        Mi[1][i1] = Ni;
        yi[1][i1][0] = 2 * (xi[0] - 0.5);
        yi[1][i1][1] = 2 * xi[1];

        Mi[2][i2] = Ni;
        yi[2][i2][0] = 2 * xi[0];
        yi[2][i2][1] = 2 * (xi[1] - 0.5);

        Mi[3][i3][0] = -Ni[0];
        Mi[3][i3][1] = -Ni[1];
        yi[3][i3][0] = 1 - 2 * xi[0];
        yi[3][i3][1] = 1 - 2 * xi[1];
      }


      void RemapSquares(const std::vector<double> &xi, const std::vector<double> &Ni,
                        const unsigned &i0, const unsigned &i1, const unsigned &i2, const unsigned &i3,
                        std::vector<std::vector<std::vector<double>>> &yi, std::vector<std::vector<std::vector<double>>> &Mi) {
        Mi[0][i0] = Ni;
        yi[0][i0][0] = -1 + 2 * (xi[0] + 1);
        yi[0][i0][1] = -1 + 2 * (xi[1] + 1);

        Mi[1][i1] = Ni;
        yi[1][i1][0] =  -1 + 2 * xi[0];
        yi[1][i1][1] =  -1 + 2 * (xi[1] + 1);

        Mi[2][i2] = Ni;
        yi[2][i2][0] = -1 + 2 * xi[0];
        yi[2][i2][1] = -1 + 2 * xi[1];

        Mi[3][i3] = Ni;
        yi[3][i3][0] = -1 + 2 * (xi[0] + 1);
        yi[3][i3][1] = -1 + 2 * xi[1];
      }

    private:

      std::vector<std::vector<double>> _xp0;
      std::vector<std::vector<double>> _Np0;
      std::vector<std::vector<double>> _Ni0;

      const std::vector<std::vector<double>> _xig = {{0., 0., 0.}, {1. / 3., 1. / 3., 1. / 3.}, {1. / 3., 1. / 3., 0.}, {0., 0.}, {1. / 3., 1. / 3.}, {0.}};

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

      std::vector<std::vector <double>> _xvig;

      std::vector<std::vector<double>> _Jac;
      std::vector<std::vector<double>> _JacI;
  };

  void AdaptiveSplit::Split(const std::vector<std::vector<double>> &xv, const unsigned ielType, const unsigned &level, const unsigned &father) {

    if(level == 0) {
      _weight1.resize(0);
      _weight2.resize(0);
      _weightI.resize(0);
      _xvig.resize(0);

      if(_map.size() < 1) _map.resize(1);
      _map[0].assign(1, std::vector<unsigned> (_xp0.size()));
      for(unsigned j = 0; j < _xp0.size(); j++) _map[0][0][j] = j;
    }


    const unsigned &n0 = _xi[level][father].size();
    const unsigned &n1 = _map[level][father].size();

    const unsigned &dim = xv.size();
    if(level == 0 || (level < 1 && n1 >= 4)) {

      unsigned nChilds = (dim == 2) ? 4 : 8;
      std::vector<int> j(nChilds, 0);
      std::vector<int> k(nChilds, n1 - 1);

      if(_xi.size() < level + 2) _xi.resize(level + 2);
      _xi[level + 1].assign(nChilds, std::vector<std::vector<double>>(n1, std::vector<double>(dim)));
      std::vector<std::vector<double>> &xi = _xi[level][father];

      if(_Ni.size() < level + 2) _Ni.resize(level + 2);
      _Ni[level + 1].assign(nChilds, std::vector<std::vector<double>>(n1, std::vector<double>(dim)));
      std::vector<std::vector<double>> &Ni = _Ni[level][father];

      if(_map.size() < level + 2) _map.resize(level + 2);
      _map[level + 1].assign(nChilds, std::vector<unsigned> (_map[level][father].size()));

      std::vector<unsigned> &map = _map[level][father];

      if(ielType == 3) {
        for(unsigned i = 0; i < n1; i++) {
          if(xi[i][0] < 0) {
            if(xi[i][1] < 0) {
              _map[level + 1][0][j[0]] = map[i];
              RemapSquares(xi[i], Ni[i], j[0]++, k[1]--, k[2]--, k[3]--, _xi[level + 1], _Ni[level + 1]);
            }
            else {
              _map[level + 1][3][j[3]] = map[i];
              RemapSquares(xi[i], Ni[i], k[0]--, k[1]--, k[2]--, j[3]++, _xi[level + 1], _Ni[level + 1]);
            }
          }
          else {
            if(xi[i][1] < 0) {
              _map[level + 1][1][j[1]] = map[i];
              RemapSquares(xi[i], Ni[i], k[0]--, j[1]++, k[2]--, k[3]--, _xi[level + 1], _Ni[level + 1]);
            }
            else {
              _map[level + 1][2][j[2]] = map[i];
              RemapSquares(xi[i], Ni[i], k[0]--, k[1]--, j[2]++, k[3]--, _xi[level + 1], _Ni[level + 1]);
            }
          }
        }
      }

      else if(ielType == 4) {
        for(unsigned i = 0; i < n1; i++) {
          if(xi[i][0] > 0.5) {
            _map[level + 1][1][j[1]] = map[i];
            RemapTriangles(xi[i], Ni[i], k[0]--, j[1]++, k[2]--, k[3]--, _xi[level + 1], _Ni[level + 1]);
          }
          else if(xi[i][1] > 0.5) {
            _map[level + 1][2][j[2]] = map[i];
            RemapTriangles(xi[i], Ni[i], k[0]--, k[1]--, j[2]++, k[3]--, _xi[level + 1], _Ni[level + 1]);
          }
          else if(xi[i][0] + xi[i][1] < 0.5) {
            _map[level + 1][0][j[0]] = map[i];
            RemapTriangles(xi[i], Ni[i], j[0]++, k[1]--, k[2]--, k[3]--, _xi[level + 1], _Ni[level + 1]);
          }
          else {
            _map[level + 1][3][j[3]] = map[i];
            RemapTriangles(xi[i], Ni[i], k[0]--, k[1]--, k[2]--, j[3]++, _xi[level + 1], _Ni[level + 1]);
          }
        }
      }

      for(unsigned l = 0; l < nChilds; l++) {
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
        this->Split(xvj, ielType, level + 1, l);
      }
    }
    else {

      std::vector <double> a;
      double d;

      unsigned nve = (ielType == 3) ? 4 : 3;

      std::vector<std::vector <double>> &xi = _xi[level][father];
      std::vector<std::vector <double>> &Ni = _Ni[level][father];

      std::vector<double> xim(dim, 0);
      std::vector<double> N(dim, 0);

      if(n0 == 1) {
        a = Ni[0];
        d = -a[0] * xi[0][0] - a[1] * xi[0][1];
      }
      else if(n1 > 0) { // it is a cut cell
        double detN = 0;
        for(unsigned k = 0; k < dim; k++) {
          for(unsigned i = 0; i < n1; i++) {
            xim[k] += xi[i][k];
            N[k] += Ni[i][k];
          }
          detN += N[k] * N[k];
        }
        detN = sqrt(detN);
        for(unsigned k = 0; k < dim; k++) {
          N[k] /= detN;
          xim[k] /= n1;
        }
        std::vector<double> d2(n0, 0);
        std::vector<double> weight(n0);
        double sigma2 = 0;
        for(unsigned i = 0; i < n0; i++) {
          for(unsigned k = 0; k < dim; k++) {
            d2[i] += (xi[i][k] - xim[k]) * (xi[i][k] - xim[k]);
          }
          sigma2 += d2[i];
        }
        sigma2 /= n0 * (4 * n1);
        for(unsigned i = 0; i < n1; i++) {
          weight[i] = 1.;
        }

        for(unsigned i = n1; i < n0; i++) {
          weight[i] = exp(-d2[i] / sigma2 );
        }
        FindBestFit(xi, weight, N, a, d);

        //std::cout << a[0] << a[1] << d << std::endl;
      }
      else {//probably, it is not a cut cell
        double d2min = 1.0e10;
        unsigned imin = 0;
        for(unsigned i = 0; i < n0; i++) {
          double d2 = 0.;
          for(unsigned k = 0; k < dim; k++) {
            d2 += (xi[i][k] - _xig[ielType][k]) * (xi[i][k] - _xig[ielType][k]);
          }
          if(d2 < d2min) {
            imin = i;
            d2min = d2;
          }
        }
        a = Ni[imin];
        d = -a[0] * xi[imin][0] - a[1] * xi[imin][1];

        //N = Ni[imin];
        //FindBestFit(xi, boost::none, N, a, d);
      }



      CutFemWeight <double, cpp_bin_float_oct> *cutElemeFem = (ielType == 3) ? _quad : _tri;
      cutElemeFem->GetWeightWithMap(0, a, d, _weightCut2);
      for(unsigned k = 0; k < dim; k++) a[k] = - a[k];
      d = -d;
      cutElemeFem->GetWeightWithMap(0, a, d, _weightCut1);
      cutElemeFem->GetWeightWithMap(-1, a, d, _weightCutI);

      const elem_type *elemFem = _femFine->GetFiniteElement(ielType, 0);
      unsigned ng = elemFem->GetGaussPointNumber();

      unsigned size0 = _weight1.size();
      _weight1.resize(size0 + ng);
      _weight2.resize(size0 + ng);
      _weightI.resize(size0 + ng);
      _xvig.resize(size0 + ng, std::vector<double>(dim, 0));


      {
        //printing intersections
        double yi1 = 0.;
        double xi1 = -d / a[0];

        double xi2 = 0.;
        double yi2 = -d / a[1];

        double xi3 = (-d - a[1]) / (a[0] - a[1]);
        double yi3 = (-d - a[0] * xi3) / a[1];

        std::vector<double> xiP;
        std::vector<double> xP;
        if(xi1 > 0 && xi1 < 1) {
          xiP.resize(2);
          xiP[0] = xi1;
          xiP[1] = yi1;
          std::vector <double> phi;
          elemFem->GetPhi(phi, xiP);
          xP.assign(2, 0);
          for(unsigned i = 0; i < phi.size(); i++) {
            for(unsigned k = 0; k < dim; k++) {
              xP[k] += phi[i] * xv[k][i];
            }
          }
          std::cout << xP[0] << " " << xP[1] << "\n";
          if(yi2 > 0 && yi2 < 1) {
            xiP[0] = xi2;
            xiP[1] = yi2;
            elemFem->GetPhi(phi, xiP);
            xP.assign(2, 0);
            for(unsigned i = 0; i < phi.size(); i++) {
              for(unsigned k = 0; k < dim; k++) {
                xP[k] += phi[i] * xv[k][i];
              }
            }
            std::cout << xP[0] << " " << xP[1] << "\n\n";
          }
          else {
            xiP[0] = xi3;
            xiP[1] = yi3;
            elemFem->GetPhi(phi, xiP);
            xP.assign(2, 0);
            for(unsigned i = 0; i < phi.size(); i++) {
              for(unsigned k = 0; k < dim; k++) {
                xP[k] += phi[i] * xv[k][i];
              }
            }
            std::cout << xP[0] << " " << xP[1] << "\n\n";
          }
          for(unsigned i = 0; i < xv[0].size(); i++) {
            std::cout << xv[0][i] << " " << xv[1][i] << "\n";
          }
          std::cout << xv[0][0] << " " << xv[1][0] << "\n\n";
        }
        else if(yi2 > 0 && yi2 < 1) {
          xiP.resize(2);
          xiP[0] = xi2;
          xiP[1] = yi2;
          std::vector <double> phi;
          elemFem->GetPhi(phi, xiP);
          xP.assign(2, 0);
          for(unsigned i = 0; i < phi.size(); i++) {
            for(unsigned k = 0; k < dim; k++) {
              xP[k] += phi[i] * xv[k][i];
            }
          }
          std::cout << xP[0] << " " << xP[1] << "\n";
          xiP[0] = xi3;
          xiP[1] = yi3;
          elemFem->GetPhi(phi, xiP);
          xP.assign(2, 0);
          for(unsigned i = 0; i < phi.size(); i++) {
            for(unsigned k = 0; k < dim; k++) {
              xP[k] += phi[i] * xv[k][i];
            }
          }
          std::cout << xP[0] << " " << xP[1] << "\n\n";


          for(unsigned i = 0; i < xv[0].size(); i++) {
            std::cout << xv[0][i] << " " << xv[1][i] << "\n";
          }
          std::cout << xv[0][0] << " " << xv[1][0] << "\n\n";
        }
      }







      double Area = 0.;
      for(unsigned ig = 0; ig < ng; ig++) {
        double weight;
        elemFem->GetJacobianMatrix(xv, ig, weight, _Jac, _JacI);
        const double *phi = elemFem->GetPhi(ig);


        for(unsigned k = 0; k < dim; k++) {
          for(unsigned i = 0; i < xv[k].size(); i++) {
            _xvig[size0 + ig][k] += phi[i] * xv[k][i];
          }
        }


        Area += weight;
        double dsN = 0.;
        std::vector <double> Np(dim, 0.);
        for(unsigned k = 0; k < dim; k++) {
          for(unsigned j = 0; j < dim; j++) {
            Np[k] += _JacI[j][k] * a[j];
          }
          dsN += Np[k] * Np[k];
        }
        dsN = sqrt(dsN);
        _weight1[size0 + ig] = weight * _weightCut1[ig];
        _weight2[size0 + ig] = weight * _weightCut2[ig];
        _weightI[size0 + ig] = weight * _weightCutI[ig] * dsN;
      }
    }
  }
}
#endif



