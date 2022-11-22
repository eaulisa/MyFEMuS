#ifndef __femus_AdaptiveSplit_hpp__
#define __femus_AdaptiveSplit_hpp__

#include "CutFemWeight.hpp"
#include "CDWeights.hpp"
#include "Fem.hpp"

namespace femus {

  unsigned FindClosestPointToParentTriangle(const std::vector<std::vector<double>> &xi);
  unsigned FindClosestPointToParentSquare(const std::vector<std::vector<double>> &xi);
  void PrintIntersectionTriangle(const std::vector<std::vector<double>> &xv, const std::vector<double> &a, const double &d, const elem_type *fem);
  void PrintIntersectionSquare(const std::vector<std::vector<double>> &xv, const std::vector<double> &a, const double &d, const elem_type *fem);

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

      void Split(const std::vector<std::vector<double>> &xv, const unsigned ielType, const unsigned &level = 0, const unsigned &father = 0);
      const std::vector<double>& GetWeight1() {
        return _weight1;
      };
      const std::vector<double>& GetWeight2() {
        return _weight2;
      };
      const std::vector<double>& GetWeightI() {
        return _weightI;
      };
      const std::vector<std::vector<double>>& GetXg() {
        return _xg;
      };

      std::vector<std::vector<double>> &GetNiFather() {
        if(_Ni.size() == 0) _Ni.resize(1, std::vector<std::vector<std::vector<double>>>(1));
        return _Ni[0][0];
      }

      std::vector<std::vector<double>> &GetXiFather() {
        if(_xi.size() == 0) _xi.resize(1, std::vector<std::vector<std::vector<double>>>(1));
        return _xi[0][0];
      }

      std::vector<double> &GetDsFather() {
        if(_ds.size() == 0) _ds.resize(1, std::vector<std::vector<double>>(1));
        return _ds[0][0];
      }

      typedef unsigned(*FindClosestPointToParentElement)(const std::vector<std::vector<double>> &xi);
      typedef void (*PrintIntersection)(const std::vector<std::vector<double>> &xv, const std::vector<double> &a, const double &d, const elem_type *fem);

      FindClosestPointToParentElement findClosestPoint;
      PrintIntersection printIntersection;

      double GetC() {
        return _volume1 / (_volume1 + _volume2);
      }

      double GetVolume1() {
        return _volume1;
      }

      double GetVolume2() {
        return _volume2;
      }

      double GetSurfaceArea() {
        return _surface;
      }

      void RemapTriangles(const std::vector<double> &xi, const std::vector<double> &Ni, const double &ds,
                          const unsigned &i0, const unsigned &i1, const unsigned &i2, const unsigned &i3,
                          std::vector<std::vector<std::vector<double>>> &yi, std::vector<std::vector<std::vector<double>>> &Mi, std::vector<std::vector<double>> &dsi) {

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

        dsi[0][i0] = ds;
        dsi[1][i1] = ds;
        dsi[2][i2] = ds;
        dsi[3][i3] = ds;
      }


      void RemapSquares(const std::vector<double> &xi, const std::vector<double> &Ni, const double &ds,
                        const unsigned &i0, const unsigned &i1, const unsigned &i2, const unsigned &i3,
                        std::vector<std::vector<std::vector<double>>> &yi, std::vector<std::vector<std::vector<double>>> &Mi, std::vector<std::vector<double>> &dsi) {
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

        dsi[0][i0] = ds;
        dsi[1][i1] = ds;
        dsi[2][i2] = ds;
        dsi[3][i3] = ds;
      }

    private:

      const std::vector<std::vector<double>> _xig = {{0., 0., 0.}, {1. / 3., 1. / 3., 1. / 3.}, {1. / 3., 1. / 3., 0.}, {0., 0.}, {1. / 3., 1. / 3.}, {0.}};

      std::vector<std::vector<std::vector<std::vector<double>>>> _xi;
      std::vector<std::vector<std::vector<std::vector<double>>>> _Ni;
      std::vector<std::vector<std::vector<double>>> _ds;

      std::vector<std::vector<unsigned>> _n0;
      std::vector<std::vector<unsigned>> _n1;


      CutFemWeight <double, cpp_bin_float_oct> *_quad;
      CutFemWeight <double, cpp_bin_float_oct> *_tri;
      CutFemWeight <double, cpp_bin_float_oct> *_cutFem;

      Fem *_femFine;
      Fem *_femCoarse;
      const elem_type *_fem;

      double _volume1, _volume2, _surface;

      std::vector <double> _weight1;
      std::vector <double> _weight2;
      std::vector <double> _weightI;

      std::vector <double> _weightCut1;
      std::vector <double> _weightCut2;
      std::vector <double> _weightCutI;

      std::vector<std::vector <double>> _xg;

      std::vector<std::vector<double>> _Jac;
      std::vector<std::vector<double>> _JacI;
  };

  void AdaptiveSplit::Split(const std::vector<std::vector<double>> &xv, const unsigned ielType, const unsigned &level, const unsigned &father) {

    if(level == 0) {
      if(ielType == 3) {
        findClosestPoint =  &FindClosestPointToParentSquare;
        printIntersection = &PrintIntersectionSquare;
      }
      else {
        findClosestPoint =  &FindClosestPointToParentTriangle;
        printIntersection = &PrintIntersectionTriangle;
      }

      _volume1 = 0;
      _volume2 = 0;
      _surface = 0;

      _weight1.resize(0);
      _weight2.resize(0);
      _weightI.resize(0);
      _xg.resize(0);

      _n0.assign(1, std::vector<unsigned>(1, _xi[0][0].size()));
      _n1.assign(1, std::vector<unsigned>(1, _xi[0][0].size()));
    }


    const unsigned &n0 = _n0[level][father];
    const unsigned &n1 = _n1[level][father];
    const unsigned &dim = xv.size();

    if(level == 0 || (level < 2 && n1 >= 4)) {

      unsigned nChilds = (dim == 2) ? 4 : 8;
      std::vector<int> k(nChilds, n1 - 1);

      if(_xi.size() < level + 2) _xi.resize(level + 2);
      _xi[level + 1].assign(nChilds, std::vector<std::vector<double>>(n1, std::vector<double>(dim)));
      std::vector<std::vector<double>> &xi = _xi[level][father];

      if(_Ni.size() < level + 2) _Ni.resize(level + 2);
      _Ni[level + 1].assign(nChilds, std::vector<std::vector<double>>(n1, std::vector<double>(dim)));
      std::vector<std::vector<double>> &Ni = _Ni[level][father];

      if(_ds.size() < level + 2) _ds.resize(level + 2);
      _ds[level + 1].assign(nChilds, std::vector<double>(n1));
      std::vector<double> &ds = _ds[level][father];


      if(_n0.size() < level + 2) _n0.resize(level + 2);
      _n0[level + 1].assign(nChilds, n1);

      if(_n1.size() < level + 2) _n1.resize(level + 2);
      _n1[level + 1].assign(nChilds, 0);

      if(ielType == 3) {
        for(unsigned i = 0; i < n1; i++) {
          if(xi[i][0] < 0) {
            if(xi[i][1] < 0) {
              RemapSquares(xi[i], Ni[i], ds[i], _n1[level + 1][0]++, k[1]--, k[2]--, k[3]--, _xi[level + 1], _Ni[level + 1], _ds[level + 1]);
            }
            else {
              RemapSquares(xi[i], Ni[i], ds[i], k[0]--, k[1]--, k[2]--, _n1[level + 1][3]++, _xi[level + 1], _Ni[level + 1], _ds[level + 1]);
            }
          }
          else {
            if(xi[i][1] < 0) {
              RemapSquares(xi[i], Ni[i], ds[i], k[0]--, _n1[level + 1][1]++, k[2]--, k[3]--, _xi[level + 1], _Ni[level + 1], _ds[level + 1]);
            }
            else {
              RemapSquares(xi[i], Ni[i], ds[i], k[0]--, k[1]--, _n1[level + 1][2]++, k[3]--, _xi[level + 1], _Ni[level + 1], _ds[level + 1]);
            }
          }
        }
      }

      else if(ielType == 4) {
        for(unsigned i = 0; i < n1; i++) {
          if(xi[i][0] > 0.5) {
            RemapTriangles(xi[i], Ni[i], ds[i], k[0]--, _n1[level + 1][1]++, k[2]--, k[3]--, _xi[level + 1], _Ni[level + 1], _ds[level + 1]);
          }
          else if(xi[i][1] > 0.5) {
            RemapTriangles(xi[i], Ni[i], ds[i], k[0]--, k[1]--, _n1[level + 1][2]++, k[3]--, _xi[level + 1], _Ni[level + 1], _ds[level + 1]);
          }
          else if(xi[i][0] + xi[i][1] < 0.5) {
            RemapTriangles(xi[i], Ni[i], ds[i], _n1[level + 1][0]++, k[1]--, k[2]--, k[3]--, _xi[level + 1], _Ni[level + 1], _ds[level + 1]);
          }
          else {
            RemapTriangles(xi[i], Ni[i], ds[i], k[0]--, k[1]--, k[2]--, _n1[level + 1][3]++, _xi[level + 1], _Ni[level + 1], _ds[level + 1]);
          }
        }
      }

      for(unsigned l = 0; l < nChilds; l++) {
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
      std::vector <double> &ds = _ds[level][father];

      std::vector<double> xim(dim, 0);
      std::vector<double> N(dim, 0);

      if(n0 == 1 || n1 == 1) { // for cut and no cut
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
          weight[i] = ds[i];
        }

        for(unsigned i = n1; i < n0; i++) {
          weight[i] = ds[i] * exp(-d2[i] / sigma2);
        }
        FindBestFit(xi, weight, N, a, d);
      }
      else {//probably, it is not a cut cell
        unsigned imin = (*findClosestPoint)(xi); //FindClosestPointToParentElement;
        a = Ni[imin];
        d = -a[0] * xi[imin][0] - a[1] * xi[imin][1];
      }


      _cutFem = (ielType == 3) ? _quad : _tri;
      (*_cutFem)(0, a, d, _weightCut2);

      unsigned ii = 1;
      while(ii < _weightCut2.size() && fabs(_weightCut2[ii] - _weightCut2[0]) < 1.0e-10) {
        ii++;
      }
      if(ii != _weightCut2.size()) { //it is a cut cell
        for(unsigned k = 0; k < dim; k++) a[k] = - a[k];
        d = -d;
        (*_cutFem)(0, a, d, _weightCut1);
        (*_cutFem)(-1, a, d, _weightCutI);

        _fem = _femFine->GetFiniteElement(ielType, 0);
        (*printIntersection)(xv, a, d, _fem);
        unsigned ng = _fem->GetGaussPointNumber();
        unsigned size0 = _weight1.size();
        _weight1.resize(size0 + ng);
        _weight2.resize(size0 + ng);
        _weightI.resize(size0 + ng);
        _xg.resize(size0 + ng, std::vector<double>(dim, 0));

        for(unsigned ig = 0; ig < ng; ig++) {
          double weight;
          _fem->GetJacobianMatrix(xv, ig, weight, _Jac, _JacI);
          const double *phi = _fem->GetPhi(ig);

          for(unsigned k = 0; k < dim; k++) {
            for(unsigned i = 0; i < xv[k].size(); i++) {
              _xg[size0 + ig][k] += phi[i] * xv[k][i];
            }
          }

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

          _volume1 += _weight1[size0 + ig];
          _volume2 += _weight2[size0 + ig];
          _surface += _weightI[size0 + ig];
        }
      }
      else { //it is full or empty

        _fem = _femCoarse->GetFiniteElement(ielType, 0);
        unsigned ng = _fem->GetGaussPointNumber();
        unsigned size0 = _weight1.size();

        double C = (_weightCut2[0] > 0.5) ? 0 : 1.;
        _weight1.resize(size0 + ng, C);
        _weight2.resize(size0 + ng, 1. - C);
        _weightI.resize(size0 + ng, 0.);
        _xg.resize(size0 + ng, std::vector<double>(dim, 0));

        std::vector <double> &weightRef = (C > 0.5) ? _weight1 : _weight2;
        double &volumeRef = (C > 0.5) ? _volume1 : _volume2;
        for(unsigned ig = 0; ig < ng; ig++) {
          double weight;
          _fem->GetJacobianMatrix(xv, ig, weight, _Jac, _JacI);
          const double *phi = _fem->GetPhi(ig);

          for(unsigned k = 0; k < dim; k++) {
            for(unsigned i = 0; i < xv[k].size(); i++) {
              _xg[size0 + ig][k] += phi[i] * xv[k][i];
            }
          }
          weightRef[size0 + ig] = weight;
          volumeRef += weight;
        }
      }



    }
  }

  unsigned FindClosestPointToParentTriangle(const std::vector<std::vector<double>> &xi) {
    double d2min = 1.e10;
    unsigned imin = 0;
    double d2;
    for(unsigned i = 0; i < xi.size(); i++) {
      const double &x = xi[i][0];
      const double &y = xi[i][1];
      if(x < 0) {
        if(y < 0) {
          d2 = x * x + y * y;
        }
        else if(y < 1) {
          d2 = x * x;
        }
        else {
          d2 = x * x + (1 - y) * (1 - y);
        }
      }
      else if(y < 0) {
        if(x < 1) {
          d2 = y * y;
        }
        else {
          d2 = (1 - x) * (1 - x) + y * y;
        }
      }
      else if(y < x - 1) {
        d2 = (1 - x) * (1 - x) + y * y;
      }
      else if(y < x + 1) {
        d2 = 0.5 * (x + y - 1) * (x + y - 1);
      }
      else {
        d2 = x * x + (1 - y) * (1 - y);
      }
      if(d2 < d2min) {
        imin = i;
        d2min = d2;
      }
    }
    return imin;
  }



  unsigned FindClosestPointToParentSquare(const std::vector<std::vector<double>> &xi) {
    double d2min = 1.e10;
    unsigned imin = 0;
    double d2;
    for(unsigned i = 0; i < xi.size(); i++) {
      const double &x = xi[i][0];
      const double &y = xi[i][1];
      if(x < -1.) {
        if(y < -1.) {
          d2 = (x + 1.) * (x + 1.) + (y + 1.) * (y + 1.);
        }
        else if(y < 1.) {
          d2 = (x + 1.) * (x + 1.);
        }
        else {
          d2 = (x + 1.) * (x + 1.) + (y - 1.) * (y - 1.);
        }
      }
      else if(x > 1.) {
        if(y < -1.) {
          d2 = (x - 1.) * (x - 1.) + (y + 1.) * (y + 1.);
        }
        else if(y < 1.) {
          d2 = (x - 1.) * (x - 1.);
        }
        else {
          d2 = (x - 1.) * (x - 1.) + (y - 1.) * (y - 1.);
        }
      }
      else if(y < -1.) {
        d2 = (y + 1.) * (y + 1.);
      }
      else {
        d2 = (y - 1.) * (y - 1.);
      }

      if(d2 < d2min) {
        imin = i;
        d2min = d2;
      }

    }
    return imin;
  }



#include <iostream>
#include <fstream>
  void PrintIntersectionTriangle(const std::vector<std::vector<double>> &xv, const std::vector<double> &a, const double &d, const elem_type *fem) {

    std::ofstream fout;
    fout.open("intersection.txt", std::ios::app);

    //printing intersections
    unsigned dim = xv.size();

    double yi1 = 0.;
    double xi1 = -d / a[0];

    double xi2 = 0.;
    double yi2 = -d / a[1];

    double xi3 = (-d - a[1]) / (a[0] - a[1]);
    double yi3 = (-d - a[0]) / (a[1] - a[0]);

    std::vector<double> xiP;
    std::vector<double> xP;
    if(xi1 > 0 && xi1 < 1) {
      xiP.resize(2);
      xiP[0] = xi1;
      xiP[1] = yi1;
      std::vector <double> phi;
      fem->GetPhi(phi, xiP);
      xP.assign(2, 0);
      for(unsigned i = 0; i < phi.size(); i++) {
        for(unsigned k = 0; k < dim; k++) {
          xP[k] += phi[i] * xv[k][i];
        }
      }
      fout << xP[0] << " " << xP[1] << "\n";
      if(yi2 > 0 && yi2 < 1) {
        xiP[0] = xi2;
        xiP[1] = yi2;
        fem->GetPhi(phi, xiP);
        xP.assign(2, 0);
        for(unsigned i = 0; i < phi.size(); i++) {
          for(unsigned k = 0; k < dim; k++) {
            xP[k] += phi[i] * xv[k][i];
          }
        }
        fout << xP[0] << " " << xP[1] << "\n\n";
      }
      else {
        xiP[0] = xi3;
        xiP[1] = yi3;
        fem->GetPhi(phi, xiP);
        xP.assign(2, 0);
        for(unsigned i = 0; i < phi.size(); i++) {
          for(unsigned k = 0; k < dim; k++) {
            xP[k] += phi[i] * xv[k][i];
          }
        }
        fout << xP[0] << " " << xP[1] << "\n\n";
      }
      for(unsigned i = 0; i < xv[0].size(); i++) {
        fout << xv[0][i] << " " << xv[1][i] << "\n";
      }
      fout << xv[0][0] << " " << xv[1][0] << "\n\n";
    }
    else if(yi2 > 0 && yi2 < 1) {
      xiP.resize(2);
      xiP[0] = xi2;
      xiP[1] = yi2;
      std::vector <double> phi;
      fem->GetPhi(phi, xiP);
      xP.assign(2, 0);
      for(unsigned i = 0; i < phi.size(); i++) {
        for(unsigned k = 0; k < dim; k++) {
          xP[k] += phi[i] * xv[k][i];
        }
      }
      fout << xP[0] << " " << xP[1] << "\n";
      xiP[0] = xi3;
      xiP[1] = yi3;
      fem->GetPhi(phi, xiP);
      xP.assign(2, 0);
      for(unsigned i = 0; i < phi.size(); i++) {
        for(unsigned k = 0; k < dim; k++) {
          xP[k] += phi[i] * xv[k][i];
        }
      }
      fout << xP[0] << " " << xP[1] << "\n\n";

      for(unsigned i = 0; i < xv[0].size(); i++) {
        fout << xv[0][i] << " " << xv[1][i] << "\n";
      }
      fout << xv[0][0] << " " << xv[1][0] << "\n\n";
    }
    fout.close();
  }


  void PrintIntersectionSquare(const std::vector<std::vector<double>> &xv, const std::vector<double> &a, const double &d, const elem_type *fem) {

    std::ofstream fout;
    fout.open("intersection.txt", std::ios::app);

    //printing intersections
    unsigned dim = xv.size();

    double yi1 = -1.;
    double xi1 = (-d + a[1]) / a[0];

    double xi2 = 1.;
    double yi2 = (-d - a[0]) / a[1];

    double yi3 = 1.;
    double xi3 = (-d - a[1]) / a[0];

    double xi4 = -1.;
    double yi4 = (-d + a[0]) / a[1];

    std::vector<double> xiP(2);
    std::vector<double> xP;
    if(xi1 > -1 && xi1 < 1) { //bottom
      xiP[0] = xi1;
      xiP[1] = yi1;
      std::vector <double> phi;
      fem->GetPhi(phi, xiP);
      xP.assign(2, 0);
      for(unsigned i = 0; i < phi.size(); i++) {
        for(unsigned k = 0; k < dim; k++) {
          xP[k] += phi[i] * xv[k][i];
        }
      }
      fout << xP[0] << " " << xP[1] << "\n";
      if(yi2 > -1 && yi2 < 1) { //bottom-right
        xiP[0] = xi2;
        xiP[1] = yi2;
        fem->GetPhi(phi, xiP);
        xP.assign(2, 0);
        for(unsigned i = 0; i < phi.size(); i++) {
          for(unsigned k = 0; k < dim; k++) {
            xP[k] += phi[i] * xv[k][i];
          }
        }
      }
      else if(xi3 > -1 && xi3 < 1) { //bottom-top
        xiP[0] = xi3;
        xiP[1] = yi3;
        fem->GetPhi(phi, xiP);
        xP.assign(2, 0);
        for(unsigned i = 0; i < phi.size(); i++) {
          for(unsigned k = 0; k < dim; k++) {
            xP[k] += phi[i] * xv[k][i];
          }
        }
      }
      else { //bottom-left
        xiP[0] = xi4;
        xiP[1] = yi4;
        fem->GetPhi(phi, xiP);
        xP.assign(2, 0);
        for(unsigned i = 0; i < phi.size(); i++) {
          for(unsigned k = 0; k < dim; k++) {
            xP[k] += phi[i] * xv[k][i];
          }
        }
      }
      fout << xP[0] << " " << xP[1] << "\n\n";
      for(unsigned i = 0; i < xv[0].size(); i++) {
        fout << xv[0][i] << " " << xv[1][i] << "\n";
      }
      fout << xv[0][0] << " " << xv[1][0] << "\n\n";
    }
    else if(yi2 > -1 && yi2 < 1) { // right
      xiP[0] = xi2;
      xiP[1] = yi2;
      std::vector <double> phi;
      fem->GetPhi(phi, xiP);
      xP.assign(2, 0);
      for(unsigned i = 0; i < phi.size(); i++) {
        for(unsigned k = 0; k < dim; k++) {
          xP[k] += phi[i] * xv[k][i];
        }
      }
      fout << xP[0] << " " << xP[1] << "\n"; // right - top
      if(xi3 > -1 && xi3 < 1) {
        xiP[0] = xi3;
        xiP[1] = yi3;
        fem->GetPhi(phi, xiP);
        xP.assign(2, 0);
        for(unsigned i = 0; i < phi.size(); i++) {
          for(unsigned k = 0; k < dim; k++) {
            xP[k] += phi[i] * xv[k][i];
          }
        }
      }
      else { // right - left
        xiP[0] = xi4;
        xiP[1] = yi4;
        fem->GetPhi(phi, xiP);
        xP.assign(2, 0);
        for(unsigned i = 0; i < phi.size(); i++) {
          for(unsigned k = 0; k < dim; k++) {
            xP[k] += phi[i] * xv[k][i];
          }
        }
      }
      fout << xP[0] << " " << xP[1] << "\n\n";
      for(unsigned i = 0; i < xv[0].size(); i++) {
        fout << xv[0][i] << " " << xv[1][i] << "\n";
      }
      fout << xv[0][0] << " " << xv[1][0] << "\n\n";
    }
    else if(xi3 > -1 && xi3 < 1) { // top
      xiP[0] = xi3;
      xiP[1] = yi3;
      std::vector <double> phi;
      fem->GetPhi(phi, xiP);
      xP.assign(2, 0);
      for(unsigned i = 0; i < phi.size(); i++) {
        for(unsigned k = 0; k < dim; k++) {
          xP[k] += phi[i] * xv[k][i];
        }
      }
      fout << xP[0] << " " << xP[1] << "\n";
      xiP[0] = xi4;  // top - left
      xiP[1] = yi4;
      fem->GetPhi(phi, xiP);
      xP.assign(2, 0);
      for(unsigned i = 0; i < phi.size(); i++) {
        for(unsigned k = 0; k < dim; k++) {
          xP[k] += phi[i] * xv[k][i];
        }
      }
      fout << xP[0] << " " << xP[1] << "\n\n";

      for(unsigned i = 0; i < xv[0].size(); i++) {
        fout << xv[0][i] << " " << xv[1][i] << "\n";
      }
      fout << xv[0][0] << " " << xv[1][0] << "\n\n";
    }
    fout.close();
  }

}
#endif



