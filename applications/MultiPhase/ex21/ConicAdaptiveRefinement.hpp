
#ifndef __Conic_Adaptive_Refinement_hpp__
#define __Conic_Adaptive_Refinement_hpp__

#include <iostream>
#include <cmath>
#include <vector>

#include "CutFemWeight.hpp"
#include "CDWeights.hpp"


typedef double TypeIO;
typedef cpp_bin_float_oct TypeA;
typedef cpp_bin_float_oct oct;

#include "Fem.hpp"

using namespace std;
using namespace femus;


class Data {
  public:
    Data(const unsigned &VType, const unsigned &PType,
         const std::vector<std::vector<double>> &V, const std::vector<double> &P1, const std::vector<double> &P2,
         std::vector<double> &res, std::vector<double> &jac,
         const std::vector<std::vector<double>> &xv, const unsigned &elType,
         const double &rho1, const double &rho2, const double &mu1, const double &mu2, const double &sigma, const double&dt) :
      _VType(VType), _PType(PType), _V(V), _P1(P1), _P2(P2), _res(res), _jac(jac), _elType(elType), _xv(xv), _rho1(rho1), _rho2(rho2), _mu1(mu1), _mu2(mu2), _sigma(sigma), _dt(dt) {}

    // DATA TO ASSEMBLE TWO PHASE NAVIER-STOKES
    const unsigned &_VType, &_PType;
    const std::vector<std::vector<double>> &_V;
    const std::vector<double> &_P1;
    const std::vector<double> &_P2;

    std::vector<double> &_res;
    std::vector<double> &_jac;

    const unsigned &_elType;
    const std::vector<std::vector <double> > &_xv;

    const double &_rho1, &_rho2, &_mu1, &_mu2, &_sigma, &_dt;
};

void AssembleNavierStokes(Data *data, const std::vector <double> &phiV, const std::vector <double> &phiV_x, std::vector <double> &phiP,
                          const double &weight, const double &weight1, const double &weight2, const double &weightI,
                          const std::vector <double> &N, const double &kappa, const double &dsN, const double &eps);


class ConicAdaptiveRefinement {
  public:
    ConicAdaptiveRefinement() {



      // coordinates in the reference element of the nodes of the target elements

      _weightCF[3] = new CutFemWeight <double, TypeA > (QUAD, 5, "legendre");
      _weightCF[4] = new CutFemWeight <double, TypeA > (TRI, 5, "legendre");

      //_quadCF = new CutFemWeight <double, TypeA > (QUAD, 5, "legendre");

      _fem1 = new Fem(3, 2);
      _fem2 = new Fem(_weightCF[3]->GetGaussQuadratureOrder(), _weightCF[3]->GetDimension());


      double dx = .025;
      double dtetha = 1.;

      _weightCD[3] = new CDWeightQUAD <TypeA> (5, dx, dtetha);
      _weightCD[4] = new CDWeightTRI <TypeA> (5, dx, dtetha);

      // dx = 0.00625;
      // dtetha = .25;
      // _quadCDI = new CDWeightQUAD <TypeA> (5, dx, dtetha, true);

    }

    ~ConicAdaptiveRefinement() {
      delete _weightCF[3];
      delete _weightCF[4];

      delete _weightCD[3];
      delete _weightCD[4];

      delete _fem1;
      delete _fem2;
    };

    const std::vector<std::vector<double>> & GetXiInParentElement() {
      return _xr[_data->_elType];
    };

    void GetXInParentElement(const std::vector<std::vector<double>> &xv, std::vector<std::vector<double>> &y) {
      if(_data->_elType == 3) {
        y = {{xv[0][0], xv[1][0]}, {xv[0][1], xv[1][1]}, {xv[0][2], xv[1][2]}, {xv[0][3], xv[1][3]} };
      }
      else if(_data->_elType == 4) {
        y = {{xv[0][0], xv[1][0]}, {xv[0][1], xv[1][1]}, {xv[0][2], xv[1][2]}};
      }
    };

    // void CalculateConicsInParentElement(const std::vector<std::vector<double>> &x, const std::vector<double> &A, std::vector<double> &A1);
    void GetConicsInTargetElement(const std::vector<std::vector<double>> &x, const std::vector<double> &A, std::vector<double> &A1);
    void ComputeJacobian(const std::vector<std::vector<double>> &x, std::vector<std::vector<double>> &J);
    void ComputeInverseJacobian(std::vector<std::vector<double>> &J, std::vector<std::vector<double>> &IJ);
    double EvaluateConic(const std::vector<double>&x, const std::vector<double>&a) {
      return a[0] * x[0] * x[0] + a[1] * x[0] * x[1] + a[2] * x[1] * x[1] + a[3] * x[0] + a[4] * x[1] + a[5];
    }

    void EvaluateConicNormal(const std::vector<double>&x, const std::vector<double>&a, std::vector<double>&n) {
      n.resize(2);
      n[0] = 2. * a[0] * x[0] + a[1] * x[1] + a[3];
      n[1] = 2. * a[2] * x[1] + a[1] * x[0] + a[4];
      double det = sqrt(n[0] * n[0] + n[1] * n[1]);
      n[0] /= det;
      n[1] /= det;
    }


    int TestIfIntesection(const std::vector<std::vector<double>>&x, const std::vector<double>&a);
    int TestIfIntesectionWithReferenceQuad(const std::vector<double>&A);

    void BestFitLinearInterpolation(const std::vector<double> &A, std::vector<double> &B);


    //double GetVolumeFraction(const std::vector<double>&a);

    std::tuple<double, double, double> AdaptiveRefinement(const unsigned & levelMax,
                                                          std::vector<std::vector<double>>&x,
                                                          std::vector<std::vector<double>>&xi,
                                                          const std::vector<double>&Ar,
                                                          const unsigned &level = 1);

    void PrintElement(const unsigned & level, const unsigned & jp, const std::vector<std::vector<double>>&x) {
      const unsigned &n = x.size();
      const unsigned &dim = x[0].size();
      for(unsigned i = 0; i < n ; i++) {
        std::cout << level << " " << jp << " ";
        for(unsigned j = 0; j < dim; j++) {
          std::cout << x[i][j] << " ";
        }
        std::cout << 0 << " ";
        std::cout << std::endl;
      }
      std::cout << level << " " << jp << " ";
      for(unsigned j = 0; j < dim; j++) {
        std::cout << x[n - 1][j] << " ";
      }
      std::cout << 0 << " ";
      std::cout << std::endl << std::endl;
    }


    bool CheckIfRootsAreInBetweenM1andP1(const double & A, const double & B, const double & C) {
      double det = sqrt(B * B - 4. * A * C);
      if(det >= 0.) {
        double t = (-B - det) / (2. * A);
        if(t >= -1. && t <= 1.) return true;
        t = (-B + det) / (2. * A);
        if(t >= -1. && t <= 1.) return true;
      }
      return false;
    }

    bool CheckIfRootsAreInBetween0andP1(const double & A, const double & B, const double & C) {
      double det = sqrt(B * B - 4. * A * C);
      if(det >= 0.) {
        double t = (-B - det) / (2. * A);
        if(t >= 0. && t <= 1.) return true;
        t = (-B + det) / (2. * A);
        if(t >= 0. && t <= 1.) return true;
      }
      return false;
    }


    void SetDataPointer(Data * data) {
      _data = data;
    }

  private:


    std::vector<std::vector<std::vector<std::vector<double>>>> _y;
    std::vector<std::vector<std::vector<std::vector<double>>>> _yi;
    Fem* _fem1, *_fem2;

    double _weight;
    std::vector <double> _phi;
    std::vector <double> _phix;

    std::vector <double> _phiV;
    std::vector <double> _phiVx;


    //std::vector<std::vector<double>> _xl0;

    std::vector<double> _xg;


    CutFemWeight <double, TypeA> *_weightCF[6];
    CDWeight <TypeA> *_weightCD[6];

    Data *_data;

    const std::vector<unsigned>_nve0 = {0, 0, 0, 4, 3};
    const std::vector<unsigned>_nve2 = {0, 0, 0, 9, 7};
    const std::vector<std::vector<std::vector<unsigned>>> _idx = {
      {}, {}, {},
      { // quad
        {0, 4, 8, 7},
        {4, 1, 5, 8},
        {8, 5, 2, 6},
        {7, 8, 6, 3}
      },
      { // triangle
        {0, 3, 5},
        {3, 1, 4},
        {5, 4, 2},
        {4, 5, 3}
      }
    };

    const std::vector<std::vector<std::vector<double>>> _xr = {
      {}, {}, {},
      {{1, -1}, {1, 1}, {-1, 1}, {-1, -1}},//quad
      {{0, 0}, {1, 0}, {0, 1}}//triange
    };
    const std::vector<std::vector<std::vector<std::vector<double>>>> _yr = {{}, {}, {},
      { {{-1, -1}, {0, -1}, {0, 0}, {-1, 0}},//quad
        {{0, -1}, {1, -1}, {1, 0}, {0, 0}},
        {{0, 0}, {1, 0}, {1, 1}, {0, 1}},
        {{-1, 0}, {0, 0}, {0, 1}, {-1, 1}}
      },
      { {{0., 0.}, {0.5, 0.}, {0., 0.5}},//triangle
        {{0.5, 0.}, {1., 0.}, {0.5, 0.5}},
        {{0., 0.5}, {0.5, 0.5}, {0., 1.}},
        {{0.5, 0.5}, {0, 0.5}, {0.5, 0.}}
      }
    };
};

//coefficients of conic in the reference system
void ConicAdaptiveRefinement::GetConicsInTargetElement(const std::vector<std::vector<double>> &x, const std::vector<double> &A, std::vector<double> &A1) {
  const double &a = A[0];
  const double &b = A[1];
  const double &c = A[2];
  const double &d = A[3];
  const double &e = A[4];
  const double &f = A[5];

  if(_data->_elType == 3) {
    const double &x1 = x[0][0];
    const double &x4 = x[3][0];

    const double &y1 = x[0][1];
    const double &y2 = x[1][1];

    double Lx = x[1][0] - x1;
    double Ly = x[3][1] - y1;

    //compact
    double Lx14 = Lx + x1 + x4;
    double Ly12 = Ly + y1 + y2;
    double y12 = y1 - y2;
    double x14 = x1 - x4;
    double Lxy = Lx * Ly;

    A1.resize(A.size());

    A1[0] = 0.25 * (a * Lx * Lx + (-b * Lx + c * y12) * y12);
    A1[1] = 0.25 * (b * Lxy - 2 * a * Lx * x14 + b * x14 * y12 - 2 * c * Ly * y12);
    A1[2] = 0.25 * (c * Ly * Ly + (-b * Ly + a * x14) * x14);
    A1[3] = 0.25 * (2 * d * Lx + b * Lxy + 2 * a * Lx * Lx14 - 2 * e * y1 - 2 * c * Ly * y1 - b * x1 * y1 - b * x4 * y1 - 2 * c * y1 * y1 + (2 * e + 2 * c * Ly + b * (2 * Lx + x1 + x4)) * y2 + 2 * c * y2 * y2);
    A1[4] = 0.25 * (2 * e * Ly - 2 * x14 * (d + a * Lx14) + 2 * c * Ly * Ly12 + b * (Lxy - x1 * (y1 + y2) + x4 * (2 * Ly + y1 + y2)));
    A1[5] = 0.25 * (4 * f + 2 * d * Lx14 + a * Lx14 * Lx14 + Ly12 * (2 * e + b * Lx14 + c * Ly12));
  }

  else if(_data->_elType == 4) {
    const double &x1 = x[0][0];
    const double &y1 = x[0][1];

    double L2x = x[1][0] - x1;
    double L2y = x[1][1] - y1;

    double L3x = x[2][0] - x1;
    double L3y = x[2][1] - y1;


    A1.resize(A.size());

    A1[0] = a * L2x * L2x + L2y * (b * L2x + c * L2y);
    A1[1] = 2. * a * L2x * L3x + b * L2y * L3x + b * L2x * L3y + 2. * c * L2y * L3y;
    A1[2] = a * L3x * L3x + L3y * (b * L3x + c * L3y);
    A1[3] = d * L2x + e * L2y + 2. * a * L2x * x1 + b * L2y * x1 + b * L2x * y1 + 2. * c * L2y * y1;
    A1[4] = d * L3x + e * L3y + 2. * a * L3x * x1 + b * L3y * x1 + b * L3x * y1 + 2. * c * L3y * y1;
    A1[5] = f + d * x1 + a * x1 * x1 + y1 * (e + b * x1 + c * y1);
  }

}


int ConicAdaptiveRefinement::TestIfIntesection(const std::vector<std::vector<double>>&x, const std::vector<double>&a) {

  double value0 = EvaluateConic(x[0], a);
  for(unsigned i = 1; i < x.size(); i++) {
    if(EvaluateConic(x[i], a) * value0 <= 0) return 0;
  }
  return (value0 > 0) ? 1 : -1;
}

int ConicAdaptiveRefinement::TestIfIntesectionWithReferenceQuad(const std::vector<double>&Ar) {

  const double& a = Ar[0];
  const double& b = Ar[1];
  const double& c = Ar[2];
  const double& d = Ar[3];
  const double& e = Ar[4];
  const double& f = Ar[5];


  double det;
  double A, B, C;

  if(_data->_elType == 3) {

    //main diagonal
    A = a + b + c;
    B = d + e;
    C = f;
    if(CheckIfRootsAreInBetweenM1andP1(A, B, C)) return 0;

    //other diagonal
    A = a - b + c;
    B = -d + e;
    C = f;
    if(CheckIfRootsAreInBetweenM1andP1(A, B, C)) return 0;

    //bottom edge
    A = a;
    B = -b + d;
    C = c - e + f;
    if(CheckIfRootsAreInBetweenM1andP1(A, B, C)) return 0;

    //top edge
    A = a;
    B = b + d;
    C = c + e + f;
    if(CheckIfRootsAreInBetweenM1andP1(A, B, C)) return 0;

    //left edge
    A = c;
    B = -b + e;
    C = a - d + f;
    if(CheckIfRootsAreInBetweenM1andP1(A, B, C)) return 0;

    //right edge
    A = c;
    B =  b + e;
    C = a + d + f;
    if(CheckIfRootsAreInBetweenM1andP1(A, B, C)) return 0;

    return (f > 0) ? 1 : -1;
  }
  else if(_data->_elType == 4) {
    //bottom edge
    A = a;
    B = d;
    C = f;
    if(CheckIfRootsAreInBetween0andP1(A, B, C)) return 0;

    //left edge
    A = c;
    B = e;
    C = f;
    if(CheckIfRootsAreInBetween0andP1(A, B, C)) return 0;

    //other edge (x,y=-x+1)
    A = a - b + c;
    B = b - 2. * c + d - e;
    C = c + e + f;
    if(CheckIfRootsAreInBetween0andP1(A, B, C)) return 0;

    //Line 1 passes through barycenter (x,y=-0.5x+0.5)
    A = a - 0.5 * b + 0.25 * c;
    B = 0.5 * b - 0.5 * c + d - 0.5 * e;
    C = 0.25 * c + 0.5 * e + f;
    if(CheckIfRootsAreInBetween0andP1(A, B, C)) return 0;

    //Line 2 passes through barycenter (x=-0.5y+0.5,y)
    A = 0.25 * a - 0.5 * b + c;
    B = -0.5 * a + 0.5 * b - 0.5 * d + e;
    C = 0.25 * a + 0.5 * d + f;
    if(CheckIfRootsAreInBetween0andP1(A, B, C)) return 0;

    return ((a + b + c) / 9. + (d + e) / 3. + f > 0) ? 1 : -1;
  }
}


std::tuple<double, double, double> ConicAdaptiveRefinement::AdaptiveRefinement(
  const unsigned &levelMax,
  std::vector<std::vector<double>>&x, // physical coordinates at the nodes of the l-mesh
  std::vector<std::vector<double>>&xil0,// l0 parent coordinates at the nodes of the l-mesh
  const std::vector<double>&Ar,
  const unsigned &level) { // myconic

  double area1 = 0.;
  double area2 = 0.;
  double arcLenght = 0.;

  if(level == 1) {
    _y.resize(levelMax);
    _yi.resize(levelMax);
  }

  if(level < levelMax) {
    int test = TestIfIntesectionWithReferenceQuad(Ar);
    if(test == 0) { // it means there is an intersection

      x.resize(_nve2[_data->_elType]);
      xil0.resize(_nve2[_data->_elType]);

      if(_data->_elType == 3) {
        x[4] = {0.5 * (x[0][0] + x[1][0]), 0.5 * (x[0][1] + x[1][1])};
        x[5] = {0.5 * (x[1][0] + x[2][0]), 0.5 * (x[1][1] + x[2][1])};
        x[6] = {0.5 * (x[2][0] + x[3][0]), 0.5 * (x[2][1] + x[3][1])};
        x[7] = {0.5 * (x[3][0] + x[0][0]), 0.5 * (x[3][1] + x[0][1])};
        x[8] = {0.5 * (x[4][0] + x[6][0]), 0.5 * (x[5][1] + x[7][1])};

        xil0[4] = {0.5 * (xil0[0][0] + xil0[1][0]), 0.5 * (xil0[0][1] + xil0[1][1])};
        xil0[5] = {0.5 * (xil0[1][0] + xil0[2][0]), 0.5 * (xil0[1][1] + xil0[2][1])};
        xil0[6] = {0.5 * (xil0[2][0] + xil0[3][0]), 0.5 * (xil0[2][1] + xil0[3][1])};
        xil0[7] = {0.5 * (xil0[3][0] + xil0[0][0]), 0.5 * (xil0[3][1] + xil0[0][1])};
        xil0[8] = {0.5 * (xil0[4][0] + xil0[6][0]), 0.5 * (xil0[5][1] + xil0[7][1])};
      }
      else if(_data->_elType == 4) {
        x[3] = {0.5 * (x[0][0] + x[1][0]), 0.5 * (x[0][1] + x[1][1])};
        x[4] = {0.5 * (x[1][0] + x[2][0]), 0.5 * (x[1][1] + x[2][1])};
        x[5] = {0.5 * (x[2][0] + x[0][0]), 0.5 * (x[2][1] + x[0][1])};
        x[6] = {(x[0][0] + 2 * x[4][0]) / 3., (x[0][1] + 2 * x[4][1]) / 3.};

        xil0[3] = {0.5 * (xil0[0][0] + xil0[1][0]), 0.5 * (xil0[0][1] + xil0[1][1])};
        xil0[4] = {0.5 * (xil0[1][0] + xil0[2][0]), 0.5 * (xil0[1][1] + xil0[2][1])};
        xil0[5] = {0.5 * (xil0[2][0] + xil0[0][0]), 0.5 * (xil0[2][1] + xil0[0][1])};
        xil0[6] = {(xil0[0][0] + 2 * xil0[4][0]) / 3., (xil0[0][1] + 2 * xil0[4][1]) / 3.};
      }

      _y[level - 1].reserve(_nve2[_data->_elType]);
      _y[level - 1].resize(4, std::vector< std::vector<double> > (_nve0[_data->_elType], std::vector<double>(2)));

      _yi[level - 1].reserve(_nve2[_data->_elType]);
      _yi[level - 1].resize(4, std::vector< std::vector<double> > (_nve0[_data->_elType], std::vector<double>(2)));
      for(unsigned i = 0; i < 4; i++) { // number of childer stays the same (in 2D is 4)
        for(unsigned j = 0; j < _nve0[_data->_elType]; j++) {
          _y[level - 1][i][j] = x[_idx[_data->_elType][i][j]];
          _yi[level - 1][i][j] = xil0[_idx[_data->_elType][i][j]];
        }
      }

      for(unsigned i = 0; i < 4; i++) {
        std::vector<double> At;                             // conic cofficient in the target element
        GetConicsInTargetElement(_yr[_data->_elType][i], Ar, At);
        std::tuple <double, double, double> a = AdaptiveRefinement(levelMax, _y[level - 1][i], _yi[level - 1][i], At, level + 1);
        area1 += std::get<0>(a);
        area2 += std::get<1>(a);
        arcLenght += std::get<2>(a);

      }
    }
    else {
      //PrintElement(level, j, x);

      if(test == -1 || test == 1) { // it means it is a full element


        const elem_type *femL = _fem1->GetFiniteElement(_data->_elType, 0);
        const elem_type *femV = _fem1->GetFiniteElement(_data->_elType, _data->_VType);
        const elem_type *femP = _fem1->GetFiniteElement(_data->_elType, _data->_PType);

        std::vector<std::vector<double>> xl; // phisical coordinates at the nodes of l-mesh
        if(_data->_elType == 3) {
          xl = {{x[0][0], x[1][0], x[2][0], x[3][0]}, {x[0][1], x[1][1], x[2][1], x[3][1]}};
        }
        else if(_data->_elType == 4) {
          xl = {{x[0][0], x[1][0], x[2][0]}, {x[0][1], x[1][1], x[2][1]}};
        }
        for(unsigned ig = 0; ig < femL->GetGaussPointNumber(); ig++) { //l-mesh gauss loop
          femL->Jacobian(xl, ig, _weight, _phi, _phix); // _phi, and _phix are the l-mesh test function and gradient at the gauss point of the l-mesh
          if(test == -1) { // inside
            unsigned dim = 2;
            std::vector <double> xil0g(2, 0);  // get the l0 parent coordinates at the gauss point l-mesh
            for(unsigned i = 0; i < _phi.size(); i++) {
              for(unsigned k = 0; k < dim; k++) {
                xil0g[k] += _phi[i] * xil0[i][k];
              }
            }
            double gaussPointWeight;
            femV->Jacobian(_data->_xv, xil0g, gaussPointWeight, _phiV, _phiVx);  //_xv are the phisical coordinates at nodes of l0-mesh
            // _phiV, and _phiVx are the l0 test function and gradient at the gauss point of the l-mesh


            // std::cout<<xil0g[0]<<" "<<xil0g[1]<<std::endl;

            std::vector <double> xg(2, 0); // get the phisical coordinate at the gauss point of the l-mesh, using only the information at the l0-mesh
            for(unsigned i = 0; i < _phiV.size(); i++) {
              for(unsigned k = 0; k < dim; k++) {
                xg[k] += _phiV[i] * _data->_xv[k][i];
              }
            }

            std::cout<<xg[0]<<" "<<xg[1]<<std::endl;

            area1 += (xg[0] * xg[0] + xg[1] * xg[1]) * _weight;
            //area1 += _weight;

            AssembleNavierStokes(_data, _phiV, _phiVx, _phiV, /*  //TODO std::vector <double> &phiP,*/
                          _weight, 1., 0., 0.,
                          {0.,0.,}, 0., 0., 0.0000000001);

          }
          else {
            area2 += _weight;


            AssembleNavierStokes(_data, _phiV, _phiVx, _phiV, /*  //TODO std::vector <double> &phiP,*/
                          _weight, 0., 1., 0.,
                          {0.,0.,}, 0., 0., 0.0000000001);


          }
        }
      }
    }
  }
  else {
    std::vector<std::vector<double>> xl; // phisical coordinates at the nodes of l-mesh
    if(_data->_elType == 3) {
      xl = {{x[0][0], x[1][0], x[2][0], x[3][0]}, {x[0][1], x[1][1], x[2][1], x[3][1]}};
    }
    else if(_data->_elType == 4) {
      xl = {{x[0][0], x[1][0], x[2][0]}, {x[0][1], x[1][1], x[2][1]}};
    }

    std::vector<double> B;
    this->BestFitLinearInterpolation(Ar, B);
    std::vector<double> weight1;
    std::vector<double> weight2;
    std::vector<double> weightI;
    //_weightCF[_data->_elType]->GetWeightWithMap(0, {-B[0], -B[1]}, -B[2], weight1);
    //_weightCF[_data->_elType]->GetWeightWithMap(0, {B[0], B[1]}, B[2], weight2);

    _weightCD[_data->_elType]->GetWeight({-B[0], -B[1]}, -B[2], weight1);
    _weightCD[_data->_elType]->GetWeight({ B[0],  B[1]},  B[2], weight2);

    _weightCF[_data->_elType]->GetWeightWithMap(-1, {B[0], B[1]}, B[2], weightI);
    //_quadCDI->GetWeight({ B[0],  B[1]},  B[2], weightI);

    const elem_type *femL = _fem2->GetFiniteElement(_data->_elType, 0);
    const elem_type *femV = _fem2->GetFiniteElement(_data->_elType, _data->_VType);
    const elem_type *femP = _fem2->GetFiniteElement(_data->_elType, _data->_PType);

    unsigned dim = 2;
    for(unsigned ig = 0; ig < femL->GetGaussPointNumber(); ig++) {

      std::vector<std::vector<double>> J, Ji;

      femL->GetJacobianMatrix(xl, ig, _weight, J, Ji);

      double dsN = 0.;
      std::vector<double> Nf(dim, 0.);
      for(unsigned k = 0; k < dim; k++) {
        for(unsigned j = 0; j < dim; j++) {
          Nf[k] += Ji[j][k] * B[j];
          //Nf[k] += Ji[j][k] * Nr[j];
        }
        dsN += Nf[k] * Nf[k];
      }
      dsN = sqrt(dsN);
      for(unsigned k = 0; k < dim; k++) {
        Nf[k] /= dsN;
      }

      femL->Jacobian(xl, ig, _weight, _phi, _phix);

      //_phip = femL->GetPhi(ig);

      std::vector <double> xil0g(2, 0);
      for(unsigned i = 0; i < _phi.size(); i++) {
        for(unsigned k = 0; k < dim; k++) {
          xil0g[k] += _phi[i] * xil0[i][k];
        }
      }
      //std::cout<<xil0g[0]<<" "<<xil0g[1]<<std::endl;

      double gaussPointWeight;
      femV->Jacobian(_data->_xv, xil0g, gaussPointWeight, _phiV, _phiVx);

      std::vector <double> xg(2, 0);
      for(unsigned i = 0; i < _phiV.size(); i++) {
        for(unsigned k = 0; k < dim; k++) {
          xg[k] += _phiV[i] * _data->_xv[k][i];
        }
      }

      //std::cout<<xg[0]<<" "<<xg[1]<<std::endl;

      arcLenght += dsN * _weight * weightI[ig];
      area1 += (xg[0] * xg[0] + xg[1] * xg[1]) * _weight * weight1[ig];
      //area1 += _weight * weight1[ig];
      area2 += _weight * weight2[ig];

       AssembleNavierStokes(_data, _phiV, _phiVx, _phiV, /*  //TODO std::vector <double> &phiP,*/
                          _weight, weight1[ig], weight2[ig], weightI[ig],
                          Nf, 0 /*curvature*/, dsN, 0);


    }

  }

  return tuple<double, double, double>(area1, area2, arcLenght);
}

void ConicAdaptiveRefinement::BestFitLinearInterpolation(const std::vector<double> &A, std::vector<double> &B) {

  std::vector<std::vector <double>> M(3, std::vector<double>(3, 0)) ;
  std::vector <double> F(3, 0);

  unsigned n = 10;
  double h, x0, y0;
  unsigned jj;
  if (_data->_elType == 3){
    h = 2. / n;
    x0 = -1;
    y0 = -1;
    jj = 0;
  }
  else if (_data->_elType == 4){
    h = 1. / n;
    x0 = 0;
    y0 = 0;
    jj = 1;
  }
  else{
    std::cout << "Error this element type has not been implemented yet\n;";
    abort();
  }

  double s2 = 0;

  double x = x0;
  for(unsigned i = 0; i <= n ; i++, x += h) {
    double y = y0;
    double z = EvaluateConic({x, y}, A);
    int nj = n;
    for(int j = 0; j <= nj ; j++, y += h, nj -= jj ) {
      s2 += z * z;
      z += h * (A[1] * x + A[2] * (2. * y + h) + A[4]); // += increase on the levelSet quadratic function for the next iteration
    }
  }
  s2 /= n;

  x = x0;
  for(unsigned i = 0; i <= n ; i++, x += h) {
    double y = y0;
    double z = EvaluateConic({x, y}, A);
    int nj = n;
    for(int j = 0; j <= nj ; j++, y += h,  nj -= jj) {

      double w = exp(-100. * z * z / s2);
      double wx = w * x;
      double wy = w * y;

      M[0][0] += wx * x;
      M[0][1] += wx * y;
      M[0][2] += wx;

      M[1][1] += wy * y;
      M[1][2] += wy;

      M[2][2] += w;

      F[0] += wx * z;
      F[1] += wy * z;
      F[2] += w * z;

      z += h * (A[1] * x + A[2] * (2. * y + h) + A[4]); // += increase on the levelSet quadratic function for the next iteration

    }
  }

  M[1][0] = M[0][1];
  M[2][0] = M[0][2];
  M[2][1] = M[1][2];

  double det = M[0][0] * (M[1][2] * M[1][2] - M[1][1] * M[2][2]) +
               M[0][1] * (M[1][0] * M[2][2] - M[2][1] * M[0][2]) +
               M[0][2] * (M[1][1] * M[0][2] - M[0][1] * M[1][2]) ;



  std::vector<std::vector <double>> Mi = {
    {(M[1][2] * M[1][2] - M[1][1] * M[2][2]) / det, 0., 0.},
    {(M[0][1] * M[2][2] - M[0][2] * M[1][2]) / det, (M[0][2] * M[0][2] - M[0][0] * M[2][2]) / det, 0.},
    {(M[0][2] * M[1][1] - M[0][1] * M[1][2]) / det, (M[0][0] * M[1][2] - M[0][1] * M[0][2]) / det, (M[0][1] * M[0][1] - M[0][0] * M[1][1]) / det}
  };
  Mi[0][1] = Mi[1][0];
  Mi[0][2] = Mi[2][0];
  Mi[1][2] = Mi[2][1];

  B.assign(3, 0);

  for(unsigned i = 0; i < 3; i++) {
    for(unsigned j = 0; j < 3; j++) {
      B[i] += Mi[i][j] * F[j];
    }
  }

  det = sqrt(B[0] * B[0] + B[1] * B[1]);
  B[0] /= det;
  B[1] /= det;
  B[2] /= det;

}





#endif





