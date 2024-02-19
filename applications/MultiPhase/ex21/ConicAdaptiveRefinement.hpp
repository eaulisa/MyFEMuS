
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

class ConicAdaptiveRefinement {
  public:
    ConicAdaptiveRefinement() {
      _xr = {{-1., -1.}, {1., 1.}, {1., -1.}, {-1., 1.},
        {-1., 0.}, {1., 0.}, {0., -1.}, { 0., 1.}, {0., 0.}
      };
      _wr = {0.0625, 0.0625, 0.0625, 0.0625, 0.125, 0.125, 0.125, 0.125, 0.25};

      _yr = {
        {{-1, -1}, {0, -1}, {0, 0}, {-1, 0}},
        {{0, -1}, {1, -1}, {1, 0}, {0, 0}},
        {{0, 0}, {1, 0}, {1, 1}, {0, 1}},
        {{-1, 0}, {0, 0}, {0, 1}, {-1, 1}}
      };
      // coordinates in the reference element of the nodes of the target elements

      _quadCF = new CutFemWeight <double, TypeA > (QUAD, 5, "legendre");

      _fem1 = new Fem(3, 2);
      _fem2 = new Fem(_quadCF->GetGaussQuadratureOrder(), _quadCF->GetDimension());
      _quad1 = _fem1->GetFiniteElement(3, 0); //quad linear fem for coarse integration
      _quad2 = _fem2->GetFiniteElement(3, 0); //quad linear fem for fine integration

      double dx = .025;
      double dtetha = 1.;

      _quadCD = new CDWeightQUAD <TypeA> (5, dx, dtetha);

      dx = 0.00625;
      dtetha = .25;
      _quadCDI = new CDWeightQUAD <TypeA> (5, dx, dtetha, true);

    }

    ~ConicAdaptiveRefinement() {
      delete _quadCF;
      delete _quadCD;
      delete _quadCDI;
      delete _fem1;
      delete _fem2;
    };

    void CalculateConicsInTargetElement(const std::vector<std::vector<double>> &x, const std::vector<double> &A, std::vector<double> &A1);
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


    double GetVolumeFraction(const std::vector<double>&a);

    std::tuple<double, double, double> AdaptiveRefinement(const unsigned &level,  const unsigned &j,
                                                          const unsigned &levelMax, const std::vector<std::vector<double>>&x,
                                                          const std::vector<std::vector<double>>&xi,
                                                          const std::vector<double>&Ar);

    void PrintElement(const unsigned &level, const unsigned &jp, const std::vector<std::vector<double>>&x) {
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

    bool CheckIfRootsAreInBetweenM1andP1(const double &A, const double &B, const double &C) {
      double det = sqrt(B * B - 4. * A * C);
      if(det >= 0.) {
        double t = (-B - det) / (2. * A);
        if(t >= -1. && t <= 1.) return true;
        t = (-B + det) / (2. * A);
        if(t >= -1. && t <= 1.) return true;
      }
      return false;
    }

  private:
    std::vector<double> _wr;
    std::vector<std::vector<double>> _xr;
    std::vector<std::vector<std::vector<double>>> _yr;
    std::vector<std::vector<std::vector<std::vector<double>>>> _y;
    std::vector<std::vector<std::vector<std::vector<double>>>> _yi;
    Fem* _fem1, *_fem2;
    const elem_type *_quad1, *_quad2;
    double _weight;
    std::vector <double> _phi;
    std::vector <double> _phix;

    double _weightl0;
    std::vector <double> _phil0;
    std::vector <double> _phil0x;


    std::vector<std::vector<double>> _xl0;

    std::vector<double> _xg;

    CutFemWeight <double, TypeA> *_quadCF;

    CDWeightQUAD <TypeA> *_quadCD;
    CDWeightQUAD <TypeA> *_quadCDI;


};


//coefficients of conic in the reference system
void ConicAdaptiveRefinement::CalculateConicsInTargetElement(const std::vector<std::vector<double>> &x, const std::vector<double> &A, std::vector<double> &A1) {

  const double &x1 = x[0][0];
  const double &x4 = x[3][0];

  const double &y1 = x[0][1];
  const double &y2 = x[1][1];

  double Lx = x[1][0] - x1;
  double Ly = x[3][1] - y1;

  const double &a = A[0];
  const double &b = A[1];
  const double &c = A[2];
  const double &d = A[3];
  const double &e = A[4];
  const double &f = A[5];

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
//Jacobian of trasformation
void ConicAdaptiveRefinement::ComputeJacobian(const std::vector<std::vector<double>> &x, std::vector<std::vector<double>> &J) {
  const double &x1 = x[0][0];
  const double &x2 = x[1][0];
  const double &x4 = x[3][0];

  const double &y1 = x[0][1];
  const double &y2 = x[1][1];
  const double &y4 = x[3][1];

  J[0][0] = 0.5 * (x2 - x1);                                //partial derivative of x wrt xi
  J[0][1] = 0.5 * (x4 - x1);                                //partial derivative of x wrt eta
  J[1][0] = 0.5 * (y2 - y1);                                //partial derivative of y wrt xi
  J[1][1] = 0.5 * (y4 - y1);                                //partial derivative of y wrt eta
}

//Inverse of Jacobian calculation
void ConicAdaptiveRefinement::ComputeInverseJacobian(std::vector<std::vector<double>> &J, std::vector<std::vector<double>> &IJ) {
  double detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];

  IJ.resize(2, std::vector<double>(2));

  if(std::abs(detJ) > 1e-10) {
    IJ[0][0] = J[1][1] / detJ;
    IJ[0][1] = -J[0][1] / detJ;
    IJ[1][0] = -J[1][0] / detJ;
    IJ[1][1] = J[0][0] / detJ;
  }
  else {
    std::cout << " Error: Jacobian Matrix is singular." << std::endl;
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

  //main diagonal
  A = a + b + c;
  B =  d + e;
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





double ConicAdaptiveRefinement::GetVolumeFraction(const std::vector<double>&a) {

  double C = 0.;
  for(unsigned i = 0; i < _xr.size(); i++) {
    if(EvaluateConic(_xr[i], a) < 0) C += _wr[i];
  }
  return C;
}


std::tuple<double, double, double> ConicAdaptiveRefinement::AdaptiveRefinement(
  const unsigned &level, // mylevel, with initial level = 1
  const unsigned &j, // son number with respect to the father, for level = 1, is only 1
  const unsigned &levelMax,
  const std::vector<std::vector<double>>&x, // physical coordinates at the nodes of the l-mesh
  const std::vector<std::vector<double>>&xil0,// l0 parent coordinates at the nodes of the l-mesh
  const std::vector<double>&Ar) { // myconic

  double area1 = 0.;
  double area2 = 0.;
  double arcLenght = 0.;

  const double &x1 = x[0][0];
  const double &x2 = x[1][0];
  const double &x3 = x[2][0];
  const double &x4 = x[3][0];

  const double &y1 = x[0][1];
  const double &y2 = x[1][1];
  const double &y3 = x[2][1];
  const double &y4 = x[3][1];

  const double &xi1 = xil0[0][0];
  const double &xi2 = xil0[1][0];
  const double &xi3 = xil0[2][0];
  const double &xi4 = xil0[3][0];

  const double &yi1 = xil0[0][1];
  const double &yi2 = xil0[1][1];
  const double &yi3 = xil0[2][1];
  const double &yi4 = xil0[3][1];

  if(level == 1) {
    _xl0 = {{x1, x2, x3, x4}, {y1, y2, y3, y4}}; //init the phisical coordinates at the nodes of the l0-mesh
    _y.resize(levelMax);
    _yi.resize(levelMax);
  }

  if(level < levelMax) {
    int test = TestIfIntesectionWithReferenceQuad(Ar);
    if(test == 0) { // it means there is an intersection

      std::vector<double> x5 = {0.5 * (x1 + x2), 0.5 * (y1 + y2)};
      std::vector<double> x6 = {0.5 * (x2 + x3), 0.5 * (y2 + y3)};
      std::vector<double> x7 = {0.5 * (x3 + x4), 0.5 * (y3 + y4)};
      std::vector<double> x8 = {0.5 * (x1 + x4), 0.5 * (y1 + y4)};
      std::vector<double> x9 = {0.5 * (x5[0] + x7[0]), 0.5 * (x5[1] + x7[1])};
      _y[level - 1] = {
        {{x1, y1}, x5, x9, x8},
        {x5, {x2, y2}, x6, x9},
        {x9, x6, {x3, y3}, x7},
        {x8, x9, x7, {x4, y4}}
      }; // physical coordinates of my children
      x5 = {0.5 * (xi1 + xi2), 0.5 * (yi1 + yi2)};
      x6 = {0.5 * (xi2 + xi3), 0.5 * (yi2 + yi3)};
      x7 = {0.5 * (xi3 + xi4), 0.5 * (yi3 + yi4)};
      x8 = {0.5 * (xi1 + xi4), 0.5 * (yi1 + yi4)};
      x9 = {0.5 * (x5[0] + x7[0]), 0.5 * (x5[1] + x7[1])};
      _yi[level - 1] = {
        {{xi1, yi1}, x5, x9, x8},
        {x5, {xi2, yi2}, x6, x9},
        {x9, x6, {xi3, yi3}, x7},
        {x8, x9, x7, {xi4, yi4}}
      };



      for(unsigned i = 0; i < 4; i++) {
        std::vector<double> At;                             // conic cofficient in the target element
        CalculateConicsInTargetElement(_yr[i], Ar, At);
        std::tuple <double, double, double> a = AdaptiveRefinement(level + 1, i + 1, levelMax, _y[level - 1][i], _yi[level - 1][i], At);
        area1 += std::get<0>(a);
        area2 += std::get<1>(a);
        arcLenght += std::get<2>(a);

      }
    }
    else {
      //PrintElement(level, j, x);

      if(test == -1 || test == 1) { // it means it is a full element

        std::vector<std::vector<double>> xl = {{x1, x2, x3, x4}, {y1, y2, y3, y4}}; // phisical coordinates at the nodes of l-mesh
        for(unsigned ig = 0; ig < _quad1->GetGaussPointNumber(); ig++) { //l-mesh gauss loop
          _quad1->Jacobian(xl, ig, _weight, _phi, _phix); // _phi, and _phix are the l-mesh test function and gradient at the gauss point of the l-mesh
          if(test == -1) { // inside
            unsigned dim = 2;
            std::vector <double> xil0g(2, 0);  // get the l0 parent coordinates at the gauss point l-mesh
            for(unsigned i = 0; i < _phi.size(); i++) {
              for(unsigned k = 0; k < dim; k++) {
                xil0g[k] += _phi[i] * xil0[i][k];
              }
            }
            _quad2->Jacobian(_xl0, xil0g, _weightl0 /*_weightl0 has no meaning*/, _phil0, _phil0x);  //_xl0 are the phisical coordinates at nodes of l0-mesh
            // _phil0, and _phil0x are the l0 test function and gradient at the gauss point of the l-mesh


            std::vector <double> xg(2, 0); // get the phisical coordinate at the gauss point of the l-mesh, using only the information at the l0-mesh
            for(unsigned i = 0; i < _phil0.size(); i++) {
              for(unsigned k = 0; k < dim; k++) {
                xg[k] += _phil0[i] * _xl0[k][i];
              }
            }
            area1 += (xg[0] * xg[0] + xg[1] * xg[1]) * _weight;
          }
          else {
            area2 += _weight;
          }
        }
      }
    }
  }
  else {

    std::vector<std::vector<double>> xl = {{x1, x2, x3, x4}, {y1, y2, y3, y4}};

    std::vector<double> B;
    this->BestFitLinearInterpolation(Ar, B);
    std::vector<double> weight1;
    std::vector<double> weight2;
    std::vector<double> weightI;
    //quadCF->GetWeightWithMap(0, {-B[0], -B[1]}, -B[2], weight1);
    //_quadCF->GetWeightWithMap(0, {B[0], B[1]}, B[2], weight2);

    _quadCD->GetWeight({-B[0], -B[1]}, -B[2], weight1);
    _quadCD->GetWeight({ B[0],  B[1]},  B[2], weight2);

    _quadCF->GetWeightWithMap(-1, {B[0], B[1]}, B[2], weightI);
    //_quadCDI->GetWeight({ B[0],  B[1]},  B[2], weightI);

    unsigned dim = 2;
    for(unsigned ig = 0; ig < _quad2->GetGaussPointNumber(); ig++) {

      std::vector<std::vector<double>> J, Ji;

      _quad2->GetJacobianMatrix(xl, ig, _weight, J, Ji);
      _quad2->Jacobian(xl, ig, _weight, _phi, _phix);

      std::vector <double> xil0g(2, 0);
      for(unsigned i = 0; i < _phi.size(); i++) {
        for(unsigned k = 0; k < dim; k++) {
          xil0g[k] += _phi[i] * xil0[i][k];
        }
      }
      _quad2->Jacobian(_xl0, xil0g, _weightl0, _phil0, _phil0x);

      std::vector <double> xg(2, 0);
      for(unsigned i = 0; i < _phil0.size(); i++) {
        for(unsigned k = 0; k < dim; k++) {
          xg[k] += _phil0[i] * _xl0[k][i];
        }
      }

      // std::vector<double> Nr;
      // this->EvaluateConicNormal(xi, Ar, Nr);

      //std::cout << ig<<" "<<Nr[0] << " " << B[0] << " " << Nr[1] <<" "<<B[1]<<std::endl;


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

      arcLenght += dsN * _weight * weightI[ig];
      area1 += (xg[0] * xg[0] + xg[1] * xg[1]) * _weight * weight1[ig];
      area2 += _weight * weight2[ig];

    }

  }

  return tuple<double, double, double>(area1, area2, arcLenght);
}


void ConicAdaptiveRefinement::BestFitLinearInterpolation(const std::vector<double> &A, std::vector<double> &B) {

  std::vector<std::vector <double>> M(3, std::vector<double>(3, 0)) ;
  std::vector <double> F(3, 0);

  unsigned n = 10;
  double h = 2. / n;

  double s2 = 0;

  double x = -1.;
  for(unsigned i = 0; i <= n ; i++, x += h) {
    double y = -1.;
    double z = EvaluateConic({x, y}, A);
    for(unsigned j = 0; j <= n ; j++, y += h) {
      s2 += z * z;
      z += h * (A[1] * x + A[2] * (2. * y + h) + A[4]); // += increase on the levelSet quadratic function for the next iteration
    }
  }
  s2 /= n;

  x = -1.;
  for(unsigned i = 0; i <= n ; i++, x += h) {
    double y = -1.;
    double z = EvaluateConic({x, y}, A);
    for(unsigned j = 0; j <= n ; j++, y += h) {

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
