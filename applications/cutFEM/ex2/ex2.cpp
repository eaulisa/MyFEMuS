
#include "../eqPoly/eqPoly.hpp"
#include <time.h>
#include "FemusInit.hpp"
#include "Elem.hpp"
#include <stdlib.h>
#include "../eqPoly/bestFit.hpp"

using namespace femus;
using namespace std;
using std::cout;

int main(int argc, char** args) {

  //FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  //SetCoefficients takes ( dim, degree of equivalent polynomial, rho, vector < a, b, c > for the dicontinuity ( point, line, or plane ), and element (3=triangle/tet, 4=square/cube
  // a*x + b*y + c*z + d = 0, and d) as inputs
  EquivalentPolynomial eqP;
  BestFit bf;

  //
   eqP.SetCoefficients(3, 2, 2, std::vector<double> {1., 2., 1., 0.}, 3);
   eqP.PrintCoefficients();
//   std::cout << eqP.GetValue(std::vector<double> {0.5, 0.5}) << " " << std::endl;

  std::vector < double >points {1.,2.,4.,5.,7.,8.,-5.,23.,12.,15.,14.,14.};
  //std::vector < double >points(2000);
  std::vector < double >normal {1., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0.};
  unsigned dim = 3;
  std::vector < double > equation(dim + 1, 0.);

  equation = bf.FindBestFit(points, normal, dim);
  std::vector < double > onepoint {1.,2.,3.};
  unsigned element = 0;
  
  for(unsigned i = 0; i < dim +1; i++) {
        std::cout << equation[i] << "  best fit" << endl;

      
  }

//   clock_t t;
//   t = clock();
// 
//   for(unsigned k = 0; k < 1000; k++) {
// 
//     for(unsigned i = 0; i < 2000; i++) {
// 
//       points[i] = rand() % 100;
// 
//     }
// 
//     
//     eqP.FindBestFit(points, normal, dim);
// 
// 
//   }
// 
//   
// 
//   t = clock() - t;
//   printf("It took me %d clicks (%f seconds).\n", t, ((float)t) / CLOCKS_PER_SEC);


  /*std::vector <double> tes(27);
  MatrixXd C;
  C.setRandom(10,3);
  JacobiSVD<MatrixXd> svd( C, ComputeThinU | ComputeThinV);
  MatrixXd Cp = svd.matrixV();
  MatrixXd sigma = svd.singularValues().asDiagonal();
  svd.matrixV().transpose();
  MatrixXd diff = Cp;
  tes[0] = Cp(1,2);

  cout << "diff:\n" << Cp.col(2) << "\n";
  cout << "diff:\n" << tes[0] << "\n";
  */


  
  std::vector<double> phi;
  std::vector<double> gradPhi;
  double weight;
  {
    std::vector<std::vector<double>> xv = {{-1., 1.}};
    double integral = 0.;
    unsigned dim = 1;
    element = 6;
    eqP.SetCoefficients(1, 3, 50, std::vector<double> {1., -0.5},  4);
    const elem_type * fe = new const elem_type_1D("line", "linear", "ninth");
    for(unsigned ig = 0; ig < fe->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      fe->Jacobian(xv, ig, weight, phi, gradPhi);
      std::vector<double> xg(dim, 0.);
      for(unsigned i = 0; i < phi.size(); i++) {
        for(unsigned k = 0; k < dim; k++) {
          xg[k] += xv[k][i] * phi[i];
        }
      }
      integral += xg[0] * xg[0] * xg[0] * eqP.GetValue(xg, element) * weight;
    }
    std::cout << "Integral = " << integral << std::endl;
    delete fe;
  }

  element = 4;
  std::vector  <double> pt {0.5, 0.5};
  eqP.SetCoefficients(2, 2, 10, std::vector<double> {2., 1., -1.},  element);
  std::cout << eqP.GetValue(pt, element) << " triangle at 0.5, 0.5 " << std::endl;
     eqP.PrintCoefficients();

     //TODO Fix triangle integration
  
  {
    std::vector<std::vector<double>> xv = {{-1., 1., 1., -1.}, {-1., -1., 1., 1.}};
    double integral = 0.;
    unsigned dim = 2;
    element = 4;
    eqP.SetCoefficients(2, 2, 10, std::vector<double> {2., 1., -1}, element);
    const elem_type * fe = new const elem_type_2D("tri", "linear", "ninth");
    for(unsigned ig = 0; ig < fe->GetGaussPointNumber(); ig++) {

      fe->Jacobian(xv, ig, weight, phi, gradPhi);
      std::vector<double> xg(dim, 0.);
      for(unsigned i = 0; i < phi.size(); i++) {
        for(unsigned k = 0; k < dim; k++) {
          xg[k] += xv[k][i] * phi[i];
        }
      }
      integral += xg[0] * xg[1] * eqP.GetValue(xg, element) * weight;
    }
    std::cout << "Integral = " << integral << std::endl;
    delete fe;
  }





  return 0;
}


