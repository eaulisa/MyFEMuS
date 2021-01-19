
#include "eqPoly.hpp"

#include "FemusInit.hpp"
#include "Elem.hpp"

//using namespace femus;
using namespace std;

int main(int argc, char** args) {

  //FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  //SetCoefficients takes ( dim, degree of equivalent polynomial, rho, vector < a, b, c > for the dicontinuity ( point, line, or plane ) 
  // a*x + b*y + c*z + d = 0, and d) as inputs
  EquivalentPolynomial eqP;
  eqP.SetCoefficients(2, 4, 20, std::vector<double> {1., 1.}, -1.);
  eqP.PrintCoefficients();
  std::cout << eqP.GetValue(std::vector<double> {0.5,0.5}) << " " << std::endl;






/*

  std::vector<double> phi;
  std::vector<double> gradPhi;
  double weight;
  {
    std::vector<std::vector<double>> xv = {{-1., 1.}};
    double integral = 0.;
    unsigned dim = 1;
    eqP.SetCoefficients(1, 3, 100, std::vector<double> {}, 0.5);
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
      integral += xg[0] * xg[0] * xg[0] * eqP.GetValue(xg) * weight;
    }
    std::cout << "Integral = " << integral << std::endl;
    delete fe;
  }

  {
    std::vector<std::vector<double>> xv = {{-1., 1., 1., -1.}, {-1., -1., 1., 1.}};
    double integral = 0.;
    unsigned dim = 2;
    eqP.SetCoefficients(2, 2, 20, std::vector<double> {1., -1.}, 1.);
    const elem_type * fe = new const elem_type_2D("quad", "linear", "ninth");
    for(unsigned ig = 0; ig < fe->GetGaussPointNumber(); ig++) {

      fe->Jacobian(xv, ig, weight, phi, gradPhi);
      std::vector<double> xg(dim, 0.);
      for(unsigned i = 0; i < phi.size(); i++) {
        for(unsigned k = 0; k < dim; k++) {
          xg[k] += xv[k][i] * phi[i];
        }
      }
      integral += xg[0] * xg[1] * eqP.GetValue(xg) * weight;
    }
    std::cout << "Integral = " << integral << std::endl;
    delete fe;
  }
*/




  return 0;
}

