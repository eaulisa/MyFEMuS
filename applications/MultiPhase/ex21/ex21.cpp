#include "FemusInit.hpp"


#include "ConicAdaptiveRefinement.hpp"

using namespace femus;


int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
  //example inputs
  // std::vector<std::vector<double>> y = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
  // std::vector<std::vector<double>> y = {{1, -1}, {1, 1}, {-1, 1},{-1, -1}};
  // std::vector<std::vector<double>> yi = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};

  unsigned VType = 2, PType = 3;
  unsigned dofsV = 7;
  unsigned dofsP = 1;
  unsigned dofsAll = 2 * dofsV + 2 * dofsP;


  std::vector<std::vector<double>> V(2,std::vector<double>(dofsV, 2.));
  std::vector<double> K1(dofsP, 1.);
  std::vector<double> P1(dofsP, 1.);
  std::vector<double> P2(dofsP, 2.);

  std::vector<double> res(dofsAll, 0.);
  std::vector<double> jac(dofsAll * dofsAll, 0.);

  // unsigned elType = 3;
  //std::vector<std::vector<double>> xv = {{-1., 1., 1., -1., 0., 1., 0., -1., 0.}, {-1., -1., 1., 1., -1., 0., 1., 0., 0.}};
  //std::vector<std::vector<double>> xv = {{1., 1., -1., -1.,  1., 0., -1., 0., 0.}, { -1., 1., 1.,-1., 0., 1., 0.,-1., 0.}};
  unsigned elType = 4;
  std::vector<std::vector<double>> xv = {{-1., 3., -1., 1., 1., -1., 1. / 3.}, {-1., -1., 3., -1., 1., 1., 1. / 3.}};
  //TODO

  double rho1 = 1., rho2 = 2., mu1 = .2, mu2 = 0.4, sigma = 1., dt = 0.01;

  //Matrix to store calculated coefficients
  //Coeficients of conics in physical system
  std::vector <double> A = {1, 0, 1., 0, 0, -0.5};
  std::vector <double> Ap;

  ConicAdaptiveRefinement cad;

  Data *data = new Data(VType, PType, V, P1, P2, K1, res, jac, xv, elType, rho1, rho2, mu1, mu2, sigma, dt, {0,0}, A, false);

  cad.SetDataPointer(data);

  std::vector<std::vector<double>> y;
  cad.GetXInParentElement(xv, y);
  std::vector<std::vector<double>> yi = cad.GetXiInParentElement();

  cad.GetConicsInTargetElement(y, A, Ap);
  tuple <double, double, double> a = cad.AdaptiveRefinement(1, y, yi, Ap);

  delete data;

  double exact;



  exact = 0.39269908169872414;

  //exact = M_PI * 0.5;
  std::cout.precision(14);
  std::cout << " Analytic Area1 = " << exact << std::endl;
  std::cout << " Computed Area1 = " << std::get<0>(a) << std::endl;
  std::cout << " Relative Error = " << fabs((exact - std::get<0>(a)) / exact) << std::endl;

  exact = 8 - M_PI * 0.5;
  std::cout << " Analytic Area2 = " << exact << std::endl;
  std::cout << " Computed Area2 = " << std::get<1>(a) << std::endl;
  std::cout << " Relative Error = " << fabs((exact - std::get<1>(a)) / exact) << std::endl;

  exact = 2 * M_PI * sqrt(2.) / 2.;
  std::cout << " Analytic perimeter = " << exact << std::endl;
  std::cout << " Computed perimeter = " << std::get<2>(a) << std::endl;
  std::cout << " Relative Error = " << fabs((exact - std::get<2>(a)) / exact) << std::endl;
  //

  return 0;
}




