#include "FemusInit.hpp"

#include "ConicAdaptiveRefinement.hpp"

using namespace femus;

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
  //example inputs
  std::vector<std::vector<double>> y = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
  std::vector<std::vector<double>> yi = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};

  //Matrix to store calculated coefficients
  //Coeficients of conics in physical system
  std::vector <double> A = {1, 0, 1., 0, 0, -0.5};
  std::vector <double> Ap;

  ConicAdaptiveRefinement cad;

  // std::vector <double> B(3);
  // cad.BestFitLinearInterpolation({-1, 3, 2, 4, 3, -1}, B);
  //
  //
  // std::cout << "a = "<< B[0] << "; b = " << B[1] << "; c = " << B[2] <<";"<< std::endl;
  //
  //
  // return 0;

  cad.CalculateConicsInTargetElement(y, A, Ap);
  tuple <double,double,double> a = cad.AdaptiveRefinement(1,  1, 6, y, yi, Ap);

  double exact;

  exact = 0.39269908169872414;

  std::cout.precision(14);
  std::cout << " Analytic Area1 = " << exact << std::endl;
  std::cout << " Computed Area1 = " << std::get<0>(a) << std::endl;
  std::cout << " Relative Error = " << fabs((exact - std::get<0>(a))/exact) << std::endl;

  exact = 4. - M_PI * 0.5;
  std::cout << " Analytic Area2 = " << exact << std::endl;
  std::cout << " Computed Area2 = " << std::get<1>(a) << std::endl;
  std::cout << " Relative Error = " << fabs((exact - std::get<1>(a))/exact) << std::endl;

  exact = 2 * M_PI * sqrt(2.)/2.;
  std::cout << " Analytic perimeter = " << exact << std::endl;
  std::cout << " Computed perimeter = " << std::get<2>(a) << std::endl;
  std::cout << " Relative Error = " << fabs((exact - std::get<2>(a))/exact) << std::endl;


  return 0;
}
