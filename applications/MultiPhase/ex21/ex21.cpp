#include "FemusInit.hpp"

#include "ConicAdaptiveRefinement.hpp"

using namespace femus;

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
  //example inputs
  std::vector<std::vector<double>> y = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};

  //Matrix to store calculated coefficients
  //Coeficients of conics in physical system
  std::vector <double> A = {4., 0, 2., -0.8, 0.8, -0.88};
  std::vector <double> Ap;

  ConicAdaptiveRefinement cad;

  cad.CalculateConicsInTargetElement(y, A, Ap);
  double area = cad.AdaptiveRefinement(1,  1, 10, y, Ap);

  std::cout.precision(14);
  std::cout << " Analytic Area = " << M_PI * 0.5 * sqrt(2.) / 2. << std::endl;
  std::cout << " Computed Area = " << area << std::endl;

  return 0;
}
