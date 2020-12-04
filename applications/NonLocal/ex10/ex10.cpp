
#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "TransientSystem.hpp"
#include "NonLinearImplicitSystem.hpp"

#include "NumericVector.hpp"
#include "adept.h"

#include "petsc.h"
#include "petscmat.h"
#include "PetscMatrix.hpp"

#include "slepceps.h"

#include "../include/nonlocal_assembly_adaptive.hpp"

//2D NONLOCAL EX : nonlocal diffusion for a body with different material properties

using namespace femus;

double InitalValueU(const std::vector < double >& x) {
  double value;

  //value =  x[0] * x[0] + x[1] * x[1]; //consistency 
  value =  x[0] * x[0] * x[0] + x[1] * x[1] * x[1]; //cubic 
  //value = x[0] * x[0] * x[0] * x[0] + x[1] * x[1] * x[1] * x[1]; //quartic




  return value;
}

void GetL2Norm(MultiLevelSolution & mlSol, MultiLevelSolution & mlSolFine);

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {

  bool dirichlet = true;
  //value = x[0] * x[0] + x[1] * x[1]; //consistency
  value = x[0] * x[0] * x[0] + x[1] * x[1] * x[1]; //cubic
  //value = x[0] * x[0] * x[0] * x[0] + x[1] * x[1] * x[1] * x[1]; //quartic


  return dirichlet;
}

//unsigned numberOfUniformLevels = 2; //consistency
unsigned numberOfUniformLevels = 1; //cubic-quartic 2->6 //cubic Marta4Quad Tri Mix 
//unsigned numberOfUniformLevels = 2; //cubic-quartic 2->4 mappa a 4->6 //cubic Marta4Fine 


unsigned numberOfUniformLevelsFine = 1;

int main(int argc, char** argv) {

  FemusInit mpinit(argc, argv, MPI_COMM_WORLD);
  std::cout << "Hello\n"<<std::flush;
  return 0;

} //end main
