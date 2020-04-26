/* This example is for quasi-conformal minimization */
/* mu controls the beltrami coefficient */

#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "TransientSystem.hpp"
#include "adept.h"
#include <cstdlib>
#include "petsc.h"
#include "petscmat.h"
#include "PetscMatrix.hpp"

unsigned counter = 0;
unsigned numberOfIterations = 1;
using namespace femus;


#include "../include/supportFunctions.hpp"
//#include "../include/updateMu1.hpp"
#include "../include/assembleConformalMinimization5.hpp"

double InitalValueCM(const std::vector < double >& x) {
//   return cos(4.* M_PI * sqrt(x[0] * x[0] + x[1] * x[1])/0.5) ;
  return cos(20 * M_PI * x[0]) + sin(20 * M_PI * x[1]);
  //return sin(28.* M_PI * x[0]) * sin(28.* M_PI * x[1]) + cos(28.* M_PI * x[0]) * cos(28.* M_PI * x[1]) ;
}
double GetTimeStep(const double t) {
  return 1;
}

bool SetBoundaryCondition(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time);

// Main program starts here.
int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);


  unsigned numberOfUniformLevels = 1;
  unsigned numberOfNonLinearSteps = 10;
  numberOfIterations = 1;
 
  // define multilevel mesh
  MultiLevelMesh mlMsh;

  // Read coarse level mesh and generate finer level meshes.
  double scalingFactor = 1.;
  unsigned numberOfSelectiveLevels = 0;
 
  mlMsh.ReadCoarseMesh("../input/hex.neu", "seventh", scalingFactor);
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels, NULL);

  // Erase all the coarse mesh levels.
  mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  const unsigned  DIM = mlMsh.GetDimension();

  // Define the multilevel solution and attach the mlMsh object to it.
  MultiLevelSolution mlSol(&mlMsh);

  // Add variables X,Y,W to mlSol.

  FEOrder feOrder = FIRST;
  mlSol.AddSolution("Dx1", LAGRANGE, feOrder, 2);
  mlSol.AddSolution("Dx2", LAGRANGE, feOrder, 2);
  mlSol.AddSolution("Dx3", LAGRANGE, feOrder, 2);

  mlSol.AddSolution("mu1", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
  mlSol.AddSolution("mu2", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
  mlSol.AddSolution("weight1", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);

  mlSol.AddSolution("mu1Edge", LAGRANGE, SECOND, 0, false);
  mlSol.AddSolution("mu2Edge", LAGRANGE, SECOND, 0, false);
  mlSol.AddSolution("cntEdge", LAGRANGE, SECOND, 0, false);

  mlSol.Initialize("All");


  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("All");

  mlSol.GenerateBdc("Dx1", "Time_dependent");
  mlSol.GenerateBdc("Dx2", "Time_dependent");
  mlSol.GenerateBdc("Dx3", "Time_dependent");
  
  MultiLevelProblem mlProb(&mlSol);

  // Add system Conformal or Shear Minimization in mlProb.
  TransientNonlinearImplicitSystem& system = mlProb.add_system < TransientNonlinearImplicitSystem > ("conformal"); //for conformal

  // Add solutions newDX, Lambda to system.
  system.AddSolutionToSystemPDE("Dx1");
  system.AddSolutionToSystemPDE("Dx2");
  system.AddSolutionToSystemPDE("Dx3");
  
  // Parameters for convergence and # of iterations.
  system.SetMaxNumberOfNonLinearIterations(numberOfNonLinearSteps);
  system.SetNonLinearConvergenceTolerance(1.e-10);

  system.AttachGetTimeIntervalFunction(GetTimeStep);
  system.init();

  mlSol.SetWriter(VTK);

  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");
  mlSol.GetWriter()->SetDebugOutput(true);
  mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);

  system.SetAssembleFunction(AssembleConformalMinimization);

  std::vector<std::string> mov_vars;
  mov_vars.push_back("Dx1");
  mov_vars.push_back("Dx2");
  mov_vars.push_back("Dx3");
  mlSol.GetWriter()->SetMovingMesh(mov_vars);

  system.CopySolutionToOldSolution();
  for(unsigned k = 0; k < numberOfIterations; k++) {
    system.MGsolve();
    //ProjectSolution(mlSol);
    mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, k + 1);
    system.CopySolutionToOldSolution();
  }
  
  return 0;
}

bool SetBoundaryCondition(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time) {

  bool dirichlet = true;
  value = 0.;
  
  if(!strcmp(solName, "Dx2") || !strcmp(solName, "Dx3") ) {
    if(faceName == 1 || faceName == 2)  
    dirichlet = false;
  }
  if(!strcmp(solName, "Dx3") || !strcmp(solName, "Dx1") ) {
    if(faceName == 3 || faceName == 4)  
    dirichlet = false;
  }
  if(!strcmp(solName, "Dx1") || !strcmp(solName, "Dx2") ) {
    if(faceName == 5 || faceName == 6)  
    dirichlet = false;
  }
  return dirichlet;
}
