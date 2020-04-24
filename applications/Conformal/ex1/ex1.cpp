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

const bool read_groups = false;                       //by default, if no argument is given, this is "true"
const bool read_boundary_groups = false;              //by default, if no argument is given, this is "true"

//unsigned muSmoothingType = 1; // the mesh should be logically structured and uniformly oriented
unsigned muSmoothingType = 2; // invariant with respect to quad orientation
unsigned conformalTriangleType = 2;
unsigned counter = 0;
const double eps = 1.e-5;

class parameter {
  public:
    parameter(const unsigned simulation1, const bool surface1,
              const bool constraintIsOn1, const unsigned numberOfUniformLevels1,
              const unsigned numberOfSmoothingSteps1,
              const bool finalSmoothIsOn1, const unsigned numberOfNonLinearSteps1,
              const unsigned numberOfIterations1, const double beltramiNorm1 = 0.
             ):
      simulation(simulation1),
      surface(surface1),
      constraintIsOn(constraintIsOn1),
      numberOfUniformLevels(numberOfUniformLevels1),
      numberOfSmoothingSteps(numberOfSmoothingSteps1),
      finalSmoothIsOn(finalSmoothIsOn1),
      numberOfNonLinearSteps(numberOfSmoothingSteps1),
      numberOfIterations(numberOfIterations1),
      beltramiNorm(beltramiNorm1) {
    }
    unsigned simulation;
    bool surface;
    bool constraintIsOn;
    unsigned numberOfUniformLevels;
    unsigned numberOfSmoothingSteps;
    bool finalSmoothIsOn;
    unsigned numberOfNonLinearSteps;
    unsigned numberOfIterations;
    double beltramiNorm;
};

const parameter squareQuad = parameter(1, false, false, 4, 1, true, 300, 1, 0.811567);
const parameter squareTri = parameter(1, false, false, 4, 1, true, 500, 1, 0.868442);
const parameter cylinderUnconstrained = parameter(2, true, false, 4, 12, false, 30, 1);
const parameter cylinderConstrained = parameter(2, true, true, 4, 12, false, 30, 3, 0.794163);
//parameter intersection = parameter();
//parameter cat = parameter();
//parameter hand = parameter();
//parameter moo = parameter();

unsigned numberOfIterations = 1;
bool surface = true;
bool constraintIsOn = true;
unsigned numberOfSmoothingSteps = 12;
bool finalSmoothIsOn = false;

using namespace femus;

#include "../include/supportFunctions.hpp"
#include "../include/updateMu3.hpp"
#include "../include/assembleConformalMinimization3.hpp"

double InitalValueCM(const std::vector < double >& x) {
//   return cos(4.* M_PI * sqrt(x[0] * x[0] + x[1] * x[1])/0.5) ;
  return cos(20 * M_PI * x[0]) + sin(20 * M_PI * x[1]);
  //return sin(28.* M_PI * x[0]) * sin(28.* M_PI * x[1]) + cos(28.* M_PI * x[0]) * cos(28.* M_PI * x[1]) ;
}
double GetTimeStep(const double t) {
  return 1;
}

bool SetBoundaryConditionZero(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time);
bool SetBoundaryConditionSquare(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time);
bool SetBoundaryConditionCylinder(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time);
bool SetBoundaryConditionIntersection(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time);


// Main program starts here.
int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  unsigned simulation;
  unsigned numberOfUniformLevels = 4;
  unsigned numberOfNonLinearSteps = 30;

  if(argc < 2) {
  simulation0: // square with quads, default
    simulation = 0;
    surface = false;
    constraintIsOn = false;
    numberOfSmoothingSteps = 1;
    finalSmoothIsOn = true;
    numberOfUniformLevels = 4;
    numberOfNonLinearSteps = 300;
  }
  else {
    if(!strcmp("1", args[1])) { // square with triangles
      simulation = 1;
      surface = false;
      constraintIsOn = false;
      numberOfSmoothingSteps = 1;
      finalSmoothIsOn = true;
      numberOfUniformLevels = 4;
      numberOfNonLinearSteps = 500;
    }
    else if(!strcmp("2", args[1])) { // cylinder
      simulation = 2;
      numberOfIterations = 3;
      numberOfUniformLevels = 4;
      surface = true;
      constraintIsOn = true;
      numberOfSmoothingSteps = 12;
      finalSmoothIsOn = false;
      numberOfNonLinearSteps = 30;
    }
    else if(!strcmp("3", args[1])) { // intersection
      simulation = 3;
      numberOfUniformLevels = 1;
      constraintIsOn = false;
    }
    else if(!strcmp("4", args[1])) { // cat
      simulation = 4;
      surface = true;
      constraintIsOn = true;
      numberOfSmoothingSteps = 12;
      finalSmoothIsOn = false;
      numberOfUniformLevels = 1;
      numberOfNonLinearSteps = 3;
    }
    else if(!strcmp("5", args[1])) { // hand
      simulation = 5;
      numberOfUniformLevels = 1;
    }
    else if(!strcmp("6", args[1])) { // moo
      simulation = 6;
      surface = true;
      constraintIsOn = true;
      numberOfSmoothingSteps = 50;
      finalSmoothIsOn = false;
      numberOfUniformLevels = 2;
      numberOfNonLinearSteps = 50;
    }
    else if(!strcmp("7", args[1])) {
      parameter par = cylinderConstrained;
      simulation = par.simulation;
      surface = par.surface;
      constraintIsOn = par.constraintIsOn;
      numberOfUniformLevels = par.numberOfUniformLevels;
      numberOfSmoothingSteps = par.numberOfSmoothingSteps;
      finalSmoothIsOn = par.finalSmoothIsOn;
      numberOfNonLinearSteps = par.numberOfNonLinearSteps;
      numberOfIterations = par.numberOfIterations;
    }
    else {
      goto simulation0;
    }
  }

  // define multilevel mesh
  MultiLevelMesh mlMsh;

  // Read coarse level mesh and generate finer level meshes.
  double scalingFactor = 1.;
  unsigned numberOfSelectiveLevels = 0;
  if(simulation == 0) {
    mlMsh.ReadCoarseMesh("../input/square.neu", "seventh", scalingFactor);
    //mlMsh.ReadCoarseMesh("../input/square1.neu", "seventh", scalingFactor);
  }
  else if(simulation == 1) {
    mlMsh.ReadCoarseMesh("../input/squareTri3D.neu", "seventh", scalingFactor);
  }
  else if(simulation == 2) {
    mlMsh.ReadCoarseMesh("../input/cylinder2.neu", "seventh", scalingFactor);
  }
  else if(simulation == 3) {
    mlMsh.ReadCoarseMesh("../input/intersection.neu", "seventh", scalingFactor);
  }
  else if(simulation == 4) {
    mlMsh.ReadCoarseMesh("../input/cat.med", "seventh", scalingFactor, read_groups, read_boundary_groups);
  }
  else if(simulation == 5) {
    mlMsh.ReadCoarseMesh("../input/hand.med", "seventh", scalingFactor);
  }
  else if(simulation == 6) {
    mlMsh.ReadCoarseMesh("../../Willmore/WillmoreSurface/input/moo.med", "seventh", scalingFactor, read_groups, read_boundary_groups);
  }

  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

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
  if(surface) mlSol.AddSolution("Dx3", LAGRANGE, feOrder, 2);

  if(constraintIsOn) mlSol.AddSolution("Lambda", LAGRANGE, feOrder, 0);

  mlSol.AddSolution("ENVN", LAGRANGE, feOrder, 0, false);

  mlSol.AddSolution("mu1", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
  mlSol.AddSolution("mu2", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
  mlSol.AddSolution("weight1", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);

  mlSol.AddSolution("mu1Edge", LAGRANGE, SECOND, 0, false);
  mlSol.AddSolution("mu2Edge", LAGRANGE, SECOND, 0, false);
  mlSol.AddSolution("cntEdge", LAGRANGE, SECOND, 0, false);

  mlSol.AddSolution("cm", LAGRANGE, SECOND, 0, false);

  // Initialize the variables and attach boundary conditions.
  mlSol.Initialize("All");

  mlSol.Initialize("cm", InitalValueCM);

  if(simulation < 2) {
    mlSol.AttachSetBoundaryConditionFunction(SetBoundaryConditionSquare);
  }
  else if(simulation == 2) {
    mlSol.AttachSetBoundaryConditionFunction(SetBoundaryConditionCylinder);
  }
  else if(simulation == 3) {
    mlSol.AttachSetBoundaryConditionFunction(SetBoundaryConditionIntersection);
  }
  else {
    mlSol.AttachSetBoundaryConditionFunction(SetBoundaryConditionZero);
  }
  mlSol.GenerateBdc("All");

  mlSol.GenerateBdc("Dx1", "Time_dependent");
  mlSol.GenerateBdc("Dx2", "Time_dependent");
  if(surface) mlSol.GenerateBdc("Dx3", "Time_dependent");

  GetElementNearVertexNumber(mlSol);
  MultiLevelProblem mlProb(&mlSol);

  // Add system Conformal or Shear Minimization in mlProb.
  TransientNonlinearImplicitSystem& system = mlProb.add_system < TransientNonlinearImplicitSystem > ("conformal"); //for conformal

  // Add solutions newDX, Lambda to system.
  system.AddSolutionToSystemPDE("Dx1");
  system.AddSolutionToSystemPDE("Dx2");
  if(surface) system.AddSolutionToSystemPDE("Dx3");
  if(constraintIsOn) system.AddSolutionToSystemPDE("Lambda");

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
  if(surface) mov_vars.push_back("Dx3");
  mlSol.GetWriter()->SetMovingMesh(mov_vars);

  for(unsigned k = 0; k < numberOfIterations; k++) {
    system.CopySolutionToOldSolution();
    system.MGsolve();
    //ProjectSolution(mlSol);
    mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, k + 1);
  }

  return 0;
}

bool SetBoundaryConditionZero(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time) {

  bool dirichlet = true;
  value = 0.;

  if(!strcmp(solName, "Lambda")) {
    dirichlet = false;
  }
  return dirichlet;
}

bool SetBoundaryConditionSquare(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time) {
  bool dirichlet = true;
  value = 0.;

  if(!strcmp(solName, "Dx1")) {
    if(1 == faceName || 3 == faceName) {
      dirichlet = false;
    }
    if(4 == faceName) {
      value = time / numberOfIterations * 0.75 * sin(x[1] / 0.5 * M_PI);
    }
  }
  else if(!strcmp(solName, "Dx2")) {
    if(2 == faceName) {
      dirichlet = false;
    }
  }

  if(!strcmp(solName, "Lambda")) {
    dirichlet = false;
  }
  return dirichlet;
}


bool SetBoundaryConditionCylinder(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time) {

  bool dirichlet = true;
  value = 0.;

  if(!strcmp(solName, "Dx1")) {
    if(1 == faceName) {
      value = time / numberOfIterations * 0.8 * sin(x[1] / 0.5 * M_PI);
    }
  }

  if(!strcmp(solName, "Lambda")) {
    dirichlet = false;
  }
  return dirichlet;
}

bool SetBoundaryConditionIntersection(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time) {

  bool dirichlet = true;
  value = 0.;

  if(1 == faceName || 2 == faceName) {
    if(!strcmp(solName, "Dx1")) {
      if(1 == faceName)
        value = time / numberOfIterations * 1;
      else {
        value = - time / numberOfIterations * 1;
      }
    }
    if(!strcmp(solName, "Dx2")) {
      value = time / numberOfIterations * 0.35 * x[1] / 0.5;
    }
    else if(!strcmp(solName, "Dx3")) {
      value = time / numberOfIterations * 0.75 * x[2] / 0.5;
    }
  }
  else if(3 == faceName || 4 == faceName) {
    if(!strcmp(solName, "Dx1")) {
      value = time / numberOfIterations * 0.45 * x[0] / 0.5;
    }
    else if(!strcmp(solName, "Dx3")) {
      value = time / numberOfIterations * 0.35 * x[2] / 0.5;
    }
  }
  else if(5 == faceName || 6 == faceName) {
    if(!strcmp(solName, "Dx1")) {
      value = time / numberOfIterations * 0.125 * x[0] / 0.5;
    }
    else if(!strcmp(solName, "Dx2")) {
      value = time / numberOfIterations * 0.65 * x[1] / 0.5;
    }
  }

  if(!strcmp(solName, "Lambda")) {
    dirichlet = false;
  }
  return dirichlet;
}



