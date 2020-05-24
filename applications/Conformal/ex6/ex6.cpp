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
const double eps = 1.e-5;
bool areaConstraint = false;

using namespace femus;

#include "../include/parameter.hpp"

Parameter parameter;

#include "../include/supportFunctions.hpp"
#include "../include/updateMu1.hpp"
#include "../include/assembleConformalMinimization6.hpp"

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

const Parameter squareQuad = Parameter("square with quads", 0, false, false, 4, 1, true, 300, 1, 0.811569);
const Parameter squareTri = Parameter("square with triangles", 1, false, false, 4, 1, true, 500, 1, 0.868445);
//const Parameter cylinderUnconstrained = Parameter("cylinder unconstrained", 2, true, false, 4, 12, false, 30, 1, 0.910958);
const Parameter cylinderUnconstrained = Parameter("cylinder unconstrained", 2, true, false, 4, 3, true, 90, 1, 0.729612);
const Parameter cylinderConstrained = Parameter("cylinder constrained", 3, true, true, 4, 12, false, 30, 3, 0.793786);
const Parameter intersection = Parameter("intersection", 4, true, false, 2, 100, true, 10, 5, 0.486729);
//const Parameter intersection = Parameter("intersection", 4, true, false, 2, 100, true, 50, 1, 0.674721);
//const Parameter intersection = Parameter("intersection", 4, true, false, 2, 12, false, 30, 1, 0.979639);
const Parameter cat = Parameter("cat", 5, true, true, 1, 12, false, 2, 1, 0.986943);
const Parameter hand = Parameter("hand", 6, true, true, 1, 12, false, 10, 1, 0.580335);
const Parameter moo = Parameter("moo", 7, true, true, 2, 50, true, 50, 1, 0.654910);

// Main program starts here.
int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  if(argc >= 2) {
    if(!strcmp("0", args[1])) { // square with triangles
      parameter = squareQuad;
    }
    else if(!strcmp("1", args[1])) { // square with triangles
      parameter = squareTri;
    }
    else if(!strcmp("2", args[1])) { // cylinder unconstraint
      parameter = cylinderUnconstrained;
    }
    else if(!strcmp("3", args[1])) { // cylinder constraint
      parameter = cylinderConstrained;
    }
    else if(!strcmp("4", args[1])) { // intersection
      parameter = intersection;
    }
    else if(!strcmp("5", args[1])) { // cat
      parameter = cat;
    }
    else if(!strcmp("6", args[1])) { // hand
      parameter = hand;
    }
    else if(!strcmp("7", args[1])) { // moo
      parameter = moo;
    }
    else {
      goto generic;
    }
  }
  else { //generic
  generic:
    unsigned simulation = 0;
    bool surface = false;
    bool constraintIsOn = false;
    unsigned numberOfUniformLevels = 4;
    unsigned numberOfSmoothingSteps = 1;
    bool finalSmoothIsOn = true;
    unsigned numberOfNonLinearSteps = 300;
    unsigned numberOfIterations = 1;
    parameter = Parameter("generic", simulation, surface, constraintIsOn, numberOfUniformLevels,
                          numberOfSmoothingSteps, finalSmoothIsOn, numberOfNonLinearSteps, numberOfIterations);
  }

  // define multilevel mesh
  MultiLevelMesh mlMsh;

  // Read coarse level mesh and generate finer level meshes.
  double scalingFactor = 1.;
  unsigned numberOfSelectiveLevels = 0;
  if(parameter.simulation == 0) {
    mlMsh.ReadCoarseMesh("../input/square.neu", "seventh", scalingFactor);
    //mlMsh.ReadCoarseMesh("../input/square1.neu", "seventh", scalingFactor);
  }
  else if(parameter.simulation == 1) {
    mlMsh.ReadCoarseMesh("../input/squareTri3D.neu", "seventh", scalingFactor);
  }
  else if(parameter.simulation == 2 || parameter.simulation == 3) {
    mlMsh.ReadCoarseMesh("../input/cylinder2.neu", "seventh", scalingFactor);
  }
  else if(parameter.simulation == 4) {
    mlMsh.ReadCoarseMesh("../input/intersection.neu", "seventh", scalingFactor);
  }
  else if(parameter.simulation == 5) {
    mlMsh.ReadCoarseMesh("../input/cat.med", "seventh", scalingFactor, false, false);
  }
  else if(parameter.simulation == 6) {
    mlMsh.ReadCoarseMesh("../input/handbndry.med", "seventh", scalingFactor, false, false);
  }
  else if(parameter.simulation == 7) {
    mlMsh.ReadCoarseMesh("../../Willmore/WillmoreSurface/input/newtorus.med", "seventh", scalingFactor, false, false);
  }
  else { //generic pick your mesh
    mlMsh.ReadCoarseMesh("../input/square.neu", "seventh", scalingFactor);
  }

  mlMsh.RefineMesh(parameter.numberOfUniformLevels , parameter.numberOfUniformLevels, NULL);

  // Erase all the coarse mesh levels.
  mlMsh.EraseCoarseLevels(parameter.numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  const unsigned  DIM = mlMsh.GetDimension();

  // Define the multilevel solution and attach the mlMsh object to it.
  MultiLevelSolution mlSol(&mlMsh);

  // Add variables X,Y,W to mlSol.

  FEOrder feOrder = FIRST;
  mlSol.AddSolution("Dx1", LAGRANGE, feOrder, 2);
  mlSol.AddSolution("Dx2", LAGRANGE, feOrder, 2);
  if(parameter.surface) mlSol.AddSolution("Dx3", LAGRANGE, feOrder, 2);

  if(parameter.constraintIsOn) mlSol.AddSolution("Lambda", LAGRANGE, feOrder, 0);

  if(areaConstraint) mlSol.AddSolution("Lambda1", DISCONTINUOUS_POLYNOMIAL, ZERO, 0);

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

  if(parameter.simulation < 2) {
    mlSol.AttachSetBoundaryConditionFunction(SetBoundaryConditionSquare);
  }
  else if(parameter.simulation < 4) {
    mlSol.AttachSetBoundaryConditionFunction(SetBoundaryConditionCylinder);
  }
  else if(parameter.simulation < 5) {
    mlSol.AttachSetBoundaryConditionFunction(SetBoundaryConditionIntersection);
  }
  else if(parameter.simulation < 8) {
    mlSol.AttachSetBoundaryConditionFunction(SetBoundaryConditionZero);
  }
  else { // generic pick your boundary
    mlSol.AttachSetBoundaryConditionFunction(SetBoundaryConditionSquare);
  }
  mlSol.GenerateBdc("All");

  mlSol.GenerateBdc("Dx1", "Time_dependent");
  mlSol.GenerateBdc("Dx2", "Time_dependent");
  if(parameter.surface) mlSol.GenerateBdc("Dx3", "Time_dependent");

  GetElementNearVertexNumber(mlSol);
  MultiLevelProblem mlProb(&mlSol);

  // Add system Conformal or Shear Minimization in mlProb.
  TransientNonlinearImplicitSystem& system = mlProb.add_system < TransientNonlinearImplicitSystem > ("conformal"); //for conformal

  // Add solutions newDX, Lambda to system.
  system.AddSolutionToSystemPDE("Dx1");
  system.AddSolutionToSystemPDE("Dx2");
  if(parameter.surface) system.AddSolutionToSystemPDE("Dx3");
  if(parameter.constraintIsOn) system.AddSolutionToSystemPDE("Lambda");
  if(areaConstraint) {
    system.AddSolutionToSystemPDE("Lambda1");
    system.SetNumberOfGlobalVariables(1);
  }

  // Parameters for convergence and # of iterations.
  system.SetMaxNumberOfNonLinearIterations(parameter.numberOfNonLinearSteps);
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
  if(parameter.surface) mov_vars.push_back("Dx3");
  mlSol.GetWriter()->SetMovingMesh(mov_vars);

  system.CopySolutionToOldSolution();
  for(unsigned k = 0; k < parameter.numberOfIterations; k++) {
    system.MGsolve();
    //ProjectSolution(mlSol);
    mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, k + 1);
    system.CopySolutionToOldSolution();
  }

  EvaluateMu(mlSol);
  mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, parameter.numberOfIterations);
  parameter.print();

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
      value = time / parameter.numberOfIterations * 0.75 * sin(x[1] / 0.5 * M_PI);
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
      value = time / parameter.numberOfIterations * 0.8 * sin(x[1] / 0.5 * M_PI);
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
        value = time / parameter.numberOfIterations * 1;
      else {
        value = - time / parameter.numberOfIterations * 1;
      }
    }
    if(!strcmp(solName, "Dx2")) {
      value = time / parameter.numberOfIterations * 0.35 * x[1] / 0.5;
    }
    else if(!strcmp(solName, "Dx3")) {
      value = time / parameter.numberOfIterations * 0.75 * x[2] / 0.5;
    }
  }
  else if(3 == faceName || 4 == faceName) {
    if(!strcmp(solName, "Dx1")) {
      value = time / parameter.numberOfIterations * 0.45 * x[0] / 0.5;
    }
    else if(!strcmp(solName, "Dx3")) {
      value = time / parameter.numberOfIterations * 0.35 * x[2] / 0.5;
    }
  }
  else if(5 == faceName || 6 == faceName) {
    if(!strcmp(solName, "Dx1")) {
      value = time / parameter.numberOfIterations * 0.125 * x[0] / 0.5;
    }
    else if(!strcmp(solName, "Dx2")) {
      value = time / parameter.numberOfIterations * 0.65 * x[1] / 0.5;
    }
  }

  if(!strcmp(solName, "Lambda")) {
    dirichlet = false;
  }
  return dirichlet;
}
