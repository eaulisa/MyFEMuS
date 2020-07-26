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

unsigned conformalType0 = 2;
unsigned conformalType = 2;

using namespace femus;

#include "../include/parameter.hpp"

Parameter parameter;

#include "../include/supportFunctions.hpp"
#include "../include/updateMu8.hpp"
#include "../include/assembleConformalMinimization9.hpp"

double InitalValueCM(const std::vector < double >& x) {
//   return cos(4.* M_PI * sqrt(x[0] * x[0] + x[1] * x[1])/0.5) ;
  return cos(20 * M_PI * x[0]) + sin(20 * M_PI * x[1]);
  //return sin(28.* M_PI * x[0]) * sin(28.* M_PI * x[1]) + cos(28.* M_PI * x[0]) * cos(28.* M_PI * x[1]) ;
}
double GetTimeStep(const double t) {
  return 1;
}


void UpdateMesh(MultiLevelSolution &mlSol);

bool SetBoundaryConditionZero(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time);
bool SetBoundaryConditionSquare(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time);
bool SetBoundaryConditionCylinder(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time);
bool SetBoundaryConditionIntersection(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time);

const Parameter squareQuad = Parameter("square with quads", 0, false, false, 4, 1, true, 300, 1, 0.811569);
const Parameter squareTri = Parameter("square with triangles", 1, false, false, 1, 1, true, 500, 1, 0.805200);
//const Parameter cylinderUnconstrained = Parameter("cylinder unconstrained", 2, true, false, 4, 12, false, 30, 1, 0.910958);
const Parameter cylinderUnconstrained = Parameter("cylinder unconstrained", 2, true, false, 4, 3, true, 250, 1, 0.746343);
//const Parameter cylinderConstrained = Parameter("cylinder constrained", 3, true, false, 4, 3, true, 100, 1, 0.730090); //areaConstraint
const Parameter cylinderConstrained = Parameter("cylinder constrained", 3, true, true, 4, 2, true, 1000, 1, 0.793786); //normal constraint
const Parameter intersection = Parameter("intersection", 4, true, false, 2, 100, true, 10, 5, 0.486729);
//const Parameter intersection = Parameter("intersection", 4, true, false, 2, 100, true, 50, 1, 0.674721);
//const Parameter intersection = Parameter("intersection", 4, true, false, 2, 12, false, 30, 1, 0.979639);
//const Parameter cat = Parameter("cat", 5, true, false, 1, 20, false, 5, 1, 0.987606); //need conformal type 0,0 // areaConstraint // no folding!?!?!
//const Parameter cat = Parameter("cat", 5, true, true, 1, 50, false, 5, 1, 0.986969); //need conformal type 0,0 // no folding?!?!?!
const Parameter cat = Parameter("cat", 5, true, true, 1, 1, true, 5, 1, 0.996086); //need conformal type 0,0
//const Parameter cat = Parameter("cat", 5, true, false, 1, 20, true, 5, 1, 0.996710); //need conformal type 0,0 // areaConstraint
//const Parameter cat = Parameter("cat", 5, true, true, 1, 25, false, 3, 1, 0.986754); //need conformal type 0,0
const Parameter hand = Parameter("hand", 6, true, true, 1, 12, false, 10, 1, 0.580335);
//const Parameter moo = Parameter("moo", 7, true, true, 1, 1, true, 40, 1, 0.654910);
const Parameter moo = Parameter("moo", 7, true, true, 2, 1, true, 50, 1, 0.602613);
const Parameter moai = Parameter("moai", 8, true, true, 1, 20, true, 20, 1, 0.888489);
const Parameter fert = Parameter("fert", 9, true, true, 1, 20, true, 3, 1, 0.995966);

// Main program starts here.
int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

//   std::vector <double> angle(3);
//   std::vector < std::vector <double> > xC(2);
// //   xC[0].resize(4);
// //   xC[1].resize(4);
// //
// //   angle[0] = angle[1] = angle[2] = 2. * M_PI / 4.;
// //   angle[3] = 2. * M_PI / 3.;
//
//   xC[0].resize(3);
//   xC[1].resize(3);
//
//   angle[2] = angle[1] = M_PI / 3.;
//   angle[0] = M_PI / 2.;
//
//
//   GetConformalStructure(angle, xC);
//
//   for(unsigned i = 0; i <= xC[0].size(); i++) {
//     unsigned ii = i % xC[0].size();
//     std::cout << xC[0][ii] << " " << xC[1][ii]   << std::endl;
//   }
//  return 1;





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
    else if(!strcmp("8", args[1])) { // moai
      parameter = moai;
    }
    else if(!strcmp("9", args[1])) { // moai
      parameter = fert;
    }
    else {
      goto generic;
    }
  }
  else { //generic
  generic:
    unsigned simulation = 8;
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
    //mlMsh.ReadCoarseMesh("../input/squareReg.neu", "seventh", scalingFactor);
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
    mlMsh.ReadCoarseMesh("../../Willmore/WillmoreSurface/input/moo.med", "seventh", scalingFactor, false, false);
  }
  else if(parameter.simulation == 8) {
    mlMsh.ReadCoarseMesh("../../Willmore/WillmoreSurface/input/DTquad.med", "seventh", scalingFactor, false, false);
  }
  else if(parameter.simulation == 9) {
    mlMsh.ReadCoarseMesh("../../Willmore/WillmoreSurface/input/stupid.med", "seventh", scalingFactor, false, false);
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

  if(areaConstraint && parameter.simulation > 3) {
    mlSol.FixSolutionAtOnePoint("Dx1");
    mlSol.FixSolutionAtOnePoint("Dx2");
    if(parameter.surface) mlSol.FixSolutionAtOnePoint("Dx3");
  }



  if(parameter.constraintIsOn) mlSol.AddSolution("Lambda", LAGRANGE, feOrder, 0);

  if(areaConstraint) mlSol.AddSolution("Lambda1", DISCONTINUOUS_POLYNOMIAL, ZERO, 0);

  mlSol.AddSolution("env", LAGRANGE, FIRST, 0, false);
  mlSol.AddSolution("vAngle", LAGRANGE, FIRST, 0);

  mlSol.AddSolution("mu1", DISCONTINUOUS_POLYNOMIAL, ZERO, 2);
  mlSol.AddSolution("mu2", DISCONTINUOUS_POLYNOMIAL, ZERO, 2);
  mlSol.AddSolution("lmu1", DISCONTINUOUS_POLYNOMIAL, ZERO, 0);
  mlSol.AddSolution("lmu2", DISCONTINUOUS_POLYNOMIAL, ZERO, 0);
//   mlSol.FixSolutionAtOnePoint("lmu1");
//   mlSol.FixSolutionAtOnePoint("lmu2");
//   mlSol.FixSolutionAtOnePoint("mu1");
//   mlSol.FixSolutionAtOnePoint("mu2");


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
  else if(parameter.simulation < 10) {
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
  system.SetNonLinearConvergenceTolerance(1.e-10);


  system.AttachGetTimeIntervalFunction(GetTimeStep);
  system.SetAssembleFunction(AssembleConformalMinimization);

  // system.init();
  // system.SetTime(-1.);
  // system.SetMaxNumberOfNonLinearIterations(1);
  // system.MGsolve();
  // UpdateMesh(mlSol);
  //
  // system.SetTime(1.);
  // mlSol.GenerateBdc("Dx1", "Time_dependent");
  // mlSol.GenerateBdc("Dx2", "Time_dependent");
  // if(parameter.surface) mlSol.GenerateBdc("Dx3", "Time_dependent");
  system.SetTime(0.);
  system.SetMaxNumberOfNonLinearIterations(parameter.numberOfNonLinearSteps);
  system.init();


  // Add system Conformal or Shear Minimization in mlProb.
  LinearImplicitSystem& systemMu = mlProb.add_system < LinearImplicitSystem > ("mu"); //for conformal

  // Add solutions newDX, Lambda to system.
  systemMu.AddSolutionToSystemPDE("mu1");
  systemMu.AddSolutionToSystemPDE("mu2");
  systemMu.AddSolutionToSystemPDE("lmu1");
  systemMu.AddSolutionToSystemPDE("lmu2");
  systemMu.SetAssembleFunction(AssembleResMu);
  systemMu.init();


  counter = 0;

  mlSol.SetWriter(VTK);

  std::vector<std::string> mov_vars;
  mov_vars.push_back("Dx1");
  mov_vars.push_back("Dx2");
  if(parameter.surface) mov_vars.push_back("Dx3");
  mlSol.GetWriter()->SetMovingMesh(mov_vars);

  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");
  mlSol.GetWriter()->SetDebugOutput(true);

  //EvaluateMu(mlSol);
  mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);
  parameter.print();

  system.CopySolutionToOldSolution();
  for(unsigned k = 0; k < parameter.numberOfIterations; k++) {
//     if(k == parameter.numberOfIterations - 1)  {
//       system.SetMaxNumberOfNonLinearIterations(1000);
//     }
    system.MGsolve();
    //ProjectSolution(mlSol);
    mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, k + 1);
    system.CopySolutionToOldSolution();
  }

  EvaluateMu(mlSol);
  mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, parameter.numberOfIterations);
  parameter.print();


  delete  PtP[0][0];
  delete  PtP[0][1];
  delete  PtP[1][0];
  delete  PtP[1][1];

  return 0;
}

bool SetBoundaryConditionZero(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time) {

  bool dirichlet = true;
  value = 0.;

  if(!strcmp(solName, "Lambda")) {
    dirichlet = false;
  }
  if(!strcmp(solName, "vAngles")) {
    value = M_PI;
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
    if(2 == faceName || (time == 0. && 4 == faceName)) {
      dirichlet = false;
    }
  }

  else if(!strcmp(solName, "Lambda")) {
    dirichlet = false;
  }

  else if(!strcmp(solName, "vAngle")) {
    value = M_PI;
    if(fabs(x[0]) > 0.49999 &&  fabs(x[1]) > 0.49999) {
      value = M_PI * 1.5;
    }
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

  else if(!strcmp(solName, "vAngle")) {
    value = M_PI;
    if(fabs(x[0]) > 0.49999 &&  fabs(x[1]) > 0.49999) {
      value = M_PI;
    }
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


void UpdateMesh(MultiLevelSolution &mlSol) {

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution* sol = mlSol.GetSolutionLevel(level);
  Mesh* msh = mlSol._mlMesh->GetLevel(level);

  NumericVector* mysol;

  mysol = NumericVector::build().release();
  mysol->init(*msh->_topology->_Sol[0]);

  unsigned dim = (parameter.surface) ?  3 : 2;

  std::string solName[3] = {"Dx1", "Dx2", "Dx3"};

  std::vector < unsigned > solDxIndex(dim);
  for(int k = 0; k < dim; k++) {
    solDxIndex[k] = mlSol.GetIndex(solName[k].c_str());
  }

  unsigned solType = mlSol.GetSolutionType(solDxIndex[0]);

  for(int k = 0; k < dim; k++) {
    mysol->matrix_mult(*sol->_Sol[solDxIndex[k]], *msh->GetQitoQjProjection(2, solType));
    *msh->_topology->_Sol[k] += *mysol;
    sol->_Sol[solDxIndex[k]]->zero();
    sol->_SolOld[solDxIndex[k]]->zero();
  }

  delete mysol;

}
