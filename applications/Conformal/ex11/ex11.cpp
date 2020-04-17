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

bool stopIterate = false;
const bool read_groups = false;                       //by default, if no argument is given, this is "true"
const bool read_boundary_groups = false;              //by default, if no argument is given, this is "true"
//unsigned muSmoothingType = 1; // the mesh should be logically structured and uniformly oriented
unsigned muSmoothingType = 2; // invariant with respect to quad orientation


unsigned conformalTriangleType = 2;

// const double normalSign = -1.;
const bool O2conformal = true;
const bool noLM = true;
unsigned counter = 0;
const double eps = 1.0e-3;// * O2conformal;

const unsigned numberOfIterations = 1;

using namespace femus;

#include "../include/supportFunctions.hpp"
#include "../include/updateMu1.hpp"
#include "../include/assembleConformalMinimization.hpp"

// Comment back in for working code
//const double mu[2] = {0.8, 0.};

void ProjectSolution(MultiLevelSolution& mlSol);

double InitalValueCM(const std::vector < double >& x) {
//   return cos(4.* M_PI * sqrt(x[0] * x[0] + x[1] * x[1])/0.5) ;
  return cos(20 * M_PI * x[0]) + sin(20 * M_PI * x[1]);
   //return sin(28.* M_PI * x[0]) * sin(28.* M_PI * x[1]) + cos(28.* M_PI * x[0]) * cos(28.* M_PI * x[1]) ;
}


double GetTimeStep (const double t) {
  return 1;
}

// IBVs.  No boundary, and IVs set to sphere (just need something).
bool SetBoundaryCondition(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time) {

    
  bool dirichlet = true;
  value = 0.;
  
  if(1 == faceName || 2 == faceName ) {
//     if(!strcmp(solName, "Dx1")) {
//       if(1 == faceName)
//        value = 1;
//       else{
//        value = -1; 
//       }
//     }  
    if(!strcmp(solName, "Dx2")) {
      value = 0.5 * x[1]/0.5;
    }
    else if(!strcmp(solName, "Dx3")) {
      value = 0.5 * x[2]/0.5;
    }
  }
  else if(3 == faceName || 4 == faceName ) {
    if(!strcmp(solName, "Dx1")) {
      value = 0.5 * x[0]/0.5;
    }
    else if(!strcmp(solName, "Dx3")) {
      value = 0.5 * x[2]/0.5;
    }
  }
  else if(5 == faceName || 6 == faceName ) {
    if(!strcmp(solName, "Dx1")) {
      value = 0.5 * x[0]/0.5;
    }
    else if(!strcmp(solName, "Dx2")) {
      value = 0.5 * x[1]/0.5;
    }
  }
  

  
  
  
  

//   bool dirichlet = false;
//   value = 0.;
// 
//   if(100 == faceName) {
//     //value = 0.6 * sin(x[1] / 0.5 * M_PI);
//     dirichlet = true;
//   }

//   if(!strcmp(solName, "Dx1")) {
//     if(3 == faceName || 3 == faceName) {
//       dirichlet = false;
//     }
//     if(4 == faceName) {
//       value = 0.75 * sin(x[1] / 0.5 * M_PI);
//     }
//   }
//   else if(!strcmp(solName, "Dx2")) {
//     if(2 == faceName) {
//       dirichlet = false;
//     }
//   }


//   if (!strcmp (solName, "Dx1")) {
//     if (1 == faceName || 3 == faceName) {
//       dirichlet = false;
//     }
//   }
//   else if (!strcmp (solName, "Dx2")) {
//     if (2 == faceName || 4 == faceName) {
//       dirichlet = false;
//     }
//   }



/*

  bool dirichlet = true;
  value = 0.;
  
  if(!strcmp(solName, "Dx1")) {
    if(1 == faceName) {
      //value = 0.04 * sin (4*(x[1] / 0.5 * acos (-1.)));
      value = time / numberOfIterations * 0.6 * sin(x[1] / 0.5 * M_PI);
      //dirichlet = false;
    }
  }*/

//   else if(!strcmp(solName, "Dx2")) {
//     if(1 == faceName) {
// 
//       //value = 0.04 * sin (4*(x[1] / 0.5 * acos (-1.)));
//       value = 0.25 * x[1];
//       //dirichlet = false;
//     }
//   }
// 
//   else if(!strcmp(solName, "Dx3")) {
//     if(1 == faceName) {
// 
//       //value = 0.04 * sin (4*(x[1] / 0.5 * acos (-1.)));
//       value = 0.25 * x[2];
//       //dirichlet = false;
//     }
//   }


  return dirichlet;
}


// Main program starts here.
int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);


  // define multilevel mesh
  unsigned maxNumberOfMeshes;
  MultiLevelMesh mlMsh;

  // Read coarse level mesh and generate finer level meshes.
  double scalingFactor = 1.;

  //mlMsh.GenerateCoarseBoxMesh(32, 32, 0, -0.5, 0.5, -0.5, 0.5, 0., 0., QUAD9, "seventh");

  //mlMsh.ReadCoarseMesh("../input/squareReg3D.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("../input/square13D.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("../input/squareTri3D.neu", "seventh", scalingFactor);
//   mlMsh.ReadCoarseMesh("../input/cylinder2.neu", "seventh", scalingFactor);
  mlMsh.ReadCoarseMesh("../input/intersection.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("../input/hand.med", "seventh", scalingFactor);

  //mlMsh.ReadCoarseMesh("../input/cat.med", "seventh", scalingFactor, read_groups, read_boundary_groups);


  unsigned numberOfUniformLevels = 2;
  unsigned numberOfSelectiveLevels = 0;
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
  mlSol.AddSolution("Dx3", LAGRANGE, feOrder, 2);

  mlSol.AddSolution("Lambda1", DISCONTINUOUS_POLYNOMIAL, ZERO, 0);

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

  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("All");
  
  mlSol.GenerateBdc("Dx1","Time_dependent");
  mlSol.GenerateBdc("Dx2","Time_dependent");
  mlSol.GenerateBdc("Dx3","Time_dependent");
 

  GetElementNearVertexNumber(mlSol);

  MultiLevelProblem mlProb(&mlSol);

  // Add system Conformal or Shear Minimization in mlProb.
  TransientNonlinearImplicitSystem& system = mlProb.add_system < TransientNonlinearImplicitSystem > ("conformal"); //for conformal

  // Add solutions newDX, Lambda1 to system.
  system.AddSolutionToSystemPDE("Dx1");
  system.AddSolutionToSystemPDE("Dx2");
  system.AddSolutionToSystemPDE("Dx3");
  system.AddSolutionToSystemPDE("Lambda1");

  // Parameters for convergence and # of iterations.
  system.SetMaxNumberOfNonLinearIterations(10);
  system.SetNonLinearConvergenceTolerance(1.e-10);

  
  system.AttachGetTimeIntervalFunction (GetTimeStep);
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
  
  for (unsigned k = 0; k < numberOfIterations; k++) {
    system.CopySolutionToOldSolution();
    system.MGsolve();
    //ProjectSolution(mlSol);
    mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, k + 1);
  }
  
  return 0;
}


void ProjectSolution(MultiLevelSolution& mlSol) {
  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1u;

  Solution* sol = mlSol.GetSolutionLevel(level);
  Mesh* msh = mlSol._mlMesh->GetLevel(level);
  elem* el = msh->el;

  //unsigned  dim = msh->GetDimension();
  unsigned dim = 2;
  unsigned DIM = 3;
  std::vector < unsigned > indexDx(DIM);
  indexDx[0] = mlSol.GetIndex("Dx1");
  indexDx[1] = mlSol.GetIndex("Dx2");
  indexDx[2] = mlSol.GetIndex("Dx3");
  unsigned solType = mlSol.GetSolutionType(indexDx[0]);

  unsigned iproc = msh->processor_id();

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned nDofs  = msh->GetElementDofNumber(iel, solType);

    for(unsigned i = 0; i < nDofs; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, solType);
      unsigned xdof = msh->GetSolutionDof(i, iel, 2);

      double Yhat = (*msh->_topology->_Sol[1])(xdof);
      double Zhat = (*msh->_topology->_Sol[2])(xdof);

      double DY  = (*sol->_Sol[indexDx[1]])(idof);
      double DZ  = (*sol->_Sol[indexDx[2]])(idof);

      double theta = atan2(Zhat + DZ, Yhat + DY);

      sol->_Sol[indexDx[1]]->set(idof, 0.5 * cos(theta) - Yhat);
      sol->_Sol[indexDx[2]]->set(idof, 0.5 * sin(theta) - Zhat);
    }
  }

  sol->_Sol[indexDx[1]]->close();
  sol->_Sol[indexDx[2]]->close();

}
