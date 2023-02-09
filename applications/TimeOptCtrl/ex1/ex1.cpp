/** \file Ex13.cpp
 *  \brief This example shows how to set and solve the weak form
 *   of the time dependent Stress-Strain Equation, for a beam with pressure
 *
 *  \nabla \cdot \sigma + \mass*acceleration = F
 *  in a Beam domain (in 2D and 3D) clamped on one end
 *
 *  \author Eugenio Aulisa
 */


#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "LinearImplicitSystem.hpp"
#include "TransientSystem.hpp"
#include "adept.h"

using namespace femus;

const unsigned alpha = 1.0e-5;
const double beta = 1.0e-5;
const double Pe = 10.;
double dt = 0.1;

unsigned iext;

double SetVariableTimeStep(const double time) {
  return dt;
}

void AssembleTimeOptimalControl(MultiLevelProblem& ml_prob);

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true;
  value = 0.;

  std::string name = SolName;
  std::string subName = name.substr(0, 1);

  if(facename == 4 &&  facename == 2) {
    if(subName == "T" || subName == "z") {
      value = 1.;
    }
  }
  else if(facename == 3) {
    if(subName == "T" || subName == "z") {
      value = 2.;
    }
  }
  else if(facename == 1) {
    if(subName == "T" || subName == "z") {
      dirichlet = false;
    }
  }
  return dirichlet;
}



int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;

  unsigned nx = 10;
  unsigned ny = 10;
  unsigned nz = 1;

  double length = .1;
  double lengthx = 1;

  mlMsh.GenerateCoarseBoxMesh(nx, ny, 0, 0., lengthx, 0., length, 0., 0., QUAD9, "seventh");

  unsigned dim = mlMsh.GetDimension();

  unsigned numberOfUniformLevels = 2;
  unsigned numberOfSelectiveLevels = 0;

  mlMsh.RefineMesh(numberOfUniformLevels, numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // erase all the coarse mesh levels
  mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  MultiLevelSolution mlSol(&mlMsh);


  unsigned numberOfIterations = 3;
  // add variables to mlSol
  for(unsigned i = 0; i < numberOfIterations; i++) {
    char TName[10];
    char lName[10];
    char zName[10];
    sprintf(TName, "T%d", i);
    sprintf(lName, "l%d", i);
    sprintf(zName, "z%d", i);
    mlSol.AddSolution(TName, LAGRANGE, SECOND, 2);
    mlSol.AddSolution(lName, LAGRANGE, SECOND, 2);
    mlSol.AddSolution(zName, LAGRANGE, SECOND, 2);
  }

  mlSol.Initialize("All");

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("All");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);


  std::vector<TransientLinearImplicitSystem*> system(numberOfIterations);
  for(unsigned i = 0; i < numberOfIterations; i++) {

    char sName[10];
    char TName[10];
    char lName[10];
    char zName[10];
    sprintf(sName, "top%d", i);
    sprintf(TName, "T%d", i);
    sprintf(lName, "l%d", i);
    sprintf(zName, "z%d", i);

    system[i] = &(mlProb.add_system < TransientLinearImplicitSystem > (sName));
    //add solution "D" to system
    system[i]->AddSolutionToSystemPDE(lName);
    system[i]->AddSolutionToSystemPDE(TName);
    system[i]->AddSolutionToSystemPDE(zName);

    // attach the assembling function to system
    system[i]->AttachGetTimeIntervalFunction(SetVariableTimeStep);

    system[i]->SetAssembleFunction(AssembleTimeOptimalControl);

    // initilaize and solve the system
    system[i]->init();
    system[i]->SetOuterSolver(PREONLY);
  }
  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(true);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);

  for(unsigned t = 0; t < 100; t++) {
    for(iext = 0; iext < numberOfIterations; iext++) {
      system[iext]->CopySolutionToOldSolution(); // Copy D, V, and A into DOld, VOld, and AOld, respectively
      system[iext]->MGsolve(); //solve for A, using DOld, VOld, and AOld
    }

    vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, t + 1);
  }
  return 0;
}



void AssembleTimeOptimalControl(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled


  char sName[10];
  char TName[10];
  char lName[10];
  char zName[10];
  sprintf(sName, "top%d", iext);
  sprintf(TName, "T%d", iext);
  sprintf(lName, "l%d", iext);
  sprintf(zName, "z%d", iext);


  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  TransientLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<TransientLinearImplicitSystem> (sName);   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual std::vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable


  unsigned TiIndex = mlSol->GetIndex(TName);
  unsigned liIndex = mlSol->GetIndex(lName);
  unsigned ziIndex = mlSol->GetIndex(zName);

  unsigned solType = mlSol->GetSolutionType(TiIndex);

  unsigned TiPdeIndex = mlPdeSys->GetSolPdeIndex(TName);
  unsigned liPdeIndex = mlPdeSys->GetSolPdeIndex(lName);
  unsigned ziPdeIndex = mlPdeSys->GetSolPdeIndex(zName);

  std::vector < adept::adouble > solTi;
  std::vector < adept::adouble > solli;
  std::vector < adept::adouble > solzi;

  std::vector < adept::adouble > mResTi;
  std::vector < adept::adouble > mResli;
  std::vector < adept::adouble > mReszi;

  std::vector < std::vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector <double> phi;  // local test function for velocity
  std::vector <double> gradPhi; // local test function first order partial derivatives
  double weight;

  std::vector< unsigned > sysDof; // local to global pdeSys dofs
  std::vector< double > Res; // local redidual std::vector
  std::vector < double > Jac;

  RES->zero(); // Set to zero all the entries of the Global Residual std::vector
  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);

    unsigned nDofs = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
    unsigned nDofsAll = 3 * nDofs;
    sysDof.resize(nDofsAll);

    for(unsigned  k = 0; k < dim; k++) {
      x[k].resize(nDofs);
    }
    solTi.resize(nDofs);
    solli.resize(nDofs);
    solzi.resize(nDofs);

    mResTi.assign(nDofs, 0.);
    mResli.assign(nDofs, 0.);
    mReszi.assign(nDofs, 0.);

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, solType);

      solTi[i] = (*sol->_Sol[TiIndex])(solDof);
      solli[i] = (*sol->_Sol[liIndex])(solDof);
      solzi[i] = (*sol->_Sol[ziIndex])(solDof);

      sysDof[i] = pdeSys->GetSystemDof(TiIndex, TiPdeIndex, i, iel);
      sysDof[nDofs + i] = pdeSys->GetSystemDof(liIndex, liPdeIndex, i, iel);
      sysDof[2 * nDofs + i] = pdeSys->GetSystemDof(ziIndex, ziPdeIndex, i, iel);
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);
      for(unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k])(xDof);
      }
    }

    std::vector<bool> nodeIsControlBoundary(nDofs, false);

    for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
      unsigned int facename = -(msh->el->GetFaceElementIndex(iel, jface) + 1);
      if(facename == 1) {
        unsigned nve = msh->GetElementFaceDofNumber(iel, jface, solType);
        const unsigned felt = msh->GetElementFaceType(iel, jface);
        for(unsigned i = 0; i < nve; i++) {
          nodeIsControlBoundary[ msh->GetLocalFaceVertexIndex(iel, jface, i)] = true;
        }
      }
    }


    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {
      msh->_finiteElement[ielGeom][solType]->Jacobian(x, ig, weight, phi, gradPhi);


      adept::adouble Tig = 0;
      adept::adouble lig = 0;
      adept::adouble zig = 0;
      std::vector < adept::adouble > gradTig(dim, 0);
      std::vector < adept::adouble > gradlig(dim, 0);
      std::vector < adept::adouble > gradzig(dim, 0);


      for(unsigned i = 0; i < nDofs; i++) {
        Tig += solTi[i] * phi[i];
        lig += solli[i] * phi[i];
        zig += solzi[i] * phi[i];
        for(unsigned  j = 0; j < dim; j++) {
          gradTig[j] += solTi[i] * gradPhi[i * dim + j];
          gradlig[j] += solli[i] * gradPhi[i * dim + j];
          gradzig[j] += solzi[i] * gradPhi[i * dim + j];
        }
      }

      double nu = 1. / Pe;

      // *** phiV_i loop ***
      for(unsigned i = 0; i < nDofs; i++) {

        adept::adouble Tieq = 0;
        adept::adouble lieq = 0;
        adept::adouble zieq = 0;

        std::vector < adept::adouble > NSV(dim, 0.);

        for(unsigned j = 0; j < dim; j++) {  // second index j in each equation
          Tieq   +=  gradPhi[i * dim + j] * (nu * gradlig[j] + alpha * gradTig[j]); // laplace
          lieq   +=  gradPhi[i * dim + j] * (nu * gradTig[j]); // laplace
          zieq   +=  gradPhi[i * dim + j] * (nu * gradzig[j]); // laplace
        }
        Tieq += beta * phi[i] * Tig;

        mResTi[i] += Tieq * weight;
        mResli[i] += lieq * weight;
        if(!nodeIsControlBoundary[i]) mReszi[i] += zieq * weight;
        else mReszi[i] += solzi[i] - solTi[i];

      } // end phiV_i loop
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube mRes in double Res and store them in RES
    Res.resize(nDofsAll);    //resize

    for(int i = 0; i < nDofs; i++) {
      Res[ i ] = -mResTi[i].value();
      Res[ nDofs + i ] = -mResli[i].value();
      Res[ 2 * nDofs + i ] = -mReszi[i].value();
    }

    RES->add_vector_blocked(Res, sysDof);

    //Extract and store the Jacobian
    Jac.resize(nDofsAll * nDofsAll);

    s.dependent(&mResTi[0], nDofs);
    s.dependent(&mResli[0], nDofs);
    s.dependent(&mReszi[0], nDofs);

    s.independent(&solTi[0], nDofs);
    s.independent(&solli[0], nDofs);
    s.independent(&solzi[0], nDofs);

    s.jacobian(&Jac[0], true);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();
  KK->close();
// ***************** END ASSEMBLY *******************

}



