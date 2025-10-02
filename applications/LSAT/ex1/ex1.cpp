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


const unsigned DIM = 2;
const double BETA = 0.25;
const double GAMMA = 0.5;
bool UseNewmarkUpdateWithD = true;
double dt = 0.025;

double SetVariableTimeStep(const double time) {
  return dt;
}


using namespace std;
using namespace femus;

double SetRegions(Solution *sol);
void SetTarget(Solution *sol, const std::string &R, const double &time);

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true;
  value = 0.;

  // if((DIM == 2 && facename == 4) || (DIM == 3 && facename == 5)) {   // left boundary condition.
  //   dirichlet = true;
  // }
  //
  // if((DIM == 2 && facename == 3) || (DIM == 3 && facename == 4)) {   // top boundary condition.
  //   dirichlet = false;
  //   value = 500.;
  // }

  return dirichlet;
}


void NewmarkUpdate(MultiLevelSolution *mlSol);
void NewmarkUpdateWithD(MultiLevelSolution *mlSol);

void AssembleResAD(MultiLevelProblem& ml_prob);

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;

  unsigned nx = 10;
  unsigned ny = 2;
  unsigned nz = 1;

  double length = 1;
  double lengthx = M_PI;

  if(DIM == 2) {
    mlMsh.GenerateCoarseBoxMesh(nx, ny, 0, 0., lengthx, 0., length, 0., 0., QUAD9, "seventh");
  }
  else if(DIM == 3) {
    nz = ny;
    mlMsh.GenerateCoarseBoxMesh(nx, ny, nz, 0., lengthx, 0., length, 0., length,  HEX27, "seventh");
  }

  unsigned dim = mlMsh.GetDimension();

  unsigned numberOfUniformLevels = 4;
  unsigned numberOfSelectiveLevels = 0;

  mlMsh.RefineMesh(numberOfUniformLevels, numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // erase all the coarse mesh levels
  mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("Zi", LAGRANGE, SECOND, 2); //2 means that we have 2 solutions sol and solOld
  mlSol.AddSolution("Xi", LAGRANGE, SECOND);
  mlSol.AddSolution("Yi", LAGRANGE, SECOND);
  mlSol.AddSolution("Ei", LAGRANGE, SECOND, false);

  mlSol.AddSolution("Z", LAGRANGE, SECOND, false);

  mlSol.AddSolution("Z0", LAGRANGE, SECOND, false); //2 means that we have 2 solutions sol and solOld
  mlSol.AddSolution("Z0Old", LAGRANGE, SECOND, false);
  mlSol.AddSolution("X0", LAGRANGE, SECOND, false);
  mlSol.AddSolution("Y0", LAGRANGE, SECOND, false);

  mlSol.AddSolution("Z1", LAGRANGE, SECOND, false);
  mlSol.AddSolution("Z1Old", LAGRANGE, SECOND, false);//2 means that we have 2 solutions sol and solOld
  mlSol.AddSolution("X1", LAGRANGE, SECOND, false);
  mlSol.AddSolution("Y1", LAGRANGE, SECOND, false);

  mlSol.AddSolution("R", LAGRANGE, SECOND, false);
  mlSol.AddSolution("E0", LAGRANGE, SECOND, false);
  mlSol.AddSolution("E1", LAGRANGE, SECOND, false);

  mlSol.AddSolution("B", DISCONTINUOUS_POLYNOMIAL, ZERO, false);
  mlSol.AddSolution("C", DISCONTINUOUS_POLYNOMIAL, ZERO, false);

  mlSol.Initialize("All");

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("All");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);

  // add system Poisson in mlProb as a Linear Implicit System
  TransientLinearImplicitSystem& system = mlProb.add_system < TransientLinearImplicitSystem > ("LSAT");

  //add solution "D" to system
  system.AddSolutionToSystemPDE("Zi");
  system.AddSolutionToSystemPDE("Xi");
  system.AddSolutionToSystemPDE("Yi");

  system.SetAssembleFunction(AssembleResAD);



  // attach the assembling function to system
  system.AttachGetTimeIntervalFunction(SetVariableTimeStep);

  // initilaize and solve the system
  system.init();
  system.SetOuterSolver(PREONLY);


  const unsigned level = mlMsh.GetNumberOfLevels() - 1;
  Solution* sol = mlSol.GetSolutionLevel(level);

  SetRegions(sol);

  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(false);


  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);

  const unsigned int n_timesteps = 150;

  sol->_Sol[(mlSol.GetIndex("Z0Old"))]->zero();
  sol->_Sol[(mlSol.GetIndex("Z1Old"))]->zero();

  for(unsigned t = 1; t <= n_timesteps; t++) {

    //system.CopySolutionToOldSolution(); // Copy D, V, and A into DOld, VOld, and AOld, respectively
    SetTarget(sol, "R", t * dt);

    sol->_Sol[mlSol.GetIndex("Z")]->zero();

    *(sol->_Sol[mlSol.GetIndex("Ei")]) = *(sol->_Sol[(mlSol.GetIndex("R"))]);
    *(sol->_Sol[mlSol.GetIndex("Zi")]) = *(sol->_Sol[(mlSol.GetIndex("Z0Old"))]);
    *(sol->_SolOld[mlSol.GetIndex("Zi")]) = *(sol->_Sol[(mlSol.GetIndex("Z0Old"))]);
    system.MGsolve(); //solve for A, using DOld, VOld, and AOld
    *(sol->_Sol[mlSol.GetIndex("Z0")]) = *(sol->_Sol[(mlSol.GetIndex("Zi"))]);
    *(sol->_Sol[mlSol.GetIndex("E0")]) = *(sol->_Sol[(mlSol.GetIndex("Ei"))]);
    *(sol->_Sol[mlSol.GetIndex("E0")]) -= *(sol->_Sol[(mlSol.GetIndex("Zi"))]);
    *(sol->_Sol[mlSol.GetIndex("Z")]) += *(sol->_Sol[(mlSol.GetIndex("Zi"))]);


    *(sol->_Sol[mlSol.GetIndex("Ei")]) = *(sol->_Sol[(mlSol.GetIndex("E0"))]);
    *(sol->_Sol[mlSol.GetIndex("Zi")]) = *(sol->_Sol[(mlSol.GetIndex("Z1Old"))]);
    *(sol->_SolOld[mlSol.GetIndex("Zi")]) = *(sol->_Sol[(mlSol.GetIndex("Z1Old"))]);
    system.MGsolve(); //solve for A, using DOld, VOld, and AOld
    *(sol->_Sol[mlSol.GetIndex("Z1")]) = *(sol->_Sol[(mlSol.GetIndex("Zi"))]);
    *(sol->_Sol[mlSol.GetIndex("E1")]) = *(sol->_Sol[(mlSol.GetIndex("Ei"))]);
    *(sol->_Sol[mlSol.GetIndex("E1")]) -= *(sol->_Sol[(mlSol.GetIndex("Zi"))]);
    *(sol->_Sol[mlSol.GetIndex("Z")]) += *(sol->_Sol[(mlSol.GetIndex("Zi"))]);

     vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, t);

     *(sol->_Sol[(mlSol.GetIndex("Z0Old"))]) = *(sol->_Sol[mlSol.GetIndex("Z0")]);
     *(sol->_Sol[(mlSol.GetIndex("Z1Old"))]) = *(sol->_Sol[mlSol.GetIndex("Z1")]);

  }



  return 0;
}


double flc4hs(double const & x, double const & eps) {

  double r = x / eps;
  if(r < -1) {
    return 0.;
  }
  else if(r < 1.) {
    double r2 = r * r;
    double r3 = r * r2;
    double r5 = r3 * r2;
    double r7 = r5 * r2;
    double r9 = r7 * r2;
    return (128. + 315. * r - 420. * r3 + 378. * r5 - 180. * r7 + 35. * r9) / 256.;
  }
  else {
    return 1.;
  }
}

double GetTargetSolution(const std::vector<double> &xv, const double &time) {
  return cos(xv[0] * xv[1]) * flc4hs(time-.5,0.5) ;
}

void SetTarget(Solution *sol, const std::string &R, const double &time) {

  Mesh* msh = sol->GetMesh();    // pointer to the mesh (level) object
  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  unsigned rIndex = sol->GetIndex(R.c_str());

  unsigned rType = sol->GetSolutionType(rIndex);

  std::vector < double > xv(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned nDofsR = msh->GetElementDofNumber(iel, rType);
    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsR; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // local to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        xv[k] = (*msh->_topology->_Sol[k])(xDof);
      }
      unsigned uDof = msh->GetSolutionDof(i, iel, rType);    // local to global mapping between solution node and solution dof
      //rotation;

      double r = GetTargetSolution(xv, time);
      sol->_Sol[rIndex]->set(uDof,r);

    }
  }
  sol->_Sol[rIndex]->close();
}


double CheckIfInsideB(const std::vector<double> &xv) {
  return (xv[0] > 0.5 && xv[0] < 1. && xv[1] > 0.25 && xv[1] < 0.75 );
}

double CheckIfInsideC(const std::vector<double> &xv) {
  return (xv[0] > 0.5 && xv[0] < 1. && xv[1] > 0.25 && xv[1] < 0.75 );
}

double SetRegions(Solution *sol) {

  Mesh* msh = sol->GetMesh();    // pointer to the mesh (level) object
  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  unsigned  IndexB = sol->GetIndex("B");
  unsigned  IndexC = sol->GetIndex("C");

  std::vector < double > xv(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned nDofs = msh->GetElementDofNumber(iel, 1);
    // local storage of global mapping and solution

    double ielIsB = 1.;
    double ielIsC = 1.;

    for(unsigned i = 0; i < nDofs; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // local to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        xv[k] = (*msh->_topology->_Sol[k])(xDof);
      }

      ielIsB = CheckIfInsideB(xv);
      ielIsC = CheckIfInsideC(xv);
    }
    sol->_Sol[IndexB]->set(iel, ielIsB);
    sol->_Sol[IndexC]->set(iel, ielIsC);
  }
  sol->_Sol[IndexB]->close();
  sol->_Sol[IndexC]->close();
}




//Assemble Residual using A to update D amd V
void AssembleResAD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled


  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  TransientLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<TransientLinearImplicitSystem> ("LSAT");   // pointer to the linear implicit system named "Beam"
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

  double dt =  mlPdeSys->GetIntervalTime();
  double time =  mlPdeSys->GetTime();


  //solution variable
  unsigned solIndexZ = mlSol->GetIndex("Zi");
  unsigned solIndexX = mlSol->GetIndex("Xi");
  unsigned solIndexY = mlSol->GetIndex("Yi");

  unsigned solIndexE = mlSol->GetIndex("Ei");

  unsigned solIndexB = mlSol->GetIndex("B");
  unsigned solIndexC = mlSol->GetIndex("C");



  unsigned solType = mlSol->GetSolutionType(solIndexZ);


  unsigned solPdeIndexZ = mlPdeSys->GetSolPdeIndex("Zi");
  unsigned solPdeIndexX = mlPdeSys->GetSolPdeIndex("Xi");
  unsigned solPdeIndexY = mlPdeSys->GetSolPdeIndex("Yi");


  std::vector < double > solZOld;    // local solution

  std::vector < adept::adouble > solZ;    // local solution
  std::vector < adept::adouble > solX;
  std::vector < adept::adouble > solY;

  std::vector < double > solE;

  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector <double> phi;  // local test function for velocity
  std::vector <double> gradPhi; // local test function first order partial derivatives
  double weight; // gauss point weight

  std::vector < unsigned > sysDof; // local to global pdeSys dofs
  std::vector < adept::adouble > aRes;
  std::vector < double > res; // local redidual std::vector
  std::vector < double > Jac;

  RES->zero(); // Set to zero all the entries of the Global Residual std::vector
  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);

    double BBs = (*sol->_Sol[solIndexB])(iel);;
    double CsC = (*sol->_Sol[solIndexC])(iel);;

    unsigned nDofs = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs

    unsigned nDofsAll = 3 * nDofs;

    sysDof.resize(nDofsAll);
    aRes.assign(nDofsAll, 0.);

    solZOld.resize(nDofs);
    solZ.resize(nDofs);
    solX.resize(nDofs);
    solY.resize(nDofs);
    solE.resize(nDofs);

    for(unsigned  k = 0; k < dim; k++) {
      coordX[k].resize(nDofs);
    }

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned iDof = msh->GetSolutionDof(i, iel, solType);

      solZOld[i] = (*sol->_SolOld[solIndexZ])(iDof);
      solZ[i] = (*sol->_Sol[solIndexZ])(iDof);
      solX[i] = (*sol->_Sol[solIndexX])(iDof);
      solY[i] = (*sol->_Sol[solIndexY])(iDof);

      solE[i] = (*sol->_Sol[solIndexE])(iDof);

      for(unsigned  k = 0; k < 3; k++) {
        unsigned solIndex = (k == 0) ? solIndexZ : ((k == 1) ? solIndexX : solIndexY);
        unsigned solPdeIndex = (k == 0) ? solPdeIndexZ : ((k == 1) ? solPdeIndexX : solPdeIndexY);
        sysDof[k * nDofs + i] = pdeSys->GetSystemDof(solIndex, solPdeIndex, i, iel);
      }
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // local to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }

    s.new_recording();

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solType]->Jacobian(coordX, ig, weight, phi, gradPhi);

      double ZOldg = 0.;
      adept::adouble Zg = 0.;
      adept::adouble Xg = 0.;
      adept::adouble Yg = 0.;
      std::vector<double> xg(dim,0.);

      double r = 0;

      std::vector < adept::adouble > gradZg(dim, 0.);
      std::vector < adept::adouble > gradXg(dim, 0.);
      std::vector < adept::adouble > gradYg(dim, 0.);

      for(unsigned i = 0; i < nDofs; i++) {

        ZOldg += solZOld[i] * phi[i];

        Zg += solZ[i] * phi[i];
        Xg += solX[i] * phi[i];
        Yg += solY[i] * phi[i];

        r += solE[i] * phi[i];

        for(unsigned j = 0; j < dim; j++) {
          gradZg[j] += solZ[i] * gradPhi[i * dim + j];
          gradXg[j] += solX[i] * gradPhi[i * dim + j];
          gradYg[j] += solY[i] * gradPhi[i * dim + j];

          xg[j] += coordX[j][i] * phi[i];
        }
      }

      double alpha = 1.0e-7;
      double beta = 0.25;

     // double r = GetTargetSolution(xg,time);
      double mu = 1.;

      // *** phiA_i loop ***
      for(unsigned i = 0; i < nDofs; i++) {

        adept::adouble aResZ = (Zg - ZOldg) / dt * phi[i] - BBs * Xg / alpha * phi[i];
        adept::adouble aResX = - CsC * (r - (1. - beta) * Zg - beta * Yg) * phi[i];
        adept::adouble aResY = - BBs * Xg / alpha * phi[i];


        for(unsigned j = 0; j < dim; j++) {  // second index j in each equation
          aResZ +=  mu * gradPhi[i * dim + j] * gradZg[j]; // diffusion
          aResX +=  mu * gradPhi[i * dim + j] * gradXg[j]; // diffusion
          aResY +=  mu * gradPhi[i * dim + j] * gradYg[j]; // diffusion
        }

        aRes[0 * nDofs + i] += aResZ * weight;
        aRes[1 * nDofs + i] += aResX * weight;
        aRes[2 * nDofs + i] += aResY * weight;

      } // end phiA_i loop
    }

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    res.resize(nDofsAll);
    //copy the value of the adept::adoube mRes in double Res and store
    for(int i = 0; i < nDofsAll; i++) {
      res[i] = -aRes[i].value();
    }
    RES->add_vector_blocked(res, sysDof);

    // define the dependent variables
    s.dependent(aRes.data(), nDofsAll);
    s.independent(solZ.data(), nDofs);
    s.independent(solX.data(), nDofs);
    s.independent(solY.data(), nDofs);

    Jac.resize(nDofsAll * nDofsAll);
    // get the jacobian matrix (ordered by column)
    s.jacobian(Jac.data(), true);

    RES->add_vector_blocked(res, sysDof);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();
  KK->close();

}
