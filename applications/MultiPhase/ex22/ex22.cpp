/** \file Ex6.cpp
 *  \brief This example shows how to set and solve the weak form
 *   of the Navier-Stokes Equation
 *
 *  \f{eqnarray*}
 *  && \mathbf{V} \cdot \nabla \mathbf{V} - \nabla \cdot \nu (\nabla \mathbf{V} +(\nabla \mathbf{V})^T)
 *  +\nabla P = 0 \\
 *  && \nabla \cdot \mathbf{V} = 0
 *  \f}
 *  in a unit box domain (in 2D and 3D) with given vertical velocity 1 on
 *  the left boundary and walls elsewhere.
 *  \author Eugenio Aulisa
 */

#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"

//#include "TransientSystem.hpp"
//#include "NonLinearImplicitSystem.hpp"
#include "LinearImplicitSystem.hpp"


#include "adept.h"

//#include "PolynomialBases.hpp"

//#include "CutFemWeight.hpp"

//#include "CDWeights.hpp"


#include "petsc.h"
#include "petscmat.h"
#include "PetscMatrix.hpp"
#include "Fem.hpp"

#include "./include/ConicAdaptiveRefinement.hpp"
void GetError(MultiLevelProblem & ml_prob, const std::vector <double> &A, std::vector <double> &globalErr);
void GetError2(MultiLevelProblem & ml_prob, const std::vector <double> &A, std::vector <double> &globalErr);
void InitCurvature(MultiLevelProblem & ml_prob, const std::vector <double> &A);
void ProjectSolution(MultiLevelSolution & mlSol, MultiLevelSolution & mlSol1, SparseMatrix * Pj[3]);
// // RT
// const double mu2 = 0.851857;
// const double mu1 = 0.567904;
// // const double mu2 = 0.009535;
// // const double mu1 = 0.006349;
// const double rho2 = 1.5;
// const double rho1 = 1.;
// const double sigma = 0.;
// const double gravity = -10.;

// // Turek 1
// const double mu1 = 1;
// const double mu2 = 10.;
// const double rho1 = 100.;
// const double rho2 = 1000;
// const double sigma = 24.5;
// const double gravity = -0.98;

// Turek 2

std::vector <double> A = {1., 0, 1., 0, 0, -1.};

const double mu1 = 1.;
const double mu2 = 1.;
const double rho1 = 1.;
const double rho2 = 1.;
const double sigma = 1;
const double gravity = -0.;
const double dt = 0.1;
std::vector <double> g = {0, 0, 0};

#include "./include/AllPenalty.hpp"

void AssembleError(Data *data, const std::vector <double> &phiV, const std::vector <double> &phiV_x, const std::vector <double> &phiP, const std::vector <double> &phiP_x,
                   const double &C, const double &weight, const double &weight1, const double &weight2, const double &weightI,
                   const std::vector <double> &N, const std::vector<double> &kappa, const double &dsN);

void AssembleError2(Data *data, const std::vector <double> &phiV, const std::vector <double> &phiV_x, const std::vector <double> &phiP, const std::vector <double> &phiP_x,
                    const double &C, const double &weight, const double &weight1, const double &weight2, const double &weightI,
                    const std::vector <double> &N, const std::vector<double> &kappa, const double &dsN);


#define RADIUS 0.25
#define XG 0.5
#define YG 0.5
#define ZG 0.
// #define RADIUS 0.25
// #define XG 0.
// #define YG -0.2
// #define ZG 0.

using namespace femus;

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet

  if(!strcmp(SolName, "U")) {  // strcmp compares two string in lexiographic sense.
    value = 0.;
  }
  else if(!strcmp(SolName, "V")) {
    //if(facename == 2 || facename == 4) dirichlet = false;
    value = 0.;
//     if(x[0] < 0. && x[1] < 0.5 && x[1] > -0.5 && x[2] < 0.5 && x[2] > -0.5) value = 1.;
  }
  else if(!strcmp(SolName, "W")) {
    value = 0.;
  }
  else if(!strcmp(SolName, "P1")) {
    dirichlet = false;
    value = 0.;
  }
  else if(!strcmp(SolName, "P2")) {
    if(x[0] > -1.999999 || x[1] > -1.999999) dirichlet = false;
    value = 0;

    // if (x[0] > -1.999 && x[0] < 1.999 && x[1] > -1.999 && x[1] < 1.999) {
    //   dirichlet = false;
    //   value = 0.;
    // }
  }
  else if(!strcmp(SolName, "K1")) {
    value = 0.;
  }


  return dirichlet;
}

void AssembleMultiphase(MultiLevelProblem& ml_prob);

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  FEOrder Vorder = SECOND;
  FEOrder Porder = SECOND;
  unsigned nx, ny;
  nx = ny = 64;

  MultiLevelMesh mlMsh;
  mlMsh.GenerateCoarseBoxMesh(nx, ny, 0, -2., 2., -2., 2., 0., 0., QUAD9, "fifth"); // Turek 1&2
  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels, numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);
  mlMsh.PrintInfo();
  unsigned dim = mlMsh.GetDimension();

  MultiLevelSolution mlSol(&mlMsh);
  mlSol.AddSolution("U", LAGRANGE, Vorder);
  mlSol.AddSolution("V", LAGRANGE, Vorder);
  if(dim == 3) mlSol.AddSolution("W", LAGRANGE, Vorder);
  mlSol.AddSolution("P1", LAGRANGE, Porder);
  mlSol.AddSolution("P2", LAGRANGE, Porder);
  mlSol.AddSolution("C0", DISCONTINUOUS_POLYNOMIAL, ZERO, 1, false);
  mlSol.AddSolution("C1", LAGRANGE, Porder, 1, false);
  mlSol.AddSolution("C2", LAGRANGE, Porder, 1, false);
  mlSol.AddSolution("cnt", LAGRANGE, Porder, 1, false);
  mlSol.AddSolution("K0", DISCONTINUOUS_POLYNOMIAL, ZERO, 1, false);
  mlSol.AddSolution("K1", LAGRANGE, Porder, 1);
  mlSol.Initialize("All");
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("All");

  MultiLevelProblem mlProb(&mlSol);
  LinearImplicitSystem& system = mlProb.add_system < LinearImplicitSystem > ("NS");
  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");
  if(dim == 3) system.AddSolutionToSystemPDE("W");
  system.AddSolutionToSystemPDE("P1");
  system.AddSolutionToSystemPDE("P2");
  system.SetSparsityPatternMinimumSize(500);
  system.SetAssembleFunction(AssembleMultiphase);
  system.init();
  system.SetOuterSolver(PREONLY);

  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");
  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(true);
  InitCurvature(mlProb, A);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);
  system.MGsolve();
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 1);

  std::vector<double> globaErr;
  GetError(mlProb, A, globaErr);

  MultiLevelMesh mlMsh1;
  mlMsh1.GenerateCoarseBoxMesh(nx, ny, 0, -2., 2., -2., 2., 0., 0., QUAD9, "fifth"); // Turek 1&2
  numberOfUniformLevels = 2;
  numberOfSelectiveLevels = 0;
  mlMsh1.RefineMesh(numberOfUniformLevels, numberOfUniformLevels + numberOfSelectiveLevels, NULL);


  Mesh* mesh1 = mlMsh1.GetLevel(1);
  SparseMatrix* Pj[3];
  Pj[0] = mesh1->GetCoarseToFineProjection(0);
  Pj[1] = mesh1->GetCoarseToFineProjection(1);
  Pj[2] = mesh1->GetCoarseToFineProjection(2);

  mlMsh1.EraseCoarseLevels(numberOfUniformLevels - 1);
  mlMsh1.PrintInfo();

  MultiLevelSolution mlSol1(&mlMsh1);
  mlSol1.AddSolution("U", LAGRANGE, Vorder);
  mlSol1.AddSolution("V", LAGRANGE, Vorder);
  if(dim == 3) mlSol.AddSolution("W", LAGRANGE, Vorder);
  mlSol1.AddSolution("P1", LAGRANGE, Porder);
  mlSol1.AddSolution("P2", LAGRANGE, Porder);
  mlSol1.AddSolution("C0", DISCONTINUOUS_POLYNOMIAL, ZERO, 1, false);
  mlSol1.AddSolution("C1", LAGRANGE, Porder, 1, false);
  mlSol1.AddSolution("C2", LAGRANGE, Porder, 1, false);
  mlSol1.AddSolution("cnt", LAGRANGE, Porder, 1, false);
  mlSol1.AddSolution("K0", DISCONTINUOUS_POLYNOMIAL, ZERO, 1, false);
  mlSol1.AddSolution("K1", LAGRANGE, Porder, 1);

  mlSol1.AddSolution("Uc", LAGRANGE, Vorder, 1, false);
  mlSol1.AddSolution("Vc", LAGRANGE, Vorder, 1, false);
  mlSol1.AddSolution("P1c", LAGRANGE, Porder, 1, false);
  mlSol1.AddSolution("P2c", LAGRANGE, Porder, 1, false);

  mlSol1.Initialize("All");
  mlSol1.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

  mlSol1.GenerateBdc("All");

  ProjectSolution(mlSol, mlSol1, Pj);

  MultiLevelProblem mlProb1(&mlSol1);
  LinearImplicitSystem& system1 = mlProb1.add_system < LinearImplicitSystem > ("NS");
  system1.AddSolutionToSystemPDE("U");
  system1.AddSolutionToSystemPDE("V");
  if(dim == 3) system1.AddSolutionToSystemPDE("W");
  system1.AddSolutionToSystemPDE("P1");
  system1.AddSolutionToSystemPDE("P2");
  system1.SetSparsityPatternMinimumSize(500);
  system1.SetAssembleFunction(AssembleMultiphase);
  system1.init();
  system1.SetOuterSolver(PREONLY);

  variablesToBePrinted.push_back("All");
  VTKWriter vtkIO1(&mlSol1);
  vtkIO1.SetDebugOutput(true);
  InitCurvature(mlProb1, A);
  vtkIO1.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);
  system1.MGsolve();
  vtkIO1.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 1);
  std::cout << "Errors using analytic solution\n";
  std::vector<double> globaErr1, globaErr2;
  GetError(mlProb1, A, globaErr1);
  std::cout << "Errors using two successive solutions\n";
  GetError2(mlProb1, A, globaErr2);

  std::cout << "nx = " << nx << " ; " << globaErr[0] << " ; " << globaErr[1] << " ; " << globaErr[2] << " ; " << globaErr[3] << " ; "
            << globaErr[4] << " ; " << globaErr[5] << " ; " << globaErr[6] << " ; " << globaErr[7] << " ; " << globaErr[8] << std::endl;
  std::cout << "conv order ; " << log(globaErr[0] / globaErr1[0]) / log(2) << " ; " << log(globaErr[1] / globaErr1[1]) / log(2) << " ; " << log(globaErr[2] / globaErr1[2]) / log(2) << " ; " << log(globaErr[3] / globaErr1[3]) / log(2) << " ; "
            << log(globaErr[4] / globaErr1[4]) / log(2) << " ; " << log(globaErr[5] / globaErr1[5]) / log(2) << " ; " << log(globaErr[6] / globaErr1[6]) / log(2) << " ; " << log(globaErr[7] / globaErr1[7]) / log(2) << " ; "
            << log(globaErr[8] / globaErr1[8]) / log(2) << std::endl;;
  std::cout << "nx = " << 2 * nx << " ; " << globaErr1[0] << " ; " << globaErr1[1] << " ; " << globaErr1[2] << " ; " << globaErr1[3] << " ; "
            << globaErr1[4] << " ; " << globaErr1[5] << " ; " << globaErr1[6] << " ; " << globaErr1[7] << " ; " << globaErr1[8] << std::endl;

  return 0;
}

void ProjectSolution(MultiLevelSolution & mlSol, MultiLevelSolution & mlSol1, SparseMatrix * Pj[3]) {
  Solution* sol0 = mlSol.GetLevel(0);
  Solution* sol1 = mlSol1.GetLevel(0);

  unsigned indexU0 = sol0->GetIndex("U");
  unsigned indexU1 = sol1->GetIndex("Uc");
  unsigned indexV0 = sol0->GetIndex("V");
  unsigned indexV1 = sol1->GetIndex("Vc");
  unsigned solTypeV = sol1->GetSolutionType(indexU1);

  sol1->_Sol[indexU1]->matrix_mult(*sol0->_Sol[indexU0], *Pj[solTypeV]);
  sol1->_Sol[indexU1]->close();

  sol1->_Sol[indexV1]->matrix_mult(*sol0->_Sol[indexV0], *Pj[solTypeV]);
  sol1->_Sol[indexV1]->close();


  unsigned indexP10 = sol0->GetIndex("P1");
  unsigned indexP11 = sol1->GetIndex("P1c");
  unsigned indexP20 = sol0->GetIndex("P2");
  unsigned indexP21 = sol1->GetIndex("P2c");
  unsigned solTypeP = sol1->GetSolutionType(indexP11);

  sol1->_Sol[indexP11]->matrix_mult(*sol0->_Sol[indexP10], *Pj[solTypeP]);
  sol1->_Sol[indexP11]->close();

  sol1->_Sol[indexP21]->matrix_mult(*sol0->_Sol[indexP20], *Pj[solTypeP]);
  sol1->_Sol[indexP21]->close();





}

void InitCurvature(MultiLevelProblem & ml_prob, const std::vector <double> &A) {

  const unsigned level = ml_prob._ml_msh->GetNumberOfLevels() - 1 ;
  Solution *sol = ml_prob._ml_sol->GetSolutionLevel(level);
  Mesh *msh = ml_prob._ml_msh->GetLevel(level);
  elem *el = msh->el;

  const unsigned  dim = msh->GetDimension();
  unsigned iproc = msh->processor_id();

  unsigned solC0Index = sol->GetIndex("C0");
  unsigned solC1Index = sol->GetIndex("C1");

  unsigned solK0Index = sol->GetIndex("K0");
  unsigned solK1Index = sol->GetIndex("K1");

  unsigned solPKCType = sol->GetSolutionType(solK1Index);

  unsigned solCntIndex = sol->GetIndex("cnt");

  unsigned xvType = 2;

  sol->_Sol[solK0Index]->zero();
  sol->_Sol[solK1Index]->zero();

  std::vector<std::vector<double>> xT;

  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    unsigned elType = msh->GetElementType(iel);
    unsigned nDofsK1 = msh->GetElementDofNumber(iel, solPKCType);

    unsigned nDofsX = msh->GetElementDofNumber(iel, 0);
    xT.assign(nDofsX, std::vector<double>(dim, 0));
    for(unsigned i = 0; i < nDofsX; i++) {
      unsigned xdof  = msh->GetSolutionDof(i, iel, xvType);
      for(unsigned k = 0; k < dim; k++) {
        xT[i][k] = (*msh->_topology->_Sol[k])(xdof);      // global extraction and local storage for the element coordinates
      }
    }
    std::vector <double> Ap;
    ConicAdaptiveRefinement::GetConicsInTargetElement(xT, A, elType, Ap);

    int testIntersection = ConicAdaptiveRefinement::TestIfIntesectionWithReferenceElement(Ap, elType);
    if(testIntersection == 0) {
      sol->_Sol[solC0Index]->set(iel, 0.5);
    }
    else if(testIntersection < 0) {
      sol->_Sol[solC0Index]->set(iel, 1.);
    }




    std::vector<double> xv(dim);
    bool atLeastOneIsZero = false;
    bool atLeastOneIsOne = false;
    double kappaIel = 0.;
    for(unsigned i = 0; i < nDofsK1; i++) {
      unsigned xdof  = msh->GetSolutionDof(i, iel, xvType);    // local to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        xv[k] = (*msh->_topology->_Sol[k])(xdof);      // global extraction and local storage for the element coordinates
      }

      unsigned idof  = msh->GetSolutionDof(i, iel, solPKCType);    // local to global mapping between coordinates node and coordinate dof
      double value = ConicAdaptiveRefinement::EvaluateConic(xv, A);
      if(value < 0.) {  //1.0e-10) {
        sol->_Sol[solC1Index]->set(idof, 1.);
        atLeastOneIsOne = true;
      }
      else { // if (value > 1.0e-10) {
        atLeastOneIsZero = true;
      }
      std::vector <double> kappa;
      ConicAdaptiveRefinement::GetConicCurvature(xv, A, kappa, false);
      if(std::isnan(kappa[0])) kappa[0] = 0.;
      else sol->_Sol[solK1Index]->set(idof, kappa[0]);
      kappaIel += kappa[0];
    }
    // if (!atLeastOneIsOne) sol->_Sol[solC0Index]->set(iel, 0.);
    // else if (!atLeastOneIsZero) sol->_Sol[solC0Index]->set(iel, 1.);
    // else sol->_Sol[solC0Index]->set(iel, 0.5);

    sol->_Sol[solK0Index]->set(iel, kappaIel / nDofsK1);
  }
  sol->_Sol[solC0Index]->close();
  sol->_Sol[solC1Index]->close();

  sol->_Sol[solK0Index]->close();
  sol->_Sol[solK1Index]->close();

  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    if((*sol->_Sol[solC0Index])(iel) == 0.5) {
      unsigned nDofsK1 = msh->GetElementDofNumber(iel, solPKCType);
      for(unsigned i = 0; i < nDofsK1; i++) {
        unsigned idof  = msh->GetSolutionDof(i, iel, solPKCType);
        double C1 = (*sol->_Sol[solC1Index])(idof);
        if(C1 < 0.2) {
          sol->_Sol[solC1Index]->set(idof, 0.25);
        }
        else if(C1 > 0.8) {
          sol->_Sol[solC1Index]->set(idof, 0.75);
        }
      }
    }
  }

  sol->_Sol[solC1Index]->close();

  unsigned numberOfSmoothings = 100;
  for(unsigned k = 0; k < numberOfSmoothings; k++) {

    //From the element to the nodes
    for(unsigned i = msh->_dofOffset[solPKCType][iproc]; i < msh->_dofOffset[solPKCType][iproc + 1]; i++) {
      double C1 = (*sol->_Sol[solC1Index])(i);
      if(C1 < 0.2 || C1 > 0.8) sol->_Sol[solK1Index]->set(i, 0.);
    }
    sol->_Sol[solK1Index]->close();
    sol->_Sol[solCntIndex]->zero();
    for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
      if((*sol->_Sol[solC0Index])(iel) > 0.9 || (*sol->_Sol[solC0Index])(iel) < 0.1) {
        double K0 = (*sol->_Sol[solK0Index])(iel);
        unsigned nDofsK1 = msh->GetElementDofNumber(iel, solPKCType);
        for(unsigned i = 0; i < nDofsK1; i++) {
          unsigned idof  = msh->GetSolutionDof(i, iel, solPKCType);
          double C1 = (*sol->_Sol[solC1Index])(idof);
          if((*sol->_Bdc[solK1Index])(idof) > 0 && (C1 > 0.8 || C1 < 0.2)) {
            sol->_Sol[solK1Index]->add(idof, K0);
            sol->_Sol[solCntIndex]->add(idof, 1.);
          }
        }
      }
    }
    sol->_Sol[solK1Index]->close();
    sol->_Sol[solCntIndex]->close();

    for(unsigned i = msh->_dofOffset[solPKCType][iproc]; i < msh->_dofOffset[solPKCType][iproc + 1]; i++) {
      double cnt = (*sol->_Sol[solCntIndex])(i);
      if(cnt > 0.1) {
        double K1 = (*sol->_Sol[solK1Index])(i);
        sol->_Sol[solK1Index]->set(i, K1 / cnt);
      }
    }
    sol->_Sol[solK1Index]->close();

    if(k < numberOfSmoothings - 1) {
      //from the nodes to the element
      for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
        if((*sol->_Sol[solC0Index])(iel) > 0.9 || (*sol->_Sol[solC0Index])(iel) < 0.1) {
          double K0 = 0.;
          unsigned nDofsK1 = msh->GetElementDofNumber(iel, solPKCType);
          for(unsigned i = 0; i < nDofsK1; i++) {
            unsigned idof  = msh->GetSolutionDof(i, iel, solPKCType);
            double K1 = (*sol->_Sol[solK1Index])(idof);
            K0 += K1;
          }
          sol->_Sol[solK0Index]->set(iel, K0 / nDofsK1);
        }
      }
      sol->_Sol[solK0Index]->close();
    }
  }
}

//Attempting to create J by hand
void AssembleMultiphase(MultiLevelProblem & ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled


  //  extract pointers to the several objects that we are going to use
  LinearImplicitSystem* mlPdeSys = &ml_prob.get_system<LinearImplicitSystem> ("NS");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual std::vector object in pdeSys (level)

  MatResetPreallocation((static_cast< PetscMatrix* >(KK))->mat());
  MatSetOption((static_cast< PetscMatrix* >(KK))->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  KK->zero();
  RES->zero();

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  std::vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("V");    // get the position of "V" in the ml_sol object
  if(dim == 3) solVIndex[2] = mlSol->GetIndex("W");       // get the position of "V" in the ml_sol object

  unsigned solVType = mlSol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

  unsigned solP1Index = mlSol->GetIndex("P1");    // get the position of "P1" in the ml_sol object
  unsigned solP2Index = mlSol->GetIndex("P2");    // get the position of "P2" in the ml_sol object
  unsigned solPKCType = mlSol->GetSolutionType(solP1Index);    // get the finite element type for "u"

  unsigned solK1Index = mlSol->GetIndex("K1");
  unsigned solC0Index = mlSol->GetIndex("C0");
  unsigned solC1Index = mlSol->GetIndex("C1");
  unsigned solC2Index = mlSol->GetIndex("C2");
  unsigned solCntIndex = mlSol->GetIndex("cnt");

  std::vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("U");    // get the position of "U" in the pdeSys object
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("V");    // get the position of "V" in the pdeSys object
  if(dim == 3) solVPdeIndex[2] = mlPdeSys->GetSolPdeIndex("W");

  unsigned solP1PdeIndex = mlPdeSys->GetSolPdeIndex("P1");    // get the position of "P" in the pdeSys object
  unsigned solP2PdeIndex = mlPdeSys->GetSolPdeIndex("P2");    // get the position of "P" in the pdeSys object

  std::vector < std::vector < double > >  solV(dim);    // local solution
  std::vector < double >  solP1; // local solution
  std::vector < double >  solP2; // local solution
  std::vector < double >  solC1; // local solution
  std::vector < double >  solK1; // local solution

  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector< unsigned > sysDof; // local to global pdeSys dofs
  std::vector< double > Res; // local redidual std::vector
  std::vector < double > Jac;

  clock_t start_time = clock();

  ConicAdaptiveRefinement cad;

  sol->_Sol[solC0Index]->zero();
  sol->_Sol[solC2Index]->zero();
  sol->_Sol[solCntIndex]->zero();

  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    unsigned elType = msh->GetElementType(iel);
    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType);
    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPKCType);    // number of solution element dofs
    unsigned nDofsVP = dim * nDofsV + 2 * nDofsP;

    // resize local arrays
    sysDof.resize(nDofsVP);
    Res.assign(nDofsVP, 0.);
    Jac.assign(nDofsVP * nDofsVP, 0.);

    for(unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofsV);
      coordX[k].resize(nDofsX);
    }
    solP1.resize(nDofsP);
    solP2.resize(nDofsP);
    solC1.resize(nDofsP);
    solK1.resize(nDofsP);

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // local to global mapping between solution node and solution dof

      for(unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);
        sysDof[k * nDofsV + i] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }

    for(unsigned i = 0; i < nDofsP; i++) {
      unsigned iDof = msh->GetSolutionDof(i, iel, solPKCType);    // local to global mapping between solution node and solution dof
      solP1[i] = (*sol->_Sol[solP1Index])(iDof);      // global extraction and local storage for the solution
      solP2[i] = (*sol->_Sol[solP2Index])(iDof);      // global extraction and local storage for the solution
      solC1[i] = (*sol->_Sol[solC1Index])(iDof);      // global extraction and local storage for the solution
      solK1[i] = (*sol->_Sol[solK1Index])(iDof);      // global extraction and local storage for the solution
      sysDof[dim * nDofsV + i ] = pdeSys->GetSystemDof(solP1Index, solP1PdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
      sysDof[dim * nDofsV + nDofsP + i ] = pdeSys->GetSystemDof(solP2Index, solP2PdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // local to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }

    std::vector <double> Ap;
    bool distributeCurvature = true;
    double eps = 1.0e-5;

    Data *data = new Data(solVType, solPKCType, solV, solP1, solP2, solK1, solC1, Res, Jac, coordX, elType, rho1, rho2, mu1, mu2, sigma, dt, g, A, distributeCurvature, eps);

    cad.SetDataPointer(data);

    std::vector<std::vector<double>> xT;
    cad.TransposeXinParentElement(coordX, elType, xT);
    std::vector<std::vector<double>> yi = cad.GetXiInParentElement();

    cad.GetConicsInTargetElement(xT, A, elType, Ap);
    tuple <double, double, double> a = cad.Assemble(AssembleNavierStokes, 5, xT, yi, Ap);

    double C = std::get<0>(a) / (std::get<0>(a) + std::get<1>(a));

    sol->_Sol[solC0Index]->set(iel, C);

    unsigned nDofs1 = msh->GetElementDofNumber(iel, solPKCType);
    for(unsigned i = 0; i < nDofs1; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, solPKCType);
      sol->_Sol[solC2Index]->add(idof, C);
      sol->_Sol[solCntIndex]->add(idof, 1);
    }

    RES->add_vector_blocked(Res, sysDof);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);
  } //end element loop for each process

  sol->_Sol[solC0Index]->close();
  sol->_Sol[solC2Index]->close();
  sol->_Sol[solCntIndex]->close();

  for(unsigned i = msh->_dofOffset[solPKCType][iproc]; i < msh->_dofOffset[solPKCType][iproc + 1]; i++) {
    double value = (*sol->_Sol[solC2Index])(i);
    double cnt = (*sol->_Sol[solCntIndex])(i);
    sol->_Sol[solC2Index]->set(i, value / cnt);
  }
  sol->_Sol[solC2Index]->close();

  std::cout << "Navier-Stokes Assembly time = " << static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC << std::endl << std::endl;

  AssembleAllPenalty(ml_prob);

  RES->close();
  KK->close();

  //KK->draw();
}

void GetError(MultiLevelProblem & ml_prob, const std::vector <double> &A, std::vector<double> &error) {

  const unsigned level = ml_prob._ml_msh->GetNumberOfLevels() - 1 ;
  Solution *sol = ml_prob._ml_sol->GetSolutionLevel(level);
  Mesh *msh = ml_prob._ml_msh->GetLevel(level);
  elem *el = msh->el;

  const unsigned  dim = msh->GetDimension();
  unsigned    iproc = msh->processor_id();

  //solution variable
  std::vector < unsigned > solVIndex(dim);
  solVIndex[0] = sol->GetIndex("U");
  solVIndex[1] = sol->GetIndex("V");
  if(dim == 3) solVIndex[2] = sol->GetIndex("W");

  unsigned solVType = sol->GetSolutionType(solVIndex[0]);

  unsigned solP1Index = sol->GetIndex("P1");
  unsigned solP2Index = sol->GetIndex("P2");
  unsigned solC2Index = sol->GetIndex("C2");
  unsigned solK1Index = sol->GetIndex("K1");

  unsigned solPKCType = sol->GetSolutionType(solP1Index);

  std::vector < std::vector < double > >  solV(dim);
  std::vector < double >  solP1;
  std::vector < double >  solP2;
  std::vector < double >  solC2;
  std::vector < double >  solK1;

  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector< double > localErr(2 * dim + 3, 0.);
  std::vector< double > globalErr(2 * dim + 3, 0.);
  std::vector < double > Jac; // support variable not used

  clock_t start_time = clock();

  ConicAdaptiveRefinement cad;

  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    unsigned elType = msh->GetElementType(iel);
    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType);
    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPKCType);    // number of solution element dofs

    unsigned nDofsVP = dim * nDofsV + 2 * nDofsP;

    // resize local arrays
    for(unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofsV);
      coordX[k].resize(nDofsX);
    }
    solP1.resize(nDofsP);
    solP2.resize(nDofsP);
    solC2.resize(nDofsP);
    solK1.resize(nDofsP);

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // local to global mapping between solution node and solution dof
      for(unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);
      }
    }

    for(unsigned i = 0; i < nDofsP; i++) {
      unsigned iDof = msh->GetSolutionDof(i, iel, solPKCType);    // local to global mapping between solution node and solution dof
      solP1[i] = (*sol->_Sol[solP1Index])(iDof);      // global extraction and local storage for the solution
      solP2[i] = (*sol->_Sol[solP2Index])(iDof);      // global extraction and local storage for the solution
      solC2[i] = (*sol->_Sol[solC2Index])(iDof);
      solK1[i] = (*sol->_Sol[solK1Index])(iDof);
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // local to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }

    bool distributeCurvature = true;
    double eps = 0.;

    Data *data = new Data(solVType, solPKCType, solV, solP1, solP2, solK1, solC2, localErr, Jac, coordX, elType, rho1, rho2, mu1, mu2, sigma, dt, g, A, distributeCurvature, eps);

    cad.SetDataPointer(data);

    std::vector<std::vector<double>> xT;
    cad.TransposeXinParentElement(coordX, elType, xT);
    std::vector<std::vector<double>> yi = cad.GetXiInParentElement();

    std::vector <double> Ap;
    cad.GetConicsInTargetElement(xT, A, elType, Ap);
    cad.Assemble(AssembleError, 5, xT, yi, Ap);

  } //end element loop for each process

  MPI_Allreduce(localErr.data(), globalErr.data(), localErr.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


  error = {sqrt(globalErr[0]), sqrt(globalErr[1]), sqrt(globalErr[2]), sqrt(globalErr[3]),
           sqrt(globalErr[0]) + sqrt(globalErr[2]), sqrt(globalErr[1]) + sqrt(globalErr[3]),
           sqrt(globalErr[4]), sqrt(globalErr[5]), sqrt(globalErr[6])
          };

  std::cout.precision(14);
  std::cout << "U-l2NormErr = "  << sqrt(globalErr[0]) << std::endl;
  std::cout << "V-l2NormErr = "  << sqrt(globalErr[1]) << std::endl;
  std::cout << "U-SemiNormErr = "  << sqrt(globalErr[2]) << std::endl;
  std::cout << "V-SeminNormErr = "  << sqrt(globalErr[3]) << std::endl;
  std::cout << "U-tripleNormErr = "  << sqrt(globalErr[0]) + sqrt(globalErr[2]) << std::endl;
  std::cout << "V-tripleNormErr = "  << sqrt(globalErr[1]) + sqrt(globalErr[3]) << std::endl;

  std::cout << "P1-l2NormErr = " << sqrt(globalErr[4]) << std::endl;
  std::cout << "P2-l2NormErr = " << sqrt(globalErr[5]) << std::endl;
  std::cout << "P-l2NormErr = "  << sqrt(globalErr[6]) << std::endl;

  std::cout << "Error Function time = " << static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC << std::endl << std::endl;
}

void GetError2(MultiLevelProblem & ml_prob, const std::vector <double> &A, std::vector<double> &error) {

  const unsigned level = ml_prob._ml_msh->GetNumberOfLevels() - 1 ;
  Solution *sol = ml_prob._ml_sol->GetSolutionLevel(level);
  Mesh *msh = ml_prob._ml_msh->GetLevel(level);
  elem *el = msh->el;

  const unsigned  dim = msh->GetDimension();
  unsigned    iproc = msh->processor_id();

  //solution variable
  std::vector < unsigned > solVIndex(dim);
  solVIndex[0] = sol->GetIndex("U");
  solVIndex[1] = sol->GetIndex("V");
  if(dim == 3) solVIndex[2] = sol->GetIndex("W");

  std::vector < unsigned > solVcIndex(dim);
  solVcIndex[0] = sol->GetIndex("Uc");
  solVcIndex[1] = sol->GetIndex("Vc");
  if(dim == 3) solVcIndex[2] = sol->GetIndex("Wc");

  unsigned solVType = sol->GetSolutionType(solVIndex[0]);

  unsigned solP1Index = sol->GetIndex("P1");
  unsigned solP1cIndex = sol->GetIndex("P1c");
  unsigned solP2Index = sol->GetIndex("P2");
  unsigned solP2cIndex = sol->GetIndex("P2c");

  unsigned solC2Index = sol->GetIndex("C2");
  unsigned solK1Index = sol->GetIndex("K1");

  unsigned solPKCType = sol->GetSolutionType(solP1Index);

  std::vector < std::vector < double > >  solV(dim);
  std::vector < double >  solP1;
  std::vector < double >  solP2;
  std::vector < double >  solC2;
  std::vector < double >  solK1;

  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector< double > localErr(2 * dim + 3, 0.);
  std::vector< double > globalErr(2 * dim + 3, 0.);

  std::vector < double > Jac; // support variable not used

  clock_t start_time = clock();

  ConicAdaptiveRefinement cad;

  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    unsigned elType = msh->GetElementType(iel);
    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType);
    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPKCType);    // number of solution element dofs

    // resize local arrays
    for(unsigned  k = 0; k < dim; k++) {
      solV[k].resize(2 * nDofsV);
      coordX[k].resize(nDofsX);
    }
    solP1.resize(2 * nDofsP);
    solP2.resize(2 * nDofsP);
    solC2.resize(nDofsP);
    solK1.resize(nDofsP);

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // local to global mapping between solution node and solution dof
      for(unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);
        solV[k][nDofsV + i] = (*sol->_Sol[solVcIndex[k]])(solVDof);
      }
    }

    for(unsigned i = 0; i < nDofsP; i++) {
      unsigned iDof = msh->GetSolutionDof(i, iel, solPKCType);    // local to global mapping between solution node and solution dof
      solP1[i] = (*sol->_Sol[solP1Index])(iDof);
      solP1[nDofsP + i] = (*sol->_Sol[solP1cIndex])(iDof);
      solP2[i] = (*sol->_Sol[solP2Index])(iDof);
      solP2[nDofsP + i] = (*sol->_Sol[solP2cIndex])(iDof); // global extraction and local storage for the solution
      solC2[i] = (*sol->_Sol[solC2Index])(iDof);
      solK1[i] = (*sol->_Sol[solK1Index])(iDof);
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // local to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }

    bool distributeCurvature = true;
    double eps = 0.;

    Data *data = new Data(solVType, solPKCType, solV, solP1, solP2, solK1, solC2, localErr, Jac, coordX, elType, rho1, rho2, mu1, mu2, sigma, dt, g, A, distributeCurvature, eps);

    cad.SetDataPointer(data);

    std::vector<std::vector<double>> xT;
    cad.TransposeXinParentElement(coordX, elType, xT);
    std::vector<std::vector<double>> yi = cad.GetXiInParentElement();

    std::vector <double> Ap;
    cad.GetConicsInTargetElement(xT, A, elType, Ap);
    cad.Assemble(AssembleError2, 5, xT, yi, Ap);

  } //end element loop for each process

  MPI_Allreduce(localErr.data(), globalErr.data(), localErr.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


  error = {sqrt(globalErr[0]), sqrt(globalErr[1]), sqrt(globalErr[2]), sqrt(globalErr[3]),
           sqrt(globalErr[0]) + sqrt(globalErr[2]), sqrt(globalErr[1]) + sqrt(globalErr[3]),
           sqrt(globalErr[4]), sqrt(globalErr[5]), sqrt(globalErr[6])
          };

  std::cout.precision(14);
  std::cout << "U-l2NormErr = "  << sqrt(globalErr[0]) << std::endl;
  std::cout << "V-l2NormErr = "  << sqrt(globalErr[1]) << std::endl;
  std::cout << "U-SemiNormErr = "  << sqrt(globalErr[2]) << std::endl;
  std::cout << "V-SeminNormErr = "  << sqrt(globalErr[3]) << std::endl;
  std::cout << "U-tripleNormErr = "  << sqrt(globalErr[0]) + sqrt(globalErr[2]) << std::endl;
  std::cout << "V-tripleNormErr = "  << sqrt(globalErr[1]) + sqrt(globalErr[3]) << std::endl;


  std::cout << "P1-l2NormErr = " << sqrt(globalErr[4]) << std::endl;
  std::cout << "P2-l2NormErr = " << sqrt(globalErr[5]) << std::endl;
  std::cout << "P-l2NormErr = "  << sqrt(globalErr[6]) << std::endl;

  std::cout << "Error Function time = " << static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC << std::endl;

}


void AssembleError(Data * data, const std::vector <double> &phiV, const std::vector <double> &phiV_x, const std::vector <double> &phiP, const std::vector <double> &phiP_x,
                   const double & C, const double & weight, const double & weight1, const double & weight2, const double & weightI,
                   const std::vector <double> &N, const std::vector<double> &kappa, const double & dsN) {


  const unsigned &dim = data->_V.size();
  const unsigned &nDofsV = phiV.size();
  const unsigned &nDofsP = phiP.size();
  unsigned nDofsVP = dim * nDofsV + 2 * nDofsP;

  std::vector < double > solVg(dim, 0);
  std::vector < std::vector < double > > solVg_x(dim, std::vector<double> (dim, 0));
  for(unsigned i = 0; i < nDofsV; i++) {
    for(unsigned  K = 0; K < dim; K++) {
      solVg[K] += data->_V[K][i] * phiV[i];
      for(unsigned J = 0; J < dim; J++) {
        solVg_x[K][J] += data->_V[K][i] * phiV_x[i * dim + J];
      }
    }
  }
  std::vector < double > solVg_exc(dim, 0);
  std::vector<std::vector < double > > solVg_x_exc(dim, std::vector<double>(dim, 0));

  double solP1g = 0;
  double solP2g = 0;
  double solPg = 0;
  double solP1g_exc = 1;
  double solP2g_exc = 0;
  double solPg_exc = 0;

  for(unsigned i = 0; i < nDofsP; i++) {
    solP1g += phiP[i] * data->_P1[i];
    solP2g += phiP[i] * data->_P2[i];
    solPg += phiP[i] * (data->_C1[i] * data->_P1[i] + (1 - data->_C1[i]) * data->_P2[i]);
    solPg_exc += phiP[i] * (data->_C1[i] * solP1g_exc + (1 - data->_C1[i]) * solP2g_exc);
  }

  for(unsigned K = 0; K < dim; K++) {
    data->_res[K] += (solVg[K] - solVg_exc[K]) * (solVg[K] - solVg_exc[K]) * weight;
    for(unsigned J = 0; J < dim; J++) {
      data->_res[dim + K] += (data->_mu1 * (K == 0) + data->_mu2 * (K == 1)) * (solVg_x[K][J] - solVg_x_exc[K][J]) * (solVg_x[K][J] - solVg_x_exc[K][J]) * weight;
    }
  }
  data->_res[2 * dim + 0] += (solP1g - solP1g_exc) * (solP1g - solP1g_exc) * weight * weight1;
  data->_res[2 * dim + 1] += (solP2g - solP2g_exc) * (solP2g - solP2g_exc) * weight * weight2;
  data->_res[2 * dim + 2] += (solPg - solPg_exc) * (solPg - solPg_exc) * weight;

}


void AssembleError2(Data * data, const std::vector <double> &phiV, const std::vector <double> &phiV_x, const std::vector <double> &phiP, const std::vector <double> &phiP_x,
                    const double & C, const double & weight, const double & weight1, const double & weight2, const double & weightI,
                    const std::vector <double> &N, const std::vector<double> &kappa, const double & dsN) {


  const unsigned &dim = data->_V.size();
  const unsigned &nDofsV = phiV.size();
  const unsigned &nDofsP = phiP.size();

  std::vector < double > solVg(dim, 0);
  std::vector < std::vector < double > > solVg_x(dim, std::vector<double> (dim, 0));
  std::vector < double > solVg_exc(dim, 0);
  std::vector < std::vector < double > > solVg_x_exc(dim, std::vector<double> (dim, 0));

  for(unsigned i = 0; i < nDofsV; i++) {
    for(unsigned  K = 0; K < dim; K++) {
      solVg[K] += data->_V[K][i] * phiV[i];
      solVg_exc[K] += data->_V[K][i + nDofsV] * phiV[i];
      for(unsigned J = 0; J < dim; J++) {
        solVg_x[K][J] += data->_V[K][i] * phiV_x[i * dim + J];
        solVg_x_exc[K][J] += data->_V[K][i + nDofsV] * phiV_x[i * dim + J];
      }
    }
  }


  double solP1g = 0;
  double solP2g = 0;
  double solPg = 0;
  double solP1g_exc = 0.;
  double solP2g_exc = 0.;
  double solPg_exc = 0.;

  for(unsigned i = 0; i < nDofsP; i++) {
    solP1g += phiP[i] * data->_P1[i];
    solP2g += phiP[i] * data->_P2[i];
    solPg += phiP[i] * (data->_C1[i] * data->_P1[i] + (1 - data->_C1[i]) * data->_P2[i]);

    solP1g_exc += phiP[i] * data->_P1[i + nDofsP];
    solP2g_exc += phiP[i] * data->_P2[i + nDofsP];
    solPg_exc += phiP[i] * (data->_C1[i] * data->_P1[i + nDofsP] + (1 - data->_C1[i]) * data->_P2[i + nDofsP]);
  }

  for(unsigned K = 0; K < dim; K++) {
    data->_res[K] += (solVg[K] - solVg_exc[K]) * (solVg[K] - solVg_exc[K]) * weight;
    for(unsigned J = 0; J < dim; J++) {
      data->_res[dim + K] += (data->_mu1 * (K == 0) + data->_mu2 * (K == 1)) * (solVg_x[K][J] - solVg_x_exc[K][J]) * (solVg_x[K][J] - solVg_x_exc[K][J]) * weight;
    }
  }
  data->_res[2 * dim + 0] += (solP1g - solP1g_exc) * (solP1g - solP1g_exc) * weight * weight1;
  data->_res[2 * dim + 1] += (solP2g - solP2g_exc) * (solP2g - solP2g_exc) * weight * weight2;
  data->_res[2 * dim + 2] += (solPg - solPg_exc) * (solPg - solPg_exc) * weight;

}


