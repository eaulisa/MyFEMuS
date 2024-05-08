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
void GetError(MultiLevelProblem & ml_prob);
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
const double mu1 = 1.;
const double mu2 = 1.;
const double rho1 = 1.;
const double rho2 = 1.;
const double sigma = 1;
const double gravity = -0.;
const double dt = 0.1;

std::vector <double> g = {0, 0, 0};

#include "./include/AllPenalty.hpp"

void AssembleError(Data *data, const std::vector <double> &phiV, const std::vector <double> &phiV_x, const std::vector <double> &phiP,
                   const double &C, const double &weight, const double &weight1, const double &weight2, const double &weightI,
                   const std::vector <double> &N, const std::vector<double> &kappa, const double &dsN, const double &eps);


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
    value = 0.;
  }

  return dirichlet;
}

void AssembleMultiphase(MultiLevelProblem& ml_prob);

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
//   mlMsh.ReadCoarseMesh("./input/cube_hex.neu", "seventh", scalingFactor);
//   mlMsh.ReadCoarseMesh("./input/square_quad.neu", "fifth", scalingFactor);
  mlMsh.GenerateCoarseBoxMesh(32, 32, 0, -2., 2., -2., 2., 0., 0., QUAD9, "fifth"); // Turek 1&2
//   mlMsh.GenerateCoarseBoxMesh(64, 256, 0, -0.5, 0.5, -2, 2, 0., 0., QUAD9, "fifth"); //RT
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
     probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();

  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels, numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // erase all the coarse mesh levels
  mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  MultiLevelSolution mlSol(&mlMsh);



  FEOrder Vorder = SECOND;
  FEOrder Porder = SECOND;

  // add variables to mlSol
  mlSol.AddSolution("U", LAGRANGE, Vorder);
  mlSol.AddSolution("V", LAGRANGE, Vorder);
  if(dim == 3) mlSol.AddSolution("W", LAGRANGE, Vorder);
  // mlSol.AddSolution("P1",  DISCONTINUOUS_POLYNOMIAL, ZERO);
  // mlSol.AddSolution("P2",  DISCONTINUOUS_POLYNOMIAL, ZERO);



  mlSol.AddSolution("P1", LAGRANGE, Porder);
  mlSol.AddSolution("P2", LAGRANGE, Porder);

  mlSol.AddSolution("C0", DISCONTINUOUS_POLYNOMIAL, ZERO, 1, false);
  mlSol.AddSolution("C1", LAGRANGE, Porder, 1, false);
  mlSol.AddSolution("cnt", LAGRANGE, Porder, 1, false);

  mlSol.AddSolution("K0", DISCONTINUOUS_POLYNOMIAL, ZERO, 1, false);
  mlSol.AddSolution("K1", LAGRANGE, Porder, 1, false);


//    //Taylor-hood
//    mlSol.AddSolution("U", LAGRANGE, SERENDIPITY);
//    mlSol.AddSolution("V", LAGRANGE, SERENDIPITY);
//    if (dim == 3) mlSol.AddSolution("W", LAGRANGE, SERENDIPITY);
//    mlSol.AddSolution("P", LAGRANGE, FIRST);

//    // Bad FEM pair - no LBB condition
//    mlSol.AddSolution("U", LAGRANGE, FIRST);
//    mlSol.AddSolution("V", LAGRANGE, FIRST);
//    if (dim == 3) mlSol.AddSolution("W", LAGRANGE, FIRST);
//    mlSol.AddSolution("P", LAGRANGE, FIRST);


  mlSol.Initialize("All");


  Solution* sol = mlSol.GetSolutionLevel(mlMsh.GetNumberOfLevels() - 1);

  //cld = new Cloud(sol);
  //cldint = new Cloud(sol);

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
//   mlSol.FixSolutionAtOnePoint("P1");
  //mlSol.FixSolutionAtOnePoint("P2");

  mlSol.GenerateBdc("All");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);

  // add system Poisson in mlProb as a Linear Implicit System
  LinearImplicitSystem& system = mlProb.add_system < LinearImplicitSystem > ("NS");

  // add solution "u" to system
  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");

  if(dim == 3) system.AddSolutionToSystemPDE("W");

  system.AddSolutionToSystemPDE("P1");
  system.AddSolutionToSystemPDE("P2");

  system.SetSparsityPatternMinimumSize(250);

  // attach the assembling function to system
  system.SetAssembleFunction(AssembleMultiphase);

  // initilaize and solve the system
  system.init();
  system.SetOuterSolver(PREONLY);

  // BEGIN Testing the class Cloud

  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");
  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(true);

  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);
  system.MGsolve();
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 1);
  system.MGsolve();
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 2);

  GetError(mlProb);

  return 0;
}

void InitCurvature(Solution* sol, const std::vector<double> &A) {
  Mesh* msh = sol->GetMesh();
  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
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

  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofsK1 = msh->GetElementDofNumber(iel, solPKCType);

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
      if(ConicAdaptiveRefinement::EvaluateConic(xv, A) < 0) {
        sol->_Sol[solC1Index]->set(idof, 1.);
        atLeastOneIsOne = true;
      }
      else {
        atLeastOneIsZero = true;
      }
      std::vector <double> kappa;
      ConicAdaptiveRefinement::GetConicCurvature(xv, A, kappa, false);
      if(std::isnan(kappa[0])) kappa[0] = 0.;
      else sol->_Sol[solK1Index]->set(idof, kappa[0]);
      kappaIel += kappa[0];
    }
    if(!atLeastOneIsOne) sol->_Sol[solC0Index]->set(iel, 0.);
    else if(!atLeastOneIsZero) sol->_Sol[solC0Index]->set(iel, 1.);
    else sol->_Sol[solC0Index]->set(iel, 0.5);

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

  unsigned numberOfSmoothings = 10;
  for(unsigned k = 0; k < numberOfSmoothings; k++) {

    //From the element to the nodes
    for(unsigned i = msh->_dofOffset[solPKCType][iproc]; i < msh->_dofOffset[solPKCType][iproc + 1]; i++) {
      double C1 = (*sol->_Sol[solC1Index])(i);
      if(C1 < 0.2 || C1 > 0.8) sol->_Sol[solK1Index]->set(i, 0.);
    }
    sol->_Sol[solK1Index]->close();
    sol->_Sol[solCntIndex]->zero();
    for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
      if((*sol->_Sol[solC0Index])(iel) > 0.9) {
        double K0 = (*sol->_Sol[solK0Index])(iel);
        unsigned nDofsK1 = msh->GetElementDofNumber(iel, solPKCType);
        for(unsigned i = 0; i < nDofsK1; i++) {
          unsigned idof  = msh->GetSolutionDof(i, iel, solPKCType);
          double C1 = (*sol->_Sol[solC1Index])(idof);
          if(C1 > 0.8) {
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
        if((*sol->_Sol[solC0Index])(iel) > 0.9) {
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
  LinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<LinearImplicitSystem> ("NS");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual std::vector object in pdeSys (level)

  MatResetPreallocation((static_cast< PetscMatrix* >(KK))->mat());
  MatSetOption((static_cast< PetscMatrix* >(KK))->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  KK->zero();
  RES->zero();

  std::vector <double> A = {1., 0, 1., 0, 0, -1.};
  InitCurvature(sol, A);

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

  std::vector < double >  solK1; // local solution

  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector< unsigned > sysDof; // local to global pdeSys dofs
  std::vector< double > Res; // local redidual std::vector
  std::vector < double > Jac;

  clock_t start_time = clock();

  ConicAdaptiveRefinement cad;

  sol->_Sol[solC0Index]->zero();
  sol->_Sol[solC1Index]->zero();
  sol->_Sol[solCntIndex]->zero();

  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPKCType);    // number of solution element dofs

    unsigned nDofsVP = dim * nDofsV + 2 * nDofsP;

    unsigned nDofsK1 = msh->GetElementDofNumber(iel, solPKCType);


    // resize local arrays
    sysDof.resize(nDofsVP);
    Res.assign(nDofsVP, 0.);
    Jac.assign(nDofsVP * nDofsVP, 0.);

    for(unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofsV);
      coordX[k].resize(nDofsV);
    }
    solP1.resize(nDofsP);
    solP2.resize(nDofsP);

    solK1.resize(nDofsK1);

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // local to global mapping between solution node and solution dof

      for(unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);
        sysDof[k * nDofsV + i] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }

    for(unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPKCType);    // local to global mapping between solution node and solution dof
      solP1[i] = (*sol->_Sol[solP1Index])(solPDof);      // global extraction and local storage for the solution
      solP2[i] = (*sol->_Sol[solP2Index])(solPDof);      // global extraction and local storage for the solution
      sysDof[dim * nDofsV + i ] = pdeSys->GetSystemDof(solP1Index, solP1PdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
      sysDof[dim * nDofsV + nDofsP + i ] = pdeSys->GetSystemDof(solP2Index, solP2PdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    for(unsigned i = 0; i < nDofsK1; i++) {
      unsigned iDof = msh->GetSolutionDof(i, iel, solPKCType);    // local to global mapping between solution node and solution dof
      solK1[i] = (*sol->_Sol[solK1Index])(iDof);
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // local to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }

    double h = 4. / 100;

    //std::vector <double> A = {0, 0, 0, 0, 1, -h/2.};

    //std::vector <double> A = {-1, 0, -1., 0, 0, 1.};

    std::vector <double> Ap;
    bool distributeCurvature = true;

    Data *data = new Data(solVType, solPKCType, solV, solP1, solP2, solK1, Res, Jac, coordX, ielGeom, rho1, rho2, mu1, mu2, sigma, dt, g, A, distributeCurvature);

    cad.SetDataPointer(data);

    std::vector<std::vector<double>> y;
    cad.GetXInParentElement(coordX, y);
    std::vector<std::vector<double>> yi = cad.GetXiInParentElement();

    cad.GetConicsInTargetElement(y, A, Ap);
    //std::cout<<Ap[0]<<" "<<Ap[1]<<" "<<Ap[2]<<" "<<Ap[3]<<" "<<Ap[4]<<" "<<Ap[5]<<std::endl;
    tuple <double, double, double> a = cad.Assemble(AssembleNavierStokes, 5, y, yi, Ap);

    double C = std::get<0>(a) / (std::get<0>(a) + std::get<1>(a));

    sol->_Sol[solC0Index]->set(iel, C);

    unsigned nDofs1 = msh->GetElementDofNumber(iel, solPKCType);
    for(unsigned i = 0; i < nDofs1; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, solPKCType);
      sol->_Sol[solC1Index]->add(idof, C);
      sol->_Sol[solCntIndex]->add(idof, 1);
    }

    RES->add_vector_blocked(Res, sysDof);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);
  } //end element loop for each process

  sol->_Sol[solC0Index]->close();
  sol->_Sol[solC1Index]->close();
  sol->_Sol[solCntIndex]->close();

  for(unsigned i = msh->_dofOffset[solPKCType][iproc]; i < msh->_dofOffset[solPKCType][iproc + 1]; i++) {
    double value = (*sol->_Sol[solC1Index])(i);
    double cnt = (*sol->_Sol[solCntIndex])(i);
    sol->_Sol[solC1Index]->set(i, value / cnt);
  }
  sol->_Sol[solC1Index]->close();

  std::cout << "Navier-Stokes Assembly time = " << static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC << std::endl;

  //AssembleAllPenalty(ml_prob);

  RES->close();
  KK->close();

}

//Attempting to create J by hand
void GetError(MultiLevelProblem & ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled


  //  extract pointers to the several objects that we are going to use
  LinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<LinearImplicitSystem> ("NS");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  std::vector <double> A = {1., 0, 1., 0, 0, -1.};

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
  unsigned solK1Index = mlSol->GetIndex("K1");

  unsigned solPKCType = mlSol->GetSolutionType(solP1Index);    // get the finite element type for "u"

  std::vector < std::vector < double > >  solV(dim);    // local solution
  std::vector < double >  solP1; // local solution
  std::vector < double >  solP2; // local solution

  std::vector < double >  solK1; // local solution

  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector< double > localErr(dim + 2, 0.);
  std::vector< double > globalErr(dim + 2, 0.);// local redidual std::vector
  std::vector < double > Jac;

  clock_t start_time = clock();

  ConicAdaptiveRefinement cad;

  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPKCType);    // number of solution element dofs

    unsigned nDofsVP = dim * nDofsV + 2 * nDofsP;

    unsigned nDofsK1 = msh->GetElementDofNumber(iel, solPKCType);

    // resize local arrays
    for(unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofsV);
      coordX[k].resize(nDofsV);
    }
    solP1.resize(nDofsP);
    solP2.resize(nDofsP);

    solK1.resize(nDofsK1);

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // local to global mapping between solution node and solution dof
      for(unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);
      }
    }

    for(unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPKCType);    // local to global mapping between solution node and solution dof
      solP1[i] = (*sol->_Sol[solP1Index])(solPDof);      // global extraction and local storage for the solution
      solP2[i] = (*sol->_Sol[solP2Index])(solPDof);      // global extraction and local storage for the solution
    }

    for(unsigned i = 0; i < nDofsK1; i++) {
      unsigned iDof = msh->GetSolutionDof(i, iel, solPKCType);    // local to global mapping between solution node and solution dof
      solK1[i] = (*sol->_Sol[solK1Index])(iDof);
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // local to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }

    bool distributeCurvature = true;
    Data *data = new Data(solVType, solPKCType, solV, solP1, solP2, solK1, localErr, Jac, coordX, ielGeom, rho1, rho2, mu1, mu2, sigma, dt, g, A, distributeCurvature);

    //std::cout<<localErr[0] << " "<<localErr[1] << " "<<localErr[2] << " "<<localErr[3] << std::endl;

    cad.SetDataPointer(data);

    std::vector<std::vector<double>> y;
    cad.GetXInParentElement(coordX, y);
    std::vector<std::vector<double>> yi = cad.GetXiInParentElement();

    std::vector <double> Ap;
    cad.GetConicsInTargetElement(y, A, Ap);
    cad.Assemble(AssembleError, 5, y, yi, Ap);

  } //end element loop for each process

  MPI_Allreduce(localErr.data(), globalErr.data(), localErr.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


  std::cout.precision(14);
  std::cout << "U-errNorm = " << sqrt(globalErr[0]) << std::endl;
  std::cout << "V-errNorm = " << sqrt(globalErr[1]) << std::endl;
  std::cout << "P1-errNorm = " << sqrt(globalErr[2]) << std::endl;
  std::cout << "P2-errNorm = " << sqrt(globalErr[3]) << std::endl;

  std::cout << "Error Function time = " << static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC << std::endl;


}


void AssembleError(Data *data, const std::vector <double> &phiV, const std::vector <double> &phiV_x, const std::vector <double> &phiP,
                   const double &C, const double &weight, const double &weight1, const double &weight2, const double &weightI,
                   const std::vector <double> &N, const std::vector<double> &kappa, const double &dsN, const double &eps) {


  const unsigned &dim = data->_V.size();
  const unsigned &nDofsV = phiV.size();
  const unsigned &nDofsP = phiP.size();
  unsigned nDofsVP = dim * nDofsV + 2 * nDofsP;

  std::vector < double > solVg(dim, 0);
  std::vector < std::vector < double > > SolVg_x(dim, std::vector<double> (dim, 0));
  for(unsigned i = 0; i < nDofsV; i++) {
    for(unsigned  K = 0; K < dim; K++) {
      solVg[K] += data->_V[K][i] * phiV[i];
      for(unsigned J = 0; J < dim; J++) {
        SolVg_x[K][J] += data->_V[K][i] * phiV_x[i * dim + J];
      }
    }
  }

  double solP1g = 0;
  double solP2g = 0;

  for(unsigned i = 0; i < nDofsP; i++) {
    solP1g += phiP[i] * data->_P1[i];
    solP2g += phiP[i] * data->_P2[i];
  }

 // if(fabs(solP1g)>0.0000001 )
  //  std::cout << solP1g <<" "<< solP2g << std::endl;

  std::vector < double > solVg_exc(dim, 0);
  double solP1g_exc = 1;
  double solP2g_exc = 0;

  double rho = data->_rho1 * weight1 + data->_rho2 * weight2;
  double mu = data->_mu1 * weight1 + data->_mu2 * weight2;
  double rhoC = rho;

  for(unsigned K = 0;K < dim; K++){
    data->_res[K] += (solVg[K] - solVg_exc[K]) * (solVg[K] - solVg_exc[K]) * weight;
  }
  data->_res[dim + 0] += (solP1g - solP1g_exc) * (solP1g - solP1g_exc) * weight * weight1;
  data->_res[dim + 1] += (solP2g - solP2g_exc) * (solP2g - solP2g_exc) * weight * weight2;

}



