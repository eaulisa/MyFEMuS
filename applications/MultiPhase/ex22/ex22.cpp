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


#include "../ex21/ConicAdaptiveRefinement.hpp"
// #include "../include/MyMarker/MyMarker.hpp"
// #include "../include/MyMarker/MyMarker.cpp"
//
// typedef double TypeIO;
// typedef cpp_bin_float_oct TypeA;
// typedef cpp_bin_float_oct oct;
//
// // CutFemWeight <double, double> quad = CutFemWeight<double, double>(QUAD, 5, "legendre");
// CutFemWeight <TypeIO, TypeA> quad  = CutFemWeight<TypeIO, TypeA >(QUAD, 5, "legendre");
// CutFemWeight <TypeIO, TypeA> tri  = CutFemWeight<TypeIO, TypeA >(TRI, 1, "legendre");
// Fem fem = Fem(quad.GetGaussQuadratureOrder(), quad.GetDimension());
//
// #include "../include/Cloud.hpp"
// Cloud *cld;
// Cloud *cldint;

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
const double sigma = 1.96;
const double gravity = -0.;
const double dt = 0.1;


std::vector <double> g={0,-1,0};

#include "./include/GhostPenalty.hpp"
#include "./include/GhostPenaltyDGP.hpp"
#include "./include/Stabilization.hpp"

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
  else if(!strcmp(SolName, "P1") || !strcmp(SolName, "P2")) {
    dirichlet = false;
    value = 0.;
  }

  return dirichlet;
}

//double TimeStepMultiphase(const double time);
void AssembleMultiphase(MultiLevelProblem& ml_prob);
//void AssembleMultiphaseAD(MultiLevelProblem& ml_prob);

//void AssembleBoussinesqAppoximation_AD(MultiLevelProblem& ml_prob);    //, unsigned level, const unsigned &levelMax, const bool &assembleMatrix );

//void GetNormalQuad(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double &R, std::vector<double> &a, double &d,  std::vector<double> &xm, std::vector<double> &b, double &db, unsigned &cut);

//void GetNormalTri(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double &R, std::vector<double> &a, double &d,  std::vector<double> &xm, std::vector<double> &b, double &db, unsigned &cut);

//void GetNormalTet(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double &R, std::vector<double> &a, double &d,  std::vector<double> &xm, std::vector<double> &b, double &db, double &vol, unsigned &cut);
//double getHeightPolyhedronSphereInt(const std::vector < std::vector <double> > &aN, const std::vector <double> &a, const std::vector <double> &xg, const double &R);
//double CurvatureQuadric(const std::vector<double> &a, const std::vector<double> &xp);
//void NormalQuadric(const std::vector<double> &a, const std::vector<double> &xp, std::vector<double> &N);

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
//   mlMsh.ReadCoarseMesh("./input/cube_hex.neu", "seventh", scalingFactor);
//   mlMsh.ReadCoarseMesh("./input/square_quad.neu", "fifth", scalingFactor);
  mlMsh.GenerateCoarseBoxMesh(100, 100, 0, -2., 2., -2., 2., 0., 0., QUAD9, "fifth"); // Turek 1&2
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

  // add variables to mlSol
  mlSol.AddSolution("U", LAGRANGE, SECOND);
  mlSol.AddSolution("V", LAGRANGE, SECOND);
  if(dim == 3) mlSol.AddSolution("W", LAGRANGE, SECOND);
  //mlSol.AddSolution("P1",  DISCONTINUOUS_POLYNOMIAL, ZERO);
  //mlSol.AddSolution("P2",  DISCONTINUOUS_POLYNOMIAL, ZERO);
  mlSol.AddSolution("P1", LAGRANGE, FIRST);
  mlSol.AddSolution("P2", LAGRANGE, FIRST);

  mlSol.AddSolution("C", DISCONTINUOUS_POLYNOMIAL, ZERO, 1, false);

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
//   mlSol.FixSolutionAtOnePoint("P2");

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
  system.MGsolve();
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 1);



  return 0;
}


//Attempting to create J by hand
void AssembleMultiphase(MultiLevelProblem& ml_prob) {
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
  unsigned solPType = mlSol->GetSolutionType(solP1Index);    // get the finite element type for "u"

  unsigned solCIndex = mlSol->GetIndex("C");

  std::vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("U");    // get the position of "U" in the pdeSys object
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("V");    // get the position of "V" in the pdeSys object
  if(dim == 3) solVPdeIndex[2] = mlPdeSys->GetSolPdeIndex("W");

  unsigned solP1PdeIndex = mlPdeSys->GetSolPdeIndex("P1");    // get the position of "P" in the pdeSys object
  unsigned solP2PdeIndex = mlPdeSys->GetSolPdeIndex("P2");    // get the position of "P" in the pdeSys object

  std::vector < std::vector < double > >  solV(dim);    // local solution
  std::vector < double >  solP1; // local solution
  std::vector < double >  solP2; // local solution

  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector< unsigned > sysDof; // local to global pdeSys dofs
  std::vector< double > Res; // local redidual std::vector
  std::vector < double > Jac;

  clock_t start_time = clock();

  ConicAdaptiveRefinement cad;

  sol->_Sol[solCIndex]->zero();

  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);    // number of solution element dofs

    unsigned nDofsVP = dim * nDofsV + 2 * nDofsP;

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


    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // local to global mapping between solution node and solution dof

      for(unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);
        sysDof[k * nDofsV + i] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }

    for(unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);    // local to global mapping between solution node and solution dof
      solP1[i] = (*sol->_Sol[solP1Index])(solPDof);      // global extraction and local storage for the solution
      solP2[i] = (*sol->_Sol[solP2Index])(solPDof);      // global extraction and local storage for the solution
      sysDof[dim * nDofsV + i ] = pdeSys->GetSystemDof(solP1Index, solP1PdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
      sysDof[dim * nDofsV + nDofsP + i ] = pdeSys->GetSystemDof(solP2Index, solP2PdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // local to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }

    double h = 4./100;

    //std::vector <double> A = {0, 0, 0, 0, 1, -h/2.};

    std::vector <double> A = {-1, 0, -1., 0, 0, 10.};
    std::vector <double> Ap;

    Data *data = new Data(solVType, solPType, solV, solP1, solP2, Res, Jac, coordX, ielGeom, rho1, rho2, mu1, mu2, sigma, dt, g, A);

    cad.SetDataPointer(data);

    std::vector<std::vector<double>> y;
    cad.GetXInParentElement(coordX, y);
    std::vector<std::vector<double>> yi = cad.GetXiInParentElement();

    cad.GetConicsInTargetElement(y, A, Ap);
    //std::cout<<Ap[0]<<" "<<Ap[1]<<" "<<Ap[2]<<" "<<Ap[3]<<" "<<Ap[4]<<" "<<Ap[5]<<std::endl;
    tuple <double, double, double> a = cad.AdaptiveRefinement(2, y, yi, Ap);

    double C = std::get<0>(a) / (std::get<0>(a) + std::get<1>(a));

    sol->_Sol[solCIndex]->set(iel,C);

    RES->add_vector_blocked(Res, sysDof);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);


  } //end element loop for each process

  sol->_Sol[solCIndex]->close();


  std::cout << "Navier-Stokes Assembly time = " << static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC << std::endl;


  //AssembleStabilizationTerms(ml_prob);
  //AssembleGhostPenalty(ml_prob);
  //AssembleGhostPenaltyDGP(ml_prob, true);
  //AssembleGhostPenaltyDGP(ml_prob, false);

  //AssembleCIPPressure(ml_prob, true);
  //AssembleCIPPressure(ml_prob, false);


  RES->close();
  KK->close();



}


