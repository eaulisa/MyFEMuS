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

#include "TransientSystem.hpp"
#include "NonLinearImplicitSystem.hpp"

#include "adept.h"

#include "PolynomialBases.hpp"

#include "CutFemWeight.hpp"

#include "CDWeights.hpp"


#include "petsc.h"
#include "petscmat.h"
#include "PetscMatrix.hpp"
#include "Fem.hpp"

#include "../include/MyMarker/MyMarker.hpp"


typedef double TypeIO;
typedef cpp_bin_float_oct TypeA;
typedef cpp_bin_float_oct oct;

// CutFemWeight <double, double> quad = CutFemWeight<double, double>(QUAD, 5, "legendre");
CutFemWeight <TypeIO, TypeA> quad  = CutFemWeight<TypeIO, TypeA >(QUAD, 5, "legendre");
CutFemWeight <TypeIO, TypeA> tri  = CutFemWeight<TypeIO, TypeA >(TRI, 1, "legendre");
Fem fem = Fem(quad.GetGaussQuadratureOrder(), quad.GetDimension());

#include "../include/Cloud.hpp"
#include "../include/AdaptiveSplit.hpp"

Cloud *cld;
Cloud *cldint;
AdaptiveSplit *asplit;


// // RT
// const double mu2 = 0.009535;
// const double mu1 = 0.006349;
// const double rho2 = 1.5;
// const double rho1 = 1.;
// const double sigma = 0.;
// const double gravity = -10.;

// // Turek 1
const double mu1 = 0.1;
const double mu2 = 10.;
const double rho1 = 1.;
const double rho2 = 1000;
const double sigma = 1.96;
const double gravity = -0.98;

std::vector <double> g;

#include "../include/GhostPenalty.hpp"
#include "../include/GhostPenaltyDGP.hpp"
#include "../include/Stabilization.hpp"

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
    if(facename == 2 || facename == 4) dirichlet = false;
    value = 0.;
//     if(x[0] < 0. && x[1] < 0.5 && x[1] > -0.5 && x[2] < 0.5 && x[2] > -0.5) value = 1.;
  }
  else if(!strcmp(SolName, "W")) {
    value = 0.;
  }
  else if(!strcmp(SolName, "P")) {
    dirichlet = false;
    value = 0.;
  }

  return dirichlet;
}

double TimeStepMultiphase(const double time);
void AssembleMultiphase(MultiLevelProblem& ml_prob);


int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  //mlMsh.ReadCoarseMesh("./input/cube_hex.neu", "seventh", scalingFactor);
//   mlMsh.ReadCoarseMesh("./input/square_quad.neu", "fifth", scalingFactor);
  mlMsh.GenerateCoarseBoxMesh(161, 321, 0, 0., 1., 0., 2., 0., 0., QUAD9, "fifth"); //turek
//   mlMsh.GenerateCoarseBoxMesh(64, 256, 0, -0.5, 0.5, -2, 2, 0., 0., QUAD9, "fifth"); //Raileigh Taylor
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
  mlSol.AddSolution("U", LAGRANGE, FIRST, 2);
  mlSol.AddSolution("V", LAGRANGE, FIRST, 2);
  if(dim == 3) mlSol.AddSolution("W", LAGRANGE, SECOND, 2);
  mlSol.AddSolution("P1",  DISCONTINUOUS_POLYNOMIAL, ZERO);
  mlSol.AddSolution("P2",  DISCONTINUOUS_POLYNOMIAL, ZERO);
//   mlSol.AddSolution("P1", LAGRANGE, FIRST);
//   mlSol.AddSolution("P2", LAGRANGE, FIRST);

  mlSol.AddSolution("C", DISCONTINUOUS_POLYNOMIAL, ZERO, false);
  mlSol.AddSolution("Cn", LAGRANGE, SECOND, false);
  mlSol.AddSolution("Q", DISCONTINUOUS_POLYNOMIAL, ZERO, false);
  mlSol.AddSolution("DIC", DISCONTINUOUS_POLYNOMIAL, ZERO, false);

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

  cld = new Cloud(sol);
  cldint = new Cloud(sol);

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
//   mlSol.FixSolutionAtOnePoint("P1");
//   mlSol.FixSolutionAtOnePoint("P2");

  mlSol.GenerateBdc("All");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);

  // add system Poisson in mlProb as a Linear Implicit System
  TransientNonlinearImplicitSystem& system = mlProb.add_system < TransientNonlinearImplicitSystem > ("NS");

  // add solution "u" to system
  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");

  if(dim == 3) system.AddSolutionToSystemPDE("W");

  system.AddSolutionToSystemPDE("P1");
  system.AddSolutionToSystemPDE("P2");

  system.SetSparsityPatternMinimumSize(250);

  // attach the assembling function to system
  system.SetAssembleFunction(AssembleMultiphase);
//   system.SetAssembleFunction(AssembleMultiphaseAD);
  system.AttachGetTimeIntervalFunction(TimeStepMultiphase);

  // initilaize and solve the system
  system.init();
  system.SetOuterSolver(PREONLY);


  // BEGIN Testing the class Cloud

  std::vector<std::string> velocity = {"U", "V"};
  std::cout << "Testing the class Cloud \n";


  unsigned nIterations = 1000;
  unsigned nMrk = 1000;


// cld->AddEllipse({XG, YG}, {RADIUS, RADIUS}, 8);
  cld->AddEllipse({XG, YG}, {RADIUS, RADIUS}, nMrk);
//   cld->AddQuadric({1.,0.,1.,-2.*XG ,-2*YG ,XG*XG+YG*YG-RADIUS*RADIUS}, 8);
//   cld->AddQuadric({1./50 ,0.,0.,0.,1., 1./256 - 0.005}, 8);
  cld->ComputeQuadraticBestFit();

  cldint->AddInteriorEllipse({XG, YG}, {RADIUS, RADIUS});
//   cldint->AddInteriorQuadric({1.,0.,1.,-2.*XG ,-2*YG ,XG*XG+YG*YG-RADIUS*RADIUS});
//   cldint->AddInteriorQuadric({1./50 ,0.,0.,0.,1., 1./256 - 0.005});
  cldint->RebuildInteriorMarkers(*cld, "C", "Cn");

  cld->PrintCSV("markerBefore", 0);
  cld->PrintCSV("marker", 0);
//   cldint->PrintCSV("markerInt", 0);


  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");
  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(true);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);

  unsigned qorder = 5;
  unsigned lmax = 1;
  asplit = new AdaptiveSplit(qorder, lmax);

  for(unsigned it = 1; it <= nIterations; it++) {

    mlSol.CopySolutionToOldSolution();

    system.MGsolve();
    double dt = system.GetIntervalTime();

    vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, it);

    cld->RKAdvection(4, velocity, dt);
    cldint->RKAdvection(4, velocity, dt);
    cld->PrintCSV("markerBefore", it);
    cld->ComputeQuadraticBestFit();
    cld->RebuildMarkers(11, 9, 10);
    cldint->RebuildInteriorMarkers(*cld, "C", "Cn");
    cld->PrintCSV("marker", it);
    //vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, it);

  }

  delete asplit;
  delete cld;
  delete cldint;

  return 0;
}

double TimeStepMultiphase(const double time) {
//   double dt =  0.005; //RT
  double dt =  0.01; //Turek
  return dt;
}

//Attempting to create J by hand
void AssembleMultiphase(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled


  //  extract pointers to the several objects that we are going to use
  TransientNonlinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<TransientNonlinearImplicitSystem> ("NS");   // pointer to the linear implicit system named "Poisson"
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

  if(dim == 2) g = {0, gravity};
  else g = {0, 0, gravity};

  AssembleGhostPenalty(ml_prob);
  AssembleGhostPenaltyDGP(ml_prob, true);
  AssembleGhostPenaltyDGP(ml_prob, false);
  AssembleStabilizationTerms(ml_prob);


  double dt =  mlPdeSys->GetIntervalTime();

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
  std::vector < std::vector < double > >  solVOld(dim);    // local solution
  std::vector < double >  solP1; // local solution
  std::vector < double >  solP2; // local solution

  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector <double> phiV;  // local test function for velocity
  std::vector <double> phiV_x; // local test function first order partial derivatives

  double* phiP; // local test function for the pressure
  double weight; // gauss point weight

  std::vector< unsigned > sysDof; // local to global pdeSys dofs
  std::vector< double > Res; // local redidual std::vector
  std::vector < double > Jac;

  /* BEGIN cutfem stuff for surface tension integration */

  double R = RADIUS;

  std::vector < std::vector < double > > x1;
  std::vector < double > xg(dim);
  xg[0] = XG;
  xg[1] = YG;
  if(dim > 2) xg[2] = ZG;

  unsigned qM = 3;
  double dx = .05;
  double dtetha = 2.;

  double eps = 0.00000001;

  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    double C = (*sol->_Sol[solCIndex])(iel);

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDof = msh->GetElementDofNumber(iel, 0);  // number of coordinate linear element dofs
    x1.resize(dim);
    for(unsigned k = 0; k < dim; k++) {
      x1[k].resize(nDof);
    }

    for(unsigned k = 0; k < dim; k++) {
      for(unsigned i = 0; i < nDof; i++) {
        unsigned xDof  = msh->GetSolutionDof(i, iel, coordXType);    // global to global mapping between coordinates node and coordinate dof
        x1[k][(i + 2) % nDof] = (*msh->_topology->_Sol[k])(xDof); // global extraction and local storage for the element coordinates
      }
    }

    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);    // number of solution element dofs

    unsigned nDofsVP = dim * nDofsV + 2 * nDofsP;

    // resize local arrays
    sysDof.resize(nDofsVP);
    Res.assign(nDofsVP, 0.);
    Jac.assign(nDofsVP * nDofsVP, 0.);

    for(unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofsV);
      solVOld[k].resize(nDofsV);
      coordX[k].resize(nDofsV);
    }
    solP1.resize(nDofsP);
    solP2.resize(nDofsP);


    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // local to global mapping between solution node and solution dof

      for(unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);
        solVOld[k][i] = (*sol->_SolOld[solVIndex[k]])(solVDof);
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


    unsigned cut = (cld->GetNumberOfMarker(iel) > 0) ? 1 : 0;

    if(cut == 1) {

      std::vector<std::vector<double>> Jacob, JacI;
      std::vector<double> a;
      std::vector<double> xm;
      double d;

      const elem_type *femV = asplit->GetFiniteElementFine(ielGeom, solVType);
      const elem_type *femP = asplit->GetFiniteElementFine(ielGeom, solPType);

      femV->GetJacobianMatrix(coordX, {0., 0.}, weight, Jacob, JacI); //TODO in all ex
      cld->GetElementQuantities(iel, Jacob, asplit->GetXiFather(), asplit->GetDsFather(),  asplit->GetNiFather());

      asplit->Split(coordX, ielGeom, true);

      const std::vector <double> &weight1 = asplit->GetWeight1();
      const std::vector <double> &weight2 = asplit->GetWeight2();
      const std::vector <double> &weightI = asplit->GetWeightI();
      const std::vector<std::vector <double>> &xg = asplit->GetXg();
      const std::vector<std::vector <double>> &xig = asplit->GetXig();

      C = asplit->GetC();
      sol->_Sol[solCIndex]->set(iel, C);

      std::vector<double> NN(dim, 0.);
      double kk = 0.;

      // *** Gauss point loop ***
      for(unsigned ig = 0; ig < weight1.size(); ig++) {
        // *** get gauss point weight, test function and test function partial derivatives ***
        femV->Jacobian(coordX, xig[ig], weight, phiV, phiV_x);
        weight = weight1[ig] + weight2[ig];

        std::vector<double> phiP;
        femP->GetPhi(phiP, xig[ig]);

        kk = cld->GetAverageCurvature(iel);
        NN = cld->GetNormal(iel, xg[ig]);

        std::vector < double > solV_gss(dim, 0);
        std::vector < double > solVOld_gss(dim, 0);
        std::vector < std::vector < double > > gradSolV_gss(dim);

        for(unsigned  k = 0; k < dim; k++) {
          gradSolV_gss[k].assign(dim, 0.);
        }

        for(unsigned i = 0; i < nDofsV; i++) {
          for(unsigned  k = 0; k < dim; k++) {
            solV_gss[k] += solV[k][i] * phiV[i];
            solVOld_gss[k] += solVOld[k][i] * phiV[i];
          }
          for(unsigned j = 0; j < dim; j++) {
            for(unsigned  k = 0; k < dim; k++) {
              gradSolV_gss[k][j] += solV[k][i] * phiV_x[i * dim + j];
            }
          }
        }

        double solP1_gss = 0;
        double solP2_gss = 0;
        for(unsigned i = 0; i < nDofsP; i++) {
          solP1_gss += phiP[i] * solP1[i];
          solP2_gss += phiP[i] * solP2[i];
        }
        double rho = rho1 * C + rho2 * (1. - C);
        double mu = mu1 * C + mu2 * (1. - C);
        double rhoC = rho1 * C + rho2 * (1. - C);

        // *** phiV_i loop ***
        for(unsigned i = 0; i < nDofsV; i++) {
          for(unsigned  I = 0; I < dim; I++) {  //momentum equation in k
            double NSV = 0.;
            for(unsigned J = 0; J < dim; J++) {  // second index j in each equation
              NSV   +=  mu * phiV_x[i * dim + J] * (gradSolV_gss[I][J] + gradSolV_gss[J][I]); // diffusion
              NSV   +=  rho * phiV[i] * (solV_gss[J] * gradSolV_gss[I][J]); // nonlinear term
            }

            NSV += rho * phiV[i] * (solV_gss[I] - solVOld_gss[I]) / dt ;
            NSV += - rhoC * phiV[i] * g[I]; // gravity term
            NSV *= weight;
            NSV += - phiV_x[i * dim + I] * (solP1_gss * weight1[ig] + solP2_gss * weight2[ig]);  // pressure gradient
            Res[I * nDofsV + i] -=  NSV;

            Res[I * nDofsV + i] += - sigma * phiV[i] * NN[I] * weightI[ig] * kk;

          }
        } // end phiV_i loop

        // *** phiP_i loop ***
        for(unsigned i = 0; i < nDofsP; i++) {
          for(int I = 0; I < dim; I++) {
            Res[dim * nDofsV + i] += - gradSolV_gss[I][I] * phiP[i]  * weight1[ig]; //continuity
            Res[dim * nDofsV + nDofsP + i] += - gradSolV_gss[I][I] * phiP[i]  * weight2[ig]; //continuity
          }
        } // end phiP_i loop
        // end gauss point loop


        //--------------------------------------------------------------------------------------------------------
        // Add the local Matrix/Vector into the global Matrix/Vector

        for(unsigned i = 0; i < nDofsV; i++) {
          for(unsigned I = 0; I < dim; I++) { //row velocity blocks or dimension
            unsigned VIrow = I * nDofsV + i;
            for(unsigned j = 0; j < nDofsV; j++) {
              unsigned VIcolumn = I * nDofsV + j;

              Jac[ VIrow * nDofsVP + VIcolumn] += rho * phiV[i] * phiV[j] * weight / dt; // inertia

              for(unsigned J = 0; J < dim ; J++) { //column velocity blocks or dimension
                unsigned VJcolumn = J * nDofsV + j;
                Jac[ VIrow * nDofsVP + VIcolumn ] += mu * phiV_x[i * dim + J] * phiV_x[j * dim + J] * weight; //diagonal diffusion
                Jac[ VIrow * nDofsVP + VJcolumn ] += mu * phiV_x[i * dim + J] * phiV_x[j * dim + I] * weight; //off-diagonal diffusion

                Jac[ VIrow * nDofsVP + VIcolumn ] += rho * phiV[i] * solV_gss[J] * phiV_x[j * dim + J] * weight; //diagonal nonlinear
                Jac[ VIrow * nDofsVP + VJcolumn ] += rho * phiV[i] * phiV[j] * gradSolV_gss[I][J] * weight; //off-diagonal nonlinear
              }
            }

            for(unsigned j = 0; j < nDofsP; j++) {
              unsigned P1column = dim * nDofsV + j;
              unsigned P2column = dim * nDofsV + nDofsP + j;
              Jac[VIrow * nDofsVP + P1column] += - phiV_x[i * dim + I] * phiP[j] * weight1[ig]; //pressure gradient
              Jac[P1column * nDofsVP + VIrow] -= - phiV_x[i * dim + I] * phiP[j] * weight1[ig]; //continuity
              Jac[VIrow * nDofsVP + P2column] += - phiV_x[i * dim + I] * phiP[j] * weight2[ig]; //pressure gradient
              Jac[P2column * nDofsVP + VIrow] -= - phiV_x[i * dim + I] * phiP[j] * weight2[ig]; //continuity
            }
          }
        }
      }
    }
    else {

      const elem_type *femV = asplit->GetFiniteElementCoarse(ielGeom, solVType);
      const elem_type *femP = asplit->GetFiniteElementCoarse(ielGeom, solPType);

      // *** Gauss point loop ***
      for(unsigned ig = 0; ig < femV->GetGaussPointNumber(); ig++) {
        // *** get gauss point weight, test function and test function partial derivatives ***
        femV->Jacobian(coordX, ig, weight, phiV, phiV_x);
        phiP = femP->GetPhi(ig);
     
        std::vector < double > solV_gss(dim, 0);
        std::vector < double > solVOld_gss(dim, 0);
        std::vector < std::vector < double > > gradSolV_gss(dim);

        for(unsigned  k = 0; k < dim; k++) {
          gradSolV_gss[k].assign(dim, 0.);
        }

        for(unsigned i = 0; i < nDofsV; i++) {
          for(unsigned  k = 0; k < dim; k++) {
            solV_gss[k] += solV[k][i] * phiV[i];
            solVOld_gss[k] += solVOld[k][i] * phiV[i];
          }
          for(unsigned j = 0; j < dim; j++) {
            for(unsigned  k = 0; k < dim; k++) {
              gradSolV_gss[k][j] += solV[k][i] * phiV_x[i * dim + j];
            }
          }
        }

        double solP1_gss = 0;
        double solP2_gss = 0;
        for(unsigned i = 0; i < nDofsP; i++) {
          solP1_gss += phiP[i] * solP1[i];
          solP2_gss += phiP[i] * solP2[i];
        }

        double rho = rho1 * C + rho2 * (1. - C);
        double mu = mu1 * C + mu2 * (1. - C);
        double rhoC = rho1 * C + rho2 * (1. - C);

        double weight1 = weight * C;
        double weight2 = weight * (1. - C);

        // *** phiV_i loop ***
        for(unsigned i = 0; i < nDofsV; i++) {
          for(unsigned  I = 0; I < dim; I++) {  //momentum equation in k
            double NSV = 0.;
            for(unsigned J = 0; J < dim; J++) {  // second index j in each equation
              NSV   +=  mu * phiV_x[i * dim + J] * (gradSolV_gss[I][J] + gradSolV_gss[J][I]); // diffusion
              NSV   +=  rho * phiV[i] * (solV_gss[J] * gradSolV_gss[I][J]); // nonlinear term
            }
            NSV += rho * phiV[i] * (solV_gss[I] - solVOld_gss[I]) / dt ;
            NSV += - rhoC * phiV[i] * g[I]; // gravity term
            NSV *= weight;
            NSV += - phiV_x[i * dim + I] * (solP1_gss * weight1 + solP2_gss * weight2);  // pressure gradient
            Res[I * nDofsV + i] -=  NSV;
          }
        } // end phiV_i loop

        // *** phiP_i loop ***
        for(unsigned i = 0; i < nDofsP; i++) {
          for(int I = 0; I < dim; I++) {
            Res[dim * nDofsV + i] += - gradSolV_gss[I][I] * phiP[i]  * weight1; //continuity
            Res[dim * nDofsV + nDofsP + i] += - gradSolV_gss[I][I] * phiP[i]  * weight2; //continuity
          }
          if(C == 0)
            Res[dim * nDofsV + i] += - solP1_gss * phiP[i]  * weight2 * eps; //penalty
          if(C == 1)
            Res[dim * nDofsV + nDofsP + i] += - solP2_gss * phiP[i]  * weight1 * eps; //penalty

        } // end phiP_i loop
        // end gauss point loop


        //--------------------------------------------------------------------------------------------------------
        // Add the local Matrix/Vector into the global Matrix/Vector

        for(unsigned i = 0; i < nDofsV; i++) {
          for(unsigned I = 0; I < dim; I++) { //row velocity blocks or dimension
            unsigned VIrow = I * nDofsV + i;
            for(unsigned j = 0; j < nDofsV; j++) {
              unsigned VIcolumn = I * nDofsV + j;

              Jac[ VIrow * nDofsVP + VIcolumn] += rho * phiV[i] * phiV[j] * weight / dt; // inertia

              for(unsigned J = 0; J < dim ; J++) { //column velocity blocks or dimension
                unsigned VJcolumn = J * nDofsV + j;
                Jac[ VIrow * nDofsVP + VIcolumn ] += mu * phiV_x[i * dim + J] * phiV_x[j * dim + J] * weight; //diagonal diffusion
                Jac[ VIrow * nDofsVP + VJcolumn ] += mu * phiV_x[i * dim + J] * phiV_x[j * dim + I] * weight; //off-diagonal diffusion

                Jac[ VIrow * nDofsVP + VIcolumn ] += rho * phiV[i] * solV_gss[J] * phiV_x[j * dim + J] * weight; //diagonal nonlinear
                Jac[ VIrow * nDofsVP + VJcolumn ] += rho * phiV[i] * phiV[j] * gradSolV_gss[I][J] * weight; //off-diagonal nonlinear
              }
            }

            for(unsigned j = 0; j < nDofsP; j++) {
              unsigned P1column = dim * nDofsV + j;
              unsigned P2column = dim * nDofsV + nDofsP + j;
              Jac[VIrow * nDofsVP + P1column] += - phiV_x[i * dim + I] * phiP[j] * weight1; //pressure gradient
              Jac[P1column * nDofsVP + VIrow] -= - phiV_x[i * dim + I] * phiP[j] * weight1; //continuity
              Jac[VIrow * nDofsVP + P2column] += - phiV_x[i * dim + I] * phiP[j] * weight2; //pressure gradient
              Jac[P2column * nDofsVP + VIrow] -= - phiV_x[i * dim + I] * phiP[j] * weight2; //continuity
            }
          }
        }
        for(unsigned i = 0; i < nDofsP; i++) {
          unsigned P1row = dim * nDofsV + i;
          unsigned P2row = dim * nDofsV + nDofsP + i;
          for(unsigned j = 0; j < nDofsP; j++) {
            unsigned P1column = dim * nDofsV + j;
            unsigned P2column = dim * nDofsV + nDofsP + j;
            if(C == 0)
              Jac[P1row * nDofsVP + P1column] += phiP[i] * phiP[j] * weight2 * eps; //continuity
            if(C == 1)
              Jac[P2row * nDofsVP + P2column] += phiP[i] * phiP[j] * weight1 * eps; //continuity
          }
        }
      }
    }

    RES->add_vector_blocked(Res, sysDof);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);


  } //end element loop for each process


  RES->close();
  KK->close();

//  VecView ( (static_cast<PetscVector*> (RES))->vec(),  PETSC_VIEWER_STDOUT_SELF);
//MatView ( (static_cast<PetscMatrix*> (KK))->mat(), PETSC_VIEWER_STDOUT_SELF);

//   PetscViewer    viewer;
//   PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL, NULL, 0, 0, 900, 900, &viewer);
//   PetscObjectSetName((PetscObject) viewer, "PWilmore matrix");
//   PetscViewerPushFormat(viewer, PETSC_VIEWER_DRAW_LG);
//   MatView((static_cast<PetscMatrix*>(KK))->mat(), viewer);
//   double a;
//   std::cin >> a;


}


