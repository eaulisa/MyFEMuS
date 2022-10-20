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

#include "BestFitPlane.hpp"
#include "Fem.hpp"
#include "GenerateTriangles.hpp"

#include "../include/MyMarker/MyMarker.hpp"
#include "../include/MyMarker/MyMarker.cpp"

typedef double TypeIO;
typedef cpp_bin_float_oct TypeA;
typedef cpp_bin_float_oct oct;

// CutFemWeight <double, double> quad = CutFemWeight<double, double>(QUAD, 5, "legendre");
CutFemWeight <TypeIO, TypeA> quad  = CutFemWeight<TypeIO, TypeA >(QUAD, 5, "legendre");
CutFemWeight <TypeIO, TypeA> tri  = CutFemWeight<TypeIO, TypeA >(TRI, 1, "legendre");
Fem fem = Fem(quad.GetGaussQuadratureOrder(), quad.GetDimension());

#include "../include/Cloud.hpp"
Cloud *cld;
Cloud *cldint;

const double mu1 = 1.;
const double mu2 = 10.;
const double rho1 = 100.;
const double rho2 = 1000.;
const double sigma = 24.5;

#include "../include/GhostPenalty.hpp"
#include "../include/GhostPenaltyDGP.hpp"

#define RADIUS 0.25
#define XG 0.5
#define YG 0.5
#define ZG 0.
bool gravity = true;

using namespace femus;


bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet

  if(!strcmp(SolName, "U")) {  // strcmp compares two string in lexiographic sense.
    value = 0.;
  }
  else if(!strcmp(SolName, "V")) {
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

void AssembleBoussinesqAppoximation_AD(MultiLevelProblem& ml_prob);    //, unsigned level, const unsigned &levelMax, const bool &assembleMatrix );

void GetNormalQuad(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double &R, std::vector<double> &a, double &d,  std::vector<double> &xm, std::vector<double> &b, double &db, unsigned &cut);

void GetNormalTri(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double &R, std::vector<double> &a, double &d,  std::vector<double> &xm, std::vector<double> &b, double &db, unsigned &cut);

void GetNormalTet(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double &R, std::vector<double> &a, double &d,  std::vector<double> &xm, std::vector<double> &b, double &db, double &vol, unsigned &cut);
double getHeightPolyhedronSphereInt(const std::vector < std::vector <double> > &aN, const std::vector <double> &a, const std::vector <double> &xg, const double &R);
double CurvatureQuadric(const std::vector<double> &a, const std::vector<double> &xp);
void NormalQuadric(const std::vector<double> &a, const std::vector<double> &xp, std::vector<double> &N);

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  //mlMsh.ReadCoarseMesh("./input/cube_hex.neu", "seventh", scalingFactor);
//   mlMsh.ReadCoarseMesh("./input/square_quad.neu", "fifth", scalingFactor);
  mlMsh.GenerateCoarseBoxMesh(4, 8, 0, 0, 1, 0, 2, 0., 0., QUAD9, "seventh"); 
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
     probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();

  unsigned numberOfUniformLevels = 6;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels, numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // erase all the coarse mesh levels
  mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("U", LAGRANGE, SECOND, 2);
  mlSol.AddSolution("V", LAGRANGE, SECOND, 2);
  if(dim == 3) mlSol.AddSolution("W", LAGRANGE, SECOND, 2);
  mlSol.AddSolution("P1",  DISCONTINUOUS_POLYNOMIAL, ZERO);
  mlSol.AddSolution("P2",  DISCONTINUOUS_POLYNOMIAL, ZERO);
//   mlSol.AddSolution("P1", LAGRANGE, FIRST);
//   mlSol.AddSolution("P2", LAGRANGE, FIRST);

  mlSol.AddSolution("C", DISCONTINUOUS_POLYNOMIAL, ZERO, false);
  mlSol.AddSolution("Cn", LAGRANGE, SECOND, false);

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
  mlSol.FixSolutionAtOnePoint("P2");

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
  system.AttachGetTimeIntervalFunction(TimeStepMultiphase);

  // initilaize and solve the system
  system.init();
  system.SetOuterSolver(PREONLY);


  // BEGIN Testing the class Cloud

  std::vector<std::string> velocity = {"U", "V"};
  std::cout << "Testing the class Cloud \n";


  unsigned nIterations = 1000;
  
  
  cld->AddEllipse({XG, YG}, {RADIUS, RADIUS}, 8);
//   cld->AddQuadric({1.,0.,1.,-2.*XG ,-2*YG ,XG*XG+YG*YG-RADIUS*RADIUS}, 8);
//   cld->AddQuadric({-1./10,0.,0.,0.,+1.,0.}, 8);
  cld->ComputeQuadraticBestFit();

  cldint->AddInteriorEllipse({XG, YG}, {RADIUS, RADIUS});
//   cldint->AddInteriorQuadric({1.,0.,1.,-2.*XG ,-2*YG ,XG*XG+YG*YG-RADIUS*RADIUS});
//   cldint->AddInteriorQuadric({-1./10,0.,0.,0.,+1.,0.});
  cldint->RebuildInteriorMarkers(*cld, "C", "Cn");

  cld->PrintCSV("markerBefore", 0);
  cld->PrintCSV("marker", 0);
  cldint->PrintCSV("markerInt", 0);


  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");
  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(true);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);


  for(unsigned it = 1; it <= nIterations; it++) {

    mlSol.CopySolutionToOldSolution();

    system.MGsolve();
    double dt = system.GetIntervalTime();

    cld->RKAdvection(4, velocity, dt);
    cldint->RKAdvection(4, velocity, dt);
    cld->PrintCSV("markerBefore", it);
    cld->ComputeQuadraticBestFit();
    cld->RebuildMarkers(10, 10, 10);
    cldint->RebuildInteriorMarkers(*cld, "C", "Cn");
    cld->PrintCSV("marker", it);
    vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, it);

  }

  delete cld;
  delete cldint;

  return 0;
}

double TimeStepMultiphase(const double time) {
  double dt =  0.01;
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


  AssembleGhostPenalty(ml_prob);
  AssembleGhostPenaltyDGP(ml_prob, true);
  AssembleGhostPenaltyDGP(ml_prob, false);


  double dt =  mlPdeSys->GetIntervalTime();

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

  std::vector<double> g(dim);
  g[0] = 0.;
  g[1] = (dim == 2) ? - gravity * 0.98 : 0;
  if(dim == 3) g[2] = - gravity * 0.98;

  CutFemWeight <TypeIO, TypeA> tet  = CutFemWeight<TypeIO, TypeA >(TET, qM, "legendre");
  CDWeightQUAD <TypeA> quadCD(qM, dx, dtetha);
  CDWeightTRI <TypeA> triCD(qM, dx, dtetha);


  /* END cutfem stuff for surface tension integration */

// cld->AddEllipse({XG, YG}, {RADIUS, RADIUS}, nMax);

//   cld.RKAdvection(4, {"U", "V"}, dtetha); // TODO dtetha sbagliato
//   cld->PrintCSV("markerBefore",it);
//  cld->ComputeQuadraticBestFit();
//   cld->RebuildMarkers(8, 12, 8);
//   cld->PrintCSV("marker",it);

  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

//       for(unsigned iel = msh->_elementOffset[msh->processor_id()]; iel < msh->_elementOffset[msh->processor_id() + 1]; iel++) {
//       std::cout << "iel = " << iel << "   ";
//       const std::vector<double> &a = cld.GetQuadraticBestFitCoefficients(iel);
//       for(unsigned i = 0; i < a.size(); i++) std::cout << a[i] << "  ";
//       std::cout << "\n";
//     }
//     std::cout << std::endl;

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

    std::vector<double> a;
    std::vector<double> xm;
    double d;
    unsigned cut = 0;
    double vol;

    std::vector<std::vector<double>> Jacob, JacI;

    const elem_type *femV = msh->_finiteElement[ielGeom][solVType];
    const elem_type *femP = msh->_finiteElement[ielGeom][solPType];

    unsigned cnt = cld->GetNumberOfMarker(iel);

    if(cnt > 0) cut = 1;

    if(cut == 1) {
      femV = fem.GetFiniteElement(ielGeom, solVType);
      femP = fem.GetFiniteElement(ielGeom, solPType);
      femV->GetJacobianMatrix(coordX, cld->GetCloudBaricenterInParentElement(iel), weight, Jacob, JacI);
      cld->GetLinearFit(iel, Jacob, a, d);
    }

//     if(ielGeom == 3) GetNormalQuad(x1, xg, R, a, d, xm, b, db, cut);
//     if(ielGeom == 3) cld->GetLinearFit(iel, Jacob, b, db);
//     if(ielGeom == 4) GetNormalTri(x1, xg, R, a, d, xm, b, db, cut);
//     else if(ielGeom == 1) GetNormalTetBF(x1, xg, R, a, d, xm, b, db, vol, cut);
//     else if(ielGeom == 0) GetNormalHexBF(x1, xg, R, a, d, xm, b, db, vol, cut, fem.GetFiniteElement(0, 0));

    std::vector <TypeIO> weightCF(quad.GetGaussQuadraturePointNumber(), 0.);
    std::vector <TypeIO> weightCFInt(quad.GetGaussQuadraturePointNumber(), 0.);
    std::vector <TypeIO> weightCFExt(quad.GetGaussQuadraturePointNumber(), 0.);

    if(cut == 1) {
      bool wMap = 1;
      if(ielGeom == 3) {
        quad.GetWeightWithMap(0, a, d, weightCFExt);
        for(unsigned k = 0; k < dim; k++) a[k] = - a[k];
        d = -d;
        quad.GetWeightWithMap(-1, a, d, weightCF);
        quad.GetWeightWithMap(0, a, d, weightCFInt);

//           quadCD.GetWeight(a, d, weightCF);
      }
//         else if(ielGeom == 4) {
//           triCD.GetWeight(b, db, weightCF);
//           const double* weightG = tri.GetGaussWeightPointer();
//         }
//         else if(ielGeom == 1) {
//           tet.GetWeightWithMap(-1, b, db, weightCF);
//           const double* weightG = tet.GetGaussWeightPointer();
//         }
    }
    else {
      for(unsigned i = 0; i < weightCFInt.size(); i++) {
        weightCFInt[i] = C;
        weightCFExt[i] = 1. - C;
      }
    }

    std::vector<double> xqp(dim);
    std::vector<double> NN(dim, 0.);
    double kk = 0.;

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < femV->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      femV->Jacobian(coordX, ig, weight, phiV, phiV_x);
      phiP = femP->GetPhi(ig);

      double dsN = 0.;
      std::vector <double> Nf(dim, 0); // unit normal in the physical element from the fluid to the solid

      if(cut == 1) {

        femV->GetJacobianMatrix(coordX, ig, weight, Jacob, JacI);

        for(unsigned k = 0; k < dim; k++) {
          for(unsigned j = 0; j < dim; j++) {
            Nf[k] += JacI[j][k] * a[j];
          }
          dsN += Nf[k] * Nf[k];
        }
        dsN = sqrt(dsN);
        for(unsigned k = 0; k < dim; k++) {
          Nf[k] /= dsN;
        }
      }



      for(unsigned k = 0; k < dim; k++) {
        xqp[k] = 0.;
      }

      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned k = 0; k < dim; k++) {
          xqp[k] += coordX[k][i] * phiV[i];
        }
      }

      if(cld->GetNumberOfMarker(iel) > 0) {
        double magN2 = 0.;
//         kk = cld->getCurvature(iel, xqp);
       kk = cld->getAverageCurvature(iel);
        NN = cld->getNormal(iel, xqp);
//       kk = CurvatureQuadric({1., 1., 0., - 2 * XG, - 2 * YG, XG * XG + YG * YG - RADIUS * RADIUS}, xqp);
//       kk = 1. / RADIUS;
//       NormalQuadric({1., 1., 0., - 2 * XG, - 2 * YG, XG * XG + YG * YG - RADIUS * RADIUS}, xqp, NN); //TODO
//       for(unsigned k = 0; k < dim; k++) magN2 += NN[k] * NN[k];
//       for(unsigned k = 0; k < dim; k++) NN[k] /= sqrt(magN2);
      }

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

//       double rho = rho1 * weightCFInt[ig] + rho1 * weightCFExt[ig];
//       double mu = mu1 * weightCFInt[ig] + mu2 * weightCFExt[ig];

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
          NSV += - phiV_x[i * dim + I] * (solP1_gss * weightCFInt[ig] + solP2_gss * weightCFExt[ig]);  // pressure gradient
          NSV += rho * phiV[i] * (solV_gss[I] - solVOld_gss[I]) / dt ;
          NSV += - rhoC * phiV[i] * g[I]; // gravity term
          Res[I * nDofsV + i] -=  NSV * weight;
          if(cut == 1) {
            Res[I * nDofsV + i] += - sigma * phiV[i] /** b[I]*/ * NN[I] * weight * weightCF[ig] * kk * dsN; 
          }
        }
      } // end phiV_i loop

      // *** phiP_i loop ***
      for(unsigned i = 0; i < nDofsP; i++) {
        for(int I = 0; I < dim; I++) {
          Res[dim * nDofsV + i] += - gradSolV_gss[I][I] * phiP[i]  * weight * weightCFInt[ig]; //continuity
          Res[dim * nDofsV + nDofsP + i] += - gradSolV_gss[I][I] * phiP[i]  * weight * weightCFExt[ig]; //continuity
        }
        if(C == 0) 
        Res[dim * nDofsV + i] += - solP1_gss * phiP[i]  * weight * (1 - C) * eps; //penalty
        if(C == 1) 
        Res[dim * nDofsV + nDofsP + i] += - solP2_gss * phiP[i]  * weight * C * eps; //penalty

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
            Jac[VIrow * nDofsVP + P1column] += - phiV_x[i * dim + I] * phiP[j] * weight * weightCFInt[ig]; //pressure gradient
            Jac[P1column * nDofsVP + VIrow] -= - phiV_x[i * dim + I] * phiP[j] * weight * weightCFInt[ig]; //continuity
            Jac[VIrow * nDofsVP + P2column] += - phiV_x[i * dim + I] * phiP[j] * weight * weightCFExt[ig]; //pressure gradient
            Jac[P2column * nDofsVP + VIrow] -= - phiV_x[i * dim + I] * phiP[j] * weight * weightCFExt[ig]; //continuity
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
          Jac[P1row * nDofsVP + P1column] += phiP[i] * phiP[j] * weight * (1 - C) * eps; //continuity
          if(C == 1) 
          Jac[P2row * nDofsVP + P2column] += phiP[i] * phiP[j] * weight * C * eps; //continuity
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


void GetNormalQuad(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double & R, std::vector<double> &a, double & d,  std::vector<double> &xm, std::vector<double> &b, double & db, unsigned & cut) {

  const unsigned &dim =  xv.size();
  const unsigned &nve =  xv[0].size();

  const double& x1 = xv[0][0];
  const double& x2 = xv[0][1];
  const double& x3 = xv[0][2];
  const double& x4 = xv[0][3];
  const double& y1 = xv[1][0];
  const double& y2 = xv[1][1];
  const double& y3 = xv[1][2];
  const double& y4 = xv[1][3];

  double hx = 0.5 * (fabs(x3 - x1) + fabs(x4 - x2));
  double hy = 0.5 * (fabs(y3 - y1) + fabs(y4 - y2));
  double h = sqrt(hx * hx + hy * hy);
  double eps = 1.0e-10 * h;

  std::vector<double> dist(nve, 0);
  std::vector<double> dist0(nve);
  unsigned cnt0 = 0;


  std::vector<double> A(2, 0.);
  std::vector<std::vector<double>> xe(2, std::vector<double>(4));
  double D = 0.;
  unsigned intMax = 2;
  unsigned nEdge = 0;
  unsigned cnt = 0;

  for(unsigned i = 0; i < nve; i++) {
    unsigned ip1 = (i + 1) % nve;
    A[0] = xv[1][ip1] - xv[1][i];
    A[1] = - xv[0][ip1] + xv[0][i];
    D = - A[0] * xv[0][i] - A[1] * xv[1][i];


    std::vector<double> inters(intMax, 0.);
    unsigned dir = (fabs(A[0]) > fabs(A[1])) ? 1 : 0 ;
    unsigned dirp1 = (dir + 1) % 2;

    double iMax = std::max(xv[dir][ip1], xv[dir][i]);
    double iMin = std::min(xv[dir][ip1], xv[dir][i]);

    double delta = ((A[0] * A[0] + A[1] * A[1]) * R * R) - (D + A[0] * xg[0] + A[1] * xg[1]) * (D + A[0] * xg[0] + A[1] * xg[1]);
    a.resize(dim);
    if(delta > 0.) {
      inters[0] = (- A[dir] * (D + A[dirp1] * xg[dirp1]) + A[dirp1] * (A[dirp1] * xg[dir] - sqrt(delta))) / (A[0] * A[0] + A[1] * A[1]);
      inters[1] = (- A[dir] * (D + A[dirp1] * xg[dirp1]) + A[dirp1] * (A[dirp1] * xg[dir] + sqrt(delta))) / (A[0] * A[0] + A[1] * A[1]);
      unsigned nInt = 0;
      unsigned jInt = 2;
      for(unsigned j = 0; j < intMax; j++) {
        if(inters[j] < iMax && inters[j] > iMin) {
          nInt++;
          jInt = j;
        }
      }
      if(nInt == 1) {
        xe[dir][cnt] = inters[jInt];
        xe[dirp1][cnt] = (- D - A[dir] * xe[dir][cnt]) / A[dirp1];
        cnt++;
      }
    }
  }
  if(cnt == 0) {
    cut = (R * R - (xv[0][0] - xg[0]) * (xv[0][0] - xg[0]) - (xv[1][0] - xg[1]) * (xv[1][0] - xg[1]) > 0) ? 0 : 2; //0 inside, 2 outside
    return;
  }
  else if(cnt == 4) { // small ball no good
    cut = 0;
    return;
  }
  else if(cnt == 2) {
    cut = 1;
    std::vector<double> theta(2);

    a[0] = xe[1][1] - xe[1][0] ;
    a[1] = - xe[0][1] + xe[0][0] ;

    xm.resize(2);
    xm[0] = 0.5 * (xe[0][0] + xe[0][1]);
    xm[1] = 0.5 * (xe[1][0] + xe[1][1]);

    double det = 0;
    for(unsigned k = 0; k < dim; k++) {
      det += a[k] * (xg[k] - xm[k]);
    }
    double sign = (det >= 0) ? 1. : -1.;

    double norm = sign * sqrt(a[0] * a[0] + a[1] * a[1]);
    a[0] /= norm;
    a[1] /= norm;

    theta[0] = atan2(xe[1][0]  - xg[1], xe[0][0] - xg[0]);
    theta[1] = atan2(xe[1][1]  - xg[1], xe[0][1] - xg[0]);

    if(theta[0] > theta[1]) {
      std::swap(theta[0], theta[1]);
    }
    double DT = theta[1] - theta[0];
    if(DT > M_PI) {
      std::swap(theta[0], theta[1]);
      theta[1] += 2. * M_PI;
      DT = theta[1] - theta[0];
    }
    xm.resize(dim);

    d = R * sqrt(0.5 * DT / tan(0.5 * DT)) ;
    a.resize(dim);
    a[0] = -cos(theta[0] + 0.5 * DT);
    a[1] = -sin(theta[0] + 0.5 * DT);

    for(unsigned k = 0; k < dim; k++) {
      xm[k] = -a[k] * d + xg[k];
    }
    d += - a[0] * xg[0] - a[1] * xg[1]; //TODO

    double d2 = sqrt(pow(xm[0] - xg[0], 2) + pow(xm[1] - xg[1], 2));
    d = d2 * tan(0.5 * DT);

//     std::cout << "xm = " << xm[0] << " " << xm[1] << std::endl;
//     std::cout << "a = " << a[0] << " b = " << a[1] << std::endl;





//   for(unsigned i = 0; i < nve; i++) {
//     for(unsigned k = 0;  k < dim; k++) {
//       dist[i] += (xv[k][i] - xg[k]) * (xv[k][i] - xg[k]);
//     }
//     dist[i] = sqrt(dist[i]) - R;
//
//     if(fabs(dist[i]) < eps) {
//       dist0[i] = (dist[i] < 0) ? -eps : eps;
//       dist[i] = 0.;
//       cnt0++;
//     }
//     else {
//       dist0[i] = dist[i];
//     }
//   }
//
//   if(cnt0 > 0) {
//     unsigned cntp = 0;
//     for(unsigned i = 0; i < nve; i++) {
//       if(dist[i] > 0) cntp++;
//       dist[i] = dist0[i];
//     }
//     if(cntp == 0) { // the element is inside the ball
//       cut = 0;
//       return;
//     }
//     else if(cntp == nve - cnt0) {  // the element in outside the ball
//       cut = 2;
//       return;
//     }
//   }
//
//   std::vector <double> theta(2);
//   unsigned cnt = 0;
//   for(unsigned e = 0; e < nve; e++) {
//     unsigned ep1 = (e + 1) % nve;
//     if(dist[e] * dist[ep1] < 0) {
//       double s = 0.5  * (1 + (dist[e] + dist[ep1]) / (dist[e] - dist[ep1]));
//       theta[cnt] = atan2((1 - s) * xv[1][e] + s * xv[1][ep1]  - xg[1], (1 - s) * xv[0][e] + s * xv[0][ep1] - xg[0]) ;
//       cnt++;
//     }
//   }
//
//   if(cnt == 0) {
//     if(dist[0] < 0) cut = 0; // cell inside the ball
//     else cut = 2; // cell outside the ball
//     return;
//   }
//   else {
//     cut = 1;
//     if(theta[0] > theta[1]) {
//       std::swap(theta[0], theta[1]);
//     }
//     double DT = theta[1] - theta[0];
//     if(DT > M_PI) {
//       std::swap(theta[0], theta[1]);
//       theta[1] += 2. * M_PI;
//       DT = theta[1] - theta[0];
//     }
//     xm.resize(dim);
//
//     d = R * sqrt(0.5 * DT / tan(0.5 * DT)) ;
//     a.resize(dim);
//     a[0] = -cos(theta[0] + 0.5 * DT);
//     a[1] = -sin(theta[0] + 0.5 * DT);
//
//     for(unsigned k = 0; k < dim; k++) {
//       xm[k] = -a[k] * d + xg[k];
//     }
//     d += - a[0] * xg[0] - a[1] * xg[1]; //TODO
//
//     double d2 = sqrt(pow(xm[0] - xg[0], 2) + pow(xm[1] - xg[1], 2));
//     d = d2 * tan(0.5 * DT);
//
//     std::cout.precision(14);
//
//     std::cout << "xm = " << xm[0] << " " << xm[1] << std::endl;
//     std::cout << "a = " << a[0] << " b = " << a[1] << " d = " << d << std::endl;

    std::vector<double> xi(dim);
    double &u = xi[0];
    double &v = xi[1];

    std::vector < std::vector < double > > J(2, std::vector<double>(2));

    double dx12 = x1 - x2;
    double dx34 = x3 - x4;
    double dy12 = y1 - y2;
    double dy34 = y3 - y4;
    double hu = dx34 * dy12 - dx12 * dy34;

    double dx14 = (x1 - x4);
    double dy23 = (y2 - y3);
    double dx23 = (x2 - x3);
    double dy14 = (y1 - y4);
    double hv = dx14 * dy23 - dx23 * dy14;

    double eps2 = 1.0e-10 * h * h;

    if(fabs(hu) > eps2) {//edges 1 and 3 are not parallel
      double gu = -x4 * y1 + x3 * y2 - x2 * y3 + x1 * y4;
      double f = xm[0] * (dy12 + dy34) - xm[1] * (dx12 + dx34);
      double fpgu = f + gu;

      double det = sqrt(hu * (- 2. * xm[0] * (dy14 + dy23)
                              + (2. * xm[1] - y3 - y4) * (x1 + x2)
                              - (2. * xm[1] - y1 - y2) * (x3 + x4))
                        + fpgu * fpgu);
      u = (fpgu + det) / hu;

      if(fabs(hv) > eps2) { //edges 2 and 4 are not parallel
        double gv = -x4 * y3 + x3 * y4 - x2 * y1 + x1 * y2;
        v = (f + gv - det) / hv;

        J[0][0] = 0.25 * ((-1. + v) * dx12 + (1. + v) * dx34);
        J[0][1] = 0.25 * ((-1. + u) * dx14 - (1. + u) * dx23);
        J[1][0] = 0.25 * ((-1. + v) * dy12 + (1. + v) * dy34);
        J[1][1] = 0.25 * ((-1. + u) * dy14 - (1. + u) * dy23);

      }
      else { //edges 2 and 4 are parallel
        //   std::cout << "2 and 4 are parallel\n";
        J[0][1] = 0.25 * ((-1. + u) * dx14 - (1. + u) * dx23);
        J[1][1] = 0.25 * ((-1. + u) * dy14 - (1. + u) * dy23);

        v = (J[0][1] > eps) ?
            (0.25 * ((-1. + u) * (x1 + x4) - (1. + u) * (x3 + x2)) + xm[0]) / J[0][1] :
            (0.25 * ((-1. + u) * (y1 + y4) - (1. + u) * (y3 + y2)) + xm[1]) / J[1][1];

        J[0][0] = 0.25 * ((-1. + v) * dx12 + (1. + v) * dx34);
        J[1][0] = 0.25 * ((-1. + v) * dy12 + (1. + v) * dy34);

      }
    }
    else if(fabs(hv) > eps2) {  //edges 1 and 3 are parallel, but edges 2 and 4 are not
      // std::cout << "1 and 3 are parallel\n";
      double f = xm[0] * (dy12 + dy34) - xm[1] * (dx12 + dx34);
      double gv = -x4 * y3 + x3 * y4 - x2 * y1 + x1 * y2;
      double fpgv = f + gv;

      double det = sqrt(hv * (- 2. * xm[0] * (dy12 - dy34)
                              + (2. * xm[1] - y2 - y3) * (x1 + x4)
                              - (2. * xm[1] - y1 - y4) * (x2 + x3))
                        +  fpgv * fpgv);

      v = (fpgv - det) / hv;

      J[0][0] = 0.25 * ((-1. + v) * dx12 + (1. + v) * dx34);
      J[1][0] = 0.25 * ((-1. + v) * dy12 + (1. + v) * dy34);

      u = (fabs(J[0][0]) > eps) ?
          (0.25 * ((-1. + v) * (x1 + x2) - (1. + v) * (x3 + x4)) + xm[0]) / J[0][0] :
          (0.25 * ((-1. + v) * (y1 + y2) - (1. + v) * (y3 + y4)) + xm[1]) / J[1][0];

      J[0][1] = 0.25 * ((-1. + u) * dx14 - (1. + u) * dx23);
      J[1][1] = 0.25 * ((-1. + u) * dy14 - (1. + u) * dy23);
    }
    else { //edges 1 and 3, and  edges 2 and 4 are parallel
      //   std::cout << "Romboid\n";
      std::vector<std::vector<unsigned> > idx = {{3, 1}, {0, 2}};

      double A[2][2] = {{ -dy14, dy23}, {dy12, dy34}};
      double B[2][2] = {{dx14, -dx23}, { -dx12, -dx34}};

      for(unsigned k = 0; k < 2; k++) {
        double d[2];
        for(unsigned j = 0 ; j < 2; j++) {
          double Ckj = - A[k][j] * xv[0][idx[k][j]] - B[k][j] * xv[1][idx[k][j]];
          d[j] = (A[k][j] * xm[0] + B[k][j] * xm[1] + Ckj) / sqrt(A[k][j] * A[k][j] + B[k][j] * B[k][j]);
        }
        xi[k] = -1. + 2. * d[0] / (d[0] + d[1]);
      }

      J[0][0] = 0.25 * ((-1. + v) * dx12 + (1. + v) * dx34);
      J[0][1] = 0.25 * ((-1. + u) * dx14 - (1. + u) * dx23);
      J[1][0] = 0.25 * ((-1. + v) * dy12 + (1. + v) * dy34);
      J[1][1] = 0.25 * ((-1. + u) * dy14 - (1. + u) * dy23);
    }

//     double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
//     std::vector < std::vector < double > > Ji = {{J[1][1] / det, -J[0][1] / det}, {-J[1][0] / det, J[0][0] / det}};

    b.assign(dim, 0);
    for(unsigned k = 0; k < dim; k++) {
      for(unsigned j = 0; j < dim; j++) {
        b[k] += J[j][k] * a[j];
      }
    }
    double bNorm = sqrt(b[0] * b[0] + b[1] * b[1]);
    b[0] /= bNorm;
    b[1] /= bNorm;
    db = - b[0] * xi[0] - b[1] * xi[1];
//   std::cout << b[0] << " " << b[1] << " " << db << " " << std::endl;

//     // Old inverse mapping for comparison
//
//     std::vector <  std::vector < std::vector <double > > > aP(1);
//     short unsigned quad = 3;
//     unsigned linear = 0;
//     bool ielIsInitialized = false;
//     if(!ielIsInitialized) {
//       ielIsInitialized = true;
//       ProjectNodalToPolynomialCoefficients(aP[0], xv, quad, linear) ;
//     }
//
//     std::vector<double> xib(dim);
//     GetClosestPointInReferenceElement(xv, xm, quad, xib);
//     bool inverseMapping = GetInverseMapping(linear, quad, aP, xm, xib, 100);
//     if(!inverseMapping) {
//       std::cout << "InverseMapping failed" << std::endl;
//     }
//
//     //std::cout << xib[0] << " " << xib[1] << std::endl;
//
//     vector < vector < double > > Jac2;
//     vector < vector < double > > JacI2;
//     finiteElementQuad->GetJacobianMatrix(xv, xib, Jac2, JacI2);
//
//
//     //std::swap(JacI[0][1],JacI[1][0]);
//     //std::cout << Jac[0][0] << " " << Jac[0][1] << " " << Jac[1][0] << " " <<Jac[1][1] << std::endl;
//
//     std::vector <double> b2(dim, 0.);
//
//     for(unsigned k = 0; k < dim; k++) {
//       for(unsigned j = 0; j < dim; j++) {
//         b2[k] += JacI2[j][k] * a[j];
//       }
//     }
//     double b2Norm = sqrt(b2[0] * b2[0] + b2[1] * b2[1]);
//     b2[0] /= b2Norm;
//     b2[1] /= b2Norm;
//     double db2 = - b2[0] * xi[0] - b2[1] * xi[1];
//     std::cout << b2[0] << " " << b2[1] << " " << db2 << " " << std::endl;

  }

}

void GetNormalTri(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double & R, std::vector<double> &a, double & d,  std::vector<double> &xm, std::vector<double> &b, double & db, unsigned & cut) {

  const unsigned &dim =  xv.size();
  const unsigned &nve =  xv[0].size();

  const double& x1 = xv[0][0];
  const double& x2 = xv[0][1];
  const double& x3 = xv[0][2];
  const double& y1 = xv[1][0];
  const double& y2 = xv[1][1];
  const double& y3 = xv[1][2];




  std::vector<double> A(2, 0.);
  std::vector<std::vector<double>> xe(2, std::vector<double>(4));
  double D = 0.;
  unsigned intMax = 2;
  unsigned nEdge = 0;
  unsigned cnt = 0;

  for(unsigned i = 0; i < nve; i++) {
    unsigned ip1 = (i + 1) % nve;
    A[0] = xv[1][ip1] - xv[1][i];
    A[1] = - xv[0][ip1] + xv[0][i];
    D = - A[0] * xv[0][i] - A[1] * xv[1][i];


    std::vector<double> inters(intMax, 0.);
    unsigned dir = (fabs(A[0]) > fabs(A[1])) ? 1 : 0 ;
    unsigned dirp1 = (dir + 1) % 2;

    double iMax = std::max(xv[dir][ip1], xv[dir][i]);
    double iMin = std::min(xv[dir][ip1], xv[dir][i]);

    double delta = ((A[0] * A[0] + A[1] * A[1]) * R * R) - (D + A[0] * xg[0] + A[1] * xg[1]) * (D + A[0] * xg[0] + A[1] * xg[1]);
    a.resize(dim);
    if(delta > 0.) {
      inters[0] = (- A[dir] * (D + A[dirp1] * xg[dirp1]) + A[dirp1] * (A[dirp1] * xg[dir] - sqrt(delta))) / (A[0] * A[0] + A[1] * A[1]);
      inters[1] = (- A[dir] * (D + A[dirp1] * xg[dirp1]) + A[dirp1] * (A[dirp1] * xg[dir] + sqrt(delta))) / (A[0] * A[0] + A[1] * A[1]);
      unsigned nInt = 0;
      unsigned jInt = 2;
      for(unsigned j = 0; j < intMax; j++) {
        if(inters[j] < iMax && inters[j] > iMin) {
          nInt++;
          jInt = j;
        }
      }
      if(nInt == 1) {
        xe[dir][cnt] = inters[jInt];
        xe[dirp1][cnt] = (- D - A[dir] * xe[dir][cnt]) / A[dirp1];
        cnt++;
      }
    }
  }
  if(cnt == 0) {
    cut = (R * R - (xv[0][0] - xg[0]) * (xv[0][0] - xg[0]) - (xv[1][0] - xg[1]) * (xv[1][0] - xg[1]) > 0) ? 0 : 2;
    return;
  }
  else if(cnt == 4) {
    cut = 0;
    return;
  }
  else if(cnt == 2) {
    cut = 1;
    std::vector<double> theta(2);

    a[0] = xe[1][1] - xe[1][0] ;
    a[1] = - xe[0][1] + xe[0][0] ;

    xm.resize(2);
    xm[0] = 0.5 * (xe[0][0] + xe[0][1]);
    xm[1] = 0.5 * (xe[1][0] + xe[1][1]);

    double det = 0;
    for(unsigned k = 0; k < dim; k++) {
      det += a[k] * (xg[k] - xm[k]);
    }
    double sign = (det >= 0) ? 1. : -1.;

    double norm = sign * sqrt(a[0] * a[0] + a[1] * a[1]);
    a[0] /= norm;
    a[1] /= norm;

    theta[0] = atan2(xe[1][0]  - xg[1], xe[0][0] - xg[0]);
    theta[1] = atan2(xe[1][1]  - xg[1], xe[0][1] - xg[0]);

    if(theta[0] > theta[1]) {
      std::swap(theta[0], theta[1]);
    }
    double DT = theta[1] - theta[0];
    if(DT > M_PI) {
      std::swap(theta[0], theta[1]);
      theta[1] += 2. * M_PI;
      DT = theta[1] - theta[0];
    }
    xm.resize(dim);

    d = R * sqrt(0.5 * DT / tan(0.5 * DT)) ;
    a.resize(dim);
    a[0] = -cos(theta[0] + 0.5 * DT);
    a[1] = -sin(theta[0] + 0.5 * DT);

    for(unsigned k = 0; k < dim; k++) {
      xm[k] = -a[k] * d + xg[k];
    }
    d += - a[0] * xg[0] - a[1] * xg[1]; //TODO

    double d2 = sqrt(pow(xm[0] - xg[0], 2) + pow(xm[1] - xg[1], 2));
    d = d2 * tan(0.5 * DT);

//     std::cout << "xm = " << xm[0] << " " << xm[1] << std::endl;
//     std::cout << "a = " << a[0] << " b = " << a[1] << std::endl;





//   double hx = (fabs(x2 - x1) + fabs(x3 - x2) + fabs(x3 - x1)) / 3.;
//   double hy = (fabs(y2 - y1) + fabs(y3 - y2) + fabs(y3 - y1)) / 3.;
//
//   double h = sqrt(hx * hx + hy * hy);
//   double eps = 1.0e-10 * h;
//
//   std::vector<double> dist(nve, 0);
//   std::vector<double> dist0(nve);
//   unsigned cnt0 = 0;
//   for(unsigned i = 0; i < nve; i++) {
//     for(unsigned k = 0;  k < dim; k++) {
//       dist[i] += (xv[k][i] - xg[k]) * (xv[k][i] - xg[k]);
//     }
//     dist[i] = sqrt(dist[i]) - R;
//
//     //std::cout << dist[i] << std::endl;
//
//     if(fabs(dist[i]) < eps) {
//       dist0[i] = (dist[i] < 0) ? -eps : eps;
//       dist[i] = 0.;
//       cnt0++;
//     }
//     else {
//       dist0[i] = dist[i];
//     }
//   }
//
//   if(cnt0 > 0) {
//     unsigned cntp = 0;
//     for(unsigned i = 0; i < nve; i++) {
//       if(dist[i] > 0) cntp++;
//       dist[i] = dist0[i];
//     }
//     if(cntp == 0) { // the element is inside the ball
//       cut = 0;
//       return;
//     }
//     else if(cntp == nve - cnt0) {  // the element in outside the ball
//       cut = 2;
//       return;
//     }
//   }
//
//   std::vector <double> theta(2);
//   unsigned cnt = 0;
//   for(unsigned e = 0; e < nve; e++) {
//     unsigned ep1 = (e + 1) % nve;
//     if(dist[e] * dist[ep1] < 0) {
//       double s = 0.5  * (1 + (dist[e] + dist[ep1]) / (dist[e] - dist[ep1]));
//       theta[cnt] = atan2((1 - s) * xv[1][e] + s * xv[1][ep1]  - xg[1], (1 - s) * xv[0][e] + s * xv[0][ep1] - xg[0]) ;
//       cnt++;
//     }
//   }
//
//   if(cnt == 0) {
//     if(dist[0] < 0) cut = 0; // cell inside the ball
//     else cut = 2; // cell outside the ball
//     return;
//   }
//   else {
//     cut = 1;
//     if(theta[0] > theta[1]) {
//       std::swap(theta[0], theta[1]);
//     }
//     double DT = theta[1] - theta[0];
//     if(DT > M_PI) {
//       std::swap(theta[0], theta[1]);
//       theta[1] += 2. * M_PI;
//       DT = theta[1] - theta[0];
//     }
//     xm.resize(dim);
//
//     d = R * sqrt(0.5 * DT / tan(0.5 * DT)) ;
//     a.resize(dim);
//     a[0] = -cos(theta[0] + 0.5 * DT);
//     a[1] = -sin(theta[0] + 0.5 * DT);
//
//     for(unsigned k = 0; k < dim; k++) {
//       xm[k] = -a[k] * d + xg[k];
//     }
//     d += - a[0] * xg[0] - a[1] * xg[1]; //TODO
//
//     //std::cout << "xm = " << xm[0] << " " << xm[1] << std::endl;
//     //std::cout << "a = " << a[0] << " b = " << a[1] << " d = " << d << std::endl;
//
//     double d2 = sqrt(pow(xm[0] - xg[0], 2) + pow(xm[1] - xg[1], 2));
//     d = d2 * tan(0.5 * DT);
//
//     std::cout.precision(14);



    std::vector<double> xi(dim);

    std::vector < std::vector < double > > J(2, std::vector<double>(2));
    J[0][0] = (-x1 + x2);
    J[0][1] = (-x1 + x3);

    J[1][0] = (-y1 + y2);
    J[1][1] = (-y1 + y3);

    double den = (x3 * y1 - x1 * y3 + x2 * J[1][1] - y2 * J[0][1]);

    xi[0] = (x3 * y1 - x1 * y3 + xm[0] * J[1][1] - xm[1] * J[0][1]) / den;
    xi[1] = (x1 * y2 - x2 * y1 - xm[0] * J[1][0] + xm[1] * J[0][0]) / den;


    //std::cout << xi[0] << " " << xi[1] << std::endl;


    b.assign(dim, 0);
    for(unsigned k = 0; k < dim; k++) {
      for(unsigned j = 0; j < dim; j++) {
        b[k] += J[j][k] * a[j];
      }
    }
    double bNorm = sqrt(b[0] * b[0] + b[1] * b[1]);
    b[0] /= bNorm;
    b[1] /= bNorm;
    db = - b[0] * xi[0] - b[1] * xi[1];


    //std::cout << b[0] << " " << b[1] << " " << db << " " << std::endl;
  }
}


void GetNormalTet(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double & R, std::vector<double> &a, double & d,  std::vector<double> &xm, std::vector<double> &a2, double & d2, double & volume,  unsigned & cut) {

  const unsigned dim =  3;
  const unsigned nve =  4;

  //std::cout<<nve<<std::endl;

  const double& x1 = xv[0][0];
  const double& x2 = xv[0][1];
  const double& x3 = xv[0][2];
  const double& x4 = xv[0][3];
  const double& y1 = xv[1][0];
  const double& y2 = xv[1][1];
  const double& y3 = xv[1][2];
  const double& y4 = xv[1][3];
  const double& z1 = xv[2][0];
  const double& z2 = xv[2][1];
  const double& z3 = xv[2][2];
  const double& z4 = xv[2][3];

  double hx = (fabs(x2 - x1) + fabs(x3 - x2) + fabs(x3 - x1) + fabs(x4 - x1) + fabs(x4 - x2) + fabs(x4 - x3)) / 6.;
  double hy = (fabs(y2 - y1) + fabs(y3 - y2) + fabs(y3 - y1) + fabs(y4 - y1) + fabs(y4 - y2) + fabs(y4 - y3)) / 6.;
  double hz = (fabs(z2 - z1) + fabs(z3 - z2) + fabs(z3 - z1) + fabs(z4 - z1) + fabs(z4 - z2) + fabs(z4 - z3)) / 6.;
  double h = sqrt(hx * hx + hy * hy + hz * hz);
  double eps = 1.0e-10 * h;




  std::vector<double> A(dim, 0.);
  std::vector < std::vector <double> > y(4, std::vector<double>(dim));
  std::vector < unsigned > i0(6);
  double D = 0.;
  unsigned intMax = 2;
  unsigned nEdge = 0;
  unsigned cnt = 0;

  for(unsigned i = 0; i < nve - 1; i++) {
    for(unsigned j = i + 1; j < nve; j++) {
      for(unsigned k = 0; k < dim; k++) {
        A[k] = xv[k][j] - xv[k][i];
      }

      std::vector<double> inters(intMax, 0.);
      unsigned dir = (fabs(A[0]) > fabs(A[1])) ? ((fabs(A[0]) > fabs(A[2])) ? 0 : 2) : ((fabs(A[1]) > fabs(A[2])) ? 1 : 2) ;
      unsigned dirp1 = (dir + 1) % dim;
      unsigned dirp2 = (dir + 2) % dim;

      double iMax = std::max(xv[dir][j], xv[dir][i]);
      double iMin = std::min(xv[dir][j], xv[dir][i]);

      double Axdi[3] = {A[1]* (xg[2] - xv[2][i]) - A[2] * (xg[1] - xv[1][i]),
                        A[2]* (xg[0] - xv[0][i]) - A[0] * (xg[2] - xv[2][i]),
                        A[0]* (xg[1] - xv[1][i]) - A[1] * (xg[0] - xv[0][i])
                       } ;

      double Addi = A[0] * (xg[0] - xv[0][i]) +
                    A[1] * (xg[1] - xv[1][i]) +
                    A[2] * (xg[2] - xv[2][i]);


      double den = A[0] * A[0] + A[1] * A[1] + A[2] * A[2];

      double delta = (- (A[dir] * A[dir]) *
                      (Axdi[0] * Axdi[0] + Axdi[1] * Axdi[1] + Axdi[2] * Axdi[2] - R * R * den));



//       double var = A[dir] * A[dir] * xg[dir] + (A[dirp1] * A[dirp1] + A[dirp2] * A[dirp2]) * xv[dir][i]
//                    + A[dir] * A[dirp1] * (xg[dirp1] - xv[dirp1][i])
//                    + A[dir] * A[dirp2] * (xg[dirp2] - xv[dirp2][i]);


      double var = den * xv[dir][i] + A[dir] * Addi;

      a.resize(dim);
      if(delta > 0.) {
        inters[0] = (var - sqrt(delta)) / den;
        inters[1] = (var + sqrt(delta)) / den;
        unsigned nInt = 0;
        unsigned jInt = 2;
        for(unsigned ii = 0; ii < intMax; ii++) {
          if(inters[ii] < iMax && inters[ii] > iMin) {
            nInt++;
            jInt = ii;
          }
        }
        if(nInt == 1) {
          y[cnt][dir] = inters[jInt];
          y[cnt][dirp1] = xv[dirp1][i] + A[dirp1] * (y[cnt][dir] - xv[dir][i]) / A[dir];
          y[cnt][dirp2] = xv[dirp2][i] + A[dirp2] * (y[cnt][dir] - xv[dir][i]) / A[dir];
          i0[cnt] = (i + j) - (i == 0);
          cnt++;
        }
      }
    }
  }
  if(cnt == 0) {
    cut = (R * R - (xv[0][0] - xg[0]) * (xv[0][0] - xg[0]) - (xv[1][0] - xg[1]) * (xv[1][0] - xg[1])  - (xv[2][0] - xg[2]) * (xv[2][0] - xg[2]) > 0) ? 0 : 2;
    return;
  }
  else if(cnt > 4) {
    cut = 0;
    return;
  }
  else if(cnt == 4 || cnt == 3) {
    cut = 1;



//   std::vector<double> dist(nve, 0);
//   std::vector<double> dist0(nve);
//   unsigned cnt0 = 0;
//   for(unsigned i = 0; i < nve; i++) {
//     for(unsigned k = 0;  k < dim; k++) {
//       dist[i] += (xv[k][i] - xg[k]) * (xv[k][i] - xg[k]);
//     }
//     dist[i] = sqrt(dist[i]) - R;
//
//     if(fabs(dist[i]) < eps) {
//       dist0[i] = (dist[i] < 0) ? -eps : eps;
//       dist[i] = 0.;
//       cnt0++;
//     }
//     else {
//       dist0[i] = dist[i];
//     }
//   }
//
//   if(cnt0 > 0) {
//     unsigned cntp = 0;
//     for(unsigned i = 0; i < nve; i++) {
//       if(dist[i] > 0) cntp++;
//       dist[i] = dist0[i];
//     }
//     if(cntp == 0) { // the element is inside the ball
//       cut = 0;
//       return;
//     }
//     else if(cntp == nve - cnt0) {  // the element in outside the ball
//       cut = 2;
//       return;
//     }
//   }
//
//   std::vector < std::vector <double> > y(4, std::vector<double>(dim));
//   std::vector < unsigned > i0(4);
//   unsigned cnt = 0;
//   for(unsigned i = 0; i < nve - 1; i++) {
//     for(unsigned j = i + 1; j < nve; j++) {
//       if(dist[i] * dist[j] < 0) {
//         double s = dist[i] / (dist[i] - dist[j]);
//         for(unsigned k = 0; k < dim; k++) {
//           y[cnt][k] = (1. - s) * xv[k][i] + s * xv[k][j];
//         }
//         i0[cnt] = (i + j) - (i == 0);
//         cnt++;
//       }
//     }
//   }
//
//   if(cnt == 0) {
//     if(dist[0] < 0) cut = 0; // cell inside the ball
//     else cut = 2; // cell outside the ball
//     return;
//   }
//   else {
//     cut = 1;

    if(cnt == 4) {

      if((i0[0] == 0 && i0[1] == 1) || (i0[0] == 1 && i0[1] == 2)) {
        std::swap(y[2], y[3]);
        std::swap(i0[2], i0[3]);
      }
      else {
        std::swap(y[1], y[3]);
        std::swap(y[1], y[2]);

        std::swap(i0[1], i0[3]);
        std::swap(i0[1], i0[2]);
      }

      //std::cout << i0[0] << " " << i0[1] << " " << i0[2] << " " << i0[3] << std::endl;
    }

    std::vector <double> yg(dim, 0);
    for(unsigned k = 0; k < dim; k++) {
      for(unsigned i = 0; i < cnt; i++) {
        yg[k] += y[i][k];
      }
      yg[k] /= cnt;
    }

    a.resize(dim);

    std::vector < std::vector <double> > b(cnt, std::vector<double>(dim));
    for(unsigned k = 0; k < dim; k++) {
      a[k] =  yg[k] - xg[k];
      for(unsigned i = 0; i < cnt; i++) {
        b[i][k] = y[i][k] - xg[k];
      }
    }
    double an = 0.;
    std::vector <double> bn(cnt, 0);
    for(unsigned k = 0; k < dim; k++) {
      an += a[k] * a[k];
      for(unsigned i = 0; i < cnt; i++) {
        bn[i] += b[i][k] * b[i][k];
      }
    }
    an = sqrt(an);
    for(unsigned i = 0; i < cnt; i++) {
      bn[i] = sqrt(bn[i]);
    }

    for(unsigned k = 0; k < dim; k++) {
      a[k] /= an;
      for(unsigned i = 0; i < cnt; i++) {
        b[i][k] /= bn[i];
      }
    }




//     double phig = 0;
//     for(unsigned i = 0; i < cnt; i++) {
//       double phii = 0;
//       for(unsigned k = 0; k < dim; k++) {
//         phii += a[k] * b[i][k];
//       }
//       phii = acos(phii / (an * bn[i]));
//       phig += phii;
//     }
//     phig /= cnt;
    //double H = R * pow(2. * (1. - cos(phig)) / (tan(phig) * tan(phig)), 1. / 3.);

    double H = getHeightPolyhedronSphereInt(b, a, xg, R);
    //else H = R * pow(2. * cos(phig) * cos(phig) / (1. + cos(phig)), 1. / 3.);;



//     if(cnt == 3) H = getHeightPolyhedronSphereInt(bn, a, R);
//     else H = R * pow(2. * cos(phig) * cos(phig) / (1. + cos(phig)), 1. / 3.);

    // H = R * pow(2. * cos(phig) * cos(phig) / (1. + cos(phig)), 1. / 3.);


    xm.resize(dim);
    for(unsigned k = 0; k < dim; k++) {
      xm[k] = xg[k] + a[k] * H;
    }

    // std::cout << "\nBBB " << H << " " << a[0] << " " << a[1] << " " << a[2] << " \n" << xm[0] << " " << xm[1] << " " << xm[2] << std::endl;


    /*
        if(theta[0] > theta[1]) {
          std::swap(theta[0], theta[1]);
        }
        double DT = theta[1] - theta[0];
        if(DT > M_PI) {
          std::swap(theta[0], theta[1]);
          theta[1] += 2. * M_PI;
          DT = theta[1] - theta[0];
        }
        xm.resize(dim);

        d = R * sqrt(0.5 * DT / tan(0.5 * DT)) ;
        a.resize(dim);
        a[0] = -cos(theta[0] + 0.5 * DT);
        a[1] = -sin(theta[0] + 0.5 * DT);

        for(unsigned k = 0; k < dim; k++) {
          xm[k] = -a[k] * d + xg[k];
        }
        d += - a[0] * xg[0] - a[1] * xg[1]; //TODO*/

    //std::cout << "xm = " << xm[0] << " " << xm[1] << std::endl;
    //std::cout << "a = " << a[0] << " b = " << a[1] << " d = " << d << std::endl;
    /*
        double d2 = sqrt(pow(xm[0] - xg[0], 2) + pow(xm[1] - xg[1], 2));
        d = d2 * tan(0.5 * DT);*/

    std::cout.precision(14);

    std::vector<double> xi(dim);

    std::vector < std::vector < double > > J(3, std::vector<double>(3));
    //std::vector < std::vector < double > > JI(3, std::vector<double>(3));
    J[0][0] = (-x1 + x2);
    J[0][1] = (-x1 + x3);
    J[0][2] = (-x1 + x4);

    J[1][0] = (-y1 + y2);
    J[1][1] = (-y1 + y3);
    J[1][2] = (-y1 + y4);

    J[2][0] = (-z1 + z2);
    J[2][1] = (-z1 + z3);
    J[2][2] = (-z1 + z4);

    double den =   J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
                   - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
                   + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

    volume = den / 6.;

    xi[0] = -(x3 * y4 * z1 - x3 * xm[1] * z1 - x1 * y4 * z3 + x1 * xm[1] * z3 - x3 * y1 * z4 + x1 * y3 * z4 - x1 * xm[1] * z4 + x3 * xm[1] * z4 +
              xm[0] * (y3 * z1 - y4 * z1 - y1 * z3 + y4 * z3 + y1 * z4 - y3 * z4) +
              x3 * y1 * xm[2] - x1 * y3 * xm[2] + x1 * y4 * xm[2] - x3 * y4 * xm[2] +
              x4 * (xm[1] * z1 + y1 * z3 - xm[1] * z3 - y1 * xm[2] + y3 * (-z1 + xm[2]))) / den;

    xi[1] = -(-(x2 * y4 * z1) + x2 * xm[1] * z1 + x1 * y4 * z2 - x1 * xm[1] * z2 + x2 * y1 * z4 - x1 * y2 * z4 + x1 * xm[1] * z4 - x2 * xm[1] * z4 +
              xm[0] * (-(y2 * z1) + y4 * z1 + y1 * z2 - y4 * z2 - y1 * z4 + y2 * z4) +
              (-(x2 * y1) + x1 * y2 - x1 * y4 + x2 * y4) * xm[2] +
              x4 * (-(xm[1] * z1) - y1 * z2 + xm[1] * z2 + y2 * (z1 - xm[2]) + y1 * xm[2])) / den;


    xi[2] = -(x2 * y3 * z1 - x2 * xm[1] * z1 - x1 * y3 * z2 + x1 * xm[1] * z2 - x2 * y1 * z3 + x1 * y2 * z3 - x1 * xm[1] * z3 + x2 * xm[1] * z3 +
              xm[0] * (y2 * z1 - y3 * z1 - y1 * z2 + y3 * z2 + y1 * z3 - y2 * z3) +
              x2 * y1 * xm[2] - x1 * y2 * xm[2] + x1 * y3 * xm[2] - x2 * y3 * xm[2] +
              x3 * (xm[1] * z1 + y1 * z2 - xm[1] * z2 - y1 * xm[2] + y2 * (-z1 + xm[2]))) / den;


    // std::cout << "AAA " << xi[0] << " " << xi[1] << " " << xi[2] << " " << den << std::endl;


    a2.assign(dim, 0);
    for(unsigned k = 0; k < dim; k++) {
      for(unsigned j = 0; j < dim; j++) {
        a2[k] -= J[j][k] * a[j]; // this normal has to point toward the center of the ball, thus -=
      }
    }
    double bNorm = sqrt(a2[0] * a2[0] + a2[1] * a2[1] + a2[2] * a2[2]);
    a2[0] /= bNorm;
    a2[1] /= bNorm;
    a2[2] /= bNorm;
    d2 = - a2[0] * xi[0] - a2[1] * xi[1] - a2[2] * xi[2];


    // std::cout << a2[0] << " " << a2[1] << " " << a2[2] << " " << d2 << " " << std::endl;
  }
}

double getHeightPolyhedronSphereInt(const std::vector < std::vector <double> > &b, const std::vector <double> &a, const std::vector <double> &xg, const double &R) {
  const unsigned& cnt = b.size();
  if(b.size() < 3) {
    abort();
  }
  const unsigned& dim = b[0].size();

  std::vector < std::vector <double> > v(cnt, std::vector<double>(dim, 0.));
  for(unsigned i = 0; i < cnt; i++) {
    unsigned ip1 = (i + 1) % cnt;
    for(unsigned k = 0; k < dim; k++) {
      v[i][k] += b[ip1][k] - b[i][k];
    }
  }


  double S = - M_PI * (cnt - 2u);
  for(unsigned i = 0; i < cnt; i++) {
    double dotf = 0.;
    double dotb = 0.;
    unsigned im1 = (cnt + i - 1u) % cnt;
    for(unsigned k = 0; k < dim; k++) {
      dotf += v[i][k] * b[i][k];
      dotb += v[im1][k] * b[i][k];
    }
    double PfdotPb = 0.;
    double normPf = 0.;
    double normPb = 0.;
    for(unsigned k = 0; k < dim; k++) {
      double pf = v[i][k] - dotf * b[i][k];
      double pb = - v[im1][k] + dotb * b[i][k];
      PfdotPb += pf * pb;
      normPf += pf * pf;
      normPb += pb * pb;
    }
    normPf = sqrt(normPf);
    normPb = sqrt(normPb);
    S += acos(PfdotPb / (normPf * normPb));

  }

  std::vector < std::vector <double> > x(cnt, xg);

  for(unsigned i = 0; i < cnt; i++) {
    double h = 0.;
    for(unsigned k = 0; k < dim; k++) {
      h += b[i][k] * a[k];
    }
    h = 1. / h;
    for(unsigned k = 0; k < dim; k++) {
      x[i][k] += h * b[i][k];
    }
  }
  for(unsigned i = 1; i < cnt; i++) {
    for(unsigned k = 0; k < dim; k++) {
      x[i][k] -= x[0][k];
    }
  }
  x[0] = {0., 0., 0.};

  double A = 0.;
  for(unsigned i = 1; i < cnt - 1; i++) {
    A += 0.5 * sqrt((x[i][1] * x[i + 1][2] - x[i][2] * x[i + 1][1]) * (x[i][1] * x[i + 1][2] - x[i][2] * x[i + 1][1]) +
                    (x[i][2] * x[i + 1][0] - x[i][0] * x[i + 1][2]) * (x[i][2] * x[i + 1][0] - x[i][0] * x[i + 1][2]) +
                    (x[i][0] * x[i + 1][1] - x[i][1] * x[i + 1][0]) * (x[i][0] * x[i + 1][1] - x[i][1] * x[i + 1][0]));
  }

  return R * pow(S / A, 1. / 3.);
}

double CurvatureQuadric(const std::vector<double> &a, const std::vector<double> &xp) {
  double k = (8 * a[0] * a[1] * a[1] * xp[1] * xp[1] + 2 * a[1] * ((a[3] + 2 * a[0] * xp[0]) * (a[3] + 2 * a[0] * xp[0]) + 4 * a[0] * (a[4] + a[2] * xp[0]) * xp[1] - a[2] * a[2] * xp[1] * xp[1]) - 2 * (a[4] + a[2] * xp[0]) * (-a[0] * a[4] + a[2] * (a[3] + a[0] * xp[0] + a[2] * xp[1]))) / pow(((a[4] + a[2] * xp[0] + 2 * a[1] * xp[1]) * (a[4] + a[2] * xp[0] + 2 * a[1] * xp[1]) + (a[3] + 2 * a[0] * xp[0] + a[2] * xp[1]) * (a[3] + 2 * a[0] * xp[0] + a[2] * xp[1])), 3. / 2.);
  return k;
}

void NormalQuadric(const std::vector<double> &a, const std::vector<double> &xp, std::vector<double> &N) {
  N.resize(xp.size());

  N[0] = ((a[3] + 2 * a[0] * xp[0] + a[2] * xp[1]) * (8 * a[0] * a[1] * a[1] * xp[1] * xp[1] +
          2 * a[1] * ((a[3] + 2 * a[0] * xp[0]) * (a[3] + 2 * a[0] * xp[0]) + 4 * a[0] * (a[4] + a[2] * xp[0]) * xp[1] - a[2] * a[2] * xp[1] * xp[1]) -
          2 * (a[4] + a[2] * xp[0]) * (-a[0] * a[4] + a[2] * (a[3] + a[0] * xp[0] + a[2] * xp[1])))) / (pow((a[4] + a[2] * xp[0] +
              2 * a[1] * xp[1]) * (a[4] + a[2] * xp[0] + 2 * a[1] * xp[1]) + (a[3] + 2 * a[0] * xp[0] + a[2] * xp[1]) * (a[3] + 2 * a[0] * xp[0] + a[2] * xp[1]), 2));

  N[1] = ((a[4] + a[2] * xp[0] + 2 * a[1] * xp[1]) * (8 * a[0] * a[1] * a[1] * xp[1] * xp[1] +
          2 * a[1] * ((a[3] + 2 * a[0] * xp[0]) * (a[3] + 2 * a[0] * xp[0]) + 4 * a[0] * (a[4] + a[2] * xp[0]) * xp[1] - a[2] * a[2] * xp[1] * xp[1]) -
          2 * (a[4] + a[2] * xp[0]) * (-a[0] * a[4] + a[2] * (a[3] + a[0] * xp[0] + a[2] * xp[1])))) / (pow((a[4] + a[2] * xp[0] +
              2 * a[1] * xp[1]) * (a[4] + a[2] * xp[0] + 2 * a[1] * xp[1]) + (a[3] + 2 * a[0] * xp[0] + a[2] * xp[1]) * (a[3] + 2 * a[0] * xp[0] + a[2] * xp[1]), 2));
}

/*

        for(unsigned k = 0; k < dim; k++) {

          laplce = 0.;
          nonLinear = 0.;


          if(i % 3  == 0) {
            if(k == 0) {
              laplce += 2 * phiV_x[i * dim + k] * phiV_x[j * dim + k] + phiV_x[i * dim + k + 1] * phiV_x[j * dim + k + 1] + phiV_x[i * dim + k + 2] * phiV_x[j * dim + k + 2];
              nonLinear += (solV_gss[0] * phiV_x[j * dim + k] + solV_gss[1] * phiV_x[j * dim + k + 1] + solV_gss[2] * phiV_x[j * dim + k + 2] + gradSolV_gss[k][k]) * phiV[i];
            }

            else if(k == 1) {
              laplce += phiV_x[i * dim + k] * phiV_x[j * dim + k];
              nonLinear += phiV[i] * phiV[j] * gradSolV_gss[0][k];

            }

            else if(k == 2) {
              laplce += phiV_x[i * dim + k] * phiV_x[j * dim + k];
              nonLinear += phiV[i] * phiV[j] * gradSolV_gss[0][k];

            }
          }

          if(i % 3  == 1) {
            if(k == 1) {
              laplce += 2 * phiV_x[i * dim + k] * phiV_x[j * dim + k] + phiV_x[i * dim + k + 1] * phiV_x[j * dim + k + 1] + phiV_x[i * dim + k - 2] * phiV_x[j * dim + k - 2];
              nonLinear += (solV_gss[0] * phiV_x[j * dim + k] + solV_gss[1] * phiV_x[j * dim + k + 1] + solV_gss[2] * phiV_x[j * dim + k + 2] + gradSolV_gss[k][k]) * phiV[i];
            }

            else if(k == 0) {
              laplce += phiV_x[i * dim + k] * phiV_x[j * dim + k + 1];
              nonLinear += phiV[i] * phiV[j] * gradSolV_gss[1][k];

            }

            else if(k == 2) {
              laplce += phiV_x[i * dim + k] * phiV_x[j * dim + k - 1];
              nonLinear += phiV[i] * phiV[j] * gradSolV_gss[1][k];

            }
          }

          if(i % 3  == 2) {
            if(k == 2) {
              laplce += 2 * phiV_x[i * dim + k] * phiV_x[j * dim + k] + phiV_x[i * dim + k - 1] * phiV_x[j * dim + k - 1] + phiV_x[i * dim + k - 2] * phiV_x[j * dim + k - 2];
              nonLinear += (solV_gss[0] * phiV_x[j * dim + k] + solV_gss[1] * phiV_x[j * dim + k + 1] + solV_gss[2] * phiV_x[j * dim + k + 2] + gradSolV_gss[k][k]) * phiV[i];
            }

            else if(k == 0) {
              laplce += phiV_x[i * dim + k] * phiV_x[j * dim + k + 2];
              nonLinear += phiV[i] * phiV[j] * gradSolV_gss[2][k];

            }

            else if(k == 1) {
              laplce += phiV_x[i * dim + k] * phiV_x[j * dim + k + 1];
              nonLinear += phiV[i] * phiV[j] * gradSolV_gss[2][k];

            }
          }

        }

        Jac[3 * i * nDofsV +  i * nDofsP + j] = nu * laplce + nonLinear;
        Jac[(3 * i + 1) * nDofsV +  i * nDofsP + j] = nu * laplce + nonLinear;
        Jac[(3 * i + 2) * nDofsV +  i * nDofsP + j] = nu * laplce + nonLinear;

        Jac[3 * nDofsV * (nDofsV + i) +  nDofsV * nDofsP + i * nDofsP + j] = nu * laplce + nonLinear;
        Jac[3 * nDofsV * (nDofsV + i) + nDofsV + nDofsV * nDofsP + i * nDofsP + j] = nu * laplce + nonLinear;
        Jac[3 * nDofsV * (nDofsV + i) + 2 * nDofsV + nDofsV * nDofsP + i * nDofsP + j] = laplce + nonLinear;

        Jac[6 * nDofsV * (nDofsV + i) + 2 * nDofsV * nDofsP + i * nDofsP + j] = nu * laplce + nonLinear;
        Jac[6 * nDofsV * (nDofsV + i) + nDofsV + 2 * nDofsV * nDofsP + i * nDofsP + j] = nu * laplce + nonLinear;
        Jac[6 * nDofsV * (nDofsV + i) + 2 * nDofsV + 2 * nDofsV * nDofsP + i * nDofsP + j] = nu * laplce + nonLinear;


      }
    }*/
//   }
// }
// double minusP = 0.;
// double minusV = 0.;
//
// for(unsigned i = 0; i < nDofsV; i++) {
//   for(unsigned j = 0; j < nDofsP; j++) {
//     minusP = 0;
//     for(unsigned  J = 0; J < dim; J++) {
//
//       minusP += phiV_x[i * dim + J];
//     }
//
//     Jac[3 * (i + 1) * nDofsV +  i * nDofsP + j] = -phiP[j] * minusP;
//     Jac[9 * (i + 1) * nDofsV + 3 * nDofsP * nDofsV + i * nDofsP + j] = -phiP[j] * minusV;
//
//   }
//
// }
// }
// }
// } //end element loop for each process





