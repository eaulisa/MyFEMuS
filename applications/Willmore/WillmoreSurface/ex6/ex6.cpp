/* This example details the full implementation of the p-Willmore flow
   algorithm, which involves three nonlinear systems.

   System0 AssembleInit computes the initial curvatures given mesh positions.
   System AssemblePWillmore solves the flow equations.
   System2 AssembleConformalMinimization "reparametrizes" the surface to
   correct the mesh. */

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

/* Vector option for P (to handle polynomials).
 ap is the coefficient in front of the power of H. */
const unsigned P[3] = {0, 3, 4};
const double ap[3] = {1, 0., 0.};

using namespace femus;

bool O2conformal = true;
bool firstTime = true;
double surface0 = 0.;
double volume0 = 0.;
double willmore0 = 0.;

unsigned counter = 0;

unsigned conformalType0 = 2;
unsigned conformalType = 2;

#include "../include/parameter.hpp"
Parameter parameter = Parameter("cow", 0, true, true, 1, 1, true, 5, 1, 0); //TODO

const double eps = 1.0e-5;

const double normalSign = -1.;

#include "../include/supportFunctions.hpp"
#include "../include/updateMu.hpp"
#include "../include/assembleConformalMinimization.hpp"
#include "../include/assembleInit.hpp"

// Penalty parameter for conformal minimization (eps).
// Trick for system0 (delta). ????
// Trick for system2 (timederiv).
const double delta = 0.00000;
const double timederiv = 0.;

// Declaration of systems.
void AssembleMCF(MultiLevelProblem&);


double dt0 = 500;

// Function to control the time stepping.
double GetTimeStep(const double t) {
  // if(time==0) return 5.0e-7;
  //return 0.0001;
  //double dt0 = .00002;
  // double dt0 = 0.67;
  //double s = 1.;
  //double n = 0.3;
  //return dt0 * pow (1. + t / pow (dt0, s), n);
  return dt0;
}

// IBVs.  No boundary, and IVs set to sphere (just need something).
bool SetBoundaryCondition(const std::vector < double >& x,
                          const char SolName[], double& value, const int facename,
                          const double time) {
  bool dirichlet = false;
  value = 0.;
  return dirichlet;
}

double InitalValueY1(const std::vector < double >& x) {
  return -2. * x[0];
}
double InitalValueY2(const std::vector < double >& x) {
  return -2. * x[1];
}
double InitalValueY3(const std::vector < double >& x) {
  return -2. * x[2];
}


// Main program starts here.
int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  unsigned maxNumberOfMeshes;
  MultiLevelMesh mlMsh;

  // Read coarse level mesh and generate finer level meshes.
  double scalingFactor = 1.; // 1 over the characteristic length

  //mlMsh.ReadCoarseMesh("../input/torus.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("../input/sphere.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("../input/ellipsoidRef3.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("../input/ellipsoidV1.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("../input/genusOne.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("../input/knot.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("../input/cube.neu", "seventh", scalingFactor);
  mlMsh.ReadCoarseMesh("../input/cylinderInBallp75.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("../input/ballTet.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("../input/horseShoe3.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("../input/tiltedTorus.neu", "seventh", scalingFactor);
  scalingFactor = 1.;
  //mlMsh.ReadCoarseMesh ("../input/dog.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("../input/virus3.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("../input/ellipsoidSphere.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("../input/CliffordTorus.neu", "seventh", scalingFactor);

  //mlMsh.ReadCoarseMesh ("../input/spot.med", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("../input/moai.med", "seventh", scalingFactor);

  const bool read_groups = false;                        //by default, if no argument is given, this is "true"
  const bool read_boundary_groups = false;              //by default, if no argument is given, this is "true"
  //mlMsh.ReadCoarseMesh ("../input/moai.med", "seventh", scalingFactor, read_groups, read_boundary_groups);

  // Set number of mesh levels.
  unsigned numberOfUniformLevels = 4;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // Erase all the coarse mesh levels.
  mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  // Define the multilevel solution and attach the mlMsh object to it.
  MultiLevelSolution mlSol(&mlMsh);

  // Add variables X,Y,W to mlSol.
  mlSol.AddSolution("Dx1", LAGRANGE, FIRST, 2);
  mlSol.AddSolution("Dx2", LAGRANGE, FIRST, 2);
  mlSol.AddSolution("Dx3", LAGRANGE, FIRST, 2);

  mlSol.AddSolution("Y1", LAGRANGE, FIRST, 2);
  mlSol.AddSolution("Y2", LAGRANGE, FIRST, 2);
  mlSol.AddSolution("Y3", LAGRANGE, FIRST, 2);

  // Add variable "Lambda" based on constraint choice.
  mlSol.AddSolution("K", LAGRANGE, FIRST, 2);


  mlSol.AddSolution("ENVN", LAGRANGE, FIRST, 0, false);

  // Add variables "nDx" and "Lambda1" for the conformal system.
  mlSol.AddSolution("nDx1", LAGRANGE, FIRST, 0);
  mlSol.AddSolution("nDx2", LAGRANGE, FIRST, 0);
  mlSol.AddSolution("nDx3", LAGRANGE, FIRST, 0);
  mlSol.AddSolution("Lambda1", DISCONTINUOUS_POLYNOMIAL, ZERO, 0);

  mlSol.AddSolution("env", LAGRANGE, FIRST, 0, false);
  mlSol.AddSolution("vAngle", LAGRANGE, FIRST, 0);

  mlSol.AddSolution("mu1", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
  mlSol.AddSolution("mu2", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
  mlSol.AddSolution("weight1", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);

  mlSol.AddSolution("mu1Edge", LAGRANGE, SECOND, 0, false);
  mlSol.AddSolution("mu2Edge", LAGRANGE, SECOND, 0, false);
  mlSol.AddSolution("cntEdge", LAGRANGE, SECOND, 0, false);

  // Initialize the variables and attach boundary conditions.
  mlSol.Initialize("All");

  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("All");

  GetElementNearVertexNumber(mlSol);

  MultiLevelProblem mlProb(&mlSol);

  LinearImplicitSystem& systemY = mlProb.add_system < LinearImplicitSystem > ("InitY");

  // Add solutions Y to systemY.
  systemY.AddSolutionToSystemPDE("Y1");
  systemY.AddSolutionToSystemPDE("Y2");
  systemY.AddSolutionToSystemPDE("Y3");

  // Add the assembling function to system0 and initialize.
  systemY.SetAssembleFunction(AssembleSystemY);
  systemY.init();

  // Add system P-Willmore in mlProb as a time-dependent system.
  TransientNonlinearImplicitSystem& system = mlProb.add_system < TransientNonlinearImplicitSystem > ("MCF");

  // Add solutions X, Y, W to P-Willmore system.
  system.AddSolutionToSystemPDE("Dx1");
  system.AddSolutionToSystemPDE("Dx2");
  system.AddSolutionToSystemPDE("Dx3");
  system.AddSolutionToSystemPDE("K");

  // Parameters for convergence and # of iterations for Willmore.
  system.SetMaxNumberOfNonLinearIterations(1);
  system.SetNonLinearConvergenceTolerance(1.e-10);

  // Attach the assembling function to P-Willmore system.
  system.SetAssembleFunction(AssembleMCF);

  // Attach time step function to P-Willmore sysyem.
  system.AttachGetTimeIntervalFunction(GetTimeStep);

  // Initialize the P-Willmore system.
  system.init();
  system.SetMgType(V_CYCLE);

  // Add system2 Conformal Minimization in mlProb.
  NonLinearImplicitSystem& system2 = mlProb.add_system < NonLinearImplicitSystem > ("conformal");

  // Add solutions newDX, Lambda1 to system2.
  system2.AddSolutionToSystemPDE("nDx1");
  system2.AddSolutionToSystemPDE("nDx2");
  system2.AddSolutionToSystemPDE("nDx3");
  system2.AddSolutionToSystemPDE("Lambda1");

  // Parameters for convergence and # of iterations.
  system2.SetMaxNumberOfNonLinearIterations(2);
  system2.SetNonLinearConvergenceTolerance(1.e-10);

  // Attach the assembling function to system2 and initialize.
  system2.SetAssembleFunction(AssembleConformalMinimization);
  system2.init();

  mlSol.SetWriter(VTK);
  std::vector<std::string> mov_vars;
  mov_vars.push_back("Dx1");
  mov_vars.push_back("Dx2");
  mov_vars.push_back("Dx3");
  mlSol.GetWriter()->SetMovingMesh(mov_vars);

  // and this?
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  // and this?
  mlSol.GetWriter()->SetDebugOutput(false);
  mlSol.GetWriter()->Write("./output1", "linear", variablesToBePrinted, 0);

  // First, solve system2 to "conformalize" the initial mesh.
  CopyDisplacement(mlSol, true);
  //system2.MGsolve();

  // Then, solve system0 to compute initial curvatures.
  CopyDisplacement(mlSol, false);
  system.CopySolutionToOldSolution();
  //systemY.MGsolve();

  mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "linear", variablesToBePrinted, 0);

  // Parameters for the main algorithm loop.
  unsigned numberOfTimeSteps = 1000u;
  unsigned printInterval = 1u;

  std::fstream fs;
  int iproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
  if(iproc == 0) {
    fs.open("Energy.txt", std::fstream::out);
  }


  // Main algorithm loop.
  for(unsigned time_step = 0; time_step < numberOfTimeSteps; time_step++) {
    system.CopySolutionToOldSolution();
    system.MGsolve();

    double energy = GetPWillmoreEnergy(mlSol);
    double dt = system.GetIntervalTime();
    double time = system.GetTime();
    if(iproc == 0) fs << dt << " " << time << " " << energy << std::endl;

    dt0 *= 1.1;

    if(time_step % 1 == 0) {
      mlSol.GetWriter()->Write("./output1", "linear", variablesToBePrinted, (time_step + 1) / printInterval);

      CopyDisplacement(mlSol, true);

      // if (time_step % 20 ==1) {
    //  system2.MGsolve();
      // }

      CopyDisplacement(mlSol, false);
      system.CopySolutionToOldSolution();
      //systemY.MGsolve();
    }

    if((time_step + 1) % printInterval == 0)
      mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "linear", variablesToBePrinted, (time_step + 1) / printInterval);
  }

  if(iproc == 0) fs.close();
  return 0;
}


// Building the P-Willmore assembly function.
void AssembleMCF(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled

  // Call the adept stack object.
  adept::Stack& s = FemusInit::_adeptStack;

  // Extract pointers to the several objects that we are going to use.
  TransientNonlinearImplicitSystem* mlPdeSys =
    &ml_prob.get_system<TransientNonlinearImplicitSystem> ("MCF");

  // Define level and time variable.
  double dt = mlPdeSys->GetIntervalTime();
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  // Point to the mesh and element objects.
  Mesh *msh = ml_prob._ml_msh->GetLevel(level);
  elem *el = msh->el;

  // Point to mlSol, solution (level), and equation (level) objects.
  MultiLevelSolution *mlSol = ml_prob._ml_sol;
  Solution *sol = ml_prob._ml_sol->GetSolutionLevel(level);
  LinearEquationSolver *pdeSys = mlPdeSys->_LinSolver[level];

  // Point to the global stiffness mtx and residual vectors in pdeSys (level).
  SparseMatrix *KK = pdeSys->_KK;
  NumericVector *RES = pdeSys->_RES;

  // Convenience variables to encode the dimension.
  const unsigned dim = 2;
  const unsigned DIM = 3;

  // Get the process_id (for parallel computation).
  unsigned iproc = msh->processor_id();

  // Extract the solution vector; get solDx positions in the ml_sol object.
  unsigned solDxIndex[DIM];
  solDxIndex[0] = mlSol->GetIndex("Dx1");
  solDxIndex[1] = mlSol->GetIndex("Dx2");
  solDxIndex[2] = mlSol->GetIndex("Dx3");

  // Extract the finite element type for solx.
  unsigned solxType;
  solxType = mlSol->GetSolutionType(solDxIndex[0]);

  // Get positions of solDx in the pdeSys object.
  unsigned solDxPdeIndex[DIM];
  solDxPdeIndex[0] = mlPdeSys->GetSolPdeIndex("Dx1");
  solDxPdeIndex[1] = mlPdeSys->GetSolPdeIndex("Dx2");
  solDxPdeIndex[2] = mlPdeSys->GetSolPdeIndex("Dx3");

  // Define solx and solxOld.
  std::vector < adept::adouble > solx[DIM];
  std::vector < double > solxOld[DIM];



  unsigned solKIndex;
  solKIndex = mlSol->GetIndex("K");
  unsigned solKType;
  solKType = mlSol->GetSolutionType(solKIndex);

  // Get positions of solDx in the pdeSys object.
  unsigned solKPdeIndex;
  solKPdeIndex = mlPdeSys->GetSolPdeIndex("K");

  // Define solx and solxOld.
  std::vector < adept::adouble > solK;

  // Get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC).
  unsigned xType = 2;

  // Local-to-global pdeSys dofs.
  std::vector < unsigned > SYSDOF;

  // Define local residual vectors.
  vector < double > Res;
  std::vector< adept::adouble > aResx[3];
  std::vector< adept::adouble > aResK;

  // Local (column-ordered) Jacobian matrix
  vector < double > Jac;

  // Optional command for debugging.
  /* MatSetOption ( ( static_cast<PetscMatrix*> ( KK ) )->mat(),
      MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE ); */

  // Zero all entries of global matrix and residual vectors
  KK->zero();
  RES->zero();

  // Initialize area, volume, P-Willmore energy.
  double surface = 0.;
  double volume = 0.;
  double willmore = 0.;

  // ELEMENT LOOP: each process loops only on the elements that it owns.
  for(int iel = msh->_elementOffset[iproc];
      iel < msh->_elementOffset[iproc + 1]; iel++) {

    // Number of solution element dofs.
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nxDofs  = msh->GetElementDofNumber(iel, solxType);
    unsigned nKDofs  = msh->GetElementDofNumber(iel, solKType);

    // Resize solution vectors.
    for(unsigned K = 0; K < DIM; K++) {
      solx[K].resize(nxDofs);
      solxOld[K].resize(nxDofs);
    }
    solK.resize(nKDofs);

    // Convenience variable for keeping track of problem size.
    unsigned sizeAll = DIM * nxDofs + nKDofs ;

    // Resize local arrays.
    SYSDOF.resize(sizeAll);
    Res.resize(sizeAll);
    for(unsigned K = 0; K < DIM; K++) {
      aResx[K].assign(nxDofs, 0.);
    }
    aResK.assign(nKDofs, 0.);

    // Loop which handles local storage of global mapping and solution X.
    for(unsigned i = 0; i < nxDofs; i++) {

      // Global-to-local mapping between solution node and solution dof.
      unsigned iDDof = msh->GetSolutionDof(i, iel, solxType);
      unsigned iXDof  = msh->GetSolutionDof(i, iel, xType);

      for(unsigned K = 0; K < DIM; K++) {
        solxOld[K][i] = (*msh->_topology->_Sol[K])(iXDof)
                        + (*sol->_SolOld[solDxIndex[K]])(iDDof);
        solx[K][i] = (*msh->_topology->_Sol[K])(iXDof)
                     + (*sol->_Sol[solDxIndex[K]])(iDDof);

        // Global-to-global mapping between solution node and pdeSys dof.
        SYSDOF[K * nxDofs + i] = pdeSys->GetSystemDof(solDxIndex[K], solDxPdeIndex[K], i, iel);
      }
    }

    for(unsigned i = 0; i < nKDofs; i++) {
      unsigned iKDof = msh->GetSolutionDof(i, iel, solKType);
      solK[i] = (*sol->_Sol[solKIndex])(iKDof);
      SYSDOF[DIM * nxDofs + i] = pdeSys->GetSystemDof(solKIndex, solKPdeIndex, i, iel);
    }


    // Start a new recording of all the operations involving adept variables.
    s.new_recording();

    for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
      int faceIndex = el->GetBoundaryIndex(iel, jface);
      // look for boundary faces
      if(faceIndex == 1) {
        const unsigned faceGeom = msh->GetElementFaceType(iel, jface); //edge
        unsigned fnxDofs = msh->GetElementFaceDofNumber(iel, jface, solxType);
        unsigned fnKDofs = msh->GetElementFaceDofNumber(iel, jface, solKType);
        std::vector  <  adept::adouble > fsolx[DIM];    // A matrix holding the face coordinates rowwise.
        std::vector  <  double > fsolxOld[DIM];    // A matrix holding the face coordinates rowwise.

        for(int I = 0; I < DIM; I++) {
          fsolx[I].resize(fnxDofs);
          fsolxOld[I].resize(fnxDofs);
        }
        for(unsigned i = 0; i < fnxDofs; i++) {
          unsigned inode = msh->GetLocalFaceVertexIndex(iel, jface, i);    // face-to-element local node mapping.
          for(unsigned I = 0; I < DIM; I++) {
            fsolx[I][i] =  solx[I][inode]; // We extract the local coordinates on the face from local coordinates on the element.
            fsolxOld[I][i] =  solxOld[I][inode];
          }
        }

        std::vector  <  adept::adouble > fsolK(fnKDofs);    // A matrix holding the face coordinates rowwise.
        for(unsigned i = 0; i < fnKDofs; i++) {
          unsigned inode = msh->GetLocalFaceVertexIndex(iel, jface, i);    // face-to-element local node mapping.
          fsolK[i] =  solK[inode]; // We extract the local coordinates on the face from local coordinates on the element.
        }


        for(unsigned ig = 0; ig  <  msh->_finiteElement[faceGeom][solxType]->GetGaussPointNumber(); ig++) {

          double fweight = msh->_finiteElement[faceGeom][solxType]->GetGaussWeight(ig);
          const double *fphix = msh->_finiteElement[faceGeom][solxType]->GetPhi(ig);
          const double *fphix_u = msh->_finiteElement[ielGeom][solxType]->GetDPhiDXi(ig);

          const double *fphiK = msh->_finiteElement[faceGeom][solKType]->GetPhi(ig);


          adept::adouble fsolxg[3] = {0., 0., 0.};
          adept::adouble fsolxOldg_u[3] = {0., 0., 0.};

          for(unsigned I=0;I<DIM;I++){
            for(unsigned i=0; i < fnxDofs; i++){
              fsolxg[I] += fphix[i] * fsolx[I][i];
              fsolxOldg_u[I] += fphix_u[i] * fsolxOld[I][i];
            }
          }

          adept::adouble fsolKg = 0.;
          for(unsigned i=0; i < fnKDofs; i++){
            fsolKg += fphiK[i] * fsolK[i];
          }

          adept::adouble length = sqrt(fsolxOldg_u[0] * fsolxOldg_u[0] + fsolxOldg_u[1] * fsolxOldg_u[1] + fsolxOldg_u[2] * fsolxOldg_u[2]) * fweight;

          for(unsigned i = 0; i < fnKDofs; i++){
            unsigned fdof = msh->GetLocalFaceVertexIndex(iel, jface, i);    // face-to-element local node mapping.
            aResK[fdof] -= fphiK[i] * ( ( fsolxg[0] * fsolxg[0] + fsolxg[1] * fsolxg[1] + fsolxg[2] * fsolxg[2] ) - 1. ) * length;
          }

          for(unsigned I = 0; I < DIM; I++){
            for(unsigned i = 0; i < fnxDofs; i++){
              unsigned fdof = msh->GetLocalFaceVertexIndex(iel, jface, i);    // face-to-element local node mapping.
              aResx[I][fdof] -= fsolKg * ( 2. * fsolxg[I] * fphix[i] ) * length;
            }
          }
        }
      }
    }




    // begin GAUSS POINT LOOP
    for(unsigned ig = 0; ig < msh->
        _finiteElement[ielGeom][solxType]->GetGaussPointNumber(); ig++) {

      const double *phix;  // local test function
      const double *phix_uv[dim]; // first order derivatives in (u,v)

      double weight; // gauss point weight

      //Extract Gauss point weight, test functions, and their partial derivatives.
      // "0" is derivative in u, "1" is derivative in v.
      weight = msh->_finiteElement[ielGeom][solxType]->GetGaussWeight(ig);

      phix = msh->_finiteElement[ielGeom][solxType]->GetPhi(ig);
      phix_uv[0] = msh->_finiteElement[ielGeom][solxType]->GetDPhiDXi(ig);
      phix_uv[1] = msh->_finiteElement[ielGeom][solxType]->GetDPhiDEta(ig);

      // Initialize quantities xNew, xOld, at the Gauss points.
      adept::adouble solxNewg[3] = {0., 0., 0.};
      double solxOldg[3] = {0., 0., 0.};

      // Initialize derivatives of x and W (new, middle, old) at the Gauss points.
      adept::adouble solxNew_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      adept::adouble solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      double solxOld_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};

      // Compute the above quantities at the Gauss points.
      for(unsigned K = 0; K < DIM; K++) {

        for(unsigned i = 0; i < nxDofs; i++) {
          solxNewg[K] += phix[i] * solx[K][i];
          solxOldg[K] += phix[i] * solxOld[K][i];
        }


        for(int j = 0; j < dim; j++) {
          for(unsigned i = 0; i < nxDofs; i++) {
            solxNew_uv[K][j]    += phix_uv[j][i] * solx[K][i];
            solx_uv[K][j]    += phix_uv[j][i] * 0.5 * (solx[K][i] + solxOld[K][i]);
            solxOld_uv[K][j] += phix_uv[j][i] * solxOld[K][i];
          }
        }
      }

      // Computing the metric, metric determinant, and area element.
      adept::adouble g[dim][dim] = {{0., 0.}, {0., 0.}};
      for(unsigned i = 0; i < dim; i++) {
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned K = 0; K < DIM; K++) {
            g[i][j] += solxOld_uv[K][i] * solxOld_uv[K][j];
          }
        }
      }
      adept::adouble detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];
      adept::adouble Area = weight * sqrt(detg);

      // Computing the unit normal vector N.
      adept::adouble normal[DIM];
      normal[0] = normalSign * (solxOld_uv[1][0] * solxOld_uv[2][1]
                                - solxOld_uv[2][0] * solxOld_uv[1][1]) / sqrt(detg);
      normal[1] = normalSign * (solxOld_uv[2][0] * solxOld_uv[0][1]
                                - solxOld_uv[0][0] * solxOld_uv[2][1]) / sqrt(detg);
      normal[2] = normalSign * (solxOld_uv[0][0] * solxOld_uv[1][1]
                                - solxOld_uv[1][0] * solxOld_uv[0][1]) / sqrt(detg);

      // Computing the metric inverse
      adept::adouble gi[dim][dim];
      gi[0][0] =  g[1][1] / detg;
      gi[0][1] = -g[0][1] / detg;
      gi[1][0] = -g[1][0] / detg;
      gi[1][1] =  g[0][0] / detg;

      // Computing the "reduced Jacobian" g^{ij}X_j .
      adept::adouble Jir[dim][DIM] = {{0., 0., 0.}, {0., 0., 0.}};
      for(unsigned i = 0; i < dim; i++) {
        for(unsigned J = 0; J < DIM; J++) {
          for(unsigned k = 0; k < dim; k++) {
            Jir[i][J] += gi[i][k] * solxOld_uv[J][k];
          }
        }
      }

      // Initializing tangential gradients of X and W (new, middle, old).
      adept::adouble solxNew_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      adept::adouble solx_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      adept::adouble solxOld_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

      // Computing tangential gradients defined above.
      for(unsigned I = 0; I < DIM; I++) {
        for(unsigned J = 0; J < DIM; J++) {
          for(unsigned i = 0; i < dim; i++) {
            //for (unsigned j = 0; j < dim; j++) {
            //   solxNew_Xtan[I][J] += gi[i][j] * solx_uv[J][j] * solxNew_uv[I][i];
            //   solx_Xtan[I][J] += gi[i][j] * solx_uv[J][j] * solx_uv[I][i];
            //   solxOld_Xtan[I][J] += gi[i][j] * solx_uv[J][j] * solxOld_uv[I][i];
            // }
            solxNew_Xtan[I][J] += solxNew_uv[I][i] * Jir[i][J];
            solx_Xtan[I][J] += solx_uv[I][i] * Jir[i][J];
            solxOld_Xtan[I][J] += solxOld_uv[I][i] * Jir[i][J];
          }
        }
      }

      // Define and compute gradients of test functions for X and W.
      std::vector < adept::adouble > phix_Xtan[DIM];
      for(unsigned J = 0; J < DIM; J++) {
        phix_Xtan[J].assign(nxDofs, 0.);

        for(unsigned inode  = 0; inode < nxDofs; inode++) {
          for(unsigned k = 0; k < dim; k++) {
            // for (unsigned j = 0; j < dim; j++) {
            //   phix_Xtan[J][inode] += g[k][j] * solx_uv[J][j] * phix_uv[k][inode];
            // }
            phix_Xtan[J][inode] += phix_uv[k][inode] * Jir[k][J];
          }

        // adept::adouble term1 = 0;
        // for(unsigned K = 0; K < DIM; K++) {
        //   term1 +=  solxNew_Xtan[J][K] * phix_Xtan[K][inode] * Area;
        // }
        }
        // willmore += term1.value();
      }

      // Implement the curvature equation Y = \Delta X .
      for(unsigned K = 0; K < DIM; K++) {

        // Implement the MCF equation.
        for(unsigned i = 0; i < nxDofs; i++) {
          adept::adouble term = 0.;

          for(unsigned J = 0; J < DIM; J++) {
            term +=  solxNew_Xtan[K][J] * phix_Xtan[J][i];
          }
          aResx[K][i] += ((solxNewg[K] - solxOldg[K]) / dt + term) * phix[i] * Area;
        }

      }

      const double *phiK = msh->_finiteElement[ielGeom][solKType]->GetPhi(ig);
      for(unsigned i = 0; i < nKDofs; i++) {
        aResK[i] += 1.0e-10 * solK[i] * phiK[i] * weight;
      }


      // Compute new surface area, volume, and P-Willmore energy.
      for(unsigned K = 0; K < DIM; K++) {
        volume += normalSign * (solxNewg[K].value() * normal[K].value()) * Area.value();
        //willmore += (Hg.value() * Hg.value()) / 4 * Area.value();
      }
      surface += Area.value();


    } // end GAUSS POINT LOOP.

    //------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector
    //copy the value of the adept::adoube aRes in double Res and store

    for(int K = 0; K < DIM; K++) {
      for(int i = 0; i < nxDofs; i++) {
        Res[ K * nxDofs + i] = -aResx[K][i].value();
      }
    }
    for(int i = 0; i < nKDofs; i++) {
      Res[ DIM * nxDofs + i] = -aResK[i].value();
    }

    RES->add_vector_blocked(Res, SYSDOF);

    // Resize Jacobian.
    Jac.resize(sizeAll * sizeAll);

    // Define the dependent variables.
    for(int K = 0; K < DIM; K++) {
      s.dependent(&aResx[K][0], nxDofs);
    }
    s.dependent(&aResK[0], nKDofs);

    // Define the independent variables.
    for(int K = 0; K < DIM; K++) {
      s.independent(&solx[K][0], nxDofs);
    }
    s.independent(&solK[0], nKDofs);

    // Get the Jacobian matrix (ordered by row).
    s.jacobian(&Jac[0], true);

    KK->add_matrix_blocked(Jac, SYSDOF, SYSDOF);

    s.clear_independents();
    s.clear_dependents();

  } // End ELEMENT LOOP for each process.

  RES->close();
  KK->close();

  //KK->draw();

// Get data from each process running in parallel.
  double surfaceAll;
  MPI_Reduce(&surface, &surfaceAll, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  std::cout << "SURFACE = " << surfaceAll << std::endl;

  double volumeAll;
  MPI_Reduce(&volume, &volumeAll, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  std::cout << "VOLUME = " << volumeAll << std::endl;

  // double willmoreAll;
  // MPI_Reduce(&willmore, &willmoreAll, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  // std::cout << "WILLMORE = " << willmoreAll << std::endl;

  firstTime = false;

} // end AssemblePWillmore.
