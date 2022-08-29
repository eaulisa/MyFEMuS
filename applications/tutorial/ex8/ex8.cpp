/** \file Ex8.cpp
 *  \brief This example shows how to set and solve the weak form
 *   of the Stress-Strain Equation, for a beam
 *
 *  \nabla \cdot \sigma = F
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
#include "adept.h"
#include "petsc.h"
#include "petscmat.h"
#include "PetscMatrix.hpp"

unsigned DIM = 3;

using namespace femus;

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = false;
  value = 0.;

  if((DIM == 2 && facename == 4) || (DIM == 3 && facename == 5)) {   // left boundary condition.
    dirichlet = true;
  }

  return dirichlet;
}


void AssembleBeamAppoximation_AD(MultiLevelProblem& ml_prob);    
void AssembleBeamAppoximation(MultiLevelProblem& ml_prob);

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;


  unsigned nx = 10; 
  unsigned ny = 2; 
  unsigned nz = 1;

  double length = .1;
  double lengthx = 1;

  if(DIM == 2) {
    mlMsh.GenerateCoarseBoxMesh(nx, ny, 0, 0., lengthx, 0., length, 0., 0., QUAD9, "seventh");
  }
  else if(DIM == 3) {
    nz = ny;
    mlMsh.GenerateCoarseBoxMesh(nx, ny, nz, 0., lengthx , 0., length, 0., length,  HEX27, "seventh");
  }

  unsigned dim = mlMsh.GetDimension();

  unsigned numberOfUniformLevels = 2;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // erase all the coarse mesh levels
  mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("DX", LAGRANGE, SECOND);
  mlSol.AddSolution("DY", LAGRANGE, SECOND);
  if(dim == 3) mlSol.AddSolution("DZ", LAGRANGE, SECOND);

  mlSol.Initialize("All");

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("All");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);

  // add system Poisson in mlProb as a Linear Implicit System
  LinearImplicitSystem& system = mlProb.add_system < LinearImplicitSystem > ("Beam");

  // add solution "u" to system
  system.AddSolutionToSystemPDE("DX");
  system.AddSolutionToSystemPDE("DY");

  if(dim == 3) system.AddSolutionToSystemPDE("DZ");

  // attach the assembling function to system
  system.SetAssembleFunction(AssembleBeamAppoximation);

  // initilaize and solve the system
  system.init();
  system.SetOuterSolver(PREONLY);
  system.MGsolve();

  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(true);

  std::vector<std::string> mov_vars;
  mov_vars.push_back("DX");
  mov_vars.push_back("DY");
  if(DIM == 3)
    mov_vars.push_back("DZ");
  vtkIO.SetMovingMesh(mov_vars);


  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  return 0;
}


void AssembleBeamAppoximation(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled
  //  extract pointers to the several objects that we are going to use

  LinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<LinearImplicitSystem> ("Beam");   // pointer to the linear implicit system named "Beam"
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
  std::vector < unsigned > solDIndex(dim);
  
  solDIndex[0] = mlSol->GetIndex("DX");    // get the position of "X" in the ml_sol object
  solDIndex[1] = mlSol->GetIndex("DY");    // get the position of "Y" in the ml_sol object
  if(dim == 3) solDIndex[2] = mlSol->GetIndex("DZ");       // get the position of "Z" in the ml_sol object

  unsigned solDType = mlSol->GetSolutionType(solDIndex[0]);    // get the finite element type for "D"
  std::vector < unsigned > solDPdeIndex(dim);
  
  solDPdeIndex[0] = mlPdeSys->GetSolPdeIndex("DX");    // get the position of "X" in the pdeSys object
  solDPdeIndex[1] = mlPdeSys->GetSolPdeIndex("DY");    // get the position of "Y" in the pdeSys object
  if(dim == 3) solDPdeIndex[2] = mlPdeSys->GetSolPdeIndex("DZ");   // get the position of "Z" in the pdeSys object

  std::vector < std::vector < double > >  solD(dim);    // local solution
  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector <double> phiD;  // local test function for displacement
  std::vector <double> phiD_x; // local test function first order partial derivatives

  double weight; // gauss point weight

  std::vector< unsigned > sysDof; // local to global pdeSys dofs
  std::vector< double > Res; // local redidual std::vector
  std::vector < double > Jac;

  RES->zero(); // Set to zero all the entries of the Global Residual std::vector
  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);

    unsigned nDofsD = msh->GetElementDofNumber(iel, solDType);    // number of solution element dofs

    Jac.assign((dim * nDofsD) * (dim * nDofsD), 0.);

    unsigned nDofsDAll = dim * nDofsD;
    
    // resize local arrays
    sysDof.resize(nDofsDAll);
    Res.assign(nDofsDAll, 0.);
    Jac.assign(nDofsDAll * nDofsDAll, 0.);

    for(unsigned  k = 0; k < dim; k++) {
      solD[k].resize(nDofsD);
      coordX[k].resize(nDofsD);
    }

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsD; i++) {
      unsigned solDDof = msh->GetSolutionDof(i, iel, solDType);    // local to global mapping between solution node and solution dof

      for(unsigned  k = 0; k < dim; k++) {
        solD[k][i] = (*sol->_Sol[solDIndex[k]])(solDDof);      // global extraction and local storage for the solution
        sysDof[k * nDofsD + i] = pdeSys->GetSystemDof(solDIndex[k], solDPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }
    
    // local storage of coordinates
    for(unsigned i = 0; i < nDofsD; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // local to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solDType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solDType]->Jacobian(coordX, ig, weight, phiD, phiD_x);

      std::vector < double > solD_gss(dim, 0);
      std::vector < std::vector < double > > gradSolD_gss(dim);

      for(unsigned  k = 0; k < dim; k++) {
        gradSolD_gss[k].assign(dim, 0.);
      }

      for(unsigned i = 0; i < nDofsD; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          solD_gss[k] += solD[k][i] * phiD[i];
        }
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            gradSolD_gss[k][j] += solD[k][i] * phiD_x[i * dim + j];
          }
        }
      }

      double E = 10.e6;
      double nu = 0.35;
      double rho = 1000.;

      double mu = E / (2. * (1. + nu));
      double lambda =  E * nu / ((1. + nu) * (1. - 2.*nu));

      // *** phiD_i loop ***
      for(unsigned i = 0; i < nDofsD; i++) {
        for(unsigned  I = 0; I < dim; I++) {  
          double NSV = 0.;
          for(unsigned J = 0; J < dim; J++) {  
            NSV +=  mu * phiD_x[i * dim + J] * (gradSolD_gss[I][J] + gradSolD_gss[J][I]); // diffusion
            NSV +=  lambda * phiD_x[i * dim + I]  * gradSolD_gss[J][J]; //div D
          }
          
          if(I == 1) {
            NSV += phiD[i] * 9.8 * rho; //Gravity only in y - direction
          }
          Res[I * nDofsD + i] -=  NSV * weight;
        }
      } // end phiD_i loop
      
      // *** Formation of Jacobian ***
      for(unsigned i = 0; i < nDofsD; i++) {
        for(unsigned I = 0; I < dim; I++) { //row displacement blocks or dimension
          unsigned VIrow = I * nDofsD + i;
          for(unsigned j = 0; j < nDofsD; j++) {
            unsigned VIcolumn = I * nDofsD + j;
            for(unsigned J = 0; J < dim ; J++) { //column displacement blocks or dimension
              unsigned VJcolumn = J * nDofsD + j;
              Jac[ VIrow * nDofsDAll + VIcolumn] += mu * phiD_x[i * dim + J] * phiD_x[j * dim + J] * weight; //diagonal diffusion
              Jac[ VIrow * nDofsDAll + VJcolumn] += mu * phiD_x[i * dim + J] * phiD_x[j * dim + I] * weight; //off-diagonal diffusion
              Jac[ VIrow * nDofsDAll + VJcolumn] += lambda * phiD_x[i * dim + I] * phiD_x[j * dim + J] * weight; //diagonal div D             
            }
          }       
        }
      }
    }

    // *** Add the local Matrix/Vector into the global Matrix/Vector ***
    RES->add_vector_blocked(Res, sysDof);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);

  } //end element loop for each process

  RES->close();
  KK->close();

  // *** Uncomment to view sparse matrix ***
  
  //VecView ( (static_cast<PetscVector*> (RES))->vec(),  PETSC_VIEWER_STDOUT_SELF);
  //MatView ( (static_cast<PetscMatrix*> (KK))->mat(), PETSC_VIEWER_STDOUT_SELF);

  //PetscViewer    viewer;
  //PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL, NULL, 0, 0, 900, 900, &viewer);
  //PetscObjectSetName((PetscObject) viewer, "PWilmore matrix");
  //PetscViewerPushFormat(viewer, PETSC_VIEWER_DRAW_LG);
  //MatView((static_cast<PetscMatrix*>(KK))->mat(), viewer);
  //double a;
  //std::cin >> a;

}
