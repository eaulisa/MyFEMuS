/** \file Ex13.cpp
 *  \brief This example shows how to set and solve the weak form
 *   of the time dependent Stress-Strain Equation, for a beam with pressure
 *
 *  \nabla \cdot \sigma + \mass*acceleration = F
 *  in a Beam domain (in 2D and 3D) clamped on one end
 *  
 *  \author Eugenio Aulisa
 */
// Idea is to set-up porous media equations to solve for u(disp of the solid matr),
// p(fluid pressure) and c(agent concentration). 

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
double dt = 0.01;

double SetVariableTimeStep(const double time) {
  return dt;
}

// u = 0, P=400;
using namespace femus;

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true;
  value = 0.;

  if(!strcmp(SolName, "P")) {   // left boundary condition.
    value = 400;
  }

  return dirichlet;
}


void NewmarkUpdate(MultiLevelSolution *mlSol);
void NewmarkUpdateWithD(MultiLevelSolution *mlSol);

void AssembleRes(MultiLevelProblem& ml_prob);
void AssembleResD(MultiLevelProblem& ml_prob);

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
  mlSol.AddSolution("DX", LAGRANGE, SECOND, 2); //2 means that we have 2 solutions sol and solOld
  mlSol.AddSolution("DY", LAGRANGE, SECOND, 2);
  if(dim == 3) mlSol.AddSolution("DZ", LAGRANGE, SECOND, 2);

  mlSol.AddSolution("VX", LAGRANGE, SECOND, 2);
  mlSol.AddSolution("VY", LAGRANGE, SECOND, 2);
  if(dim == 3) mlSol.AddSolution("VZ", LAGRANGE, SECOND, 2);

  mlSol.AddSolution("AX", LAGRANGE, SECOND, 2);
  mlSol.AddSolution("AY", LAGRANGE, SECOND, 2);
  if(dim == 3) mlSol.AddSolution("AZ", LAGRANGE, SECOND, 2);

  mlSol.AddSolution("P",  DISCONTINUOUS_POLYNOMIAL, FIRST);

  mlSol.Initialize("All");

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("All");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);

  // add system Poisson in mlProb as a Linear Implicit System
  TransientLinearImplicitSystem& system = mlProb.add_system < TransientLinearImplicitSystem > ("Beam");

  if(UseNewmarkUpdateWithD){
    //add solution "D" to system  
    system.AddSolutionToSystemPDE("DX");
    system.AddSolutionToSystemPDE("DY");
    if(dim == 3) system.AddSolutionToSystemPDE("DZ");
    system.AddSolutionToSystemPDE("P");
    system.SetAssembleFunction(AssembleResD);
  }
  
  else{
    //add solution "A" to system
    system.AddSolutionToSystemPDE("AX");
    system.AddSolutionToSystemPDE("AY");
    if(dim == 3) system.AddSolutionToSystemPDE("AZ");
    system.AddSolutionToSystemPDE("P");
    system.SetAssembleFunction(AssembleRes);
  }

  // attach the assembling function to system
  system.AttachGetTimeIntervalFunction(SetVariableTimeStep);

  // initilaize and solve the system
  system.init();
  system.SetOuterSolver(PREONLY);

  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(true);

  std::vector<std::string> mov_vars;
  mov_vars.push_back("DX");
  mov_vars.push_back("DY");
  if(DIM == 3) mov_vars.push_back("DZ");
    
  vtkIO.SetMovingMesh(mov_vars);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);

  const unsigned int n_timesteps = 150;
  if(UseNewmarkUpdateWithD){
      
    for(unsigned t = 0; t < n_timesteps; t++) {
        system.CopySolutionToOldSolution(); // Copy D, V, and A into DOld, VOld, and AOld, respectively
        system.MGsolve(); //solve for A, using DOld, VOld, and AOld
        NewmarkUpdateWithD(&mlSol); // update for A and V, using DOld, VOld, AOld and D

        vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, t + 1);
    }
  }
  
  else{
      
      for(unsigned t = 0; t < n_timesteps; t++) {
        system.CopySolutionToOldSolution(); // Copy D, V, and A into DOld, VOld, and AOld, respectively
        system.MGsolve(); //solve for A, using DOld, VOld, and AOld
        NewmarkUpdate(&mlSol); // update for D an V, using DOld, VOld, AOld and A
        
        vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, t + 1);
    }
  }
  return 0;
}


//Assemble Residual using A to update D amd V
void AssembleRes(MultiLevelProblem& ml_prob) {
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
  if(dim == 3) solDIndex[2] = mlSol->GetIndex("DZ");       // get the position of "Y" in the ml_sol object

  std::vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("VX");    // get the position of "X" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("VY");    // get the position of "Y" in the ml_sol object
  if(dim == 3) solVIndex[2] = mlSol->GetIndex("VZ");       // get the position of "Y" in the ml_sol object

  //solution variable
  std::vector < unsigned > solAIndex(dim);
  solAIndex[0] = mlSol->GetIndex("AX");    // get the position of "X" in the ml_sol object
  solAIndex[1] = mlSol->GetIndex("AY");    // get the position of "Y" in the ml_sol object
  if(dim == 3) solAIndex[2] = mlSol->GetIndex("AZ");       // get the position of "Y" in the ml_sol object

  unsigned solAType = mlSol->GetSolutionType(solAIndex[0]);    // get the finite element type for "u"

  unsigned solPIndex;
  solPIndex = mlSol->GetIndex("P");    // get the position of "P" in the ml_sol object
  unsigned solPType = mlSol->GetSolutionType(solPIndex);    // get the finite element type for "u"

  std::vector < unsigned > solAPdeIndex(dim);
  solAPdeIndex[0] = mlPdeSys->GetSolPdeIndex("AX");    // get the position of "X" in the pdeSys object
  solAPdeIndex[1] = mlPdeSys->GetSolPdeIndex("AY");    // get the position of "Y" in the pdeSys object
  if(dim == 3) solAPdeIndex[2] = mlPdeSys->GetSolPdeIndex("AZ");

  unsigned solPPdeIndex;
  solPPdeIndex = mlPdeSys->GetSolPdeIndex("P");    // get the position of "P" in the pdeSys object

  std::vector < std::vector < double > >  solDOld(dim);    // local solution
  std::vector < std::vector < double > >  solVOld(dim);    // local solution
  std::vector < std::vector < double > >  solAOld(dim);    // local solution

  std::vector < std::vector < double > >  solA(dim);    // local solution
  std::vector < double >  solP; // local solution

  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector <double> phiA;  // local test function for velocity
  std::vector <double> phiA_x; // local test function first order partial derivatives

  double* phiP; // local test function for the pressure
  double weight; // gauss point weight

  std::vector< unsigned > sysDof; // local to global pdeSys dofs
  std::vector< double > Res; // local redidual std::vector
  std::vector < double > Jac;

  RES->zero(); // Set to zero all the entries of the Global Residual std::vector
  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);

    unsigned nDofsA = msh->GetElementDofNumber(iel, solAType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);    // number of solution element dofs

    unsigned nDofsAll = dim * nDofsA + nDofsP;

    sysDof.resize(nDofsAll);
    Res.assign(nDofsAll, 0.);
    Jac.assign(nDofsAll * nDofsAll, 0.);

    for(unsigned  k = 0; k < dim; k++) {
      solDOld[k].resize(nDofsA);
      solVOld[k].resize(nDofsA);
      solAOld[k].resize(nDofsA);
      solA[k].resize(nDofsA);
      coordX[k].resize(nDofsA);
    }
    solP.resize(nDofsP);


    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsA; i++) {
      unsigned iDof = msh->GetSolutionDof(i, iel, solAType);    // local to global mapping between solution node and solution dof

      for(unsigned  k = 0; k < dim; k++) {
        solDOld[k][i] = (*sol->_SolOld[solDIndex[k]])(iDof);
        solVOld[k][i] = (*sol->_SolOld[solVIndex[k]])(iDof);
        solAOld[k][i] = (*sol->_SolOld[solAIndex[k]])(iDof);

        solA[k][i] = (*sol->_Sol[solAIndex[k]])(iDof);      // global extraction and local storage for the solution
        sysDof[k * nDofsA + i] = pdeSys->GetSystemDof(solAIndex[k], solAPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }

    for(unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);    // local to global mapping between solution node and solution dof
      solP[i] = (*sol->_Sol[solPIndex])(solPDof);      // global extraction and local storage for the solution
      sysDof[dim * nDofsA + i ] = pdeSys->GetSystemDof(solPIndex, solPPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofsA; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // local to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }

// *** Face Gauss point loop (boundary Integral) ***
    for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
      int faceIndex = el->GetBoundaryIndex(iel, jface);
      // look for boundary faces
      if((dim == 2 && faceIndex == 3) || (dim == 3 && faceIndex == 4)) {
        const unsigned faceGeom = msh->GetElementFaceType(iel, jface);
        unsigned faceDofs = msh->GetElementFaceDofNumber(iel, jface, solAType);
        std::vector  < std::vector  <  double> > faceCoordinates(dim);    // A matrix holding the face coordinates rowwise.
        for(int k = 0; k < dim; k++) {
          faceCoordinates[k].resize(faceDofs);
        }
        for(unsigned i = 0; i < faceDofs; i++) {
          unsigned inode = msh->GetLocalFaceVertexIndex(iel, jface, i);    // face-to-element local node mapping.
          for(unsigned k = 0; k < dim; k++) {
            faceCoordinates[k][i] =  coordX[k][inode]; // We extract the local coordinates on the face from local coordinates on the element.
          }
        }
        for(unsigned ig = 0; ig  <  msh->_finiteElement[faceGeom][solAType]->GetGaussPointNumber(); ig++) {
          // We call the method GetGaussPointNumber from the object finiteElement in the mesh object msh.
          std::vector < double > normal;
          msh->_finiteElement[faceGeom][solAType]->JacobianSur(faceCoordinates, ig, weight, phiA, phiA_x, normal);

          std::vector< double > xg(dim, 0.);
          for(unsigned i = 0; i < faceDofs; i++) {
            for(unsigned k = 0; k < dim; k++) {
              xg[k] += phiA[i] * faceCoordinates[k][i]; 
            }
          }
          double tau; // sigma (dot) normal (dot) v
          SetBoundaryCondition(xg, "D", tau, faceIndex, 0.);   // return tau
          // *** phiD loop ***
          for(unsigned i = 0; i < faceDofs; i++) {
            unsigned inode = msh->GetLocalFaceVertexIndex(iel, jface, i);
            for(unsigned k = 0; k < dim; k++) {
              Res[k * nDofsA + inode] -=  phiA[i] * tau * normal[k] * weight;
            }
          }
        }
      }
    }

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solAType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solAType]->Jacobian(coordX, ig, weight, phiA, phiA_x);
      phiP = msh->_finiteElement[ielGeom][solPType]->GetPhi(ig);

      std::vector < double > solA_gss(dim, 0);
      std::vector < std::vector < double > > gradSolD_gss(dim);

      for(unsigned  k = 0; k < dim; k++) {
        gradSolD_gss[k].assign(dim, 0.);
      }

      for(unsigned i = 0; i < nDofsA; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          solA_gss[k] += solA[k][i] * phiA[i];
        }
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            double solD = solDOld[k][i] + dt * solVOld[k][i] + 0.5 * dt * dt * ((1. - 2. * BETA) * solAOld[k][i] + 2. * BETA * solA[k][i]);
            gradSolD_gss[k][j] += solD * phiA_x[i * dim + j];
          }
        }
      }

      double solP_gss = 0;
      for(unsigned i = 0; i < nDofsP; i++) {
        solP_gss += phiP[i] * solP[i];
      }

      double E = 10.e6;
      double nu = 0.4;
      double rho = 1000.;

      double mu = E / (2. * (1. + nu));
      double lambdaI = ((1. + nu) * (1. - 2.*nu)) / (E * nu);

      // *** phiA_i loop ***
      for(unsigned i = 0; i < nDofsA; i++) {
        for(unsigned  I = 0; I < dim; I++) {  //momentum equation in k
          double term = 0.;
          for(unsigned J = 0; J < dim; J++) {  // second index j in each equation
            term +=  mu * phiA_x[i * dim + J] * (gradSolD_gss[I][J] + gradSolD_gss[J][I]); // diffusion
          }
          term +=  -phiA_x[i * dim + I]  * solP_gss; //div D

          if(I == 1) {
            term += phiA[i] * 9.8 * rho; //Gravity only in y - direction
          }
          Res[I * nDofsA + i] -= (rho * phiA[i] * solA_gss[I] + term) * weight;
        }
      } // end phiA_i loop


      double divD = 0.;
      for(int I = 0; I < dim; I++) {
        divD += gradSolD_gss[I][I];
      }
      // *** phiP_i loop ***
      for(unsigned i = 0; i < nDofsP; i++) {
        Res[dim * nDofsA + i] -= phiP[i]  * (-divD - lambdaI * solP_gss) / (dt * dt * BETA) *  weight ; //continuity
      } // end phiP_i loop
      // end gauss point loop


      //--------------------------------------------------------------------------------------------------------
      // Add the local Matrix/Vector into the global Matrix/Vector

      for(unsigned i = 0; i < nDofsA; i++) {
        for(unsigned I = 0; I < dim; I++) { //row velocity blocks or dimension
          unsigned VIrow = I * nDofsA + i;
          for(unsigned j = 0; j < nDofsA; j++) {
            unsigned VIcolumn = I * nDofsA + j;

            Jac[ VIrow * nDofsAll + VIcolumn] += rho * phiA[i] * phiA[j] * weight; // inertia

            for(unsigned J = 0; J < dim ; J++) { //column velocity blocks or dimension
              unsigned VJcolumn = J * nDofsA + j;
              Jac[ VIrow * nDofsAll + VIcolumn] += dt * dt * BETA * mu * phiA_x[i * dim + J] * phiA_x[j * dim + J] * weight; //diagonal diffusion
              Jac[ VIrow * nDofsAll + VJcolumn] += dt * dt * BETA * mu * phiA_x[i * dim + J] * phiA_x[j * dim + I] * weight; //off-diagonal diffusion
            }
          }

          for(unsigned j = 0; j < nDofsP; j++) {
            unsigned Pcolumn = dim * nDofsA + j;
            Jac[VIrow * nDofsAll + Pcolumn] += - phiA_x[i * dim + I] * phiP[j] * weight; //pressure gradient
            Jac[Pcolumn * nDofsAll + VIrow] += - phiA_x[i * dim + I] * phiP[j] * weight; //continuity
          }
        }
      }

      for(unsigned i = 0; i < nDofsP; i++) {
        unsigned Prow = dim * nDofsA + i;
        for(unsigned j = 0; j < nDofsP; j++) {
          unsigned Pcolumn = dim * nDofsA + j;
          Jac[Prow * nDofsAll + Pcolumn] += - lambdaI / (dt * dt * BETA) * phiP[i]  * phiP[j] *  weight;   //continuity
        }
      } // end phiP_i loop
    }// end gauss loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    RES->add_vector_blocked(Res, sysDof);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);

  } //end element loop for each process


  RES->close();
  KK->close();

}


void AssembleResD(MultiLevelProblem& ml_prob) {
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

  std::vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("VX");    // get the position of "X" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("VY");    // get the position of "Y" in the ml_sol object
  if(dim == 3) solVIndex[2] = mlSol->GetIndex("VZ");       // get the position of "Z" in the ml_sol object

  //solution variable
  std::vector < unsigned > solAIndex(dim);
  solAIndex[0] = mlSol->GetIndex("AX");    // get the position of "X" in the ml_sol object
  solAIndex[1] = mlSol->GetIndex("AY");    // get the position of "Y" in the ml_sol object
  if(dim == 3) solAIndex[2] = mlSol->GetIndex("AZ");       // get the position of "Z" in the ml_sol object

  unsigned solDType = mlSol->GetSolutionType(solDIndex[0]);    // get the finite element type for "A"

  unsigned solPIndex;
  solPIndex = mlSol->GetIndex("P");    // get the position of "P" in the ml_sol object
  unsigned solPType = mlSol->GetSolutionType(solPIndex);    // get the finite element type for "P"

  std::vector < unsigned > solDPdeIndex(dim);
  solDPdeIndex[0] = mlPdeSys->GetSolPdeIndex("DX");    // get the position of "X" in the pdeSys object
  solDPdeIndex[1] = mlPdeSys->GetSolPdeIndex("DY");    // get the position of "Y" in the pdeSys object
  if(dim == 3) solDPdeIndex[2] = mlPdeSys->GetSolPdeIndex("DZ");   // get the position of "Z" in the ml_sol object

  unsigned solPPdeIndex;
  solPPdeIndex = mlPdeSys->GetSolPdeIndex("P");    // get the position of "P" in the pdeSys object

  std::vector < std::vector < double > >  solDOld(dim);    // local solution
  std::vector < std::vector < double > >  solVOld(dim);    // local solution
  std::vector < std::vector < double > >  solAOld(dim);    // local solution

  std::vector < std::vector < double > >  solD(dim);    // local solution
  std::vector < double >  solP; // local solution

  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector <double> phiD;  // local test function for velocity
  std::vector <double> phiD_x; // local test function first order partial derivatives

  double* phiP; // local test function for the pressure
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
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);    // number of solution element dofs

    unsigned nDofsAll = dim * nDofsD + nDofsP;

    sysDof.resize(nDofsAll);
    Res.assign(nDofsAll, 0.);
    Jac.assign(nDofsAll * nDofsAll, 0.);

    for(unsigned  k = 0; k < dim; k++) {
      solDOld[k].resize(nDofsD);
      solVOld[k].resize(nDofsD);
      solAOld[k].resize(nDofsD);
      solD[k].resize(nDofsD);
      coordX[k].resize(nDofsD);
    }
    solP.resize(nDofsP);


    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsD; i++) {
      unsigned iDof = msh->GetSolutionDof(i, iel, solDType);    // local to global mapping between solution node and solution dof

      for(unsigned  k = 0; k < dim; k++) {
        solAOld[k][i] = (*sol->_SolOld[solAIndex[k]])(iDof);
        solVOld[k][i] = (*sol->_SolOld[solVIndex[k]])(iDof);
        solDOld[k][i] = (*sol->_SolOld[solDIndex[k]])(iDof);

        solD[k][i] = (*sol->_Sol[solDIndex[k]])(iDof);      // global extraction and local storage for the solution
        sysDof[k * nDofsD + i] = pdeSys->GetSystemDof(solDIndex[k], solDPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }

    for(unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);    // local to global mapping between solution node and solution dof
      solP[i] = (*sol->_Sol[solPIndex])(solPDof);      // global extraction and local storage for the solution
      sysDof[dim * nDofsD + i ] = pdeSys->GetSystemDof(solPIndex, solPPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofsD; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // local to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }

// *** Face Gauss point loop (boundary Integral) ***
    for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
      int faceIndex = el->GetBoundaryIndex(iel, jface);
      // look for boundary faces
      if((dim == 2 && faceIndex == 3) || (dim == 3 && faceIndex == 4)) {
        const unsigned faceGeom = msh->GetElementFaceType(iel, jface);
        unsigned faceDofs = msh->GetElementFaceDofNumber(iel, jface, solDType);
        std::vector  < std::vector  <  double> > faceCoordinates(dim);    // A matrix holding the face coordinates rowwise.
        for(int k = 0; k < dim; k++) {
          faceCoordinates[k].resize(faceDofs);
        }
        for(unsigned i = 0; i < faceDofs; i++) {
          unsigned inode = msh->GetLocalFaceVertexIndex(iel, jface, i);    // face-to-element local node mapping.
          for(unsigned k = 0; k < dim; k++) {
            faceCoordinates[k][i] =  coordX[k][inode]; // We extract the local coordinates on the face from local coordinates on the element.
          }
        }
        for(unsigned ig = 0; ig  <  msh->_finiteElement[faceGeom][solDType]->GetGaussPointNumber(); ig++) {
          // We call the method GetGaussPointNumber from the object finiteElement in the mesh object msh.
          std::vector < double > normal;
          msh->_finiteElement[faceGeom][solDType]->JacobianSur(faceCoordinates, ig, weight, phiD, phiD_x, normal);

          std::vector< double > xg(dim, 0.);
          for(unsigned i = 0; i < faceDofs; i++) {
            for(unsigned k = 0; k < dim; k++) {
              xg[k] += phiD[i] * faceCoordinates[k][i]; // xg(ig)= \sum_{i=0}^faceDofs phi[i](xig) facecoordinates[i]
            }
          }
          double tau; // a(u)*grad_u\cdot normal
          SetBoundaryCondition(xg, "D", tau, faceIndex, 0.);   // return tau
          // *** phiD loop ***
          for(unsigned i = 0; i < faceDofs; i++) {
            unsigned inode = msh->GetLocalFaceVertexIndex(iel, jface, i);
            for(unsigned k = 0; k < dim; k++) {
              Res[k * nDofsD + inode] -=  phiD[i] * tau * normal[k] * weight;
            }
          }
        }
      }
    }

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solDType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solDType]->Jacobian(coordX, ig, weight, phiD, phiD_x);
      phiP = msh->_finiteElement[ielGeom][solPType]->GetPhi(ig);

      std::vector < double > solA_gss(dim, 0);
      std::vector < std::vector < double > > gradSolD_gss(dim);

      for(unsigned  k = 0; k < dim; k++) {
        gradSolD_gss[k].assign(dim, 0.);
      }

      for(unsigned i = 0; i < nDofsD; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          solA_gss[k] += ((solD[k][i] - solDOld[k][i]) / (BETA * dt * dt) - solVOld[k][i] / (BETA * dt) + ((BETA - 0.5) * solAOld[k][i]) / BETA) * phiD[i];
        }
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            gradSolD_gss[k][j] += solD[k][i] * phiD_x[i * dim + j];
          }
        }
      }

      double solP_gss = 0;
      for(unsigned i = 0; i < nDofsP; i++) {
        solP_gss += phiP[i] * solP[i];
      }

      double E = 10.e6;
      double nu = 0.4;
      double rho = 1000.;

      double mu = E / (2. * (1. + nu));
      double lambdaI = ((1. + nu) * (1. - 2.*nu)) / (E * nu);

      // *** phiD_i loop ***
      for(unsigned i = 0; i < nDofsD; i++) {
        for(unsigned  I = 0; I < dim; I++) {  //momentum equation in k
          double term = 0.;
          for(unsigned J = 0; J < dim; J++) {  // second index j in each equation
            term +=  mu * phiD_x[i * dim + J] * (gradSolD_gss[I][J] + gradSolD_gss[J][I]); // diffusion
          }
          term +=  -phiD_x[i * dim + I]  * solP_gss; //div D

          if(I == 1) {
            term += phiD[i] * 9.8 * rho; //Gravity only in y - direction
          }
          Res[I * nDofsD + i] -= (rho * phiD[i] * solA_gss[I] + term) * weight;
        }
      } // end phiD_i loop


      double divD = 0.;
      for(int I = 0; I < dim; I++) {
        divD += gradSolD_gss[I][I];
      }
      // *** phiP_i loop ***
      for(unsigned i = 0; i < nDofsP; i++) {
        Res[dim * nDofsD + i] -= phiP[i]  * (-divD - lambdaI * solP_gss) * weight ; //Continuity with Pressure
      } // end phiP_i loop
      // end gauss point loop

      //--------------------------------------------------------------------------------------------------------
      // Add the local Matrix/Vector into the global Matrix/Vector

      for(unsigned i = 0; i < nDofsD; i++) {
        for(unsigned I = 0; I < dim; I++) { //row velocity blocks or dimension
          unsigned VIrow = I * nDofsD + i;
          for(unsigned j = 0; j < nDofsD; j++) {
            unsigned VIcolumn = I * nDofsD + j;

            Jac[ VIrow * nDofsAll + VIcolumn] += rho / (BETA * dt * dt) * phiD[i] * phiD[j] * weight;   // inertia

            for(unsigned J = 0; J < dim ; J++) { //column velocity blocks or dimension
              unsigned VJcolumn = J * nDofsD + j;
              Jac[ VIrow * nDofsAll + VIcolumn] += mu * phiD_x[i * dim + J] * phiD_x[j * dim + J] * weight; //diagonal diffusion 
              Jac[ VIrow * nDofsAll + VJcolumn] += mu * phiD_x[i * dim + J] * phiD_x[j * dim + I] * weight; //off-diagonal diffusion 
            }
          }

          for(unsigned j = 0; j < nDofsP; j++) {
            unsigned Pcolumn = dim * nDofsD + j;
            Jac[VIrow * nDofsAll + Pcolumn] += - phiD_x[i * dim + I] * phiP[j] * weight; //pressure gradient
            Jac[Pcolumn * nDofsAll + VIrow] += - phiD_x[i * dim + I] * phiP[j] * weight; //continuity
          }
        }
      }

      for(unsigned i = 0; i < nDofsP; i++) {
        unsigned Prow = dim * nDofsD + i;
        for(unsigned j = 0; j < nDofsP; j++) {
          unsigned Pcolumn = dim * nDofsD + j;
          Jac[Prow * nDofsAll + Pcolumn] += - lambdaI * phiP[i]  * phiP[j] *  weight;   //continuity 
        }
      } // end phiP_i loop
    }// end gauss loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    RES->add_vector_blocked(Res, sysDof);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);

  } //end element loop for each process


  RES->close();
  KK->close();
}

//This updates D and V using A
void NewmarkUpdate(MultiLevelSolution *mlSol) {
  unsigned level = mlSol->_mlMesh->GetNumberOfLevels() - 1u;
  //  extract pointers to the several objects that we are going to use
  Mesh*     msh = mlSol->_mlMesh->GetLevel(level);    // pointer to the mesh (level) object
  elem*     el  = msh->el;  // pointer to the elem object in msh (level)
  Solution* sol = mlSol->GetSolutionLevel(level);    // pointer to the solution (level) object

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  std::vector < unsigned > solDIndex(dim);
  solDIndex[0] = mlSol->GetIndex("DX");    // get the position of "X" in the ml_sol object
  solDIndex[1] = mlSol->GetIndex("DY");    // get the position of "Y" in the ml_sol object
  if(dim == 3) solDIndex[2] = mlSol->GetIndex("DZ");       // get the position of "Y" in the ml_sol object

  std::vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("VX");    // get the position of "X" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("VY");    // get the position of "Y" in the ml_sol object
  if(dim == 3) solVIndex[2] = mlSol->GetIndex("VZ");       // get the position of "Y" in the ml_sol object

  std::vector < unsigned > solAIndex(dim);
  solAIndex[0] = mlSol->GetIndex("AX");    // get the position of "X" in the ml_sol object
  solAIndex[1] = mlSol->GetIndex("AY");    // get the position of "Y" in the ml_sol object
  if(dim == 3) solAIndex[2] = mlSol->GetIndex("AZ");       // get the position of "Y" in the ml_sol object

  unsigned solType = mlSol->GetSolutionType(solDIndex[0]);

  //Newmark Update
  for(unsigned k = 0; k < dim; k++) {
    for(unsigned i = msh->_dofOffset[solType][iproc]; i < msh->_dofOffset[solType][iproc + 1]; i++) {
      double Anew = (*sol->_Sol[solAIndex[k]])(i);
      double Aold = (*sol->_SolOld[solAIndex[k]])(i);
      double Vold = (*sol->_SolOld[solVIndex[k]])(i);
      double Dold = (*sol->_SolOld[solDIndex[k]])(i);

      double Vnew = Vold + dt * ((1. - GAMMA) * Aold + GAMMA * Anew);
      double Dnew = Dold + dt * Vold + 0.5 * dt * dt * ((1. - 2. * BETA) * Aold + 2. * BETA * Anew);

      sol->_Sol[solVIndex[k]]->set(i, Vnew);
      sol->_Sol[solDIndex[k]]->set(i, Dnew);
    }
    sol->_Sol[solVIndex[k]]->close();
    sol->_Sol[solDIndex[k]]->close();
  }
}

//Update A and V using D
void NewmarkUpdateWithD(MultiLevelSolution *mlSol) {
  unsigned level = mlSol->_mlMesh->GetNumberOfLevels() - 1u;
  //  extract pointers to the several objects that we are going to use
  Mesh*     msh = mlSol->_mlMesh->GetLevel(level);    // pointer to the mesh (level) object
  elem*     el  = msh->el;  // pointer to the elem object in msh (level)
  Solution* sol = mlSol->GetSolutionLevel(level);    // pointer to the solution (level) object

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  std::vector < unsigned > solDIndex(dim);
  solDIndex[0] = mlSol->GetIndex("DX");    // get the position of "X" in the ml_sol object
  solDIndex[1] = mlSol->GetIndex("DY");    // get the position of "Y" in the ml_sol object
  if(dim == 3) solDIndex[2] = mlSol->GetIndex("DZ");       // get the position of "Z" in the ml_sol object

  std::vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("VX");    // get the position of "X" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("VY");    // get the position of "Y" in the ml_sol object
  if(dim == 3) solVIndex[2] = mlSol->GetIndex("VZ");       // get the position of "Z" in the ml_sol object

  std::vector < unsigned > solAIndex(dim);
  solAIndex[0] = mlSol->GetIndex("AX");    // get the position of "X" in the ml_sol object
  solAIndex[1] = mlSol->GetIndex("AY");    // get the position of "Y" in the ml_sol object
  if(dim == 3) solAIndex[2] = mlSol->GetIndex("AZ");       // get the position of "Z" in the ml_sol object

  unsigned solType = mlSol->GetSolutionType(solDIndex[0]);

  //Newmark Update
  for(unsigned k = 0; k < dim; k++) {
    for(unsigned i = msh->_dofOffset[solType][iproc]; i < msh->_dofOffset[solType][iproc + 1]; i++) {
      double Dnew = (*sol->_Sol[solDIndex[k]])(i);
      double Aold = (*sol->_SolOld[solAIndex[k]])(i);
      double Vold = (*sol->_SolOld[solVIndex[k]])(i);
      double Dold = (*sol->_SolOld[solDIndex[k]])(i);

      double Anew = (Dnew - Dold) / (BETA * dt * dt) - Vold / (BETA * dt) + ((BETA - 0.5) * Aold) / BETA;
      double Vnew = Vold + (1 - GAMMA) * dt * Aold + GAMMA * dt * Anew;

      sol->_Sol[solAIndex[k]]->set(i, Anew);
      sol->_Sol[solVIndex[k]]->set(i, Vnew);
    }
    sol->_Sol[solAIndex[k]]->close();
    sol->_Sol[solVIndex[k]]->close();
  }
}







