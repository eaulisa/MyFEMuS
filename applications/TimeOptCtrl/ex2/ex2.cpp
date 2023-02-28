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
#include "NonLinearImplicitSystem.hpp"
#include "adept.h"


#include "petsc.h"
#include "petscmat.h"
#include "PetscMatrix.hpp"

double beta = 0.00001;
double alpha = 0.00001;

using namespace femus;

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0;

  if(!strcmp(SolName, "bV")) {
    if(x[0] < 0. && x[1] < 0.5 && x[1] > -0.5 && x[2] < 0.5 && x[2] > -0.5) dirichlet = false;
  }
  if(!strcmp(SolName, "Vc")) {
    if(x[0] < 0. && x[1] < 0.5 && x[1] > -0.5 && x[2] < 0.5 && x[2] > -0.5) value = 1;
  }
  else if(!strcmp(SolName, "bP") || !strcmp(SolName, "lP") || !strcmp(SolName, "Pc")) {
    dirichlet = false;
    value = 0.;
  }

  return dirichlet;
}


void AssembleSteadyStateControl(MultiLevelProblem& ml_prob);    //, unsigned level, const unsigned &levelMax, const bool &assembleMatrix );
void AssembleManifactureSolution(MultiLevelProblem& ml_prob);    //, unsigned level, const unsigned &levelMax, const bool &assembleMatrix );

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  //mlMsh.ReadCoarseMesh("./input/cube_hex.neu", "seventh", scalingFactor);
  mlMsh.ReadCoarseMesh("./input/square_quad.neu", "seventh", scalingFactor);
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
     probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();

  unsigned numberOfUniformLevels = 5;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels, numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // erase all the coarse mesh levels
  mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("Uc", LAGRANGE, SECOND);
  mlSol.AddSolution("Vc", LAGRANGE, SECOND);
  mlSol.AddSolution("Pc",  DISCONTINUOUS_POLYNOMIAL, FIRST);


  // add variables to mlSol
  mlSol.AddSolution("bU", LAGRANGE, SECOND);
  mlSol.AddSolution("bV", LAGRANGE, SECOND);
  mlSol.AddSolution("bP",  DISCONTINUOUS_POLYNOMIAL, FIRST);


  mlSol.AddSolution("lU", LAGRANGE, SECOND);
  mlSol.AddSolution("lV", LAGRANGE, SECOND);
  mlSol.AddSolution("lP",  DISCONTINUOUS_POLYNOMIAL, FIRST);

  mlSol.Initialize("All");

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.FixSolutionAtOnePoint("Pc");
  mlSol.FixSolutionAtOnePoint("bP");
  mlSol.FixSolutionAtOnePoint("lP");

  mlSol.GenerateBdc("All");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);

  // add system Poisson in mlProb as a Linear Implicit System
  LinearImplicitSystem& systemC = mlProb.add_system < LinearImplicitSystem > ("ManSol");
  
  systemC.AddSolutionToSystemPDE("Uc");
  systemC.AddSolutionToSystemPDE("Vc");
  systemC.AddSolutionToSystemPDE("Pc");

  // attach the assembling function to system
  systemC.SetAssembleFunction(AssembleManifactureSolution);

  // initilaize and solve the system
  systemC.init();
  systemC.SetOuterSolver(PREONLY);
  systemC.MGsolve();
  
  
  LinearImplicitSystem& system = mlProb.add_system < LinearImplicitSystem > ("NS");

  // add solution "u" to system

  system.AddSolutionToSystemPDE("bU");
  system.AddSolutionToSystemPDE("bV");
  system.AddSolutionToSystemPDE("bP");


  system.AddSolutionToSystemPDE("lU");
  system.AddSolutionToSystemPDE("lV");
  system.AddSolutionToSystemPDE("lP");

  // attach the assembling function to system
  system.SetAssembleFunction(AssembleSteadyStateControl);

  // initilaize and solve the system
  system.init();
  system.SetOuterSolver(PREONLY);
  system.MGsolve();

  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(true);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 1);

  return 0;
}


void AssembleSteadyStateControl(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

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

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  
  std::vector < unsigned > solVcIndex(dim);
  solVcIndex[0] = mlSol->GetIndex("Uc");
  solVcIndex[1] = mlSol->GetIndex("Vc");
  
  //solution variable
  std::vector < unsigned > solbVIndex(dim);
  std::vector < unsigned > sollVIndex(dim);
  solbVIndex[0] = mlSol->GetIndex("bU");
  solbVIndex[1] = mlSol->GetIndex("bV");
  sollVIndex[0] = mlSol->GetIndex("lU");
  sollVIndex[1] = mlSol->GetIndex("lV");

  unsigned solVType = mlSol->GetSolutionType(solbVIndex[0]);

  unsigned solbPIndex;
  unsigned sollPIndex;
  solbPIndex = mlSol->GetIndex("bP");
  sollPIndex = mlSol->GetIndex("lP");
  unsigned solPType = mlSol->GetSolutionType(solbPIndex);

  std::vector < unsigned > solbVPdeIndex(dim);
  std::vector < unsigned > sollVPdeIndex(dim);
  solbVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("bU");
  solbVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("bV");
  sollVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("lU");
  sollVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("lV");


  unsigned solbPPdeIndex;
  unsigned sollPPdeIndex;
  solbPPdeIndex = mlPdeSys->GetSolPdeIndex("bP");
  sollPPdeIndex = mlPdeSys->GetSolPdeIndex("lP");

  
  std::vector < std::vector < double > >  solVc(dim);
  
  std::vector < std::vector < adept::adouble > >  solbV(dim);
  std::vector < std::vector < adept::adouble > >  sollV(dim);
  std::vector < adept::adouble >  solbP;
  std::vector < adept::adouble >  sollP;

  std::vector< std::vector < adept::adouble > > mResbV(dim);
  std::vector< std::vector < adept::adouble > > mReslV(dim);
  std::vector< adept::adouble > mResbP;
  std::vector< adept::adouble > mReslP;

  std::vector < std::vector < double > > x(dim);
  unsigned xType = 2;

  std::vector <double> phiV;
  std::vector <double> phiVx;

  double* phiP;
  double weight;

  std::vector< unsigned > sysDof;
  std::vector< double > Res;
  std::vector < double > Jac;

  RES->zero();
  KK->zero();

  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielType = msh->GetElementType(iel);

    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);

    unsigned nDofsVP = dim * nDofsV + nDofsP;
    sysDof.resize(2 * nDofsVP);

    for(unsigned  k = 0; k < dim; k++) {
      solVc[k].resize(nDofsV);  
      solbV[k].resize(nDofsV);
      sollV[k].resize(nDofsV);
      x[k].resize(nDofsV);
    }
    solbP.resize(nDofsP);
    sollP.resize(nDofsP);

    for(unsigned  k = 0; k < dim; k++) {
      mResbV[k].assign(nDofsV, 0.);
      mReslV[k].assign(nDofsV, 0.);
    }
    mResbP.assign(nDofsP, 0.);
    mReslP.assign(nDofsP, 0.);

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);

      for(unsigned  k = 0; k < dim; k++) {
        
        solVc[k][i] = (*sol->_Sol[solVcIndex[k]])(solVDof);  
          
        solbV[k][i] = (*sol->_Sol[solbVIndex[k]])(solVDof);
        sollV[k][i] = (*sol->_Sol[sollVIndex[k]])(solVDof);

        sysDof[k * nDofsV + i] = pdeSys->GetSystemDof(solbVIndex[k], solbVPdeIndex[k], i, iel);
        sysDof[nDofsVP + k * nDofsV + i] = pdeSys->GetSystemDof(sollVIndex[k], sollVPdeIndex[k], i, iel);
      }
    }

    for(unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);
      solbP[i] = (*sol->_Sol[solbPIndex])(solPDof);
      sollP[i] = (*sol->_Sol[sollPIndex])(solPDof);
      sysDof[dim * nDofsV + i ] = pdeSys->GetSystemDof(solbPIndex, solbPPdeIndex, i, iel);
      sysDof[nDofsVP + dim * nDofsV + i ] = pdeSys->GetSystemDof(sollPIndex, sollPPdeIndex, i, iel);
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // local to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k])(xDof);      // global extraction and local storage for the element coordinates
      }
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielType][solVType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielType][solVType]->Jacobian(x, ig, weight, phiV, phiVx);
      phiP = msh->_finiteElement[ielType][solPType]->GetPhi(ig);

      std::vector < double > solVcg(dim, 0);
      
      std::vector < adept::adouble > solbVg(dim, 0);
      std::vector < adept::adouble > sollVg(dim, 0);
      std::vector < double > xg(dim, 0);

      std::vector < std::vector < adept::adouble > > solbVxg(dim, std::vector < adept::adouble >(dim, 0.));
      std::vector < std::vector < adept::adouble > > sollVxg(dim, std::vector < adept::adouble >(dim, 0.));

      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned  k = 0; k < dim; k++) {
            
          solVcg[k] += solVc[k][i] * phiV[i];  
          solbVg[k] += solbV[k][i] * phiV[i];
          sollVg[k] += sollV[k][i] * phiV[i];
          xg[k] += x[k][i] * phiV[i];
        }
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            solbVxg[k][j] += solbV[k][i] * phiVx[i * dim + j];
            sollVxg[k][j] += sollV[k][i] * phiVx[i * dim + j];
          }
        }
      }

      //std::vector < double > Vc = {xg[1], -xg[0]};
     
      std::vector < double > Vc = solVcg;


      adept::adouble solbPg = 0;
      adept::adouble sollPg = 0;
      for(unsigned i = 0; i < nDofsP; i++) {
        solbPg += phiP[i] * solbP[i];
        sollPg += phiP[i] * sollP[i];
      }
      double nu = 1.;

      // *** phiV_i loop ***
      for(unsigned i = 0; i < nDofsV; i++) {
        std::vector < adept::adouble > NSVb(dim, 0.);
        std::vector < adept::adouble > NSVl(dim, 0.);

        for(unsigned  k = 0; k < dim; k++) {  //momentum equation in k
          for(unsigned j = 0; j < dim; j++) {  // second index j in each equation
            NSVb[k]   +=  nu * phiVx[i * dim + j] * (solbVxg[k][j] + solbVxg[j][k]);
            NSVl[k]   +=  nu * phiVx[i * dim + j] * (sollVxg[k][j] + sollVxg[j][k]) + beta * phiVx[i * dim + j] * (solbVxg[k][j] + solbVxg[j][k]);
          }
          NSVb[k] += -solbPg * phiVx[i * dim + k];
          NSVl[k] += -sollPg * phiVx[i * dim + k]  + alpha * solbVg[k] * phiV[i] + (solbVg[k] - Vc[k]) * phiV[i];
        }
        for(unsigned  k = 0; k < dim; k++) {
          mResbV[k][i] += - NSVb[k] * weight;
          mReslV[k][i] += - NSVl[k] * weight;
        }
      } // end phiV_i loop

      // *** phiP_i loop ***
      for(unsigned i = 0; i < nDofsP; i++) {
        for(int k = 0; k < dim; k++) {
          mResbP[i] += - (-solbVxg[k][k]) * phiP[i]  * weight;
          mReslP[i] += - (-sollVxg[k][k]) * phiP[i]  * weight;
        }
      } // end phiP_i loop
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube mRes in double Res and store them in RES
    Res.resize(2 * nDofsVP);    //resize

    for(int i = 0; i < nDofsV; i++) {
      for(unsigned  k = 0; k < dim; k++) {
        Res[ k * nDofsV + i ] = -mReslV[k][i].value();
        Res[nDofsVP +  k * nDofsV + i ] = -mResbV[k][i].value();
      }
    }

    for(int i = 0; i < nDofsP; i++) {
      Res[ dim * nDofsV + i ] = -mReslP[i].value();
      Res[nDofsVP + dim * nDofsV + i ] = -mResbP[i].value();
    }

    RES->add_vector_blocked(Res, sysDof);

    //Extract and store the Jacobian

    Jac.resize(2 * nDofsVP * 2 * nDofsVP);
    // define the dependent variables


    for(unsigned  k = 0; k < dim; k++) {
      s.dependent(&mReslV[k][0], nDofsV);
    }
    s.dependent(&mReslP[0], nDofsP);

    for(unsigned  k = 0; k < dim; k++) {
      s.dependent(&mResbV[k][0], nDofsV);
    }
    s.dependent(&mResbP[0], nDofsP);






    // define the independent variables
    for(unsigned  k = 0; k < dim; k++) {
      s.independent(&solbV[k][0], nDofsV);
    }
    s.independent(&solbP[0], nDofsP);
    // define the independent variables
    for(unsigned  k = 0; k < dim; k++) {
      s.independent(&sollV[k][0], nDofsV);
    }
    s.independent(&sollP[0], nDofsP);


    // get the and store jacobian matrix (row-major)
    s.jacobian(&Jac[0], true);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();
  KK->close();

// KK->draw();
  // ***************** END ASSEMBLY *******************

}


void AssembleManifactureSolution(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  LinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<LinearImplicitSystem> ("ManSol");   // pointer to the linear implicit system named "Poisson"
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
  std::vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("Uc");
  solVIndex[1] = mlSol->GetIndex("Vc");

  unsigned solVType = mlSol->GetSolutionType(solVIndex[0]);

  unsigned solPIndex;
  solPIndex = mlSol->GetIndex("Pc");
  unsigned solPType = mlSol->GetSolutionType(solPIndex);

  std::vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("Uc");
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("Vc");
  

  unsigned solPPdeIndex;
  solPPdeIndex = mlPdeSys->GetSolPdeIndex("Pc");

  std::vector < std::vector < adept::adouble > >  solV(dim);
  std::vector < adept::adouble >  solP;
  
  std::vector< std::vector < adept::adouble > > mResV(dim);
  std::vector< adept::adouble > mResP;

  std::vector < std::vector < double > > x(dim);
  unsigned xType = 2;

  std::vector <double> phiV;
  std::vector <double> phiVx;

  double* phiP;
  double weight;

  std::vector< unsigned > sysDof;
  std::vector< double > Res;
  std::vector < double > Jac;

  RES->zero();
  KK->zero();

  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielType = msh->GetElementType(iel);

    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);

    unsigned nDofsVP = dim * nDofsV + nDofsP;
    sysDof.resize(nDofsVP);

    for(unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofsV);;
      x[k].resize(nDofsV);
    }
    solP.resize(nDofsP);

    for(unsigned  k = 0; k < dim; k++) {
      mResV[k].assign(nDofsV, 0.);
    }
    mResP.assign(nDofsP, 0.);
    
    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);

      for(unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);

        sysDof[k * nDofsV + i] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);
        
      }
    }

    for(unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);
      solP[i] = (*sol->_Sol[solPIndex])(solPDof);
      sysDof[dim * nDofsV + i ] = pdeSys->GetSystemDof(solPIndex, solPPdeIndex, i, iel);
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // local to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k])(xDof);      // global extraction and local storage for the element coordinates
      }
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielType][solVType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielType][solVType]->Jacobian(x, ig, weight, phiV, phiVx);
      phiP = msh->_finiteElement[ielType][solPType]->GetPhi(ig);

      std::vector < adept::adouble > solVg(dim, 0);

      std::vector < std::vector < adept::adouble > > solVxg(dim, std::vector < adept::adouble >(dim, 0.));

      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          solVg[k] += solV[k][i] * phiV[i];
        }
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            solVxg[k][j] += solV[k][i] * phiVx[i * dim + j];
          }
        }
      }     

      adept::adouble solPg = 0;
      for(unsigned i = 0; i < nDofsP; i++) {
        solPg += phiP[i] * solP[i];
      }
      double nu = 1.;

      // *** phiV_i loop ***
      for(unsigned i = 0; i < nDofsV; i++) {
        std::vector < adept::adouble > NSV(dim, 0.);

        for(unsigned  k = 0; k < dim; k++) {  //momentum equation in k
          for(unsigned j = 0; j < dim; j++) {  // second index j in each equation
            NSV[k]   +=  nu * phiVx[i * dim + j] * (solVxg[k][j] + solVxg[j][k]);
          }
          NSV[k] += -solPg * phiVx[i * dim + k];
        }
        for(unsigned  k = 0; k < dim; k++) {
          mResV[k][i] += - NSV[k] * weight;
        }
      } // end phiV_i loop

      // *** phiP_i loop ***
      for(unsigned i = 0; i < nDofsP; i++) {
        for(int k = 0; k < dim; k++) {
          mResP[i] += - (-solVxg[k][k]) * phiP[i]  * weight;
        }
      } // end phiP_i loop
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube mRes in double Res and store them in RES
    Res.resize(nDofsVP);    //resize

    for(int i = 0; i < nDofsV; i++) {
      for(unsigned  k = 0; k < dim; k++) {
        Res[ k * nDofsV + i ] = -mResV[k][i].value();
      }
    }

    for(int i = 0; i < nDofsP; i++) {
      Res[ dim * nDofsV + i ] = -mResP[i].value();
    }

    RES->add_vector_blocked(Res, sysDof);

    //Extract and store the Jacobian

    Jac.resize(nDofsVP * nDofsVP);
    // define the dependent variables


    for(unsigned  k = 0; k < dim; k++) {
      s.dependent(&mResV[k][0], nDofsV);
    }
    s.dependent(&mResP[0], nDofsP);

    // define the independent variables
    for(unsigned  k = 0; k < dim; k++) {
      s.independent(&solV[k][0], nDofsV);
    }
    s.independent(&solP[0], nDofsP);


    // get the and store jacobian matrix (row-major)
    s.jacobian(&Jac[0], true);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();
  KK->close();

// KK->draw();
  // ***************** END ASSEMBLY *******************

}




