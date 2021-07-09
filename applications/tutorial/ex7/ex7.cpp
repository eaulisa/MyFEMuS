/** \file Ex7.cpp
 *  \brief This example shows how to set and solve the weak form
 *   of the Boussinesq appoximation of the Navier-Stokes Equation
 *
 *  \f{eqnarray*}
 *  && \mathbf{V} \cdot \nabla T - \nabla \cdot\alpha \nabla T = 0 \\
 *  && \mathbf{V} \cdot \nabla \mathbf{V} - \nabla \cdot \nu (\nabla \mathbf{V} +(\nabla \mathbf{V})^T)
 *  +\nabla P = \beta T \mathbf{j} \\
 *  && \nabla \cdot \mathbf{V} = 0
 *  \f}
 *  in a unit box domain (in 2D and 3D) with given temperature 0 and 1 on
 *  the left and right walls, respectively, and insulated walls elsewhere.
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


using namespace femus;

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0.;

  if(!strcmp(SolName, "T")) {
    if(facename == 2) {
      value = 1.;
    }
    else if(facename == 3) {
      dirichlet = false; //Neumann
    }
  }
  else if(!strcmp(SolName, "P")) {
    dirichlet = false;
  }

  return dirichlet;
}


bool SetRefinementFlag(const std::vector < double >& x, const int& elemgroupnumber, const int& level) {

  bool refine = 0;

  if(elemgroupnumber == 6 && level < 4) refine = 1;

  if(elemgroupnumber == 7 && level < 5) refine = 1;

  if(elemgroupnumber == 8 && level < 6) refine = 1;

  return refine;

}



void AssembleBoussinesqAppoximation_AD(MultiLevelProblem& ml_prob);    //, unsigned level, const unsigned &levelMax, const bool &assembleMatrix );
void AssembleBoussinesqAppoximation(MultiLevelProblem& ml_prob);

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  //mlMsh.ReadCoarseMesh("./input/cube_hex.neu","seventh",scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/square_quad.neu", "seventh", scalingFactor);
  mlMsh.ReadCoarseMesh("./input/quadAMR.neu", "seventh", scalingFactor);
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
     probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();

//   unsigned numberOfUniformLevels = 3;
//   unsigned numberOfSelectiveLevels = 0;
//   mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  unsigned numberOfUniformLevels = 2;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , SetRefinementFlag);


  // erase all the coarse mesh levels
  //mlMsh.EraseCoarseLevels(numberOfUniformLevels - 3);

  // print mesh info
  mlMsh.PrintInfo();

  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("T", LAGRANGE, FIRST);

  mlSol.AddSolution("U", LAGRANGE, SECOND);
  mlSol.AddSolution("V", LAGRANGE, SECOND);

  if(dim == 3) mlSol.AddSolution("W", LAGRANGE, SECOND);

  //mlSol.AddSolution("P", LAGRANGE, FIRST);
  mlSol.AddSolution("P",  DISCONTINUOUS_POLYNOMIAL, FIRST);

  mlSol.AssociatePropertyToSolution("P", "Pressure"); //Maybe needed for adaptivity
  mlSol.Initialize("All");

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.FixSolutionAtOnePoint("P");
  mlSol.GenerateBdc("All");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);

  // add system Poisson in mlProb as a Linear Implicit System
  NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("Boussinesq");

  // add solution "u" to system
  system.AddSolutionToSystemPDE("T");
  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");

  if(dim == 3) system.AddSolutionToSystemPDE("W");

  system.AddSolutionToSystemPDE("P");

  system.SetLinearEquationSolverType(FEMuS_DEFAULT); // GMRES precontioned with geometric multigrid
  //system.SetLinearEquationSolverType(FEMuS_ASM); // Additive Swartz Method
  // attach the assembling function to system
  system.SetAssembleFunction(AssembleBoussinesqAppoximation_AD);

  system.SetMaxNumberOfNonLinearIterations(20);
  system.SetNonLinearConvergenceTolerance(1.e-12);

  system.SetMaxNumberOfLinearIterations(1);
  system.SetAbsoluteLinearConvergenceTolerance(1.e-12);

  system.SetMgType(F_CYCLE); //for steady state

  system.SetNumberPreSmoothingStep(1);
  system.SetNumberPostSmoothingStep(1);
  // initilaize and solve the system
  system.init();

  system.SetSolverFineGrids(RICHARDSON); //level Solver
  system.SetPreconditionerFineGrids(ILU_PRECOND); //level Preconditioner
  //system.SetTolerances(1.e-20, 1.e-20, 1.e+50, 40);
  system.SetTolerances(1.e-3, 1.e-20, 1.e+50, 5); //relative, absolute, breaking, number of iterations for the linear solver


// system.ClearVariablesToBeSolved(); //option for the FEMuS_ASM
  //   system.AddVariableToBeSolved("All");
  //   system.SetNumberOfSchurVariables(1);
  //   system.SetElementBlockNumber(4);
  //system.SetDirichletBCsHandling(ELIMINATION);
  //system.solve();
  system.MGsolve();

  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  GMVWriter gmvIO(&mlSol);
  variablesToBePrinted.push_back("all");
  gmvIO.SetDebugOutput(true);
  gmvIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);


  return 0;
}


void AssembleBoussinesqAppoximation_AD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<NonLinearImplicitSystem> ("Boussinesq");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*           msh         = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*           el          = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*   mlSol         = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*   sol         = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level];  // pointer to the equation (level) object
  SparseMatrix*   KK          = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*  RES         = pdeSys->_RES; // pointer to the global residual std::vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)


  //solution variable
  unsigned solTIndex;
  solTIndex = mlSol->GetIndex("T");    // get the position of "T" in the ml_sol object
  unsigned solTType = mlSol->GetSolutionType(solTIndex);    // get the finite element type for "T"

  std::vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("V");    // get the position of "V" in the ml_sol object

  if(dim == 3) solVIndex[2] = mlSol->GetIndex("W");       // get the position of "V" in the ml_sol object

  unsigned solVType = mlSol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

  unsigned solPIndex;
  solPIndex = mlSol->GetIndex("P");    // get the position of "P" in the ml_sol object
  unsigned solPType = mlSol->GetSolutionType(solPIndex);    // get the finite element type for "u"

  unsigned solTPdeIndex;
  solTPdeIndex = mlPdeSys->GetSolPdeIndex("T");    // get the position of "T" in the pdeSys object

  std::vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("U");    // get the position of "U" in the pdeSys object
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("V");    // get the position of "V" in the pdeSys object

  if(dim == 3) solVPdeIndex[2] = mlPdeSys->GetSolPdeIndex("W");

  unsigned solPPdeIndex;
  solPPdeIndex = mlPdeSys->GetSolPdeIndex("P");    // get the position of "P" in the pdeSys object

  std::vector < adept::adouble >  solT; // local solution
  std::vector < std::vector < adept::adouble > >  solV(dim);    // local solution
  std::vector < adept::adouble >  solP; // local solution

  std::vector< adept::adouble > mResT; // local redidual std::vector
  std::vector< std::vector < adept::adouble > > mResV(dim);    // local redidual std::vector
  std::vector< adept::adouble > mResP; // local redidual std::vector

  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector <double> phiV;  // local test function
  std::vector <double> phiV_x; // local test function first order partial derivatives


  std::vector <double> phiT;  // local test function
  std::vector <double> phiT_x; // local test function first order partial derivatives

  double* phiP;
  double weight; // gauss point weight

  std::vector< int > sysDof; // local to global pdeSys dofs
  std::vector< double > Res; // local redidual std::vector
  std::vector < double > Jac;


  KK->zero(); // Set to zero all the entries of the Global Matrix
  RES->zero();

  // element loop: each process loops only on the elements that owns
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);

    unsigned nDofsT = msh->GetElementDofNumber(iel, solTType);   // number of solution element dofs
    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);   // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);   // number of solution element dofs
    unsigned nDofsX = (nDofsT > nDofsV) ? nDofsT : nDofsV;       // number of coordinate element dofs

    unsigned nDofsTVP = nDofsT + dim * nDofsV + nDofsP;
    // resize local arrays
    sysDof.resize(nDofsTVP);

    solT.resize(nDofsT);

    for(unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofsV);
      coordX[k].resize(nDofsX);
    }

    solP.resize(nDofsP);

    mResT.assign(nDofsV, 0.);

    for(unsigned  k = 0; k < dim; k++) {
      mResV[k].assign(nDofsV, 0.);
    }

    mResP.assign(nDofsP, 0.);

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsT; i++) {
      unsigned solTDof = msh->GetSolutionDof(i, iel, solTType);    // global to global mapping between solution node and solution dof
      solT[i] = (*sol->_Sol[solTIndex])(solTDof);      // global extraction and local storage for the solution
      sysDof[i] = pdeSys->GetSystemDof(solTIndex, solTPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dofs
    }

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // global to global mapping between solution node and solution dof

      for(unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);      // global extraction and local storage for the solution
        sysDof[nDofsT + (k * nDofsV) + i ] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }

    for(unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);    // global to global mapping between solution node and solution dof
      solP[i] = (*sol->_Sol[solPIndex])(solPDof);      // global extraction and local storage for the solution
      sysDof[nDofsT + (dim * nDofsV) + i ] = pdeSys->GetSystemDof(solPIndex, solPPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // global to global mapping between coordinates node and coordinate dof

      for(unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solVType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solTType]->Jacobian(coordX, ig, weight, phiT, phiT_x);
      msh->_finiteElement[ielGeom][solVType]->Jacobian(coordX, ig, weight, phiV, phiV_x);

      phiP = msh->_finiteElement[ielGeom][solPType]->GetPhi(ig); // very fast, the gradient is expensive

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      adept::adouble solT_gss = 0;
      std::vector < adept::adouble > gradSolT_gss(dim, 0.);

      for(unsigned i = 0; i < nDofsT; i++) {
        solT_gss += phiT[i] * solT[i];

        for(unsigned j = 0; j < dim; j++) {
          gradSolT_gss[j] += phiT_x[i * dim + j] * solT[i];
        }
      }

      std::vector < adept::adouble > solV_gss(dim, 0);
      std::vector < std::vector < adept::adouble > > gradSolV_gss(dim);

      for(unsigned  k = 0; k < dim; k++) {
        gradSolV_gss[k].assign(dim,0.);
      }

      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          solV_gss[k] += phiV[i] * solV[k][i];
        }

        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            gradSolV_gss[k][j] += phiV_x[i * dim + j] * solV[k][i];
          }
        }
      }

      adept::adouble solP_gss = 0;

      for(unsigned i = 0; i < nDofsP; i++) {
        solP_gss += phiP[i] * solP[i];
      }

      double nu = 1.;
      double alpha = .5;
      double beta = 2000.;

      // *** phiT_i loop ***
      for(unsigned i = 0; i < nDofsT; i++) {
        adept::adouble Temp = 0.;

        for(unsigned k = 0; k < dim; k++) {
          Temp +=  alpha * phiT_x[i * dim + k] * gradSolT_gss[k];
          Temp +=  phiT[i] * (solV_gss[k] * gradSolT_gss[k]);
        }

        mResT[i] -= - Temp * weight;
      } // end phiT_i loop


      // *** phiV_i loop ***
      for(unsigned i = 0; i < nDofsV; i++) {
        std::vector < adept::adouble > NSV(dim, 0.);

        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            NSV[k]   +=  nu * phiV_x[i * dim + j] * (gradSolV_gss[k][j] + gradSolV_gss[j][k]);
            NSV[k]   +=  phiV[i] * (solV_gss[j] * gradSolV_gss[k][j]);
          }
        }

        for(unsigned  k = 0; k < dim; k++) {
          NSV[k] += -solP_gss * phiV_x[i * dim + k];
        }

        NSV[1] += -beta * solT_gss * phiV[i];

        for(unsigned  k = 0; k < dim; k++) {
          mResV[k][i] -= - NSV[k] * weight;
        }
      } // end phiV_i loop

      // *** phiP_i loop ***
      for(unsigned i = 0; i < nDofsP; i++) {
        for(int k = 0; k < dim; k++) {
          mResP[i] -= gradSolV_gss[k][k] * phiP[i]  * weight;
        }
      } // end phiP_i loop

    } // end gauss point loop

    // } // endif single element not refined or fine grid loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube mRes in double Res and store them in RES
    Res.resize(nDofsTVP);    //resize

    for(int i = 0; i < nDofsT; i++) {
      Res[i] = -mResT[i].value();
    }

    for(int i = 0; i < nDofsV; i++) {
      for(unsigned  k = 0; k < dim; k++) {
        Res[ i + nDofsT + k * nDofsV ] = -mResV[k][i].value();
      }
    }

    for(int i = 0; i < nDofsP; i++) {
      Res[ i + nDofsT + dim * nDofsV ] = -mResP[i].value();
    }

    RES->add_vector_blocked(Res, sysDof);

    //Extarct and store the Jacobian
    Jac.resize(nDofsTVP * nDofsTVP);
    // define the dependent variables
    s.dependent(&mResT[0], nDofsT);

    for(unsigned  k = 0; k < dim; k++) {
      s.dependent(&mResV[k][0], nDofsV);
    }

    s.dependent(&mResP[0], nDofsP);

    // define the independent variables
    s.independent(&solT[0], nDofsT);

    for(unsigned  k = 0; k < dim; k++) {
      s.independent(&solV[k][0], nDofsV);
    }

    s.independent(&solP[0], nDofsP);

    // get the and store jacobian matrix (row-major)
    s.jacobian(&Jac[0] , true);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();

  KK->close();

  // ***************** END ASSEMBLY *******************
}



//Attempting to create J by hand
void AssembleBoussinesqAppoximation(MultiLevelProblem& ml_prob) {
//  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled


  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<NonLinearImplicitSystem> ("Boussinesq");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*           msh         = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*           el          = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*   mlSol         = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*   sol         = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level];  // pointer to the equation (level) object
  SparseMatrix*   KK          = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*  RES         = pdeSys->_RES; // pointer to the global residual std::vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)


  //solution variable
  unsigned solTIndex;
  solTIndex = mlSol->GetIndex("T");    // get the position of "T" in the ml_sol object
  unsigned solTType = mlSol->GetSolutionType(solTIndex);    // get the finite element type for "T"

  std::vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("V");    // get the position of "V" in the ml_sol object

  if(dim == 3) solVIndex[2] = mlSol->GetIndex("W");       // get the position of "V" in the ml_sol object

  unsigned solVType = mlSol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

  unsigned solPIndex;
  solPIndex = mlSol->GetIndex("P");    // get the position of "P" in the ml_sol object
  unsigned solPType = mlSol->GetSolutionType(solPIndex);    // get the finite element type for "u"

  unsigned solTPdeIndex;
  solTPdeIndex = mlPdeSys->GetSolPdeIndex("T");    // get the position of "T" in the pdeSys object

  std::vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("U");    // get the position of "U" in the pdeSys object
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("V");    // get the position of "V" in the pdeSys object

  if(dim == 3) solVPdeIndex[2] = mlPdeSys->GetSolPdeIndex("W");

  unsigned solPPdeIndex;
  solPPdeIndex = mlPdeSys->GetSolPdeIndex("P");    // get the position of "P" in the pdeSys object

  std::vector < double >  solT; // local solution
  std::vector < std::vector < double > >  solV(dim);    // local solution
  std::vector < double >  solP; // local solution

  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector <double> phiV;  // local test function
  std::vector <double> phiV_x; // local test function first order partial derivatives


  std::vector <double> phiT;  // local test function
  std::vector <double> phiT_x; // local test function first order partial derivatives

  double* phiP;
  double weight; // gauss point weight

  std::vector< int > sysDof; // local to global pdeSys dofs
  std::vector< double > Res; // local redidual std::vector
  std::vector < double > Jac;


  KK->zero(); // Set to zero all the entries of the Global Matrix
  RES->zero();

  // element loop: each process loops only on the elements that owns


  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);

    unsigned nDofsT = msh->GetElementDofNumber(iel, solTType);   // number of solution element dofs
    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);   // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);   // number of solution element dofs
    unsigned nDofsX = (nDofsT > nDofsV) ? nDofsT : nDofsV;       // number of coordinate element dofs

    unsigned nDofsTVP = nDofsT + dim * nDofsV + nDofsP;
    Jac.assign(nDofsTVP * nDofsTVP, 0.);

    // resize local arrays
    sysDof.resize(nDofsTVP);
    Res.assign(nDofsTVP, 0.);
    Jac.assign(nDofsTVP * nDofsTVP, 0.);


    solT.resize(nDofsT);

    for(unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofsV);
      coordX[k].resize(nDofsX);
    }

    solP.resize(nDofsP);

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsT; i++) {
      unsigned solTDof = msh->GetSolutionDof(i, iel, solTType);    // global to global mapping between solution node and solution dof
      solT[i] = (*sol->_Sol[solTIndex])(solTDof);      // global extraction and local storage for the solution
      sysDof[i] = pdeSys->GetSystemDof(solTIndex, solTPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dofs
    }

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // global to global mapping between solution node and solution dof

      for(unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);      // global extraction and local storage for the solution
        sysDof[nDofsT + (k * nDofsV) + i ] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }

    for(unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);    // global to global mapping between solution node and solution dof
      solP[i] = (*sol->_Sol[solPIndex])(solPDof);      // global extraction and local storage for the solution
      sysDof[nDofsT + (dim * nDofsV) + i ] = pdeSys->GetSystemDof(solPIndex, solPPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // global to global mapping between coordinates node and coordinate dof

      for(unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solVType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solTType]->Jacobian(coordX, ig, weight, phiT, phiT_x);
      msh->_finiteElement[ielGeom][solVType]->Jacobian(coordX, ig, weight, phiV, phiV_x);

      phiP = msh->_finiteElement[ielGeom][solPType]->GetPhi(ig); // very fast, the gradient is expensive

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      double solT_gss = 0.;
      std::vector < double > gradSolT_gss(dim, 0.);

      for(unsigned i = 0; i < nDofsT; i++) {
        solT_gss += phiT[i] * solT[i];

        for(unsigned j = 0; j < dim; j++) {
          gradSolT_gss[j] += phiT_x[i * dim + j] * solT[i];
        }
      }

      std::vector < double > solV_gss(dim, 0.);
      std::vector < std::vector < double > > gradSolV_gss(dim);

      for(unsigned  k = 0; k < dim; k++) {
        gradSolV_gss[k].assign(dim,0.);
      }

      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          solV_gss[k] += phiV[i] * solV[k][i];
        }

        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            gradSolV_gss[k][j] += phiV_x[i * dim + j] * solV[k][i];
          }
        }
      }

      double solP_gss = 0.;

      for(unsigned i = 0; i < nDofsP; i++) {
        solP_gss += phiP[i] * solP[i];
      }

      double nu = 1.;
      double alpha = .5;
      double beta = 2000.;

      // *** phiT_i loop ***
      for(unsigned i = 0; i < nDofsT; i++) {
        double Temp = 0.;

        for(unsigned k = 0; k < dim; k++) {
          Temp +=  alpha * phiT_x[i * dim + k] * gradSolT_gss[k];
          Temp +=  phiT[i] * (solV_gss[k] * gradSolT_gss[k]);
        }

        Res[i] -= Temp * weight;
      } // end phiT_i loop


      // *** phiV_i loop ***
      for(unsigned i = 0; i < nDofsV; i++) {
        std::vector < double > NSV(dim, 0.);

        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            NSV[k]   +=  nu * phiV_x[i * dim + j] * (gradSolV_gss[k][j] + gradSolV_gss[j][k]);
            NSV[k]   +=  phiV[i] * (solV_gss[j] * gradSolV_gss[k][j]);
          }
        }

        for(unsigned  k = 0; k < dim; k++) {
          NSV[k] += -solP_gss * phiV_x[i * dim + k];
        }

        //Temperature only for y - axis
        NSV[1] += -beta * solT_gss * phiV[i];

        for(unsigned  k = 0; k < dim; k++) {
          Res[nDofsT + k * nDofsV + i] -= NSV[k] * weight;
        }
      } // end phiV_i loop

      // *** phiP_i loop ***
      for(unsigned i = 0; i < nDofsP; i++) {
        for(int k = 0; k < dim; k++) {
          Res[nDofsT + dim * nDofsV + i] -= -gradSolV_gss[k][k] * phiP[i]  * weight;
        }
      } // end phiP_i loop


      for(unsigned i = 0; i < nDofsT; i++) { //Jacobian for the temperature equation
        unsigned Trow = i * nDofsTVP;
        for(unsigned j = 0; j < nDofsT; j++) {
          unsigned Tcolumn = j;
          for(unsigned k = 0; k < dim; k++) {
            Jac[ Trow + Tcolumn] += alpha * phiT_x[i * dim + k] * phiT_x[j * dim + k] * weight; //Laplace of temperature
            Jac[ Trow + Tcolumn] += phiT[i] * solV_gss[k] * phiT_x[j * dim + k] * weight;  //nonlinear temperature
          }
        }
        for(unsigned j = 0; j < nDofsV; j++) {
          for(unsigned J = 0; J < dim ; J++) { //column velocity blocks or dimension
            unsigned VJcolumn = nDofsT + J * nDofsV + j;
            Jac[ Trow + VJcolumn] += phiT[i] * phiV[j] * gradSolT_gss[J] * weight;  //nonlinear temperature
          }
        }
      }
      


      for(unsigned i = 0; i < nDofsV; i++) { //Jacobian for the Navier-Stokes equation
        for(unsigned I = 0; I < dim; I++) { //row velocity blocks or dimension
          unsigned VIrow = nDofsT + I * nDofsV + i;

          if(I == 1) { // buoyancy 
            for(unsigned j = 0; j < nDofsT; j++) {
              unsigned Tcolumn = j;
              Jac[ VIrow * nDofsTVP + Tcolumn] += -beta * phiV[i] * phiT[j] * weight; //beta * tempature
            }
          }

          for(unsigned j = 0; j < nDofsV; j++) {
            unsigned VIcolumn = nDofsT + I * nDofsV + j;
            for(unsigned J = 0; J < dim ; J++) { //column velocity blocks or dimension
              unsigned VJcolumn = nDofsT + J * nDofsV + j;
              Jac[ VIrow * nDofsTVP + VIcolumn] += nu * phiV_x[i * dim + J] * phiV_x[j * dim + J] * weight; //diagonal diffusion
              Jac[ VIrow * nDofsTVP + VJcolumn] += nu * phiV_x[i * dim + J] * phiV_x[j * dim + I] * weight; //off-diagonal diffusion

              Jac[ VIrow * nDofsTVP + VIcolumn] += phiV[i] * solV_gss[J] * phiV_x[j * dim + J] * weight; //diagonal nonlinear
              Jac[ VIrow * nDofsTVP + VJcolumn] += phiV[i] * phiV[j] * gradSolV_gss[I][J] * weight; //off-diagonal nonlinear

            }
          }

          for(unsigned j = 0; j < nDofsP; j++) {
            unsigned Pcolumn = nDofsT + dim * nDofsV + j;
            Jac[ VIrow * nDofsTVP + Pcolumn] += - phiV_x[i * dim + I] * phiP[j] * weight; //pressure gradient
            Jac[ Pcolumn * nDofsTVP + VIrow] += - phiV_x[i * dim + I] * phiP[j] * weight; //continuity
          }

        }
      }

      
    }
    // end gauss point loop

    RES->add_vector_blocked(Res, sysDof);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);

  } //end element loop for each process

  RES->close();
  KK->close();

 // VecView((static_cast<PetscVector*>(RES))->vec(),  PETSC_VIEWER_STDOUT_SELF);
 // MatView((static_cast<PetscMatrix*>(KK))->mat(), PETSC_VIEWER_STDOUT_SELF);

//   PetscViewer    viewer;
//   PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL, NULL, 0, 0, 900, 900, &viewer);
//   PetscObjectSetName((PetscObject) viewer, "Boussinesq matrix");
//   PetscViewerPushFormat(viewer, PETSC_VIEWER_DRAW_LG);
//   MatView((static_cast<PetscMatrix*>(KK))->mat(), viewer);
//   double a;
//   std::cin >> a;

// ***************** END ASSEMBLY *******************
}










