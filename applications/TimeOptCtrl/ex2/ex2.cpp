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
#include "LinearImplicitSystem.hpp"
#include "TransientSystem.hpp"
#include "adept.h"


#include "petsc.h"
#include "petscmat.h"
#include "PetscMatrix.hpp"

const double beta = 0.00001;
const double alpha = 0.00001;
const double t0 = 1.;
double Re = 100.;

const unsigned numberOfIterations = 3;
unsigned iext;

using namespace femus;

double flc4hs(double const &x, double const &eps);

double SetVariableTimeStep(const double time) {
  return 0.05;
}

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0;

  if(!strcmp(SolName, "bV") || !strcmp(SolName, "Vi")) {
    if(x[0] < 0. && x[1] < 0.5 && x[1] > -0.5 && x[2] < 0.5 && x[2] > -0.5) dirichlet = false;
  }
  if(!strcmp(SolName, "Vc")) {
    if(x[0] < 0. && x[1] < 0.5 && x[1] > -0.5 && x[2] < 0.5 && x[2] > -0.5) value = sin(time) * flc4hs(time - t0, t0);
  }
  else if(!strcmp(SolName, "bP") || !strcmp(SolName, "lP") || !strcmp(SolName, "Pc")) {
    dirichlet = false;
    value = 0.;
  }

  return dirichlet;
}

void AssembleSteadyStateControl(MultiLevelProblem& ml_prob);
void AssembleSystemZi(MultiLevelProblem& ml_prob);
void AssembleManifactureSolution(MultiLevelProblem& ml_prob);
void GetError(MultiLevelProblem& ml_prob);

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
  mlSol.AddSolution("Uc", LAGRANGE, SECOND, 2);
  mlSol.AddSolution("Vc", LAGRANGE, SECOND, 2);
  mlSol.AddSolution("Pc",  DISCONTINUOUS_POLYNOMIAL, FIRST);

  // add variables to mlSol
  mlSol.AddSolution("bU", LAGRANGE, SECOND);
  mlSol.AddSolution("bV", LAGRANGE, SECOND);
  mlSol.AddSolution("bP",  DISCONTINUOUS_POLYNOMIAL, FIRST);

  mlSol.AddSolution("lU", LAGRANGE, SECOND);
  mlSol.AddSolution("lV", LAGRANGE, SECOND);
  mlSol.AddSolution("lP",  DISCONTINUOUS_POLYNOMIAL, FIRST);

  mlSol.AddSolution("Ui", LAGRANGE, SECOND, 2);
  mlSol.AddSolution("Vi", LAGRANGE, SECOND, 2);
  mlSol.AddSolution("Pi",  DISCONTINUOUS_POLYNOMIAL, FIRST);

  char uName[10];
  char vName[10];
  char pName[10];
  const unsigned level = 0;
  Solution* sol = mlSol.GetSolutionLevel(level);

  for(unsigned i = 0; i < numberOfIterations; i++) {
    sprintf(uName, "U%d", i); //cascade solution
    mlSol.AddSolution(uName, LAGRANGE, SECOND, 2); //
    sprintf(vName, "V%d", i); //cascade solution
    mlSol.AddSolution(vName, LAGRANGE, SECOND, 2); //
    sprintf(pName, "P%d", i); //cascade solution
    mlSol.AddSolution(pName,  DISCONTINUOUS_POLYNOMIAL, FIRST);
  }
  mlSol.AddSolution("U0Older", LAGRANGE, SECOND, false);
  mlSol.AddSolution("V0Older", LAGRANGE, SECOND, false);



  mlSol.Initialize("All");

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.FixSolutionAtOnePoint("Pc");
  mlSol.FixSolutionAtOnePoint("bP");
  mlSol.FixSolutionAtOnePoint("lP");
  mlSol.FixSolutionAtOnePoint("Pi");

  mlSol.GenerateBdc("All");
  mlSol.GenerateBdc("Vc", "Time_dependent");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);

  // add system Poisson in mlProb as a Linear Implicit System
  TransientLinearImplicitSystem& systemC = mlProb.add_system < TransientLinearImplicitSystem > ("ManSol");

  systemC.AddSolutionToSystemPDE("Uc");
  systemC.AddSolutionToSystemPDE("Vc");
  systemC.AddSolutionToSystemPDE("Pc");

  // attach the assembling function to system
  systemC.SetAssembleFunction(AssembleManifactureSolution);
  systemC.AttachGetTimeIntervalFunction(SetVariableTimeStep);

  // initilaize and solve the system
  systemC.init();
  systemC.SetOuterSolver(PREONLY);

  TransientLinearImplicitSystem& system = mlProb.add_system < TransientLinearImplicitSystem > ("systembZ");

  // add solution "u" to system

  system.AddSolutionToSystemPDE("bU");
  system.AddSolutionToSystemPDE("bV");
  system.AddSolutionToSystemPDE("bP");


  system.AddSolutionToSystemPDE("lU");
  system.AddSolutionToSystemPDE("lV");
  system.AddSolutionToSystemPDE("lP");

  // attach the assembling function to system
  system.SetAssembleFunction(AssembleSteadyStateControl);
  system.AttachGetTimeIntervalFunction(SetVariableTimeStep);


  // initilaize and solve the system
  system.init();
  system.SetOuterSolver(PREONLY);


  TransientLinearImplicitSystem& systemi = mlProb.add_system < TransientLinearImplicitSystem > ("systemZi");

  systemi.AddSolutionToSystemPDE("Ui");
  systemi.AddSolutionToSystemPDE("Vi");
  systemi.AddSolutionToSystemPDE("Pi");

  // attach the assembling function to system
  systemi.SetAssembleFunction(AssembleSystemZi);
  systemi.AttachGetTimeIntervalFunction(SetVariableTimeStep);

  // initilaize and solve the system
  systemi.init();
  systemi.SetOuterSolver(PREONLY);

  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");
  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(false);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);

  *(sol->_SolOld[mlProb._ml_sol->GetIndex("U0")]) = *(sol->_Sol[mlProb._ml_sol->GetIndex("U0")]);
  *(sol->_SolOld[mlProb._ml_sol->GetIndex("V0")]) = *(sol->_Sol[mlProb._ml_sol->GetIndex("V0")]);

  for(unsigned t = 0; t < 100; t++) {

    *(sol->_Sol[mlProb._ml_sol->GetIndex("U0Older")]) = *(sol->_SolOld[mlProb._ml_sol->GetIndex("U0")]);
    *(sol->_Sol[mlProb._ml_sol->GetIndex("V0Older")]) = *(sol->_SolOld[mlProb._ml_sol->GetIndex("V0")]);

    mlSol.CopySolutionToOldSolution();
    systemC.MGsolve();

    double time =  system.GetTime();

    for(iext = 0; iext < numberOfIterations; iext++) {
      sprintf(uName, "U%d", iext);
      sprintf(vName, "V%d", iext);
      sprintf(pName, "P%d", iext);

      *(sol->_Sol[mlProb._ml_sol->GetIndex("Ui")]) = *(sol->_Sol[mlProb._ml_sol->GetIndex(uName)]);
      *(sol->_SolOld[mlProb._ml_sol->GetIndex("Ui")]) = *(sol->_SolOld[mlProb._ml_sol->GetIndex(uName)]);
      *(sol->_Sol[mlProb._ml_sol->GetIndex("Vi")]) = *(sol->_Sol[mlProb._ml_sol->GetIndex(vName)]);
      *(sol->_SolOld[mlProb._ml_sol->GetIndex("Vi")]) = *(sol->_SolOld[mlProb._ml_sol->GetIndex(vName)]);
      *(sol->_Sol[mlProb._ml_sol->GetIndex("Pi")]) = *(sol->_Sol[mlProb._ml_sol->GetIndex(pName)]);

      system.SetTime(time);
      systemi.SetTime(time);

      system.MGsolve();
      systemi.MGsolve();

      *(sol->_Sol[mlProb._ml_sol->GetIndex(uName)]) = *(sol->_Sol[mlProb._ml_sol->GetIndex("Ui")]);
      *(sol->_SolOld[mlProb._ml_sol->GetIndex(uName)]) = *(sol->_SolOld[mlProb._ml_sol->GetIndex("Ui")]);
      *(sol->_Sol[mlProb._ml_sol->GetIndex(vName)]) = *(sol->_Sol[mlProb._ml_sol->GetIndex("Vi")]);
      *(sol->_SolOld[mlProb._ml_sol->GetIndex(vName)]) = *(sol->_SolOld[mlProb._ml_sol->GetIndex("Vi")]);
      *(sol->_Sol[mlProb._ml_sol->GetIndex(pName)]) = *(sol->_Sol[mlProb._ml_sol->GetIndex("Pi")]);


      GetError(mlProb);

    }
    vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, t + 1);
  }

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
  TransientLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<TransientLinearImplicitSystem> ("systembZ");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh    = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el     = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual std::vector object in pdeSys (level)

  double dt =  mlPdeSys->GetIntervalTime();
  double time =  mlPdeSys->GetTime();

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  std::vector < unsigned > solVcIndex(dim);
  solVcIndex[0] = mlSol->GetIndex("Uc");
  solVcIndex[1] = mlSol->GetIndex("Vc");

  char Um1[10], Vm1[10];
  if(iext > 0) sprintf(Um1, "U%d", iext - 1);
  else sprintf(Um1, "U%d", 0);
  if(iext > 0) sprintf(Vm1, "V%d", iext - 1);
  else sprintf(Vm1, "V%d", 0);

  std::vector < unsigned > solVm1Index(dim);
  solVm1Index[0] = mlSol->GetIndex(Um1);
  solVm1Index[1] = mlSol->GetIndex(Vm1);



  std::vector < unsigned > solV0Index(dim);
  solV0Index[0] = mlSol->GetIndex("U0");
  solV0Index[1] = mlSol->GetIndex("V0");

  std::vector < unsigned > solV0OlderIndex(dim);
  solV0OlderIndex[0] = mlSol->GetIndex("U0Older");
  solV0OlderIndex[1] = mlSol->GetIndex("V0Older");

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

  std::vector < std::vector < double > >  solVm1(dim);
  std::vector < std::vector < double > >  solVm1Old(dim);

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
      solVm1[k].resize(nDofsV);
      solVm1Old[k].resize(nDofsV);
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

        if(iext != 0) {
          solVm1[k][i] = (*sol->_Sol[solVm1Index[k]])(solVDof);
          solVm1Old[k][i] = (*sol->_SolOld[solVm1Index[k]])(solVDof);
        }
        else {
          solVm1[k][i] = (*sol->_SolOld[solV0Index[k]])(solVDof);
          solVm1Old[k][i] = (*sol->_Sol[solV0OlderIndex[k]])(solVDof);
        }

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
      std::vector < double > solVm1g(dim, 0);
      std::vector < double > solVm1Oldg(dim, 0);

      std::vector < adept::adouble > solbVg(dim, 0);
      std::vector < adept::adouble > sollVg(dim, 0);
      std::vector < double > xg(dim, 0);

      std::vector < std::vector < adept::adouble > > solbVxg(dim, std::vector < adept::adouble >(dim, 0.));
      std::vector < std::vector < adept::adouble > > sollVxg(dim, std::vector < adept::adouble >(dim, 0.));

      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned  k = 0; k < dim; k++) {

          solVcg[k] += solVc[k][i] * phiV[i];
          solVm1g[k] += solVm1[k][i] * phiV[i];
          solVm1Oldg[k] += solVm1Old[k][i] * phiV[i];

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

      double iRe = 1. / Re;


      adept::adouble solbPg = 0;
      adept::adouble sollPg = 0;
      for(unsigned i = 0; i < nDofsP; i++) {
        solbPg += phiP[i] * solbP[i];
        sollPg += phiP[i] * sollP[i];
      }


      // *** phiV_i loop ***
      for(unsigned i = 0; i < nDofsV; i++) {
        std::vector < adept::adouble > NSVb(dim, 0.);
        std::vector < adept::adouble > NSVl(dim, 0.);

        for(unsigned  k = 0; k < dim; k++) {  //momentum equation in k
          for(unsigned j = 0; j < dim; j++) {  // second index j in each equation
            NSVb[k]   +=  iRe * phiVx[i * dim + j] * (solbVxg[k][j] + solbVxg[j][k]);
            NSVl[k]   +=  iRe * phiVx[i * dim + j] * (sollVxg[k][j] + sollVxg[j][k]) + beta * phiVx[i * dim + j] * (solbVxg[k][j] + solbVxg[j][k]);
          }
          NSVb[k] += (1 - 0.1 * (iext == 0)) * (solVm1g[k] - solVm1Oldg[k]) / dt * phiV[i] - solbPg * phiVx[i * dim + k];
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


void AssembleSystemZi(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  TransientLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<TransientLinearImplicitSystem> ("systemZi");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  double dt =  mlPdeSys->GetIntervalTime();
  double time =  mlPdeSys->GetTime();

  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual std::vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  std::vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("Ui");
  solVIndex[1] = mlSol->GetIndex("Vi");

  std::vector < unsigned > solbVIndex(dim);
  solbVIndex[0] = mlSol->GetIndex("bU");
  solbVIndex[1] = mlSol->GetIndex("bV");

  unsigned solVType = mlSol->GetSolutionType(solVIndex[0]);

  unsigned solPIndex;
  solPIndex = mlSol->GetIndex("Pi");
  unsigned solPType = mlSol->GetSolutionType(solPIndex);

  std::vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("Ui");
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("Vi");


  unsigned solPPdeIndex;
  solPPdeIndex = mlPdeSys->GetSolPdeIndex("Pi");

  std::vector < std::vector < double > >  solbV(dim);
  std::vector < std::vector < adept::adouble > >  solV(dim);
  std::vector < std::vector < double > >  solVOld(dim);
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
      solbV[k].resize(nDofsV);;
      solV[k].resize(nDofsV);
      solVOld[k].resize(nDofsV);
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
        solbV[k][i] = (*sol->_Sol[solbVIndex[k]])(solVDof);

        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);
        solVOld[k][i] = (*sol->_SolOld[solVIndex[k]])(solVDof);

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


    std::vector<bool> nodeIsControlBoundary(nDofsV, false);

    for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
      unsigned int facename = -(msh->el->GetFaceElementIndex(iel, jface) + 1);
      if(facename == 1) {
        unsigned nve = msh->GetElementFaceDofNumber(iel, jface, solVType);
        const unsigned felt = msh->GetElementFaceType(iel, jface);
        for(unsigned i = 0; i < nve; i++) {
          nodeIsControlBoundary[ msh->GetLocalFaceVertexIndex(iel, jface, i)] = true;
        }
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
      std::vector < double > solVOldg(dim, 0);


      std::vector < std::vector < adept::adouble > > solVxg(dim, std::vector < adept::adouble >(dim, 0.));

      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          solVg[k] += solV[k][i] * phiV[i];
          solVOldg[k] += solVOld[k][i] * phiV[i];
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
      double iRe = 1. / Re;

      // *** phiV_i loop ***
      for(unsigned i = 0; i < nDofsV; i++) {
        std::vector < adept::adouble > NSV(dim, 0.);

        for(unsigned  k = 0; k < dim; k++) {  //momentum equation in k
          for(unsigned j = 0; j < dim; j++) {  // second index j in each equation
            NSV[k]   +=  iRe * phiVx[i * dim + j] * (solVxg[k][j] + solVxg[j][k]);
          }
          NSV[k] += (solVg[k] - solVOldg[k]) / dt * phiV[i] - solPg * phiVx[i * dim + k];
        }
        for(unsigned  k = 0; k < dim; k++) {
          if(k != 1 || !nodeIsControlBoundary[i]) {
            mResV[k][i] += - NSV[k] * weight;
          }
          else {
            mResV[k][i] += solV[k][i] - solbV[k][i];
          }


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

void AssembleManifactureSolution(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  TransientLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<TransientLinearImplicitSystem> ("ManSol");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  double dt =  mlPdeSys->GetIntervalTime();
  double time =  mlPdeSys->GetTime();


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
  std::vector < std::vector < double > >  solVOld(dim);
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
      solV[k].resize(nDofsV);
      solVOld[k].resize(nDofsV);;
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
        solVOld[k][i] = (*sol->_SolOld[solVIndex[k]])(solVDof);

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
      std::vector < double > solVOldg(dim, 0);

      std::vector < std::vector < adept::adouble > > solVxg(dim, std::vector < adept::adouble >(dim, 0.));

      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          solVg[k] += solV[k][i] * phiV[i];
          solVOldg[k] += solVOld[k][i] * phiV[i];
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
      double iRe = 1. / Re;

      // *** phiV_i loop ***
      for(unsigned i = 0; i < nDofsV; i++) {
        std::vector < adept::adouble > NSV(dim, 0.);

        for(unsigned  k = 0; k < dim; k++) {  //momentum equation in k
          for(unsigned j = 0; j < dim; j++) {  // second index j in each equation
            NSV[k]   +=  iRe * phiVx[i * dim + j] * (solVxg[k][j] + solVxg[j][k]);
          }
          NSV[k] += (solVg[k] - solVOldg[k]) / dt * phiV[i] - solPg * phiVx[i * dim + k];
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



void GetError(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled


  char uName[10];
  char vName[10];
  sprintf(uName, "U%d", iext);
  sprintf(vName, "V%d", iext);

  //char sName[10];
  //sprintf(sName, "top%d", iext);

  //  extract pointers to the several objects that we are going to use
  //TransientLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<TransientLinearImplicitSystem> (sName);   // pointer to the linear implicit system named "Poisson"
  TransientLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<TransientLinearImplicitSystem> ("systemZi");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable

  double time =  mlPdeSys->GetTime();

  std::vector<unsigned> VcIndex(dim);
  VcIndex[0] = mlSol->GetIndex("Uc");
  VcIndex[1] = mlSol->GetIndex("Vc");

  std::vector<unsigned> VIndex(dim);
  VIndex[0] = mlSol->GetIndex(uName);
  VIndex[1] = mlSol->GetIndex(vName);

  unsigned solType = mlSol->GetSolutionType(VcIndex[0]);

  std::vector<std::vector < double >> solVc(dim);
  std::vector<std::vector < double >> solV(dim);

  std::vector < std::vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector <double> phi, gradPhi;  // local test function for velocity
  double weight;

  double  localIntegral[5] = {0, 0, 0, 0, 0};

  
  
  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

      
      
    unsigned group = msh->GetElementGroup(iel);
    if(group == group) {

      short unsigned ielGeom = msh->GetElementType(iel);
      unsigned nDofs = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
      for(unsigned  k = 0; k < dim; k++) {
        x[k].resize(nDofs);
      }
      for(unsigned k = 0; k < dim; k++) {
        solV[k].resize(nDofs);
        solVc[k].resize(nDofs);
      }
      
      // local storage of global mapping and solution
      for(unsigned i = 0; i < nDofs; i++) {
        unsigned solDof = msh->GetSolutionDof(i, iel, solType);
        for(unsigned k = 0; k < dim; k++) {
          solVc[k][i] = (*sol->_Sol[VcIndex[k]])(solDof);
          solV[k][i] = (*sol->_Sol[VIndex[k]])(solDof);
        }
      }

      // local storage of coordinates
      for(unsigned i = 0; i < nDofs; i++) {
        unsigned xDof  = msh->GetSolutionDof(i, iel, xType);
        for(unsigned k = 0; k < dim; k++) {
          x[k][i] = (*msh->_topology->_Sol[k])(xDof);
        }
      }



      // *** Gauss point loop ***
      for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {
        msh->_finiteElement[ielGeom][solType]->Jacobian(x, ig, weight, phi, gradPhi);

        std::vector<double> Vg(dim, 0.);
        std::vector<double> Vcg(dim, 0.);

        for(unsigned i = 0; i < nDofs; i++) {
          for(unsigned k = 0; k < dim; k++) {
            Vg[k] += solV[k][i] * phi[i];
            Vcg[k] += solVc[k][i] * phi[i];
          }
        }

        localIntegral[0] += Vcg[0] * weight;
        localIntegral[1] += Vcg[1] * weight;
        localIntegral[2] += Vg[0] * weight;
        localIntegral[3] += Vg[1] * weight;
        localIntegral[4] += ((Vg[0] - Vcg[0]) * (Vg[0] - Vcg[0]) + (Vg[1] - Vcg[1]) * (Vg[1] - Vcg[1])) * weight;

      } // end gauss point loop

    }
  } //end element loop for each process

  double globalIntegral[5] = {0., 0., 0., 0., 0.};

  MPI_Reduce(localIntegral, globalIntegral, 5, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


  if(iproc == 0) {

    std::string filename = "SolutionIntegral" + std::to_string(iext) + ".txt";
    std::ofstream fout;
    if(fabs(time) < 1.0e-10) {
      fout.open(filename.c_str(), std::ios::out);
    }
    else {
      fout.open(filename.c_str(), std::ios::app);
    }
    fout << time << " " << globalIntegral[0] << " " << globalIntegral[1] << " " << globalIntegral[2] << " " << globalIntegral[3] << " " << globalIntegral[4] << std::endl;
    fout.close();
  }


// ***************** END ASSEMBLY *******************

}



















double flc4hs(double const &x, double const &eps) {

  double r = x / eps;
  if(r < -1) {
    return 0.;
  }
  else if(r < 1.) {
    double r2 = r * r;
    double r3 = r * r2;
    double r5 = r3 * r2;
    double r7 = r5 * r2;
    double r9 = r7 * r2;
    return (128. + 315. * r - 420. * r3 + 378. * r5 - 180. * r7 + 35. * r9) / 256.;
  }
  else {
    return 1.;
  }
}


