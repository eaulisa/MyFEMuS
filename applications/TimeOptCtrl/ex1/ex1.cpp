/** \file Ex13.cpp
 *  \brief This example shows how to set and solve the weak form
 *   of the time dependent Stress-Strain Equation, for a beam with pressure
 *
 *  \nabla \cdot \sigma + \mass*acceleration = F
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
#include "TransientSystem.hpp"
#include "adept.h"
#include <fstream>



using namespace femus;

const unsigned alpha = 1.0e-5;
const double beta = 1.0e-5;
const double Pe = 5.;
const double t0 = 3.;
const std::vector <double> v = {1., 0.2};

bool steady = false;
bool manifacturedSolution = true;
unsigned iext;
const unsigned numberOfIterations = 6;

void AssembleTimeOptimalControl(MultiLevelProblem& ml_prob);
void AssembleManifacturedSolution(MultiLevelProblem& ml_prob);
void GetError(MultiLevelProblem& ml_prob);
void StoreOldSolutions(MultiLevelProblem& ml_prob);
double flc4hs(double const &x, double const &eps);

double SetVariableTimeStep(const double time) {
  return 0.05;
}
double GetTargetSolution(const double &t) {
  return 0.9 * sin(t) * flc4hs(t - t0, t0);
}

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true;
  value = 0.;

  std::string name = SolName;
  std::string subName = name.substr(0, 1);

  double Td1 = 1. + sin(2.* time) * flc4hs(time - t0, t0);
  double Td2 = 2. + sin(.5 * time) * flc4hs(time - t0, t0);
  double Tc = Td1 + (x[0] - 0.25) * (0.75 - x[0]) * 3. * sin(2.* time) * flc4hs(time - t0, t0);

  if(facename == 4 ||  facename == 2) {
    if(subName == "T" || subName == "z") {
      value =  Td1;
    }
  }
  else if(facename == 3) {
    if(subName == "T" || subName == "z") {
      value = Td2;
    }
  }
  else if(facename == 1) {
    if(name == "Tc") {
      value = Tc;
    }
    else if(subName == "T" || subName == "z") {
      dirichlet = false;
    }
  }
  return dirichlet;
}



int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
//   MultiLevelMesh mlMsh;
//
//   unsigned nx = 10;
//   unsigned ny = 20;
//   unsigned nz = 1;
//
//   double length =  2;
//   double lengthx = 1;
//
//   mlMsh.GenerateCoarseBoxMesh(nx, ny, 0, 0., lengthx, 0., length, 0., 0., QUAD9, "seventh");

  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;

  MultiLevelMesh mlMsh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels, "./input/rectangle.neu", "fifth", 1., NULL);

  unsigned dim = mlMsh.GetDimension();

  //mlMsh.RefineMesh(numberOfUniformLevels, numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // erase all the coarse mesh levels
  mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  MultiLevelSolution mlSol(&mlMsh);

  mlSol.AddSolution("Tc", LAGRANGE, SECOND, 2); // manifactured solution
  // add variables to mlSol
  for(unsigned i = 0; i < numberOfIterations; i++) {
    char zName[10];
    sprintf(zName, "z%d", i); //cascade solution
    mlSol.AddSolution(zName, LAGRANGE, SECOND, 2); //
  }
  //support variables
  mlSol.AddSolution("Ti", LAGRANGE, SECOND);
  mlSol.AddSolution("li", LAGRANGE, SECOND);
  mlSol.AddSolution("zi", LAGRANGE, SECOND, 2);
  mlSol.AddSolution("z0Older", LAGRANGE, SECOND, false);

  mlSol.Initialize("All");

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

  mlSol.GenerateBdc("Tc", "Time_dependent");
  for(unsigned i = 0; i < numberOfIterations; i++) {
    char zName[10];
    sprintf(zName, "z%d", i);
    mlSol.GenerateBdc(zName, "Time_dependent");
  }
  mlSol.GenerateBdc("zi", "Time_dependent");
  mlSol.GenerateBdc("Ti", "Time_dependent");
  mlSol.GenerateBdc("li", "Steady");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);

  TransientLinearImplicitSystem *systemTc = &(mlProb.add_system<TransientLinearImplicitSystem> ("sysTc"));
  //add solution "D" to system
  systemTc->AddSolutionToSystemPDE("Tc");
  systemTc->AttachGetTimeIntervalFunction(SetVariableTimeStep);
  systemTc->SetAssembleFunction(AssembleManifacturedSolution);
  systemTc->init();
  systemTc->SetOuterSolver(PREONLY);


  TransientLinearImplicitSystem *system = &(mlProb.add_system<TransientLinearImplicitSystem> ("TmOptCtr"));

  system->AddSolutionToSystemPDE("Ti");
  system->AddSolutionToSystemPDE("li");
  system->AddSolutionToSystemPDE("zi");

  // attach the assembling function to system
  system->AttachGetTimeIntervalFunction(SetVariableTimeStep);
  system->SetAssembleFunction(AssembleTimeOptimalControl);

  // initilaize and solve the system
  system->init();
  system->SetOuterSolver(PREONLY);


//   std::vector<TransientLinearImplicitSystem*> system(numberOfIterations);
//   for(unsigned i = 0; i < numberOfIterations; i++) {
//     char sName[10];
//     sprintf(sName, "top%d", i);
//     system[i] = &(mlProb.add_system < TransientLinearImplicitSystem > (sName));
//     //add solution "D" to system
//     system[i]->AddSolutionToSystemPDE("Ti");
//     system[i]->AddSolutionToSystemPDE("li");
//     system[i]->AddSolutionToSystemPDE("zi");
//
//     // attach the assembling function to system
//     system[i]->AttachGetTimeIntervalFunction(SetVariableTimeStep);
//     system[i]->SetAssembleFunction(AssembleTimeOptimalControl);
//
//     // initilaize and solve the system
//     system[i]->init();
//     system[i]->SetOuterSolver(PREONLY);
//   }
  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");


  //Get Initial Data
  steady = true;
  {
    double dt =  systemTc->GetIntervalTime();
    systemTc->SetTime(-dt);
    systemTc->MGsolve(); //solve for A, using DOld, VOld, and AOld
  }
  for(iext = 0; iext < numberOfIterations; iext++) {

    char zName[10];
    sprintf(zName, "z%d", iext);
    const unsigned level = systemTc->GetLevelToAssemble();
    Solution* sol = mlProb._ml_sol->GetSolutionLevel(level);

    //double dt =  system[iext]->GetIntervalTime();
    //system[iext]->SetTime(-dt);

    double time = system->GetTime();

    std::cout<<"AAAAAAAAAA "<<time;

    double dt =  system->GetIntervalTime();
    system->SetTime(-dt);

    *(sol->_Sol[mlProb._ml_sol->GetIndex("zi")]) = *(sol->_Sol[mlProb._ml_sol->GetIndex(zName)]);
    //system[iext]->MGsolve(); //solve for A, using DOld, VOld, and AOld
    system->MGsolve(); //solve for A, using DOld, VOld, and AOld
    *(sol->_Sol[mlProb._ml_sol->GetIndex(zName)]) = *(sol->_Sol[mlProb._ml_sol->GetIndex("zi")]);

  }

  {
    //Old and Older solution: inizialization
    //const unsigned level = system[0]->GetLevelToAssemble();
    const unsigned level = system->GetLevelToAssemble();
    Solution* sol = mlProb._ml_sol->GetSolutionLevel(level);
    *(sol->_SolOld[mlProb._ml_sol->GetIndex("z0")]) = *(sol->_Sol[mlProb._ml_sol->GetIndex("z0")]);
    *(sol->_Sol[mlProb._ml_sol->GetIndex("z0Older")]) = *(sol->_SolOld[mlProb._ml_sol->GetIndex("z0")]);
  }

  steady = false;

  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(false);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);

  for(iext = 0; iext < numberOfIterations; iext++) {
    GetError(mlProb);
  }

  for(unsigned t = 0; t < 200; t++) {
    {
      //store Older solution
      //const unsigned level = system[0]->GetLevelToAssemble();
      const unsigned level = system->GetLevelToAssemble();
      Solution*  sol = mlProb._ml_sol->GetSolutionLevel(level);
      *(sol->_Sol[mlProb._ml_sol->GetIndex("z0Older")]) = *(sol->_SolOld[mlProb._ml_sol->GetIndex("z0")]);
    }

    mlSol.CopySolutionToOldSolution();
    systemTc->MGsolve();

    //double time =  system[0]->GetTime();
    double time =  system->GetTime();

    for(iext = 0; iext < numberOfIterations; iext++) {

      char zName[10];
      sprintf(zName, "z%d", iext);
      //const unsigned level = system[iext]->GetLevelToAssemble();
      const unsigned level = system->GetLevelToAssemble();
      Solution* sol = mlProb._ml_sol->GetSolutionLevel(level);

      *(sol->_Sol[mlProb._ml_sol->GetIndex("zi")]) = *(sol->_Sol[mlProb._ml_sol->GetIndex(zName)]);
      *(sol->_SolOld[mlProb._ml_sol->GetIndex("zi")]) = *(sol->_SolOld[mlProb._ml_sol->GetIndex(zName)]);

      //system[iext]->SetTime(time);
      //system[iext]->MGsolve(); //solve for A, using DOld, VOld, and AOld

      system->SetTime(time);
      system->MGsolve(); //solve for A, using DOld, VOld, and AOld

      *(sol->_Sol[mlProb._ml_sol->GetIndex(zName)]) = *(sol->_Sol[mlProb._ml_sol->GetIndex("zi")]);
      *(sol->_SolOld[mlProb._ml_sol->GetIndex(zName)]) = *(sol->_SolOld[mlProb._ml_sol->GetIndex("zi")]);

      GetError(mlProb);
    }
    vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, t + 1);
  }
  return 0;
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




void AssembleTimeOptimalControl(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled


  //char sName[10];
  char zNameM1[10];
  //sprintf(sName, "top%d", iext);

  //sprintf(zName, "z%d", iext);
  if(iext > 0) sprintf(zNameM1, "z%d", iext - 1);
  else sprintf(zNameM1, "z%d", 0);

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  //TransientLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<TransientLinearImplicitSystem> (sName);   // pointer to the linear implicit system named "Poisson"

  TransientLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<TransientLinearImplicitSystem> ("TmOptCtr");   // pointer to the linear implicit

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

  double dt =  mlPdeSys->GetIntervalTime();
  double time =  mlPdeSys->GetTime();

  unsigned TcIndex = mlSol->GetIndex("Tc");

  unsigned TiIndex = mlSol->GetIndex("Ti");
  unsigned liIndex = mlSol->GetIndex("li");
  unsigned ziIndex = mlSol->GetIndex("zi");
  unsigned ziOlderIndex = mlSol->GetIndex("z0Older");

  unsigned ziM1Index = mlSol->GetIndex(zNameM1);

  unsigned solType = mlSol->GetSolutionType(TiIndex);

  unsigned TiPdeIndex = mlPdeSys->GetSolPdeIndex("Ti");
  unsigned liPdeIndex = mlPdeSys->GetSolPdeIndex("li");

  //unsigned ziiIndex = mlSol->GetIndex("zi");
  unsigned ziPdeIndex = mlPdeSys->GetSolPdeIndex("zi");

  std::vector < adept::adouble > solTi;
  std::vector < adept::adouble > solli;
  std::vector < adept::adouble > solzi;
  std::vector < double > solziOld;

  std::vector < double > solTc;
  std::vector < double > solziOlder;

  std::vector < double > solziM1;
  std::vector < double > solziM1Old;

  std::vector < adept::adouble > mResTi;
  std::vector < adept::adouble > mResli;
  std::vector < adept::adouble > mReszi;

  std::vector < std::vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector <double> phi;  // local test function for velocity
  std::vector <double> gradPhi; // local test function first order partial derivatives
  double weight;

  std::vector< unsigned > sysDof; // local to global pdeSys dofs
  std::vector< double > Res; // local redidual std::vector
  std::vector < double > Jac;

  RES->zero(); // Set to zero all the entries of the Global Residual std::vector
  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);

    unsigned group = msh->GetElementGroup(iel);

    unsigned nDofs = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
    unsigned nDofsAll = 3 * nDofs;
    sysDof.resize(nDofsAll);

    for(unsigned  k = 0; k < dim; k++) {
      x[k].resize(nDofs);
    }
    solTc.resize(nDofs);



    solTi.resize(nDofs);
    solli.resize(nDofs);
    solzi.resize(nDofs);
    solziOld.resize(nDofs);
    solziOlder.resize(nDofs);

    solziM1.resize(nDofs);
    solziM1Old.resize(nDofs);

    mResTi.assign(nDofs, 0.);
    mResli.assign(nDofs, 0.);
    mReszi.assign(nDofs, 0.);

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, solType);

      solTc[i] = (*sol->_Sol[TcIndex])(solDof);

      solTi[i] = (*sol->_Sol[TiIndex])(solDof);
      solli[i] = (*sol->_Sol[liIndex])(solDof);

      solzi[i] = (*sol->_Sol[ziIndex])(solDof);
      solziOld[i] = (*sol->_SolOld[ziIndex])(solDof);
      solziOlder[i] = (*sol->_Sol[ziOlderIndex])(solDof);

      solziM1[i] = (*sol->_Sol[ziM1Index])(solDof);
      solziM1Old[i] = (*sol->_SolOld[ziM1Index])(solDof);

      sysDof[i] = pdeSys->GetSystemDof(TiIndex, TiPdeIndex, i, iel);
      sysDof[nDofs + i] = pdeSys->GetSystemDof(liIndex, liPdeIndex, i, iel);
      sysDof[2 * nDofs + i] = pdeSys->GetSystemDof(ziIndex, ziPdeIndex, i, iel);
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);
      for(unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k])(xDof);
      }
    }

    std::vector<bool> nodeIsControlBoundary(nDofs, false);

    for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
      unsigned int facename = -(msh->el->GetFaceElementIndex(iel, jface) + 1);
      if(facename == 1) {
        unsigned nve = msh->GetElementFaceDofNumber(iel, jface, solType);
        const unsigned felt = msh->GetElementFaceType(iel, jface);
        for(unsigned i = 0; i < nve; i++) {
          nodeIsControlBoundary[ msh->GetLocalFaceVertexIndex(iel, jface, i)] = true;
        }
      }
    }


    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {
      msh->_finiteElement[ielGeom][solType]->Jacobian(x, ig, weight, phi, gradPhi);


      adept::adouble Tig = 0;
      adept::adouble lig = 0;
      adept::adouble zig = 0;
      double ziOldg = 0;
      double ziOlderg = 0;
      double ziM1g = 0;
      double ziM1Oldg = 0;
      double Tc = 0;

      std::vector < adept::adouble > gradTig(dim, 0);
      std::vector < adept::adouble > gradlig(dim, 0);
      std::vector < adept::adouble > gradzig(dim, 0);
      std::vector < double > gradziOldg(dim, 0);

      for(unsigned i = 0; i < nDofs; i++) {

        Tig += solTi[i] * phi[i];
        lig += solli[i] * phi[i];
        zig += solzi[i] * phi[i];
        ziOldg += solziOld[i] * phi[i];
        ziOlderg += solziOlder[i] * phi[i];
        ziM1g += solziM1[i] * phi[i];
        ziM1Oldg += solziM1Old[i] * phi[i];
        Tc += solTc[i] * phi[i];

        for(unsigned  j = 0; j < dim; j++) {
          gradTig[j] += solTi[i] * gradPhi[i * dim + j];
          gradlig[j] += solli[i] * gradPhi[i * dim + j];
          gradzig[j] += solzi[i] * gradPhi[i * dim + j];
          gradziOldg[j] += solziOld[i] * gradPhi[i * dim + j];
        }
      }

      if(!manifacturedSolution) Tc = GetTargetSolution(time);

      // std::cout.precision(14);
      // std::cout << ziOldg <<" ";

      double nu = 1. / Pe;

      // *** phiV_i loop ***
      for(unsigned i = 0; i < nDofs; i++) {

        adept::adouble Tieq = 0.;
        adept::adouble lieq = 0.;
        adept::adouble zieq = 0.;

        for(unsigned j = 0; j < dim; j++) {  // second index j in each equation
          Tieq   +=  gradPhi[i * dim + j] * (nu * gradlig[j] + beta * gradTig[j]); // weak laplace
          Tieq   +=  gradPhi[i * dim + j] * v[j] * lig; // weak convection

          lieq   +=  gradPhi[i * dim + j] * (nu * gradTig[j]); // weal laplace
          lieq   +=  phi[i] * v[j] * gradTig[j]; // strong convection

          zieq   +=  gradPhi[i * dim + j] * nu * (gradzig[j] + 0. * gradziOldg[j]);  // weak laplace
          zieq   +=  phi[i] * v[j] * gradzig[j]; // strong convection
        }
        Tieq += alpha * phi[i] * Tig;
        if(group == 5) Tieq += phi[i] * (Tig - Tc);

        if(!steady) {
          if(iext > 0) lieq += phi[i] * (ziM1g - ziM1Oldg) / dt;
          else lieq += 0.9 * phi[i] * (ziOldg - ziOlderg) / dt;
          zieq += phi[i] * (zig - ziOldg) / dt;
        }


        mResTi[i] += Tieq * weight;
        mResli[i] += lieq * weight;
        if(!nodeIsControlBoundary[i])  mReszi[i] += zieq * weight;
        else mReszi[i] += solzi[i] - solTi[i];

      } // end phiV_i loop
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube mRes in double Res and store them in RES
    Res.resize(nDofsAll);    //resize

    for(int i = 0; i < nDofs; i++) {
      Res[ i ] = -mResTi[i].value();
      Res[ nDofs + i ] = -mResli[i].value();
      Res[ 2 * nDofs + i ] = -mReszi[i].value();
    }

    RES->add_vector_blocked(Res, sysDof);

    //Extract and store the Jacobian
    Jac.resize(nDofsAll * nDofsAll);

    s.dependent(&mResTi[0], nDofs);
    s.dependent(&mResli[0], nDofs);
    s.dependent(&mReszi[0], nDofs);

    s.independent(&solTi[0], nDofs);
    s.independent(&solli[0], nDofs);
    s.independent(&solzi[0], nDofs);

    s.jacobian(&Jac[0], true);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();
  KK->close();

  //KK->draw();
// ***************** END ASSEMBLY *******************

}

void AssembleManifacturedSolution(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled



  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  TransientLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<TransientLinearImplicitSystem> ("sysTc");   // pointer to the linear implicit system named "Poisson"
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

  double dt =  mlPdeSys->GetIntervalTime();
  double time =  mlPdeSys->GetTime();

  unsigned TcIndex = mlSol->GetIndex("Tc");
  unsigned solType = mlSol->GetSolutionType(TcIndex);
  unsigned TcPdeIndex = mlPdeSys->GetSolPdeIndex("Tc");

  std::vector < adept::adouble > solTc;
  std::vector < double > solTcOld;

  std::vector < adept::adouble > mResTc;

  std::vector < std::vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector <double> phi;  // local test function for velocity
  std::vector <double> gradPhi; // local test function first order partial derivatives
  double weight;

  std::vector< unsigned > sysDof; // local to global pdeSys dofs
  std::vector< double > Res; // local redidual std::vector
  std::vector < double > Jac;

  RES->zero(); // Set to zero all the entries of the Global Residual std::vector
  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);

    unsigned nDofs = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs

    sysDof.resize(nDofs);

    for(unsigned  k = 0; k < dim; k++) {
      x[k].resize(nDofs);
    }
    solTc.resize(nDofs);
    solTcOld.resize(nDofs);
    mResTc.assign(nDofs, 0.);

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, solType);
      solTc[i] = (*sol->_Sol[TcIndex])(solDof);
      solTcOld[i] = (*sol->_SolOld[TcIndex])(solDof);
      sysDof[i] = pdeSys->GetSystemDof(TcIndex, TcPdeIndex, i, iel);
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);
      for(unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k])(xDof);
      }
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {
      msh->_finiteElement[ielGeom][solType]->Jacobian(x, ig, weight, phi, gradPhi);


      adept::adouble Tcg = 0;
      double TcOldg = 0;

      std::vector < adept::adouble > gradTcg(dim, 0);
      std::vector < double > gradTcOldg(dim, 0);

      for(unsigned i = 0; i < nDofs; i++) {
        Tcg += solTc[i] * phi[i];
        TcOldg += solTcOld[i] * phi[i];
        for(unsigned  j = 0; j < dim; j++) {
          gradTcg[j] += solTc[i] * gradPhi[i * dim + j];
          gradTcOldg[j] += solTcOld[i] * gradPhi[i * dim + j];
        }
      }

      // std::cout.precision(14);
      // std::cout << ziOldg <<" ";

      double nu = 1. / Pe;

      // *** phiV_i loop ***
      for(unsigned i = 0; i < nDofs; i++) {

        adept::adouble Tceq = 0.;

        for(unsigned j = 0; j < dim; j++) {  // second index j in each equation
          Tceq   +=  gradPhi[i * dim + j] * nu * (gradTcg[j] + 0. * gradTcOldg[j]);  // weak laplace
          Tceq   +=  phi[i] * v[j] * gradTcg[j]; // strong convection
        }

        if(!steady) {
          Tceq += phi[i] * (Tcg - TcOldg) / dt;
        }

        mResTc[i] += Tceq * weight;

      } // end phiV_i loop
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube mRes in double Res and store them in RES
    Res.resize(nDofs);    //resize
    for(int i = 0; i < nDofs; i++) {
      Res[ i ] = -mResTc[i].value();
    }

    RES->add_vector_blocked(Res, sysDof);

    //Extract and store the Jacobian
    Jac.resize(nDofs * nDofs);

    s.dependent(&mResTc[0], nDofs);

    s.independent(&solTc[0], nDofs);

    s.jacobian(&Jac[0], true);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();
  KK->close();

  //KK->draw();
// ***************** END ASSEMBLY *******************

}


void GetError(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled


  char zName[10];
  sprintf(zName, "z%d", iext);

  //char sName[10];
  //sprintf(sName, "top%d", iext);

  //  extract pointers to the several objects that we are going to use
  //TransientLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<TransientLinearImplicitSystem> (sName);   // pointer to the linear implicit system named "Poisson"
  TransientLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<TransientLinearImplicitSystem> ("TmOptCtr");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable

  double time =  mlPdeSys->GetTime();

  unsigned TcIndex = mlSol->GetIndex("Tc");
  unsigned ziIndex = mlSol->GetIndex(zName);
  unsigned solType = mlSol->GetSolutionType(ziIndex);

  std::vector < double > solTc;
  std::vector < double > solzi;

  std::vector < std::vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector <double> phi, gradPhi;  // local test function for velocity
  double weight;

  double  localIntegral[3] = {0, 0, 0};

  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned group = msh->GetElementGroup(iel);
    if(group == 5) {

      short unsigned ielGeom = msh->GetElementType(iel);
      unsigned nDofs = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
      for(unsigned  k = 0; k < dim; k++) {
        x[k].resize(nDofs);
      }
      solzi.resize(nDofs);
      solTc.resize(nDofs);

      // local storage of global mapping and solution
      for(unsigned i = 0; i < nDofs; i++) {
        unsigned solDof = msh->GetSolutionDof(i, iel, solType);
        solTc[i] = (*sol->_Sol[TcIndex])(solDof);
        solzi[i] = (*sol->_Sol[ziIndex])(solDof);
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

        double zig = 0;
        double Tc = 0;

        for(unsigned i = 0; i < nDofs; i++) {
          zig += solzi[i] * phi[i];
          Tc += solTc[i] * phi[i];
        }
        if(!manifacturedSolution)Tc = GetTargetSolution(time);

        localIntegral[0] += Tc * weight;
        localIntegral[1] += zig * weight;
        localIntegral[2] += (zig - Tc) * (zig - Tc) * weight;

      } // end gauss point loop

    }
  } //end element loop for each process

  double globalIntegral[3] = {0., 0., 0.};

  MPI_Reduce(localIntegral, globalIntegral, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


  if(iproc == 0) {

    std::string filename = "SolutionIntegral" + std::to_string(iext) + ".txt";
    std::ofstream fout;
    if(fabs(time) < 1.0e-10) {
      fout.open(filename.c_str(), std::ios::out);
    }
    else {
      fout.open(filename.c_str(), std::ios::app);
    }
    fout << time << " " << globalIntegral[0] << " " << globalIntegral[1] << " " << globalIntegral[2] << " " << std::endl;
    fout.close();
  }


// ***************** END ASSEMBLY *******************

}







