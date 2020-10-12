
#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "TransientSystem.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "Marker.hpp"
#include "Line.hpp"

#include "Fluid.hpp"
#include "Solid.hpp"
#include "Parameter.hpp"

#include "NumericVector.hpp"
#include "adept.h"

using namespace femus;
Line *bulk;
Line *lineI;
void BuildFlag(MultiLevelSolution & mlSol);
double eps;

#include "../../Nitsche/support/particleInit.hpp"
#include "../../Nitsche/support/sharedFunctions.hpp"
#include "../include/mpmFsi8.hpp"
using namespace femus;

double dt = 1.;

double SetVariableTimeStep(const double time) {
  if(time < 2.5) dt = 0.05;
  else dt = 0.0125;

  return dt;
}

void BuildFlagSolidRegion(MultiLevelSolution & mlSol);
void ProjectNewmarkDisplacemenet(MultiLevelSolution & mlSol);

bool SetBoundaryCondition(const std::vector < double >&x, const char name[], double &value, const int facename, const double t) {
  bool dirichlet = true;
  value = 0.;

  const double Ubar = 1.0;    // guess?
  const double L = 0.41;
  const double H = 2.5;


  if(!strcmp(name, "UY")) {
    if(1 == facename) {     //inflow
      if(t < 2.0) {
        value = 1.5 * Ubar * 4.0 / 0.1681 * (x[0] + 0.2) * (-x[0] + 0.21) * 0.5 * (1. - cos(0.5 * M_PI * t));
      }
      else {
        value = 1.5 * Ubar * 4.0 / 0.1681 * (x[0] + 0.2) * (-x[0] + 0.21);
      }
    }
    else if(2 == facename || 6 == facename) {     //outflow and porous media Neumann
      dirichlet = false;
      value = 0.;
    }
  }

  else if(!strcmp(name, "UX")) {
    if(6 == facename) {     //porous media Neumann
      dirichlet = false;
      value = 0.;
    }
  }

  else if(!strcmp(name, "P")) {
    dirichlet = 0;
    value = 0;
  }

  return dirichlet;

}


int main(int argc, char **args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;
  unsigned numberOfUniformLevels = 1; //for refinement in 3D
  //unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0 ;

  double Lref = 1.;
  double Uref = 1.;
//   double rhos = 7850;
//   double rhof = 1000;
//   double nu = 0.3;
//   double E = 2.e05;
//   double muf = 1.0e-3;



  double rhof = 1000.;
  double muf = 1.;
  double rhos = 10000.;
  double nu = 0.4;
  double E = 1400000;


  beta = 0.3;
  Gamma = 0.5;


  Parameter par(Lref, Uref);

  // Generate Solid Object
  Solid solid(par, E, nu, rhos, "Neo-Hookean");
  Fluid fluid(par, muf, rhof, "Newtonian");

  mlMsh.ReadCoarseMesh("../input/turek2Db.neu", "fifth", scalingFactor);
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels,
                   numberOfUniformLevels, NULL);

  mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);
  numberOfUniformLevels = 1;

  unsigned dim = mlMsh.GetDimension();

  FEOrder femOrder = SECOND;

  MultiLevelSolution mlSol(&mlMsh);
  // add variables to mlSol
  mlSol.AddSolution("UX", LAGRANGE, femOrder, 2);
  mlSol.AddSolution("UY", LAGRANGE, femOrder, 2);
  if(dim > 2) mlSol.AddSolution("UZ", LAGRANGE, femOrder, 2);

  //mlSol.AddSolution ("P", DISCONTINUOUS_POLYNOMIAL, ZERO, 2);
  //mlSol.AddSolution ("P", DISCONTINUOUS_POLYNOMIAL, FIRST, 2);
  mlSol.AddSolution("P", LAGRANGE, FIRST, 2);

  mlSol.AddSolution("eflag", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
  mlSol.AddSolution("nflag", LAGRANGE, femOrder, 0, false);

  mlSol.AddSolution("DX", LAGRANGE, femOrder, 0, false);
  mlSol.AddSolution("DY", LAGRANGE, femOrder, 0, false);
  if(dim > 2) mlSol.AddSolution("DZ", LAGRANGE, femOrder, 0, false);

  mlSol.AddSolution("VX", LAGRANGE, femOrder, 0, false);
  mlSol.AddSolution("VY", LAGRANGE, femOrder, 0, false);
  if(dim > 2) mlSol.AddSolution("VZ", LAGRANGE, femOrder, 0, false);

  mlSol.AddSolution("AX", LAGRANGE, femOrder, 0, false);
  mlSol.AddSolution("AY", LAGRANGE, femOrder, 0, false);
  if(dim > 2) mlSol.AddSolution("AZ", LAGRANGE, femOrder, 0, false);

  
  
  //mlSol.SetIfFSI(true);

  mlSol.Initialize("All");

  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  //mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition1);

  // ******* Set boundary conditions *******
  mlSol.GenerateBdc("UX", "Steady");
  mlSol.GenerateBdc("UY", "Time_dependent");
  if(dim > 2)
    mlSol.GenerateBdc("UZ", "Steady");

  mlSol.GenerateBdc("P", "Steady");

  MultiLevelProblem ml_prob(&mlSol);

  ml_prob.parameters.set < Solid > ("SolidMPM") = solid;
  ml_prob.parameters.set < Solid > ("SolidFEM") = solid;
  ml_prob.parameters.set < Fluid > ("FluidFEM") = fluid;

  // ******* Add MPM system to the MultiLevel problem *******
  TransientNonlinearImplicitSystem & system =
    ml_prob.add_system < TransientNonlinearImplicitSystem > ("MPM_FSI");
  system.AddSolutionToSystemPDE("UX");
  system.AddSolutionToSystemPDE("UY");
  if(dim > 2)
    system.AddSolutionToSystemPDE("UZ");
  system.AddSolutionToSystemPDE("P");

  system.SetSparsityPatternMinimumSize(1000);
//   if(dim > 1) system.SetSparsityPatternMinimumSize (500, "VY");
//   if(dim > 2) system.SetSparsityPatternMinimumSize (500, "VZ");
//

  // ******* System MPM-FSI Assembly *******
  system.SetAssembleFunction(AssembleMPMSys);
  //system.SetAssembleFunction (AssembleMPMSysOld);
  //system.SetAssembleFunction(AssembleFEM);
  // ******* set MG-Solver *******
  system.SetMgType(V_CYCLE);

  system.SetAbsoluteLinearConvergenceTolerance(1.0e-10);
  system.SetMaxNumberOfLinearIterations(1);
  system.SetNonLinearConvergenceTolerance(1.e-9);
  system.SetMaxNumberOfNonLinearIterations(2);

  system.SetNumberPreSmoothingStep(1);
  system.SetNumberPostSmoothingStep(1);

  // ******* Set Preconditioner *******
  system.SetLinearEquationSolverType(FEMuS_DEFAULT);
  system.init();

  // ******* Set Smoother *******
  system.SetSolverFineGrids(GMRES);
  system.SetPreconditionerFineGrids(ILU_PRECOND);
  system.SetTolerances(1.e-10, 1.e-15, 1.e+50, 40, 40);

  BuildFlagSolidRegion(mlSol);

  mlSol.SetWriter(VTK);

  std::vector < std::string > mov_vars;
  mov_vars.push_back("DX");
  mov_vars.push_back("DY");
  //mov_vars.push_back("DZ");
  mlSol.GetWriter()->SetMovingMesh(mov_vars);

  std::vector < std::string > print_vars;
  print_vars.push_back("All");

  mlSol.GetWriter()->SetDebugOutput(true);

  ProjectNewmarkDisplacemenet(mlSol);
  mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, 0);
  
  system.AttachGetTimeIntervalFunction(SetVariableTimeStep);
  unsigned n_timesteps = 10000;
  for(unsigned time_step = 1; time_step <= n_timesteps; time_step++) {
       
    system.CopySolutionToOldSolution();
    system.MGsolve();
    ProjectNewmarkDisplacemenet(mlSol);
    mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, time_step);
  }

  return 0;

}



void BuildFlagSolidRegion(MultiLevelSolution & mlSol) {

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution *sol = mlSol.GetSolutionLevel(level);
  Mesh *msh = mlSol._mlMesh->GetLevel(level);
  unsigned iproc = msh->processor_id();

  unsigned eflagIndex = mlSol.GetIndex("eflag");
  unsigned nflagIndex = mlSol.GetIndex("nflag");

  unsigned nflagType = mlSol.GetSolutionType(nflagIndex);

  sol->_Sol[eflagIndex]->zero();
  sol->_Sol[nflagIndex]->zero();

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    unsigned flag_mat = msh->GetElementMaterial(iel);
    sol->_Sol[eflagIndex]->set(iel, flag_mat);
    if(flag_mat == 4) {
      unsigned nDofs = msh->GetElementDofNumber(iel, nflagType);
      for(unsigned  i = 0; i < nDofs; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, nflagType);
        sol->_Sol[nflagIndex]->set(idof, 6);
      }
    }
  }
  sol->_Sol[eflagIndex]->close();
  sol->_Sol[nflagIndex]->close();

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    unsigned flag_mat = msh->GetElementMaterial(iel);
    sol->_Sol[eflagIndex]->set(iel, flag_mat);
    if(flag_mat == 3) {
      unsigned nDofs = msh->GetElementDofNumber(iel, nflagType);
      for(unsigned  i = 0; i < nDofs; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, nflagType);

        unsigned nodeFlag = (*sol->_Sol[nflagIndex])(idof);
        if(nodeFlag == 0) {
          sol->_Sol[nflagIndex]->set(idof, 4);
        }
        else if(nodeFlag == 6) {
          sol->_Sol[nflagIndex]->set(idof, 5);
        }


      }
    }
  }
  sol->_Sol[nflagIndex]->close();
}



void ProjectNewmarkDisplacemenet(MultiLevelSolution & mlSol) {

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution *sol = mlSol.GetSolutionLevel(level);
  Mesh *msh = mlSol._mlMesh->GetLevel(level);
  unsigned iproc = msh->processor_id();

  const unsigned dim = msh->GetDimension();

  unsigned nflagIndex = mlSol.GetIndex("nflag");

  const char varname[12][5] = {"UX", "UY", "UZ", "DX", "DY", "DZ", "VX", "VY", "VZ", "AX", "AY", "AZ"};

  vector <unsigned> indexSolU(dim);
  vector <unsigned> indexSolD(dim);
  vector <unsigned> indexSolV(dim);
  vector <unsigned> indexSolA(dim);
  for(unsigned k = 0; k < dim; k++) {
    indexSolU[k] = mlSol.GetIndex(&varname[k][0]);
    indexSolD[k] = mlSol.GetIndex(&varname[k + 3][0]);
    indexSolV[k] = mlSol.GetIndex(&varname[k + 6][0]);
    indexSolA[k] = mlSol.GetIndex(&varname[k + 9][0]);
  }
  unsigned solType = mlSol.GetSolutionType(&varname[0][0]);

  for(int i = msh->_dofOffset[solType][iproc]; i < msh->_dofOffset[solType][iproc + 1]; i++) {
    if((*sol->_Sol[nflagIndex])(i) >= 4) {
      for(unsigned k = 0 ; k < dim ; k++) {
          
        double Dnew = (*sol->_Sol[indexSolU[k]])(i);
        double Dold = (*sol->_SolOld[indexSolU[k]])(i);
        
        double Vold = (*sol->_Sol[indexSolV[k]])(i);
        double Aold = (*sol->_Sol[indexSolA[k]])(i);
        double Anew = (Dnew - Dold) / (beta * dt * dt) - Vold / (beta * dt) + ((beta - 0.5) * Aold) / beta;
        double Vnew = Vold + (1 - Gamma) * dt * Aold + Gamma * dt * Anew;
                
        sol->_Sol[indexSolD[k]]->set(i, Dnew);
        sol->_Sol[indexSolV[k]]->set(i, Vnew);          
        sol->_Sol[indexSolA[k]]->set(i, Anew);
          
      }
    }
  }

  for(unsigned k = 0 ; k < dim ; k++) {
    sol->_Sol[indexSolD[k]]->close();
    sol->_Sol[indexSolV[k]]->close();
    sol->_Sol[indexSolA[k]]->close();
  }
}



void BuildFlag(MultiLevelSolution & mlSol) {

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;
  Solution *sol = mlSol.GetSolutionLevel(level);
  Mesh *msh = mlSol._mlMesh->GetLevel(level);
  unsigned iproc = msh->processor_id();

  unsigned eflagIndex = mlSol.GetIndex("eflag");
  unsigned nflagIndex = mlSol.GetIndex("nflag");

  unsigned nflagType = mlSol.GetSolutionType(nflagIndex);

  std::vector < Marker * >particle3 = bulk->GetParticles();
  std::vector < unsigned >markerOffset3 = bulk->GetMarkerOffset();
  unsigned imarker3 = markerOffset3[iproc];

  sol->_Sol[eflagIndex]->zero();
  sol->_Sol[nflagIndex]->zero();

  for(int iel = msh->_elementOffset[iproc];
      iel < msh->_elementOffset[iproc + 1]; iel++) {

    bool inside = false;
    bool outside = false;
    bool interface = false;

    while(imarker3 < markerOffset3[iproc + 1]
          && iel > particle3[imarker3]->GetMarkerElement()) {
      imarker3++;
    }
    while(imarker3 < markerOffset3[iproc + 1]
          && iel == particle3[imarker3]->GetMarkerElement()) {
      double dg1 = particle3[imarker3]->GetMarkerDistance();
      if(dg1 < -1.0e-10) {
        outside = true;
        if(inside) {
          interface = true;
          break;
        }
      }
      else if(dg1 > 1.0e-10) {
        inside = true;
        if(outside) {
          interface = true;
          break;
        }
      }
      else {
        interface = true;
        break;
      }
      imarker3++;
    }
    if(interface) {
      sol->_Sol[eflagIndex]->set(iel, 1.);
      unsigned nDofu = msh->GetElementDofNumber(iel, nflagType);      // number of solution element dofs
      for(unsigned i = 0; i < nDofu; i++) {
        unsigned iDof = msh->GetSolutionDof(i, iel, nflagType);
        sol->_Sol[nflagIndex]->set(iDof, 1.);
      }
    }
    else if(inside) {
      sol->_Sol[eflagIndex]->set(iel, 2.);
    }
  }
  sol->_Sol[eflagIndex]->close();
  sol->_Sol[nflagIndex]->close();
}


