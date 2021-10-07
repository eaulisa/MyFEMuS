// This is a benchmark problem taken from 2007 Klaus-Jürgen Bathe
// 2D and 3D FSI patch tests. Solid is placed on the top of the fluid
//constant traction is applied from the bottom and slip on the other walls.
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
Line* bulk;
Line* lineI;
Line* intPoints;
void BuildFlag(MultiLevelSolution& mlSol);


void BuidProjection(MultiLevelProblem& ml_prob);
void ProjectSolutionIntoGradient(MultiLevelProblem& ml_prob);

bool NeoHookean = false;
bool particleSmoothingIsOn = false;

double eps = 0.00;

double gravity[3] = {0., 0., 0.};
bool weakP = false;

double theta = .8;
double af = 1. - theta;
double am = af - 0.1;
double beta = 0.25 + 0.5 * (af - am);
double Gamma = 0.5 + (af - am);
// double beta = .5;
// double Gamma = 1.;

double DTMIN = 0.001;

double factor = 1.;

double gammacF = factor * 0.05;
double gammacS = factor * 0.05;
double gammap = factor * 0.05;
double gammau = 0.05 * gammacF;

double GAMMA = factor * 45;//45;   // 10, 45 in the paper.

#include "./include/mpmFsi15.hpp"
using namespace femus;

double dt = 0.001;
double SetVariableTimeStep(const double time) {
    
//   double dt0 = 1. / 360;
//   double dt1 = 1. / 360; //FSI3
// 
//   double T = 2.;
//   if(time < T)
//     dt = (0.5 * (1. + cos(M_PI * time / T))) * dt0    + (0.5 * (1. - cos(M_PI * time / T))) * dt1;
//   else
//     dt = dt1;

  return dt;
}


bool SetBoundaryCondition(const std::vector < double >&x, const char name[], double &value, const int facename, const double t) {
  bool test = 1;      //dirichlet
  value = 0.;      

if(!strcmp(name, "DY")) {
    if(3 == facename) {
      test = 0;   // neumann do-nothing
    }
  }
  else if(!strcmp(name, "DX")) {
    if(1 == facename || 2 == facename) {
      test = 0; // neumann do-nothing
    }
  }
  else if(!strcmp(name, "VX")) {
    if(1 == facename || 2 == facename) {
      test = 0; // neumann do-nothing
    }
  }
  else if(!strcmp(name, "VY")) {
    if(3 == facename) {
      test = 0; // neumann do-nothing
    }
  }
  else if(!strcmp(name, "P")) {
    if(facename == 1){
        value = 0.1;
    }
    else{
        test = 0.;
    }
  }

  return test;

}

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  MultiLevelMesh mlMsh;

  double scalingFactor = 1.;
  unsigned numberOfUniformLevels = 4; //for refinement in 3D
  unsigned numberOfUniformLevelsStart = numberOfUniformLevels;
  //unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;


  double Lref = 1.;
  double Uref = 1.;
  double rhof = 1000.;
  double muf = 0.001;
  double rhos = 800.;
  double nu = 0.4999;
  double E = 200; //FSI3


  Parameter par(Lref, Uref);

  // Generate Solid Object
  Solid solid(par, E, nu, rhos, "Neo-Hookean"); // Mooney-Rivlin?
  Fluid fluid(par, muf, rhof, "Newtonian");

  mlMsh.ReadCoarseMesh("../input/benchmarkFSI.neu", "fifth", 1); // 
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels, NULL);

  mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);
  numberOfUniformLevels = 1;

  unsigned dim = mlMsh.GetDimension();

  FEOrder femOrder = FIRST;

  MultiLevelSolution mlSol(&mlMsh);
  // add variables to mlSol
  mlSol.AddSolution("DX", LAGRANGE, femOrder, 2);
  if(dim > 1) mlSol.AddSolution("DY", LAGRANGE, femOrder, 2);
  if(dim > 2) mlSol.AddSolution("DZ", LAGRANGE, femOrder, 2);

  mlSol.AddSolution("VX", LAGRANGE, femOrder, 2);
  if(dim > 1) mlSol.AddSolution("VY", LAGRANGE, femOrder, 2);
  if(dim > 2) mlSol.AddSolution("VZ", LAGRANGE, femOrder, 2);

  //mlSol.AddSolution ("P", DISCONTINUOUS_POLYNOMIAL, ZERO, 2);
  //mlSol.AddSolution ("P", DISCONTINUOUS_POLYNOMIAL, FIRST, 2);
  mlSol.AddSolution("P", LAGRANGE, FIRST, 2);

  mlSol.AddSolution("eflag", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
  mlSol.AddSolution("nflag", LAGRANGE, femOrder, 0, false);


  std::string Uxname[3][3] = {"DXx", "DXy", "DXz", "DYx", "DYy", "DYz", "DZx", "DZy", "DZz"};
  for(unsigned j = 0; j < dim; j++) {
    for(unsigned k = 0; k < dim; k++) {
      mlSol.AddSolution(Uxname[j][k].c_str(), LAGRANGE, femOrder, 0, false);
    }
  }
  mlSol.AddSolution("weight", LAGRANGE, femOrder, 0, false);
  mlSol.SetIfFSI(true);
  mlSol.Initialize("All");
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

  // ******* Set boundary conditions *******
  mlSol.GenerateBdc("DX", "Time_dependent");
  if(dim > 1) mlSol.GenerateBdc("DY", "Time_dependent");
  if(dim > 2) mlSol.GenerateBdc("DZ", "Steady");
  mlSol.GenerateBdc("VX", "Time_dependent");
  if(dim > 1) mlSol.GenerateBdc("VY", "Time_dependent");
  if(dim > 2) mlSol.GenerateBdc("VZ", "Steady");
  mlSol.GenerateBdc("P", "Steady");


  ////////////////////////////////////////////////////////////////


  MultiLevelProblem ml_prob(&mlSol);


  LinearImplicitSystem* systemP[3];
  std::string Pname[3] = {"Px", "Py", "Pz"};
  for(unsigned k = 0; k < dim; k++) {
    systemP[k] = &ml_prob.add_system < LinearImplicitSystem > (Pname[k]);
    systemP[k]->AddSolutionToSystemPDE("DX");
    systemP[k]->init();
  }

  BuidProjection(ml_prob);

  ml_prob.parameters.set<Solid> ("SolidMPM") = solid;
  ml_prob.parameters.set<Solid> ("SolidFEM") = solid;
  ml_prob.parameters.set<Fluid> ("FluidFEM") = fluid;

  // ******* Add MPM system to the MultiLevel problem *******
  TransientNonlinearImplicitSystem& system = ml_prob.add_system < TransientNonlinearImplicitSystem > ("MPM_FSI");
  system.AddSolutionToSystemPDE("DX");
  if(dim > 1) system.AddSolutionToSystemPDE("DY");
  if(dim > 2) system.AddSolutionToSystemPDE("DZ");
  system.AddSolutionToSystemPDE("VX");
  if(dim > 1) system.AddSolutionToSystemPDE("VY");
  if(dim > 2) system.AddSolutionToSystemPDE("VZ");
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
  system.SetMaxNumberOfNonLinearIterations(5);

  system.SetNumberPreSmoothingStep(1);
  system.SetNumberPostSmoothingStep(1);

  // ******* Set Preconditioner *******
  system.SetLinearEquationSolverType(FEMuS_DEFAULT);

  system.init();

  // ******* Set Smoother *******
  system.SetSolverFineGrids(GMRES);

  system.SetPreconditionerFineGrids(ILU_PRECOND);

  system.SetTolerances(1.e-10, 1.e-15, 1.e+50, 40, 40);


  //BuildIntegrationPoints(mlSol);

  std::ifstream fin;
  std::ostringstream fileName;
  std::ostringstream level_number;
  fileName << "../input/benchmarkFSISolid";
  level_number << 6;

  fileName << level_number.str();

  //BEGIN bulk reading

  std::string bulkfile = fileName.str();
  bulkfile += ".bulk.txt";

  std::vector < std::vector <double> > xp;
  std::vector <double> wp;
  std::vector <double> dist;
  std::vector < MarkerType > markerType;

  fin.open(bulkfile);



  unsigned size;
  fin >> dim >> size;



  //std::cout << dim << " " << size << std::endl << std::flush;
  xp.resize(size);
  wp.resize(size);
  dist.resize(size);
  markerType.assign(size, VOLUME);
  for(unsigned ip = 0; ip < xp.size(); ip++) {
    xp[ip].resize(dim);
    for(unsigned k = 0; k < dim; k++) {
      fin >> xp[ip][k];
    }
    fin >> wp[ip];
    fin >> dist[ip];
  }
  fin.close();


  double delta_max = 0.75 / scalingFactor * 1.5 / numberOfUniformLevelsStart;
  //double delta_max = 0.005 / (numberOfUniformLevelsStart - 3);

  for(int i = 0; i < xp.size(); i++) {
    if(dist[i] < -delta_max) {
      xp.erase(xp.begin() + i);
      wp.erase(wp.begin() + i);
      dist.erase(dist.begin() + i);
      markerType.erase(markerType.begin() + i);
      i--;
    }
  }

  double shift = .0e-4;
  for(int i = 0; i < xp.size(); i++) {
    xp[i][0] += shift;
  }



  bulk = new Line(xp, wp, dist, markerType, mlSol.GetLevel(numberOfUniformLevels - 1), 2);

  std::vector < std::vector < std::vector < double > > >  bulkPoints(1);
  bulk->GetLine(bulkPoints[0]);
  PrintLine(DEFAULT_OUTPUTDIR, "bulk", bulkPoints, 0);

  //END bulk reading

  //BEGIN interface reading

  std::string interfacefile = fileName.str();
  interfacefile += ".interface.txt";

  std::vector < std::vector < std::vector < double > > > T;
  fin.open(interfacefile);
  fin >> dim >> size;
  xp.resize(size);
  T.resize(size);
  markerType.assign(size, INTERFACE);
  for(unsigned ip = 0; ip < xp.size(); ip++) {
    xp[ip].resize(dim);
    for(unsigned k = 0; k < dim; k++) {
      fin >> xp[ip][k];
    }
    T[ip].resize(dim - 1);
    for(unsigned l = 0; l < dim - 1; l++) {
      T[ip][l].resize(dim);
      for(unsigned k = 0; k < dim; k++) {
        fin >> T[ip][l][k];
      }
    }
  }
  fin.close();

  for(int i = 0; i < xp.size(); i++) {
    xp[i][0] += shift;
  }
  lineI = new Line(xp, T, markerType, mlSol.GetLevel(numberOfUniformLevels - 1), 2);

  std::vector < std::vector < std::vector < double > > > lineIPoints(1);
  lineI->GetLine(lineIPoints[0]);

  PrintLine(DEFAULT_OUTPUTDIR, "interfaceMarkers", lineIPoints, 0);



//END interface reading

  unsigned iproc = mlMsh.GetLevel(0)->processor_id();

  lineI->GetParticlesToGridMaterial(false);
  bulk->GetParticlesToGridMaterial(false);

  BuildFlag(mlSol);

// ******* Print solution *******
  mlSol.SetWriter(VTK);

  std::vector<std::string> mov_vars;
  mov_vars.push_back("DX");
  mov_vars.push_back("DY");
  if(dim == 3) mov_vars.push_back("DZ");
  mlSol.GetWriter()->SetMovingMesh(mov_vars);

  std::vector<std::string> print_vars;
  print_vars.push_back("All");

  mlSol.GetWriter()->SetDebugOutput(true);
  mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, 0);
  mlSol.GetWriter()->Write("./output1", "biquadratic", print_vars, 0);


  system.AttachGetTimeIntervalFunction(SetVariableTimeStep);
  unsigned n_timesteps = 2000;
  unsigned printTimeInterval = 1;
  for(unsigned time_step = 1; time_step <= n_timesteps; time_step++) {

    system.CopySolutionToOldSolution();

    mlSol.GetWriter()->Write("output1", "biquadratic", print_vars, time_step / printTimeInterval);
    bulk->GetLine(bulkPoints[0]);
    PrintLine("output1", "bulk", bulkPoints, time_step);
    lineI->GetLine(lineIPoints[0]);
    PrintLine("output1", "interfaceMarkers", lineIPoints, time_step / printTimeInterval);

    system.MGsolve();
    ProjectSolutionIntoGradient(ml_prob);

    double time = system.GetTime();

    // GetDragAndLift(ml_prob, time, pfile);

    mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, time_step / printTimeInterval);
    GridToParticlesProjection(ml_prob, *bulk, *lineI);

    bulk->GetLine(bulkPoints[0]);
    PrintLine(DEFAULT_OUTPUTDIR, "bulk", bulkPoints, time_step / printTimeInterval);
    lineI->GetLine(lineIPoints[0]);
    PrintLine(DEFAULT_OUTPUTDIR, "interfaceMarkers", lineIPoints, time_step / printTimeInterval);

  }


  delete bulk;
  delete lineI;

  return 0;

} //end main





void BuildFlag(MultiLevelSolution & mlSol) {

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution *sol  = mlSol.GetSolutionLevel(level);
  Mesh     *msh   = mlSol._mlMesh->GetLevel(level);
  unsigned iproc  = msh->processor_id();

  unsigned eflagIndex = mlSol.GetIndex("eflag");
  unsigned nflagIndex = mlSol.GetIndex("nflag");

  unsigned nflagType = mlSol.GetSolutionType(nflagIndex);

  std::vector<Marker*> particle3 = bulk->GetParticles();
  std::vector<unsigned> markerOffset3 = bulk->GetMarkerOffset();
  unsigned imarker3 = markerOffset3[iproc];

  sol->_Sol[eflagIndex]->zero();
  sol->_Sol[nflagIndex]->zero();

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    bool inside = false;
    bool outside = false;
    bool interface = false;

    while(imarker3 < markerOffset3[iproc + 1] && iel > particle3[imarker3]->GetMarkerElement()) {
      imarker3++;
    }
    while(imarker3 < markerOffset3[iproc + 1] && iel == particle3[imarker3]->GetMarkerElement()) {
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
      unsigned nDofu  = msh->GetElementDofNumber(iel, nflagType);  // number of solution element dofs
      for(unsigned i = 0; i < nDofu; i++) {
        unsigned iDof = msh->GetSolutionDof(i, iel, nflagType);
        sol->_Sol[nflagIndex]->set(iDof, 1.);
      }
    }
    else if(inside) {
      sol->_Sol[eflagIndex]->set(iel, 2.);
    }
    else if(outside) {
      sol->_Sol[eflagIndex]->set(iel, .5);
    }
  }
  sol->_Sol[eflagIndex]->close();
  sol->_Sol[nflagIndex]->close();
}

void BuidProjection(MultiLevelProblem& ml_prob) {

  adept::Stack& s = FemusInit::_adeptStack;

  MultiLevelSolution*  mlSol = ml_prob._ml_sol;
  unsigned level = mlSol->_mlMesh->GetNumberOfLevels() - 1u;

  Solution* sol = ml_prob._ml_sol->GetSolutionLevel(level);
  Mesh* msh = ml_prob._ml_msh->GetLevel(level);
  elem* el = msh->el;

  unsigned  dim = msh->GetDimension();

  std::vector < LinearImplicitSystem* > mlSysP(dim);
  std::vector < LinearEquationSolver* > sysP(dim);
  std::vector < SparseMatrix*> P(dim);

  std::string Pname[3] = {"Px", "Py", "Pz"};
  std::string Uname = "DX";
  unsigned soluIndex = mlSol->GetIndex(Uname.c_str());
  for(unsigned k = 0; k < dim; k++) {
    mlSysP[k] =  &ml_prob.get_system< LinearImplicitSystem > (Pname[k]);
    sysP[k] = mlSysP[k]->_LinSolver[level];
    P[k] = sysP[k]->_KK;
    P[k]->zero();
  }

  vector < vector < double > > x(dim);
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector< int > sysDof;
  vector <double> phi;
  double* phi2;
  vector <double> phi_x;
  double weight;

  unsigned solwIndex = mlSol->GetIndex("weight");
  unsigned solType = mlSol->GetSolutionType(solwIndex);

  sol->_Sol[solwIndex]->zero();

  unsigned    iproc = msh->processor_id();

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofs  = msh->GetElementDofNumber(iel, solType);
    sysDof.resize(nDofs);
    for(int k = 0; k < dim; k++) {
      x[k].resize(nDofs);
    }
    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, solType);
      sysDof[i] = msh->GetSolutionDof(i, iel, solType);
    }
    // local storage of coordinates
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);
      for(unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k])(xDof);
      }
    }

    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {
      msh->_finiteElement[ielGeom][solType]->Jacobian(x, ig, weight, phi, phi_x);
      phi2 = (solType != 1) ? msh->_finiteElement[ielGeom][solType]->GetPhi(ig) : msh->_finiteElement[ielGeom][2]->GetPhi(ig);
      // *** phi_i loop ***
      for(unsigned i = 0; i < nDofs; i++) {
        sol->_Sol[solwIndex]->add(sysDof[i], phi2[i] * weight);
      } // end phi_i loop
    } // end gauss point loop

  } //end element loop for each process*/
  sol->_Sol[solwIndex]->close();


  //solution variable


  {
    unsigned soluType = mlSol->GetSolutionType(soluIndex);
    if(soluType != solType) {
      std::cout << "error weight and u should be of the same type\n";
      abort();
    }
  }
  unsigned soluPdeIndex = mlSysP[0]->GetSolPdeIndex("DX");
  std::vector < adept::adouble > solu;

  std::vector <double> Jac;
  std::vector < std::vector< adept::adouble > > aRes(dim);  // local redidual vector

  std::vector < double > solw;

  //BEGIN element loop
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofs  = msh->GetElementDofNumber(iel, solType);
    solu.resize(nDofs);
    solw.resize(nDofs);
    sysDof.resize(nDofs);
    for(int k = 0; k < dim; k++) {
      x[k].resize(nDofs);
      aRes[k].assign(nDofs, 0.);
    }
    Jac.resize(nDofs * nDofs);
    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, solType);
      solu[i] = (*sol->_Sol[soluIndex])(solDof);
      solw[i] = (*sol->_Sol[solwIndex])(solDof);
      sysDof[i] = sysP[0]->GetSystemDof(soluIndex, soluPdeIndex, i, iel);
    }
    // local storage of coordinates
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);
      for(unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k])(xDof);
      }
    }

    s.new_recording();
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {

      msh->_finiteElement[ielGeom][solType]->Jacobian(x, ig, weight, phi, phi_x);
      phi2 = (solType != 1) ? msh->_finiteElement[ielGeom][solType]->GetPhi(ig) : msh->_finiteElement[ielGeom][2]->GetPhi(ig);
      std::vector < adept::adouble > solux_g(dim, 0.);
      for(unsigned i = 0; i < nDofs; i++) {
        for(unsigned k = 0; k < dim; k++) {
          solux_g[k] += phi_x[i * dim + k] * solu[i];
        }
      }

      // *** phi_i loop ***
      for(unsigned i = 0; i < nDofs; i++) {
        for(unsigned k = 0; k < dim; k++) {
          aRes[k][i] += solux_g[k] * phi2[i] * weight / solw[i]; //smoothed gradient part.
        }
      } // end phi_i loop
    } // end gauss point loop

    s.independent(&solu[0], nDofs);
    for(unsigned k = 0; k < dim; k++) {
      s.dependent(&aRes[k][0], nDofs);
      s.jacobian(&Jac[0], true);
      P[k]->add_matrix_blocked(Jac, sysDof, sysDof);
      s.clear_dependents();
    }
    s.clear_independents();
  } //end element loop for each process*/

  for(unsigned k = 0; k < dim; k++) {
    P[k]->close();
  }

}


void ProjectSolutionIntoGradient(MultiLevelProblem& ml_prob) {

  MultiLevelSolution*  mlSol = ml_prob._ml_sol;
  unsigned level = mlSol->_mlMesh->GetNumberOfLevels() - 1u;

  Solution* sol = ml_prob._ml_sol->GetSolutionLevel(level);
  Mesh* msh = ml_prob._ml_msh->GetLevel(level);
  elem* el = msh->el;

  unsigned  dim = msh->GetDimension();

  std::vector < LinearImplicitSystem* > mlSysP(dim + 1);
  std::vector < LinearEquationSolver* > sysP(dim + 1);
  std::vector < SparseMatrix*> P(dim + 1);

  std::string Pname[3] = {"Px", "Py", "Pz"};
  std::string Uxname[3][3] = {"DXx", "DXy", "DXz", "DYx", "DYy", "DYz", "DZx", "DZy", "DZz"};

  std::vector < std::vector < unsigned > > solIndex(dim);
  for(unsigned j = 0; j < dim; j++) {
    solIndex[j].resize(dim);
    for(unsigned k = 0; k < dim; k++) {
      mlSysP[k] =  &ml_prob.get_system< LinearImplicitSystem > (Pname[k]);
      sysP[k] = mlSysP[k]->_LinSolver[level];
      solIndex[j][k] = mlSol->GetIndex(Uxname[j][k].c_str());
      P[k] = sysP[k]->_KK;
    }
  }

  std::string Uname[3] = {"DX", "DY", "DZ"};

  std::vector < unsigned > soluIndex(dim);
  for(unsigned k = 0; k < dim; k++) {
    soluIndex[k] = mlSol->GetIndex(Uname[k].c_str());
  }

  for(unsigned j = 0; j < dim; j++) {
    for(unsigned k = 0; k < dim; k++) {
      (*sol->_Sol[solIndex[j][k]]).matrix_mult((*sol->_Sol[soluIndex[j]]), (*P[k]));
    }
  }

}


