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
void BuildFlag(MultiLevelSolution& mlSol);
double eps = 0.;

double gravity[3] = {0., 0., 0.};
bool weakP = false;

double theta = .75;
double af = 1. - theta;
double am = af - 0.25;
double beta = 0.25 + 0.5 * (af - am);
double Gamma = 0.5 + (af - am);

double gammac = 0.05;
double gammap = 0.01;
double gammau = 0.01 * gammac;

double GAMMA = 10;   // 10, 45 in the paper.

#include "../ex12/include/mpmFsi10.hpp"
using namespace femus;

double SetVariableTimeStep(const double time) {
  double dt = 1.;
  if(time < 2.) dt = 0.05;
  else dt = 0.01;

  return dt;
}



bool SetBoundaryCondition(const std::vector < double >&x, const char name[], double &value, const int facename, const double t) {
  bool test = 1;      //dirichlet
  value = 0.;

  const double Ubar = 1.0;    // guess?
  const double L = 0.41;
  const double H = 2.5;

  if(!strcmp(name, "DY")) {
    if(1 == facename || 2 == facename) {
      test = 0;
      value = 0.;
    }
  }

  else if(!strcmp(name, "DX")) {
    if(3 == facename) {    //fluid wall
      test = 0;
      value = 0.;
    }
  }

  else if(!strcmp(name, "VY")) {
    if(2 == facename) {     //outflow
      test = 0;
      value = 0.;
    }
  }

  else if(!strcmp(name, "VX")) {
    if(1 == facename) {     //inflow
      test = 1;
      if(t < 2.0) {
        value = 1.5 * Ubar * 4.0 / 0.1681 * (x[1] + 0.21) * (-x[1] + 0.2) * 0.5 * (1. - cos(0.5 * M_PI * t));
      }
      else {
        value = 1.5 * Ubar * 4.0 / 0.1681 * (x[1] + 0.21) * (-x[1] + 0.2);
      }
    }
    else if(2 == facename) {    //outflow
      test = 0;
      value = 0.;
    }
  }

  else if(!strcmp(name, "P")) {
    if(weakP || 2 != facename) {
      test = 0;
    }
    value = 0;
  }

  return test;

}

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  MultiLevelMesh mlMsh;

  double scalingFactor = 1.;
  unsigned numberOfUniformLevels = 5; //for refinement in 3D
  unsigned numberOfUniformLevelsStart = numberOfUniformLevels;
  //unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;


  double Lref = 1.;
  double Uref = 1.;
  double rhof = 1000.;
  double muf = 1.;
  double rhos = 10000.;
  double nu = 0.4;
  double E = 1400000;


  Parameter par(Lref, Uref);

  // Generate Solid Object
  Solid solid(par, E, nu, rhos, "Neo-Hookean");
  Fluid fluid(par, muf, rhof, "Newtonian");

  mlMsh.ReadCoarseMesh("../input/turek2D.neu", "fifth", scalingFactor);
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

  mlSol.SetIfFSI(true);

  mlSol.Initialize("All");

  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

  // ******* Set boundary conditions *******
  mlSol.GenerateBdc("DX", "Steady");
  if(dim > 1) mlSol.GenerateBdc("DY", "Steady");
  if(dim > 2) mlSol.GenerateBdc("DZ", "Steady");
  mlSol.GenerateBdc("VX", "Time_dependent");
  if(dim > 1) mlSol.GenerateBdc("VY", "Steady");
  if(dim > 2) mlSol.GenerateBdc("VZ", "Steady");
  mlSol.GenerateBdc("P", "Steady");

  MultiLevelProblem ml_prob(&mlSol);

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



  // ******* Add MPM system to the MultiLevel problem *******
  NonLinearImplicitSystem& system2 = ml_prob.add_system < NonLinearImplicitSystem > ("DISP");
  system2.AddSolutionToSystemPDE("DX");
  if(dim > 1) system2.AddSolutionToSystemPDE("DY");
  if(dim > 2) system2.AddSolutionToSystemPDE("DZ");

  // ******* System MPM Assembly *******
//     system2.SetAssembleFunction(AssembleSolidDisp);
  //system2.SetAssembleFunction(AssembleFEM);
  // ******* set MG-Solver *******
  system2.SetMgType(V_CYCLE);


  system2.SetAbsoluteLinearConvergenceTolerance(1.0e-10);
  system2.SetMaxNumberOfLinearIterations(1);
  system2.SetNonLinearConvergenceTolerance(1.e-9);
  system2.SetMaxNumberOfNonLinearIterations(1);

  system2.SetNumberPreSmoothingStep(1);
  system2.SetNumberPostSmoothingStep(1);

  // ******* Set Preconditioner *******
  system2.SetLinearEquationSolverType(FEMuS_DEFAULT);

  system2.init();

  // ******* Set Smoother *******
  system2.SetSolverFineGrids(GMRES);

  system2.SetPreconditionerFineGrids(ILU_PRECOND);

  system2.SetTolerances(1.e-10, 1.e-15, 1.e+50, 2, 2);


  std::ifstream fin;
  std::ostringstream level_number;
  level_number << 2;

  //BEGIN bulk reading

  std::string bulkfile = "../input/";
  bulkfile += "turekBeam2D";
  bulkfile += level_number.str();
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


  double delta_max = 0.013 / (numberOfUniformLevelsStart - 4);

  for(int i = 0; i < xp.size(); i++) {
    if(dist[i] < -delta_max) {
      xp.erase(xp.begin() + i);
      wp.erase(wp.begin() + i);
      dist.erase(dist.begin() + i);
      markerType.erase(markerType.begin() + i);
      i--;
    }
  }


  double shift = 1.0e-4;
  for(int i = 0; i < xp.size(); i++) {
    xp[i][0] += shift;
  }


  bulk = new Line(xp, wp, dist, markerType, mlSol.GetLevel(numberOfUniformLevels - 1), 2);

  std::vector < std::vector < std::vector < double > > >  bulkPoints(1);
  bulk->GetLine(bulkPoints[0]);
  PrintLine(DEFAULT_OUTPUTDIR, "bulk", bulkPoints, 0);

  //END bulk reading

  //BEGIN interface reading



  std::string interfacefile = "../input/";
  interfacefile += "turekBeam2D";
  interfacefile += level_number.str();
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

  double xmax = -1.0e10;
  double ymax =  1.0e10;
  unsigned imax = lineIPoints[0].size();
  for(unsigned i = 0; i < lineIPoints[0].size(); i++) {
    if(lineIPoints[0][i][0] > xmax) {
      imax = i;
      xmax = lineIPoints[0][i][0];
      ymax = fabs(lineIPoints[0][i][1]);
    }
    else if(lineIPoints[0][i][0] = xmax) {
      if(fabs(lineIPoints[0][i][1]) < ymax) {
        imax = i;
        xmax = lineIPoints[0][i][0];
        ymax = fabs(lineIPoints[0][i][1]);
      }
    }
  }

  std::cout << "imax = " << imax << " xmax = " << xmax << " ymax = " << ymax << std::endl;
  PrintLine(DEFAULT_OUTPUTDIR, "interfaceMarkers", lineIPoints, 0);
//END interface reading


  unsigned iproc = mlMsh.GetLevel(0)->processor_id();

  std::ofstream fout;
  std::ofstream pout;

  level_number << numberOfUniformLevelsStart;

  std::string ofile = "./save/";
  ofile += "beamTipPositionLevel";
  ofile += level_number.str();
  ofile += ".turekFSI2.txt";

  std::string pfile = "./save/";
  pfile += "pressureLevel";
  pfile += level_number.str();
  pfile += ".turekFSI2.txt";

  if(iproc == 0) {
    fout.open(ofile);
    fout << 0 << " " << lineIPoints[0][imax][0] << " " << lineIPoints[0][imax][1] << std::endl;

    pout.open(pfile);
    pout.close();
  }

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
  unsigned n_timesteps = 3000;
  for(unsigned time_step = 1; time_step <= n_timesteps; time_step++) {

    system.CopySolutionToOldSolution();

    mlSol.GetWriter()->Write("output1", "biquadratic", print_vars, time_step);
    bulk->GetLine(bulkPoints[0]);
    PrintLine("output1", "bulk", bulkPoints, time_step);
    lineI->GetLine(lineIPoints[0]);
    PrintLine("output1", "interfaceMarkers", lineIPoints, time_step);

    system.MGsolve();

    double time = system.GetTime();

    GetPressureDragAndLift(ml_prob, time, imax, pfile);

    mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, time_step);
    GridToParticlesProjection(ml_prob, *bulk, *lineI);
    bulk->GetLine(bulkPoints[0]);
    PrintLine(DEFAULT_OUTPUTDIR, "bulk", bulkPoints, time_step);
    lineI->GetLine(lineIPoints[0]);
    PrintLine(DEFAULT_OUTPUTDIR, "interfaceMarkers", lineIPoints, time_step);

    if(iproc == 0) fout << time << " " << lineIPoints[0][imax][0] << " " << lineIPoints[0][imax][1] << std::endl;
  }

  if(iproc == 0) fout.close();

  delete bulk;
  delete lineI;

  return 0;

} //end main

void BuildFlag(MultiLevelSolution& mlSol) {

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


