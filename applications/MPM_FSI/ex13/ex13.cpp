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
void BuildIntegrationPoints(MultiLevelSolution& mlSol);
void GetDragAndLift(MultiLevelProblem& ml_prob, const double & time, const std::string &pfile);

bool NeoHookean = false;
bool particleSmoothingIsOn = false;

double eps = 0.00;

double gravity[3] = {0., 0., 0.};
bool weakP = false;

double theta = 1.;
double af = 1. - theta;
double am = af - 0.;
double beta = 0.25 + 0.5 * (af - am);
double Gamma = 0.5 + (af - am);
// double beta = .5;
// double Gamma = 1.;

double DTMIN = 0.001;

double factor = 2.;

double gammacF = factor * 0.05;
double gammacS = factor * 0.05;
double gammap = factor * 0.05;
double gammau = 0.05 * gammacF;

double GAMMA = 45;//45;   // 10, 45 in the paper.

#include "../ex12/include/mpmFsi10.hpp"
using namespace femus;

double SetVariableTimeStep(const double time) {
  //double dt = 1.; //FSI1
  double dt = 0.001; //FSI3

  if(time < 2.) dt = .05;
  else dt = 0.001;

  return dt;
}



bool SetBoundaryCondition(const std::vector < double >&x, const char name[], double &value, const int facename, const double t) {
  bool test = 1;      //dirichlet
  value = 0.;

  //const double Ubar = 0.2;    // FSI1
  const double Ubar = 2;    // FSI3
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
  double rhos = 1000.;
  double nu = 0.4;
  //double E = 1400000; //FSI1
  double E = 4 * 1400000; //FSI3


  Parameter par(Lref, Uref);

  // Generate Solid Object
  Solid solid(par, E, nu, rhos, "Neo-Hookean");
  Fluid fluid(par, muf, rhof, "Newtonian");

  mlMsh.ReadCoarseMesh("../input/turek2D.neu", "fifth", scalingFactor); // FSI_mesh1
  //mlMsh.ReadCoarseMesh("../input/turek2DNew.neu", "fifth", scalingFactor); // FSI_mesh2
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

  BuildIntegrationPoints(mlSol);

  std::ifstream fin;
  std::ostringstream fileName;
  std::ostringstream level_number;
  fileName << "../input/turekBeam2D";
  level_number << 2;
  //  fileName<<"../input/turekBeam2DNew"; level_number << 4;

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


  double delta_max = 0.013 / (numberOfUniformLevelsStart - 4);
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


  double xmax = -1.0e10;
  double ymax =  1.0e10;
  unsigned imax = 0;

  for(unsigned i = 0; i < lineIPoints[0].size(); i++) {
    if(lineIPoints[0][i][0] > xmax) {
      imax = i;
      xmax = lineIPoints[0][i][0];
      ymax = fabs(lineIPoints[0][i][1]);
    }
    else if(lineIPoints[0][i][0] == xmax) {
      if(fabs(lineIPoints[0][i][1]) < ymax) {
        imax = i;
        xmax = lineIPoints[0][i][0];
        ymax = fabs(lineIPoints[0][i][1]);
      }
    }
  }
  std::cout << "imax = " << imax << " xmax = " << xmax << " ymax = " << ymax << std::endl;

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

  intPoints->GetParticlesToGridMaterial(false);
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
  unsigned n_timesteps = 20000;
  unsigned printTimeInterval = 50;
  for(unsigned time_step = 1; time_step <= n_timesteps; time_step++) {

    system.CopySolutionToOldSolution();

    mlSol.GetWriter()->Write("output1", "biquadratic", print_vars, time_step / printTimeInterval);
    bulk->GetLine(bulkPoints[0]);
    PrintLine("output1", "bulk", bulkPoints, time_step);
    lineI->GetLine(lineIPoints[0]);
    PrintLine("output1", "interfaceMarkers", lineIPoints, time_step / printTimeInterval);

    system.MGsolve();

    double time = system.GetTime();

    GetDragAndLift(ml_prob, time, pfile);

    mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, time_step / printTimeInterval);
    GridToParticlesProjection(ml_prob, *bulk, *lineI);
    bulk->GetLine(bulkPoints[0]);
    PrintLine(DEFAULT_OUTPUTDIR, "bulk", bulkPoints, time_step / printTimeInterval);
    lineI->GetLine(lineIPoints[0]);
    PrintLine(DEFAULT_OUTPUTDIR, "interfaceMarkers", lineIPoints, time_step / printTimeInterval);

    if(iproc == 0) fout << time << " " << lineIPoints[0][imax][0] << " " << lineIPoints[0][imax][1] << std::endl;
  }

  if(iproc == 0) fout.close();

  delete bulk;
  delete lineI;
  delete intPoints;

  return 0;

} //end main




void BuildIntegrationPoints(MultiLevelSolution& mlSol) {
  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution *sol  = mlSol.GetSolutionLevel(level);
  Mesh     *msh   = mlSol._mlMesh->GetLevel(level);
  elem* el = msh->el;   // pointer to the elem object in msh (level)
  unsigned iproc  = msh->processor_id();

  const unsigned dim = msh->GetDimension();
  std::vector <std::vector < double> > vx(dim);

  std::vector < double> normal;
  std::vector < double> phi;
  std::vector < double> phi_x;
  double weight;

  std::vector < std::vector <double> > x(dim);
  std::vector < std::vector <double> > N(2);

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned nDofs = msh->GetElementDofNumber(iel, 2);
    for(int k = 0; k < dim; k++) {
      vx[k].resize(nDofs);
    }
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned idofX = msh->GetSolutionDof(i, iel, 2);
      for(unsigned  k = 0; k < dim; k++) {
        vx[k][i] = (*msh->_topology->_Sol[k])(idofX);
      }
    }
    for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
      if(el->GetBoundaryIndex(iel, jface) == 4) {

        const unsigned typef = msh->GetElementFaceType(iel, jface);
        unsigned nDofsf = msh->GetElementFaceDofNumber(iel, jface, 2);

        std::vector  < std::vector  <  double> > vxf(dim);    // A matrix holding the face coordinates rowwise.
        for(int k = 0; k < dim; k++) {
          vxf[k].resize(nDofsf);
        }

        bool solid = true;
        for(unsigned i = 0; i < nDofsf; i++) {
          unsigned inode = msh->GetLocalFaceVertexIndex(iel, jface, i);    // face-to-element local node mapping.
          for(unsigned k = 0; k < dim; k++) {
            vxf[k][i] =  vx[k][inode]; // We extract the local coordinates on the face from local coordinates on the element.
          }
          if(vxf[0][i] < 0. || fabs(vxf[1][i]) > 0.01) solid = false;

        }
        if(solid) {
          for(unsigned ig = 0; ig  <  msh->_finiteElement[typef][2]->GetGaussPointNumber(); ig++) {
            msh->_finiteElement[typef][2]->JacobianSur(vxf, ig, weight, phi, phi_x, normal);
            std::vector<double> xg(dim, 0);

            unsigned n = x[0].size();

            for(unsigned k = 0; k < dim; k++) {
              x[k].resize(n + 1);
            }
            for(unsigned i = 0; i < nDofsf; i++) {
              for(unsigned k = 0; k < dim; k++) {
                x[k][n] += phi[i] * vxf[k][i];
              }
            }

            N[0].resize(n + 1);
            N[1].resize(n + 1);
            N[0][n] = weight * normal[0];
            N[1][n] = weight * normal[1];

          }
        }
      }
    }
  }

  unsigned sizeLocal = x[0].size();
  unsigned sizeAll;

  MPI_Allreduce(&sizeLocal, &sizeAll, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);

  std::vector < std::vector <double> > xr(dim);
  std::vector < std::vector <double> > Nr(2);

  xr[0].resize(sizeAll);
  xr[1].resize(sizeAll);
  Nr[0].resize(sizeAll);
  Nr[1].resize(sizeAll);

  unsigned nprocs = msh->n_processors();
  std::vector <int> sizeLocals(nprocs);

  MPI_Allgather(&sizeLocal, 1, MPI_UNSIGNED, sizeLocals.data(), 1, MPI_UNSIGNED, PETSC_COMM_WORLD);

  std::vector <int> relPos(nprocs);
  relPos[0] = 0.;
  for(unsigned i = 1; i < relPos.size(); i++) {
    relPos[i] = relPos[i - 1] + sizeLocals[i - 1];
  }
    
  MPI_Allgatherv(x[0].data(), x[0].size(), MPI_DOUBLE, xr[0].data(), sizeLocals.data(), relPos.data(), MPI_DOUBLE, PETSC_COMM_WORLD);
  MPI_Allgatherv(x[1].data(), x[1].size(), MPI_DOUBLE, xr[1].data(), sizeLocals.data(), relPos.data(), MPI_DOUBLE, PETSC_COMM_WORLD);
  MPI_Allgatherv(N[0].data(), N[0].size(), MPI_DOUBLE, Nr[0].data(), sizeLocals.data(), relPos.data(), MPI_DOUBLE, PETSC_COMM_WORLD);
  MPI_Allgatherv(N[1].data(), N[1].size(), MPI_DOUBLE, Nr[1].data(), sizeLocals.data(), relPos.data(), MPI_DOUBLE, PETSC_COMM_WORLD);

  std::vector < std::vector <double> > xp(sizeAll);
  std::vector < std::vector < std::vector <double> > > Tp(sizeAll);
  std::vector < MarkerType > markerType(sizeAll, INTERFACE);

  for(unsigned i = 0; i < sizeAll; i++) {
    xp[i].resize(dim);
    Tp[i].resize(1);
    Tp[i][0].resize(dim);
    for(unsigned k = 0; k < dim; k++) {
      xp[i][k] = xr[k][i];
    }

    Tp[i][0][0] = -Nr[1][i];
    Tp[i][0][1] = Nr[0][i];
  }

  intPoints = new Line(xp, Tp, markerType, mlSol.GetLevel(level), 2);

  std::vector < std::vector < std::vector < double > > > integrationPoints(1);
  intPoints->GetLine(integrationPoints[0]);

  PrintLine(DEFAULT_OUTPUTDIR, "integrationPoints", integrationPoints, 0);

}



void GetDragAndLift(MultiLevelProblem& ml_prob, const double & time, const std::string &pfile) {

  //pointers and references

  TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("MPM_FSI");
  const unsigned  level = my_nnlin_impl_sys.GetLevelToAssemble();
  MultiLevelSolution* mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* mysolution = mlSol->GetSolutionLevel(level);     // pointer to the solution (level) object

  Mesh* msh = ml_prob._ml_msh->GetLevel(level);     // pointer to the mesh (level) object
  elem* el = msh->el;   // pointer to the elem object in msh (level)

  const unsigned dim = msh->GetDimension();

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  // data
  unsigned iproc  = msh->processor_id();

  std::ofstream pout;
  if(iproc == 0) {
    pout.open(pfile,  std::ios_base::app);
    pout << time;
    pout.close();
  }


  vector< vector< double > > solDTld(dim);      // as in the paper
  vector< vector< double > > solV(dim);         // velocity at n+1
  vector< double > solP;

  vector< vector< double > > solDOld(dim);      // displacement at n

  vector < double > phi;
  vector < double > phiP;
  double weight;

  double lenght = 0.;

  vector <vector < double> > vx(dim); // background mesh configuration at n + 1
  vector <vector < double> > vxOld(dim); // background mesh configuration at n

  vector < double> gradPhi;  // gradient with respect to vx
  vector < double > gradOldPhi; // gradient with respect to vxOld

  double dragF = 0.;
  double liftF = 0.;

  double dragS = 0.;
  double liftS = 0.;

  //reading parameters for fluid FEM domain

  double muFluid = ml_prob.parameters.get<Fluid> ("FluidFEM").get_viscosity();
  double muMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_lame_shear_modulus();
  double lambdaMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_lame_lambda();

  std::cout.precision(10);

  //variable-name handling
  const char varname[10][5] = {"DX", "DY", "DZ", "VX", "VY", "VZ"};

  vector <unsigned> indexSolD(dim);
  vector <unsigned> indexSolV(dim);
  for(unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolD[ivar] = mlSol->GetIndex(&varname[ivar][0]);
    indexSolV[ivar] = mlSol->GetIndex(&varname[ivar + 3][0]);
  }
  unsigned solType = mlSol->GetSolutionType(&varname[0][0]);

  unsigned indexSolP = mlSol->GetIndex("P");
  unsigned solTypeP = mlSol->GetSolutionType("P");

  unsigned eflagIndex = mlSol->GetIndex("eflag");

  std::vector < std::vector < std::vector <double > > > aP(3);

  std::vector<Marker*> particleI = intPoints->GetParticles();
  std::vector<unsigned> markerOffsetI = intPoints->GetMarkerOffset();
  unsigned imarkerI = markerOffsetI[iproc];

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielt = msh->GetElementType(iel);

    unsigned nDofs = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solTypeP);  // number of pressure dofs

    unsigned eFlag = static_cast <unsigned>(floor((*mysolution->_Sol[eflagIndex])(iel) + 0.25));

    for(unsigned  k = 0; k < dim; k++) {
      solDTld[k].resize(nDofs);
      solDOld[k].resize(nDofs);

      solV[k].resize(nDofs);

      vx[k].resize(nDofs);
      vxOld[k].resize(nDofs);
    }
    solP.resize(nDofsP);


    for(unsigned i = 0; i < nDofs; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, solType);

      for(unsigned  k = 0; k < dim; k++) {
        solDOld[k][i] = (*mysolution->_SolOld[indexSolD[k]])(idof);//t_n
        solDTld[k][i] = (*mysolution->_Sol[indexSolD[k]])(idof) - solDOld[k][i]; //t_{n+1} -t_n
        solV[k][i] = (*mysolution->_Sol[indexSolV[k]])(idof);
      }
    }

    for(unsigned i = 0; i < nDofsP; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, solTypeP);
      solP[i] = (*mysolution->_Sol[indexSolP])(idof);
    }

    for(unsigned i = 0; i < nDofs; i++) {
      unsigned idofX = msh->GetSolutionDof(i, iel, 2);
      for(unsigned  k = 0; k < dim; k++) {
        vxOld[k][i] = (*msh->_topology->_Sol[k])(idofX) + solDOld[k][i]; // deformed reference configuration
        vx[k][i]  = vxOld[k][i] + solDTld[k][i]; // also equal to vxHat + SolD
      }
    }

    for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
      if(el->GetBoundaryIndex(iel, jface) == 4) {

        const unsigned typef = msh->GetElementFaceType(iel, jface);
        unsigned nDofsf = msh->GetElementFaceDofNumber(iel, jface, solType);

        std::vector  < std::vector  <  double> > vxf(dim);    // A matrix holding the face coordinates rowwise.
        for(int k = 0; k < dim; k++) {
          vxf[k].resize(nDofsf);
        }

        bool solid = true;
        for(unsigned i = 0; i < nDofsf; i++) {
          unsigned inode = msh->GetLocalFaceVertexIndex(iel, jface, i);    // face-to-element local node mapping.
          for(unsigned k = 0; k < dim; k++) {
            vxf[k][i] =  vx[k][inode]; // We extract the local coordinates on the face from local coordinates on the element.
          }
          if(vxf[0][i] < 0. || fabs(vxf[1][i]) > 0.01) solid = false;

        }
        if(!solid) {

          for(unsigned jtype = 0; jtype < solType + 1; jtype++) {
            ProjectNodalToPolynomialCoefficients(aP[jtype], vx, ielt, jtype);
          }

          for(unsigned ig = 0; ig  <  msh->_finiteElement[typef][solType]->GetGaussPointNumber(); ig++) {

            std::vector<double> normal(dim);
            double area;
            msh->_finiteElement[typef][solType]->JacobianSur(vxf, ig, area, phi, gradPhi, normal);
            std::vector<double> xg(dim, 0);

            for(unsigned i = 0; i < nDofsf; i++) {
              for(unsigned k = 0; k < dim; k++) {
                xg[k] += phi[i] * vxf[k][i];
              }
            }

            std::vector <double> xi;//local coordinates of the face gauss point with respect to iel
            GetClosestPointInReferenceElement(vx, xg, ielt, xi);
            bool inverseMapping = GetInverseMapping(solType, ielt, aP, xg, xi, 100);

            msh->_finiteElement[ielt][solType]->Jacobian(vx, xi, weight, phi, gradPhi);

            std::vector < double > tauF(dim, 0.);
            for(unsigned k = 0; k < dim; k++) {
              for(unsigned i = 0; i < nDofs; i++) {
                for(unsigned j = 0; j < dim; j++) {
                  tauF[k] += muFluid * (solV[k][i] * gradPhi[i * dim + j] + solV[j][i] * gradPhi[i * dim + k]) * normal[j];
                }
              }
            }

            msh->_finiteElement[ielt][solTypeP]->Jacobian(vx, xi, weight, phi, gradPhi);
            for(unsigned k = 0; k < dim; k++) {
              for(unsigned i = 0; i < nDofsP; i++) {
                tauF[k] -= solP[i] * phi[i] * normal[k];
              }
            }

            lenght += area;

            if(dim == 2) {
              dragF += tauF[0] * area;
              liftF += tauF[1] * area;
            }
            else if(dim == 3) {
              liftF += (tauF[0] * normal[0] + tauF[1] * normal[1] + tauF[2] * normal[2]) * area;
            }
          }
        }
      }
    }


    //BEGIN INTERFACE PARTICLE

    if(eFlag >= 1) {  //interface markers

      while(imarkerI < markerOffsetI[iproc + 1] && iel > particleI[imarkerI]->GetMarkerElement()) {
        imarkerI++;
      }

      while(imarkerI < markerOffsetI[iproc + 1] && iel == particleI[imarkerI]->GetMarkerElement()) {

        //std::cout << imarkerI << " " << particleI[imarkerI]->GetMarkerElement() << " " << iel << std::endl;

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particleI[imarkerI]->GetMarkerLocalCoordinates();
        msh->_finiteElement[ielt][solType]->Jacobian(vxOld, xi, weight, phi, gradOldPhi);
        msh->_finiteElement[ielt][solType]->Jacobian(vx, xi, weight, phi, gradPhi);

        vector<vector < double > > gradOldSolDTld(dim); // $\widetilde{\nabla} \widetilde{\bm D}$
        for(int j = 0; j < dim; j++) {
          gradOldSolDTld[j].assign(dim, 0.);
        }
        for(int j = 0; j < dim; j++) {
          for(unsigned i = 0; i < nDofs; i++) {
            for(int k = 0; k < dim; k++) {
              gradOldSolDTld[k][j] += solDTld[k][i] * gradOldPhi[i * dim + j];
            }
          }
        }

        std::vector <std::vector < double > > T;
        particleI[imarkerI]->GetMarkerTangent(T);

        std::vector < std::vector < double > > FHatOld;
        FHatOld = particleI[imarkerI]->GetDeformationGradient(); //extraction of the deformation gradient

        double FTld[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
        double FHat[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
        std::vector<std::vector <double>> svFHat(dim);

        double Cauchy[3][3];
        double Id2th[3][3] = {{ 1., 0., 0.}, { 0., 1., 0.}, { 0., 0., 1.}};

        for(unsigned j = 0; j < dim; j++) {
          for(unsigned k = 0; k < dim; k++) {
            FTld[j][k] += gradOldSolDTld[j][k];
          }
        }

        for(unsigned i = 0; i < dim; i++) {
          svFHat[i].resize(dim);
          for(unsigned j = 0; j < dim; j++) {
            svFHat[i][j] = 0.;
            for(unsigned k = 0; k < dim; k++) {
              FHat[i][j] += FTld[i][k] * FHatOld[k][j];
              svFHat[i][j] += FTld[i][k] * FHatOld[k][j];
            }
          }
        }

        particleI[imarkerI]->SetDeformationGradient(svFHat);


        if(dim == 2) FHat[2][2] = 1.;
        double JHat =  FHat[0][0] * FHat[1][1] * FHat[2][2] + FHat[0][1] * FHat[1][2] * FHat[2][0] + FHat[0][2] * FHat[1][0] * FHat[2][1]
                       - FHat[2][0] * FHat[1][1] * FHat[0][2] - FHat[2][1] * FHat[1][2] * FHat[0][0] - FHat[2][2] * FHat[1][0] * FHat[0][1];


        //std::cout << JHat << " ";

        if(NeoHookean) {
          double B[3][3];
          for(unsigned i = 0; i < 3; i++) {
            for(int j = 0; j < 3; j++) {
              B[i][j] = 0.;
              for(unsigned k = 0; k < 3; k++) {
                //left Cauchy-Green deformation tensor or Finger tensor (B = F*F^T)
                B[i][j] += FHat[i][k] * FHat[j][k];
              }
            }
          }

          for(unsigned j = 0; j < 3; j++) {
            for(unsigned k = 0; k < 3; k++) {
              Cauchy[j][k] = lambdaMpm * log(JHat) / JHat * Id2th[j][k] + muMpm / JHat * (B[j][k] - Id2th[j][k]);     //alternative formulation
            }
          }
        }
        else {
          double E[3][3];
          double S[3][3];

          for(unsigned i = 0; i < 3; i++) { //E = 0.5(F^T F - I)
            for(unsigned j = 0; j < 3; j++) {
              E[i][j] = 0.;
              for(unsigned k = 0; k < 3; k++) {
                E[i][j] += FHat[k][i] * FHat[k][j];
              }
              E[i][j] = 0.5 * (E[i][j] - Id2th[i][j]);
            }
          }

          double traceE = E[0][0] + E[1][1] + E[2][2];

          for(unsigned i = 0; i < 3; i++) { // S = lambda Tr(E) +  2 mu E
            for(unsigned j = 0; j < 3; j++) {
              S[i][j] = lambdaMpm * traceE * Id2th[i][j] + 2. * muMpm * E[i][j];     //alternative formulation
            }
          }

          double SFt[3][3];
          for(unsigned i = 0; i < 3; i++) { // S F^t
            for(unsigned j = 0; j < 3; j++) {
              SFt[i][j] = 0.;
              for(unsigned k = 0; k < 3; k++) {
                SFt[i][j] += S[i][k] * FHat[j][k];
              }
            }
          }

          for(unsigned i = 0; i < 3; i++) { // 1./J F S F^t
            for(unsigned j = 0; j < 3; j++) {
              Cauchy[i][j] = 0.;
              for(unsigned k = 0; k < 3; k++) {
                Cauchy[i][j] += FHat[i][k] * SFt[k][j] / JHat;
              }
            }
          }
        }


        std::vector < double > N(dim); // normal pointing toward the outside
        if(dim == 2) {
          N[0] =  T[0][1];
          N[1] = -T[0][0];
          weight = sqrt(N[0] * N[0] + N[1] * N[1]);
          N[0] /= weight;
          N[1] /= weight;
        }
        else {
          N[0] = T[0][1] * T[1][2] - T[0][2] * T[1][1];
          N[1] = T[0][2] * T[1][0] - T[0][0] * T[1][2];
          N[2] = T[0][0] * T[1][1] - T[0][1] * T[1][0];
          weight = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
          N[0] /= weight;
          N[1] /= weight;
          N[2] /= weight;
        }

        lenght += weight;

        std::vector < double > tauS(dim, 0.);
        for(unsigned i = 0; i < dim; i++) {
          for(unsigned j = 0; j < dim; j++) {
            tauS[i] += Cauchy[i][j] * N[j];
          }
        }

        if(dim == 2) {
          dragS += tauS[0] * weight;
          liftS += tauS[1] * weight;
        }
        else if(dim == 3) {
          liftS += (tauS[0] * N[0] + tauS[1] * N[1] + tauS[2] * N[2]) * weight;
        }
        imarkerI++;
      }
    }
  }

  std::cout << "arclength = " << lenght << " " << 2.*M_PI * 0.05 << std::endl;

  double dragFAll, liftFAll;
  double dragSAll, liftSAll;
  MPI_Reduce(&dragF, &dragFAll, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  MPI_Reduce(&liftF, &liftFAll, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);

  MPI_Reduce(&dragS, &dragSAll, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  MPI_Reduce(&liftS, &liftSAll, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);

  if(iproc == 0) {
    pout.open(pfile,  std::ios_base::app);
    pout << " " << -dragFAll << " " << -liftFAll << " " << -dragSAll << " " << -liftSAll
         << " " << -dragFAll - dragSAll << " " << -liftFAll - liftSAll << std::endl;
    pout.close();
  }
}








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


