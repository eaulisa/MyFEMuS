#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "MultiLevelSolution.hpp"
#include "MeshRefinement.hpp"

#include "TransientSystem.hpp"
#include "NonLinearImplicitSystem.hpp"

#include "Fluid.hpp"
#include "Solid.hpp"
#include "Parameter.hpp"

#include "PetscMatrix.hpp"

bool weakP = false;
double theta = 1.;
double af = 1. - theta;
double am = af;
double Beta = 0.25 + 0.5 * (af - am);
double Gamma = 0.5 + (af - am);

using namespace femus;

#include "marker.hpp"
#include "background.hpp"
#include "projection.hpp"

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  MultiLevelMesh mlMshM;
  mlMshM.ReadCoarseMesh("../input/benchmarkFSISolid.neu", "fifth", 1);
  mlMshM.RefineMesh(3, 3, NULL); //uniform refinement, this goes with the background mesh refinement. For COMSOL we use 8 = 3 turekBeam2D

  unsigned numberOfRefinementM = 3;
  for(unsigned i = 0; i < numberOfRefinementM; i++) {
    FlagElements(mlMshM, 2);
    mlMshM.AddAMRMeshLevel();
  }
  mlMshM.EraseCoarseLevels(numberOfRefinementM);

  MultiLevelSolution mlSolM(&mlMshM);
  InitializeMarkerVariables(mlSolM);

  UpdateMeshQuantities(&mlSolM);
  
  unsigned numberOfRefinementB = numberOfRefinementM;
  MultiLevelMesh mlMshB;
  //mlMshB.ReadCoarseMesh("../input/benchmarkFSI.neu", "fifth", 1);
  mlMshB.ReadCoarseMesh("../input/turek2D.neu", "fifth", 1); // FSI_mesh2
  mlMshB.RefineMesh(numberOfRefinementB, numberOfRefinementB, NULL);
  mlMshB.EraseCoarseLevels(numberOfRefinementB - 1u);

  MultiLevelSolution mlSolB(&mlMshB);
  InitializeBackgroundVariables(mlSolB);


  unsigned dim = mlMshB.GetDimension();
  // ******* Set boundary conditions *******
  mlSolB.GenerateBdc("DX", "Steady");
  if(dim > 1) mlSolB.GenerateBdc("DY", "Steady");
  if(dim > 2) mlSolB.GenerateBdc("DZ", "Steady");
  mlSolB.GenerateBdc("VX", "Time_dependent");
  if(dim > 1) mlSolB.GenerateBdc("VY", "Steady");
  if(dim > 2) mlSolB.GenerateBdc("VZ", "Steady");
  mlSolB.GenerateBdc("P", "Steady");


  double Lref = 1.;
  double Uref = 1.;
  double rhof = 1000.;
  double muf = 1.;
  double rhos = 1000.;
  double nu = 0.4;
  double E = 1400000; //FSI1
  //double E = 4 * 1400000; //FSI3
  
  Parameter par(Lref, Uref);
  // Generate Solid Object
  Solid solid(par, E, nu, rhos, "Neo-Hookean");
  Fluid fluid(par, muf, rhof, "Newtonian");
  MultiLevelProblem ml_prob(&mlSolB);
  ml_prob.parameters.set<Solid> ("SolidMPM") = solid;
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

  // ******* System MPM-FSI Assembly *******
  system.SetAssembleFunction(AssembleMPMSys);
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

  system.AttachGetTimeIntervalFunction(SetVariableTimeStepB);

  //******* Print solution *******
  mlSolM.SetWriter(VTK);
  mlSolB.SetWriter(VTK);
  
  std::vector<std::string> mov_vars;
  mov_vars.push_back("DX");
  mov_vars.push_back("DY");
  if(dim == 3) mov_vars.push_back("DZ");
  mlSolM.GetWriter()->SetMovingMesh(mov_vars);
  //mlSolB.GetWriter()->SetMovingMesh(mov_vars);
  mlSolM.GetWriter()->SetDebugOutput(false);
  mlSolB.GetWriter()->SetDebugOutput(true);
  std::vector<std::string> print_vars;
  print_vars.push_back("All");

  mlSolM.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, 0);
  mlSolB.GetWriter()->Write(DEFAULT_OUTPUTDIR, "linear", print_vars, 0);

  Projection projection(&mlSolM, &mlSolB);
 
  for(unsigned t = 1; t <= 20; t++) {
    projection.SetNewmarkParameters(Beta, Gamma, 1.);  
    clock_t time = clock();
    //projection.FromMarkerToBackground();

    system.MGsolve();
    
    //projection.FromBackgroundToMarker();
    std::cout << "time" << t << " = " << static_cast<double>((clock() - time)) / CLOCKS_PER_SEC << std::endl;
    mlSolM.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, t);
    mlSolB.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, t);
  }

  return 0;

} //end main

double SetVariableTimeStepB(const double time) {
  double dt; 
  double dt0 = 0.1;
  double dt1 = 0.1; //FSI3
  //double dt1 = 0.15;  // Steady state-FSI1

  double T = 6.;
  if(time < T)
    dt = (0.5 * (1. + cos(M_PI * time / T))) * dt0    + (0.5 * (1. - cos(M_PI * time / T))) * dt1;
    //dt = dt0; // steady state-FSI1
  else
    dt = dt1;

  return dt;
}


bool SetBoundaryConditionB(const std::vector < double >&x, const char name[], double &value, const int facename, const double t) {
  bool test = 1;      //dirichlet
  value = 0.;

  const double Ubar = 0.2;    // FSI1
  //const double Ubar = 2;    // FSI3
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
  }

  return test;

}
