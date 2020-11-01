
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

double beta = 0.25;
double Gamma = 0.5;
double theta = 0.5;
double af = theta;
double pInf = (1. + af) / (2. - af);
double am = pInf / (1. + pInf);

const elem_type *fem[2][6][5];

//#include "../../Nitsche/support/particleInit.hpp"
//#include "../../Nitsche/support/sharedFunctions.hpp"
#include "./include/assemblySolid.hpp"
#include "./include/assemblyBoundaryLayer.hpp"
#include "./include/assemblyFluid.hpp"
using namespace femus;

double dt = 1.;

double SetVariableTimeStep(const double time) {
  if(time < 2.5) dt = 0.05;
  else dt = 0.0125;

  return dt;
}

void Assemble(MultiLevelProblem& ml_prob);

void InitPElement(MultiLevelSolution & mlSol, const std::vector <double> &x0);
void BuildFlagSolidRegion(MultiLevelSolution & mlSol);
void ProjectNewmarkDisplacemenet(MultiLevelSolution & mlSol);

void NewFem(const char* GaussOrderCoarse, const char* GaussOrderFine);
void DeleteFem();

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

  NewFem("fifth", "twenty fifth");

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
  mlSol.AddSolution("pElem", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
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
  system.SetAssembleFunction(Assemble);
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
  InitPElement(mlSol, std::vector <double> {0., 0.4, 0.});

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
  unsigned n_timesteps = 1000;
  for(unsigned time_step = 1; time_step <= n_timesteps; time_step++) {

    system.CopySolutionToOldSolution();
    system.MGsolve();
    ProjectNewmarkDisplacemenet(mlSol);
    mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, time_step);
  }


  DeleteFem();
  return 0;

}

void Assemble(MultiLevelProblem& ml_prob) {

  // ml_prob is the global object from/to where get/set all the data
  // level is the level of the PDE system to be assembled
  // levelMax is the Maximum level of the MultiLevelProblem
  // assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  clock_t AssemblyTime = 0;
  clock_t start_time, end_time;

  //pointers and references

  TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("MPM_FSI");
  const unsigned  level = my_nnlin_impl_sys.GetLevelToAssemble();
  MultiLevelSolution* mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* mysolution = mlSol->GetSolutionLevel(level);     // pointer to the solution (level) object

  //Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* myLinEqSolver = my_nnlin_impl_sys._LinSolver[level];  // pointer to the equation (level) object

  Mesh* msh = ml_prob._ml_msh->GetLevel(level);     // pointer to the mesh (level) object
  elem* el = msh->el;   // pointer to the elem object in msh (level)
  SparseMatrix* myKK = myLinEqSolver->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* myRES =  myLinEqSolver->_RES;  // pointer to the global residual vector object in pdeSys (level)

  MatSetOption((static_cast<PetscMatrix*>(myKK))->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  myKK->zero();
  myRES->zero();

  AssembleSolid(ml_prob);
  AssembleSolidInterface(ml_prob);
  AssembleBoundaryLayer(ml_prob);
  AssembleBoundaryLayerProjection(ml_prob);
  AssembleFluid(ml_prob);
  AssembleGhostPenalty(ml_prob);

  myKK->close();
  myRES->close();

  //double tolerance = 1.0e-12 * myKK->linfty_norm();
  //myKK->RemoveZeroEntries(tolerance);

  end_time = clock();
  AssemblyTime += (end_time - start_time);


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


void InitPElement(MultiLevelSolution & mlSol, const std::vector <double> &x0) {

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution *sol = mlSol.GetSolutionLevel(level);
  Mesh *msh = mlSol._mlMesh->GetLevel(level);
  const unsigned dim = msh->GetDimension();

  unsigned iproc  = msh->processor_id();
  unsigned nprocs = msh->n_processors();

  unsigned pElemIndex = mlSol.GetIndex("pElem");

  struct {
    double dist2;
    int elem;
  } in, out;

  in.dist2 = 10e10;

  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    unsigned flag_mat = msh->GetElementMaterial(iel);
    if(flag_mat == 2) {
      unsigned nDofs = msh->GetElementDofNumber(iel, 2);
      unsigned idofX = msh->GetSolutionDof(nDofs - 1, iel, 2); //global dof for mesh coordinates
      double dist2 = 0.;
      for(unsigned k = 0; k < dim; k++) {
        dist2 += (x0[k] - (*msh->_topology->_Sol[k])(idofX)) * (x0[k] - (*msh->_topology->_Sol[k])(idofX));
      }
      if(dist2 < in.dist2) {
        in.dist2 = dist2;
        in.elem = iel;
      }
    }
  }

  MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, PETSC_COMM_WORLD);

  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    unsigned flag_mat = msh->GetElementMaterial(iel);
    if(flag_mat != 2) {
      sol->_Sol[pElemIndex]->set(iel, out.elem);
    }
    else {
      sol->_Sol[pElemIndex]->set(iel, -1.);
    }
  }
  sol->_Sol[pElemIndex]->close();


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

  for(unsigned k = 0; k < dim; k++) {
    *sol->_Sol[indexSolV[k]] = *sol->_Sol[indexSolU[k]];
  }

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



void NewFem(const char* GaussOrderCoarse, const char* GaussOrderFine) {

  fem[0][0][0] = new const elem_type_3D("hex", "linear", GaussOrderCoarse);
  fem[0][0][1] = new const elem_type_3D("hex", "quadratic", GaussOrderCoarse);
  fem[0][0][2] = new const elem_type_3D("hex", "biquadratic", GaussOrderCoarse);
  fem[0][0][3] = new const elem_type_3D("hex", "constant", GaussOrderCoarse);
  fem[0][0][4] = new const elem_type_3D("hex", "disc_linear", GaussOrderCoarse);

  fem[0][1][0] = new const elem_type_3D("tet", "linear", GaussOrderCoarse);
  fem[0][1][1] = new const elem_type_3D("tet", "quadratic", GaussOrderCoarse);
  fem[0][1][2] = new const elem_type_3D("tet", "biquadratic", GaussOrderCoarse);
  fem[0][1][3] = new const elem_type_3D("tet", "constant", GaussOrderCoarse);
  fem[0][1][4] = new const elem_type_3D("tet", "disc_linear", GaussOrderCoarse);

  fem[0][2][0] = new const elem_type_3D("wedge", "linear", GaussOrderCoarse);
  fem[0][2][1] = new const elem_type_3D("wedge", "quadratic", GaussOrderCoarse);
  fem[0][2][2] = new const elem_type_3D("wedge", "biquadratic", GaussOrderCoarse);
  fem[0][2][3] = new const elem_type_3D("wedge", "constant", GaussOrderCoarse);
  fem[0][2][4] = new const elem_type_3D("wedge", "disc_linear", GaussOrderCoarse);

  fem[0][3][0] = new const elem_type_2D("quad", "linear", GaussOrderCoarse);
  fem[0][3][1] = new const elem_type_2D("quad", "quadratic", GaussOrderCoarse);
  fem[0][3][2] = new const elem_type_2D("quad", "biquadratic", GaussOrderCoarse);
  fem[0][3][3] = new const elem_type_2D("quad", "constant", GaussOrderCoarse);
  fem[0][3][4] = new const elem_type_2D("quad", "disc_linear", GaussOrderCoarse);

  fem[0][4][0] = new const elem_type_2D("tri", "linear", GaussOrderCoarse);
  fem[0][4][1] = new const elem_type_2D("tri", "quadratic", GaussOrderCoarse);
  fem[0][4][2] = new const elem_type_2D("tri", "biquadratic", GaussOrderCoarse);
  fem[0][4][3] = new const elem_type_2D("tri", "constant", GaussOrderCoarse);
  fem[0][4][4] = new const elem_type_2D("tri", "disc_linear", GaussOrderCoarse);

  fem[0][5][0] = new const elem_type_1D("line", "linear", GaussOrderCoarse);
  fem[0][5][1] = new const elem_type_1D("line", "quadratic", GaussOrderCoarse);
  fem[0][5][2] = new const elem_type_1D("line", "biquadratic", GaussOrderCoarse);
  fem[0][5][3] = new const elem_type_1D("line", "constant", GaussOrderCoarse);
  fem[0][5][4] = new const elem_type_1D("line", "disc_linear", GaussOrderCoarse);


  fem[1][0][0] = new const elem_type_3D("hex", "linear", GaussOrderFine);
  fem[1][0][1] = new const elem_type_3D("hex", "quadratic", GaussOrderFine);
  fem[1][0][2] = new const elem_type_3D("hex", "biquadratic", GaussOrderFine);
  fem[1][0][3] = new const elem_type_3D("hex", "constant", GaussOrderFine);
  fem[1][0][4] = new const elem_type_3D("hex", "disc_linear", GaussOrderFine);

  fem[1][1][0] = new const elem_type_3D("tet", "linear", GaussOrderFine);
  fem[1][1][1] = new const elem_type_3D("tet", "quadratic", GaussOrderFine);
  fem[1][1][2] = new const elem_type_3D("tet", "biquadratic", GaussOrderFine);
  fem[1][1][3] = new const elem_type_3D("tet", "constant", GaussOrderFine);
  fem[1][1][4] = new const elem_type_3D("tet", "disc_linear", GaussOrderFine);

  fem[1][2][0] = new const elem_type_3D("wedge", "linear", GaussOrderFine);
  fem[1][2][1] = new const elem_type_3D("wedge", "quadratic", GaussOrderFine);
  fem[1][2][2] = new const elem_type_3D("wedge", "biquadratic", GaussOrderFine);
  fem[1][2][3] = new const elem_type_3D("wedge", "constant", GaussOrderFine);
  fem[1][2][4] = new const elem_type_3D("wedge", "disc_linear", GaussOrderFine);

  fem[1][3][0] = new const elem_type_2D("quad", "linear", GaussOrderFine);
  fem[1][3][1] = new const elem_type_2D("quad", "quadratic", GaussOrderFine);
  fem[1][3][2] = new const elem_type_2D("quad", "biquadratic", GaussOrderFine);
  fem[1][3][3] = new const elem_type_2D("quad", "constant", GaussOrderFine);
  fem[1][3][4] = new const elem_type_2D("quad", "disc_linear", GaussOrderFine);

  fem[1][4][0] = new const elem_type_2D("tri", "linear", GaussOrderFine);
  fem[1][4][1] = new const elem_type_2D("tri", "quadratic", GaussOrderFine);
  fem[1][4][2] = new const elem_type_2D("tri", "biquadratic", GaussOrderFine);
  fem[1][4][3] = new const elem_type_2D("tri", "constant", GaussOrderFine);
  fem[1][4][4] = new const elem_type_2D("tri", "disc_linear", GaussOrderFine);

  fem[1][5][0] = new const elem_type_1D("line", "linear", GaussOrderFine);
  fem[1][5][1] = new const elem_type_1D("line", "quadratic", GaussOrderFine);
  fem[1][5][2] = new const elem_type_1D("line", "biquadratic", GaussOrderFine);
  fem[1][5][3] = new const elem_type_1D("line", "constant", GaussOrderFine);
  fem[1][5][4] = new const elem_type_1D("line", "disc_linear", GaussOrderFine);

}

void DeleteFem() {
  for(unsigned i = 0; i < 2; i++) {
    for(unsigned j = 0; j < 6; j++) {
      for(unsigned k = 0; k < 5; k++) {
        delete fem[i][j][k];
      }
    }
  }
}
