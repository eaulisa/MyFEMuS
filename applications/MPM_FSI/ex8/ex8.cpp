
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


#include "./include/RefineElement.hpp"



double theta = 1.;
double af = 1. - theta;
// double pInf = (1. + af) / (2. - af);
double am = af - 0.1; //pInf / (1. + pInf);
double beta = 0.25 + 0.5 * (af - am);
double Gamma = 0.5 + (af - am);

const elem_type *fem[2][6][5];
RefineElement *refinedFem[6][3];

//#include "../../Nitsche/support/particleInit.hpp"
//#include "../../Nitsche/support/sharedFunctions.hpp"
#include "./include/assemblySolid.hpp"
#include "./include/assemblyBoundaryLayer.hpp"
#include "./include/assemblyFluid.hpp"


using namespace femus;



double dt = 1.;

double SetVariableTimeStep(const double time) {
  if(time < 20.) dt = 0.05;
  else dt = 0.01;

  return dt;
}

void Assemble(MultiLevelProblem& ml_prob);

void InitPElement(MultiLevelSolution & mlSol);
void BuildFlagSolidRegion(MultiLevelSolution & mlSol);
void BuildFlagFluidRegion(MultiLevelSolution & mlSol);
void ProjectNewmarkDisplacemenet(MultiLevelSolution & mlSol);

void NewFem(const char* GaussOrderCoarse, const char* GaussOrderFine, const unsigned &lmax);
void DeleteFem();

bool SetBoundaryCondition(const std::vector < double >&x, const char name[], double &value, const int facename, const double t) {
  bool dirichlet = true;
  value = 0.;

  const double Ubar = 1.0;    // guess?
  const double L = 0.41;
  const double H = 2.5;


  if(!strcmp(name, "UX")) {
    if(1 == facename) {     //inflow
      if(t < 2.0) {
        value = 1.5 * Ubar * 4.0 / 0.1681 * (x[1] + 0.2) * (-x[1] + 0.21) * 0.5 * (1. - cos(0.5 * M_PI * t));
      }
      else {
        value = 1.5 * Ubar * 4.0 / 0.1681 * (x[1] + 0.2) * (-x[1] + 0.21);
      }
    }
    else if(2 == facename || 6 == facename) {     //outflow and porous media Neumann
      dirichlet = false;
      value = 0.;
    }
  }

  else if(!strcmp(name, "UY")) {
    if(6 == facename) {     //porous media Neumann
      dirichlet = false;
      value = 0.;
    }
  }

  else if(!strcmp(name, "P")) {
    if(2 != facename)  {
      dirichlet = 0;
      value = 0;
    }
  }

  return dirichlet;

}


int main(int argc, char **args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;
  unsigned numberOfUniformLevels = 3; //for refinement in 3D
  //unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0 ;

  double Lref = 1.;
  double Uref = 1.;
//   double rhos = 7850;
//   double rhof = 1000;
//   double nu = 0.3;
//   double E = 2.e05;
//   double muf = 1.0e-3;

  unsigned lmax = 2;

  NewFem("fifth", "twenty fifth", lmax);

//   std::vector <std::vector <double> > xv = {{ -1., 1., 1., -1.}, { -1., -1., 1., 1.}};
//   std::vector <double>  il(lmax);
//
//   refinedFem[3][0]->InitElement1(xv);
//
//   double area = 0.;
//
//   int l = 0;
//   for(il[l] = 0; il[l] < refinedFem[3][0]->GetLevelSize(l); il[l]++) {
//     if(l < lmax - 1 && rand() > RAND_MAX/2) { //refine
//     //  std::cout << l <<" "<< il[l] <<std::endl;
//       refinedFem[3][0]->BuildElement1Prolongation(l, il[l]);  // this creates the next level for the il[l] element
//       l++;// move to the next level
//       il[l] = -1; // this reset the il[l++] index to zero (notice that the il[l]++, will augment it by one)
//     }
//     else{ //integrate
//
//      // std::cout << l <<" "<< il[l] <<std::endl;
//       const std::vector<std::vector<double>> & xvl = refinedFem[3][0]->GetElement1NodeCoordinates(l, il[l]);
//
//       area += (xvl[0][2] - xvl[0][0]) * (xvl[1][2] - xvl[1][0]);
//
//       for(unsigned i = 0; i < xvl[0].size();i ++){
//         for(unsigned k = 0; k < xvl.size();k++){
//           std::cout << xvl[k][i] << " ";
//         }
//         std::cout << std::endl;
//       }
//       std::cout << std::endl;
//
//       while(l > 0 && il[l] + 1 == refinedFem[3][0]->GetLevelSize(l)) {
//         //once all the elements have been integrated move to the previous level
//         l--;
//       }
//     }
//   }
//
//   std::cout << area <<std::endl;
//
//   DeleteFem();
//   return 1;


  double rhof = 1000.;
  double muf = 1.;
  double rhos = 10000.;
  double nu = 0.4;
  double E = 1400000;

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

  mlSol.Initialize("All");

  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

  // ******* Set boundary conditions *******
  mlSol.GenerateBdc("UX", "Time_dependent");
  mlSol.GenerateBdc("UY", "Steady");
  if(dim > 2) mlSol.GenerateBdc("UZ", "Steady");

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

  InitPElement(mlSol);

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
  unsigned n_timesteps = 2000;
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

  BuildFlagFluidRegion(*mlSol);

  AssembleFluid(ml_prob);
  AssembleGhostPenalty(ml_prob);

  myKK->close();
  myRES->close();

  //double tolerance = 1.0e-12 * myKK->linfty_norm();
  //myKK->RemoveZeroEntries(tolerance);

//   PetscInt row = 16133 * 0 + 462;
//   PetscInt ncols;
//   const PetscInt *cols;
//   const PetscScalar *vals;
//   MatGetRow((static_cast<PetscMatrix*>(myKK))->mat(),row,&ncols,&cols,&vals);

//   for(unsigned i = 0; i<ncols; i++){
//     unsigned minus = 0;
//     if(cols[i] > 16133 * 2) minus = 2;
//     else if(cols[i] > 16133) minus = 1;
//     std::cout << cols[i] - minus * 16133 << " " << vals[i]<< std::endl;
//   }

// MatRestoreRow((static_cast<PetscMatrix*>(myKK))->mat(),row,&ncols,&cols,&vals);

  //myKK->draw();

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



void BuildFlagFluidRegion(MultiLevelSolution & mlSol) {

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution *sol = mlSol.GetSolutionLevel(level);
  Mesh *msh = mlSol._mlMesh->GetLevel(level);
  unsigned iproc = msh->processor_id();

  unsigned eflagIndex = mlSol.GetIndex("eflag");
  unsigned nflagIndex = mlSol.GetIndex("nflag");

  unsigned nflagType = mlSol.GetSolutionType(nflagIndex);

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned eFlag = static_cast <unsigned>(floor((*sol->_Sol[eflagIndex])(iel) + 0.5));

    if(eFlag < 3) {
      unsigned nDofs = msh->GetElementDofNumber(iel, nflagType);
      for(unsigned  i = 0; i < nDofs; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, nflagType);
        sol->_Sol[nflagIndex]->set(idof, 2);
      }
    }
  }
  sol->_Sol[nflagIndex]->close();




  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned eFlag = static_cast <unsigned>(floor((*sol->_Sol[eflagIndex])(iel) + 0.5));

    if(eFlag == 0) {
      unsigned nDofs = msh->GetElementDofNumber(iel, nflagType);
      for(unsigned  i = 0; i < nDofs; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, nflagType);
        sol->_Sol[nflagIndex]->set(idof, 0);
      }
    }
  }
  sol->_Sol[nflagIndex]->close();


  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned eFlag = static_cast <unsigned>(floor((*sol->_Sol[eflagIndex])(iel) + 0.5));

    if(eFlag == 1) {
      unsigned nDofs = msh->GetElementDofNumber(iel, nflagType);
      for(unsigned  i = 0; i < nDofs; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, nflagType);
        sol->_Sol[nflagIndex]->set(idof, 1);
      }
    }
  }
  sol->_Sol[nflagIndex]->close();


  unsigned counterAll = 1;
  while(counterAll > 0) {
    unsigned iprocCounter = 0;

    for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

      unsigned eFlag = static_cast <unsigned>(floor((*sol->_Sol[eflagIndex])(iel) + 0.5));

      if(eFlag == 2) {
        unsigned nDofs = msh->GetElementDofNumber(iel, nflagType);
        for(unsigned  i = 0; i < nDofs; i++) {
          unsigned idof = msh->GetSolutionDof(i, iel, nflagType);
          unsigned inodeFlag = (*sol->_Sol[nflagIndex])(idof);
          if(inodeFlag == 1 || inodeFlag == 10) {
            sol->_Sol[eflagIndex]->set(iel, 10);
            break;
          }
        }
      }
    }
    sol->_Sol[eflagIndex]->close();

    for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
      unsigned eFlag = static_cast <unsigned>(floor((*sol->_Sol[eflagIndex])(iel) + 0.5));
      if(eFlag == 10) {
        unsigned nDofs = msh->GetElementDofNumber(iel, nflagType);
        for(unsigned  i = 0; i < nDofs; i++) {
          unsigned idof = msh->GetSolutionDof(i, iel, nflagType);
          unsigned inodeFlag = (*sol->_Sol[nflagIndex])(idof);
          if(inodeFlag == 2) {
            iprocCounter++;
            sol->_Sol[nflagIndex]->set(idof, 10);
          }
        }
      }
    }
    sol->_Sol[nflagIndex]->close();
    MPI_Allreduce(&iprocCounter, &counterAll, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);
  }
}




void InitPElement(MultiLevelSolution & mlSol) {

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution *sol = mlSol.GetSolutionLevel(level);
  Mesh *msh = mlSol._mlMesh->GetLevel(level);
  const unsigned dim = msh->GetDimension();

  unsigned iproc  = msh->processor_id();
  unsigned nprocs = msh->n_processors();

  unsigned pElemIndex = mlSol.GetIndex("pElem");

  std::vector <double> xc(dim);

  struct {
    double dist2;
    int elem;
  } in, out;

  for(unsigned kproc = 0; kproc < nprocs; kproc++) {
    for(unsigned iel = msh->_elementOffset[kproc]; iel < msh->_elementOffset[kproc + 1]; iel++) {
      unsigned flag_mat;
      if(iproc == kproc) {
        flag_mat = msh->GetElementMaterial(iel);
        if(flag_mat > 2) {
          unsigned nDofs = msh->GetElementDofNumber(iel, 2);
          unsigned idofX = msh->GetSolutionDof(nDofs - 1, iel, 2); //center of the element dof (coordinates)

          for(unsigned k = 0; k < dim; k++) {
            xc[k] =  (*msh->_topology->_Sol[k])(idofX); // x-y-z coordinates
          }
        }
        else {
          sol->_Sol[pElemIndex]->set(iel, -1.);
        }
      }


      MPI_Bcast(&flag_mat, 1, MPI_UNSIGNED, kproc, PETSC_COMM_WORLD);

      if(flag_mat > 2) {
        MPI_Bcast(xc.data(), xc.size(), MPI_DOUBLE, kproc, PETSC_COMM_WORLD);

        in.dist2 = 1000;
        in.elem = INT_MAX;
        for(unsigned jel = msh->_elementOffset[iproc]; jel < msh->_elementOffset[iproc + 1]; jel++) {

          unsigned jflag_mat = msh->GetElementMaterial(jel);
          if(jflag_mat == 2) {

            unsigned nDofs = msh->GetElementDofNumber(jel, 2);
            unsigned jdofX = msh->GetSolutionDof(nDofs - 1, jel, 2); //center of the element dof (coordinates)

            double dist2 = 0.;
            for(unsigned k = 0; k < dim; k++) {
              dist2 += (xc[k] - (*msh->_topology->_Sol[k])(jdofX)) * (xc[k] - (*msh->_topology->_Sol[k])(jdofX));
            }
            if(dist2 < in.dist2) {
              in.dist2 = dist2;
              in.elem = jel;
            }
          }
        }

        MPI_Reduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, kproc, PETSC_COMM_WORLD);

        if(iproc == kproc) {
          sol->_Sol[pElemIndex]->set(iel, out.elem);
        }
      }
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

  for(int i = msh->_dofOffset[solType][iproc]; i < msh->_dofOffset[solType][iproc + 1]; i++) {
    unsigned nodeFlag = (*sol->_Sol[nflagIndex])(i);  
    if(nodeFlag >= 4 && nodeFlag < 8) {
      for(unsigned k = 0 ; k < dim ; k++) {

        double Dnew = (*sol->_Sol[indexSolU[k]])(i);
        double Dold = (*sol->_SolOld[indexSolU[k]])(i);

        double Vold = (*sol->_Sol[indexSolV[k]])(i);
        double Aold = (*sol->_Sol[indexSolA[k]])(i);
        double Anew = (Dnew - Dold) / (beta * dt * dt) - Vold / (beta * dt) + Aold * (beta - 0.5) / beta;
        double Vnew = Vold + (1 - Gamma) * dt * Aold + Gamma * dt * Anew;

        sol->_Sol[indexSolD[k]]->set(i, (1. - af) * Dnew + af * Dold);
        sol->_Sol[indexSolV[k]]->set(i, Vnew);
        sol->_Sol[indexSolA[k]]->set(i, Anew);

      }
    }
    else {
      for(unsigned k = 0 ; k < dim ; k++) {
        double Vnew = (*sol->_Sol[indexSolU[k]])(i);
        sol->_Sol[indexSolV[k]]->set(i, Vnew);
      }
    }
  }

  for(unsigned k = 0 ; k < dim ; k++) {
    sol->_Sol[indexSolD[k]]->close();
    sol->_Sol[indexSolV[k]]->close();
    sol->_Sol[indexSolA[k]]->close();
  }
}



void NewFem(const char* GaussOrderCoarse, const char* GaussOrderFine, const unsigned &lmax) {

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


  refinedFem[0][0] = new RefineElement(lmax, "hex", "linear", GaussOrderCoarse, "zero", "lobatto");
  refinedFem[0][1] = new RefineElement(lmax, "hex", "quadratic", GaussOrderCoarse, "zero", "lobatto");
  refinedFem[0][2] = new RefineElement(lmax, "hex", "biquadratic", GaussOrderCoarse, "zero", "lobatto");

  refinedFem[1][0] = new RefineElement(lmax, "tet", "linear", GaussOrderCoarse, "zero", "lobatto");
  refinedFem[1][1] = new RefineElement(lmax, "tet", "quadratic", GaussOrderCoarse, "zero", "lobatto");
  refinedFem[1][2] = new RefineElement(lmax, "tet", "biquadratic", GaussOrderCoarse, "zero", "lobatto");

  refinedFem[2][0] = new RefineElement(lmax, "wedge", "linear", GaussOrderCoarse, "zero", "lobatto");
  refinedFem[2][1] = new RefineElement(lmax, "wedge", "quadratic", GaussOrderCoarse, "zero", "lobatto");
  refinedFem[2][2] = new RefineElement(lmax, "wedge", "biquadratic", GaussOrderCoarse, "zero", "lobatto");

  refinedFem[3][0] = new RefineElement(lmax, "quad", "linear", GaussOrderCoarse, "zero", "lobatto");
  refinedFem[3][1] = new RefineElement(lmax, "quad", "quadratic", GaussOrderCoarse, "zero", "lobatto");
  refinedFem[3][2] = new RefineElement(lmax, "quad", "biquadratic", GaussOrderCoarse, "zero", "lobatto");

  refinedFem[4][0] = new RefineElement(lmax, "tri", "linear", GaussOrderCoarse, "zero", "lobatto");
  refinedFem[4][1] = new RefineElement(lmax, "tri", "quadratic", GaussOrderCoarse, "zero", "lobatto");
  refinedFem[4][2] = new RefineElement(lmax, "tri", "biquadratic", GaussOrderCoarse, "zero", "lobatto");

  refinedFem[5][0] = new RefineElement(lmax, "line", "linear", GaussOrderCoarse, "zero", "lobatto");
  refinedFem[5][1] = new RefineElement(lmax, "line", "quadratic", GaussOrderCoarse, "zero", "lobatto");
  refinedFem[5][2] = new RefineElement(lmax, "line", "biquadratic", GaussOrderCoarse, "zero", "lobatto");
}

void DeleteFem() {
  for(unsigned i = 0; i < 2; i++) {
    for(unsigned j = 0; j < 6; j++) {
      for(unsigned k = 0; k < 5; k++) {
        delete fem[i][j][k];
      }
    }
  }

  for(unsigned j = 0; j < 6; j++) {
    for(unsigned k = 0; k < 3; k++) {
      delete refinedFem[j][k];
    }
  }


}
