
#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "TransientSystem.hpp"
#include "NonLinearImplicitSystem.hpp"

#include "NumericVector.hpp"
#include "adept.h"

#include "petsc.h"
#include "petscmat.h"
#include "PetscMatrix.hpp"

#include "slepceps.h"

#include "../include/nonlocal_assembly_adaptive.hpp"

//2D NONLOCAL EX : nonlocal diffusion for a body with different material properties

using namespace femus;

double InitalValueU(const std::vector < double >& x) {
  double value;

  //value =  x[0] * x[0] + x[1] * x[1]; //consistency 
  value =  x[0] * x[0] * x[0] + x[1] * x[1] * x[1]; //cubic 
  //value = x[0] * x[0] * x[0] * x[0] + x[1] * x[1] * x[1] * x[1]; //quartic




  return value;
}

void GetL2Norm(MultiLevelSolution & mlSol, MultiLevelSolution & mlSolFine);

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {

  bool dirichlet = true;
  //value = x[0] * x[0] + x[1] * x[1]; //consistency
  value = x[0] * x[0] * x[0] + x[1] * x[1] * x[1]; //cubic
  //value = x[0] * x[0] * x[0] * x[0] + x[1] * x[1] * x[1] * x[1]; //quartic


  return dirichlet;
}

//unsigned numberOfUniformLevels = 2; //consistency
unsigned numberOfUniformLevels = 1; //cubic-quartic 2->6 //cubic Marta4Quad Tri Mix 
//unsigned numberOfUniformLevels = 2; //cubic-quartic 2->4 mappa a 4->6 //cubic Marta4Fine 


unsigned numberOfUniformLevelsFine = 1;

int main(int argc, char** argv) {

  clock_t total_time = clock();

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, argv, MPI_COMM_WORLD);

  MultiLevelMesh mlMsh;
  MultiLevelMesh mlMshFine;

  double scalingFactor = 1.;
  unsigned numberOfSelectiveLevels = 0;



  //char fileName[100] = "../input/martaTest4.neu"; // good form 2->6 in serial but in parallel use martaTest4Fine
  //char fileName[100] = "../input/martaTest4Fine.neu"; // works till 144 nprocs +2
  char fileName[100] = "../input/martaTest4Finer.neu"; // works till 144 nprocs +4
  //char fileName[100] = "../input/martaTest4Tri.neu";
  //char fileName[100] = "../input/martaTest4Unstr.neu"; // works till 144 nprocs
    
  mlMsh.ReadCoarseMesh(fileName, "fifth", scalingFactor);
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , NULL);

  //consistency
  mlMshFine.ReadCoarseMesh(fileName, "fifth", scalingFactor);
  mlMshFine.RefineMesh(numberOfUniformLevelsFine + numberOfSelectiveLevels, numberOfUniformLevelsFine , NULL);

//   mlMsh.EraseCoarseLevels (numberOfUniformLevels - 1);
//   mlMshFine.EraseCoarseLevels (numberOfUniformLevelsFine - 1);

  unsigned dim = mlMsh.GetDimension();

  MultiLevelSolution mlSol(&mlMsh);
  MultiLevelSolution mlSolFine(&mlMshFine);

  // add variables to mlSol
  FEOrder femType = SERENDIPITY;
  //FEOrder femType = FIRST;
  
  std::vector < std::string > femTypeName = {"zero", "linear", "quadratic", "biquadratic"};

  mlSol.AddSolution("u", LAGRANGE,  femType, 0);
  mlSol.AddSolution("u_local", LAGRANGE,  femType, 0);
  mlSol.AddSolution("u_exact", LAGRANGE,  femType, 0, false);

  mlSol.AddSolution("cnt", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);


  mlSolFine.AddSolution("u", LAGRANGE,  femType, 0);
  mlSolFine.AddSolution("up", LAGRANGE, femType, 0, false);
  mlSolFine.AddSolution("cnt", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);

  mlSol.Initialize("All");
  mlSolFine.Initialize("All");

  mlSol.Initialize("u_exact", InitalValueU);
  mlSol.Initialize("u", InitalValueU);

  // ******* Set boundary conditions *******
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSolFine.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("All");
  mlSolFine.GenerateBdc("All");

  // ******* Set volume constraints for the nonlocal *******
  std::vector<unsigned> volumeConstraintFlags(3);
  volumeConstraintFlags[0] = 5;
  volumeConstraintFlags[1] = 6;
  volumeConstraintFlags[2] = 9;

  unsigned soluIndex = mlSol.GetIndex("u");
  mlSol.GenerateBdcOnVolumeConstraint(volumeConstraintFlags, soluIndex, 0);

  unsigned soluLocalIndex = mlSol.GetIndex("u_local");
  mlSol.GenerateBdcOnVolumeConstraint(volumeConstraintFlags, soluLocalIndex, 0);

  unsigned soluIndexFine = mlSolFine.GetIndex("u");
  mlSolFine.GenerateBdcOnVolumeConstraint(volumeConstraintFlags, soluIndexFine, 0);

  //BEGIN assemble and solve nonlocal problem
  MultiLevelProblem ml_prob(&mlSol);

  // ******* Add FEM system to the MultiLevel problem *******
  LinearImplicitSystem& system = ml_prob.add_system < LinearImplicitSystem > ("NonLocal");
  system.AddSolutionToSystemPDE("u");

  // ******* System FEM Assembly *******
  system.SetAssembleFunction(AssembleNonLocalRefined);
  //system.SetAssembleFunction(AssembleNonLocalSys);
  system.SetMaxNumberOfLinearIterations(1);
  // ******* set MG-Solver *******
  system.SetMgType(V_CYCLE);

  system.SetAbsoluteLinearConvergenceTolerance(1.e-50);
  //   system.SetNonLinearConvergenceTolerance(1.e-9);
//   system.SetMaxNumberOfNonLinearIterations(20);

  system.SetNumberPreSmoothingStep(1);
  system.SetNumberPostSmoothingStep(1);

  // ******* Set Preconditioner *******
  system.SetLinearEquationSolverType(FEMuS_DEFAULT);

  system.SetSparsityPatternMinimumSize(40000u);    //TODO tune

  system.init();

  // ******* Set Smoother *******
  system.SetSolverFineGrids(RICHARDSON);
//   system.SetRichardsonScaleFactor(0.7);

  system.SetPreconditionerFineGrids(ILU_PRECOND);

  system.SetTolerances(1.e-40, 1.e-40, 1.e+50, 100);

// ******* Solution *******

  system.MGsolve();

  //END assemble and solve nonlocal problem

  //BEGIN assemble and solve local problem
  MultiLevelProblem ml_prob2(&mlSol);

  // ******* Add FEM system to the MultiLevel problem *******
  LinearImplicitSystem& system2 = ml_prob2.add_system < LinearImplicitSystem > ("Local");
  system2.AddSolutionToSystemPDE("u_local");

  // ******* System FEM Assembly *******
  system2.SetAssembleFunction(AssembleLocalSys);
  system2.SetMaxNumberOfLinearIterations(1);
  // ******* set MG-Solver *******
  system2.SetMgType(V_CYCLE);

  system2.SetAbsoluteLinearConvergenceTolerance(1.e-50);

  system2.SetNumberPreSmoothingStep(1);
  system2.SetNumberPostSmoothingStep(1);

  system2.SetLinearEquationSolverType(FEMuS_DEFAULT);

  system2.init();

  // ******* Set Smoother *******
  system2.SetSolverFineGrids(RICHARDSON);

  system2.SetPreconditionerFineGrids(ILU_PRECOND);

  system2.SetTolerances(1.e-20, 1.e-20, 1.e+50, 100);

// ******* Solution *******

  system2.MGsolve();

  //END assemble and solve local problem

  //BEGIN assemble and solve fine nonlocal problem
  MultiLevelProblem ml_probFine(&mlSolFine);

  // ******* Add FEM system to the MultiLevel problem *******
  LinearImplicitSystem& systemFine = ml_probFine.add_system < LinearImplicitSystem > ("NonLocal");
  systemFine.AddSolutionToSystemPDE("u");

  // ******* System FEM Assembly *******
  //systemFine.SetAssembleFunction(AssembleNonLocalSys);
  systemFine.SetAssembleFunction(AssembleNonLocalRefined);
  systemFine.SetMaxNumberOfLinearIterations(1);
  // ******* set MG-Solver *******
  systemFine.SetMgType(V_CYCLE);

  systemFine.SetAbsoluteLinearConvergenceTolerance(1.e-50);
  //   systemFine.SetNonLinearConvergenceTolerance(1.e-9);
  //   systemFine.SetMaxNumberOfNonLinearIterations(20);

  systemFine.SetNumberPreSmoothingStep(1);
  systemFine.SetNumberPostSmoothingStep(1);

  // ******* Set Preconditioner *******
  systemFine.SetLinearEquationSolverType(FEMuS_DEFAULT);

  systemFine.SetSparsityPatternMinimumSize(5000u);    //TODO tune


  systemFine.init();

  // ******* Set Smoother *******
  systemFine.SetSolverFineGrids(RICHARDSON);
  // systemFine.SetRichardsonScaleFactor(0.7);

  systemFine.SetPreconditionerFineGrids(ILU_PRECOND);

  systemFine.SetTolerances(1.e-20, 1.e-20, 1.e+50, 100);

// ******* Solution *******

  systemFine.MGsolve(); //TODO

  //END assemble and solve nonlocal problem


  //BEGIN compute errors
  GetL2Norm(mlSol, mlSolFine);
  //END compute errors

  // ******* Print solution *******
  mlSol.SetWriter(VTK);
  std::vector<std::string> print_vars;
  print_vars.push_back("All");
  mlSol.GetWriter()->SetDebugOutput(true);
  mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, femTypeName[femType].c_str(), print_vars, 0);

  mlSolFine.SetWriter(VTK);
  std::vector<std::string> print_vars2;
  print_vars2.push_back("All");
  mlSolFine.GetWriter()->SetDebugOutput(true);
  mlSolFine.GetWriter()->Write(DEFAULT_OUTPUTDIR, femTypeName[femType].c_str(), print_vars2, 1);

  std::cout << std::endl << " total CPU time : " << std::setw(11) << std::setprecision(6) << std::fixed
            << static_cast<double>((clock() - total_time)) / CLOCKS_PER_SEC << " s" << std::endl;

  return 0;

} //end main


void GetL2Norm(MultiLevelSolution & mlSol, MultiLevelSolution & mlSolFine) {

  const unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;
  Mesh* msh = mlSol._mlMesh->GetLevel(level);
  Solution* sol  = mlSol.GetSolutionLevel(level);


  const unsigned  dim = msh->GetDimension();

  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  double error_solExact_norm2 = 0.;

  double error_solExact_local_norm2 = 0.;

  double error_solLocal_norm2 = 0.;

  double solNonlocal_norm2 = 0.;

  double solLocal_norm2 = 0.;

  double sol_exact_norm2 = 0.;

  unsigned soluIndex;
  soluIndex = mlSol.GetIndex("u");
  unsigned soluType = mlSol.GetSolutionType(soluIndex);

  unsigned soluIndexLocal;
  soluIndexLocal = mlSol.GetIndex("u_local");

  unsigned    iproc = msh->processor_id();
  unsigned    nprocs = msh->n_processors();

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);

    vector < vector < double > > x1(dim);

    for(int i = 0; i < dim; i++) {
      x1[i].resize(nDofx);
    }

    vector < double >  soluNonLoc(nDofu);
    vector < double >  soluLoc(nDofu);

    for(unsigned i = 0; i < nDofu; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);
      soluNonLoc[i] = (*sol->_Sol[soluIndex])(solDof);
      soluLoc[i] = (*sol->_Sol[soluIndexLocal])(solDof);
    }

    for(unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);

      for(unsigned jdim = 0; jdim < dim; jdim++) {
        x1[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);
      }
    }

    vector <double> phi;  // local test function
    vector <double> phi_x; // local test function first order partial derivatives
    double weight; // gauss point weight

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian(x1, ig, weight, phi, phi_x);
      double soluNonLoc_gss = 0.;
      double soluLoc_gss = 0.;
      double soluExact_gss = 0.;
      double x_gss = 0.;
      double y_gss = 0.;

      for(unsigned i = 0; i < nDofu; i++) {
        soluNonLoc_gss += phi[i] * soluNonLoc[i];
        soluLoc_gss += phi[i] * soluLoc[i];
        x_gss += phi[i] * x1[0][i]; // this is x at the Gauss point
        y_gss += phi[i] * x1[1][i]; // this is y at the Gauss point
      }




      
      //soluExact_gss = x_gss * x_gss + y_gss * y_gss; //consistency
      soluExact_gss = x_gss * x_gss * x_gss + y_gss * y_gss * y_gss; // cubic
      //soluExact_gss = x_gss * x_gss * x_gss * x_gss + y_gss * y_gss * y_gss * y_gss; // quartic

      error_solExact_norm2 += (soluNonLoc_gss - soluExact_gss) * (soluNonLoc_gss - soluExact_gss) * weight;

      error_solExact_local_norm2 += (soluLoc_gss - soluExact_gss) * (soluLoc_gss - soluExact_gss) * weight;

      error_solLocal_norm2 += (soluNonLoc_gss - soluLoc_gss) * (soluNonLoc_gss - soluLoc_gss) * weight;

      solNonlocal_norm2 += soluNonLoc_gss * soluNonLoc_gss * weight;

      solLocal_norm2 += soluLoc_gss * soluLoc_gss * weight;

      sol_exact_norm2 += soluExact_gss * soluExact_gss * weight;
    }
  }

  double norm2 = 0.;
  MPI_Allreduce(&error_solExact_norm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  double norm = sqrt(norm2);
  std::cout.precision(14);
  std::cout << "L2 norm of ERROR: Nonlocal - exact = " << norm << std::endl;

  norm2 = 0.;
  MPI_Allreduce(&error_solExact_local_norm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  norm = sqrt(norm2);
  std::cout.precision(14);
  std::cout << "L2 norm of ERROR: Local - exact = " << norm << std::endl;

  norm2 = 0.;
  MPI_Allreduce(&error_solLocal_norm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  norm = sqrt(norm2);
  std::cout.precision(14);
  std::cout << "L2 norm of ERROR: Nonlocal - local = " << norm << std::endl;

  norm2 = 0.;
  MPI_Allreduce(&solNonlocal_norm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  norm = sqrt(norm2);
  std::cout.precision(14);
  std::cout << "L2 norm of NONLOCAL soln = " << norm << std::endl;

  norm2 = 0.;
  MPI_Allreduce(&solLocal_norm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  norm = sqrt(norm2);
  std::cout.precision(14);
  std::cout << "L2 norm of LOCAL soln = " << norm << std::endl;

  norm2 = 0.;
  MPI_Allreduce(&sol_exact_norm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  norm = sqrt(norm2);
  std::cout.precision(14);
  std::cout << "L2 norm of EXACT soln = " << norm << std::endl;


  //BEGIN nonlocal fine - coarse L2 norm on fine grid

  std::cout << "------------------------------------- " << std::endl;

  const unsigned levelFine = mlSolFine._mlMesh->GetNumberOfLevels() - 1;
  if(levelFine >= level) {
    Mesh* mshFine = mlSolFine._mlMesh->GetLevel(levelFine);
    Solution* solFine  = mlSolFine.GetSolutionLevel(levelFine);

    unsigned soluIndexFine;
    soluIndexFine = mlSolFine.GetIndex("u");

    unsigned solupIndexFine;
    solupIndexFine = mlSolFine.GetIndex("up");

    Solution* solFineL = mlSolFine.GetSolutionLevel(level);
    *solFineL->_Sol[solupIndexFine] = *sol->_Sol[soluIndex];
    for(unsigned l = level + 1; l <= levelFine; l++) {
      Mesh* mshFineL = mlSolFine._mlMesh->GetLevel(l);
      Solution* solFineL = mlSolFine.GetSolutionLevel(l);
      Solution* solFineLm1 = mlSolFine.GetSolutionLevel(l - 1);
      solFineL->_Sol[solupIndexFine]->matrix_mult(*solFineLm1->_Sol[solupIndexFine], *mshFineL->GetCoarseToFineProjection(soluType));
    }

    double error_NonLocCoarse_NonLocFine_norm2 = 0.;

    double solNonlocalFine_norm2 = 0.;


    for(int iel = solFine->GetMesh()->_elementOffset[iproc]; iel < solFine->GetMesh()->_elementOffset[iproc + 1]; iel ++) {

      short unsigned ielGeom = mshFine->GetElementType(iel);
      unsigned nDofs  = mshFine->GetElementDofNumber(iel, soluType);

      vector < double >  solu(nDofs);
      vector < double >  soluP(nDofs);

      std::vector < std::vector <double> > x(dim);

      for(int k = 0; k < dim; k++) {
        x[k].assign(nDofs, 0.);
      }

      for(unsigned i = 0; i < nDofs; i++) {
        unsigned solDof = mshFine->GetSolutionDof(i, iel, soluType);
        solu[i] = (*solFine->_Sol[soluIndexFine])(solDof);
        soluP[i] = (*solFine->_Sol[solupIndexFine])(solDof);
        unsigned xDof  = mshFine->GetSolutionDof(i, iel, xType);
        for(unsigned k = 0; k < dim; k++) {
          x[k][i] = (*mshFine->_topology->_Sol[k])(xDof);
        }
      }

      const double* phi;  // local test function
      double weight; // gauss point weight
      unsigned igNumber = mshFine->_finiteElement[ielGeom][soluType]->GetGaussPointNumber();

      // *** Gauss point loop ***
      for(unsigned ig = 0; ig < igNumber; ig++) {

        mshFine->_finiteElement[ielGeom][soluType]->GetGaussQuantities(x, ig, weight, phi);
        double soluP_gss = 0.;
        double solu_gss = 0.;
        for(unsigned i = 0; i < nDofs; i++) {
          solu_gss += phi[i] * solu[i];
          soluP_gss += phi[i] * soluP[i];
        }

        error_NonLocCoarse_NonLocFine_norm2 += (solu_gss - soluP_gss) * (solu_gss - soluP_gss) * weight;
        solNonlocalFine_norm2 += solu_gss * solu_gss * weight;

      }

    }


    norm2 = 0.;
    MPI_Allreduce(&error_NonLocCoarse_NonLocFine_norm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    norm = sqrt(norm2);
    std::cout.precision(14);
    std::cout << "L2 norm of ERROR: Nonlocal - Nonlocal Fine = " << norm << std::endl;

    norm2 = 0.;
    MPI_Allreduce(&solNonlocalFine_norm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    norm = sqrt(norm2);
    std::cout.precision(14);
    std::cout << "L2 norm of NONLOCAL FINE soln = " << norm << std::endl;
  }
  //END

}







