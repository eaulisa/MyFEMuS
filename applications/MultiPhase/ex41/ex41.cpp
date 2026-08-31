
#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"

#include "TransientSystem.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "adept.h"
#include "PolynomialBases.hpp"
#include "CutFemWeight.hpp"
#include "CDWeights.hpp"

#include "petsc.h"
#include "petscmat.h"
#include "PetscMatrix.hpp"
#include "Fem.hpp"

#include "MyMatrix.hpp"
#include "VTKWriter.hpp"

#include <gperftools/profiler.h>

using namespace femus;
#include "../includeLS/Shape.hpp"
#include "../includeLS/Boundary.hpp"
#include "../includeLS/BBoxToIel.hpp"
#include "../includeLS/ElementTopology.hpp"
#include "../includeLS/GradientApproximation.hpp"
#include "../includeLS/LevelSetMarkers.hpp"
#include "../includeLS/Mollifier.hpp"
#include "../includeLS/Psi.hpp"
#include "../includeLS/Reinit.hpp"
#include "../includeLS/RungeKutta.hpp"

const double mu1 = 0.1;
const double mu2 = 10.;
const double rho1 = 1.;
const double rho2 = 1000;
const double sigma = 1.96;
const double gravity = 0.;//-0.98;

// //Parasitic Test
// const double mu1 = 0.1; // TODO Sandro put mu_1 = mu_2 = 0.005
// const double mu2 = 0.1;
// const double rho1 = 100.;
// const double rho2 = 100.;
// const double sigma = 3.; // ???
// const double gravity = 0.;

std::vector <double> g;

#define RADIUS 0.15
#define XG 0.
#define YG 0.
#define ZG 0.

typedef double TypeIO;
typedef cpp_bin_float_oct TypeA;
typedef cpp_bin_float_oct oct;

// CutFemWeight <double, double> quad = CutFemWeight<double, double>(QUAD, 5, "legendre");
CutFemWeight <TypeIO, TypeA> quad  = CutFemWeight<TypeIO, TypeA >(QUAD, 5, "legendre");
CutFemWeight <TypeIO, TypeA> tri  = CutFemWeight<TypeIO, TypeA >(TRI, 5, "legendre");

const std::vector< CutFemWeight <TypeIO, TypeA> *> cfw = {&quad, &quad, &quad, &quad, &tri};

Fem fem = Fem(quad.GetGaussQuadratureOrder(), quad.GetDimension());

#include "../includeLS/Stabilization.hpp"
#include "../includeLS/GhostPenalty.hpp"
#include "../includeLS/GhostPenaltyDGP.hpp"

const RungeKutta::VelKind velocityType = RungeKutta::VelKind::Zero;
//const RungeKutta::VelKind velocityType = RungeKutta::VelKind::Translation;
//const RungeKutta::VelKind velocityType = RungeKutta::VelKind::Rotation;
//const RungeKutta::VelKind velocityType = RungeKutta::VelKind::Vortex;

#include "../includeLS/Utils.hpp"



bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet

  if(!strcmp(SolName, "U")) {  // strcmp compares two string in lexiographic sense.
    //if(facename == 1 || facename == 3) dirichlet = false;
    value = 0.;
  }
  else if(!strcmp(SolName, "V")) {
    // if(facename == 2 || facename == 4) dirichlet = false;
    value = 0.;
//     if(x[0] < 0. && x[1] < 0.5 && x[1] > -0.5 && x[2] < 0.5 && x[2] > -0.5) value = 1.;
  }
  else if(!strcmp(SolName, "W")) {
    value = 0.;
  }
  else if(!strcmp(SolName, "P1") || !strcmp(SolName, "P2") ) {
    dirichlet = false;
    value = 0.;
  }

  return dirichlet;
}



void AssembleMultiphase(MultiLevelProblem& ml_prob);
double TimeStepMultiphase(const double time);
void BestFitLinearInterpolation(std::vector<const double *>& xg, std::vector<double>& psi, std::vector<double> &B);

int main(int argc, char **argv) {

  // Initialize PETSc/MPI
  FemusInit mpinit(argc, argv, MPI_COMM_WORLD);

  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if (nprocs == 1)
    ProfilerStart("profiling.prof");

  int iproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);

  MultiLevelMesh mlMsh0;

  const double scalingFactor = 1.0;
  const unsigned numberOfUniformLevels = 2u;
  const unsigned numberOfSelectiveLevels = 4u;

  std::string meshName = "./input/tri.neu";

  // Load coarse mesh and build uniform refinement levels
  mlMsh0.ReadCoarseMesh(meshName.c_str(), "seventh", scalingFactor);
  mlMsh0.RefineMesh(numberOfUniformLevels, numberOfUniformLevels, nullptr);

  unsigned dim = mlMsh0.GetDimension();

  // Parameters for selective AMR (ball centered at xc with radius r)
  const double r = RADIUS;
  std::vector<double> xc = {XG, YG, ZG};
  xc.resize(dim);

  // Iteratively flag and create new AMR levels
  for (unsigned k = 0; k < numberOfSelectiveLevels; ++k) {
    FlagFinestMeshLevel(mlMsh0, r, xc);
    mlMsh0.AddAMRMeshLevel(false);
  }

  mlMsh0.PrintInfo();
  BBoxToIel bbox(mlMsh0, 0, 3);

  // Define solution on the multilevel mesh
  MultiLevelSolution mlSol0(&mlMsh0);
  std::string psiName = "Psi";
  //std::vector<std::string> dPsiName = {"Psi_x", "Psi_y", "Psi_z"};
  //dPsiName.resize(dim);

  mlSol0.AddSolution(psiName.c_str(), LAGRANGE, SECOND, 0, false);
  //for(unsigned d = 0; d < dim; d++) mlSol0.AddSolution(dPsiName[d].c_str(), LAGRANGE, SECOND, 2);
  //mlSol0.AddSolution("Gamma", LAGRANGE, SECOND);

  std::vector<std::string> vName = {"U", "V", "W"};
  vName.resize(dim);

  for(unsigned d = 0; d < dim; d++) mlSol0.AddSolution(vName[d].c_str(), LAGRANGE, SECOND, 2, false);


  mlSol0.AddSolution("P1",  DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
  mlSol0.AddSolution("P2",  DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);

  //mlSol0.AddSolution("P1", LAGRANGE, FIRST);
  //mlSol0.AddSolution("P2", LAGRANGE, FIRST);

  std::string cName = "C";
  mlSol0.AddSolution(cName.c_str(), DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);

  mlSol0.Initialize("All");

  InitSol(mlSol0, vName, 0, 1.);


  // mlSol0.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  //
  // //mlSol0.FixSolutionAtOnePoint("P1");
  // //mlSol0.FixSolutionAtOnePoint("P2");
  // mlSol0.GenerateBdc("All");

  unsigned sigmoidType = 1;
  double eps = 0.25; //(dim == 2) ? 1. / pow(2, std::max(levelN - 7u, 1u))
  Mollifier m = Mollifier(eps, sigmoidType);

  std::vector<double> xc_1 = xc;
  double r_0 = r;
  Circle c1(xc_1, r_0);

  Boundary zero_bd;

  PsiBall psi2D(xc, r, m);
  InitLevelSet(mlSol0, psiName, psi2D);
  UpdateColorFunction(mlSol0, psiName, cName);


  // Export solution to VTK (selected levels)
  std::vector<std::string> variablesToBePrinted = {"All"};
  VTKWriter vtkIO(&mlSol0);
  vtkIO.SetDebugOutput(true);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);

  MultiLevelMesh mlMsh1;
  MultiLevelSolution mlSol1;

  MultiLevelMesh *mlmsh0 = &mlMsh0;
  MultiLevelMesh *mlmsh1 = &mlMsh1;

  MultiLevelSolution *mlsol0 = &mlSol0;
  MultiLevelSolution *mlsol1 = &mlSol1;

  // Initialize markers object
  LevelSetMarkers markers(psiName, dim);

  // Load coarse mesh and build uniform refinement levels
  mlmsh1->ReadCoarseMesh(meshName.c_str(), "seventh", scalingFactor);
  mlmsh1->RefineMesh(numberOfUniformLevels, numberOfUniformLevels, nullptr);

  double period =
    (velocityType == RungeKutta::VelKind::Vortex) ? 2 : 2.0 * M_PI;
  unsigned nSteps = 320;
  double dt = period / nSteps;

  if (iproc == 0) {

    std::ofstream out("area.dat", std::ios::app);

    if (!out) {
      throw std::runtime_error(
        "computeArea: cannot open output file area.dat");
    }

    out << "#" << std::setw(20) << "Time"
        << std::setw(25) << "Area"
        << "\n";

    out.close();
  }

  for (unsigned t = 1; t <= 0 + 1 * nSteps; t++) {

    double time = t * dt;

    mlsol0->CopySolutionToOldSolution();

    unsigned level = numberOfUniformLevels + numberOfSelectiveLevels - 1u;

    MultiLevelMesh mlMsh2(*mlmsh0, level, "seventh");
    auto msh = mlMsh2.GetLevel(0);
    level = msh->GetLevel();
    msh->SetLevel(0);

    MultiLevelSolution mlSol2(&mlMsh2);

    mlSol2.AddSolution(psiName.c_str(), LAGRANGE, SECOND, 0, false);
    for(unsigned d = 0; d < dim; d++) mlSol2.AddSolution(vName[d].c_str(), LAGRANGE, SECOND, 2);
    mlSol2.AddSolution("P1",  DISCONTINUOUS_POLYNOMIAL, ZERO);
    mlSol2.AddSolution("P2",  DISCONTINUOUS_POLYNOMIAL, ZERO);
    mlSol2.AddSolution(cName.c_str(), DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
    mlSol2.Initialize("All");

    mlSol2.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
    mlSol2.FixSolutionAtOnePoint("P1");
    mlSol2.FixSolutionAtOnePoint("P2");
    mlSol2.GenerateBdc("All");

    auto &Sol0 = (mlsol0->GetSolutionLevel(level))->_Sol;
    auto &Sol2 = (mlSol2.GetSolutionLevel(0))->_Sol;
    auto &Sol0Old = (mlsol0->GetSolutionLevel(level))->_SolOld;
    auto &Sol2Old = (mlSol2.GetSolutionLevel(0))->_SolOld;


    for(unsigned i = 0; i < Sol0.size(); i++) {
      *(Sol2[i]) = *(Sol0[i]);
      if((mlsol0->GetSolutionLevel(level))->GetSolutionTimeOrder(i) == 2) {
        *(Sol2Old[i]) = *(Sol0Old[i]);
      }
    }

    MultiLevelProblem mlProb(&mlSol2);

    // add system Navier-Stokes in mlProb as a Linear Implicit System
    TransientNonlinearImplicitSystem& system = mlProb.add_system < TransientNonlinearImplicitSystem > ("NS");

    // add velocity to system
    for(unsigned d = 0; d < dim; d++) system.AddSolutionToSystemPDE(vName[d].c_str());

    //add pressure
    system.AddSolutionToSystemPDE("P1");
    system.AddSolutionToSystemPDE("P2");

    system.SetSparsityPatternMinimumSize(250);

    // attach the assembling function to system
    system.SetAssembleFunction(AssembleMultiphase);
    system.AttachGetTimeIntervalFunction(TimeStepMultiphase);
    //
    // initilaize and solve the system
    system.init();
    system.SetOuterSolver(PREONLY);
    //
    system.MGsolve();

    VTKWriter vtkIO1(&mlSol2);
    vtkIO1.SetDebugOutput(true);
    if (t % 1 == 0)
      vtkIO1.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, t / 1);

    for(unsigned i = 0; i < Sol0.size(); i++) {
      *(Sol0[i]) = *(Sol2[i]);
      if((mlsol0->GetSolutionLevel(level))->GetSolutionTimeOrder(i) == 2) {
        *(Sol0Old[i]) = *(Sol2Old[i]);
      }
    }
    msh->SetLevel(level);


    //InitSol(*mlsol0, vName, time, period);// TODO with Navier-Stokes solve

    bbox.SetMesh(mlmsh0->GetLevel(0));

    RungeKutta rk(time, dt, period, velocityType);

    //RungeKuttaN rk(time, dt, mlmsh0, bbox);

    unsigned nLevels = numberOfUniformLevels + numberOfSelectiveLevels;
    std::vector<MyVector<double>> X0;
    MyVector<unsigned> X0Iel;
    GetCutElementPoints(*mlsol0, psiName, X0, X0Iel);
    // markers.GetCutElementPoints(*mlsol0, X0, X0Iel);
    //
    // if (t % 10 == 0) {
    //   Reinit reinit(psiName, *mlsol0, eps);
    //
    //   reinit.farFieldReinit(X0);
    //   reinit.interfaceFieldReinit(bbox);
    //   reinit.updateSolution();
    // }



    if (t == 1)
      WritePointsVTK("./output/points.0.vtk", X0);

    RungeKutta4(X0, *mlsol0, bbox, vName, dt);
    //rk.rkForward(X0);

    if (t % 1 == 0)
      WritePointsVTK("./output/points." + std::to_string(t / 1) + ".vtk", X0);


    // std::vector<MyVector<double>> field = X0;
    LevelMarkers l0;
    std::vector<LevelMarkers> lX(nLevels);

    bbox.GetInverseMappingOnCoarseLevel(X0, l0, lX[0]);

    for (unsigned k = 1; k < numberOfUniformLevels; k++) {
      bbox.Project(*mlmsh1, lX[k - 1], lX[k]);
    }

    for (unsigned k = numberOfUniformLevels; k < nLevels; ++k) {
      FlagFinestMeshLevel(*mlmsh1, lX[k - 1].GetElements());
      mlmsh1->AddAMRMeshLevel(
        false); // false -> it does not re-evaluate the AMR flag vector
      bbox.Project(*mlmsh1, lX[k - 1], lX[k]);
    }

    mlsol1->Build(mlmsh1);
    mlsol1->AddSolution(psiName.c_str(), LAGRANGE, SECOND, 0, false);
    for(unsigned d = 0; d < dim; d++) mlsol1->AddSolution(vName[d].c_str(), LAGRANGE, SECOND, 2, false);
    mlsol1->AddSolution("P1", DISCONTINUOUS_POLYNOMIAL, ZERO, false);
    mlsol1->AddSolution("P2", DISCONTINUOUS_POLYNOMIAL, ZERO, false);
    mlsol1->AddSolution(cName.c_str(), DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);

    mlsol1->Initialize("All");

    ProjectSolution(*mlsol0, *mlsol1, bbox, rk, {psiName}, vName, zero_bd, -dt, time, period);
    ProjectSolution(*mlsol0, *mlsol1, bbox, rk, vName);

    UpdateColorFunction(*mlsol1, psiName, cName);

    // Export solution to VTK (selected levels)



    std::swap(mlsol0, mlsol1);
    std::swap(mlmsh0, mlmsh1);

     VTKWriter vtkIO2(mlsol0);
    if (t % 1 == 0)
      vtkIO2.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted,
                   t / 1);

    mlsol1->clear();
    mlmsh1->resize(numberOfUniformLevels);


    // mlsol0->AttachSetBoundaryConditionFunction(SetBoundaryCondition);
    // mlsol0->FixSolutionAtOnePoint("P1");
    // mlsol0->FixSolutionAtOnePoint("P2");
    // mlsol0->GenerateBdc("All");


    double area = ComputeArea(*mlsol0, psiName);

    if (iproc == 0) {

      std::ofstream out("area.dat", std::ios::app);

      if (!out) {
        throw std::runtime_error(
          "computeArea: cannot open output file area.dat");
      }

      out << std::setprecision(16)
          << std::setw(20) << time
          << std::setw(25) << area
          << "\n";

      out.close();
    }
  }

  if (nprocs == 1)
    ProfilerStop();
  return 0;
}

double TimeStepMultiphase(const double time) {
  // double dt =  0.005; //RT
  // double dt =  0.001; //RT
  double dt =  0.015; //Turek
  // double sigma = 3;
  // double rho = 100.;
  // // double totalT = sqrt(rho*0.4*0.4*0.4) / sqrt(sigma);
  // // double dt =  totalT/800; //Parasitic Test
  //
  // double dt =   0.001 * sqrt(rho * 0.4 * 0.4 * 0.4 / sigma);
  // // double dt =  0.0001; //TODO if you use the 320x320 you have to change this
  return dt;
}

//Attempting to create J by hand
void AssembleMultiphase(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled


  //  extract pointers to the several objects that we are going to use
  TransientNonlinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<TransientNonlinearImplicitSystem> ("NS");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual std::vector object in pdeSys (level)

  MatResetPreallocation((static_cast< PetscMatrix* >(KK))->mat());
  MatSetOption((static_cast< PetscMatrix* >(KK))->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  KK->zero();
  RES->zero();

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  if(dim == 2) g = {0, gravity};
  else g = {0, 0, gravity};

  //AssembleGhostPenalty(ml_prob);
  //AssembleGhostPenaltyDGP(ml_prob, true);
  //AssembleGhostPenaltyDGP(ml_prob, false);
  AssembleStabilizationTerms(ml_prob);

  double dt =  mlPdeSys->GetIntervalTime();

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  std::vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("V");    // get the position of "V" in the ml_sol object
  if(dim == 3) solVIndex[2] = mlSol->GetIndex("W");       // get the position of "V" in the ml_sol object

  unsigned solVType = mlSol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

  unsigned solP1Index = mlSol->GetIndex("P1");    // get the position of "P1" in the ml_sol object
  unsigned solP2Index = mlSol->GetIndex("P2");    // get the position of "P2" in the ml_sol object
  unsigned solPType = mlSol->GetSolutionType(solP1Index);    // get the finite element type for "u"

  unsigned psiIndex = mlSol->GetIndex("Psi");
  unsigned psiType = mlSol->GetSolutionType("Psi");
  unsigned cIndex = mlSol->GetIndex("C");

  std::vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("U");    // get the position of "U" in the pdeSys object
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("V");    // get the position of "V" in the pdeSys object
  if(dim == 3) solVPdeIndex[2] = mlPdeSys->GetSolPdeIndex("W");

  unsigned solP1PdeIndex = mlPdeSys->GetSolPdeIndex("P1");    // get the position of "P" in the pdeSys object
  unsigned solP2PdeIndex = mlPdeSys->GetSolPdeIndex("P2");    // get the position of "P" in the pdeSys object

  std::vector < std::vector < double > >  solV(dim);    // local solution
  std::vector < std::vector < double > >  solVOld(dim);    // local solution
  std::vector < double >  solP1; // local solution
  std::vector < double >  solP2; // local solution

  std::vector < double >  psi; // local solution

  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector <double> phiV;  // local test function for velocity
  std::vector <double> phiV_x; // local test function first order partial derivatives

  double* phiP; // local test function for the pressure
  double weight; // gauss point weight

  std::vector< unsigned > sysDof; // local to global pdeSys dofs
  std::vector< double > Res; // local redidual std::vector
  std::vector < double > Jac;

  double eps = 0.00000001;

  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    double C = (*sol->_Sol[cIndex])(iel);

    short unsigned ielGeom = msh->GetElementType(iel);

    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);



    unsigned nDofsVP = dim * nDofsV + 2 * nDofsP;

    // resize local arrays
    sysDof.resize(nDofsVP);
    Res.assign(nDofsVP, 0.);
    Jac.assign(nDofsVP * nDofsVP, 0.);

    for(unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofsV);
      solVOld[k].resize(nDofsV);
      coordX[k].resize(nDofsV);
    }
    solP1.resize(nDofsP);
    solP2.resize(nDofsP);


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
      solP1[i] = (*sol->_Sol[solP1Index])(solPDof);
      solP2[i] = (*sol->_Sol[solP2Index])(solPDof);
      sysDof[dim * nDofsV + i ] = pdeSys->GetSystemDof(solP1Index, solP1PdeIndex, i, iel);
      sysDof[dim * nDofsV + nDofsP + i ] = pdeSys->GetSystemDof(solP2Index, solP2PdeIndex, i, iel);
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);
      for(unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);
      }
    }

    std::vector<double> a;
    std::vector<double> xm;
    double d;
    unsigned cut = 0;
    double vol;

    std::vector<std::vector<double>> Jacob, JacI;

    const elem_type *femV = msh->_finiteElement[ielGeom][solVType];
    const elem_type *femP = msh->_finiteElement[ielGeom][solPType];

    //unsigned cnt =

    cut = (fabs(C - 0.5) < 0.1 ) ? 1 : 0;

    if(cut == 1) {
      femV = fem.GetFiniteElement(ielGeom, solVType);
      femP = fem.GetFiniteElement(ielGeom, solPType);

      unsigned nDofsPsi = msh->GetElementDofNumber(iel, psiType);
      psi.resize(nDofsPsi);
      for(unsigned i = 0; i < nDofsPsi; i++) {
        unsigned psiDof = msh->GetSolutionDof(i, iel, psiType);
        psi[i] = (*sol->_Sol[psiIndex])(psiDof);
      }
      const elem_type *femPsi = fem.GetFiniteElement(ielGeom, psiType);

      unsigned ng = femPsi->GetGaussPointNumber();
      std::vector<double> psig(ng, 0.);
      for(unsigned ig = 0; ig < ng; ig++) {
        double *phi = femPsi->GetPhi(ig);
        for(unsigned i = 0; i < nDofsPsi; i++) {
          psig[ig] += psi[i] * phi[i];
        }
      }
      std::vector<const double *> xg(dim);
      for(unsigned d = 0; d < dim; d++) xg[d] = (femPsi->GetGaussRule()).GetGaussCoordinatePointer(d);

      // femV->GetJacobianMatrix(coordX, {0., 0.}/*cld->GetCloudBaricenterInParentElement(iel) //TODO*/, weight, Jacob, JacI);
      //cld->GetLinearFit(iel, Jacob, a, d); //TODO
      BestFitLinearInterpolation(xg, psig, a);
      d = a[dim];
      a.resize(dim);
    }

    std::vector <TypeIO> weightCF(cfw[ielGeom]->GetGaussQuadraturePointNumber(), 0.);
    std::vector <TypeIO> weightCFInt(cfw[ielGeom]->GetGaussQuadraturePointNumber(), 0.);
    std::vector <TypeIO> weightCFExt(cfw[ielGeom]->GetGaussQuadraturePointNumber(), 0.);

    if(cut == 1) {
      cfw[ielGeom]->GetWeightWithMap(0, a, d, weightCFInt);
      for(unsigned k = 0; k < dim; k++) a[k] = - a[k];
      d = -d;
      cfw[ielGeom]->GetWeightWithMap(-1, a, d, weightCF);
      cfw[ielGeom]->GetWeightWithMap(0, a, d, weightCFExt);
    }
    else {
      for(unsigned i = 0; i < weightCFInt.size(); i++) {
        weightCFInt[i] = C;
        weightCFExt[i] = 1. - C;
      }
    }

    std::vector<double> xg(dim);
    std::vector<double> NN(dim, 0.);
    double kk = 0.;

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < femV->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      femV->Jacobian(coordX, ig, weight, phiV, phiV_x);
      phiP = femP->GetPhi(ig);

      double dsN = 0.;
      std::vector <double> Nf(dim, 0); // unit normal in the physical element from the fluid to the solid

      if(cut == 1) {

        femV->GetJacobianMatrix(coordX, ig, weight, Jacob, JacI);

        for(unsigned k = 0; k < dim; k++) {
          for(unsigned j = 0; j < dim; j++) {
            Nf[k] += JacI[j][k] * a[j];
          }
          dsN += Nf[k] * Nf[k];
        }
        dsN = sqrt(dsN);
        for(unsigned k = 0; k < dim; k++) {
          Nf[k] /= dsN;
        }


        xg.assign(dim, 0);

        for(unsigned i = 0; i < nDofsV; i++) {
          for(unsigned k = 0; k < dim; k++) {
            xg[k] += coordX[k][i] * phiV[i];
          }
        }

        double det = 0.;
        for(unsigned k = 0; k < dim; k++) det += xg[k] * xg[k];
        det = std::sqrt(det);
        for(unsigned k = 0; k < dim; k++) NN[k] = xg[k] / det;

        kk = 1. / RADIUS;


      }





//       if(true ) { // TODO cld->GetNumberOfMarker(iel) > 0) {
//         double magN2 = 0.;
// //         kk = cld->getCurvature(iel, xqp);
//         //kk = cld->GetAverageCurvature(iel); // TODO get curvature
//         //NN = cld->GetNormal(iel, xqp); //TODO get normal
// //       kk = CurvatureQuadric({1., 1., 0., - 2 * XG, - 2 * YG, XG * XG + YG * YG - RADIUS * RADIUS}, xqp);
// //       kk = 1. / RADIUS;
// //       NormalQuadric({1., 1., 0., - 2 * XG, - 2 * YG, XG * XG + YG * YG - RADIUS * RADIUS}, xqp, NN); //TODO
// //       for(unsigned k = 0; k < dim; k++) magN2 += NN[k] * NN[k];
// //       for(unsigned k = 0; k < dim; k++) NN[k] /= sqrt(magN2);
//       }

      std::vector < double > solV_gss(dim, 0);
      std::vector < double > solVOld_gss(dim, 0);
      std::vector < std::vector < double > > gradSolV_gss(dim);

      for(unsigned  k = 0; k < dim; k++) {
        gradSolV_gss[k].assign(dim, 0.);
      }

      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          solV_gss[k] += solV[k][i] * phiV[i];
          solVOld_gss[k] += solVOld[k][i] * phiV[i];
        }
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            gradSolV_gss[k][j] += solV[k][i] * phiV_x[i * dim + j];
          }
        }
      }

      double solP1_gss = 0;
      double solP2_gss = 0;
      for(unsigned i = 0; i < nDofsP; i++) {
        solP1_gss += phiP[i] * solP1[i];
        solP2_gss += phiP[i] * solP2[i];
      }

//       double rho = rho1 * weightCFInt[ig] + rho2 * weightCFExt[ig];
//       double mu = mu1 * weightCFInt[ig] + mu2 * weightCFExt[ig];

      double rho = rho1 * C + rho2 * (1. - C);
      double mu = mu1 * C + mu2 * (1. - C);

      double rhoC = rho1 * C + rho2 * (1. - C);

      // *** phiV_i loop ***
      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned  I = 0; I < dim; I++) {  //momentum equation in k
          double NSV = 0.;
          for(unsigned J = 0; J < dim; J++) {  // second index j in each equation
            NSV   +=  mu * phiV_x[i * dim + J] * (gradSolV_gss[I][J] + gradSolV_gss[J][I]); // diffusion
            NSV   +=  rho * phiV[i] * (solV_gss[J] * gradSolV_gss[I][J]); // nonlinear term
          }
          NSV += - phiV_x[i * dim + I] * (solP1_gss * weightCFInt[ig] + solP2_gss * weightCFExt[ig]);  // pressure gradient
          NSV += rho * phiV[i] * (solV_gss[I] - solVOld_gss[I]) / dt ;
          NSV += - rhoC * phiV[i] * g[I]; // gravity term
          Res[I * nDofsV + i] -=  NSV * weight;
          if(cut == 1) {
            //std::cout << - sigma * phiV[i] * NN[I] * weight * weightCF[ig] * kk * dsN << " ";
            Res[I * nDofsV + i] += - sigma * phiV[i] * NN[I] * weight * weightCF[ig] * kk * dsN;
          }
        }
      } // end phiV_i loop

      // *** phiP_i loop ***
      for(unsigned i = 0; i < nDofsP; i++) {
        for(int I = 0; I < dim; I++) {
          Res[dim * nDofsV + i] += - gradSolV_gss[I][I] * phiP[i]  * weight * weightCFInt[ig]; //continuity
          Res[dim * nDofsV + nDofsP + i] += - gradSolV_gss[I][I] * phiP[i]  * weight * weightCFExt[ig]; //continuity
        }
        if(C == 0)
          Res[dim * nDofsV + i] += - solP1_gss * phiP[i]  * weight * (1 - C) * eps; //penalty
        if(C == 1)
          Res[dim * nDofsV + nDofsP + i] += - solP2_gss * phiP[i]  * weight * C * eps; //penalty

      } // end phiP_i loop
      // end gauss point loop


      //--------------------------------------------------------------------------------------------------------
      // Add the local Matrix/Vector into the global Matrix/Vector

      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned I = 0; I < dim; I++) { //row velocity blocks or dimension
          unsigned VIrow = I * nDofsV + i;
          for(unsigned j = 0; j < nDofsV; j++) {
            unsigned VIcolumn = I * nDofsV + j;

            Jac[ VIrow * nDofsVP + VIcolumn] += rho * phiV[i] * phiV[j] * weight / dt; // inertia

            for(unsigned J = 0; J < dim ; J++) { //column velocity blocks or dimension
              unsigned VJcolumn = J * nDofsV + j;
              Jac[ VIrow * nDofsVP + VIcolumn ] += mu * phiV_x[i * dim + J] * phiV_x[j * dim + J] * weight; //diagonal diffusion
              Jac[ VIrow * nDofsVP + VJcolumn ] += mu * phiV_x[i * dim + J] * phiV_x[j * dim + I] * weight; //off-diagonal diffusion

              Jac[ VIrow * nDofsVP + VIcolumn ] += rho * phiV[i] * solV_gss[J] * phiV_x[j * dim + J] * weight; //diagonal nonlinear
              Jac[ VIrow * nDofsVP + VJcolumn ] += rho * phiV[i] * phiV[j] * gradSolV_gss[I][J] * weight; //off-diagonal nonlinear
            }
          }

          for(unsigned j = 0; j < nDofsP; j++) {
            unsigned P1column = dim * nDofsV + j;
            unsigned P2column = dim * nDofsV + nDofsP + j;
            Jac[VIrow * nDofsVP + P1column] += - phiV_x[i * dim + I] * phiP[j] * weight * weightCFInt[ig]; //pressure gradient
            Jac[P1column * nDofsVP + VIrow] -= - phiV_x[i * dim + I] * phiP[j] * weight * weightCFInt[ig]; //continuity
            Jac[VIrow * nDofsVP + P2column] += - phiV_x[i * dim + I] * phiP[j] * weight * weightCFExt[ig]; //pressure gradient
            Jac[P2column * nDofsVP + VIrow] -= - phiV_x[i * dim + I] * phiP[j] * weight * weightCFExt[ig]; //continuity
          }
        }
      }
      for(unsigned i = 0; i < nDofsP; i++) {
        unsigned P1row = dim * nDofsV + i;
        unsigned P2row = dim * nDofsV + nDofsP + i;
        for(unsigned j = 0; j < nDofsP; j++) {
          unsigned P1column = dim * nDofsV + j;
          unsigned P2column = dim * nDofsV + nDofsP + j;
          if(C == 0)
            Jac[P1row * nDofsVP + P1column] += phiP[i] * phiP[j] * weight * (1 - C) * eps; //continuity
          if(C == 1)
            Jac[P2row * nDofsVP + P2column] += phiP[i] * phiP[j] * weight * C * eps; //continuity
        }
      }


    }

    RES->add_vector_blocked(Res, sysDof);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);


  } //end element loop for each process


  RES->close();
  KK->close();

//  VecView ( (static_cast<PetscVector*> (RES))->vec(),  PETSC_VIEWER_STDOUT_SELF);
//MatView ( (static_cast<PetscMatrix*> (KK))->mat(), PETSC_VIEWER_STDOUT_SELF);

//   PetscViewer    viewer;
//   PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL, NULL, 0, 0, 900, 900, &viewer);
//   PetscObjectSetName((PetscObject) viewer, "PWilmore matrix");
//   PetscViewerPushFormat(viewer, PETSC_VIEWER_DRAW_LG);
//   MatView((static_cast<PetscMatrix*>(KK))->mat(), viewer);
//   double a;
//   std::cin >> a;


}


void AssembleMultiphaseAD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;



  //  extract pointers to the several objects that we are going to use
  TransientNonlinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<TransientNonlinearImplicitSystem> ("NS");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual std::vector object in pdeSys (level)

  MatResetPreallocation((static_cast< PetscMatrix* >(KK))->mat());
  MatSetOption((static_cast< PetscMatrix* >(KK))->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  KK->zero();
  RES->zero();

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  if(dim == 2) g = {0, gravity};
  else g = {0, 0, gravity};

  AssembleGhostPenalty(ml_prob);
  AssembleGhostPenaltyDGP(ml_prob, true);
  AssembleGhostPenaltyDGP(ml_prob, false);
  AssembleStabilizationTerms(ml_prob);

  double dt =  mlPdeSys->GetIntervalTime();
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  std::vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("V");    // get the position of "V" in the ml_sol object
  if(dim == 3) solVIndex[2] = mlSol->GetIndex("W");       // get the position of "V" in the ml_sol object

  unsigned solVType = mlSol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

  unsigned solP1Index = mlSol->GetIndex("P1");    // get the position of "P1" in the ml_sol object
  unsigned solP2Index = mlSol->GetIndex("P2");    // get the position of "P2" in the ml_sol object
  unsigned solPType = mlSol->GetSolutionType(solP1Index);    // get the finite element type for "u"

  unsigned solCIndex = mlSol->GetIndex("C");

  std::vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("U");    // get the position of "U" in the pdeSys object
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("V");    // get the position of "V" in the pdeSys object
  if(dim == 3) solVPdeIndex[2] = mlPdeSys->GetSolPdeIndex("W");

  unsigned solP1PdeIndex = mlPdeSys->GetSolPdeIndex("P1");    // get the position of "P" in the pdeSys object
  unsigned solP2PdeIndex = mlPdeSys->GetSolPdeIndex("P2");    // get the position of "P" in the pdeSys object

  std::vector < std::vector < adept::adouble > >  solV(dim);    // local solution
  std::vector < std::vector < double > >  solVOld(dim);    // local solution
  std::vector < adept::adouble >  solP1; // local solution
  std::vector < adept::adouble >  solP2; // local solution

  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector <double> phiV;  // local test function for velocity
  std::vector <double> phiV_x; // local test function first order partial derivatives

  double* phiP; // local test function for the pressure
  double weight; // gauss point weight

  std::vector< unsigned > sysDof; // local to global pdeSys dofs
  std::vector< adept::adouble > Res; // local redidual std::vector
  std::vector < double > Jac;

  /* BEGIN cutfem stuff for surface tension integration */

  double R = RADIUS;

  std::vector < std::vector < double > > x1;
  std::vector < double > xg(dim);
  xg[0] = XG;
  xg[1] = YG;
  if(dim > 2) xg[2] = ZG;

  unsigned qM = 3;
  double dx = .05;
  double dtetha = 2.;

  double eps = 0.00000001;

  CutFemWeight <TypeIO, TypeA> tet  = CutFemWeight<TypeIO, TypeA >(TET, qM, "legendre");
  CDWeightQUAD <TypeA> quadCD(qM, dx, dtetha);
  CDWeightTRI <TypeA> triCD(qM, dx, dtetha);


  /* END cutfem stuff for surface tension integration */

// cld->AddEllipse({XG, YG}, {RADIUS, RADIUS}, nMax);

//   cld.RKAdvection(4, {"U", "V"}, dtetha); // TODO dtetha sbagliato
//   cld->PrintCSV("markerBefore",it);
//  cld->ComputeQuadraticBestFit();
//   cld->RebuildMarkers(8, 12, 8);
//   cld->PrintCSV("marker",it);

  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

//       for(unsigned iel = msh->_elementOffset[msh->processor_id()]; iel < msh->_elementOffset[msh->processor_id() + 1]; iel++) {
//       std::cout << "iel = " << iel << "   ";
//       const std::vector<double> &a = cld.GetQuadraticBestFitCoefficients(iel);
//       for(unsigned i = 0; i < a.size(); i++) std::cout << a[i] << "  ";
//       std::cout << "\n";
//     }
//     std::cout << std::endl;

    double C = (*sol->_Sol[solCIndex])(iel);

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDof = msh->GetElementDofNumber(iel, 0);  // number of coordinate linear element dofs
    x1.resize(dim);
    for(unsigned k = 0; k < dim; k++) {
      x1[k].resize(nDof);
    }

    for(unsigned k = 0; k < dim; k++) {
      for(unsigned i = 0; i < nDof; i++) {
        unsigned xDof  = msh->GetSolutionDof(i, iel, coordXType);    // global to global mapping between coordinates node and coordinate dof
        x1[k][(i + 2) % nDof] = (*msh->_topology->_Sol[k])(xDof); // global extraction and local storage for the element coordinates
      }
    }

    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);    // number of solution element dofs

    unsigned nDofsVP = dim * nDofsV + 2 * nDofsP;

    // resize local arrays
    sysDof.resize(nDofsVP);
    Res.assign(nDofsVP, 0.);
    Jac.assign(nDofsVP * nDofsVP, 0.);

    for(unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofsV);
      solVOld[k].resize(nDofsV);
      coordX[k].resize(nDofsV);
    }
    solP1.resize(nDofsP);
    solP2.resize(nDofsP);


    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // local to global mapping between solution node and solution dof

      for(unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);
        solVOld[k][i] = (*sol->_SolOld[solVIndex[k]])(solVDof);
        sysDof[k * nDofsV + i] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }

    for(unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);    // local to global mapping between solution node and solution dof
      solP1[i] = (*sol->_Sol[solP1Index])(solPDof);      // global extraction and local storage for the solution
      solP2[i] = (*sol->_Sol[solP2Index])(solPDof);      // global extraction and local storage for the solution
      sysDof[dim * nDofsV + i ] = pdeSys->GetSystemDof(solP1Index, solP1PdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
      sysDof[dim * nDofsV + nDofsP + i ] = pdeSys->GetSystemDof(solP2Index, solP2PdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // local to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }

    std::vector<double> a;
    std::vector<double> xm;
    double d;
    unsigned cut = 0;
    double vol;

    std::vector<std::vector<double>> Jacob, JacI;

    const elem_type *femV = msh->_finiteElement[ielGeom][solVType];
    const elem_type *femP = msh->_finiteElement[ielGeom][solPType];

    unsigned cnt = 0;//TODO cld->GetNumberOfMarker(iel);

    if(cnt > 0) cut = 1;

    if(cut == 1) {
      femV = fem.GetFiniteElement(ielGeom, solVType);
      femP = fem.GetFiniteElement(ielGeom, solPType);
      femV->GetJacobianMatrix(coordX, {0., 0.}/*cld->GetCloudBaricenterInParentElement(iel) TODO*/, weight, Jacob, JacI);
      //cld->GetLinearFit(iel, Jacob, a, d); //TODO
    }

//     if(ielGeom == 3) GetNormalQuad(x1, xg, R, a, d, xm, b, db, cut);
//     if(ielGeom == 3) cld->GetLinearFit(iel, Jacob, b, db);
//     if(ielGeom == 4) GetNormalTri(x1, xg, R, a, d, xm, b, db, cut);
//     else if(ielGeom == 1) GetNormalTetBF(x1, xg, R, a, d, xm, b, db, vol, cut);
//     else if(ielGeom == 0) GetNormalHexBF(x1, xg, R, a, d, xm, b, db, vol, cut, fem.GetFiniteElement(0, 0));

    std::vector <TypeIO> weightCF(quad.GetGaussQuadraturePointNumber(), 0.);
    std::vector <TypeIO> weightCFInt(quad.GetGaussQuadraturePointNumber(), 0.);
    std::vector <TypeIO> weightCFExt(quad.GetGaussQuadraturePointNumber(), 0.);

    if(cut == 1) {
      bool wMap = 1;
      if(ielGeom == 3) {
        quad.GetWeightWithMap(0, a, d, weightCFExt);
        for(unsigned k = 0; k < dim; k++) a[k] = - a[k];
        d = -d;
        quad.GetWeightWithMap(-1, a, d, weightCF);
        quad.GetWeightWithMap(0, a, d, weightCFInt);

//           quadCD.GetWeight(a, d, weightCF);
      }
//         else if(ielGeom == 4) {
//           triCD.GetWeight(b, db, weightCF);
//           const double* weightG = tri.GetGaussWeightPointer();
//         }
//         else if(ielGeom == 1) {
//           tet.GetWeightWithMap(-1, b, db, weightCF);
//           const double* weightG = tet.GetGaussWeightPointer();
//         }
    }
    else {
      for(unsigned i = 0; i < weightCFInt.size(); i++) {
        weightCFInt[i] = C;
        weightCFExt[i] = 1. - C;
      }
    }

    std::vector<double> xqp(dim);
    std::vector<double> NN(dim, 0.);
    double kk = 0.;

    s.new_recording();

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < femV->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      femV->Jacobian(coordX, ig, weight, phiV, phiV_x);
      phiP = femP->GetPhi(ig);

      double dsN = 0.;
      std::vector <double> Nf(dim, 0); // unit normal in the physical element from the fluid to the solid

      if(cut == 1) {

        femV->GetJacobianMatrix(coordX, ig, weight, Jacob, JacI);

        for(unsigned k = 0; k < dim; k++) {
          for(unsigned j = 0; j < dim; j++) {
            Nf[k] += JacI[j][k] * a[j];
          }
          dsN += Nf[k] * Nf[k];
        }
        dsN = sqrt(dsN);
        for(unsigned k = 0; k < dim; k++) {
          Nf[k] /= dsN;
        }
      }



      for(unsigned k = 0; k < dim; k++) {
        xqp[k] = 0.;
      }

      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned k = 0; k < dim; k++) {
          xqp[k] += coordX[k][i] * phiV[i];
        }
      }

      if(true) { //TODO cld->GetNumberOfMarker(iel) > 0) {
        double magN2 = 0.;
//         kk = cld->GetCurvature(iel, xqp);
        //kk = cld->GetAverageCurvature(iel); //TODO get curvature
        // NN = cld->GetNormal(iel, xqp); //TODO get normal
//       kk = CurvatureQuadric({1., 1., 0., - 2 * XG, - 2 * YG, XG * XG + YG * YG - RADIUS * RADIUS}, xqp);
//       kk = 1. / RADIUS;
//       NormalQuadric({1., 1., 0., - 2 * XG, - 2 * YG, XG * XG + YG * YG - RADIUS * RADIUS}, xqp, NN); //TODO
//       for(unsigned k = 0; k < dim; k++) magN2 += NN[k] * NN[k];
//       for(unsigned k = 0; k < dim; k++) NN[k] /= sqrt(magN2);
      }

      std::vector < adept::adouble > solV_gss(dim, 0);
      std::vector < double > solVOld_gss(dim, 0);
      std::vector < std::vector < adept::adouble > > gradSolV_gss(dim);

      for(unsigned  k = 0; k < dim; k++) {
        gradSolV_gss[k].assign(dim, 0.);
      }

      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          solV_gss[k] += solV[k][i] * phiV[i];
          solVOld_gss[k] += solVOld[k][i] * phiV[i];
        }
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            gradSolV_gss[k][j] += solV[k][i] * phiV_x[i * dim + j];
          }
        }
      }

      adept::adouble  solP1_gss = 0;
      adept::adouble  solP2_gss = 0;
      for(unsigned i = 0; i < nDofsP; i++) {
        solP1_gss += phiP[i] * solP1[i];
        solP2_gss += phiP[i] * solP2[i];
      }

//       double rho = rho1 * weightCFInt[ig] + rho1 * weightCFExt[ig];
//       double mu = mu1 * weightCFInt[ig] + mu2 * weightCFExt[ig];

      double rho = rho1 * C + rho2 * (1. - C);
      double mu = mu1 * C + mu2 * (1. - C);

      double rhoC = rho1 * C + rho2 * (1. - C);

      // *** phiV_i loop ***
      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned  I = 0; I < dim; I++) {  //momentum equation in k
          adept::adouble  NSV = 0.;
          for(unsigned J = 0; J < dim; J++) {  // second index j in each equation
            NSV   +=  mu * phiV_x[i * dim + J] * (gradSolV_gss[I][J] + gradSolV_gss[J][I]); // diffusion
            NSV   +=  rho * phiV[i] * (solV_gss[J] * gradSolV_gss[I][J]); // nonlinear term
          }
          NSV += - phiV_x[i * dim + I] * (solP1_gss * weightCFInt[ig] + solP2_gss * weightCFExt[ig]);  // pressure gradient
          NSV += rho * phiV[i] * (solV_gss[I] - solVOld_gss[I]) / dt ;
          NSV += - rhoC * phiV[i] * g[I]; // gravity term
          Res[I * nDofsV + i] +=  NSV * weight;
          if(cut == 1) {
            Res[I * nDofsV + i] -= - sigma * phiV[i] /** b[I]*/ * NN[I] * weight * weightCF[ig] * kk * dsN;
          }
        }
      } // end phiV_i loop

      // *** phiP_i loop ***
      for(unsigned i = 0; i < nDofsP; i++) {
        for(int I = 0; I < dim; I++) {
          Res[dim * nDofsV + i] -= - gradSolV_gss[I][I] * phiP[i]  * weight * weightCFInt[ig]; //continuity
          Res[dim * nDofsV + nDofsP + i] -= - gradSolV_gss[I][I] * phiP[i]  * weight * weightCFExt[ig]; //continuity
        }
        if(C == 0)
          Res[dim * nDofsV + i] -= - solP1_gss * phiP[i]  * weight * (1 - C) * eps; //penalty
        if(C == 1)
          Res[dim * nDofsV + nDofsP + i] -= - solP2_gss * phiP[i]  * weight * C * eps; //penalty

      } // end phiP_i loop
      // end gauss point loop


      //--------------------------------------------------------------------------------------------------------
      // Add the local Matrix/Vector into the global Matrix/Vector

//       for(unsigned i = 0; i < nDofsV; i++) {
//         for(unsigned I = 0; I < dim; I++) { //row velocity blocks or dimension
//           unsigned VIrow = I * nDofsV + i;
//           for(unsigned j = 0; j < nDofsV; j++) {
//             unsigned VIcolumn = I * nDofsV + j;
//
//             Jac[ VIrow * nDofsVP + VIcolumn] += rho * phiV[i] * phiV[j] * weight / dt; // inertia
//
//             for(unsigned J = 0; J < dim ; J++) { //column velocity blocks or dimension
//               unsigned VJcolumn = J * nDofsV + j;
//               Jac[ VIrow * nDofsVP + VIcolumn ] += mu * phiV_x[i * dim + J] * phiV_x[j * dim + J] * weight; //diagonal diffusion
//               Jac[ VIrow * nDofsVP + VJcolumn ] += mu * phiV_x[i * dim + J] * phiV_x[j * dim + I] * weight; //off-diagonal diffusion
//
//               Jac[ VIrow * nDofsVP + VIcolumn ] += rho * phiV[i] * solV_gss[J] * phiV_x[j * dim + J] * weight; //diagonal nonlinear
//               Jac[ VIrow * nDofsVP + VJcolumn ] += rho * phiV[i] * phiV[j] * gradSolV_gss[I][J] * weight; //off-diagonal nonlinear
//             }
//           }
//
//           for(unsigned j = 0; j < nDofsP; j++) {
//             unsigned P1column = dim * nDofsV + j;
//             unsigned P2column = dim * nDofsV + nDofsP + j;
//             Jac[VIrow * nDofsVP + P1column] += - phiV_x[i * dim + I] * phiP[j] * weight * weightCFInt[ig]; //pressure gradient
//             Jac[P1column * nDofsVP + VIrow] -= - phiV_x[i * dim + I] * phiP[j] * weight * weightCFInt[ig]; //continuity
//             Jac[VIrow * nDofsVP + P2column] += - phiV_x[i * dim + I] * phiP[j] * weight * weightCFExt[ig]; //pressure gradient
//             Jac[P2column * nDofsVP + VIrow] -= - phiV_x[i * dim + I] * phiP[j] * weight * weightCFExt[ig]; //continuity
//           }
//         }
//       }
//       for(unsigned i = 0; i < nDofsP; i++) {
//         unsigned P1row = dim * nDofsV + i;
//         unsigned P2row = dim * nDofsV + nDofsP + i;
//         for(unsigned j = 0; j < nDofsP; j++) {
//           unsigned P1column = dim * nDofsV + j;
//           unsigned P2column = dim * nDofsV + nDofsP + j;
//           if(C == 0)
//           Jac[P1row * nDofsVP + P1column] += phiP[i] * phiP[j] * weight * (1 - C) * eps; //continuity
//           if(C == 1)
//           Jac[P2row * nDofsVP + P2column] += phiP[i] * phiP[j] * weight * C * eps; //continuity
//         }
//       }


    }
    std::vector< double > rhs;
    rhs.resize(nDofsVP);   //resize
    for(int i = 0; i < nDofsVP; i++) {
      rhs[ i ] = -Res[i].value();
    }
    RES->add_vector_blocked(rhs, sysDof);

    s.dependent(&Res[0], nDofsVP);


    // define the independent variables J11
    for(unsigned  k = 0; k < dim; k++) {
      s.independent(&solV[k][0], nDofsV);
    }

    s.independent(&solP1[0], nDofsP);
    s.independent(&solP2[0], nDofsP);

    Jac.assign(nDofsVP * nDofsVP, 0);
    // get the and store jacobian matrix (row-major)
    s.jacobian(&Jac[0], true);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);

    s.clear_independents();

    s.clear_dependents(); // for J21 and J22

  } //end element loop for each process


  RES->close();
  KK->close();

//  VecView ( (static_cast<PetscVector*> (RES))->vec(),  PETSC_VIEWER_STDOUT_SELF);
//MatView ( (static_cast<PetscMatrix*> (KK))->mat(), PETSC_VIEWER_STDOUT_SELF);

//   PetscViewer    viewer;
//   PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL, NULL, 0, 0, 900, 900, &viewer);
//   PetscObjectSetName((PetscObject) viewer, "PWilmore matrix");
//   PetscViewerPushFormat(viewer, PETSC_VIEWER_DRAW_LG);
//   MatView((static_cast<PetscMatrix*>(KK))->mat(), viewer);
//   double a;
//   std::cin >> a;

}

void BestFitLinearInterpolation(std::vector<const double*>& xg,
                                std::vector<double>& psi,
                                std::vector<double>& B) {

  const unsigned n = psi.size();
  const unsigned dim = xg.size();
  const unsigned m = dim + 1;

  double s2 = 0.0;

  for (unsigned i = 0; i < n; i++) {
    s2 += psi[i] * psi[i];
  }

  s2 /= n;

  if (s2 < 1.e-14) {
    throw std::runtime_error("uniform  zero level-set in cutcell");
    return;
  }

  std::vector<std::vector<double>> M(m, std::vector<double>(m, 0.0));
  std::vector<double> F(m, 0.0);

  for (unsigned i = 0; i < n; i++) {

    const double f = psi[i];

    const double w = std::exp(- 100 * f * f / s2);

    for (unsigned d = 0; d < dim; d++) {

      F[d] += w * xg[d][i] * f;

      for (unsigned e = 0; e < dim; e++) {
        M[d][e] += w * xg[d][i] * xg[e][i];
      }

      M[d][dim] += w * xg[d][i];
      M[dim][d] += w * xg[d][i];
    }

    M[dim][dim] += w;
    F[dim] += w * f;
  }

  B = F;

  for (unsigned k = 0; k < m; k++) {

    unsigned pivot = k;

    for (unsigned i = k + 1; i < m; i++) {
      if (std::fabs(M[i][k]) > std::fabs(M[pivot][k])) {
        pivot = i;
      }
    }

    if (std::fabs(M[pivot][k]) < 1.e-14) {
      B.assign(m, 0.0);
      return;
    }

    if (pivot != k) {
      std::swap(M[k], M[pivot]);
      std::swap(B[k], B[pivot]);
    }

    for (unsigned i = k + 1; i < m; i++) {

      const double factor = M[i][k] / M[k][k];

      for (unsigned j = k; j < m; j++) {
        M[i][j] -= factor * M[k][j];
      }

      B[i] -= factor * B[k];
    }
  }

  for (int i = static_cast<int>(m) - 1; i >= 0; i--) {

    for (unsigned j = i + 1; j < m; j++) {
      B[i] -= M[i][j] * B[j];
    }

    B[i] /= M[i][i];
  }

  double norm = 0.0;

  for (unsigned d = 0; d < dim; d++) {
    norm += B[d] * B[d];
  }

  norm = std::sqrt(norm);

  if (norm < 1.e-14) {
    B.assign(m, 0.0);
    return;
  }

  for (unsigned d = 0; d < m; d++) {
    B[d] /= norm;
  }

  double min_psi = psi[0];
  double max_psi = psi[0];

  unsigned i_min = 0;
  unsigned i_max = 0;

  for (unsigned i = 1; i < n; i++) {

    if (psi[i] > max_psi) {
      max_psi = psi[i];
      i_max = i;
    }

    if (psi[i] < min_psi) {
      min_psi = psi[i];
      i_min = i;
    }
  }

  double test_value_max = B[dim];
  double test_value_min = B[dim];

  for (unsigned d = 0; d < dim; d++) {
    test_value_max += B[d] * xg[d][i_max];
    test_value_min += B[d] * xg[d][i_min];
  }

  const bool maxWrong = (max_psi * test_value_max < 0.0);
  const bool minWrong = (min_psi * test_value_min < 0.0);

  if (maxWrong && minWrong) {

    for (unsigned d = 0; d < dim + 1; d++) {
      B[d] *= -1.0;
    }

  }
  else if (maxWrong || minWrong) {

    std::cout << "Warning: incoherent linear approximation" << std::endl;
  }

}



