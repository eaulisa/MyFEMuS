
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
  else if(!strcmp(SolName, "NX") || !strcmp(SolName, "NY") || !strcmp(SolName, "NZ") || !strcmp(SolName, "K") ) {
    dirichlet = false;
    value = 0.;
  }

  return dirichlet;
}


//Global object to be used in the assembly as an external ml_problem
MultiLevelProblem* ml_prob0;
unsigned levelF, levelC;


#include "../includeLS/Stabilization.hpp"
#include "../includeLS/GhostPenalty.hpp"
#include "../includeLS/GhostPenaltyDGP.hpp"

void BestFitLinearInterpolation(std::vector<const double*>& xg,
                                std::vector<double>& psi,
                                std::vector<double>& B);

void AssembleMultiphase(MultiLevelProblem& ml_prob);
void AssembleNormal(MultiLevelProblem& ml_prob);
void AssembleCurvature(MultiLevelProblem& ml_prob);
double TimeStepMultiphase(const double time);


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
  for(unsigned d = 0; d < dim; d++) mlSol0.AddSolution(vName[d].c_str(), LAGRANGE, SECOND, 2);

  std::vector<std::string> pName = {"P1", "P2"};
  for(unsigned d = 0; d < pName.size(); d++) mlSol0.AddSolution(pName[d].c_str(), DISCONTINUOUS_POLYNOMIAL, ZERO, 0);

  //mlSol0.AddSolution("P1", LAGRANGE, FIRST);
  //mlSol0.AddSolution("P2", LAGRANGE, FIRST);

  std::string cName = "C";
  mlSol0.AddSolution(cName.c_str(), DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);

  mlSol0.AddSolution("NX", LAGRANGE, SECOND, 0);
  mlSol0.AddSolution("NY", LAGRANGE, SECOND, 0);
  if (dim == 3) mlSol0.AddSolution("NZ", LAGRANGE, SECOND, 0);

  mlSol0.AddSolution("K", LAGRANGE, SECOND, 0);

  mlSol0.Initialize("All");

  InitSol(mlSol0, vName, 0, 1.);


  // mlSol0.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

  //mlSol0.FixSolutionAtOnePoint("P1");
  //mlSol0.FixSolutionAtOnePoint("P2");
  // mlSol0.GenerateBdc("All");

  unsigned sigmoidType = 0;
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
  unsigned nSteps = 640;
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

  unsigned levelN = numberOfUniformLevels + numberOfSelectiveLevels;
  levelF = levelN - 1u; //fine level associated for mlmsh0 and mlmsh1
  levelC = levelN - 1u; //coarse level associated to mlmsh2, but existing also mlmsh0 and mlmsh1


  for (unsigned t = 1; t <= 0 + 1 * nSteps; t++) {

    double time = t * dt;

    mlsol0->CopySolutionToOldSolution();

    UpdateColorFunction(*mlsol0, psiName, cName);
    if(levelC < levelF) RestrictPWDCField(*mlsol0, cName, levelC, levelF);

    // mlProb0 is used to assemble and solve for N and K and assemble only the Navier-Stokes block and its stabilization terms
    mlsol0->AttachSetBoundaryConditionFunction(SetBoundaryCondition);
    mlsol0->GenerateBdc("All");
    MultiLevelProblem mlProb0(mlsol0);

    // add system Normal in mlProb as a Linear Implicit System
    LinearImplicitSystem& system0_N = mlProb0.add_system < LinearImplicitSystem > ("N");
    system0_N.AddSolutionToSystemPDE("NX");
    system0_N.AddSolutionToSystemPDE("NY");
    if (dim == 3) system0_N.AddSolutionToSystemPDE("NZ");
    system0_N.SetAssembleFunction(AssembleNormal);
    // initilaize and solve the system
    system0_N.SetMgType(V_CYCLE);
    system0_N.init();
    // system0_N.SetOuterSolver(PREONLY);
    system0_N.MGsolve();

    LinearImplicitSystem& system0_K = mlProb0.add_system < LinearImplicitSystem > ("K");
    system0_K.AddSolutionToSystemPDE("K");
    // system0_K.SetSparsityPatternMinimumSize(250);
    system0_K.SetAssembleFunction(AssembleCurvature);
    // initilaize and solve the system
    system0_K.SetMgType(V_CYCLE);
    system0_K.init();
    // system0_K.SetOuterSolver(PREONLY);
    system0_K.MGsolve();

    // mlProb2 is used to assemble the Ghost Penalty and solve the full Navier-Stokes + Ghost penalty

    MultiLevelMesh mlMsh2(*mlmsh0, levelC, levelC + 1, "seventh");
    auto msh = mlMsh2.GetLevel(0);

    MultiLevelSolution mlSol2(&mlMsh2);
    mlSol2.AddSolution(psiName.c_str(), LAGRANGE, SECOND, 0, false);
    for(unsigned d = 0; d < dim; d++) mlSol2.AddSolution(vName[d].c_str(), LAGRANGE, SECOND, 2);
    mlSol2.AddSolution("P1",  DISCONTINUOUS_POLYNOMIAL, ZERO, 0); //TODO
    mlSol2.AddSolution("P2",  DISCONTINUOUS_POLYNOMIAL, ZERO, 0);
    mlSol2.AddSolution(cName.c_str(), DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
    mlSol2.AddSolution("NX", LAGRANGE, SECOND, 0); //TODO
    mlSol2.AddSolution("NY", LAGRANGE, SECOND, 0);
    if (dim == 3) mlSol2.AddSolution("NZ", LAGRANGE, SECOND, 0);
    mlSol2.AddSolution("K", LAGRANGE, SECOND, 0);

    mlSol2.Initialize("All");

    mlSol2.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
    mlSol2.FixSolutionAtOnePoint("P1"); //TODO
    mlSol2.FixSolutionAtOnePoint("P2");
    mlSol2.GenerateBdc("All");


    auto &Sol0_C = (mlsol0->GetSolutionLevel(levelC))->_Sol;
    auto &Sol0Old_C = (mlsol0->GetSolutionLevel(levelC))->_SolOld;

    auto &Sol2 = (mlSol2.GetSolutionLevel(0))->_Sol; //careful here sol2(0) corresponds to the same level of sol0(l0)
    auto &Sol2Old = (mlSol2.GetSolutionLevel(0))->_SolOld;

    for(unsigned i = 0; i < Sol0_C.size(); i++) {
      //copy velocity
      for(unsigned d = 0; d < dim; d++) {
        unsigned vel0index = mlsol0->GetIndex(vName[d].c_str());
        unsigned vel2index = mlSol2.GetIndex(vName[d].c_str());

        *(Sol2[vel2index]) = *(Sol0_C[vel0index]);
        if((mlsol0->GetSolutionLevel(levelC))->GetSolutionTimeOrder(vel0index) == 2) {
          *(Sol2Old[vel2index]) = *(Sol0Old_C[vel0index]);
        }
      }

      unsigned c0index = mlsol0->GetIndex(cName.c_str());
      unsigned c2index = mlSol2.GetIndex(cName.c_str());
      *(Sol2[c2index]) = *(Sol0_C[c0index]);

      for(unsigned d = 0; d < pName.size(); d++) {
        unsigned Pindex = mlSol2.GetIndex(pName[d].c_str());
        Sol2[Pindex]->zero();
      }
    }

    MultiLevelProblem mlProb2(&mlSol2);

    // add system Navier-Stokes in mlProb as a Linear Implicit System
    TransientNonlinearImplicitSystem& system2 = mlProb2.add_system < TransientNonlinearImplicitSystem > ("NS");
    // add velocity to system
    for(unsigned d = 0; d < dim; d++) system2.AddSolutionToSystemPDE(vName[d].c_str());
    //add pressure
    system2.AddSolutionToSystemPDE("P1");
    system2.AddSolutionToSystemPDE("P2");
    system2.SetSparsityPatternMinimumSize(250);
    // attach the assembling function to system
    system2.SetAssembleFunction(AssembleMultiphase);
    system2.AttachGetTimeIntervalFunction(TimeStepMultiphase);
    // initilaize and solve the system
    system2.init();
    system2.SetOuterSolver(PREONLY);

    // add system Navier-Stokes in mlProb as a Linear Implicit System
    TransientNonlinearImplicitSystem& system0 = mlProb0.add_system < TransientNonlinearImplicitSystem > ("NS");
    // add velocity to system
    for(unsigned d = 0; d < dim; d++) system0.AddSolutionToSystemPDE(vName[d].c_str());
    //add pressure
    system0.AddSolutionToSystemPDE("P1");
    system0.AddSolutionToSystemPDE("P2");
    system0.SetSparsityPatternMinimumSize(250);
    // attach the assembling function to system
    system0.SetAssembleFunction(AssembleMultiphase);
    system0.AttachGetTimeIntervalFunction(TimeStepMultiphase);
    // initilaize and solve the system
    system0.init();
    system0.SetOuterSolver(PREONLY);

    ml_prob0 = &mlProb0;

    msh->SetLevel(0);
    system2.MGsolve();

    VTKWriter vtkIO1(&mlSol2);
    vtkIO1.SetDebugOutput(true);
    if (t % 1 == 0)
      vtkIO1.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, t / 1);

    msh->SetLevel(levelC);


    for(unsigned i = 0; i < Sol0_C.size(); i++) {
      //copy velocity
      for(unsigned d = 0; d < dim; d++) {
        unsigned vel0index = mlsol0->GetIndex(vName[d].c_str());
        unsigned vel2index = mlSol2.GetIndex(vName[d].c_str());

        *(Sol0_C[vel0index]) = *(Sol2[vel2index]);
        if((mlsol0->GetSolutionLevel(levelC))->GetSolutionTimeOrder(vel0index) == 2) {
          *(Sol0Old_C[vel0index]) = *(Sol2Old[vel2index]);
        }
      }

      for(unsigned d = 0; d < pName.size(); d++) {
        unsigned P0index = mlsol0->GetIndex(pName[d].c_str());
        unsigned P2index = mlSol2.GetIndex(pName[d].c_str());
        *(Sol0_C[P0index]) = *(Sol2[P2index]);
      }
    }

    bbox.SetMesh(mlmsh0->GetLevel(0));

    RungeKutta rk(time, dt, period, velocityType); // TODO to remove all together

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

    RungeKutta4(X0, *mlsol0, bbox, vName, levelC, dt); // move the interface points forward in time using the velocity mls0(lC)

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
    for(unsigned d = 0; d < dim; d++) mlsol1->AddSolution(vName[d].c_str(), LAGRANGE, SECOND, 2);
    mlsol1->AddSolution("P1", DISCONTINUOUS_POLYNOMIAL, ZERO, 0);
    mlsol1->AddSolution("P2", DISCONTINUOUS_POLYNOMIAL, ZERO, 0);
    mlsol1->AddSolution(cName.c_str(), DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
    mlsol1->AddSolution("NX", LAGRANGE, SECOND, 0);
    mlsol1->AddSolution("NY", LAGRANGE, SECOND, 0);
    if (dim == 3) mlsol1->AddSolution("NZ", LAGRANGE, SECOND, 0);
    mlsol1->AddSolution("K", LAGRANGE, SECOND, 0);


    mlsol1->Initialize("All");

    ProjectSolution(*mlsol0, *mlsol1, bbox, rk, {psiName}, levelF, levelF, vName, levelC, zero_bd, -dt, time, period);
    ProjectSolution(*mlsol0, *mlsol1, bbox, rk, vName, levelC, levelC);


    std::swap(mlsol0, mlsol1);
    std::swap(mlmsh0, mlmsh1);

    // Export solution to VTK (selected levels)
    VTKWriter vtkIO2(mlsol0);
    if (t % 1 == 0) {
      vtkIO2.Write(levelF, DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted,
                   t / 1);
      vtkIO2.Write(levelC, DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted,
                   t / 1);
    }

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

void AssembleNormal(MultiLevelProblem& ml_prob) {
  LinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<LinearImplicitSystem> ("N");
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh* msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem* el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix* KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* RES = pdeSys->_RES; // pointer to the global residual std::vector object in pdeSys (level)

  MatSetOption((static_cast< PetscMatrix* >(KK))->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  unsigned psiIndex = mlSol->GetIndex("Psi");
  unsigned psiType = mlSol->GetSolutionType("Psi");

  std::vector < unsigned > solNIndex(dim);
  solNIndex[0] = mlSol->GetIndex("NX");    // get the position of "U" in the ml_sol object
  solNIndex[1] = mlSol->GetIndex("NY");    // get the position of "V" in the ml_sol object
  if(dim == 3) solNIndex[2] = mlSol->GetIndex("NZ");       // get the position of "V" in the ml_sol object

  std::vector < unsigned > solNPdeIndex(dim);
  solNPdeIndex[0] = mlPdeSys->GetSolPdeIndex("NX");    // get the position of "U" in the pdeSys object
  solNPdeIndex[1] = mlPdeSys->GetSolPdeIndex("NY");    // get the position of "V" in the pdeSys object
  if(dim == 3) solNPdeIndex[2] = mlPdeSys->GetSolPdeIndex("NZ");

  unsigned solNType = mlSol->GetSolutionType(solNIndex[0]);

  std::vector < double >  psi; // local solution

  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector <double> phiN;  // local test function for velocity
  std::vector <double> phiN_x; // local test function first order partial derivatives

  std::vector <double> phiPsi;
  std::vector <double>  phiPsi_x;

  std::vector<std::vector<double>> N(dim);
  double weight; // gauss point weight
  double weightPsi;

  std::vector< unsigned > sysDof; // local to global pdeSys dofs
  std::vector< double > Res; // local redidual std::vector
  std::vector < double > Jac;

  KK->zero();
  RES->zero();

  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);

    unsigned nDofsN = msh->GetElementDofNumber(iel, solNType);
    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType);

    unsigned nDofs =  dim * nDofsN;

    // resize local arrays
    sysDof.resize(nDofs);
    Res.assign(nDofs, 0.);
    Jac.assign(nDofs * nDofs, 0.);

    for(unsigned  k = 0; k < dim; k++) {
      coordX[k].resize(nDofsX);
      N[k].resize(nDofsN);
    }

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsN; i++) {
      unsigned solNDof = msh->GetSolutionDof(i, iel, solNType);
      for(unsigned  d = 0; d < dim; d++) {
        N[d][i] = (*sol->_Sol[solNIndex[d]])(solNDof);
        sysDof[d * nDofsN + i] = pdeSys->GetSystemDof(solNIndex[d], solNPdeIndex[d], i, iel);
      }
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);
      for(unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);
      }
    }

    unsigned nDofsPsi;

    nDofsPsi = msh->GetElementDofNumber(iel, psiType);
    psi.resize(nDofsPsi);
    for(unsigned i = 0; i < nDofsPsi; i++) {
      unsigned psiDof = msh->GetSolutionDof(i, iel, psiType);
      psi[i] = (*sol->_Sol[psiIndex])(psiDof);
    }

    const elem_type *femPsi = msh->_finiteElement[ielGeom][psiType];
    const elem_type *femN = msh->_finiteElement[ielGeom][solNType];

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < femN->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      femN->Jacobian(coordX, ig, weight, phiN, phiN_x);
      femPsi->Jacobian(coordX, ig, weightPsi, phiPsi, phiPsi_x);


      std::vector<double> NN(dim, 0.);
      for (unsigned i = 0; i < nDofsPsi; i++) {
        for(unsigned d = 0; d < dim; d++) {
          NN[d] -= psi[i] * phiPsi_x[i * dim + d];
        }
      }
      double det = 0;
      for (unsigned d = 0; d < dim; d++) {
        det += NN[d] * NN[d];
      }
      det = std::sqrt(det + 1.e-10);
      for (unsigned d = 0; d < dim; d++) {
        NN[d] /= det;
      }

      std::vector<double> N_g(dim, 0.);
      for (unsigned i = 0; i < nDofsN; i++) {
        for(unsigned d = 0; d < dim; d++) {
          N_g[d] += N[d][i] * phiN[i];
        }
      }

      // *** phiV_i loop ***
      for(unsigned i = 0; i < nDofsN; i++) {
        for(unsigned  d = 0; d < dim; d++) {  //momentum equation in k
          double rhs = 0.;
          rhs += phiN[i] * (NN[d] - N_g[d]);
          Res[d * nDofsN + i] +=  rhs * weight;
        }
      } // end phiV_i loop


      //--------------------------------------------------------------------------------------------------------
      // Add the local Matrix/Vector into the global Matrix/Vector

      for(unsigned i = 0; i < nDofsN; i++) {
        for(unsigned d = 0; d < dim; d++) { //row velocity blocks or dimension
          unsigned VIrow = d * nDofsN + i;
          for(unsigned j = 0; j < nDofsN; j++) {
            unsigned VIcolumn = d * nDofsN + j;

            Jac[ VIrow * nDofs + VIcolumn] += phiN[i] * phiN[j] * weight ; // inertia


          }
        }
      }
    }

    RES->add_vector_blocked(Res, sysDof);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);


  } //end element loop for each process

  RES->close();
  KK->close();

}

void AssembleCurvature(MultiLevelProblem& ml_prob) {

  LinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<LinearImplicitSystem> ("K");
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh* msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem* el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix* KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* RES = pdeSys->_RES; // pointer to the global residual std::vector object in pdeSys (level)

  MatSetOption((static_cast< PetscMatrix* >(KK))->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  std::vector<unsigned> solNIndex(dim);
  solNIndex[0] = mlSol->GetIndex("NX");
  solNIndex[1] = mlSol->GetIndex("NY");
  if(dim == 3) solNIndex[2] = mlSol->GetIndex("NZ");
  unsigned solNType = mlSol->GetSolutionType("NX");

  unsigned  solKIndex;
  solKIndex = mlSol->GetIndex("K");    // get the position of "U" in the ml_sol object

  unsigned  solKPdeIndex;
  solKPdeIndex = mlPdeSys->GetSolPdeIndex("K");    // get the position of "U" in the pdeSys object


  unsigned solKType = mlSol->GetSolutionType(solKIndex);

  // std::vector < double >  psi; // local solution

  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector < std::vector < double > > normal(dim);
  std::vector < double >  K;

  std::vector <double> phi;  // local test function for velocity
  std::vector <double> phi_x; // local test function first order partial derivatives
  std::vector <double> bdphi;  // local test function for velocity
  std::vector <double> bdphi_x;

  std::vector <double> phiN;
  std::vector <double> phiN_x;
  std::vector <double> bdphiN;  // local test function for velocity
  std::vector <double> bdphiN_x;

  std::vector < double> normal_face;
  std::vector < double> normal_faceN;
  double weight_face = 0.;
  double weight_faceN = 0.;


  double weight; // gauss point weight
  double weightN;

  std::vector< unsigned > sysDof; // local to global pdeSys dofs
  std::vector< double > Res; // local redidual std::vector
  std::vector < double > Jac;

  KK->zero();
  RES->zero();

  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);

    unsigned nDofs = msh->GetElementDofNumber(iel, solKType);
    unsigned nDofsN = msh->GetElementDofNumber(iel, solNType);
    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType);

    // resize local arrays
    sysDof.resize(nDofs);
    Res.assign(nDofs, 0.);
    Jac.assign(nDofs * nDofs, 0.);

    K.resize(nDofs);
    for(unsigned  d = 0; d < dim; d++) {
      normal[d].resize(nDofsN);
      coordX[d].resize(nDofsX);
    }


    for(unsigned i = 0; i < nDofs; i++) {
      unsigned KDof  = msh->GetSolutionDof(i, iel, solKType);
      K[i] = (*sol->_Sol[solKIndex])(KDof);
      sysDof[i] = pdeSys->GetSystemDof(solKIndex, solKPdeIndex, i, iel);
    }

    for(unsigned i = 0; i < nDofsN; i++) {
      unsigned normalDof = msh->GetSolutionDof(i, iel, solNType);
      for (unsigned d = 0; d < dim; d++) {
        normal[d][i] = (*sol->_Sol[solNIndex[d]])(normalDof);
      }
    }

    for(unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);
      for(unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);
      }
    }

    const elem_type *femK = msh->_finiteElement[ielGeom][solKType];
    const elem_type *femN = msh->_finiteElement[ielGeom][solNType];

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < femN->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      femK->Jacobian(coordX, ig, weight, phi, phi_x);
      femN->Jacobian(coordX, ig, weightN, phiN, phiN_x);

      double K_g = 0.;
      for (unsigned j = 0; j < nDofs; j++) {
        K_g += K[j] * phi[j];
      }

      std::vector<double> normal_g(dim, 0.);
      for(unsigned d = 0; d < dim; d++) {
        for (unsigned j = 0; j < nDofsN; j++) {
          normal_g[d] += normal[d][j] * phiN[j];
        }
      }

      // *** phiV_i loop ***
      for(unsigned i = 0; i < nDofs; i++) {
        double rhs = 0.;
        for(unsigned  d = 0; d < dim; d++) {  //momentum equation in k
          rhs -= phi_x[i * dim + d] * normal_g[d];
        }
        rhs -= K_g * phi[i];
        Res[i] += rhs * weight;
      } // end phiV_i loop


      //--------------------------------------------------------------------------------------------------------
      // Add the local Matrix/Vector into the global Matrix/Vector

      for(unsigned i = 0; i < nDofs; i++) {
        // for(unsigned I = 0; I < dim; I++) { //row velocity blocks or dimension
        unsigned VIrow = i * nDofs;
        for(unsigned j = 0; j < nDofs; j++) {
          unsigned VIcolumn = j;

          Jac[ VIrow + VIcolumn] += phi[i] * phi[j] * weight ; // inertia


        }
        // }
      }
    }

    // *** Face Gauss point loop (boundary Integral) ***
    for ( unsigned jface = 0; jface < msh->GetElementFaceNumber ( iel ); jface++ ) {
      int faceIndex = el->GetBoundaryIndex(iel, jface);
      // look for boundary faces

      if ( faceIndex > 0 ) {
        const unsigned faceGeom = msh->GetElementFaceType ( iel, jface );
        unsigned faceDofs = msh->GetElementFaceDofNumber (iel, jface, solKType);
        unsigned faceDofsN = msh->GetElementFaceDofNumber (iel, jface, solNType);
        unsigned faceDofsX = msh->GetElementFaceDofNumber (iel, jface, coordXType);
        std::vector  < std::vector  <  double> > faceCoordinates ( dim ); // A matrix holding the face coordinates rowwise.
        for ( int k = 0; k < dim; k++ ) {
          faceCoordinates[k].resize (faceDofsX);
        }
        for ( unsigned i = 0; i < faceDofsX; i++ ) {
          unsigned inode = msh->GetLocalFaceVertexIndex ( iel, jface, i ); // face-to-element local node mapping.
          for ( unsigned k = 0; k < dim; k++ ) {
            faceCoordinates[k][i] =  coordX[k][inode]; // We extract the local coordinates on the face from local coordinates on the element.
          }
        }
        for ( unsigned ig = 0; ig  <  msh->_finiteElement[faceGeom][solKType]->GetGaussPointNumber(); ig++ ) {
          // We call the method GetGaussPointNumber from the object finiteElement in the mesh object msh.

          msh->_finiteElement[faceGeom][solKType]->JacobianSur ( faceCoordinates, ig, weight_face, bdphi, bdphi_x, normal_face );
          msh->_finiteElement[faceGeom][solNType]->JacobianSur ( faceCoordinates, ig, weight_faceN, bdphiN, bdphiN_x, normal_faceN );

          std::vector<double> normal_g(dim, 0.);
          for(unsigned d = 0; d < dim; d++) {
            for (unsigned j = 0; j < faceDofsN; j++) {
              unsigned jnode = msh->GetLocalFaceVertexIndex (iel, jface, j );
              normal_g[d] += normal[d][jnode] * bdphiN[j];
            }
          }

          // *** phi_i loop ***
          for ( unsigned i = 0; i < faceDofs; i++ ) {
            double rhs_bd = 0;
            unsigned inode = msh->GetLocalFaceVertexIndex (iel, jface, i );
            for( unsigned d = 0; d < dim; d++) {
              rhs_bd +=  bdphi[i] * normal_face[d] * normal_g[d];
            }
            Res[inode] += rhs_bd * weight_face;
          }
        }
      }
    }

    RES->add_vector_blocked(Res, sysDof);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);


  } //end element loop for each process

  RES->close();
  KK->close();

}


//Attempting to create J by hand
void AssembleMultiphase(MultiLevelProblem& ml_prob2) {

  //  extract pointers to the several objects that we are going to use
  TransientNonlinearImplicitSystem* mlPdeSys2   = &ml_prob2.get_system<TransientNonlinearImplicitSystem> ("NS");
  const unsigned level2 = mlPdeSys2->GetLevelToAssemble();

  Mesh* msh2 = ml_prob2._ml_msh->GetLevel(level2);    // pointer to the mesh (level) object
  elem* el2 = msh2->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol2        = ml_prob2._ml_sol;  // pointer to the multilevel solution object
  Solution* sol2 = ml_prob2._ml_sol->GetSolutionLevel(level2);    // pointer to the solution (level) object


  LinearEquationSolver* pdeSys2        = mlPdeSys2->_LinSolver[level2]; // pointer to the equation (level) object
  SparseMatrix* KK2 = pdeSys2->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* RES2 = pdeSys2->_RES; // pointer to the global residual std::vector object in pdeSys (level)

  //MatResetPreallocation((static_cast< PetscMatrix* >(KK2))->mat());
  MatSetOption((static_cast< PetscMatrix* >(KK2))->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  KK2->zero();
  RES2->zero();

  AssembleGhostPenalty(ml_prob2);
  AssembleGhostPenaltyDGP(ml_prob2, true);
  AssembleGhostPenaltyDGP(ml_prob2, false);

  RES2->close();
  KK2->close();


  //  extract pointers to the several objects that we are going to use
  TransientNonlinearImplicitSystem* mlPdeSys   = &ml_prob0->get_system<TransientNonlinearImplicitSystem> ("NS");


  Mesh*          msh          = ml_prob0->_ml_msh->GetLevel(levelF);    // pointer to the mesh (levelF) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (levelF)

  MultiLevelSolution*  mlSol        = ml_prob0->_ml_sol;  // pointer to the multilevelF solution object
  Solution*    sol        = ml_prob0->_ml_sol->GetSolutionLevel(levelF);    // pointer to the solution (levelF) object


  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[levelF]; // pointer to the equation (levelF) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (levelF)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual std::vector object in pdeSys (levelF)


  //MatResetPreallocation((static_cast< PetscMatrix* >(KK))->mat());
  MatSetOption((static_cast< PetscMatrix* >(KK))->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);



  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  if(dim == 2) g = {0, gravity};
  else g = {0, 0, gravity};

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

  unsigned cIndex = mlSol->GetIndex("C");

  std::vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("U");    // get the position of "U" in the pdeSys object
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("V");    // get the position of "V" in the pdeSys object
  if(dim == 3) solVPdeIndex[2] = mlPdeSys->GetSolPdeIndex("W");

  unsigned solP1PdeIndex = mlPdeSys->GetSolPdeIndex("P1");    // get the position of "P" in the pdeSys object
  unsigned solP2PdeIndex = mlPdeSys->GetSolPdeIndex("P2");    // get the position of "P" in the pdeSys object

  unsigned psiIndex = mlSol->GetIndex("Psi");
  unsigned psiType = mlSol->GetSolutionType("Psi");

  std::vector<unsigned> solNIndex(dim);
  solNIndex[0] = mlSol->GetIndex("NX");
  solNIndex[1] = mlSol->GetIndex("NY");
  if(dim == 3) solNIndex[2] = mlSol->GetIndex("NZ");
  unsigned solNType = mlSol->GetSolutionType("NX");

  unsigned solKIndex = mlSol->GetIndex("K");
  unsigned solKType = mlSol->GetSolutionType("K");

  std::vector < std::vector < double > >  solV(dim);    // local solution
  std::vector < std::vector < double > >  solVOld(dim);    // local solution
  std::vector < double >  solP1; // local solution
  std::vector < double >  solP2; // local solution

  std::vector < double >  psi; // local solution
  std::vector < double >  k;
  std::vector < std::vector < double > >  n(dim);

  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector <double> phiV;  // local test function for velocity
  std::vector <double> phiV_x; // local test function first order partial derivatives

  std::vector <double> phiPsi;
  std::vector <double>  phiPsi_x;
  std::vector <double> phiPsi_xx;

  std::vector <double> phiN;
  std::vector <double> phiN_x;

  std::vector <double> phiK;
  std::vector <double> phiK_x;

  unsigned dim2 = 3 * (dim - 1);

  double* phiP; // local test function for the pressure
  double weight; // gauss point weight
  double weightPsi;

  std::vector< unsigned > sysDof; // local to global pdeSys dofs
  std::vector< double > Res; // local redidual std::vector
  std::vector < double > Jac;

  double eps = 0.00000001;


  {
    Solution*    solC        = ml_prob0->_ml_sol->GetSolutionLevel(levelC);

    for(unsigned d = 0; d < dim; d++) {
      *(solC->_Sol[solVIndex[d]]) = *(sol2->_Sol[solVIndex[d]]); //TODO prolongation of sol2 into sol0^ln
    }
    *(solC->_Sol[solP1Index]) = *(sol2->_Sol[solP1Index]); //TODO prolongation of sol2 into sol0^ln
    *(solC->_Sol[solP2Index]) = *(sol2->_Sol[solP2Index]); //TODO prolongation of sol2 into sol0^ln

    for(unsigned level = levelC; level < levelF; level++) {
      solC = ml_prob0->_ml_sol->GetSolutionLevel(level);
      Solution*    solF = ml_prob0->_ml_sol->GetSolutionLevel(level + 1);
      Mesh*        mshF = ml_prob0->_ml_msh->GetLevel(level + 1);
      for(unsigned d = 0; d < dim; d++) {
        solF->_Sol[solVIndex[d]]->matrix_mult(*(solC->_Sol[solVIndex[d]]), *(mshF->GetCoarseToFineProjection(solVType)));
      }
      solF->_Sol[solP1Index]->matrix_mult(*(solC->_Sol[solP1Index]), *(mshF->GetCoarseToFineProjection(solPType)));
      solF->_Sol[solP2Index]->matrix_mult(*(solC->_Sol[solP2Index]), *(mshF->GetCoarseToFineProjection(solPType)));
    }
  }


  KK->zero();
  RES->zero();

  AssembleStabilizationTerms(*ml_prob0);


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
    const elem_type *femPsi = msh->_finiteElement[ielGeom][psiType];
    const elem_type *femN = msh->_finiteElement[ielGeom][solNType];
    const elem_type *femK = msh->_finiteElement[ielGeom][solKType];
    //unsigned cnt =

    unsigned nDofsPsi;
    unsigned nDofsN;
    unsigned nDofsK;

    cut = (fabs(C - 0.5) < 0.1 ) ? 1 : 0;

    if(cut == 1) {
      femV = fem.GetFiniteElement(ielGeom, solVType);
      femP = fem.GetFiniteElement(ielGeom, solPType);
      femPsi = fem.GetFiniteElement(ielGeom, psiType);
      femN = fem.GetFiniteElement(ielGeom, solNType);
      femK = fem.GetFiniteElement(ielGeom, solKType);

      nDofsPsi = msh->GetElementDofNumber(iel, psiType);
      psi.resize(nDofsPsi);
      for(unsigned i = 0; i < nDofsPsi; i++) {
        unsigned psiDof = msh->GetSolutionDof(i, iel, psiType);
        psi[i] = (*sol->_Sol[psiIndex])(psiDof);
      }

      nDofsN = msh->GetElementDofNumber(iel, solNType);
      for (unsigned d = 0; d < dim; d++) {
        n[d].resize(nDofsN);

        for(unsigned i = 0; i < nDofsN; i++) {
          unsigned solNDof = msh->GetSolutionDof(i, iel, solNType);
          n[d][i] = (*sol->_Sol[solNIndex[d]])(solNDof);
        }
      }

      nDofsK = msh->GetElementDofNumber(iel, solKType);
      k.resize(nDofsK);
      for(unsigned i = 0; i < nDofsK; i++) {
        unsigned solKDof = msh->GetSolutionDof(i, iel, solKType);
        k[i] = (*sol->_Sol[solKIndex])(solKDof);
      }


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

    double NN_exact[2];
    double kk = 0.;



    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < femV->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      femV->Jacobian(coordX, ig, weight, phiV, phiV_x);
      femN->Jacobian(coordX, ig, weightPsi, phiN, phiN_x);
      femK->Jacobian(coordX, ig, weightPsi, phiK, phiK_x);
      femPsi->Jacobian(coordX, ig, weightPsi, phiPsi, phiPsi_x, phiPsi_xx);
      phiP = femP->GetPhi(ig);

      double dsN = 0.;
      std::vector <double> Nf(dim, 0); // unit normal in the physical element from the fluid to the solid

      std::vector<double> Ng(dim, 0.);
      double Kg = 0.;

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


        double det_exact = 0.;
        for(unsigned k = 0; k < dim; k++) det_exact += xg[k] * xg[k];
        det_exact = std::sqrt(det_exact);
        for(unsigned k = 0; k < dim; k++) NN_exact[k] = xg[k] / det_exact;

        kk = 1. / RADIUS;

        for(unsigned d = 0; d < dim; d++)
          NN[d] = 0;
        std::vector<double> hess(dim2);
        for (unsigned i = 0; i < nDofsPsi; i++) {
          for(unsigned d = 0; d < dim; d++) {
            NN[d] -= psi[i] * phiPsi_x[i * dim + d];
          }
          for (unsigned d = 0; d < dim2; d++) {
            hess[d] += psi[i] * phiPsi_xx[i * dim2 + d];
          }
        }
        double det = 0;
        for (unsigned j = 0; j < dim; j++) {
          det += NN[j] * NN[j];
        }
        det = sqrt(det);
        for (unsigned j = 0; j < dim; j++) {
          NN[j] /= det;
        }
        double H = 0;
        for(unsigned J = 0; J < dim; J++) {
          for(unsigned K = 0; K < dim; K++) {
            //2D xx, yy, xy
            //3D xx, yy, zz, xy, yz ,zx
            unsigned L;
            if(J == K) L = J;
            else if(1 == J + K) L = dim;     // xy
            else if(2 == J + K) L = dim + 2; // xz
            else if(3 == J + K) L = dim + 1; // yz
            H += NN[J] * hess[L] * NN[K];
          }
          H -= hess[J];
        }
        H /= (dim - 1) * det;

        kk = H;


        //===========================================================================================

        for (unsigned i = 0; i < nDofsN; i++) {
          for(unsigned d = 0; d < dim; d++) {
            Ng[d] += n[d][i] * phiN[i];
          }
        }
        for (unsigned i = 0; i < nDofsK; i++) {
          Kg += k[i] * phiK[i];
        }

        // std::cerr<<"H "<<H<<" NN "<<NN[0]<<" "<<NN[1] <<" vs kk "<<kk<<" n "<<NN_exact[0]<<" "<<NN_exact[1]<<std::endl;

      }


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
          for(unsigned k = 0; k < dim; k++) {
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

      double rho = rho1 * weightCFInt[ig] + rho2 * weightCFExt[ig];
      double mu = mu1 * weightCFInt[ig] + mu2 * weightCFExt[ig];

      //double rho = rho1 * C + rho2 * (1. - C);
      //double mu = mu1 * C + mu2 * (1. - C);

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
          NSV += - rho * phiV[i] * g[I]; // gravity term
          Res[I * nDofsV + i] -=  NSV * weight;
          if(cut == 1) {
            //std::cout << - sigma * phiV[i] * NN[I] * weight * weightCF[ig] * kk * dsN << " ";
            //Res[I * nDofsV + i] += - sigma * phiV[i] * NN[I] * weight * weightCF[ig] * kk * dsN;
            Res[I * nDofsV + i] += - sigma * phiV[i] * Ng[I] * weight * weightCF[ig] * Kg * dsN;

            // std::vector<std::vector<double>> P (dim);
            // for (int d = 0; d < dim; d ++)
            //   P[d].resize(dim);
            //
            // for(int i = 0; i < dim; i++){
            //   for(int j = 0; j < dim; j++){
            //     if(i==j) P[i][j] += 1.;
            //     P[i][j] -= Ng[i] * Ng[j];
            //   }
            // }
            //
            // for (int d = 0; d < dim; d++) {
            //   Res[I * nDofsV + i] += - sigma  * P[I][d] * phiV_x[i * dim + d] * weight * weightCF[ig] * dsN;
            // }
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


  //TODO restrict KK into KK2 space
  KK2->matrix_add (1., *KK, "different_nonzero_pattern");
  *RES2 += *RES;


  // Mat A = (static_cast<PetscMatrix*>(KK2))->mat();
  // Mat B = (static_cast<PetscMatrix*>(KK))->mat();
  //
  // PetscInt AM, AN, BM, BN;
  // PetscInt Am, An, Bm, Bn;
  // PetscInt ars, are, brs, bre;
  // PetscInt acs, ace, bcs, bce;
  //
  // MatGetSize(A, &AM, &AN);
  // MatGetSize(B, &BM, &BN);
  //
  // MatGetLocalSize(A, &Am, &An);
  // MatGetLocalSize(B, &Bm, &Bn);
  //
  // MatGetOwnershipRange(A, &ars, &are);
  // MatGetOwnershipRange(B, &brs, &bre);
  //
  // MatGetOwnershipRangeColumn(A, &acs, &ace);
  // MatGetOwnershipRangeColumn(B, &bcs, &bce);
  //
  // PetscPrintf(PETSC_COMM_WORLD,
  //             "KK2 global %d x %d, local %d x %d\n",
  //             AM, AN, Am, An);
  //
  // PetscPrintf(PETSC_COMM_WORLD,
  //             "KK  global %d x %d, local %d x %d\n",
  //             BM, BN, Bm, Bn);
  //
  // PetscSynchronizedPrintf(PETSC_COMM_WORLD,
  //                         "KK2 rows [%d,%d), cols [%d,%d)\n",
  //                         ars, are, acs, ace);
  //
  // PetscSynchronizedPrintf(PETSC_COMM_WORLD,
  //                         "KK  rows [%d,%d), cols [%d,%d)\n",
  //                         brs, bre, bcs, bce);
  //
  // PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);








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

    const double w = std::exp(- 10 * f * f / s2);

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
