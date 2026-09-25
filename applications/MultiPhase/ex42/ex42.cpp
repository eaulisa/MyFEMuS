#include <sys/stat.h>
#include <sys/types.h>

#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"

#include "TransientSystem.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "adept.h"
#include "PolynomialBases.hpp"

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
#include "../includeLS/Stabilization.hpp"
#include "../includeLS/GhostPenalty.hpp"
#include "../includeLS/GhostPenaltyDGP.hpp"

bool printdb = false;

const RungeKutta::VelKind velocityType = RungeKutta::VelKind::Zero;
//const RungeKutta::VelKind velocityType = RungeKutta::VelKind::Translation;
//const RungeKutta::VelKind velocityType = RungeKutta::VelKind::Rotation;
//const RungeKutta::VelKind velocityType = RungeKutta::VelKind::Vortex;

#include "../includeLS/Utils.hpp"

struct TaylorHood {
  struct Variable {
    FEFamily family;
    FEOrder order;
  };

  Variable U;
  Variable p;

  TaylorHood(const unsigned type) {
    if (type == 1) {
      U.family = LAGRANGE;
      U.order = FIRST;

      p.family = DISCONTINUOUS_POLYNOMIAL;
      p.order = ZERO;
    }
    else if (type == 2) {
      U.family = LAGRANGE;
      U.order = SECOND;

      p.family = LAGRANGE;
      p.order = FIRST;
    }
    else if (type == 3) {
      U.family = LAGRANGE;
      U.order = FIRST;

      p.family = LAGRANGE;
      p.order = FIRST;
    }
    else if (type == 4) {
      U.family = LAGRANGE;
      U.order = SECOND;

      p.family = LAGRANGE;
      p.order = SECOND;
    }
    else if (type == 5) {
      U.family = LAGRANGE;
      U.order = SECOND;

      p.family = DISCONTINUOUS_POLYNOMIAL;
      p.order = ZERO;
    }
    else {
      throw std::runtime_error("TaylorHood: type must be 1 or 2");
    }
  }
};

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet

  if(!strcmp(SolName, "U")) {  // strcmp compares two string in lexiographic sense.
    //if(facename == 1 || facename == 3) dirichlet = false;
    value = 0.;
  }
  else if(!strcmp(SolName, "V")) {
    // if(facename == 2 || facename == 4) dirichlet = false;
    if ((x[0] < 1.e-10 || x[0] > 1. - 1.e-10) && (x[1] > 1.e-10 && x[1] < 2 - 1.e-10)) {
      dirichlet = false;
    }
    value = 0.;
  }
  else if(!strcmp(SolName, "W")) {
    value = 0.;
  }
  else if(!strcmp(SolName, "P1") || !strcmp(SolName, "P2") ) {
    dirichlet = false;
    value = 0.;
  }
  else if(!strcmp(SolName, "NX") || !strcmp(SolName, "NY") || !strcmp(SolName, "NZ") || !strcmp(SolName, "K") || !strcmp(SolName, "AuxPsi")) {
    dirichlet = false;
    value = 0.;
  }

  return dirichlet;
}

void AssembleMultiphase(MultiLevelProblem& ml_prob);
void AssembleNormalLumped(MultiLevelProblem& ml_prob);
void AssembleCurvatureLumped(MultiLevelProblem& ml_prob);
void AssembleNormal(MultiLevelProblem& ml_prob);
void AssembleCurvature(MultiLevelProblem& ml_prob);
void AssembleSmoothLevelSet(MultiLevelProblem& ml_prob);
// double TimeStepMultiphase(const double time);

int main(int argc, char **argv) {

  // Initialize PETSc/MPI
  FemusInit mpinit(argc, argv, MPI_COMM_WORLD);

  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if (nprocs == 1)
    ProfilerStart("profiling.prof");

  int iproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);

  SimulationArgs args;

  try {
    args = ParseSimulationArgs(argc, argv);
  }
  catch (const std::exception& e) {

    if (iproc == 0) {
      std::cerr
          << "Command-line error: "
          << e.what()
          << '\n';
    }

    return 1;
  }

  const unsigned numberOfUniformLevels = args.uniformLevels;
  const std::vector<unsigned>& adaptiveLevelsList = args.adaptiveLevels;
  const std::vector<unsigned>& nStepsList = args.nSteps;
  const unsigned levelOffset = args.levelOffset;

  //domain settings
  double xmin, xmax, ymin, ymax, zmin, zmax = 0.;
  xmax = 1.;
  ymax = 2.;
  int nx, ny, nz = 0;
  nx = 8;
  ny = 16;

  // time settings
  double period = 0.5;
  static double dt = 0.0;

  // interface settings
  const double r = 0.25;
  std::vector<double> xc = {0.5, 0.5, 0.5};

  // variables name
  std::string psiName = "Psi";
  std::vector<std::string> vName = {"U", "V", "W"};
  std::vector<std::string> pName = {"P1", "P2"};
  std::string cName = "C";

  for (const unsigned numberOfSelectiveLevels : adaptiveLevelsList) {

    const double scalingFactor = 1.0;

    unsigned levelN = numberOfUniformLevels + numberOfSelectiveLevels;
    const unsigned levelF = levelN - 1u; //fine level associated for mlmsh0 and mlmsh1
    const unsigned levelC = levelF - levelOffset; //coarse level associated to mlmsh2, but existing also mlmsh0 and mlmsh1
    const unsigned level0 = 0;//levelC;

    // create uniform fine mesh for storing solutions to compare
    MultiLevelMesh mlMshReference;
    mlMshReference.GenerateCoarseBoxMesh(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, QUAD9, "seventh");
    mlMshReference.RefineMesh(levelF + 1, levelF + 1, nullptr);

    // prepare uniform solution vector
    std::vector<TemporalSnapshot> temporalSnapshots;
    temporalSnapshots.reserve(nStepsList.size());

    for (const unsigned nSteps : nStepsList) {

      MultiLevelMesh mlMsh0;
      // std::string meshName = "./input/tri.neu";
      // mlMsh0.ReadCoarseMesh(meshName.c_str(), "seventh", scalingFactor);
      mlMsh0.GenerateCoarseBoxMesh(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, QUAD9, "seventh"); // Turek 1&2

      dt = period / nSteps;

      auto TimeStepMultiphase = [](const double time) {
        return dt;
      };
      auto firstTimeStepMultiphase = [](const double time) {
        return 0.5 * dt;
      };

      // const unsigned nPrint  = std::max(1u, nSteps / 1u);
      // const unsigned nReinit = std::max(1u, nSteps / 100u);

      const unsigned nPrint = 1;
      const unsigned nReinit = 100000;

      mlMsh0.RefineMesh(numberOfUniformLevels, numberOfUniformLevels, nullptr);

      unsigned dim = mlMsh0.GetDimension();

      xc.resize(dim);

      unsigned sigmoidType = 0;
      double eps = 0.25; //(dim == 2) ? 1. / pow(2, std::max(levelN - 7u, 1u))
      Mollifier m = Mollifier(eps, sigmoidType);

      PsiBall psi2D(xc, r, m);

      // Iteratively flag and create new AMR levels
      for (unsigned k = 0; k < numberOfSelectiveLevels; ++k) {
        FlagFinestMeshLevel(mlMsh0, psi2D);
        mlMsh0.AddAMRMeshLevel(false);
      }

      mlMsh0.PrintInfo();
      BBoxToIel bbox(mlMsh0, 0, 3);

      // Define solution on the multilevel mesh
      MultiLevelSolution mlSol0(&mlMsh0);

      // FEOrder velOrder = SECOND;
      TaylorHood TH(2);
      std::cout << "Velocity discretization : " << TH.U.family << " " << TH.U.order << std::endl;
      std::cout << "Pressure discretization : " << TH.p.family << " " << TH.p.order << std::endl;

      vName.resize(dim);

      // final uniform solutions
      auto mlSolReference = std::make_unique<MultiLevelSolution>(&mlMshReference);
      mlSolReference->AddSolution(psiName.c_str(), LAGRANGE, SECOND, 0, false);
      for(unsigned d = 0; d < dim; ++d) mlSolReference->AddSolution(vName[d].c_str(), TH.U.family, TH.U.order, 0, false);

      mlSolReference->Initialize("All");

      // system solutions
      mlSol0.AddSolution(psiName.c_str(), LAGRANGE, SECOND, 0, false);
      for(unsigned d = 0; d < dim; d++) mlSol0.AddSolution(vName[d].c_str(), TH.U.family, TH.U.order, 2);
      for(unsigned d = 0; d < pName.size(); d++) mlSol0.AddSolution(pName[d].c_str(), TH.p.family, TH.p.order, 2);
      mlSol0.AddSolution(cName.c_str(), DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);

      mlSol0.Initialize("All");

      InitSol(mlSol0, vName, 0, 1.);

      Mesh* mshC = mlMsh0.GetLevel(levelC);
      const double hC = GetMaxElementH(mshC, levelC);

      // std::vector<double> xc_1 = xc;
      // double r_0 = r;
      // Circle c1(xc_1, r_0);

      Boundary zero_bd;

      InitLevelSet(mlSol0, psiName, psi2D);
      // UpdateColorFunction(mlSol0, psiName, cName);

      // Export solution to VTK (selected levels)
      std::vector<std::string> variablesToBePrinted = {"All"};

      MultiLevelMesh mlMsh1;
      MultiLevelSolution mlSol1;

      MultiLevelMesh *mlmsh0 = &mlMsh0;
      MultiLevelMesh *mlmsh1 = &mlMsh1;

      MultiLevelSolution *mlsol0 = &mlSol0;
      MultiLevelSolution *mlsol1 = &mlSol1;

      // Initialize markers object
      LevelSetMarkers markers(psiName, dim);

      // Load coarse mesh and build uniform refinement levels
      mlmsh1->GenerateCoarseBoxMesh(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, QUAD9, "seventh"); // Turek 1&2
      // mlmsh1->ReadCoarseMesh(meshName.c_str(), "seventh", scalingFactor);
      mlmsh1->RefineMesh(numberOfUniformLevels, numberOfUniformLevels, nullptr);

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

      MultiphasePhysicalProperties properties;
      properties.mu1 = 0.1;
      properties.mu2 = 10.;
      properties.rho1 = 1.;
      properties.rho2 = 1000;
      properties.sigma = 0.0;//1.96;
      properties.gravity = -0.98;

      UpdateColorFunction(*mlsol0, psiName, cName);
      if(levelC < levelF) RestrictPWDCField(*mlsol0, cName, levelC, levelF);

      mlsol0->AttachSetBoundaryConditionFunction(SetBoundaryCondition);
      mlsol0->GenerateBdc("All");

      LevelSetDiagnostics diagnostics = ComputeLevelSetDiagnostics(*mlsol0, psiName, simulation_type::rising_bubble, vName, levelC);

      const unsigned w = 24;

      if(iproc == 0) mkdir("levelsetdiagnostic", 0755);

      const std::string diagnosticsFile =
        "levelsetdiagnostic/diagnostics_LC" + std::to_string(levelC) +
        "_LF" + std::to_string(levelF) + "_NS" + std::to_string(nSteps) + ".dat";

      if(iproc == 0) {
        std::ofstream out(diagnosticsFile, std::ios::trunc);

        out << std::setw(w) << "time"
            << std::setw(w) << "innerArea"
            << std::setw(w) << "outerArea"
            << std::setw(w) << "totalArea"
            << std::setw(w) << "interfaceLength";

        for(unsigned d = 0; d < dim; ++d)
          out << std::setw(w) << ("barycenter" + std::to_string(d));

        for(unsigned d = 0; d < dim; ++d)
          out << std::setw(w) << ("meanVelocity" + std::to_string(d));

        out << std::setw(w) << "circularity"
            << '\n';
      }

      PrintLevelSetDiagnostics(diagnostics, iproc, 0.0, diagnosticsFile);

      const std::string ls_outputdir = "output_ls";
      const std::string vel_outputdir = "output_vel";

      VTKWriter vtkIO(&mlSol0);
      vtkIO.Write(levelF, ls_outputdir, "biquadratic", variablesToBePrinted, 0);

      LevelSetDiagnostics final_diagnostics;

      for (unsigned t = 1; t <= 0 + 1 * nSteps + 1; t++) {

        TimeDiscretization td = (t == 1) ? TimeDiscretization::BackEuler : TimeDiscretization::CrankNicholson;

        double time = t * dt;

        mlsol0->CopySolutionToOldSolution();

        // mlProb2 is used to assemble the Ghost Penalty and solve the full Navier-Stokes + Ghost penalty

        MultiLevelMesh mlMsh2(*mlmsh0, level0, levelC + 1, "seventh");
        std::vector<Mesh*> msh2(levelC + 1 - level0);
        for(unsigned l = 0; l < msh2.size(); l++)
          msh2[l] = mlMsh2.GetLevel(l);

        MultiLevelSolution mlSol2(&mlMsh2);
        mlSol2.AddSolution(psiName.c_str(), LAGRANGE, SECOND, 0, false);
        for(unsigned d = 0; d < dim; d++) mlSol2.AddSolution(vName[d].c_str(), TH.U.family, TH.U.order, 2);
        for(unsigned d = 0; d < pName.size(); d++) mlSol2.AddSolution(pName[d].c_str(), TH.p.family, TH.p.order, 2);
        mlSol2.AddSolution(cName.c_str(), DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);

        mlSol2.Initialize("All");

        mlSol2.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
        bool allLevel = true;
        mlSol2.GenerateBdc("All");

        auto &Sol0_C = (mlsol0->GetSolutionLevel(levelC))->_Sol;
        auto &Sol0Old_C = (mlsol0->GetSolutionLevel(levelC))->_SolOld;

        auto &Sol2 = (mlSol2.GetSolutionLevel(levelC - level0))->_Sol; //careful here sol2(0) corresponds to the same level of sol0(l0)
        auto &Sol2Old = (mlSol2.GetSolutionLevel(levelC - level0))->_SolOld;

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

        if(level0 < levelC) RestrictPWDCField(mlSol2, cName, 0, levelC - level0);

        std::vector<double> xtarget = {xmin, ymin, zmin};
        xtarget.resize(dim);
        SetUnphysicalPressureDofs(mlSol2, cName, pName, 0, levelC - level0, xtarget, true);

        MultiLevelProblem mlProb2(&mlSol2);

        // add system Navier-Stokes in mlProb as a Linear Implicit System
        TransientNonlinearImplicitSystem& system2 = mlProb2.add_system < TransientNonlinearImplicitSystem > ("NS");

        // add velocity to system
        for(unsigned d = 0; d < dim; d++) system2.AddSolutionToSystemPDE(vName[d].c_str());
        //add pressure
        for(unsigned d = 0; d < pName.size(); d++) system2.AddSolutionToSystemPDE(pName[d].c_str());

        unsigned sparsity_pattern_size = dim * std::pow(5, dim) + 2 * std::pow(3, dim) + 4 * dim * std::pow(5, dim - 1); // only Q2-Q1
        system2.SetSparsityPatternMinimumSize(sparsity_pattern_size);
        // attach the assembling function to system
        system2.SetAssembleFunction(AssembleMultiphase);
        system2.AttachGetTimeIntervalFunction(
          (t == 1) ? firstTimeStepMultiphase : TimeStepMultiphase
        );
        // initilaize and solve the system

        //system2.SetOuterSolver(PREONLY);
        system2.SetMaxNumberOfNonLinearIterations(10);

        MultiLevelProblem mlProb0(mlsol0);
        // add system Navier-Stokes in mlProb as a Linear Implicit System
        TransientNonlinearImplicitSystem& system0 = mlProb0.add_system < TransientNonlinearImplicitSystem > ("NS");
        // add velocity to system
        for(unsigned d = 0; d < dim; d++) system0.AddSolutionToSystemPDE(vName[d].c_str());
        //add pressure
        for(unsigned d = 0; d < pName.size(); d++) system0.AddSolutionToSystemPDE(pName[d].c_str());
        //system0.SetSparsityPatternMinimumSize(250);
        // attach the assembling function to system
        system0.SetAssembleFunction(AssembleMultiphase);
        system0.AttachGetTimeIntervalFunction(
          (t == 1) ? firstTimeStepMultiphase : TimeStepMultiphase
        );
        // initilaize and solve the system
        system0.init();

        mlProb2.SetMultiphaseParams(&mlProb0, levelF, levelC, level0, hC, properties, td);
        mlProb0.SetMultiphaseParams(nullptr, levelF, levelC, level0, hC, properties, td);

        for(unsigned l = 0; l < msh2.size(); l++)
          msh2[l]->SetLevel(l);

        if (t == 1)
          system2.SetMgType(V_CYCLE);
        else
          system2.SetMgType(V_CYCLE);

        system2.SetLinearEquationSolverType(FEMuS_ASM);

        system2.init();

        for (unsigned l = 0; l < levelC + 1 - level0; l++) {
          LinearEquationSolver* pdeSys2_l  = system2._LinSolver[l];
          // const  std::vector<NumericVector*> *
          pdeSys2_l->SetSolution(&mlSol2.GetSolutionLevel(l)->_Sol);
          pdeSys2_l->MergeNullSpaceBases(true);

        }

        // ******* Set Smoother *******

        // system2.SetSolverFineGrids(GMRES);

        system2.SetSolverFineGrids(RICHARDSON);
        system2.SetRichardsonScaleFactor(.8);

        system2.SetNumberPreSmoothingStep(4);
        system2.SetNumberPostSmoothingStep(4);
        // system2.SetTolerances(1.e-20, 1.e-20, 1.e+50, 50, 30);

        system2.SetPreconditionerFineGrids(MLU_PRECOND);
        system2.SetTolerances(1.e-10, 1.e-12, 1.e+50, 40, 40);

        system2.SetNumberOfSchurVariables(2);
        system2.SetElementBlockNumber(3);

        //system2.SetPreconditionerFineGrids(ILU_PRECOND);
        system2.MGsolve();
        for(unsigned l = 0; l < msh2.size(); l++)
          msh2[l]->SetLevel(l + level0);

        VTKWriter vtkIO2(&mlSol2);
        vtkIO2.SetDebugOutput(true);
        if (t % nPrint == 0) {
          vtkIO2.Write(levelC - level0, vel_outputdir, "biquadratic", variablesToBePrinted, t / 1);
        }

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

        if(t == nSteps + 1) {
          Solution* solC = mlsol0->GetSolutionLevel(levelC);

          for(unsigned d = 0; d < dim; ++d) {
            const unsigned velIndex = mlsol0->GetIndex(vName[d].c_str());

            NumericVector* uNew = solC->_Sol[velIndex];
            NumericVector* uOld = solC->_SolOld[velIndex];

            const unsigned first = uNew->first_local_index();
            const unsigned last = uNew->last_local_index();

            for(unsigned i = first; i < last; ++i) {
              const double value = 0.5 * ((*uNew)(i) + (*uOld)(i));
              uNew->set(i, value);
            }
            uNew->close();
          }

          bbox.SetMesh(mlmsh0->GetLevel(0));
          ProjectSolution(*mlsol0, *mlSolReference, bbox, vName, levelC, levelF);

          final_diagnostics = ComputeLevelSetDiagnostics(*mlSolReference, psiName, simulation_type::rising_bubble, vName, levelF);

          break;
        }

        bbox.SetMesh(mlmsh0->GetLevel(0));

        unsigned nLevels = numberOfUniformLevels + numberOfSelectiveLevels;
        std::vector<MyVector<double>> X0;
        MyVector<int> X0Iel;
        // GetCutElementPoints(*mlsol0, psiName, X0, X0Iel);

        std::vector<std::vector<std::vector<double>>> inflow_markers0(0);

        std::vector<std::vector<double>> inflow_markers(dim);
        // zero_bd.updateMarkers(inflow_markers0, inflow_markers, time-dt, period, dt);
        //
        // std::vector<MyVector<double>> IX(dim);
        // for (unsigned k = 0; k < dim; ++k) {
        //   IX[k].buildFromLocal(inflow_markers[k]);
        // }
        //
        // {
        //   LevelMarkers l0;
        //   const unsigned bboxLevels = nLevels - bbox.GetLevel();
        //
        //   std::vector<LevelMarkers> lX(bboxLevels);
        //
        //   bbox.GetInverseMappingOnCoarseLevel(IX, l0, lX[0]);
        //
        //   const std::vector<bool> &isInsideDomain = l0.GetPointInsideDomain();
        //
        //   for (unsigned d = 0; d < dim; d++)
        //     inflow_markers[d].clear();
        //
        //   for (unsigned d = 0; d < dim; d++) {
        //     unsigned offset = IX[d].begin();
        //     for (unsigned i = IX[d].begin(); i < IX[d].end(); ++i) {
        //       if (!isInsideDomain[i - offset]) {
        //           inflow_markers[d].push_back(IX[d][i]);
        //       }
        //     }
        //   }
        // }

        markers.GetCutElementPoints(*mlsol0, X0, X0Iel, inflow_markers);

        // if (t % nReinit == 0) {
        //   Reinit reinit(psiName, *mlsol0, m);

        //   reinit.farFieldReinit(X0);
        //   reinit.interfaceFieldReinit(bbox);
        //   reinit.updateSolution();
        // }

        // if (t == 1)
        //   WritePointsVTK("./output/points.0.vtk", X0);

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
        for(unsigned d = 0; d < dim; d++) mlsol1->AddSolution(vName[d].c_str(), TH.U.family, TH.U.order, 2);
        for(unsigned d = 0; d < pName.size(); d++) mlsol1->AddSolution(pName[d].c_str(), TH.p.family, TH.p.order, 2);
        mlsol1->AddSolution(cName.c_str(), DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);

        mlsol1->Initialize("All");
        mlsol1->AttachSetBoundaryConditionFunction(SetBoundaryCondition);
        mlsol1->GenerateBdc("All");

        ProjectSolution(*mlsol0, *mlsol1, bbox, {psiName}, levelF, levelF, vName, levelC, zero_bd, -dt, time, period);
        ProjectSolution(*mlsol0, *mlsol1, bbox, vName, levelC, levelC);

        UpdateColorFunction(*mlsol1, psiName, cName);
        if(levelC < levelF) RestrictPWDCField(*mlsol1, cName, levelC, levelF);

        for(unsigned i = 0; i < cfw.size(); i++) cfw[i]->ClearMap();

        LevelSetDiagnostics diagnostics = ComputeLevelSetDiagnostics(*mlsol1, psiName, simulation_type::rising_bubble, vName, levelC);
        PrintLevelSetDiagnostics(diagnostics, iproc, time, diagnosticsFile);

        // Export solution to VTK (selected levels)
        VTKWriter vtkIO1(mlsol1);
        if (t % nPrint == 0) {
          vtkIO1.Write(levelF, ls_outputdir, "biquadratic", variablesToBePrinted, t / 1);
        }

        std::swap(mlsol0, mlsol1);
        std::swap(mlmsh0, mlmsh1);

        mlsol1->clear();
        mlmsh1->resize(numberOfUniformLevels);

        for(unsigned i = 0; i < cfw.size(); i++) cfw[i]->clear();

        if(t == nSteps) {
          bbox.SetMesh(mlmsh0->GetLevel(0));
          ProjectSolution(*mlsol0, *mlSolReference, bbox, {psiName}, levelF, levelF);
        }

      }

      temporalSnapshots.push_back( {
        nSteps,
        dt,
        final_diagnostics,
        std::move(mlSolReference)
      }
                                 );

    }

    PrintTemporalConvergence(temporalSnapshots, psiName, vName, levelF, iproc);

  }

  if (nprocs == 1)
    ProfilerStop();
  return 0;
}

// double TimeStepMultiphase(const double time) {
//   // double dt =  0.005; //RT
//   // double dt =  0.001; //RT
//   double dt =  0.0025; //Turek
//   // double sigma = 3;
//   // double rho = 100.;
//   // // double totalT = sqrt(rho*0.4*0.4*0.4) / sqrt(sigma);
//   // // double dt =  totalT/800; //Parasitic Test
//   //
//   // double dt =   0.001 * sqrt(rho * 0.4 * 0.4 * 0.4 / sigma);
//   // // double dt =  0.0001; //TODO if you use the 320x320 you have to change this
//   return dt;
// }

void AssembleMultiphase(MultiLevelProblem& ml_prob2) {

  //  extract pointers to the several objects that we are going to use
  TransientNonlinearImplicitSystem* mlPdeSys2   = &ml_prob2.get_system<TransientNonlinearImplicitSystem> ("NS");
  const unsigned level2 = mlPdeSys2->GetLevelToAssemble();

  Mesh* msh2 = ml_prob2._ml_msh->GetLevel(level2);    // pointer to the mesh (level) object
  elem* el2 = msh2->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol2        = ml_prob2._ml_sol;  // pointer to the multilevel solution object
  Solution* sol2 = ml_prob2._ml_sol->GetSolutionLevel(level2);    // pointer to the solution (level) object

  const unsigned  dim = msh2->GetDimension(); // get the domain dimension of the problem
  std::vector < unsigned > sol2VIndex(dim);
  sol2VIndex[0] = mlSol2->GetIndex("U");    // get the position of "U" in the ml_sol object
  sol2VIndex[1] = mlSol2->GetIndex("V");    // get the position of "V" in the ml_sol object
  if(dim == 3) sol2VIndex[2] = mlSol2->GetIndex("W");

  unsigned sol2P1Index = mlSol2->GetIndex("P1");    // get the position of "P1" in the ml_sol object
  unsigned sol2P2Index = mlSol2->GetIndex("P2");    // get the position of "P2" in the ml_sol object
  unsigned sol2PType = mlSol2->GetSolutionType(sol2P1Index);

  LinearEquationSolver* pdeSys2        = mlPdeSys2->_LinSolver[level2]; // pointer to the equation (level) object
  SparseMatrix* KK2 = pdeSys2->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* RES2 = pdeSys2->_RES; // pointer to the global residual std::vector object in pdeSys (level)

  //MatResetPreallocation((static_cast< PetscMatrix* >(KK2))->mat());
  MatSetOption((static_cast< PetscMatrix* >(KK2))->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  KK2->zero();
  RES2->zero();

  if (sol2PType != mlSol2->GetSolutionType(sol2P2Index)) {
    throw std::runtime_error("Pressure type mismatch");
  }

  AssembleGhostPenaltyVelocity(ml_prob2);

  if (sol2PType == 3) {
    AssembleGhostPenaltyDGP(ml_prob2, true);
    AssembleGhostPenaltyDGP(ml_prob2, false);
  }
  else if (sol2PType == 0) {
    AssembleGhostPenaltyLinearPressure(ml_prob2, true);
    AssembleGhostPenaltyLinearPressure(ml_prob2, false);
  }
  else {
    throw std::runtime_error("Pressure type must be 3 or 0");
  }

  AssembleStabilizationTerms(ml_prob2);

  RES2->close();
  KK2->close();

  MultiphaseParams mParam = ml_prob2.GetMultiphaseParams();
  MultiLevelProblem *ml_prob0 = mParam.mlProbF;
  const unsigned levelF = mParam.levelF;
  const unsigned levelC = mParam.levelC;
  const unsigned level0 = mParam.level0;

  const TimeDiscretization td = mParam.td;
  const double cold = (td == TimeDiscretization::CrankNicholson) ? 0.5 : 0.;
  const double cnew = (td == TimeDiscretization::CrankNicholson) ? 0.5 : 1.;

  std::cout << "levelC = " << levelC << " levelF = " << levelF << std::endl;
  std::cout << "level to assemble = " << level2  << " mapping level to assemble to levelC = " << level0 + level2 << std::endl;

  MultiphasePhysicalProperties properties = mParam.properties;

  const bool firstNonlinearIt = (mlPdeSys2->GetNonlinearIt() == 0);

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
  //MatSetOption((static_cast< PetscMatrix* >(KK))->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  double mu1 = properties.mu1;
  double mu2 = properties.mu2;
  double rho1 = properties.rho1;
  double rho2 = properties.rho2;
  double sigma = properties.sigma;
  double gravity = properties.gravity;
  std::vector <double> g;
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

  std::vector < std::vector < double > >  solV(dim);    // local solution
  std::vector < std::vector < double > >  solVOld(dim);    // local solution
  std::vector < double >  solP1; // local solution
  std::vector < double >  solP2; // local solution

  std::vector < double >  psi; // local solution
  std::vector < double >  k;
  std::vector < std::vector < double > >  n(dim);

  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  unsigned solXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

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

  double eps = 0.;//1.e-14;

  if(printdb) std::cout << "Before Solution Projection\n" << std::flush;

  {

    Solution* sol2_l  = ml_prob2._ml_sol->GetSolutionLevel(level2);    // pointer to the solution (level) object
    Solution* sol0_l  = ml_prob0->_ml_sol->GetSolutionLevel(level0 + level2);

    for(unsigned d = 0; d < dim; d++) {
      *(sol0_l->_Sol[solVIndex[d]]) = *(sol2_l->_Sol[sol2VIndex[d]]); //TODO prolongation of sol2 into sol0^ln
      *(sol0_l->_SolOld[solVIndex[d]]) = *(sol2_l->_SolOld[sol2VIndex[d]]);
    }
    *(sol0_l->_Sol[solP1Index]) = *(sol2_l->_Sol[solP1Index]); //TODO prolongation of sol2 into sol0^ln
    *(sol0_l->_Sol[solP2Index]) = *(sol2_l->_Sol[solP2Index]); //TODO prolongation of sol2 into sol0^ln

    for(unsigned level = level0 + level2; level < levelF; level++) {
      if(printdb) std::cout << "\t Before Velocity Solution Projection " << level << std::endl << std::flush;
      sol0_l = ml_prob0->_ml_sol->GetSolutionLevel(level);
      Solution*  sol0_lp1 = ml_prob0->_ml_sol->GetSolutionLevel(level + 1);
      Mesh*      msh0_lp1 = ml_prob0->_ml_msh->GetLevel(level + 1);
      for(unsigned d = 0; d < dim; d++) {
        sol0_lp1->_Sol[solVIndex[d]]->matrix_mult(*(sol0_l->_Sol[solVIndex[d]]), *(msh0_lp1->GetCoarseToFineProjection(solVType)));
        sol0_lp1->_SolOld[solVIndex[d]]->matrix_mult(*(sol0_l->_SolOld[solVIndex[d]]), *(msh0_lp1->GetCoarseToFineProjection(solVType)));
      }
      if(printdb) std::cout << "\t Afer Velocity Solution Projection " << level << std::endl << std::flush;
      if(printdb) std::cout << "\t Before Pressure Solution Projection " << level << std::endl << std::flush;
      sol0_lp1->_Sol[solP1Index]->matrix_mult(*(sol0_l->_Sol[solP1Index]), *(msh0_lp1->GetCoarseToFineProjection(solPType)));
      sol0_lp1->_Sol[solP2Index]->matrix_mult(*(sol0_l->_Sol[solP2Index]), *(msh0_lp1->GetCoarseToFineProjection(solPType)));
      if(printdb) std::cout << "\t After Pressure Solution Projection " << level << std::endl << std::flush;
    }
  }

  KK->zero();
  RES->zero();

  if(printdb) std::cout << "After Solution Projection\n" << std::flush;

  // AssembleStabilizationTerms(*ml_prob0);

  double epsilonp = 0;

  clock_t start_time = clock();

  // element loop: each process loops only on the elements that owns

  if(printdb) std::cout << "Before KK assembly\n" << std::flush;
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    double C = (*sol->_Sol[cIndex])(iel);

    short unsigned ielGeom = msh->GetElementType(iel);

    unsigned nDofsX = msh->GetElementDofNumber(iel, solXType);
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
      coordX[k].resize(nDofsX);
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
    for(unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, solXType);
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
    //unsigned cnt =

    unsigned nDofsPsi;
    unsigned nDofsN;
    unsigned nDofsK;

    cut = (fabs(C - 0.5) < 0.1 ) ? 1 : 0;

    if(cut == 1) {
      femV = fem.GetFiniteElement(ielGeom, solVType);
      femP = fem.GetFiniteElement(ielGeom, solPType);
      femPsi = fem.GetFiniteElement(ielGeom, psiType);

      nDofsPsi = msh->GetElementDofNumber(iel, psiType);
      psi.resize(nDofsPsi);
      for(unsigned i = 0; i < nDofsPsi; i++) {
        unsigned psiDof = msh->GetSolutionDof(i, iel, psiType);
        psi[i] = (*sol->_Sol[psiIndex])(psiDof);
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

      // (*cfw[ielGeom])(0, a, d, weightCFInt);
      // for(unsigned k = 0; k < dim; k++) a[k] = - a[k];
      // d = -d;
      // (*cfw[ielGeom])(-1, a, d, weightCF);
      // (*cfw[ielGeom])(0, a, d, weightCFExt);

      cfCDw0[ielGeom]->GetWeight(a, d, weightCFInt);
      for(unsigned i = 0; i < weightCFInt.size(); i++) weightCFExt[i] = 1. - weightCFInt[i];
      //(*cfw[ielGeom])(-1, a, d, weightCF);
      cfw[ielGeom]->GetWeightWithMap(-1, a, d, weightCF);
      /*
      for(unsigned k = 0; k < dim; k++) a[k] = - a[k];
      d = -d;*/

      //cfCDw0[ielGeom]->GetWeight(a, d, weightCFExt);
      //weightCFExt.resize(weightCFInt.size());

    }
    else {
      for(unsigned i = 0; i < weightCFInt.size(); i++) {
        weightCFInt[i] = C;
        weightCFExt[i] = 1. - C;
      }
    }

    std::vector<double> xg(dim);

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < femV->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      femV->Jacobian(coordX, ig, weight, phiV, phiV_x);
      femPsi->Jacobian(coordX, ig, weightPsi, phiPsi, phiPsi_x, phiPsi_xx);
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
        dsN = sqrt(dsN) + 1.e-20;
        for(unsigned k = 0; k < dim; k++) {
          Nf[k] /= dsN;
        }

        xg.assign(dim, 0);

        for(unsigned i = 0; i < nDofsV; i++) {
          for(unsigned k = 0; k < dim; k++) {
            xg[k] += coordX[k][i] * phiV[i];
          }
        }

      }

      std::vector < double > solV_gss(dim, 0);
      std::vector < double > solVOld_gss(dim, 0);
      std::vector < std::vector < double > > gradSolV_gss(dim);
      std::vector < std::vector < double > > gradSolVOld_gss(dim);

      for(unsigned  k = 0; k < dim; k++) {
        gradSolV_gss[k].assign(dim, 0.);
        gradSolVOld_gss[k].assign(dim, 0.);
      }

      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          solV_gss[k] += solV[k][i] * phiV[i];
          solVOld_gss[k] += solVOld[k][i] * phiV[i];
        }
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned k = 0; k < dim; k++) {
            gradSolV_gss[k][j] += solV[k][i] * phiV_x[i * dim + j];
            gradSolVOld_gss[k][j] += solVOld[k][i] * phiV_x[i * dim + j];
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

      double rhoC = rho1 * C + rho2 * (1. - C);

      // *** phiV_i loop ***
      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned  I = 0; I < dim; I++) {  //momentum equation in k
          double NSV = 0.;
          for(unsigned J = 0; J < dim; J++) {  // second index j in each equation
            // residual terms
            NSV   +=  cnew * mu * phiV_x[i * dim + J] * (gradSolV_gss[I][J] + gradSolV_gss[J][I]); // diffusion
            NSV   +=  cnew * rho * phiV[i] * (solV_gss[J] * gradSolV_gss[I][J]); // nonlinear term
            // crank-nicholson old terms
            NSV   +=  cold * mu * phiV_x[i * dim + J] * (gradSolVOld_gss[I][J] + gradSolVOld_gss[J][I]); // diffusion
            NSV   +=  cold * rho * phiV[i] * (solVOld_gss[J] * gradSolVOld_gss[I][J]); // nonlinear term
          }
          NSV += - phiV_x[i * dim + I] * (solP1_gss * weightCFInt[ig] + solP2_gss * weightCFExt[ig]);  // pressure gradient
          NSV += rho * phiV[i] * (solV_gss[I] - solVOld_gss[I]) / dt ;
          NSV += - rho * phiV[i] * g[I]; // gravity term
          Res[I * nDofsV + i] -=  NSV * weight;

          // surface tension stabilization -- rhs contribution
          if(cut == 1) {

            std::vector<std::vector<double>> P (dim);
            for (int d = 0; d < dim; d ++)
              P[d].resize(dim);

            for(int i = 0; i < dim; i++) {
              for(int j = 0; j < dim; j++) {
                if(i == j) P[i][j] += 1.;
                P[i][j] -= Nf[i] * Nf[j];
              }
            }

            for (int d = 0; d < dim; d++) {
              Res[I * nDofsV + i] += - sigma  * P[I][d] * phiV_x[i * dim + d] * weight * weightCF[ig] * dsN;
            }

            double stabSF = 0.0;

            for (unsigned a = 0; a < dim; ++a) {

              double gradDeltaU_tg = 0.0;
              double gradPhi_tg    = 0.0;

              for (unsigned b = 0; b < dim; ++b) {

                gradDeltaU_tg += P[a][b] * (cnew * gradSolV_gss[I][b] - cold * gradSolVOld_gss[I][b]);
                gradPhi_tg += P[a][b] * phiV_x[i * dim + b];

              }

              stabSF += gradDeltaU_tg * gradPhi_tg;
            }

            Res[I * nDofsV + i] += - cnew * sigma * dt * stabSF * weight * weightCF[ig] * dsN;

          }
        }
      } // end phiV_i loop

      // *** phiP_i loop ***
      for(unsigned i = 0; i < nDofsP; i++) {
        for(int I = 0; I < dim; I++) {
          Res[dim * nDofsV + i] += - gradSolV_gss[I][I] * phiP[i]  * weight * weightCFInt[ig]; //continuity
          Res[dim * nDofsV + nDofsP + i] += - gradSolV_gss[I][I] * phiP[i]  * weight * weightCFExt[ig]; //continuity

          Res[dim * nDofsV + i] += - 2 * cold * gradSolVOld_gss[I][I] * phiP[i]  * weight * weightCFInt[ig]; //continuity
          Res[dim * nDofsV + nDofsP + i] += - 2 * cold * gradSolVOld_gss[I][I] * phiP[i]  * weight * weightCFExt[ig]; //continuity
        }
        if(C == 0)
          Res[dim * nDofsV + i] += - solP1_gss * phiP[i]  * weight * (1 - C) * eps; //penalty
        if(C == 1)
          Res[dim * nDofsV + nDofsP + i] += - solP2_gss * phiP[i]  * weight * C * eps; //penalty

        if(C > 0.1)
          Res[dim * nDofsV + i] += - solP1_gss * phiP[i]  * weight * epsilonp; //penalty
        if(C < 0.9)
          Res[dim * nDofsV + nDofsP + i] += - solP2_gss * phiP[i] * weight * epsilonp; //penalty

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
              Jac[ VIrow * nDofsVP + VIcolumn ] += cnew * mu * phiV_x[i * dim + J] * phiV_x[j * dim + J] * weight; //diagonal diffusion
              Jac[ VIrow * nDofsVP + VJcolumn ] += cnew * mu * phiV_x[i * dim + J] * phiV_x[j * dim + I] * weight; //off-diagonal diffusion

              Jac[ VIrow * nDofsVP + VIcolumn ] += cnew * rho * phiV[i] * solV_gss[J] * phiV_x[j * dim + J] * weight; //diagonal nonlinear
              Jac[ VIrow * nDofsVP + VJcolumn ] += cnew * rho * phiV[i] * phiV[j] * gradSolV_gss[I][J] * weight; //off-diagonal nonlinear
            }

            // surface tension stabilization -- matrix contribution
            if(cut == 1) {

              double stabSFJac = 0.0;

              std::vector<std::vector<double>> P (dim);
              for (int d = 0; d < dim; d ++)
                P[d].resize(dim);

              for(int i = 0; i < dim; i++) {
                for(int j = 0; j < dim; j++) {
                  if(i == j) P[i][j] += 1.;
                  P[i][j] -= Nf[i] * Nf[j];
                }
              }

              for(unsigned a = 0; a < dim; ++a) {
                for(unsigned b = 0; b < dim; ++b) {

                  stabSFJac +=
                    cnew
                    * phiV_x[j * dim + a]
                    * P[a][b]
                    * phiV_x[i * dim + b];
                }
              }

              Jac[VIrow * nDofsVP + VIcolumn] += cnew * sigma * dt * stabSFJac * weight * weightCF[ig] * dsN;

            }
          }

          for(unsigned j = 0; j < nDofsP; j++) {
            unsigned P1column = dim * nDofsV + j;
            unsigned P2column = dim * nDofsV + nDofsP + j;
            Jac[VIrow * nDofsVP + P1column] += - phiV_x[i * dim + I] * phiP[j] * weight * weightCFInt[ig]; //pressure gradient
            Jac[VIrow * nDofsVP + P2column] += - phiV_x[i * dim + I] * phiP[j] * weight * weightCFExt[ig]; //pressure gradient

            Jac[P1column * nDofsVP + VIrow] -= - phiV_x[i * dim + I] * phiP[j] * weight * weightCFInt[ig]; //continuity
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

          if(C > 0.1)
            Jac[P1row * nDofsVP + P1column] += phiP[i] * phiP[j] * weight  * epsilonp; //continuity
          if(C < 0.9)
            Jac[P2row * nDofsVP + P2column] += phiP[i] * phiP[j] * weight  * epsilonp; //continuity
        }
      }

    }

    RES->add_vector_blocked(Res, sysDof);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);

  } //end element loop for each process

  RES->close();
  KK->close();

  if(printdb) std::cout << "After KK assembly\n" << std::flush;

  std::cout << "Matrix Assembly time        = " << static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC << std::flush << std::endl;
  start_time = clock();

  vector < SparseMatrix* > PP, RR, PPamr, RRamr;
  PP = mlPdeSys->GetProjectionMatrix();
  RR = mlPdeSys->GetRestrictionMatrix();
  PPamr = mlPdeSys->GetAMRProjectionMatrix();
  RRamr = mlPdeSys->GetAMRRestrictionMatrix();

  vector < LinearEquationSolver*> LinSolver = mlPdeSys->GetLinearSolver();

  MultiLevelMesh * mlmsh0 = ml_prob0->_ml_msh;
  for(unsigned level = levelF; level > level0 + level2; level--) {
    if(!mlmsh0->GetLevel(level)->GetIfHomogeneous() && level == levelF) { //AMR RESTRICTION
      if(printdb) std::cout << "Before KK amr restriction\n" << std::flush;

      if(!RRamr[level]) {
        (LinSolver[level]->_RESC)->matrix_mult_transpose(*LinSolver[level]->_RES, *PPamr[level]);
        *(LinSolver[level]->_RES) = *(LinSolver[level]->_RESC);
        LinSolver[level]->SwapMatrices();
        LinSolver[level]->_KK->matrix_PtAP(*PPamr[level], *LinSolver[level]->_KKamr, false); // cannot use !firstNonlinearIt here
      }
      else {
        (LinSolver[level]->_RESC)->matrix_mult(*LinSolver[level]->_RES, *RRamr[level]);
        *(LinSolver[level]->_RES) = * (LinSolver[level]->_RESC);
        LinSolver[level]->SwapMatrices();
        LinSolver[level]->_KK->matrix_ABC(*RRamr[level], *LinSolver[level]->_KKamr, *PPamr[level], false); // cannot use !firstNonlinearIt here
      }
      if(printdb) std::cout << "After KK amr restriction\n" << std::flush;
    }

    if(printdb) std::cout << "Before KK level" << level << " restriction\n" << std::flush;
    if(!RR[level]) { //Multilevel Restriction
      (LinSolver[level - 1u]->_RES)->matrix_mult_transpose(*LinSolver[level]->_RES, *PP[level]); // Resc = Pt Resf
      LinSolver[level - 1u]->_KK->matrix_PtAP(*PP[level], *LinSolver[level]->_KK, !firstNonlinearIt); // Kc = Pt Kf P // mat_reuse works only with level == levelF above
    }
    else {
      (LinSolver[level - 1u]->_RES)->matrix_mult(*LinSolver[level]->_RES, *RR[level]); // Resc = R Resf
      LinSolver[level - 1u]->_KK->matrix_ABC(*RR[level], *LinSolver[level]->_KK, *PP[level], !firstNonlinearIt); // Kc = R Kf P // mat_reuse works only with level == levelF above
    }
    if(printdb) std::cout << "After KK level" << level << " restriction\n" << std::flush;
  }

  std::cout << "Matrix Restriction time     = " << static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC << std::endl << std::flush;

  start_time = clock();

  if(printdb) std::cout << "Before KK sum \n" << std::flush;
  KK2->matrix_add (1., *LinSolver[level0 + level2]->_KK, "different_nonzero_pattern");
  *RES2 += *LinSolver[level0 + level2]->_RES;
  if(printdb) std::cout << "After KK sum \n" << std::flush;

  double tolerance = 0.;
  KK2->RemoveZeroEntries(tolerance);

  std::cout << "Matrix Clean Entry time     = " << static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC << std::endl << std::flush;

}
