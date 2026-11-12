#include <sys/stat.h>
#include <sys/types.h>

#include "FemusInit.hpp"

#include "LinearEquationSolverPetscFieldSplit.hpp"
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

// #include <gperftools/profiler.h>

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
#include "../includeLS/GhostPenaltySparsityPattern.hpp"

bool printdb = false;

// const RungeKutta::VelKind velocityType = RungeKutta::VelKind::Zero;
//const RungeKutta::VelKind velocityType = RungeKutta::VelKind::Translation;
//const RungeKutta::VelKind velocityType = RungeKutta::VelKind::Rotation;
//const RungeKutta::VelKind velocityType = RungeKutta::VelKind::Vortex;

#include "../includeLS/Utils.hpp"

std::vector<std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > > PMatrix(6); //projection matrix between elements;

CutFemIntegration cutFemIntegration;

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

void AssembleMultiphase(MultiLevelProblem& ml_prob);
void AssemblePsi(MultiLevelProblem& ml_prob);
// double TimeStepMultiphase(const double time);

int main(int argc, char **argv) {

  // Initialize PETSc/MPI
  FemusInit mpinit(argc, argv, MPI_COMM_WORLD);

  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  int iproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);

  // if (nprocs == 1)
  //   ProfilerStart("profiling.prof");

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
  const unsigned levelOffsetPsi = args.levelOffsetPsi;

  //simulation setting
  static const SimulationCase simulation = args.simulation;
  static const unsigned requestedDim = args.dim;
  const unsigned symmetric = args.symmetric;
  const unsigned axisymmetric = args.axisymmetric;

  cutFemIntegration = args.cutFemIntegration;

  const ElemType elementType = (requestedDim == 2) ? QUAD9 : HEX27;

  static const SimulationConfig cfg = GetSimulationConfig(simulation, requestedDim, symmetric, axisymmetric);

  MultiphasePhysicalProperties properties = cfg.properties;

  static std::vector<double> x_min;
  static std::vector<double> x_max;

  x_min = cfg.x_min;
  x_max = cfg.x_max;

  const double axisCoordinate = cfg.x_min[0];
  const std::vector<double> symmetryCoordinates = {cfg.x_min[0], cfg.x_min[1]};

  const double xmin = x_min[0];
  const double xmax = x_max[0];
  const double ymin = x_min[1];
  const double ymax = x_max[1];
  const double zmin = x_min[2];
  const double zmax = x_max[2];

  const unsigned nx = cfg.n[0];
  const unsigned ny = cfg.n[1];
  const unsigned nz = cfg.n[2];

  auto SetBoundaryCondition =
  [](const std::vector<double>& x, const char SolName[], double & value, const int /*facename*/, const double /*time*/) -> bool {

    const double tol = 1.e-10;
    const unsigned dim = cfg.boundary.min.size();

    int component = -1;
    if(!strcmp(SolName, "U")) component = 0;
    else if(!strcmp(SolName, "V")) component = 1;
    else if(!strcmp(SolName, "W")) component = 2;

    const bool pressure = !strcmp(SolName, "P1") || !strcmp(SolName, "P2") || !strcmp(SolName, "Psi");

    value = 0.;

    if(component < 0 && !pressure) return true;

    bool dirichlet = false;

    auto ApplyBC = [&](const BoundaryCondition & bc, const unsigned normalDir) {

      if(component >= 0) {

        switch(bc.type) {

          case BoundaryType::DIRICHLET_S:
            dirichlet = true;
            value = bc.value;
            break;

          case BoundaryType::NEUMANN_S:
            break;

          case BoundaryType::WALL:
            dirichlet = true;
            value = 0.;
            break;

          case BoundaryType::SLIP:
            if(static_cast<unsigned>(component) == normalDir) {
              dirichlet = true;
              value = 0.;
            }
            break;

          case BoundaryType::OUTFLOW:
            break;
        }

        return;
      }

      switch(bc.type) {

        case BoundaryType::DIRICHLET_S:
          dirichlet = true;
          value = bc.value;
          break;

        case BoundaryType::NEUMANN_S:
          value = bc.value;
          break;

        case BoundaryType::OUTFLOW:
          dirichlet = true;
          value = bc.value;
          break;

        case BoundaryType::WALL:
        case BoundaryType::SLIP:
          break;
      }
    };

    for(unsigned d = 0; d < dim; ++d) {

      if(x[d] <= cfg.x_min[d] + tol) ApplyBC(cfg.boundary.min[d], d);
      if(x[d] >= cfg.x_max[d] - tol) ApplyBC(cfg.boundary.max[d], d);
    }

    return dirichlet;
  };

  const double period = cfg.period;
  static double dt = 0.0;

  // interface settings
  const double r = cfg.psi.radius;
  std::vector<double> xc = cfg.psi.center;

  // variables name
  std::string psiName = "Psi";
  std::vector<std::string> vName = {"U", "V", "W"};
  std::vector<std::string> pName = {"P1", "P2"};
  std::string cName = "C";

  for (const unsigned numberOfSelectiveLevels : adaptiveLevelsList) {

    const double scalingFactor = 1.0;

    unsigned levelN = numberOfUniformLevels + numberOfSelectiveLevels;
    const unsigned levelF = levelN - 1u; //fine level associated for mlmsh0 and mlmsh1
    const unsigned levelPsi = levelF - levelOffsetPsi;
    const unsigned levelC = levelF - levelOffset; //coarse level associated to mlmsh2, but existing also mlmsh0 and mlmsh1
    const unsigned level0 = 0;

    for (const unsigned nSteps : nStepsList) {

      MultiLevelMesh mlMsh0;
      std::string meshName = "./input/cube_01_unstruct.neu";
      //mlMsh0.ReadCoarseMesh(meshName.c_str(), "seventh", scalingFactor);
      mlMsh0.GenerateCoarseBoxMesh(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, elementType, "seventh"); // Turek 1&2

      dt = period / nSteps;

      auto TimeStepMultiphase = [](const double time) {
        return dt;
      };
      auto firstTimeStepMultiphase = [](const double time) {
        return 0.5 * dt;
      };

      if(args.nPrint == 0u) {
        throw std::runtime_error("--nprint must be > 0");
      }

      if(args.nReinit == 0u) {
        throw std::runtime_error("--nreinit must be > 0");
      }

      const unsigned nPrint = std::max(1u, nSteps / args.nPrint);
      const unsigned nReinit = std::max(1u, nSteps / args.nReinit);

      mlMsh0.RefineMesh(numberOfUniformLevels, numberOfUniformLevels, nullptr);

      unsigned dim = mlMsh0.GetDimension();

      if(dim == 3) {
        auto refineElement = new RefineElement(1, "hex", "linear", "fifth", "fifth", "legendre");
        PMatrix[0] = refineElement->GetProlongationMatrix();
        delete refineElement;

        refineElement = new RefineElement(1, "tet", "linear", "fifth", "fifth", "legendre");
        PMatrix[1] = refineElement->GetProlongationMatrix();
        delete refineElement;

        refineElement = new RefineElement(1, "wedge", "linear", "fifth", "fifth", "legendre");
        PMatrix[2] = refineElement->GetProlongationMatrix();
        delete refineElement;
      }
      else if (dim == 2) {
        auto refineElement = new RefineElement(1, "quad", "linear", "fifth", "fifth", "legendre");
        PMatrix[3] = refineElement->GetProlongationMatrix();
        delete refineElement;

        refineElement = new RefineElement(1, "tri", "linear", "fifth", "fifth", "legendre");
        PMatrix[4] = refineElement->GetProlongationMatrix();
        delete refineElement;
      }

      xc.resize(dim);

      unsigned sigmoidType = 0;
      double eps = 0.25; //(dim == 2) ? 1. / pow(2, std::max(levelN - 7u, 1u))
      Mollifier m = Mollifier(eps, sigmoidType);

      Psi psi(cfg.psi);

      // Iteratively flag and create new AMR levels
      for (unsigned k = 0; k < numberOfSelectiveLevels; ++k) {
        FlagFinestMeshLevel(mlMsh0, psi);
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

      InitLevelSet(mlSol0, psiName, psi);
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
      mlmsh1->GenerateCoarseBoxMesh(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, elementType, "seventh"); // Turek 1&2
      // mlmsh1->ReadCoarseMesh(meshName.c_str(), "seventh", scalingFactor);
      mlmsh1->RefineMesh(numberOfUniformLevels, numberOfUniformLevels, nullptr);

      UpdateColorFunction(*mlsol0, psiName, cName);
      if(levelC < levelF) RestrictPWDCField(*mlsol0, cName, levelC, levelF);

      mlsol0->AttachSetBoundaryConditionFunction(SetBoundaryCondition);
      mlsol0->GenerateBdc("All");

      LevelSetDiagnostics diagnostics = ComputeLevelSetDiagnostics(*mlsol0, psiName, vName, pName, levelC, simulation, cfg, symmetric, axisymmetric, axisCoordinate, symmetryCoordinates);

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

        if(simulation == SimulationCase::RB1 || simulation == SimulationCase::RB2) {

          for(unsigned d = 0; d < dim; ++d) {
            out << std::setw(w) << ("barycenter" + std::to_string(d));
          }

          for(unsigned d = 0; d < dim; ++d) {
            out << std::setw(w) << ("meanVelocity" + std::to_string(d));
          }

          out << std::setw(w) << "circularity";
        }

        if(simulation == SimulationCase::RT) {
          out << std::setw(w) << "pb"
              << std::setw(w) << "ps"
              << std::setw(w) << "psad";
        }

        if(simulation == SimulationCase::PI) {
          out << std::setw(w) << "fluid1UpperVolume"
              << std::setw(w) << "fluid1UpperRatio"
              << std::setw(w) << "fluid1VolumeRatio"
              << std::setw(w) << "potentialEnergy1"
              << std::setw(w) << "potentialEnergy2"
              << std::setw(w) << "kineticEnergy1"
              << std::setw(w) << "kineticEnergy2"
              << std::setw(w) << "enstrophy1"
              << std::setw(w) << "enstrophy2"
              << std::setw(w) << "mechanicalEnergy"
              << std::setw(w) << "potentialEnergy1Star"
              << std::setw(w) << "potentialEnergy2Star"
              << std::setw(w) << "kineticEnergy1Star"
              << std::setw(w) << "kineticEnergy2Star"
              << std::setw(w) << "interfaceAreaStar";
        }

        out << '\n';
      }

      // PrintLevelSetDiagnostics(diagnostics, iproc, 0.0, diagnosticsFile);

      const std::string ls_outputdir = "output_ls";
      const std::string vel_outputdir = "output_vel";

      VTKWriter vtkIO(&mlSol0);
      vtkIO.Write(levelF, ls_outputdir, "biquadratic", variablesToBePrinted, 0);

      LevelSetDiagnostics final_diagnostics;

      BurstingBubbleState bbState;
      BurstingBubbleOptions bbOptions;
      bbOptions.liquidPositive = false;
      bbOptions.minRadiusOverH = 2.;
      bbOptions.minEquivalentRadius = 0.;
      bbOptions.minDropletVolume = 0.;
      bbOptions.graphLiquidFractionTolerance = 1.e-12;
      bbOptions.persistenceSteps = 2u;
      bbOptions.freeSurfaceZ = cfg.psi.freeSurfaceHeight;
      bbOptions.axisTolerance = -1.;
      bbOptions.psiZeroTolerance = 1.e-12;
      bbOptions.dropletTrackingTolerance = -1.;
      const double bbRadius = 1.;
      const double bbRhoLiquid = properties.rho2;
      const double bbMuLiquid = properties.mu2;
      const double bbSigma = properties.sigma;
      const double bbHMin = GetMaxElementH(mlmsh0->GetLevel(levelF), levelF);
      const std::string bbOutputPrefix = "levelsetdiagnostic/bursting_bubble_LC" + std::to_string(levelC) + "_LF" + std::to_string(levelF) + "_NS" + std::to_string(nSteps);
      if(iproc == 0) {
        std::ofstream(bbOutputPrefix + "_jet.dat", std::ios::trunc).close();
        std::ofstream(bbOutputPrefix + "_droplets.dat", std::ios::trunc).close();
      }

      for (unsigned t = 1; t <= 1 + nSteps; t++) {

        std::cout << "==============================================\n";
        std::cout << "  ITERATION " << t << std::endl;
        std::cout << "==============================================\n";

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
            unsigned P0index = mlsol0->GetIndex(pName[d].c_str());
            unsigned P2index = mlSol2.GetIndex(pName[d].c_str());
            *(Sol2[P2index]) = *(Sol0_C[P0index]);
          }
        }

        // mlProb3 is used to assemble the auxiliary level-set
        MultiLevelMesh mlMsh3(*mlmsh0, level0, levelPsi + 1, "seventh");
        std::vector<Mesh*> msh3(levelPsi + 1 - level0);
        for(unsigned l = 0; l < msh3.size(); l++)
          msh3[l] = mlMsh3.GetLevel(l);
        MultiLevelSolution mlSol3(&mlMsh3);
        mlSol3.AddSolution(psiName.c_str(), LAGRANGE, SECOND, 0, true);

        mlSol3.Initialize("All");

        mlSol3.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
        mlSol3.GenerateBdc("All");

        auto &Sol0_Psi = (mlsol0->GetSolutionLevel(levelPsi))->_Sol;
        auto &Sol3 = (mlSol3.GetSolutionLevel(levelPsi - level0))->_Sol;

        const unsigned c0Index = mlsol0->GetIndex(cName.c_str());
        const unsigned c2Index = mlSol2.GetIndex(cName.c_str());

        auto originalC = NumericVector::build();

        originalC->init(*Sol0_C[c0Index]);

        *originalC = *Sol0_C[c0Index];

        ReclassifyUnsupportedPressureCells(mlSol2, cName, levelC - level0, levelC);

        if(level0 < levelC) {
          RestrictPWDCField(mlSol2, cName, 0, levelC - level0);
        }

        *Sol0_C[c0Index] = *Sol2[c2Index];

        std::vector<double> xtarget = (simulation == SimulationCase::PI) ?
                                      std::vector<double> {0.075, 0.075, 0.} :
                                      std::vector<double> {xmin, ymin, zmin};
        xtarget.resize(dim);
        SetUnphysicalPressureDofs(mlSol2, cName, pName, 0, levelC - level0, xtarget, true);

        MultiLevelProblem mlProb2(&mlSol2);
        MultiLevelProblem mlProb3(&mlSol3);

        // add system Navier-Stokes in mlProb as a Linear Implicit System
        TransientNonlinearImplicitSystem& system2 = mlProb2.add_system < TransientNonlinearImplicitSystem > ("NS");
        LinearImplicitSystem& system3 = mlProb3.add_system < LinearImplicitSystem > ("PSI");

        // add velocity to system
        for(unsigned d = 0; d < dim; d++) system2.AddSolutionToSystemPDE(vName[d].c_str());
        //add pressure
        for(unsigned d = 0; d < pName.size(); d++) system2.AddSolutionToSystemPDE(pName[d].c_str());

        system3.AddSolutionToSystemPDE(psiName.c_str());

        auto extraSparsityFunction = [&](
                                       LinearEquation & linEq,
                                       vector < std::map < int, bool > >& BlgToMe_d,
                                       vector < std::map < int, bool > >& BlgToMe_o,
                                       std::map < int, std::map <int, bool > >& DnBlgToMe_o,
                                       std::map < int, std::map <int, bool > >& DnBlgToMe_d
        ) {
          if (linEq._msh->GetLevel() != levelC - level0) {
            return;
          }

          for (unsigned idim = 0; idim < dim; idim++) {
            GetGPSparsityPattern (mlProb2, cName, vName[idim], levelC - level0, levelC, 0, BlgToMe_d, BlgToMe_o, DnBlgToMe_o, DnBlgToMe_d);
          }
          for (unsigned idim = 0; idim < pName.size(); idim++) {
            GetGPSparsityPattern (mlProb2, cName, pName[idim], levelC - level0, levelC, idim + 1, BlgToMe_d, BlgToMe_o, DnBlgToMe_o, DnBlgToMe_d);
          }

        };

        system2.SetExtraSparsityFunction(extraSparsityFunction);

        // attach the assembling function to system
        system2.SetAssembleFunction(AssembleMultiphase);
        system2.AttachGetTimeIntervalFunction(
          (t == 1) ? firstTimeStepMultiphase : TimeStepMultiphase
        );
        // initilaize and solve the system

        system3.SetAssembleFunction(AssemblePsi);

        //system2.SetOuterSolver(PREONLY);
        //system2.SetMaxNumberOfNonLinearIterations(10);

        bool linearScheme = true;

        mlProb2.SetMultiphaseParams(mlsol0, levelPsi, levelC, level0, hC, properties, td, linearScheme, axisymmetric, axisCoordinate);
        mlProb3.SetMultiphaseParams(mlsol0, levelF, levelPsi, level0, hC, properties, td, linearScheme, axisymmetric, axisCoordinate);

        for(unsigned l = 0; l < msh2.size(); l++)
          msh2[l]->SetLevel(l);

        for(unsigned l = 0; l < msh3.size(); l++)
          msh3[l]->SetLevel(l);

        std::vector < unsigned > fieldV(dim);
        for(unsigned d = 0; d < dim; d++) fieldV[d] = system2.GetSolPdeIndex(vName[d].c_str());

        FieldSplitTree FS_V(PREONLY, MULTIGRID_PRECOND, fieldV, "Velocity");
        FS_V.MgSetSolverFineGrids(CHEBYSHEV);
        FS_V.MgSetPreconditionerFineGrids (JACOBI_PRECOND);
        FS_V.MgSetParameters (4, 4, 0);

        std::vector < unsigned > fieldP(pName.size());
        for(unsigned d = 0; d < pName.size(); d++) fieldP[d] = system2.GetSolPdeIndex(pName[d].c_str());

        FieldSplitTree FS_P(PREONLY, MULTIGRID_PRECOND, fieldP, "Pressure");
        FS_P.MgSetSolverFineGrids(CHEBYSHEV);
        FS_P.MgSetPreconditionerFineGrids (JACOBI_PRECOND);
        FS_P.MgSetParameters (4, 4, 0);

        std::vector < FieldSplitTree *> FS1;
        FS1.reserve(2);
        FS1.push_back(&FS_V);
        FS1.push_back(&FS_P);

        FieldSplitTree FS_NS(GMRES, FIELDSPLIT_SCHUR_PRECOND, FS1, "Navier-Stokes");
        FS_NS.SetTolerances(1.e-8, 1.e-12, 1.e+50, 10);

        CoarseLevelInclude includeCoarseLevel = INCLUDE_COARSE_LEVEL_TRUE;

        system2.SetLinearEquationSolverType(FEMuS_FIELDSPLIT, includeCoarseLevel);

        // attach the assembling function to system
        system2.SetMaxNumberOfNonLinearIterations((linearScheme) ? 1 : 3);
        system2.SetMaxNumberOfLinearIterations(10);
        system2.SetMgType(V_CYCLE);

        // initilaize and solve the system
        system2.init();

        system2.SetOuterSolver(FGMRES);
        system2.SetFieldSplitTree(&FS_NS);
        system2.SetTolerances(1.e-8, 1.e-12, 1.e+50, 50);
        system2.SetAbsoluteLinearConvergenceTolerance(1.e-8);

        system3.SetMaxNumberOfLinearIterations(10);
        system3.SetMgType(V_CYCLE);

        // initilaize and solve the system
        if (levelPsi < levelF) {
          system3.init();

          system3.SetOuterSolver(GMRES);
          system3.SetTolerances(1.e-8, 1.e-12, 1.e+50, 50);
          system3.SetAbsoluteLinearConvergenceTolerance(1.e-8);

          PetscOptions psiOptions;
          char* originalOptions = nullptr;

          PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsCreate(&psiOptions));
          PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsGetAll(nullptr, &originalOptions));
          PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsInsertString(psiOptions, originalOptions));
          PetscCallAbort(PETSC_COMM_WORLD, PetscFree(originalOptions));

          // PSI outer solver
          PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsSetValue(psiOptions, "-ksp_type", "gmres"));
          PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsSetValue(psiOptions, "-ksp_rtol", "1.e-8"));
          PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsSetValue(psiOptions, "-ksp_atol", "1.e-12"));
          PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsSetValue(psiOptions, "-ksp_max_it", "50"));

          // PSI multigrid smoothers
          for(unsigned l = 1; l <= levelPsi - level0; ++l) {
            const std::string option = "-level-" + std::to_string(l) + "ksp_max_it";
            PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsSetValue(psiOptions, option.c_str(), "1"));
          }

          // Activate PSI-specific PETSc options
          PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsPush(psiOptions));

          system3.MGsolve();

          // Restore the original PETSc options for NS
          PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsPop());
          PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsDestroy(&psiOptions));

          unsigned psi0index = mlsol0->GetIndex(psiName.c_str());
          unsigned psi3index = mlSol3.GetIndex(psiName.c_str());
          *(Sol0_Psi[psi0index]) = *(Sol3[psi3index]);
        }

        for (unsigned l = 0; l < levelC + 1 - level0; l++) {
          LinearEquationSolver* pdeSys2_base  = system2._LinSolver[l];
          pdeSys2_base->SetSolution(&mlSol2.GetSolutionLevel(l)->_Sol);
          // pdeSys2_base->MergeNullSpaceBases(true);
        }

        //system2.MGsolve();
        system2.Solve();

        *Sol0_C[c0Index] = *originalC;

        for(unsigned l = 0; l < msh2.size(); l++)
          msh2[l]->SetLevel(l + level0);

        for(unsigned l = 0; l < msh3.size(); l++)
          msh3[l]->SetLevel(l + level0);

        VTKWriter vtkIO2(&mlSol2);
        VTKWriter vtkIO3(&mlSol3);
        vtkIO2.SetDebugOutput(true);
        if (t % nPrint == 0) {
          vtkIO2.Write(levelC - level0, vel_outputdir, "biquadratic", variablesToBePrinted, t / nPrint);
          vtkIO3.Write(levelPsi - level0, ls_outputdir, "biquadratic", variablesToBePrinted, t / nPrint);
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

          final_diagnostics = ComputeLevelSetDiagnostics(*mlsol0, psiName, vName, pName, levelC, simulation, cfg, symmetric, axisymmetric, axisCoordinate, symmetryCoordinates);

          VTKWriter vtkIO(mlsol0);
          vtkIO.Write(levelF, ls_outputdir, "biquadratic", variablesToBePrinted, t / nPrint);

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

        if (t % nReinit == 0) {
          Reinit reinit(psiName, *mlsol0, m);

          reinit.farFieldReinit(X0);
          reinit.interfaceFieldReinit(bbox);
          reinit.updateSolution();
        }

        // if (t == 1)
        //   WritePointsVTK("./output/points.0.vtk", X0);

        RungeKutta4(X0, *mlsol0, bbox, x_min, x_max, vName, levelC, dt); // move the interface points forward in time using the velocity mls0(lC)

        // if (t % 1 == 0)
        //   WritePointsVTK("./output/points." + std::to_string(t / 1) + ".vtk", X0);

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

        ProjectSolution(*mlsol0, *mlsol1, bbox, x_min, x_max, {psiName}, levelF, levelF, vName, levelC, zero_bd, -dt, time, period);
        ProjectSolution(*mlsol0, *mlsol1, bbox, x_min, x_max, vName, levelC, levelC);
        ProjectSolution(*mlsol0, *mlsol1, bbox, x_min, x_max, pName, levelC, levelC);

        UpdateColorFunction(*mlsol1, psiName, cName);
        if(levelC < levelF) RestrictPWDCField(*mlsol1, cName, levelC, levelF);

        for(unsigned i = 0; i < cfw.size(); i++) cfw[i]->ClearMap();
        if(cutFemIntegration == CutFemIntegration::Quadratic) {
          for(auto* rule : cfwLS) {
            if(rule != nullptr) rule->ClearCache();
          }
        }

        LevelSetDiagnostics diagnostics = ComputeLevelSetDiagnostics(*mlsol1, psiName, vName, pName, levelC, simulation, cfg, symmetric, axisymmetric, axisCoordinate, symmetryCoordinates);
        PrintLevelSetDiagnostics(diagnostics, iproc, time, diagnosticsFile);

        if(simulation == SimulationCase::BB) {
          const BurstingBubbleMetrics bbMetrics = ComputeBurstingBubbleMetrics(*mlsol1, psiName, vName, levelC, time, bbHMin, bbRadius, bbRhoLiquid, bbMuLiquid, bbSigma, axisCoordinate, bbState, bbOptions);
          PrintBurstingBubbleMetrics(bbMetrics, iproc, bbOutputPrefix);
        }

        // Export solution to VTK (selected levels)
        VTKWriter vtkIO1(mlsol1);
        if (t % nPrint == 0) {
          vtkIO1.Write(levelF, ls_outputdir, "biquadratic", variablesToBePrinted, t / nPrint);
        }

        std::swap(mlsol0, mlsol1);
        std::swap(mlmsh0, mlmsh1);

        mlsol1->clear();
        mlmsh1->resize(numberOfUniformLevels);

        for(unsigned i = 0; i < cfw.size(); i++) cfw[i]->clear();

      }

    }

  }

  // if (nprocs == 1)
  //   ProfilerStop();
  return 0;
}

void postRecv(void* data_recv, unsigned count_recv, MPI_Datatype datatype, unsigned source, int tag, std::vector<MPI_Request>& requests) {
  MPI_Request request = MPI_REQUEST_NULL;
  MPI_Irecv(data_recv, static_cast<int>(count_recv), datatype, static_cast<int>(source), tag, MPI_COMM_WORLD, &request);
  requests.push_back(request);
}

void postSend(const void* data_send, unsigned count_send, MPI_Datatype datatype, unsigned destination, int tag, std::vector<MPI_Request>& requests) {
  if (count_send == 0) {
    return;
  }
  MPI_Request request = MPI_REQUEST_NULL;
  MPI_Isend(data_send, static_cast<int>(count_send), datatype, static_cast<int>(destination), tag, MPI_COMM_WORLD, &request);
  requests.push_back(request);
}

void waitAll(std::vector<MPI_Request>& requests) {
  if (!requests.empty()) {
    MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
    requests.clear();
  }
}

void AssembleMultiphase(MultiLevelProblem& ml_prob2) {

  //  extract pointers to the several objects that we are going to use
  TransientNonlinearImplicitSystem* mlPdeSys   = &ml_prob2.get_system<TransientNonlinearImplicitSystem> ("NS");
  const unsigned levelToAssemble = mlPdeSys->GetLevelToAssemble();

  MultiLevelSolution*  mlSol    = ml_prob2._ml_sol;
  LinearEquationSolver* pdeSys  = mlPdeSys->_LinSolver[levelToAssemble];
  SparseMatrix*    KK           = pdeSys->_KK;
  NumericVector*   RES          = pdeSys->_RES;

  Mesh* msh = ml_prob2._ml_msh->GetLevel(levelToAssemble);    // pointer to the mesh (level) object
  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  // Stuff to be taken/used on level2 from system2 and mlsol
  std::vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("V");    // get the position of "V" in the ml_sol object
  if(dim == 3) solVIndex[2] = mlSol->GetIndex("W");       // get the position of "V" in the ml_sol object
  unsigned solVType = mlSol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

  unsigned solP1Index = mlSol->GetIndex("P1");    // get the position of "P1" in the ml_sol object
  unsigned solP2Index = mlSol->GetIndex("P2");    // get the position of "P2" in the ml_sol object
  unsigned solPType = mlSol->GetSolutionType(solP1Index);    // get the finite element type for "u"

  std::vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("U");    // get the position of "U" in the pdeSys object
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("V");    // get the position of "V" in the pdeSys object
  if(dim == 3) solVPdeIndex[2] = mlPdeSys->GetSolPdeIndex("W");

  unsigned solP1PdeIndex = mlPdeSys->GetSolPdeIndex("P1");    // get the position of "P" in the pdeSys object
  unsigned solP2PdeIndex = mlPdeSys->GetSolPdeIndex("P2");    // get the position of "P" in the pdeSys object

  KK->zero();
  RES->zero();

  if (solPType != mlSol->GetSolutionType(solP2Index)) {
    throw std::runtime_error("Pressure type mismatch");
  }

  // AssembleGhostPenaltyVelocity(ml_prob2);

  if (solPType == 3) {
    AssembleGhostPenaltyDGP(ml_prob2, true);
    AssembleGhostPenaltyDGP(ml_prob2, false);
  }
  else if (solPType == 0) {
    AssembleGhostPenaltyLinearPressure(ml_prob2, true);
    AssembleGhostPenaltyLinearPressure(ml_prob2, false);
  }
  else {
    throw std::runtime_error("Pressure type must be 3 or 0");
  }

  AssembleStabilizationTerms(ml_prob2);

  //RES->close();
  //KK->flush();

  MultiphaseParams mParam = ml_prob2.GetMultiphaseParams();
  const unsigned levelF = mParam.levelF;
  const unsigned levelC = mParam.levelC;
  const unsigned level0 = mParam.level0;

  const TimeDiscretization td = mParam.td;
  const double cold = (td == TimeDiscretization::CrankNicholson) ? 0.5 : 0.;
  const double cnew = (td == TimeDiscretization::CrankNicholson) ? 0.5 : 1.;

  std::cout << "levelC = " << levelC << " levelF = " << levelF << std::endl;
  std::cout << "level to assemble = " << levelToAssemble  << " mapping level to assemble to levelC = " << level0 + levelToAssemble << std::endl;

  MultiphasePhysicalProperties properties = mParam.properties;
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
  bool linear = mParam.linear;
  bool axisymmetric = mParam.axisymmetric;
  double axisCoordinate = mParam.axisCoordinate;

  const bool firstNonlinearIt = (mlPdeSys->GetNonlinearIt() == 0);

  // Stuff to be taken/used from levelC to levelF from mlsol0
  MultiLevelSolution*  mlSol0  = static_cast<MultiLevelSolution*>(mParam.extProblem);  // pointer to the multilevelF solution object
  unsigned cIndex = mlSol0->GetIndex("C");
  unsigned psiIndex = mlSol0->GetIndex("Psi");
  unsigned psiType = mlSol0->GetSolutionType("Psi");

  std::vector < std::vector < double > >  solV(dim);    // local solution
  std::vector < std::vector < double > >  solVOld(dim);    // local solution
  std::vector < double >  solP1; // local solution
  std::vector < double >  solP2; // local solution

  std::vector < double >  psi; // local solution
  std::vector < double >  k;
  std::vector < std::vector < double > >  n(dim);

  std::vector < std::vector < double > > coordXi(dim);
  std::vector < std::vector < double > > coordX(dim); // local coordinates
  std::vector < std::vector < double > > coordXgf(dim);    // local coordinates
  unsigned solXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
  unsigned solXiType = 0; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector <double> phiV;  // local test function for velocity
  std::vector <double> phiV_x; // local test function first order partial derivatives

  std::vector <double> phiPsi;
  std::vector <double> phiPsi_x;

  std::vector <double> phiXi;
  std::vector <double> phiXi_x;

  unsigned dim2 = 3 * (dim - 1);

  std::vector<double> phiP; // local test function for the pressure
  double weight; // gauss point weight
  double weightFake;

  std::vector< unsigned > sysDof; // local to global pdeSys dofs
  std::vector< double > Res; // local redidual std::vector
  std::vector < double > Jac;

  double eps = 0.;

  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  // if (nprocs == 1)
  //   ProfilerStart("profiling.prof");

  if(printdb) std::cout << "After Solution Projection\n" << std::flush;

  // AssembleStabilizationTerms(*ml_prob0);

  double epsilonp = 0;

  clock_t start_time = clock();

  // element loop: each process loops only on the elements that owns

  if(printdb) std::cout << "Before KK assembly\n" << std::flush;

  double P[3][3] = {};
  double solV_gss[3] = {};
  double solVOld_gss[3] = {};
  double gradSolV_gss[3][3] = {};
  double gradSolVOld_gss[3][3] = {};

  const unsigned nElementTypes = 6;
  const unsigned nSolutionTypes = 5;
  // const unsigned nChildren = 1u << dim; // we dont assume that if C==0.5 then the element is refined. we check in the next loop

  const int tagCounts = 100;
  const int tagV = 101;
  const int tagP = 102;
  const int tagX = 103;
  const int tagXi = tagX + static_cast<int>(dim);
  const int tagElCounts = tagXi + static_cast<int>(dim);
  const int tagEl = tagElCounts + 1;
  // const int tagElGF = tagEl + 1;
  const int tagSolV = tagEl + 1;
  const int tagSolVOld = tagSolV + static_cast<int>(dim);
  const int tagSolP1 = tagSolVOld + static_cast<int>(dim);
  const int tagSolP2 = tagSolP1 + 1;

  std::vector<MPI_Request> requests;
  std::vector<unsigned> el_recv;
  // std::vector<unsigned> el_gf_recv;
  std::vector<unsigned> dofV_recv;
  std::vector<unsigned> dofP_recv;
  std::vector<std::vector<double>> X_recv;
  std::vector<std::vector<double>> Xi_recv;
  std::vector<std::vector<double>> solV_recv;
  std::vector<std::vector<double>> solVOld_recv;
  std::vector<double> solP1_recv;
  std::vector<double> solP2_recv;

  std::vector<unsigned> offsetV;
  std::vector<unsigned> offsetP;
  std::vector<unsigned> offsetX;
  std::vector<unsigned> offsetXi;

  Solution* sol = mlSol->GetLevel(levelToAssemble);

  start_time = clock();

  double A = 0;
  double Per = 0;

  for (unsigned l = levelC; l <= levelF; ++l) {
    Solution* sol0 = mlSol0->GetLevel(l);

    Mesh* msh = sol0->GetMesh();
    elem* el = msh->GetElementArray();

    const unsigned nChildrenRefined = msh->GetRefIndex();

    auto GetNumberOfChildren = [&](const unsigned iel) -> unsigned {
      return msh->GetRefinedElementIndex(iel) ? nChildrenRefined : 1u;
    };

    auto solC = sol0->_Sol[cIndex];
    auto solPsi = sol0->_Sol[psiIndex];

    unsigned    iproc = msh->processor_id();

    const unsigned nElementsToLoop = (l == levelC) ? msh->_elementOffset[iproc + 1] - msh->_elementOffset[iproc] : static_cast<unsigned>(el_recv.size());

    std::vector<unsigned> el_to_assemble(nElementsToLoop, 1);
    offsetV.resize(nElementsToLoop + 1, 0);
    offsetP.resize(nElementsToLoop + 1, 0);
    offsetX.resize(nElementsToLoop + 1, 0);
    offsetXi.resize(nElementsToLoop + 1, 0);

    for (unsigned k = 0; k < nElementsToLoop; ++k) {
      const unsigned iel = (l == levelC) ? msh->_elementOffset[iproc] + k : el_recv[k];

      if (l < levelF) {
        const double C = (*solC)(iel);
        el_to_assemble[k] = !(C > 0.25 && C < 0.75);
      }

      offsetV[k + 1] = offsetV[k] + msh->GetElementDofNumber(iel, solVType);
      offsetP[k + 1] = offsetP[k] + msh->GetElementDofNumber(iel, solPType);
      offsetX[k + 1] = offsetX[k] + msh->GetElementDofNumber(iel, solXType);
      offsetXi[k + 1] = offsetXi[k] + msh->GetElementDofNumber(iel, solXiType);
    }

    for (unsigned k = 0; k < nElementsToLoop; ++k) {
      if (el_to_assemble[k] == 0) {
        continue;
      }

      const unsigned iel = (l == levelC) ? msh->_elementOffset[iproc] + k : el_recv[k];
      // const unsigned iel_gf = (l == levelC) ? iel : el_gf_recv[k];
      const unsigned iel_type = msh->GetElementType(iel);
      const unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);
      const unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);
      const unsigned nDofsX = msh->GetElementDofNumber(iel, solXType);
      const unsigned nDofsXi = msh->GetElementDofNumber(iel, solXiType);

      // Existing element assembly goes here; use local fields at levelC and received fields above levelC.
    }

    Solution* sol0p = (l < levelF) ? mlSol0->GetLevel(l + 1) : nullptr;
    Mesh* mshp = (l < levelF) ? sol0p->GetMesh() : nullptr;

    std::vector<unsigned> n_el_send(nprocs, 0);
    std::vector<std::vector<unsigned>> n_el_to_interpolate_send(nprocs, std::vector<unsigned>(nElementTypes, 0));

    // for (unsigned k = 0; k < nElementsToLoop; ++k) {
    //   if (el_to_assemble[k] != 0) {
    //     continue;
    //   }

    //   const unsigned iel = (l == levelC) ? msh->_elementOffset[iproc] + k : el_recv[k];
    //   const unsigned iel_type = msh->GetElementType(iel);

    //   for (unsigned j = 0; j < nChildren; ++j) {
    //     const unsigned jel = el->GetChild(iel, j);
    //     const unsigned jproc = mshp->IsdomBisectionSearch(jel, 3);
    //     ++n_el_to_interpolate_send[jproc][iel_type];
    //     ++n_el_send[jproc];
    //   }
    // }

    for(unsigned k = 0; k < nElementsToLoop; ++k) {

      if(el_to_assemble[k] != 0) {
        continue;
      }

      const unsigned iel = (l == levelC) ? msh->_elementOffset[iproc] + k : el_recv[k];

      const unsigned iel_type = msh->GetElementType(iel);
      const unsigned nChildrenCurrent = GetNumberOfChildren(iel);

      for(unsigned j = 0; j < nChildrenCurrent; ++j) {

        const unsigned jel = el->GetChild(iel, j);
        const unsigned jproc = mshp->IsdomBisectionSearch(jel, 3);

        ++n_el_to_interpolate_send[jproc][iel_type];
        ++n_el_send[jproc];
      }

    }

    std::vector<std::vector<unsigned>> n_dof_send(nprocs, std::vector<unsigned>(nSolutionTypes, 0));

    for (unsigned jproc = 0; jproc < nprocs; ++jproc) {
      for (unsigned soltype = 0; soltype < nSolutionTypes; ++soltype) {
        for (unsigned ieltype = 0; ieltype < nElementTypes; ++ieltype) {
          n_dof_send[jproc][soltype] += n_el_to_interpolate_send[jproc][ieltype] * NVE[ieltype][soltype];
        }
      }
    }

    std::vector<std::vector<unsigned>> el_send(nprocs);
    // std::vector<std::vector<unsigned>> el_gf_send(nprocs);
    std::vector<std::vector<unsigned>> dofV_send(nprocs);
    std::vector<std::vector<unsigned>> dofP_send(nprocs);
    std::vector<std::vector<double>> X_send(dim * nprocs);
    std::vector<std::vector<double>> Xi_send(dim * nprocs);
    std::vector<std::vector<double>> solV_send(dim * nprocs);
    std::vector<std::vector<double>> solVOld_send(dim * nprocs);
    std::vector<std::vector<double>> solP1_send(nprocs);
    std::vector<std::vector<double>> solP2_send(nprocs);

    for (unsigned jproc = 0; jproc < nprocs; ++jproc) {
      el_send[jproc].resize(n_el_send[jproc]);
      // el_gf_send[jproc].resize(n_el_send[jproc]);
      dofV_send[jproc].resize(n_dof_send[jproc][solVType]);
      dofP_send[jproc].resize(n_dof_send[jproc][solPType]);
      solP1_send[jproc].resize(n_dof_send[jproc][solPType]);
      solP2_send[jproc].resize(n_dof_send[jproc][solPType]);

      for (unsigned d = 0; d < dim; ++d) {
        X_send[d * nprocs + jproc].resize(n_dof_send[jproc][solXType]);
        Xi_send[d * nprocs + jproc].resize(n_dof_send[jproc][solXiType]);
        solV_send[d * nprocs + jproc].resize(n_dof_send[jproc][solVType]);
        solVOld_send[d * nprocs + jproc].resize(n_dof_send[jproc][solVType]);
      }
    }

    std::vector<unsigned> posEl_send(nprocs, 0);
    std::vector<unsigned> posV_send(nprocs, 0);
    std::vector<unsigned> posP_send(nprocs, 0);
    std::vector<unsigned> posX_send(nprocs, 0);
    std::vector<unsigned> posXi_send(nprocs, 0);

    for (unsigned k = 0; k < nElementsToLoop; ++k) {

      const unsigned iel = (l == levelC) ? msh->_elementOffset[iproc] + k : el_recv[k];
      // const unsigned iel_gf = (l == levelC) ? iel : el_gf_recv[k];
      const unsigned iel_type = msh->GetElementType(iel);
      const unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);
      const unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);
      const unsigned nDofsX = msh->GetElementDofNumber(iel, solXType);
      const unsigned nDofsXi = msh->GetElementDofNumber(iel, solXiType);

      for(unsigned  k = 0; k < dim; k++) {
        coordXi[k].resize(nDofsXi);
      }
      for (unsigned i = 0; i < nDofsXi; ++i) {
        if (l == levelC) {
          const double* xi_pt = msh->_finiteElement[iel_type][solXiType]->GetBasis()->GetXcoarse(i);
          for (unsigned d = 0; d < dim; ++d) {
            coordXi[d][i] = xi_pt[d];
          }
        }
        else {
          for (unsigned d = 0; d < dim; ++d) {
            coordXi[d][i] = Xi_recv[d][offsetXi[k] + i];
          }
        }
      }

      if (el_to_assemble[k] == 1) {

        //here is where the assembly happens
        double C = (*solC)(iel);
        short unsigned ielGeom = iel_type;
        unsigned nDofsVP = dim * nDofsV + 2 * nDofsP;

        //resize local arrays
        sysDof.resize(nDofsVP);
        Res.assign(nDofsVP, 0.);
        Jac.assign(nDofsVP * nDofsVP, 0.);

        for(unsigned  d = 0; d < dim; d++) {
          solV[d].resize(nDofsV);
          solVOld[d].resize(nDofsV);
          coordX[d].resize(nDofsX);
          coordXgf[d].resize(nDofsX);
        }
        solP1.resize(nDofsP);
        solP2.resize(nDofsP);

        //local storage of global mapping and solution
        if (l == levelC) {
          for(unsigned i = 0; i < nDofsV; i++) {
            unsigned solVDof =  msh->GetSolutionDof(i, iel, solVType);
            for(unsigned  d = 0; d < dim; d++) {
              solV[d][i] = (*sol->_Sol[solVIndex[d]])(solVDof);
              solVOld[d][i] = (*sol->_SolOld[solVIndex[d]])(solVDof);
              sysDof[d * nDofsV + i] = pdeSys->GetSystemDof(solVIndex[d], solVPdeIndex[d], i, iel);
            }
          }
          for(unsigned i = 0; i < nDofsP; i++) {
            unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);
            solP1[i] = (*sol->_Sol[solP1Index])(solPDof);
            solP2[i] = (*sol->_Sol[solP2Index])(solPDof);
            sysDof[dim * nDofsV + i ] = pdeSys->GetSystemDof(solP1Index, solP1PdeIndex, i, iel);
            sysDof[dim * nDofsV + nDofsP + i ] = pdeSys->GetSystemDof(solP2Index, solP2PdeIndex, i, iel);
          }
        }
        else {
          for(unsigned i = 0; i < nDofsV; i++) {
            unsigned solVDof = dofV_recv[offsetV[k] + i];
            for(unsigned  d = 0; d < dim; d++) {
              solV[d][i] = solV_recv[d][offsetV[k] + i];
              solVOld[d][i] = solVOld_recv[d][offsetV[k] + i];
              sysDof[d * nDofsV + i] = pdeSys->GetSystemDof(solVIndex[d], solVPdeIndex[d], solVDof);
            }
          }
          for(unsigned i = 0; i < nDofsP; i++) {
            unsigned solPDof = dofP_recv[offsetP[k] + i];
            solP1[i] = solP1_recv[offsetP[k] + i];
            solP2[i] = solP2_recv[offsetP[k] + i];;
            sysDof[dim * nDofsV + i ] = pdeSys->GetSystemDof(solP1Index, solP1PdeIndex, solPDof);
            sysDof[dim * nDofsV + nDofsP + i ] = pdeSys->GetSystemDof(solP2Index, solP2PdeIndex, solPDof);
          }
        }

        // local storage of coordinates
        for(unsigned i = 0; i < nDofsX; i++) {
          unsigned coordXDof  = msh->GetSolutionDof(i, iel, solXType);
          for(unsigned d = 0; d < dim; d++) {
            coordX[d][i] = (*msh->_topology->_Sol[d])(coordXDof);
            coordXgf[d][i] = (l == levelC) ? coordX[d][i] : X_recv[d][offsetX[k] + i];
          }
        }

        std::vector<double> a(dim);
        double d;

        std::vector<std::vector<double>> Jacob, JacI;

        const elem_type *femV = msh->_finiteElement[ielGeom][solVType];
        const elem_type *femP = msh->_finiteElement[ielGeom][solPType];
        const elem_type *femXi = msh->_finiteElement[ielGeom][solXiType];
        const elem_type *femPsi = nullptr;

        unsigned nDofsPsi;
        unsigned nDofsN;
        unsigned nDofsK;

        unsigned cut = (fabs(C - 0.5) < 0.1 ) ? 1 : 0;

        if(cut == 1) { // this can happen only on l == levelF
          if (l != levelF) throw std::runtime_error("CUT CELL SHOULD BE ONLY AT LEVELF\n");
          femV = fem.GetFiniteElement(ielGeom, solVType);
          femP = fem.GetFiniteElement(ielGeom, solPType);
          femXi = fem.GetFiniteElement(ielGeom, solXiType);
          femPsi = fem.GetFiniteElement(ielGeom, psiType);

          nDofsPsi = msh->GetElementDofNumber(iel, psiType);
          psi.resize(nDofsPsi);
          for(unsigned i = 0; i < nDofsPsi; i++) {
            unsigned psiDof = msh->GetSolutionDof(i, iel, psiType);
            psi[i] = (*solPsi)(psiDof);
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

          if(cutFemIntegration == CutFemIntegration::Linear) {
            BestFitLinearInterpolation(xg, psig, a);
            d = a[dim];
            a.resize(dim);
          }
        }

        
        const unsigned nGauss = femV->GetGaussPointNumber();

        std::vector<TypeIO> weightCF(nGauss, 0.);
        std::vector<TypeIO> weightCFInt(nGauss, 0.);
        std::vector<TypeIO> weightCFExt(nGauss, 0.);

        if(cut == 1) {

          if(cutFemIntegration == CutFemIntegration::Linear) {

            if(cfw[ielGeom]->GetGaussQuadraturePointNumber() != nGauss) {
              throw std::runtime_error("Linear CutFEM Gauss rule mismatch");
            }

            cfw[ielGeom]->GetWeightWithMap(0, a, d, weightCFInt);
            cfw[ielGeom]->GetWeightWithMap(-1, a, d, weightCF);

          }
          else {

            auto* rule = cfwLS[ielGeom];

            if(rule == nullptr) {
              throw std::runtime_error("Quadratic CutFEM does not support WEDGE");
            }

            if(rule->GetGaussQuadraturePointNumber() != nGauss) {
              throw std::runtime_error("Quadratic CutFEM Gauss rule mismatch");
            }

            rule->SetIntegrationDirection(0);
            rule->SetDirectionValidation(0);

            rule->GetPhysicalSurfaceWeightsFromNodalGeometry(psi, coordX, weightCFInt, weightCF, 1.e-8, 1.e-8);
          }

          for(unsigned i = 0; i < nGauss; ++i) {
            weightCFExt[i] = 1. - weightCFInt[i];
          }
        }
        else {

          for(unsigned i = 0; i < nGauss; ++i) {
            weightCFInt[i] = C;
            weightCFExt[i] = 1. - C;
          }
        }


        // *** Gauss point loop ***
        std::vector<double> Nf(dim);
        std::vector<double> Nls(dim);
        std::vector<double> gradPsi(dim);
        std::vector<double> normalDotGradPhi(nDofsV);
        for(unsigned ig = 0; ig < femV->GetGaussPointNumber(); ig++) {
          // *** get gauss point weight, test function and test function partial derivatives ***

          if(l > levelC) {
            femXi->Jacobian(coordX, ig, weight, phiXi, phiXi_x);
            std::vector <double> xi_gf(dim, 0.);
            for(unsigned i = 0; i < nDofsXi; i++) {
              for(unsigned d = 0; d < dim; d++) {
                xi_gf[d] += coordXi[d][i] * phiXi[i];
              }
            }
            femV->Jacobian(coordXgf, xi_gf, weightFake, phiV, phiV_x);
            femP->GetPhi(phiP, xi_gf);
          }
          else {
            femV->Jacobian(coordX, ig, weight, phiV, phiV_x);
            double *phiPtemp = femP->GetPhi(ig);
            phiP.resize(nDofsP);
            for(unsigned i = 0; i < nDofsP; i++) phiP[i] = phiPtemp[i];
          }

          double x_gss = 0.;

          if (nDofsX != nDofsV) throw std::runtime_error("ADD FEMX");

          for(unsigned i = 0; i < nDofsX; ++i) {
            x_gss += coordXgf[0][i] * phiV[i];
          }

          const double r = axisymmetric ? x_gss - axisCoordinate : 1.;
          const double axisWeight = axisymmetric ? r : 1.;

          if(axisymmetric && dim != 2) {
            throw std::runtime_error("Axisymmetric assembly requires dim == 2");
          }

          if(axisymmetric && r < -1.e-12) {
            throw std::runtime_error("Negative axisymmetric radius");
          }

          double dsN = 0.;
          if (cut == 1) {
            femV->GetJacobianMatrix(coordX, ig, weightFake, Jacob, JacI);
            dsN = 0.0;
            for (unsigned d = 0; d < dim; ++d) {
              double nd = 0.0;
              for (unsigned j = 0; j < dim; ++j)
                nd += JacI[j][d] * a[j];
              Nf[d] = nd;
              dsN += nd * nd;
            }
            dsN = std::sqrt(dsN) + 1.e-20;
            for (unsigned d = 0; d < dim; ++d) {
              Nf[d] /= dsN;
            }
            // for (unsigned i = 0; i < dim; ++i) {
            //   const double ni = Nf[i];
            //   auto& Pi = P[i];
            //   for (unsigned j = 0; j < dim; ++j)
            //     Pi[j] = (i == j ? 1.0 : 0.0) - ni * Nf[j];
            // }
            // for (unsigned j = 0; j < nDofsV; j++) {
            //   const unsigned jOffset = j * dim;
            //   normalDotGradPhi[j] = 0.;
            //   for (unsigned J = 0; J < dim; J++) {
            //     const double gradPhiJ = phiV_x[jOffset + J];
            //     normalDotGradPhi[j] += Nf[J] * gradPhiJ;
            //   }
            // }
            std::vector<double> Nls(dim);
            std::vector<double> gradPsi(dim);
            std::vector<double> xiPsi(dim);

            for(unsigned d = 0; d < dim; ++d) {
              xiPsi[d] = femV->GetGaussRule().GetGaussCoordinatePointer(d)[ig];
            }

            double weightPsi = 0.;

            femPsi->Jacobian(coordX, xiPsi, weightPsi, phiPsi, phiPsi_x);

            std::fill(gradPsi.begin(), gradPsi.end(), 0.);

            for(unsigned i = 0; i < nDofsPsi; ++i) {
              for(unsigned d = 0; d < dim; ++d) {
                gradPsi[d] += psi[i] * phiPsi_x[i * dim + d];
              }
            }

            double normGradPsi = 0.;

            for(unsigned d = 0; d < dim; ++d) {
              normGradPsi += gradPsi[d] * gradPsi[d];
            }

            normGradPsi = std::sqrt(normGradPsi);

            if(normGradPsi > 1.e-12) {
              for(unsigned d = 0; d < dim; ++d) {
                Nls[d] = gradPsi[d] / normGradPsi;
              }
            }
            else {
              for(unsigned d = 0; d < dim; ++d) {
                Nls[d] = Nf[d];
              }
            }

            double normalAlignment = 0.;

            for(unsigned d = 0; d < dim; ++d) {
              normalAlignment += Nls[d] * Nf[d];
            }

            if(normalAlignment < 0.) {
              for(unsigned d = 0; d < dim; ++d) {
                Nls[d] *= -1.;
              }
            }

            for(unsigned i = 0; i < dim; ++i) {
              for(unsigned j = 0; j < dim; ++j) {
                P[i][j] = (i == j ? 1. : 0.) - Nls[i] * Nls[j];
              }
            }

            for(unsigned j = 0; j < nDofsV; ++j) {
              normalDotGradPhi[j] = 0.;
              for(unsigned J = 0; J < dim; ++J) {
                normalDotGradPhi[j] += Nls[J] * phiV_x[j * dim + J];
              }
            }
          }

          double surfaceWeight = 0.;

          if(cut == 1) {

            if(cutFemIntegration == CutFemIntegration::Quadratic) {
              surfaceWeight = cfwLS[ielGeom]->GetGaussWeightPointer()[ig] * weightCF[ig];
            }
            else {
              surfaceWeight = weight * weightCF[ig] * dsN;
            }
          }

          A += weight * weightCFInt[ig];
          Per += surfaceWeight;

          for (unsigned d = 0; d < dim; ++d) {
            solV_gss[d] = 0.0;
            solVOld_gss[d] = 0.0;
            for (unsigned j = 0; j < dim; ++j) {
              gradSolV_gss[d][j] = 0.0;
              gradSolVOld_gss[d][j] = 0.0;
            }
          }

          for (unsigned i = 0; i < nDofsV; ++i) {
            const double phi = phiV[i];
            const unsigned offset = i * dim;
            for (unsigned d = 0; d < dim; ++d) {
              const double u = solV[d][i];
              const double uOld = solVOld[d][i];
              solV_gss[d] += u * phi;
              solVOld_gss[d] += uOld * phi;
              for (unsigned j = 0; j < dim; ++j) {
                const double dphi = phiV_x[offset + j];
                gradSolV_gss[d][j] += u * dphi;
                gradSolVOld_gss[d][j] += uOld * dphi;
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

          const double massCoefficient = rho * weight * axisWeight / dt;

          const double diffusionCoefficientNew = cnew * mu * weight * axisWeight;
          const double diffusionCoefficientOld = cold * mu * weight * axisWeight;

          const double convectionCoefficientNew = cnew * rho * weight * axisWeight;
          const double convectionCoefficientOld = cold * rho * weight * axisWeight;

          // *** phiV_i loop ***
          for(unsigned i = 0; i < nDofsV; i++) {
            for(unsigned  I = 0; I < dim; I++) {  //momentum equation in k
              double NSV = 0.;
              NSV += massCoefficient * phiV[i] * (solV_gss[I] - solVOld_gss[I]);
              for(unsigned J = 0; J < dim; J++) {  // second index j in each equation
                // residual terms
                NSV   +=  diffusionCoefficientNew * phiV_x[i * dim + J] * (gradSolV_gss[I][J] + gradSolV_gss[J][I]); // diffusion
                NSV   +=  diffusionCoefficientOld * phiV_x[i * dim + J] * (gradSolVOld_gss[I][J] + gradSolVOld_gss[J][I]); // diffusion
                if(linear) { // Picard/Hasan linear convective term
                  NSV   +=  phiV[i] * (convectionCoefficientNew * solVOld_gss[J] * gradSolV_gss[I][J] +
                                       convectionCoefficientOld * solV_gss[J] * gradSolVOld_gss[I][J]);
                }
                else { // BackEuler/CrankNicholson nonlinear convective term
                  NSV   +=  convectionCoefficientNew * phiV[i] * (solV_gss[J] * gradSolV_gss[I][J]);
                  NSV   +=  convectionCoefficientOld * phiV[i] * (solVOld_gss[J] * gradSolVOld_gss[I][J]);
                }

              }

              if(axisymmetric && I == 0) {
                NSV += 2. * cnew * mu * solV_gss[0] * phiV[i] * weight / r;
                NSV += 2. * cold * mu * solVOld_gss[0] * phiV[i] * weight / r;
              }

              double weightedDivPhi = axisWeight * phiV_x[i * dim + I];

              if(axisymmetric && I == 0) {
                weightedDivPhi += phiV[i];
              }

              NSV += - weightedDivPhi * (solP1_gss * weightCFInt[ig] + solP2_gss * weightCFExt[ig]) * weight;  // pressure gradient
              NSV += - rho * phiV[i] * g[I] * weight * axisWeight; // gravity term

              Res[I * nDofsV + i] -=  NSV;

              // surface tension stabilization -- rhs contribution
              if(cut == 1) {
                for (int d = 0; d < dim; d++) {
                  Res[I * nDofsV + i] += - sigma * axisWeight * P[I][d] * phiV_x[i * dim + d] * surfaceWeight;
                }

                if(axisymmetric && I == 0) {
                  Res[I * nDofsV + i] += -sigma * phiV[i] * surfaceWeight;
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

                if(axisymmetric && I == 0) {
                  stabSF += (cnew * solV_gss[0] - cold * solVOld_gss[0]) * phiV[i] / (r * r);
                }

                Res[I * nDofsV + i] += -cnew * sigma * dt * stabSF * surfaceWeight * axisWeight;
              }
            }
          } // end phiV_i loop

          // *** phiP_i loop ***
          for(unsigned i = 0; i < nDofsP; i++) {

            double weightedDivU = 0.;

            for(unsigned I = 0; I < dim; ++I) {
              weightedDivU += axisWeight * gradSolV_gss[I][I];
            }

            if(axisymmetric) {
              weightedDivU += solV_gss[0];
            }

            Res[dim * nDofsV + i] += - weightedDivU * phiP[i]  * weight * weightCFInt[ig]; //continuity
            Res[dim * nDofsV + nDofsP + i] += - weightedDivU * phiP[i]  * weight * weightCFExt[ig]; //continuity

            // if(C == 0)
            //   Res[dim * nDofsV + i] += - solP1_gss * phiP[i]  * weight * (1 - C) * eps; //penalty
            // if(C == 1)
            //   Res[dim * nDofsV + nDofsP + i] += - solP2_gss * phiP[i]  * weight * C * eps; //penalty

            // if(C > 0.1)
            //   Res[dim * nDofsV + i] += - solP1_gss * phiP[i]  * weight * epsilonp; //penalty
            // if(C < 0.9)
            //   Res[dim * nDofsV + nDofsP + i] += - solP2_gss * phiP[i] * weight * epsilonp; //penalty

          } // end phiP_i loop
          // end gauss point loop

          //--------------------------------------------------------------------------------------------------------
          // Add the local Matrix/Vector into the global Matrix/Vector
          std::vector<double> solVdotGradphi_gss(nDofsV, 0.);
          if(linear) {
            for(unsigned j = 0; j < nDofsV; j++) {
              for(unsigned J = 0; J < dim ; J++) {
                solVdotGradphi_gss[j] += solVOld_gss[J] * phiV_x[j * dim + J];
              }
            }
          }
          else {
            for(unsigned j = 0; j < nDofsV; j++) {
              for(unsigned J = 0; J < dim ; J++) {
                solVdotGradphi_gss[j] += solV_gss[J] * phiV_x[j * dim + J];
              }
            }
          }

          double surfaceTensionCoefficient = cnew * cnew * sigma * dt * surfaceWeight * axisWeight * 0;

          for(unsigned i = 0; i < nDofsV; i++) {
            for(unsigned j = 0; j < nDofsV; j++) {

              const double phiProduct = phiV[i] * phiV[j];

              double gradPhi_dot_gradPhi = 0.0;
              for(unsigned J = 0; J < dim ; J++) {
                gradPhi_dot_gradPhi += phiV_x[i * dim + J] * phiV_x[j * dim + J]; //diagonal diffusion
              }

              double Phi_dot_solVGradPhi = solVdotGradphi_gss[j] * phiV[i];

              double diagonalContribution = massCoefficient * phiProduct +
                                            diffusionCoefficientNew * gradPhi_dot_gradPhi +
                                            convectionCoefficientNew * Phi_dot_solVGradPhi ;

              if (cut == 1) {
                const double tangentialGradDot = gradPhi_dot_gradPhi - normalDotGradPhi[i] * normalDotGradPhi[j];
                diagonalContribution += surfaceTensionCoefficient * tangentialGradDot;
              }

              for(unsigned I = 0; I < dim; I++) { //row velocity blocks or dimension
                unsigned VIrow = (I * nDofsV + i) * nDofsVP;
                unsigned VIcolumn = I * nDofsV + j;

                Jac[ VIrow + VIcolumn] += diagonalContribution; //diagonal

                if(axisymmetric && cut == 1 && I == 0) {
                  Jac[VIrow + VIcolumn] += surfaceTensionCoefficient * phiProduct / (r * r);
                }

                if(axisymmetric && I == 0) {
                  Jac[VIrow + VIcolumn] += 2. * cnew * mu * phiV[i] * phiV[j] * weight / r;
                }

                for(unsigned J = 0; J < dim ; J++) { //column velocity blocks or dimension
                  unsigned VJcolumn = J * nDofsV + j;
                  Jac[ VIrow + VJcolumn ] += diffusionCoefficientNew * phiV_x[i * dim + J] * phiV_x[j * dim + I]; //off-diagonal diffusion
                  Jac[ VIrow + VJcolumn ] += (linear) ? convectionCoefficientOld * phiProduct * gradSolVOld_gss[I][J] //off-diagonal convection
                                             : convectionCoefficientNew * phiProduct * gradSolV_gss[I][J]; //off-diagonal convection
                }
              }
            }
          }

          for(unsigned i = 0; i < nDofsV; i++) {
            for(unsigned I = 0; I < dim; I++) { //row velocity blocks or dimension
              unsigned VIrow = I * nDofsV + i;

              double weightedDivPhi = axisWeight * phiV_x[i * dim + I];

              if(axisymmetric && I == 0) {
                weightedDivPhi += phiV[i];
              }

              for(unsigned j = 0; j < nDofsP; j++) {
                unsigned P1column = dim * nDofsV + j;
                unsigned P2column = dim * nDofsV + nDofsP + j;
                Jac[VIrow * nDofsVP + P1column] += - weightedDivPhi * phiP[j] * weight * weightCFInt[ig]; //pressure gradient
                Jac[VIrow * nDofsVP + P2column] += - weightedDivPhi * phiP[j] * weight * weightCFExt[ig]; //pressure gradient

                Jac[P1column * nDofsVP + VIrow] -= - weightedDivPhi * phiP[j] * weight * weightCFInt[ig]; //continuity
                Jac[P2column * nDofsVP + VIrow] -= - weightedDivPhi * phiP[j] * weight * weightCFExt[ig]; //continuity
              }
            }
          }
          // for(unsigned i = 0; i < nDofsP; i++) {
          //   unsigned P1row = dim * nDofsV + i;
          //   unsigned P2row = dim * nDofsV + nDofsP + i;
          //   for(unsigned j = 0; j < nDofsP; j++) {
          //     unsigned P1column = dim * nDofsV + j;
          //     unsigned P2column = dim * nDofsV + nDofsP + j;
          //     if(C == 0)
          //       Jac[P1row * nDofsVP + P1column] += phiP[i] * phiP[j] * weight * (1 - C) * eps; //continuity
          //     if(C == 1)
          //       Jac[P2row * nDofsVP + P2column] += phiP[i] * phiP[j] * weight * C * eps; //continuity

          //     if(C > 0.1)
          //       Jac[P1row * nDofsVP + P1column] += phiP[i] * phiP[j] * weight  * epsilonp; //continuity
          //     if(C < 0.9)
          //       Jac[P2row * nDofsVP + P2column] += phiP[i] * phiP[j] * weight  * epsilonp; //continuity
          //   }
          // }

        }

        RES->add_vector_blocked(Res, sysDof);
        KK->add_matrix_blocked(Jac, sysDof, sysDof);
      }
      // else {

      //   std::vector<std::vector<double>> Xichild(dim * nChildren);

      //   for (unsigned d = 0; d < dim; ++d) {
      //     std::vector<std::vector<double>> XichildComponent;
      //     ProlongateFieldToChildren(XichildComponent, PMatrix[iel_type], coordXi[d]);
      //     for (unsigned j = 0; j < nChildren; ++j) {
      //       Xichild[d * nChildren + j].swap(XichildComponent[j]);
      //     }
      //   }

      //   for (unsigned j = 0; j < nChildren; ++j) {
      //     const unsigned jel = el->GetChild(iel, j);
      //     const unsigned jproc = mshp->IsdomBisectionSearch(jel, 3);
      //     el_send[jproc][posEl_send[jproc]] = jel;
      //     // el_gf_send[jproc][posEl_send[jproc]] = iel_gf;
      //     ++posEl_send[jproc];

      //     for (unsigned i = 0; i < nDofsV; ++i) {
      //       if (l == levelC) {
      //         const unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);
      //         dofV_send[jproc][posV_send[jproc]] = solVDof;
      //         for (unsigned d = 0; d < dim; ++d) {
      //           solV_send[d * nprocs + jproc][posV_send[jproc]] = (*sol->_Sol[solVIndex[d]])(solVDof);
      //           solVOld_send[d * nprocs + jproc][posV_send[jproc]] = (*sol->_SolOld[solVIndex[d]])(solVDof);
      //         }
      //       }
      //       else {
      //         dofV_send[jproc][posV_send[jproc]] = dofV_recv[offsetV[k] + i];
      //         for (unsigned d = 0; d < dim; ++d) {
      //           solV_send[d * nprocs + jproc][posV_send[jproc]] = solV_recv[d][offsetV[k] + i];
      //           solVOld_send[d * nprocs + jproc][posV_send[jproc]] = solVOld_recv[d][offsetV[k] + i];
      //         }
      //       }
      //       ++posV_send[jproc];
      //     }

      //     for (unsigned i = 0; i < nDofsP; ++i) {
      //       if (l == levelC) {
      //         const unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);
      //         dofP_send[jproc][posP_send[jproc]] = solPDof;
      //         solP1_send[jproc][posP_send[jproc]] = (*sol->_Sol[solP1Index])(solPDof);
      //         solP2_send[jproc][posP_send[jproc]] = (*sol->_Sol[solP2Index])(solPDof);
      //       }
      //       else {
      //         dofP_send[jproc][posP_send[jproc]] = dofP_recv[offsetP[k] + i];
      //         solP1_send[jproc][posP_send[jproc]] = solP1_recv[offsetP[k] + i];
      //         solP2_send[jproc][posP_send[jproc]] = solP2_recv[offsetP[k] + i];
      //       }
      //       ++posP_send[jproc];
      //     }

      //     for (unsigned i = 0; i < nDofsX; ++i) {
      //       if (l == levelC) {
      //         const unsigned XDof = msh->GetSolutionDof(i, iel, solXType);
      //         for (unsigned d = 0; d < dim; ++d) {
      //           X_send[d * nprocs + jproc][posX_send[jproc]] = (*msh->_topology->_Sol[d])(XDof);
      //         }
      //       }
      //       else {
      //         for (unsigned d = 0; d < dim; ++d) {
      //           X_send[d * nprocs + jproc][posX_send[jproc]] = X_recv[d][offsetX[k] + i];
      //         }
      //       }
      //       ++posX_send[jproc];
      //     }

      //     for (unsigned i = 0; i < nDofsXi; ++i) {
      //       for (unsigned d = 0; d < dim; ++d) {
      //         Xi_send[d * nprocs + jproc][posXi_send[jproc]] = Xichild[d * nChildren + j][i];
      //       }
      //       ++posXi_send[jproc];
      //     }
      //   }
      // }
      else {

        const bool refined = (msh->GetRefinedElementIndex(iel) != 0u);
        const unsigned nChildrenCurrent = GetNumberOfChildren(iel);

        std::vector<std::vector<double>> Xichild;

        if(refined) {

          Xichild.resize(dim * nChildrenCurrent);

          for(unsigned d = 0; d < dim; ++d) {
            std::vector<std::vector<double>> XichildComponent;
            ProlongateFieldToChildren(XichildComponent, PMatrix[iel_type], coordXi[d]);

            for(unsigned j = 0; j < nChildrenCurrent; ++j) {
              Xichild[d * nChildrenCurrent + j].swap(XichildComponent[j]);
            }
          }

        }

        for(unsigned j = 0; j < nChildrenCurrent; ++j) {

          const unsigned jel = el->GetChild(iel, j);
          const unsigned jproc = mshp->IsdomBisectionSearch(jel, 3);

          el_send[jproc][posEl_send[jproc]] = jel;
          ++posEl_send[jproc];

          for(unsigned i = 0; i < nDofsV; ++i) {

            if(l == levelC) {

              const unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);
              dofV_send[jproc][posV_send[jproc]] = solVDof;

              for(unsigned d = 0; d < dim; ++d) {
                solV_send[d * nprocs + jproc][posV_send[jproc]] = (*sol->_Sol[solVIndex[d]])(solVDof);
                solVOld_send[d * nprocs + jproc][posV_send[jproc]] = (*sol->_SolOld[solVIndex[d]])(solVDof);
              }
            }
            else {

              dofV_send[jproc][posV_send[jproc]] = dofV_recv[offsetV[k] + i];

              for(unsigned d = 0; d < dim; ++d) {
                solV_send[d * nprocs + jproc][posV_send[jproc]] = solV_recv[d][offsetV[k] + i];
                solVOld_send[d * nprocs + jproc][posV_send[jproc]] = solVOld_recv[d][offsetV[k] + i];
              }
            }

            ++posV_send[jproc];
          }

          for(unsigned i = 0; i < nDofsP; ++i) {

            if(l == levelC) {
              const unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);

              dofP_send[jproc][posP_send[jproc]] = solPDof;
              solP1_send[jproc][posP_send[jproc]] = (*sol->_Sol[solP1Index])(solPDof);
              solP2_send[jproc][posP_send[jproc]] = (*sol->_Sol[solP2Index])(solPDof);
            }
            else {
              dofP_send[jproc][posP_send[jproc]] = dofP_recv[offsetP[k] + i];
              solP1_send[jproc][posP_send[jproc]] = solP1_recv[offsetP[k] + i];
              solP2_send[jproc][posP_send[jproc]] = solP2_recv[offsetP[k] + i];
            }

            ++posP_send[jproc];
          }

          for(unsigned i = 0; i < nDofsX; ++i) {

            if(l == levelC) {

              const unsigned XDof = msh->GetSolutionDof(i, iel, solXType);

              for(unsigned d = 0; d < dim; ++d) {
                X_send[d * nprocs + jproc][posX_send[jproc]] = (*msh->_topology->_Sol[d])(XDof);
              }
            }
            else {

              for(unsigned d = 0; d < dim; ++d) {
                X_send[d * nprocs + jproc][posX_send[jproc]] = X_recv[d][offsetX[k] + i];
              }
            }

            ++posX_send[jproc];
          }

          for(unsigned i = 0; i < nDofsXi; ++i) {

            for(unsigned d = 0; d < dim; ++d) {
              Xi_send[d * nprocs + jproc][posXi_send[jproc]] = refined ? Xichild[d * nChildrenCurrent + j][i] : coordXi[d][i];
            }

            ++posXi_send[jproc];
          }
        }
      }
    }

    if(l < levelF) {
      std::vector<unsigned> n_el_recv(nprocs, 0);
      std::vector<std::vector<unsigned>> n_dof_recv(nprocs, std::vector<unsigned>(nSolutionTypes, 0));
      n_el_recv[iproc] = n_el_send[iproc];
      n_dof_recv[iproc] = n_dof_send[iproc];

      for (unsigned p = 0; p < nprocs; ++p) {
        if (p != iproc) {
          postRecv(&n_el_recv[p], 1, MPI_UNSIGNED, p, tagElCounts, requests);
          postRecv(n_dof_recv[p].data(), nSolutionTypes, MPI_UNSIGNED, p, tagCounts, requests);
        }
      }

      for (unsigned p = 0; p < nprocs; ++p) {
        if (p != iproc) {
          postSend(&n_el_send[p], 1, MPI_UNSIGNED, p, tagElCounts, requests);
          postSend(n_dof_send[p].data(), nSolutionTypes, MPI_UNSIGNED, p, tagCounts, requests);
        }
      }

      waitAll(requests);

      std::vector<unsigned> offsetEl_recv(nprocs + 1, 0);
      std::vector<unsigned> offsetV_recv(nprocs + 1, 0);
      std::vector<unsigned> offsetP_recv(nprocs + 1, 0);
      std::vector<unsigned> offsetX_recv(nprocs + 1, 0);
      std::vector<unsigned> offsetXi_recv(nprocs + 1, 0);

      for (unsigned p = 0; p < nprocs; ++p) {
        offsetEl_recv[p + 1] = offsetEl_recv[p] + n_el_recv[p];
        offsetV_recv[p + 1] = offsetV_recv[p] + n_dof_recv[p][solVType];
        offsetP_recv[p + 1] = offsetP_recv[p] + n_dof_recv[p][solPType];
        offsetX_recv[p + 1] = offsetX_recv[p] + n_dof_recv[p][solXType];
        offsetXi_recv[p + 1] = offsetXi_recv[p] + n_dof_recv[p][solXiType];
      }

      el_recv.resize(offsetEl_recv[nprocs]);
      // el_gf_recv.resize(offsetEl_recv[nprocs]);
      dofV_recv.resize(offsetV_recv[nprocs]);
      dofP_recv.resize(offsetP_recv[nprocs]);
      solP1_recv.resize(offsetP_recv[nprocs]);
      solP2_recv.resize(offsetP_recv[nprocs]);
      X_recv.resize(dim);
      Xi_recv.resize(dim);
      solV_recv.resize(dim);
      solVOld_recv.resize(dim);

      if (n_el_send[iproc] != 0) {
        std::copy(el_send[iproc].begin(), el_send[iproc].end(), el_recv.data() + offsetEl_recv[iproc]);
        // std::copy(el_gf_send[iproc].begin(), el_gf_send[iproc].end(), el_gf_recv.data() + offsetEl_recv[iproc]);
      }

      if (n_dof_send[iproc][solVType] != 0) {
        std::copy(dofV_send[iproc].begin(), dofV_send[iproc].end(), dofV_recv.data() + offsetV_recv[iproc]);
      }

      if (n_dof_send[iproc][solPType] != 0) {
        std::copy(dofP_send[iproc].begin(), dofP_send[iproc].end(), dofP_recv.data() + offsetP_recv[iproc]);
        std::copy(solP1_send[iproc].begin(), solP1_send[iproc].end(), solP1_recv.data() + offsetP_recv[iproc]);
        std::copy(solP2_send[iproc].begin(), solP2_send[iproc].end(), solP2_recv.data() + offsetP_recv[iproc]);
      }

      for (unsigned d = 0; d < dim; ++d) {
        X_recv[d].resize(offsetX_recv[nprocs]);
        Xi_recv[d].resize(offsetXi_recv[nprocs]);
        solV_recv[d].resize(offsetV_recv[nprocs]);
        solVOld_recv[d].resize(offsetV_recv[nprocs]);

        if (n_dof_send[iproc][solXType] != 0) {
          std::copy(X_send[d * nprocs + iproc].begin(), X_send[d * nprocs + iproc].end(), X_recv[d].data() + offsetX_recv[iproc]);
        }
        if (n_dof_send[iproc][solXiType] != 0) {
          std::copy(Xi_send[d * nprocs + iproc].begin(), Xi_send[d * nprocs + iproc].end(), Xi_recv[d].data() + offsetXi_recv[iproc]);
        }
        if (n_dof_send[iproc][solVType] != 0) {
          std::copy(solV_send[d * nprocs + iproc].begin(), solV_send[d * nprocs + iproc].end(), solV_recv[d].data() + offsetV_recv[iproc]);
          std::copy(solVOld_send[d * nprocs + iproc].begin(), solVOld_send[d * nprocs + iproc].end(), solVOld_recv[d].data() + offsetV_recv[iproc]);
        }
      }

      for (unsigned p = 0; p < nprocs; ++p) {
        if (p == iproc) {
          continue;
        }

        if (n_el_recv[p] != 0) {
          postRecv(el_recv.data() + offsetEl_recv[p], n_el_recv[p], MPI_UNSIGNED, p, tagEl, requests);
          // postRecv(el_gf_recv.data() + offsetEl_recv[p], n_el_recv[p], MPI_UNSIGNED, p, tagElGF, requests);
        }
        if (n_dof_recv[p][solVType] != 0) {
          postRecv(dofV_recv.data() + offsetV_recv[p], n_dof_recv[p][solVType], MPI_UNSIGNED, p, tagV, requests);
        }
        if (n_dof_recv[p][solPType] != 0) {
          postRecv(dofP_recv.data() + offsetP_recv[p], n_dof_recv[p][solPType], MPI_UNSIGNED, p, tagP, requests);
          postRecv(solP1_recv.data() + offsetP_recv[p], n_dof_recv[p][solPType], MPI_DOUBLE, p, tagSolP1, requests);
          postRecv(solP2_recv.data() + offsetP_recv[p], n_dof_recv[p][solPType], MPI_DOUBLE, p, tagSolP2, requests);
        }

        for (unsigned d = 0; d < dim; ++d) {
          if (n_dof_recv[p][solXType] != 0) {
            postRecv(X_recv[d].data() + offsetX_recv[p], n_dof_recv[p][solXType], MPI_DOUBLE, p, tagX + static_cast<int>(d), requests);
          }
          if (n_dof_recv[p][solXiType] != 0) {
            postRecv(Xi_recv[d].data() + offsetXi_recv[p], n_dof_recv[p][solXiType], MPI_DOUBLE, p, tagXi + static_cast<int>(d), requests);
          }
          if (n_dof_recv[p][solVType] != 0) {
            postRecv(solV_recv[d].data() + offsetV_recv[p], n_dof_recv[p][solVType], MPI_DOUBLE, p, tagSolV + static_cast<int>(d), requests);
            postRecv(solVOld_recv[d].data() + offsetV_recv[p], n_dof_recv[p][solVType], MPI_DOUBLE, p, tagSolVOld + static_cast<int>(d), requests);
          }
        }
      }

      for (unsigned p = 0; p < nprocs; ++p) {
        if (p == iproc) {
          continue;
        }

        postSend(el_send[p].data(), n_el_send[p], MPI_UNSIGNED, p, tagEl, requests);
        // postSend(el_gf_send[p].data(), n_el_send[p], MPI_UNSIGNED, p, tagElGF, requests);
        postSend(dofV_send[p].data(), n_dof_send[p][solVType], MPI_UNSIGNED, p, tagV, requests);
        postSend(dofP_send[p].data(), n_dof_send[p][solPType], MPI_UNSIGNED, p, tagP, requests);
        postSend(solP1_send[p].data(), n_dof_send[p][solPType], MPI_DOUBLE, p, tagSolP1, requests);
        postSend(solP2_send[p].data(), n_dof_send[p][solPType], MPI_DOUBLE, p, tagSolP2, requests);

        for (unsigned d = 0; d < dim; ++d) {
          postSend(X_send[d * nprocs + p].data(), n_dof_send[p][solXType], MPI_DOUBLE, p, tagX + static_cast<int>(d), requests);
          postSend(Xi_send[d * nprocs + p].data(), n_dof_send[p][solXiType], MPI_DOUBLE, p, tagXi + static_cast<int>(d), requests);
          postSend(solV_send[d * nprocs + p].data(), n_dof_send[p][solVType], MPI_DOUBLE, p, tagSolV + static_cast<int>(d), requests);
          postSend(solVOld_send[d * nprocs + p].data(), n_dof_send[p][solVType], MPI_DOUBLE, p, tagSolVOld + static_cast<int>(d), requests);
        }
      }

      waitAll(requests);
    }

    RES->close();
    KK->flush();
  }

  //END parallel exchange between levels

  RES->close();
  KK->close();

  std::cout << "Matrix Assembly time        = " << static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC << std::flush << std::endl;

  MPI_Allreduce(MPI_IN_PLACE, &A, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &Per, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


  double radius  = 0.25;
  double radius2 = radius * radius;
  double radius3 = radius2 * radius;
  double per_th = (dim == 2) ? 2 * M_PI * radius : 4. * M_PI * radius2;
  double A_th = (dim == 2) ? M_PI * radius2 : 4. / 3. * M_PI * radius3;
  std::cout << " Area error " << std::abs(A - A_th)/A_th << std::endl;
  std::cout << " Perimetro error " << std::abs(Per - per_th)/per_th << std::endl;

}

void AssemblePsi(MultiLevelProblem& ml_prob3) {
  LinearImplicitSystem* mlPdeSys = &ml_prob3.get_system<LinearImplicitSystem>("PSI");
  const unsigned levelToAssemble = mlPdeSys->GetLevelToAssemble();
  MultiLevelSolution* mlSol3 = ml_prob3._ml_sol;
  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[levelToAssemble];
  SparseMatrix* KK = pdeSys->_KK;
  NumericVector* RES = pdeSys->_RES;

  const MultiphaseParams mParam = ml_prob3.GetMultiphaseParams();
  MultiLevelSolution* mlSol0 = static_cast<MultiLevelSolution*>(mParam.extProblem);

  const unsigned levelF = mParam.levelF;
  const unsigned levelPsi = mParam.levelC;
  const unsigned level0 = mParam.level0;

  const unsigned levelStart = level0 + levelToAssemble;
  if (levelStart > levelPsi || levelPsi > levelF) {
    throw std::runtime_error("AssemblePsi: invalid source/target level range");
  }

  const unsigned solPsiIndex = mlSol3->GetIndex("Psi");
  const unsigned solPsiPdeIndex = mlPdeSys->GetSolPdeIndex("Psi");
  const unsigned psiType = mlSol3->GetSolutionType(solPsiIndex);
  const unsigned sourcePsiIndex = mlSol0->GetIndex("Psi");
  const unsigned sourcePsiType = mlSol0->GetSolutionType(sourcePsiIndex);

  const unsigned solXiType = 0;
  const unsigned solXType = 2;
  if (psiType != sourcePsiType) {
    throw std::runtime_error("AssemblePsi: source and target Psi FE types differ");
  }

  const bool axisymmetric = mParam.axisymmetric;
  const double axisCoordinate = mParam.axisCoordinate;

  Mesh* targetMesh = ml_prob3._ml_msh->GetLevel(levelToAssemble);
  const unsigned dim = targetMesh->GetDimension();
  if (axisymmetric && dim != 2) {
    throw std::runtime_error("AssemblePsi: axisymmetry requires a 2D mesh");
  }

  int nprocs = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  const unsigned nProc = static_cast<unsigned>(nprocs);

  KK->zero();
  RES->zero();

  const int tagCounts = 210;
  const int tagElements = 211;
  const int tagPsiDofs = 212;
  const int tagPsiValues = 213;
  const int tagXi = 214;

  std::vector<MPI_Request> requests;
  std::vector<unsigned> el_recv;
  std::vector<unsigned> dofPsi_recv;
  std::vector<double> psiTarget_recv;
  std::vector<std::vector<double>> Xi_recv(dim);

  Solution* sol3 = mlSol3->GetLevel(levelToAssemble);

  for (unsigned l = levelStart; l <= levelF; ++l) {
    Solution* sol0 = mlSol0->GetLevel(l);
    Mesh* msh = sol0->GetMesh();
    elem* el = msh->GetElementArray();
    const unsigned iproc = msh->processor_id();

    const unsigned nElementsToLoop = (l == levelStart) ? msh->_elementOffset[iproc + 1] - msh->_elementOffset[iproc] : static_cast<unsigned>(el_recv.size());

    std::vector<unsigned> offsetPsi(nElementsToLoop + 1, 0);
    std::vector<unsigned> offsetXi(nElementsToLoop + 1, 0);

    for (unsigned k = 0; k < nElementsToLoop; ++k) {
      const unsigned iel = (l == levelStart) ? msh->_elementOffset[iproc] + k : el_recv[k];
      offsetPsi[k + 1] = offsetPsi[k] + msh->GetElementDofNumber(iel, psiType);
      offsetXi[k + 1] = offsetXi[k] + msh->GetElementDofNumber(iel, solXiType);
    }

    Mesh* mshp = (l < levelF) ? mlSol0->GetLevel(l + 1)->GetMesh() : nullptr;
    const unsigned nChildrenRefined = msh->GetRefIndex();

    std::vector<std::vector<unsigned>> sendCounts(nProc, std::vector<unsigned>(3, 0));

    if (l < levelF) {
      for (unsigned k = 0; k < nElementsToLoop; ++k) {
        const unsigned iel = (l == levelStart) ? msh->_elementOffset[iproc] + k : el_recv[k];

        const bool refined = (msh->GetRefinedElementIndex(iel) != 0u);
        const unsigned nChildren = refined ? nChildrenRefined : 1u;
        const unsigned nDofsPsi = msh->GetElementDofNumber(iel, psiType);

        const unsigned nDofsXi = msh->GetElementDofNumber(iel, solXiType);

        for (unsigned j = 0; j < nChildren; ++j) {
          const unsigned jel = el->GetChild(iel, j);
          const unsigned jproc = mshp->IsdomBisectionSearch(jel, 3);
          ++sendCounts[jproc][0];
          sendCounts[jproc][1] += nDofsPsi;
          sendCounts[jproc][2] += nDofsXi;
        }
      }
    }

    std::vector<std::vector<unsigned>> el_send(nProc);
    std::vector<std::vector<unsigned>> dofPsi_send(nProc);
    std::vector<std::vector<double>> psiTarget_send(nProc);
    std::vector<std::vector<double>> Xi_send(dim * nProc);

    for (unsigned p = 0; p < nProc; ++p) {
      el_send[p].resize(sendCounts[p][0]);
      dofPsi_send[p].resize(sendCounts[p][1]);
      psiTarget_send[p].resize(sendCounts[p][1]);
      for (unsigned d = 0; d < dim; ++d) {
        Xi_send[d * nProc + p].resize(sendCounts[p][2]);
      }
    }

    std::vector<unsigned> posEl(nProc, 0);
    std::vector<unsigned> posPsi(nProc, 0);
    std::vector<unsigned> posXi(nProc, 0);

    for (unsigned k = 0; k < nElementsToLoop; ++k) {
      const unsigned iel = (l == levelStart) ? msh->_elementOffset[iproc] + k : el_recv[k];
      const unsigned ielType = msh->GetElementType(iel);
      const unsigned nDofsPsi = msh->GetElementDofNumber(iel, psiType);
      const unsigned nDofsXi = msh->GetElementDofNumber(iel, solXiType);

      std::vector<std::vector<double>> coordXi(dim, std::vector<double>(nDofsXi));
      for (unsigned i = 0; i < nDofsXi; ++i) {
        if (l == levelStart) {
          const double* xi_pt = msh->_finiteElement[ielType][solXiType]->GetBasis()->GetXcoarse(i);
          for (unsigned d = 0; d < dim; ++d) {
            coordXi[d][i] = xi_pt[d];
          }
        }
        else {
          for (unsigned d = 0; d < dim; ++d) {
            coordXi[d][i] = Xi_recv[d][offsetXi[k] + i];
          }
        }
      }

      if (l == levelF) {
        std::vector<unsigned> sysDof(nDofsPsi);
        std::vector<double> targetPsi(nDofsPsi);

        for (unsigned i = 0; i < nDofsPsi; ++i) {
          if (l == levelStart) {
            const unsigned psiDof = msh->GetSolutionDof(i, iel, psiType);
            sysDof[i] = pdeSys->GetSystemDof(solPsiIndex, solPsiPdeIndex, i, iel);
            targetPsi[i] = (*sol3->_Sol[solPsiIndex])(psiDof);
          }
          else {
            const unsigned psiDof = dofPsi_recv[offsetPsi[k] + i];
            sysDof[i] = pdeSys->GetSystemDof(solPsiIndex, solPsiPdeIndex, psiDof);
            targetPsi[i] = psiTarget_recv[offsetPsi[k] + i];
          }
        }

        const unsigned nDofsSource = msh->GetElementDofNumber(iel, sourcePsiType);
        std::vector<double> sourcePsi(nDofsSource);
        for (unsigned i = 0; i < nDofsSource; ++i) {
          const unsigned sourceDof = msh->GetSolutionDof(i, iel, sourcePsiType);
          sourcePsi[i] = (*sol0->_Sol[sourcePsiIndex])(sourceDof);
        }

        const unsigned nDofsX = msh->GetElementDofNumber(iel, solXType);
        std::vector<std::vector<double>> coordX(dim, std::vector<double>(nDofsX));
        for (unsigned i = 0; i < nDofsX; ++i) {
          const unsigned XDof = msh->GetSolutionDof(i, iel, solXType);
          for (unsigned d = 0; d < dim; ++d) {
            coordX[d][i] = (*msh->_topology->_Sol[d])(XDof);
          }
        }

        const elem_type* femSource = msh->_finiteElement[ielType][sourcePsiType];
        const elem_type* femTarget = msh->_finiteElement[ielType][psiType];
        const elem_type* femXi = msh->_finiteElement[ielType][solXiType];

        // Weighted H1 projection: (u,v) + h1GradientWeight*(grad u,grad v).
        // Use 1.0 for the standard H1 inner product; 0.0 reproduces L2.
        constexpr double h1GradientWeight = 1.0;

        std::vector<double> phiSource, phiSource_x, phiTarget, phiTarget_x;
        std::vector<double> Res(nDofsPsi, 0.);
        std::vector<double> Jac(nDofsPsi * nDofsPsi, 0.);

        for (unsigned ig = 0; ig < femSource->GetGaussPointNumber(); ++ig) {
          double weight = 0.;
          femSource->Jacobian(coordX, ig, weight, phiSource, phiSource_x);

          if (l == levelStart) {
            phiTarget = phiSource;
            phiTarget_x = phiSource_x;
          }
          else {
            // eta is the reference coordinate on the CURRENT (fine) element.
            // xi  is the reference coordinate on the TARGET (coarse) element.
            std::vector<double> eta(dim, 0.), xi(dim, 0.);
            const double* phiXi = femXi->GetPhi(ig);
            for (unsigned d = 0; d < dim; ++d) {
              eta[d] = femSource->GetGaussRule().GetGaussCoordinatePointer(d)[ig];
            }
            for (unsigned i = 0; i < nDofsXi; ++i) {
              for (unsigned d = 0; d < dim; ++d) {
                xi[d] += coordXi[d][i] * phiXi[i];
              }
            }
            femTarget->GetPhi(phiTarget, xi);

            // dN_target/dxi: target basis derivatives in coarse reference coordinates.
            std::vector<std::vector<double>> dPhiTarget(dim);
            femTarget->GetDPhiDXi(dPhiTarget[0], xi);
            femTarget->GetDPhiDEta(dPhiTarget[1], xi);
            if (dim == 3) femTarget->GetDPhiDZeta(dPhiTarget[2], xi);

            // dN_Xi/deta: mapping derivatives from current reference element
            // to the coarse reference element (stored in coordXi).
            std::vector<std::vector<double>> dPhiXi(dim);
            femXi->GetDPhiDXi(dPhiXi[0], eta);
            femXi->GetDPhiDEta(dPhiXi[1], eta);
            if (dim == 3) femXi->GetDPhiDZeta(dPhiXi[2], eta);

            // dxi_m/deta_q.
            std::vector<std::vector<double>> mapJac(dim, std::vector<double>(dim, 0.));
            for (unsigned m = 0; m < dim; ++m) {
              for (unsigned q = 0; q < dim; ++q) {
                for (unsigned kXi = 0; kXi < nDofsXi; ++kXi) {
                  mapJac[m][q] += coordXi[m][kXi] * dPhiXi[q][kXi];
                }
              }
            }

            // deta_q/dx_d: inverse Jacobian of the actual (possibly curved)
            // fine-element geometry. Do NOT use coordXi as physical geometry.
            double unusedWeight = 0.;
            std::vector<std::vector<double>> geometricJac, geometricJacInv;
            femSource->GetJacobianMatrix(coordX, ig, unusedWeight,
                                         geometricJac, geometricJacInv);

            // Chain rule: dN_target/dx_d =
            // sum_{m,q} dN_target/dxi_m * dxi_m/deta_q * deta_q/dx_d.
            phiTarget_x.assign(nDofsPsi * dim, 0.);
            for (unsigned i = 0; i < nDofsPsi; ++i) {
              for (unsigned d = 0; d < dim; ++d) {
                for (unsigned m = 0; m < dim; ++m) {
                  for (unsigned q = 0; q < dim; ++q) {
                    phiTarget_x[i * dim + d] +=
                      dPhiTarget[m][i] * mapJac[m][q] * geometricJacInv[q][d];
                  }
                }
              }
            }
          }

          double psiFine = 0.;
          double psiCoarse = 0.;
          std::vector<double> gradPsiFine(dim, 0.), gradPsiCoarse(dim, 0.);

          for (unsigned j = 0; j < nDofsSource; ++j) {
            psiFine += sourcePsi[j] * phiSource[j];
            for (unsigned d = 0; d < dim; ++d) {
              gradPsiFine[d] += sourcePsi[j] * phiSource_x[j * dim + d];
            }
          }
          for (unsigned j = 0; j < nDofsPsi; ++j) {
            psiCoarse += targetPsi[j] * phiTarget[j];
            for (unsigned d = 0; d < dim; ++d) {
              gradPsiCoarse[d] += targetPsi[j] * phiTarget_x[j * dim + d];
            }
          }

          if (axisymmetric) {
            if (nDofsX != nDofsSource) {
              throw std::runtime_error("AssemblePsi: expected Q2 geometry");
            }
            double xGauss = 0.;
            for (unsigned j = 0; j < nDofsX; ++j) {
              xGauss += coordX[0][j] * phiSource[j];
            }
            const double r = xGauss - axisCoordinate;
            if (r < -1.e-12) {
              throw std::runtime_error("AssemblePsi: negative radial coordinate");
            }
            weight *= std::max(r, 0.);
          }

          for (unsigned i = 0; i < nDofsPsi; ++i) {
            double gradientResidual = 0.;
            for (unsigned d = 0; d < dim; ++d) {
              gradientResidual += phiTarget_x[i * dim + d] *
                                  (gradPsiFine[d] - gradPsiCoarse[d]);
            }
            Res[i] += weight * (phiTarget[i] * (psiFine - psiCoarse)
                                + h1GradientWeight * gradientResidual);

            for (unsigned j = 0; j < nDofsPsi; ++j) {
              double gradientMatrix = 0.;
              for (unsigned d = 0; d < dim; ++d) {
                gradientMatrix += phiTarget_x[i * dim + d] * phiTarget_x[j * dim + d];
              }
              Jac[i * nDofsPsi + j] += weight *
                (phiTarget[i] * phiTarget[j] + h1GradientWeight * gradientMatrix);
            }
          }
        }

        RES->add_vector_blocked(Res, sysDof);
        KK->add_matrix_blocked(Jac, sysDof, sysDof);
      }
      else {
        const bool refined = (msh->GetRefinedElementIndex(iel) != 0u);
        const unsigned nChildren = refined ? nChildrenRefined : 1u;
        std::vector<std::vector<double>> Xichild;

        if (refined) {
          Xichild.resize(dim * nChildren);
          for (unsigned d = 0; d < dim; ++d) {
            std::vector<std::vector<double>> XichildComponent;
            ProlongateFieldToChildren(XichildComponent, PMatrix[ielType], coordXi[d]);
            for (unsigned j = 0; j < nChildren; ++j) {
              Xichild[d * nChildren + j].swap(XichildComponent[j]);
            }
          }
        }

        for (unsigned j = 0; j < nChildren; ++j) {
          const unsigned jel = el->GetChild(iel, j);
          const unsigned jproc = mshp->IsdomBisectionSearch(jel, 3);
          el_send[jproc][posEl[jproc]++] = jel;

          for (unsigned i = 0; i < nDofsPsi; ++i) {
            if (l == levelStart) {
              const unsigned psiDof = msh->GetSolutionDof(i, iel, psiType);
              dofPsi_send[jproc][posPsi[jproc]] = psiDof;
              psiTarget_send[jproc][posPsi[jproc]] = (*sol3->_Sol[solPsiIndex])(psiDof);
            }
            else {
              dofPsi_send[jproc][posPsi[jproc]] = dofPsi_recv[offsetPsi[k] + i];
              psiTarget_send[jproc][posPsi[jproc]] = psiTarget_recv[offsetPsi[k] + i];
            }
            ++posPsi[jproc];
          }

          for (unsigned i = 0; i < nDofsXi; ++i) {
            for (unsigned d = 0; d < dim; ++d) {
              Xi_send[d * nProc + jproc][posXi[jproc]] = refined ? Xichild[d * nChildren + j][i] : coordXi[d][i];
            }
            ++posXi[jproc];
          }
        }
      }
    }

    if (l == levelF) {
      continue;
    }

    std::vector<std::vector<unsigned>> recvCounts(nProc, std::vector<unsigned>(3, 0));
    recvCounts[iproc] = sendCounts[iproc];

    for (unsigned p = 0; p < nProc; ++p) {
      if (p != iproc) {
        postRecv(recvCounts[p].data(), 3, MPI_UNSIGNED, p, tagCounts, requests);
      }
    }

    for (unsigned p = 0; p < nProc; ++p) {
      if (p != iproc) {
        postSend(sendCounts[p].data(), 3, MPI_UNSIGNED, p, tagCounts, requests);
      }
    }

    waitAll(requests);

    std::vector<unsigned> offsetElRecv(nProc + 1, 0);

    std::vector<unsigned> offsetPsiRecv(nProc + 1, 0);
    std::vector<unsigned> offsetXiRecv(nProc + 1, 0);

    for (unsigned p = 0; p < nProc; ++p) {
      offsetElRecv[p + 1] = offsetElRecv[p] + recvCounts[p][0];
      offsetPsiRecv[p + 1] = offsetPsiRecv[p] + recvCounts[p][1];
      offsetXiRecv[p + 1] = offsetXiRecv[p] + recvCounts[p][2];
    }

    el_recv.resize(offsetElRecv[nProc]);
    dofPsi_recv.resize(offsetPsiRecv[nProc]);
    psiTarget_recv.resize(offsetPsiRecv[nProc]);
    for (unsigned d = 0; d < dim; ++d) {
      Xi_recv[d].resize(offsetXiRecv[nProc]);
    }

    if (sendCounts[iproc][0]) {
      std::copy(el_send[iproc].begin(), el_send[iproc].end(), el_recv.data() + offsetElRecv[iproc]);
    }

    if (sendCounts[iproc][1]) {
      std::copy(dofPsi_send[iproc].begin(), dofPsi_send[iproc].end(), dofPsi_recv.data() + offsetPsiRecv[iproc]);
      std::copy(psiTarget_send[iproc].begin(), psiTarget_send[iproc].end(), psiTarget_recv.data() + offsetPsiRecv[iproc]);
    }

    if (sendCounts[iproc][2]) {
      for (unsigned d = 0; d < dim; ++d) {
        std::copy(Xi_send[d * nProc + iproc].begin(), Xi_send[d * nProc + iproc].end(), Xi_recv[d].data() + offsetXiRecv[iproc]);
      }
    }

    for (unsigned p = 0; p < nProc; ++p) {
      if (p == iproc) continue;
      if (recvCounts[p][0]) {
        postRecv(el_recv.data() + offsetElRecv[p], recvCounts[p][0], MPI_UNSIGNED, p, tagElements, requests);
      }
      if (recvCounts[p][1]) {
        postRecv(dofPsi_recv.data() + offsetPsiRecv[p], recvCounts[p][1], MPI_UNSIGNED, p, tagPsiDofs, requests);
        postRecv(psiTarget_recv.data() + offsetPsiRecv[p], recvCounts[p][1], MPI_DOUBLE, p, tagPsiValues, requests);
      }
      if (recvCounts[p][2]) {
        for (unsigned d = 0; d < dim; ++d) {
          postRecv(Xi_recv[d].data() + offsetXiRecv[p], recvCounts[p][2], MPI_DOUBLE, p, tagXi + static_cast<int>(d), requests);
        }
      }
    }

    for (unsigned p = 0; p < nProc; ++p) {
      if (p == iproc) continue;
      postSend(el_send[p].data(), sendCounts[p][0], MPI_UNSIGNED, p, tagElements, requests);
      postSend(dofPsi_send[p].data(), sendCounts[p][1], MPI_UNSIGNED, p, tagPsiDofs, requests);
      postSend(psiTarget_send[p].data(), sendCounts[p][1], MPI_DOUBLE, p, tagPsiValues, requests);
      for (unsigned d = 0; d < dim; ++d) {
        postSend(Xi_send[d * nProc + p].data(), sendCounts[p][2], MPI_DOUBLE, p, tagXi + static_cast<int>(d), requests);
      }
    }

    waitAll(requests);
  }

  RES->close();
  KK->close();
}
