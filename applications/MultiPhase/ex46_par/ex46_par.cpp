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
#include "../includeLS/TangentialSurfaceFactors.hpp"
#include <memory>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#include "petsc.h"
#include "petscmat.h"
#include "PetscMatrix.hpp"
#include "Fem.hpp"

#include "MyMatrix.hpp"
#include "VTKWriter.hpp"

#include "SurfacePressureRangeDiagnostic.hpp"

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

// Dedicated static-droplet / parasitic-current benchmark.
bool printdb = false;
bool parVelocityGhostPenalty = false;  // matches the supplied ex46.cpp
bool parPressureGhostPenalty = true;   // matches the supplied ex46.cpp
bool parVolumetricStabilization = false; // matches the supplied ex46.cpp
double parSphereRadius = 0.002;


// const RungeKutta::VelKind velocityType = RungeKutta::VelKind::Zero;
//const RungeKutta::VelKind velocityType = RungeKutta::VelKind::Translation;
//const RungeKutta::VelKind velocityType = RungeKutta::VelKind::Rotation;
//const RungeKutta::VelKind velocityType = RungeKutta::VelKind::Vortex;

#include "../includeLS/Utils.hpp"

std::vector<std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > > PMatrix(6); //projection matrix between elements;

CutFemIntegration cutFemIntegration;
double surfaceMomentTolerance = 1.e-8;

namespace {

SurfacePressureRangeDiagnostic* activeSurfaceRange = nullptr;

class SurfaceRangeCaptureGuard {
public:
  explicit SurfaceRangeCaptureGuard(
    SurfacePressureRangeDiagnostic& diagnostic) {

    if (activeSurfaceRange != nullptr) {
      throw std::logic_error(
        "Surface range diagnostic capture already active");
    }

    activeSurfaceRange = &diagnostic;
  }

  ~SurfaceRangeCaptureGuard() {
    activeSurfaceRange = nullptr;
  }

  SurfaceRangeCaptureGuard(const SurfaceRangeCaptureGuard&) = delete;

  SurfaceRangeCaptureGuard& operator=(
    const SurfaceRangeCaptureGuard&) = delete;
};

void ExportSurfacePressureDofs(
    SurfacePressureRangeDiagnostic& diagnostic,
    MultiLevelSolution& mlSol,
    TransientNonlinearImplicitSystem& system,
    unsigned level,
    unsigned dim) {

  Solution* solution = mlSol.GetSolutionLevel(level);

  LinearEquationSolver* linear = system._LinSolver[level];

  const char* velocityNames[] = {"U", "V", "W"};
  const char* pressureNames[] = {"P1", "P2"};

  for (unsigned field = 0; field < dim + 2; ++field) {

    const char* name =
      field < dim
        ? velocityNames[field]
        : pressureNames[field-dim];

    const unsigned solIndex = mlSol.GetIndex(name);

    const unsigned pdeIndex = system.GetSolPdeIndex(name);

    NumericVector* vector = solution->_Sol[solIndex];

    NumericVector* bdcVector = solution->_Bdc[solIndex];

    for (unsigned dof = vector->first_local_index();
         dof < vector->last_local_index();
         ++dof) {

      const unsigned systemDof =
        linear->GetSystemDof(solIndex, pdeIndex, dof);

      // const double bdc = (*bdcVector)(dof);
      const bool isPressure = (field >= dim);

      const double bdc = isPressure
                 ? 2.0
                 : (*bdcVector)(dof);

      // FEMuS: 2 = free, 1 = AMR constrained, 0 = Dirichlet.
      const double fixedValue =
        bdc > 1.5 || (bdc > 0.5 && bdc < 1.5)
          ? 0.0
          : (*vector)(dof);

      diagnostic.AddDof(
        systemDof,
        name,
        bdc,
        fixedValue
      );
    }
  }
}

} // namespace

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
// No level-set transport, reinitialization, or PSI solve in this benchmark.
// double TimeStepMultiphase(const double time);

void SurfaceReferenceGradients(const elem_type* element, const std::vector<double>& point, std::vector<double>& gradient) {
  const unsigned dimension = point.size();
  std::vector<double> derivative;
  element->GetDPhiDXi(derivative,point);
  gradient.assign(derivative.size()*dimension,0.);
  for(unsigned d = 0; d < dimension; ++d) {
    if(d == 1) element->GetDPhiDEta(derivative,point);
    if(d == 2) element->GetDPhiDZeta(derivative,point);
    for(unsigned i = 0; i < derivative.size(); ++i) gradient[i*dimension+d] = derivative[i];
  }
}

void SurfaceReferenceGradients(const elem_type* element, unsigned gauss, std::vector<double>& gradient) {
  const unsigned dimension = element->GetDim();
  const unsigned count = element->GetNDofs();
  gradient.resize(count*dimension);
  const double* derivative[3] = {element->GetDPhiDXi(gauss),dimension > 1 ? element->GetDPhiDEta(gauss) : nullptr,dimension > 2 ? element->GetDPhiDZeta(gauss) : nullptr};
  for(unsigned i = 0; i < count; ++i) for(unsigned d = 0; d < dimension; ++d) gradient[i*dimension+d] = derivative[d][i];
}

// Parasitic-current metrics on the actual NS velocity space.
struct ParVelocityNorms {
  double l2 = 0.;
  double maxGauss = 0.;
  double maxNodal = 0.;
};

ParVelocityNorms MeasureParVelocity(MultiLevelSolution& mlSol, unsigned level) {
  Mesh* mesh = mlSol.GetMultilevelMesh()->GetLevel(level);
  Solution* sol = mlSol.GetSolutionLevel(level);
  const unsigned dim = mesh->GetDimension();
  const unsigned rank = mesh->processor_id();
  const unsigned velocityType = mlSol.GetSolutionType("U");
  const char* names[3] = {"U", "V", "W"};
  std::vector<NumericVector*> vectors(dim);
  for(unsigned d = 0; d < dim; ++d) {
    vectors[d] = sol->_Sol[mlSol.GetIndex(names[d])];
  }
  double localIntegral = 0.;
  double localGaussMax = 0.;
  double localNodalMax = 0.;
  const unsigned begin = vectors[0]->first_local_index();
  const unsigned end = vectors[0]->last_local_index();
  for(unsigned dof = begin; dof < end; ++dof) {
    double squared = 0.;
    for(unsigned d = 0; d < dim; ++d) {
      const double v = (*vectors[d])(dof);
      squared += v*v;
    }
    localNodalMax = std::max(localNodalMax,std::sqrt(squared));
  }
  for(unsigned iel = mesh->_elementOffset[rank];
      iel < mesh->_elementOffset[rank+1U]; ++iel) {
    const unsigned type = mesh->GetElementType(iel);
    const unsigned count = mesh->GetElementDofNumber(iel,velocityType);
    std::vector<std::vector<double>> x(dim,std::vector<double>(count));
    std::vector<std::vector<double>> u(dim,std::vector<double>(count));
    for(unsigned j = 0; j < count; ++j) {
      const unsigned xdof = mesh->GetSolutionDof(j,iel,2);
      const unsigned vdof = mesh->GetSolutionDof(j,iel,velocityType);
      for(unsigned d = 0; d < dim; ++d) {
        x[d][j] = (*mesh->_topology->_Sol[d])(xdof);
        u[d][j] = (*vectors[d])(vdof);
      }
    }
    const elem_type* fe = mesh->_finiteElement[type][velocityType];
    std::vector<double> phi,deriv;
    for(unsigned ig = 0; ig < fe->GetGaussPointNumber(); ++ig) {
      double weight = 0.;
      fe->Jacobian(x,ig,weight,phi,deriv);
      double speed2 = 0.;
      for(unsigned d = 0; d < dim; ++d) {
        double ug = 0.;
        for(unsigned j = 0; j < count; ++j) ug += phi[j]*u[d][j];
        speed2 += ug*ug;
      }
      localIntegral += weight*speed2;
      localGaussMax = std::max(localGaussMax,std::sqrt(speed2));
    }
  }
  double integral=0.,maxGauss=0.,maxNodal=0.;
  MPI_Allreduce(&localIntegral,&integral,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&localGaussMax,&maxGauss,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(&localNodalMax,&maxNodal,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  ParVelocityNorms result;
  result.l2 = std::sqrt(std::max(0.,integral));
  result.maxGauss = maxGauss;
  result.maxNodal = maxNodal;
  return result;
}

struct ParOptions {
  unsigned dim = 3;
  unsigned uniformLevels = 1;
  unsigned nStepsPerPeriod = 1000; // dt=period/nStepsPerPeriod
  unsigned solveSteps = 1;         // number of actual NS steps
  unsigned nx = 8,ny = 8,nz = 8;
  unsigned momentOrder = 10;
  unsigned surfaceWeights = 3;
  double momentTolerance = 1.e-10;
  double meshScale = 100.;
  std::string mesh = "imported";  // generated | imported
  std::string meshFile;
  std::string initialPressure = "young-laplace"; // young-laplace | zero
  std::string outputPrefix;
  std::string pressurePrefix;
  bool surfacePressureDiagnostic = false;
  bool pressureGP = true;
  bool velocityGP = false;
  bool stabilization = false;
};

ParOptions ParseParOptions(int argc, char** argv) {
  ParOptions opt;
  auto need = [&](int& i,const std::string& name) -> std::string {
    if(i+1 >= argc) throw std::invalid_argument("Missing value for "+name);
    return argv[++i];
  };
  auto yesno = [&](int& i,const std::string& name) -> bool {
    const std::string v = need(i,name);
    if(v != "0" && v != "1") throw std::invalid_argument(name+" must be 0 or 1");
    return v=="1";
  };
  for(int i=1;i<argc;++i) {
    const std::string s = argv[i];
    if(s=="--dim") opt.dim = std::stoul(need(i,s));
    else if(s=="--uniform-levels") opt.uniformLevels = std::stoul(need(i,s));
    else if(s=="--nsteps") opt.nStepsPerPeriod = std::stoul(need(i,s));
    else if(s=="--solve-steps") opt.solveSteps = std::stoul(need(i,s));
    else if(s=="--nx") opt.nx = std::stoul(need(i,s));
    else if(s=="--ny") opt.ny = std::stoul(need(i,s));
    else if(s=="--nz") opt.nz = std::stoul(need(i,s));
    else if(s=="--mesh") opt.mesh = need(i,s);
    else if(s=="--mesh-file") opt.meshFile = need(i,s);
    else if(s=="--mesh-scale") opt.meshScale = std::stod(need(i,s));
    else if(s=="--initial-pressure") opt.initialPressure = need(i,s);
    else if(s=="--surface-moment-order") opt.momentOrder = std::stoul(need(i,s));
    else if(s=="--surface-moment-tolerance") opt.momentTolerance = std::stod(need(i,s));
    else if(s=="--surface-weights") opt.surfaceWeights = std::stoul(need(i,s));
    else if(s=="--surface-pressure-diagnostic") opt.surfacePressureDiagnostic = yesno(i,s);
    else if(s=="--surface-pressure-prefix") opt.pressurePrefix = need(i,s);
    else if(s=="--output-prefix") opt.outputPrefix = need(i,s);
    else if(s=="--pressure-gp") opt.pressureGP = yesno(i,s);
    else if(s=="--velocity-gp") opt.velocityGP = yesno(i,s);
    else if(s=="--stabilization") opt.stabilization = yesno(i,s);
    else if(s=="--gk-regular-2d"||s=="--gk-tangent-2d"||
            s=="--gk-regular-3d"||s=="--gk-tangent-3d") { (void) need(i,s); }
    else if(s=="-options_file"||s=="-matptap_via") { (void) need(i,s); }
    else if(s=="--help"||s=="-h") {
      std::cout << "PAR stationary sphere only (Cartesian, Q2 velocity / Q1 pressure).\n"
        << "  --mesh generated|imported  --mesh-file PATH  --mesh-scale S\n"
        << "  --dim 2|3  --uniform-levels N  --nx N --ny N --nz N\n"
        << "  --nsteps N (sets dt=period/N)  --solve-steps N (actual solves)\n"
        << "  --initial-pressure zero|young-laplace\n"
        << "  --pressure-gp 0|1 --velocity-gp 0|1 --stabilization 0|1\n"
        << "  --surface-moment-order Q --surface-moment-tolerance T --surface-weights 1|3\n"
        << "  --surface-pressure-diagnostic 0|1 --surface-pressure-prefix PREFIX\n"
        << "  --output-prefix PREFIX\n"
        << FELevelSetMoments::GaussKronrodRules::Help();
      std::exit(0);
    }
    // PETSc-specific command-line options are intentionally ignored here.
  }
  if(opt.dim!=2 && opt.dim!=3) throw std::invalid_argument("--dim must be 2 or 3");
  if(opt.mesh!="generated"&&opt.mesh!="imported") throw std::invalid_argument("--mesh must be generated or imported");
  if(opt.initialPressure!="zero"&&opt.initialPressure!="young-laplace") throw std::invalid_argument("--initial-pressure must be zero or young-laplace");
  if(!opt.uniformLevels||!opt.nStepsPerPeriod||!opt.solveSteps||!opt.nx||!opt.ny||!opt.nz) throw std::invalid_argument("Level/step/mesh counts must be positive");
  if(opt.momentOrder<1||opt.momentOrder>18) throw std::invalid_argument("Moment order must be 1..18");
  if(opt.surfaceWeights!=1&&opt.surfaceWeights!=3) throw std::invalid_argument("Surface weights must be 1 or 3");
  if(!(opt.momentTolerance>0.)||!std::isfinite(opt.momentTolerance)) throw std::invalid_argument("Invalid moment tolerance");
  if(!(opt.meshScale>0.)||!std::isfinite(opt.meshScale)) throw std::invalid_argument("Invalid mesh scale");
  if(opt.meshFile.empty()) opt.meshFile = opt.dim==2 ?
    "./input/unstructured_PC_coarse_Q1.neu" : "./input/unstructured_cube_hex27_Q1.neu";
  if(opt.outputPrefix.empty()) opt.outputPrefix =
    "levelsetdiagnostic/par_" + opt.mesh + "_D" + std::to_string(opt.dim)
    + "_U" + std::to_string(opt.uniformLevels) + "_NS" + std::to_string(opt.nStepsPerPeriod);
  if(opt.pressurePrefix.empty()) opt.pressurePrefix = opt.outputPrefix + "_surface_pressure";
  return opt;
}

int main(int argc,char** argv) {
  FemusInit mpinit(argc,argv,MPI_COMM_WORLD);
  int rank = 0,nprocs = 1;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

  ParOptions opt;
  try {
    opt = ParseParOptions(argc,argv);
    surfaceQuadratureRules = FELevelSetMoments::OuterQuadratureRules::ParseArguments(argc,argv);
  }
  catch(const std::exception& e) {
    if(rank==0) std::cerr << "PAR option error: " << e.what() << '\n';
    return 1;
  }
  // This is a static interface benchmark. No AMR, level offsets,
  // surface advection, reinitialization or other simulation cases.
  const unsigned dim = opt.dim;
  const unsigned levelC = opt.uniformLevels-1U;
  const unsigned levelF = levelC;
  const unsigned levelPsi = levelF;
  const unsigned level0 = 0;
  const SimulationCase simulation = SimulationCase::PAR;
  static const SimulationConfig cfg = GetSimulationConfig(
    SimulationCase::PAR,dim,false,false);
  const MultiphasePhysicalProperties properties = cfg.properties;
  parSphereRadius = cfg.psi.radius;
  parPressureGhostPenalty = opt.pressureGP;
  parVelocityGhostPenalty = opt.velocityGP;
  parVolumetricStabilization = opt.stabilization;
  surfaceMomentOrder = opt.momentOrder;
  surfaceMomentTolerance = opt.momentTolerance;
  surfaceWeightCount = opt.surfaceWeights;
  cutFemIntegration = CutFemIntegration::Quadratic;
  for(auto* rule:cfwLS) {
    rule->SetIntegrationDirection(0);
    rule->SetDirectionValidation(0);
    rule->SetOuterQuadratureRules(surfaceQuadratureRules);
  }

  const double xmin=cfg.x_min[0], xmax=cfg.x_max[0];
  const double ymin=cfg.x_min[1], ymax=cfg.x_max[1];
  const double zmin=cfg.x_min[2], zmax=cfg.x_max[2];
  const ElemType elementType = dim==2 ? QUAD9 : HEX27;
  if(rank==0) {
    std::cout << "PAR static-interface benchmark: dim=" << dim
              << ", mesh=" << opt.mesh << ", ranks=" << nprocs
              << ", levels=" << opt.uniformLevels
              << ", actual solve steps=" << opt.solveSteps << '\n';
    std::cout << "  initial pressure=" << opt.initialPressure
              << ", GP pressure=" << opt.pressureGP
              << ", GP velocity=" << opt.velocityGP
              << ", stabilization=" << opt.stabilization << '\n';
    std::cout << "  CutFEM q=" << surfaceMomentOrder
              << ", tolerance=" << std::setprecision(16) << surfaceMomentTolerance
              << ", surface weights=" << surfaceWeightCount << '\n';
    if(opt.mesh=="imported") std::cout << "  input=" << opt.meshFile
              << ", scale=" << opt.meshScale << '\n';
    else std::cout << "  generated cells=" << opt.nx << 'x' << opt.ny
                   << 'x' << (dim==3?opt.nz:1) << '\n';
  }

  // Keep the same free-slip boundary conditions as the PAR case in ex46.
  auto SetBoundaryCondition = [](
    const std::vector<double>& x,const char SolName[],double& value,
    const int /*face*/,const double /*time*/) -> bool {
    value = 0.;
    int component = -1;
    if(!strcmp(SolName,"U")) component=0;
    else if(!strcmp(SolName,"V")) component=1;
    else if(!strcmp(SolName,"W")) component=2;
    if(component<0) return false; // pressure and Psi: no external Dirichlet
    const double tol=1.e-10;
    const unsigned bcDim=cfg.boundary.min.size();
    for(unsigned d=0;d<bcDim;++d) {
      if(x[d]<=cfg.x_min[d]+tol &&
         cfg.boundary.min[d].type==BoundaryType::SLIP &&
         static_cast<unsigned>(component)==d) return true;
      if(x[d]>=cfg.x_max[d]-tol &&
         cfg.boundary.max[d].type==BoundaryType::SLIP &&
         static_cast<unsigned>(component)==d) return true;
    }
    return false;
  };

  MultiLevelMesh mlMsh0;
  if(opt.mesh=="generated") {
    mlMsh0.GenerateCoarseBoxMesh(opt.nx,opt.ny,dim==3?opt.nz:0,
      xmin,xmax,ymin,ymax,zmin,zmax,elementType,"seventh");
  }
  else {
    mlMsh0.ReadCoarseMesh(opt.meshFile.c_str(),"seventh",opt.meshScale);
  }
  if(mlMsh0.GetDimension()!=dim) throw std::runtime_error("Mesh dimension mismatch");
  InitializeSurfaceQuadrature(*mlMsh0.GetLevel(0));
  mlMsh0.RefineMesh(opt.uniformLevels,opt.uniformLevels,nullptr);
  mlMsh0.PrintInfo();

  BBoxToIel bbox(mlMsh0, 0, 3);

  MultiLevelSolution mlSol0(&mlMsh0);
  const char* velocities[] = {"U","V","W"};
  const char* pressures[] = {"P1","P2"};
  std::vector<std::string> vName;
  for(unsigned d=0;d<dim;++d) vName.emplace_back(velocities[d]);
  const std::vector<std::string> pName = {"P1","P2"};
  const std::string psiName="Psi",cName="C";
  mlSol0.AddSolution("Psi",LAGRANGE,SECOND,0,false);
  for(unsigned d=0;d<dim;++d) mlSol0.AddSolution(velocities[d],LAGRANGE,SECOND,2);
  for(unsigned d=0;d<2;++d) mlSol0.AddSolution(pressures[d],LAGRANGE,FIRST,2);
  mlSol0.AddSolution("C",DISCONTINUOUS_POLYNOMIAL,ZERO,0,false);
  mlSol0.Initialize("All");
  // Stationary Young--Laplace state: enforce u_h(x,0)=0 independently
  // of the general-purpose prescribed-velocity function in RungeKutta.hpp.
  for(unsigned l=0;l<=levelC;++l) {
    Solution* initial = mlSol0.GetSolutionLevel(l);
    for(unsigned d=0;d<dim;++d) {
      NumericVector* u = initial->_Sol[mlSol0.GetIndex(velocities[d])];
      u->zero();
      u->close();
    }
  }
  Psi psi(cfg.psi);
  InitLevelSet(mlSol0,psiName,psi);
  UpdateColorFunction(mlSol0,psiName,cName);
  mlSol0.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol0.GenerateBdc("All");

  // Build coarse-level phase classifications from the finest level.
    if(levelC > 0)
    RestrictPWDCField(mlSol0, cName, 0, levelC);

    // Temporary solution receiving the transported level set.
    // It shares the same fixed multilevel mesh.
    MultiLevelSolution transportedPsi(&mlMsh0);

    transportedPsi.AddSolution(
    psiName.c_str(), LAGRANGE, SECOND, 0, false
    );

    transportedPsi.Initialize("All");

  const double hC=GetMaxElementH(mlMsh0.GetLevel(levelC),levelC);
  const double dt=cfg.period/opt.nStepsPerPeriod;
  static double parDt=0.;
  parDt=dt;
  auto timeStep = [](const double) -> double { return parDt; };
  auto firstTimeStep = [](const double) -> double { return 0.5*parDt; };

  // Write an initial snapshot to compare the two mesh families.
  const std::vector<std::string> allVariables = {"All"};
  {
    VTKWriter writer(&mlSol0);
    writer.Write(levelF,"output_par_ls","biquadratic",allVariables,0);
  }
  if(rank==0) {
    mkdir("levelsetdiagnostic",0755);
    const std::string metricsFile=opt.outputPrefix+"_currents.dat";
    std::ofstream out(metricsFile,std::ios::trunc);
    if(!out) throw std::runtime_error("Cannot open "+metricsFile);
    out << "# time L2_u max_gauss_u max_nodal_u\n";
  }
  MPI_Barrier(MPI_COMM_WORLD);

  InitializeParasiticCurrents(levelC);

  for(unsigned step=1;step<=opt.solveSteps;++step) {
    mlSol0.CopySolutionToOldSolution();
    const TimeDiscretization td=step==1 ?
      TimeDiscretization::BackEuler : TimeDiscretization::CrankNicholson;
    MultiLevelMesh mlMsh2(mlMsh0,level0,levelC+1U,"seventh");
    MultiLevelSolution mlSol2(&mlMsh2);
    mlSol2.AddSolution("Psi",LAGRANGE,SECOND,0,false);
    for(unsigned d=0;d<dim;++d) mlSol2.AddSolution(velocities[d],LAGRANGE,SECOND,2);
    for(unsigned d=0;d<2;++d) mlSol2.AddSolution(pressures[d],LAGRANGE,FIRST,2);
    mlSol2.AddSolution("C",DISCONTINUOUS_POLYNOMIAL,ZERO,0,false);
    mlSol2.Initialize("All");
    mlSol2.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
    mlSol2.GenerateBdc("All");
    Solution* src=mlSol0.GetSolutionLevel(levelC);
    Solution* dst=mlSol2.GetSolutionLevel(levelC);
    for(unsigned d=0;d<dim;++d) {
      const unsigned i0=mlSol0.GetIndex(velocities[d]);
      const unsigned i2=mlSol2.GetIndex(velocities[d]);
      *(dst->_Sol[i2])=*(src->_Sol[i0]);
      *(dst->_SolOld[i2])=*(src->_SolOld[i0]);
    }
    for(const char* pname:pressures) {
      *(dst->_Sol[mlSol2.GetIndex(pname)])=*(src->_Sol[mlSol0.GetIndex(pname)]);
    }

    *(dst->_Sol[mlSol2.GetIndex("Psi")]) = *(src->_Sol[mlSol0.GetIndex("Psi")]);

    const unsigned c0Index = mlSol0.GetIndex(cName.c_str());
    const unsigned c2Index = mlSol2.GetIndex(cName.c_str());

    for(unsigned l = 0; l <= levelC; ++l) *(mlSol2.GetSolutionLevel(l)->_Sol[c2Index]) = *(mlSol0.GetSolutionLevel(l)->_Sol[c0Index]);

    // Production pressure constraints are essential even when the range
    // diagnostic exports its pressure space as unconstrained.
    std::vector<double> gaugePoint={xmin,ymin,zmin};
    gaugePoint.resize(dim);
    SetUnphysicalPressureDofs(mlSol2,cName,pName,0,levelC,gaugePoint,true);
    if(step==1 && opt.initialPressure=="young-laplace") {
      const unsigned iP1=mlSol2.GetIndex("P1");
      NumericVector* p1=dst->_Sol[iP1];
      NumericVector* bdc=dst->_Bdc[iP1];
      const double jump=(dim==3?2.:1.)*properties.sigma/cfg.psi.radius;
      for(unsigned dof=p1->first_local_index();dof<p1->last_local_index();++dof) {
        if((*bdc)(dof)>1.5) p1->set(dof,jump);
      }
      p1->close();
    }

    MultiLevelProblem problem(&mlSol2);
    TransientNonlinearImplicitSystem& system=
      problem.add_system<TransientNonlinearImplicitSystem>("NS");
    for(const std::string& name:vName) system.AddSolutionToSystemPDE(name.c_str());
    for(const std::string& name:pName) system.AddSolutionToSystemPDE(name.c_str());

    auto extraSparsityFunction = [&] (
        LinearEquation& linEq,
        std::vector<std::map<int,bool>>& d,
        std::vector<std::map<int,bool>>& o,
        std::map<int,std::map<int,bool>>& od,
        std::map<int,std::map<int,bool>>& dd) {
      if(linEq._msh->GetLevel()!=levelC) return;
      for(unsigned idim=0;idim<dim;++idim)
        GetGPSparsityPattern(problem,cName,vName[idim],levelC,levelC,0,d,o,od,dd);
      for(unsigned idim=0;idim<pName.size();++idim)
        GetGPSparsityPattern(problem,cName,pName[idim],levelC,levelC,idim+1,d,o,od,dd);
    };
    system.SetExtraSparsityFunction(extraSparsityFunction);
    system.SetAssembleFunction(AssembleMultiphase);
    system.AttachGetTimeIntervalFunction(step==1?firstTimeStep:timeStep);
    problem.SetMultiphaseParams(&mlSol0,levelPsi,levelC,level0,hC,
      properties,td,true,false,cfg.x_min[0]);
    for(unsigned l=0;l<=levelC;++l) mlMsh2.GetLevel(l)->SetLevel(l);

    std::vector<unsigned> velocityFields(dim),pressureFields(pName.size());
    for(unsigned d=0;d<dim;++d)
      velocityFields[d]=system.GetSolPdeIndex(vName[d].c_str());
    for(unsigned d=0;d<pName.size();++d)
      pressureFields[d]=system.GetSolPdeIndex(pName[d].c_str());
    FieldSplitTree FS_V(PREONLY,MULTIGRID_PRECOND,velocityFields,"Velocity");
    FS_V.MgSetSolverFineGrids(CHEBYSHEV);
    FS_V.MgSetPreconditionerFineGrids(JACOBI_PRECOND);
    FS_V.MgSetParameters(4,4,0);
    FieldSplitTree FS_P(PREONLY,MULTIGRID_PRECOND,pressureFields,"Pressure");
    FS_P.MgSetSolverFineGrids(CHEBYSHEV);
    FS_P.MgSetPreconditionerFineGrids(JACOBI_PRECOND);
    FS_P.MgSetParameters(4,4,0);
    std::vector<FieldSplitTree*> splits={&FS_V,&FS_P};
    FieldSplitTree FS_NS(GMRES,FIELDSPLIT_SCHUR_PRECOND,splits,"Navier-Stokes");
    FS_NS.SetTolerances(1.e-14,1.e-14,1.e+50,50);
    system.SetLinearEquationSolverType(FEMuS_FIELDSPLIT,INCLUDE_COARSE_LEVEL_TRUE);
    system.SetMaxNumberOfNonLinearIterations(1);
    system.SetMaxNumberOfLinearIterations(10);
    system.SetMgType(V_CYCLE);
    system.init();
    system.SetOuterSolver(FGMRES);
    system.SetFieldSplitTree(&FS_NS);
    system.SetTolerances(1.e-14,1.e-14,1.e+50,200);
    system.SetAbsoluteLinearConvergenceTolerance(1.e-14);
    for(unsigned l=0;l<=levelC;++l)
      system._LinSolver[l]->SetSolution(&mlSol2.GetSolutionLevel(l)->_Sol);

    if(opt.surfacePressureDiagnostic && step==1) {
      const std::string prefix=opt.pressurePrefix;
      SurfacePressureRangeDiagnostic diagnostic;
      diagnostic.Begin(prefix,rank,dim,
        (dim==3?2.:1.)*properties.sigma/cfg.psi.radius,surfaceWeightCount);
      ExportSurfacePressureDofs(diagnostic,mlSol2,system,levelC,dim);
      {
        SurfaceRangeCaptureGuard guard(diagnostic);
        system.Solve();
      }
      if(!diagnostic.IsComplete())
        throw std::runtime_error("Surface-pressure diagnostic capture incomplete");
    }
    else system.Solve();

    const ParVelocityNorms norms=MeasureParVelocity(mlSol2,levelC);
    const double time = (static_cast<double>(step) - 0.5) * dt;

    PrintParasiticCurrents(&mlSol2, levelC, levelC, time, vName, cfg,cutFemIntegration);

    if(rank==0) {
      std::cout << "PAR step " << step << '/' << opt.solveSteps
                << ", t=" << std::setprecision(15) << time
                << ": L2(u)=" << norms.l2
                << ", max Gauss |u|=" << norms.maxGauss
                << ", max nodal |u|=" << norms.maxNodal << '\n';
      std::ofstream out(opt.outputPrefix+"_currents.dat",std::ios::app);
      out << std::setprecision(17) << time << ' '
          << norms.l2 << ' ' << norms.maxGauss << ' ' << norms.maxNodal << '\n';
    }
    // The pressure and velocity are advanced, but Psi and C remain fixed.
    // for(unsigned d=0;d<dim;++d) {
    //   const unsigned i0=mlSol0.GetIndex(velocities[d]);
    //   const unsigned i2=mlSol2.GetIndex(velocities[d]);
    //   *(src->_Sol[i0])=*(dst->_Sol[i2]);
    // }
    // for(const char* pname:pressures)
    //   *(src->_Sol[mlSol0.GetIndex(pname)])=
    //     *(dst->_Sol[mlSol2.GetIndex(pname)]);
    // {
    //   VTKWriter writer(&mlSol2);
    //   writer.Write(levelC,"output_par_vel","biquadratic",allVariables,step);
    // }


    for(unsigned d = 0; d < dim; ++d) {
        const unsigned i0 = mlSol0.GetIndex(velocities[d]);
        const unsigned i2 = mlSol2.GetIndex(velocities[d]);

        *(src->_Sol[i0]) = *(dst->_Sol[i2]);
    }

    for(const char* pname : pressures) *(src->_Sol[mlSol0.GetIndex(pname)]) = *(dst->_Sol[mlSol2.GetIndex(pname)]);

    const double dtTransport = (step == 1) ? 0.5 * dt : dt;

    ProjectSolution(mlSol0,transportedPsi,bbox,cfg.x_min,cfg.x_max,{psiName},levelC,levelC,vName,levelC,Boundary(),-dtTransport,time,cfg.period);

    const unsigned psi0Index = mlSol0.GetIndex(psiName.c_str());

    const unsigned psiTransportIndex = transportedPsi.GetIndex(psiName.c_str());

    *(src->_Sol[psi0Index]) = *(transportedPsi.GetSolutionLevel(levelC)->_Sol[psiTransportIndex]);

    UpdateColorFunction(mlSol0, psiName, cName);

    if(levelC > 0)
    RestrictPWDCField(mlSol0, cName, 0, levelC);

    for(auto* rule : cfwLS) if(rule != nullptr) rule->ClearCache();

    // {
    //     VTKWriter writer(&mlSol0);
    //     writer.Write(levelC,"output_par_vel","biquadratic",allVariables,step);
    // }
  }
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

  if(parVelocityGhostPenalty) AssembleGhostPenaltyVelocity(ml_prob2);

  if (parPressureGhostPenalty && solPType == 3) {
    AssembleGhostPenaltyDGP(ml_prob2, true);
    AssembleGhostPenaltyDGP(ml_prob2, false);
  }
  else if (parPressureGhostPenalty && solPType == 0) {
    AssembleGhostPenaltyLinearPressure(ml_prob2, true);
    AssembleGhostPenaltyLinearPressure(ml_prob2, false);
  }
  else if (solPType != 3 && solPType != 0) {
    throw std::runtime_error("Pressure type must be 3 or 0");
  }

  if(parVolumetricStabilization) AssembleStabilizationTerms(ml_prob2);

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

  const bool captureSurfaceRange = activeSurfaceRange != nullptr && activeSurfaceRange->IsCapturing() && levelToAssemble == levelC - level0 && firstNonlinearIt;

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

        std::vector<double> surfaceResidual;
        std::vector<double> pressureActionYL;

        if (captureSurfaceRange) {
          surfaceResidual.assign(dim * nDofsV, 0.0);
          pressureActionYL.assign(dim * nDofsV, 0.0);
        }

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
        const bool tangentialSurfaceWeights = cutFemIntegration == CutFemIntegration::Quadratic && !axisymmetric;
        // CutFemWeight_LS<double,long double>* elementSurfaceRule = cfwLS[ielGeom];

        CutFemWeight_LS<double,long double>* elementSurfaceRule = cfwLS[ielGeom];
        CutFemWeight_LS<TypeIO,SurfaceFitPrecision>* elementSurfaceQuadRule = nullptr;
        const double* elementSurfaceGaussWeights = elementSurfaceRule->GetGaussWeightPointer();

        if(cut == 1) { // this can happen only on l == levelF
          if (l != levelF) throw std::runtime_error("CUT CELL SHOULD BE ONLY AT LEVELF\n");
          Fem& integrationFem = tangentialSurfaceWeights ? *surfaceFem : fem;
          femV = integrationFem.GetFiniteElement(ielGeom,solVType);
          femP = integrationFem.GetFiniteElement(ielGeom,solPType);
          femXi = integrationFem.GetFiniteElement(ielGeom,solXiType);
          femPsi = integrationFem.GetFiniteElement(ielGeom,psiType);
          // if(tangentialSurfaceWeights) {
          //   elementSurfaceRule = surfaceRules.at(ielGeom).get();
          // }
          if(tangentialSurfaceWeights) {
            elementSurfaceQuadRule = surfaceRules.at(ielGeom).get();
            elementSurfaceGaussWeights = elementSurfaceQuadRule->GetGaussWeightPointer();
          }

          nDofsPsi = msh->GetElementDofNumber(iel, psiType);
          psi.resize(nDofsPsi);
          for(unsigned i = 0; i < nDofsPsi; i++) {
            unsigned psiDof = msh->GetSolutionDof(i, iel, psiType);
            psi[i] = (*solPsi)(psiDof);
          }

          if(cutFemIntegration == CutFemIntegration::Linear) {
            const unsigned ng = femPsi->GetGaussPointNumber();
            std::vector<double> psig(ng,0.);
            for(unsigned ig = 0; ig < ng; ++ig) for(unsigned i = 0; i < nDofsPsi; ++i) psig[ig] += psi[i]*femPsi->GetPhi(ig)[i];
            std::vector<const double*> xg(dim);
            for(unsigned axis = 0; axis < dim; ++axis) xg[axis] = femPsi->GetGaussRule().GetGaussCoordinatePointer(axis);
            BestFitLinearInterpolation(xg, psig, a);
            d = a[dim];
            a.resize(dim);
          }
        }

        
        const unsigned nGauss = femV->GetGaussPointNumber();

        std::vector<TypeIO> weightCF(nGauss, 0.);
        std::vector<std::vector<TypeIO>> weightCFNormal;
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
          // else {

          //   CutFemWeight_LS<double,long double>* rule = elementSurfaceRule;

          //   if(rule == nullptr) {
          //     throw std::runtime_error("Unsupported quadratic CutFEM element.");
          //   }

          //   if(rule->GetGaussQuadraturePointNumber() != nGauss) {
          //     throw std::runtime_error("Quadratic CutFEM Gauss rule mismatch");
          //   }

          //   const double tol = surfaceMomentTolerance;
          //   if(tangentialSurfaceWeights) {
          //     rule->GetTangentialSurfaceWeightsFromNodalGeometry(psi,coordX,weightCFInt,weightCF,tol,tol);
          //     if(dim == 3 && surfaceWeightCount == 3) rule->GetNormalSurfaceWeightsFromNodalGeometry(psi,coordX,weightCFInt,weightCFNormal,tol,tol);
          //   }
          //   else rule->GetPhysicalSurfaceWeightsFromNodalGeometry(psi,coordX,weightCFInt,weightCF,tol,tol);
          // }
          else {

            const double tol = surfaceMomentTolerance;
            const CutFemLevelSetData levelSet(static_cast<GeomElType>(ielGeom),psi);

            if(tangentialSurfaceWeights) {
              auto* rule = elementSurfaceQuadRule;
              if(rule == nullptr) throw std::runtime_error("Unsupported quadratic CutFEM element.");
              if(rule->GetGaussQuadraturePointNumber() != nGauss) throw std::runtime_error("Quadratic CutFEM Gauss rule mismatch");
              if(dim == 3 && surfaceWeightCount == 3) {
                try { rule->GetTangentialAndNormalSurfaceWeightsFromLevelSetData(levelSet,coordX,weightCFInt,weightCF,weightCFNormal,tol,tol); }
                catch(const std::exception& e) { std::cerr << "CUT MOMENT FAILURE rank=" << iproc << " level=" << l << " element=" << iel << " stage=surface error=" << e.what() << std::endl; throw; }
              }
              else rule->GetTangentialSurfaceWeightsFromLevelSetData(levelSet,coordX,weightCFInt,weightCF,tol,tol);
            }
            else {
              auto* rule = elementSurfaceRule;
              if(rule->GetGaussQuadraturePointNumber() != nGauss) throw std::runtime_error("Quadratic CutFEM Gauss rule mismatch");
              rule->GetPhysicalSurfaceWeightsFromLevelSetData(levelSet,coordX,weightCFInt,weightCF,tol,tol);
            }
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
        std::vector<double> surfaceForce, surfaceTangent, surfaceShapeGradient, levelSetGradient(dim), levelSetShapeGradient, mapShapeGradient, parentGradient, xi_gf(dim);
        std::vector<std::vector<double>> mapGradient(dim,std::vector<double>(dim));
        const unsigned surfaceTangentSize = dim == 2 ? 1U : 3U;
        for(unsigned ig = 0; ig < femV->GetGaussPointNumber(); ig++) {
          // *** get gauss point weight, test function and test function partial derivatives ***

          if(l > levelC) {
            femXi->Jacobian(coordX, ig, weight, phiXi, phiXi_x);
            std::fill(xi_gf.begin(),xi_gf.end(),0.);
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
          double surfaceMetricSquared = 0.;
          if(cut == 1 && tangentialSurfaceWeights) {
            std::fill(levelSetGradient.begin(),levelSetGradient.end(),0.);
            SurfaceReferenceGradients(femPsi,ig,levelSetShapeGradient);
            for(unsigned i = 0; i < nDofsPsi; ++i) for(unsigned d = 0; d < dim; ++d) levelSetGradient[d] += psi[i]*levelSetShapeGradient[i*dim+d];
            femPsi->GetJacobianMatrix(coordX,ig,weightFake,Jacob,JacI);
            if(l > levelC) SurfaceReferenceGradients(femV,xi_gf,surfaceShapeGradient);
            else SurfaceReferenceGradients(femV,ig,surfaceShapeGradient);
            if(l > levelC) {
              SurfaceReferenceGradients(femXi,ig,mapShapeGradient);
              for(unsigned a = 0; a < dim; ++a) std::fill(mapGradient[a].begin(),mapGradient[a].end(),0.);
              for(unsigned a = 0; a < dim; ++a) for(unsigned b = 0; b < dim; ++b) for(unsigned j = 0; j < nDofsXi; ++j) mapGradient[a][b] += coordXi[a][j]*mapShapeGradient[j*dim+b];
              parentGradient = surfaceShapeGradient;
              std::fill(surfaceShapeGradient.begin(),surfaceShapeGradient.end(),0.);
              for(unsigned i = 0; i < nDofsV; ++i) for(unsigned b = 0; b < dim; ++b) for(unsigned a = 0; a < dim; ++a) surfaceShapeGradient[i*dim+b] += parentGradient[i*dim+a]*mapGradient[a][b];
            }
            surfaceMetricSquared = FELevelSetMoments::TangentialSurfaceFactors(Jacob,levelSetGradient,surfaceShapeGradient,surfaceForce,surfaceTangent);
          }
          if(cut == 1 && !tangentialSurfaceWeights) {
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
              // surfaceWeight = elementSurfaceRule->GetGaussWeightPointer()[ig] * weightCF[ig];
              surfaceWeight = elementSurfaceGaussWeights[ig] * weightCF[ig];
            }
            else {
              surfaceWeight = weight * weightCF[ig] * dsN;
            }
          }

          A += weight * weightCFInt[ig];
          Per += surfaceWeight*(tangentialSurfaceWeights ? surfaceMetricSquared : 1.);

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

              if (captureSurfaceRange) {
                pressureActionYL[I*nDofsV+i] += activeSurfaceRange->PressureJump() * weightedDivPhi * weight * weightCFInt[ig];
              }

              NSV += - weightedDivPhi * (solP1_gss * weightCFInt[ig] + solP2_gss * weightCFExt[ig]) * weight;  // pressure gradient
              NSV += - rho * phiV[i] * g[I] * weight * axisWeight; // gravity term

              Res[I * nDofsV + i] -=  NSV;

              if (cut == 1) {

                double fSigma = 0.0;

                if (tangentialSurfaceWeights && dim == 3 && surfaceWeightCount == 3) {
                  // fSigma = -sigma * elementSurfaceRule->GetGaussWeightPointer()[ig] *
                  fSigma = -sigma * elementSurfaceGaussWeights[ig] *
                    (
                      weightCFNormal[(I+1U)%3U][ig] *
                      surfaceTangent[3U*i+(I+2U)%3U]
                      -
                      weightCFNormal[(I+2U)%3U][ig] *
                      surfaceTangent[3U*i+(I+1U)%3U]
                    );
                }
                else if (tangentialSurfaceWeights) {
                  fSigma = -sigma * surfaceForce[i*dim+I] * surfaceWeight;
                }
                else {
                  for (unsigned d = 0; d < dim; ++d) {
                    fSigma -= sigma * axisWeight * P[I][d] * phiV_x[i*dim+d] * surfaceWeight;
                  }
                }

                if (axisymmetric && I == 0) {
                  fSigma -= sigma * phiV[i] * surfaceWeight;
                }

                const unsigned row = I*nDofsV+i;

                // Original production residual.
                Res[row] += fSigma;

                // Independent diagnostic accumulation.
                if (captureSurfaceRange) {
                  surfaceResidual[row] += fSigma;
                }

                // Keep the existing stabilization calculations below.
              }

              // surface tension stabilization -- rhs contribution
              if(cut == 1) {
                // if(tangentialSurfaceWeights && dim == 3 && surfaceWeightCount == 3) Res[I*nDofsV+i] -= sigma*elementSurfaceRule->GetGaussWeightPointer()[ig]*(weightCFNormal[(I+1U)%3U][ig]*surfaceTangent[3U*i+(I+2U)%3U]-weightCFNormal[(I+2U)%3U][ig]*surfaceTangent[3U*i+(I+1U)%3U]);
                // else if(tangentialSurfaceWeights) Res[I*nDofsV+i] -= sigma*surfaceForce[i*dim+I]*surfaceWeight;
                // else for(unsigned d = 0; d < dim; ++d) Res[I*nDofsV+i] -= sigma*axisWeight*P[I][d]*phiV_x[i*dim+d]*surfaceWeight;

                // if(axisymmetric && I == 0) {
                //   Res[I * nDofsV + i] += -sigma * phiV[i] * surfaceWeight;
                // }

                double stabSF = 0.0;

                if(tangentialSurfaceWeights) {
                  for(unsigned a = 0; a < surfaceTangentSize; ++a) {
                    double increment = 0.;
                    for(unsigned j = 0; j < nDofsV; ++j) increment += (cnew*solV[I][j]-cold*solVOld[I][j])*surfaceTangent[j*surfaceTangentSize+a];
                    stabSF += increment*surfaceTangent[i*surfaceTangentSize+a];
                  }
                }
                else {
                  for (unsigned a = 0; a < dim; ++a) {

                    double gradDeltaU_tg = 0.0;
                    double gradPhi_tg    = 0.0;

                    for (unsigned b = 0; b < dim; ++b) {

                      gradDeltaU_tg += P[a][b] * (cnew * gradSolV_gss[I][b] - cold * gradSolVOld_gss[I][b]);
                      gradPhi_tg += P[a][b] * phiV_x[i * dim + b];

                    }
                    stabSF += gradDeltaU_tg * gradPhi_tg;
                  }
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

          double surfaceTensionCoefficient = cnew*cnew*sigma*dt*surfaceWeight*axisWeight*(tangentialSurfaceWeights ? 1. : 0.);

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
                double tangentialGradDot = 0.;
                if(tangentialSurfaceWeights) for(unsigned a = 0; a < surfaceTangentSize; ++a) tangentialGradDot += surfaceTangent[i*surfaceTangentSize+a]*surfaceTangent[j*surfaceTangentSize+a];
                else tangentialGradDot = gradPhi_dot_gradPhi-normalDotGradPhi[i]*normalDotGradPhi[j];
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

        if (captureSurfaceRange) {

          const unsigned nV = dim * nDofsV;
          const unsigned nP = 2 * nDofsP;

          std::vector<double> pressureOperator(nV * nP, 0.0);

          for (unsigned row = 0; row < nV; ++row) {
            for (unsigned col = 0; col < nP; ++col) {
              pressureOperator[row*nP+col] = -Jac[row*nDofsVP + nV + col];
            }
          }

          activeSurfaceRange->AddElement(sysDof, nV, nP, surfaceResidual, pressureOperator, pressureActionYL);
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

  if (captureSurfaceRange) {
    activeSurfaceRange->FinalizeAssembly();
  }

  std::cout << "Matrix Assembly time        = " << static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC << std::flush << std::endl;

  MPI_Allreduce(MPI_IN_PLACE, &A, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &Per, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


  const double radius = parSphereRadius;
  double radius2 = radius * radius;
  double radius3 = radius2 * radius;
  double per_th = (dim == 2) ? 2 * M_PI * radius : 4. * M_PI * radius2;
  double A_th = (dim == 2) ? M_PI * radius2 : 4. / 3. * M_PI * radius3;
  if(msh->processor_id() == 0) {
    std::cout << " PAR phase volume: " << std::setprecision(16) << A
              << " (relative error " << std::abs(A-A_th)/A_th << ")\n"
              << " PAR interface measure: " << Per
              << " (relative error " << std::abs(Per-per_th)/per_th << ")" << std::endl;
  }

}

