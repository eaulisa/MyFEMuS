
#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
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

const RungeKutta::VelKind velocityType = RungeKutta::VelKind::Translation;
//const RungeKutta::VelKind velocityType = RungeKutta::VelKind::Rotation;
//const RungeKutta::VelKind velocityType = RungeKutta::VelKind::Vortex;

#include "../includeLS/Utils.hpp"


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
  const unsigned numberOfUniformLevels = 1u;
  const unsigned numberOfSelectiveLevels = 3u;

  std::string meshName = "./input/tri.neu";

  // Load coarse mesh and build uniform refinement levels
  mlMsh0.ReadCoarseMesh(meshName.c_str(), "seventh", scalingFactor);
  mlMsh0.RefineMesh(numberOfUniformLevels, numberOfUniformLevels, nullptr);

  unsigned dim = mlMsh0.GetDimension();

  // Parameters for selective AMR (ball centered at xc with radius r)
  const double r = 0.125;
  std::vector<double> xc = {0.0, (velocityType == RungeKutta::VelKind::Translation) ? 0.751 : 0.25};

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

  mlSol0.AddSolution(psiName.c_str(), LAGRANGE, SECOND);
  //for(unsigned d = 0; d < dim; d++) mlSol0.AddSolution(dPsiName[d].c_str(), LAGRANGE, SECOND, 2);
  mlSol0.AddSolution("Gamma", LAGRANGE, SECOND);

  std::vector<std::string> vName = {"U", "V", "W"};
  vName.resize(dim);

  for(unsigned d = 0; d < dim; d++) mlSol0.AddSolution(vName[d].c_str(), LAGRANGE, SECOND, 2);
  mlSol0.Initialize("All");

  InitSol(mlSol0, vName, 0, 1.);


  //for(unsigned d = 0; d < dim; d++) mlSol0.Initialize(vName[d].c_str(), Initvel[d]);

  unsigned sigmoidType = 1;
  double eps = 0.25; //(dim == 2) ? 1. / pow(2, std::max(levelN - 7u, 1u))
  Mollifier m = Mollifier(eps, sigmoidType);

  std::vector<double> xc_1 = xc;
  std::vector<double> xc_2 = (dim == 2) ? std::vector<double>({xc[0]-0.35, xc[1] + 0.5}) : std::vector<double>({xc_1[0], xc_1[1], xc_1[2] + 0.5});
  std::vector<double> xc_3 = (dim == 2) ? std::vector<double>({xc[0]+0.35, xc[1] + 0.25}) : std::vector<double>({xc_2[0], xc_2[1], xc_2[2] + 0.5});
  double r_0 = r;
  Circle c1(xc_1,r_0);
  Circle c2(xc_2, r_0);
  Circle c3(xc_3, r_0);
  // Shape* shape = &c1;

  std::vector<Shape*> shape = {&c1};
  std::vector<double> timeOffset = {1.};

  Inflow inflow(InflowVelocity, shape, timeOffset, m);

  std::vector<std::vector<std::vector<double>>> inflow_markers0(shape.size());
  
  for (unsigned s = 0; s < shape.size(); s++) {
    inflow_markers0[s].resize(dim);
    if (iproc == 0) {
      shape[s]->GetMarkers(inflow_markers0[s], 1000);
    }
  }

  Boundary* inflow_bd = &inflow;

  PsiBall psi2D(xc, r, m);
  InitLevelSet(mlSol0, "Psi", psi2D);

  //GetSolutionGradient(mlSol0, psiName, dPsiName);

  // Export solution to VTK (selected levels)
  std::vector<std::string> variablesToBePrinted = {"All"};
  VTKWriter vtkIO(&mlSol0);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);

  MultiLevelMesh mlMsh1;
  MultiLevelSolution mlSol1;

  MultiLevelMesh *mlmsh0 = &mlMsh0;
  MultiLevelMesh *mlmsh1 = &mlMsh1;

  MultiLevelSolution *mlsol0 = &mlSol0;
  MultiLevelSolution *mlsol1 = &mlSol1;

  // Initialize markers object
  LevelSetMarkers markers("Psi", dim);

  // Load coarse mesh and build uniform refinement levels
  mlmsh1->ReadCoarseMesh(meshName.c_str(), "seventh", scalingFactor);
  mlmsh1->RefineMesh(numberOfUniformLevels, numberOfUniformLevels, nullptr);

  // RungeKutta::VelKind velocityType = RungeKutta::VelKind::Rotation;
  //RungeKutta::VelKind velocityType = RungeKutta::VelKind::Vortex;
  // RungeKutta::VelKind velocityType = RungeKutta::VelKind::Rotation;

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
    InitSol(*mlsol0, vName, time, period);// TODO with Navier-Stokes solve

    bbox.SetMesh(mlmsh0->GetLevel(0));

    RungeKutta rk(time, dt, period, velocityType);

    //RungeKuttaN rk(time, dt, mlmsh0, bbox);

    unsigned nLevels = numberOfUniformLevels + numberOfSelectiveLevels;
    std::vector<MyVector<double>> X0;
    MyVector<int> X0Iel;
    // GetCutElementPoints(*mlsol0, "Psi", X0, X0Iel);

    std::vector<std::vector<double>> inflow_markers(dim);
    inflow_bd->updateMarkers(inflow_markers0, inflow_markers, time-dt, period, dt);

    std::vector<MyVector<double>> IX(dim);
    for (unsigned k = 0; k < dim; ++k) {
      IX[k].buildFromLocal(inflow_markers[k]);
    }

    {
      LevelMarkers l0;
      const unsigned bboxLevels = nLevels - bbox.GetLevel();

      std::vector<LevelMarkers> lX(bboxLevels);

      bbox.GetInverseMappingOnCoarseLevel(IX, l0, lX[0]);

      const std::vector<bool> &isInsideDomain = l0.GetPointInsideDomain();

      for (unsigned d = 0; d < dim; d++)
        inflow_markers[d].clear();

      for (unsigned d = 0; d < dim; d++) {
        unsigned offset = IX[d].begin();
        for (unsigned i = IX[d].begin(); i < IX[d].end(); ++i) {
          if (!isInsideDomain[i - offset]) {
              inflow_markers[d].push_back(IX[d][i]);
          }
        }
      }
    }
    
    markers.GetCutElementPoints(*mlsol0, X0, X0Iel, inflow_markers);
    
    if (t % 10 == 0) {
      Reinit reinit("Psi", *mlsol0, m);
    
      reinit.farFieldReinit(X0);
      reinit.interfaceFieldReinit(bbox);
      reinit.updateSolution();
    }

    if (t == 1)
      WritePointsVTK("./output/points.0.vtk", X0);

    RungeKutta4(X0, *mlsol0, bbox, vName, dt);
    //rk.rkForward(X0);

    if (t % 10 == 0)
      WritePointsVTK("./output/points." + std::to_string(t / 10) + ".vtk", X0);

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

    // TestProjections(l0, lX);

    mlsol1->Build(mlmsh1);
    mlsol1->AddSolution("Psi", LAGRANGE, SECOND);

    for(unsigned d = 0; d < dim; d++) mlsol1->AddSolution(vName[d].c_str(), LAGRANGE, SECOND, 2);

    mlsol1->Initialize("All");
    //InitSol(*mlsol1, vName, time, period);
    //for(unsigned d = 0; d < dim; d++) mlsol1->Initialize(vName[d].c_str(), Initvel[d]);


    ProjectSolution(*mlsol0, *mlsol1, bbox, rk, {"Psi"}, vName, *inflow_bd, -dt, time, period);
    ProjectSolution(*mlsol0, *mlsol1, bbox, rk, vName);

    // Export solution to VTK (selected levels)

    std::swap(mlsol0, mlsol1);
    std::swap(mlmsh0, mlmsh1);
    mlsol1->clear();
    mlmsh1->resize(numberOfUniformLevels);

    VTKWriter vtkIO1(mlsol0);
    if (t % 10 == 0)
      vtkIO1.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted,
                   t / 10);

    double area = ComputeArea(*mlsol0, psiName);

    if (iproc == 0) {

      std::ofstream out("area.dat", std::ios::app);

      if (!out) {
        throw std::runtime_error(
            "computeArea: cannot open output file area.dat");
      }

      out << std::setprecision(16)
          <<std::setw(20) << time
          << std::setw(25) << area
          << "\n";

      out.close();
    }
  }

  if (nprocs == 1)
    ProfilerStop();
  return 0;
}
