#include <iostream>
#include <vector>
#include <cmath>
#include <memory>

#include <gperftools/profiler.h>

#include "Mesh.hpp"
#include "Field.hpp"
#include "FemProjection.hpp"   // FemProjection, Quad9Projection, Tri7Projection
#include "RefineMesh.hpp"
#include "VtkOutput.hpp"
#include "VtuOutput.hpp"
#include "Mollifier.hpp"
#include "Psi.hpp"
#include "RungeKutta.hpp"
#include "PointLocator.hpp"
#include "Errors.hpp"
#include "MeshSeed.hpp"


#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"

#include "SimConfig.hpp"

#include "Reinit.hpp"



int main(int argc, char** args) {


  ProfilerStart("profiling.prof");


  SimConfig cfg = parseArgs(argc, args);

  std::cout << "dim=" << cfg.dim
            << " vel=" << (int)cfg.velocityType
            << " meshSeed=" << (int)cfg.meshSeed
            << " shape=" << (int)cfg.meshShape
            << " period=" << cfg.period
            << " steps=" << cfg.nSteps
            << " uniformRefinementLevel=" <<cfg.uniformRefinementLevel
            << "\n";


  std::array<std::unique_ptr<FemProjection>, 6> elProj;
  elProj[0].reset(new Hex27Projection());
  elProj[1].reset(new Tet15Projection());
  elProj[2].reset(new Wedge21Projection());
  elProj[3].reset(new Quad9Projection());
  elProj[4].reset(new Tri7Projection());
  elProj[5].reset(new Line3Projection());

  RungeKutta::VelKind velocityType = cfg.velocityType;
  unsigned uniformRefinementLevel = cfg.uniformRefinementLevel;
  std::vector<double> xc = cfg.xc;
  double r = cfg.r;
  const double period = cfg.period;
  const unsigned print_step = cfg.print_step;
  const unsigned reinit_step = cfg.reinit_step;
  const unsigned levelNstart = cfg.levelNstart;
  unsigned delta_depth = cfg.delta_depth;
  unsigned nSteps = cfg.nSteps;
  const bool advect_markers = cfg.advect_markers;
  const bool reinit_adaptive = cfg.reinit_adaptive;
  const double reinit_tol = cfg.reinit_tol;


  std::vector<std::tuple<double, double, double >> Er;
  Er.reserve(delta_depth);

  std::vector<clock_t> Time;
  Time.reserve(delta_depth);

  unsigned reinit_cnt[delta_depth];


  //BEGIN LEVEL LOOP
  for (unsigned levelN = levelNstart; levelN < levelNstart + delta_depth; levelN++) {

    clock_t start_time =  clock();

    std::vector<std::vector<std::vector<double>>> X0(levelN + 1);
    std::vector<std::vector<unsigned>> elType0(levelN + 1);
    std::vector<std::vector<std::vector<unsigned>>> elTplgy0(levelN + 1);
    std::vector<std::vector<unsigned>> elLevel0(levelN + 1);
    std::vector<std::vector<int>> AMR0(levelN + 1);
    std::vector<std::vector<unsigned>> elFather0(levelN + 1);
    std::vector<std::vector<std::vector<unsigned>>> elChildern0(levelN + 1);
    std::vector<std::vector<std::vector<unsigned>>> nodeElements0(levelN + 1);
    std::vector<std::vector<std::vector<unsigned>>> neighbors0(levelN + 1);


    std::vector<std::vector<std::vector<double>>> X1(levelN + 1);
    std::vector<std::vector<unsigned>> elType1(levelN + 1);
    std::vector<std::vector<std::vector<unsigned>>> elTplgy1(levelN + 1);
    std::vector<std::vector<unsigned>> elLevel1(levelN + 1);
    std::vector<std::vector<int>> AMR1(levelN + 1);
    std::vector<std::vector<unsigned>> elFather1(levelN + 1);
    std::vector<std::vector<std::vector<unsigned>>> elChildern1(levelN + 1);
    std::vector<std::vector<std::vector<unsigned>>> nodeElements1(levelN + 1);
    std::vector<std::vector<std::vector<unsigned>>> neighbors1(levelN + 1);


    std::vector<std::vector<std::vector<double>>> X2(levelN + 1);
    std::vector<std::vector<unsigned>> elType2(levelN + 1);
    std::vector<std::vector<std::vector<unsigned>>> elTplgy2(levelN + 1);
    std::vector<std::vector<unsigned>> elLevel2(levelN + 1);
    std::vector<std::vector<int>> AMR2(levelN + 1);
    std::vector<std::vector<unsigned>> elFather2(levelN + 1);
    std::vector<std::vector<std::vector<unsigned>>> elChildern2(levelN + 1);
    std::vector<std::vector<std::vector<unsigned>>> nodeElements2(levelN + 1);
    std::vector<std::vector<std::vector<unsigned>>> neighbors2(levelN + 1);




    std::vector<Mesh> mesh0, mesh1, mesh2;
    std::vector<Field> field0, field1, field2;
    mesh0.reserve(levelN + 1);
    field0.reserve(levelN + 1);

    mesh1.reserve(levelN + 1);
    field1.reserve(levelN + 1);

    mesh2.reserve(levelN + 1);
    field2.reserve(levelN + 1);

    for (unsigned l = 0; l <= levelN; l++) {
      mesh0.emplace_back(elTplgy0[l], elType0[l], elLevel0[l], X0[l], AMR0[l], elFather0[l], elChildern0[l], nodeElements0[l], neighbors0[l]);
      field0.emplace_back(mesh0[l]);

      mesh1.emplace_back(elTplgy1[l], elType1[l], elLevel1[l], X1[l], AMR1[l], elFather1[l], elChildern1[l], nodeElements1[l], neighbors1[l]);
      field1.emplace_back(mesh1[l]);

      mesh2.emplace_back(elTplgy2[l], elType2[l], elLevel2[l], X2[l], AMR2[l], elFather2[l], elChildern2[l], nodeElements2[l], neighbors2[l]);
      field2.emplace_back(mesh2[l]);
    }


    //BEGIN MESH AND FIELD INITIALIZATION
    MeshSeed seed = MeshSeedFactory::make(cfg.meshSeed, cfg.meshShape, cfg.shift);

//    auto seed = MeshSeedFactory::make(meshSeed, meshShape, shift);

    elLevel0[0] = seed.elLevel;
    elType0[0]  = seed.elType;
    elTplgy0[0] = seed.elTplgy;
    X0[0]       = seed.X;
    unsigned dim = X0[0].size();

    std::vector<unsigned> candidateIndices(mesh0[0].numNodes());
    std::iota(candidateIndices.begin(), candidateIndices.end(), 0u);
    dedupNodesHash(mesh0[0], candidateIndices);

    mesh0[0].resetAllFathersToNoFather();
    mesh0[0].buildNodeToElementAdjacency();
    mesh0[0].buildFaceNeighborsFromNodeToElement();
    field0[0].clear();

    std::string filename = "./output/refined_mesh" + std::to_string(dim) + "D_level" + std::to_string(levelN) + ".";

    double eps = (dim == 2) ? 1. / pow(2, std::max(levelN - 7u, 1u)) : 1. / pow(2, std::max(levelN - 4u, 1u));

    //init multilevel mesh
    unsigned neighMode = 3; // 0=no-ring, 1=vertices, 2=faces, 3=hybrid
    for (unsigned l = 1; l <= levelN; l++) {
      if (l <= uniformRefinementLevel) {
        mesh0[l - 1].setUniformRefinement();
      }
      else {
        mesh0[l - 1].setRefinementFromBallLevelSetCrossing_OneRing(xc, r, neighMode);
        mesh0[l - 1].adjustAMRForOneLevelDiscontinuity();
      }
      refineAndProjectMesh(elProj, mesh0[l - 1], mesh0[l]);
      field0[l].clear();
    }

    // init field at the top level
    const unsigned topLevel = levelN;
    PsiBall psi2D(xc, r, eps);
    field0[topLevel].addField("Psi", Field::Location::Nodal, 1.);
    const unsigned psiId0 = field0[topLevel].id("Psi");
    auto& Psi0 = field0[topLevel].getById(psiId0);
    for (std::size_t k = 0; k < Psi0.size(); ++k) {
      const std::vector<double> x = field0[topLevel].dofCoordById(psiId0, k); // nodal => mesh node coord
      Psi0[k] = psi2D(x);
    }

    std::cout << "Iteration = " << 0 << " Number of Points = " << mesh0[topLevel].X()[0].size() << std::endl;
    writeMeshFieldVTU(filename + std::to_string(0) + ".vtu", field0[topLevel]);

    for (unsigned l = 0; l <= levelN; l++) field2[l] = field0[l];
    //END MESH AND FIELD INITIALIZATION


    //BEGIN TIME LOOP
    std::vector<std::vector<double>> Xp;
    PointLocator pl = PointLocator(mesh0[0], .1);
    std::vector<PointLocatorResult> out, in;

    const unsigned nIter = 320;
    const double dt = period / static_cast<double>(nIter);
    reinit_cnt[levelN-levelNstart] = 0;

    mesh1[0] = mesh0[0];

    // BEGIN MARKERS PLACEMENT
    std::vector<std::vector<double>> markers;
    Reinit reinit(elProj, field0, psiId0, psi2D._m);
    reinit.computeMarkersAdvection(markers);
    // END MARKERS PLACEMENT
    
    for (unsigned k = 1; k <= nSteps; ++k) {
      const double time = k * dt;
      bool reinit_flag = (k % reinit_step == 0);

      out.clear();

      if (advect_markers) { 
        RungeKutta::rkForward(markers, time, dt, velocityType);
        pl.locateAll(out, markers);
      } else {
        field0[topLevel].extractInterfaceVerticesAndCentersByName("Psi", Xp, topLevel, levelN + 1);
        RungeKutta::rkForward(Xp, time, dt, velocityType);
        pl.locateAll(out, Xp);
      }

      mesh1[0].setRefinementFromLocatedPoints_OneRing(out, neighMode);
      mesh1[0].adjustAMRForOneLevelDiscontinuity();
      field1[0].clear();

      for (unsigned l = 1; l <= levelN; l++) {
        mesh1[l].clearAllData();
        refineAndProjectMesh(elProj, mesh1[l - 1], mesh1[l]);

        std::swap(in, out);
        mesh1[l - 1].projectPointLocatorResultsToNextLevel(mesh1[l], in, out);
        mesh1[l].setRefinementFromLocatedPoints_OneRing(out, neighMode);
        mesh1[l].adjustAMRForOneLevelDiscontinuity();

        field1[l].clear();
      }

      //move new multilevel mesh1 nodes backward to get the new topLevel field1 from the old filed0
      Xp = mesh1[topLevel].X();
      RungeKutta::rkBackward(Xp, time, dt, velocityType);
      out.clear();
      pl.locateAll(out, Xp);
      for (unsigned l = 1; l <= levelN; l++) {
        std::swap(in, out);
        mesh0[l - 1].projectPointLocatorResultsToNextLevel(mesh0[l], in, out);
      }

      field1[topLevel].addField("Psi", Field::Location::Nodal, 1.);
      const unsigned psiId1 = field1[topLevel].id("Psi");
      auto& Psi1 = field1[topLevel].getById(psiId1);
      const unsigned psiId0 = field0[topLevel].id("Psi");
      field0[topLevel].evalNodalAtLocatedPointsById(psiId0, out, elProj, Psi1, -1.0);

      //swap for the next iteration
      for (unsigned l = 0; l <= levelN; l++) swap(field0[l], field1[l]);

      // BEGIN PSI reinit
      {
        Reinit reinit(elProj, field0, psiId0, psi2D._m);

        // markers for advection and/or field distorsion evaluation
        if (advect_markers || reinit_adaptive)
          reinit.computeMarkersAdvection(markers);

        // reinitialization step
        if (reinit_flag) {
          reinit.reinitializeSignedDistance();
          reinit_cnt[levelN-levelNstart]++;
        } else if (reinit_adaptive) {
          double error = reinit.fieldDistortion(markers, false);
          if (error > reinit_tol) {
            reinit.reinitializeSignedDistance();
            reinit_cnt[levelN-levelNstart]++;
          }
        }
      }
      // END PSI reinit

      if (k % print_step == 0) std::cout << "\x1b[1A" << "\x1b[2K";   // cursor up 1, and erase entire line
      std::cout << "\x1b[1A" << "\x1b[2K"                        // cursor up 1, and erase entire line
                << "\r"                                          // return to column 1
                << std::flush;

      std::cout << "Iteration = " << k << " Number of Points = " << mesh1[topLevel].X()[0].size() << std::endl;
      if (k % print_step == 0) writeMeshFieldVTU(filename + std::to_string(k / print_step) + ".vtu", field0[topLevel]);

      
    }
    //END TIME LOOP

    //BEGIN MERGE INITIAL AND FINAL SOLUTIONS
    {
      for (unsigned l = 0; l <= levelN; l++) {
        swap(field1[l], field2[l]);
      }

      std::vector<std::vector<double>> Yp;
      field0[topLevel].extractInterfaceVerticesAndCentersByName("Psi", Xp, topLevel, levelN + 1);
      field1[topLevel].extractInterfaceVerticesAndCentersByName("Psi", Yp, topLevel, levelN + 1);

      if (Yp.size() != Xp.size()) {
        throw std::runtime_error("Xp and Yp have different spatial dimensions");
      }

      for (std::size_t d = 0; d < Xp.size(); ++d) {
        Xp[d].insert(Xp[d].end(), Yp[d].begin(), Yp[d].end());
      }
      pl.locateAll(out, Xp);


      mesh2[0].setRefinementFromLocatedPoints_OneRing(out, neighMode);
      mesh2[0].adjustAMRForOneLevelDiscontinuity();
      field2[0].clear();
      for (unsigned l = 1; l <= levelN; l++) {
        mesh2[l].clearAllData();
        refineAndProjectMesh(elProj, mesh2[l - 1], mesh2[l]);

        std::swap(in, out);
        mesh2[l - 1].projectPointLocatorResultsToNextLevel(mesh2[l], in, out);
        mesh2[l].setRefinementFromLocatedPoints_OneRing(out, neighMode);
        mesh2[l].adjustAMRForOneLevelDiscontinuity();
        field2[l].clear();
      }

      const auto& Xp = mesh2[topLevel].X();
      out.clear();
      pl.locateAll(out, Xp);

      for (unsigned l = 1; l <= levelN; l++) {
        std::swap(in, out);
        mesh0[l - 1].projectPointLocatorResultsToNextLevel(mesh0[l], in, out);
      }

      field2[topLevel].addField("PsiE", Field::Location::Nodal, 1.);
      const unsigned psiEId2 = field2[topLevel].id("PsiE");
      auto& PsiE = field2[topLevel].getById(psiEId2);

      const unsigned psiId0 = field0[topLevel].id("Psi");
      field0[topLevel].evalNodalAtLocatedPointsById(psiId0, out, elProj, PsiE, -1.0);

      out.clear();
      pl.locateAll(out, Xp);
      for (unsigned j = 1; j <= levelN; j++) {
        std::swap(in, out);
        mesh1[j - 1].projectPointLocatorResultsToNextLevel(mesh1[j], in, out);
      }

      field2[topLevel].addField("PsiS", Field::Location::Nodal, 1.);
      const unsigned psiSId2 = field2[topLevel].id("PsiS");
      auto& PsiS = field2[topLevel].getById(psiSId2);

      const unsigned psiId1 = field1[topLevel].id("Psi");
      field1[topLevel].evalNodalAtLocatedPointsById(psiId1, out, elProj, PsiS, -1.0);

      writeMeshFieldVTU(filename + "StartAndEnd.0" + ".vtu", field2[topLevel]);
      writeMeshFieldVTU(filename + "StartAndEnd.1" + ".vtu", field2[topLevel]);

      Er.push_back(computeMassAndGeometricError("PsiS", "PsiE", field2[topLevel], elProj));

      Time.push_back(clock() - start_time);
    }
    //END MERGE INITIAL AND FINAL SOLUTIONS
  }
  //END LEVEL LOOP


  std::cout << "Max_Depth\tMass_Error\tGeom_Error\tScaled_Geom_Error\tCompt_Time(s)\tReinit_Counter" << std::endl;
  for (unsigned i = 0; i < delta_depth; i++) {
    std::cout << levelNstart + i << "\t\t"
              << std::get<0>(Er[i]) << "\t"
              << std::get<1>(Er[i]) << "\t"
              << std::get<2>(Er[i]) << "\t\t"
              << static_cast<double>(Time[i]) / CLOCKS_PER_SEC << "\t\t"
              << reinit_cnt[i] << std::endl;
    if (i + 1 < delta_depth) {
      std::cout << "conv.\t\t"
                << log(std::get<0>(Er[i]) / std::get<0>(Er[i + 1])) / log(2.) << "\t\t"
                << log(std::get<1>(Er[i]) / std::get<1>(Er[i + 1])) / log(2.) << "\t\t"
                << log(std::get<2>(Er[i]) / std::get<2>(Er[i + 1])) / log(2.) << "\t\t"
                << log((double)Time[i + 1] / (double)Time[i]) / log(2.) << std::endl;
    }
  }
  if (delta_depth > 1) {
    std::cout << "\naver. conv. \t"
              << log(std::get<0>(Er[0]) / std::get<0>(Er[delta_depth - 1])) / ((delta_depth - 1) * log(2.)) << "\t\t"
              << log(std::get<1>(Er[0]) / std::get<1>(Er[delta_depth - 1])) / ((delta_depth - 1) * log(2.)) << "\t\t"
              << log(std::get<2>(Er[0]) / std::get<2>(Er[delta_depth - 1])) / ((delta_depth - 1) * log(2.)) << "\t\t"
              << log((double)Time[delta_depth - 1] / (double)Time[0]) / ((delta_depth - 1) * log(2.)) << std::endl;
  }
  return 0;


  ProfilerStop();

  return 1;
}
