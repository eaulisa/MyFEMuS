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

int main() {


  ProfilerStart("profiling.prof");

  std::array<std::unique_ptr<FemProjection>, 6> elProj;
  elProj[0].reset(new Hex27Projection());
  elProj[1].reset(new Tet15Projection());
  elProj[2].reset(new Wedge21Projection());
  elProj[3].reset(new Quad9Projection());
  elProj[4].reset(new Tri7Projection());
  elProj[5].reset(new Line3Projection());

  //std::cout<<elProj[0]->GetProjection()<<std::endl;
  //std::cout << elProj[1]->GetProjection() << std::endl;
  //std::cout << elProj[2]->GetProjection() << std::endl;


  const unsigned levelNstart = 6;
  unsigned delta_depth = 4;

  std::vector<std::pair<double, double>> Er;
  Er.reserve(delta_depth);

  std::vector<clock_t> Time;
  Time.reserve(delta_depth);


  //MeshSeedFactory::Type meshSeed = MeshSeedFactory::Type::SquareTri7;
  //MeshSeedFactory::Type meshSeed = MeshSeedFactory::Type::SquareQuad9;

  //std::vector<double> xc = {0, 0.25};
  //double r = 0.15;
  // const double period = 8.;

  MeshSeedFactory::Type meshSeed = MeshSeedFactory::Type::CubeHex27;

  std::vector<double> xc = {0, 0, 0.25};
  double r = 0.15;
  const double period = 4.;

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


    /* elLevel0[0] = elLevel1[0] = {0, 0};
     elType0[0] = elType1[0] = {3, 4};
     elTplgy0[0] = elTplgy1[0] = {
       {0, 1, 2, 3, 5, 6, 7, 8, 11},
       {1, 4, 2, 9, 10, 6, 12}
     };
     X0[0] = X1[0] = {
       { -0.5,  0.5,  0.5, -0.5,  1.5,  0.0,  0.5,  0.0, -0.5,  1.0,  1.0,  0.0,  0.8333333333333333 },
       { -0.5, -0.5,  0.5,  0.5, -0.5, -0.5,  0.0,  0.5,  0.0, -0.5,  0.0,  0.0, -0.16666666666666669 }
     };
    */


    auto seed = MeshSeedFactory::make(meshSeed);

    elLevel0[0] = seed.elLevel;
    elType0[0]  = seed.elType;
    elTplgy0[0] = seed.elTplgy;
    X0[0]       = seed.X;
    unsigned dim = X0[0].size();


    // //square box
    // elLevel0[0] = {0};
    // elType0[0] = {3};
    // elTplgy0[0] = {
    //   {0, 1, 2, 3, 4, 5, 6, 7, 8}
    // };
    //
    // X0[0] = {
    //   { -0.5,  0.5,  0.5, -0.5,  0.0,  0.5,  0.0, -0.5,  0.0},
    //   { -0.5, -0.5,  0.5,  0.5, -0.5,  0.0,  0.5,  0.0,  0.0}
    // };
    //
    //
    //
    // //triangle box
    // elLevel0[0] = {0, 0};
    // elType0[0] = {4, 4};
    // elTplgy0[0] = {
    //   {0, 1, 3, 4, 8, 7, 9},
    //   {1, 2, 3, 5, 6, 8, 10}
    // };
    //
    // X0[0] = {
    //   { -0.5,  0.5,  0.5, -0.5,  0.0,  0.5,  0.0, -0.5,  0.0,  -1. / 6., 1. / 6. },
    //   { -0.5, -0.5,  0.5,  0.5, -0.5,  0.0,  0.5,  0.0,  0.0,  -1. / 6., 1. / 6. }
    // };



// X0[0] = X1[0] = {
//    {0, 1, 1, 0., 2., 0.5, 1., 0.5, 0, 1.5, 1.5, 0.5, 1 + 1. / 3.},
//    {0, 0, 1, 1, 0, 0, 0.5, 1, 0.5, 0, 0.5, 0.5, 1. / 3}
//  };

    std::vector<unsigned> candidateIndices(mesh0[0].numNodes());
    std::iota(candidateIndices.begin(), candidateIndices.end(), 0u);
    dedupNodesHash(mesh0[0], candidateIndices);

    mesh0[0].resetAllFathersToNoFather();
    mesh0[0].buildNodeToElementAdjacency();
    mesh0[0].buildFaceNeighborsFromNodeToElement();
    field0[0].clear();

    std::string filename = "./output/refined_mesh2D_level" + std::to_string(levelN) + ".";

    double eps = (dim == 2) ? 1. / pow(2, std::max(levelN - 7u, 1u)) : 1. / pow(2, std::max(levelN - 4u, 1u));

    //init multilevel mesh
    unsigned neighMode = 3; // 0=no-ring, 1=vertices, 2=faces, 3=hybrid
    for (unsigned l = 1; l <= levelN; l++) {
      if (l == 1) {
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

    mesh1[0] = mesh0[0];
    for (unsigned k = 1; k <= nIter; ++k) {
      const double time = k * dt;

      //move the old topLevel mesh0 nodes forward and use them to build the new multilevel mesh1 recursively
      field0[topLevel].extractInterfaceVerticesAndCentersByName("Psi", Xp, topLevel, levelN + 1);
      RungeKutta::rkForward(Xp, time, dt, RungeKutta::VelKind::Vortex);

      out.clear();
      pl.locateAll(out, Xp);
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
      RungeKutta::rkBackward(Xp, time, dt, RungeKutta::VelKind::Vortex);
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

      std::cout << "\x1b[1A" << "\x1b[2K"   // cursor up 1, and erase entire line
                << "\x1b[1A" << "\x1b[2K"   // cursor up 1, and erase entire line
                << "\r"        // return to column 1
                << std::flush;
      std::cout << "Iteration = " << k << " Number of Points = " << mesh1[topLevel].X()[0].size() << std::endl;
      writeMeshFieldVTU(filename + std::to_string(k) + ".vtu", field1[topLevel]);

      //swap for the next iteration
      for (unsigned l = 0; l <= levelN; l++) swap(field0[l], field1[l]);
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



  std::cout << "Max_Depth\tMass_Error\tGeometric_Error\tCompt_Time(s)" << std::endl;
  for (unsigned i = 0; i < delta_depth; i++) {
    std::cout << levelNstart + i << "\t\t" << Er[i].first << "\t" << Er[i].second << "\t"
              /*      */ << static_cast<double>(Time[i]) / CLOCKS_PER_SEC << std::endl;
    if (i + 1 < delta_depth) {
      std::cout << "conv.\t\t" << log(Er[i].first / Er[i + 1].first) / log(2.) << "\t\t"
                /*      */ << log(Er[i].second / Er[i + 1].second) / log(2.) << "\t\t"
                /*      */ << log((double)Time[i + 1] / (double)Time[i]) / log(2.) << std::endl;
    }
  }
  if (delta_depth > 1) {
    std::cout << "\naver. conv. \t" << log(Er[0].first / Er[delta_depth - 1].first) / ((delta_depth - 1) * log(2.))
              << "\t\t" << log(Er[0].second / Er[delta_depth - 1].second) / ((delta_depth - 1) * log(2.))
              << "\t\t" << log((double)Time[delta_depth - 1] / (double)Time[0]) / ((delta_depth - 1) * log(2.)) << std::endl;
  }
  return 0;


  ProfilerStop();

  return 1;
}

// mesh0[0].clearAllData();
//
//
// elType0[0] = {0, 2, 1};
// elLevel0[0] = {0, 0, 0};
// elTplgy0[0] = {
//   {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26},
//   {27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47},
//   {48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62}
// };
//
// X0[0] = {{
//     //hex
//     0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
//     0.5, 1.0, 0.5, 0.0, 0.5, 1.0, 0.5, 0.0, 0.0, 1.0, 1.0, 0.0,
//     0.5, 1.0, 0.5, 0.0, 0.5, 0.5,
//     0.5,
//     //wedge
//     1.0, 2.0, 1.0, 1.0, 2.0, 1.0,
//     1.5, 1.5, 1.0, 1.5, 1.5, 1.0, 1.0, 2.0, 1.0,
//     1.5, 1.5, 1.0, 1.3333333333333333, 1.3333333333333333,
//     1.3333333333333333,
//     //tet
//     1.0, 2.0, 1.0, 1.0,
//     1.5, 1.5, 1.0, 1.0, 1.5, 1.0,
//     1.3333333333333333, 1.3333333333333333,
//     1.3333333333333333, 1.0,
//     1.25
//   }, {
//     //hex
//     0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0,
//     0.0, 0.5, 1.0, 0.5, 0.0, 0.5, 1.0, 0.5, 0.0, 0.0, 1.0, 1.0,
//     0.0, 0.5, 1.0, 0.5, 0.5, 0.5,
//     0.5,
//     //wedge
//     0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
//     0.0, 0.5, 0.5, 0.0, 0.5, 0.5, 0.0, 0.0, 1.0,
//     0.0, 0.5, 0.5, 1.0 / 3.0, 1.0 / 3.0,
//     1.0 / 3.0,
//     //tet
//     0.0, 0.0, 1.0, 0.0,
//     0.0, 0.5, 0.5, 0.0, 0.0, 0.5,
//     0.3333333333333333, 0.0,
//     0.3333333333333333, 0.3333333333333333,
//     0.25
//   }, {
//     //hex
//     0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,
//     0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5,
//     0.5, 0.5, 0.5, 0.5, 0.0, 1.0,
//     0.5,
//     //wedge
//     0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
//     0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.5, 0.5, 0.5,
//     0.5, 0.5, 0.5, 0.0, 1.0,
//     0.5,
//     //tet
//     1.0, 1.0, 1.0, 2.0,
//     1.0, 1.0, 1.0, 1.5, 1.5, 1.5,
//     1.0, 1.3333333333333333,
//     1.3333333333333333, 1.3333333333333333,
//     1.25
//   }
// };
//
// candidateIndices.resize(mesh0[0].numNodes());
// std::iota(candidateIndices.begin(), candidateIndices.end(), 0u);
// dedupNodesHash(mesh0[0], candidateIndices);
//
// std::cout << "Number of nodes after deduplication = " << mesh0[0].numNodes() << std::endl;
// mesh0[0].resetAllFathersToNoFather();
// mesh0[0].buildNodeToElementAdjacency();
// mesh0[0].buildFaceNeighborsFromNodeToElement();
//
// //field0[0].rebindMeshAndResize(mesh0[0]);
//
// PsiBall psi3D(std::vector<double> {1., 0, 1.}, 0.125, 0.001);
// auto psiId = field0[0].id("Psi");
// auto Psi = field0[0].getById(psiId);
// for (std::size_t k = 0; k < Psi.size(); ++k) {
//   const std::vector<double> x = field0[0].dofCoordById(psiId, k); // nodal => mesh node coord
//   Psi[k] = psi3D(x);
// }
//
// filename = "./output/refined_mesh3D.";
// writeMeshFieldVTK(filename + "0.vtk", field0[0]);
//
//
// for (unsigned l = 1; l < levelN; l++) {
//   mesh0[l - 1].setRefinementFromBallLevelSetCrossing_OneRing({1., 0., 1.}, 0.125, neighMode);
//   mesh0[l - 1].adjustAMRForOneLevelDiscontinuity();
//   refineAndProjectMesh(elProj, mesh0[l - 1], mesh0[l]);
//
//   //field0[l].rebindMeshAndResize(mesh0[l]);
//   const unsigned psiId = field0[l].id("Psi");
//   auto& Psi = field0[l].getById(psiId);
//   for (std::size_t k = 0; k < Psi.size(); ++k) {
//     const std::vector<double> x = field0[l].dofCoordById(psiId, k); // nodal => mesh node coord
//     Psi[k] = psi3D(x);
//   }
//   writeMeshFieldVTK(filename + std::to_string(l) + ".vtk", field0[l]);
// }
//
// return 0;
//}
























/*

#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "adept.h"
#include "MeshRefinement.hpp"

#include "projection.hpp"

using namespace femus;

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0.;
  return dirichlet;
}


bool SetRefinementFlag(const std::vector < double >& x, const int& elemgroupnumber, const int& level) {

  bool refine = 0;

  if (elemgroupnumber == 6 && level < 4) refine = 1;

  if (elemgroupnumber == 7 && level < 5) refine = 1;

  if (elemgroupnumber == 8 && level < 6) refine = 1;

//   if (elemgroupnumber==6 && level<1) refine=1;
//   if (elemgroupnumber==7 && level<2) refine=1;
//   if (elemgroupnumber==8 && level<3) refine=1;

  return refine;

}



void AssemblePoisson_AD(MultiLevelProblem& ml_prob);    //, unsigned level, const unsigned &levelMax, const bool &assembleMatrix );


int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  //mlMsh.ReadCoarseMesh("./input/cube_hex.neu","seventh",scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/square_quad.neu", "seventh", scalingFactor);
  mlMsh.GenerateCoarseBoxMesh(2, 1, 0, 0., 1., 0, 1, 0., 0., QUAD9, "seventh");
  //mlMsh.ReadCoarseMesh("./input/quadAMR.neu", "seventh", scalingFactor);

  unsigned dim = mlMsh.GetDimension();

  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels, numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // unsigned numberOfUniformLevels = 4;
  // unsigned numberOfSelectiveLevels = 3;
  // mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , SetRefinementFlag);

  // erase all the coarse mesh levels
  //mlMsh.EraseCoarseLevels(numberOfUniformLevels - 3);

  // print mesh info
  mlMsh.PrintInfo();
  MultiLevelSolution mlSol(&mlMsh);
  mlSol.AddSolution("PHI", LAGRANGE, SECOND);
  mlSol.Initialize("All");
  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("All");

  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);
  for (unsigned k = 0; k < 2; k++) {

    MeshRefinement meshcoarser(*(mlMsh.GetLevel(mlMsh.GetNumberOfLevels() - 1)));
    if (k == 0) meshcoarser.FlagOnlyEvenElementsToBeRefined();
    else meshcoarser.FlagAllElementsToBeRefined();

    mlMsh.AddAMRMeshLevel();
    mlSol.AddSolutionLevel();

    vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  }


  return 0;
}
*/



