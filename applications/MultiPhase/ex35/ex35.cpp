#include <iostream>
#include <vector>
#include <cmath>
#include <memory>
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

int main() {

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

  unsigned dim = 2;

  unsigned levelN = 12;
  std::vector<std::vector<std::vector<double>>> X0(levelN);
  std::vector<std::vector<unsigned>> elType0(levelN);
  std::vector<std::vector<std::vector<unsigned>>> elTplgy0(levelN);
  std::vector<std::vector<unsigned>> elLevel0(levelN);
  std::vector<std::vector<int>> AMR0(levelN);
  std::vector<std::vector<unsigned>> elFather0(levelN);
  std::vector<std::vector<std::vector<unsigned>>> elChildern0(levelN);
  std::vector<std::vector<std::vector<unsigned>>> nodeElements0(levelN);
  std::vector<std::vector<std::vector<unsigned>>> neighbors0(levelN);


  std::vector<std::vector<std::vector<double>>> X1(levelN);
  std::vector<std::vector<unsigned>> elType1(levelN);
  std::vector<std::vector<std::vector<unsigned>>> elTplgy1(levelN);
  std::vector<std::vector<unsigned>> elLevel1(levelN);
  std::vector<std::vector<int>> AMR1(levelN);
  std::vector<std::vector<unsigned>> elFather1(levelN);
  std::vector<std::vector<std::vector<unsigned>>> elChildern1(levelN);
  std::vector<std::vector<std::vector<unsigned>>> nodeElements1(levelN);
  std::vector<std::vector<std::vector<unsigned>>> neighbors1(levelN);



  std::vector<Mesh> mesh0, mesh1;
  std::vector<Field> field0, field1;
  mesh0.reserve(levelN);
  field0.reserve(levelN);

  mesh1.reserve(levelN);
  field1.reserve(levelN);

  for (unsigned l = 0; l < levelN; l++) {
    mesh0.emplace_back(elTplgy0[l], elType0[l], elLevel0[l], X0[l], AMR0[l], elFather0[l], elChildern0[l], nodeElements0[l], neighbors0[l]);
    field0.emplace_back(mesh0[l]);

    mesh1.emplace_back(elTplgy1[l], elType1[l], elLevel1[l], X1[l], AMR1[l], elFather1[l], elChildern1[l], nodeElements1[l], neighbors1[l]);
    field1.emplace_back(mesh1[l]);
  }


  elLevel0[0] = elLevel1[0] = {0, 0};
  elType0[0] = elType1[0] = {3, 4};
  elTplgy0[0] = elTplgy1[0] = {
    {0, 1, 2, 3, 5, 6, 7, 8, 11},
    {1, 4, 2, 9, 10, 6, 12}
  };
  X0[0] = X1[0] = {
    { -0.5,  0.5,  0.5, -0.5,  1.5,  0.0,  0.5,  0.0, -0.5,  1.0,  1.0,  0.0,  0.8333333333333333 },
    { -0.5, -0.5,  0.5,  0.5, -0.5, -0.5,  0.0,  0.5,  0.0, -0.5,  0.0,  0.0, -0.16666666666666669 }
  };

// X0[0] = X1[0] = {
//    {0, 1, 1, 0., 2., 0.5, 1., 0.5, 0, 1.5, 1.5, 0.5, 1 + 1. / 3.},
//    {0, 0, 1, 1, 0, 0, 0.5, 1, 0.5, 0, 0.5, 0.5, 1. / 3}
//  };
  mesh0[0].resetAllFathersToNoFather();
  mesh0[0].buildNodeToElementAdjacency();
  mesh0[0].buildFaceNeighborsFromNodeToElement();

  mesh1[0].resetAllFathersToNoFather();
  mesh1[0].buildNodeToElementAdjacency();
  mesh1[0].buildFaceNeighborsFromNodeToElement();

  field0[0].addField("Psi", Field::Location::Nodal, 1.);

  PsiBall psi2D(std::vector<double> {.0, 0.25}, 0.15, 0.001);
  unsigned psiId = field0[0].id("Psi");
  auto& Psi = field0[0].getById(psiId);
  for (std::size_t k = 0; k < Psi.size(); ++k) {
    const std::vector<double> x = field0[0].dofCoordById(psiId, k); // nodal => mesh node coord
    Psi[k] = psi2D(x);
  }

  std::string filename = "./output/refined_mesh2D.";

  unsigned neighMode = 3; // 0=no-ring, 1=vertices, 2=faces, 3=hybrid
  for (unsigned l = 1; l < levelN; l++) {
    if (l == 1) mesh0[l - 1].setUniformRefinement();
    else {
      mesh0[l - 1].setRefinementFromBallLevelSetCrossing_OneRing({0., 0.25}, 0.15, neighMode);
      mesh0[l - 1].adjustAMRForOneLevelDiscontinuity();
    }
    refineAndProjectMesh(elProj, mesh0[l - 1], mesh0[l]);
    field0[l].addField("Psi", Field::Location::Nodal, 1.);

    const unsigned psiId = field0[l].id("Psi");
    auto& Psi = field0[l].getById(psiId);
    for (std::size_t k = 0; k < Psi.size(); ++k) {
      const std::vector<double> x = field0[l].dofCoordById(psiId, k); // nodal => mesh node coord
      Psi[k] = psi2D(x);
    }
  }

  writeMeshFieldVTU(filename + std::to_string(0) + ".vtu", field0[levelN - 1]);

  std::vector<std::vector<double>> Xp;
  PointLocator pl0 = PointLocator(mesh0[0], .1);
  PointLocator pl1 = PointLocator(mesh1[0], .1);
  std::vector<PointLocatorResult> out, in;
  //writePointsVTK("./output/points2D.0.vtk", Xp);

  const unsigned nIter = 320;
  const double period = 8.;
  const double dt = period / static_cast<double>(nIter);

  for (unsigned k = 1; k <= nIter; ++k) {
    const double time = k * dt;

    field0[levelN - 1].extractInterfaceVerticesAndCentersByName("Psi", Xp, levelN - 1, levelN);
    RungeKutta::rkForward(Xp, time, dt, RungeKutta::VelKind::Vortex);
    pl1.locateAll(out, Xp);

    mesh1[0].setRefinementFromLocatedPoints_OneRing(out, neighMode);
    mesh1[0].adjustAMRForOneLevelDiscontinuity();

    for (unsigned l = 1; l < levelN; l++) {
      mesh1[l].clearAllData();
      refineAndProjectMesh(elProj, mesh1[l - 1], mesh1[l]);

      std::swap(in, out);
      mesh1[l - 1].projectPointLocatorResultsToNextLevel(mesh1[l], in, out);
      mesh1[l].setRefinementFromLocatedPoints_OneRing(out, neighMode);
      mesh1[l].adjustAMRForOneLevelDiscontinuity();

      field1[l].clear();
      field1[l].addField("Psi", Field::Location::Nodal, 1.);
      const unsigned psiId = field1[l].id("Psi");

      if (l < levelN - 1) {
        auto& Psi = field1[l].getById(psiId);
        for (std::size_t k = 0; k < Psi.size(); ++k) {
          Psi[k] = 0.;
        }
      }
      else {
        Xp = mesh1[levelN - 1].X();
        RungeKutta::rkBackward(Xp, time, dt, RungeKutta::VelKind::Vortex);
        pl0.locateAll(out, Xp);
        for (unsigned j = 1; j < levelN; j++) {
          std::swap(in, out);
          mesh0[j - 1].projectPointLocatorResultsToNextLevel(mesh0[j], in, out);
        }
        auto& Psi = field1[l].getById(psiId);
        const unsigned psiId0 = field0[levelN - 1].id("Psi");
        field0[levelN - 1].evalNodalAtLocatedPointsById(psiId0, out, elProj, Psi, -1.0);
      }
    }
    for (unsigned l = 0; l < levelN; l++) {
      swap(mesh0[l], mesh1[l]);
      swap(field0[l], field1[l]);
    }
    writeMeshFieldVTU(filename + std::to_string(k) + ".vtu", field0[levelN - 1]);
  }

  return 1;

  //RungeKutta::rkBackward(Xp, time, dt, RungeKutta::VelKind::Vortex);

  //writePointsVTK("./output/points2D.2.vtk", Xp);


  mesh0[0].clearAllData();


  elType0[0] = {0, 2, 1};
  elLevel0[0] = {0, 0, 0};
  elTplgy0[0] = {
    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26},
    {27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47},
    {48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62}
  };

  X0[0] = {{
      //hex
      0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
      0.5, 1.0, 0.5, 0.0, 0.5, 1.0, 0.5, 0.0, 0.0, 1.0, 1.0, 0.0,
      0.5, 1.0, 0.5, 0.0, 0.5, 0.5,
      0.5,
      //wedge
      1.0, 2.0, 1.0, 1.0, 2.0, 1.0,
      1.5, 1.5, 1.0, 1.5, 1.5, 1.0, 1.0, 2.0, 1.0,
      1.5, 1.5, 1.0, 1.3333333333333333, 1.3333333333333333,
      1.3333333333333333,
      //tet
      1.0, 2.0, 1.0, 1.0,
      1.5, 1.5, 1.0, 1.0, 1.5, 1.0,
      1.3333333333333333, 1.3333333333333333,
      1.3333333333333333, 1.0,
      1.25
    }, {
      //hex
      0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0,
      0.0, 0.5, 1.0, 0.5, 0.0, 0.5, 1.0, 0.5, 0.0, 0.0, 1.0, 1.0,
      0.0, 0.5, 1.0, 0.5, 0.5, 0.5,
      0.5,
      //wedge
      0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
      0.0, 0.5, 0.5, 0.0, 0.5, 0.5, 0.0, 0.0, 1.0,
      0.0, 0.5, 0.5, 1.0 / 3.0, 1.0 / 3.0,
      1.0 / 3.0,
      //tet
      0.0, 0.0, 1.0, 0.0,
      0.0, 0.5, 0.5, 0.0, 0.0, 0.5,
      0.3333333333333333, 0.0,
      0.3333333333333333, 0.3333333333333333,
      0.25
    }, {
      //hex
      0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,
      0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5,
      0.5, 0.5, 0.5, 0.5, 0.0, 1.0,
      0.5,
      //wedge
      0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
      0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.5, 0.5, 0.5,
      0.5, 0.5, 0.5, 0.0, 1.0,
      0.5,
      //tet
      1.0, 1.0, 1.0, 2.0,
      1.0, 1.0, 1.0, 1.5, 1.5, 1.5,
      1.0, 1.3333333333333333,
      1.3333333333333333, 1.3333333333333333,
      1.25
    }
  };

  std::vector<unsigned> candidateIndices(mesh0[0].numNodes());
  for (unsigned i = 0; i < candidateIndices.size(); i++) candidateIndices[i] = i;
  dedupNodesHash(mesh0[0], candidateIndices);
  std::cout << "Number of nodes after deduplication = " << mesh0[0].numNodes() << std::endl;
  mesh0[0].resetAllFathersToNoFather();
  mesh0[0].buildNodeToElementAdjacency();
  mesh0[0].buildFaceNeighborsFromNodeToElement();

  field0[0].rebindMeshAndResize(mesh0[0]);

  PsiBall psi3D(std::vector<double> {1., 0, 1.}, 0.125, 0.001);
  psiId = field0[0].id("Psi");
  Psi = field0[0].getById(psiId);
  for (std::size_t k = 0; k < Psi.size(); ++k) {
    const std::vector<double> x = field0[0].dofCoordById(psiId, k); // nodal => mesh node coord
    Psi[k] = psi3D(x);
  }

  filename = "./output/refined_mesh3D.";
  writeMeshFieldVTK(filename + "0.vtk", field0[0]);


  for (unsigned l = 1; l < levelN; l++) {
    mesh0[l - 1].setRefinementFromBallLevelSetCrossing_OneRing({1., 0., 1.}, 0.125, neighMode);
    mesh0[l - 1].adjustAMRForOneLevelDiscontinuity();
    refineAndProjectMesh(elProj, mesh0[l - 1], mesh0[l]);

    field0[l].rebindMeshAndResize(mesh0[l]);
    const unsigned psiId = field0[l].id("Psi");
    auto& Psi = field0[l].getById(psiId);
    for (std::size_t k = 0; k < Psi.size(); ++k) {
      const std::vector<double> x = field0[l].dofCoordById(psiId, k); // nodal => mesh node coord
      Psi[k] = psi3D(x);
    }
    writeMeshFieldVTK(filename + std::to_string(l) + ".vtk", field0[l]);
  }

  return 0;
}
























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


