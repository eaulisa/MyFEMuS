#include <iostream>
#include <vector>
#include <cmath>
#include <memory>
#include "Mesh.hpp"
#include "FemProjection.hpp"   // FemProjection, Quad9Projection, Tri7Projection
#include "RefineMesh.hpp"
#include "VtkOutput.hpp"

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

  unsigned levelN = 5;
  std::vector<std::vector<std::vector<double>>> X(levelN);
  // for (unsigned l = 0; l < levelN; l++) {
  //   X[l].resize(dim);
  // }
  std::vector<std::vector<unsigned>> elType(levelN);
  std::vector<std::vector<std::vector<unsigned>>> elTplgy(levelN);
  std::vector<std::vector<unsigned>> elLevel(levelN);
  std::vector<std::vector<int>> AMR(levelN);
  std::vector<std::vector<unsigned>> elFather(levelN);
  std::vector<std::vector<std::vector<unsigned>>> elChildern(levelN);
  std::vector<std::vector<std::vector<unsigned>>> nodeElements(levelN);

  std::vector<Mesh> mesh;;
  mesh.reserve(levelN);

  for (unsigned l = 0; l < levelN; l++) {
    mesh.emplace_back(elTplgy[l], elType[l], elLevel[l], X[l], AMR[l], elFather[l], elChildern[l], nodeElements[l]);
  }


  elLevel[0] = {0, 0};
  elType[0] = {3, 4};
  elTplgy[0] = {
    {0, 1, 2, 3, 5, 6, 7, 8, 11},
    {1, 4, 2, 9, 10, 6, 12}
  };
  X[0] = {
    {0, 1, 1, 0., 2., 0.5, 1., 0.5, 0, 1.5, 1.5, 0.5, 1 + 1. / 3.},
    {0, 0, 1, 1, 0, 0, 0.5, 1, 0.5, 0, 0.5, 0.5, 1. / 3}
  };
  mesh[0].resetAllFathersToNoFather();
  mesh[0].buildNodeToElementAdjacency();

  std::string filename = "./output/refined_mesh2D.";
  writeMeshVTK(filename + "0.vtk", mesh[0]);

  AMR[0] = {1, 0};
  refineAndProjectMesh(elProj, mesh[0], mesh[1]);
  writeMeshVTK(filename + "1.vtk", mesh[1]);

  for (unsigned l = 2; l < levelN; l++) {
    mesh[l - 1].setRefineEvenElements();
    mesh[l - 1].adjustAMRForOneLevelDiscontinuity();
    refineAndProjectMesh(elProj, mesh[l - 1], mesh[l]);
    writeMeshVTK(filename + std::to_string(l) + ".vtk", mesh[l]);
  }

  elType[0] = {0, 2, 1};
  elLevel[0] = {0, 0, 0};
  elTplgy[0] = {
    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26},
    {27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47},
    {48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62}
  };

  X[0] = {{
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

  std::vector<unsigned> candidateIndices(mesh[0].numNodes());
  for (unsigned i = 0; i < candidateIndices.size(); i++) candidateIndices[i] = i;
  dedupNodesHash(mesh[0], candidateIndices);
  std::cout << "Number of nodes after deduplication = " << mesh[0].numNodes() << std::endl;
  mesh[0].resetAllFathersToNoFather();
  mesh[0].buildNodeToElementAdjacency();



  filename = "./output/refined_mesh3D.";
  writeMeshVTK(filename + "0.vtk", mesh[0]);

  AMR[0] = {1, 0, 1};
  refineAndProjectMesh(elProj, mesh[0], mesh[1]);
  writeMeshVTK(filename + "1.vtk", mesh[1]);

  for (unsigned l = 2; l < levelN; l++) {
    mesh[l - 1].setRefineEvenElements();
    mesh[l - 1].adjustAMRForOneLevelDiscontinuity();
    refineAndProjectMesh(elProj, mesh[l - 1], mesh[l]);
    writeMeshVTK(filename + std::to_string(l) + ".vtk", mesh[l]);
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

