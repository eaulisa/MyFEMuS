#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "MultiLevelSolution.hpp"

#include "MeshRefinement.hpp"




using namespace femus;


#include "marker.hpp"
#include "background.hpp"
#include "projection.hpp"

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  MultiLevelMesh mlMshM;
  mlMshM.ReadCoarseMesh("../input/benchmarkFSISolid.neu", "fifth", 1);
  mlMshM.RefineMesh(3, 3, NULL); //uniform refinement, this goes with the background mesh refinement. For COMSOL we use 8 = 3 turekBeam2D

  unsigned numberOfRefinement = 3;
  for(unsigned i = 0; i < numberOfRefinement; i++) {
    FlagElements(mlMshM, 2);
    mlMshM.AddAMRMeshLevel();
  }
  mlMshM.EraseCoarseLevels(numberOfRefinement);

  MultiLevelSolution mlSolM(&mlMshM);
  InitializeMarkerVariables(mlSolM);

  UpdateMeshQuantities(mlSolM);

  MultiLevelMesh mlMshB;
  mlMshB.ReadCoarseMesh("../input/benchmarkFSI.neu", "fifth", 1);
  mlMshB.RefineMesh(3, 3, NULL);

  MultiLevelSolution mlSolB(&mlMshB);
  InitializeBackgroundVariables(mlSolB);

  //******* Print solution *******
  mlSolM.SetWriter(VTK);
  mlSolB.SetWriter(VTK);


  unsigned dim = mlMshM.GetDimension();
  std::vector<std::string> mov_vars;
  mov_vars.push_back("DX");
  mov_vars.push_back("DY");
  if(dim == 3) mov_vars.push_back("DZ");
  mlSolM.GetWriter()->SetMovingMesh(mov_vars);
  mlSolB.GetWriter()->SetMovingMesh(mov_vars);
  mlSolM.GetWriter()->SetDebugOutput(false);
  mlSolB.GetWriter()->SetDebugOutput(false);
  std::vector<std::string> print_vars;
  print_vars.push_back("All");

  mlSolM.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, 0);
  mlSolB.GetWriter()->Write(DEFAULT_OUTPUTDIR, "linear", print_vars, 0);

  Projection projection(&mlSolM, &mlSolB);
  for(unsigned t = 1; t <= 10; t++) {
    clock_t time = clock();
    projection.Init();
    projection.FromMarkerToBackground();
    projection.FakeMovement();
    projection.FromBackgroundToMarker();
    std::cout << "time" << t << " = "<< static_cast<double>((clock() - time)) / CLOCKS_PER_SEC << std::endl;
    
    UpdateMeshQuantities(mlSolM);

    mlSolM.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, t);
    mlSolB.GetWriter()->Write(DEFAULT_OUTPUTDIR, "linear", print_vars, t);
  }

  return 0;

} //end main

