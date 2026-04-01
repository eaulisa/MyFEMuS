/** tutorial/Ex1
 * This example shows how to:
 * initialize a femus application;
 * define the multilevel-mesh object mlMsh;
 * read from the file ./input/square.neu the coarse-level mesh and associate it to mlMsh;
 * add in mlMsh uniform refined level-meshes;
 * define the multilevel-solution object mlSol associated to mlMsh;
 * add in mlSol different types of finite element solution variables; +
 * initialize the solution variables;
 * define vtk and gmv writer objects associated to mlSol;
 * print vtk and gmv binary-format files in ./output directory.
 **/

#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"

using namespace femus;

double InitalValueU(const std::vector < double >& x) {
  return x[0] + x[1];
}

double InitalValueP(const std::vector < double >& x) {
  return x[0];
}

double InitalValueT(const std::vector < double >& x) {
  return x[1];
}


void FlagFinestMeshLevel(MultiLevelMesh &mlMsh, const double &r, const std::vector<double> &xc);

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh; // Consider the "mlMsh" object of "MultiLevel" class.
  double scalingFactor = 1.;
  // read coarse level mesh and generate finers level meshes
  mlMsh.ReadCoarseMesh("./input/square.neu", "seventh", scalingFactor); // Let mlMsh read the coarse mesh.
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
      probably in the furure it is not going to be an argument of this function   */
  unsigned numberOfUniformLevels = 2; // We apply uniform refinement.
  unsigned numberOfSelectiveLevels = 0; // we may want to see solutions on different level of meshes.
  mlMsh.RefineMesh(numberOfUniformLevels, numberOfUniformLevels + numberOfSelectiveLevels, NULL);  // Add those refined meshed in mlMsh object.

  double r = 0.125;
  std::vector<double> xc = {0.0, 0.25};

  for(unsigned k = 0; k < 5; k++) {
    FlagFinestMeshLevel(mlMsh, r, xc);
    mlMsh.AddAMRMeshLevel(false);
  }
  mlMsh.PrintInfo();


  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol(&mlMsh); // Define "mlSol" object with argument mlMsh of "MultiLevelSolution" class.
  // add variables to mlSol
  mlSol.AddSolution("U", LAGRANGE, FIRST);
  mlSol.AddSolution("V", LAGRANGE, SERENDIPITY);
  mlSol.AddSolution("W", LAGRANGE, SECOND);
  mlSol.AddSolution("P", DISCONTINUOUS_POLYNOMIAL, ZERO);
  mlSol.AddSolution("T", DISCONTINUOUS_POLYNOMIAL, FIRST);

  mlSol.Initialize("All");    // initialize all varaibles to zero

  mlSol.Initialize("U", InitalValueU);
  mlSol.Initialize("V", InitalValueP);
  mlSol.Initialize("W", InitalValueT);

  mlSol.Initialize("P", InitalValueP);
  mlSol.Initialize("T", InitalValueT);  // note that this initialization is the same as piecewise constant element

  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");
//   variablesToBePrinted.push_back("U");
//   variablesToBePrinted.push_back("P");
//   variablesToBePrinted.push_back("T");
//   variablesToBePrinted.push_back("V");
//   variablesToBePrinted.push_back("W");

  VTKWriter vtkIO(&mlSol);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

//   GMVWriter gmvIO(&mlSol);
//   variablesToBePrinted.push_back("all");
//   gmvIO.SetDebugOutput(false);
//   gmvIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  //   From what I understand, you define the type of mesh and how to refine it and what integration scheme to use,
  //then you associate the mesh and integration scheme to a specific type of problem, with specifed order, then solve and print.

  return 0;
}

void FlagFinestMeshLevel(MultiLevelMesh &mlMsh, const double &r, const std::vector<double> &xc) {

  const unsigned level    = mlMsh.GetNumberOfLevels() - 1;
  Mesh &msh               = *mlMsh.GetLevel(level);
  const unsigned iproc    = msh.processor_id();
  const unsigned dim      = msh.GetDimension();
  const unsigned xType    = 2;
  const unsigned amrIndex = msh.GetAmrIndex();

  const double r2 = r * r;
  unsigned localRefined = 0;

  auto& xv = msh._topology->_Sol;
  auto* amrFlag = msh._topology->_Sol[amrIndex];

  amrFlag->zero();

  for (unsigned iel = msh._elementOffset[iproc]; iel < msh._elementOffset[iproc + 1]; ++iel) {

    const unsigned nDof = msh.GetElementDofNumber(iel, xType);

    const unsigned xDof0 = msh.GetSolutionDof(0, iel, xType);
    double d2_0 = r2;
    for (unsigned k = 0; k < dim; ++k) {
      const double d = (*xv[k])(xDof0) - xc[k];
      d2_0 -= d * d;
    }
    const int sign0 = (d2_0 > 0.) ? 1 : -1;

    bool signChange = false;
    for (unsigned i = 1; i < nDof; ++i) {
      const unsigned xDof = msh.GetSolutionDof(i, iel, xType);
      double d2 = r2;

      for (unsigned k = 0; k < dim; ++k) {
        const double d = (*xv[k])(xDof) - xc[k];
        d2 -= d * d;
      }

      const int signi = (d2 > 0.) ? 1 : -1;
      if (signi != sign0) {
        signChange = true;
        break;
      }
    }

    if (signChange) {
      amrFlag->set(iel, 1.);
      ++localRefined;
    }
  }

  amrFlag->close();

  unsigned globalRefined = 0;
  MPI_Allreduce(&localRefined, &globalRefined, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

  msh.el->SetRefinedElementNumber(globalRefined);
}


