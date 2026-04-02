
#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "VTKWriter.hpp"

#include "include/Mollifier.hpp"
#include "include/Psi.hpp"

using namespace femus;

void FlagFinestMeshLevel(MultiLevelMesh &mlMsh, const double &r, const std::vector<double> &xc);
void Init(MultiLevelSolution &mlSol, const std::string &name, const PsiBall &psi2D);

int main(int argc, char** argv) {

  // Initialize PETSc/MPI
  FemusInit mpinit(argc, argv, MPI_COMM_WORLD);

  MultiLevelMesh mlMsh;

  const double scalingFactor = 1.0;
  const unsigned numberOfUniformLevels   = 2u;
  const unsigned numberOfSelectiveLevels = 4u;

  // Load coarse mesh and build uniform refinement levels
  mlMsh.ReadCoarseMesh("./input/square.neu", "seventh", scalingFactor);
  mlMsh.RefineMesh(numberOfUniformLevels, numberOfUniformLevels, nullptr);

  // Parameters for selective AMR (ball centered at xc with radius r)
  const double r = 0.125;
  std::vector<double> xc = {0.0, 0.25};
  double eps = 0.01;//(dim == 2) ? 1. / pow(2, std::max(levelN - 7u, 1u))

  // Iteratively flag and create new AMR levels
  for (unsigned k = 0; k < numberOfSelectiveLevels; ++k) {
    FlagFinestMeshLevel(mlMsh, r, xc);
    mlMsh.AddAMRMeshLevel(false); // false -> it does not re-evaluate the AMR flag vector
  }

  mlMsh.PrintInfo();

  // Define solution on the multilevel mesh
  MultiLevelSolution mlSol(&mlMsh);
  mlSol.AddSolution("Psi", LAGRANGE, SECOND);
  mlSol.Initialize("All");

  PsiBall psi2D(xc, r, eps);
  Init(mlSol, "Psi", psi2D);


  // Export solution to VTK (selected levels)
  std::vector<std::string> variablesToBePrinted = {"All"};
  VTKWriter vtkIO(&mlSol);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 1);

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

  auto& xv = msh._topology->_Sol;
  auto* amrFlag = msh._topology->_Sol[amrIndex];

  amrFlag->zero();

  unsigned offset = msh._elementOffset[iproc];
  unsigned offsetp1 = msh._elementOffset[iproc + 1];

  for (unsigned iel = offset; iel < offsetp1; ++iel) {

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
      amrFlag->set(iel, 2);
    }
  }
  amrFlag->close();


  for (unsigned iel = offset; iel < offsetp1; ++iel) {
    if( (*amrFlag)(iel) == 2) {
      for(unsigned j = 0; j < msh.el->GetElementNearElementSize(iel, 1); j++) {
        unsigned jel = msh.el->GetElementNearElement(iel, j);
        if(offset <= jel && jel < offsetp1) {
          if((*amrFlag)(jel) < 0.5) amrFlag->set(jel, 1); //this is on spot since jel belongs to iproc
        }
        else {
          amrFlag->add(jel, 1); // this is buffered since jel does not belong to iproc
        }
      }
    }
  }
  amrFlag->close();

  unsigned localRefined = 0;
  for (unsigned iel = offset; iel < offsetp1; ++iel) {
    if( (*amrFlag)(iel) > 0.5) {
      amrFlag->set(iel, 1);
      ++localRefined;
    }
  }
  amrFlag->close();

  unsigned globalRefined = 0;
  MPI_Allreduce(&localRefined, &globalRefined, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

  msh.el->SetRefinedElementNumber(globalRefined);

}


void Init(MultiLevelSolution &mlSol, const std::string &name, const PsiBall &psi2D) {

  MultiLevelMesh &mlMsh = *mlSol.GetMultilevelMesh();
  const unsigned level  = mlMsh.GetNumberOfLevels() - 1u;
  Mesh &msh             = *mlMsh.GetLevel(level);
  Solution &sol         = *mlSol.GetLevel(level);
  const unsigned iproc  = msh.processor_id();
  const unsigned dim    = msh.GetDimension();

  unsigned solIndex = mlSol.GetIndex(name.c_str());
  unsigned solType  = mlSol.GetSolutionType(name.c_str());

  const unsigned xType = 2u; // coordinate field type

  auto& xv = msh._topology->_Sol;
  std::vector<double> x(dim);

  const unsigned offset   = msh._elementOffset[iproc];
  const unsigned offsetp1 = msh._elementOffset[iproc + 1];

  auto& solVec = sol._Sol[solIndex];

  solVec->zero();

  // Loop over local elements and interpolate psi2D at solution DoFs
  for (unsigned iel = offset; iel < offsetp1; ++iel) {

    const unsigned nDof = msh.GetElementDofNumber(iel, solType);

    for (unsigned i = 0; i < nDof; ++i) {

      // Get physical coordinates of the current DoF
      const unsigned xDof = msh.GetSolutionDof(i, iel, xType);
      for (unsigned k = 0; k < dim; ++k) {
        x[k] = (*xv[k])(xDof);
      }

      // Evaluate and assign field value
      const unsigned solDof = msh.GetSolutionDof(i, iel, solType);
      solVec->set(solDof, psi2D(x));
    }
  }

  solVec->close();
}

