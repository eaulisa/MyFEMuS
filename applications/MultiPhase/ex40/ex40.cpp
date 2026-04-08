
#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "VTKWriter.hpp"
#include "MyMatrix.hpp"


using namespace femus;

#include "include/Mollifier.hpp"
#include "include/Psi.hpp"
#include "include/GradientApproximation.hpp"
#include "include/BBoxToIel.hpp"



void FlagFinestMeshLevel(MultiLevelMesh &mlMsh, const double &r, const std::vector<double> &xc);
void Init(MultiLevelSolution &mlSol, const std::string &name, const PsiBall &psi2D);
void GetCutElementPoints(MultiLevelSolution &mlSol, const std::string &name, std::vector<MyVector<double>> &X, MyVector<unsigned> &Xiel);
static void WritePointsVTK(const std::string& filename, const std::vector<MyVector<double>>& X);

int main(int argc, char** argv) {

  // Initialize PETSc/MPI
  FemusInit mpinit(argc, argv, MPI_COMM_WORLD);

  MultiLevelMesh mlMsh;

  const double scalingFactor = 1.0;
  const unsigned numberOfUniformLevels   = 1u;
  const unsigned numberOfSelectiveLevels = 0u;

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
  BBoxToIel bboxInfo(mlMsh, 3, 0);


  // Define solution on the multilevel mesh
  MultiLevelSolution mlSol(&mlMsh);
  mlSol.AddSolution("Psi", LAGRANGE, SECOND);
  mlSol.Initialize("All");

  PsiBall psi2D(xc, r, eps);
  Init(mlSol, "Psi", psi2D);

  std::vector<MyVector<double>> X;
  MyVector<unsigned> Xiel;
  GetCutElementPoints(mlSol, "Psi", X, Xiel);

  WritePointsVTK("./output/points.0.vtk", X);
  WritePointsVTK("./output/points.1.vtk", X);

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


void GetCutElementPoints(MultiLevelSolution &mlSol,
                         const std::string &name,
                         std::vector<MyVector<double>> &X,
                         MyVector<unsigned> &Xiel) {

  MultiLevelMesh &mlMsh = *mlSol.GetMultilevelMesh();
  const unsigned level  = mlMsh.GetNumberOfLevels() - 1u;
  Mesh &msh             = *mlMsh.GetLevel(level);
  Solution &sol         = *mlSol.GetLevel(level);
  const unsigned iproc  = msh.processor_id();
  const unsigned dim    = msh.GetDimension();

  const unsigned solIndex = mlSol.GetIndex(name.c_str());
  const unsigned solType  = mlSol.GetSolutionType(name.c_str());
  const unsigned xType    = 2u;

  const double c1 = 2. / 3., c2 = 1. / 3.;

  auto& xv     = msh._topology->_Sol;
  auto& solVec = sol._Sol[solIndex];

  const unsigned offset   = msh._elementOffset[iproc];
  const unsigned offsetp1 = msh._elementOffset[iproc + 1];

  unsigned maxPtsPerEl = 1u;
  for(unsigned d = 0; d < dim; ++d) maxPtsPerEl *= 2u;
  ++maxPtsPerEl;

  std::vector<std::vector<double>> x(dim);
  for(unsigned k = 0; k < dim; ++k) x[k].reserve(maxPtsPerEl);

  std::vector<double> phi;
  std::vector<double> gradPhi(dim);

  std::vector<std::vector<double>> Y(dim);
  std::vector<unsigned> Yiel;

  for(unsigned k = 0; k < dim; ++k) {
    Y[k].reserve((offsetp1 - offset) * maxPtsPerEl);
  }
  Yiel.reserve((offsetp1 - offset) * maxPtsPerEl);

  // Loop over local elements and collect interior points associated with cut elements
  for(unsigned iel = offset; iel < offsetp1; ++iel) {

    const unsigned nDof = msh.GetElementDofNumber(iel, solType);

    const unsigned solDof0 = msh.GetSolutionDof(0, iel, solType);
    const double val0 = (*solVec)(solDof0);
    const int sign0 = (val0 > 0.) - (val0 < 0.);

    for(unsigned i = 1; i < nDof; ++i) {
      const unsigned solDofi = msh.GetSolutionDof(i, iel, solType);
      const double vali = (*solVec)(solDofi);
      const int signi = (vali > 0.) - (vali < 0.);

      if(signi != sign0) {

        const unsigned nDof0 = msh.GetElementDofNumber(iel, 0);
        phi.resize(nDof0 + 1u);
        for(unsigned k = 0; k < dim; ++k) {
          x[k].resize(nDof0 + 1u);
        }
        for(unsigned j = 0; j < nDof0; ++j) {
          const unsigned solDofj = msh.GetSolutionDof(j, iel, solType);
          phi[j] = (*solVec)(solDofj);
          const unsigned xDof = msh.GetSolutionDof(j, iel, xType);
          for(unsigned k = 0; k < dim; ++k) {
            x[k][j] = (*xv[k])(xDof);
          }
        }
        const unsigned nDof2 = msh.GetElementDofNumber(iel, 2);
        const unsigned xDofc = msh.GetSolutionDof(nDof2 - 1u, iel, xType);
        if(solType == xType) {
          phi[nDof0] = (*solVec)(xDofc);

        }
        else {
          phi[nDof0] = 0.;
          for(unsigned  j = 0; j < nDof0; ++j) {
            phi[nDof0] += phi[j];
          }
          phi[nDof0] /= nDof0;
        }
        for(unsigned k = 0; k < dim; ++k) {
          x[k][nDof0] = (*xv[k])(xDofc);
        }

        computeElementGradientFromLocalData(x, phi, gradPhi);

        //shift points
        for(unsigned j = 0; j < nDof0; ++j) {
          phi[j] = c1 * phi[j] + c2 * phi[nDof0];
          for(unsigned k = 0; k < dim; ++k) {
            x[k][j] = c1 * x[k][j] + c2 * x[k][nDof0];
          }
        }

        double gradNorm2 = 0.;
        for(unsigned k = 0; k < dim; ++k) {
          gradNorm2 += gradPhi[k] * gradPhi[k];
        }

        if(gradNorm2 < 1.e-20) {
          for(unsigned j = 0; j <= nDof0; ++j) {
            for(unsigned k = 0; k < dim; ++k) {
              Y[k].push_back(x[k][j]);
            }
            Yiel.push_back(iel);
          }
          break; //use current points
        }

        const double invGradNorm2 = 1. / gradNorm2;
        for(unsigned j = 0; j <= nDof0; ++j) {
          for(unsigned k = 0; k < dim; ++k) {
            Y[k].push_back(x[k][j] - phi[j] * gradPhi[k] * invGradNorm2);
          }
          Yiel.push_back(iel);
        }
        break;
      }
    }
  }

  X.resize(dim);
  for(unsigned k = 0; k < dim; ++k) {
    X[k].buildFromLocal(Y[k]);
  }
  Xiel.buildFromLocal(Yiel);
}


static void WritePointsVTK(const std::string& filename,
                           const std::vector<MyVector<double>>& X) {

  const unsigned dim = X.size();
  if (dim == 0) {
    throw std::runtime_error("writePointsVTK: X.size()==0");
  }
  if (dim > 3) {
    throw std::runtime_error("writePointsVTK: dim > 3 not supported");
  }

  // Localize all components
  std::vector<std::vector<double>> Xp(dim);
  for (unsigned k = 0; k < dim; ++k) {
    X[k].localize(Xp[k]);
  }

  const std::size_t nPts = Xp[0].size();
  for (unsigned k = 1; k < dim; ++k) {
    if (Xp[k].size() != nPts) {
      throw std::runtime_error("writePointsVTK: inconsistent sizes across components");
    }
  }

  int iproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);

  // Only rank 0 writes
  if (iproc != 0) return;

  std::ofstream out(filename);
  if (!out) {
    throw std::runtime_error("writePointsVTK: cannot open file");
  }

  out << "# vtk DataFile Version 3.0\n";
  out << "Point cloud\n";
  out << "ASCII\n";
  out << "DATASET POLYDATA\n";
  out << "POINTS " << nPts << " double\n";

  for (std::size_t i = 0; i < nPts; ++i) {
    const double x = (dim >= 1) ? Xp[0][i] : 0.0;
    const double y = (dim >= 2) ? Xp[1][i] : 0.0;
    const double z = (dim >= 3) ? Xp[2][i] : 0.0;
    out << x << " " << y << " " << z << "\n";
  }

  out << "VERTICES " << nPts << " " << (2 * nPts) << "\n";
  for (std::size_t i = 0; i < nPts; ++i) {
    out << "1 " << i << "\n";
  }
}

// void BuildBoundedBoxCartesianMesh(const MultiLevelMesh& mlMsh,
//                                   MyMatrix<double>& BBox,
//                                   const unsigned nPartition,
//                                   const unsigned level) {
//   if(nPartition == 0u) {
//     throw std::runtime_error("BuildBoundedBoxCartesianMesh: nPartition must be > 0");
//   }
//
//   const Mesh& msh = *mlMsh.GetLevel(level);
//   const unsigned iproc = msh.processor_id();
//   const unsigned dim   = msh.GetDimension();
//   const unsigned xType = 2u;
//
//   auto& xv = msh._topology->_Sol;
//
//   const unsigned offset   = msh._elementOffset[iproc];
//   const unsigned offsetp1 = msh._elementOffset[iproc + 1];
//
//   if(offset == offsetp1) {
//     throw std::runtime_error("BuildBoundedBoxCartesianMesh: empty local partition");
//   }
//
//   std::vector<double> xMin(dim,  DBL_MAX);
//   std::vector<double> xMax(dim, -DBL_MAX);
//   std::vector<double> hCellMin(dim, DBL_MAX);
//
//   std::vector<double> xMinIel(dim);
//   std::vector<double> xMaxIel(dim);
//
//   for(unsigned iel = offset; iel < offsetp1; ++iel) {
//     const unsigned nDof = msh.GetElementDofNumber(iel, xType);
//
//     const unsigned xDof0 = msh.GetSolutionDof(0, iel, xType);
//     for(unsigned k = 0; k < dim; ++k) {
//       const double x0 = (*xv[k])(xDof0);
//       xMinIel[k] = x0;
//       xMaxIel[k] = x0;
//     }
//
//     for(unsigned i = 1; i < nDof; ++i) {
//       const unsigned xDof = msh.GetSolutionDof(i, iel, xType);
//       for(unsigned k = 0; k < dim; ++k) {
//         const double xval = (*xv[k])(xDof);
//         if(xMinIel[k] > xval) xMinIel[k] = xval;
//         if(xMaxIel[k] < xval) xMaxIel[k] = xval;
//       }
//     }
//
//     for(unsigned k = 0; k < dim; ++k) {
//       if(xMin[k] > xMinIel[k]) xMin[k] = xMinIel[k];
//       if(xMax[k] < xMaxIel[k]) xMax[k] = xMaxIel[k];
//
//       const double hIel = xMaxIel[k] - xMinIel[k];
//       if(hCellMin[k] > hIel) hCellMin[k] = hIel;
//     }
//   }
//
//   std::vector<unsigned> N(dim);
//   std::vector<double> hBox(dim);
//
//   for(unsigned k = 0; k < dim; ++k) {
//     if(hCellMin[k] <= 0.0) {
//       throw std::runtime_error("BuildBoundedBoxCartesianMesh: nonpositive local element size");
//     }
//
//     hBox[k] = hCellMin[k] / static_cast<double>(nPartition);
//     xMin[k] -= hBox[k] / 10.0;
//     xMax[k] += hBox[k] / 10.0;
//
//     N[k] = static_cast<unsigned>(std::ceil((xMax[k] - xMin[k]) / hBox[k]));
//     hBox[k] = (xMax[k] - xMin[k]) / static_cast<double>(N[k]);
//
//     std::cout << iproc << " " << xMin[k] << " " << xMax[k] << " " << N[k] << " " << hBox[k] << std::endl;
//   }
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
// }
/*
inline unsigned FlattenBoxIndex2D(const unsigned i,
                                  const unsigned j,
                                  const unsigned Nx) {
  return i + Nx * j;
}

inline unsigned FlattenBoxIndex3D(const unsigned i,
                                  const unsigned j,
                                  const unsigned k,
                                  const unsigned Nx,
                                  const unsigned Ny) {
  return i + Nx * (j + Ny * k);
}*/

/*
void BuildBoundedBoxCartesianMesh(const MultiLevelMesh& mlMsh,
                                  std::vector<std::vector<unsigned>> &bboxToIel,
                                  const unsigned nPartition,
                                  const unsigned level) {
  if(nPartition == 0u) {
    throw std::runtime_error("BuildBoundedBoxCartesianMesh: nPartition must be > 0");
  }

  const Mesh& msh    = *mlMsh.GetLevel(level);
  const unsigned iproc = msh.processor_id();
  const unsigned dim   = msh.GetDimension();
  const unsigned xType = 2u;

  auto& xv = msh._topology->_Sol;

  const unsigned offset   = msh._elementOffset[iproc];
  const unsigned offsetp1 = msh._elementOffset[iproc + 1];

  if(dim < 1u || dim > 3u) {
    throw std::runtime_error("BuildBoundedBoxCartesianMesh: dim must be 1, 2, or 3");
  }

  if(offset == offsetp1) {
    return;
  }

  std::vector<double> xMin(dim,  DBL_MAX);
  std::vector<double> xMax(dim, -DBL_MAX);
  std::vector<double> hCellMin(dim, DBL_MAX);

  std::vector<double> xMinIel(dim);
  std::vector<double> xMaxIel(dim);

  // --------------------------------------------------
  // First pass: local partition bounding box and min h
  // --------------------------------------------------
  for(unsigned iel = offset; iel < offsetp1; ++iel) {

    const unsigned nDof = msh.GetElementDofNumber(iel, xType);

    const unsigned xDof0 = msh.GetSolutionDof(0, iel, xType);
    for(unsigned k = 0; k < dim; ++k) {
      const double x0 = (*xv[k])(xDof0);
      xMinIel[k] = x0;
      xMaxIel[k] = x0;
    }

    for(unsigned i = 1; i < nDof; ++i) {
      const unsigned xDof = msh.GetSolutionDof(i, iel, xType);
      for(unsigned k = 0; k < dim; ++k) {
        const double xval = (*xv[k])(xDof);
        if(xMinIel[k] > xval) xMinIel[k] = xval;
        if(xMaxIel[k] < xval) xMaxIel[k] = xval;
      }
    }

    for(unsigned k = 0; k < dim; ++k) {
      if(xMin[k] > xMinIel[k]) xMin[k] = xMinIel[k];
      if(xMax[k] < xMaxIel[k]) xMax[k] = xMaxIel[k];

      const double hIel = xMaxIel[k] - xMinIel[k];
      if(hCellMin[k] > hIel) hCellMin[k] = hIel;
    }
  }

  std::vector<unsigned> N(dim, 0u);
  std::vector<double> hBox(dim, 0.0);

  for(unsigned k = 0; k < dim; ++k) {
    if(hCellMin[k] <= 0.0) {
      throw std::runtime_error("BuildBoundedBoxCartesianMesh: nonpositive local element size");
    }

    hBox[k] = hCellMin[k] / static_cast<double>(nPartition);
    xMin[k] -= hBox[k] / 10.0;
    xMax[k] += hBox[k] / 10.0;

    N[k] = static_cast<unsigned>(std::ceil((xMax[k] - xMin[k]) / hBox[k]));
    if(N[k] == 0u) N[k] = 1u;

    hBox[k] = (xMax[k] - xMin[k]) / static_cast<double>(N[k]);

    std::cout << iproc << " " << xMin[k] << " " << xMax[k] << " " << N[k] << " " << hBox[k] << std::endl;
  }


  // --------------------------------------------------
  // Allocate candidate list: one list per BBox cell
  // --------------------------------------------------
  unsigned nBoxElem = 1u;
  for(unsigned k = 0; k < dim; ++k) {
    nBoxElem *= N[k];
  }

  bboxToIel.clear();
  bboxToIel.resize(nBoxElem);
  for(unsigned i = 0; i < bboxToIel.size(); ++i) bboxToIel[i].reserve(10);

  // Small tolerance against roundoff at box boundaries
  std::vector<double> eps(dim, 0.0);
  for(unsigned k = 0; k < dim; ++k) {
    eps[k] = 1.e-12 * hBox[k];
  }

  // --------------------------------------------------
  // Second pass: for each unstructured element,
  // find BBox index range in each dimension,
  // then append iel to all covered BBox cells
  // --------------------------------------------------
  for(unsigned iel = offset; iel < offsetp1; ++iel) {

    const unsigned nDof = msh.GetElementDofNumber(iel, xType);

    const unsigned xDof0 = msh.GetSolutionDof(0, iel, xType);
    for(unsigned k = 0; k < dim; ++k) {
      const double x0 = (*xv[k])(xDof0);
      xMinIel[k] = x0;
      xMaxIel[k] = x0;
    }

    for(unsigned i = 1; i < nDof; ++i) {
      const unsigned xDof = msh.GetSolutionDof(i, iel, xType);
      for(unsigned k = 0; k < dim; ++k) {
        const double xval = (*xv[k])(xDof);
        if(xMinIel[k] > xval) xMinIel[k] = xval;
        if(xMaxIel[k] < xval) xMaxIel[k] = xval;
      }
    }

    std::vector<unsigned> iMin(dim, 0u);
    std::vector<unsigned> iMax(dim, 0u);

    for(unsigned k = 0; k < dim; ++k) {
      int ikMin = static_cast<int>(std::floor((xMinIel[k] - xMin[k] - eps[k]) / hBox[k]));
      int ikMax = static_cast<int>(std::floor((xMaxIel[k] - xMin[k] + eps[k]) / hBox[k]));

      if(ikMin < 0) ikMin = 0;
      if(ikMax < 0) ikMax = 0;
      if(ikMin >= static_cast<int>(N[k])) ikMin = static_cast<int>(N[k]) - 1;
      if(ikMax >= static_cast<int>(N[k])) ikMax = static_cast<int>(N[k]) - 1;

      if(ikMax < ikMin) std::swap(ikMin, ikMax);

      iMin[k] = static_cast<unsigned>(ikMin);
      iMax[k] = static_cast<unsigned>(ikMax);
    }

    if(dim == 1u) {
      for(unsigned i = iMin[0]; i <= iMax[0]; ++i) {
        const unsigned ibox = i;
        bboxToIel[ibox].push_back(iel);
      }
    }
    else if(dim == 2u) {
      for(unsigned j = iMin[1]; j <= iMax[1]; ++j) {
        for(unsigned i = iMin[0]; i <= iMax[0]; ++i) {
          const unsigned ibox = FlattenBoxIndex2D(i, j, N[0]);
          bboxToIel[ibox].push_back(iel);
        }
      }
    }
    else {
      for(unsigned k = iMin[2]; k <= iMax[2]; ++k) {
        for(unsigned j = iMin[1]; j <= iMax[1]; ++j) {
          for(unsigned i = iMin[0]; i <= iMax[0]; ++i) {
            const unsigned ibox = FlattenBoxIndex3D(i, j, k, N[0], N[1]);
            bboxToIel[ibox].push_back(iel);
          }
        }
      }
    }
  }

  for(unsigned ibox = 0; ibox < bboxToIel.size(); ++ibox) {
    std::cout << "iproc = " << iproc
              << ", ibox = " << ibox
              << ", nCandidates = " << bboxToIel[ibox].size() << std::endl;
  }

  // --------------------------------------------------
  // Here you would still build BBox itself if needed.
  // This code only builds bboxToIel.
  // --------------------------------------------------
}*/
