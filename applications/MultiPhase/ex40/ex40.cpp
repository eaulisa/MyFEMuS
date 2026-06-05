
#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MyMatrix.hpp"
#include "VTKWriter.hpp"

#include <gperftools/profiler.h>

using namespace femus;

#include "include/BBoxToIel.hpp"
#include "include/ElementTopology.hpp"
#include "include/GradientApproximation.hpp"
#include "include/LevelSetMarkers.hpp"
#include "include/Mollifier.hpp"
#include "include/Psi.hpp"
#include "include/Reinit.hpp"
#include "include/RungeKutta.hpp"

void FlagFinestMeshLevel(MultiLevelMesh &mlMsh, const double &r,
                         const std::vector<double> &xc);
void FlagFinestMeshLevel(MultiLevelMesh &mlMsh, MyVector<unsigned> &XIel);
void Init(MultiLevelSolution &mlSol, const std::string &name,
          const PsiBall &psi2D);
void GetCutElementPoints(MultiLevelSolution &mlSol, const std::string &name,
                         std::vector<MyVector<double>> &X,
                         MyVector<unsigned> &Xiel);
static void WritePointsVTK(const std::string &filename,
                           const std::vector<MyVector<double>> &X);
void Shift(std::vector<MyVector<double>> &X, const std::vector<double> &dx);

void ProjectSolution(MultiLevelSolution &mlSol0, MultiLevelSolution &mlSol1,
                     BBoxToIel &bbox, RungeKutta &rk);
void TestProjections(LevelMarkers &l0, std::vector<LevelMarkers> &lX);

int main(int argc, char **argv) {

  // Initialize PETSc/MPI
  FemusInit mpinit(argc, argv, MPI_COMM_WORLD);

  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if (nprocs == 1)
    ProfilerStart("profiling.prof");

  MultiLevelMesh mlMsh0;

  const double scalingFactor = 1.0;
  const unsigned numberOfUniformLevels = 1u;
  const unsigned numberOfSelectiveLevels = 3u;

  std::string meshName = "./input/tri.neu";

  // Load coarse mesh and build uniform refinement levels
  mlMsh0.ReadCoarseMesh(meshName.c_str(), "seventh", scalingFactor);
  mlMsh0.RefineMesh(numberOfUniformLevels, numberOfUniformLevels, nullptr);

  // Parameters for selective AMR (ball centered at xc with radius r)
  const double r = 0.125;
  std::vector<double> xc = {0.0, 0.25};
  double eps = 0.25; //(dim == 2) ? 1. / pow(2, std::max(levelN - 7u, 1u))

  // Iteratively flag and create new AMR levels
  for (unsigned k = 0; k < numberOfSelectiveLevels; ++k) {
    FlagFinestMeshLevel(mlMsh0, r, xc);
    mlMsh0.AddAMRMeshLevel(
        false); // false -> it does not re-evaluate the AMR flag vector
  }

  mlMsh0.PrintInfo();
  BBoxToIel bbox(mlMsh0, 0, 3);

  // Define solution on the multilevel mesh
  MultiLevelSolution mlSol0(&mlMsh0);
  mlSol0.AddSolution("Psi", LAGRANGE, SECOND);
  mlSol0.Initialize("All");

  PsiBall psi2D(xc, r, eps);
  Init(mlSol0, "Psi", psi2D);

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
  LevelSetMarkers markers("Psi", 2);

  // Load coarse mesh and build uniform refinement levels
  mlmsh1->ReadCoarseMesh(meshName.c_str(), "seventh", scalingFactor);
  mlmsh1->RefineMesh(numberOfUniformLevels, numberOfUniformLevels, nullptr);

  // RungeKutta::VelKind velocityType = RungeKutta::VelKind::Rotation;
  RungeKutta::VelKind velocityType = RungeKutta::VelKind::Vortex;
  // RungeKutta::VelKind velocityType = RungeKutta::VelKind::Rotation;

  double period =
      (velocityType == RungeKutta::VelKind::Vortex) ? 2 : 2.0 * M_PI;
  unsigned nSteps = 320;
  double dt = period / nSteps;

  for (unsigned t = 1; t <= 0 + 1 * nSteps; t++) {

    double time = t * dt;

    bbox.SetMesh(mlmsh0->GetLevel(0));

    RungeKutta rk(time, dt, period, velocityType);

    unsigned nLevels = numberOfUniformLevels + numberOfSelectiveLevels;
    std::vector<MyVector<double>> X0;
    MyVector<unsigned> X0Iel;
    // GetCutElementPoints(*mlsol0, "Psi", X0, X0Iel);
    markers.GetCutElementPoints(*mlsol0, X0, X0Iel);

    if (t % 10 == 0) {
      Reinit reinit("Psi", *mlsol0, eps);

      reinit.farFieldReinit(X0);
      reinit.interfaceFieldReinit(bbox);
      reinit.updateSolution();
    }

    if (t == 1)
      WritePointsVTK("./output/points.0.vtk", X0);

    rk.rkForward(X0);

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
    mlsol1->Initialize("All");

    ProjectSolution(*mlsol0, *mlsol1, bbox, rk);

    // Export solution to VTK (selected levels)

    std::swap(mlsol0, mlsol1);
    std::swap(mlmsh0, mlmsh1);
    mlsol1->clear();
    mlmsh1->resize(numberOfUniformLevels);

    VTKWriter vtkIO1(mlsol0);
    if (t % 10 == 0)
      vtkIO1.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted,
                   t / 10);
  }

  if (nprocs == 1)
    ProfilerStop();
  return 0;
}

void Shift(std::vector<MyVector<double>> &X, const std::vector<double> &dx) {
  for (unsigned k = 0; k < X.size(); k++) {
    for (unsigned i = X[k].begin(); i < X[k].end(); i++) {
      X[k][i] += dx[k];
    }
  }
}

void FlagFinestMeshLevel(MultiLevelMesh &mlMsh, const double &r,
                         const std::vector<double> &xc) {

  const unsigned level = mlMsh.GetNumberOfLevels() - 1;
  Mesh &msh = *mlMsh.GetLevel(level);
  const unsigned iproc = msh.processor_id();
  const unsigned dim = msh.GetDimension();
  const unsigned xType = 2;
  const unsigned amrIndex = msh.GetAmrIndex();

  const double r2 = r * r;

  auto &xv = msh._topology->_Sol;
  auto *amrFlag = msh._topology->_Sol[amrIndex];

  amrFlag->zero();

  unsigned offset = msh._elementOffset[iproc];
  unsigned offsetp1 = msh._elementOffset[iproc + 1];

  for (unsigned iel = offset; iel < offsetp1; ++iel) {

    const unsigned nDof = msh.GetElementDofNumber(iel, xType);

    const unsigned xDof0 = msh.GetSolutionDof(0, iel, xType);
    double d2_0 = r2;
    for (unsigned k = 0; k < dim; ++k) {
      const double d = (*xv[k])(xDof0)-xc[k];
      d2_0 -= d * d;
    }
    const int sign0 = (d2_0 > 0.) ? 1 : -1;

    bool signChange = false;
    for (unsigned i = 1; i < nDof; ++i) {
      const unsigned xDof = msh.GetSolutionDof(i, iel, xType);
      double d2 = r2;

      for (unsigned k = 0; k < dim; ++k) {
        const double d = (*xv[k])(xDof)-xc[k];
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
    if ((*amrFlag)(iel) == 2) {
      for (unsigned j = 0; j < msh.el->GetElementNearElementSize(iel, 1); j++) {
        unsigned jel = msh.el->GetElementNearElement(iel, j);
        if (offset <= jel && jel < offsetp1) {
          if ((*amrFlag)(jel) < 0.5)
            amrFlag->set(jel, 1); // this is on spot since jel belongs to iproc
        } else {
          amrFlag->add(
              jel, 1); // this is buffered since jel does not belong to iproc
        }
      }
    }
  }
  amrFlag->close();

  unsigned localRefined = 0;
  for (unsigned iel = offset; iel < offsetp1; ++iel) {
    if ((*amrFlag)(iel) > 0.5) {
      amrFlag->set(iel, 1);
      ++localRefined;
    }
  }
  amrFlag->close();

  unsigned globalRefined = 0;
  MPI_Allreduce(&localRefined, &globalRefined, 1, MPI_UNSIGNED, MPI_SUM,
                MPI_COMM_WORLD);

  msh.el->SetRefinedElementNumber(globalRefined);
}

void FlagFinestMeshLevel(MultiLevelMesh &mlMsh, MyVector<unsigned> &XIel) {

  const unsigned level = mlMsh.GetNumberOfLevels() - 1;
  Mesh &msh = *mlMsh.GetLevel(level);
  const unsigned iproc = msh.processor_id();
  const unsigned dim = msh.GetDimension();
  const unsigned amrIndex = msh.GetAmrIndex();

  auto *amrFlag = msh._topology->_Sol[amrIndex];

  amrFlag->zero();

  unsigned offset = msh._elementOffset[iproc];
  unsigned offsetp1 = msh._elementOffset[iproc + 1];

  for (unsigned i = XIel.begin(); i < XIel.end(); i++) {
    unsigned iel = XIel[i];
    amrFlag->set(iel, 2);
  }
  amrFlag->close();

  for (unsigned iel = offset; iel < offsetp1; ++iel) {
    if ((*amrFlag)(iel) == 2) {
      for (unsigned j = 0; j < msh.el->GetElementNearElementSize(iel, 1); j++) {
        unsigned jel = msh.el->GetElementNearElement(iel, j);
        if (offset <= jel && jel < offsetp1) {
          if ((*amrFlag)(jel) < 0.5)
            amrFlag->set(jel, 1); // this is on spot since jel belongs to iproc
        } else {
          amrFlag->add(
              jel, 1); // this is buffered since jel does not belong to iproc
        }
      }
    }
  }
  amrFlag->close();

  unsigned localRefined = 0;
  for (unsigned iel = offset; iel < offsetp1; ++iel) {
    if ((*amrFlag)(iel) > 0.5) {
      amrFlag->set(iel, 1);
      ++localRefined;
    }
  }
  amrFlag->close();

  unsigned globalRefined = 0;
  MPI_Allreduce(&localRefined, &globalRefined, 1, MPI_UNSIGNED, MPI_SUM,
                MPI_COMM_WORLD);

  msh.el->SetRefinedElementNumber(globalRefined);
}

void Init(MultiLevelSolution &mlSol, const std::string &name,
          const PsiBall &psi2D) {

  MultiLevelMesh &mlMsh = *mlSol.GetMultilevelMesh();
  const unsigned level = mlMsh.GetNumberOfLevels() - 1u;
  Mesh &msh = *mlMsh.GetLevel(level);
  Solution &sol = *mlSol.GetLevel(level);
  const unsigned iproc = msh.processor_id();
  const unsigned dim = msh.GetDimension();

  unsigned solIndex = mlSol.GetIndex(name.c_str());
  unsigned solType = mlSol.GetSolutionType(name.c_str());

  const unsigned xType = 2u; // coordinate field type

  auto &xv = msh._topology->_Sol;
  std::vector<double> x(dim);

  const unsigned offset = msh._elementOffset[iproc];
  const unsigned offsetp1 = msh._elementOffset[iproc + 1];

  auto &solVec = sol._Sol[solIndex];

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

void GetCutElementPoints(MultiLevelSolution &mlSol, const std::string &name,
                         std::vector<MyVector<double>> &X,
                         MyVector<unsigned> &Xiel) {

  MultiLevelMesh &mlMsh = *mlSol.GetMultilevelMesh();
  const unsigned level = mlMsh.GetNumberOfLevels() - 1u;
  Mesh &msh = *mlMsh.GetLevel(level);
  Solution &sol = *mlSol.GetLevel(level);
  const unsigned iproc = msh.processor_id();
  const unsigned dim = msh.GetDimension();

  const unsigned solIndex = mlSol.GetIndex(name.c_str());
  const unsigned solType = mlSol.GetSolutionType(name.c_str());
  const unsigned xType = 2u;

  const double c1 = 2. / 3., c2 = 1. / 3.;

  auto &xv = msh._topology->_Sol;
  auto &solVec = sol._Sol[solIndex];

  const unsigned offset = msh._elementOffset[iproc];
  const unsigned offsetp1 = msh._elementOffset[iproc + 1];

  unsigned maxPtsPerEl = 1u;
  for (unsigned d = 0; d < dim; ++d)
    maxPtsPerEl *= 2u;
  ++maxPtsPerEl;

  std::vector<std::vector<double>> x(dim);
  for (unsigned k = 0; k < dim; ++k)
    x[k].reserve(maxPtsPerEl);

  std::vector<double> phi;
  std::vector<double> gradPhi(dim);

  std::vector<std::vector<double>> Y(dim);
  std::vector<unsigned> Yiel;

  for (unsigned k = 0; k < dim; ++k) {
    Y[k].reserve((offsetp1 - offset) * maxPtsPerEl);
  }
  Yiel.reserve((offsetp1 - offset) * maxPtsPerEl);

  // Loop over local elements and collect interior points associated with cut
  // elements
  for (unsigned iel = offset; iel < offsetp1; ++iel) {

    const unsigned nDof = msh.GetElementDofNumber(iel, solType);

    const unsigned solDof0 = msh.GetSolutionDof(0, iel, solType);
    const double val0 = (*solVec)(solDof0);
    const int sign0 = (val0 > 0.) - (val0 < 0.);

    for (unsigned i = 1; i < nDof; ++i) {
      const unsigned solDofi = msh.GetSolutionDof(i, iel, solType);
      const double vali = (*solVec)(solDofi);
      const int signi = (vali > 0.) - (vali < 0.);

      if (signi != sign0) {

        const unsigned nDof0 = msh.GetElementDofNumber(iel, 0);
        phi.resize(nDof0 + 1u);
        for (unsigned k = 0; k < dim; ++k) {
          x[k].resize(nDof0 + 1u);
        }
        for (unsigned j = 0; j < nDof0; ++j) {
          const unsigned solDofj = msh.GetSolutionDof(j, iel, solType);
          phi[j] = (*solVec)(solDofj);
          const unsigned xDof = msh.GetSolutionDof(j, iel, xType);
          for (unsigned k = 0; k < dim; ++k) {
            x[k][j] = (*xv[k])(xDof);
          }
        }
        const unsigned nDof2 = msh.GetElementDofNumber(iel, 2);
        const unsigned xDofc = msh.GetSolutionDof(nDof2 - 1u, iel, xType);
        if (solType == xType) {
          phi[nDof0] = (*solVec)(xDofc);

        } else {
          phi[nDof0] = 0.;
          for (unsigned j = 0; j < nDof0; ++j) {
            phi[nDof0] += phi[j];
          }
          phi[nDof0] /= nDof0;
        }
        for (unsigned k = 0; k < dim; ++k) {
          x[k][nDof0] = (*xv[k])(xDofc);
        }

        computeElementGradientFromLocalData(x, phi, gradPhi);

        // shift points
        for (unsigned j = 0; j < nDof0; ++j) {
          phi[j] = c1 * phi[j] + c2 * phi[nDof0];
          for (unsigned k = 0; k < dim; ++k) {
            x[k][j] = c1 * x[k][j] + c2 * x[k][nDof0];
          }
        }

        double gradNorm2 = 0.;
        for (unsigned k = 0; k < dim; ++k) {
          gradNorm2 += gradPhi[k] * gradPhi[k];
        }

        if (gradNorm2 < 1.e-20) {
          for (unsigned j = 0; j <= nDof0; ++j) {
            for (unsigned k = 0; k < dim; ++k) {
              Y[k].push_back(x[k][j]);
            }
            Yiel.push_back(iel);
          }
          break; // use current points
        }

        const double invGradNorm2 = 1. / gradNorm2;
        for (unsigned j = 0; j <= nDof0; ++j) {
          for (unsigned k = 0; k < dim; ++k) {
            Y[k].push_back(x[k][j] - phi[j] * gradPhi[k] * invGradNorm2);
          }
          Yiel.push_back(iel);
        }
        break;
      }
    }
  }

  X.resize(dim);
  for (unsigned k = 0; k < dim; ++k) {
    X[k].buildFromLocal(Y[k]);
  }
  Xiel.buildFromLocal(Yiel);
}

static void WritePointsVTK(const std::string &filename,
                           const std::vector<MyVector<double>> &X) {

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
      throw std::runtime_error(
          "writePointsVTK: inconsistent sizes across components");
    }
  }

  int iproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);

  // Only rank 0 writes
  if (iproc != 0)
    return;

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

void TestProjections(LevelMarkers &l0, std::vector<LevelMarkers> &lX) {

  std::vector<MyVector<double>> field0Copy = l0.GetFields();
  std::vector<MyVector<double>> &field0 =
      l0.GetFields(); // assume we have a valid field in lX

  const unsigned nFields = field0.size();
  if (nFields > 0) {

    std::vector<std::vector<double>> Wfield_r, Wfield_s;

    const bool backward = true;
    const bool forward = !backward;
    const bool checkIfInsideProcDomain = true;

    l0.RebuildLocalFromField(Wfield_s, nFields, forward);
    l0.SendLocalField(Wfield_s, Wfield_r);
    lX[0].RebuildFieldFromLocal(Wfield_r, nFields, forward);

    for (unsigned l = 0; l < lX.size() - 1u; l++) {
      lX[l].RebuildLocalFromField(Wfield_s, nFields, forward);
      lX[l].SendLocalField(Wfield_s, Wfield_r);
      lX[l + 1].RebuildFieldFromLocal(Wfield_r, nFields, forward);

      lX[l + 1].RebuildLocalFromField(Wfield_r, nFields, backward);
      lX[l].SendLocalField(Wfield_r, Wfield_s);
      lX[l].RebuildFieldFromLocal(Wfield_s, nFields, backward);
    }

    lX[0].RebuildLocalFromField(Wfield_r, nFields, backward); // Probably okay
    l0.SendLocalField(Wfield_r, Wfield_s);
    l0.RebuildFieldFromLocal(Wfield_s, nFields, backward);

    const std::vector<bool> &ZisInside = l0.GetPointInsideDomain();

    unsigned offset0 = field0[0].begin();
    unsigned offset1 = field0[0].end();

    for (unsigned i = offset0; i < offset1; i++) {
      if (ZisInside[i - offset0] == false) {
        for (unsigned k = 0; k < nFields; k++) {
          field0[k][i] = 0.; // TODO add BC
        }
      }
    }

    for (unsigned i = offset0; i < offset1; i++) {
      if (ZisInside[i - offset0] == true) {
        for (unsigned k = 0; k < nFields; k++) {
          if (fabs(field0Copy[k][i] - field0[k][i]) > 1.0e-12)
            std::cerr << "error ";
        }
      }
    }
  }
}

void GetAllSolutionPoints(MultiLevelSolution &mlSol, const std::string &name,
                          std::vector<MyVector<double>> &X) {

  MultiLevelMesh &mlMsh = *mlSol.GetMultilevelMesh();
  const unsigned level = mlMsh.GetNumberOfLevels() - 1u;
  Mesh &msh = *mlMsh.GetLevel(level);
  const unsigned dim = msh.GetDimension();

  const unsigned solType = mlSol.GetSolutionType(name.c_str());
  const unsigned xType = 2u;

  if (solType > xType) {
    throw std::runtime_error("GetAllSolutionPoints: coordinate FE space is too "
                             "low-order for solType");
  }

  const unsigned iproc = msh.processor_id();

  auto &xv = msh._topology->_Sol;

  const unsigned solOffset = msh._dofOffset[solType][iproc];
  const unsigned solOffsetp1 = msh._dofOffset[solType][iproc + 1];

  std::vector<std::vector<double>> Xloc(dim);
  for (unsigned k = 0; k < dim; k++) {
    Xloc[k].resize(solOffsetp1 - solOffset);
  }

  const unsigned elOffset = msh._elementOffset[iproc];
  const unsigned elOffsetp1 = msh._elementOffset[iproc + 1];

  for (unsigned iel = elOffset; iel < elOffsetp1; ++iel) {
    const unsigned nDofSol = msh.GetElementDofNumber(iel, solType);
    for (unsigned i = 0; i < nDofSol; ++i) {
      const unsigned sdof = msh.GetSolutionDof(i, iel, solType);
      if (solOffset <= sdof && sdof < solOffsetp1) {
        const unsigned xdof = msh.GetSolutionDof(i, iel, xType);
        for (unsigned k = 0; k < dim; k++) {
          Xloc[k][sdof - solOffset] = (*xv[k])(xdof);
        }
      }
    }
  }

  X.resize(dim);
  for (unsigned k = 0; k < dim; ++k) {
    X[k].buildFromLocal(Xloc[k]);
  }
}

void ProjectSolution(MultiLevelSolution &mlSol0, MultiLevelSolution &mlSol1,
                     BBoxToIel &bbox, RungeKutta &rk) {

  const std::string name = "Psi";

  MultiLevelMesh &mlMsh0 = *mlSol0.GetMultilevelMesh();

  const unsigned nLevels = mlMsh0.GetNumberOfLevels();

  // Extract all Psi grid points on the finest level of mlSol1
  std::vector<MyVector<double>> X1;
  GetAllSolutionPoints(mlSol1, name, X1);

  // Advect the points backward in time
  rk.rkBackward(X1);

  // const unsigned bbxLevel = bbox.GetLevel();
  const unsigned bboxLevels = nLevels - bbox.GetLevel();
  // Build l0 and lX forward projection using the advected points
  LevelMarkers l0;
  std::vector<LevelMarkers> lX(bboxLevels);

  bbox.GetInverseMappingOnCoarseLevel(X1, l0, lX[0]);

  for (unsigned k = 1; k < bboxLevels; ++k) {
    bbox.Project(mlMsh0, lX[k - 1], lX[k]);
  }

  // Evaluate Psi field on finest level of mlSol0 using top-level iel and xi
  const unsigned level0 = nLevels - 1u;

  Mesh &msh0 = *mlMsh0.GetLevel(level0);
  Solution &sol0 = *mlSol0.GetLevel(level0);

  const unsigned solIndex0 = mlSol0.GetIndex(name.c_str());
  const unsigned solType0 = mlSol0.GetSolutionType(name.c_str());

  auto &solVec0 = sol0._Sol[solIndex0];

  LevelMarkers &lTop = lX.back();

  std::vector<MyVector<double>> &Xi = lTop.GetLocalCoordinates();
  MyVector<unsigned> &Iel = lTop.GetElements();

  const unsigned dim = Xi.size();

  std::vector<double> psiLocal;
  psiLocal.resize(Iel.end() - Iel.begin(), 0.0);

  std::vector<double> xi(dim);
  std::vector<double> phi;

  for (unsigned ip = Iel.begin(); ip < Iel.end(); ++ip) {

    const unsigned iel = Iel[ip];
    short unsigned ielType = msh0.GetElementType(iel);

    for (unsigned k = 0; k < dim; ++k) {
      xi[k] = Xi[k][ip];
    }

    const unsigned nDof = msh0.GetElementDofNumber(iel, solType0);

    phi.resize(nDof);

    msh0._finiteElement[ielType][solType0]->GetPhi(phi, xi);

    double value = 0.0;
    for (unsigned j = 0; j < nDof; ++j) {
      const unsigned solDof = msh0.GetSolutionDof(j, iel, solType0);
      value += phi[j] * (*solVec0)(solDof);
    }

    psiLocal[ip - Iel.begin()] = value;
  }

  lTop.GetFields().resize(1);
  lTop.GetFields()[0].buildFromLocal(psiLocal);

  // Project Psi backward through the marker hierarchy
  std::vector<std::vector<double>> Wfield_r;
  std::vector<std::vector<double>> Wfield_s;

  const unsigned nFields = 1u;
  const bool backward = true;

  for (int l = static_cast<int>(bboxLevels) - 1; l >= 1; --l) {
    lX[l].RebuildLocalFromField(Wfield_r, nFields, backward);
    lX[l - 1].SendLocalField(Wfield_r, Wfield_s);
    lX[l - 1].RebuildFieldFromLocal(Wfield_s, nFields, backward);
  }

  lX[0].RebuildLocalFromField(Wfield_r, nFields, backward);
  l0.SendLocalField(Wfield_r, Wfield_s);
  l0.RebuildFieldFromLocal(Wfield_s, nFields, backward);

  // Update Psi on mlSol1 with the backward-projected values

  MultiLevelMesh &mlMsh1 = *mlSol1.GetMultilevelMesh();
  const unsigned level1 = mlMsh1.GetNumberOfLevels() - 1u;

  Solution &sol1 = *mlSol1.GetLevel(level1);
  const unsigned solIndex1 = mlSol1.GetIndex(name.c_str());

  auto &solVec1 = sol1._Sol[solIndex1];

  const MyVector<double> &psiProjected = l0.GetFields()[0];
  const std::vector<bool> &isInsideDomain = l0.GetPointInsideDomain();

  solVec1->zero();
  unsigned offset = psiProjected.begin();
  for (unsigned i = psiProjected.begin(); i < psiProjected.end(); ++i) {
    if (isInsideDomain[i - offset]) {
      solVec1->set(i, psiProjected[i]);
    } else { // TODO add boundarycondition for psi
      solVec1->set(i, -1.);
    }
  }

  solVec1->close();
}
