
#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "VTKWriter.hpp"
#include "MyMatrix.hpp"


using namespace femus;

class LevelMarkers {
  public:
    LevelMarkers() = default;

    std::vector<MyVector<double>>& GetFields() {
      return _field;
    }

    const std::vector<MyVector<double>>& GetFields() const {
      return _field;
    }

    std::vector<MyVector<double>>& GetLocalCoordinates() {
      return _Xi;
    }

    const std::vector<MyVector<double>>& GetLocalCoordinates() const {
      return _Xi;
    }

    MyVector<unsigned>& GetElements() {
      return _XIel;
    }

    const MyVector<unsigned>& GetElements() const {
      return _XIel;
    }


    unsigned& GetLevel() {
      return _level;
    }

    unsigned GetLevel() const {
      return _level;
    }

    std::vector<std::vector<unsigned>> & GetMap_s() {
      return _map_s;
    }

    const std::vector<std::vector<unsigned>> & GetMap_s() const {
      return _map_s;
    }

    std::vector<std::vector<unsigned>> & GetMap_r() {
      return _map_r;
    }

    const std::vector<std::vector<unsigned>> & GetMap_r() const {
      return _map_r;
    }



    void SetLevel(unsigned level) {
      _level = level;
    }

    // void RebuildLocalReceivers(std::vector<std::vector<double>> &W_r,
    //                            std::vector<std::vector<double>> &Wi_r,
    //                            std::vector<std::vector<unsigned>> &WIel_r) {
    //   const unsigned nprocs = _offset.size() - 1u;
    //   const unsigned dim    = _Xi.size();
    //   const unsigned nFields  = _field.size();
    //
    //   assert(_X.size() == _Xi.size());
    //   assert(_offset.back() == _XIel.end() - _XIel.begin());
    //
    //   W_r.resize(nprocs);
    //   Wi_r.resize(nprocs);
    //   WIel_r.resize(nprocs);
    //
    //   for (unsigned p = 0; p < nprocs; ++p) {
    //     assert(_offset[p] <= _offset[p + 1]);
    //     const unsigned nloc = _offset[p + 1] - _offset[p];
    //     WIel_r[p].resize(nloc);
    //     W_r[p].resize(nloc * dim);
    //     Wi_r[p].resize(nloc * dim);
    //   }
    //
    //   unsigned kproc = 0u;
    //
    //   for (unsigned i = _XIel.begin(); i < _XIel.end(); ++i) {
    //     const unsigned i0 = i - _XIel.begin();
    //
    //     while (kproc + 1 < nprocs && i0 >= _offset[kproc + 1]) {
    //       ++kproc;
    //     }
    //
    //     const unsigned jloc = i0 - _offset[kproc];
    //     const unsigned base = jloc * dim;
    //
    //     // point-major packing: [x0,y0,(z0), x1,y1,(z1), ...]
    //     for (unsigned k = 0; k < dim; ++k) {
    //       W_r[kproc][base + k]  = _field[k][i];
    //       Wi_r[kproc][base + k] = _Xi[k][i];
    //     }
    //
    //     WIel_r[kproc][jloc] = _XIel[i];
    //   }
    // }


    void RebuildReceiverFromField(std::vector<std::vector<double>> &field_r, const unsigned nFields) {

      const unsigned nprocs = _map_r.size();
      field_r.resize(nprocs);

      assert(nFields  == _field.size());
      if(nFields > 0) {

        const unsigned offset1 = (nFields == 0) ? 0 : _field[0].begin();

        for (unsigned kproc = 0; kproc < nprocs; ++kproc) {
          field_r[kproc].resize(_map_r[kproc].size() * nFields);


          for (unsigned i = 0; i < _map_r[kproc].size(); ++i) {
            for (unsigned j = 0; j < nFields; ++j) {
              field_r[kproc][i * nFields + j] = _field[j][offset1 + _map_r[kproc][i]];
            }
          }
        }
      }
    }



    void RebuildFieldFromSender(const std::vector<std::vector<double>> &field_s, const unsigned nFields) {
      unsigned nprocs = field_s.size();
      for (unsigned p = 0; p < nprocs; ++p) {
        assert(field_s[p].size() % nFields == 0);
        assert(_map[p].size() == field_s[p].size());
      }
      assert(nFields == _field.size());

      unsigned sizeAll = 0u;
      for(unsigned p = 0; p < nprocs; p++) {
        sizeAll += field_s[p].size();
      }
      std::vector<std::vector<double>> locField(nFields);
      for(unsigned k = 0; k < nFields; k++) locField[k].resize(sizeAll / nFields);
      for(unsigned p = 0; p < nprocs; p++) {
        for(unsigned i = 0; i < field_s[p].size(); i++ ) {
          locField[i % nFields][_map_s[p][i / nFields]] = field_s[p][i];
        }
      }
      _field.resize(nFields);
      for(unsigned k = 0; k < nFields; k++) _field[k].buildFromLocal(locField[k]);
    }

  private:
    std::vector<MyVector<double>> _field;
    std::vector<MyVector<double>> _Xi;
    MyVector<unsigned> _XIel;

    std::vector<std::vector<unsigned>> _map_s;
    std::vector<std::vector<unsigned>> _map_r;

    unsigned _level = 0;
};



#include "include/Mollifier.hpp"
#include "include/Psi.hpp"
#include "include/GradientApproximation.hpp"
#include "include/BBoxToIel.hpp"

void FlagFinestMeshLevel(MultiLevelMesh &mlMsh, const double &r, const std::vector<double> &xc);
void FlagFinestMeshLevel(MultiLevelMesh &mlMsh, MyVector<unsigned> &XIel);
void Init(MultiLevelSolution &mlSol, const std::string &name, const PsiBall &psi2D);
void GetCutElementPoints(MultiLevelSolution &mlSol, const std::string &name, std::vector<MyVector<double>> &X, MyVector<unsigned> &Xiel);
static void WritePointsVTK(const std::string& filename, const std::vector<MyVector<double>>& X);
void Shift( std::vector<MyVector<double>> &X, const std::vector<double>&dx);


int main(int argc, char** argv) {

  // Initialize PETSc/MPI
  FemusInit mpinit(argc, argv, MPI_COMM_WORLD);

  MultiLevelMesh mlMsh0;

  const double scalingFactor = 1.0;
  const unsigned numberOfUniformLevels   = 3u;
  const unsigned numberOfSelectiveLevels = 3u;

  // Load coarse mesh and build uniform refinement levels
  mlMsh0.ReadCoarseMesh("./input/square.neu", "seventh", scalingFactor);
  mlMsh0.RefineMesh(numberOfUniformLevels, numberOfUniformLevels, nullptr);

  // Parameters for selective AMR (ball centered at xc with radius r)
  const double r = 0.125;
  std::vector<double> xc = {0.0, 0.25};
  double eps = 0.01;//(dim == 2) ? 1. / pow(2, std::max(levelN - 7u, 1u))

  // Iteratively flag and create new AMR levels
  for (unsigned k = 0; k < numberOfSelectiveLevels; ++k) {
    FlagFinestMeshLevel(mlMsh0, r, xc);
    mlMsh0.AddAMRMeshLevel(false); // false -> it does not re-evaluate the AMR flag vector
  }

  mlMsh0.PrintInfo();
  BBoxToIel bbox0(mlMsh0, 3, 0);


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


  std::vector<MyVector<double>> X0;
  MyVector<unsigned> X0Iel;
  GetCutElementPoints(mlSol0, "Psi", X0, X0Iel);
  WritePointsVTK("./output/points.0.vtk", X0);
  Shift(X0, {0.1, 0});
  WritePointsVTK("./output/points.1.vtk", X0);

  MultiLevelMesh mlMsh1;


  // Load coarse mesh and build uniform refinement levels
  mlMsh1.ReadCoarseMesh("./input/square.neu", "seventh", scalingFactor);
  mlMsh1.RefineMesh(numberOfUniformLevels, numberOfUniformLevels, nullptr);

  unsigned nLevels = numberOfUniformLevels + numberOfSelectiveLevels;

  BBoxToIel bbox1(mlMsh1, 3, 0);


  std::vector<MyVector<double>> field = X0;
  field.resize(1);
  //field.clear();

  std::vector<LevelMarkers> lX(nLevels);
  bbox1.GetInverseMappingOnCoarseLevel(X0, field, lX[0]);

  for(unsigned k = 1; k < numberOfUniformLevels; k++ ) {
    bbox1.Project(lX[k - 1], lX[k]);
  }

  for (unsigned k = numberOfUniformLevels; k < nLevels; ++k) {
    FlagFinestMeshLevel(mlMsh1, lX[k - 1].GetElements());
    mlMsh1.AddAMRMeshLevel(false); // false -> it does not re-evaluate the AMR flag vector
    bbox1.Project(lX[k - 1], lX[k]);
  }

  MultiLevelSolution mlSol1(&mlMsh1);
  mlSol1.AddSolution("Psi", LAGRANGE, SECOND);
  mlSol1.Initialize("All");

  // Export solution to VTK (selected levels)

  VTKWriter vtkIO1(&mlSol1);
  //vtkIO1.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);
  vtkIO1.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 1);

  return 0;
}


void Shift( std::vector<MyVector<double>> &X, const std::vector<double>&dx) {
  for(unsigned k = 0; k < X.size(); k++) {
    for(unsigned i = X[k].begin(); i < X[k].end(); i++) {
      X[k][i] += dx[k];
    }
  }
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


void FlagFinestMeshLevel(MultiLevelMesh &mlMsh, MyVector<unsigned> &XIel) {

  const unsigned level    = mlMsh.GetNumberOfLevels() - 1;
  Mesh &msh               = *mlMsh.GetLevel(level);
  const unsigned iproc    = msh.processor_id();
  const unsigned dim      = msh.GetDimension();
  const unsigned amrIndex = msh.GetAmrIndex();

  auto* amrFlag = msh._topology->_Sol[amrIndex];

  amrFlag->zero();

  unsigned offset = msh._elementOffset[iproc];
  unsigned offsetp1 = msh._elementOffset[iproc + 1];

  for(unsigned i = XIel.begin(); i < XIel.end(); i++) {
    unsigned iel = XIel[i];
    amrFlag->set(iel, 2);
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

