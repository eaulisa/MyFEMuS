#pragma once

void RungeKutta4(std::vector<MyVector<double>> &X,
                 MultiLevelSolution & mlSol,
                 BBoxToIel & bbox,
                 const std::vector<std::string> &velName,
                 const double dt);

void rkStep(MultiLevelSolution & mlSol,
            BBoxToIel & bbox,
            const std::vector<MyVector<double>> &X,
            std::vector<std::vector<MyVector<double>>> &K,
            const unsigned rkStep,
            const std::vector<std::string> &velName,
            const double dt,
            const double c,
            const std::vector<double> &a);

void InterpolateSolution(LevelMarkers &l0,
                         MultiLevelSolution &mlSol0,
                         BBoxToIel &bbox,
                         std::vector<MyVector<double>> &X,
                         const std::vector< std::string >solName,
                         const double c = 1.
                        );



void ProjectSolution(MultiLevelSolution &mlSol0 /* target */, MultiLevelSolution &mlSol1 /* source */,
                     BBoxToIel &bbox, RungeKutta &rk,
                     const std::vector<std::string> solName,
                     const std::vector<std::string> vName = {},
                     const Boundary& bd_inflow = Boundary(),
                     const double dt = 0.,
                     const double time = 0.,
                     const double period = 0.);




void RestrictPWDCField(MultiLevelSolution &mlSol,
                       const std::string &CName,
                       const unsigned level0,
                       const unsigned level1) {

  MultiLevelMesh &mlMsh = *mlSol.GetMultilevelMesh();

  const unsigned solIndex = mlSol.GetIndex(CName.c_str());
  const unsigned solType = mlSol.GetSolutionType(CName.c_str());

  if(solType != 3) {
    std::cout << "Error! The C Field is not PWC\n" << std::endl;
    abort();
  }

  for(int l = level1; l > level0; l--) {

    Mesh &msh_l = *mlMsh.GetLevel(l);
    Mesh &msh_lm1 = *mlMsh.GetLevel(l - 1);

    const unsigned iproc = msh_l.processor_id();
    const unsigned dim = msh_l.GetDimension();
    const unsigned maxNumberOfChildren = 1u << dim;

    auto &solC_l = (mlSol.GetSolutionLevel(l))->_Sol[solIndex];
    auto &solC_lm1 = (mlSol.GetSolutionLevel(l - 1))->_Sol[solIndex];

    NumericVector *father = NumericVector::build().release();
    father->init(*solC_l);
    father->zero();

    for(unsigned iel_lm1 = msh_lm1._elementOffset[iproc];
        iel_lm1 < msh_lm1._elementOffset[iproc + 1];
        iel_lm1++) {

      const unsigned numberOfChildren =
          msh_lm1.GetRefinedElementIndex(iel_lm1) ?
          maxNumberOfChildren : 1;

      for(unsigned j = 0; j < numberOfChildren; j++) {
        const unsigned iel_l =
            msh_lm1.el->GetChildElement(iel_lm1, j);

        father->set(iel_l, iel_lm1);
      }
    }

    father->close();

    solC_lm1->zero();

    for(unsigned iel_l = msh_l._elementOffset[iproc];
        iel_l < msh_l._elementOffset[iproc + 1];
        iel_l++) {

      const unsigned iel_lm1 =
          static_cast<unsigned>((*father)(iel_l));

      solC_lm1->add(iel_lm1, (*solC_l)(iel_l));
    }

    solC_lm1->close();

    const double tol = 1.e-10;

    for(unsigned iel_lm1 = msh_lm1._elementOffset[iproc];
        iel_lm1 < msh_lm1._elementOffset[iproc + 1];
        iel_lm1++) {

      if(msh_lm1.GetRefinedElementIndex(iel_lm1)) {

        const double value = (*solC_lm1)(iel_lm1);

        if(value <= tol) {
          solC_lm1->set(iel_lm1, 0.);
        }
        else if(value >= maxNumberOfChildren - tol) {
          solC_lm1->set(iel_lm1, 1.);
        }
        else {
          solC_lm1->set(iel_lm1, 0.5);
        }
      }
    }

    solC_lm1->close();

    delete father;
  }
}


inline std::vector<double>
Velocity(std::vector<double> &xp, const double time, double period) noexcept {
  double u = 0.0, v = 0.0;

  switch(velocityType) {
    case RungeKutta::VelKind::Vortex: {
      const double T = period;
      const double x = xp[0] + 0.5;
      const double y = xp[1] + 0.5;

      const double sx = std::sin(M_PI * x);
      const double cx = std::cos(M_PI * x);
      const double sy = std::sin(M_PI * y);
      const double cy = std::cos(M_PI * y);
      const double cosT = std::cos(M_PI * time / T);

      u = -2.0 * sx * sx * sy * cy * cosT;
      v =  2.0 * sx * cx * sy * sy * cosT;

      break;
    }
    case RungeKutta::VelKind::Rotation: {

      u = xp[1];
      v = -xp[0];

      break;
    }
    case RungeKutta::VelKind::Translation: {

      u = 0;
      v = -0.3 * 4 * (0.5 - xp[0]) * (0.5 + xp[0]);
      // v = -0.3;

      break;
    }
    case RungeKutta::VelKind::Zero: {

      u = 0;
      v = 0;


      break;
    }

  }

  return {u, v};
}




void Shift(std::vector<MyVector<double>> &X, const std::vector<double> &dx) {
  for (unsigned k = 0; k < X.size(); k++) {
    for (unsigned i = X[k].begin(); i < X[k].end(); i++) {
      X[k][i] += dx[k];
    }
  }
}

void FlagFinestMeshLevel(MultiLevelMesh & mlMsh, const double & r,
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
    if ((*amrFlag)(iel) == 2) {
      for (unsigned j = 0; j < msh.el->GetElementNearElementSize(iel, 1); j++) {
        unsigned jel = msh.el->GetElementNearElement(iel, j);
        if (offset <= jel && jel < offsetp1) {
          if ((*amrFlag)(jel) < 0.5)
            amrFlag->set(jel, 1); // this is on spot since jel belongs to iproc
        }
        else {
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

void FlagFinestMeshLevel(MultiLevelMesh & mlMsh, MyVector<unsigned> &XIel) {

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
        }
        else {
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

void InitLevelSet(MultiLevelSolution & mlSol, const std::string & name,
                  const PsiBall & psi2D) {

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


void UpdateColorFunction(MultiLevelSolution & mlSol, const std::string & psiName, const std::string & cName) {

  MultiLevelMesh &mlMsh = *mlSol.GetMultilevelMesh();
  const unsigned level = mlMsh.GetNumberOfLevels() - 1u;
  Mesh &msh = *mlMsh.GetLevel(level);
  Solution &sol = *mlSol.GetLevel(level);
  const unsigned iproc = msh.processor_id();
  const unsigned dim = msh.GetDimension();

  unsigned psiIndex = mlSol.GetIndex(psiName.c_str());
  unsigned psiType = mlSol.GetSolutionType(psiName.c_str());

  unsigned cIndex = mlSol.GetIndex(cName.c_str());

  const unsigned offset = msh._elementOffset[iproc];
  const unsigned offsetp1 = msh._elementOffset[iproc + 1];

  auto &psiVec = sol._Sol[psiIndex];
  auto &cVec = sol._Sol[cIndex];

  cVec->zero();

  // Loop over local elements and interpolate psi2D at solution DoFs
  for (unsigned iel = offset; iel < offsetp1; ++iel) {
    const unsigned nDof = msh.GetElementDofNumber(iel, psiType);
    double value0 = (*psiVec)(msh.GetSolutionDof(0, iel, psiType));
    bool signChanged = false;
    for (unsigned i = 1; i < nDof; ++i) {
      double value = (*psiVec)(msh.GetSolutionDof(i, iel, psiType));
      if(value0 * value <= 0.) {
        cVec->set(iel, 0.5);
        signChanged = true;
        break;
      }
    }
    if(!signChanged && value0 > 0.) {
      cVec->set(iel, 1.);
    }
  }
  cVec->close();
}


void InitSol(MultiLevelSolution & mlSol, const std::vector<std::string> &solName, const double time, const double period) {

  MultiLevelMesh &mlMsh = *mlSol.GetMultilevelMesh();
  const unsigned level = mlMsh.GetNumberOfLevels() - 1u;
  Mesh &msh = *mlMsh.GetLevel(level);
  Solution &sol = *mlSol.GetLevel(level);
  const unsigned iproc = msh.processor_id();
  const unsigned dim = msh.GetDimension();

  assert(dim == solName.size());


  std::vector<unsigned> solIndex(dim);
  for (unsigned d = 0; d < dim; d++) solIndex[d] = mlSol.GetIndex(solName[d].c_str());

  unsigned solType = mlSol.GetSolutionType(solName[0].c_str());

  const unsigned xType = 2u; // coordinate field type

  auto &xv = msh._topology->_Sol;
  std::vector<double> x(dim);

  const unsigned offset = msh._elementOffset[iproc];
  const unsigned offsetp1 = msh._elementOffset[iproc + 1];

  std::vector<NumericVector*> solVec(dim);
  for(unsigned d = 0; d < dim; d++) {
    solVec[d] = sol._Sol[solIndex[d]];
    solVec[d]->zero();
  }

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

      auto vel = Velocity(x, time, period);
      for(unsigned d = 0; d < dim; d++) solVec[d]->set(solDof, vel[d]);
    }
  }

  for(unsigned d = 0; d < dim; d++) solVec[d]->close();
}





void GetCutElementPoints(MultiLevelSolution & mlSol, const std::string & name,
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

        }
        else {
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

static void WritePointsVTK(const std::string & filename,
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

void GetAllSolutionPoints(MultiLevelSolution & mlSol, const std::string & name,
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

void ProjectSolution(MultiLevelSolution & mlSol0 /* marker receive */,
                     MultiLevelSolution & mlSol1 /* marker send */,
                     BBoxToIel & bbox, RungeKutta & rk,
                     const std::vector<std::string> solName,
                     const std::vector<std::string> vName,
                     const Boundary & bd,
                     const double dt,
                     const double time,
                     const double period) {

  assert(!solName.empty());
  const unsigned nFields = solName.size();

  const unsigned solType =
    mlSol0.GetSolutionType(solName[0].c_str());

  assert(solType <= 2);

  for (unsigned k = 0; k < solName.size(); ++k) {

    const unsigned solType0 =
      mlSol0.GetSolutionType(solName[k].c_str());

    const unsigned solType1 =
      mlSol1.GetSolutionType(solName[k].c_str());

    assert(solType0 == solType);
    assert(solType1 == solType);
  }

  MultiLevelMesh &mlMsh0 = *mlSol0.GetMultilevelMesh();

  const unsigned nLevels = mlMsh0.GetNumberOfLevels();

  // Extract all Psi grid points on the finest level of mlSol1
  std::vector<MyVector<double>> X1;
  GetAllSolutionPoints(mlSol1, solName[0], X1);

  unsigned dim = X1.size();

  if(fabs(dt) > 1.0e-10) RungeKutta4(X1, mlSol0, bbox, vName, dt);

  // rk.rkBackward(X1);

  LevelMarkers l0;
  InterpolateSolution(l0, mlSol0, bbox, X1, solName);

  MultiLevelMesh &mlMsh1 = *mlSol1.GetMultilevelMesh();
  const unsigned level1 = mlMsh1.GetNumberOfLevels() - 1u;

  Solution &sol1 = *mlSol1.GetLevel(level1);

  for(unsigned k = 0; k < nFields; k++) {
    const unsigned solIndex1 = mlSol1.GetIndex(solName[k].c_str());

    auto &solVec1 = sol1._Sol[solIndex1];

    const MyVector<double> &psiProjected = l0.GetFields()[k];
    const std::vector<bool> &isInsideDomain = l0.GetPointInsideDomain();

    solVec1->zero();
    unsigned offset = psiProjected.begin();
    for (unsigned i = psiProjected.begin(); i < psiProjected.end(); ++i) {
      if (isInsideDomain[i - offset]) {
        solVec1->set(i, psiProjected[i]);
      }
      else { // TODO add boundarycondition for psi
        std::vector<double> x1 = (dim == 1) ? std::vector<double>({X1[0][i]})
                                 : (dim == 2) ? std::vector<double>({X1[0][i], X1[1][i]})
                                 : std::vector<double>({X1[0][i], X1[1][i], X1[2][i]});

        double value = bd.getValue(x1, time + dt, period, dt);
        solVec1->set(i, value);

      }
    }
    solVec1->close();
  }
}



void RungeKutta4(std::vector<MyVector<double>> &X,
                 MultiLevelSolution & mlSol,
                 BBoxToIel & bbox,
                 const std::vector<std::string> &velName,
                 const double dt) {
  const unsigned &dim = X.size();
  const unsigned rk_nsteps = 4;
  const std::vector <double> c_forward = {0., 0.5, 0.5, 1.};
  const std::vector <double> c_backward = {1., 0.5, 0.5, 0.};
  const std::vector <double> c = (dt > 0) ? c_forward : c_backward;
  const std::vector<std::vector <double> > a = {{}, {0.5}, {0, 0.5}, {0., 0., 1.}};
  const std::vector <double> b = {1. / 6., 1. / 3., 1. / 3., 1. / 6.} ;
  std::vector<std::vector<MyVector<double>>> K;
  for(unsigned rk = 0; rk < rk_nsteps; rk++) {
    rkStep(mlSol, bbox, X, K, rk, velName,  dt, c[rk], a[rk]);
  }
  for(unsigned rk = 0; rk < rk_nsteps; rk++) {
    for(unsigned d = 0; d < dim; d++) {
      for(unsigned i = X[d].begin(); i < X[d].end(); i++) {
        X[d][i] += dt * K[rk][d][i] * b[rk];
      }
    }
  }
}

void rkStep(MultiLevelSolution & mlSol,
            BBoxToIel & bbox,
            const std::vector<MyVector<double>> &X,
            std::vector<std::vector<MyVector<double>>> &K,
            const unsigned rkStep,
            const std::vector<std::string> &velName,
            const double dt,
            const double c,
            const std::vector<double> &a) {

  assert (a.size() == rkStep);
  assert(!velName.empty());
  const unsigned nFields = velName.size();

  assert(X.size() == nFields);

  assert(K.size() == rkStep);
  for (unsigned j = 0; j < rkStep; ++j) {
    assert(K[j].size() == nFields);
    for (unsigned d = 0; d < nFields; ++d) {
      assert(K[j][d].begin() == X[d].begin());
      assert(K[j][d].end()   == X[d].end());
    }
  }


  const unsigned velType =
    mlSol.GetSolutionType(velName[0].c_str());

  assert(velType <= 2);

  for (unsigned k = 0; k < velName.size(); ++k) {
    const unsigned velTypek =
      mlSol.GetSolutionType(velName[k].c_str());
    assert(velTypek == velType);
  }

  MultiLevelMesh &mlMsh = *mlSol.GetMultilevelMesh();

  const unsigned nLevels = mlMsh.GetNumberOfLevels();

  // Extract all Psi grid points on the finest level of mlSol1
  std::vector<MyVector<double>> Xk = X;
  for(unsigned d = 0; d < nFields; d++) {
    for (unsigned i = Xk[d].begin(); i < Xk[d].end(); ++i) {
      for(unsigned j = 0; j < a.size(); j++) {
        Xk[d][i] += a[j] * K[j][d][i] * dt;
      }
    }
  }

  LevelMarkers l0;

  InterpolateSolution(l0, mlSol, bbox, Xk, velName, c);

  K.resize(rkStep + 1);
  K[rkStep].resize(nFields);
  for(unsigned d = 0; d < nFields; d++) {
    K[rkStep][d] = l0.GetFields()[d];
  }
  if(rkStep > 0) { // to check if marker went out the domain
    const std::vector<bool> &isInsideDomain = l0.GetPointInsideDomain();
    const unsigned offset = K[rkStep][0].begin();
    const unsigned offsetp1 = K[rkStep][0].end();
    for (unsigned i = offset; i < offsetp1; ++i) {
      if (!isInsideDomain[i - offset]) {
        for(unsigned d = 0; d < nFields; d++) {
          K[rkStep][d][i] = K[rkStep - 1][d][i];
        }
      }
    }
  }
}


void InterpolateSolution(LevelMarkers & l0,
                         MultiLevelSolution & mlSol0,
                         BBoxToIel & bbox,
                         std::vector<MyVector<double>> &X,
                         const std::vector< std::string >solName,
                         const double c) {

  const unsigned &nFields = solName.size();
  MultiLevelMesh &mlMsh0 = *mlSol0.GetMultilevelMesh();
  const unsigned nLevels = mlMsh0.GetNumberOfLevels();
  assert(bbox.GetLevel() < nLevels);
  const unsigned bboxLevels = nLevels - bbox.GetLevel();

  std::vector<LevelMarkers> lX(bboxLevels);

  bbox.GetInverseMappingOnCoarseLevel(X, l0, lX[0]);

  for (unsigned k = 1; k < bboxLevels; ++k) {
    bbox.Project(mlMsh0, lX[k - 1], lX[k]);
  }

  const unsigned level0 = nLevels - 1u;

  Mesh &msh0 = *mlMsh0.GetLevel(level0);
  Solution &sol0 = *mlSol0.GetLevel(level0);

  LevelMarkers &lTop = lX.back();
  lTop.GetFields().resize(nFields);

  std::vector<MyVector<double>> &Xi = lTop.GetLocalCoordinates();
  MyVector<unsigned> &Iel = lTop.GetElements();

  const unsigned dim = Xi.size();

  for(unsigned d = 0; d < nFields; d++) {
    const unsigned solIndex0 = mlSol0.GetIndex(solName[d].c_str());

    const unsigned solType = mlSol0.GetSolutionType(solName[d].c_str());

    auto &solNew = sol0._Sol[solIndex0];
    auto &solOld = ( fabs(c - 1.) < 1e-5) ? sol0._Sol[solIndex0] : sol0._SolOld[solIndex0];

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

      const unsigned nDof = msh0.GetElementDofNumber(iel, solType);

      phi.resize(nDof);

      msh0._finiteElement[ielType][solType]->GetPhi(phi, xi);

      double value = 0.0;
      for (unsigned j = 0; j < nDof; ++j) {
        const unsigned solDof = msh0.GetSolutionDof(j, iel, solType);
        value += phi[j] * ((1. - c) * (*solOld)(solDof) + c * (*solNew)(solDof));
      }

      psiLocal[ip - Iel.begin()] = value;
    }
    lTop.GetFields()[d].buildFromLocal(psiLocal);
  }

  // Project Psi backward through the marker hierarchy
  std::vector<std::vector<double>> Wfield_r;
  std::vector<std::vector<double>> Wfield_s;


  const bool backward = true;

  for (int l = static_cast<int>(bboxLevels) - 1; l >= 1; --l) {
    lX[l].RebuildLocalFromField(Wfield_r, nFields, backward);
    lX[l - 1].SendLocalField(Wfield_r, Wfield_s);
    lX[l - 1].RebuildFieldFromLocal(Wfield_s, nFields, backward);
  }

  lX[0].RebuildLocalFromField(Wfield_r, nFields, backward);
  l0.SendLocalField(Wfield_r, Wfield_s);
  l0.RebuildFieldFromLocal(Wfield_s, nFields, backward);
}

void GetSolutionGradient(MultiLevelSolution & mlSol, const std::string & solName, std::vector<std::string> &gradSolName) {

  MultiLevelMesh &mlMsh = *mlSol.GetMultilevelMesh();
  const unsigned level = mlMsh.GetNumberOfLevels() - 1u;
  Mesh &msh = *mlMsh.GetLevel(level);
  Solution &sol = *mlSol.GetLevel(level);
  const unsigned iproc = msh.processor_id();
  const unsigned dim = msh.GetDimension();

  unsigned solIndex = mlSol.GetIndex(solName.c_str());
  unsigned gammaIndex = mlSol.GetIndex("Gamma");

  std::vector<unsigned> gradSolIndex(dim);
  for(unsigned d = 0; d < dim; d++) gradSolIndex[d] = mlSol.GetIndex(gradSolName[d].c_str());

  unsigned solType = mlSol.GetSolutionType(solName.c_str());

  const unsigned xType = 2u; // coordinate field type

  auto &xv = msh._topology->_Sol;
  std::vector<std::vector<double>> x(dim);
  std::vector<double> isol;
  std::vector<double> phi;
  std::vector<double> phi_x;
  double weight;
  std::vector<unsigned> idof;


  auto &gammaVec = sol._Sol[gammaIndex];
  auto &solVec = sol._Sol[solIndex];

  std::vector<NumericVector*> gradSolVec(dim);


  gammaVec->zero();
  for(unsigned d = 0; d < dim; d++)  {
    gradSolVec[d] = sol._Sol[gradSolIndex[d]];
    gradSolVec[d]->zero();
  }


  unsigned offset = msh._elementOffset[iproc];
  unsigned offsetp1 = msh._elementOffset[iproc + 1];
  // Loop over local elements and interpolate psi2D at solution DoFs
  for (unsigned iel = offset; iel < offsetp1; ++iel) {

    unsigned ielType = msh.GetElementType(iel);

    const unsigned nDof = msh.GetElementDofNumber(iel, solType);

    isol.resize(nDof);
    idof.resize(nDof);
    for(unsigned d = 0; d < dim; d++)  {
      x[d].resize(nDof);
    }

    for (unsigned i = 0; i < nDof; ++i) {

      const unsigned iDof = msh.GetSolutionDof(i, iel, solType);
      idof[i] = iDof;
      isol[i] = (*solVec)(iDof);
      // Get physical coordinates of the current DoF
      const unsigned xDof = msh.GetSolutionDof(i, iel, xType);
      for (unsigned d = 0; d < dim; ++d) {
        x[d][i] = (*xv[d])(xDof);
      }
    }

    for(unsigned ig = 0; ig < msh._finiteElement[ielType][solType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh._finiteElement[ielType][solType]->Jacobian(x, ig, weight, phi, phi_x);

      std::vector<double> grad(dim, 0.);

      for(unsigned i = 0; i < nDof; i++) {
        for(unsigned d = 0; d < dim; d++) {
          grad[d] += isol[i] * phi_x[d + i * dim] * weight;
        }
      }

      for(unsigned i = 0; i < nDof; i++) {
        for(unsigned d = 0; d < dim; d++) {
          (*gradSolVec[d]).add(idof[i], grad[d] * phi[i] * weight);
        }
        (*gammaVec).add(idof[i], phi[i] * weight);
      }


    }
  }

  gammaVec->close();
  for(unsigned d = 0; d < dim; d++)  {
    gradSolVec[d]->close();
  }

  offset = msh._dofOffset[solType][iproc];
  offsetp1 = msh._dofOffset[solType][iproc + 1];

  for(unsigned i = offset; i < offsetp1; i++) {
    for(unsigned d = 0; d < dim; d++) {
      double value = (*gradSolVec[d])(i);
      (*gradSolVec[d]).set(i, value / (*gammaVec)(i) );
    }
  }
  for(unsigned d = 0; d < dim; d++)  {
    gradSolVec[d]->close();
  }

}

double ComputeArea(MultiLevelSolution & mlSol,
                   const std::string & solName) {

  MultiLevelMesh &mlMsh = *mlSol.GetMultilevelMesh();
  const unsigned level = mlMsh.GetNumberOfLevels() - 1u;

  Mesh &msh = *mlMsh.GetLevel(level);
  Solution &sol = *mlSol.GetLevel(level);

  const unsigned iproc = msh.processor_id();
  const unsigned dim   = msh.GetDimension();

  const unsigned solIndex = mlSol.GetIndex(solName.c_str());

  const unsigned solType = mlSol.GetSolutionType(solName.c_str());

  const unsigned xType = 2u; // coordinate field type

  auto &xv = msh._topology->_Sol;
  auto &solVec = sol._Sol[solIndex];

  std::vector<std::vector<double>> x(dim);
  std::vector<double> isol;

  std::vector<double> phi;
  std::vector<double> phi_x;

  double weight = 0.0;

  auto Hpos = [](const double psi) -> double {
    return (psi > 0.0) ? 1.0 : 0.0;
  };

  double localArea = 0.0;

  const unsigned offset = msh._elementOffset[iproc];

  const unsigned offsetp1 = msh._elementOffset[iproc + 1];

  for (unsigned iel = offset; iel < offsetp1; ++iel) {

    const unsigned ielType = msh.GetElementType(iel);

    const unsigned nDof = msh.GetElementDofNumber(iel, solType);

    isol.resize(nDof);

    for (unsigned d = 0; d < dim; ++d) {
      x[d].resize(nDof);
    }

    for (unsigned i = 0; i < nDof; ++i) {

      const unsigned iDof = msh.GetSolutionDof(i, iel, solType);

      isol[i] = (*solVec)(iDof);

      const unsigned xDof = msh.GetSolutionDof(i, iel, xType);

      for (unsigned d = 0; d < dim; ++d) {
        x[d][i] = (*xv[d])(xDof);
      }
    }

    double A = 0.0;
    double Ap = 0.0;

    const unsigned nGauss = msh._finiteElement[ielType][solType]->GetGaussPointNumber();

    for (unsigned ig = 0; ig < nGauss; ++ig) {

      phi.clear();
      phi_x.clear();

      msh._finiteElement[ielType][solType]->Jacobian(x, ig, weight, phi, phi_x);

      A += weight;

      double psi_q = 0.0;

      for (unsigned i = 0; i < nDof; ++i) {
        psi_q += isol[i] * phi[i];
      }

      Ap += weight * Hpos(psi_q);
    }

    if (A > 0.0) {
      localArea += Ap;
    }
  }

  double area = 0.0;

  MPI_Allreduce(&localArea,
                &area,
                1,
                MPI_DOUBLE,
                MPI_SUM,
                MPI_COMM_WORLD);

  if (iproc == 0) {
    std::cout << "=======================================" << std::endl;
    std::cout
        << "Internal area(" << solName << ") = " << area
        << std::endl;
    std::cout << "=======================================" << std::endl;
  }

  return area;
}



