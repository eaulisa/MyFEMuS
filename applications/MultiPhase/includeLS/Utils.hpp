#pragma once

#include "CutFemWeight.hpp"
#include "CDWeights.hpp"
typedef double TypeIO;
typedef cpp_bin_float_oct TypeA;
typedef cpp_bin_float_oct oct;

// CutFemWeight <double, double> quad = CutFemWeight<double, double>(QUAD, 5, "legendre");
CutFemWeight <TypeIO, TypeA> quad  = CutFemWeight<TypeIO, TypeA >(QUAD, 5, "legendre");
CutFemWeight <TypeIO, TypeA> tri  = CutFemWeight<TypeIO, TypeA >(TRI, 5, "legendre");
CutFemWeight <TypeIO, TypeA> hex  = CutFemWeight<TypeIO, TypeA >(HEX, 5, "legendre");
CutFemWeight <TypeIO, TypeA> tet  = CutFemWeight<TypeIO, TypeA >(TET, 5, "legendre");
CutFemWeight <TypeIO, TypeA> wedge  = CutFemWeight<TypeIO, TypeA >(WEDGE, 5, "legendre");

const std::vector< CutFemWeight <TypeIO, TypeA> *> cfw = {&hex, &tet, &wedge, &quad, &tri};

Fem fem = Fem(hex.GetGaussQuadratureOrder(), hex.GetDimension());

void RungeKutta4(std::vector<MyVector<double>> &X,
                 MultiLevelSolution & mlSol,
                 BBoxToIel & bbox,
                 const std::vector<std::string> &velName,
                 const unsigned vlevel,
                 const double dt);

void rkStep(MultiLevelSolution & mlSol,
            BBoxToIel & bbox,
            const std::vector<MyVector<double>> &X,
            std::vector<std::vector<MyVector<double>>> &K,
            const unsigned rkStep,
            const std::vector<std::string> &velName,
            const unsigned vlevel,
            const double dt,
            const double c,
            const std::vector<double> &a);

void InterpolateSolution(LevelMarkers &l0,
                         MultiLevelSolution &mlSol0,
                         BBoxToIel &bbox,
                         std::vector<MyVector<double>> &X,
                         const std::vector< std::string >solName,
                         const unsigned level,
                         const double c
                        );

void ProjectSolution(MultiLevelSolution &mlSol0 /* target */, MultiLevelSolution &mlSol1 /* source */,
                     BBoxToIel &bbox,
                     const std::vector<std::string> solName,
                     const unsigned s0Level,
                     const unsigned s1Level,
                     const std::vector<std::string> vName = {},
                     const unsigned vLevel = UINT_MAX,
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

    auto el_lm1 = msh_lm1.el;
    auto el_l = msh_l.el;

    for(unsigned iel_l = msh_l._elementOffset[iproc];
        iel_l < msh_l._elementOffset[iproc + 1];
        iel_l++) {

      const unsigned iel_lm1 =
        static_cast<unsigned>((*father)(iel_l));

      double value_l = (*solC_l)(iel_l);

      solC_lm1->add(iel_lm1, value_l);

      // if (value_l == 0.) {
      //   el_l->SetElementMaterial(iel_l, 4);
      // }
      // else if (value_l < 0.9) {
      //   el_l->SetElementMaterial(iel_l, 3);
      // }
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
          // el_lm1->SetElementMaterial(iel_lm1, 4);
        }
        else if(value >= maxNumberOfChildren - tol) {
          solC_lm1->set(iel_lm1, 1.);
        }
        else {
          solC_lm1->set(iel_lm1, 0.5);
          // el_lm1->SetElementMaterial(iel_lm1, 3);
        }
      }
    }

    solC_lm1->close();

    delete father;
  }
}

void SetUnphysicalPressureDofs(MultiLevelSolution& mlSol, const std::string& CName, const std::vector<std::string>& PName,
                               const unsigned level0, const unsigned level1, const std::vector<double>& xtarget, const bool fixPressureAtOnePoint) {
  if (PName.size() != 2) {
    std::cout << "Error! Expected two pressure fields.\n";
    abort();
  }

  MultiLevelMesh& mlMsh = *mlSol.GetMultilevelMesh();
  const unsigned solCIndex = mlSol.GetIndex(CName.c_str());
  const unsigned solCType = mlSol.GetSolutionType(CName.c_str());

  if (solCType != 3) {
    std::cout << "Error! The C field must be PWC.\n";
    abort();
  }

  std::vector<unsigned> solPIndex(2);
  std::vector<unsigned> solPType(2);

  for (unsigned n = 0; n < 2; ++n) {
    solPIndex[n] = mlSol.GetIndex(PName[n].c_str());
    solPType[n] = mlSol.GetSolutionType(PName[n].c_str());
  }

  if (solPType[0] != solPType[1]) {
    std::cout << "Error! Pressure fields must have the same FE type.\n";
    abort();
  }

  for (unsigned l = level0; l <= level1; ++l) {
    Mesh& msh = *mlMsh.GetLevel(l);
    const unsigned iproc = msh.processor_id();
    const unsigned dim = msh.GetDimension();

    if (xtarget.size() < dim) {
      std::cout << "Error! xtarget has wrong dimension.\n";
      abort();
    }

    auto& solC = (mlSol.GetSolutionLevel(l))->_Sol[solCIndex];
    auto& solP1 = (mlSol.GetSolutionLevel(l))->_Sol[solPIndex[0]];
    auto& solP2 = (mlSol.GetSolutionLevel(l))->_Sol[solPIndex[1]];
    auto& solP1Bdc = (mlSol.GetSolutionLevel(l))->_Bdc[solPIndex[0]];
    auto& solP2Bdc = (mlSol.GetSolutionLevel(l))->_Bdc[solPIndex[1]];

    NumericVector* supportP1 = NumericVector::build().release();
    NumericVector* supportP2 = NumericVector::build().release();

    supportP1->init(*solP1);
    supportP2->init(*solP2);
    supportP1->zero();
    supportP2->zero();

    for (unsigned iel = msh._elementOffset[iproc]; iel < msh._elementOffset[iproc + 1]; ++iel) {
      const double C = (*solC)(iel);
      const unsigned nDofP1 = msh.GetElementDofNumber(iel, solPType[0]);
      const unsigned nDofP2 = msh.GetElementDofNumber(iel, solPType[1]);

      if (C >= 0.1) {
        for (unsigned i = 0; i < nDofP1; ++i) {
          const unsigned dof = msh.GetSolutionDof(i, iel, solPType[0]);
          supportP1->add(dof, 1.0);
        }
      }

      if (C <= 0.9) {
        for (unsigned i = 0; i < nDofP2; ++i) {
          const unsigned dof = msh.GetSolutionDof(i, iel, solPType[1]);
          supportP2->add(dof, 1.0);
        }
      }
    }

    supportP1->close();
    supportP2->close();

    const unsigned firstDof = solP1->first_local_index();
    const unsigned lastDof = solP1->last_local_index();

    for (unsigned dof = firstDof; dof < lastDof; ++dof) {
      if ((*supportP1)(dof) < 0.5) {
        solP1->set(dof, 0.0);
        solP1Bdc->set(dof, 0.0);
      }
    }

    for (unsigned dof = firstDof; dof < lastDof; ++dof) {
      if ((*supportP2)(dof) < 0.5) {
        solP2->set(dof, 0.0);
        solP2Bdc->set(dof, 0.0);
      }
    }

    if(fixPressureAtOnePoint) {

      double minDistance = std::numeric_limits<double>::max();
      int minDof = -1;

      for (unsigned iel = msh._elementOffset[iproc]; iel < msh._elementOffset[iproc + 1]; ++iel) {
        const unsigned nDofP = msh.GetElementDofNumber(iel, solPType[0]);

        if (solPType[0] == 3) {
          const unsigned dof = msh.GetSolutionDof(0, iel, solPType[0]);

          if (dof < firstDof || dof >= lastDof) continue;

          const unsigned nDofX = msh.GetElementDofNumber(iel, 2);
          const unsigned xDof = msh.GetSolutionDof(nDofX - 1, iel, 2);
          double distance2 = 0.0;

          for (unsigned k = 0; k < dim; ++k) {
            const double dx = (*msh._topology->_Sol[k])(xDof) - xtarget[k];
            distance2 += dx * dx;
          }

          if (((*supportP1)(dof) > 0.5 || (*supportP2)(dof) > 0.5) && distance2 < minDistance) {
            minDistance = distance2;
            minDof = static_cast<int>(dof);
          }
        }
        else {
          for (unsigned i = 0; i < nDofP; ++i) {
            const unsigned dof = msh.GetSolutionDof(i, iel, solPType[0]);

            if (dof < firstDof || dof >= lastDof) continue;

            const unsigned xDof = msh.GetSolutionDof(i, iel, 2);
            double distance2 = 0.0;

            for (unsigned k = 0; k < dim; ++k) {
              const double dx = (*msh._topology->_Sol[k])(xDof) - xtarget[k];
              distance2 += dx * dx;
            }

            if (((*supportP1)(dof) > 0.5 || (*supportP2)(dof) > 0.5) && distance2 < minDistance) {
              minDistance = distance2;
              minDof = static_cast<int>(dof);
            }
          }
        }
      }

      struct {
        double value;
        int index;
      } localMin, globalMin;

      localMin.value = minDistance;
      localMin.index = minDof;

      MPI_Allreduce(&localMin, &globalMin, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

      if (globalMin.index >= 0) {
        const unsigned dof = static_cast<unsigned>(globalMin.index);

        if (dof >= firstDof && dof < lastDof) {
          if ((*supportP1)(dof) > 0.5) {
            solP1->set(dof, 0.0);
            solP1Bdc->set(dof, 0.0);
          }
          else {
            solP2->set(dof, 0.0);
            solP2Bdc->set(dof, 0.0);
          }
        }
      }
    }

    solP1->close();
    solP2->close();
    solP1Bdc->close();
    solP2Bdc->close();

    delete supportP1;
    delete supportP2;
  }
}

inline std::vector<double>
Velocity(std::vector<double> &xp, const double time, double period) noexcept {
  double u = 0.0, v = 0.0, w = 0.0;

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
      w = 0.;

      break;
    }
    case RungeKutta::VelKind::Rotation: {

      u = xp[1];
      v = -xp[0];
      w = 0.;

      break;
    }
    case RungeKutta::VelKind::Translation: {

      u = 0;
      v = -0.3 * 4 * (0.5 - xp[0]) * (0.5 + xp[0]);
      // v = -0.3;
      w = 0.;

      break;
    }
    case RungeKutta::VelKind::Zero: {

      u = 0;
      v = 0;
      w = 0.;

      break;
    }

  }

  if (xp.size() == 2) return {u, v};
  else if (xp.size() == 3) return {u, v, w};
}

void Shift(std::vector<MyVector<double>> &X, const std::vector<double> &dx) {
  for (unsigned k = 0; k < X.size(); k++) {
    for (unsigned i = X[k].begin(); i < X[k].end(); i++) {
      X[k][i] += dx[k];
    }
  }
}

template <class PsiType>
void FlagFinestMeshLevel(MultiLevelMesh& mlMsh, const PsiType& psi) {

  const unsigned level = mlMsh.GetNumberOfLevels() - 1;
  Mesh& msh = *mlMsh.GetLevel(level);

  const unsigned iproc = msh.processor_id();
  const unsigned dim = msh.GetDimension();
  const unsigned xType = 2;
  const unsigned amrIndex = msh.GetAmrIndex();

  auto& xv = msh._topology->_Sol;
  auto* amrFlag = msh._topology->_Sol[amrIndex];

  amrFlag->zero();

  const unsigned offset = msh._elementOffset[iproc];
  const unsigned offsetp1 = msh._elementOffset[iproc + 1];

  std::vector<double> x(dim);

  for(unsigned iel = offset; iel < offsetp1; ++iel) {

    const unsigned nDof = msh.GetElementDofNumber(iel, xType);

    const unsigned xDof0 = msh.GetSolutionDof(0, iel, xType);

    for(unsigned k = 0; k < dim; ++k)
      x[k] = (*xv[k])(xDof0);

    const double psi0 = psi(x);
    const int sign0 = (psi0 > 0.) ? 1 : -1;

    bool signChange = false;

    for(unsigned i = 1; i < nDof; ++i) {

      const unsigned xDof = msh.GetSolutionDof(i, iel, xType);

      for(unsigned k = 0; k < dim; ++k)
        x[k] = (*xv[k])(xDof);

      const double psii = psi(x);
      const int signi = (psii > 0.) ? 1 : -1;

      if(signi != sign0) {
        signChange = true;
        break;
      }
    }

    if(signChange)
      amrFlag->set(iel, 2);
  }

  amrFlag->close();

  for(unsigned iel = offset; iel < offsetp1; ++iel) {

    if((*amrFlag)(iel) == 2) {

      for(unsigned j = 0; j < msh.el->GetElementNearElementSize(iel, 1); ++j) {

        const unsigned jel = msh.el->GetElementNearElement(iel, j);

        if(offset <= jel && jel < offsetp1) {

          if((*amrFlag)(jel) < 0.5)
            amrFlag->set(jel, 1);
        }
        else {
          amrFlag->add(jel, 1);
        }
      }
    }
  }

  amrFlag->close();

  unsigned localRefined = 0;

  for(unsigned iel = offset; iel < offsetp1; ++iel) {

    if((*amrFlag)(iel) > 0.5) {

      amrFlag->set(iel, 1);

      ++localRefined;
    }
  }

  amrFlag->close();

  unsigned globalRefined = 0;

  MPI_Allreduce(
    &localRefined,
    &globalRefined,
    1,
    MPI_UNSIGNED,
    MPI_SUM,
    MPI_COMM_WORLD);

  msh.el->SetRefinedElementNumber(globalRefined);
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

template <class PsiType>
void InitLevelSet(MultiLevelSolution & mlSol, const std::string & name,
                  const PsiType & psi2D) {

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

void GetAllSolutionPoints(MultiLevelSolution & mlSol, const std::string & name, const unsigned s0Level,
                          std::vector<MyVector<double>> &X) {

  MultiLevelMesh &mlMsh = *mlSol.GetMultilevelMesh();
  const unsigned level = s0Level;//mlMsh.GetNumberOfLevels() - 1u;
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
                     BBoxToIel & bbox,
                     const std::vector<std::string> solName,
                     const unsigned s0Level,
                     const unsigned s1Level,
                     const std::vector<std::string> vName,
                     const unsigned vLevel,
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

  // Extract all Psi grid points on the finest level of mlSol1
  std::vector<MyVector<double>> X1;
  GetAllSolutionPoints(mlSol1, solName[0], s1Level, X1);

  unsigned dim = X1.size();

  if(fabs(dt) > 1.0e-10) RungeKutta4(X1, mlSol0, bbox, vName, vLevel, dt);

  LevelMarkers l0;
  double useSol = 1.; // rather than solOld = 0.
  InterpolateSolution(l0, mlSol0, bbox, X1, solName, s0Level, useSol);

  MultiLevelMesh &mlMsh1 = *mlSol1.GetMultilevelMesh();
  const unsigned level1 = s1Level;

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
        throw std::runtime_error("POINT OUTSIDE THE DOMAIN");
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
                 const unsigned vLevel,
                 const double dt) {
  const unsigned &dim = X.size();
  const unsigned rk_nsteps = 4;
  // const std::vector <double> c_forward = {0., 0.5, 0.5, 1.};
  // const std::vector <double> c_backward = {1., 0.5, 0.5, 0.};
  // const std::vector <double> c = (dt > 0) ? c_forward : c_backward;
  const std::vector <double> c = {1., 1., 1., 1.};
  const std::vector<std::vector <double> > a = {{}, {0.5}, {0, 0.5}, {0., 0., 1.}};
  const std::vector <double> b = {1. / 6., 1. / 3., 1. / 3., 1. / 6.} ;
  std::vector<std::vector<MyVector<double>>> K;
  for(unsigned rk = 0; rk < rk_nsteps; rk++) {
    rkStep(mlSol, bbox, X, K, rk, velName, vLevel, dt, c[rk], a[rk]);
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
            const std::vector<std::string> &vName,
            const unsigned vLevel,
            const double dt,
            const double c,
            const std::vector<double> &a) {

  assert (a.size() == rkStep);
  assert(!vName.empty());
  const unsigned nFields = vName.size();

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
    mlSol.GetSolutionType(vName[0].c_str());

  assert(velType <= 2);

  for (unsigned k = 0; k < vName.size(); ++k) {
    const unsigned velTypek =
      mlSol.GetSolutionType(vName[k].c_str());
    assert(velTypek == velType);
  }

  MultiLevelMesh &mlMsh = *mlSol.GetMultilevelMesh();

  //const unsigned nLevels = mlMsh.GetNumberOfLevels();

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

  InterpolateSolution(l0, mlSol, bbox, Xk, vName, vLevel,  c);

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
                         const unsigned level,
                         const double c) {

  const unsigned &nFields = solName.size();
  MultiLevelMesh &mlMsh0 = *mlSol0.GetMultilevelMesh();
  const unsigned nLevels = level + 1u; //mlMsh0.GetNumberOfLevels();
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

double GetMaxElementH(Mesh* msh, const unsigned level) {
  const unsigned dim = msh->GetDimension();
  const unsigned iproc = msh->processor_id();
  const unsigned solXType = 2;

  double hLocalMax = 0.;

  std::vector<std::vector<double>> coordX(dim);
  std::vector<double> phiX;
  std::vector<double> phiX_x;

  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; ++iel) {
    if(msh->el->GetElementLevel(iel) != static_cast<int>(level)) continue;

    const unsigned ielGeom = msh->GetElementType(iel);
    const unsigned nDofsX = msh->GetElementDofNumber(iel, solXType);

    for(unsigned d = 0; d < dim; ++d) coordX[d].resize(nDofsX);

    for(unsigned i = 0; i < nDofsX; ++i) {
      const unsigned xDof = msh->GetSolutionDof(i, iel, solXType);

      for(unsigned d = 0; d < dim; ++d) coordX[d][i] = (*msh->_topology->_Sol[d])(xDof);
    }

    const elem_type* femX = msh->_finiteElement[ielGeom][solXType];

    double cellMeasure = 0.;

    for(unsigned ig = 0; ig < femX->GetGaussPointNumber(); ++ig) {
      double weightX = 0.;

      femX->Jacobian(coordX, ig, weightX, phiX, phiX_x);

      cellMeasure += weightX;
    }

    const double h = std::pow(cellMeasure, 1.0 / static_cast<double>(dim));

    hLocalMax = std::max(hLocalMax, h);
  }

  double hGlobalMax = 0.;

  MPI_Allreduce(&hLocalMax, &hGlobalMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  return hGlobalMax;
}

void BestFitLinearInterpolation(std::vector<const double*>& xg,
                                std::vector<double>& psi,
                                std::vector<double>& B) {

  const unsigned n = psi.size();
  const unsigned dim = xg.size();
  const unsigned m = dim + 1;

  if(n == 0) {
    throw std::runtime_error("BestFitLinearInterpolation: empty psi vector");
  }

  for(unsigned d = 0; d < dim; d++) {
    if(xg[d] == nullptr) {
      throw std::runtime_error("BestFitLinearInterpolation: null Gauss coordinate pointer");
    }
  }

  double s2 = 0.;
  double maxAbsPsi = 0.;

  for(unsigned i = 0; i < n; i++) {
    s2 += psi[i] * psi[i];
    maxAbsPsi = std::max(maxAbsPsi, std::fabs(psi[i]));
  }

  s2 /= static_cast<double>(n);

  if(maxAbsPsi <= std::numeric_limits<double>::min() ||
      s2 <= std::numeric_limits<double>::min()) {
    throw std::runtime_error("BestFitLinearInterpolation: degenerate level-set values");
  }

  const double beta = 10.;

  std::vector<double> weights(n, 0.);
  double maxWeight = 0.;

  for(unsigned i = 0; i < n; i++) {
    const double q = psi[i] * psi[i] / s2;
    weights[i] = std::exp(-beta * q);
    maxWeight = std::max(maxWeight, weights[i]);
  }

  if(maxWeight <= std::numeric_limits<double>::min()) {
    throw std::runtime_error("BestFitLinearInterpolation: all weights vanished");
  }

  for(unsigned i = 0; i < n; i++) {
    weights[i] /= maxWeight;
  }

  std::vector<std::vector<double>> M(m, std::vector<double>(m, 0.));
  std::vector<double> F(m, 0.);

  for(unsigned i = 0; i < n; i++) {

    const double f = psi[i];
    const double w = weights[i];

    for(unsigned d = 0; d < dim; d++) {

      F[d] += w * xg[d][i] * f;

      for(unsigned e = 0; e < dim; e++) {
        M[d][e] += w * xg[d][i] * xg[e][i];
      }

      M[d][dim] += w * xg[d][i];
      M[dim][d] += w * xg[d][i];
    }

    M[dim][dim] += w;
    F[dim] += w * f;
  }

  double matrixScale = 0.;

  for(unsigned i = 0; i < m; i++) {
    for(unsigned j = 0; j < m; j++) {
      matrixScale = std::max(matrixScale, std::fabs(M[i][j]));
    }
  }

  if(matrixScale <= std::numeric_limits<double>::min()) {
    throw std::runtime_error("BestFitLinearInterpolation: zero least-squares matrix");
  }

  B = F;

  const double pivotTolerance = 1.e-12;

  for(unsigned k = 0; k < m; k++) {

    unsigned pivot = k;

    for(unsigned i = k + 1; i < m; i++) {
      if(std::fabs(M[i][k]) > std::fabs(M[pivot][k])) {
        pivot = i;
      }
    }

    if(std::fabs(M[pivot][k]) <= pivotTolerance * matrixScale) {
      throw std::runtime_error("BestFitLinearInterpolation: singular least-squares matrix");
    }

    if(pivot != k) {
      std::swap(M[k], M[pivot]);
      std::swap(B[k], B[pivot]);
    }

    for(unsigned i = k + 1; i < m; i++) {

      const double factor = M[i][k] / M[k][k];

      M[i][k] = 0.;

      for(unsigned j = k + 1; j < m; j++) {
        M[i][j] -= factor * M[k][j];
      }

      B[i] -= factor * B[k];
    }
  }

  for(int ii = static_cast<int>(m) - 1; ii >= 0; ii--) {

    const unsigned i = static_cast<unsigned>(ii);

    double rhs = B[i];

    for(unsigned j = i + 1; j < m; j++) {
      rhs -= M[i][j] * B[j];
    }

    if(std::fabs(M[i][i]) <= pivotTolerance * matrixScale) {
      throw std::runtime_error("BestFitLinearInterpolation: singular back substitution");
    }

    B[i] = rhs / M[i][i];
  }

  double normGradient = 0.;

  for(unsigned d = 0; d < dim; d++) {
    normGradient += B[d] * B[d];
  }

  normGradient = std::sqrt(normGradient);

  if(normGradient <= 1.e-12) {
    throw std::runtime_error("BestFitLinearInterpolation: degenerate fitted gradient");
  }

  for(unsigned d = 0; d < m; d++) {
    B[d] /= normGradient;
  }

  double correlation = 0.;

  for(unsigned i = 0; i < n; i++) {

    double fittedPsi = B[dim];

    for(unsigned d = 0; d < dim; d++) {
      fittedPsi += B[d] * xg[d][i];
    }

    correlation += weights[i] * psi[i] * fittedPsi;
  }

  if(correlation < 0.) {
    for(unsigned d = 0; d < m; d++) {
      B[d] *= -1.;
    }
  }
}

enum class SimulationCase {
  PAR,
  RB1,
  RB2,
  SAYE1,
  SAYE2
};

std::string SimulationCaseName(const SimulationCase s) {
  switch(s) {
    case SimulationCase::PAR:
      return "par";
    case SimulationCase::RB1:
      return "rb1";
    case SimulationCase::RB2:
      return "rb2";
    case SimulationCase::SAYE1:
      return "saye1";
    case SimulationCase::SAYE2:
      return "saye2";
  }

  throw std::runtime_error("Invalid SimulationCase");
}

SimulationCase ParseSimulationCase(const std::string& s) {

  if(s == "par")  return SimulationCase::PAR;
  if(s == "rb1")  return SimulationCase::RB1;
  if(s == "rb2")  return SimulationCase::RB2;
  if(s == "saye1") return SimulationCase::SAYE1;
  if(s == "saye2") return SimulationCase::SAYE2;

  throw std::runtime_error(
    "Unknown simulation '" + s +
    "'. Valid values are: par, rb1, rb2, saye"
  );
}

struct LevelSetDiagnostics {
  double innerArea = 0.;
  double outerArea = 0.;
  double totalArea = 0.;
  double interfaceLength = 0.;
  std::vector<double> barycenter;
  std::vector<double> meanVelocity;
  double circularity = 0.;
};

LevelSetDiagnostics ComputeLevelSetDiagnostics(MultiLevelSolution& mlSol, const std::string& psiName,
    const std::vector<std::string>& velocityName, const std::vector<std::string>& pressureName, const unsigned velocityLevel) {
  MultiLevelMesh& mlMsh = *mlSol.GetMultilevelMesh();
  const unsigned level = mlMsh.GetNumberOfLevels() - 1u;
  Mesh* msh = mlMsh.GetLevel(level);
  Solution* sol = mlSol.GetSolutionLevel(level);
  const unsigned dim = msh->GetDimension();
  const unsigned iproc = msh->processor_id();
  const unsigned psiIndex = mlSol.GetIndex(psiName.c_str());
  const unsigned psiType = mlSol.GetSolutionType(psiIndex);
  const unsigned xType = 2;

  std::vector<unsigned> velocityIndex;
  unsigned velocityType = 0;

  std::vector<unsigned> pressureIndex;
  unsigned pressureType = 0;

  if(dim != 2 && dim != 3) {
    throw std::runtime_error("ComputeLevelSetDiagnostics: rising_bubble requires dim = 2 or 3");
  }

  if(velocityName.size() != dim) {
    throw std::runtime_error("ComputeLevelSetDiagnostics: velocityName must contain dim components");
  }

  velocityIndex.resize(dim);

  for(unsigned d = 0; d < dim; ++d) {
    velocityIndex[d] = mlSol.GetIndex(velocityName[d].c_str());
  }

  velocityType = mlSol.GetSolutionType(velocityIndex[0]);

  for(unsigned d = 1; d < dim; ++d) {
    if(mlSol.GetSolutionType(velocityIndex[d]) != velocityType) {
      throw std::runtime_error("ComputeLevelSetDiagnostics: velocity components must have the same solution type");
    }
  }

  if(pressureName.size() != 2) {
    throw std::runtime_error(
      "ComputeLevelSetDiagnostics: pressureName must contain P1 and P2"
    );
  }

  pressureIndex.resize(2);

  for(unsigned p = 0; p < 2; ++p) {
    pressureIndex[p] = mlSol.GetIndex(pressureName[p].c_str());
  }

  pressureType = mlSol.GetSolutionType(pressureIndex[0]);

  if(mlSol.GetSolutionType(pressureIndex[1]) != pressureType) {
    throw std::runtime_error(
      "ComputeLevelSetDiagnostics: pressure fields must have the same solution type"
    );
  }

  if(velocityLevel > level) {
    throw std::runtime_error(
      "ComputeLevelSetDiagnostics: velocityLevel cannot be finer than diagnostics level");
  }

  for(unsigned l = velocityLevel; l < level; ++l) {

    Solution* sol_l = mlSol.GetSolutionLevel(l);
    Solution* sol_lp1 = mlSol.GetSolutionLevel(l + 1u);
    Mesh* msh_lp1 = mlMsh.GetLevel(l + 1u);
    SparseMatrix* Pv = msh_lp1->GetCoarseToFineProjection(velocityType);
    SparseMatrix* Pp = msh_lp1->GetCoarseToFineProjection(pressureType);

    for(unsigned d = 0; d < dim; ++d) {
      sol_lp1->_Sol[velocityIndex[d]]->matrix_mult(
        *(sol_l->_Sol[velocityIndex[d]]),
        *Pv);
    }

    for(unsigned p = 0; p < 2; ++p) {
      sol_lp1->_Sol[pressureIndex[p]]->matrix_mult(
        *(sol_l->_Sol[pressureIndex[p]]),
        *Pp);
    }

  }

  double innerAreaLocal = 0.0;
  double outerAreaLocal = 0.0;
  double totalAreaLocal = 0.0;
  double interfaceLengthLocal = 0.0;

  std::vector<double> barycenterIntegralLocal(dim, 0.0);
  std::vector<double> velocityIntegralLocal(dim, 0.0);

  std::vector<std::vector<double>> coordX(dim);
  std::vector<double> psi;
  std::vector<std::vector<double>> velocity(dim);

  std::vector<double> phiPsi;
  std::vector<double> phiPsi_x;
  std::vector<double> phiX;
  std::vector<double> phiX_x;
  std::vector<double> phiVelocity;
  std::vector<double> phiVelocity_x;

  std::vector<std::vector<double>> Jacob;
  std::vector<std::vector<double>> JacI;

  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; ++iel) {
    const unsigned ielGeom = msh->GetElementType(iel);
    const unsigned nDofsPsi = msh->GetElementDofNumber(iel, psiType);
    const unsigned nDofsX = msh->GetElementDofNumber(iel, xType);

    unsigned nDofsVelocity = 0;

    nDofsVelocity = msh->GetElementDofNumber(iel, velocityType);

    psi.resize(nDofsPsi);

    for(unsigned d = 0; d < dim; ++d) {
      coordX[d].resize(nDofsX);
    }

    for(unsigned d = 0; d < dim; ++d) {
      velocity[d].resize(nDofsVelocity);
    }

    for(unsigned i = 0; i < nDofsPsi; ++i) {
      const unsigned dof = msh->GetSolutionDof(i, iel, psiType);
      psi[i] = (*sol->_Sol[psiIndex])(dof);
    }

    for(unsigned i = 0; i < nDofsVelocity; ++i) {
      const unsigned dof = msh->GetSolutionDof(i, iel, velocityType);

      for(unsigned d = 0; d < dim; ++d) {
        velocity[d][i] = (*sol->_Sol[velocityIndex[d]])(dof);
      }
    }

    for(unsigned i = 0; i < nDofsX; ++i) {
      const unsigned dof = msh->GetSolutionDof(i, iel, xType);

      for(unsigned d = 0; d < dim; ++d) {
        coordX[d][i] = (*msh->_topology->_Sol[d])(dof);
      }
    }

    const elem_type* femPsi = fem.GetFiniteElement(ielGeom, psiType);
    const elem_type* femX = fem.GetFiniteElement(ielGeom, xType);
    const elem_type* femVelocity = nullptr;

    femVelocity = fem.GetFiniteElement(ielGeom, velocityType);

    const unsigned nGauss = femPsi->GetGaussPointNumber();

    if(nGauss != cfw[ielGeom]->GetGaussQuadraturePointNumber()) {
      throw std::runtime_error("ComputeLevelSetDiagnostics: incompatible quadrature");
    }

    if(femVelocity->GetGaussPointNumber() != nGauss) {
      throw std::runtime_error("ComputeLevelSetDiagnostics: incompatible velocity quadrature");
    }

    std::vector<double> psiG(nGauss, 0.0);

    double psiMin = std::numeric_limits<double>::max();
    double psiMax = -std::numeric_limits<double>::max();

    for(unsigned i = 0; i < nDofsPsi; ++i) {
      psiMin = std::min(psiMin, psi[i]);
      psiMax = std::max(psiMax, psi[i]);
    }

    for(unsigned ig = 0; ig < nGauss; ++ig) {
      double* phi = femPsi->GetPhi(ig);

      for(unsigned i = 0; i < nDofsPsi; ++i) {
        psiG[ig] += psi[i] * phi[i];
      }

      psiMin = std::min(psiMin, psiG[ig]);
      psiMax = std::max(psiMax, psiG[ig]);
    }

    const bool cut = psiMin <= 0.0 && psiMax >= 0.0 && psiMax - psiMin > 1.e-14;

    if(!cut) {
      double elementArea = 0.0;

      for(unsigned ig = 0; ig < nGauss; ++ig) {
        double weight = 0.0;

        femPsi->Jacobian(coordX, ig, weight, phiPsi, phiPsi_x);

        elementArea += weight;

        if(psiMax > 0.0) {
          double weightX = 0.0;
          double weightVelocity = 0.0;

          femX->Jacobian(coordX, ig, weightX, phiX, phiX_x);
          femVelocity->Jacobian(coordX, ig, weightVelocity, phiVelocity, phiVelocity_x);

          std::vector<double> x(dim, 0.0);
          std::vector<double> u(dim, 0.0);

          for(unsigned i = 0; i < nDofsX; ++i) {
            for(unsigned d = 0; d < dim; ++d) {
              x[d] += coordX[d][i] * phiX[i];
            }
          }

          for(unsigned i = 0; i < nDofsVelocity; ++i) {
            for(unsigned d = 0; d < dim; ++d) {
              u[d] += velocity[d][i] * phiVelocity[i];
            }
          }

          for(unsigned d = 0; d < dim; ++d) {
            barycenterIntegralLocal[d] += x[d] * weight;
            velocityIntegralLocal[d] += u[d] * weight;
          }
        }
      }

      totalAreaLocal += elementArea;

      if(psiMax > 0.0) {
        innerAreaLocal += elementArea;
      }
      else {
        outerAreaLocal += elementArea;
      }

      continue;
    }

    std::vector<const double*> xg(dim);

    for(unsigned d = 0; d < dim; ++d) {
      xg[d] = femPsi->GetGaussRule().GetGaussCoordinatePointer(d);
    }

    std::vector<double> a;

    BestFitLinearInterpolation(xg, psiG, a);

    double interfaceConstant = a[dim];

    a.resize(dim);

    std::vector<TypeIO> weightInner(cfw[ielGeom]->GetGaussQuadraturePointNumber(), 0.0);
    std::vector<TypeIO> weightOuter(cfw[ielGeom]->GetGaussQuadraturePointNumber(), 0.0);
    std::vector<TypeIO> weightInterface(cfw[ielGeom]->GetGaussQuadraturePointNumber(), 0.0);

    (*cfw[ielGeom])(0, a, interfaceConstant, weightInner);

    for(unsigned k = 0; k < dim; ++k) {
      a[k] = -a[k];
    }

    interfaceConstant = -interfaceConstant;

    (*cfw[ielGeom])(-1, a, interfaceConstant, weightInterface);
    (*cfw[ielGeom])(0, a, interfaceConstant, weightOuter);

    for(unsigned ig = 0; ig < nGauss; ++ig) {
      double weight = 0.0;
      double weightAux = 0.0;
      double weightX = 0.0;

      femPsi->Jacobian(coordX, ig, weight, phiPsi, phiPsi_x);
      femX->Jacobian(coordX, ig, weightX, phiX, phiX_x);

      const double innerWeight = weight * static_cast<double>(weightInner[ig]);
      const double outerWeight = weight * static_cast<double>(weightOuter[ig]);

      totalAreaLocal += weight;
      innerAreaLocal += innerWeight;
      outerAreaLocal += outerWeight;

      double weightVelocity = 0.0;

      femVelocity->Jacobian(coordX, ig, weightVelocity, phiVelocity, phiVelocity_x);

      std::vector<double> x(dim, 0.0);
      std::vector<double> u(dim, 0.0);

      for(unsigned i = 0; i < nDofsX; ++i) {
        for(unsigned d = 0; d < dim; ++d) {
          x[d] += coordX[d][i] * phiX[i];
        }
      }

      for(unsigned i = 0; i < nDofsVelocity; ++i) {
        for(unsigned d = 0; d < dim; ++d) {
          u[d] += velocity[d][i] * phiVelocity[i];
        }
      }

      for(unsigned d = 0; d < dim; ++d) {
        barycenterIntegralLocal[d] += x[d] * innerWeight;
        velocityIntegralLocal[d] += u[d] * innerWeight;
      }

      femPsi->GetJacobianMatrix(coordX, ig, weight, Jacob, JacI);

      std::vector<double> Nref(dim, 0.0);

      double dsN2 = 0.0;

      for(unsigned k = 0; k < dim; ++k) {
        for(unsigned j = 0; j < dim; ++j) {
          Nref[k] += JacI[j][k] * a[j];
        }

        dsN2 += Nref[k] * Nref[k];
      }

      const double dsN = std::sqrt(dsN2);
      const double dGamma = weight * static_cast<double>(weightInterface[ig]) * dsN;

      interfaceLengthLocal += dGamma;

    }
  }

  const unsigned nGlobalQuantities = 4u + 2u * dim;

  std::vector<double> localSum(nGlobalQuantities, 0.0);
  std::vector<double> globalSum(nGlobalQuantities, 0.0);

  localSum[0] = innerAreaLocal;
  localSum[1] = outerAreaLocal;
  localSum[2] = totalAreaLocal;
  localSum[3] = interfaceLengthLocal;

  for(unsigned d = 0; d < dim; ++d) {
    localSum[4 + d] = barycenterIntegralLocal[d];
    localSum[4 + dim + d] = velocityIntegralLocal[d];
  }

  MPI_Allreduce(localSum.data(), globalSum.data(), static_cast<int>(nGlobalQuantities), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  LevelSetDiagnostics result;

  result.innerArea = globalSum[0];
  result.outerArea = globalSum[1];
  result.totalArea = globalSum[2];
  result.interfaceLength = globalSum[3];

  result.barycenter.assign(dim, 0.0);
  result.meanVelocity.assign(dim, 0.0);

  if(result.innerArea > 1.e-14) {
    for(unsigned d = 0; d < dim; ++d) {
      result.barycenter[d] = globalSum[4 + d] / result.innerArea;
      result.meanVelocity[d] = globalSum[4 + dim + d] / result.innerArea;
    }

    const double pi = std::acos(-1.0);

    if(dim == 2 && result.interfaceLength > 1.e-14) {
      result.circularity = 2.0 * std::sqrt(pi * result.innerArea) / result.interfaceLength;
    }

    if(dim == 3 && result.interfaceLength > 1.e-14) {
      result.circularity = std::cbrt(36.0 * pi * result.innerArea * result.innerArea) / result.interfaceLength;
    }
  }

  return result;
}

void PrintLevelSetDiagnostics(
  const LevelSetDiagnostics& diagnostics,
  const unsigned iproc,
  const double time,
  const std::string& fileName = "") {
  if(iproc != 0) {
    return;
  }

  // ============================================================
  // Terminal output
  // ============================================================

  if(fileName.empty()) {

    std::cout << std::setprecision(16);

    std::cout << "Time                = "
              << time << std::endl;

    std::cout << "Inner area          = "
              << diagnostics.innerArea << std::endl;

    std::cout << "Outer area          = "
              << diagnostics.outerArea << std::endl;

    std::cout << "Total area          = "
              << diagnostics.totalArea << std::endl;

    std::cout << "Interface length    = "
              << diagnostics.interfaceLength << std::endl;

    for(unsigned d = 0;
        d < diagnostics.barycenter.size();
        ++d) {

      std::cout << "Barycenter[" << d << "]      = "
                << diagnostics.barycenter[d]
                << std::endl;
    }

    for(unsigned d = 0;
        d < diagnostics.meanVelocity.size();
        ++d) {

      std::cout << "Mean velocity[" << d << "]   = "
                << diagnostics.meanVelocity[d]
                << std::endl;
    }

    std::cout << "Circularity         = "
              << diagnostics.circularity
              << std::endl;

    return;
  }

  // ============================================================
  // File output
  // ============================================================

  std::ofstream out(fileName, std::ios::app);

  if(!out) {
    throw std::runtime_error(
      "PrintLevelSetDiagnostics: cannot open file " + fileName);
  }

  const unsigned w = 24;

  out << std::scientific
      << std::setprecision(16);

  out << std::setw(w) << time
      << std::setw(w) << diagnostics.innerArea
      << std::setw(w) << diagnostics.outerArea
      << std::setw(w) << diagnostics.totalArea
      << std::setw(w) << diagnostics.interfaceLength;

  for(unsigned d = 0;
      d < diagnostics.barycenter.size();
      ++d) {

    out << std::setw(w)
        << diagnostics.barycenter[d];
  }

  for(unsigned d = 0;
      d < diagnostics.meanVelocity.size();
      ++d) {

    out << std::setw(w)
        << diagnostics.meanVelocity[d];
  }

  out << std::setw(w)
      << diagnostics.circularity;

  out << '\n';
}

struct SimulationArgs {
  unsigned uniformLevels = 0u;
  unsigned levelOffset   = 2u;

  unsigned dim = 2u;

  SimulationCase simulation = SimulationCase::RB2;

  std::vector<unsigned> adaptiveLevels;
  std::vector<unsigned> nSteps;
};

SimulationArgs ParseSimulationArgs(const int argc, char** argv) {

  SimulationArgs args;

  for (int i = 1; i < argc; ++i) {

    const std::string option = argv[i];

    if (option == "--uniform-levels") {

      if (i + 1 >= argc) {
        throw std::runtime_error(
          "Missing value after --uniform-levels"
        );
      }

      args.uniformLevels = std::stoul(argv[++i]);
    }

    else if (option == "--adaptive-levels") {

      while (i + 1 < argc) {

        const std::string next = argv[i + 1];

        if (next.rfind("--", 0) == 0)
          break;

        args.adaptiveLevels.push_back(
          std::stoul(argv[++i])
        );
      }
    }

    else if (option == "--nsteps") {

      while (i + 1 < argc) {

        const std::string next = argv[i + 1];

        if (next.rfind("--", 0) == 0)
          break;

        args.nSteps.push_back(
          std::stoul(argv[++i])
        );
      }
    }

    else if (option == "--level-offset") {

      if (i + 1 >= argc) {
        throw std::runtime_error(
          "Missing value after --level-offset"
        );
      }

      args.levelOffset = std::stoul(argv[++i]);
    }

    else if(option == "--simulation") {

      if(i + 1 >= argc) {
        throw std::runtime_error(
          "Missing value after --simulation"
        );
      }

      args.simulation =
        ParseSimulationCase(argv[++i]);
    }

    else if(option == "--dim") {

      if(i + 1 >= argc) {
        throw std::runtime_error(
          "Missing value after --dim"
        );
      }

      args.dim = std::stoul(argv[++i]);

      if(args.dim != 2u && args.dim != 3u) {
        throw std::runtime_error(
          "--dim must be 2 or 3"
        );
      }
    }

    else if (option == "--help" || option == "-h") {

      std::cout
          << "Usage:\n"
          << "  " << argv[0]
          << " --uniform-levels N"
          << " --adaptive-levels A1 [A2 ...]"
          << " --nsteps N1 [N2 ...]"
          << " [--level-offset K]\n\n"

          << "Example:\n"
          << "  " << argv[0]
          << " --uniform-levels 2"
          << " --adaptive-levels 4 5 6"
          << " --nsteps 600 1200 2400"
          << " --level-offset 2\n";

      std::exit(0);
    }

    // else {

    //   throw std::runtime_error(
    //     "Unknown command-line option: " + option
    //   );
    // }
  }

  if (args.uniformLevels == 0u) {
    throw std::runtime_error(
      "--uniform-levels must be specified and > 0"
    );
  }

  if (args.adaptiveLevels.empty()) {
    throw std::runtime_error(
      "--adaptive-levels must contain at least one value"
    );
  }

  if (args.nSteps.empty()) {
    throw std::runtime_error(
      "--nsteps must contain at least one value"
    );
  }

  for (const unsigned n : args.nSteps) {
    if (n == 0u) {
      throw std::runtime_error(
        "All --nsteps values must be > 0"
      );
    }
  }

  return args;
}

double ComputeL2Error(MultiLevelSolution& solA, MultiLevelSolution& solB, const std::vector<std::string>& fieldNames, const unsigned level) {
  MultiLevelMesh* mlMsh = solA.GetMultilevelMesh();
  Mesh* msh = mlMsh->GetLevel(level);
  Solution* A = solA.GetSolutionLevel(level);
  Solution* B = solB.GetSolutionLevel(level);
  const unsigned dim = msh->GetDimension();
  const unsigned iproc = msh->processor_id();
  const unsigned nFields = fieldNames.size();

  std::vector<unsigned> indexA(nFields);
  std::vector<unsigned> indexB(nFields);

  for(unsigned k = 0; k < nFields; ++k) {
    indexA[k] = solA.GetIndex(fieldNames[k].c_str());
    indexB[k] = solB.GetIndex(fieldNames[k].c_str());
  }

  const unsigned solType = solA.GetSolutionType(fieldNames[0].c_str());

  for(unsigned k = 0; k < nFields; ++k) {
    if(solA.GetSolutionType(fieldNames[k].c_str()) != solType) abort();
    if(solB.GetSolutionType(fieldNames[k].c_str()) != solType) abort();
  }

  double localError2 = 0.0;

  std::vector<std::vector<double>> vx(dim);
  std::vector<double> phi;
  std::vector<double> gradPhi;

  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; ++iel) {
    const unsigned ielt = msh->GetElementType(iel);
    const unsigned nDofs = msh->GetElementDofNumber(iel, solType);
    const unsigned nDofsX = msh->GetElementDofNumber(iel, 2);

    for(unsigned d = 0; d < dim; ++d) {
      vx[d].resize(nDofsX);

      for(unsigned i = 0; i < nDofsX; ++i) {
        const unsigned xDof = msh->GetSolutionDof(i, iel, 2);
        vx[d][i] = (*msh->_topology->_Sol[d])(xDof);
      }
    }

    std::vector<std::vector<double>> a(nFields, std::vector<double>(nDofs, 0.0));
    std::vector<std::vector<double>> b(nFields, std::vector<double>(nDofs, 0.0));

    for(unsigned k = 0; k < nFields; ++k) {
      for(unsigned i = 0; i < nDofs; ++i) {
        const unsigned dof = msh->GetSolutionDof(i, iel, solType);
        a[k][i] = (*A->_Sol[indexA[k]])(dof);
        b[k][i] = (*B->_Sol[indexB[k]])(dof);
      }
    }

    const unsigned nGauss = msh->_finiteElement[ielt][solType]->GetGaussPointNumber();

    for(unsigned ig = 0; ig < nGauss; ++ig) {
      double weight = 0.0;

      msh->_finiteElement[ielt][solType]->Jacobian(vx, ig, weight, phi, gradPhi);

      double pointError2 = 0.0;

      for(unsigned k = 0; k < nFields; ++k) {
        double diff = 0.0;

        for(unsigned i = 0; i < nDofs; ++i) {
          diff += phi[i] * (a[k][i] - b[k][i]);
        }

        pointError2 += diff * diff;
      }

      localError2 += pointError2 * weight;
    }
  }

  double globalError2 = 0.0;

  MPI_Allreduce(&localError2, &globalError2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return std::sqrt(globalError2);
}

enum class ConvergenceType {
  Temporal,
  Spatial
};

struct ConvergenceSnapshot {
  unsigned resolution;   // nSteps for temporal, levelC for spatial
  double step;           // dt for temporal, h for spatial

  LevelSetDiagnostics diagnostics;
  std::unique_ptr<MultiLevelSolution> solution;
};

void PrintConvergence(
  std::vector<ConvergenceSnapshot>& snapshots,
  const std::string& psiName,
  const std::vector<std::string>& vName,
  const unsigned referenceLevel,
  const unsigned iproc,
  const ConvergenceType type,
  const SimulationCase simulation) {

  const unsigned nSnapshots = snapshots.size();

  if(nSnapshots < 3) {
    if(iproc == 0) {
      std::cerr << "Need at least three solutions for convergence analysis." << std::endl;
    }
    return;
  }

  std::vector<ConvergenceSnapshot*> ordered(nSnapshots);

  for(unsigned k = 0; k < nSnapshots; ++k) ordered[k] = &snapshots[k];

  std::sort(ordered.begin(), ordered.end(), [](const ConvergenceSnapshot * a, const ConvergenceSnapshot * b) {
    return a->step > b->step;
  });

  const bool hasRisingBubbleMetrics = (simulation == SimulationCase::RB1) || (simulation == SimulationCase::RB2);

  unsigned verticalDirection = 0;

  if(hasRisingBubbleMetrics) {
    if(ordered[0]->diagnostics.barycenter.empty()) throw std::runtime_error("PrintConvergence: empty barycenter");
    verticalDirection = ordered[0]->diagnostics.barycenter.size() - 1u;
  }

  // -------------------------------------------------------
  // Errors
  // -------------------------------------------------------

  std::vector<double> errorPsi(nSnapshots - 1u, 0.);
  std::vector<double> errorVel(nSnapshots - 1u, 0.);

  std::vector<double> errorBarycenter;
  std::vector<double> errorRisingVelocity;
  std::vector<double> errorCircularity;

  if(hasRisingBubbleMetrics) {
    errorBarycenter.resize(nSnapshots - 1u, 0.);
    errorRisingVelocity.resize(nSnapshots - 1u, 0.);
    errorCircularity.resize(nSnapshots - 1u, 0.);
  }

  // -------------------------------------------------------
  // Orders
  // -------------------------------------------------------

  const double nan = std::numeric_limits<double>::quiet_NaN();

  std::vector<double> orderPsi(nSnapshots - 1u, nan);
  std::vector<double> orderVel(nSnapshots - 1u, nan);

  std::vector<double> orderBarycenter;
  std::vector<double> orderRisingVelocity;
  std::vector<double> orderCircularity;

  if(hasRisingBubbleMetrics) {
    orderBarycenter.resize(nSnapshots - 1u, nan);
    orderRisingVelocity.resize(nSnapshots - 1u, nan);
    orderCircularity.resize(nSnapshots - 1u, nan);
  }

  // -------------------------------------------------------
  // Compute differences between consecutive resolutions
  // -------------------------------------------------------

  for(unsigned k = 0; k + 1u < nSnapshots; ++k) {

    errorPsi[k] = ComputeL2Error(*ordered[k]->solution, *ordered[k + 1u]->solution, {psiName}, referenceLevel);

    errorVel[k] = ComputeL2Error(*ordered[k]->solution, *ordered[k + 1u]->solution, vName, referenceLevel);

    if(hasRisingBubbleMetrics) {

      const double barycenterA = ordered[k]->diagnostics.barycenter[verticalDirection];
      const double barycenterB = ordered[k + 1u]->diagnostics.barycenter[verticalDirection];

      const double risingVelocityA = ordered[k]->diagnostics.meanVelocity[verticalDirection];
      const double risingVelocityB = ordered[k + 1u]->diagnostics.meanVelocity[verticalDirection];

      const double circularityA = ordered[k]->diagnostics.circularity;
      const double circularityB = ordered[k + 1u]->diagnostics.circularity;

      errorBarycenter[k] = std::abs(barycenterA - barycenterB);
      errorRisingVelocity[k] = std::abs(risingVelocityA - risingVelocityB);
      errorCircularity[k] = std::abs(circularityA - circularityB);
    }
  }

  // -------------------------------------------------------
  // Orders
  // -------------------------------------------------------

  auto ComputeOrder = [](const double e0, const double e1, const double ratio) {
    if(e0 <= 0. || e1 <= 0. || ratio <= 0. || std::abs(std::log(ratio)) < 1.e-14) return std::numeric_limits<double>::quiet_NaN();
    return std::log(e0 / e1) / std::log(ratio);
  };

  for(unsigned k = 0; k + 2u < nSnapshots; ++k) {

    const double ratio = ordered[k]->step / ordered[k + 1u]->step;

    orderPsi[k] = ComputeOrder(errorPsi[k], errorPsi[k + 1u], ratio);
    orderVel[k] = ComputeOrder(errorVel[k], errorVel[k + 1u], ratio);

    if(hasRisingBubbleMetrics) {
      orderBarycenter[k] = ComputeOrder(errorBarycenter[k], errorBarycenter[k + 1u], ratio);
      orderRisingVelocity[k] = ComputeOrder(errorRisingVelocity[k], errorRisingVelocity[k + 1u], ratio);
      orderCircularity[k] = ComputeOrder(errorCircularity[k], errorCircularity[k + 1u], ratio);
    }
  }

  if(iproc != 0) return;

  // -------------------------------------------------------
  // Print
  // -------------------------------------------------------

  const std::string title = (type == ConvergenceType::Temporal) ? "TEMPORAL SELF-CONVERGENCE" : "SPATIAL SELF-CONVERGENCE";

  const std::string pairName = (type == ConvergenceType::Temporal) ? "NS pair" : "LC pair";

  const std::string stepName = (type == ConvergenceType::Temporal) ? "dt" : "h";

  const unsigned wPair = 14;
  const unsigned wStep = 14;
  const unsigned wError = 18;
  const unsigned wOrder = 10;

  std::cout << '\n';

  std::cout << "========================================================================================";

  if(hasRisingBubbleMetrics) std::cout << "==============================================================";

  std::cout << '\n';

  std::cout << std::setw(60) << title << '\n';

  std::cout << "========================================================================================";

  if(hasRisingBubbleMetrics) std::cout << "==============================================================";

  std::cout << '\n';

  std::cout << std::setw(wPair) << pairName
            << std::setw(wStep) << stepName
            << std::setw(wError) << "L2(Psi)"
            << std::setw(wOrder) << "p(Psi)"
            << std::setw(wError) << "L2(V)"
            << std::setw(wOrder) << "p(V)";

  if(hasRisingBubbleMetrics) {
    std::cout << std::setw(wError) << "dBarycenter"
              << std::setw(wOrder) << "p(B)"
              << std::setw(wError) << "dRisingVel"
              << std::setw(wOrder) << "p(RV)"
              << std::setw(wError) << "dCircularity"
              << std::setw(wOrder) << "p(C)";
  }

  std::cout << '\n';

  for(unsigned k = 0; k + 1u < nSnapshots; ++k) {

    const std::string pair = std::to_string(ordered[k]->resolution) + "-" + std::to_string(ordered[k + 1u]->resolution);

    std::cout << std::setw(wPair) << pair;

    std::cout << std::setw(wStep) << std::scientific << std::setprecision(5) << ordered[k]->step;

    auto PrintErrorAndOrder = [&](const double error, const double order) {
      std::cout << std::setw(wError) << std::scientific << std::setprecision(8) << error;

      if(std::isfinite(order)) {
        std::cout << std::setw(wOrder) << std::fixed << std::setprecision(4) << order;
      }
      else {
        std::cout << std::setw(wOrder) << "-";
      }
    };

    PrintErrorAndOrder(errorPsi[k], orderPsi[k]);
    PrintErrorAndOrder(errorVel[k], orderVel[k]);

    if(hasRisingBubbleMetrics) {
      PrintErrorAndOrder(errorBarycenter[k], orderBarycenter[k]);
      PrintErrorAndOrder(errorRisingVelocity[k], orderRisingVelocity[k]);
      PrintErrorAndOrder(errorCircularity[k], orderCircularity[k]);
    }

    std::cout << '\n';
  }

  std::cout << "========================================================================================";

  if(hasRisingBubbleMetrics) std::cout << "==============================================================";

  std::cout << '\n';
}

struct SimulationConfig {

  std::vector<double> x_min;
  std::vector<double> x_max;
  std::vector<unsigned> n;

  double period = 0.;

  std::vector<double> center;
  double radius = 0.;

  MultiphasePhysicalProperties properties;
};

SimulationConfig GetSimulationConfig(
  const SimulationCase simulation,
  const unsigned dim) {

  SimulationConfig cfg;

  cfg.x_min.resize(3);
  cfg.x_max.resize(3);
  cfg.n.resize(3);

  switch(simulation) {

    case SimulationCase::PAR: {

      cfg.x_min[0] = -0.5;
      cfg.x_max[0] =  0.5;
      cfg.x_min[1] = -0.5;
      cfg.x_max[1] =  0.5;
      cfg.x_min[2] = (dim == 3) ? -0.5 : 0.;
      cfg.x_max[2] = (dim == 3) ? 0.5 : 0.;

      cfg.n[0] = 8;
      cfg.n[1] = 8;
      cfg.n[2] = (dim == 3) ? 8 : 0;

      cfg.period = 10.;

      cfg.center = {0.0, 0.0, 0.0};
      cfg.center.resize(dim);
      cfg.radius = 0.25;

      cfg.properties.rho1 = 1000.;
      cfg.properties.rho2 = 1.;
      cfg.properties.mu1  = 10.;
      cfg.properties.mu2  = 0.1;
      cfg.properties.sigma = 1.;
      cfg.properties.gravity = 0.;

      break;
    }

    case SimulationCase::SAYE1: {

      cfg.x_min[0] = -0.5;
      cfg.x_max[0] =  0.5;
      cfg.x_min[1] = -0.5;
      cfg.x_max[1] =  0.5;
      cfg.x_min[2] = (dim == 3) ? -0.5 : 0.;
      cfg.x_max[2] = (dim == 3) ? 0.5 : 0.;

      cfg.n[0] = 8;
      cfg.n[1] = 8;
      cfg.n[2] = (dim == 3) ? 8 : 0;

      cfg.period = 0.01;

      cfg.properties.mu1 = 1.;
      cfg.properties.mu2 = 1.;
      cfg.properties.rho1 = 1.;
      cfg.properties.rho2 = 1.;
      cfg.properties.sigma = 10.;
      cfg.properties.gravity = 0.;

      break;
    }

    case SimulationCase::SAYE2: {

      cfg.x_min[0] = -0.5;
      cfg.x_max[0] =  0.5;
      cfg.x_min[1] = -0.5;
      cfg.x_max[1] =  0.5;
      cfg.x_min[2] = (dim == 3) ? -0.5 : 0.;
      cfg.x_max[2] = (dim == 3) ?  0.5 : 0.;

      cfg.n[0] = 8;
      cfg.n[1] = 8;
      cfg.n[2] = (dim == 3) ? 8 : 0;

      cfg.period = 0.01;

      cfg.properties.mu1 = 1.;
      cfg.properties.mu2 = 1.;
      cfg.properties.rho1 = 1.;
      cfg.properties.rho2 = 1.;
      cfg.properties.sigma = 10.;
      cfg.properties.gravity = 0.;

      break;
    }

    case SimulationCase::RB1: {

      cfg.x_min[0] = 0.;
      cfg.x_max[0] = 1.;
      cfg.x_min[1] = 0.;
      cfg.x_max[1] = (dim == 3) ? 1. : 2.;
      cfg.x_min[2] = 0.0;
      cfg.x_max[2] = (dim == 3) ? 2. : 0.;

      cfg.n[0] = 8;
      cfg.n[1] = (dim == 3) ? 8 : 16;
      cfg.n[2] = (dim == 3) ? 16 : 0;

      cfg.period = 3.;

      cfg.center = {0.5, 0.5, 0.5};
      cfg.center.resize(dim);
      cfg.radius = 0.25;

      cfg.properties.rho2 = 1000.;
      cfg.properties.rho1 = 100.;
      cfg.properties.mu2  = 10.;
      cfg.properties.mu1  = 1.;
      cfg.properties.sigma = 24.5;
      cfg.properties.gravity = -0.98;

      break;
    }

    case SimulationCase::RB2: {

      cfg.x_min[0] = 0.;
      cfg.x_max[0] = 1.;
      cfg.x_min[1] = 0.;
      cfg.x_max[1] = (dim == 3) ? 1. : 2.;
      cfg.x_min[2] = 0.0;
      cfg.x_max[2] = (dim == 3) ? 2. : 0.;

      cfg.n[0] = 8;
      cfg.n[1] = (dim == 3) ? 8 : 16;
      cfg.n[2] = (dim == 3) ? 16 : 0;

      cfg.period = 3.;

      cfg.center = {0.5, 0.5, 0.5};
      cfg.center.resize(dim);
      cfg.radius = 0.25;

      cfg.properties.rho2 = 1000.;
      cfg.properties.rho1 = 1.;
      cfg.properties.mu2  = 10.;
      cfg.properties.mu1  = 0.1;
      cfg.properties.sigma = 1.96;
      cfg.properties.gravity = -0.98;

      break;
    }
  }

  return cfg;
}