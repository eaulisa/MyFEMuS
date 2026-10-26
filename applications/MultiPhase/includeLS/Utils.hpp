#pragma once

#include "CutFemWeightNew.hpp"
#include "CDWeights.hpp"
#include "MultiLevelProblem.hpp"
typedef double TypeIO;
typedef cpp_bin_float_oct TypeA;
//typedef double TypeA;
typedef cpp_bin_float_oct oct;

CutFemWeight <TypeIO, TypeA> hex  = CutFemWeight<TypeIO, TypeA >(HEX, 5, "legendre");
CutFemWeight <TypeIO, TypeA> tet  = CutFemWeight<TypeIO, TypeA >(TET, 5, "legendre");
CutFemWeight <TypeIO, TypeA> wedge = CutFemWeight<TypeIO, TypeA >(WEDGE, 5, "legendre");
CutFemWeight <TypeIO, TypeA> quad = CutFemWeight<TypeIO, TypeA >(QUAD, 5, "legendre");
CutFemWeight <TypeIO, TypeA> tri = CutFemWeight<TypeIO, TypeA >(TRI, 5, "legendre");

const std::vector< CutFemWeight <TypeIO, TypeA> *> cfw = {&hex, &tet, &wedge, &quad, &tri};

// unsigned qM = 5;
// double dx = .5;
// double dtetha = 30.;
// CDWeightHEX <TypeA> hexCD0(qM, dx, dtetha);
// CDWeightQUAD <TypeA> quadCD0(qM, dx, dtetha);
// CDWeightTRI <TypeA> triCD0(qM, dx, dtetha);
//
// const std::vector< CDWeight <TypeA> *> cfCDw0 = {&hexCD0, &quadCD0, &quadCD0, &quadCD0, &triCD0};

Fem fem = Fem(hex.GetGaussQuadratureOrder(), hex.GetDimension());

void RungeKutta4(std::vector<MyVector<double>> &X,
                 MultiLevelSolution & mlSol,
                 BBoxToIel & bbox,
                 const std::vector<double> x_min,
                 const std::vector<double> x_max,
                 const std::vector<std::string> &velName,
                 const unsigned vlevel,
                 const double dt);

void rkStep(MultiLevelSolution & mlSol,
            BBoxToIel & bbox,
            const std::vector<double> x_min,
            const std::vector<double> x_max,
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
                     const std::vector<double> x_min,
                     const std::vector<double> x_max,
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
                     const std::vector<double> x_min,
                     const std::vector<double> x_max,
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

  const unsigned dim = X1.size();
  const double tol = 1.e-12;

  for(unsigned i = X1[0].begin(); i < X1[0].end(); ++i) {
    for(unsigned d = 0; d < dim; ++d) {
      if(X1[d][i] <= x_min[d] + tol && X1[d][i] >= x_min[d] - tol) X1[d][i] = x_min[d];
      else if(X1[d][i] >= x_max[d] - tol && X1[d][i] <= x_max[d] + tol) X1[d][i] = x_max[d];
    }
  }

  if(fabs(dt) > 1.0e-10) RungeKutta4(X1, mlSol0, bbox, x_min, x_max, vName, vLevel, dt);

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

    // solVec1->zero();
    unsigned offset = psiProjected.begin();
    for (unsigned i = psiProjected.begin(); i < psiProjected.end(); ++i) {
      if (isInsideDomain[i - offset]) {
        solVec1->set(i, psiProjected[i]);
      }
      else { // TODO add boundarycondition for psi
        // throw std::runtime_error("POINT OUTSIDE THE DOMAIN");
        std::cerr << std::setprecision(17);
        std::cerr << "POINT OUTSIDE THE DOMAIN" << std::endl;
        std::cerr << X1[0][i] << " " << X1[1][i] << std::endl;
        abort();
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
                 const std::vector<double> x_min,
                 const std::vector<double> x_max,
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
  std::vector<MyVector<double>> Xold = X;
  std::vector<std::vector<MyVector<double>>> K;
  for(unsigned rk = 0; rk < rk_nsteps; rk++) {
    rkStep(mlSol, bbox, x_min, x_max, X, K, rk, velName, vLevel, dt, c[rk], a[rk]);
  }
  for(unsigned rk = 0; rk < rk_nsteps; rk++) {
    for(unsigned d = 0; d < dim; d++) {
      for(unsigned i = X[d].begin(); i < X[d].end(); i++) {
        X[d][i] += dt * K[rk][d][i] * b[rk];
      }
    }
  }
  const double boundaryTol = 1.e-12;

  for(unsigned i = X[0].begin();
      i < X[0].end();
      ++i) {

    for(unsigned d = 0; d < dim; ++d) {
      if(std::abs(X[d][i] - x_min[d]) < boundaryTol) {

        X[d][i] = x_min[d];
      }

      else if(std::abs(X[d][i] - x_max[d]) < boundaryTol) {

        X[d][i] = x_max[d];
      }
    }
  }

}

void rkStep(MultiLevelSolution & mlSol,
            BBoxToIel & bbox,
            const std::vector<double> x_min,
            const std::vector<double> x_max,
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

      const double boundaryTol = 1.e-12;

      for (unsigned k = 0; k < X.size(); ++k) {

        if (std::abs(Xk[k][i] - x_min[k]) < boundaryTol) {

          Xk[k][i] = x_min[k];
        }

        else if (std::abs(Xk[k][i] - x_max[k]) < boundaryTol) {

          Xk[k][i] = x_max[k];
        }
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

double FindLevelSetZeroOnNodalLine(
  Mesh* msh,
  Solution* sol,
  const unsigned psiIndex,
  const unsigned psiType,
  const std::vector<double>& lineCoordinates,
  const double coordTol = 1.e-12,
  const double psiTol = 1.e-14) {

  const unsigned dim = msh->GetDimension();
  const unsigned iproc = msh->processor_id();
  const unsigned xType = 2;

  if(lineCoordinates.size() != dim - 1u) {
    throw std::runtime_error("FindLevelSetZeroOnNodalLine: lineCoordinates must have dim - 1 entries");
  }

  struct LineNode {
    unsigned dof;
    double xVertical;
    double psi;
  };

  std::vector<unsigned> localDof;
  std::vector<double> localX;
  std::vector<double> localPsi;

  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1u]; ++iel) {

    const unsigned nDofsPsi = msh->GetElementDofNumber(iel, psiType);
    const unsigned nDofsX = msh->GetElementDofNumber(iel, xType);

    if(nDofsPsi != nDofsX) {
      throw std::runtime_error("FindLevelSetZeroOnNodalLine: psi and coordinate spaces must have the same nodal structure");
    }

    for(unsigned i = 0; i < nDofsPsi; ++i) {

      const unsigned psiDof = msh->GetSolutionDof(i, iel, psiType);
      const unsigned xDof = msh->GetSolutionDof(i, iel, xType);

      bool onLine = true;

      for(unsigned d = 0; d < dim - 1u; ++d) {

        const double xd = (*msh->_topology->_Sol[d])(xDof);

        if(std::abs(xd - lineCoordinates[d]) > coordTol) {
          onLine = false;
        }
      }

      if(!onLine) {
        continue;
      }

      const double xVertical = (*msh->_topology->_Sol[dim - 1u])(xDof);
      const double psiValue = (*sol->_Sol[psiIndex])(psiDof);

      localDof.push_back(psiDof);
      localX.push_back(xVertical);
      localPsi.push_back(psiValue);
    }
  }

  int nprocs = 1;

  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  const int localSize = static_cast<int>(localDof.size());

  std::vector<int> counts(nprocs, 0);
  std::vector<int> displs(nprocs, 0);

  MPI_Allgather(&localSize, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

  int globalSize = 0;

  for(int p = 0; p < nprocs; ++p) {
    displs[p] = globalSize;
    globalSize += counts[p];
  }

  std::vector<unsigned> globalDof(globalSize);
  std::vector<double> globalX(globalSize);
  std::vector<double> globalPsi(globalSize);

  MPI_Allgatherv(localDof.data(), localSize, MPI_UNSIGNED, globalDof.data(), counts.data(), displs.data(), MPI_UNSIGNED, MPI_COMM_WORLD);

  MPI_Allgatherv(localX.data(), localSize, MPI_DOUBLE, globalX.data(), counts.data(), displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);

  MPI_Allgatherv(localPsi.data(), localSize, MPI_DOUBLE, globalPsi.data(), counts.data(), displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);

  std::vector<LineNode> nodes(globalSize);

  for(int i = 0; i < globalSize; ++i) {
    nodes[i].dof = globalDof[i];
    nodes[i].xVertical = globalX[i];
    nodes[i].psi = globalPsi[i];
  }

  std::sort(nodes.begin(), nodes.end(), [](const LineNode & a, const LineNode & b) {
    return a.dof < b.dof;
  });

  nodes.erase(std::unique(nodes.begin(), nodes.end(), [](const LineNode & a, const LineNode & b) {
    return a.dof == b.dof;
  }), nodes.end());

  std::sort(nodes.begin(), nodes.end(), [](const LineNode & a, const LineNode & b) {
    return a.xVertical < b.xVertical;
  });

  if(nodes.size() < 2u) {
    throw std::runtime_error("FindLevelSetZeroOnNodalLine: fewer than two nodes found on line");
  }

  std::vector<double> roots;

  for(const LineNode& node : nodes) {
    if(std::abs(node.psi) <= psiTol) {
      roots.push_back(node.xVertical);
    }
  }

  for(unsigned i = 0; i + 1u < nodes.size(); ++i) {

    const double psi0 = nodes[i].psi;
    const double psi1 = nodes[i + 1u].psi;

    if(psi0 * psi1 < 0.0) {

      const double x0 = nodes[i].xVertical;
      const double x1 = nodes[i + 1u].xVertical;

      const double root = x0 - psi0 * (x1 - x0) / (psi1 - psi0);

      roots.push_back(root);
    }
  }

  std::sort(roots.begin(), roots.end());

  roots.erase(std::unique(roots.begin(), roots.end(), [coordTol](const double a, const double b) {
    return std::abs(a - b) <= coordTol;
  }), roots.end());

  if(roots.empty()) {
    throw std::runtime_error("FindLevelSetZeroOnNodalLine: no zero crossing found");
  }

  if(roots.size() != 1u) {
    throw std::runtime_error("FindLevelSetZeroOnNodalLine: multiple zero crossings found");
  }

  return roots[0];
}

enum class SimulationCase {
  PAR,
  RB1,
  RB2,
  SAYE1,
  SAYE2,
  RT,
  DB,
  BB,
  MB,
  PI
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
    case SimulationCase::RT:
      return "rt";
    case SimulationCase::DB:
      return "rt";
    case SimulationCase::BB:
      return "bb";
    case SimulationCase::MB:
      return "mb";
    case SimulationCase::PI:
      return "pi";
  }

  throw std::runtime_error("Invalid SimulationCase");
}

SimulationCase ParseSimulationCase(const std::string& s) {

  if(s == "par")  return SimulationCase::PAR;
  if(s == "rb1")  return SimulationCase::RB1;
  if(s == "rb2")  return SimulationCase::RB2;
  if(s == "saye1") return SimulationCase::SAYE1;
  if(s == "saye2") return SimulationCase::SAYE2;
  if(s == "rt") return SimulationCase::RT;
  if(s == "db") return SimulationCase::DB;
  if(s == "bb") return SimulationCase::BB;
  if(s == "mb") return SimulationCase::MB;
  if(s == "pi") return SimulationCase::PI;

  throw std::runtime_error(
    "Unknown simulation '" + s +
    "'. Valid values are: par, rb1, rb2, saye"
  );
}

enum class BoundaryType : unsigned {
  DIRICHLET_S,
  NEUMANN_S,
  WALL,
  SLIP,
  OUTFLOW
};

struct BoundaryCondition {

  BoundaryType type = BoundaryType::NEUMANN_S;
  double value = 0.;
};

struct BoundaryConfig {

  std::vector<BoundaryCondition> min;
  std::vector<BoundaryCondition> max;

  void Resize(const unsigned dim) {

    min.resize(dim);
    max.resize(dim);
  }
};

struct SimulationConfig {

  std::vector<double> x_min;
  std::vector<double> x_max;
  std::vector<unsigned> n;

  double period = 0.;

  PsiConfig psi;

  MultiphasePhysicalProperties properties;

  BoundaryConfig boundary;

  bool symmetric = false;
  bool axisymmetric = false;
};

struct LevelSetDiagnostics {

  // ============================================================
  // Common metrics
  // ============================================================

  double innerArea = 0.;
  double outerArea = 0.;
  double totalArea = 0.;
  double interfaceLength = 0.;

  // ============================================================
  // Rising bubble
  // ============================================================

  bool hasRisingBubbleMetrics = false;
  std::vector<double> barycenter;
  std::vector<double> meanVelocity;
  double circularity = 0.;

  // ============================================================
  // Rayleigh-Taylor
  // ============================================================

  bool hasRayleighTaylorMetrics = false;
  double pb = 0.;
  double ps = 0.;
  double psad = 0.;

  // ============================================================
  // Phase inversion
  // ============================================================

  bool hasPhaseInversionMetrics = false;

  double fluid1UpperVolume = 0.;
  double fluid1UpperRatio = 0.;
  double fluid1VolumeRatio = 0.;

  double potentialEnergy1 = 0.;
  double potentialEnergy2 = 0.;

  double kineticEnergy1 = 0.;
  double kineticEnergy2 = 0.;

  double enstrophy1 = 0.;
  double enstrophy2 = 0.;

  double mechanicalEnergy = 0.;

  double potentialEnergy1Star = 0.;
  double potentialEnergy2Star = 0.;

  double kineticEnergy1Star = 0.;
  double kineticEnergy2Star = 0.;

  double interfaceAreaStar = 0.;
};

LevelSetDiagnostics ComputeLevelSetDiagnostics(MultiLevelSolution& mlSol, const std::string& psiName, const std::vector<std::string>& velocityName, const std::vector<std::string>& pressureName, const unsigned velocityLevel, const SimulationCase simulation, const SimulationConfig& cfg, const bool symmetric, const bool axisymmetric, const double axisCoordinate, const std::vector<double>& symmetryCoordinates) {

  MultiLevelMesh& mlMsh = *mlSol.GetMultilevelMesh();
  const unsigned level = mlMsh.GetNumberOfLevels() - 1u;
  Mesh* msh = mlMsh.GetLevel(level);
  Solution* sol = mlSol.GetSolutionLevel(level);
  const unsigned dim = msh->GetDimension();
  const unsigned iproc = msh->processor_id();
  const unsigned psiIndex = mlSol.GetIndex(psiName.c_str());
  const unsigned psiType = mlSol.GetSolutionType(psiIndex);
  const unsigned xType = 2;
  const double pi = std::acos(-1.0);

  if(dim != 2 && dim != 3) {
    throw std::runtime_error("ComputeLevelSetDiagnostics: requires dim = 2 or 3");
  }

  if(symmetric && axisymmetric) {
    throw std::runtime_error("ComputeLevelSetDiagnostics: symmetric and axisymmetric cannot both be true");
  }

  if(axisymmetric && dim != 2) {
    throw std::runtime_error("ComputeLevelSetDiagnostics: axisymmetric diagnostics require dim = 2");
  }

  const double symmetryFactor = symmetric ? ((dim == 2) ? 2.0 : 4.0) : 1.0;
  const bool computeRisingBubbleMetrics = simulation == SimulationCase::RB1 || simulation == SimulationCase::RB2;
  const bool computeRayleighTaylorMetrics = simulation == SimulationCase::RT;
  const bool computePhaseInversionMetrics = simulation == SimulationCase::PI;

  if(computePhaseInversionMetrics && dim != 3) {
    throw std::runtime_error("ComputeLevelSetDiagnostics: PI requires dim = 3");
  }

  if(computePhaseInversionMetrics && (symmetric || axisymmetric)) {
    throw std::runtime_error("ComputeLevelSetDiagnostics: PI diagnostics do not support symmetry reduction");
  }

  std::vector<unsigned> velocityIndex;
  unsigned velocityType = 0;
  std::vector<unsigned> pressureIndex;
  unsigned pressureType = 0;

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
    throw std::runtime_error("ComputeLevelSetDiagnostics: pressureName must contain P1 and P2");
  }

  pressureIndex.resize(2);

  for(unsigned p = 0; p < 2; ++p) {
    pressureIndex[p] = mlSol.GetIndex(pressureName[p].c_str());
  }

  pressureType = mlSol.GetSolutionType(pressureIndex[0]);

  if(mlSol.GetSolutionType(pressureIndex[1]) != pressureType) {
    throw std::runtime_error("ComputeLevelSetDiagnostics: pressure fields must have the same solution type");
  }

  if(velocityLevel > level) {
    throw std::runtime_error("ComputeLevelSetDiagnostics: velocityLevel cannot be finer than diagnostics level");
  }

  for(unsigned l = velocityLevel; l < level; ++l) {
    Solution* sol_l = mlSol.GetSolutionLevel(l);
    Solution* sol_lp1 = mlSol.GetSolutionLevel(l + 1u);
    Mesh* msh_lp1 = mlMsh.GetLevel(l + 1u);
    SparseMatrix* Pv = msh_lp1->GetCoarseToFineProjection(velocityType);
    SparseMatrix* Pp = msh_lp1->GetCoarseToFineProjection(pressureType);

    for(unsigned d = 0; d < dim; ++d) {
      sol_lp1->_Sol[velocityIndex[d]]->matrix_mult(*(sol_l->_Sol[velocityIndex[d]]), *Pv);
    }

    for(unsigned p = 0; p < 2; ++p) {
      sol_lp1->_Sol[pressureIndex[p]]->matrix_mult(*(sol_l->_Sol[pressureIndex[p]]), *Pp);
    }
  }

  const unsigned verticalDir = (dim == 2) ? 1u : 2u;

  double phaseInversionH = 0.;
  double phaseInversionTopMin = 0.;
  double phaseInversionInitialVolume = 0.;
  double phaseInversionInitialInterfaceArea = 0.;
  double rho1 = 0.;
  double rho2 = 0.;
  double gravity = 0.;
  double potentialEnergy1Final = 0.;
  double potentialEnergy2Final = 0.;
  double kineticEnergy1Scale = 0.;
  double kineticEnergy2Scale = 0.;

  if(computePhaseInversionMetrics) {
    phaseInversionH = cfg.x_max[0] - cfg.x_min[0];

    for(unsigned d = 1; d < 3; ++d) {
      const double Hd = cfg.x_max[d] - cfg.x_min[d];

      if(std::abs(Hd - phaseInversionH) > 1.e-12 * phaseInversionH) {
        throw std::runtime_error("ComputeLevelSetDiagnostics: PI requires a cubic domain");
      }
    }

    rho1 = cfg.properties.rho1;
    rho2 = cfg.properties.rho2;
    gravity = std::abs(cfg.properties.gravity);

    if(rho1 <= 0. || rho2 <= 0. || gravity <= 0.) {
      throw std::runtime_error("ComputeLevelSetDiagnostics: invalid PI physical properties");
    }

    if(rho2 <= rho1) {
      throw std::runtime_error("ComputeLevelSetDiagnostics: PI requires rho2 > rho1");
    }

    phaseInversionTopMin = cfg.x_max[verticalDir] - phaseInversionH / 8.;
    phaseInversionInitialVolume = phaseInversionH * phaseInversionH * phaseInversionH / 8.;
    phaseInversionInitialInterfaceArea = 3. * phaseInversionH * phaseInversionH / 4.;
    const double H4 = phaseInversionH * phaseInversionH * phaseInversionH * phaseInversionH;
    potentialEnergy1Final = 15. * rho1 * gravity * H4 / 128.;
    potentialEnergy2Final = 49. * rho2 * gravity * H4 / 128.;
    const double Ug = (rho2 - rho1) / rho1 * std::sqrt(phaseInversionH * gravity / 2.);
    const double H3 = phaseInversionH * phaseInversionH * phaseInversionH;
    kineticEnergy1Scale = rho1 * Ug * Ug * H3 / 16.;
    kineticEnergy2Scale = rho2 * Ug * Ug * H3 / 16.;
  }

  double innerAreaLocal = 0.;
  double outerAreaLocal = 0.;
  double totalAreaLocal = 0.;
  double interfaceLengthLocal = 0.;
  std::vector<double> barycenterIntegralLocal(dim, 0.);
  std::vector<double> velocityIntegralLocal(dim, 0.);
  double fluid1UpperVolumeLocal = 0.;
  double potentialEnergy1Local = 0.;
  double potentialEnergy2Local = 0.;
  double kineticEnergy1Local = 0.;
  double kineticEnergy2Local = 0.;
  double enstrophy1Local = 0.;
  double enstrophy2Local = 0.;

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

  auto GetPhysicalMeasureFactor = [&](const std::vector<double>& x) -> double {
    if(axisymmetric) {
    const double r = x[0] - axisCoordinate;

      if(r < -1.e-12) {
        throw std::runtime_error("ComputeLevelSetDiagnostics: negative axisymmetric radius");
      }

      const double radius = (r > 0.) ? r : 0.;
      return 2. * pi * radius;
    }

    if(symmetric) {
    return symmetryFactor;
  }

  return 1.;
};

for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; ++iel) {
    const unsigned ielGeom = msh->GetElementType(iel);
    const unsigned nDofsPsi = msh->GetElementDofNumber(iel, psiType);
    const unsigned nDofsX = msh->GetElementDofNumber(iel, xType);
    const unsigned nDofsVelocity = msh->GetElementDofNumber(iel, velocityType);

    psi.resize(nDofsPsi);

    for(unsigned d = 0; d < dim; ++d) {
      coordX[d].resize(nDofsX);
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
    const elem_type* femVelocity = fem.GetFiniteElement(ielGeom, velocityType);
    const unsigned nGauss = femPsi->GetGaussPointNumber();

    if(nGauss != cfw[ielGeom]->GetGaussQuadraturePointNumber()) {
      throw std::runtime_error("ComputeLevelSetDiagnostics: incompatible quadrature");
    }

    if(femVelocity->GetGaussPointNumber() != nGauss) {
      throw std::runtime_error("ComputeLevelSetDiagnostics: incompatible velocity quadrature");
    }

    std::vector<double> psiG(nGauss, 0.);
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

    const bool cut = psiMin <= 0. && psiMax >= 0. && psiMax - psiMin > 1.e-14;

    auto EvaluateVelocity = [&](const unsigned ig, std::vector<double>& u, double& omega2) {
      double weightVelocity = 0.;
      femVelocity->Jacobian(coordX, ig, weightVelocity, phiVelocity, phiVelocity_x);
      u.assign(dim, 0.);
      std::vector<std::vector<double>> gradU(dim, std::vector<double>(dim, 0.));

      for(unsigned i = 0; i < nDofsVelocity; ++i) {
        for(unsigned d = 0; d < dim; ++d) {
          u[d] += velocity[d][i] * phiVelocity[i];

          for(unsigned j = 0; j < dim; ++j) {
            gradU[d][j] += velocity[d][i] * phiVelocity_x[i * dim + j];
          }
        }
      }

      omega2 = 0.;

      if(dim == 2) {
        const double omega = gradU[1][0] - gradU[0][1];
        omega2 = omega * omega;
      }
      else {
        const double omegaX = gradU[2][1] - gradU[1][2];
        const double omegaY = gradU[0][2] - gradU[2][0];
        const double omegaZ = gradU[1][0] - gradU[0][1];
        omega2 = omegaX * omegaX + omegaY * omegaY + omegaZ * omegaZ;
      }
    };

    if(!cut) {
      for(unsigned ig = 0; ig < nGauss; ++ig) {
        double weight = 0.;
        double weightX = 0.;
        femPsi->Jacobian(coordX, ig, weight, phiPsi, phiPsi_x);
        femX->Jacobian(coordX, ig, weightX, phiX, phiX_x);
        std::vector<double> x(dim, 0.);

        for(unsigned i = 0; i < nDofsX; ++i) {
          for(unsigned d = 0; d < dim; ++d) {
            x[d] += coordX[d][i] * phiX[i];
          }
        }

        const double measureFactor = GetPhysicalMeasureFactor(x);
        const double physicalWeight = weight * measureFactor;
        const bool inside = psiMax > 0.;

        totalAreaLocal += physicalWeight;

        if(inside) {
          innerAreaLocal += physicalWeight;
        }
        else {
          outerAreaLocal += physicalWeight;
        }

        if(computeRisingBubbleMetrics && inside) {
          std::vector<double> u;
          double omega2 = 0.;
          EvaluateVelocity(ig, u, omega2);

          for(unsigned d = 0; d < dim; ++d) {
            barycenterIntegralLocal[d] += x[d] * physicalWeight;
            velocityIntegralLocal[d] += u[d] * physicalWeight;
          }
        }

        if(computePhaseInversionMetrics) {
          std::vector<double> u;
          double omega2 = 0.;
          EvaluateVelocity(ig, u, omega2);
          double u2 = 0.;

          for(unsigned d = 0; d < dim; ++d) {
            u2 += u[d] * u[d];
          }

          const double height = x[verticalDir] - cfg.x_min[verticalDir];

          if(inside) {
            potentialEnergy1Local += rho1 * gravity * height * physicalWeight;
            kineticEnergy1Local += 0.5 * rho1 * u2 * physicalWeight;
            enstrophy1Local += 0.5 * omega2 * physicalWeight;

            if(x[verticalDir] >= phaseInversionTopMin) {
              fluid1UpperVolumeLocal += physicalWeight;
            }
          }
          else {
            potentialEnergy2Local += rho2 * gravity * height * physicalWeight;
            kineticEnergy2Local += 0.5 * rho2 * u2 * physicalWeight;
            enstrophy2Local += 0.5 * omega2 * physicalWeight;
          }
        }
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

    std::vector<TypeIO> weightInner(cfw[ielGeom]->GetGaussQuadraturePointNumber(), 0.);
    std::vector<TypeIO> weightOuter(cfw[ielGeom]->GetGaussQuadraturePointNumber(), 0.);
    std::vector<TypeIO> weightInterface(cfw[ielGeom]->GetGaussQuadraturePointNumber(), 0.);

    cfw[ielGeom]->GetWeightWithMap(0, a, interfaceConstant, weightInner);

    for(unsigned i = 0; i < weightInner.size(); ++i) {
      weightOuter[i] = 1. - weightInner[i];
    }

    cfw[ielGeom]->GetWeightWithMap(-1, a, interfaceConstant, weightInterface);

    for(unsigned ig = 0; ig < nGauss; ++ig) {
      double weight = 0.;
      double weightX = 0.;
      femPsi->Jacobian(coordX, ig, weight, phiPsi, phiPsi_x);
      femX->Jacobian(coordX, ig, weightX, phiX, phiX_x);
      std::vector<double> x(dim, 0.);

      for(unsigned i = 0; i < nDofsX; ++i) {
        for(unsigned d = 0; d < dim; ++d) {
          x[d] += coordX[d][i] * phiX[i];
        }
      }

      const double measureFactor = GetPhysicalMeasureFactor(x);
      const double physicalWeight = weight * measureFactor;
      const double innerPhysicalWeight = weight * static_cast<double>(weightInner[ig]) * measureFactor;
      const double outerPhysicalWeight = weight * static_cast<double>(weightOuter[ig]) * measureFactor;

      totalAreaLocal += physicalWeight;
      innerAreaLocal += innerPhysicalWeight;
      outerAreaLocal += outerPhysicalWeight;

      if(computeRisingBubbleMetrics) {
        std::vector<double> u;
        double omega2 = 0.;
        EvaluateVelocity(ig, u, omega2);

        for(unsigned d = 0; d < dim; ++d) {
          barycenterIntegralLocal[d] += x[d] * innerPhysicalWeight;
          velocityIntegralLocal[d] += u[d] * innerPhysicalWeight;
        }
      }

      if(computePhaseInversionMetrics) {
        std::vector<double> u;
        double omega2 = 0.;
        EvaluateVelocity(ig, u, omega2);
        double u2 = 0.;

        for(unsigned d = 0; d < dim; ++d) {
          u2 += u[d] * u[d];
        }

        const double height = x[verticalDir] - cfg.x_min[verticalDir];

        potentialEnergy1Local += rho1 * gravity * height * innerPhysicalWeight;
        potentialEnergy2Local += rho2 * gravity * height * outerPhysicalWeight;
        kineticEnergy1Local += 0.5 * rho1 * u2 * innerPhysicalWeight;
        kineticEnergy2Local += 0.5 * rho2 * u2 * outerPhysicalWeight;
        enstrophy1Local += 0.5 * omega2 * innerPhysicalWeight;
        enstrophy2Local += 0.5 * omega2 * outerPhysicalWeight;

        if(x[verticalDir] >= phaseInversionTopMin) {
          fluid1UpperVolumeLocal += innerPhysicalWeight;
        }
      }

      femPsi->GetJacobianMatrix(coordX, ig, weight, Jacob, JacI);
      std::vector<double> Nref(dim, 0.);
      double dsN2 = 0.;

      for(unsigned k = 0; k < dim; ++k) {
        for(unsigned j = 0; j < dim; ++j) {
          Nref[k] += JacI[j][k] * a[j];
        }

        dsN2 += Nref[k] * Nref[k];
      }

      const double dsN = std::sqrt(dsN2);
      const double dGamma = weight * static_cast<double>(weightInterface[ig]) * dsN * measureFactor;
      interfaceLengthLocal += dGamma;
    }
  }

  const unsigned basePhaseInversion = 4u + 2u * dim;
  const unsigned nGlobalQuantities = basePhaseInversion + 7u;
  std::vector<double> localSum(nGlobalQuantities, 0.);
  std::vector<double> globalSum(nGlobalQuantities, 0.);

  localSum[0] = innerAreaLocal;
  localSum[1] = outerAreaLocal;
  localSum[2] = totalAreaLocal;
  localSum[3] = interfaceLengthLocal;

  for(unsigned d = 0; d < dim; ++d) {
    localSum[4u + d] = barycenterIntegralLocal[d];
    localSum[4u + dim + d] = velocityIntegralLocal[d];
  }

  localSum[basePhaseInversion + 0u] = fluid1UpperVolumeLocal;
  localSum[basePhaseInversion + 1u] = potentialEnergy1Local;
  localSum[basePhaseInversion + 2u] = potentialEnergy2Local;
  localSum[basePhaseInversion + 3u] = kineticEnergy1Local;
  localSum[basePhaseInversion + 4u] = kineticEnergy2Local;
  localSum[basePhaseInversion + 5u] = enstrophy1Local;
  localSum[basePhaseInversion + 6u] = enstrophy2Local;

  MPI_Allreduce(localSum.data(), globalSum.data(), static_cast<int>(nGlobalQuantities), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  LevelSetDiagnostics result;

  result.innerArea = globalSum[0];
  result.outerArea = globalSum[1];
  result.totalArea = globalSum[2];
  result.interfaceLength = globalSum[3];

  result.hasRisingBubbleMetrics = computeRisingBubbleMetrics;

  if(computeRisingBubbleMetrics) {
    result.barycenter.assign(dim, 0.);
    result.meanVelocity.assign(dim, 0.);

    if(result.innerArea > 1.e-14) {
      for(unsigned d = 0; d < dim; ++d) {
        result.barycenter[d] = globalSum[4u + d] / result.innerArea;
        result.meanVelocity[d] = globalSum[4u + dim + d] / result.innerArea;
      }

      if(symmetric) {
        for(unsigned d = 0; d < dim - 1u; ++d) {
          result.barycenter[d] = symmetryCoordinates[d];
          result.meanVelocity[d] = 0.;
        }
      }

      if(axisymmetric) {
        result.barycenter[0] = axisCoordinate;
        result.meanVelocity[0] = 0.;
      }

      if(result.interfaceLength > 1.e-14) {
        if(axisymmetric || dim == 3) {
          result.circularity = std::cbrt(36. * pi * result.innerArea * result.innerArea) / result.interfaceLength;
        }
        else {
          result.circularity = 2. * std::sqrt(pi * result.innerArea) / result.interfaceLength;
        }
      }
    }
  }

  result.hasRayleighTaylorMetrics = computeRayleighTaylorMetrics;

  if(computeRayleighTaylorMetrics) {
    if(dim == 2) {
      result.ps = FindLevelSetZeroOnNodalLine(msh, sol, psiIndex, psiType, {0.0});
      result.pb = FindLevelSetZeroOnNodalLine(msh, sol, psiIndex, psiType, {0.5});
    }
    else {
      result.ps = FindLevelSetZeroOnNodalLine(msh, sol, psiIndex, psiType, {0.0, 0.0});
      result.pb = FindLevelSetZeroOnNodalLine(msh, sol, psiIndex, psiType, {0.5, 0.5});
      const double saddle1 = FindLevelSetZeroOnNodalLine(msh, sol, psiIndex, psiType, {0.5, 0.0});
      const double saddle2 = FindLevelSetZeroOnNodalLine(msh, sol, psiIndex, psiType, {0.0, 0.5});
      result.psad = 0.5 * (saddle1 + saddle2);
    }
  }

  result.hasPhaseInversionMetrics = computePhaseInversionMetrics;

  if(computePhaseInversionMetrics) {
    result.fluid1UpperVolume = globalSum[basePhaseInversion + 0u];
    result.potentialEnergy1 = globalSum[basePhaseInversion + 1u];
    result.potentialEnergy2 = globalSum[basePhaseInversion + 2u];
    result.kineticEnergy1 = globalSum[basePhaseInversion + 3u];
    result.kineticEnergy2 = globalSum[basePhaseInversion + 4u];
    result.enstrophy1 = globalSum[basePhaseInversion + 5u];
    result.enstrophy2 = globalSum[basePhaseInversion + 6u];

    if(phaseInversionInitialVolume > 1.e-14) {
      result.fluid1UpperRatio = result.fluid1UpperVolume / phaseInversionInitialVolume;
      result.fluid1VolumeRatio = result.innerArea / phaseInversionInitialVolume;
    }

    result.mechanicalEnergy = result.potentialEnergy1 + result.potentialEnergy2 + result.kineticEnergy1 + result.kineticEnergy2;

    if(potentialEnergy1Final > 1.e-14) {
      result.potentialEnergy1Star = result.potentialEnergy1 / potentialEnergy1Final;
    }

    if(potentialEnergy2Final > 1.e-14) {
      result.potentialEnergy2Star = result.potentialEnergy2 / potentialEnergy2Final;
    }

    if(kineticEnergy1Scale > 1.e-14) {
      result.kineticEnergy1Star = result.kineticEnergy1 / kineticEnergy1Scale;
    }

    if(kineticEnergy2Scale > 1.e-14) {
      result.kineticEnergy2Star = result.kineticEnergy2 / kineticEnergy2Scale;
    }

    if(phaseInversionInitialInterfaceArea > 1.e-14) {
      result.interfaceAreaStar = result.interfaceLength / phaseInversionInitialInterfaceArea;
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

    std::cout << "Time                = " << time << std::endl;
    std::cout << "Inner area          = " << diagnostics.innerArea << std::endl;
    std::cout << "Outer area          = " << diagnostics.outerArea << std::endl;
    std::cout << "Total area          = " << diagnostics.totalArea << std::endl;
    std::cout << "Interface length    = " << diagnostics.interfaceLength << std::endl;

    // ==========================================================
    // Rising bubble
    // ==========================================================

    if(diagnostics.hasRisingBubbleMetrics) {

      for(unsigned d = 0; d < diagnostics.barycenter.size(); ++d) {
        std::cout << "Barycenter[" << d << "]      = " << diagnostics.barycenter[d] << std::endl;
      }

      for(unsigned d = 0; d < diagnostics.meanVelocity.size(); ++d) {
        std::cout << "Mean velocity[" << d << "]   = " << diagnostics.meanVelocity[d] << std::endl;
      }

      std::cout << "Circularity         = " << diagnostics.circularity << std::endl;
    }

    // ==========================================================
    // Rayleigh-Taylor
    // ==========================================================

    if(diagnostics.hasRayleighTaylorMetrics) {

      std::cout << "Bubble position     = " << diagnostics.pb << std::endl;
      std::cout << "Spike position      = " << diagnostics.ps << std::endl;
      std::cout << "Saddle position     = " << diagnostics.psad << std::endl;
    }

    // ==========================================================
    // Phase inversion
    // ==========================================================

    if(diagnostics.hasPhaseInversionMetrics) {

      std::cout << "Fluid 1 upper volume = " << diagnostics.fluid1UpperVolume << std::endl;
      std::cout << "Fluid 1 upper ratio  = " << diagnostics.fluid1UpperRatio << std::endl;
      std::cout << "Fluid 1 volume ratio = " << diagnostics.fluid1VolumeRatio << std::endl;

      std::cout << "Potential energy 1   = " << diagnostics.potentialEnergy1 << std::endl;
      std::cout << "Potential energy 2   = " << diagnostics.potentialEnergy2 << std::endl;

      std::cout << "Kinetic energy 1     = " << diagnostics.kineticEnergy1 << std::endl;
      std::cout << "Kinetic energy 2     = " << diagnostics.kineticEnergy2 << std::endl;

      std::cout << "Enstrophy 1          = " << diagnostics.enstrophy1 << std::endl;
      std::cout << "Enstrophy 2          = " << diagnostics.enstrophy2 << std::endl;

      std::cout << "Mechanical energy    = " << diagnostics.mechanicalEnergy << std::endl;

      std::cout << "Potential energy 1*  = " << diagnostics.potentialEnergy1Star << std::endl;
      std::cout << "Potential energy 2*  = " << diagnostics.potentialEnergy2Star << std::endl;

      std::cout << "Kinetic energy 1*    = " << diagnostics.kineticEnergy1Star << std::endl;
      std::cout << "Kinetic energy 2*    = " << diagnostics.kineticEnergy2Star << std::endl;

      std::cout << "Interface area*      = " << diagnostics.interfaceAreaStar << std::endl;
    }

    return;
  }

  // ============================================================
  // File output
  // ============================================================

  std::ofstream out(fileName, std::ios::app);

  if(!out) {
    throw std::runtime_error("PrintLevelSetDiagnostics: cannot open file " + fileName);
  }

  const unsigned w = 24;

  out << std::scientific;
  out << std::setprecision(16);

  // ============================================================
  // Common quantities
  // ============================================================

  out << std::setw(w) << time;
  out << std::setw(w) << diagnostics.innerArea;
  out << std::setw(w) << diagnostics.outerArea;
  out << std::setw(w) << diagnostics.totalArea;
  out << std::setw(w) << diagnostics.interfaceLength;

  // ============================================================
  // Rising bubble
  // ============================================================

  if(diagnostics.hasRisingBubbleMetrics) {

    for(unsigned d = 0; d < diagnostics.barycenter.size(); ++d) {
      out << std::setw(w) << diagnostics.barycenter[d];
    }

    for(unsigned d = 0; d < diagnostics.meanVelocity.size(); ++d) {
      out << std::setw(w) << diagnostics.meanVelocity[d];
    }

    out << std::setw(w) << diagnostics.circularity;
  }

  // ============================================================
  // Rayleigh-Taylor
  // ============================================================

  if(diagnostics.hasRayleighTaylorMetrics) {

    out << std::setw(w) << diagnostics.pb;
    out << std::setw(w) << diagnostics.ps;
    out << std::setw(w) << diagnostics.psad;
  }

  // ============================================================
  // Phase inversion
  // ============================================================

  if(diagnostics.hasPhaseInversionMetrics) {

    out << std::setw(w) << diagnostics.fluid1UpperVolume;
    out << std::setw(w) << diagnostics.fluid1UpperRatio;
    out << std::setw(w) << diagnostics.fluid1VolumeRatio;

    out << std::setw(w) << diagnostics.potentialEnergy1;
    out << std::setw(w) << diagnostics.potentialEnergy2;

    out << std::setw(w) << diagnostics.kineticEnergy1;
    out << std::setw(w) << diagnostics.kineticEnergy2;

    out << std::setw(w) << diagnostics.enstrophy1;
    out << std::setw(w) << diagnostics.enstrophy2;

    out << std::setw(w) << diagnostics.mechanicalEnergy;

    out << std::setw(w) << diagnostics.potentialEnergy1Star;
    out << std::setw(w) << diagnostics.potentialEnergy2Star;

    out << std::setw(w) << diagnostics.kineticEnergy1Star;
    out << std::setw(w) << diagnostics.kineticEnergy2Star;

    out << std::setw(w) << diagnostics.interfaceAreaStar;
  }

  out << '\n';
}
struct SimulationArgs {
  unsigned uniformLevels = 1u;
  unsigned levelOffset   = 1u;

  unsigned dim = 2u;

  SimulationCase simulation = SimulationCase::RB2;

  std::vector<unsigned> adaptiveLevels = {2};
  std::vector<unsigned> nSteps = {1200};
  double period = 1.288;

  bool symmetric = false;
  bool axisymmetric = false;
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

      args.adaptiveLevels.clear();

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

      args.nSteps.clear();

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

    else if(option == "--symmetric") {

      if(i + 1 >= argc) {
        throw std::runtime_error(
          "Missing value after --symmetric"
        );
      }

      args.symmetric = std::stoul(argv[++i]);
    }

    else if(option == "--axisymmetric") {

      if(i + 1 >= argc) {
        throw std::runtime_error(
          "Missing value after --axisymmetric"
        );
      }

      args.axisymmetric = std::stoul(argv[++i]);
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

  // if (args.uniformLevels == 0u) {
  //   throw std::runtime_error(
  //     "--uniform-levels must be specified and > 0"
  //   );
  // }
  //
  // if (args.adaptiveLevels.empty()) {
  //   throw std::runtime_error(
  //     "--adaptive-levels must contain at least one value"
  //   );
  // }

  // if (args.nSteps.empty()) {
  //   throw std::runtime_error(
  //     "--nsteps must contain at least one value"
  //   );
  // }

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

SimulationConfig GetSimulationConfig(const SimulationCase simulation, const unsigned dim, const bool symmetric, const bool axisymmetric) {

  if(dim != 2 && dim != 3) {
    throw std::runtime_error("GetSimulationConfig: dim must be either 2 or 3");
  }

  if(symmetric && axisymmetric) {
    throw std::runtime_error("GetSimulationConfig: symmetric and axisymmetric cannot both be enabled");
  }

  if(axisymmetric && dim != 2) {
    throw std::runtime_error("GetSimulationConfig: axisymmetric formulation is available only in 2D");
  }

  SimulationConfig cfg;

  cfg.symmetric = symmetric;
  cfg.axisymmetric = axisymmetric;

  cfg.x_min.resize(3);
  cfg.x_max.resize(3);
  cfg.n.resize(3);
  cfg.boundary.Resize(dim);

  const double pi = std::acos(-1.0);

  auto SetAllNeumann = [&cfg, dim]() {
    for(unsigned d = 0; d < dim; ++d) {
      cfg.boundary.min[d] = {BoundaryType::NEUMANN_S, 0.};
      cfg.boundary.max[d] = {BoundaryType::NEUMANN_S, 0.};
    }
  };

  auto SetAllSlip = [&cfg, dim]() {
    for(unsigned d = 0; d < dim; ++d) {
      cfg.boundary.min[d] = {BoundaryType::SLIP, 0.};
      cfg.boundary.max[d] = {BoundaryType::SLIP, 0.};
    }
  };

  auto SetSlipSidesWallEnds = [&cfg, dim]() {
    for(unsigned d = 0; d < dim; ++d) {
      cfg.boundary.min[d] = {BoundaryType::SLIP, 0.};
      cfg.boundary.max[d] = {BoundaryType::SLIP, 0.};
    }

    const unsigned verticalDir = (dim == 2) ? 1 : 2;

    cfg.boundary.min[verticalDir] = {BoundaryType::WALL, 0.};
    cfg.boundary.max[verticalDir] = {BoundaryType::WALL, 0.};
  };

  auto HalveDomainDirection = [&cfg](const unsigned d) {
    if(cfg.n[d] == 0) {
      throw std::runtime_error("GetSimulationConfig: cannot halve a direction with zero elements");
    }

    if(cfg.n[d] % 2 != 0) {
      throw std::runtime_error("GetSimulationConfig: number of elements must be even in a symmetry direction");
    }

    cfg.x_min[d] = 0.5 * (cfg.x_min[d] + cfg.x_max[d]);
    cfg.n[d] /= 2;
    cfg.boundary.min[d] = {BoundaryType::SLIP, 0.};
  };

  switch(simulation) {

    case SimulationCase::SAYE1:
    case SimulationCase::SAYE2: {

      cfg.x_min[0] = -0.5;
      cfg.x_max[0] = 0.5;
      cfg.x_min[1] = -0.5;
      cfg.x_max[1] = 0.5;
      cfg.x_min[2] = (dim == 3) ? -0.5 : 0.;
      cfg.x_max[2] = (dim == 3) ? 0.5 : 0.;

      cfg.n[0] = 8;
      cfg.n[1] = 8;
      cfg.n[2] = (dim == 3) ? 8 : 0;

      cfg.period = 0.01;

      cfg.psi.type = PsiType::SAYE;
      cfg.psi.amplitude = 0.025;
      cfg.psi.waveNumber = 2.0 * pi;

      cfg.properties.rho1 = 1.;
      cfg.properties.rho2 = 1.;
      cfg.properties.mu1 = 1.;
      cfg.properties.mu2 = 1.;
      cfg.properties.sigma = 10.;
      cfg.properties.gravity = 0.;

      SetSlipSidesWallEnds();

      break;
    }

    case SimulationCase::RB1:
    case SimulationCase::RB2: {

      const bool rb1 = simulation == SimulationCase::RB1;

      cfg.x_min[0] = 0.;
      cfg.x_max[0] = 1.;
      cfg.x_min[1] = 0.;
      cfg.x_max[1] = (dim == 3) ? 1. : 2.;
      cfg.x_min[2] = 0.;
      cfg.x_max[2] = (dim == 3) ? 2. : 0.;

      cfg.n[0] = (dim == 3) ? 8 : 8;
      cfg.n[1] = (dim == 3) ? 8 : 16;
      cfg.n[2] = (dim == 3) ? 16 : 0;

      cfg.period = 3.;

      cfg.psi.type = PsiType::BALL;
      cfg.psi.center = {0.5, 0.5, 0.5};
      cfg.psi.center.resize(dim);
      cfg.psi.radius = 0.25;

      cfg.properties.rho1 = rb1 ? 100. : 1.;
      cfg.properties.rho2 = 1000.;
      cfg.properties.mu1 = rb1 ? 1. : 0.1;
      cfg.properties.mu2 = 10.;
      cfg.properties.sigma = rb1 ? 24.5 : 1.96;
      cfg.properties.gravity = -0.98;

      SetSlipSidesWallEnds();

      break;
    }

    case SimulationCase::BB: {
      if(dim != 2) throw std::runtime_error("GetSimulationConfig: BB requires dim = 2");
      const double La = 6.7e4;
      const double Bo = 0.1172;
      const double Oh = 1. / std::sqrt(La);

      cfg.x_min[0] = -10.;
      cfg.x_max[0] = 10.;
      cfg.x_min[1] = -5.;
      cfg.x_max[1] = 5.;
      cfg.x_min[2] = 0.;
      cfg.x_max[2] = 0.;

      cfg.n[0] = 16;
      cfg.n[1] = 8;
      cfg.n[2] = 0;

      cfg.period = 1.2;

      cfg.psi.type = PsiType::BURSTING_BUBBLE;
      cfg.psi.bondNumber = Bo;
      cfg.psi.burstingPointLess = 400u;
      cfg.psi.burstingOuterRadius = 10.;
      cfg.psi.freeSurfaceHeight = 0.;

      cfg.properties.rho1 = 1. / 998.;
      cfg.properties.rho2 = 1.;
      cfg.properties.mu1 = Oh / 55.;
      cfg.properties.mu2 = Oh;
      cfg.properties.sigma = 1.;
      cfg.properties.gravity = -Bo;
      SetAllSlip();
      break;
    }

    case SimulationCase::MB: {

      cfg.x_min[0] = 0;
      cfg.x_max[0] = 1.;
      cfg.x_min[1] = 0.;
      cfg.x_max[1] = (dim == 2) ? 1.5 : 1.;
      cfg.x_min[2] = 0.;
      cfg.x_max[2] = (dim == 2) ? 0. : 1.;

      cfg.n[0] = 8;
      cfg.n[1] = 8;
      cfg.n[2] = (dim == 2) ? 0 : 8;

      cfg.period = (dim == 2) ? 0.35 : 0.25;

      cfg.psi.type = PsiType::TWO_BUBBLES;
      if (dim == 2) {
        cfg.psi.bubbleCenters = {{0.5, 0.65}, {0.5, 0.35}};
      }
      else {
        cfg.psi.bubbleCenters = {{0.5, 0.5, 0.58}, {0.5, 0.5, 0.3}};
      }
      cfg.psi.bubbleRadii = {0.15, 0.1};

      cfg.properties.rho1 = (dim == 2) ? 100. : 1.;
      cfg.properties.rho2 = (dim == 2) ? 1000. : 10.;
      cfg.properties.mu1 = (dim == 2) ? 0.078 : 2.5e-4;
      cfg.properties.mu2 = (dim == 2) ? 0.156 : 5.e-4;
      cfg.properties.sigma = (dim == 2) ? 1.226 : 0;
      cfg.properties.gravity = -9.81;

      SetSlipSidesWallEnds();

      break;
    }

    case SimulationCase::RT: {

      cfg.x_min[0] = -0.5;
      cfg.x_max[0] = 0.5;
      cfg.x_min[1] = (dim == 2) ? -2.0 : -0.5;
      cfg.x_max[1] = (dim == 2) ? 2.0 : 0.5;
      cfg.x_min[2] = (dim == 2) ? 0. : -2.;
      cfg.x_max[2] = (dim == 2) ? 0. : 2.;

      cfg.n[0] = 8;
      cfg.n[1] = (dim == 2) ? 32 : 8;
      cfg.n[2] = (dim == 2) ? 0 : 32;

      cfg.period = (dim == 2) ? 3. : 4.5;

      cfg.psi.type = PsiType::RAYLEIGH_TAYLOR;
      cfg.psi.amplitude = (dim == 2) ? 0.005 : 0.05;
      cfg.psi.waveNumber = 2.0 * pi;

      cfg.properties.rho1 = (dim == 2) ? 1.5 : 1.;
      cfg.properties.rho2 = (dim == 2) ? 1. : 1. / 3.;
      cfg.properties.mu1 = (dim == 2) ? 0.009535 : 1. / 1024.;
      cfg.properties.mu2 = (dim == 2) ? 0.006349 : 1. / 1024.;
      cfg.properties.sigma = (dim == 2) ? 0 : 1.e-5;
      cfg.properties.gravity = (dim == 2) ? -10. : -1.;

      (dim == 2) ? SetSlipSidesWallEnds() : SetAllSlip();

      break;
    }

    case SimulationCase::DB: {

      if(dim != 2) {
        throw std::runtime_error("GetSimulationConfig: DB is currently implemented only in 2D");
      }

      const double a = 0.05715;

      cfg.x_min[0] = 0.;
      cfg.x_max[0] = 5.0 * a;
      cfg.x_min[1] = 0.;
      cfg.x_max[1] = 3.0 * a;
      cfg.x_min[2] = 0.;
      cfg.x_max[2] = 0.;

      cfg.n[0] = 10;
      cfg.n[1] = 6;
      cfg.n[2] = 0;

      cfg.period = 1.;

      cfg.psi.type = PsiType::DAM_BREAK;
      cfg.psi.damWidth = a;
      cfg.psi.damHeight = 2.0 * a;

      cfg.properties.rho1 = 1000.;
      cfg.properties.rho2 = 1.2;
      cfg.properties.mu1 = 1.e-3;
      cfg.properties.mu2 = 1.8e-5;
      cfg.properties.sigma = 0.;
      cfg.properties.gravity = -9.81;

      SetAllSlip();

      break;
    }

    case SimulationCase::PI: {

      if(dim != 3) {
        throw std::runtime_error("GetSimulationConfig: PI requires dim = 3");
      }

      if(symmetric || axisymmetric) {
        throw std::runtime_error("GetSimulationConfig: PI does not support symmetric or axisymmetric reduction");
      }

      const double H = 0.1;

      cfg.x_min[0] = 0.;
      cfg.x_max[0] = H;
      cfg.x_min[1] = 0.;
      cfg.x_max[1] = H;
      cfg.x_min[2] = 0.;
      cfg.x_max[2] = H;

      cfg.n[0] = 8;
      cfg.n[1] = 8;
      cfg.n[2] = 8;

      cfg.period = 12.86;

      cfg.psi.type = PsiType::PHASE_INVERSION;
      cfg.psi.phaseInversionSize = 0.5 * H;

      cfg.properties.rho1 = 900.;
      cfg.properties.rho2 = 1000.;
      cfg.properties.mu1 = 0.1;
      cfg.properties.mu2 = 0.001;
      cfg.properties.sigma = 0.045;
      cfg.properties.gravity = -9.81;

      SetAllSlip();

      break;
    }

    default: {
      throw std::runtime_error("GetSimulationConfig: unknown SimulationCase");
    }
  }

  if(axisymmetric && dim == 2) {
    HalveDomainDirection(0);
  }
  else if(symmetric) {
    HalveDomainDirection(0);

    if(dim == 3) {
      HalveDomainDirection(1);
    }
  }

  return cfg;
}

#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

struct BBDropletMetrics {
  unsigned id = 0u;
  double volume = 0.;
  double radius = 0.;
  double radiusOverRb = 0.;
  double zCm = 0.;
  double zCmOverRb = 0.;
  double verticalVelocity = 0.;
  double capillaryNumber = 0.;
};

struct BurstingBubbleMetrics {
  double time = 0.;
  double timeOverTc = 0.;

  bool hasT0 = false;
  double t0 = std::numeric_limits<double>::quiet_NaN();
  double shiftedTimeOverTc = std::numeric_limits<double>::quiet_NaN();

  bool hasJetTip = false;
  double jetTipZ = std::numeric_limits<double>::quiet_NaN();
  double jetTipZOverRb = std::numeric_limits<double>::quiet_NaN();
  double jetTipVelocity = std::numeric_limits<double>::quiet_NaN();
  double jetTipVelocityFD = std::numeric_limits<double>::quiet_NaN();
  double jetTipCa = std::numeric_limits<double>::quiet_NaN();

  double bulkLiquidVolume = 0.;
  unsigned nLiquidComponents = 0u;
  unsigned nDroplets = 0u;
  std::vector<BBDropletMetrics> droplets;

  bool firstDropletDetected = false;
  double firstDropletTime = std::numeric_limits<double>::quiet_NaN();
  double firstDropletTimeOverTc = std::numeric_limits<double>::quiet_NaN();
  double firstDropletRadiusOverRb = std::numeric_limits<double>::quiet_NaN();
  double firstDropletVelocity = std::numeric_limits<double>::quiet_NaN();
  double firstDropletCa = std::numeric_limits<double>::quiet_NaN();
};

struct BurstingBubbleOptions {
  // In BB2 psi > 0 is the gas phase and psi < 0 is the liquid phase.
  bool liquidPositive = false;

  // A detached component is considered a physical droplet only if
  // R_eq >= minRadiusOverH * hMin.
  double minRadiusOverH = 2.;
  double minEquivalentRadius = 0.;
  double minDropletVolume = 0.;

  // Connectivity threshold. An element belongs to the liquid graph only if
  // V_liquid / V_element is larger than this value. This deliberately avoids
  // connecting components only through a numerically non-negligible liquid amount.
  double graphLiquidFractionTolerance = 1.e-12;

  unsigned persistenceSteps = 2u;
  double freeSurfaceZ = 0.;
  double axisTolerance = -1.;
  double psiZeroTolerance = 1.e-12;
  double dropletTrackingTolerance = -1.;
};

struct BurstingBubbleState {
  struct TrackedDroplet {
    unsigned id = 0u;
    double zCm = 0.;
    double verticalVelocity = 0.;
  };

  bool hasPreviousCall = false;
  double previousCallTime = 0.;

  bool hasPreviousTip = false;
  double previousTipTime = 0.;
  double previousJetTipZ = 0.;

  bool hasT0 = false;
  double t0 = std::numeric_limits<double>::quiet_NaN();

  bool firstDropletDetected = false;
  unsigned pendingFirstDropSteps = 0u;
  double pendingFirstDropTime = std::numeric_limits<double>::quiet_NaN();
  double pendingFirstDropRadius = std::numeric_limits<double>::quiet_NaN();
  double pendingFirstDropVelocity = std::numeric_limits<double>::quiet_NaN();
  double pendingFirstDropCa = std::numeric_limits<double>::quiet_NaN();

  double firstDropletTime = std::numeric_limits<double>::quiet_NaN();
  double firstDropletRadius = std::numeric_limits<double>::quiet_NaN();
  double firstDropletVelocity = std::numeric_limits<double>::quiet_NaN();
  double firstDropletCa = std::numeric_limits<double>::quiet_NaN();

  unsigned nextDropletId = 1u;
  std::vector<TrackedDroplet> previousDroplets;

  void ResetPendingFirstDrop() {
    pendingFirstDropSteps = 0u;
    pendingFirstDropTime = std::numeric_limits<double>::quiet_NaN();
    pendingFirstDropRadius = std::numeric_limits<double>::quiet_NaN();
    pendingFirstDropVelocity = std::numeric_limits<double>::quiet_NaN();
    pendingFirstDropCa = std::numeric_limits<double>::quiet_NaN();
  }
};

namespace BurstingBubbleDetail {

  struct UnionFind {
    explicit UnionFind(const unsigned n) : parent(n), rank(n, 0u) {
      std::iota(parent.begin(), parent.end(), 0u);
    }

    unsigned Find(unsigned i) {
      if(parent[i] != i) parent[i] = Find(parent[i]);
      return parent[i];
    }

    void Unite(unsigned a, unsigned b) {
      a = Find(a);
      b = Find(b);

      if(a == b) return;

      if(rank[a] < rank[b]) std::swap(a, b);

      parent[b] = a;

      if(rank[a] == rank[b]) ++rank[a];
    }

    std::vector<unsigned> parent;
    std::vector<unsigned> rank;
  };

  template <class T>
  std::vector<T> AllGatherVector(const std::vector<T>& local, MPI_Datatype mpiType) {
    int nprocs = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    const int nLocal = static_cast<int>(local.size());

    std::vector<int> counts(nprocs, 0);
    std::vector<int> displs(nprocs, 0);

    MPI_Allgather(
      &nLocal,
      1,
      MPI_INT,
      counts.data(),
      1,
      MPI_INT,
      MPI_COMM_WORLD);

    int nGlobal = 0;

    for(int p = 0; p < nprocs; ++p) {
      displs[p] = nGlobal;
      nGlobal += counts[p];
    }

    std::vector<T> global(nGlobal);

    MPI_Allgatherv(
      local.empty() ? nullptr : local.data(),
      nLocal,
      mpiType,
      global.empty() ? nullptr : global.data(),
      counts.data(),
      displs.data(),
      mpiType,
      MPI_COMM_WORLD);

    return global;
  }

  inline bool FileEmptyOrMissing(const std::string& fileName) {
    std::ifstream in(fileName, std::ios::binary | std::ios::ate);
    return !in || in.tellg() == std::streampos(0);
  }

  struct AxisRoot {
    bool found = false;
    double z = std::numeric_limits<double>::quiet_NaN();
    double uz = std::numeric_limits<double>::quiet_NaN();
  };

  AxisRoot FindLowestAxisRoot(
    Mesh* msh,
    Solution* sol,
    const unsigned psiIndex,
    const unsigned psiType,
    const unsigned uzIndex,
    const unsigned velocityType,
    const double axisCoordinate,
    const double coordTol,
    const double psiTol) {

    const unsigned dim = msh->GetDimension();
    const unsigned iproc = msh->processor_id();
    const unsigned xType = 2u;

    if(dim != 2u) {
      throw std::runtime_error("FindLowestAxisRoot: requires a 2D mesh");
    }

    struct AxisNode {
      unsigned psiDof = 0u;
      double z = 0.;
      double psi = 0.;
      double uz = 0.;
    };

    std::vector<unsigned> localDof;
    std::vector<double> localZ;
    std::vector<double> localPsi;
    std::vector<double> localUz;

    for(unsigned iel = msh->_elementOffset[iproc];
        iel < msh->_elementOffset[iproc + 1u];
        ++iel) {

      const unsigned nPsi = msh->GetElementDofNumber(iel, psiType);
      const unsigned nVel = msh->GetElementDofNumber(iel, velocityType);
      const unsigned nX = msh->GetElementDofNumber(iel, xType);

      if(nPsi != nX) {
        throw std::runtime_error(
          "FindLowestAxisRoot: psi and coordinates must have the same nodal structure");
      }

      if(nVel != nPsi) {
        throw std::runtime_error(
          "FindLowestAxisRoot: velocity and psi must have the same nodal structure");
      }

      for(unsigned i = 0; i < nPsi; ++i) {
        const unsigned xDof = msh->GetSolutionDof(i, iel, xType);
        const double r = (*msh->_topology->_Sol[0])(xDof);

        if(std::abs(r - axisCoordinate) > coordTol) continue;

        const unsigned psiDof = msh->GetSolutionDof(i, iel, psiType);
        const unsigned velDof = msh->GetSolutionDof(i, iel, velocityType);

        localDof.push_back(psiDof);
        localZ.push_back((*msh->_topology->_Sol[1])(xDof));
        localPsi.push_back((*sol->_Sol[psiIndex])(psiDof));
        localUz.push_back((*sol->_Sol[uzIndex])(velDof));
      }
    }

    const std::vector<unsigned> globalDof =
      AllGatherVector(localDof, MPI_UNSIGNED);

    const std::vector<double> globalZ =
      AllGatherVector(localZ, MPI_DOUBLE);

    const std::vector<double> globalPsi =
      AllGatherVector(localPsi, MPI_DOUBLE);

    const std::vector<double> globalUz =
      AllGatherVector(localUz, MPI_DOUBLE);

    if(globalDof.size() != globalZ.size() ||
        globalDof.size() != globalPsi.size() ||
        globalDof.size() != globalUz.size()) {
      throw std::runtime_error("FindLowestAxisRoot: inconsistent gathered nodal data");
    }

    std::vector<AxisNode> nodes(globalDof.size());

    for(unsigned i = 0; i < nodes.size(); ++i) {
      nodes[i].psiDof = globalDof[i];
      nodes[i].z = globalZ[i];
      nodes[i].psi = globalPsi[i];
      nodes[i].uz = globalUz[i];
    }

    std::sort(
      nodes.begin(),
      nodes.end(),
    [](const AxisNode & a, const AxisNode & b) {
      return a.psiDof < b.psiDof;
    });

    nodes.erase(
      std::unique(
        nodes.begin(),
        nodes.end(),
    [](const AxisNode & a, const AxisNode & b) {
      return a.psiDof == b.psiDof;
    }),
    nodes.end());

    std::sort(
      nodes.begin(),
      nodes.end(),
    [](const AxisNode & a, const AxisNode & b) {
      return a.z < b.z;
    });

    AxisRoot result;

    auto UpdateRoot = [&](const double z, const double uz) {
      if(!result.found || z < result.z) {
        result.found = true;
        result.z = z;
        result.uz = uz;
      }
    };

    for(const AxisNode& node : nodes) {
      if(std::abs(node.psi) <= psiTol) {
        UpdateRoot(node.z, node.uz);
      }
    }

    for(unsigned i = 0; i + 1u < nodes.size(); ++i) {
      const double psi0 = nodes[i].psi;
      const double psi1 = nodes[i + 1u].psi;

      if(psi0 * psi1 < 0.) {
        const double alpha = psi0 / (psi0 - psi1);
        const double z = nodes[i].z + alpha * (nodes[i + 1u].z - nodes[i].z);
        const double uz = nodes[i].uz + alpha * (nodes[i + 1u].uz - nodes[i].uz);
        UpdateRoot(z, uz);
      }
    }

    return result;
  }

} // namespace BurstingBubbleDetail

BurstingBubbleMetrics ComputeBurstingBubbleMetrics(
  MultiLevelSolution& mlSol,
  const std::string& psiName,
  const std::vector<std::string>& velocityName,
  const unsigned velocityLevel,
  const double time,
  const double hMin,
  const double bubbleRadius,
  const double rhoLiquid,
  const double muLiquid,
  const double sigma,
  const double axisCoordinate,
  BurstingBubbleState& state,
  const BurstingBubbleOptions& options = BurstingBubbleOptions()) {

  using namespace BurstingBubbleDetail;

  MultiLevelMesh& mlMsh = *mlSol.GetMultilevelMesh();
  const unsigned level = mlMsh.GetNumberOfLevels() - 1u;
  Mesh* msh = mlMsh.GetLevel(level);
  Solution* sol = mlSol.GetSolutionLevel(level);

  const unsigned dim = msh->GetDimension();
  const unsigned iproc = msh->processor_id();
  const unsigned psiIndex = mlSol.GetIndex(psiName.c_str());
  const unsigned psiType = mlSol.GetSolutionType(psiIndex);
  const unsigned xType = 2u;
  const double pi = std::acos(-1.0);

  if(dim != 2u) {
    throw std::runtime_error(
      "ComputeBurstingBubbleMetrics: requires a 2D axisymmetric mesh");
  }

  if(velocityName.size() != 2u) {
    throw std::runtime_error(
      "ComputeBurstingBubbleMetrics: velocityName must contain two components");
  }

  if(hMin <= 0. ||
      bubbleRadius <= 0. ||
      rhoLiquid <= 0. ||
      muLiquid <= 0. ||
      sigma <= 0.) {

    std::ostringstream oss;
    oss << "ComputeBurstingBubbleMetrics: invalid parameters:"
        << " hMin=" << hMin
        << " Rb=" << bubbleRadius
        << " rhoL=" << rhoLiquid
        << " muL=" << muLiquid
        << " sigma=" << sigma;

    throw std::runtime_error(oss.str());
  }

  if(options.graphLiquidFractionTolerance <= 0. ||
      options.graphLiquidFractionTolerance >= 1.) {
    throw std::runtime_error(
      "ComputeBurstingBubbleMetrics: graphLiquidFractionTolerance must be in (0,1)");
  }

  if(velocityLevel > level) {
    throw std::runtime_error(
      "ComputeBurstingBubbleMetrics: velocityLevel > diagnostics level");
  }

  const unsigned urIndex = mlSol.GetIndex(velocityName[0].c_str());
  const unsigned uzIndex = mlSol.GetIndex(velocityName[1].c_str());
  const unsigned velocityType = mlSol.GetSolutionType(urIndex);

  if(mlSol.GetSolutionType(uzIndex) != velocityType) {
    throw std::runtime_error(
      "ComputeBurstingBubbleMetrics: velocity FE types differ");
  }

  for(unsigned l = velocityLevel; l < level; ++l) {
    Solution* solL = mlSol.GetSolutionLevel(l);
    Solution* solLp1 = mlSol.GetSolutionLevel(l + 1u);
    Mesh* mshLp1 = mlMsh.GetLevel(l + 1u);
    SparseMatrix* Pv = mshLp1->GetCoarseToFineProjection(velocityType);

    solLp1->_Sol[urIndex]->matrix_mult(*(solL->_Sol[urIndex]), *Pv);
    solLp1->_Sol[uzIndex]->matrix_mult(*(solL->_Sol[uzIndex]), *Pv);
  }

  const double tc =
    std::sqrt(
      rhoLiquid *
      bubbleRadius * bubbleRadius * bubbleRadius /
      sigma);

  const double axisTol =
    options.axisTolerance > 0.
    ? options.axisTolerance
    : 1.e-10 * std::max(1., bubbleRadius);

  const double minRadius =
    std::max(
      options.minEquivalentRadius,
      options.minRadiusOverH * hMin);

  const double minVolume =
    std::max(
      options.minDropletVolume,
      4. * pi * minRadius * minRadius * minRadius / 3.);

  auto LiquidPhi = [&](const double value) {
    return options.liquidPositive ? value : -value;
  };

  // ------------------------------------------------------------
  // Build one record for every element containing a non-negligible
  // amount of liquid. Connectivity uses the SAME mesh-neighbour API
  // already used by FlagFinestMeshLevel():
  //
  //   msh->el->GetElementNearElementSize(iel, 1)
  //   msh->el->GetElementNearElement(iel, j)
  //
  // This is robust across MPI partitions and AMR transitions.
  // ------------------------------------------------------------

  std::vector<unsigned> localIds;
  std::vector<double> localData;
  std::vector<unsigned> localNeighbourPairs;

  std::vector<std::vector<double>> coordX(dim);
  std::vector<std::vector<double>> velocity(dim);
  std::vector<double> psi;
  std::vector<double> phiPsi;
  std::vector<double> phiPsi_x;
  std::vector<double> phiX;
  std::vector<double> phiX_x;
  std::vector<double> phiVelocity;
  std::vector<double> phiVelocity_x;

  for(unsigned iel = msh->_elementOffset[iproc];
      iel < msh->_elementOffset[iproc + 1u];
      ++iel) {

    const unsigned geom = msh->GetElementType(iel);
    const unsigned nPsi = msh->GetElementDofNumber(iel, psiType);
    const unsigned nX = msh->GetElementDofNumber(iel, xType);
    const unsigned nV = msh->GetElementDofNumber(iel, velocityType);

    psi.resize(nPsi);
    coordX.assign(dim, std::vector<double>(nX, 0.));
    velocity.assign(dim, std::vector<double>(nV, 0.));

    for(unsigned i = 0; i < nPsi; ++i) {
      const unsigned dof = msh->GetSolutionDof(i, iel, psiType);
      psi[i] = (*sol->_Sol[psiIndex])(dof);
    }

    for(unsigned i = 0; i < nV; ++i) {
      const unsigned dof = msh->GetSolutionDof(i, iel, velocityType);
      velocity[0][i] = (*sol->_Sol[urIndex])(dof);
      velocity[1][i] = (*sol->_Sol[uzIndex])(dof);
    }

    for(unsigned i = 0; i < nX; ++i) {
      const unsigned dof = msh->GetSolutionDof(i, iel, xType);

      for(unsigned d = 0; d < dim; ++d) {
        coordX[d][i] = (*msh->_topology->_Sol[d])(dof);
      }
    }

    const elem_type* femPsi = fem.GetFiniteElement(geom, psiType);
    const elem_type* femX = fem.GetFiniteElement(geom, xType);
    const elem_type* femVelocity = fem.GetFiniteElement(geom, velocityType);

    const unsigned nGauss = femPsi->GetGaussPointNumber();

    if(nGauss != cfw[geom]->GetGaussQuadraturePointNumber()) {
      throw std::runtime_error(
        "ComputeBurstingBubbleMetrics: incompatible quadrature");
    }

    if(femVelocity->GetGaussPointNumber() != nGauss) {
      throw std::runtime_error(
        "ComputeBurstingBubbleMetrics: incompatible velocity quadrature");
    }

    std::vector<double> psiG(nGauss, 0.);

    double psiMin = std::numeric_limits<double>::max();
    double psiMax = -std::numeric_limits<double>::max();

    for(const double value : psi) {
      psiMin = std::min(psiMin, value);
      psiMax = std::max(psiMax, value);
    }

    for(unsigned ig = 0; ig < nGauss; ++ig) {
      double* phi = femPsi->GetPhi(ig);

      for(unsigned i = 0; i < nPsi; ++i) {
        psiG[ig] += psi[i] * phi[i];
      }

      psiMin = std::min(psiMin, psiG[ig]);
      psiMax = std::max(psiMax, psiG[ig]);
    }

    const bool cut =
      psiMin <= 0. &&
      psiMax >= 0. &&
      psiMax - psiMin > 1.e-14;

    std::vector<TypeIO> weightPositive(nGauss, 0.);
    std::vector<TypeIO> weightLiquid(nGauss, 0.);

    if(cut) {
      std::vector<const double*> xg(dim);

      for(unsigned d = 0; d < dim; ++d) {
        xg[d] = femPsi->GetGaussRule().GetGaussCoordinatePointer(d);
      }

      std::vector<double> a;
      BestFitLinearInterpolation(xg, psiG, a);

      double interfaceConstant = a[dim];
      a.resize(dim);

      cfw[geom]->GetWeightWithMap(
        0,
        a,
        interfaceConstant,
        weightPositive);

      for(unsigned ig = 0; ig < nGauss; ++ig) {
        weightLiquid[ig] =
          options.liquidPositive
          ? weightPositive[ig]
          : static_cast<TypeIO>(1. - weightPositive[ig]);
      }
    }
    else {
      const bool liquidElement =
        options.liquidPositive
        ? psiMax > 0.
        : psiMax <= 0.;

      if(liquidElement) {
        std::fill(
          weightLiquid.begin(),
          weightLiquid.end(),
          static_cast<TypeIO>(1.));
      }
    }

    double elementVolume = 0.;
    double liquidVolume = 0.;
    double zMoment = 0.;
    double uzMoment = 0.;

    for(unsigned ig = 0; ig < nGauss; ++ig) {
      double weight = 0.;
      double weightX = 0.;
      double weightV = 0.;

      femPsi->Jacobian(
        coordX,
        ig,
        weight,
        phiPsi,
        phiPsi_x);

      femX->Jacobian(
        coordX,
        ig,
        weightX,
        phiX,
        phiX_x);

      femVelocity->Jacobian(
        coordX,
        ig,
        weightV,
        phiVelocity,
        phiVelocity_x);

      std::vector<double> x(dim, 0.);

      for(unsigned i = 0; i < nX; ++i) {
        for(unsigned d = 0; d < dim; ++d) {
          x[d] += coordX[d][i] * phiX[i];
        }
      }

      const double r = x[0] - axisCoordinate;

      if(r < -1.e-12) {
        throw std::runtime_error(
          "ComputeBurstingBubbleMetrics: negative axisymmetric radius");
      }

      const double dV =
        weight *
        2. * pi * std::max(r, 0.);

      const double dVL =
        dV * static_cast<double>(weightLiquid[ig]);

      elementVolume += dV;

      if(dVL <= 0.) continue;

      double uz = 0.;

      for(unsigned i = 0; i < nV; ++i) {
        uz += velocity[1][i] * phiVelocity[i];
      }

      liquidVolume += dVL;
      zMoment += x[1] * dVL;
      uzMoment += uz * dVL;
    }

    if(elementVolume <= 0.) continue;

    const double liquidFraction =
      liquidVolume / elementVolume;

    if(liquidFraction <= options.graphLiquidFractionTolerance) {
      continue;
    }

    localIds.push_back(iel);
    localData.push_back(liquidVolume);
    localData.push_back(zMoment);
    localData.push_back(uzMoment);
    localData.push_back(liquidFraction);

    const unsigned nNear =
      msh->el->GetElementNearElementSize(iel, 1);

    for(unsigned j = 0; j < nNear; ++j) {
      const unsigned jel =
        msh->el->GetElementNearElement(iel, j);

      if(jel == iel) continue;

      localNeighbourPairs.push_back(iel);
      localNeighbourPairs.push_back(jel);
    }
  }

  const std::vector<unsigned> ids =
    AllGatherVector(localIds, MPI_UNSIGNED);

  const std::vector<double> data =
    AllGatherVector(localData, MPI_DOUBLE);

  const std::vector<unsigned> neighbourPairs =
    AllGatherVector(localNeighbourPairs, MPI_UNSIGNED);

  if(data.size() != 4u * ids.size()) {
    throw std::runtime_error(
      "ComputeBurstingBubbleMetrics: inconsistent gathered element data");
  }

  if(neighbourPairs.size() % 2u != 0u) {
    throw std::runtime_error(
      "ComputeBurstingBubbleMetrics: inconsistent gathered neighbour data");
  }

  BurstingBubbleMetrics result;
  result.time = time;
  result.timeOverTc = time / tc;

  if(ids.empty()) {
    state.previousDroplets.clear();
    state.ResetPendingFirstDrop();
    state.hasPreviousCall = true;
    state.previousCallTime = time;
    return result;
  }

  std::unordered_map<unsigned, unsigned> idToPos;
  idToPos.reserve(2u * ids.size());

  for(unsigned i = 0; i < ids.size(); ++i) {
    idToPos[ids[i]] = i;
  }

  UnionFind uf(static_cast<unsigned>(ids.size()));

  // A pair is connected only if BOTH elements survived the liquid-fraction
  // threshold. Therefore a psi=0 contact or an infinitesimal liquid sliver
  // cannot keep two macroscopic components connected.
  for(unsigned k = 0; k < neighbourPairs.size(); k += 2u) {
    const auto ia = idToPos.find(neighbourPairs[k]);
    const auto ib = idToPos.find(neighbourPairs[k + 1u]);

    if(ia != idToPos.end() && ib != idToPos.end()) {
      uf.Unite(ia->second, ib->second);
    }
  }

  struct Component {
    double volume = 0.;
    double zMoment = 0.;
    double uzMoment = 0.;
  };

  std::unordered_map<unsigned, Component> components;

  for(unsigned i = 0; i < ids.size(); ++i) {
    Component& c = components[uf.Find(i)];
    c.volume += data[4u * i];
    c.zMoment += data[4u * i + 1u];
    c.uzMoment += data[4u * i + 2u];
  }

  result.nLiquidComponents =
    static_cast<unsigned>(components.size());

  unsigned bulkRoot = components.begin()->first;
  double bulkVolume = components.begin()->second.volume;

  for(const auto& item : components) {
    if(item.second.volume > bulkVolume) {
      bulkRoot = item.first;
      bulkVolume = item.second.volume;
    }
  }

  result.bulkLiquidVolume = bulkVolume;

  struct CurrentDrop {
    unsigned id = 0u;
    double volume = 0.;
    double radius = 0.;
    double z = 0.;
    double vz = 0.;
  };

  std::vector<CurrentDrop> current;

  for(const auto& item : components) {
    if(item.first == bulkRoot) continue;
    if(item.second.volume < minVolume) continue;

    CurrentDrop d;
    d.volume = item.second.volume;
    d.radius = std::cbrt(3. * d.volume / (4. * pi));
    d.z = item.second.zMoment / d.volume - options.freeSurfaceZ;
    d.vz = item.second.uzMoment / d.volume;

    current.push_back(d);
  }

  // ------------------------------------------------------------
  // Persistent droplet IDs
  // ------------------------------------------------------------

  const double dtTrack =
    state.hasPreviousCall
    ? std::max(time - state.previousCallTime, 0.)
    : 0.;

  std::vector<bool> used(
    state.previousDroplets.size(),
    false);

  for(CurrentDrop& d : current) {
    double best = std::numeric_limits<double>::max();
    int bestIndex = -1;

    for(unsigned j = 0; j < state.previousDroplets.size(); ++j) {
      if(used[j]) continue;

      const auto& old = state.previousDroplets[j];

      const double predicted =
        old.zCm +
        old.verticalVelocity * dtTrack;

      const double distance =
        std::abs(d.z - predicted);

      if(distance < best) {
        best = distance;
        bestIndex = static_cast<int>(j);
      }
    }

    const double trackingTolerance =
      options.dropletTrackingTolerance > 0.
      ? options.dropletTrackingTolerance
      : std::max(8. * hMin, 0.05 * bubbleRadius);

    if(bestIndex >= 0 && best <= trackingTolerance) {
      d.id =
        state.previousDroplets[
          static_cast<unsigned>(bestIndex)].id;

      used[static_cast<unsigned>(bestIndex)] = true;
    }
    else {
      d.id = state.nextDropletId++;
    }
  }

  std::sort(
    current.begin(),
    current.end(),
  [](const CurrentDrop & a, const CurrentDrop & b) {
    return a.id < b.id;
  });

  state.previousDroplets.clear();

  for(const CurrentDrop& d : current) {
    BBDropletMetrics m;
    m.id = d.id;
    m.volume = d.volume;
    m.radius = d.radius;
    m.radiusOverRb = d.radius / bubbleRadius;
    m.zCm = d.z;
    m.zCmOverRb = d.z / bubbleRadius;
    m.verticalVelocity = d.vz;
    m.capillaryNumber = muLiquid * d.vz / sigma;

    result.droplets.push_back(m);

    BurstingBubbleState::TrackedDroplet tracked;
    tracked.id = d.id;
    tracked.zCm = d.z;
    tracked.verticalVelocity = d.vz;

    state.previousDroplets.push_back(tracked);
  }

  result.nDroplets =
    static_cast<unsigned>(result.droplets.size());

  // ------------------------------------------------------------
  // First-drop event
  // ------------------------------------------------------------

  if(!state.firstDropletDetected) {
    if(!current.empty()) {
      const auto first =
        std::min_element(
          current.begin(),
          current.end(),
      [](const CurrentDrop & a, const CurrentDrop & b) {
        return a.id < b.id;
      });

      if(state.pendingFirstDropSteps == 0u) {
        state.pendingFirstDropTime = time;
        state.pendingFirstDropRadius = first->radius;
        state.pendingFirstDropVelocity = first->vz;
        state.pendingFirstDropCa = muLiquid * first->vz / sigma;
      }

      ++state.pendingFirstDropSteps;

      if(state.pendingFirstDropSteps >=
          std::max(1u, options.persistenceSteps)) {

        state.firstDropletDetected = true;
        state.firstDropletTime = state.pendingFirstDropTime;
        state.firstDropletRadius = state.pendingFirstDropRadius;
        state.firstDropletVelocity = state.pendingFirstDropVelocity;
        state.firstDropletCa = state.pendingFirstDropCa;
      }
    }
    else {
      state.ResetPendingFirstDrop();
    }
  }

  // ------------------------------------------------------------
  // Central cavity apex / jet tip.
  //
  // Do NOT use the connectivity graph here. Along the axis the physical jet
  // tip is the LOWEST psi=0 root. Detached jet drops generate additional roots
  // above it, so choosing the lowest root keeps tracking the bulk jet.
  // ------------------------------------------------------------

  const AxisRoot axisRoot =
    FindLowestAxisRoot(
      msh,
      sol,
      psiIndex,
      psiType,
      uzIndex,
      velocityType,
      axisCoordinate,
      axisTol,
      options.psiZeroTolerance);

  if(axisRoot.found) {
    result.hasJetTip = true;
    result.jetTipZ = axisRoot.z - options.freeSurfaceZ;
    result.jetTipZOverRb = result.jetTipZ / bubbleRadius;
    result.jetTipVelocity = axisRoot.uz;
    result.jetTipCa = muLiquid * result.jetTipVelocity / sigma;

    if(state.hasPreviousTip && time > state.previousTipTime) {
      result.jetTipVelocityFD =
        (result.jetTipZ - state.previousJetTipZ) /
        (time - state.previousTipTime);
    }

    if(!state.hasT0 &&
        state.hasPreviousTip &&
        state.previousJetTipZ <= 0. &&
        result.jetTipZ > 0.) {

      const double dz =
        result.jetTipZ -
        state.previousJetTipZ;

      const double alpha =
        std::abs(dz) > 1.e-14
        ? -state.previousJetTipZ / dz
        : 1.;

      state.t0 =
        state.previousTipTime +
        alpha *
        (time - state.previousTipTime);

      state.hasT0 = true;
    }

    state.previousTipTime = time;
    state.previousJetTipZ = result.jetTipZ;
    state.hasPreviousTip = true;
  }

  state.hasPreviousCall = true;
  state.previousCallTime = time;

  result.hasT0 = state.hasT0;

  if(state.hasT0) {
    result.t0 = state.t0;
    result.shiftedTimeOverTc =
      (time - state.t0) / tc;
  }

  result.firstDropletDetected =
    state.firstDropletDetected;

  if(state.firstDropletDetected) {
    result.firstDropletTime = state.firstDropletTime;
    result.firstDropletTimeOverTc = state.firstDropletTime / tc;
    result.firstDropletRadiusOverRb = state.firstDropletRadius / bubbleRadius;
    result.firstDropletVelocity = state.firstDropletVelocity;
    result.firstDropletCa = state.firstDropletCa;
  }

  return result;
}

void PrintBurstingBubbleMetrics(
  const BurstingBubbleMetrics& m,
  const unsigned iproc,
  const std::string& prefix) {

  using namespace BurstingBubbleDetail;

  if(iproc != 0u) return;

  const unsigned w = 24u;
  const std::string jetFile = prefix + "_jet.dat";
  const std::string dropFile = prefix + "_droplets.dat";

  const bool jetHeader = FileEmptyOrMissing(jetFile);

  std::ofstream jet(jetFile, std::ios::app);

  if(!jet) {
    throw std::runtime_error(
      "PrintBurstingBubbleMetrics: cannot open " + jetFile);
  }

  if(jetHeader) {
    jet
        << "#"
        << std::setw(w - 1) << "t"
        << std::setw(w) << "t/tc"
        << std::setw(w) << "(t-t0)/tc"
        << std::setw(w) << "z_tip"
        << std::setw(w) << "z_tip/Rb"
        << std::setw(w) << "V_tip"
        << std::setw(w) << "Ca_tip"
        << std::setw(w) << "V_tip_FD"
        << std::setw(w) << "N_components"
        << std::setw(w) << "N_drops"
        << std::setw(w) << "t_d1"
        << std::setw(w) << "t_d1/tc"
        << std::setw(w) << "R_d1/Rb"
        << std::setw(w) << "V_d1"
        << std::setw(w) << "Ca_d1"
        << '\n';
  }

  jet
      << std::scientific
      << std::setprecision(16)
      << std::setw(w) << m.time
      << std::setw(w) << m.timeOverTc
      << std::setw(w) << m.shiftedTimeOverTc
      << std::setw(w) << m.jetTipZ
      << std::setw(w) << m.jetTipZOverRb
      << std::setw(w) << m.jetTipVelocity
      << std::setw(w) << m.jetTipCa
      << std::setw(w) << m.jetTipVelocityFD
      << std::setw(w) << m.nLiquidComponents
      << std::setw(w) << m.nDroplets
      << std::setw(w) << m.firstDropletTime
      << std::setw(w) << m.firstDropletTimeOverTc
      << std::setw(w) << m.firstDropletRadiusOverRb
      << std::setw(w) << m.firstDropletVelocity
      << std::setw(w) << m.firstDropletCa
      << '\n';

  const bool dropHeader = FileEmptyOrMissing(dropFile);

  std::ofstream drop(dropFile, std::ios::app);

  if(!drop) {
    throw std::runtime_error(
      "PrintBurstingBubbleMetrics: cannot open " + dropFile);
  }

  if(dropHeader) {
    drop
        << "#"
        << std::setw(w - 1) << "t"
        << std::setw(w) << "t/tc"
        << std::setw(w) << "drop_id"
        << std::setw(w) << "z_cm"
        << std::setw(w) << "z_cm/Rb"
        << std::setw(w) << "V_z"
        << std::setw(w) << "Ca_d"
        << std::setw(w) << "volume"
        << std::setw(w) << "R_eq"
        << std::setw(w) << "R_eq/Rb"
        << '\n';
  }

  for(const BBDropletMetrics& d : m.droplets) {
    drop
        << std::scientific
        << std::setprecision(16)
        << std::setw(w) << m.time
        << std::setw(w) << m.timeOverTc
        << std::setw(w) << d.id
        << std::setw(w) << d.zCm
        << std::setw(w) << d.zCmOverRb
        << std::setw(w) << d.verticalVelocity
        << std::setw(w) << d.capillaryNumber
        << std::setw(w) << d.volume
        << std::setw(w) << d.radius
        << std::setw(w) << d.radiusOverRb
        << '\n';
  }
}
