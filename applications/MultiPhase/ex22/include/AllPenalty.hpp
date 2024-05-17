#ifndef ALLPENALTY_HPP_INCLUDED
#define ALLPENALTY_HPP_INCLUDED

void BuildAllPenalty(const double &mu1, const double &mu2, const double &rho1, const double &rho2, const double &dt, const double &Ci, const double &Cj,
                     const double &h, const double &h2, const double &h3, const double &weight, const std::vector<double> &normal, const bool &ghostPenalty,
                     const std::vector<std::vector<adept::adouble>> &solVi, const std::vector<std::vector<adept::adouble>> &solVj,
                     const std::vector<double> &gradPhiVi, std::vector<double> &gradPhiVj, const std::vector<double> &hessPhiVi, std::vector<double> &hessPhiVj,
                     const std::vector<std::vector<adept::adouble>> &solPi, const std::vector<std::vector<adept::adouble>> &solPj,
                     const std::vector<double> &gradPhiPi, std::vector<double> &gradPhiPj, const std::vector<double> &hessPhiPi, std::vector<double> &hessPhiPj,
                     std::vector<adept::adouble> &aResi, std::vector<adept::adouble> &aResj) {

  double C2 = 1;
  double PHIT1 = mu1 + C2 * rho1 * h2 / dt;
  double PHIT2 = mu2 + C2 * rho2 * h2 / dt;

  double PHITi = Ci * PHIT1 + (1. - Ci) * PHIT2;
  double PHITj = Cj * PHIT1 + (1. - Cj) * PHIT2;

  double gammaV = 0.1;
  double PHIT = gammaV * 0.5 * (PHITi + PHITj);

  double gammaP = 10;
  double PHIP[2] = {gammaP * h2 / PHIT1,gammaP * h2 / PHIT2};

  if (Ci < 1.0e-10   || Cj < 1.0e-10) PHIP[0] = 0.;
  if (Ci > 1. - 1.0e-10 || Cj > 1. - 1.0e-10) PHIP[1] = 0.;

  const unsigned &dim = solVi.size();
  const unsigned dim2 = 3 * (dim - 1);
  const unsigned &nDofsVi = solVi[0].size();
  const unsigned &nDofsVj = solVj[0].size();
  const unsigned &nDofsPi = solPi[0].size();
  const unsigned &nDofsPj = solPj[0].size();

  std::vector<unsigned> offsetVi(dim);
  std::vector<unsigned> offsetVj(dim);
  for (unsigned J = 0; J < dim; J++) {
    offsetVi[J] = J * nDofsVi;
    offsetVj[J] = J * nDofsVj;
  }

  std::vector<unsigned> offsetPi = {dim * nDofsVi, dim * nDofsVi + nDofsPi};
  std::vector<unsigned> offsetPj = {dim * nDofsVj, dim * nDofsVj + nDofsPj};

  std::vector<std::vector<adept::adouble>> gradSolVi(dim, std::vector<adept::adouble>(dim, 0));
  std::vector<std::vector<adept::adouble>> gradSolVj(dim, std::vector<adept::adouble>(dim, 0));

  for (unsigned J = 0; J < dim; J++) {
    for (unsigned K = 0; K < dim; K++) {
      for (unsigned i = 0; i < nDofsVi; i++) {
        gradSolVi[J][K] += solVi[J][i] * gradPhiVi[i * dim + K];
      }
      for (unsigned j = 0; j < nDofsVj; j++) {
        gradSolVj[J][K] += solVj[J][j] * gradPhiVj[j * dim + K];
      }
    }
  }
  for (unsigned J = 0; J < dim; J++) {
    for (unsigned K = 0; K < dim; K++) {
      for (unsigned i = 0; i < nDofsVi; i++) {
        aResi[offsetVi[J] + i] +=  PHIT * h * (gradSolVi[K][K] -  gradSolVj[K][K])  * gradPhiVi[i * dim + J] * weight;
      }
      for (unsigned j = 0; j < nDofsVj; j++) {
        aResj[offsetVj[J] + j] -=  PHIT * h * (gradSolVi[K][K] -  gradSolVj[K][K])  * gradPhiVj[j * dim + J] * weight;
      }
    }
  }

  std::vector<std::vector<adept::adouble>> gradSolPi(2, std::vector<adept::adouble>(dim, 0));
  std::vector<std::vector<adept::adouble>> gradSolPj(2, std::vector<adept::adouble>(dim, 0));

  for (unsigned  K = 0; K < 2; K++) gradSolPi[K].assign(dim, 0.);
  for (unsigned  K = 0; K < 2; K++) gradSolPj[K].assign(dim, 0.);

  for (unsigned J = 0; J < 2; J++) {
    for (unsigned K = 0; K < dim; K++) {
      for (unsigned i = 0; i < nDofsPi; i++) {
        gradSolPi[J][K] += solPi[J][i] * gradPhiPi[i * dim + K];
      }
      for (unsigned j = 0; j < nDofsPj; j++) {
        gradSolPj[J][K] += solPj[J][j] * gradPhiPj[j * dim + K];
      }
    }
  }
  for (unsigned J = 0; J < 2; J++) {
    for (unsigned K = 0; K < dim; K++) {
      for (unsigned i = 0; i < nDofsPi; i++) {
        aResi[offsetPi[J] + i] +=  PHIP[J] * h * (gradSolPi[J][K] -  gradSolPj[J][K])  * gradPhiPi[i * dim + K] * weight;
      }
      for (unsigned j = 0; j < nDofsPj; j++) {
        aResj[offsetPj[J] + j] -=  PHIP[J] * h * (gradSolPi[J][K] -  gradSolPj[J][K])  * gradPhiPj[j * dim + K] * weight;
      }
    }
  }


  if (ghostPenalty) {

    std::vector<adept::adouble> gradSolVdotNi(dim, 0.);
    std::vector<adept::adouble> gradSolVdotNj(dim, 0.);

    for (unsigned J = 0; J < dim; J++) {
      for (unsigned K = 0; K < dim; K++) {
        gradSolVdotNi[J] += gradSolVi[J][K] * normal[K];
        gradSolVdotNj[J] += gradSolVj[J][K] * normal[K];
      }
    }

    for (unsigned J = 0; J < dim; J++) {
      for (unsigned K = 0; K < dim; K++) {
        for (unsigned i = 0; i < nDofsVi; i++) {
          aResi[offsetVi[J] + i] +=  PHIT * h * (gradSolVdotNi[J] - gradSolVdotNj[J])  * gradPhiVi[i * dim + K] * normal[K] * weight;
        }
        for (unsigned j = 0; j < nDofsVj; j++) {
          aResj[offsetVj[J] + j] -=  PHIT * h * (gradSolVdotNi[J] - gradSolVdotNj[J])  * gradPhiVj[j * dim + K] * normal[K] * weight;
        }
      }
    }

    std::vector<adept::adouble> gradSolPdotNi(2, 0.);
    std::vector<adept::adouble> gradSolPdotNj(2, 0.);

    for (unsigned J = 0; J < 2; J++) {
      for (unsigned K = 0; K < dim; K++) {
        gradSolPdotNi[J] += gradSolPi[J][K] * normal[K];
        gradSolPdotNj[J] += gradSolPj[J][K] * normal[K];
      }
    }


    //std::cout << Ci << " " << Cj << " " << PHIT << " " << PHIT << std::endl;

    for (unsigned J = 0; J < 2; J++) {
      for (unsigned K = 0; K < dim; K++) {
        for (unsigned i = 0; i < nDofsPi; i++) {
          aResi[offsetPi[J] + i] +=  PHIP[J] * h * (gradSolPdotNi[J] - gradSolPdotNj[J]) * gradPhiPi[i * dim + K] * normal[K] * weight;
        }
        for (unsigned j = 0; j < nDofsPj; j++) {
          aResj[offsetPj[J] + j] -=  PHIP[J] * h * (gradSolPdotNi[J] - gradSolPdotNj[J]) * gradPhiPj[j * dim + K] * normal[K] * weight;
        }
      }
    }


    std::vector<std::vector<std::vector<adept::adouble>>> HessSolVi(dim, std::vector<std::vector<adept::adouble>>(dim, std::vector<adept::adouble>(dim, 0.)));
    std::vector<std::vector<std::vector<adept::adouble>>> HessSolVj(dim, std::vector<std::vector<adept::adouble>>(dim, std::vector<adept::adouble>(dim, 0.)));

    for (unsigned I = 0; I < dim; I++) {
      for (unsigned J = 0; J < dim; J++) {
        for (unsigned K = 0; K < dim; K++) {
          unsigned L;
          if (J == K) L = J;
          else if (1 == J + K) L = dim;    // xy
          else if (2 == J + K) L = dim + 2; // xz
          else if (3 == J + K) L = dim + 1; // yz

          for (unsigned i = 0; i < nDofsVi; i++) {
            HessSolVi[I][J][K] += solVi[I][i] * hessPhiVi[i * dim2 + L];
          }
          for (unsigned j = 0; j < nDofsVj; j++) {
            HessSolVj[I][J][K] += solVj[I][j] * hessPhiVj[j * dim2 + L];
          }
        }
      }
    }

    std::vector<std::vector<adept::adouble>> HessSolVdotNi(dim, std::vector<adept::adouble>(dim, 0.));
    std::vector<std::vector<adept::adouble>> HessSolVdotNj(dim, std::vector<adept::adouble>(dim, 0.));

    for (unsigned I = 0; I < dim; I++) {
      for (unsigned J = 0; J < dim; J++) {
        for (unsigned K = 0; K < dim; K++) {
          HessSolVdotNi[I][J] += HessSolVi[I][J][K] * normal[K];
          HessSolVdotNj[I][J] += HessSolVj[I][J][K] * normal[K];
        }
      }
    }

    for (unsigned I = 0; I < dim; I++) {
      unsigned J = I;
      for (unsigned K = 0; K < dim; K++) {
        unsigned L;
        if (J == K) L = J;
        else if (1 == J + K) L = dim;    // xy
        else if (2 == J + K) L = dim + 2; // xz
        else if (3 == J + K) L = dim + 1; // yz
        for (unsigned i = 0; i < nDofsVi; i++) {
          aResi[offsetVi[I] + i] +=  0.00 * PHIT * h3 * (HessSolVdotNi[I][J] - HessSolVdotNj[I][J]) * hessPhiVi[i * dim2 + L] * normal[K] * weight;
        }
        for (unsigned j = 0; j < nDofsVj; j++) {
          aResj[offsetVj[I] + j] -=  0.00 * PHIT * h3 * (HessSolVdotNi[I][J] - HessSolVdotNj[I][J]) * hessPhiVj[j * dim2 + L] * normal[K] * weight;
        }
      }
    }

    for (unsigned I = 0; I < dim; I++) {
      for (unsigned J = 0; J < dim; J++) {
        for (unsigned K = 0; K < dim; K++) {
          unsigned L;
          if (J == K) L = J;
          else if (1 == J + K) L = dim;    // xy
          else if (2 == J + K) L = dim + 2; // xz
          else if (3 == J + K) L = dim + 1; // yz
          for (unsigned i = 0; i < nDofsVi; i++) {
            aResi[offsetVi[I] + i] +=  0.00 * PHIT * h3 * (HessSolVdotNi[I][J] - HessSolVdotNj[I][J]) * hessPhiVi[i * dim2 + L] * normal[K] * weight;
          }
          for (unsigned j = 0; j < nDofsVj; j++) {
            aResj[offsetVj[I] + j] -=  0.00 * PHIT * h3 * (HessSolVdotNi[I][J] - HessSolVdotNj[I][J]) * hessPhiVj[j * dim2 + L] * normal[K] * weight;
          }
        }
      }
    }

    std::vector<std::vector<std::vector<adept::adouble>>> HessSolPi(2, std::vector<std::vector<adept::adouble>>(dim, std::vector<adept::adouble>(dim, 0.)));
    std::vector<std::vector<std::vector<adept::adouble>>> HessSolPj(2, std::vector<std::vector<adept::adouble>>(dim, std::vector<adept::adouble>(dim, 0.)));

    for (unsigned I = 0; I < 2; I++) {
      for (unsigned J = 0; J < dim; J++) {
        for (unsigned K = 0; K < dim; K++) {
          unsigned L;
          if (J == K) L = J;
          else if (1 == J + K) L = dim;    // xy
          else if (2 == J + K) L = dim + 2; // xz
          else if (3 == J + K) L = dim + 1; // yz

          for (unsigned i = 0; i < nDofsPi; i++) {
            HessSolPi[I][J][K] += solPi[I][i] * hessPhiPi[i * dim2 + L];
          }
          for (unsigned j = 0; j < nDofsPj; j++) {
            HessSolPj[I][J][K] += solPj[I][j] * hessPhiPj[j * dim2 + L];
          }
        }
      }
    }

    std::vector<std::vector<adept::adouble>> HessSolPdotNi(2, std::vector<adept::adouble>(dim, 0.));
    std::vector<std::vector<adept::adouble>> HessSolPdotNj(2, std::vector<adept::adouble>(dim, 0.));

    for (unsigned I = 0; I < 2; I++) {
      for (unsigned J = 0; J < dim; J++) {
        for (unsigned K = 0; K < dim; K++) {
          HessSolPdotNi[I][J] += HessSolPi[I][J][K] * normal[K];
          HessSolPdotNj[I][J] += HessSolPj[I][J][K] * normal[K];
        }
      }
    }

    for (unsigned I = 0; I < 2; I++) {
      for (unsigned J = 0; J < dim; J++) {
        for (unsigned K = 0; K < dim; K++) {
          unsigned L;
          if (J == K) L = J;
          else if (1 == J + K) L = dim;    // xy
          else if (2 == J + K) L = dim + 2; // xz
          else if (3 == J + K) L = dim + 1; // yz

          for (unsigned i = 0; i < nDofsPi; i++) {
            aResi[offsetPi[I] + i] +=  PHIP[I] * h3 * (HessSolPdotNi[I][J] - HessSolPdotNj[I][J]) * hessPhiPi[i * dim2 + L] * normal[K] * weight;
          }
          for (unsigned j = 0; j < nDofsPj; j++) {
            aResj[offsetPj[I] + j] -=  PHIP[I] * h3 * (HessSolPdotNi[I][J] - HessSolPdotNj[I][J]) * hessPhiPj[j * dim2 + L] * normal[K] * weight;
          }
        }
      }
    }
  }
}


void AssembleAllPenalty(MultiLevelProblem& ml_prob) {

  //this function works both for fluid and solid ghost penalty, the boolean fluid switches between the two

  clock_t start_time;

  //pointers and references

  LinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<LinearImplicitSystem> ("NS");
  const unsigned  level = my_nnlin_impl_sys.GetLevelToAssemble();
  MultiLevelSolution* mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* mysolution = mlSol->GetSolutionLevel(level);     // pointer to the solution (level) object

  LinearEquationSolver* myLinEqSolver = my_nnlin_impl_sys._LinSolver[level];  // pointer to the equation (level) object

  Mesh* msh = ml_prob._ml_msh->GetLevel(level);     // pointer to the mesh (level) object
  elem* el = msh->el;   // pointer to the elem object in msh (level)

  SparseMatrix* myKK = myLinEqSolver->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* myRES =  myLinEqSolver->_RES;  // pointer to the global residual vector object in pdeSys (level)

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  const unsigned dim = msh->GetDimension();
  const unsigned dim2 = 3 * (dim - 1);


  // data
  unsigned iproc  = msh->processor_id();
  unsigned nprocs  = msh->n_processors();

  std::vector<std::vector<adept::adouble>> solVi(dim);
  std::vector<std::vector<adept::adouble>> solVj(dim);
  std::vector<std::vector<double>> solVjD(dim);

  std::vector<std::vector<adept::adouble>> solPi(2);
  std::vector<std::vector<adept::adouble>> solPj(2);
  std::vector<std::vector<double>> solPjD(2);

  std::vector < std::vector < adept::adouble > > gradSolVi(dim);
  std::vector < std::vector < adept::adouble > > gradSolVj(dim);

  std::vector < std::vector < adept::adouble > > gradSolPi(2);
  std::vector < std::vector < adept::adouble > > gradSolPj(2);

  std::vector<adept::adouble> aResi;
  std::vector<adept::adouble> aResj;

  std::vector<double> rhsi;
  std::vector<double> rhsj;

  vector < double > Jac;

  std::vector <unsigned> sysDofsi;
  std::vector <unsigned> sysDofsj;

  double weight;
  std::vector < double > phi;
  std::vector < double> gradPhi;

  double weightPi;
  std::vector<double> phiPi, gradPhiPi, hessPhiPi;

  double weightPj;
  std::vector<double> phiPj, gradPhiPj, hessPhiPj;

  double weightVi;
  std::vector<double> phiVi, gradPhiVi, hessPhiVi;

  double weightVj;
  std::vector<double> phiVj, gradPhiVj, hessPhiVj;


  std::vector <std::vector < double> > xi(dim);
  std::vector <std::vector < double> > xj(dim);
  std::vector <std::vector < double > > faceVx(dim);


  //variable-name handling

  std::vector<unsigned> indexSolV(dim);
  std::vector<unsigned> indexPdeV(dim);
  indexSolV[0] = mlSol->GetIndex("U");
  indexSolV[1] = mlSol->GetIndex("V");
  indexPdeV[0] = my_nnlin_impl_sys.GetSolPdeIndex("U");
  indexPdeV[1] = my_nnlin_impl_sys.GetSolPdeIndex("V");
  unsigned solTypeV = mlSol->GetSolutionType("U");

  std::vector<unsigned> indexSolP(2);
  indexSolP[0] = mlSol->GetIndex("P1");
  indexSolP[1] = mlSol->GetIndex("P2");
  std::vector<unsigned> indexPdeP(2);
  indexPdeP[0] = my_nnlin_impl_sys.GetSolPdeIndex("P1");
  indexPdeP[1] = my_nnlin_impl_sys.GetSolPdeIndex("P2");
  unsigned solTypeP = mlSol->GetSolutionType("P1");

  if (solTypeP > 2 || solTypeV > 2) {
    std::cout << "error the AssembleAllPenalty works only for Lagrangian Solutions\n";
    abort();
  }

  unsigned indexSolC = mlSol->GetIndex("C0");

  unsigned solTypeX = 2;
  start_time = clock();

  std::vector < std::vector < std::vector <double > > > aPi(3);
  std::vector < std::vector < std::vector <double > > > aPj(3);

  for (unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielType = msh->GetElementType(iel);
    double Ci = (*mysolution->_Sol[indexSolC])(iel);

    unsigned nDofsVi = msh->GetElementDofNumber(iel, solTypeV);
    unsigned nDofsPi = msh->GetElementDofNumber(iel, solTypeP);

    for (unsigned K = 0; K < dim; K++) solVi[K].resize(nDofsVi);
    for (unsigned K = 0; K < 2; K++) solPi[K].resize(nDofsPi);
    sysDofsi.resize(dim * nDofsVi + 2 * nDofsPi);

    for (unsigned i = 0; i < nDofsVi; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, solTypeV);
      for (unsigned K = 0; K < dim; K++) {
        solVi[K][i] = (*mysolution->_Sol[indexSolV[K]])(idof);
        sysDofsi[K * nDofsVi + i] = myLinEqSolver->GetSystemDof(indexSolV[K], indexPdeV[K], i, iel);
      }
    }

    for (unsigned i = 0; i < nDofsPi; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, solTypeP);
      for (unsigned K = 0; K < 2; K++) {
        solPi[K][i] = (*mysolution->_Sol[indexSolP[K]])(idof);
        sysDofsi[dim * nDofsVi + K * nDofsPi + i] = myLinEqSolver->GetSystemDof(indexSolP[K], indexPdeP[K], i, iel);
      }
    }

    unsigned nDofsXi = msh->GetElementDofNumber(iel, solTypeX);    // number of solution element dofs
    for (unsigned  K = 0; K < dim; K++) xi[K].resize(nDofsXi);
    for (unsigned i = 0; i < nDofsXi; i++) {
      unsigned idofX = msh->GetSolutionDof(i, iel, solTypeX);
      for (unsigned  K = 0; K < dim; K++) {
        xi[K][i] = (*msh->_topology->_Sol[K])(idofX);
      }
    }
    double h11 = (xi[0][2] - xi[0][0]);
    double h12 = (xi[1][2] - xi[1][0]);


    for (unsigned ktype = 0; ktype < solTypeP + 1; ktype++) {
      ProjectNodalToPolynomialCoefficients(aPi[ktype], xi, ielType, ktype);
    }

    for (unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
      int jel = el->GetFaceElementIndex(iel, iface) - 1;

      if (jel >= 0) { // iface is not a boundary of the domain
        unsigned jproc = msh->IsdomBisectionSearch(jel, 3);
        if (jproc == iproc) {

          if (jel > iel) {

            short unsigned jelType = msh->GetElementType(jel);
            double Cj = (*mysolution->_Sol[indexSolC])(jel);

            bool ghostPenalty = ((Ci > 0 && Ci < 1) || (Cj > 0 && Cj < 1)) ? true : false;

            unsigned nDofsVj = msh->GetElementDofNumber(jel, solTypeV);
            unsigned nDofsPj = msh->GetElementDofNumber(jel, solTypeP);

            for (unsigned K = 0; K < dim; K++) solVj[K].resize(nDofsVj);
            for (unsigned K = 0; K < 2; K++) solPj[K].resize(nDofsPj);
            sysDofsj.resize(dim * nDofsVj + 2 * nDofsPj);

            for (unsigned j = 0; j < nDofsVj; j++) {
              unsigned jdof = msh->GetSolutionDof(j, jel, solTypeV);
              for (unsigned K = 0; K < dim; K++) {
                solVj[K][j] = (*mysolution->_Sol[indexSolV[K]])(jdof);
                sysDofsj[K * nDofsVj + j] = myLinEqSolver->GetSystemDof(indexSolV[K], indexPdeV[K], j, jel);
              }
            }

            for (unsigned j = 0; j < nDofsPj; j++) {
              unsigned jdof = msh->GetSolutionDof(j, jel, solTypeP);
              for (unsigned K = 0; K < 2; K++) {
                solPj[K][j] = (*mysolution->_Sol[indexSolP[K]])(jdof);
                sysDofsj[dim * nDofsVj + K * nDofsPj + j] = myLinEqSolver->GetSystemDof(indexSolP[K], indexPdeP[K], j, jel);
              }
            }

            unsigned nDofsXj = msh->GetElementDofNumber(jel, solTypeX);    // number of solution element dofs
            for (unsigned  k = 0; k < dim; k++) {
              xj[k].resize(nDofsXj);
            }
            for (unsigned j = 0; j < nDofsXj; j++) {
              unsigned jdofX = msh->GetSolutionDof(j, jel, solTypeX);
              for (unsigned  K = 0; K < dim; K++) xj[K][j] = (*msh->_topology->_Sol[K])(jdofX);
            }
            double h21 = (xj[0][2] - xj[0][0]);
            double h22 = (xj[1][2] - xj[1][0]);

            aResi.assign(dim * nDofsVi + 2 * nDofsPi, 0);
            aResj.assign(dim * nDofsVj + 2 * nDofsPj, 0);

            s.new_recording();

            const unsigned faceGeom = msh->GetElementFaceType(iel, iface);
            unsigned faceDofs = msh->GetElementFaceDofNumber(iel, iface, solTypeX);
            for (unsigned  K = 0; K < dim; K++) faceVx[K].resize(faceDofs);
            for (unsigned i = 0; i < faceDofs; i++) {
              unsigned inode = msh->GetLocalFaceVertexIndex(iel, iface, i);    // face-to-element local node mapping.
              for (unsigned K = 0; K < dim; K++) {
                faceVx[K][i] =  xi[K][inode]; // We extract the local coordinates on the face from local coordinates on the element.
              }
            }

            for (unsigned ktype = 0; ktype < solTypeP + 1; ktype++) {
              ProjectNodalToPolynomialCoefficients(aPj[ktype], xj, jelType, ktype);
            }

            for (unsigned ig = 0; ig  <  msh->_finiteElement[faceGeom][solTypeP]->GetGaussPointNumber(); ig++) {

              std::vector < double> normal;
              msh->_finiteElement[faceGeom][solTypeX]->JacobianSur(faceVx, ig, weight, phi, gradPhi, normal);

              double h = 0.5 * (fabs(h11 * normal[0] + h12 * normal[1]) + fabs(h21 * normal[0] + h22 * normal[1])); //characteristic lenght in normal direction
              double h2 = h * h;
              double h3 = h2 * h;

              std::vector< double > xg(dim, 0.); // physical coordinates of the face Gauss point
              for (unsigned i = 0; i < faceDofs; i++) {
                for (unsigned K = 0; K < dim; K++) {
                  xg[K] += phi[i] * faceVx[K][i];
                }
              }

              std::vector <double> etai;//local coordinates of the face gauss point with respect to iel
              GetClosestPointInReferenceElement(xi, xg, ielType, etai);

              bool inverseMapping = GetInverseMapping(solTypeP, ielType, aPi, xg, etai, 100);
              if (!inverseMapping) {
                std::cout << "InverseMapping1 failed at " << iel << " " << jel << " " << iface << std::endl;
              }

              std::vector <double> etaj;//local coordinates of the face gauss point with respect to jel
              GetClosestPointInReferenceElement(xj, xg, jelType, etaj);

              inverseMapping = GetInverseMapping(solTypeP, jelType, aPj, xg, etaj, 100);
              if (!inverseMapping) {
                std::cout << "InverseMapping2 failed at " << iel << " " << jel << " " << iface << std::endl;
              }

              msh->_finiteElement[ielType][solTypeP]->Jacobian(xi, etai, weightPi, phiPi, gradPhiPi, hessPhiPi);
              msh->_finiteElement[jelType][solTypeP]->Jacobian(xj, etaj, weightPj, phiPj, gradPhiPj, hessPhiPj);

              msh->_finiteElement[ielType][solTypeV]->Jacobian(xi, etai, weightVi, phiVi, gradPhiVi, hessPhiVi);
              msh->_finiteElement[jelType][solTypeV]->Jacobian(xj, etaj, weightVj, phiVj, gradPhiVj, hessPhiVj);

              BuildAllPenalty(mu1, mu2, rho1, rho2, dt, Ci, Cj, h, h2, h3, weight, normal, ghostPenalty,
                              solVi, solVj, gradPhiVi, gradPhiVj, hessPhiVi, hessPhiVj,
                              solPi, solPj, gradPhiPi, gradPhiPj, hessPhiPi, hessPhiPj,
                              aResi, aResj);
            }

            rhsi.resize(aResi.size());   //resize
            for (int i = 0; i < aResi.size(); i++) {
              rhsi[i] = -aResi[i].value();
            }
            myRES->add_vector_blocked(rhsi, sysDofsi);

            rhsj.resize(aResj.size());   //resize
            for (int i = 0; i < aResj.size(); i++) {
              rhsj[i] = -aResj[i].value();
            }
            myRES->add_vector_blocked(rhsj, sysDofsj);

            std::vector<unsigned> sysDofs(0);
            sysDofs.reserve(sysDofsi.size() + sysDofsj.size());
            sysDofs.insert(sysDofs.end(), sysDofsi.begin(), sysDofsi.end());
            sysDofs.insert(sysDofs.end(), sysDofsj.begin(), sysDofsj.end());

            s.clear_dependents();
            s.clear_independents();

            s.dependent(&aResi[0], aResi.size());
            s.dependent(&aResj[0], aResi.size());
            for (unsigned K = 0; K < dim; K++) s.independent(solVi[K].data(), solVi[K].size());
            for (unsigned K = 0; K < 2; K++) s.independent(solPi[K].data(), solPi[K].size());
            for (unsigned K = 0; K < dim; K++) s.independent(solVj[K].data(), solVj[K].size());
            for (unsigned K = 0; K < 2; K++) s.independent(solPj[K].data(), solPj[K].size());

            Jac.resize(sysDofs.size() * sysDofs.size()); //J22
            s.jacobian(&Jac[0], true);
            myKK->add_matrix_blocked(Jac, sysDofs, sysDofs);
            s.clear_dependents();
            s.clear_independents();

          }
        }
      }
    }
  }



  if (nprocs > 1) {
    for (unsigned kproc = 0; kproc < nprocs; kproc++) {
      for (unsigned iel = msh->_elementOffset[kproc]; iel < msh->_elementOffset[kproc + 1]; iel++) {
        unsigned nFaces;
        if (iproc == kproc) {
          nFaces = msh->GetElementFaceNumber(iel);
        }
        MPI_Bcast(&nFaces, 1, MPI_UNSIGNED, kproc, PETSC_COMM_WORLD);

        for (unsigned iface = 0; iface < nFaces; iface++) {

          int jel;
          if (iproc == kproc) {
            jel = el->GetFaceElementIndex(iel, iface) - 1;
          }
          MPI_Bcast(&jel, 1, MPI_INT, kproc, PETSC_COMM_WORLD);

          if (jel >= 0) { // iface is not a boundary of the domain
            unsigned jproc = msh->IsdomBisectionSearch(jel, 3);
            if (jproc != kproc  && (iproc == kproc || iproc == jproc)) {
              if (jel > iel) {
                std::vector < MPI_Request > reqs(2 * dim + 3);
                if (iproc == jproc) {

                  unsigned short jelType = msh->GetElementType(jel);
                  double Cj = (*mysolution->_Sol[indexSolC])(jel);

                  unsigned nDofsVj = msh->GetElementDofNumber(jel, solTypeV);
                  unsigned nDofsPj = msh->GetElementDofNumber(jel, solTypeP);

                  for (unsigned K = 0; K < dim; K++) solVjD[K].resize(nDofsVj);
                  for (unsigned K = 0; K < 2; K++) solPjD[K].resize(nDofsPj);
                  sysDofsj.resize(dim * nDofsVj + 2 * nDofsPj);

                  for (unsigned j = 0; j < nDofsVj; j++) {
                    unsigned jdof = msh->GetSolutionDof(j, jel, solTypeV);
                    for (unsigned K = 0; K < dim; K++) {
                      solVjD[K][j] = (*mysolution->_Sol[indexSolV[K]])(jdof);
                      sysDofsj[K * nDofsVj + j] = myLinEqSolver->GetSystemDof(indexSolV[K], indexPdeV[K], j, jel);
                    }
                  }

                  for (unsigned j = 0; j < nDofsPj; j++) {
                    unsigned jdof = msh->GetSolutionDof(j, jel, solTypeP);
                    for (unsigned K = 0; K < 2; K++) {
                      solPjD[K][j] = (*mysolution->_Sol[indexSolP[K]])(jdof);
                      sysDofsj[dim * nDofsVj + K * nDofsPj + j] = myLinEqSolver->GetSystemDof(indexSolP[K], indexPdeP[K], j, jel);
                    }
                  }

                  unsigned nDofsXj = msh->GetElementDofNumber(jel, solTypeX);    // number of solution element dofs
                  for (unsigned  K = 0; K < dim; K++) xj[K].resize(nDofsXj);

                  for (unsigned j = 0; j < nDofsXj; j++) {
                    unsigned jdofX = msh->GetSolutionDof(j, jel, solTypeX);
                    for (unsigned  K = 0; K < dim; K++) xj[K][j] = (*msh->_topology->_Sol[K])(jdofX);
                  }

                  //blocking send
                  MPI_Send(&jelType, 1, MPI_UNSIGNED_SHORT, kproc, 0, PETSC_COMM_WORLD);
                  MPI_Send(&Cj, 1, MPI_DOUBLE, kproc, 1, PETSC_COMM_WORLD);
                  MPI_Send(&nDofsVj, 1, MPI_UNSIGNED, kproc, 2, PETSC_COMM_WORLD);
                  MPI_Send(&nDofsPj, 1, MPI_UNSIGNED, kproc, 3, PETSC_COMM_WORLD);
                  MPI_Send(&nDofsXj, 1, MPI_UNSIGNED, kproc, 4, PETSC_COMM_WORLD);

                  // //non-blocking send
                  // for(unsigned K = 0; K < dim; K++) MPI_Isend(solVjD[K].data(), solVjD[K].size(), MPI_DOUBLE, kproc, 5 + K, PETSC_COMM_WORLD, &reqs[K]);
                  // for(unsigned K = 0; K < 2; K++) MPI_Isend(solPjD[K].data(), solPjD[K].size(), MPI_DOUBLE, kproc, 5 + dim + K, PETSC_COMM_WORLD, &reqs[dim + K]);
                  // MPI_Isend(sysDofsj.data(), sysDofsj.size(), MPI_UNSIGNED, kproc, 7 + dim, PETSC_COMM_WORLD, &reqs[2 + dim]);
                  // for(unsigned K = 0; K < dim; K++) MPI_Isend(xj[K].data(), xj[K].size(), MPI_DOUBLE, kproc, 8 + dim + K, PETSC_COMM_WORLD, &reqs[3 + dim + K]);
                  //
                  // MPI_Status status;
                  // for(unsigned m = 0; m < reqs.size(); m++) {
                  //   MPI_Wait(&reqs[m], &status);
                  // }


                  for (unsigned K = 0; K < dim; K++) MPI_Send(solVjD[K].data(), solVjD[K].size(), MPI_DOUBLE, kproc, 5 + K, PETSC_COMM_WORLD);
                  for (unsigned K = 0; K < 2; K++) MPI_Send(solPjD[K].data(), solPjD[K].size(), MPI_DOUBLE, kproc, 5 + dim + K, PETSC_COMM_WORLD);
                  MPI_Send(sysDofsj.data(), sysDofsj.size(), MPI_UNSIGNED, kproc, 7 + dim, PETSC_COMM_WORLD);
                  for (unsigned K = 0; K < dim; K++) MPI_Send(xj[K].data(), xj[K].size(), MPI_DOUBLE, kproc, 8 + dim + K, PETSC_COMM_WORLD);


                  //  if(false && jel == 9882) {
                  //
                  //   std::cerr << "AAAAAAAAAAAAA " << iel << " " << jel << std::endl;
                  //   std::cerr << "jType=" << jelType <<  " Cj = " << Cj << " nDofsVj=" << nDofsVj << " nDofsPj=" << nDofsVj << " nDofsXj=" << nDofsXj << std::endl;
                  //
                  //   std::cerr << solVj[0].size() << " "<<solVj[1].size()<<" "<<solPj[0].size()<<" "<<solPj[1].size()<<std::endl;
                  //
                  //
                  //
                  //   std::cerr << solVj[0][0] << " " << solVj[0][1] << " " << solVj[0][2] << " " << solVj[0][3] << std::endl;
                  //   std::cerr << solVj[1][0] << " " << solVj[1][1] << " " << solVj[1][2] << " " << solVj[1][3] << std::endl;
                  //
                  //   std::cerr << solPj[0][0] << " " << solPj[0][1] << " " << solPj[0][2] << " " << solPj[0][3] << std::endl;
                  //   std::cerr << solPj[1][0] << " " << solPj[1][1] << " " << solPj[1][2] << " " << solPj[1][3] << std::endl;
                  //
                  //   std::cerr << xj[0][0] << " " << xj[0][1] << " " << xj[0][2] << " " << xj[0][3] << std::endl;
                  //   std::cerr << xj[1][0] << " " << xj[1][1] << " " << xj[1][2] << " " << xj[1][3] << std::endl;
                  //
                  // }



                }

                if (iproc == kproc) {
                  unsigned short ielType = msh->GetElementType(iel);
                  double Ci = (*mysolution->_Sol[indexSolC])(iel);

                  unsigned nDofsVi = msh->GetElementDofNumber(iel, solTypeV);
                  unsigned nDofsPi = msh->GetElementDofNumber(iel, solTypeP);

                  for (unsigned K = 0; K < dim; K++) solVi[K].resize(nDofsVi);
                  for (unsigned K = 0; K < 2; K++) solPi[K].resize(nDofsPi);
                  sysDofsi.resize(dim * nDofsVi + 2 * nDofsPi);

                  for (unsigned i = 0; i < nDofsVi; i++) {
                    unsigned idof = msh->GetSolutionDof(i, iel, solTypeV);
                    for (unsigned K = 0; K < dim; K++) {
                      solVi[K][i] = (*mysolution->_Sol[indexSolV[K]])(idof);
                      sysDofsi[K * nDofsVi + i] = myLinEqSolver->GetSystemDof(indexSolV[K], indexPdeV[K], i, iel);
                    }
                  }

                  for (unsigned i = 0; i < nDofsPi; i++) {
                    unsigned idof = msh->GetSolutionDof(i, iel, solTypeP);
                    for (unsigned K = 0; K < 2; K++) {
                      solPi[K][i] = (*mysolution->_Sol[indexSolP[K]])(idof);
                      sysDofsi[dim * nDofsVi + K * nDofsPi + i] = myLinEqSolver->GetSystemDof(indexSolP[K], indexPdeP[K], i, iel);
                    }
                  }

                  unsigned nDofsXi = msh->GetElementDofNumber(iel, solTypeX);    // number of solution element dofs
                  for (unsigned  K = 0; K < dim; K++) xi[K].resize(nDofsXi);
                  for (unsigned i = 0; i < nDofsXi; i++) {
                    unsigned idofX = msh->GetSolutionDof(i, iel, solTypeX);
                    for (unsigned  K = 0; K < dim; K++) xi[K][i] = (*msh->_topology->_Sol[K])(idofX);
                  }

                  unsigned short jelType;
                  double Cj;
                  unsigned nDofsVj, nDofsPj, nDofsXj;
                  //blocking recv
                  MPI_Recv(&jelType, 1, MPI_UNSIGNED_SHORT, jproc, 0, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                  MPI_Recv(&Cj, 1, MPI_DOUBLE, jproc, 1, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                  MPI_Recv(&nDofsVj, 1, MPI_UNSIGNED, jproc, 2, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                  MPI_Recv(&nDofsPj, 1, MPI_UNSIGNED, jproc, 3, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                  MPI_Recv(&nDofsXj, 1, MPI_UNSIGNED, jproc, 4, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);


                  for (unsigned K = 0; K < dim; K++) solVj[K].resize(nDofsVj);
                  for (unsigned K = 0; K < dim; K++) solVjD[K].resize(nDofsVj);
                  for (unsigned K = 0; K < 2; K++) solPj[K].resize(nDofsPj);
                  for (unsigned K = 0; K < 2; K++) solPjD[K].resize(nDofsPj);
                  sysDofsj.resize(dim * nDofsVj + 2 * nDofsPj);
                  for (unsigned  K = 0; K < dim; K++) xj[K].resize(nDofsXj);

                  // // non-blocking recv
                  // for(unsigned K = 0; K < dim; K++) MPI_Irecv(solVj[K].data(), solVj[K].size(), MPI_DOUBLE, jproc, 5 + K, PETSC_COMM_WORLD, &reqs[K]);
                  // for(unsigned K = 0; K < 2; K++) MPI_Irecv(solPj[K].data(), solPj[K].size(), MPI_DOUBLE, jproc, 5 + dim + K, PETSC_COMM_WORLD, &reqs[dim + K]);
                  // MPI_Irecv(sysDofsj.data(), sysDofsj.size(), MPI_UNSIGNED, jproc, 7 + dim, PETSC_COMM_WORLD, &reqs[2 + dim]);
                  // for(unsigned K = 0; K < dim; K++) MPI_Irecv(xj[K].data(), xj[K].size(), MPI_DOUBLE, jproc, 8 + dim + K, PETSC_COMM_WORLD, &reqs[3 + dim + K]);
                  //
                  // MPI_Status status;
                  // for(unsigned m = 0; m < reqs.size(); m++) {
                  //   MPI_Wait(&reqs[m], &status);
                  // }


                  for (unsigned K = 0; K < dim; K++) MPI_Recv(solVjD[K].data(), solVjD[K].size(), MPI_DOUBLE, jproc, 5 + K, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                  for (unsigned K = 0; K < 2; K++) MPI_Recv(solPjD[K].data(), solPjD[K].size(), MPI_DOUBLE, jproc, 5 + dim + K, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                  MPI_Recv(sysDofsj.data(), sysDofsj.size(), MPI_UNSIGNED, jproc, 7 + dim, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                  for (unsigned K = 0; K < dim; K++) MPI_Recv(xj[K].data(), xj[K].size(), MPI_DOUBLE, jproc, 8 + dim + K, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);


                  for (unsigned K = 0; K < dim; K++) {
                    for (unsigned i = 0; i < solVj[K].size(); i++) {
                      solVj[K][i] = solVjD[K][i];
                    }
                  }

                  for (unsigned K = 0; K < 2; K++) {
                    for (unsigned i = 0; i < solPj[K].size(); i++) {
                      solPj[K][i] = solPjD[K][i];
                    }
                  }

                  bool ghostPenalty = ((Ci > 0 && Ci < 1) || (Cj > 0 && Cj < 1)) ? true : false;

                  double h11 = (xi[0][2] - xi[0][0]);
                  double h12 = (xi[1][2] - xi[1][0]);

                  double h21 = (xj[0][2] - xj[0][0]);
                  double h22 = (xj[1][2] - xj[1][0]);

                  for (unsigned ktype = 0; ktype < solTypeP + 1; ktype++) {
                    ProjectNodalToPolynomialCoefficients(aPi[ktype], xi, ielType, ktype);
                    ProjectNodalToPolynomialCoefficients(aPj[ktype], xj, jelType, ktype);
                  }

                  aResi.assign(dim * nDofsVi + 2 * nDofsPi, 0);
                  aResj.assign(dim * nDofsVj + 2 * nDofsPj, 0);

                  s.new_recording();

                  const unsigned faceGeom = msh->GetElementFaceType(iel, iface);
                  unsigned faceDofs = msh->GetElementFaceDofNumber(iel, iface, solTypeX);
                  for (unsigned  K = 0; K < dim; K++) faceVx[K].resize(faceDofs);
                  for (unsigned i = 0; i < faceDofs; i++) {
                    unsigned inode = msh->GetLocalFaceVertexIndex(iel, iface, i);    // face-to-element local node mapping.
                    for (unsigned K = 0; K < dim; K++) {
                      faceVx[K][i] =  xi[K][inode]; // We extract the local coordinates on the face from local coordinates on the element.
                    }
                  }

                  for (unsigned ig = 0; ig  <  msh->_finiteElement[faceGeom][solTypeP]->GetGaussPointNumber(); ig++) {

                    std::vector < double> normal;
                    msh->_finiteElement[faceGeom][solTypeX]->JacobianSur(faceVx, ig, weight, phi, gradPhi, normal);

                    double h = 0.5 * (fabs(h11 * normal[0] + h12 * normal[1]) + fabs(h21 * normal[0] + h22 * normal[1])); //characteristic lenght in normal direction
                    double h2 = h * h;
                    double h3 = h2 * h;

                    std::vector< double > xg(dim, 0.); // physical coordinates of the face Gauss point
                    for (unsigned i = 0; i < faceDofs; i++) {
                      for (unsigned K = 0; K < dim; K++) {
                        xg[K] += phi[i] * faceVx[K][i];
                      }
                    }

                    std::vector <double> etai;//local coordinates of the face gauss point with respect to iel
                    GetClosestPointInReferenceElement(xi, xg, ielType, etai);

                    bool inverseMapping = GetInverseMapping(solTypeP, ielType, aPi, xg, etai, 100);
                    if (!inverseMapping) {
                      std::cout << "InverseMapping1 failed at " << iel << " " << jel << " " << iface << std::endl;
                    }

                    std::vector <double> etaj;//local coordinates of the face gauss point with respect to jel
                    GetClosestPointInReferenceElement(xj, xg, jelType, etaj);

                    inverseMapping = GetInverseMapping(solTypeP, jelType, aPj, xg, etaj, 100);
                    if (!inverseMapping) {
                      std::cout << "InverseMapping2 failed at " << iel << " " << jel << " " << iface << std::endl;
                    }

                    msh->_finiteElement[ielType][solTypeP]->Jacobian(xi, etai, weightPi, phiPi, gradPhiPi);
                    msh->_finiteElement[jelType][solTypeP]->Jacobian(xj, etaj, weightPj, phiPj, gradPhiPj);

                    msh->_finiteElement[ielType][solTypeV]->Jacobian(xi, etai, weightVi, phiVi, gradPhiVi, hessPhiVi);
                    msh->_finiteElement[jelType][solTypeV]->Jacobian(xj, etaj, weightVj, phiVj, gradPhiVj, hessPhiVj);

                    BuildAllPenalty(mu1, mu2, rho1, rho2, dt, Ci, Cj, h, h2, h3, weight, normal, ghostPenalty,
                                    solVi, solVj, gradPhiVi, gradPhiVj, hessPhiVi, hessPhiVj,
                                    solPi, solPj, gradPhiPi, gradPhiPj, hessPhiPi, hessPhiPj,
                                    aResi, aResj);
                  }

                  rhsi.resize(aResi.size());   //resize
                  for (int i = 0; i < aResi.size(); i++) {
                    rhsi[i] = -aResi[i].value();
                  }
                  myRES->add_vector_blocked(rhsi, sysDofsi);

                  rhsj.resize(aResj.size());   //resize
                  for (int i = 0; i < aResj.size(); i++) {
                    rhsj[i] = -aResj[i].value();
                  }
                  myRES->add_vector_blocked(rhsj, sysDofsj);

                  std::vector<unsigned> sysDofs(0);
                  sysDofs.reserve(sysDofsi.size() + sysDofsj.size());
                  sysDofs.insert(sysDofs.end(), sysDofsi.begin(), sysDofsi.end());
                  sysDofs.insert(sysDofs.end(), sysDofsj.begin(), sysDofsj.end());

                  s.clear_dependents();
                  s.clear_independents();

                  s.dependent(&aResi[0], aResi.size());
                  s.dependent(&aResj[0], aResi.size());
                  for (unsigned K = 0; K < dim; K++) s.independent(solVi[K].data(), solVi[K].size());
                  for (unsigned K = 0; K < 2; K++) s.independent(solPi[K].data(), solPi[K].size());
                  for (unsigned K = 0; K < dim; K++) s.independent(solVj[K].data(), solVj[K].size());
                  for (unsigned K = 0; K < 2; K++) s.independent(solPj[K].data(), solPj[K].size());

                  Jac.resize(sysDofs.size() * sysDofs.size()); //J22
                  s.jacobian(&Jac[0], true);
                  myKK->add_matrix_blocked(Jac, sysDofs, sysDofs);
                  s.clear_dependents();
                  s.clear_independents();
                }
              }
            }
          }
        }
      }
    }
  }




  std::cout << "All Penalty Assembly time = " << static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC << std::endl;

}

#endif // ALLPENALTY_HPP_INCLUDED
