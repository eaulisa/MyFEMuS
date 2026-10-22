#include <cmath>
#include <ctime>
#include <iostream>
#include <vector>

void AssembleStabilizationTerms(MultiLevelProblem& ml_prob) {
  const std::clock_t start_time = std::clock();

  TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem>("NS");
  const unsigned level = my_nnlin_impl_sys.GetLevelToAssemble();
  MultiLevelSolution* mlSol = ml_prob._ml_sol;
  Solution* mysolution = mlSol->GetSolutionLevel(level);
  LinearEquationSolver* myLinEqSolver = my_nnlin_impl_sys._LinSolver[level];
  Mesh* msh = ml_prob._ml_msh->GetLevel(level);
  SparseMatrix* myKK = myLinEqSolver->_KK;
  NumericVector* myRES = myLinEqSolver->_RES;

  const unsigned dim = msh->GetDimension();
  const unsigned dim2 = 3 * (dim - 1);
  const unsigned iproc = msh->processor_id();
  const double dt = my_nnlin_impl_sys.GetIntervalTime();
  const bool linear = ml_prob.GetMultiphaseParams().linear;
  const TimeDiscretization td = ml_prob.GetMultiphaseParams().td;
  const bool axisymmetric = ml_prob.GetMultiphaseParams().axisymmetric;
  const double axisCoordinate = ml_prob.GetMultiphaseParams().axisCoordinate;
  const double cold = (td == TimeDiscretization::CrankNicholson) ? 0.5 : 0.;
  const double cnew = (td == TimeDiscretization::CrankNicholson) ? 0.5 : 1.;
  const double derivativeW = linear ? 0. : cnew;
  const double advectionOldWeight = linear ? cold - cnew : cold;

  const MultiphasePhysicalProperties properties = ml_prob.GetMultiphaseParams().properties;
  std::vector<double> g(dim, 0.);
  g[dim - 1] = properties.gravity;

  const char varname[3][5] = {"U", "V", "W"};
  std::vector<unsigned> indexSolV(dim);
  std::vector<unsigned> indexPdeV(dim);
  for(unsigned k = 0; k < dim; ++k) {
    indexSolV[k] = mlSol->GetIndex(&varname[k][0]);
    indexPdeV[k] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[k][0]);
  }

  const unsigned solTypeV = mlSol->GetSolutionType(&varname[0][0]);
  const unsigned indexSolP1 = mlSol->GetIndex("P1");
  const unsigned indexPdeP1 = my_nnlin_impl_sys.GetSolPdeIndex("P1");
  const unsigned indexSolP2 = mlSol->GetIndex("P2");
  const unsigned indexPdeP2 = my_nnlin_impl_sys.GetSolPdeIndex("P2");
  const unsigned indexSolC = mlSol->GetIndex("C");
  const unsigned solTypeP = mlSol->GetSolutionType("P1");

  std::vector<std::vector<unsigned>> secondDerivativeIndex(dim, std::vector<unsigned>(dim, 0));
  for(unsigned k = 0; k < dim; ++k) {
    for(unsigned l = 0; l < dim; ++l) {
      if(k == l) secondDerivativeIndex[k][l] = k;
      else if(k + l == 1) secondDerivativeIndex[k][l] = dim;
      else if(k + l == 2) secondDerivativeIndex[k][l] = dim + 2;
      else secondDerivativeIndex[k][l] = dim + 1;
    }
  }

  std::vector<std::vector<double>> solV(dim);
  std::vector<std::vector<double>> solVOld(dim);
  std::vector<std::vector<double>> vx(dim);
  std::vector<double> solP;
  std::vector<unsigned> sysDofsAll;
  std::vector<double> rhs;
  std::vector<double> Jac;
  std::vector<double> phi;
  std::vector<double> gradPhi;
  std::vector<double> nablaPhi;
  std::vector<double> phiP;
  std::vector<double> gradPhiP;
  double weight = 0.;

  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; ++iel) {
    const double C = (*mysolution->_Sol[indexSolC])(iel);
    if(C != 0. && C != 1.) continue;

    const unsigned ielt = msh->GetElementType(iel);
    const unsigned nDofsV = msh->GetElementDofNumber(iel, solTypeV);
    const unsigned nDofsP = msh->GetElementDofNumber(iel, solTypeP);
    const unsigned nDofsVelocity = dim * nDofsV;
    const unsigned nDofsAll = nDofsVelocity + nDofsP;
    const unsigned indexSolP = (C == 1.) ? indexSolP1 : indexSolP2;
    const unsigned indexPdeP = (C == 1.) ? indexPdeP1 : indexPdeP2;
    const double rho = properties.rho1 * C + properties.rho2 * (1. - C);
    const double mu = properties.mu1 * C + properties.mu2 * (1. - C);

    sysDofsAll.resize(nDofsAll);
    rhs.assign(nDofsAll, 0.);
    Jac.assign(nDofsAll * nDofsAll, 0.);
    solP.resize(nDofsP);
    for(unsigned k = 0; k < dim; ++k) {
      solV[k].resize(nDofsV);
      solVOld[k].resize(nDofsV);
      vx[k].resize(nDofsV);
    }

    for(unsigned i = 0; i < nDofsV; ++i) {
      const unsigned idof = msh->GetSolutionDof(i, iel, solTypeV);
      const unsigned idofX = msh->GetSolutionDof(i, iel, 2);
      for(unsigned k = 0; k < dim; ++k) {
        solV[k][i] = (*mysolution->_Sol[indexSolV[k]])(idof);
        solVOld[k][i] = (*mysolution->_SolOld[indexSolV[k]])(idof);
        sysDofsAll[k * nDofsV + i] = myLinEqSolver->GetSystemDof(indexSolV[k], indexPdeV[k], i, iel);
        vx[k][i] = (*msh->_topology->_Sol[k])(idofX);
      }
    }

    for(unsigned i = 0; i < nDofsP; ++i) {
      const unsigned idof = msh->GetSolutionDof(i, iel, solTypeP);
      solP[i] = (*mysolution->_Sol[indexSolP])(idof);
      sysDofsAll[nDofsVelocity + i] = myLinEqSolver->GetSystemDof(indexSolP, indexPdeP, i, iel);
    }

    for(unsigned ig = 0; ig < msh->_finiteElement[ielt][solTypeV]->GetGaussPointNumber(); ++ig) {
      msh->_finiteElement[ielt][solTypeP]->Jacobian(vx, ig, weight, phiP, gradPhiP);
      msh->_finiteElement[ielt][solTypeV]->Jacobian(vx, ig, weight, phi, gradPhi, nablaPhi);

      std::vector<std::vector<double>> jacobianMatrix;
      std::vector<std::vector<double>> jacobianInverse;
      msh->_finiteElement[ielt][solTypeV]->GetJacobianMatrix(vx, ig, weight, jacobianMatrix, jacobianInverse);

      std::vector<double> solVg(dim, 0.);
      std::vector<double> solVgOld(dim, 0.);
      std::vector<std::vector<double>> gradSolVg(dim, std::vector<double>(dim, 0.));
      std::vector<std::vector<double>> gradSolVgOld(dim, std::vector<double>(dim, 0.));
      std::vector<std::vector<double>> DeltaSolVg(dim, std::vector<double>(dim2, 0.));
      std::vector<std::vector<double>> DeltaSolVgOld(dim, std::vector<double>(dim2, 0.));

      double x_gss = 0.;

      for(unsigned i = 0; i < nDofsV; ++i) {
        x_gss += phi[i] * vx[0][i];
      }

      const double r = axisymmetric ? x_gss - axisCoordinate : 1.;
      const double axisWeight = axisymmetric ? r : 1.;
      const double volumeWeight = weight * axisWeight;

      if(axisymmetric && r < -1.e-12) {
        throw std::runtime_error("Negative axisymmetric radius in stabilization");
      }

      for(unsigned i = 0; i < nDofsV; ++i) {
        for(unsigned k = 0; k < dim; ++k) {
          solVg[k] += phi[i] * solV[k][i];
          solVgOld[k] += phi[i] * solVOld[k][i];
          for(unsigned j = 0; j < dim; ++j) {
            gradSolVg[k][j] += gradPhi[i * dim + j] * solV[k][i];
            gradSolVgOld[k][j] += gradPhi[i * dim + j] * solVOld[k][i];
          }
          for(unsigned j = 0; j < dim2; ++j) {
            DeltaSolVg[k][j] += nablaPhi[i * dim2 + j] * solV[k][i];
            DeltaSolVgOld[k][j] += nablaPhi[i * dim2 + j] * solVOld[k][i];
          }
        }
      }

      std::vector<double> w(dim, 0.);
      for(unsigned k = 0; k < dim; ++k) {
        w[k] = linear ? solVgOld[k] : cnew * solVg[k] + cold * solVgOld[k];
      }

      std::vector<double> gradSolPg(dim, 0.);
      if(solTypeP != 3) {
        for(unsigned i = 0; i < nDofsP; ++i) {
          for(unsigned k = 0; k < dim; ++k) {
            gradSolPg[k] += solP[i] * gradPhiP[i * dim + k];
          }
        }
      }

      std::vector<std::vector<double>> G(dim, std::vector<double>(dim, 0.));
      for(unsigned i = 0; i < dim; ++i) {
        for(unsigned j = 0; j < dim; ++j) {
          for(unsigned k = 0; k < dim; ++k) {
            G[i][j] += jacobianInverse[k][i] * jacobianInverse[k][j];
          }
        }
      }

      double denom = std::pow(2. * rho / dt, 2.);
      double traceG = 0.;
      std::vector<double> metricVelocity(dim, 0.);
      for(unsigned i = 0; i < dim; ++i) {
        traceG += G[i][i];
        for(unsigned j = 0; j < dim; ++j) {
          denom += rho * rho * w[i] * G[i][j] * w[j] + 36. * mu * mu * G[i][j] * G[i][j];
          metricVelocity[i] += G[i][j] * w[j];
        }
      }

      const double tauM = 1. / std::sqrt(denom);
      const double tauC = 1. / (tauM * traceG);
      std::vector<double> supgPhi(nDofsV, 0.);
      std::vector<double> tauMsupgPhi(nDofsV, 0.);
      std::vector<double> laplacianPhi(nDofsV, 0.);
      std::vector<double> advectionPhi(nDofsV, 0.);
      for(unsigned i = 0; i < nDofsV; ++i) {
        for(unsigned k = 0; k < dim; ++k) {
          supgPhi[i] += rho * w[k] * gradPhi[i * dim + k];
          laplacianPhi[i] += nablaPhi[i * dim2 + k];
          advectionPhi[i] += (linear ? solVgOld[k] : solVg[k]) * gradPhi[i * dim + k];
        }
        tauMsupgPhi[i] = tauM * supgPhi[i];
      }

      double divVg = 0.;
      for(unsigned k = 0; k < dim; ++k) {
        divVg += gradSolVg[k][k];
      }

      if(axisymmetric) {
        divVg += solVg[0] / r;
      }

      std::vector<double> rM(dim, 0.);
      for(unsigned k = 0; k < dim; ++k) {
        double advectionNew = 0.;
        double advectionOld = 0.;
        double diffusionNew = 0.;
        double diffusionOld = 0.;
        for(unsigned j = 0; j < dim; ++j) {
          if(linear) advectionNew += rho * (solVgOld[j] * gradSolVg[k][j] + solVg[j] * gradSolVgOld[k][j]);
          else advectionNew += rho * solVg[j] * gradSolVg[k][j];
          advectionOld += rho * solVgOld[j] * gradSolVgOld[k][j];
          diffusionNew -= mu * (DeltaSolVg[k][j] + DeltaSolVg[j][secondDerivativeIndex[k][j]]);
          diffusionOld -= mu * (DeltaSolVgOld[k][j] + DeltaSolVgOld[j][secondDerivativeIndex[k][j]]);
        }

        if(axisymmetric) {
          if(k == 0) {
            diffusionNew -= mu * (2. * gradSolVg[0][0] / r - 2. * solVg[0] / (r * r));
            diffusionOld -= mu * (2. * gradSolVgOld[0][0] / r - 2. * solVgOld[0] / (r * r));
          }

          if(k == 1) {
            diffusionNew -= mu * (gradSolVg[1][0] + gradSolVg[0][1]) / r;
            diffusionOld -= mu * (gradSolVgOld[1][0] + gradSolVgOld[0][1]) / r;
          }
        }

        rM[k] = rho * (solVg[k] - solVgOld[k]) / dt + cnew * advectionNew + advectionOldWeight * advectionOld + cnew * diffusionNew + cold * diffusionOld + gradSolPg[k] - rho * g[k];
      }

      std::vector<std::vector<double>> derivativeRM(dim, std::vector<double>(nDofsVelocity, 0.));
      std::vector<double> derivativeTauM(nDofsVelocity, 0.);
      std::vector<double> derivativeTauC(nDofsVelocity, 0.);
      for(unsigned l = 0; l < dim; ++l) {
        for(unsigned j = 0; j < nDofsV; ++j) {
          const unsigned column = l * nDofsV + j;
          derivativeTauM[column] = -derivativeW * rho * rho * tauM * tauM * tauM * phi[j] * metricVelocity[l];
          derivativeTauC[column] = -tauC * derivativeTauM[column] / tauM;
          for(unsigned k = 0; k < dim; ++k) {
            derivativeRM[k][column] = cnew * rho * phi[j] * (linear ? gradSolVgOld[k][l] : gradSolVg[k][l]) - cnew * mu * nablaPhi[j * dim2 + secondDerivativeIndex[k][l]];
            if(k == l) derivativeRM[k][column] += rho * phi[j] / dt + cnew * rho * advectionPhi[j] - cnew * mu * laplacianPhi[j];

            if(axisymmetric && k == 0 && l == 0) {
              derivativeRM[k][column] -= cnew * mu * (2. * gradPhi[j * dim + 0] / r - 2. * phi[j] / (r * r));
            }
            if(axisymmetric && k == 1 && l == 1) {
              derivativeRM[k][column] -= cnew * mu * gradPhi[j * dim + 0] / r;
            }
            if(axisymmetric && k == 1 && l == 0) {
              derivativeRM[k][column] -= cnew * mu * gradPhi[j * dim + 1] / r;
            }
          }
        }
      }

      for(unsigned i = 0; i < nDofsV; ++i) {
        for(unsigned k = 0; k < dim; ++k) {
          double divTest = gradPhi[i * dim + k];

          if(axisymmetric && k == 0) {
            divTest += phi[i] / r;
          }

          const unsigned row = k * nDofsV + i;
          rhs[row] -= (rM[k] * tauMsupgPhi[i] + tauC * divVg * divTest) * volumeWeight;
          for(unsigned l = 0; l < dim; ++l) {
            for(unsigned j = 0; j < nDofsV; ++j) {
              const unsigned column = l * nDofsV + j;
              const double derivativeSupg = derivativeTauM[column] * supgPhi[i] + tauM * rho * derivativeW * phi[j] * gradPhi[i * dim + l];
              double derivativeDiv = gradPhi[j * dim + l];

              if(axisymmetric && l == 0) {
                derivativeDiv += phi[j] / r;
              }

              Jac[row * nDofsAll + column] += (derivativeRM[k][column] * tauMsupgPhi[i] + rM[k] * derivativeSupg + (derivativeTauC[column] * divVg + tauC * derivativeDiv) * divTest) * volumeWeight;
            }
          }
          if(solTypeP != 3) {
            for(unsigned j = 0; j < nDofsP; ++j) {
              Jac[row * nDofsAll + nDofsVelocity + j] += tauMsupgPhi[i] * gradPhiP[j * dim + k] * volumeWeight;
            }
          }
        }
      }

      if(solTypeP != 3) {
        for(unsigned i = 0; i < nDofsP; ++i) {
          const unsigned row = nDofsVelocity + i;
          double pressureResidual = 0.;
          for(unsigned k = 0; k < dim; ++k) {
            pressureResidual += gradPhiP[i * dim + k] * rM[k];
          }
          rhs[row] -= tauM * pressureResidual * volumeWeight;
          for(unsigned column = 0; column < nDofsVelocity; ++column) {
            double pressureDerivative = 0.;
            for(unsigned k = 0; k < dim; ++k) {
              pressureDerivative += gradPhiP[i * dim + k] * derivativeRM[k][column];
            }
            Jac[row * nDofsAll + column] += (derivativeTauM[column] * pressureResidual + tauM * pressureDerivative) * volumeWeight;
          }
          for(unsigned j = 0; j < nDofsP; ++j) {
            double pressureDerivative = 0.;
            for(unsigned k = 0; k < dim; ++k) {
              pressureDerivative += gradPhiP[i * dim + k] * gradPhiP[j * dim + k];
            }
            Jac[row * nDofsAll + nDofsVelocity + j] += tauM * pressureDerivative * volumeWeight;
          }
        }
      }
    }

    myRES->add_vector_blocked(rhs, sysDofsAll);
    myKK->add_matrix_blocked(Jac, sysDofsAll, sysDofsAll);
  }

  std::cout << "Stabilization Assembly time = " << static_cast<double>(std::clock() - start_time) / CLOCKS_PER_SEC << std::flush << std::endl;
}