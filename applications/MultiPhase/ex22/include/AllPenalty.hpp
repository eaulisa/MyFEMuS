

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

  std::vector<std::vector<adept::adouble>> solPi(2);
  std::vector<std::vector<adept::adouble>> solPj(2);

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
  std::vector<double> phiPi, gradPhiPi;

  double weightPj;
  std::vector<double> phiPj, gradPhiPj;

  double weightVi;
  std::vector<double> phiVi, gradPhiVi, hessPhiVi;

  double weightVj;
  std::vector<double> phiVj, gradPhiVj, hessPhiVj;


  std::vector <std::vector < double> > xi(dim);
  std::vector <std::vector < double> > xj(dim);
  std::vector <std::vector < double > > faceVx(dim);

  std::cout.precision(10);

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

  if(solTypeP > 2 || solTypeV > 2) {
    std::cout << "error the AssembleAllPenalty works only for Lagrangian Solutions\n";
    abort();
  }

  unsigned indexSolC = mlSol->GetIndex("C0");

  unsigned solTypeX = 2;
  start_time = clock();

  std::vector < std::vector < std::vector <double > > > aPi(3);
  std::vector < std::vector < std::vector <double > > > aPj(3);

  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielType = msh->GetElementType(iel);
    double Ci = (*mysolution->_Sol[indexSolC])(iel);

    unsigned nDofsVi = msh->GetElementDofNumber(iel, solTypeV);
    unsigned nDofsPi = msh->GetElementDofNumber(iel, solTypeP);

    for(unsigned K = 0; K < dim; K++) solVi[K].resize(nDofsVi);
    for(unsigned K = 0; K < 2; K++) solPi[K].resize(nDofsPi);
    sysDofsi.resize(dim * nDofsVi + 2 * nDofsPi);

    for(unsigned i = 0; i < nDofsVi; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, solTypeV);
      for(unsigned K = 0; K < dim; K++) {
        solVi[K][i] = (*mysolution->_Sol[indexSolV[K]])(idof);
        sysDofsi[K * nDofsVi + i] = myLinEqSolver->GetSystemDof(indexSolV[K], indexPdeV[K], i, iel);
      }
    }

    for(unsigned i = 0; i < nDofsPi; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, solTypeP);
      for(unsigned K = 0; K < 2; K++) {
        solPi[K][i] = (*mysolution->_Sol[indexSolP[K]])(idof);
        sysDofsi[dim * nDofsVi + K * nDofsPi + i] = myLinEqSolver->GetSystemDof(indexSolP[K], indexPdeP[K], i, iel);
      }
    }

    unsigned nDofsXi = msh->GetElementDofNumber(iel, solTypeX);    // number of solution element dofs
    for(unsigned  K = 0; K < dim; K++) xi[K].resize(nDofsXi);
    for(unsigned i = 0; i < nDofsXi; i++) {
      unsigned idofX = msh->GetSolutionDof(i, iel, solTypeX);
      for(unsigned  K = 0; K < dim; K++) {
        xi[K][i] = (*msh->_topology->_Sol[K])(idofX);
      }
    }
    double h11 = (xi[0][2] - xi[0][0]);
    double h12 = (xi[1][2] - xi[1][0]);


    for(unsigned ktype = 0; ktype < solTypeP + 1; ktype++) {
      ProjectNodalToPolynomialCoefficients(aPi[ktype], xi, ielType, ktype);
    }

    for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
      int jel = el->GetFaceElementIndex(iel, iface) - 1;

      if(jel >= 0) { // iface is not a boundary of the domain
        unsigned jproc = msh->IsdomBisectionSearch(jel, 3);
        if(jproc == iproc) {

          if(jel > iel) {

            short unsigned jelType = msh->GetElementType(jel);
            double Cj = (*mysolution->_Sol[indexSolC])(jel);

            bool ghostPenalty = ((Ci > 0 && Ci < 1) || (Cj > 0 && Cj < 1)) ? true : false;

            unsigned nDofsVj = msh->GetElementDofNumber(jel, solTypeV);
            unsigned nDofsPj = msh->GetElementDofNumber(jel, solTypeP);

            for(unsigned K = 0; K < dim; K++) solVj[K].resize(nDofsVj);
            for(unsigned K = 0; K < 2; K++) solPj[K].resize(nDofsPj);
            sysDofsj.resize(dim * nDofsVj + 2 * nDofsPj);

            for(unsigned j = 0; j < nDofsVj; j++) {
              unsigned jdof = msh->GetSolutionDof(j, jel, solTypeV);
              for(unsigned K = 0; K < dim; K++) {
                solVj[K][j] = (*mysolution->_Sol[indexSolV[K]])(jdof);
                sysDofsj[K * nDofsVj + j] = myLinEqSolver->GetSystemDof(indexSolV[K], indexPdeV[K], j, jel);
              }
            }

            for(unsigned j = 0; j < nDofsPj; j++) {
              unsigned jdof = msh->GetSolutionDof(j, jel, solTypeP);
              for(unsigned K = 0; K < 2; K++) {
                solPj[K][j] = (*mysolution->_Sol[indexSolP[K]])(jdof);
                sysDofsj[dim * nDofsVj + K * nDofsPj + j] = myLinEqSolver->GetSystemDof(indexSolP[K], indexPdeP[K], j, jel);
              }
            }

            unsigned nDofsXj = msh->GetElementDofNumber(jel, solTypeX);    // number of solution element dofs
            for(unsigned  k = 0; k < dim; k++) {
              xj[k].resize(nDofsXj);
            }
            for(unsigned j = 0; j < nDofsXj; j++) {
              unsigned jdofX = msh->GetSolutionDof(j, jel, solTypeX);
              for(unsigned  K = 0; K < dim; K++) xj[K][j] = (*msh->_topology->_Sol[K])(jdofX);
            }
            double h21 = (xj[0][2] - xj[0][0]);
            double h22 = (xj[1][2] - xj[1][0]);

            aResi.assign(dim * nDofsVi + 2 * nDofsPi, 0);
            aResj.assign(dim * nDofsVj + 2 * nDofsPj, 0);

            s.new_recording();

            const unsigned faceGeom = msh->GetElementFaceType(iel, iface);
            unsigned faceDofs = msh->GetElementFaceDofNumber(iel, iface, solTypeX);
            for(unsigned  K = 0; K < dim; K++) faceVx[K].resize(faceDofs);
            for(unsigned i = 0; i < faceDofs; i++) {
              unsigned inode = msh->GetLocalFaceVertexIndex(iel, iface, i);    // face-to-element local node mapping.
              for(unsigned K = 0; K < dim; K++) {
                faceVx[K][i] =  xi[K][inode]; // We extract the local coordinates on the face from local coordinates on the element.
              }
            }

            for(unsigned ktype = 0; ktype < solTypeP + 1; ktype++) {
              ProjectNodalToPolynomialCoefficients(aPj[ktype], xj, jelType, ktype);
            }

            for(unsigned ig = 0; ig  <  msh->_finiteElement[faceGeom][solTypeP]->GetGaussPointNumber(); ig++) {

              std::vector < double> normal;
              msh->_finiteElement[faceGeom][solTypeP]->JacobianSur(faceVx, ig, weight, phi, gradPhi, normal);

              double h = 0.5 * (fabs(h11 * normal[0] + h12 * normal[1]) + fabs(h21 * normal[0] + h22 * normal[1])); //characteristic lenght in normal direction
              double h2 = h * h;
              double h3 = h2 * h;

              std::vector< double > xg(dim, 0.); // physical coordinates of the face Gauss point
              for(unsigned i = 0; i < faceDofs; i++) {
                for(unsigned K = 0; K < dim; K++) {
                  xg[K] += phi[i] * faceVx[K][i];
                }
              }

              std::vector <double> etai;//local coordinates of the face gauss point with respect to iel
              GetClosestPointInReferenceElement(xi, xg, ielType, etai);

              bool inverseMapping = GetInverseMapping(solTypeP, ielType, aPi, xg, etai, 100);
              if(!inverseMapping) {
                std::cout << "InverseMapping1 failed at " << iel << " " << jel << " " << iface << std::endl;
              }

              std::vector <double> etaj;//local coordinates of the face gauss point with respect to jel
              GetClosestPointInReferenceElement(xj, xg, jelType, etaj);

              inverseMapping = GetInverseMapping(solTypeP, jelType, aPj, xg, etaj, 100);
              if(!inverseMapping) {
                std::cout << "InverseMapping2 failed at " << iel << " " << jel << " " << iface << std::endl;
              }

              msh->_finiteElement[ielType][solTypeP]->Jacobian(xi, etai, weightPi, phiPi, gradPhiPi);
              msh->_finiteElement[jelType][solTypeP]->Jacobian(xj, etaj, weightPj, phiPj, gradPhiPj);

              msh->_finiteElement[ielType][solTypeV]->Jacobian(xi, etai, weightVi, phiVi, gradPhiVi, hessPhiVi);
              msh->_finiteElement[jelType][solTypeV]->Jacobian(xj, etaj, weightVj, phiVj, gradPhiVj, hessPhiVj);

              double C2 = 1;
              double PHIT1 = mu1 + C2 * rho1 * h2 / dt;
              double PHIT2 = mu2 + C2 * rho2 * h2 / dt;

              double PHITi = Ci * PHIT1 + (1. - Ci) * PHIT2;
              double PHITj = Cj * PHIT1 + (1. - Cj) * PHIT2;
              double PHIT = 0.5 * (PHITi + PHITj);

              double PHIP[2] = {0.05 * h2 / PHIT1, 0.05 * h2 / PHIT2};

              for(unsigned  K = 0; K < dim; K++) gradSolVi[K].assign(dim, 0.);
              for(unsigned  K = 0; K < dim; K++) gradSolVj[K].assign(dim, 0.);

              for(unsigned J = 0; J < dim; J++) {
                for(unsigned K = 0; K < dim; K++) {
                  for(unsigned i = 0; i < nDofsVi; i++) {
                    gradSolVi[J][K] += solVi[J][i] * gradPhiVi[i * dim + K];
                  }
                  for(unsigned i = 0; i < nDofsVj; i++) {
                    gradSolVj[J][K] += solVj[J][i] * gradPhiVj[i * dim + K];
                  }
                }
              }
              for(unsigned J = 0; J < dim; J++) {
                for(unsigned K = 0; K < dim; K++) {
                  for(unsigned i = 0; i < nDofsVi; i++) {
                    aResi[J * nDofsVi + i] +=  PHIT * h * (gradSolVi[K][K] -  gradSolVj[K][K])  * gradPhiVi[i * dim + J] * weight;
                  }
                  for(unsigned i = 0; i < nDofsVj; i++) {
                    aResj[J * nDofsVj + i] -=  PHIT * h * (gradSolVi[K][K] -  gradSolVj[K][K])  * gradPhiVj[i * dim + J] * weight;
                  }
                }
              }

              for(unsigned  K = 0; K < 2; K++) gradSolPi[K].assign(dim, 0.);
              for(unsigned  K = 0; K < 2; K++) gradSolPj[K].assign(dim, 0.);

              for(unsigned J = 0; J < 2; J++) {
                for(unsigned K = 0; K < dim; K++) {
                  for(unsigned i = 0; i < nDofsPi; i++) {
                    gradSolPi[J][K] += solPi[J][i] * gradPhiPi[i * dim + K];
                  }
                  for(unsigned i = 0; i < nDofsPj; i++) {
                    gradSolPj[J][K] += solPj[J][i] * gradPhiPj[i * dim + K];
                  }
                }
              }
              for(unsigned J = 0; J < 2; J++) {
                for(unsigned K = 0; K < dim; K++) {
                  for(unsigned i = 0; i < nDofsPi; i++) {
                    aResi[dim * nDofsVi + J * nDofsPi + i] +=  PHIP[J] * h * (gradSolPi[J][K] -  gradSolPj[J][K])  * gradPhiPi[i * dim + K] * weight;
                  }
                  for(unsigned i = 0; i < nDofsPj; i++) {
                    aResj[dim * nDofsVj + J * nDofsPj + i] -=  PHIP[J] * h * (gradSolPi[J][K] -  gradSolPj[J][K])  * gradPhiPj[i * dim + K] * weight;
                  }
                }
              }


              if(ghostPenalty) {

                std::vector<adept::adouble> gradSolVdotNi(dim, 0.);
                std::vector<adept::adouble> gradSolVdotNj(dim, 0.);





                for(unsigned J = 0; J < dim; J++) {
                  for(unsigned K = 0; K < dim; K++) {
                    gradSolVdotNi[J] = gradSolVi[J][K] * normal[K];
                    gradSolVdotNj[J] = gradSolVj[J][K] * normal[K];
                  }
                }

                for(unsigned J = 0; J < dim; J++) {
                  for(unsigned K = 0; K < dim; K++) {
                    for(unsigned i = 0; i < nDofsVi; i++) {
                      aResi[J * nDofsVi + i] +=  PHIT * h * (gradSolVdotNi[J] - gradSolVdotNj[J])  * gradPhiVi[i * dim + K] * normal[K] * weight;
                    }
                    for(unsigned i = 0; i < nDofsVj; i++) {
                      aResj[J * nDofsVj + i] -=  PHIT * h * (gradSolVdotNi[J] - gradSolVdotNj[J])  * gradPhiVj[i * dim + K] * normal[K] * weight;
                    }
                  }
                }

                std::vector<adept::adouble> gradSolPdotNi(2, 0.);
                std::vector<adept::adouble> gradSolPdotNj(2, 0.);

                for(unsigned J = 0; J < 2; J++) {
                  for(unsigned K = 0; K < dim; K++) {
                    gradSolPdotNi[J] = gradSolPi[J][K] * normal[K];
                    gradSolPdotNj[J] = gradSolPj[J][K] * normal[K];
                  }
                }

                for(unsigned J = 0; J < 2; J++) {
                  for(unsigned K = 0; K < dim; K++) {
                    for(unsigned i = 0; i < nDofsPi; i++) {
                      aResi[dim * nDofsVi + J * nDofsPi + i] +=  PHIP[J] * h * (gradSolPdotNi[J] - gradSolPdotNj[J]) * gradPhiPi[i * dim + K] * normal[K] * weight;
                    }
                    for(unsigned i = 0; i < nDofsPj; i++) {
                      aResj[dim * nDofsVj + J * nDofsPj + i] -=  PHIP[J] * h * (gradSolPdotNi[J] - gradSolPdotNj[J]) * gradPhiPj[i * dim + K] * normal[K] * weight;
                    }
                  }
                }


                std::vector<std::vector<std::vector<adept::adouble>>> HessSolVi(dim, std::vector<std::vector<adept::adouble>>(dim, std::vector<adept::adouble>(dim, 0.)));
                std::vector<std::vector<std::vector<adept::adouble>>> HessSolVj(dim, std::vector<std::vector<adept::adouble>>(dim, std::vector<adept::adouble>(dim, 0.)));

                for(unsigned I = 0; I < dim; I++) {
                  unsigned J = I;
                  //for(unsigned J = 0; J < dim; J++) {
                  for(unsigned K = 0; K < dim; K++) {
                    unsigned L;
                    if(J == K) L = J;
                    else if(1 == J + K) L = dim;     // xy
                    else if(2 == J + K) L = dim + 2; // xz
                    else if(3 == J + K) L = dim + 1; // yz

                    for(unsigned i = 0; i < nDofsVi; i++) {
                      HessSolVi[I][J][K] += solVi[I][i] * hessPhiVi[i * dim2 + L];
                    }
                    for(unsigned i = 0; i < nDofsVj; i++) {
                      HessSolVj[I][J][K] += solVj[I][i] * hessPhiVj[i * dim2 + L];
                    }
                  }
                  //}
                }

                std::vector<std::vector<adept::adouble>> HessSolVdotNi(dim, std::vector<adept::adouble>(dim, 0.));
                std::vector<std::vector<adept::adouble>> HessSolVdotNj(dim, std::vector<adept::adouble>(dim, 0.));

                for(unsigned I = 0; I < dim; I++) {
                  unsigned J = I;
                  //for(unsigned J = 0; J < dim; J++) {
                  for(unsigned K = 0; K < dim; K++) {
                    HessSolVdotNi[I][J] += HessSolVi[I][J][K] * normal[K];
                    HessSolVdotNj[I][J] += HessSolVj[I][J][K] * normal[K];
                  }
                  //}
                }

                for(unsigned I = 0; I < dim; I++) {
                  unsigned J = I;
                  //for(unsigned J = 0; J < dim; J++) {
                  for(unsigned K = 0; K < dim; K++) {

                    unsigned L;
                    if(J == K) L = J;
                    else if(1 == J + K) L = dim;     // xy
                    else if(2 == J + K) L = dim + 2; // xz
                    else if(3 == J + K) L = dim + 1; // yz

                    for(unsigned i = 0; i < nDofsVi; i++) {
                      aResi[I * nDofsVi + i] +=  PHIT * h3 * (HessSolVdotNi[I][J] - HessSolVdotNj[I][J]) * hessPhiVi[i * dim2 + L] * normal[K] * weight;
                    }

                    for(unsigned i = 0; i < nDofsVj; i++) {
                      aResi[I * nDofsVj + i] -=  PHIT * h3 * (HessSolVdotNi[I][J] - HessSolVdotNj[I][J])  * hessPhiVj[i * dim2 + L] * normal[K] * weight;
                    }
                  }
                }
              }
            }

            rhsi.resize(aResi.size());   //resize
            for(int i = 0; i < aResi.size(); i++) {
              rhsi[i] = -aResi[i].value();
            }
            myRES->add_vector_blocked(rhsi, sysDofsi);

            rhsj.resize(aResj.size());   //resize
            for(int i = 0; i < aResj.size(); i++) {
              rhsj[i] = -aResj[i].value();
            }
            myRES->add_vector_blocked(rhsj, sysDofsj);

            s.dependent(&aResi[0], aResi.size());
            for(unsigned K = 0; K < dim; K++) s.independent(solVi[K].data(), solVi[K].size());
            for(unsigned K = 0; K < 2; K++) s.independent(solPi[K].data(), solPi[K].size());

            Jac.resize(sysDofsi.size() * sysDofsi.size()); //J11
            s.jacobian(&Jac[0], true);
            myKK->add_matrix_blocked(Jac, sysDofsi, sysDofsi);

            s.clear_independents();
            for(unsigned K = 0; K < dim; K++) s.independent(solVj[K].data(), solVj[K].size());
            for(unsigned K = 0; K < 2; K++) s.independent(solPj[K].data(), solPj[K].size());

            Jac.resize(sysDofsi.size() * sysDofsj.size()); //J12
            s.jacobian(&Jac[0], true);
            myKK->add_matrix_blocked(Jac, sysDofsi, sysDofsj);

            s.clear_dependents();
            s.clear_independents();
            s.dependent(&aResj[0], aResj.size());
            for(unsigned K = 0; K < dim; K++) s.independent(solVi[K].data(), solVi[K].size());
            for(unsigned K = 0; K < 2; K++) s.independent(solPi[K].data(), solPi[K].size());

            Jac.resize(sysDofsj.size() * sysDofsi.size()); //J21
            s.jacobian(&Jac[0], true);
            myKK->add_matrix_blocked(Jac, sysDofsj, sysDofsi);

            s.clear_independents();
            for(unsigned K = 0; K < dim; K++) s.independent(solVj[K].data(), solVj[K].size());
            for(unsigned K = 0; K < 2; K++) s.independent(solPj[K].data(), solPj[K].size());

            Jac.resize(sysDofsj.size() * sysDofsj.size()); //J22
            s.jacobian(&Jac[0], true);
            myKK->add_matrix_blocked(Jac, sysDofsj, sysDofsj);
            s.clear_independents();
          }
        }
      }
    }
  }
  std::cout << "All Penalty Assembly time = " << static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC << std::endl;

}





















