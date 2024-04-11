

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

  // data
  unsigned iproc  = msh->processor_id();
  unsigned nprocs  = msh->n_processors();

  //quantities for iel will have index1
  //quantities for jel will have index2


  std::vector<adept::adouble> solVi;
  std::vector<adept::adouble> solVj;

  std::vector<adept::adouble> solPi;
  std::vector<adept::adouble> solPj;

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
  std::vector<double> phiVi, gradPhiVi;

  double weightVj;
  std::vector<double> phiVj, gradPhiVj;


  vector <vector < double> > xi(dim);
  vector <vector < double> > xj(dim);

  std::cout.precision(10);

  //variable-name handling

  std::vector<unsigned> indexSolV(dim);
  std::vector<unsigned> indexPdeV(dim);
  indexSolV[0] = mlSol->GetIndex("U");
  indexSolV[1] = mlSol->GetIndex("V");
  indexPdeV[0] = my_nnlin_impl_sys.GetSolPdeIndex("U");
  indexPdeV[1] = my_nnlin_impl_sys.GetSolPdeIndex("V");
  unsigned solTypeV = mlSol->GetSolutionType("U");

  unsigned indexSolP1 = mlSol->GetIndex("P1");
  unsigned indexSolP2 = mlSol->GetIndex("P2");
  unsigned indexPdeP1 = my_nnlin_impl_sys.GetSolPdeIndex("P1");
  unsigned indexPdeP2 = my_nnlin_impl_sys.GetSolPdeIndex("P2");
  unsigned solTypeP = mlSol->GetSolutionType("P1");

  if(solTypeP > 2 || solTypeV > 2) {
    std::cout << "error the AssembleAllPenalty works only for Lagrangian Solutions\n";
    abort();
  }

  unsigned solTypeX = 2;
  start_time = clock();

  std::vector < std::vector < std::vector <double > > > aPi(3);
  std::vector < std::vector < std::vector <double > > > aPj(3);

  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielType = msh->GetElementType(iel);

    unsigned nDofsVi = msh->GetElementDofNumber(iel, solTypeV);
    unsigned nDofsPi = msh->GetElementDofNumber(iel, solTypeP);

    solVi.resize(dim * nDofsVi);
    solPi.resize(2 * nDofsPi);
    sysDofsi.resize(dim * nDofsVi + 2 * nDofsPi);

    for(unsigned i = 0; i < nDofsVi; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, solTypeV);
      for(unsigned k = 0; k < dim; k++) {
        solVi[i + k * nDofsVi] = (*mysolution->_Sol[indexSolV[k]])(idof);
        sysDofsi[i + k * nDofsVi] = myLinEqSolver->GetSystemDof(indexSolV[k], indexPdeV[k], i, iel);
      }
    }

    for(unsigned i = 0; i < nDofsPi; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, solTypeP);
      solPi[i] = (*mysolution->_Sol[indexSolP1])(idof);
      solPi[i + nDofsPi] = (*mysolution->_Sol[indexSolP2])(idof);
      sysDofsi[i + dim * nDofsVi] = myLinEqSolver->GetSystemDof(indexSolP1, indexPdeP1, i, iel);
      sysDofsi[i + dim * nDofsVi + nDofsPi] = myLinEqSolver->GetSystemDof(indexSolP2, indexPdeP2, i, iel);
    }

    unsigned nDofsXi = msh->GetElementDofNumber(iel, solTypeX);    // number of solution element dofs
    for(unsigned  k = 0; k < dim; k++) {
      xi[k].resize(nDofsXi);
    }
    for(unsigned i = 0; i < nDofsXi; i++) {
      unsigned idofX = msh->GetSolutionDof(i, iel, solTypeX);
      for(unsigned  k = 0; k < dim; k++) {
        xi[k][i] = (*msh->_topology->_Sol[k])(idofX);
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

            unsigned nDofsVj = msh->GetElementDofNumber(jel, solTypeV);
            unsigned nDofsPj = msh->GetElementDofNumber(jel, solTypeP);

            solVj.resize(dim * nDofsVj);
            solPj.resize(2 * nDofsPj);
            sysDofsj.resize(dim * nDofsVj + 2 * nDofsPj);

            for(unsigned j = 0; j < nDofsVj; j++) {
              unsigned jdof = msh->GetSolutionDof(j, jel, solTypeV);
              for(unsigned k = 0; k < dim; k++) {
                solVj[j + k * nDofsVj] = (*mysolution->_Sol[indexSolV[k]])(jdof);
                sysDofsj[j + k * nDofsVj] = myLinEqSolver->GetSystemDof(indexSolV[k], indexPdeV[k], j, jel);
              }
            }

            for(unsigned j = 0; j < nDofsPj; j++) {
              unsigned jdof = msh->GetSolutionDof(j, jel, solTypeP);
              solPj[j] = (*mysolution->_Sol[indexSolP1])(jdof);
              solPj[j + nDofsPj] = (*mysolution->_Sol[indexSolP2])(jdof);
              sysDofsj[j + dim * nDofsVj] = myLinEqSolver->GetSystemDof(indexSolP1, indexPdeP1, j, jel);
              sysDofsj[j + dim * nDofsVj + nDofsPj] = myLinEqSolver->GetSystemDof(indexSolP2, indexPdeP2, j, jel);
            }

            unsigned nDofsXj = msh->GetElementDofNumber(jel, solTypeX);    // number of solution element dofs
            for(unsigned  k = 0; k < dim; k++) {
              xj[k].resize(nDofsXj);
            }
            for(unsigned j = 0; j < nDofsXj; j++) {
              unsigned jdofX = msh->GetSolutionDof(j, jel, solTypeX);
              for(unsigned  k = 0; k < dim; k++) {
                xj[k][j] = (*msh->_topology->_Sol[k])(jdofX);
              }
            }
            double h21 = (xj[0][2] - xj[0][0]);
            double h22 = (xj[1][2] - xj[1][0]);

            aResi.assign(dim * nDofsVi + 2 * nDofsPi, 0);
            aResj.assign(dim * nDofsVj + 2 * nDofsPj, 0);

            s.new_recording();

            const unsigned faceGeom = msh->GetElementFaceType(iel, iface);
            unsigned faceDofs = msh->GetElementFaceDofNumber(iel, iface, solTypeX);
            std::vector  < std::vector  <  double > > faceVx(dim);    // A matrix holding the face coordinates rowwise.
            for(int k = 0; k < dim; k++) {
              faceVx[k].resize(faceDofs);
            }
            for(unsigned i = 0; i < faceDofs; i++) {
              unsigned inode = msh->GetLocalFaceVertexIndex(iel, iface, i);    // face-to-element local node mapping.
              for(unsigned k = 0; k < dim; k++) {
                faceVx[k][i] =  xi[k][inode]; // We extract the local coordinates on the face from local coordinates on the element.
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

              std::vector< double > xg(dim, 0.); // physical coordinates of the face Gauss point
              for(unsigned i = 0; i < faceDofs; i++) {
                for(unsigned k = 0; k < dim; k++) {
                  xg[k] += phi[i] * faceVx[k][i];
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

              msh->_finiteElement[ielType][solTypeV]->Jacobian(xi, etai, weightVi, phiVi, gradPhiVi);
              msh->_finiteElement[jelType][solTypeV]->Jacobian(xj, etaj, weightVj, phiVj, gradPhiVj);

              double C2 = 1;
              double PHIT1 = mu1 + C2 * rho1 * h2 / dt;
              double PHI1 = 0.05 * h2 / PHIT1;

              double PHIT2 = mu2 + C2 * rho2 * h2 / dt;
              double PHI2 = 0.05 * h2 / PHIT2;



              std::vector < std::vector < adept::adouble > > gradSolVi(dim, std::vector < adept::adouble >(dim, 0.));
              std::vector < std::vector < adept::adouble > > gradSolVj(dim, std::vector < adept::adouble >(dim, 0.));

              for(unsigned j = 0; j < dim; j++) {
                for(unsigned k = 0; k < dim; k++) {
                  for(unsigned i = 0; i < nDofsVi; i++) {
                    gradSolVi[j][k] += solVi[i + j * nDofsVi] * gradPhiVi[i * dim + k];
                  }
                  for(unsigned i = 0; i < nDofsVj; i++) {
                    gradSolVj[j][k] += solVj[i + j * nDofsVj] * gradPhiVj[i * dim + k];
                  }
                }
              }
              for(unsigned j = 0; j < dim; j++) {
                for(unsigned k = 0; k < dim; k++) {
                  for(unsigned i = 0; i < nDofsVi; i++) {
                    aResi[i + j * nDofsVi] +=  PHIT1 * h * (gradSolVi[k][k] -  gradSolVj[k][k])  * gradPhiVi[i * dim + j] * weightVi;

                  }
                  for(unsigned i = 0; i < nDofsVj; i++) {
                    aResj[i + j * nDofsVj] -=  PHIT2 * h * (gradSolVi[k][k] -  gradSolVj[k][k])  * gradPhiVj[i * dim + j] * weightVj;

                  }
                }
              }



              std::vector < adept::adouble > gradSolP1i(dim, 0.);
              std::vector < adept::adouble > gradSolP2i(dim, 0.);

              std::vector < adept::adouble > gradSolP1j(dim, 0.);
              std::vector < adept::adouble > gradSolP2j(dim, 0.);

              for(unsigned k = 0; k < dim; k++) {
                for(unsigned i = 0; i < nDofsPi; i++) {
                  gradSolP1i[k] += solPi[i] * gradPhiPi[i * dim + k];
                  gradSolP2i[k] += solPi[i + nDofsPi] * gradPhiPi[i * dim + k];
                }
                for(unsigned i = 0; i < nDofsPj; i++) {
                  gradSolP1j[k] += solPj[i] * gradPhiPj[i * dim + k];
                  gradSolP2j[k] += solPj[i + nDofsPj] * gradPhiPj[i * dim + k];
                }
              }

              for(unsigned k = 0; k < dim; k++) {
                for(unsigned i = 0; i < nDofsPi; i++) {
                  aResi[i + dim * nDofsVi] +=  PHI1 * h * (gradSolP1i[k] -  gradSolP1j[k])  * gradPhiPi[i * dim + k] * weightPi;
                  aResi[i + dim * nDofsVi + nDofsPi] +=  PHI1 * h * (gradSolP2i[k] -  gradSolP2j[k])  * gradPhiPi[i * dim + k] * weightPi;
                }
                for(unsigned i = 0; i < nDofsPj; i++) {
                  aResj[i + dim * nDofsVj] -=  PHI1 * h * (gradSolP1i[k] -  gradSolP1j[k])  * gradPhiPj[i * dim + k] * weightPj;
                  aResj[i + dim * nDofsVj + nDofsPj] -=  PHI2 * h * (gradSolP2i[k] -  gradSolP2j[k])  * gradPhiPj[i * dim + k] * weightPj;
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
            s.independent(&solVi[0], solVi.size());
            s.independent(&solPi[0], solPi.size());

            Jac.resize(sysDofsi.size() * sysDofsi.size()); //J11
            s.jacobian(&Jac[0], true);
            myKK->add_matrix_blocked(Jac, sysDofsi, sysDofsi);

            s.clear_independents();
            s.independent(&solVj[0], solVj.size());
            s.independent(&solPj[0], solPj.size());

            Jac.resize(sysDofsi.size() * sysDofsj.size()); //J12
            s.jacobian(&Jac[0], true);
            myKK->add_matrix_blocked(Jac, sysDofsi, sysDofsj);

            s.clear_dependents();
            s.clear_independents();
            s.dependent(&aResj[0], aResj.size());
            s.independent(&solVi[0], solVi.size());
            s.independent(&solPi[0], solPi.size());

            Jac.resize(sysDofsj.size() * sysDofsi.size()); //J21
            s.jacobian(&Jac[0], true);
            myKK->add_matrix_blocked(Jac, sysDofsj, sysDofsi);

            s.clear_independents();
            s.independent(&solVj[0], solVj.size());
            s.independent(&solPj[0], solPj.size());

            Jac.resize(sysDofsj.size() * sysDofsj.size()); //J22
            s.jacobian(&Jac[0], true);
            myKK->add_matrix_blocked(Jac, sysDofsj, sysDofsj);
            s.clear_independents();
          }
        }
      }
    }
  }

  std::cout << "CIP  Pressure Assembly time = " << static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC << std::endl;
}




















