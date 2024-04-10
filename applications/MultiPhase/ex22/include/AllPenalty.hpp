

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

  std::vector<adept::adouble> solPi;
  std::vector<adept::adouble> solPj;

  std::vector<adept::adouble> aResPi;
  std::vector<adept::adouble> aResPj;

  std::vector<double> rhsPi;
  std::vector<double> rhsPj;

   vector < double > Jac;

  std::vector <unsigned> sysDofsPi;
  std::vector <unsigned> sysDofsPj;

  double weight;
  std::vector < double > phi;
  std::vector < double> gradPhi;

  double weightPi;
  std::vector<double> phiPi, gradPhiPi;

  double weightPj;
  std::vector<double> phiPj, gradPhiPj;


  vector <vector < double> > xi(dim);
  vector <vector < double> > xj(dim);

  std::cout.precision(10);

  //variable-name handling
  unsigned indexSolP1 = mlSol->GetIndex("P1");
  unsigned indexSolP2 = mlSol->GetIndex("P2");
  unsigned indexPdeP1 = my_nnlin_impl_sys.GetSolPdeIndex("P1");
  unsigned indexPdeP2 = my_nnlin_impl_sys.GetSolPdeIndex("P2");
  unsigned solTypeP = mlSol->GetSolutionType("P1");

  if(solTypeP > 2) {
    std::cout << "error this function works only for Lagrangian Solutions\n";
    abort();
  }

  unsigned solTypeX = 2;
  start_time = clock();

  std::vector < std::vector < std::vector <double > > > aPi(3);
  std::vector < std::vector < std::vector <double > > > aPj(3);

  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielType = msh->GetElementType(iel);

    unsigned nDofsPi = msh->GetElementDofNumber(iel, solTypeP);
    solPi.resize(2 * nDofsPi);
    sysDofsPi.resize(2 * nDofsPi);

    for(unsigned i = 0; i < nDofsPi; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, solTypeP);
      solPi[i] = (*mysolution->_Sol[indexSolP1])(idof);
      solPi[i + nDofsPi] = (*mysolution->_Sol[indexSolP2])(idof);
      sysDofsPi[i] = myLinEqSolver->GetSystemDof(indexSolP1, indexPdeP1, i, iel);
      sysDofsPi[i + nDofsPi] = myLinEqSolver->GetSystemDof(indexSolP2, indexPdeP2, i, iel);
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


    for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
      int jel = el->GetFaceElementIndex(iel, iface) - 1;

      if(jel >= 0) { // iface is not a boundary of the domain
        unsigned jproc = msh->IsdomBisectionSearch(jel, 3);
        if(jproc == iproc) {

          if(jel > iel) {

            short unsigned jelType = msh->GetElementType(jel);

            unsigned nDofsPj = msh->GetElementDofNumber(jel, solTypeP);
            solPj.resize(2 * nDofsPj);
            sysDofsPj.resize(2 * nDofsPj);

            for(unsigned j = 0; j < nDofsPj; j++) {
              unsigned jdof = msh->GetSolutionDof(j, jel, solTypeP);
              solPj[j] = (*mysolution->_Sol[indexSolP1])(jdof);
              solPj[j + nDofsPj] = (*mysolution->_Sol[indexSolP2])(jdof);
              sysDofsPj[j] = myLinEqSolver->GetSystemDof(indexSolP1, indexPdeP1, j, jel);
              sysDofsPj[j + nDofsPj] = myLinEqSolver->GetSystemDof(indexSolP2, indexPdeP2, j, jel);
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

            aResPi.assign(2 * nDofsPi, 0);
            aResPj.assign(2 * nDofsPj, 0);

            s.new_recording();

            const unsigned faceGeom = msh->GetElementFaceType(iel, iface);
            unsigned faceDofs = msh->GetElementFaceDofNumber(iel, iface, solTypeP);
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
              ProjectNodalToPolynomialCoefficients(aPi[ktype], xi, ielType, ktype);
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

              double C2 = 1;
              double PHIT1 = mu1 + C2 * rho1 *h2 / dt;
              double PHI1 = 0.05 * h2 /PHIT1;

              double PHIT2 = mu2 + C2 * rho2 *h2 / dt;
              double PHI2 = 0.05 * h2 /PHIT2;

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
                  aResPi[i] +=  PHI1 * h * (gradSolP1i[k] -  gradSolP1j[k])  * gradPhiPi[i * dim + k] * weightPi;
                  aResPi[i + nDofsPi] +=  PHI2 * h * (gradSolP2i[k] -  gradSolP2j[k])  * gradPhiPi[i * dim + k] * weightPi;
                }
                for(unsigned i = 0; i < nDofsPj; i++) {
                  aResPj[i] -=  PHI1 * h * (gradSolP1i[k] -  gradSolP1j[k])  * gradPhiPj[i * dim + k] * weightPj;
                  aResPj[i + nDofsPj] -=  PHI2 * h * (gradSolP2i[k] -  gradSolP2j[k])  * gradPhiPj[i * dim + k] * weightPj;
                }
              }
            }


            rhsPi.resize(sysDofsPi.size());   //resize
            for(int i = 0; i < sysDofsPi.size(); i++) {
              rhsPi[i] = -aResPi[i].value();
            }
            myRES->add_vector_blocked(rhsPi, sysDofsPi);

            rhsPj.resize(sysDofsPj.size());   //resize
            for(int i = 0; i < sysDofsPj.size(); i++) {
              rhsPj[i] = -aResPj[i].value();
            }
            myRES->add_vector_blocked(rhsPj, sysDofsPj);




            s.dependent(&aResPi[0], aResPi.size());
            s.independent(&solPi[0], solPi.size());

            Jac.resize(sysDofsPi.size() * sysDofsPi.size()); //J11
            s.jacobian(&Jac[0], true);
            myKK->add_matrix_blocked(Jac, sysDofsPi, sysDofsPi);

            s.clear_independents();
            s.independent(&solPj[0], solPj.size());

            Jac.resize( sysDofsPi.size() * sysDofsPj.size());//J12
            s.jacobian(&Jac[0], true);
            myKK->add_matrix_blocked(Jac, sysDofsPi, sysDofsPj);

            s.clear_dependents();
            s.clear_independents();
            s.dependent(&aResPj[0], aResPj.size());
            s.independent(&solPi[0], solPi.size());

            Jac.resize(sysDofsPj.size() * sysDofsPi.size()); //J21
            s.jacobian(&Jac[0], true);
            myKK->add_matrix_blocked(Jac, sysDofsPj, sysDofsPi);

            s.clear_independents();
            s.independent(&solPj[0], solPj.size());

            Jac.resize(sysDofsPj.size() * sysDofsPj.size()); //J22
            s.jacobian(&Jac[0], true);
            myKK->add_matrix_blocked(Jac, sysDofsPj, sysDofsPj);
            s.clear_independents();
          }
        }
      }
    }
  }

  std::cout << "CIP  Pressure Assembly time = " << static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC << std::endl;
}




















