
void AssembleGhostPenaltyDGP(MultiLevelProblem& ml_prob, const bool &P1) {

  //this function works both for fluid and solid ghost penalty, the boolean fluid switches between the two

  clock_t start_time;

  //pointers and references

  TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("NS");
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

  adept::adouble solP1;
  adept::adouble solP2;

  adept::adouble aRes1;     // local redidual vector
  adept::adouble aRes2;     // local redidual vector

  std::vector < double >  rhs1(1); // local redidual vector
  std::vector < double >  rhs2(1); // local redidual vector
  vector < double > Jac;

  std::vector <unsigned> sysDofsP1(1);
  std::vector <unsigned> sysDofsP2(1);

  double weight;
  std::vector < double > phi;
  std::vector < double> gradPhi;

  vector <vector < double> > vx(dim);

  double mu = 2. * mu1 * mu2 / (mu1 + mu2);
  double rho = 2. * rho1 * rho2 / (rho1 + rho2);

  double dt =  my_nnlin_impl_sys.GetIntervalTime();

  std::cout.precision(10);

  //variable-name handling

  unsigned indexSol = (P1) ? mlSol->GetIndex("P1") : mlSol->GetIndex("P2");
  unsigned indexPde = (P1) ? my_nnlin_impl_sys.GetSolPdeIndex("P1") : my_nnlin_impl_sys.GetSolPdeIndex("P2");
  //unsigned solType = mlSol->GetSolutionType("P1");

  unsigned indexSolC = mlSol->GetIndex("C");
  unsigned solTypeX = 2;

  start_time = clock();

  //flagmark
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    double Ciel = (*mysolution->_Sol[indexSolC])(iel);

    if(Ciel > 0 && Ciel < 1) {

      short unsigned ielt1 = msh->GetElementType(iel);
      unsigned nDofsX = msh->GetElementDofNumber(iel, solTypeX);    // number of solution element dofs
      for(unsigned  k = 0; k < dim; k++) {
        vx[k].resize(nDofsX);
      }

      solP1 = (*mysolution->_Sol[indexSol])(iel);
      sysDofsP1[0] = myLinEqSolver->GetSystemDof(indexSol, indexPde, 0, iel);


      for(unsigned i = 0; i < nDofsX; i++) {
        unsigned idofX = msh->GetSolutionDof(i, iel, solTypeX);
        for(unsigned  k = 0; k < dim; k++) {
          vx[k][i] = (*msh->_topology->_Sol[k])(idofX);
        }
      }



      for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
        int jel = el->GetFaceElementIndex(iel, iface) - 1;

        if(jel >= 0) { // iface is not a boundary of the domain
          unsigned jproc = msh->IsdomBisectionSearch(jel, 3);
          if(jproc == iproc) {



            double Cjel = (*mysolution->_Sol[indexSolC])(jel);

            if((Cjel > 0 && Cjel < 1 && jel > iel) || Cjel == P1) {

              short unsigned ielt2 = msh->GetElementType(jel);
              solP2 = (*mysolution->_Sol[indexSol])(jel);
              sysDofsP2[0] = myLinEqSolver->GetSystemDof(indexSol, indexPde, 0, jel);


              aRes1 = 0;
              aRes2 = 0;

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
                  faceVx[k][i] =  vx[k][inode]; // We extract the local coordinates on the face from local coordinates on the element.
                }
              }

              double h11 = (vx[0][2] - vx[0][0]);
              double h12 = (vx[1][2] - vx[1][0]);

              //double h21 = (vx2[0][2] - vx2[0][0]);
              //double h22 = (vx2[1][2] - vx2[1][0]);

              for(unsigned ig = 0; ig  <  msh->_finiteElement[faceGeom][solTypeX]->GetGaussPointNumber(); ig++) {

                std::vector < double> normal;
                msh->_finiteElement[faceGeom][solTypeX]->JacobianSur(faceVx, ig, weight, phi, gradPhi, normal);

                double h = fabs(h11 * normal[0] + h12 * normal[1]);// + fabs(h21 * normal[0] + h22 * normal[1])); //characteristic lenght in normal direction

                double C1 = 0.05 * h / mu; // h * h / (sigma * dt); dt /(rho * h);
                aRes1 +=  C1 * (solP1 - solP2) * weight;
                aRes2 -=  C1 * (solP1 - solP2) * weight;

              }

              //copy the value of the adept::adoube aRes in double Res and store them in RES
              rhs1[0] = -aRes1.value();
              myRES->add_vector_blocked(rhs1, sysDofsP1);

              rhs2[0] = -aRes2.value();
              myRES->add_vector_blocked(rhs2, sysDofsP2);

              std::vector<unsigned> sysDofsAll = {sysDofsP1[0], sysDofsP2[0]};

              // define the dependent variables J11 and J12
              s.dependent(&aRes1, 1);
              s.dependent(&aRes2, 1);
              s.independent(&solP1, 1);
              s.independent(&solP2, 1);

              Jac.resize(2 * 2);
              // get the and store jacobian matrix (row-major)
              s.jacobian(&Jac[0], true);
              myKK->add_matrix_blocked(Jac, sysDofsAll, sysDofsAll);
              s.clear_independents();
              s.clear_dependents(); // for J21 and J22
            }
          }
        }
      }
    }
  }
  // *************************************
  std::cout << "Ghost Penalty Assembly time = " << static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC << std::endl;

}
