
void AssembleGhostPenalty(MultiLevelProblem& ml_prob, const bool &fluid) {

  // ml_prob is the global object from/to where get/set all the data
  // level is the level of the PDE system to be assembled
  // levelMax is the Maximum level of the MultiLevelProblem
  // assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  clock_t start_time;

  //pointers and references

  TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("MPM_FSI");
  const unsigned  level = my_nnlin_impl_sys.GetLevelToAssemble();
  MultiLevelSolution* mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* mysolution = mlSol->GetSolutionLevel(level);     // pointer to the solution (level) object

  //Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* myLinEqSolver = my_nnlin_impl_sys._LinSolver[level];  // pointer to the equation (level) object

  Mesh* msh = ml_prob._ml_msh->GetLevel(level);     // pointer to the mesh (level) object
  elem* el = msh->el;   // pointer to the elem object in msh (level)
  SparseMatrix* myKK = myLinEqSolver->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* myRES =  myLinEqSolver->_RES;  // pointer to the global residual vector object in pdeSys (level)

  MatSetOption((static_cast< PetscMatrix* >(myKK))->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);


  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  const unsigned dim = msh->GetDimension();
  const unsigned dim2 = 3 * (dim - 1);

  // data
  unsigned iproc  = msh->processor_id();
  unsigned nprocs  = msh->n_processors();

  //quantities for iel will have index1
  //quantities for jel will have index2


  vector< vector< double > > solV1Old(dim);
  vector< vector< double > > solV2Old(dim);

  vector< vector< adept::adouble > > solV1(dim); // local solution (velocity)
  vector< adept::adouble > solP1; // local solution (velocity)

  vector< vector< adept::adouble > > solV2(dim); // local solution (velocity)
  vector< vector< double > > solV2d(dim); // local solution (velocity)
  vector< adept::adouble > solP2; // local solution (velocity)
  vector< double > solP2d; // local solution (velocity)

  vector< vector< adept::adouble > > aResV1(dim);     // local redidual vector
  vector< adept::adouble > aResP1; // local redidual vector

  vector< vector< adept::adouble > > aResV2(dim);     // local redidual vector
  vector< adept::adouble > aResP2; // local redidual vector

  vector< double > rhs1; // local redidual vector
  vector< double > rhs2; // local redidual vector
  vector < double > Jac;

  std::vector <unsigned> sysDofs1;
  std::vector <unsigned> sysDofs2;
  double weight;
  std::vector < double > phi;
  std::vector < double> gradPhi;

  double weight1;
  std::vector < double > phi1;
  std::vector < double> gradPhi1;
  std::vector < double> nablaPhi1;

  double weight2;
  std::vector < double > phi2;
  std::vector < double> gradPhi2;
  std::vector < double> nablaPhi2;


  vector <vector < double> > vx1(dim);
  vector <vector < double> > vx2(dim);

  //reading parameters for fluid FEM domain
  double rhoFluid = (fluid) ? ml_prob.parameters.get<Fluid> ("FluidFEM").get_density() :
                    ml_prob.parameters.get<Solid> ("SolidMPM").get_density();
  double muFluid = (fluid) ? ml_prob.parameters.get<Fluid> ("FluidFEM").get_viscosity() :
                   ml_prob.parameters.get<Solid> ("SolidMPM").get_lame_shear_modulus();

  double dt =  my_nnlin_impl_sys.GetIntervalTime();


  std::cout.precision(10);

  //variable-name handling
  const char varname[10][5] = {"DX", "DY", "DZ", "VX", "VY", "VZ"};

  vector <unsigned> indexSolV(dim);
  vector <unsigned> indexPdeV(dim);
  for(unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolV[ivar] = mlSol->GetIndex(&varname[ivar + 3 * fluid][0]);
    indexPdeV[ivar] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar + 3 * fluid][0]);
  }
  unsigned solTypeV = mlSol->GetSolutionType(&varname[0][0]);

  unsigned indexSolP = mlSol->GetIndex("P");
  unsigned indexPdeP = my_nnlin_impl_sys.GetSolPdeIndex("P");
  unsigned solTypeP = mlSol->GetSolutionType("P");

  unsigned eflagIndex = mlSol->GetIndex("eflag");

  start_time = clock();

  std::vector < std::vector < std::vector <double > > > aP1(3);
  std::vector < std::vector < std::vector <double > > > aP2(3);


  //flagmark
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned eFlag1 = static_cast <unsigned>(floor((*mysolution->_Sol[eflagIndex])(iel) + 0.5));
    if(eFlag1 == 1) {

      short unsigned ielt1 = msh->GetElementType(iel);

      unsigned nDofsV1 = msh->GetElementDofNumber(iel, solTypeV);    // number of solution element dofs
      unsigned nDofsP1 = msh->GetElementDofNumber(iel, solTypeP);    // number of solution element dofs
      unsigned nDofs1 = dim * nDofsV1 + nDofsP1;

      // resize local arrays
      sysDofs1.resize(nDofs1);
      //nodeFlag1.resize(nDofsV1);

      for(unsigned  k = 0; k < dim; k++) {
        solV1[k].resize(nDofsV1);
        solV1Old[k].resize(nDofsV1);
        vx1[k].resize(nDofsV1);
      }
      solP1.resize(nDofsP1);

      for(unsigned i = 0; i < nDofsV1; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, solTypeV);

        //nodeFlag1[i] = (*mysolution->_Sol[nflagIndex])(idof);

        for(unsigned  k = 0; k < dim; k++) {
          solV1[k][i] = (*mysolution->_Sol[indexSolV[k]])(idof);
          solV1Old[k][i] = (*mysolution->_SolOld[indexSolV[k]])(idof);
          sysDofs1[k * nDofsV1 + i] = myLinEqSolver->GetSystemDof(indexSolV[k], indexPdeV[k], i, iel);
        }
      }

      for(unsigned i = 0; i < nDofsP1; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, solTypeP);
        solP1[i] = (*mysolution->_Sol[indexSolP])(idof);
        sysDofs1[dim * nDofsV1 + i] = myLinEqSolver->GetSystemDof(indexSolP, indexPdeP, i, iel);
      }

      for(unsigned i = 0; i < nDofsV1; i++) {
        unsigned idofX = msh->GetSolutionDof(i, iel, 2);
        for(unsigned  k = 0; k < dim; k++) {
          vx1[k][i] = (*msh->_topology->_Sol[k])(idofX);
        }
      }

      bool aP1IsInitialized = false;

      for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
        int jel = el->GetFaceElementIndex(iel, iface) - 1;
        if(jel >= 0) { // iface is not a boundary of the domain
          unsigned jproc = msh->IsdomBisectionSearch(jel , 3);
          if(jproc == iproc) {
            unsigned eFlag2 = static_cast <unsigned>(floor((*mysolution->_Sol[eflagIndex])(jel) + 0.5));

            if(eFlag2 == 0 + !fluid * 2 || (eFlag2 == 1 && jel > iel)) {

              short unsigned ielt2 = msh->GetElementType(jel);

              unsigned nDofsV2 = msh->GetElementDofNumber(jel, solTypeV);    // number of solution element dofs
              unsigned nDofsP2 = msh->GetElementDofNumber(jel, solTypeP);    // number of solution element dofs
              unsigned nDofs2 = dim * nDofsV2 + nDofsP2;

              // resize local arrays
              sysDofs2.resize(nDofs2);
              //nodeFlag2.resize(nDofsV2);

              for(unsigned  k = 0; k < dim; k++) {
                solV2[k].resize(nDofsV2);
                solV2Old[k].resize(nDofsV2);
                vx2[k].resize(nDofsV2);
              }
              solP2.resize(nDofsP2);


              for(unsigned i = 0; i < nDofsV2; i++) {
                unsigned idof = msh->GetSolutionDof(i, jel, solTypeV);
                //nodeFlag2[i] = (*mysolution->_Sol[nflagIndex])(idof);
                for(unsigned  k = 0; k < dim; k++) {
                  solV2[k][i] = (*mysolution->_Sol[indexSolV[k]])(idof);
                  solV2Old[k][i] = (*mysolution->_SolOld[indexSolV[k]])(idof);
                  sysDofs2[k * nDofsV2 + i] = myLinEqSolver->GetSystemDof(indexSolV[k], indexPdeV[k], i, jel);
                }
              }

              for(unsigned i = 0; i < nDofsP2; i++) {
                unsigned idof = msh->GetSolutionDof(i, jel, solTypeP);
                solP2[i] = (*mysolution->_Sol[indexSolP])(idof);
                sysDofs2[dim * nDofsV2 + i] = myLinEqSolver->GetSystemDof(indexSolP, indexPdeP, i, jel);
              }

              for(unsigned i = 0; i < nDofsV2; i++) {
                unsigned idofX = msh->GetSolutionDof(i, jel, 2);
                for(unsigned  k = 0; k < dim; k++) {
                  vx2[k][i] = (*msh->_topology->_Sol[k])(idofX);
                }
              }

              for(unsigned  k = 0; k < dim; k++) {
                aResV1[k].assign(nDofsV1, 0.);
                aResV2[k].assign(nDofsV2, 0.);
              }
              aResP1.assign(nDofsP1, 0.);
              aResP2.assign(nDofsP2, 0.);

              s.new_recording();


              const unsigned faceGeom = msh->GetElementFaceType(iel, iface);
              unsigned faceDofs = msh->GetElementFaceDofNumber(iel, iface, solTypeV);
              std::vector  < std::vector  <  double > > faceVx(dim);    // A matrix holding the face coordinates rowwise.
              for(int k = 0; k < dim; k++) {
                faceVx[k].resize(faceDofs);
              }
              for(unsigned i = 0; i < faceDofs; i++) {
                unsigned inode = msh->GetLocalFaceVertexIndex(iel, iface, i);    // face-to-element local node mapping.
                for(unsigned k = 0; k < dim; k++) {
                  faceVx[k][i] =  vx1[k][inode]; // We extract the local coordinates on the face from local coordinates on the element.
                }
              }

              double h = sqrt((faceVx[0][0] - faceVx[0][1]) * (faceVx[0][0] - faceVx[0][1]) +
                              (faceVx[1][0] - faceVx[1][1]) * (faceVx[1][0] - faceVx[1][1])) ;
              double h2 = h * h;
              double h3 = h * h * h;


              if(!aP1IsInitialized) { //build the basis 1,x,y,z... corresponding to the solution type
                aP1IsInitialized = true;
                for(unsigned jtype = 0; jtype < solTypeV + 1; jtype++) {
                  ProjectNodalToPolynomialCoefficients(aP1[jtype], vx1, ielt1, jtype);
                }
              }

              for(unsigned jtype = 0; jtype < solTypeV + 1; jtype++) {
                ProjectNodalToPolynomialCoefficients(aP2[jtype], vx2, ielt2, jtype);
              }

              for(unsigned ig = 0; ig  <  msh->_finiteElement[faceGeom][solTypeV]->GetGaussPointNumber(); ig++) {

                std::vector < double> normal;
                msh->_finiteElement[faceGeom][solTypeV]->JacobianSur(faceVx, ig, weight, phi, gradPhi, normal);

                std::vector< double > xg(dim, 0.); // physical coordinates of the face Gauss point
                for(unsigned i = 0; i < faceDofs; i++) {
                  for(unsigned k = 0; k < dim; k++) {
                    xg[k] += phi[i] * faceVx[k][i];
                  }
                }

                std::vector <double> xi1;//local coordinates of the face gauss point with respect to iel
                GetClosestPointInReferenceElement(vx1, xg, ielt1, xi1);

                bool inverseMapping = GetInverseMapping(solTypeV, ielt1, aP1, xg, xi1, 100);
                if(!inverseMapping) {
                  std::cout << "InverseMapping1 failed at " << iel << " " << jel << " " << iface << std::endl;
                }

                std::vector <double> xi2;//local coordinates of the face gauss point with respect to jel
                GetClosestPointInReferenceElement(vx2, xg, ielt2, xi2);

                inverseMapping = GetInverseMapping(solTypeV, ielt2, aP2, xg, xi2, 100);
                if(!inverseMapping) {
                  std::cout << "InverseMapping2 failed at " << iel << " " << jel << " " << iface << std::endl;
                }

                msh->_finiteElement[ielt1][solTypeV]->Jacobian(vx1, xi1, weight1, phi1, gradPhi1, nablaPhi1);
                msh->_finiteElement[ielt2][solTypeV]->Jacobian(vx2, xi2, weight2, phi2, gradPhi2, nablaPhi2);


                std::vector < double > solV1gOld(dim, 0.);
                std::vector < double > solV2gOld(dim, 0.);

                std::vector < adept::adouble > solV1g(dim, 0.);
                std::vector < adept::adouble > solV2g(dim, 0.);

                std::vector < adept::adouble > gradSolV1DotN(dim, 0.);
                std::vector < adept::adouble > gradSolV2DotN(dim, 0.);

                std::vector < adept::adouble >  hessSolV1DotN(dim, 0.);
                std::vector < adept::adouble >  hessSolV2DotN(dim, 0.);

                for(unsigned I = 0; I < dim; I++) {
                  for(unsigned i = 0; i < nDofsV1; i++) {
                    solV1g[I] += phi1[i] * solV1[I][i];
                    solV1gOld[I] += phi1[i] * solV1Old[I][i];
                    for(unsigned J = 0; J < dim; J++) {
                      gradSolV1DotN[I] += solV1[I][i] * gradPhi1[i * dim + J] * normal[J];
                      for(unsigned K = 0; K < dim; K++) {
                        //2D xx, yy, xy
                        //3D xx, yy, zz, xy, yz ,zx
                        unsigned L;
                        if(J == K) L = J;
                        else if(1 == J + K) L = dim;     // xy
                        else if(2 == J + K) L = dim + 2; // xz
                        else if(3 == J + K) L = dim + 1; // yz
                        hessSolV1DotN[I] += solV1[I][i] * normal[J] * nablaPhi1[i * dim2 + L] * normal[K];
                      }
                    }
                  }
                }

                for(unsigned I = 0; I < dim; I++) {
                  for(unsigned i = 0; i < nDofsV2; i++) {
                    solV2g[I] += phi2[i] * solV2[I][i];
                    solV2gOld[I] += phi2[i] * solV2Old[I][i];
                    for(unsigned J = 0; J < dim; J++) {
                      gradSolV2DotN[I] += solV2[I][i] * gradPhi2[i * dim + J] * normal[J];
                      for(unsigned K = 0; K < dim; K++) {
                        //2D xx, yy, xy
                        //3D xx, yy, zz, xy, yz ,zx
                        unsigned L;
                        if(J == K) L = J;
                        else if(1 == J + K) L = dim;     // xy
                        else if(2 == J + K) L = dim + 2; // xz
                        else if(3 == J + K) L = dim + 1; // yz
                        hessSolV2DotN[I] += solV2[I][i] * normal[J] * nablaPhi2[i * dim2 + L] * normal[K];
                      }
                    }
                  }
                }


                double V1NormL2 = 0.;
                double V2NormL2 = 0.;

                for(unsigned k = 0; k < dim; k++) {
                  V1NormL2 += solV1gOld[k] * solV1gOld[k];
                  V2NormL2 += solV2gOld[k] * solV2gOld[k];
                }

                V1NormL2 = sqrt(V1NormL2);
                V2NormL2 = sqrt(V2NormL2);


                double psiT1 = (muFluid / rhoFluid + (1. / 6.) * V1NormL2 * h + (1. / 12.) * h * h / (theta * dt));
                double psiT2 = (muFluid / rhoFluid + (1. / 6.) * V2NormL2 * h + (1. / 12.) * h * h / (theta * dt));
                double psiC = 0.5 * h * h * (1. / psiT1 + 1. / psiT2);


                double C1 = gammac * (muFluid + rhoFluid * psiC * 0.5 * (V1NormL2 * V1NormL2 + V2NormL2 * V2NormL2) + rhoFluid * h2 / (theta * dt));


                for(unsigned I = 0; I < dim; I++) {
                  for(unsigned i = 0; i < nDofsV1; i++) {
                    for(unsigned J = 0; J < dim; J++) {
                      aResV1[I][i] +=  C1 * h * gradPhi1[i * dim + J] * normal[J] * (gradSolV1DotN[I] - gradSolV2DotN[I]) * weight;
                      for(unsigned K = 0; K < dim; K++) {

                        unsigned L;
                        if(J == K) L = J;
                        else if(1 == J + K) L = dim;     // xy
                        else if(2 == J + K) L = dim + 2; // xz
                        else if(3 == J + K) L = dim + 1; // yz

                        aResV1[I][i] += C1 * h3 * normal[J] * nablaPhi1[i * dim2 + L] * normal[K] * (hessSolV1DotN[I] - hessSolV2DotN[I]) * weight;
                      }
                    }
                  }

                  for(unsigned i = 0; i < nDofsV2; i++) {
                    for(unsigned J = 0; J < dim; J++) {
                      aResV2[I][i] +=  -C1 * h * gradPhi2[i * dim + J] * normal[J] * (gradSolV1DotN[I] - gradSolV2DotN[I]) * weight;

                      for(unsigned K = 0; K < dim; K++) {

                        unsigned L;
                        if(J == K) L = J;
                        else if(1 == J + K) L = dim;     // xy
                        else if(2 == J + K) L = dim + 2; // xz
                        else if(3 == J + K) L = dim + 1; // yz

                        aResV2[I][i] +=  -C1 * h3 * normal[J] * nablaPhi2[i * dim2 + L] * normal[K] * (hessSolV1DotN[I] - hessSolV2DotN[I]) * weight;
                      }
                    }
                  }
                }

                if(fluid && solTypeP < 3) {
                  msh->_finiteElement[ielt1][solTypeP]->Jacobian(vx1, xi1, weight1, phi1, gradPhi1);
                  msh->_finiteElement[ielt2][solTypeP]->Jacobian(vx2, xi2, weight2, phi2, gradPhi2);

                  adept::adouble gradSolP1DotN = 0.;
                  adept::adouble gradSolP2DotN = 0.;


                  for(unsigned i = 0; i < nDofsP1; i++) {
                    for(unsigned J = 0; J < dim; J++) {
                      gradSolP1DotN += solP1[i] * gradPhi1[i * dim + J] * normal[J];
                    }
                  }

                  for(unsigned i = 0; i < nDofsP2; i++) {
                    for(unsigned J = 0; J < dim; J++) {
                      gradSolP2DotN += solP2[i] * gradPhi2[i * dim + J] * normal[J];
                    }
                  }

                  for(unsigned i = 0; i < nDofsP1; i++) {
                    for(unsigned J = 0; J < dim; J++) {
                      aResP1[i] +=  gammap * psiC / rhoFluid * h *  gradPhi1[i * dim + J] * normal[J] * (gradSolP1DotN - gradSolP2DotN) * weight;
                    }
                  }

                  for(unsigned i = 0; i < nDofsP2; i++) {
                    for(unsigned J = 0; J < dim; J++) {
                      aResP2[i] +=  - gammap * psiC / rhoFluid * h * gradPhi2[i * dim + J] * normal[J] * (gradSolP1DotN - gradSolP2DotN) * weight;
                    }
                  }

                }

              }

              //copy the value of the adept::adoube aRes in double Res and store them in RES
              rhs1.resize(nDofs1);   //resize
              for(int i = 0; i < nDofsV1; i++) {
                for(unsigned  k = 0; k < dim; k++) {
                  rhs1[ k * nDofsV1 + i ] = -aResV1[k][i].value();
                }
              }
              for(int i = 0; i < nDofsP1; i++) {
                rhs1[ dim * nDofsV1 + i] = -aResP1[i].value();
              }
              myRES->add_vector_blocked(rhs1, sysDofs1);


              rhs2.resize(nDofs2);   //resize
              for(int i = 0; i < nDofsV2; i++) {
                for(unsigned  k = 0; k < dim; k++) {
                  rhs2[ k * nDofsV2 + i ] = -aResV2[k][i].value();
                }
              }
              for(int i = 0; i < nDofsP2; i++) {
                rhs2[ dim * nDofsV2 + i] = -aResP2[i].value();
              }
              myRES->add_vector_blocked(rhs2, sysDofs2);


              // define the dependent variables J11 and J12
              for(unsigned  k = 0; k < dim; k++) {
                s.dependent(&aResV1[k][0], nDofsV1);
              }
              s.dependent(&aResP1[0], nDofsP1);


              // define the independent variables J11
              for(unsigned  k = 0; k < dim; k++) {
                s.independent(&solV1[k][0], nDofsV1);
              }
              s.independent(&solP1[0], nDofsP1);
              Jac.resize(nDofs1 * nDofs1);
              // get the and store jacobian matrix (row-major)
              s.jacobian(&Jac[0], true);
              myKK->add_matrix_blocked(Jac, sysDofs1, sysDofs1);
              s.clear_independents();


              // define the independent variables J12
              for(unsigned  k = 0; k < dim; k++) {
                s.independent(&solV2[k][0], nDofsV2);
              }
              s.independent(&solP2[0], nDofsP2);
              Jac.resize(nDofs1 * nDofs2);
              // get the and store jacobian matrix (row-major)
              s.jacobian(&Jac[0], true);
              myKK->add_matrix_blocked(Jac, sysDofs1, sysDofs2);
              s.clear_independents();

              s.clear_dependents(); // for J11 and J12


              // define the dependent variables J21 and J22
              for(unsigned  k = 0; k < dim; k++) {
                s.dependent(&aResV2[k][0], nDofsV2);
              }
              s.dependent(&aResP2[0], nDofsP2);


              // define the independent variables J21
              for(unsigned  k = 0; k < dim; k++) {
                s.independent(&solV1[k][0], nDofsV1);
              }
              s.independent(&solP1[0], nDofsP1);
              Jac.resize(nDofs2 * nDofs1);
              // get the and store jacobian matrix (row-major)
              s.jacobian(&Jac[0], true);
              myKK->add_matrix_blocked(Jac, sysDofs2, sysDofs1);
              s.clear_independents();


              // define the independent variables J22
              for(unsigned  k = 0; k < dim; k++) {
                s.independent(&solV2[k][0], nDofsV2);
              }
              s.independent(&solP2[0], nDofsP2);
              Jac.resize(nDofs2 * nDofs2);
              // get the and store jacobian matrix (row-major)
              s.jacobian(&Jac[0], true);
              myKK->add_matrix_blocked(Jac, sysDofs2, sysDofs2);
              s.clear_independents();

              s.clear_dependents(); // for J21 and J22
            }
          }
        }
      }
    }
  }

  //flagmark
  if(nprocs > 1) {
    for(unsigned kproc = 0; kproc < nprocs; kproc++) {
      for(int iel = msh->_elementOffset[kproc]; iel < msh->_elementOffset[kproc + 1]; iel++) {

        unsigned eFlag1;

        if(iproc == kproc) {
          eFlag1 = static_cast <unsigned>(floor((*mysolution->_Sol[eflagIndex])(iel) + 0.5));
        }
        MPI_Bcast(&eFlag1, 1, MPI_UNSIGNED, kproc, PETSC_COMM_WORLD);

        if(eFlag1 == 1) {
          unsigned nFaces;
          if(iproc == kproc) {
            nFaces = msh->GetElementFaceNumber(iel);
          }
          MPI_Bcast(&nFaces, 1, MPI_UNSIGNED, kproc, PETSC_COMM_WORLD);

          for(unsigned iface = 0; iface < nFaces; iface++) {

            int jel;
            if(iproc == kproc) {
              jel = el->GetFaceElementIndex(iel, iface) - 1;
            }
            MPI_Bcast(&jel, 1, MPI_INT, kproc, PETSC_COMM_WORLD);

            if(jel >= 0) { // iface is not a boundary of the domain
              unsigned jproc = msh->IsdomBisectionSearch(jel , 3); // return  jproc for piece-wise constant discontinuous type (3)
              if(jproc != kproc && (iproc == kproc || iproc == jproc) ) {

                unsigned eFlag2;
                if(iproc == jproc) {
                  eFlag2 = static_cast <unsigned>(floor((*mysolution->_Sol[eflagIndex])(jel) + 0.5));
                  MPI_Send(&eFlag2, 1, MPI_UNSIGNED, kproc, 0, PETSC_COMM_WORLD);
                }
                else if(iproc == kproc) {
                  MPI_Recv(&eFlag2, 1, MPI_UNSIGNED, jproc, 0, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                }

                if(eFlag2 == 0 + !fluid * 2 || (eFlag2 == 1 && jel > iel)) {
                  //std::cout << "I am " << iel << " on " << kproc << " talking with " << jel << " on " << jproc << std::endl;

                  short unsigned ielt1;
                  short unsigned ielt2;

                  if(iproc == kproc) {
                    ielt1 = msh->GetElementType(iel);
                    MPI_Recv(&ielt2, 1, MPI_UNSIGNED_SHORT, jproc, 0, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                  }
                  else if (iproc == jproc) {
                    ielt2 = msh->GetElementType(jel);
                    MPI_Send(&ielt2, 1, MPI_UNSIGNED_SHORT, kproc, 0, PETSC_COMM_WORLD);
                  }

                  unsigned nDofsV1;
                  unsigned nDofsP1;
                  unsigned nDofs1;

                  unsigned nDofsV2 = el->GetNVE(ielt2, solTypeV);
                  unsigned nDofsP2 = el->GetNVE(ielt2, solTypeP);
                  unsigned nDofs2 = dim * nDofsV2 + nDofsP2;

                  sysDofs2.resize(nDofs2);
                  for(unsigned  k = 0; k < dim; k++) {
                    solV2d[k].resize(nDofsV2);
                    solV2Old[k].resize(nDofsV2);
                    vx2[k].resize(nDofsV2);
                  }
                  solP2d.resize(nDofsP2);
                  std::vector < MPI_Request > reqs(3 * dim + 2);
                  if(iproc == kproc) {

                    nDofsV1 = el->GetNVE(ielt1, solTypeV);
                    nDofsP1 = el->GetNVE(ielt1, solTypeP);
                    nDofs1 = dim * nDofsV1 + nDofsP1;

                    sysDofs1.resize(nDofs1);

                    for(unsigned  k = 0; k < dim; k++) {
                      solV1[k].resize(nDofsV1);
                      solV1Old[k].resize(nDofsV1);
                      vx1[k].resize(nDofsV1);
                    }
                    solP1.resize(nDofsP1);

                    for(unsigned i = 0; i < nDofsV1; i++) {
                      unsigned idof = msh->GetSolutionDof(i, iel, solTypeV);

                      for(unsigned  k = 0; k < dim; k++) {
                        solV1[k][i] = (*mysolution->_Sol[indexSolV[k]])(idof);
                        solV1Old[k][i] = (*mysolution->_SolOld[indexSolV[k]])(idof);
                        sysDofs1[k * nDofsV1 + i] = myLinEqSolver->GetSystemDof(indexSolV[k], indexPdeV[k], i, iel);
                      }
                    }

                    for(unsigned i = 0; i < nDofsP1; i++) {
                      unsigned idof = msh->GetSolutionDof(i, iel, solTypeP);
                      solP1[i] = (*mysolution->_Sol[indexSolP])(idof);
                      sysDofs1[dim * nDofsV1 + i] = myLinEqSolver->GetSystemDof(indexSolP, indexPdeP, i, iel);
                    }

                    for(unsigned i = 0; i < nDofsV1; i++) {
                      unsigned idofX = msh->GetSolutionDof(i, iel, 2);
                      for(unsigned  k = 0; k < dim; k++) {
                        vx1[k][i] = (*msh->_topology->_Sol[k])(idofX);
                      }
                    }
                    for(unsigned k = 0; k < dim; k++) {
                      MPI_Irecv(solV2d[k].data(), solV2d[k].size(), MPI_DOUBLE, jproc, k, PETSC_COMM_WORLD, &reqs[k]);
                      MPI_Irecv(solV2Old[k].data(), solV2Old[k].size(), MPI_DOUBLE, jproc, k + dim, PETSC_COMM_WORLD, &reqs[k + dim]);
                      MPI_Irecv(vx2[k].data(), vx2[k].size(), MPI_DOUBLE, jproc, k + 2 * dim, PETSC_COMM_WORLD, &reqs[k + 2 * dim]);
                    }
                    MPI_Irecv(solP2d.data(), solP2d.size(), MPI_DOUBLE, jproc, 3 * dim, PETSC_COMM_WORLD,  &reqs[3 * dim]);
                    MPI_Irecv(sysDofs2.data(), sysDofs2.size(), MPI_UNSIGNED, jproc, 3 * dim + 1, PETSC_COMM_WORLD,  &reqs[3 * dim + 1]);
                  }
                  else if(iproc == jproc) {
                    for(unsigned i = 0; i < nDofsV2; i++) {
                      unsigned idof = msh->GetSolutionDof(i, jel, solTypeV);

                      for(unsigned  k = 0; k < dim; k++) {
                        solV2d[k][i] = (*mysolution->_Sol[indexSolV[k]])(idof);
                        solV2Old[k][i] = (*mysolution->_SolOld[indexSolV[k]])(idof);
                        sysDofs2[k * nDofsV2 + i] = myLinEqSolver->GetSystemDof(indexSolV[k], indexPdeV[k], i, jel);
                      }
                    }

                    for(unsigned i = 0; i < nDofsP2; i++) {
                      unsigned idof = msh->GetSolutionDof(i, jel, solTypeP);
                      solP2d[i] = (*mysolution->_Sol[indexSolP])(idof);
                      sysDofs2[dim * nDofsV2 + i] = myLinEqSolver->GetSystemDof(indexSolP, indexPdeP, i, jel);
                    }

                    for(unsigned i = 0; i < nDofsV2; i++) {
                      unsigned idofX = msh->GetSolutionDof(i, jel, 2);
                      for(unsigned  k = 0; k < dim; k++) {
                        vx2[k][i] = (*msh->_topology->_Sol[k])(idofX);
                      }
                    }

                    for(unsigned k = 0; k < dim; k++) {
                      MPI_Isend(solV2d[k].data(), solV2d[k].size(), MPI_DOUBLE, kproc, k, PETSC_COMM_WORLD, &reqs[k]);
                      MPI_Isend(solV2Old[k].data(), solV2Old[k].size(), MPI_DOUBLE, kproc, k + dim, PETSC_COMM_WORLD, &reqs[k + dim]);
                      MPI_Isend(vx2[k].data(), vx2[k].size(), MPI_DOUBLE, kproc, k + 2 * dim, PETSC_COMM_WORLD, &reqs[k + 2 * dim]);
                    }
                    MPI_Isend(solP2d.data(), solP2d.size(), MPI_DOUBLE, kproc, 3 * dim, PETSC_COMM_WORLD, &reqs[3 * dim]);
                    MPI_Isend(sysDofs2.data(), sysDofs2.size(), MPI_UNSIGNED, kproc, 3 * dim + 1, PETSC_COMM_WORLD, &reqs[3 * dim + 1]);
                  }

                  MPI_Status status;
                  for(unsigned m = 0; m < 3 * dim + 2; m++) {
                    MPI_Wait(&reqs[m], &status);
                  }

                  if(iproc == kproc) {

                    for(unsigned  k = 0; k < dim; k++) {
                      solV2[k].resize(nDofsV2);
                      for(unsigned i = 0; i < nDofsV2; i++) {
                        solV2[k][i] = solV2d[k][i];
                      }
                    }
                    solP2.resize(nDofsP2);
                    for(unsigned i = 0; i < nDofsP2; i++) {
                      solP2[i] = solP2d[i];
                    }


                    for(unsigned  k = 0; k < dim; k++) {
                      aResV1[k].assign(nDofsV1, 0.);
                      aResV2[k].assign(nDofsV2, 0.);
                    }
                    aResP1.assign(nDofsP1, 0.);
                    aResP2.assign(nDofsP2, 0.);

                    s.new_recording();

                    const unsigned faceGeom = el->GetFaceType(ielt1, iface);
                    unsigned faceDofs = el->GetNFACENODES(ielt1, iface, solTypeV);

                    std::vector  < std::vector  <  double > > faceVx(dim);    // A matrix holding the face coordinates rowwise.
                    for(int k = 0; k < dim; k++) {
                      faceVx[k].resize(faceDofs);
                    }
                    for(unsigned i = 0; i < faceDofs; i++) {
                      unsigned inode = msh->GetLocalFaceVertexIndex(iel, iface, i);    // face-to-element local node mapping.
                      for(unsigned k = 0; k < dim; k++) {
                        faceVx[k][i] =  vx1[k][inode]; // We extract the local coordinates on the face from local coordinates on the element.
                      }
                    }

                    double h = sqrt((faceVx[0][0] - faceVx[0][1]) * (faceVx[0][0] - faceVx[0][1]) +
                                    (faceVx[1][0] - faceVx[1][1]) * (faceVx[1][0] - faceVx[1][1])) ;
                    double h2 = h * h;
                    double h3 = h * h * h;

                    //std::cout << " h = " << h << " ";

                    for(unsigned jtype = 0; jtype < solTypeV + 1; jtype++) {
                      ProjectNodalToPolynomialCoefficients(aP1[jtype], vx1, ielt1, jtype);
                      ProjectNodalToPolynomialCoefficients(aP2[jtype], vx2, ielt2, jtype);
                    }

                    for(unsigned ig = 0; ig  <  msh->_finiteElement[faceGeom][solTypeV]->GetGaussPointNumber(); ig++) {

                      std::vector < double> normal;
                      msh->_finiteElement[faceGeom][solTypeV]->JacobianSur(faceVx, ig, weight, phi, gradPhi, normal);

                      std::vector< double > xg(dim, 0.); // physical coordinates of the face Gauss point
                      for(unsigned i = 0; i < faceDofs; i++) {
                        for(unsigned k = 0; k < dim; k++) {
                          xg[k] += phi[i] * faceVx[k][i];
                        }
                      }

                      std::vector <double> xi1;//local coordinates of the face gauss point with respect to iel
                      GetClosestPointInReferenceElement(vx1, xg, ielt1, xi1);

                      bool inverseMapping = GetInverseMapping(solTypeV, ielt1, aP1, xg, xi1, 100);
                      if(!inverseMapping) {
                        std::cout << "InverseMapping1 failed at " << iel << " " << jel << " " << iface << std::endl;
                      }

                      std::vector <double> xi2;//local coordinates of the face gauss point with respect to jel
                      GetClosestPointInReferenceElement(vx2, xg, ielt2, xi2);
                      
                      inverseMapping = GetInverseMapping(solTypeV, ielt2, aP2, xg, xi2, 100);
                      if(!inverseMapping) {
                        std::cout << "InverseMapping2 failed at " << iel << " " << jel << " " << iface << std::endl;
                      }

                      msh->_finiteElement[ielt1][solTypeV]->Jacobian(vx1, xi1, weight1, phi1, gradPhi1, nablaPhi1);
                      msh->_finiteElement[ielt2][solTypeV]->Jacobian(vx2, xi2, weight2, phi2, gradPhi2, nablaPhi2);

                      std::vector < double > solV1gOld(dim, 0.);
                      std::vector < double > solV2gOld(dim, 0.);

                      std::vector < adept::adouble > solV1g(dim, 0.);
                      std::vector < adept::adouble > solV2g(dim, 0.);

                      std::vector < adept::adouble > gradSolV1DotN(dim, 0.);
                      std::vector < adept::adouble > gradSolV2DotN(dim, 0.);

                      std::vector < adept::adouble >  hessSolV1DotN(dim, 0.);
                      std::vector < adept::adouble >  hessSolV2DotN(dim, 0.);

                      for(unsigned I = 0; I < dim; I++) {
                        for(unsigned i = 0; i < nDofsV1; i++) {
                          solV1g[I] += phi1[i] * solV1[I][i];
                          solV1gOld[I] += phi1[i] * solV1Old[I][i];
                          for(unsigned J = 0; J < dim; J++) {
                            gradSolV1DotN[I] += solV1[I][i] * gradPhi1[i * dim + J] * normal[J];
                            for(unsigned K = 0; K < dim; K++) {
                              //2D xx, yy, xy
                              //3D xx, yy, zz, xy, yz ,zx
                              unsigned L;
                              if(J == K) L = J;
                              else if(1 == J + K) L = dim;     // xy
                              else if(2 == J + K) L = dim + 2; // xz
                              else if(3 == J + K) L = dim + 1; // yz
                              hessSolV1DotN[I] += solV1[I][i] * normal[J] * nablaPhi1[i * dim2 + L] * normal[K];
                            }
                          }
                        }
                      }

                      for(unsigned I = 0; I < dim; I++) {
                        for(unsigned i = 0; i < nDofsV2; i++) {
                          solV2g[I] += phi2[i] * solV2[I][i];
                          solV2gOld[I] += phi2[i] * solV2Old[I][i];
                          for(unsigned J = 0; J < dim; J++) {
                            gradSolV2DotN[I] += solV2[I][i] * gradPhi2[i * dim + J] * normal[J];
                            for(unsigned K = 0; K < dim; K++) {
                              //2D xx, yy, xy
                              //3D xx, yy, zz, xy, yz ,zx
                              unsigned L;
                              if(J == K) L = J;
                              else if(1 == J + K) L = dim;     // xy
                              else if(2 == J + K) L = dim + 2; // xz
                              else if(3 == J + K) L = dim + 1; // yz
                              hessSolV2DotN[I] += solV2[I][i] * normal[J] * nablaPhi2[i * dim2 + L] * normal[K];
                            }
                          }
                        }
                      }

                      double V1NormL2 = 0.;
                      double V2NormL2 = 0.;

                      for(unsigned k = 0; k < dim; k++) {
                        V1NormL2 += solV1gOld[k] * solV1gOld[k];
                        V2NormL2 += solV2gOld[k] * solV2gOld[k];
                      }

                      V1NormL2 = sqrt(V1NormL2);
                      V2NormL2 = sqrt(V2NormL2);


                      double psiT1 = (muFluid / rhoFluid + (1. / 6.) * V1NormL2 * h + (1. / 12.) * h * h / (theta * dt));
                      double psiT2 = (muFluid / rhoFluid + (1. / 6.) * V2NormL2 * h + (1. / 12.) * h * h / (theta * dt));
                      double psiC = 0.5 * h * h * (1. / psiT1 + 1. / psiT2);

                      double C1 = gammac * (muFluid + rhoFluid * psiC * 0.5 * (V1NormL2 * V1NormL2 + V2NormL2 * V2NormL2) + rhoFluid * h2 / (theta * dt));

                      for(unsigned I = 0; I < dim; I++) {
                        for(unsigned i = 0; i < nDofsV1; i++) {
                          for(unsigned J = 0; J < dim; J++) {
                            aResV1[I][i] +=  C1 * h * gradPhi1[i * dim + J] * normal[J] * (gradSolV1DotN[I] - gradSolV2DotN[I]) * weight;
                            for(unsigned K = 0; K < dim; K++) {

                              unsigned L;
                              if(J == K) L = J;
                              else if(1 == J + K) L = dim;     // xy
                              else if(2 == J + K) L = dim + 2; // xz
                              else if(3 == J + K) L = dim + 1; // yz

                              aResV1[I][i] += C1 * h3 * normal[J] * nablaPhi1[i * dim2 + L] * normal[K] * (hessSolV1DotN[I] - hessSolV2DotN[I]) * weight;
                            }
                          }
                        }

                        for(unsigned i = 0; i < nDofsV2; i++) {
                          for(unsigned J = 0; J < dim; J++) {
                            aResV2[I][i] +=  -C1 * h * gradPhi2[i * dim + J] * normal[J] * (gradSolV1DotN[I] - gradSolV2DotN[I]) * weight;

                            for(unsigned K = 0; K < dim; K++) {

                              unsigned L;
                              if(J == K) L = J;
                              else if(1 == J + K) L = dim;     // xy
                              else if(2 == J + K) L = dim + 2; // xz
                              else if(3 == J + K) L = dim + 1; // yz

                              aResV2[I][i] +=  -C1 * h3 * normal[J] * nablaPhi2[i * dim2 + L] * normal[K] * (hessSolV1DotN[I] - hessSolV2DotN[I]) * weight;
                            }
                          }
                        }
                      }

                      if(fluid && solTypeP < 3) {
                        msh->_finiteElement[ielt1][solTypeP]->Jacobian(vx1, xi1, weight1, phi1, gradPhi1);
                        msh->_finiteElement[ielt2][solTypeP]->Jacobian(vx2, xi2, weight2, phi2, gradPhi2);

                        adept::adouble gradSolP1DotN = 0.;
                        adept::adouble gradSolP2DotN = 0.;


                        for(unsigned i = 0; i < nDofsP1; i++) {
                          for(unsigned J = 0; J < dim; J++) {
                            gradSolP1DotN += solP1[i] * gradPhi1[i * dim + J] * normal[J];
                          }
                        }

                        for(unsigned i = 0; i < nDofsP2; i++) {
                          for(unsigned J = 0; J < dim; J++) {
                            gradSolP2DotN += solP2[i] * gradPhi2[i * dim + J] * normal[J];
                          }
                        }

                        for(unsigned i = 0; i < nDofsP1; i++) {
                          for(unsigned J = 0; J < dim; J++) {
                            aResP1[i] +=  gammap * psiC / rhoFluid * h *  gradPhi1[i * dim + J] * normal[J] * (gradSolP1DotN - gradSolP2DotN) * weight;
                          }
                        }

                        for(unsigned i = 0; i < nDofsP2; i++) {
                          for(unsigned J = 0; J < dim; J++) {
                            aResP2[i] +=  - gammap * psiC / rhoFluid * h * gradPhi2[i * dim + J] * normal[J] * (gradSolP1DotN - gradSolP2DotN) * weight;
                          }
                        }

                      }

                    }

                    //copy the value of the adept::adoube aRes in double Res and store them in RES
                    rhs1.resize(nDofs1);   //resize
                    for(int i = 0; i < nDofsV1; i++) {
                      for(unsigned  k = 0; k < dim; k++) {
                        rhs1[ k * nDofsV1 + i ] = -aResV1[k][i].value();
                      }
                    }
                    for(int i = 0; i < nDofsP1; i++) {
                      rhs1[ dim * nDofsV1 + i] = -aResP1[i].value();
                    }
                    myRES->add_vector_blocked(rhs1, sysDofs1);

                    rhs2.resize(nDofs2);   //resize
                    for(int i = 0; i < nDofsV2; i++) {
                      for(unsigned  k = 0; k < dim; k++) {
                        rhs2[ k * nDofsV2 + i ] = -aResV2[k][i].value();
                      }
                    }
                    for(int i = 0; i < nDofsP2; i++) {
                      rhs2[ dim * nDofsV2 + i] = -aResP2[i].value();
                    }
                    myRES->add_vector_blocked(rhs2, sysDofs2);

                    // define the dependent variables J11 and J12
                    for(unsigned  k = 0; k < dim; k++) {
                      s.dependent(&aResV1[k][0], nDofsV1);
                    }
                    s.dependent(&aResP1[0], nDofsP1);

                    // define the independent variables J11
                    for(unsigned  k = 0; k < dim; k++) {
                      s.independent(&solV1[k][0], nDofsV1);
                    }
                    s.independent(&solP1[0], nDofsP1);
                    Jac.resize(nDofs1 * nDofs1);
                    // get the and store jacobian matrix (row-major)
                    s.jacobian(&Jac[0], true);
                    myKK->add_matrix_blocked(Jac, sysDofs1, sysDofs1);
                    s.clear_independents();

                    // define the independent variables J12
                    for(unsigned  k = 0; k < dim; k++) {
                      s.independent(&solV2[k][0], nDofsV2);
                    }
                    s.independent(&solP2[0], nDofsP2);
                    Jac.resize(nDofs1 * nDofs2);
                    // get the and store jacobian matrix (row-major)
                    s.jacobian(&Jac[0], true);
                    myKK->add_matrix_blocked(Jac, sysDofs1, sysDofs2);
                    s.clear_independents();

                    s.clear_dependents(); // for J11 and J12

                    // define the dependent variables J21 and J22
                    for(unsigned  k = 0; k < dim; k++) {
                      s.dependent(&aResV2[k][0], nDofsV2);
                    }
                    s.dependent(&aResP2[0], nDofsP2);

                    // define the independent variables J21
                    for(unsigned  k = 0; k < dim; k++) {
                      s.independent(&solV1[k][0], nDofsV1);
                    }
                    s.independent(&solP1[0], nDofsP1);
                    Jac.resize(nDofs2 * nDofs1);
                    // get the and store jacobian matrix (row-major)
                    s.jacobian(&Jac[0], true);
                    myKK->add_matrix_blocked(Jac, sysDofs2, sysDofs1);
                    s.clear_independents();

                    // define the independent variables J22
                    for(unsigned  k = 0; k < dim; k++) {
                      s.independent(&solV2[k][0], nDofsV2);
                    }
                    s.independent(&solP2[0], nDofsP2);
                    Jac.resize(nDofs2 * nDofs2);
                    // get the and store jacobian matrix (row-major)
                    s.jacobian(&Jac[0], true);
                    myKK->add_matrix_blocked(Jac, sysDofs2, sysDofs2);
                    s.clear_independents();

                    s.clear_dependents(); // for J21 and J22
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  // *************************************
  std::cout << "Ghost Penalty Assembly time = " << static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC << std::endl;

}

