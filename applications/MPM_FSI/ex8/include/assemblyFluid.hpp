#include "MultiLevelSolution.hpp"
#include "PetscMatrix.hpp"

using namespace femus;


void AssembleGhostPenalty(MultiLevelProblem& ml_prob) {

  // ml_prob is the global object from/to where get/set all the data
  // level is the level of the PDE system to be assembled
  // levelMax is the Maximum level of the MultiLevelProblem
  // assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  clock_t AssemblyTime = 0;
  clock_t start_time, end_time;

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

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  const unsigned dim = msh->GetDimension();
  const unsigned dim2 = 3 * (dim - 1);

  // data
  unsigned iproc  = msh->processor_id();

  //quantities for iel will have index1
  //quantities for jel will have index2
  vector< vector< double > > solV1Old(dim);      // local solution (displacement)
  vector< vector< double > > solV2Old(dim);      // local solution (velocity)

  //vector< vector< double > > solD1Old(dim);      // local solution (displacement)
  //vector< vector< double > > solD2Old(dim);



  vector< vector< adept::adouble > > solV1(dim); // local solution (velocity)
  vector< adept::adouble > solP1; // local solution (velocity)

  vector< vector< adept::adouble > > solV2(dim); // local solution (velocity)
  vector< adept::adouble > solP2; // local solution (velocity)

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
  double rhoFluid = ml_prob.parameters.get<Fluid> ("FluidFEM").get_density();
  double muFluid = ml_prob.parameters.get<Fluid> ("FluidFEM").get_viscosity();

  double dt =  my_nnlin_impl_sys.GetIntervalTime();

  double gammac = 50;
  double gammap = 50;

  std::cout.precision(10);

  //variable-name handling
  const char varname[3][5] = {"UX", "UY", "UZ"};

  vector <unsigned> indexSolV(dim);
  //vector <unsigned> indexSolD(dim);
  vector <unsigned> indexPdeV(dim);
  for(unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolV[ivar] = mlSol->GetIndex(&varname[ivar][0]);
    indexPdeV[ivar] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar][0]);
  }
  unsigned solTypeV = mlSol->GetSolutionType(&varname[0][0]);

  unsigned indexSolP = mlSol->GetIndex("P");
  unsigned indexPdeP = my_nnlin_impl_sys.GetSolPdeIndex("P");
  unsigned solTypeP = mlSol->GetSolutionType("P");

  unsigned eflagIndex = mlSol->GetIndex("eflag");
//

  start_time = clock();

  std::vector < std::vector < std::vector <double > > > aP1(3);
  std::vector < std::vector < std::vector <double > > > aP2(3);

  //BEGIN loop on elements (to initialize the "soft" stiffness matrix)
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
        //solD1Old[k].resize(nDofsV1);
        vx1[k].resize(nDofsV1);
      }
      solP1.resize(nDofsP1);

      for(unsigned i = 0; i < nDofsV1; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, solTypeV);

        //nodeFlag1[i] = (*mysolution->_Sol[nflagIndex])(idof);

        for(unsigned  k = 0; k < dim; k++) {
          solV1[k][i] = (*mysolution->_Sol[indexSolV[k]])(idof);
          solV1Old[k][i] = (*mysolution->_SolOld[indexSolV[k]])(idof);
          //solD1Old[k][i] = (*mysolution->_SolOld[indexSolD[k]])(idof);
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

      double h = sqrt((vx1[0][0] - vx1[0][2]) * (vx1[0][0] - vx1[0][2]) +
                      (vx1[1][0] - vx1[1][2]) * (vx1[1][0] - vx1[1][2])) ;
      double h2 = h * h;
      double h3 = h * h * h;

      bool aP1IsInitialized = false;

      //std::cout << iel<< " ";
      for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {

        int jel = el->GetFaceElementIndex(iel, iface) - 1;

        unsigned eFlag2 = (jel >= 0) ? static_cast <unsigned>(floor((*mysolution->_Sol[eflagIndex])(jel) + 0.5)) : 0;

        if(eFlag2 == 2 || (eFlag2 == 1 && jel > iel)) {

          //std::cout << jel << " ";

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
            //solD2Old[k].resize(nDofsV2);
            vx2[k].resize(nDofsV2);
          }
          solP2.resize(nDofsP2);


          for(unsigned i = 0; i < nDofsV2; i++) {
            unsigned idof = msh->GetSolutionDof(i, jel, solTypeV);
            //nodeFlag2[i] = (*mysolution->_Sol[nflagIndex])(idof);
            for(unsigned  k = 0; k < dim; k++) {
              solV2[k][i] = (*mysolution->_Sol[indexSolV[k]])(idof);
              solV2Old[k][i] = (*mysolution->_SolOld[indexSolV[k]])(idof);
              //solD2Old[k][i] = (*mysolution->_SolOld[indexSolD[k]])(idof);
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

            std::vector < double> N; //N12
            msh->_finiteElement[faceGeom][solTypeV]->JacobianSur(faceVx, ig, weight, phi, gradPhi, N);

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

            std::vector < adept::adouble > solV1g(dim, 0.);
            std::vector < adept::adouble > solV2g(dim, 0.);

            std::vector < double > solV1gOld(dim, 0.);
            std::vector < double > solV2gOld(dim, 0.);

            //std::vector < double > solD1gOld(dim, 0.);
            //std::vector < double > solD2gOld(dim, 0.);

            std::vector < adept::adouble > gradSolV1DotN(dim, 0.);
            std::vector < adept::adouble > gradSolV2DotN(dim, 0.);

            std::vector < adept::adouble >  hessSolV1DotN(dim, 0.);
            std::vector < adept::adouble >  hessSolV2DotN(dim, 0.);

            for(unsigned I = 0; I < dim; I++) {
              for(unsigned i = 0; i < nDofsV1; i++) {
                solV1g[I] += phi1[i] * solV1[I][i];
                solV1gOld[I] += phi1[i] * solV1Old[I][i];
                //solD1gOld[I] += phi1[i] * solD1Old[I][i];
                for(unsigned J = 0; J < dim; J++) {
                  gradSolV1DotN[I] += solV1[I][i] * gradPhi1[i * dim + J] * N[J];
                  for(unsigned K = 0; K < dim; K++) {
                    //2D xx, yy, xy
                    //3D xx, yy, zz, xy, yz ,zx
                    unsigned L;
                    if(J == K) L = J;
                    else if(1 == J + K) L = dim;     // xy
                    else if(2 == J + K) L = dim + 2; // xz
                    else if(3 == J + K) L = dim + 1; // yz
                    hessSolV1DotN[I] += solV1[I][i] * N[J] * nablaPhi1[i * dim2 + L] * N[K];
                  }
                }
              }
            }




            for(unsigned I = 0; I < dim; I++) {
              for(unsigned i = 0; i < nDofsV2; i++) {
                solV2g[I] += phi2[i] * solV2[I][i];
                solV2gOld[I] += phi2[i] * solV2Old[I][i];
                //solD2gOld[I] += phi2[i] * solD2Old[I][i];
                for(unsigned J = 0; J < dim; J++) {
                  gradSolV2DotN[I] += solV2[I][i] * gradPhi2[i * dim + J] * N[J];
                  for(unsigned K = 0; K < dim; K++) {
                    //2D xx, yy, xy
                    //3D xx, yy, zz, xy, yz ,zx
                    unsigned L;
                    if(J == K) L = J;
                    else if(1 == J + K) L = dim;     // xy
                    else if(2 == J + K) L = dim + 2; // xz
                    else if(3 == J + K) L = dim + 1; // yz
                    hessSolV2DotN[I] += solV2[I][i] * N[J] * nablaPhi2[i * dim2 + L] * N[K];
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

            
            double C1 = gammac * ( muFluid + rhoFluid * psiC * V1NormL2 * V1NormL2 + rhoFluid * h2 / (theta * dt) );
            
            
            for(unsigned I = 0; I < dim; I++) {
              for(unsigned i = 0; i < nDofsV1; i++) {
                for(unsigned J = 0; J < dim; J++) {
                  aResV1[I][i] += C1 * h * gradPhi1[i * dim + J] * N[J] * (gradSolV1DotN[I] - gradSolV2DotN[I]) * weight;
                  for(unsigned K = 0; K < dim; K++) {

                    unsigned L;
                    if(J == K) L = J;
                    else if(1 == J + K) L = dim;     // xy
                    else if(2 == J + K) L = dim + 2; // xz
                    else if(3 == J + K) L = dim + 1; // yz

                    aResV1[I][i] += C1 * h3 * N[J] * nablaPhi1[i * dim2 + L] * N[K] * (hessSolV1DotN[I] - hessSolV2DotN[I]) * weight;
                  }
                }
              }

              for(unsigned i = 0; i < nDofsV2; i++) {
                for(unsigned J = 0; J < dim; J++) {
                  aResV2[I][i] +=  - C1 * h * gradPhi2[i * dim + J] * N[J] * (gradSolV1DotN[I] - gradSolV2DotN[I]) * weight;

                  for(unsigned K = 0; K < dim; K++) {

                    unsigned L;
                    if(J == K) L = J;
                    else if(1 == J + K) L = dim;     // xy
                    else if(2 == J + K) L = dim + 2; // xz
                    else if(3 == J + K) L = dim + 1; // yz

                    aResV2[I][i] += - C1 * h3 * N[J] * nablaPhi2[i * dim2 + L] * N[K] * (hessSolV1DotN[I] - hessSolV2DotN[I]) * weight;
                  }
                }
              }
            }

            if(solTypeP < 3) {

              double psiP = psiC;
              msh->_finiteElement[ielt1][solTypeP]->Jacobian(vx1, xi1, weight1, phi1, gradPhi1);
              msh->_finiteElement[ielt2][solTypeP]->Jacobian(vx2, xi2, weight2, phi2, gradPhi2);

              adept::adouble gradSolP1DotN = 0.;
              adept::adouble gradSolP2DotN = 0.;


              for(unsigned i = 0; i < nDofsP1; i++) {
                for(unsigned J = 0; J < dim; J++) {
                  gradSolP1DotN += solP1[i] * gradPhi1[i * dim + J] * N[J];
                }
              }

              for(unsigned i = 0; i < nDofsP2; i++) {
                for(unsigned J = 0; J < dim; J++) {
                  gradSolP2DotN += solP2[i] * gradPhi2[i * dim + J] * N[J];
                }
              }

              for(unsigned i = 0; i < nDofsP1; i++) {
                for(unsigned J = 0; J < dim; J++) {
                  aResP1[i] += gammap * psiC / rhoFluid * h * gradPhi1[i * dim + J] * N[J] * (gradSolP1DotN - gradSolP2DotN) * weight;
                }
              }

              for(unsigned i = 0; i < nDofsP2; i++) {
                for(unsigned J = 0; J < dim; J++) {
                  aResP2[i] +=  - gammap * psiC / rhoFluid * h * gradPhi2[i * dim + J] * N[J] * (gradSolP1DotN - gradSolP2DotN) * weight;
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
      //std::cout << std::endl;
    }
  }

// *************************************
  end_time = clock();
  AssemblyTime += (end_time - start_time);


}





void AssembleFluid(MultiLevelProblem& ml_prob) {

  // ml_prob is the global object from/to where get/set all the data
  // level is the level of the PDE system to be assembled
  // levelMax is the Maximum level of the MultiLevelProblem
  // assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  clock_t AssemblyTime = 0;
  clock_t start_time, end_time;

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

  // Remove this for no-marker simulation

  // AssembleGhostPenalty(ml_prob, true);
  // AssembleGhostPenalty(ml_prob, false);

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  const unsigned dim = msh->GetDimension();

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  // data
  unsigned iproc  = msh->processor_id();


  vector< vector< adept::adouble > > solV(dim);      // local solution (velocity)
  vector< adept::adouble > solP;
  vector< vector< double > > solVOld(dim);

  vector< double > rhs;    // local redidual vector
  vector< vector< adept::adouble > > aResV(dim);     // local redidual vector
  vector< adept::adouble > aResP;    // local redidual vector
  vector < double > Jac;

  std::vector <unsigned> sysDofsAll;

  vector < double > phiV;
  vector < double > phiP;

  vector < double > gradPhiP;
  vector < double > gradPhiV;

  vector < double> nablaphiP;
  vector < double> nablaphiV;

  unsigned dim2 = 3 * (dim - 1);

  vector <vector < double> > vxHat(dim);

  double weightP;
  double weightV;

  //reading parameters for fluid FEM domain
  double rhoFluid = ml_prob.parameters.get<Fluid> ("FluidFEM").get_density();
  double muFluid = ml_prob.parameters.get<Fluid> ("FluidFEM").get_viscosity();

  double dt =  my_nnlin_impl_sys.GetIntervalTime();

  std::cout.precision(10);

  //variable-name handling
  const char varname[12][5] = {"UX", "UY", "UZ"};

  vector <unsigned> indexSolV(dim);
  vector <unsigned> indexPdeD(dim);
  vector <unsigned> indexPdeV(dim);
  for(unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolV[ivar] = mlSol->GetIndex(&varname[ivar][0]);
    indexPdeV[ivar] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar][0]);
  }
  unsigned solType = mlSol->GetSolutionType(&varname[0][0]);

  unsigned indexSolP = mlSol->GetIndex("P");
  unsigned indexPdeP = my_nnlin_impl_sys.GetSolPdeIndex("P");
  unsigned solTypeP = mlSol->GetSolutionType("P");

  unsigned eflagIndex = mlSol->GetIndex("eflag");
  unsigned nflagIndex = mlSol->GetIndex("nflag");
  std::vector < unsigned >  nodeFlag; // local solution

  start_time = clock();

  //BEGIN loop on elements (to initialize the "soft" stiffness matrix)
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielt = msh->GetElementType(iel);

    unsigned nDofs = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solTypeP);    // number of solution element dofs
    unsigned nDofsAll = dim * nDofs + nDofsP;

    // resize local arrays
    sysDofsAll.resize(nDofsAll);
    nodeFlag.resize(nDofs);

    unsigned eFlag = static_cast <unsigned>(floor((*mysolution->_Sol[eflagIndex])(iel) + 0.5));

    if(eFlag == 2) {

      for(unsigned  k = 0; k < dim; k++) {
        solV[k].resize(nDofs);
        solVOld[k].resize(nDofs);
        aResV[k].assign(nDofs, 0.);
        vxHat[k].resize(nDofs);
      }
      solP.resize(nDofsP);
      aResP.assign(nDofsP, 0.);

      for(unsigned i = 0; i < nDofs; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, solType);

        nodeFlag[i] = (*mysolution->_Sol[nflagIndex])(idof); // set it to 0 for no-marker

        for(unsigned  k = 0; k < dim; k++) {
          solV[k][i] = (*mysolution->_Sol[indexSolV[k]])(idof);
          solVOld[k][i] = (*mysolution->_SolOld[indexSolV[k]])(idof);
          sysDofsAll[k * nDofs + i] = myLinEqSolver->GetSystemDof(indexSolV[k], indexPdeV[k], i, iel);
        }
      }

      for(unsigned i = 0; i < nDofsP; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, solTypeP);
        solP[i] = (*mysolution->_Sol[indexSolP])(idof);
        sysDofsAll[dim * nDofs + i] = myLinEqSolver->GetSystemDof(indexSolP, indexPdeP, i, iel);
      }

      // start a new recording of all the operations involving adept::adouble variables
      s.new_recording();

      for(unsigned i = 0; i < nDofs; i++) {
        unsigned idofX = msh->GetSolutionDof(i, iel, 2);
        for(unsigned  k = 0; k < dim; k++) {
          vxHat[k][i] = (*msh->_topology->_Sol[k])(idofX);
        }
      }

      // *** Gauss point loop ***
      for(unsigned ig = 0; ig < msh->_finiteElement[ielt][solType]->GetGaussPointNumber(); ig++) {

        msh->_finiteElement[ielt][solTypeP]->Jacobian(vxHat, ig, weightP, phiP, gradPhiP);
        msh->_finiteElement[ielt][solType]->Jacobian(vxHat, ig, weightV, phiV, gradPhiV, nablaphiV);

        vector < adept::adouble > solVg(dim, 0.);
        vector < double > solVgOld(dim, 0.);
        vector < adept::adouble > solVgTheta(dim, 0.);

        vector < vector < adept::adouble > > gradSolVgTheta(dim);
        vector<vector<adept::adouble> > DeltaSolVgTheta(dim);

        for(unsigned  k = 0; k < dim; k++) {
          gradSolVgTheta[k].assign(dim, 0.);
          DeltaSolVgTheta[k].resize(dim2, 0.);
        }

        for(unsigned i = 0; i < nDofs; i++) {
          for(unsigned j = 0; j < dim; j++) {
            solVg[j] += phiV[i] * solV[j][i];
            solVgOld[j] += phiV[i] * solVOld[j][i];
            solVgTheta[j] += phiV[i] * (theta * solV[j][i] + (1. - theta) * solVOld[j][i]);
            for(unsigned  k = 0; k < dim; k++) {
              gradSolVgTheta[k][j] += gradPhiV[i * dim + j] * (theta * solV[k][i] + (1. - theta) * solVOld[k][i]);
              DeltaSolVgTheta[k][j]   += nablaphiV[i * dim2 + j] * (theta * solV[k][i] + (1. - theta) * solVOld[k][i]) ;
            }
          }
        }

        adept::adouble solPg = 0.;
        vector<adept::adouble> gradSolPg(dim, 0.);

        for(unsigned i = 0; i < nDofsP; i++) {
          solPg += phiP[i] * solP[i];
          for(unsigned k = 0; k < dim; k++) {
            gradSolPg[k] += solP[i] * gradPhiP[i * dim + k];
          }
        }


        //start SUPG paramters, tauM, tauC, G to get tauM_SupgPhi
        std::vector <std::vector <double> > JacMatrix;
        msh->_finiteElement[ielt][solType]->GetJacobian(vxHat, ig, weightV, JacMatrix);


        std::vector <std::vector <adept::adouble> > G(dim); // J^T . J
        for(unsigned i = 0; i < dim; i++) {
          G[i].assign(dim, 0.);
          for(unsigned j = 0; j < dim; j++) {
            for(unsigned k = 0; k < dim; k++) {
              G[i][j] += JacMatrix[k][i] * JacMatrix[k][j];
            }
          }
        }

        adept::adouble tauM = 0.;
        double CI = 36.;
        adept::adouble denom = pow(2 * rhoFluid / dt, 2.);
        for(unsigned i = 0; i < dim; i++) {
          for(unsigned j = 0; j < dim; j++) {
            denom += rhoFluid * solVgTheta[i] * G[i][j] * rhoFluid * solVgTheta[j]
                     + CI * muFluid * muFluid * G[i][j] * G[i][j];
          }
        }
        tauM += 1. / sqrt(denom);

        adept::adouble tauMtrG = 0.;
        for(unsigned k = 0; k < dim; k++) {
          tauMtrG += G[k][k];
        }
        tauMtrG *= tauM;
        adept::adouble tauC = 1. / tauMtrG;

        //end SUPG parameters
        //tauM = tauC = 0.;

        std::vector < adept::adouble > tauM_SupgPhi(nDofs, 0.);
        for(unsigned i = 0; i < nDofs; i++) {
          for(unsigned j = 0; j < dim; j++) {
            tauM_SupgPhi[i] += tauM * (rhoFluid * solVgTheta[j] * gradPhiV[i * dim + j]);
          }
        }

        adept::adouble divVgTheta = 0.;
        for(unsigned k = 0; k < dim; k++) {
          divVgTheta += gradSolVgTheta[k][k];
        }

        for(unsigned i = 0; i < nDofs; i++) {
          for(unsigned k = 0; k < dim; k++) {

            adept::adouble wlaplace = 0.;
            adept::adouble SupgLaplace = 0.;
            adept::adouble advection = 0.;



            for(unsigned j = 0; j < dim; j++) {
              wlaplace += muFluid * gradPhiV[i * dim + j] * (gradSolVgTheta[k][j] + gradSolVgTheta[j][k]);

              unsigned kdim;
              if(k == j) kdim = j;
              else if(1 == k + j) kdim = dim;        // xy
              else if(2 == k + j) kdim = dim + 2;    // xz
              else if(3 == k + j) kdim = dim + 1;    // yz
              SupgLaplace += -muFluid * (DeltaSolVgTheta[k][j] + DeltaSolVgTheta[j][kdim]) * tauM_SupgPhi[i];  // SUPG laplace

              advection += rhoFluid * solVgTheta[j] * gradSolVgTheta[k][j] * (phiV[i] + tauM_SupgPhi[i]);

            }

            adept::adouble SupgDiv = tauC * divVgTheta * gradPhiV[i * dim + k];
            adept::adouble SupgPressure = gradSolPg[k] * tauM_SupgPhi[i];

            aResV[k][i] += (rhoFluid * (solVg[k] - solVgOld[k]) / dt * (phiV[i] + tauM_SupgPhi[i])
                            + advection
                            + wlaplace +  SupgLaplace
                            - gradPhiV[i * dim + k] * solPg + SupgPressure
                            + SupgDiv
                           ) * weightV;
          }
        }


        //continuity block
        for(unsigned i = 0; i < nDofsP; i++) {
          for(unsigned k = 0; k < dim; k++) {

            adept::adouble sLaplace = 0.;
            adept::adouble advection = 0.;


            for(unsigned j = 0; j < dim; j++) {
              unsigned kdim;

              if(k == j) kdim = j;
              else if(1 == k + j) kdim = dim;       // xy
              else if(2 == k + j) kdim = dim + 2;   // xz
              else if(3 == k + j) kdim = dim + 1;   // yz

              sLaplace += (- muFluid * (DeltaSolVgTheta[k][j] + DeltaSolVgTheta[j][kdim]));
              advection += rhoFluid * solVgTheta[j]  * gradSolVgTheta[k][j];

            }


            aResP[i] += (phiP[i] * gradSolVgTheta[k][k] +
                         (rhoFluid * (solVg[k] - solVgOld[k]) / dt + advection +
                          sLaplace +  gradSolPg[k]) * tauM * gradPhiP[i * dim + k]) * weightV;

          }
        }

      } // end gauss point loop

      //copy the value of the adept::adoube aRes in double Res and store them in RES
      rhs.resize(nDofsAll);   //resize

      for(int i = 0; i < nDofs; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          rhs[k * nDofs + i] = -aResV[k][i].value();
        }
      }
      for(int i = 0; i < nDofsP; i++) {
        rhs[ dim * nDofs  + i] = -aResP[i].value();
      }

      myRES->add_vector_blocked(rhs, sysDofsAll);


      Jac.resize(nDofsAll * nDofsAll);
      // define the dependent variables

      for(unsigned  k = 0; k < dim; k++) {
        s.dependent(&aResV[k][0], nDofs);
      }
      s.dependent(&aResP[0], nDofsP);

      for(unsigned  k = 0; k < dim; k++) {
        s.independent(&solV[k][0], nDofs);
      }
      s.independent(&solP[0], nDofsP);

      // get the and store jacobian matrix (row-major)
      s.jacobian(&Jac[0], true);
      myKK->add_matrix_blocked(Jac, sysDofsAll, sysDofsAll);

      s.clear_independents();
      s.clear_dependents();

    }





  }

// *************************************
  end_time = clock();
  AssemblyTime += (end_time - start_time);
// ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

}



