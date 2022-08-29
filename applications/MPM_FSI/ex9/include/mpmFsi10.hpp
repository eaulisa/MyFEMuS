#include "MultiLevelSolution.hpp"
#include "PetscMatrix.hpp"

using namespace femus;


//void GetParticlesToNodeFlag(MultiLevelSolution &mlSol, Line & solidLine, Line & fluidLine);
//void GetPressureNeighbor(MultiLevelSolution &mlSol, Line & solidLine, Line & fluidLine);

void ProjectGridVelocity(MultiLevelSolution &mlSol);

double GetSmoothStepFunction(const double &dg1, const double &eps) {

  double a0 = 0.5; // 128./256.;
  double a1 = pow(eps, -1.) * 1.23046875; // 315/256.;
  double a3 = -pow(eps, -3.) * 1.640625; //420./256.;
  double a5 = pow(eps, -5.) * 1.4765625; // 378./256.;
  double a7 = -pow(eps, -7.) * 0.703125; // 180./256.;
  double a9 = pow(eps, -9.) * 0.13671875; // 35./256.;

  double dg2 = dg1 * dg1;
  double xi;
  if(dg1 < -eps)
    xi = 0.;
  else if(dg1 > eps) {
    xi = 1.;
  }
  else {
    xi = (a0 + dg1 * (a1 + dg2 * (a3 + dg2 * (a5 + dg2 * (a7 + dg2 * a9)))));
  }

  return xi;
}


double clamp(double x, double lowerlimit, double upperlimit) {
  if(x < lowerlimit)
    x = lowerlimit;
  if(x > upperlimit)
    x = upperlimit;
  return x;
}

double smoothstep(double edge0, double edge1, double x) {
  // Scale, bias and saturate x to 0..1 range
  x = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
  // Evaluate polynomial
  return x * x * (3 - 2 * x);
}



void AssembleGhostPenalty(MultiLevelProblem& ml_prob, const bool &fluid) {

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

  MatSetOption((static_cast< PetscMatrix* >(myKK))->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);


  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  const unsigned dim = msh->GetDimension();
  const unsigned dim2 = 3 * (dim - 1);

  // data
  unsigned iproc  = msh->processor_id();

  //quantities for iel will have index1
  //quantities for jel will have index2


  vector< vector< double > > solV1Old(dim);
  vector< vector< double > > solV2Old(dim);

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

      double h = sqrt((vx1[0][0] - vx1[0][2]) * (vx1[0][0] - vx1[0][2]) +
                      (vx1[1][0] - vx1[1][2]) * (vx1[1][0] - vx1[1][2])) ;
      double h2 = h * h;
      double h3 = h * h * h;

      bool aP1IsInitialized = false;

      for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
        int jel = el->GetFaceElementIndex(iel, iface) - 1;

        unsigned eFlag2 = (jel >= 0) ? static_cast <unsigned>(floor((*mysolution->_Sol[eflagIndex])(jel) + 0.5)) : 3;

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

// *************************************
  end_time = clock();
  AssemblyTime += (end_time - start_time);


}





void AssembleMPMSys(MultiLevelProblem& ml_prob) {

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

  myKK->zero();
  myRES->zero();

  AssembleGhostPenalty(ml_prob, true);
  AssembleGhostPenalty(ml_prob, false);


  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  const unsigned dim = msh->GetDimension();

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  // data
  unsigned iproc  = msh->processor_id();

  vector< vector< adept::adouble > > solD(dim);      // local solution (displacement)
  vector< vector< adept::adouble > > solV(dim);      // local solution (velocity)
  vector< adept::adouble > solP;

  vector< vector< double > > solDOld(dim);      // local solution (displacement)
  vector< vector< double > > solVOld(dim);

  vector< double > rhs;    // local redidual vector
  vector< vector< adept::adouble > > aResD(dim);     // local redidual vector
  vector< vector< adept::adouble > > aResV(dim);     // local redidual vector
  vector< adept::adouble > aResP;    // local redidual vector
  vector < double > Jac;

  std::vector <unsigned> sysDofsAll;
  vector < double > phi;
  vector < double > phiHat;
  vector < double > phiP;
  vector < adept::adouble> gradPhiP;  // phiP_x

  vector < adept::adouble> gradPhi;  // phi_x
  vector < adept::adouble> nablaphi; //phi_xx
  vector < double > gradPhiHat;



  unsigned dim2 = 3 * (dim - 1);




  vector <vector < adept::adouble> > vx(dim);
  vector <vector < double> > vxHat(dim);

  adept::adouble weight;
  double weightHat;


  //reading parameters for MPM body
  double rhoMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_density();
  double EMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_young_module();
  double muMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_lame_shear_modulus();
  double nuMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_poisson_coeff();
  double lambdaMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_lame_lambda();
  double KMmp = EMpm / (3.* (1. - 2. * nuMpm));     //bulk modulus

  //reading parameters for fluid FEM domain
  double rhoFluid = ml_prob.parameters.get<Fluid> ("FluidFEM").get_density();
  double muFluid = ml_prob.parameters.get<Fluid> ("FluidFEM").get_viscosity();

  double dt =  my_nnlin_impl_sys.GetIntervalTime();

  std::cout.precision(10);

  //variable-name handling
  const char varname[10][5] = {"DX", "DY", "DZ", "VX", "VY", "VZ"};

  vector <unsigned> indexSolD(dim);
  vector <unsigned> indexSolV(dim);
  vector <unsigned> indexPdeD(dim);
  vector <unsigned> indexPdeV(dim);
  for(unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolD[ivar] = mlSol->GetIndex(&varname[ivar][0]);
    indexSolV[ivar] = mlSol->GetIndex(&varname[ivar + 3][0]);
    indexPdeD[ivar] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar][0]);       // DX, DY, DZ
    indexPdeV[ivar] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar + 3][0]);   //VX, VY, VZ
  }
  unsigned solType = mlSol->GetSolutionType(&varname[0][0]);

  unsigned indexSolP = mlSol->GetIndex("P");
  unsigned indexPdeP = my_nnlin_impl_sys.GetSolPdeIndex("P");
  unsigned solTypeP = mlSol->GetSolutionType("P");

  unsigned eflagIndex = mlSol->GetIndex("eflag");
  unsigned nflagIndex = mlSol->GetIndex("nflag");
  std::vector < unsigned >  nodeFlag; // local solution

  start_time = clock();


  std::vector<Marker*> particlesBulk = bulk->GetParticles();
  std::vector<unsigned> markerOffsetBulk = bulk->GetMarkerOffset();
  unsigned iBmarker = markerOffsetBulk[iproc];


  std::vector<Marker*> particleI = lineI->GetParticles();
  std::vector<unsigned> markerOffsetI = lineI->GetMarkerOffset();
  unsigned imarkerI = markerOffsetI[iproc];



  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielt = msh->GetElementType(iel);

    unsigned nDofs = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solTypeP);  // number of pressure dofs

    unsigned nDofsAll = 2 * dim * nDofs + nDofsP;


    // resize local arrays
    sysDofsAll.resize(nDofsAll);

    nodeFlag.resize(nDofs);
    unsigned eFlag = static_cast <unsigned>(floor((*mysolution->_Sol[eflagIndex])(iel) + 0.5));

    for(unsigned  k = 0; k < dim; k++) {
      solD[k].resize(nDofs);
      solDOld[k].resize(nDofs);

      solV[k].resize(nDofs);
      solVOld[k].resize(nDofs);

      aResD[k].assign(nDofs, 0.);
      aResV[k].assign(nDofs, 0.);

      vx[k].resize(nDofs);
      vxHat[k].resize(nDofs);
    }
    solP.resize(nDofsP);
    aResP.assign(nDofsP, 0.);


    for(unsigned i = 0; i < nDofs; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, solType);

      nodeFlag[i] = (*mysolution->_Sol[nflagIndex])(idof);

      for(unsigned  k = 0; k < dim; k++) {
        solD[k][i] = (*mysolution->_Sol[indexSolD[k]])(idof);
        solDOld[k][i] = (*mysolution->_SolOld[indexSolD[k]])(idof);

        solV[k][i] = (*mysolution->_Sol[indexSolV[k]])(idof);
        solVOld[k][i] = (*mysolution->_SolOld[indexSolV[k]])(idof);

        sysDofsAll[i + k * nDofs] = myLinEqSolver->GetSystemDof(indexSolD[k], indexPdeD[k], i, iel);
        sysDofsAll[i + (k + dim) * nDofs] = myLinEqSolver->GetSystemDof(indexSolV[k], indexPdeV[k], i, iel);
      }
    }

    for(unsigned i = 0; i < nDofsP; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, solTypeP);
      solP[i] = (*mysolution->_Sol[indexSolP])(idof);
      sysDofsAll[i + (2 * dim) * nDofs] = myLinEqSolver->GetSystemDof(indexSolP, indexPdeP, i, iel);
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();


    for(unsigned i = 0; i < nDofs; i++) {
      unsigned idofX = msh->GetSolutionDof(i, iel, 2);
      for(unsigned  k = 0; k < dim; k++) {
        vxHat[k][i] = (*msh->_topology->_Sol[k])(idofX) + solDOld[k][i]; // reference configuration
        //vx[k][i] = vxHat[k][i] + solD[k][i]; // current configuration
        //vx[k][i] = (*msh->_topology->_Sol[k])(idofX) +  af * solD[k][i] + (1. - af) * solDOld[k][i]; // af configuration
        //vx[k][i] = vxHat[k][i] +  af * (solD[k][i] - solDOld[k][i]); // af configuration
        vx[k][i]  = vxHat[k][i] + (1. - af) * (solD[k][i] - solDOld[k][i]);
      }
    }

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielt][solType]->GetGaussPointNumber(); ig++) {

      msh->_finiteElement[ielt][solTypeP]->Jacobian(vx, ig, weight, phiP, gradPhiP);
      msh->_finiteElement[ielt][solType]->Jacobian(vx, ig, weight, phi, gradPhi, nablaphi);

      msh->_finiteElement[ielt][solType]->Jacobian(vxHat, ig, weightHat, phiHat, gradPhiHat);



      vector < adept::adouble > solVg(dim, 0.);
      vector < adept::adouble > solVgOld(dim, 0.);
      vector < adept::adouble > solVgTheta(dim, 0.);


      vector < adept::adouble > solDg(dim, 0.);
      vector < adept::adouble > solDgOld(dim, 0.);

      vector < vector < adept::adouble > > gradSolVgTheta(dim);
      vector<vector<adept::adouble> > DeltaSolVgTheta(dim);


      vector < vector < adept::adouble > > gradSolDgHat(dim);
      vector < vector < adept::adouble > > gradSolVg(dim);
      vector<vector<adept::adouble> > DeltaSolVg(dim); // DeltaSol = [ [uh0_xx, uh0_yy, uh0_xy], [uh1_xx, uh1_yy, uh1_xy] ]



      for(unsigned  k = 0; k < dim; k++) {
        gradSolDgHat[k].assign(dim, 0);
        gradSolVg[k].assign(dim, 0);
        DeltaSolVg[k].resize(dim2, 0.);
        gradSolVgTheta[k].assign(dim, 0.);
        DeltaSolVgTheta[k].resize(dim2, 0.);
      }

      for(unsigned i = 0; i < nDofs; i++) {
        for(unsigned j = 0; j < dim; j++) {
          solVg[j] += phi[i] * solV[j][i];
          solVgOld[j] += phi[i] * solVOld[j][i];
          solVgTheta[j] += phi[i] * (theta * solV[j][i] + (1. - theta) * solVOld[j][i]);

          solDgOld[j] += phi[i] * solDOld[j][i];
          solDg[j] += phi[i] * solD[j][i];
          for(unsigned  k = 0; k < dim; k++) {
            gradSolDgHat[k][j] += gradPhiHat[i * dim + j] * solD[k][i];
            gradSolVg[k][j] += gradPhi[i * dim + j] * solV[k][i];
            gradSolVgTheta[k][j] += gradPhi[i * dim + j] * (theta * solV[k][i] + (1. - theta) * solVOld[k][i]);
          }
        }

        for(unsigned j = 0; j < dim2; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            DeltaSolVgTheta[k][j]   += nablaphi[i * dim2 + j] * (theta * solV[k][i] + (1. - theta) * solVOld[k][i]) ;
          }
        }
      }


      adept::adouble solPg = 0.;
      for(unsigned i = 0; i < nDofsP; i++) {
        solPg += phiP[i] * solP[i];
      }





      if(eFlag == 0) { //bulk fluid: smooth extension nabla D = 0

        adept::adouble F[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned k = 0; k < dim; k++) {
            F[j][k] += gradSolDgHat[j][k];
          }
        }

        adept::adouble J_hat =  F[0][0] * F[1][1] * F[2][2] + F[0][1] * F[1][2] * F[2][0] + F[0][2] * F[1][0] * F[2][1]
                                - F[2][0] * F[1][1] * F[0][2] - F[2][1] * F[1][2] * F[0][0] - F[2][2] * F[1][0] * F[0][1];


        adept::adouble B[3][3];
        for(unsigned i = 0; i < 3; i++) {
          for(int j = 0; j < 3; j++) {
            B[i][j] = 0.;
            for(unsigned k = 0; k < 3; k++) {
              //left Cauchy-Green deformation tensor or Finger tensor (B = F*F^T)
              B[i][j] += F[i][k] * F[j][k];
            }
          }
        }

        adept::adouble I1_B = B[0][0] + B[1][1] + B[2][2];
        adept::adouble Id2th[3][3] = {{ 1., 0., 0.}, { 0., 1., 0.}, { 0., 0., 1.}};
        adept::adouble sigma[3][3];

        double E = 1;
        double nu = 0.4;

        double mu = E / (2. * (1. + nu));
        double lambda = (E * nu) / ((1. + nu) * (1. - 2.*nu));

        for(unsigned j = 0; j < 3; j++) {
          for(unsigned k = 0; k < 3; k++) {
            sigma[j][k] = lambda * log(J_hat) / J_hat * Id2th[j][k] + mu / J_hat * (B[j][k] - Id2th[j][k]);    // alternative formulation
          }
        }
        //END computation of the Cauchy Stress
        for(unsigned i = 0; i < nDofs; i++) {//Auxiliary Equations
          if(nodeFlag[i] == 0) {
            for(unsigned k = 0.; k < dim; k++) {
              adept::adouble cauchy = 0.;
              for(unsigned j = 0.; j < dim; j++) {
                cauchy += sigma[k][j] * gradPhi[i * dim + j] ;
              }
              aResD[k][i] += cauchy * weight;
            }
          }
        }
      }


//               adept::adouble wlaplaceD  = 0.;
//               for(unsigned  j = 0; j < dim; j++) {
//                 wlaplaceD +=  gradPhiHat[i * dim + j] * (gradSolDgHat[k][j] + gradSolDgHat[j][k]);
//               }
//               aResD[k][i] +=  wlaplaceD * weightHat;



      if(eFlag == 0) {   // only fluid cells

        //start SUPG paramters, tauM, tauC, G to get tauM_SupgPhi
        std::vector <std::vector <adept::adouble> > JacI;
        std::vector <std::vector <adept::adouble> > Jac;
        std::vector <std::vector <double> > JacIHat;
        msh->_finiteElement[ielt][solType]->GetJacobianMatrix(vx, ig, weight, Jac, JacI);


        if(solTypeP == 4) { //discontinuous pressure <1,\xi,\eta> bases
          for(unsigned j = 0; j < dim; j++) {
            gradPhiP[0 * dim + j]  = 0.;
          }
          for(unsigned i = 0; i < dim; i++) {
            for(unsigned j = 0; j < dim; j++) {
              gradPhiP[(i + 1) * dim + j] = JacI[i][j];
            }
          }
        }

        vector<adept::adouble> gradSolPg(dim, 0.);
        for(unsigned i = 0; i < nDofsP; i++) {
          for(unsigned k = 0; k < dim; k++) {
            gradSolPg[k] += solP[i] * gradPhiP[i * dim + k];
          }
        }



        std::vector <std::vector <adept::adouble> > G(dim); // J^T . J
        for(unsigned i = 0; i < dim; i++) {
          G[i].assign(dim, 0.);
          for(unsigned j = 0; j < dim; j++) {
            for(unsigned k = 0; k < dim; k++) {
              G[i][j] += JacI[k][i] * JacI[k][j];
            }
          }
        }

        adept::adouble tauM = 0.;
        double CI = 36.;
        adept::adouble denom = pow(2 * rhoFluid / dt, 2.);
        for(unsigned i = 0; i < dim; i++) {
          for(unsigned j = 0; j < dim; j++) {
            denom += rhoFluid * (solVgTheta[i] - (solDg[i] - solDgOld[i]) / dt) * G[i][j] * rhoFluid * (solVgTheta[j] - (solDg[j] - solDgOld[j]) / dt) // we can get the mesh velocity at af
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

        std::vector < adept::adouble > tauM_SupgPhi(nDofs, 0.);
        for(unsigned i = 0; i < nDofs; i++) {
          for(unsigned j = 0; j < dim; j++) {
            tauM_SupgPhi[i] += tauM * (rhoFluid * (solVgTheta[j] - (solDg[j] - solDgOld[j]) / dt) * gradPhi[i * dim + j]);
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
              wlaplace += muFluid * gradPhi[i * dim + j] * (gradSolVgTheta[k][j] + gradSolVgTheta[j][k]);

              unsigned kdim;
              if(k == j) kdim = j;
              else if(1 == k + j) kdim = dim;       // xy
              else if(2 == k + j) kdim = dim + 2;    // xz
              else if(3 == k + j) kdim = dim + 1;    // yz
              SupgLaplace += (-muFluid * (DeltaSolVgTheta[k][j] + DeltaSolVgTheta[j][kdim])) * tauM_SupgPhi[i];  // SUPG laplace

              advection  +=  rhoFluid * (solVgTheta[j] - (solDg[j] - solDgOld[j]) / dt) * gradSolVgTheta[k][j] * (phi[i] + tauM_SupgPhi[i]);  // ALE advection + SUPG Advection

            }

            adept::adouble SupgDiv = tauC * divVgTheta * gradPhi[i * dim + k];
            adept::adouble SupgPressure = gradSolPg[k] * tauM_SupgPhi[i];


            aResV[k][i] += (rhoFluid * (solVg[k] - solVgOld[k]) / dt * (phi[i] + tauM_SupgPhi[i])
                            + advection
                            + wlaplace + SupgLaplace
                            - gradPhi[i * dim + k] * solPg
                            //+ phi[i] * gradSolPg[k]
                            + SupgPressure
                            + SupgDiv
                           ) * weight;
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
              advection += rhoFluid * (solVgTheta[j] - (solDg[j] - solDgOld[j]) / dt) * gradSolVgTheta[k][j];

            }

            aResP[i] += (phiP[i] * gradSolVgTheta[k][k] +
                         (rhoFluid * (solVg[k] - solVgOld[k]) / dt + advection +
                          sLaplace +  gradSolPg[k]) * tauM * gradPhiP[i * dim + k]
                        ) * weight;

          }


        }

      }



      else if(eFlag == 2) {   // only solid cells: fake pressure //TODO
        for(unsigned i = 0; i < nDofsP; i++) {
          aResP[i] += 1.0e-10 * phiP[i] * solPg * weight;
        }
      }
    } // end gauss point loop




    //BEGIN BULK PARTICLE
    if(eFlag > 0) {   // interface and solid
      while(iBmarker < markerOffsetBulk[iproc + 1] && iel > particlesBulk[iBmarker]->GetMarkerElement()) {
        iBmarker++;
      }
      while(iBmarker < markerOffsetBulk[iproc + 1] && iel == particlesBulk[iBmarker]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particlesBulk[iBmarker]->GetMarkerLocalCoordinates();
        double dg1 = particlesBulk[iBmarker]->GetMarkerDistance();

        double U = GetSmoothStepFunction(dg1, eps);

        double area = particlesBulk[iBmarker]->GetMarkerMass();

        msh->_finiteElement[ielt][solType]->Jacobian(vx, xi, weight, phi, gradPhi);
        msh->_finiteElement[ielt][solTypeP]->Jacobian(vx, xi, weight, phiP, gradPhiP);
        msh->_finiteElement[ielt][solType]->Jacobian(vxHat, xi, weightHat, phiHat, gradPhiHat);


        // BEGIN EVALUATION Quantities at the particles
        std::vector <double> solVpOld(dim);
        particlesBulk[iBmarker]->GetMarkerVelocity(solVpOld);

        std::vector <double> solApOld(dim);
        particlesBulk[iBmarker]->GetMarkerAcceleration(solApOld);

        vector<adept::adouble> solDp(dim, 0.);
        vector<double> solDpOld(dim, 0.);
        vector<adept::adouble> solVp(dim, 0.);
        vector<adept::adouble> solVpTheta(dim, 0.);
        //vector<adept::adouble> solVgOldp(dim, 0.);

        //vector<vector < adept::adouble > > gradSolVp(dim);
        vector<vector < adept::adouble > > gradSolVpTheta(dim);
        //vector<vector < adept::adouble > > gradSolDp(dim);
        vector<vector < adept::adouble > > gradSolDpHat(dim);

        for(int j = 0; j < dim; j++) {
          //gradSolVp[j].assign(dim, 0.);
          gradSolVpTheta[j].assign(dim, 0.);
          //gradSolDp[j].assign(dim, 0.);
          gradSolDpHat[j].assign(dim, 0.);
        }

        for(int j = 0; j < dim; j++) {
          for(unsigned i = 0; i < nDofs; i++) {
            solDp[j] += phi[i] * solD[j][i];
            solDpOld[j] += phi[i] * solDOld[j][i];
            solVp[j] += phi[i] * solV[j][i]; //
            //solVgOldp[j] += phi[i] * solVOld[j][i]; //TODO
            for(int k = 0; k < dim; k++) {
              //gradSolDp[j][k] +=  gradPhi[i * dim + k] * solD[j][i];
              //gradSolVp[j][k] +=  gradPhi[i * dim + k] * solV[j][i];
              gradSolVpTheta[j][k] +=  gradPhi[i * dim + k] * (theta * solV[j][i] + (1. - theta) * solVOld[j][i]);
              //gradSolDpHat[j][k] +=  gradPhiHat[i * dim + k] * solD[j][i];
              gradSolDpHat[k][j] += ((1. - af) * solD[k][i] + af * solDOld[k][i]) * gradPhiHat[i * dim + j];
            }
          }
        }

        std::vector<adept::adouble> solAp(dim);
        std::vector < adept::adouble > solApAm(dim);
        std::vector < adept::adouble > solVpS(dim);

        for(int j = 0; j < dim; j++) {
          solVpTheta[j] = theta * solVp[j] + (1. - theta) * solVpOld[j]; //TODO this comes from the particle, try the old mesh velocity
          solAp[j] = (solDp[j] - solDpOld[j]) / (beta * dt * dt) - solVpOld[j] / (beta * dt) + solApOld[j] * (beta - 0.5) / beta ;   //NEWMARK ACCELERATION
          solApAm[j] = (1. - am) * solAp[j] + am * solApOld[j];
          solVpS[j] = solVpOld[j] + dt * (Gamma * solAp[j] + (1. - Gamma) * solApOld[j]);   //velocity from the solid at xp, gamma configuration
        }

//         adept::adouble solPp = 0.;
//         msh->_finiteElement[ielt][solTypeP]->GetPhi(phiP, xi);
//         for(unsigned i = 0; i < nDofsP; i++) {
//           solPp += phiP[i] * solP[i];
//         }

        adept::adouble solPg = 0.;
        vector<adept::adouble> gradSolPg(dim, 0.);
        for(unsigned i = 0; i < nDofsP; i++) {
          solPg += phiP[i] * solP[i];
          for(unsigned k = 0; k < dim; k++) {
            gradSolPg[k] += solP[i] * gradPhiP[i * dim + k];
          }
        }




        //BEGIN computation of the Cauchy Stress
        std::vector < std::vector < double > > FpOld;
        FpOld = particlesBulk[iBmarker]->GetDeformationGradient(); //extraction of the deformation gradient

        adept::adouble FpNew[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
        adept::adouble F[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
        adept::adouble B[3][3];
        adept::adouble Id2th[3][3] = {{ 1., 0., 0.}, { 0., 1., 0.}, { 0., 0., 1.}};
        adept::adouble Cauchy[3][3];

        for(unsigned j = 0; j < dim; j++) {
          for(unsigned k = 0; k < dim; k++) {
            FpNew[j][k] += gradSolDpHat[j][k];
          }
        }

        for(unsigned i = 0; i < dim; i++) {
          for(unsigned j = 0; j < dim; j++) {
            for(unsigned k = 0; k < dim; k++) {
              F[i][j] += FpNew[i][k] * FpOld[k][j];
            }
          }
        }

        if(dim == 2) F[2][2] = 1.;

        adept::adouble J_hat =  F[0][0] * F[1][1] * F[2][2] + F[0][1] * F[1][2] * F[2][0] + F[0][2] * F[1][0] * F[2][1]
                                - F[2][0] * F[1][1] * F[0][2] - F[2][1] * F[1][2] * F[0][0] - F[2][2] * F[1][0] * F[0][1];

        for(unsigned i = 0; i < 3; i++) {
          for(int j = 0; j < 3; j++) {
            B[i][j] = 0.;
            for(unsigned k = 0; k < 3; k++) {
              //left Cauchy-Green deformation tensor or Finger tensor (B = F*F^T)
              B[i][j] += F[i][k] * F[j][k];
            }
          }
        }

        adept::adouble I1_B = B[0][0] + B[1][1] + B[2][2];

        for(unsigned j = 0; j < 3; j++) {
          for(unsigned k = 0; k < 3; k++) {
            Cauchy[j][k] = lambdaMpm * log(J_hat) / J_hat * Id2th[j][k] + muMpm / J_hat * (B[j][k] - Id2th[j][k]);     //alternative formulation
          }
        }
        //END computation of the Cauchy Stress



        //BEGIN Navier-Stokes in the bulk interface cells (integration is on the particles in \Omega_f)
        if((1. - U) > 0 && eFlag == 1) {

          double dM = (1 - U) * area * rhoFluid;
          for(unsigned i = 0; i < nDofs; i++) {
            for(unsigned k = 0; k < dim; k++) {
              adept::adouble Vlaplace = 0.;
              adept::adouble advection = 0.;
              for(unsigned j = 0; j < dim; j++) {
                Vlaplace  +=  gradPhi[i * dim + j] * (gradSolVpTheta[k][j] + gradSolVpTheta[j][k]);
                advection +=  phi[i] * (solVpTheta[j] - (solDp[j] - solDpOld[j]) / dt) * gradSolVpTheta[k][j]; //ALE
              }

              aResV[k][i] += (phi[i] * (solVp[k] - solVpOld[k]) / dt + advection +
                              muFluid / rhoFluid * Vlaplace
                              - gradPhi[i * dim + k] * solPg / rhoFluid
//                               + phi[i] * gradSolPg[k] / rhoFluid
                             ) * dM;
            }
          }

          for(unsigned i = 0; i < nDofsP; i++) {
            if(eFlag == 1) {
              for(unsigned  k = 0; k < dim; k++) {
                aResP[i] += phiP[i] *  gradSolVpTheta[k][k] * (1 - U) * area;
              }
            }
          }
        }

        //BEGIN solid Momentum in the bulk interface and solid cells (integration is on the particles in \Omega_s)
        if(U > 0) {
          double dM = (eFlag == 1) ? U * area * rhoMpm : area * rhoMpm;
          for(unsigned i = 0; i < nDofs; i++) {
            adept::adouble CauchyDIR[3] = {0., 0., 0.};
            for(unsigned j = 0.; j < dim; j++) {
              for(unsigned k = 0.; k < dim; k++) {
                CauchyDIR[j] += gradPhi[i * dim + k] * Cauchy[j][k];
              }
            }



            for(unsigned k = 0; k < dim; k++) {

              aResD[k][i] += (phi[i] * solApAm[k] + J_hat * CauchyDIR[k] / rhoMpm)  * dM;

              if(nodeFlag[i] == 0) { //bulk solid nodes: kinematic: v - dD/dt = 0

                aResV[k][i] += -phiHat[i] * (solVp[k] - solVpS[k]) * area;
              }

            }
          }
        }
        iBmarker++;
      }
    }
    //END BULK PARTICLE

    if(true) {

      //BEGIN INTERFACE PARTICLE

      if(eFlag == 1) {  //interface markers

        double h = sqrt((vxHat[0][0] - vxHat[0][2]) * (vxHat[0][0] - vxHat[0][2]) +
                        (vxHat[1][0] - vxHat[1][2]) * (vxHat[1][0] - vxHat[1][2])) ;


        double thetaM = GAMMA * muFluid / h;
        double thetaL = GAMMA * rhoFluid / (theta * dt) * h  + thetaM;

        while(imarkerI < markerOffsetI[iproc + 1] && iel != particleI[imarkerI]->GetMarkerElement()) {
          imarkerI++;
        }

        while(imarkerI < markerOffsetI[iproc + 1] && iel == particleI[imarkerI]->GetMarkerElement()) {

          // the local coordinates of the particles are the Gauss points in this context
          std::vector <double> xi = particleI[imarkerI]->GetMarkerLocalCoordinates();
          msh->_finiteElement[ielt][solType]->Jacobian(vxHat, xi, weightHat, phiHat, gradPhiHat);
          msh->_finiteElement[ielt][solType]->Jacobian(vx, xi, weight, phi, gradPhi);

          std::vector <std::vector < double > > THat;
          particleI[imarkerI]->GetMarkerTangent(THat);

          std::vector < std::vector < double > > FpOld;
          FpOld = particleI[imarkerI]->GetDeformationGradient(); //extraction of the deformation gradient

          std::vector <std::vector < double > > T;
          T.resize(THat.size());

          for(unsigned k = 0; k < T.size(); k++) {
            T[k].assign(dim, 0.);
            for(unsigned i = 0; i < dim; i++) {
              for(unsigned j = 0; j < dim; j++) {
                T[k][i] += FpOld[i][j] * THat[k][j]; // can be improved
              }
            }
          }

          double weight;
          std::vector < double > N(dim);
          if(dim == 2) {
            N[0] =  T[0][1];
            N[1] = -T[0][0];
            weight = sqrt(N[0] * N[0] + N[1] * N[1]);
            N[0] /= weight;
            N[1] /= weight;
          }
          else {
            N[0] = T[0][1] * T[1][2] - T[0][2] * T[1][1];
            N[1] = T[0][2] * T[1][0] - T[0][0] * T[1][2];
            N[2] = T[0][0] * T[1][1] - T[0][1] * T[1][0];
            weight = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
            N[0] /= weight;
            N[1] /= weight;
            N[2] /= weight;
          }

          for(unsigned k = 0; k < dim; k++) {
            N[k] *= -1.; // fluid to solid normal
          }

          msh->_finiteElement[ielt][solTypeP]->GetPhi(phiP, xi);
          adept::adouble solPp = 0.;
          for(unsigned i = 0; i < nDofsP; i++) {
            solPp += phiP[i] * solP[i];
          }

          std::vector < adept::adouble > v1(dim, 0.); //fluid velocity in theta configuration
          std::vector < adept::adouble > v2(dim, 0.); //solid velocity in af configuration

          std::vector < adept::adouble > tau(dim, 0.);

          std::vector <double> solVpOld(dim);
          particleI[imarkerI]->GetMarkerVelocity(solVpOld);

          std::vector <double> solApOld(dim);
          particleI[imarkerI]->GetMarkerAcceleration(solApOld);

          std::vector <adept::adouble> solDp(dim, 0.);
          std::vector <double> solDpOld(dim, 0.);
          std::vector <adept::adouble> solAp(dim);
          //std::vector < adept::adouble > solApAm(dim);

          //update displacement and acceleration
          for(int k = 0; k < dim; k++) {
            for(unsigned i = 0; i < nDofs; i++) {
              v1[k] += phi[i] * (theta * solV[k][i] + (1. - theta) * solVOld[k][i]);
              solDp[k] += phi[i] * solD[k][i];
              solDpOld[k] += phi[i] * solDOld[k][i];
            }
            solAp[k] = (solDp[k] - solDpOld[k]) / (beta * dt * dt) - solVpOld[k] / (beta * dt) + (beta - 0.5) / beta * solApOld[k]; // Newmark acceleration
            //solApAm[k] = (1. - am) * solAp[k] + am * solApOld[k]; // generalized alpha acceleration
            v2[k] = solVpOld[k] + (1. - af) * (dt * (Gamma * solAp[k] + (1. - Gamma) * solApOld[k]));

          }

          for(unsigned k = 0; k < dim; k++) {
            tau[k] += solPp * N[k];
            for(unsigned i = 0; i < nDofs; i++) {
              for(unsigned j = 0; j < dim; j++) {
                //tau[k] += - muFluid * (solV[k][i] * gradPhi[i * dim + j] + solV[j][i] * gradPhi[i * dim + k]) * N[j];
                tau[k] += -muFluid * ((theta * solV[k][i] + (1. - theta) * solVOld[k][i]) * gradPhi[i * dim + j] +
                                      (theta * solV[j][i] + (1. - theta) * solVOld[j][i]) * gradPhi[i * dim + k]) * N[j];
              }
            }
          }

          // *** phi_i loop ***
          for(unsigned i = 0; i < nDofs; i++) {
            for(unsigned k = 0; k < dim; k++) {
              aResV[k][i] += (tau[k] - 0. * solPp * N[k]) * phi[i] * weight;  // correct sign due to the normal
              aResV[k][i] += thetaM * (v1[k] - v2[k]) * phi[i] * weight;

              aResD[k][i] += -tau[k] * phi[i] * weight; // correct sign due to the normal
              aResD[k][i] +=  thetaM * (v1[k] - v2[k]) * (-phi[i]) * weight;

              for(unsigned j = 0; j < dim; j++) {
                aResV[k][i] += - muFluid * gradPhi[i * dim + j] * N[j] * (v1[k] - v2[k]) * weight;
                aResV[k][i] += - muFluid * gradPhi[i * dim + j] * N[k] * (v1[j] - v2[j]) * weight;

                aResV[k][i] += thetaL * (v1[j] - v2[j]) * N[j] * phi[i] * N[k] * weight;
                aResD[k][i] += thetaL * (v1[j] - v2[j]) * N[j] * (-phi[i]) * N[k] * weight;

              }
            }
          } // end phi_i loop


          for(unsigned i = 0; i < nDofsP; i++) {
            for(unsigned k = 0; k < dim; k++) {
              aResP[i] += - (phiP[i]  * (v1[k] - v2[k]) * N[k]) * weight;
            }
          }
          imarkerI++;
        }
      }
    }
    //END INTERFACE PARTICLES

    //copy the value of the adept::adoube aRes in double Res and store them in RES
    rhs.resize(nDofsAll);   //resize

    for(int i = 0; i < nDofs; i++) {
      for(unsigned  k = 0; k < dim; k++) {
        rhs[ i +  k * nDofs ] = -aResD[k][i].value();
        rhs[ i + (k + dim) * nDofs ] = -aResV[k][i].value();
      }
    }
    for(int i = 0; i < nDofsP; i++) {
      rhs[ i + (2 * dim) * nDofs ] = -aResP[i].value();
    }

    myRES->add_vector_blocked(rhs, sysDofsAll);


    Jac.resize(nDofsAll * nDofsAll);
    // define the dependent variables

    for(unsigned  k = 0; k < dim; k++) {
      s.dependent(&aResD[k][0], nDofs);
    }
    for(unsigned  k = 0; k < dim; k++) {
      s.dependent(&aResV[k][0], nDofs);
    }
    s.dependent(&aResP[0], nDofsP);

    // define the independent variables
    for(unsigned  k = 0; k < dim; k++) {
      s.independent(&solD[k][0], nDofs);
    }
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



  myRES->close();
  myKK->close();


//   PetscViewer    viewer1;
//   PetscViewerDrawOpen (PETSC_COMM_WORLD, NULL, NULL, 0, 0, 1800, 1800, &viewer1);
//   PetscObjectSetName ( (PetscObject) viewer1, "FSI matrix");
//   PetscViewerPushFormat (viewer1, PETSC_VIEWER_DRAW_LG);
//   MatView ( (static_cast< PetscMatrix* > (myKK))->mat(), viewer1);
//
//   double a;
//   std::cin >> a;

// *************************************
  end_time = clock();
  AssemblyTime += (end_time - start_time);
// ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

}





void GridToParticlesProjection(MultiLevelProblem & ml_prob,
//                                Line & solidLine, Line & fluidLine, Line & interfaceLine,
                               Line & bulk, Line & lineI) {

  //pointers and references

  TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("MPM_FSI");

  const unsigned  level = my_nnlin_impl_sys.GetLevelToAssemble();
  MultiLevelSolution* mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* mysolution = mlSol->GetSolutionLevel(level);     // pointer to the solution (level) object

  Mesh* msh = ml_prob._ml_msh->GetLevel(level);     // pointer to the mesh (level) object
  elem* el = msh->el;   // pointer to the elem object in msh (level)

  double dt =  my_nnlin_impl_sys.GetIntervalTime();
  const unsigned dim = msh->GetDimension();

  // data
  unsigned iproc  = msh->processor_id();

  // local objects
  vector< vector < double > > solD(dim);
  vector< vector < double > > solDOld(dim);
  vector< vector < double > > gradSolDHat(dim);

  for(int k = 0; k < dim; k++) {
    gradSolDHat[k].resize(dim);
  }

  vector < double > phiHat;
  vector < double > gradPhiHat;

  vector <vector < double> > vxHat(dim);   //vx is coordX in assembly of ex30

  double weightHat;

  //variable-name handling
  const char varname[9][3] = {"DX", "DY", "DZ"};
  vector <unsigned> indexSolD(dim);
  unsigned solType = mlSol->GetSolutionType(&varname[0][0]);

  for(unsigned k = 0; k < dim; k++) {
    indexSolD[k] = mlSol->GetIndex(&varname[k][0]);
  }

  //BEGIN loop on bulk particles

  std::vector<Marker*> particles =  bulk.GetParticles();
  std::vector<unsigned> markerOffset = bulk.GetMarkerOffset();

  unsigned ielOld = UINT_MAX;
  for(unsigned iMarker = markerOffset[iproc]; iMarker < markerOffset[iproc + 1]; iMarker++) {
    unsigned iel = particles[iMarker]->GetMarkerElement();
    if(iel != UINT_MAX) {
      short unsigned ielt;
      unsigned nDofs;
      //update element related quantities only if we are in a different element
      if(iel != ielOld) {


        ielt = msh->GetElementType(iel);
        nDofs = msh->GetElementDofNumber(iel, solType);
        for(int i = 0; i < dim; i++) {
          solD[i].resize(nDofs);
          solDOld[i].resize(nDofs);
          vxHat[i].resize(nDofs);
        }

        for(unsigned inode = 0; inode < nDofs; inode++) {
          unsigned idof = msh->GetSolutionDof(inode, iel, solType);   //local 2 global solution
          unsigned idofX = msh->GetSolutionDof(inode, iel, 2);   //local 2 global solution
          for(int i = 0; i < dim; i++) {
            solDOld[i][inode] = (*mysolution->_SolOld[indexSolD[i]])(idof);
            solD[i][inode] = (*mysolution->_Sol[indexSolD[i]])(idof) - solDOld[i][inode];
            //moving domain
            vxHat[i][inode] = (*msh->_topology->_Sol[i])(idofX) + solDOld[i][inode];
          }
        }

      }
      std::vector <double> xi = particles[iMarker]->GetMarkerLocalCoordinates();

      msh->_finiteElement[ielt][solType]->Jacobian(vxHat, xi, weightHat, phiHat, gradPhiHat);

      std::vector <double> solVpOld(dim);
      particles[iMarker]->GetMarkerVelocity(solVpOld);

      std::vector <double> solApOld(dim);
      particles[iMarker]->GetMarkerAcceleration(solApOld);

      std::vector <double> solDp(dim, 0.);
      //update displacement and acceleration
      for(int i = 0; i < dim; i++) {
        for(unsigned inode = 0; inode < nDofs; inode++) {
          solDp[i] += phiHat[inode] * solD[i][inode];
        }
      }

      particles[iMarker]->SetMarkerDisplacement(solDp);
      particles[iMarker]->UpdateParticleCoordinates();

      std::vector <double> solAp(dim);
      std::vector <double> solVp(dim);
      for(unsigned i = 0; i < dim; i++) {
        solAp[i] = 1. / (beta * dt * dt) * solDp[i] - 1. / (beta * dt) * solVpOld[i] - (1. - 2.* beta) / (2. * beta) * solApOld[i];
        solVp[i] = solVpOld[i] + dt * ((1. - Gamma) * solApOld[i] + Gamma * solAp[i]);
      }

      particles[iMarker]->SetMarkerVelocity(solVp);
      particles[iMarker]->SetMarkerAcceleration(solAp);

      //   update the deformation gradient
      for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
          gradSolDHat[i][j] = 0.;
          for(unsigned inode = 0; inode < nDofs; inode++) {
            gradSolDHat[i][j] +=  gradPhiHat[inode * dim + j] * solD[i][inode];
          }
        }
      }
      std::vector < std::vector < double > > FpOld;
      FpOld = particles[iMarker]->GetDeformationGradient(); //extraction of the deformation gradient
      double FpNew[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
      std::vector < std::vector < double > > Fp(dim);
      for(unsigned i = 0; i < dim; i++) {
        for(unsigned j = 0; j < dim; j++) {
          FpNew[i][j] +=  gradSolDHat[i][j];
        }
      }
      for(unsigned i = 0; i < dim; i++) {
        Fp[i].resize(dim);
        for(unsigned j = 0; j < dim; j++) {
          Fp[i][j] = 0.;
          for(unsigned k = 0; k < dim; k++) {
            Fp[i][j] += FpNew[i][k] * FpOld[k][j];
          }
        }
      }
      particles[iMarker]->SetDeformationGradient(Fp);
      ielOld = iel;
    }
    else {
      break;
    }
  }
  //END loop on bulk particles

  //BEGIN loop on interface particles
  particles = lineI.GetParticles();
  markerOffset = lineI.GetMarkerOffset();
  ielOld = UINT_MAX;
  for(unsigned iMarker = markerOffset[iproc]; iMarker < markerOffset[iproc + 1]; iMarker++) {
    unsigned iel = particles[iMarker]->GetMarkerElement();
    if(iel != UINT_MAX) {
      short unsigned ielt;
      unsigned nDofs;
      //update element related quantities only if we are in a different element
      if(iel != ielOld) {
        ielt = msh->GetElementType(iel);
        nDofs = msh->GetElementDofNumber(iel, solType);
        for(int i = 0; i < dim; i++) {
          solD[i].resize(nDofs);
          solDOld[i].resize(nDofs);
          vxHat[i].resize(nDofs);
        }

        for(unsigned inode = 0; inode < nDofs; inode++) {
          unsigned idof = msh->GetSolutionDof(inode, iel, solType);   //local 2 global solution
          unsigned idofX = msh->GetSolutionDof(inode, iel, 2);   //local 2 global solution
          for(int i = 0; i < dim; i++) {
            solDOld[i][inode] = (*mysolution->_SolOld[indexSolD[i]])(idof);
            solD[i][inode] = (*mysolution->_Sol[indexSolD[i]])(idof) - solDOld[i][inode];
            //moving domain
            vxHat[i][inode] = (*msh->_topology->_Sol[i])(idofX) + solDOld[i][inode];
          }
        }

      }
      std::vector <double> xi = particles[iMarker]->GetMarkerLocalCoordinates();

      msh->_finiteElement[ielt][solType]->Jacobian(vxHat, xi, weightHat, phiHat, gradPhiHat);

      std::vector <double> solVpOld(dim);
      particles[iMarker]->GetMarkerVelocity(solVpOld);

      std::vector <double> solApOld(dim);
      particles[iMarker]->GetMarkerAcceleration(solApOld);

      std::vector <double> solDp(dim, 0.);
      //update displacement and acceleration
      for(int i = 0; i < dim; i++) {
        for(unsigned inode = 0; inode < nDofs; inode++) {
          solDp[i] += phiHat[inode] * solD[i][inode];
        }
      }

      particles[iMarker]->SetMarkerDisplacement(solDp);
      particles[iMarker]->UpdateParticleCoordinates();

      std::vector <double> solAp(dim);
      std::vector <double> solVp(dim);
      for(unsigned i = 0; i < dim; i++) {
        solAp[i] = 1. / (beta * dt * dt) * solDp[i] - 1. / (beta * dt) * solVpOld[i] - (1. - 2.* beta) / (2. * beta) * solApOld[i];
        solVp[i] = solVpOld[i] + dt * ((1. - Gamma) * solApOld[i] + Gamma * solAp[i]);
      }

      particles[iMarker]->SetMarkerVelocity(solVp);
      particles[iMarker]->SetMarkerAcceleration(solAp);

      //   update the deformation gradient
      for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
          gradSolDHat[i][j] = 0.;
          for(unsigned inode = 0; inode < nDofs; inode++) {
            gradSolDHat[i][j] +=  gradPhiHat[inode * dim + j] * solD[i][inode];
          }
        }
      }
      std::vector < std::vector < double > > FpOld;
      FpOld = particles[iMarker]->GetDeformationGradient(); //extraction of the deformation gradient
      double FpNew[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
      std::vector < std::vector < double > > Fp(dim);
      for(unsigned i = 0; i < dim; i++) {
        for(unsigned j = 0; j < dim; j++) {
          FpNew[i][j] += gradSolDHat[i][j];
        }
      }
      for(unsigned i = 0; i < dim; i++) {
        Fp[i].resize(dim);
        for(unsigned j = 0; j < dim; j++) {
          Fp[i][j] = 0.;
          for(unsigned k = 0; k < dim; k++) {
            Fp[i][j] += FpNew[i][k] * FpOld[k][j];
          }
        }
      }
      particles[iMarker]->SetDeformationGradient(Fp);
      ielOld = iel;
    }
    else {
      break;
    }
  }
  //END loop on interface particle

  ProjectGridVelocity(*mlSol);

  //BEGIN loop on elements to update grid velocity and acceleration
  for(unsigned idof = msh->_dofOffset[solType][iproc]; idof < msh->_dofOffset[solType][iproc + 1]; idof++) {
    for(int i = 0; i < dim; i++) {
      mysolution->_Sol[indexSolD[i]]->set(idof, 0.);
    }
  }

  for(int i = 0; i < dim; i++) {
    mysolution->_Sol[indexSolD[i]]->close();
  }
  //END loop on elements to update grid velocity and acceleration

  lineI.UpdateLineMPM();
  bulk.UpdateLineMPM();

  bool updateMat = false;
  lineI.GetParticlesToGridMaterial(updateMat);
  bulk.GetParticlesToGridMaterial(updateMat);


  BuildFlag(*mlSol);
// GetParticleWeights(*mlSol, &bulk);

}

unsigned getNumberOfLayers(const double & a, const double & fac, const bool inverse = true) {

  double fac1  = (inverse) ? fac : 1. / fac;
  double da = 1. / fac1;
  double b =  da;
  unsigned n = 1;

  while(b < a) {
    da /= fac1;
    b += da;
    n++;
    if(n >= 100) {
      std::cout << "Error: number of layer is unbounded, try with a smaller factor\n";
      abort();
    }
  }
  return n;
}

void ProjectGridVelocity(MultiLevelSolution &mlSol) {

  const unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;
  Mesh* msh = mlSol._mlMesh->GetLevel(level);
  Solution* sol  = mlSol.GetSolutionLevel(level);
  const unsigned dim = msh->GetDimension();

  unsigned iproc  = msh->processor_id();
  unsigned nprocs = msh->n_processors();

  vector< vector< double > > solV(dim);      // local solution (velocity)

  vector< bool > nodeFlag;     // local solution (velocity)
  vector< unsigned > idof;     // local solution (velocity)

  vector < double > phi;
  vector <vector < double> > vx(dim);   //vx is coordX in assembly of ex30
  vector <vector < double> > xp;


  //variable-name handling
  const char varname[10][5] = {"DX", "DY", "DZ", "VX", "VY", "VZ"};


  vector <unsigned> indexSolD(dim);
  vector <unsigned> indexSolV(dim);

  for(unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolD[ivar] = mlSol.GetIndex(&varname[ivar][0]);
    indexSolV[ivar] = mlSol.GetIndex(&varname[ivar + 3][0]);
  }
//   unsigned indexNodeFlag =  mlSol.GetIndex("NodeFlag");
  unsigned indexNodeFlag =  mlSol.GetIndex("nflag");

  unsigned solType = mlSol.GetSolutionType(&varname[0][0]);

  sol->_Sol[indexNodeFlag]->zero();

  for(unsigned k = 0; k < dim; k++) {
    (*sol->_SolOld[indexSolV[k]]) = (*sol->_Sol[indexSolV[k]]);
    sol->_Sol[indexSolV[k]]->zero();
  }

  unsigned counter = 0;

  std::vector < std::vector < std::vector <double > > > aP(3);

  //BEGIN loop on elements
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielType = msh->GetElementType(iel);

    unsigned nDofs = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs

    for(unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofs);
      vx[k].resize(nDofs);
    }
    nodeFlag.resize(nDofs);
    idof.resize(nDofs);
    xp.resize(nDofs);
    for(unsigned i = 0; i < nDofs; i++) {
      xp[i].resize(dim);
    }

    for(unsigned i = 0; i < nDofs; i++) {
      idof[i] = msh->GetSolutionDof(i, iel, solType);
      unsigned idofX = msh->GetSolutionDof(i, iel, 2);

      for(unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_SolOld[indexSolV[k]])(idof[i]);
        xp[i][k] = (*msh->_topology->_Sol[k])(idofX);     // coordinates of the reference configuration;
        vx[k][i] = xp[i][k] + (*sol->_Sol[indexSolD[k]])(idof[i]);     // coordinates of the deformed configuration
      }
      nodeFlag[i] = ((*sol->_Sol[indexNodeFlag])(idof[i]) > 0.5) ? true : false;
    }

    bool aPIsInitialized = false;

    double r;
    std::vector <double> xc;
    GetConvexHullSphere(vx, xc, r, 0.0001);  // get the ball that circumscribe the element
    double r2 = r * r;

    std::vector < std::vector< double > > xe; // get the box that encloses the element
    GetBoundingBox(vx, xe, 0.0001);

    for(unsigned i = 0; i < nDofs; i++) {  // loop on the nodes of the reference elements now considered as independent points
      if(!nodeFlag[i]) {
        double d2 = 0.;
        for(int k = 0; k < dim; k++) {
          d2 += (xp[i][k] - xc[k]) * (xp[i][k] - xc[k]);
        }
        bool insideHull = true;
        if(d2 > r2) {
          insideHull = false;
        }
        for(unsigned k = 0; k < dim; k++) {
          if(xp[i][k] < xe[k][0] || xp[i][k] > xe[k][1]) {
            insideHull = false;
          }
        }

        if(insideHull) {  //rough test
          if(!aPIsInitialized) {
            aPIsInitialized = true;
            std::vector < std::vector <double> > x1(dim);
            for(unsigned jtype = 0; jtype < solType + 1; jtype++) {
              ProjectNodalToPolynomialCoefficients(aP[jtype], vx, ielType, jtype) ;
            }
          }
          std::vector <double> xi;
          GetClosestPointInReferenceElement(vx, xp[i], ielType, xi);

          bool inverseMapping = GetInverseMapping(solType, ielType, aP, xp[i], xi, 100);
          if(!inverseMapping) {
            std::cout << "InverseMapping failed at " << iel << " " << idof[i] << std::endl;
          }


          bool insideDomain = CheckIfPointIsInsideReferenceDomain(xi, ielType, 1.e-3);  // fine testing

          if(inverseMapping && insideDomain) {
            sol->_Sol[indexNodeFlag]->add(idof[i], 1.);
            msh->_finiteElement[ielType][solType]->GetPhi(phi, xi);
            //std::cout << iel << " " << i << "  ";
            counter++;
            for(unsigned k = 0; k < dim; k++) {
              double solVk = 0.;
              for(unsigned j = 0; j < nDofs; j++)    {
                solVk += phi[j] * solV[k][j];
              }
              sol->_Sol[indexSolV[k]]->add(idof[i], solVk);
            }
          }
        }
      }
    }
  }

  sol->_Sol[indexNodeFlag]->close();

  for(unsigned k = 0; k < dim; k++) {
    sol->_Sol[indexSolV[k]]->close();
  }

  unsigned c0 = 0;
  for(unsigned i = msh->_dofOffset[solType][iproc]; i < msh->_dofOffset[solType][iproc + 1]; i++) {
    unsigned cnt = static_cast < unsigned >(floor((*sol->_Sol[indexNodeFlag])(i) + 0.5));
    if(cnt == 0) {
      c0++;
    }
    else if(cnt > 1) {
      counter -= (cnt - 1);
      for(unsigned k = 0; k < dim; k++) {
        double velk = (*sol->_Sol[indexSolV[k]])(i) / cnt;
        sol->_Sol[indexSolV[k]]->set(i, velk);
      }
    }
  }

  unsigned counterAll;
  MPI_Reduce(&counter, &counterAll, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  std::cout << "COUNTER = " << counterAll << " " << msh->GetTotalNumberOfDofs(solType) << std::endl;

  idof.resize(c0);
  vector< double > xp0(c0 * dim);

  unsigned c1 = 0;
  for(unsigned i = msh->_dofOffset[solType][iproc]; i < msh->_dofOffset[solType][iproc + 1]; i++) {
    if(static_cast < unsigned >(floor((*sol->_Sol[indexNodeFlag])(i) + 0.5)) == 0) {
      idof[c1] = i;
      for(unsigned k = 0; k < dim; k++) {
        xp0[c1 * dim + k] = (*msh->_topology->_Sol[k])(i);
      }
      c1++;
      if(c1 == c0) break;
    }
  }

  vector< double > xp1;
  for(unsigned jproc = 0; jproc < nprocs; jproc++) {
    c1 = c0;
    MPI_Bcast(&c1, 1, MPI_UNSIGNED, jproc, PETSC_COMM_WORLD);
    if(c1) {
      xp1.resize(c1 * dim);
      if(iproc == jproc) {
        for(unsigned i = 0; i < c1; i++) {
          for(unsigned k = 0; k < dim; k++) {
            xp1[i * dim + k] = xp0[i * dim + k];
          }
        }
      }
      MPI_Bcast(&xp1[0], c1 * dim, MPI_DOUBLE, jproc, PETSC_COMM_WORLD);

      for(unsigned i = 0; i < c1; i++) {
        std::vector < double > xp(dim);
        for(unsigned k = 0; k < dim; k++) {
          xp[k] = xp1[i * dim + k];
        }
        Marker p(xp, 1, VOLUME, sol, solType, true, 1.);
        unsigned mproc = p.GetMarkerProc(sol);
        if(iproc == mproc) {
          unsigned jel = p.GetMarkerElement();
          short unsigned jelType = msh->GetElementType(jel);
          unsigned nDofs = msh->GetElementDofNumber(jel, solType);    // number of solution element dofs

          for(unsigned  k = 0; k < dim; k++) {
            solV[k].resize(nDofs);
            vx[k].resize(nDofs);
          }

          for(unsigned j = 0; j < nDofs; j++) {
            unsigned jdof = msh->GetSolutionDof(j, jel, solType);
            unsigned jdofX = msh->GetSolutionDof(j, jel, 2);
            for(unsigned  k = 0; k < dim; k++) {
              solV[k][j] = (*sol->_SolOld[indexSolV[k]])(jdof);     //velocity to be projected
              vx[k][j] = (*msh->_topology->_Sol[k])(jdofX) + (*sol->_Sol[indexSolD[k]])(jdof);         // coordinates of the deformed configuration
            }
          }
          std::vector <double> xi = p.GetMarkerLocalCoordinates();
          msh->_finiteElement[jelType][solType]->GetPhi(phi, xi);

          std::vector < double > vel(dim, 0.);
          for(unsigned k = 0; k < dim; k++) {
            for(unsigned j = 0; j < nDofs; j++)    {
              vel[k] += phi[j] * solV[k][j];
            }
          }
          MPI_Send(&vel[0], dim, MPI_DOUBLE, jproc, 1, PETSC_COMM_WORLD);

        }
        if(iproc == jproc) {
          std::vector < double > vel(dim);
          MPI_Recv(&vel[0], dim, MPI_DOUBLE, mproc, 1, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
          for(unsigned k = 0; k < dim; k++) {
            sol->_Sol[indexSolV[k]]->set(idof[i], vel[k]);
          }
          counter++;
        }
      }
    }
  }

  for(unsigned k = 0; k < dim; k++) {
    sol->_Sol[indexSolV[k]]->close();
  }

  MPI_Reduce(&counter, &counterAll, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  std::cout << "COUNTER = " << counterAll << " " << msh->GetTotalNumberOfDofs(solType) << std::endl;

}




















