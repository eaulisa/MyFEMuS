

double surfaceArea0;

// Building the Conformal Minimization system.
void AssembleConformalMinimization(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled

  //  Extract pointers to the several objects that we are going to use.
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system< NonLinearImplicitSystem> ("conformal");   // pointer to the linear implicit system named "Poisson"

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  // Pointers to the mesh (level) object and elem object in mesh (level).
  Mesh *msh = ml_prob._ml_msh->GetLevel(level);
  elem *el = msh->el;

  // Pointers to the multilevel solution, solution (level) and equation (level).
  MultiLevelSolution *mlSol = ml_prob._ml_sol;

  Solution *sol = ml_prob._ml_sol->GetSolutionLevel(level);
  LinearEquationSolver *pdeSys = mlPdeSys->_LinSolver[level];

  // Pointers to global stiffness matrix and residual vector in pdeSys (level).
  SparseMatrix *KK = pdeSys->_KK;
  NumericVector *RES = pdeSys->_RES;

  // Convenience variables to keep track of the dimension.
  const unsigned  dim = 2;
  unsigned  DIM = (parameter.surface) ? 3 : 2;

  // Get the process_id (for parallel computation).
  unsigned iproc = msh->processor_id();

  // Setting the reference elements to be equilateral triangles.
  std::vector < std::vector < double > > cX(2);

  std::vector< double > phi;
  std::vector< double > dphidu;
  std::vector <double> phix_uv[dim];


  // Get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC).
  unsigned xType = 2;

  std::vector < unsigned > solDxIndex(dim);
  solDxIndex[0] = mlSol->GetIndex("Dx1");
  solDxIndex[1] = mlSol->GetIndex("Dx2");
  if(parameter.surface) solDxIndex[2] = mlSol->GetIndex("Dx3");
  unsigned solType = mlSol->GetSolutionType(solDxIndex[0]);
  std::vector < unsigned > solDxPdeIndex(DIM);
  solDxPdeIndex[0] = mlPdeSys->GetSolPdeIndex("Dx1");
  solDxPdeIndex[1] = mlPdeSys->GetSolPdeIndex("Dx2");
  if(parameter.surface) solDxPdeIndex[2] = mlPdeSys->GetSolPdeIndex("Dx3");

  std::vector < std::vector < double > > xhat(DIM);
  std::vector < std::vector < double > > solDx(DIM);


  // Get Lambda info
  unsigned solLIndex;
  unsigned solLType;
  unsigned solLPdeIndex;
  std::vector < double > solL;
  if(parameter.constraintIsOn) {
    solLIndex = mlSol->GetIndex("Lambda");
    solLType = mlSol->GetSolutionType(solLIndex);
    solLPdeIndex = mlPdeSys->GetSolPdeIndex("Lambda");
  }

  // Get mu info
  std::vector < unsigned > solMuIndex(dim);
  solMuIndex[0] = mlSol->GetIndex("mu1");
  solMuIndex[1] = mlSol->GetIndex("mu2");
  unsigned solTypeMu = mlSol->GetSolutionType(solMuIndex[0]);
  std::vector < std::vector < double > > solMu(dim);

  if(counter > 0) {
    UpdateMu(*mlSol);
  }
  else {
    sol->_Sol[solMuIndex[0]]->zero();
    sol->_Sol[solMuIndex[1]]->zero();
  }

  unsigned vAngleIndex = mlSol->GetIndex("vAngle");
  unsigned vAngleType = mlSol->GetSolutionType(vAngleIndex);
  std::vector <double> vAngle;

  // Local-to-global pdeSys dofs.
  std::vector < int > SYSDOF;
  // Local residual vectors.
  vector< double > Res;
  // Local Jacobian matrix (ordered by column).
  vector < double > Jac;

  KK->zero();  // Zero all the entries of the Global Matrix
  RES->zero(); // Zero all the entries of the Global Residual

  double surfaceArea = 0.;

  // Setting up solLambda1 (vol) and solLambda2 (area).
  unsigned solLambda1PdeIndex;

  double solLambda1 = 0.;
  unsigned lambda1PdeDof;

  if(areaConstraint) {
    unsigned solLambda1Index;
    solLambda1Index = mlSol->GetIndex("Lambda1");
    solLambda1PdeIndex = mlPdeSys->GetSolPdeIndex("Lambda1");

    if(areaConstraint) {
      double lambda1;
      if(iproc == 0) {
        lambda1 = (*sol->_Sol[solLambda1Index])(0);  // global to local solution
        lambda1PdeDof = pdeSys->GetSystemDof(solLambda1Index, solLambda1PdeIndex, 0, 0);
      }
      MPI_Bcast(&lambda1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&lambda1PdeDof, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
      solLambda1 = lambda1;
    }
    std::vector < double > value(2);
    std::vector < int > row(1);
    std::vector < int > columns(2);
    value[0] = 1;
    value[1] = -1;
    columns[1] = lambda1PdeDof;

    // For equations other than Lagrange multiplier:
    for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
      if(iel > 0) {
        row[0] = pdeSys->GetSystemDof(solLambda1Index, solLambda1PdeIndex, 0, iel);
        columns[0] = row[0];
        KK->add_matrix_blocked(value, row, columns);
      }
    }
  }



  // ELEMENT LOOP: each process loops only on the elements that it owns.
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    // Numer of solution element dofs.
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nxDofs  = msh->GetElementDofNumber(iel, solType);
    unsigned nLDofs = (parameter.constraintIsOn) ? msh->GetElementDofNumber(iel, solLType) : 0;

    // Resize local arrays.
    for(unsigned K = 0; K < DIM; K++) {
      xhat[K].resize(nxDofs);
      solDx[K].resize(nxDofs);
      solL.resize(nLDofs);
    }

    // Resize local arrays
    unsigned sizeAll = DIM * nxDofs + nLDofs + areaConstraint;

    SYSDOF.resize(sizeAll);
    Res.assign(sizeAll, 0.);
    Jac.assign(sizeAll * sizeAll, 0.);

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nxDofs; i++) {

      // Global-to-local mapping between X solution node and solution dof.
      unsigned iDDof = msh->GetSolutionDof(i, iel, solType);
      unsigned iXDof  = msh->GetSolutionDof(i, iel, xType);

      for(unsigned K = 0; K < DIM; K++) {
        xhat[K][i] = (*msh->_topology->_Sol[K])(iXDof) + (*sol->_SolOld[solDxIndex[K]])(iDDof);
        solDx[K][i] = (*sol->_Sol[solDxIndex[K]])(iDDof) - (*sol->_SolOld[solDxIndex[K]])(iDDof);
        // Global-to-global mapping between NDx solution node and pdeSys dof.
        SYSDOF[ K * nxDofs + i] = pdeSys->GetSystemDof(solDxIndex[K], solDxPdeIndex[K], i, iel);
      }
    }

    unsigned nvAngle = msh->GetElementDofNumber(iel, vAngleType);
    vAngle.resize(nvAngle);
    for(unsigned i = 0; i < nvAngle; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, vAngleType);
      vAngle[i] = (*sol->_Sol[vAngleIndex])(idof);
    }

    if(counter == 0) {
      GetConformalCoordinates(msh, conformalType0, iel, solType, vAngle, cX);
    }
    else {
      GetConformalCoordinates(msh, conformalType, iel, solType, vAngle, cX);
    }

    // Local storage of global mapping and solution.
    if(parameter.constraintIsOn) {
      for(unsigned i = 0; i < nLDofs; i++) {

        // Global-to-local mapping between Lambda solution node and solution dof.
        unsigned iLDof = msh->GetSolutionDof(i, iel, solLType);
        solL[i] = (*sol->_Sol[solLIndex])(iLDof);

        // Global-to-global mapping between Lambda solution node and pdeSys dof.
        SYSDOF[DIM * nxDofs + i] = pdeSys->GetSystemDof(solLIndex, solLPdeIndex, i, iel);
      }
    }

    if(areaConstraint) {
      SYSDOF[sizeAll - 1u ] = lambda1PdeDof;
    }


    unsigned nDofsMu  = msh->GetElementDofNumber(iel, solTypeMu);
    for(unsigned k = 0; k < dim; k++) {
      solMu[k].resize(nDofsMu);
    }

    for(unsigned i = 0; i < nDofsMu; i++) {
      unsigned iDof = msh->GetSolutionDof(i, iel, solTypeMu);
      for(unsigned k = 0; k < dim; k++) {
        solMu[k][i] = (*sol->_Sol[solMuIndex[k]])(iDof);
      }
    }

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {


      double weight; // gauss point weight

      msh->_finiteElement[ielGeom][solType]->Jacobian(cX, ig, weight, phi, dphidu);

      phix_uv[0].resize(nxDofs);
      phix_uv[1].resize(nxDofs);
      for(unsigned i = 0; i < nxDofs; i++) {
        phix_uv[0][i] = dphidu[i * dim];
        phix_uv[1][i] = dphidu[i * dim + 1];
      }

      double xhat_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      for(unsigned K = 0; K < DIM; K++) {
        for(int j = 0; j < dim; j++) {
          for(unsigned i = 0; i < nxDofs; i++) {
            xhat_uv[K][j]    += phix_uv[j][i] * xhat[K][i];
          }
        }
      }

      // Compute the metric, metric determinant, and area element.
      double g[dim][dim] = {{0., 0.}, {0., 0.}};
      for(unsigned i = 0; i < dim; i++) {
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned K = 0; K < DIM; K++) {
            g[i][j] += xhat_uv[K][i] * xhat_uv[K][j];
          }
        }
      }
      double detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];
      double Area = weight * sqrt(detg);
      double Area2 = weight; // Trick to give equal weight to each element.

      // Compute the metric inverse.
      double gi[dim][dim];
      gi[0][0] =  g[1][1] / detg;
      gi[0][1] = -g[0][1] / detg;
      gi[1][0] = -g[1][0] / detg;
      gi[1][1] =  g[0][0] / detg;

      // Compute components of the unit normal N to the reference surface
      double normal[3] = {0., 0., 1.};
      if(parameter.surface) {
        normal[0] = (xhat_uv[1][0] * xhat_uv[2][1] - xhat_uv[2][0] * xhat_uv[1][1]) / sqrt(detg);
        normal[1] = (xhat_uv[2][0] * xhat_uv[0][1] - xhat_uv[0][0] * xhat_uv[2][1]) / sqrt(detg);
        normal[2] = (xhat_uv[0][0] * xhat_uv[1][1] - xhat_uv[1][0] * xhat_uv[0][1]) / sqrt(detg);
      }

      double mu[2] = {0., 0.};
      const double *phiMu = msh->_finiteElement[ielGeom][solTypeMu]->GetPhi(ig);  // local test function
      for(unsigned i = 0; i < nDofsMu; i++) {
        for(unsigned k = 0; k < 2; k++) {
          mu[k] += phiMu[i] * solMu[k][i];
        }
      }

      boost::math::quaternion <double> N(0, normal[0], normal[1], normal[2]);
      boost::math::quaternion <double> e1(0., 1., 0., 0.);
      boost::math::quaternion <double> e2(0., 0., 1., 0.);
      boost::math::quaternion <double> e3(0., 0., 0., 1.);
      boost::math::quaternion <double> MU(mu[0], mu[1] * normal[0], mu[1] * normal[1], mu[1] * normal[2]);

      boost::math::quaternion <double> a[6];
      a[0] = e1 - MU * e1;
      a[2] = e2 - MU * e2;
      a[4] = e3 - MU * e3;

      a[1] = N * e1 + MU * N * e1;
      a[3] = N * e2 + MU * N * e2;
      a[5] = N * e3 + MU * N * e3;

      boost::math::quaternion <double> b[6];
      b[0] = -N * e1 + MU * N * e1;
      b[2] = -N * e2 + MU * N * e2;
      b[4] = -N * e3 + MU * N * e3;

      b[1] = e1 + MU * e1;
      b[3] = e2 + MU * e2;
      b[5] = e3 + MU * e3;

      double D[6][6];
      for(unsigned i = 0; i < 6; i++) {
        for(unsigned j = 0; j < 6; j++) {
          D[i][j] = gi[0][0] * (a[i] % a[j]) + gi[0][1] * (a[i] % b[j]) +
                    gi[1][0] * (b[i] % a[j]) + gi[1][1] * (b[i] % b[j]);
        }
      }

      // Quasi-Conformal Minimization Residual and Jacobian.
      for(unsigned I = 0; I < DIM; I++) {
        for(unsigned i = 0; i < nxDofs; i++) {
          unsigned irow = I * nxDofs + i;
          unsigned istart = irow * sizeAll;
          for(unsigned J = 0; J < DIM; J++) {
            for(unsigned j = 0; j < nxDofs; j++) {
              double term = 0.;
              for(unsigned k = 0; k < dim; k++) {
                for(unsigned l = 0; l < dim; l++) {
                  term += phix_uv[k][i] * D[I * dim + k][J * dim + l] * phix_uv[l][j];
                }
              }
              Jac[istart + J * nxDofs + j] += term * Area;
              Res[I * nxDofs + i] -= term * Area * (xhat[J][j] + solDx[J][j]);
            }
          }
        }
      }

      // With normal motion constraint.
      if(parameter.constraintIsOn) {

        const double *phiL = msh->_finiteElement[ielGeom][solLType]->GetPhi(ig);
        double solLg = 0.;
        for(unsigned i = 0; i < nLDofs; i++) {
          solLg += phiL[i] * solL[i];
        }

        // Penalty Residual and Jacobian.
        for(unsigned i = 0; i < nLDofs; i++) {
          unsigned irow = DIM * nxDofs + i;
          Res[irow] -= phiL[i] * (-eps * solLg * Area);
          unsigned istart = irow * sizeAll;
          for(unsigned j = 0; j < nLDofs; j++) {
            Jac[istart + DIM * nxDofs + j] += -eps * phiL[i] * phiL[j] *  Area;
          }
        }

        // Preliminary quantities for constraint equation.
        double solDxg[3] = {0., 0., 0.};
        double solNx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
        for(unsigned K = 0; K < DIM; K++) {
          for(unsigned i = 0; i < nxDofs; i++) {
            solDxg[K] += phi[i] * solDx[K][i];
          }
          for(int j = 0; j < dim; j++) {
            for(unsigned i = 0; i < nxDofs; i++) {
              solNx_uv[K][j]   += phix_uv[j][i] * (xhat[K][i] + solDx[K][i]);
            }
          }
        }

        double gN[dim][dim] = {{0., 0.}, {0., 0.}};
        for(unsigned i = 0; i < dim; i++) {
          for(unsigned j = 0; j < dim; j++) {
            for(unsigned K = 0; K < DIM; K++) {
              gN[i][j] += solNx_uv[K][i] * solNx_uv[K][j];
            }
          }
        }
        double detgN = gN[0][0] * gN[1][1] - gN[0][1] * gN[1][0];
        double AreaN = weight * sqrt(detgN);

        double gNi[dim][dim];
        gNi[0][0] =  gN[1][1] / detgN;
        gNi[0][1] = -gN[0][1] / detgN;
        gNi[1][0] = -gN[1][0] / detgN;
        gNi[1][1] =  gN[0][0] / detgN;

        double normalN[DIM];
        normalN[0] = (solNx_uv[1][0] * solNx_uv[2][1] - solNx_uv[2][0] * solNx_uv[1][1]) / sqrt(detgN);
        normalN[1] = (solNx_uv[2][0] * solNx_uv[0][1] - solNx_uv[0][0] * solNx_uv[2][1]) / sqrt(detgN);
        normalN[2] = (solNx_uv[0][0] * solNx_uv[1][1] - solNx_uv[1][0] * solNx_uv[0][1]) / sqrt(detgN);

        double varXdotN = 0.;
        for(unsigned K = 0; K < DIM; K++) {
          varXdotN += solDxg[K] * (normalN[K] + normal[K]);
        }

        double NvarX[3][3]; //N \otimes DX
        double NdX[2][3][3]; //X_a \otimes N
        double dXdotVarX[2] = {0., 0.}; //c_a
        for(unsigned a = 0; a < 2; a++) {
          for(unsigned I = 0; I < 3; I++) {
            for(unsigned J = 0; J < 3; J++) {
              NvarX[I][J] = normalN[I] * solDxg[J];
              NdX[a][I][J] = normalN[I] * solNx_uv[J][a];
            }
            dXdotVarX[a] += solNx_uv[I][a] * solDxg[I];
          }
        }

        double tensor1[2][2][3][3] = {{
            {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}},
            {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}
          }, {
            {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}},
            {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}
          }
        };

        double tensor2[2][2][3][3] = {{
            {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}},
            {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}
          }, {
            {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}},
            {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}
          }
        };

        for(unsigned a = 0; a < 2; a++) {
          for(unsigned b = 0; b < 2; b++) {
            for(unsigned l = 0; l < 2; l++) {
                for(unsigned k = 0; k < 2; k++) {
                  for(unsigned I = 0; I < 3; I++) {
                    for(unsigned J = 0; J < 3; J++) {
                      tensor1[a][b][I][J] += gNi[a][b] * gNi[k][l] * dXdotVarX[l] * NdX[k][I][J];
                      tensor2[a][b][I][J] += gNi[a][l] * gNi[b][k] * dXdotVarX[l] * NdX[k][J][I];
                    }
                  }
                }
              }
            }
          }

        // First Residual and Linear Jacobian (1')
        for(unsigned i = 0; i < nLDofs; i++) {
          unsigned irow = DIM * nxDofs + i;
          Res[irow] -= phiL[i] * varXdotN * Area;

          unsigned istart = irow * sizeAll;
          for(unsigned K = 0; K < DIM; K++) {
            for(unsigned j = 0; j < nxDofs; j++) {
              Jac[istart + K * nxDofs + j] += phiL[i] * phi[j] * (normalN[K] + normal[K]) * Area;
            }
          }
        }

        // Second Residual and Linear Jacobian (1' transpose)
        for(unsigned K = 0; K < DIM; K++) {
          for(unsigned i = 0; i < nxDofs; i++) {
            unsigned irow = K * nxDofs + i;
            Res[irow] -= solLg * phi[i] * (normalN[K] + normal[K]) * Area;

            unsigned istart = irow * sizeAll;
            for(unsigned j = 0; j < nLDofs; j++) {
              Jac[istart + DIM * nxDofs + j] += phiL[j] * phi[i] * (normalN[K] + normal[K]) * Area;
            }
          }
        }

        // Third Residual and Nonlinear Constraint Jacobian (2'' and transpose)
        for(unsigned I = 0; I < DIM; I++) {
          for(unsigned i = 0; i < nxDofs; i++) {
            unsigned irow = I * nxDofs + i;

            double term1 = 0.;
            for(unsigned a = 0; a < 2; a++) {
              for(unsigned b = 0; b < 2; b++) {
              //term += mGidDXDXN[b][I] * phix_uv[b][i];
                term1 += -gNi[a][b] * dXdotVarX[a] * normalN[I] * phix_uv[b][i];
              }
            }
          //Res[irow] -= solLg * term1 * Area;

            unsigned istart = irow * sizeAll;
            for(unsigned J = 0; J < DIM; J++) {
              for(unsigned j = 0; j < nxDofs; j++) {

                unsigned jrow =  J * nxDofs + j;
                unsigned jstart = jrow * sizeAll;

                double term2 = 0.;
                for(unsigned a = 0; a < 2; a++) {
                  for(unsigned b = 0; b < 2; b++) {
                  //term += mGidDXN[b][I][J] * phix_uv[b][j];
                    term2 += -gNi[a][b] * NdX[a][J][I] * phix_uv[b][j];
                  }
                }
                Jac[istart + J * nxDofs + j] +=  phi[i] * solLg * term2 * Area;
                Jac[jstart + I * nxDofs + i] +=  phi[i] * solLg * term2 * Area;
              }
            }
          }
        }

        // Jacobian (1'' and transpose)
        for(unsigned i = 0; i < nLDofs; i++) {
          unsigned irow = DIM * nxDofs + i;
          unsigned istart = irow * sizeAll;
          for(unsigned J = 0; J < DIM; J++) {
            for(unsigned j = 0; j < nxDofs; j++) {

              unsigned jrow =  J * nxDofs + j;
              unsigned jstart = jrow * sizeAll;

              double term = 0.;
              for(unsigned b = 0; b < 2; b++) {
                for(unsigned a = 0; a < 2; a++) {
                //term += mGidDXDXN[b][J] * phix_uv[b][j];
                    term += -gNi[a][b] * dXdotVarX[a] * normalN[J] * phix_uv[b][j];
                }
              }
              Jac[istart + J * nxDofs + j] += phiL[i] * term * Area; //1''
              Jac[jstart + DIM * nxDofs + i] += phiL[i] * term * Area; //1'' transpose
            }
          }
        }

        // Nonlinear Constraint Jacobian 3'', 3''', 3'
        for(unsigned I = 0; I < DIM; I++) {
          for(unsigned i = 0; i < nxDofs; i++) {
            unsigned irow = I * nxDofs + i;
            unsigned istart = irow * sizeAll;
            for(unsigned J = 0; J < DIM; J++) {
              for(unsigned j = 0; j < nxDofs; j++) {

                double term = 0.;
                for(unsigned a = 0; a < 2; a++) {
                  for(unsigned b = 0; b < 2; b++) {
                     //term += phix_uv[b][i] * (gIgINdDXDXdDX[b][l][I][J] + gIgIdDXdDXDXN[b][l][I][J] - gNi[b][l] * DXN[J][I]) * phix_uv[l][j];
                      //term += phix_uv[a][i] * (2. * tensor1[a][b][I][J] + tensor2[a][b][I][J] - gNi[a][b] * NvarX[I][J]) * phix_uv[b][j];
                  }
                }
                Jac[istart + J * nxDofs + j] +=  solLg * term * Area;
              }
            }
          }
        }

        // double mGidDXDXN[2][3] = {{0., 0., 0.}, {0., 0., 0.}};
        // double mGidDXN[2][3][3] = {{{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}, {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}};
        // for(unsigned b = 0; b < 2; b++) {
        //   for(unsigned a = 0; a < 2; a++) {
        //     for(unsigned I = 0; I < 3; I++) {
        //       for(unsigned J = 0; J < 3; J++) {
        //         mGidDXDXN[b][J] += -gNi[a][b] * solNx_uv[I][a] * DXN[I][J];
        //         mGidDXN[b][I][J] += -gNi[a][b] * solNx_uv[I][a] * normalN[J];
        //       }
        //     }
        //   }
        // }

        // double gIgINdDXDXdDX[2][2][3][3] = {{
        //     {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}},
        //     {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}
        //   }, {
        //     {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}},
        //     {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}
        //   }
        // };
        //
        // double gIgIdDXdDXDXN[2][2][3][3] = {{
        //     {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}},
        //     {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}
        //   }, {
        //     {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}},
        //     {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}
        //   }
        // };
        //
        // for(unsigned b = 0; b < 2; b++) {
        //   for(unsigned l = 0; l < 2; l++) {
        //     for(unsigned a = 0; a < 2; a++) {
        //       for(unsigned k = 0; k < 2; k++) {
        //         for(unsigned I = 0; I < 3; I++) {
        //           for(unsigned J = 0; J < 3; J++) {
        //             gIgINdDXDXdDX[b][l][I][J] += 2. * gNi[b][l] * gNi[a][k] * dDXDX[a] * normalN[I] * solNx_uv[J][k];
        //             gIgIdDXdDXDXN[b][l][I][J] += gNi[a][b] * gNi[k][l] * dDXDX[a] * solNx_uv[I][k] * normalN[J];
        //           }
        //         }
        //       }
        //     }
        //   }
        // }

      }

      // With surface area constraint.
      if(areaConstraint) {

        unsigned irow = sizeAll - 1;
        unsigned istart = irow * sizeAll;
        for(unsigned J = 0; J < DIM; J++) {
          for(unsigned j = 0; j < nxDofs; j++) {
            double term = 0.;
            for(unsigned ii = 0; ii < dim; ii++) {
              for(unsigned jj = 0; jj < dim; jj++) {
                term +=  gi[ii][jj] * xhat_uv[J][ii]  *  phix_uv[jj][j]  ;
              }
            }
            Jac[istart + J * nxDofs + j] += term * Area;
            Res[irow] -= term * solDx[J][j] * Area;
          }
        }


        for(unsigned I = 0; I < DIM; I++) {
          for(unsigned i = 0; i < nxDofs; i++) {
            unsigned irow = I * nxDofs + i;
            unsigned istart = irow * sizeAll;
            double term = 0.;
            for(unsigned ii = 0; ii < dim; ii++) {
              for(unsigned jj = 0; jj < dim; jj++) {
                term +=  gi[ii][jj] * xhat_uv[I][ii]  *  phix_uv[jj][i]  ;
              }
            }
            Jac[istart + DIM * nxDofs] += term * Area;
            Res[irow] -= term * solLambda1 * Area;
          }
        }

        surfaceArea += Area;

      }

    } // end GAUSS POINT LOOP

//------------------------------------------------------------------------
// Add the local Matrix/Vector into the global Matrix/Vector
//copy the value of the adept::adoube aRes in double Res and store

    RES->add_vector_blocked(Res, SYSDOF);
    KK->add_matrix_blocked(Jac, SYSDOF, SYSDOF);


  } //end ELEMENT LOOP for each process.

  RES->close();
  KK->close();

  double surfaceAreaAll;
  MPI_Reduce(&surfaceArea, &surfaceAreaAll, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(counter == 0) surfaceArea0 = surfaceAreaAll;
  std::cout << "SURFACE = " << surfaceAreaAll << " SURFACE0 = " << surfaceArea0
            <<  " error = " << (surfaceArea0 - surfaceAreaAll) / surfaceArea0 << std::endl;



  counter++;

//   VecView ( (static_cast<PetscVector*> (RES))->vec(),  PETSC_VIEWER_STDOUT_SELF);
//     MatView ( (static_cast<PetscMatrix*> (KK))->mat(), PETSC_VIEWER_STDOUT_SELF);
//
//     PetscViewer    viewer;
//     PetscViewerDrawOpen (PETSC_COMM_WORLD, NULL, NULL, 0, 0, 900, 900, &viewer);
//     PetscObjectSetName ( (PetscObject) viewer, "PWilmore matrix");
//     PetscViewerPushFormat (viewer, PETSC_VIEWER_DRAW_LG);
//     MatView ( (static_cast<PetscMatrix*> (KK))->mat(), viewer);
//     double a;
//     std::cin >> a;

} // end AssembleO2ConformalMinimization.
