

double surfaceArea0;

// Building the Conformal Minimization system.
void AssembleConformalMinimization(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled

  adept::Stack& s = FemusInit::_adeptStack;

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

  std::vector < unsigned > solDxIndex(DIM);
  solDxIndex[0] = mlSol->GetIndex("Dx1");
  solDxIndex[1] = mlSol->GetIndex("Dx2");
  if(parameter.surface) solDxIndex[2] = mlSol->GetIndex("Dx3");
  unsigned solType = mlSol->GetSolutionType(solDxIndex[0]);
  std::vector < unsigned > solDxPdeIndex(DIM);
  solDxPdeIndex[0] = mlPdeSys->GetSolPdeIndex("Dx1");
  solDxPdeIndex[1] = mlPdeSys->GetSolPdeIndex("Dx2");
  if(parameter.surface) solDxPdeIndex[2] = mlPdeSys->GetSolPdeIndex("Dx3");

  std::vector < std::vector < double > > xhat(DIM);
  std::vector < std::vector < adept::adouble > > x(DIM);
  std::vector < std::vector < adept::adouble > > solDx(DIM);


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

  std::vector < unsigned > solMuPdeIndex(dim);
  solMuPdeIndex[0] = mlPdeSys->GetSolPdeIndex("mu1");
  solMuPdeIndex[1] = mlPdeSys->GetSolPdeIndex("mu2");

  unsigned solTypeMu = mlSol->GetSolutionType(solMuIndex[0]);
  std::vector < std::vector < adept::adouble > > solMu(dim);

//   if(counter > 0) {
//     //LinearImplicitSystem* mlPdeSysMu   = &ml_prob.get_system< LinearImplicitSystem> ("mu");
//     //mlPdeSysMu->MGsolve();
//     //GetFinalMu(*mlSol);
//     //if(counter == 1 ) BuildPMatrix(*mlSol);
//     UpdateMu(*mlSol);
//   }
//   else {
  if(counter == 0) {
    //*(sol->_Sol[solMuIndex[0]]) = 0.1;
    //*(sol->_Sol[solMuIndex[1]]) = 0.2;
    sol->_Sol[solMuIndex[0]]->zero();
    sol->_Sol[solMuIndex[1]]->zero();
  }

  unsigned vAngleIndex = mlSol->GetIndex("vAngle");
  unsigned vAngleType = mlSol->GetSolutionType(vAngleIndex);
  std::vector <double> vAngle;

  // Local-to-global pdeSys dofs.
  std::vector < int > SYSDOF;
  // Local residual vectors.
  vector< adept::adouble > maRes;
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
    unsigned nMuDofs = msh->GetElementDofNumber(iel, solTypeMu);


    // Resize local arrays.
    for(unsigned K = 0; K < DIM; K++) {
      xhat[K].resize(nxDofs);
      x[K].resize(nxDofs);
      solDx[K].resize(nxDofs);
      solL.resize(nLDofs);
    }

    for(unsigned k = 0; k < dim; k++) {
      solMu[k].resize(nMuDofs);
    }

    // Resize local arrays
    unsigned sizeAll = DIM * nxDofs + dim * nMuDofs + nLDofs + areaConstraint;

    SYSDOF.resize(sizeAll);
    maRes.assign(sizeAll, 0.);
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

    for(unsigned i = 0; i < nMuDofs; i++) {

      unsigned iMuDof  = msh->GetSolutionDof(i, iel, solTypeMu);

      for(unsigned k = 0; k < dim; k++) {
        solMu[k][i] = (*sol->_Sol[solMuIndex[k]])(iMuDof);
        // Global-to-global mapping between NDx solution node and pdeSys dof.
        SYSDOF[ DIM * nxDofs + k * nMuDofs + i] = pdeSys->GetSystemDof(solMuIndex[k], solMuPdeIndex[k], i, iel);
      }
    }

    // Local storage of global mapping and solution.
    if(parameter.constraintIsOn) {
      for(unsigned i = 0; i < nLDofs; i++) {

        // Global-to-local mapping between Lambda solution node and solution dof.
        unsigned iLDof = msh->GetSolutionDof(i, iel, solLType);
        solL[i] = (*sol->_Sol[solLIndex])(iLDof);

        // Global-to-global mapping between Lambda solution node and pdeSys dof.
        SYSDOF[DIM * nxDofs + dim * nMuDofs + i] = pdeSys->GetSystemDof(solLIndex, solLPdeIndex, i, iel);
      }
    }
    if(areaConstraint) {
      SYSDOF[sizeAll - 1u ] = lambda1PdeDof;
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

    s.new_recording();

    for(unsigned K = 0; K < DIM; K++) {
      for(unsigned i = 0; i < nxDofs; i++) {
        x[K][i] = xhat[K][i] + solDx[K][i];
      }
    }



    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {


      double weight; // gauss point weight

      msh->_finiteElement[ielGeom][solType]->Jacobian(cX, ig, weight, phi, dphidu);

      adept::adouble mu[2] = {0., 0.};
      const double *phiMu = msh->_finiteElement[ielGeom][solTypeMu]->GetPhi(ig);  // local test function
      for(unsigned i = 0; i < nMuDofs; i++) {
        for(unsigned k = 0; k < 2; k++) {
          mu[k] += phiMu[i] * solMu[k][i];
        }
      }

      //std::cout<<mu[0] << " " <<mu[1]<<" ";

      // Compute the metric, metric determinant, and area element.
      adept::adouble g[dim][dim];
      g[0][0] = (1. + 2.* mu[0] + mu[0] * mu[0]) + mu[1] * mu[1];
      g[1][0] = 2. * mu[1];
      g[0][1] = 2. * mu[1];
      g[1][1] = (1. - 2. * mu[0] + mu[0] * mu[0]) + mu[1] * mu[1];

      adept::adouble detg = g[0][0] * g[1][1] - g[0][1] * g[1][0]; //(1-|mu|^2)^2
      adept::adouble Area = weight * sqrt(detg); // (1-|mu|^2) dA



      // Compute the metric inverse.
      adept::adouble gi[dim][dim];
      gi[0][0] =  g[1][1] / detg;
      gi[0][1] = -g[0][1] / detg;
      gi[1][0] = -g[1][0] / detg;
      gi[1][1] =  g[0][0] / detg;

      //std::cout<<gi[0][0] <<" " << gi[1][0] <<" "<<gi[0][1] <<" " << gi[1][1];

      std::vector<std::vector <adept::adouble> > gradSolX(DIM, std::vector<adept::adouble>(dim, 0.));
      for(unsigned I = 0; I < DIM; I++) {
        for(unsigned i = 0; i < nxDofs; i++) {
          for(unsigned j = 0; j < dim; j++) {
            gradSolX[I][j] += dphidu[i * dim + j] * x[I][i];
          }
        }
      }

      for(unsigned I = 0; I < DIM; I++) {
        for(unsigned i = 0; i < nxDofs; i++) {
          adept::adouble term = 0.;
          for(unsigned j = 0; j < dim; j++) {
            for(unsigned k = 0; k < dim; k++) {
              term +=  dphidu[i * dim + j] * gi[j][k] * gradSolX[I][k];
            }
          }
          maRes[I * nxDofs + i] -= term * Area;
        }
      }

      // adept::adouble fz[DIM][dim];
      // adept::adouble fzb[DIM][dim];

      // for(unsigned I = 0; I < DIM; I++) {
      //   fz[I][0] = 0.5 * ((1. - mu[0]) * gradSolDx[I][0] - mu[1]  * gradSolDx[I][1]);
      //   fz[I][1] = 0.5 * (mu[1] * gradSolDx[I][0] - (1 + mu[1]) * gradSolDx[I][1]);
      //   // fzb[I][0] = 0.5 * ((1. - mu[0]) * gradSolDx[I][0] + mu[1] * gradSolDx[I][1]);
      //   // fzb[I][1] = 0.5 * (-mu[1] * gradSolDx[I][0] + (1 + mu[1]) * gradSolDx[I][1]);
      //   fzb[I][0] = 0.5 * ((1. + mu[0]) * gradSolDx[I][0] + mu[1] * gradSolDx[I][1]);
      //   fzb[I][1] = 0.5 * (mu[1] * gradSolDx[I][0] + (1 - mu[1]) * gradSolDx[I][1]);
      // }

      adept::adouble normal[3] = {0., 0., 1.};

      normal[0] = gradSolX[1][0] * gradSolX[2][1] - gradSolX[2][0] * gradSolX[1][1];
      normal[1] = gradSolX[2][0] * gradSolX[0][1] - gradSolX[0][0] * gradSolX[2][1];
      normal[2] = gradSolX[0][0] * gradSolX[1][1] - gradSolX[1][0] * gradSolX[0][1];

      adept::adouble normN = 0.;
      for(unsigned k = 0; k < DIM; k++) {
        normN += normal[k] * normal[k];
      }
      normN = sqrt(normN);
      for(unsigned k = 0; k < DIM; k++) {
        normal[k] /= normN;
      }

      boost::math::quaternion <double> N(0, normal[0].value(), normal[1].value(), normal[2].value());
      boost::math::quaternion <double> MU(mu[0].value(), mu[1].value() * normal[0].value(), mu[1].value() * normal[1].value(), mu[1].value() * normal[2].value());
      boost::math::quaternion <double> DX1(0, gradSolX[0][0].value(), gradSolX[1][0].value(),gradSolX[2][0].value());
      boost::math::quaternion <double> DX2(0, gradSolX[0][1].value(), gradSolX[1][1].value(),gradSolX[2][1].value());

      boost::math::quaternion <double> DXp = 0.5 * (DX1 + N * DX2);
      boost::math::quaternion <double> DXm = 0.5 * (DX1 - N * DX2);
      boost::math::quaternion <double> Q = DXp * conj(DXm) +
                                           norm(MU)*DXm*conj(DXp) -
                                           (MU*norm(DXm)+conj(MU)*norm(DXp));

      adept::adouble dXm[DIM];
      adept::adouble dXp[DIM];

      adept::adouble NCrossX2[DIM];
      NCrossX2[0] = normal[1] * gradSolX[2][1] - normal[2] * gradSolX[1][1];
      NCrossX2[1] = normal[2] * gradSolX[0][1] - normal[0] * gradSolX[2][1];
      NCrossX2[2] = normal[0] * gradSolX[1][1] - normal[1] * gradSolX[0][1];

      for(unsigned I = 0; I < DIM; I++) {
        dXm[I] = 0.5 * (gradSolX[I][0] - NCrossX2[I]);
        dXp[I] = 0.5 * (gradSolX[I][0] + NCrossX2[I]);
      }

      adept::adouble dXpCrossdXm[DIM];
      dXpCrossdXm[0] = dXp[1] * dXm[2] - dXp[2] * dXm[1];
      dXpCrossdXm[1] = dXp[2] * dXm[0] - dXp[0] * dXm[2];
      dXpCrossdXm[2] = dXp[0] * dXm[1] - dXp[1] * dXm[0];

      adept::adouble dXpDotdXm = 0.;
      adept::adouble dXmNorm2 = 0.;
      adept::adouble dXpNorm2 = 0.;
      // adept::adouble dXpCrossdXmNorm2 = 0.
      for(unsigned I = 0; I < DIM; I++) {
        dXpDotdXm += dXp[I] * dXm[I];
        dXmNorm2  += dXm[I] * dXm[I];
        dXpNorm2  += dXp[I] * dXp[I];
        // dXpCrossdXmNorm2  += dXpCrossdXm[I] * dXpCrossdXm[I];
      }

      adept::adouble muNorm2 = mu[0] * mu[0] + mu[1] * mu[1];

      adept::adouble term[dim];    // (dXp * conj(dXm) + |mu|^2 dXm * conj(dXp)) (real part and N part)
      term[0] = (1 + muNorm2) * dXpDotdXm;
      for(unsigned I = 0; I < DIM; I++) {
        term[1] -= (1 - muNorm2) * dXpCrossdXm[I] * normal[I];
      }

      adept::adouble hopfD[dim];  // real part and N part of Hopf differential
      hopfD[0] = term[0] - mu[0] * (dXmNorm2 + dXpNorm2);
      hopfD[1] = term[1] - mu[1] * (dXmNorm2 - dXpNorm2);


//       std::cout << hopfD[0].value() <<" "<< Q.real() << " aaa ";
//      std::cout << hopfD[1].value() <<" "<< Q.R_component_2() * normal[0].value() + Q.R_component_3() * normal[1].value() + Q.R_component_4() * normal[2].value() << " aaa ";

      // adept::adouble fzfzbb[2];

      // // fzfzbb[0] = fz[0][0] * fzb[0][0] + fz[1][0] * fzb[1][0] + fz[2][0] * fzb[2][0] +
      // //             fz[0][1] * fzb[0][1] + fz[1][1] * fzb[1][1] + fz[2][1] * fzb[2][1];
      // // fzfzbb[1] = fz[0][1] * fzb[0][0] + fz[1][1] * fzb[1][0] + fz[2][1] * fzb[2][0] -
      // //             fz[0][0] * fzb[0][1] - fz[1][0] * fzb[1][1] - fz[2][0] * fzb[2][1];

      // fzfzbb[0] = fz[0][0] * fz[0][0] + fz[1][0] * fz[1][0] + fz[2][0] * fz[2][0] -
      //             fz[0][1] * fz[0][1] - fz[1][1] * fz[1][1] - fz[2][1] * fz[2][1];
      // fzfzbb[1] = 2 * (fz[0][0] * fz[0][1] + fz[1][0] * fz[1][1] + fz[2][0] * fz[2][1]);

      for(unsigned I = 0; I < 2; I++) {
        for(unsigned i = 0; i < nMuDofs; i++) {
          // maRes[DIM * nxDofs + I * nMuDofs  + i] -= phiMu[i] * fzfzbb[I] * weight / sqrt(detg);
          // maRes[DIM * nxDofs + I * nMuDofs  + i] -= 0.001 * mu[I] * phiMu[i] * Area;
          maRes[DIM * nxDofs + I * nMuDofs  + i] -= phiMu[i] * hopfD[I] * Area;
        }
      }


      /*
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

              double solDxg[3] = {0., 0., 0.};
              double solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};

              for(unsigned K = 0; K < DIM; K++) {
                for(unsigned i = 0; i < nxDofs; i++) {
                  solDxg[K] += phi[i] * solDx[K][i];
                }
                for(int j = 0; j < dim; j++) {
                  for(unsigned i = 0; i < nxDofs; i++) {
                    solx_uv[K][j]   += phix_uv[j][i] * (xhat[K][i] + 0.5 * solDx[K][i]);
                  }
                }
              }

              double NArea[DIM];
              NArea[0] = (solx_uv[1][0] * solx_uv[2][1] - solx_uv[2][0] * solx_uv[1][1]);
              NArea[1] = (solx_uv[2][0] * solx_uv[0][1] - solx_uv[0][0] * solx_uv[2][1]);
              NArea[2] = (solx_uv[0][0] * solx_uv[1][1] - solx_uv[1][0] * solx_uv[0][1]);

              // Compute new X minus old X dot N, for "reparametrization".
              double varXdotNArea = 0.;
              for(unsigned K = 0; K < DIM; K++) {
                varXdotNArea += solDxg[K] * NArea[K];
              }

              double delNArea[2][3][3] = {{{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}
                , {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}
              };

              delNArea[0][0][0] =  0.;
              delNArea[0][0][1] =  0.5 * solx_uv[2][1];
              delNArea[0][0][2] = -0.5 * solx_uv[1][1];
              delNArea[0][1][0] = -0.5 * solx_uv[2][1];
              delNArea[0][1][1] =  0.;
              delNArea[0][1][2] =  0.5 * solx_uv[0][1];
              delNArea[0][2][0] =  0.5 * solx_uv[1][1];
              delNArea[0][2][1] = -0.5 * solx_uv[0][1];
              delNArea[0][2][2] =  0.;

              delNArea[1][0][0] =  0.;
              delNArea[1][0][1] = -0.5 * solx_uv[2][0];
              delNArea[1][0][2] =  0.5 * solx_uv[1][0];
              delNArea[1][1][0] =  0.5 * solx_uv[2][0];
              delNArea[1][1][1] =  0.;
              delNArea[1][1][2] = -0.5 * solx_uv[0][0];
              delNArea[1][2][0] = -0.5 * solx_uv[1][0];
              delNArea[1][2][1] =  0.5 * solx_uv[0][0];
              delNArea[1][2][2] =  0.;

              double deldelNArea[2][3][3][2][3] = {
                {
                  { {{0., 0., 0.}, { 0., 0., 0. }},
                    {{0., 0., 0.}, { 0., 0., 0.5}},
                    {{0., 0., 0.}, { 0., -0.5, 0. }}
                  },
                  { {{0., 0., 0.}, { 0., 0., -0.5 }},
                    {{0., 0., 0.}, { 0., 0., 0.}},
                    {{0., 0., 0.}, { 0.5, 0., 0. }}
                  },
                  { {{0., 0., 0.}, { 0., 0.5, 0. }},
                    {{0., 0., 0.}, { -0.5, 0., 0.}},
                    {{0., 0., 0.}, { 0., 0., 0. }}
                  }
                },

                { { {{ 0., 0., 0.}, {0., 0., 0.}},
                    {{ 0., 0., -0.5}, {0., 0., 0.}},
                    {{ 0., 0.5, 0.}, {0., 0., 0.}}
                  },
                  { {{ 0., 0., 0.5}, {0., 0., 0.}},
                    {{ 0., 0., 0.}, {0., 0., 0.}},
                    {{ -0.5, 0., 0.}, {0., 0., 0.}}
                  },
                  { {{ 0., -0.5, 0.}, {0., 0., 0.}},
                    {{ 0.5, 0., 0.}, {0., 0., 0.}},
                    {{ 0., 0., 0.}, {0., 0., 0.}}
                  }
                }
              };

              // Residual term 1 and Linear Jacobian
              for(unsigned i = 0; i < nLDofs; i++) {
                unsigned irow = DIM * nxDofs + dim * nMuDofs + i;
                Res[irow] -= phiL[i] * varXdotNArea * Area2;
                unsigned istart = irow * sizeAll;
                for(unsigned J = 0; J < DIM; J++) {
                  for(unsigned j = 0; j < nxDofs; j++) {
                    double term = 0.;
                    for(unsigned K = 0; K < DIM; K++) {
                      for(unsigned a = 0; a < dim; a++) {
                        term += solDxg[K] * delNArea[a][K][J] * phix_uv[a][j];
                      }
                    }
                    Jac[istart + J * nxDofs + j] += phiL[i] * (phi[j] * NArea[J] + term) * Area2;

                    unsigned jrow = J * nxDofs + j;
                    unsigned jstart = jrow * sizeAll;
                    Jac[jstart + DIM * nxDofs + dim * nMuDofs + i] += phiL[i] * term * Area2;

                  }
                }
              }

              // Residual term 2 and Linear Jacobian.
              for(unsigned I = 0; I < DIM; I++) {
                for(unsigned i = 0; i < nxDofs; i++) {
                  unsigned irow = I * nxDofs + i;
                  Res[irow] -= solLg * phi[i] * NArea[I] * Area2;
                  unsigned istart = irow * sizeAll;
                  for(unsigned j = 0; j < nLDofs; j++) {
                    Jac[istart + DIM * nxDofs + dim * nMuDofs + j] += phiL[j] * phi[i] * NArea[I] * Area2;
                  }

      //             for(unsigned J = 0; J < DIM; J++) {
      //               for(unsigned j = 0; j < nxDofs; j++) {
      //                 double term = 0;
      //                 for(unsigned a = 0; a < dim; a++) {
      //                   term += delNArea[a][I][J] * phix_uv[a][j];
      //                 }
      //                 //Jac[istart + J * nxDofs + j] += solLg * phi[i] * term * Area2;
      //
      //                 unsigned jrow = J * nxDofs + j;
      //                 unsigned jstart = jrow * sizeAll;
      //                 //Jac[jstart + I * nxDofs + i] += solLg * phi[i] * term * Area2;
      //
      //               }
      //             }
                }
              }

              // Residual term 3 and Linear Jacobian.
              for(unsigned I = 0; I < DIM; I++) {
                for(unsigned i = 0; i < nxDofs; i++) {
                  unsigned irow = I * nxDofs + i;
                  double term = 0;
                  for(unsigned a = 0; a < dim; a++) {
                    for(unsigned K = 0; K < DIM; K++) {
                      term += solDxg[K] * delNArea[a][K][I] * phix_uv[a][i];
                    }
                  }
                  Res[irow] -= solLg * term * Area2;

                  // unsigned istart = irow * sizeAll;
                  // for(unsigned J = 0; J < DIM; J++) {
                  //   for(unsigned j = 0; j < nxDofs; j++) {
                  //     double term = 0;
                  //     for(unsigned a = 0; a < dim; a++) {
                  //       for(unsigned K = 0; K < DIM; K++) {
                  //         for(unsigned b = 0; b < dim; b++) {
                  //           term += phix_uv[a][i] * solDxg[K] * deldelNArea[a][K][I][b][J] * phix_uv[b][j];
                  //         }
                  //       }
                  //     }
                  //     Jac[istart + J * nxDofs + j] += solLg * term * Area2;
                  //   }
                  // }
                }
              }
            }


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

            }*/

    } // end GAUSS POINT LOOP

    //------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector
    //copy the value of the adept::adoube aRes in double Res and store


    std::vector < double > res (sizeAll);

    for (int i = 0; i < sizeAll; i++) {
      res[i] = -maRes[i].value();
    }

    RES->add_vector_blocked(res, SYSDOF);

    Jac.resize( sizeAll * sizeAll);

    // define the dependent variables
    s.dependent(&maRes[0], sizeAll);

    // define the independent variables
    for(unsigned I = 0; I < DIM; I++) {
      s.independent(&solDx[I][0], nxDofs);
    }
    for(unsigned I = 0; I < dim; I++) {
      s.independent(&solMu[I][0], nMuDofs);
    }

    // get the jacobian matrix (ordered by column)
    s.jacobian(&Jac[0], true);

    KK->add_matrix_blocked(Jac, SYSDOF, SYSDOF);

    s.clear_independents();
    s.clear_dependents();

    //RES->add_vector_blocked(Res, SYSDOF);
    //KK->add_matrix_blocked(Jac, SYSDOF, SYSDOF);


  } //end ELEMENT LOOP for each process.

  RES->close();
  KK->close();

//  KK->draw();

  if(counter == 0) KK->print_matlab("matrixA", "ascii");
//   if(areaConstraint) {
//     double surfaceAreaAll;
//     MPI_Reduce(&surfaceArea, &surfaceAreaAll, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//     if(counter == 0) surfaceArea0 = surfaceAreaAll;
//     std::cout << "SURFACE = " << surfaceAreaAll << " SURFACE0 = " << surfaceArea0
//               <<  " error = " << (surfaceArea0 - surfaceAreaAll) / surfaceArea0 << std::endl;
//   }


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
