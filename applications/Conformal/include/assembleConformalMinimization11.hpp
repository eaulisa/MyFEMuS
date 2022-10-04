

double surfaceArea0;



void GetConformalCoordinates(Mesh *msh, const unsigned &conformalType, const unsigned &iel, const unsigned &solType,
                             std::vector<double> &vAngle, std::vector<std::vector<double>> &cX) {
  //this works only for DIM == 2, if DIM == 3 we need to project the mesh coordinates on the tangent plane
  cX.resize(2);
  unsigned nDofs = msh->GetElementDofNumber(iel, solType);
  short unsigned ielGeom = msh->GetElementType(iel);
  if(conformalType == 0) {
  conformal_default:
    if(ielGeom == QUAD) {
      cX[0] = { -1., 1., 1., -1., 0., 1., 0., -1., 0.};
      cX[1] = { -1., -1., 1., 1., -1., 0., 1., 0., 0.};
    }
    else {
      cX[0] = { -0.5, 0.5, 0., 0., 0.25, -0.25, 0. };
      cX[1] = {0., 0., sqrt(3.) / 2., 0., sqrt(3.) / 4., sqrt(3.) / 4., sqrt(3.) / 6.};
    }
    cX[0].resize(nDofs);
    cX[1].resize(nDofs);
  }
  else if(conformalType == 1) {
    cX[0].resize(nDofs);
    cX[1].resize(nDofs);
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned iXDof  = msh->GetSolutionDof(i, iel, 2);
      for(unsigned K = 0; K < 2; K++) {
        cX[K][i] = (*msh->_topology->_Sol[K])(iXDof);
      }
    }
    double scale;
    if(ielGeom == QUAD) {
      scale = sqrt(2. / (sqrt((cX[0][1] * cX[1][2] - cX[0][2] * cX[1][1]) * (cX[0][1] * cX[1][2] - cX[0][2] * cX[1][1]) +
                              (cX[0][2] * cX[1][0] - cX[0][0] * cX[1][2]) * (cX[0][2] * cX[1][0] - cX[0][0] * cX[1][2]) +
                              (cX[0][0] * cX[1][1] - cX[0][1] * cX[1][0]) * (cX[0][0] * cX[1][1] - cX[0][1] * cX[1][0])) +
                         sqrt((cX[0][2] * cX[1][3] - cX[0][3] * cX[1][2]) * (cX[0][2] * cX[1][3] - cX[0][3] * cX[1][2]) +
                              (cX[0][3] * cX[1][1] - cX[0][1] * cX[1][3]) * (cX[0][3] * cX[1][1] - cX[0][1] * cX[1][3]) +
                              (cX[0][1] * cX[1][2] - cX[0][2] * cX[1][1]) * (cX[0][1] * cX[1][2] - cX[0][2] * cX[1][1]))));
    }
    else {
      scale = sqrt((sqrt(3.) / 2.) / sqrt((cX[0][1] * cX[1][2] - cX[0][2] * cX[1][1]) * (cX[0][1] * cX[1][2] - cX[0][2] * cX[1][1]) +
                                          (cX[0][2] * cX[1][0] - cX[0][0] * cX[1][2]) * (cX[0][2] * cX[1][0] - cX[0][0] * cX[1][2]) +
                                          (cX[0][0] * cX[1][1] - cX[0][1] * cX[1][0]) * (cX[0][0] * cX[1][1] - cX[0][1] * cX[1][0])));
    }
    for(unsigned i = 0; i < nDofs; i++) {
      for(unsigned K = 0; K < 2; K++) {
        cX[K][i] *= scale;
      }
    }
  }
  else if(conformalType == 2) {
    cX[0].resize(nDofs);
    cX[1].resize(nDofs);
    GetConformalStructure(vAngle, cX);
  }
  else {
    goto conformal_default;
  }
}

void GetConformalAngles(Mesh *msh, const unsigned &conformalType, const unsigned &iel,
                        std::vector<double> &vAngle, std::vector<double> &eAngle) {

  if(conformalType == 0) {
  conformal_default:
    short unsigned ielGeom = msh->GetElementType(iel);
    if(ielGeom == QUAD) {
      eAngle = {0., 0.5 * M_PI, M_PI, 1.5 * M_PI}; //for square
    }
    else {
      eAngle = {0., 2. / 3. * M_PI, 4. / 3. * M_PI}; // for equilateral triangle
    }
  }
  else if(conformalType == 1 || conformalType == 2) {
    unsigned solType = 0;
    unsigned nDofs = msh->GetElementDofNumber(iel, solType);
    std::vector<std::vector<double>> cX;

    GetConformalCoordinates(msh, conformalType, iel, solType, vAngle,  cX);
    eAngle.resize(nDofs);

    for(unsigned i = 0; i < nDofs; i++) {
      unsigned ip = (i + 1) % nDofs;
      eAngle[i] =  atan2(cX[1][ip] - cX[1][i], cX[0][ip] - cX[0][i]);
    }
  }
  else {
    goto conformal_default;
  }
}



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

  std::vector < unsigned > solDxIndex(DIM);
  solDxIndex[0] = mlSol->GetIndex("Dx1");
  solDxIndex[1] = mlSol->GetIndex("Dx2");
  if(parameter.surface) solDxIndex[2] = mlSol->GetIndex("Dx3");
  unsigned solType = mlSol->GetSolutionType(solDxIndex[0]);
  std::vector < unsigned > solDxPdeIndex(DIM);
  solDxPdeIndex[0] = mlPdeSys->GetSolPdeIndex("Dx1");
  solDxPdeIndex[1] = mlPdeSys->GetSolPdeIndex("Dx2");
  if(parameter.surface) solDxPdeIndex[2] = mlPdeSys->GetSolPdeIndex("Dx3");

  std::vector < unsigned > solMuIndex(dim);
  solMuIndex[0] = mlSol->GetIndex("mu1");
  solMuIndex[1] = mlSol->GetIndex("mu2");
  unsigned solTypeMu = mlSol->GetSolutionType(solMuIndex[0]);
  std::vector < std::vector < double > > solMu(dim);
  std::vector < unsigned > solMuPdeIndex(dim);
  solMuPdeIndex[0] = mlPdeSys->GetSolPdeIndex("mu1");
  solMuPdeIndex[1] = mlPdeSys->GetSolPdeIndex("mu2");

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


  boost::math::quaternion <double> e[3];
  e[0] = boost::math::quaternion <double> (0., 1., 0., 0.);
  e[1] = boost::math::quaternion <double> (0., 0., 1., 0.);
  e[2] = boost::math::quaternion <double> (0., 0., 0., 1.);

  std::vector < boost::math::quaternion <double> > dphicombo[2][3];
  boost::math::quaternion <double> dphiI[2];
  std::vector < boost::math::quaternion <double> > dphim[2][3];
  std::vector < boost::math::quaternion <double> > dphip[2][3];

  // ELEMENT LOOP: each process loops only on the elements that it owns.
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    // Numer of solution element dofs.
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nxDofs  = msh->GetElementDofNumber(iel, solType);
    unsigned nMuDofs  = msh->GetElementDofNumber(iel, solTypeMu);
    unsigned nLDofs = (parameter.constraintIsOn) ? msh->GetElementDofNumber(iel, solLType) : 0;

    // Resize local arrays.
    for(unsigned K = 0; K < DIM; K++) {
      xhat[K].resize(nxDofs);
      solDx[K].resize(nxDofs);
      solL.resize(nLDofs);
    }
    for(unsigned k = 0; k < dim; k++) {
      solMu[k].resize(nMuDofs);
    }

    // Resize local arrays
    unsigned sizeAll = DIM * nxDofs + dim * nMuDofs + nLDofs + areaConstraint;

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

    for(unsigned i = 0; i < nMuDofs; i++) {
      unsigned iDof = msh->GetSolutionDof(i, iel, solTypeMu);
      for(unsigned k = 0; k < dim; k++) {
        solMu[k][i] = (*sol->_Sol[solMuIndex[k]])(iDof);
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










    for(unsigned I = 0; I < DIM; I++) {
      dphim[0][I].resize(nxDofs);
      dphim[1][I].resize(nxDofs);
      dphip[0][I].resize(nxDofs);
      dphip[1][I].resize(nxDofs);

      dphicombo[0][I].resize(nxDofs);
      dphicombo[1][I].resize(nxDofs);
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
      double solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};

      for(unsigned K = 0; K < DIM; K++) {
        for(int j = 0; j < dim; j++) {
          for(unsigned i = 0; i < nxDofs; i++) {
            xhat_uv[K][j]    += phix_uv[j][i] * xhat[K][i];
            solx_uv[K][j]    += phix_uv[j][i] * (xhat[K][i] + solDx[K][i]);
          }
        }
      }

      double mu[2] = {0., 0.};
      const double *phiMu = msh->_finiteElement[ielGeom][solTypeMu]->GetPhi(ig);  // local test function
      for(unsigned i = 0; i < nMuDofs; i++) {
        for(unsigned k = 0; k < 2; k++) {
          mu[k] += phiMu[i] * solMu[k][i];
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

      boost::math::quaternion <double> N(0, normal[0], normal[1], normal[2]);
//       boost::math::quaternion <double> e1(0., 1., 0., 0.);
//       boost::math::quaternion <double> e2(0., 0., 1., 0.);
//       boost::math::quaternion <double> e3(0., 0., 0., 1.);

      boost::math::quaternion <double> MU[2];
      
      MU[0] = boost::math::quaternion <double> (mu[0], mu[1] * normal[0], mu[1] * normal[1], mu[1] * normal[2]);
      MU[1] = -MU[0];

      //std::cout << MUU <<"\n";

      boost::math::quaternion <double> df[2];
      df[0] = boost::math::quaternion <double>(0, solx_uv[0][0], solx_uv[1][0], solx_uv[2][0]);
      df[1] = boost::math::quaternion <double>(0, solx_uv[0][1], solx_uv[1][1], solx_uv[2][1]);

      boost::math::quaternion <double> dfp[2];
      dfp[0] = df[0] - N * df[1];
      dfp[1] = df[1] + N * df[0];

      boost::math::quaternion <double> dfm[2];
      dfm[0] = df[0] + N * df[1];
      dfm[1] = df[1] - N * df[0];

      boost::math::quaternion <double> dfcombo[2];
      dfcombo[0] = dfm[0] - MU[0] * dfp[0];
      dfcombo[1] = dfm[1] - MU[1] * dfp[1];

      
      // Quasi-Conformal Minimization Residual and Jacobian.

      for(unsigned I = 0; I < DIM; I++) {

        for(unsigned i = 0; i < nxDofs; i++) {
          unsigned irow = I * nxDofs + i;
          unsigned istart = irow * sizeAll;

          dphiI[0] = phix_uv[0][i] * e[I];
          dphiI[1] = phix_uv[1][i] * e[I];

          dphim[0][I][i] = dphiI[0] + N * dphiI[1];
          dphim[1][I][i] = dphiI[1] - N * dphiI[0];

          dphip[0][I][i] = dphiI[0] - N * dphiI[1];
          dphip[1][I][i] = dphiI[1] + N * dphiI[0];

          dphicombo[0][I][i] = dphim[0][I][i] - MU[0] * dphip[0][I][i];
          dphicombo[1][I][i] = dphim[1][I][i] - MU[1] * dphip[1][I][i];
        }
      }

      std::vector < boost::math::quaternion <double> > PHIMu[2][2];
      PHIMu[0][0].resize(nMuDofs);
      PHIMu[0][1].resize(nMuDofs);
      PHIMu[1][0].resize(nMuDofs);
      PHIMu[1][1].resize(nMuDofs);
      for(unsigned i = 0; i < nMuDofs; i++) {
        PHIMu[0][0][i] = boost::math::quaternion <double> (phiMu[i], 0, 0, 0);
        PHIMu[0][1][i] = boost::math::quaternion <double> (0, phiMu[i] * normal[0], phiMu[i] * normal[1], phiMu[i] * normal[2]);
        PHIMu[1][0][i] = -PHIMu[0][0][i];
        PHIMu[1][1][i] = -PHIMu[0][1][i];
      }


      // Quasi-Conformal Minimization Residual and Jacobian.
      for(unsigned I = 0; I < DIM; I++) {
        for(unsigned i = 0; i < nxDofs; i++) {
          unsigned irow = I * nxDofs + i;
          unsigned istart = irow * sizeAll;

          for(unsigned J = 0; J < DIM; J++) {
            for(unsigned j = 0; j < nxDofs; j++) {

              double term = gi[0][0] * dphicombo[0][I][i] % dphicombo[0][J][j] +
                            gi[1][0] * dphicombo[1][I][i] % dphicombo[0][J][j] +
                            gi[0][1] * dphicombo[0][I][i] % dphicombo[1][J][j] +
                            gi[1][1] * dphicombo[1][I][i] % dphicombo[1][J][j];

              Jac[istart + J * nxDofs + j] += term * Area;
              Res[I * nxDofs + i] -= (term * (xhat[J][j] + solDx[J][j])) * Area;
            }
          }

          for(unsigned J = 0; J < dim; J++) {
            for(unsigned j = 0; j < nMuDofs; j++) {

              double term = gi[0][0] * dphicombo[0][I][i] % (- PHIMu[0][J][j] * dfp[0]) +
                            gi[1][0] * dphicombo[1][I][i] % (- PHIMu[0][J][j] * dfp[0]) +
                            gi[0][1] * dphicombo[0][I][i] % (- PHIMu[1][J][j] * dfp[1]) +
                            gi[1][1] * dphicombo[1][I][i] % (- PHIMu[1][J][j] * dfp[1]);

              term += gi[0][0] * (- PHIMu[0][J][j] * dphip[0][I][i]) % dfcombo[0] +
                      gi[1][0] * (- PHIMu[1][J][j] * dphip[1][I][i]) % dfcombo[0] +
                      gi[0][1] * (- PHIMu[0][J][j] * dphip[0][I][i]) % dfcombo[1] +
                      gi[1][1] * (- PHIMu[1][J][j] * dphip[1][I][i]) % dfcombo[1];

              Jac[istart + DIM * nxDofs + J * nMuDofs + j] += term * Area;
            }
          }
        }
      }


      for(unsigned I = 0; I < dim; I++) {
        for(unsigned i = 0; i < nMuDofs; i++) {
          double term = gi[0][0] * (-PHIMu[0][I][i] * dfp[0]) % dfcombo[0] +
                        gi[1][0] * (-PHIMu[1][I][i] * dfp[1]) % dfcombo[0] +
                        gi[0][1] * (-PHIMu[0][I][i] * dfp[0]) % dfcombo[1] +
                        gi[1][1] * (-PHIMu[1][I][i] * dfp[1]) % dfcombo[1];



          unsigned irow = DIM * nxDofs + I * nMuDofs + i;
          Res[irow] -= term * Area;
          unsigned istart = irow * sizeAll;
          for(unsigned J = 0; J < dim; J++) {
            for(unsigned j = 0; j < nMuDofs; j++) {

              double term = gi[0][0] * (-PHIMu[0][I][i] * dfp[0]) % (-PHIMu[0][J][j] * dfp[0]) +
                            gi[1][0] * (-PHIMu[1][I][i] * dfp[1]) % (-PHIMu[0][J][j] * dfp[0]) +
                            gi[0][1] * (-PHIMu[0][I][i] * dfp[0]) % (-PHIMu[1][J][j] * dfp[1]) +
                            gi[1][1] * (-PHIMu[1][I][i] * dfp[1]) % (-PHIMu[1][J][j] * dfp[1]);

              Jac[istart + DIM * nxDofs + J * nMuDofs + j] += term *  Area;
            }
          }

          for(unsigned J = 0; J < DIM; J++) {
            for(unsigned j = 0; j < nxDofs; j++) {

              double term = gi[0][0] * (-PHIMu[0][I][i] * dphip[0][J][j]) % dfcombo[0] +
                            gi[1][0] * (-PHIMu[1][I][i] * dphip[1][J][j]) % dfcombo[0] +
                            gi[0][1] * (-PHIMu[0][I][i] * dphip[0][J][j]) % dfcombo[1] +
                            gi[1][1] * (-PHIMu[1][I][i] * dphip[1][J][j]) % dfcombo[1];

              term += gi[0][0] * (-PHIMu[0][I][i] * dfp[0]) % dphicombo[0][J][j] +
                      gi[1][0] * (-PHIMu[1][I][i] * dfp[1]) % dphicombo[0][J][j] +
                      gi[0][1] * (-PHIMu[0][I][i] * dfp[0]) % dphicombo[1][J][j] +
                      gi[1][1] * (-PHIMu[1][I][i] * dfp[1]) % dphicombo[1][J][j];

              Jac[istart + J * nxDofs + j] += term *  Area;
            }
          }



        }
      }


      if(parameter.constraintIsOn) {

        const double *phiL = msh->_finiteElement[ielGeom][solLType]->GetPhi(ig);
        double solLg = 0.;
        for(unsigned i = 0; i < nLDofs; i++) {
          solLg += phiL[i] * solL[i];
        }

        // Penalty Residual and Jacobian.
        for(unsigned i = 0; i < nLDofs; i++) {
          unsigned irow = DIM * nxDofs + dim * nMuDofs + i;
          Res[irow] -= phiL[i] * (-eps * solLg * Area);
          unsigned istart = irow * sizeAll;
          for(unsigned j = 0; j < nLDofs; j++) {
            Jac[istart + DIM * nxDofs + dim * nMuDofs + j] += -eps * phiL[i] * phiL[j] *  Area;
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
              solx_uv[K][j] += phix_uv[j][i] * (xhat[K][i] + 0.5 * solDx[K][i]);
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
            { {{0., 0., 0.}, { 0. , 0. , 0. }},
              {{0., 0., 0.}, { 0. , 0. , 0.5}},
              {{0. , 0., 0.}, { 0. , -0.5, 0. }}
            },
            { {{0., 0., 0.}, { 0. , 0. , -0.5 }},
              {{0., 0., 0.}, { 0. , 0. , 0.}},
              {{0. , 0., 0.}, { 0.5, 0. , 0. }}
            },
            { {{0., 0., 0.}, { 0. , 0.5 , 0. }},
              {{0., 0., 0.}, { -0.5, 0. , 0.}},
              {{0. , 0., 0.}, { 0. , 0. , 0. }}
            }
          },

          { { {{ 0. , 0. , 0.}, {0., 0., 0.}},
              {{ 0. , 0. , -0.5}, {0., 0., 0.}},
              {{ 0. , 0.5, 0.}, {0., 0., 0.}}
            },
            { {{ 0. , 0. , 0.5}, {0., 0., 0.}},
              {{ 0. , 0. , 0.}, {0., 0., 0.}},
              {{ -0.5, 0. , 0.}, {0., 0., 0.}}
            },
            { {{ 0. , -0.5, 0.}, {0., 0., 0.}},
              {{ 0.5, 0. , 0.}, {0., 0., 0.}},
              {{ 0. , 0. , 0.}, {0., 0., 0.}}
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
            Jac[istart + sizeAll - 1] += term * Area;
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

  //KK->draw();

  double surfaceAreaAll;
  MPI_Reduce(&surfaceArea, &surfaceAreaAll, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(counter == 0) surfaceArea0 = surfaceAreaAll;
  std::cout << "SURFACE = " << surfaceAreaAll << " SURFACE0 = " << surfaceArea0
            <<  " error = " << (surfaceArea0 - surfaceAreaAll) / surfaceArea0 << std::endl;

  counter++;

  //VecView((static_cast<PetscVector*>(RES))->vec(),  PETSC_VIEWER_STDOUT_SELF);

} // end AssembleO2ConformalMinimization.


