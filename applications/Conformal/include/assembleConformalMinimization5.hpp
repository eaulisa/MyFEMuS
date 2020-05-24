
// Building the Conformal Minimization system.
void AssembleConformalMinimization(MultiLevelProblem& ml_prob) {



  // Pointers to the femus objects
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system< NonLinearImplicitSystem> ("conformal");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();
  Mesh *msh = ml_prob._ml_msh->GetLevel(level);
  elem *el = msh->el;
  MultiLevelSolution *mlSol = ml_prob._ml_sol;

//   if(counter > 0) {
//     UpdateMu(*mlSol);
//   }

  Solution *sol = ml_prob._ml_sol->GetSolutionLevel(level);
  LinearEquationSolver *pdeSys = mlPdeSys->_LinSolver[level];

  // Pointers to global stiffness matrix and residual vector in pdeSys (level).
  SparseMatrix *KK = pdeSys->_KK;
  NumericVector *RES = pdeSys->_RES;

  // Convenience variables to keep track of the dimension.

  const unsigned  DIM = 3;
  const unsigned  dim = 2;


  // Get the process_id (for parallel computation).
  unsigned iproc = msh->processor_id();

  // Setting the reference elements to be equilateral triangles.


  // Get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC).
  unsigned xType = 2;

  std::vector < unsigned > solDxIndex(DIM);
  solDxIndex[0] = mlSol->GetIndex("Dx1");
  solDxIndex[1] = mlSol->GetIndex("Dx2");
  solDxIndex[2] = mlSol->GetIndex("Dx3");
  unsigned solType = mlSol->GetSolutionType(solDxIndex[0]);
  std::vector < unsigned > solDxPdeIndex(DIM);
  solDxPdeIndex[0] = mlPdeSys->GetSolPdeIndex("Dx1");
  solDxPdeIndex[1] = mlPdeSys->GetSolPdeIndex("Dx2");
  solDxPdeIndex[2] = mlPdeSys->GetSolPdeIndex("Dx3");

  std::vector < std::vector < double > > xhat(DIM);
  std::vector < std::vector < double > > solDx(DIM);

  // Get mu info
  std::vector < unsigned > solMuIndex(dim);
  solMuIndex[0] = mlSol->GetIndex("mu1");
  solMuIndex[1] = mlSol->GetIndex("mu2");
  unsigned solTypeMu = mlSol->GetSolutionType(solMuIndex[0]);
  std::vector < std::vector < double > > solMu(dim);

  // Local-to-global pdeSys dofs.
  std::vector < int > SYSDOF;
  // Local residual vectors.
  vector< double > Res;
  // Local Jacobian matrix (ordered by column).
  vector < double > Jac;

  KK->zero();  // Zero all the entries of the Global Matrix
  RES->zero(); // Zero all the entries of the Global Residual

  // ELEMENT LOOP: each process loops only on the elements that it owns.
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    // Numer of solution element dofs.
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nxDofs  = msh->GetElementDofNumber(iel, solType);

    // Resize local arrays.
    for(unsigned K = 0; K < DIM; K++) {
      xhat[K].resize(nxDofs);
      solDx[K].resize(nxDofs);
    }

    // Resize local arrays
    unsigned sizeAll = DIM * nxDofs;

    SYSDOF.resize(sizeAll);
    Res.assign(sizeAll, 0.);
    Jac.assign(sizeAll * sizeAll, 0.);

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nxDofs; i++) {

      // Global-to-local mapping between X solution node and solution dof.
      unsigned iDDof = msh->GetSolutionDof(i, iel, solType);
      unsigned iXDof  = msh->GetSolutionDof(i, iel, xType);

      for(unsigned K = 0; K < DIM; K++) {
//         xhat[K][i] = (*msh->_topology->_Sol[K])(iXDof) + (*sol->_SolOld[solDxIndex[K]])(iDDof);
//         solDx[K][i] = (*sol->_Sol[solDxIndex[K]])(iDDof) - (*sol->_Sol[solDxIndex[K]])(iDDof);

        xhat[K][i] = (*msh->_topology->_Sol[K])(iXDof) + (*sol->_SolOld[solDxIndex[K]])(iDDof);
        solDx[K][i] = (*sol->_Sol[solDxIndex[K]])(iDDof) - (*sol->_Sol[solDxIndex[K]])(iDDof);

        SYSDOF[ K * nxDofs + i] = pdeSys->GetSystemDof(solDxIndex[K], solDxPdeIndex[K], i, iel);
      }
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

      const double *phix;  // local test function

      const double *phix_uv[DIM]; // local test function first order partial derivatives

      double weight; // gauss point weight

      // Get Gauss point weight, test function, and first order derivatives.
      if(ielGeom == HEX) {
        phix = msh->_finiteElement[ielGeom][solType]->GetPhi(ig);
        phix_uv[0] = msh->_finiteElement[ielGeom][solType]->GetDPhiDXi(ig);
        phix_uv[1] = msh->_finiteElement[ielGeom][solType]->GetDPhiDEta(ig);
        phix_uv[2] = msh->_finiteElement[ielGeom][solType]->GetDPhiDZeta(ig);

        std::vector< double > stdVectorPhi;
        std::vector< double > stdVectorPhi_uv;

        msh->_finiteElement[ielGeom][solType]->Jacobian(xhat, ig, weight, stdVectorPhi, stdVectorPhi_uv);
        //weight = msh->_finiteElement[ielGeom][solType]->GetGaussWeight(ig);
      }


      double xhat_uv[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      for(unsigned K = 0; K < DIM; K++) {
        for(int J = 0; J < DIM; J++) {
          for(unsigned i = 0; i < nxDofs; i++) {
            xhat_uv[K][J]    +=  phix_uv[J][i] * xhat[K][i];
          }
        }
      }

      for(unsigned S = 0; S < DIM; S++) { // loop on the shift
        // Compute the shifted metric, metric determinant, and area element.
        double g[dim][dim] = {{0., 0.}, {0., 0.}};
        for(unsigned i = 0; i < dim; i++) { //loop on the derivatives, shift
          unsigned ii = (S + i) % DIM;
          for(unsigned j = 0; j < dim; j++) { //loop on the derivatives, shift
            unsigned jj = (S + j) % DIM;
            for(unsigned K = 0; K < DIM; K++) { // dot product, no shift
              g[i][j] += xhat_uv[K][ii] * xhat_uv[K][jj];
            }
          }
        }
        double detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];

        double Area = weight;// * sqrt(detg);
        double Area2 = weight; // Trick to give equal weight to each element.

        // Compute the metric inverse.
        double gi[dim][dim];
        gi[0][0] =  g[1][1] / detg;
        gi[0][1] = -g[0][1] / detg;
        gi[1][0] = -g[1][0] / detg;
        gi[1][1] =  g[0][0] / detg;

        // Compute components of the unit normal N to the reference surface, shifted
        double normal[3];
        normal[0] = (xhat_uv[(S + 1) % DIM][(S + 0) % DIM] * xhat_uv[(S + 2) % DIM][(S + 1) % DIM] -
                     xhat_uv[(S + 2) % DIM][(S + 0) % DIM] * xhat_uv[(S + 1) % DIM][(S + 1) % DIM]) / sqrt(detg);
        normal[1] = (xhat_uv[(S + 2) % DIM][(S + 0) % DIM] * xhat_uv[(S + 0) % DIM][(S + 1) % DIM] -
                     xhat_uv[(S + 0) % DIM][(S + 0) % DIM] * xhat_uv[(S + 2) % DIM][(S + 1) % DIM]) / sqrt(detg);
        normal[2] = (xhat_uv[(S + 0) % DIM][(S + 0) % DIM] * xhat_uv[(S + 1) % DIM][(S + 1) % DIM] -
                     xhat_uv[(S + 1) % DIM][(S + 0) % DIM] * xhat_uv[(S + 0) % DIM][(S + 1) % DIM]) / sqrt(detg);


        double mu[dim] = {0., 0.};
//         const double *phiMu = msh->_finiteElement[ielGeom][solTypeMu]->GetPhi(ig);  // local test function
//         for(unsigned i = 0; i < nDofsMu; i++) {
//           for(unsigned k = 0; k < 2; k++) {
//             mu[k] += phiMu[i] * solMu[k][i];
//           }
//         }

        //with respect to the shifted variables and normals the quaternions and D remain the same
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
        for(unsigned I = 0; I < DIM; I++) { //loop on the variables, shift
          unsigned II = (S + I) % DIM;
          for(unsigned i = 0; i < nxDofs; i++) { //loop on the dofs, no shift
            unsigned irow = II * nxDofs + i;
            unsigned istart = irow * sizeAll;
            for(unsigned J = 0; J < DIM; J++) { //loop on the variables, shift
              unsigned JJ = (S + J) % DIM;
              for(unsigned j = 0; j < nxDofs; j++) { // loop on the dofs, no shift
                double term = 0.;
                for(unsigned k = 0; k < dim; k++) { //loop on the derivatives, shift
                  unsigned kk = (S + k) % DIM;
                  for(unsigned l = 0; l < dim; l++) { //loop on the derivatives, shift
                    unsigned ll = (S + l) % DIM;
                    term += phix_uv[kk][i] * D[I * 2 + k][J * 2 + l] * phix_uv[ll][j];
                  }
                }
                Jac[istart + JJ * nxDofs + j] += term * Area;
                Res[II * nxDofs + i] -= term * Area * (xhat[JJ][j] + solDx[JJ][j]);
              }
            }
          }
        }
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

  counter++;

  //VecView ( (static_cast<PetscVector*> (RES))->vec(),  PETSC_VIEWER_STDOUT_SELF);
//   MatView((static_cast<PetscMatrix*>(KK))->mat(), PETSC_VIEWER_STDOUT_SELF);
//
//   PetscViewer    viewer;
//   PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL, NULL, 0, 0, 900, 900, &viewer);
//   PetscObjectSetName((PetscObject) viewer, "PWilmore matrix");
//   PetscViewerPushFormat(viewer, PETSC_VIEWER_DRAW_LG);
//   MatView((static_cast<PetscMatrix*>(KK))->mat(), viewer);
//   double a;
//   std::cin >> a;
//
//   std::cout << std::flush;

} // end AssembleO2ConformalMinimization.
