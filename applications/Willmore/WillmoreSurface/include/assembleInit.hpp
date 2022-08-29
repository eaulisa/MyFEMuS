
//BEGIN Assemble SystemY
void AssembleSystemY(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled

  // Call the adept stack object.
  // Extract pointers to the several objects that we are going to use.
  LinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<LinearImplicitSystem> ("InitY");   // pointer to the linear implicit system named "Poisson"

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  // Point to the mesh and element objects.z
  Mesh *msh = ml_prob._ml_msh->GetLevel(level);
  elem *el = msh->el;

  // Point to mlSol, solution (level), and equation (level) objects.
  MultiLevelSolution *mlSol = ml_prob._ml_sol;
  Solution *sol = ml_prob._ml_sol->GetSolutionLevel(level);
  LinearEquationSolver *pdeSys = mlPdeSys->_LinSolver[level];

  // Point to the global stiffness mtx and residual vectors in pdeSys (level).
  SparseMatrix *KK = pdeSys->_KK;
  NumericVector *RES = pdeSys->_RES;

  // Convenience variables to encode the dimension.
  const unsigned dim = 2;
  const unsigned DIM = 3;

  // Get the process_id (for parallel computation).
  unsigned iproc = msh->processor_id();

  // Extract the solution vector; get solDx positions in the ml_sol object.
  unsigned solDxIndex[DIM];
  solDxIndex[0] = mlSol->GetIndex("Dx1");
  solDxIndex[1] = mlSol->GetIndex("Dx2");
  solDxIndex[2] = mlSol->GetIndex("Dx3");

  // Extract the finite element type for solx.
  unsigned solxType;
  solxType = mlSol->GetSolutionType(solDxIndex[0]);

  // Define solx and solxOld.
  std::vector < double > solx[DIM];

  // Get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC).
  unsigned xType = 2;

  // Get positions of Y in the ml_sol object.
  unsigned solYIndex[DIM];
  solYIndex[0] = mlSol->GetIndex("Y1");
  solYIndex[1] = mlSol->GetIndex("Y2");
  solYIndex[2] = mlSol->GetIndex("Y3");

  // Extract the finite element type for Y.
  unsigned solYType;
  solYType = mlSol->GetSolutionType(solYIndex[0]);

  // Get positions of Y in the pdeSys object.
  unsigned solYPdeIndex[DIM];
  solYPdeIndex[0] = mlPdeSys->GetSolPdeIndex("Y1");
  solYPdeIndex[1] = mlPdeSys->GetSolPdeIndex("Y2");
  solYPdeIndex[2] = mlPdeSys->GetSolPdeIndex("Y3");

  // Define solY and solYOld.
  std::vector < double > solY[DIM];

  // Local-to-global pdeSys dofs.
  std::vector< unsigned > SYSDOF;

  // Define local residual vectors.
  vector < double > Res;
  std::vector< double > aResY[3];

  // Local (column-ordered) Jacobian matrix
  vector < double > Jac;

  KK->zero();  // Set to zero all the entries of the Global Matrix
  RES->zero(); // Set to zero all the entries of the Global Residual

  // Initialize area, volume, P-Willmore energy.
  double surface = 0.;
  double volume = 0.;
  double energy = 0.;

  // ELEMENT LOOP: each process loops only on the elements that it owns.
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    // Number of solution element dofs.
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nxDofs  = msh->GetElementDofNumber(iel, solxType);
    unsigned nYDofs  = msh->GetElementDofNumber(iel, solYType);

    // Resize solution vectors.
    for(unsigned K = 0; K < DIM; K++) {
      solx[K].resize(nxDofs);
      solY[K].resize(nYDofs);
    }

    // Convenience variable for keeping track of problem size.
    unsigned sizeAll = DIM * nYDofs;

    // Resize local arrays.
    SYSDOF.resize(sizeAll);

    Res.assign(sizeAll, 0.);
    Jac.assign(sizeAll * sizeAll, 0.);

    for(unsigned K = 0; K < DIM; K++) {
      aResY[K].assign(nYDofs, 0.);   //resize and set to zero
    }

    // Loop which handles local storage of global mapping and solution X.
    for(unsigned i = 0; i < nxDofs; i++) {

      // Global-to-local mapping between solution node and solution dof.
      unsigned iDDof = msh->GetSolutionDof(i, iel, solxType);
      unsigned iXDof  = msh->GetSolutionDof(i, iel, xType);

      for(unsigned K = 0; K < DIM; K++) {
        solx[K][i] = (*msh->_topology->_Sol[K])(iXDof) + (*sol->_Sol[solDxIndex[K]])(iDDof);
      }
    }

    // Loop which handles local storage of global mapping and solution Y.
    for(unsigned i = 0; i < nYDofs; i++) {

      // Global-to-local mapping between solution node and solution dof.
      unsigned iYDof = msh->GetSolutionDof(i, iel, solYType);
      for(unsigned K = 0; K < DIM; K++) {

        // Global-to-local solutions.
        solY[K][i] = (*sol->_Sol[solYIndex[K]])(iYDof);

        // Global-to-global mapping between solution node and pdeSys dof.
        SYSDOF[K * nYDofs + i] =
          pdeSys->GetSystemDof(solYIndex[K], solYPdeIndex[K], i, iel);
      }
    }


//     for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
//       if(el->GetBoundaryIndex(iel, jface) == 1) {
// 
//           
//         std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAA\n";  
//           
//         const unsigned faceGeom = msh->GetElementFaceType(iel, jface); //edge
//         unsigned nxDofs = msh->GetElementFaceDofNumber(iel, jface, solxType);
// 
//         for(unsigned ig = 0; ig  <  msh->_finiteElement[faceGeom][solxType]->GetGaussPointNumber(); ig++) {
// 
//           double weight = msh->_finiteElement[faceGeom][solxType]->GetGaussWeight(ig);
//           const double *phi = msh->_finiteElement[faceGeom][solxType]->GetPhi(ig);
//           const double *phi_u = msh->_finiteElement[faceGeom][solxType]->GetDPhiDXi(ig);
// 
//           double solxg[3] = {0., 0., 0.};
//           double solxg_u[3] = {0., 0., 0.};
// 
//           for(unsigned I = 0; I < DIM; I++) {
//             for(unsigned i = 0; i < nxDofs; i++) {
//               unsigned inode = msh->GetLocalFaceVertexIndex(iel, jface, i);
//               solxg[I] += phi[i] * solx[I][inode];
//               solxg_u[I] += phi_u[i] * solx[I][inode];
//             }
//           }
// 
//           double length = sqrt(solxg_u[0] * solxg_u[0] + solxg_u[1] * solxg_u[1] + solxg_u[2] * solxg_u[2]) * weight;
// 
// 
//           for(unsigned i = 0; i < nxDofs; i++) {
//             for(unsigned I = 0; I < DIM; I++) {
//               unsigned idof = msh->GetLocalFaceVertexIndex(iel, jface, i);    // face-to-element local node mapping.
//               Res[I * nYDofs + idof] += +solxg[I] * phi[i] * length;
//             }
//           }
//         }
//       }
//     }




    // begin GAUSS POINT LOOP
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solxType]->GetGaussPointNumber(); ig++) {

      const double *phix_uv[dim]; // first order derivatives in (u,v)

      const double *phix;  // local test function
      const double *phiY;  // local test function

      double weight; // gauss point weight

      //Extract Gauss point weight, test functions, and their partial derivatives.
      // "0" is derivative in u, "1" is derivative in v.
      weight = msh->_finiteElement[ielGeom][solxType]->GetGaussWeight(ig);

      phix_uv[0] = msh->_finiteElement[ielGeom][solxType]->GetDPhiDXi(ig);
      phix_uv[1] = msh->_finiteElement[ielGeom][solxType]->GetDPhiDEta(ig);

      phix = msh->_finiteElement[ielGeom][solxType]->GetPhi(ig);
      phiY = msh->_finiteElement[ielGeom][solYType]->GetPhi(ig);

      double solxg[3] = {0., 0., 0.};
      double solYg[3] = {0., 0., 0.};

      // Initialize derivatives of x and W (new, middle, old) at the Gauss points.
      double solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      double solY_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};

      for(unsigned K = 0; K < DIM; K++) {
        for(unsigned i = 0; i < nYDofs; i++) {
          solxg[K] += phix[i] * solx[K][i];
        }
        for(unsigned i = 0; i < nYDofs; i++) {
          solYg[K] += phiY[i] * solY[K][i];
        }
        for(int j = 0; j < dim; j++) {
          for(unsigned i = 0; i < nxDofs; i++) {
            solx_uv[K][j] += phix_uv[j][i] * solx[K][i];
            solY_uv[K][j] += phix_uv[j][i] * solY[K][i];
          }
        }
      }

      // Computing the metric, metric determinant, and area element.
      double g[dim][dim] = {{0., 0.}, {0., 0.}};
      for(unsigned i = 0; i < dim; i++) {
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned K = 0; K < DIM; K++) {
            g[i][j] += solx_uv[K][i] * solx_uv[K][j];
          }
        }
      }
      double detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];
      double Area = weight * sqrt(detg);

      // Computing the unit normal vector N.
      double normal[DIM];
      normal[0] = normalSign * (solx_uv[1][0] * solx_uv[2][1]
                                - solx_uv[2][0] * solx_uv[1][1]) / sqrt(detg);
      normal[1] = normalSign * (solx_uv[2][0] * solx_uv[0][1]
                                - solx_uv[0][0] * solx_uv[2][1]) / sqrt(detg);
      normal[2] = normalSign * (solx_uv[0][0] * solx_uv[1][1]
                                - solx_uv[1][0] * solx_uv[0][1]) / sqrt(detg);


      boost::math::quaternion <double> dfu(0, solx_uv[0][0], solx_uv[1][0], solx_uv[2][0]);
      boost::math::quaternion <double> dfv(0, solx_uv[0][1], solx_uv[1][1], solx_uv[2][1]);
      boost::math::quaternion <double> N(0, normal[0], normal[1], normal[2]);

      boost::math::quaternion <double> dYu(0, solY_uv[0][0], solY_uv[1][0], solY_uv[2][0]);
      boost::math::quaternion <double> dYv(0, solY_uv[0][1], solY_uv[1][1], solY_uv[2][1]);

      // Computing the metric inverse
      double gi[dim][dim];
      gi[0][0] =  g[1][1] / detg;
      gi[0][1] = -g[0][1] / detg;
      gi[1][0] = -g[1][0] / detg;
      gi[1][1] =  g[0][0] / detg;

      // Implement the curvature equation Y = \Delta X .
      for(unsigned K = 0; K < DIM; K++) {
        for(unsigned i = 0; i < nYDofs; i++) {
          unsigned irow = K * nYDofs + i;
          unsigned istart = irow * sizeAll;

          boost::math::quaternion <double> dphiu(0, (K == 0) * phix_uv[0][i], (K == 1) * phix_uv[0][i], (K == 2) * phix_uv[0][i]);
          boost::math::quaternion <double> dphiv(0, (K == 0) * phix_uv[1][i], (K == 1) * phix_uv[1][i], (K == 2) * phix_uv[1][i]);

          double term2 =  (gi[0][0] * (dYu % dphiu) +
                           gi[0][1] * (dYu % dphiv) +
                           gi[1][0] * (dYv % dphiu) +
                           gi[1][1] * (dYv % dphiv));

          double term3 =        (gi[0][0] * (dfu % dphiu) +
                                 gi[0][1] * (dfu % dphiv) +
                                 gi[1][0] * (dfv % dphiu) +
                                 gi[1][1] * (dfv % dphiv));

          // double term4 =        (gi[0][0] * (dfu * conj(dphiu)) +
          //                        gi[0][1] * (dfu * conj(dphiv)) +
          //                        gi[1][0] * (dfv * conj(dphiu)) +
          //                        gi[1][1] * (dfv * conj(dphiv))).R_component_1();

          // double df2 = -0.5 * (N * (dfu * dfv - dfv * dfu)).R_component_1(); //quaternion sqrt(detg)

          // double df3 = -0.5 * (dfu * dfu + dfv * dfv).R_component_1(); //conformal area


          Res[irow] += (solYg[K] * phiY[i] + 0.0005 * term2 + term3) * Area;

          unsigned jstart = istart + K * nYDofs;
          for(unsigned j = 0; j < nYDofs; j++) {

            boost::math::quaternion <double> dphiuJ(0, (K == 0) * phix_uv[0][j], (K == 1) * phix_uv[0][j], (K == 2) * phix_uv[0][j]);
            boost::math::quaternion <double> dphivJ(0, (K == 0) * phix_uv[1][j], (K == 1) * phix_uv[1][j], (K == 2) * phix_uv[1][j]);

            double term2Jac =  (gi[0][0] * (dphiu % dphiuJ) +
                                gi[0][1] * (dphiu % dphivJ) +
                                gi[1][0] * (dphiv % dphiuJ) +
                                gi[1][1] * (dphiv % dphivJ));

            Jac [jstart + j] -= (phiY[i] * phiY[j] + 0.0005 * term2Jac) * Area;
          }
        }
      }

      // Compute new surface area, volume, and P-Willmore energy.
      surface += Area;


      // Computing Y.N and |Y|^2, which are essentially 2H and 4H^2.
      double YdotN = 0.;
      double YdotY = 0.;
      for(unsigned K = 0; K < DIM; K++) {
        YdotN += solYg[K] * normal[K];
        YdotY += solYg[K] * solYg[K];
      }
      double signYdotN = (YdotN >= 0.) ? 1. : -1.;

      // Some necessary quantities when working with polynomials.
      double sumP3 = 0.;
      for(unsigned p = 0; p < 3; p++) {
        double signP = (P[p] % 2u == 0) ? 1. : signYdotN;
        sumP3 += signP * ap[p] * pow(YdotY, P[p] / 2.);
      }

      for(unsigned K = 0; K < DIM; K++) {
        volume += normalSign * (solxg[K] * normal[K]) * Area;
      }
      energy += sumP3 * Area;

    } // end GAUSS POINT LOOP.

    //------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector
    RES->add_vector_blocked(Res, SYSDOF);
    KK->add_matrix_blocked(Jac, SYSDOF, SYSDOF);
  } // End ELEMENT LOOP for each process.

  RES->close();
  KK->close();

  // Get data from each process running in parallel.
  double surfaceAll;
  MPI_Reduce(&surface, &surfaceAll, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(firstTime) surface0 = surfaceAll;
  std::cout << "SURFACE = " << surfaceAll << " SURFACE0 = " << surface0 <<  " error = " << (surface0 - surfaceAll) / surface0 << std::endl;

  double volumeAll;
  MPI_Reduce(&volume, &volumeAll, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(firstTime) volume0 = volumeAll;
  std::cout << "VOLUME = " << volumeAll << " VOLUME0 = " << volume0 <<  " error = " << (volume0 - volumeAll) / volume0 << std::endl;

  double energyAll;
  MPI_Reduce(&energy, &energyAll, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  std::cout << "ENERGY = " << energyAll << std::endl;


  firstTime = false;
}
//END Assemble SystemY


//BEGIN Assemble SystemW
void AssembleSystemW(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled

  // Call the adept stack object.
  // Extract pointers to the several objects that we are going to use.
  LinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<LinearImplicitSystem> ("InitW");   // pointer to the linear implicit system named "Poisson"

  // Define level and time variable.
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  // Point to the mesh and element objects.
  Mesh *msh = ml_prob._ml_msh->GetLevel(level);
  elem *el = msh->el;

  // Point to mlSol, solution (level), and equation (level) objects.
  MultiLevelSolution *mlSol = ml_prob._ml_sol;
  Solution *sol = ml_prob._ml_sol->GetSolutionLevel(level);
  LinearEquationSolver *pdeSys = mlPdeSys->_LinSolver[level];

  // Point to the global stiffness mtx and residual vectors in pdeSys (level).
  SparseMatrix *KK = pdeSys->_KK;
  NumericVector *RES = pdeSys->_RES;

  // Convenience variables to encode the dimension.
  const unsigned dim = 2;
  const unsigned DIM = 3;

  // Get the process_id (for parallel computation).
  unsigned iproc = msh->processor_id();

  // Extract the solution vector; get solDx positions in the ml_sol object.
  unsigned solDxIndex[DIM];
  solDxIndex[0] = mlSol->GetIndex("Dx1");
  solDxIndex[1] = mlSol->GetIndex("Dx2");
  solDxIndex[2] = mlSol->GetIndex("Dx3");

  // Extract the finite element type for solx.
  unsigned solxType;
  solxType = mlSol->GetSolutionType(solDxIndex[0]);

  // Define solx and solxOld.
  std::vector < double > solx[DIM];

  // Get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC).
  unsigned xType = 2;

  // Get positions of Y in the ml_sol object.
  unsigned solYIndex[DIM];
  solYIndex[0] = mlSol->GetIndex("Y1");
  solYIndex[1] = mlSol->GetIndex("Y2");
  solYIndex[2] = mlSol->GetIndex("Y3");

  // Extract the finite element type for Y.
  unsigned solYType;
  solYType = mlSol->GetSolutionType(solYIndex[0]);

  // Define solY and solYOld.
  std::vector < double > solY[DIM];

  // Get positions of W in the ml_sol object.
  unsigned solWIndex[DIM];
  solWIndex[0] = mlSol->GetIndex("W1");
  solWIndex[1] = mlSol->GetIndex("W2");
  solWIndex[2] = mlSol->GetIndex("W3");

  // Extract the finite element type for W.
  unsigned solWType;
  solWType = mlSol->GetSolutionType(solWIndex[0]);

  // Get positions of W in the pdeSys object.
  unsigned solWPdeIndex[DIM];
  solWPdeIndex[0] = mlPdeSys->GetSolPdeIndex("W1");
  solWPdeIndex[1] = mlPdeSys->GetSolPdeIndex("W2");
  solWPdeIndex[2] = mlPdeSys->GetSolPdeIndex("W3");

  // Define local W, WOld solutions.
  std::vector < double > solW[DIM];

  // Local-to-global pdeSys dofs.
  std::vector< unsigned > SYSDOF;

  // Define local residual vectors.
  vector < double > Res;
  std::vector< double > aResW[3];

  // Local (column-ordered) Jacobian matrix
  vector < double > Jac;

  KK->zero();  // Set to zero all the entries of the Global Matrix
  RES->zero(); // Set to zero all the entries of the Global Residual

  // ELEMENT LOOP: each process loops only on the elements that it owns.
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    // Number of solution element dofs.
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nxDofs  = msh->GetElementDofNumber(iel, solxType);
    unsigned nYDofs  = msh->GetElementDofNumber(iel, solYType);
    unsigned nWDofs  = msh->GetElementDofNumber(iel, solWType);

    // Resize solution vectors.
    for(unsigned K = 0; K < DIM; K++) {
      solx[K].resize(nxDofs);
      solY[K].resize(nYDofs);
      solW[K].resize(nWDofs);
    }

    // Convenience variable for keeping track of problem size.
    unsigned sizeAll = DIM * nWDofs;

    // Resize local arrays.
    SYSDOF.resize(sizeAll);

    Res.assign(sizeAll, 0.);
    Jac.assign(sizeAll * sizeAll, 0.);

    for(unsigned K = 0; K < DIM; K++) {
      aResW[K].assign(nWDofs, 0.);   //resize and zet to zero
    }

    // Loop which handles local storage of global mapping and solution X.
    for(unsigned i = 0; i < nxDofs; i++) {
      // Global-to-local mapping between solution node and solution dof.
      unsigned iDDof = msh->GetSolutionDof(i, iel, solxType);
      unsigned iXDof  = msh->GetSolutionDof(i, iel, xType);
      for(unsigned K = 0; K < DIM; K++) {
        solx[K][i] = (*msh->_topology->_Sol[K])(iXDof) + (*sol->_Sol[solDxIndex[K]])(iDDof);
      }
    }

    // Loop which handles local storage of global mapping and solution Y.
    for(unsigned i = 0; i < nYDofs; i++) {
      // Global-to-local mapping between solution node and solution dof.
      unsigned iYDof = msh->GetSolutionDof(i, iel, solYType);
      for(unsigned K = 0; K < DIM; K++) {
        // Global-to-local solutions.
        solY[K][i] = (*sol->_Sol[solYIndex[K]])(iYDof);
      }
    }

    // Loop which handles local storage of global mapping and solution W.
    for(unsigned i = 0; i < nWDofs; i++) {
      // Global-to-local mapping between solution node and solution dof.
      unsigned iWDof = msh->GetSolutionDof(i, iel, solWType);
      for(unsigned K = 0; K < DIM; K++) {
        // Global-to-local solutions.
        solW[K][i] = (*sol->_Sol[solWIndex[K]])(iWDof);
        // Global-to-global mapping between solution node and pdeSys dof.
        SYSDOF[K * nWDofs + i] =
          pdeSys->GetSystemDof(solWIndex[K], solWPdeIndex[K], i, iel);
      }
    }

    // begin GAUSS POINT LOOP
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solxType]->GetGaussPointNumber(); ig++) {

      const double *phix_uv[dim]; // first order derivatives in (u,v)

      const double *phiY;  // local test function
      const double *phiW;  // local test function

      double weight; // gauss point weight

      //Extract Gauss point weight, test functions, and their partial derivatives.
      // "0" is derivative in u, "1" is derivative in v.
      weight = msh->_finiteElement[ielGeom][solxType]->GetGaussWeight(ig);

      phix_uv[0] = msh->_finiteElement[ielGeom][solxType]->GetDPhiDXi(ig);
      phix_uv[1] = msh->_finiteElement[ielGeom][solxType]->GetDPhiDEta(ig);

      phiY = msh->_finiteElement[ielGeom][solYType]->GetPhi(ig);
      phiW = msh->_finiteElement[ielGeom][solWType]->GetPhi(ig);

      // Initialize quantities Y, W at the Gauss points.
      double solYg[3] = {0., 0., 0.};
      double solWg[3] = {0., 0., 0.};

      // Initialize derivatives of x at the Gauss points.
      double solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};

      for(unsigned K = 0; K < DIM; K++) {
        for(unsigned i = 0; i < nYDofs; i++) {
          solYg[K] += phiY[i] * solY[K][i];
        }
        for(unsigned i = 0; i < nWDofs; i++) {
          solWg[K] += phiW[i] * solW[K][i];
        }

        for(int j = 0; j < dim; j++) {
          for(unsigned i = 0; i < nxDofs; i++) {
            solx_uv[K][j] += phix_uv[j][i] * solx[K][i];
          }
        }
      }

      // Computing the metric, metric determinant, and area element.
      double g[dim][dim] = {{0., 0.}, {0., 0.}};
      for(unsigned i = 0; i < dim; i++) {
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned K = 0; K < DIM; K++) {
            g[i][j] += solx_uv[K][i] * solx_uv[K][j];
          }
        }
      }
      double detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];
      double Area = weight * sqrt(detg);

      // Computing the unit normal vector N.
      double normal[DIM];
      normal[0] = normalSign * (solx_uv[1][0] * solx_uv[2][1]
                                - solx_uv[2][0] * solx_uv[1][1]) / sqrt(detg);
      normal[1] = normalSign * (solx_uv[2][0] * solx_uv[0][1]
                                - solx_uv[0][0] * solx_uv[2][1]) / sqrt(detg);
      normal[2] = normalSign * (solx_uv[0][0] * solx_uv[1][1]
                                - solx_uv[1][0] * solx_uv[0][1]) / sqrt(detg);

      // Computing Y.N and |Y|^2, which are essentially 2H and 4H^2.
      double YdotN = 0.;
      double YdotY = 0.;
      for(unsigned K = 0; K < DIM; K++) {
        YdotN += solYg[K] * normal[K];
        YdotY += solYg[K] * solYg[K];
      }
      double signYdotN = (YdotN >= 0.) ? 1. : -1.;

      // Some necessary quantities when working with polynomials.
      double sumP1 = 0.;
      for(unsigned p = 0; p < 3; p++) {
        double signP = (P[p] % 2u == 0) ? 1. : signYdotN;
        sumP1 += signP * ap[p] * P[p] * pow(YdotY, (P[p] - 2.) / 2.);
      }

      // Implement the curvature equation Y = \Delta X .
      for(unsigned K = 0; K < DIM; K++) {

        // Implement the equation relating Y and W.
        for(unsigned i = 0; i < nWDofs; i++) {
          unsigned irow = K * nWDofs + i;
          unsigned istart = irow * sizeAll;

          Res[irow] -= (solWg[K] - sumP1 * solYg[K]) * phiW[i] * Area;

          unsigned jstart = istart + K * nWDofs;
          for(unsigned j = 0; j < nWDofs; j++) {
            Jac[jstart + j] += phiW[i] * phiW[j] * Area;
          }
        }
      }
    } // end GAUSS POINT LOOP.

    //------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector
    RES->add_vector_blocked(Res, SYSDOF);
    KK->add_matrix_blocked(Jac, SYSDOF, SYSDOF);
  } // End ELEMENT LOOP for each process.

  RES->close();
  KK->close();
}
//END Assemble SystemW


