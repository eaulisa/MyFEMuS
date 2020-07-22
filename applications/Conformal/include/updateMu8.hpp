
SparseMatrix* PtP[2][2];

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

double EvaluateMu(MultiLevelSolution & mlSol) {

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1u;

  Solution* sol = mlSol.GetSolutionLevel(level);
  Mesh* msh = mlSol._mlMesh->GetLevel(level);
  elem* el = msh->el;

  //unsigned  dim = msh->GetDimension();
  unsigned dim = 2;
  unsigned DIM = (parameter.surface) ? 3 : 2;

  std::vector < unsigned > indexDx(DIM);
  indexDx[0] = mlSol.GetIndex("Dx1");
  indexDx[1] = mlSol.GetIndex("Dx2");
  if(parameter.surface) indexDx[2] = mlSol.GetIndex("Dx3");
  unsigned solTypeDx = mlSol.GetSolutionType(indexDx[0]);

  std::vector < unsigned > indexMu(dim);
  indexMu[0] = mlSol.GetIndex("mu1");
  indexMu[1] = mlSol.GetIndex("mu2");

  unsigned indexW1 = mlSol.GetIndex("weight1");
  unsigned solType1 = mlSol.GetSolutionType(indexMu[0]);

  std::vector< double > dof1;

  std::vector < std::vector < double > > xhat(DIM);
  std::vector < std::vector < double > > solx(DIM);

  for(unsigned k = 0; k < dim; k++) {
    sol->_Sol[indexMu[k]]->zero();
  }
  sol->_Sol[indexW1]->zero();

  double weight; // gauss point weight

  unsigned iproc = msh->processor_id();
  unsigned nprocs = msh->n_processors();


  std::vector < std::vector < double > > cX(2);

  unsigned vAngleIndex = mlSol.GetIndex("vAngle");
  unsigned vAngleType = mlSol.GetSolutionType(vAngleIndex);
  std::vector <double> vAngle;

  std::vector<double> phi_uv0;
  std::vector<double> phi_uv1;

  std::vector< double > phi;
  std::vector< double > dphidu;

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofs1  = msh->GetElementDofNumber(iel, solType1);
    unsigned nDofsDx  = msh->GetElementDofNumber(iel, solTypeDx);

    dof1.resize(nDofs1);

    for(int K = 0; K < DIM; K++) {
      xhat[K].resize(nDofsDx);
      solx[K].resize(nDofsDx);
    }

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofs1; i++) {
      dof1[i] = msh->GetSolutionDof(i, iel, solType1);
    }
    // local storage of coordinates
    for(unsigned i = 0; i < nDofsDx; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, solTypeDx);
      unsigned xDof  = msh->GetSolutionDof(i, iel, 2);
      for(unsigned K = 0; K < DIM; K++) {
        xhat[K][i] = (*msh->_topology->_Sol[K])(xDof) + (*sol->_SolOld[indexDx[K]])(idof);
        solx[K][i] = (*msh->_topology->_Sol[K])(xDof) + (*sol->_Sol[indexDx[K]])(idof);
      }
    }

    unsigned nvAngle = msh->GetElementDofNumber(iel, vAngleType);
    vAngle.resize(nvAngle);
    for(unsigned i = 0; i < nvAngle; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, vAngleType);
      vAngle[i] = (*sol->_Sol[vAngleIndex])(idof);
    }

    GetConformalCoordinates(msh, conformalType, iel, solTypeDx, vAngle, cX);


// *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solTypeDx]->GetGaussPointNumber(); ig++) {

      double weight; // gauss point weight
      msh->_finiteElement[ielGeom][solTypeDx]->Jacobian(cX, ig, weight, phi, dphidu);
      const double *phi1 = msh->_finiteElement[ielGeom][solType1]->GetPhi(ig);  // local test function

      // Initialize and compute fields at the Gauss points.
      double xhat_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      double solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};

      for(unsigned K = 0; K < DIM; K++) {
        for(int j = 0; j < dim; j++) {
          for(unsigned i = 0; i < nDofsDx; i++) {
            xhat_uv[K][j] += dphidu[i * dim + j] * xhat[K][i];
            solx_uv[K][j] += dphidu[i * dim + j] * solx[K][i];
          }
        }
      }

      // Compute the metric, metric determinant, and area element.
      double g[2][2] = {{0., 0.}, {0., 0.}};
      for(unsigned i = 0; i < dim; i++) {
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned K = 0; K < DIM; K++) {
            g[i][j] += xhat_uv[K][i] * xhat_uv[K][j];
          }
        }
      }
      double detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];

      double normal[3] = {0., 0., 1.};

      if(parameter.surface) {
        normal[0] = (xhat_uv[1][0] * xhat_uv[2][1] - xhat_uv[2][0] * xhat_uv[1][1]) / sqrt(detg);
        normal[1] = (xhat_uv[2][0] * xhat_uv[0][1] - xhat_uv[0][0] * xhat_uv[2][1]) / sqrt(detg);
        normal[2] = (xhat_uv[0][0] * xhat_uv[1][1] - xhat_uv[1][0] * xhat_uv[0][1]) / sqrt(detg);
      }

      boost::math::quaternion <double> N(0, normal[0], normal[1], normal[2]);

      boost::math::quaternion <double> du(0, solx_uv[0][0], solx_uv[1][0], solx_uv[2][0]);
      boost::math::quaternion <double> dv(0, solx_uv[0][1], solx_uv[1][1], solx_uv[2][1]);

      boost::math::quaternion <double> dup = du - N * dv;
      boost::math::quaternion <double> dum = du + N * dv;

      boost::math::quaternion <double> MU = (dum * conj(dup)) / norm(dup);

      double mu[2];

      mu[0] = MU.R_component_1();
      mu[1] = (MU * conj(N)).R_component_1();

      for(unsigned i = 0; i < nDofs1; i++) {
        sol->_Sol[indexW1]->add(dof1[i], phi1[i] * weight);
        for(unsigned k = 0; k < dim; k++) {
          sol->_Sol[indexMu[k]]->add(dof1[i], mu[k] * phi1[i] * weight);
        }
      } // end phi_i loop

    } // end gauss point loop
  } //end element loop for each process*/



  for(unsigned k = 0; k < dim; k++) {
    sol->_Sol[indexMu[k]]->close();
  }
  sol->_Sol[indexW1]->close();

  for(unsigned i = msh->_dofOffset[solType1][iproc]; i < msh->_dofOffset[solType1][iproc + 1]; i++) {

    double weight = (*sol->_Sol[indexW1])(i);

    double mu[2];
    for(unsigned k = 0; k < dim; k++) {
      mu[k] = (*sol->_Sol[indexMu[k]])(i) / weight;
      sol->_Sol[indexMu[k]]->set(i, mu[k]);
    }
  }
  for(unsigned k = 0; k < dim; k++) {
    sol->_Sol[indexMu[k]]->close();
  }

  double MuNormAverage;

  double MuNormLocalSum = 0.;
  double muNormLocalMax = 0.;
  for(unsigned i = msh->_dofOffset[solType1][iproc]; i < msh->_dofOffset[solType1][iproc + 1]; i++) {
    double muNorm = sqrt(pow((*sol->_Sol[indexMu[0]])(i), 2) + pow((*sol->_Sol[indexMu[1]])(i), 2));
    MuNormLocalSum += muNorm;

    muNormLocalMax = (muNorm > muNormLocalMax) ? muNorm : muNormLocalMax;
  }
  double muNormMax;
  MPI_Allreduce(&muNormLocalMax, &muNormMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  std::cout << "\nun-smoothed mu infinity norm = " << muNormMax << std::endl;

  MPI_Allreduce(&MuNormLocalSum, &MuNormAverage, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MuNormAverage /= msh->_dofOffset[solType1][nprocs];

  std::cout << "un-smoothed mu average norm = " << MuNormAverage << std::endl;
  std::cout << "relative difference = " << (muNormMax - MuNormAverage) / MuNormAverage << std::endl;

  return MuNormAverage;
}




void UpdateMu(MultiLevelSolution & mlSol) {

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1u;

  Solution* sol = mlSol.GetSolutionLevel(level);
  Mesh* msh = mlSol._mlMesh->GetLevel(level);

  unsigned dim = 2;

  std::vector < unsigned > indexMu(dim);
  indexMu[0] = mlSol.GetIndex("mu1");
  indexMu[1] = mlSol.GetIndex("mu2");

  *(sol->_SolOld[indexMu[0]]) = *(sol->_Sol[indexMu[0]]);
  *(sol->_SolOld[indexMu[1]]) = *(sol->_Sol[indexMu[1]]);

  double MuNormAverageBefore = EvaluateMu(mlSol);

  unsigned solType1 = mlSol.GetSolutionType(indexMu[0]);

  unsigned iproc = msh->processor_id();
  unsigned nprocs = msh->n_processors();

  NumericVector  *mu1 = NumericVector::build().release();
  mu1->init(*sol->_Sol[indexMu[0]]);
  NumericVector  *mu2 = NumericVector::build().release();
  mu2->init(*sol->_Sol[indexMu[1]]);

  for(unsigned ismooth = 0; ismooth < parameter.numberOfSmoothingSteps; ismooth++) {

    *mu1 = *(sol->_Sol[indexMu[0]]);
    *mu2 = *(sol->_Sol[indexMu[1]]);

    sol->_Sol[indexMu[0]] -> matrix_mult(*mu1, *PtP[0][0]);
    sol->_Sol[indexMu[0]] -> add_vector(*mu2, *PtP[0][1]);

    sol->_Sol[indexMu[1]] -> matrix_mult(*mu1, *PtP[1][0]);
    sol->_Sol[indexMu[1]] -> add_vector(*mu2, *PtP[1][1]);
  }

  delete mu1;
  delete mu2;


  double MuNormAverageAfter;
  {
    double MuNormLocalSum = 0.;
    double muNormLocalMax = 0.;
    for(unsigned i = msh->_dofOffset[solType1][iproc]; i < msh->_dofOffset[solType1][iproc + 1]; i++) {
      double muNorm = sqrt(pow((*sol->_Sol[indexMu[0]])(i), 2) + pow((*sol->_Sol[indexMu[1]])(i), 2));
      MuNormLocalSum += muNorm;

      muNormLocalMax = (muNorm > muNormLocalMax) ? muNorm : muNormLocalMax;
    }
    double muNormMax;
    MPI_Allreduce(&muNormLocalMax, &muNormMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    std::cout << "smoothed mu infinity norm = " << muNormMax << std::endl;

    MPI_Allreduce(&MuNormLocalSum, &MuNormAverageAfter, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MuNormAverageAfter /= msh->_dofOffset[solType1][nprocs];

    std::cout << "smoothed mu average norm = " << MuNormAverageAfter << std::endl;
    std::cout << "relative difference = " << (muNormMax - MuNormAverageAfter) / MuNormAverageAfter << std::endl;

  }


  if(parameter.finalSmoothIsOn) {
    for(unsigned i = msh->_dofOffset[solType1][iproc]; i < msh->_dofOffset[solType1][iproc + 1]; i++) {

      double mu[2];
      for(unsigned k = 0; k < 2; k++) {
        mu[k] = (*sol->_Sol[indexMu[k]])(i);
      }

      double norm = sqrt(mu[0] * mu[0] + mu[1] * mu[1]);
      double cosTheta = mu[0] / norm;
      double sinTheta = mu[1] / norm;

      sol->_Sol[indexMu[0]]->set(i, MuNormAverageBefore * cosTheta);
      sol->_Sol[indexMu[1]]->set(i, MuNormAverageBefore * sinTheta);

    }
    for(unsigned k = 0; k < 2; k++) {
      sol->_Sol[indexMu[k]]->close();
    }
  }


//   //start line search algorithm
//   unsigned DIM = (parameter.surface) ? 3 : 2;
//   std::vector < unsigned > indexDx(DIM);
//   indexDx[0] = mlSol.GetIndex("Dx1");
//   indexDx[1] = mlSol.GetIndex("Dx2");
//   if(parameter.surface) indexDx[2] = mlSol.GetIndex("Dx3");
//   unsigned solTypeDx = mlSol.GetSolutionType(indexDx[0]);
// 
//   std::vector< double > dof1;
// 
//   std::vector < std::vector < double > > xhat(DIM);
//   std::vector < std::vector < double > > solx(DIM);
// 
//   std::vector < std::vector < double > > cX(2);
// 
//   unsigned vAngleIndex = mlSol.GetIndex("vAngle");
//   unsigned vAngleType = mlSol.GetSolutionType(vAngleIndex);
//   std::vector <double> vAngle;
// 
//   std::vector<double> phi_uv0;
//   std::vector<double> phi_uv1;
// 
//   std::vector< double > phi;
//   std::vector< double > dphidu;
// 
// 
//   double num = 0.;
//   double den = 0.;
// 
//   for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
// 
//     short unsigned ielGeom = msh->GetElementType(iel);
//     unsigned nDofs1  = msh->GetElementDofNumber(iel, solType1);
//     unsigned nDofsDx  = msh->GetElementDofNumber(iel, solTypeDx);
// 
//     double mu1i = (*sol->_Sol[indexMu[0]])(iel);
//     double mu2i = (*sol->_Sol[indexMu[1]])(iel);
// 
//     double mu1im1 = (*sol->_SolOld[indexMu[0]])(iel);
//     double mu2im1 = (*sol->_SolOld[indexMu[1]])(iel);
// 
//     for(int K = 0; K < DIM; K++) {
//       xhat[K].resize(nDofsDx);
//       solx[K].resize(nDofsDx);
//     }
// 
//     // local storage of coordinates
//     for(unsigned i = 0; i < nDofsDx; i++) {
//       unsigned idof = msh->GetSolutionDof(i, iel, solTypeDx);
//       unsigned xDof  = msh->GetSolutionDof(i, iel, 2);
//       for(unsigned K = 0; K < DIM; K++) {
//         xhat[K][i] = (*msh->_topology->_Sol[K])(xDof) + (*sol->_SolOld[indexDx[K]])(idof);
//         solx[K][i] = (*msh->_topology->_Sol[K])(xDof) + (*sol->_Sol[indexDx[K]])(idof);
//       }
//     }
// 
//     unsigned nvAngle = msh->GetElementDofNumber(iel, vAngleType);
//     vAngle.resize(nvAngle);
//     for(unsigned i = 0; i < nvAngle; i++) {
//       unsigned idof = msh->GetSolutionDof(i, iel, vAngleType);
//       vAngle[i] = (*sol->_Sol[vAngleIndex])(idof);
//     }
// 
//     GetConformalCoordinates(msh, conformalType, iel, solTypeDx, vAngle, cX);
// 
// 
// // *** Gauss point loop ***
//     for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solTypeDx]->GetGaussPointNumber(); ig++) {
// 
//       double weight; // gauss point weight
//       msh->_finiteElement[ielGeom][solTypeDx]->Jacobian(cX, ig, weight, phi, dphidu);
//       const double *phi1 = msh->_finiteElement[ielGeom][solType1]->GetPhi(ig);  // local test function
// 
//       // Initialize and compute fields at the Gauss points.
//       double xhat_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
//       double solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
// 
//       for(unsigned K = 0; K < DIM; K++) {
//         for(int j = 0; j < dim; j++) {
//           for(unsigned i = 0; i < nDofsDx; i++) {
//             xhat_uv[K][j] += dphidu[i * dim + j] * xhat[K][i];
//             solx_uv[K][j] += dphidu[i * dim + j] * solx[K][i];
//           }
//         }
//       }
// 
//       // Compute the metric, metric determinant, and area element.
//       double g[2][2] = {{0., 0.}, {0., 0.}};
//       for(unsigned i = 0; i < dim; i++) {
//         for(unsigned j = 0; j < dim; j++) {
//           for(unsigned K = 0; K < DIM; K++) {
//             g[i][j] += xhat_uv[K][i] * xhat_uv[K][j];
//           }
//         }
//       }
//       double detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];
// 
//       double normal[3] = {0., 0., 1.};
// 
//       if(parameter.surface) {
//         normal[0] = (xhat_uv[1][0] * xhat_uv[2][1] - xhat_uv[2][0] * xhat_uv[1][1]) / sqrt(detg);
//         normal[1] = (xhat_uv[2][0] * xhat_uv[0][1] - xhat_uv[0][0] * xhat_uv[2][1]) / sqrt(detg);
//         normal[2] = (xhat_uv[0][0] * xhat_uv[1][1] - xhat_uv[1][0] * xhat_uv[0][1]) / sqrt(detg);
//       }
// 
//       boost::math::quaternion <double> N(0, normal[0], normal[1], normal[2]);
// 
//       boost::math::quaternion <double> du(0, solx_uv[0][0], solx_uv[1][0], solx_uv[2][0]);
//       boost::math::quaternion <double> dv(0, solx_uv[0][1], solx_uv[1][1], solx_uv[2][1]);
// 
//       boost::math::quaternion <double> dup = du - N * dv;
//       boost::math::quaternion <double> dum = du + N * dv;
// 
//       boost::math::quaternion <double> MUi(mu1i, mu2i * normal[0], mu2i * normal[1], mu2i * normal[2]);
//       boost::math::quaternion <double> MUim1(mu1im1, mu2im1 * normal[0], mu2im1 * normal[1], mu2im1 * normal[2]);
// 
//       num += (dum - MUi * dup) % ((MUim1 - MUi) * dup);
//       den += ((MUim1 - MUi) * dup) % ((MUim1 - MUi) * dup);
// 
//     }
//   }
// 
//   double t = num / den;
// 
//   sol->_Sol[indexMu[0]]->scale(1. - t);
//   sol->_Sol[indexMu[1]]->scale(1. - t);
// 
//   sol->_Sol[indexMu[0]]->add(t, *(sol->_SolOld[indexMu[0]]));
//   sol->_Sol[indexMu[1]]->add(t, *(sol->_SolOld[indexMu[1]]));




}




void BuildMuMatrix(MultiLevelProblem& ml_prob) {
  LinearImplicitSystem* mlPdeSys   = &ml_prob.get_system< LinearImplicitSystem> ("mu");   // pointer to the linear implicit system named "Poisson"

  const unsigned level = mlPdeSys->GetLevelToAssemble();
  LinearEquationSolver *pdeSys = mlPdeSys->_LinSolver[level];
  SparseMatrix *KK = pdeSys->_KK;
  // Pointers to the mesh (level) object and elem object in mesh (level).
  Mesh *msh = ml_prob._ml_msh->GetLevel(level);
  elem *el = msh->el;

  // Pointers to the multilevel solution, solution (level) and equation (level).
  MultiLevelSolution *mlSol = ml_prob._ml_sol;
  Solution* sol = mlSol->GetSolutionLevel(level);

  unsigned dim = 2;

  unsigned iproc = msh->processor_id();
  unsigned nprocs = msh->n_processors();

  unsigned vAngleIndex = mlSol->GetIndex("vAngle");
  unsigned vAngleType = mlSol->GetSolutionType(vAngleIndex);
  std::vector <double> vAngle;
  std::vector <double> eAngle;

  unsigned indexCntEdge = mlSol->GetIndex("cntEdge");

  unsigned elType = 3;
  unsigned faceType = 2;

  unsigned nel = msh->_dofOffset[elType][nprocs];
  unsigned nel_loc = msh->_dofOffset[elType][iproc + 1] - msh->_dofOffset[elType][iproc];

  unsigned nface = msh->_dofOffset[faceType][nprocs];
  unsigned nface_loc = msh->_dofOffset[faceType][iproc + 1] - msh->_dofOffset[faceType][iproc];


  std::vector < SparseMatrix* > PIJ(4);

  for(unsigned k = 0; k < 4; k++) {
    PIJ[k] = SparseMatrix::build().release();
    PIJ[k]->init(nface, nel, nface_loc, nel_loc, 10, 10);
  }

  std::vector < double > PIJl[4];

  std::vector< unsigned > irow;//loval to global mapping
  std::vector< unsigned > icolumn;//loval to global mapping

  for(int iel = msh->_dofOffset[elType][iproc]; iel < msh->_dofOffset[elType][iproc + 1]; iel++) {

    unsigned nvAngle = msh->GetElementDofNumber(iel, vAngleType);
    vAngle.resize(nvAngle);
    for(unsigned i = 0; i < nvAngle; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, vAngleType);
      vAngle[i] = (*sol->_Sol[vAngleIndex])(idof);
    }

    GetConformalAngles(msh, conformalType, iel, vAngle, eAngle);

    unsigned localDofOffset = msh->GetElementDofNumber(iel, 0);
    unsigned nFaces = msh->GetElementFaceNumber(iel);
    for(unsigned k = 0; k < 4; k++) {
      PIJl[k].resize(nFaces);
    }
    irow.resize(nFaces);
    icolumn.assign(1, iel);
    for(unsigned iface = 0; iface < nFaces; iface++) {

      irow[iface] = msh->GetSolutionDof(localDofOffset + iface, iel, faceType);

      double a = cos(eAngle[iface]);
      double b = sin(eAngle[iface]);

      PIJl[0][iface] = (a * a - b * b);
      PIJl[1][iface] = 2. * a * b;

      PIJl[2][iface] = - 2. * a * b;
      PIJl[3][iface] = (a * a - b * b);

      sol->_Sol[indexCntEdge]->add(irow[iface], 1);
    }
    for(unsigned k = 0; k < 4; k++) {
      PIJ[k]->add_matrix_blocked(PIJl[k], irow, icolumn);
    }
  }
  for(unsigned k = 0; k < 4; k++) {
    PIJ[k]->close();
  }
  sol->_Sol[indexCntEdge]->close();

  for(int iface = msh->_dofOffset[faceType][iproc]; iface < msh->_dofOffset[faceType][iproc + 1]; iface++) {
    double value = (*sol->_Sol[indexCntEdge])(iface);
    if(value > 0.5) {
      sol->_Sol[indexCntEdge]->set(iface, 1. / value);
    }
  }
  sol->_Sol[indexCntEdge]->close();
  for(unsigned k = 0; k < 4; k++) {
    MatDiagonalScale((static_cast<PetscMatrix*>(PIJ[k]))->mat(), (static_cast<PetscVector*>(sol->_Sol[indexCntEdge]))->vec(), NULL);
  }

  SparseMatrix *P = SparseMatrix::build().release();
  P->init(2, 2, PIJ);

  std::vector < std::vector < double > > Rl[2][2];
  std::vector < double > Ql[2][2];
  std::vector < double > Sl[2][2];

  std::vector < SparseMatrix* >  SIJ(4);

  for(unsigned k = 0; k < 4; k++) {
    SIJ[k] = SparseMatrix::build().release();
    SIJ[k]->init(nel, nface, nel_loc, nface_loc, 10, 10);
  }

  for(int iel = msh->_dofOffset[elType][iproc]; iel < msh->_dofOffset[elType][iproc + 1]; iel++) {

    unsigned nvAngle = msh->GetElementDofNumber(iel, vAngleType);
    vAngle.resize(nvAngle);
    for(unsigned i = 0; i < nvAngle; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, vAngleType);
      vAngle[i] = (*sol->_Sol[vAngleIndex])(idof);
    }

    GetConformalAngles(msh, conformalType, iel, vAngle, eAngle);

    unsigned localDofOffset = msh->GetElementDofNumber(iel, 0);
    unsigned nFaces = msh->GetElementFaceNumber(iel);
    for(unsigned I = 0; I < 2; I++) {
      for(unsigned J = 0; J < 2; J++) {
        Rl[I][J].resize(nFaces);
        Ql[I][J].assign(nFaces, 0);
        Sl[I][J].assign(nFaces, 0);
        for(unsigned i = 0; i < nFaces; i++) {
          Rl[I][J][i].resize(nFaces);
        }
      }
    }


    icolumn.resize(nFaces);
    irow.assign(1, iel);

    for(unsigned iface = 0; iface < nFaces; iface++) {

      icolumn[iface] = msh->GetSolutionDof(localDofOffset + iface, iel, faceType);

      double a = cos(eAngle[iface]);
      double b = sin(eAngle[iface]);

      Rl[0][0][iface][iface] = (a * a - b * b);
      Rl[0][1][iface][iface] = - 2. * a * b;
      Rl[1][0][iface][iface] =  2. * a * b;
      Rl[1][1][iface][iface] = (a * a - b * b);

    }

//     Ql[0][0][1] = 1;
//     Ql[0][0][3] = -1;
//     Ql[0][1][2] = -1;
//     Ql[0][1][0] = 1;
//
//     Ql[1][0][2] =  1;
//     Ql[1][0][0] = -1;
//     Ql[1][1][1] =  1;
//     Ql[1][1][3] = -1;



    Ql[0][0][1] = -0.25;
    Ql[0][0][3] = -0.25;
    Ql[0][0][2] = -0.25;
    Ql[0][0][0] = -0.25;

    Ql[1][1][2] = -0.25;
    Ql[1][1][0] = -0.25;
    Ql[1][1][1] = -0.25;
    Ql[1][1][3] = -0.25;


    for(unsigned I = 0; I < 2; I++) {
      for(unsigned J = 0; J < 2; J++) {
        for(unsigned K = 0; K < 2; K++) {
          for(unsigned j = 0; j < nFaces; j++) {
            for(unsigned k = 0; k < nFaces; k++) {
              Sl[I][J][j] +=  Ql[I][K][k]   *  Rl[K][J][k][j];
            }
          }
        }
      }
    }


    for(unsigned k = 0; k < 4; k++) {
      SIJ[k]->add_matrix_blocked(Sl[k / 2][k % 2], irow, icolumn);
    }
  }
  for(unsigned k = 0; k < 4; k++) {
    SIJ[k]->close();
  }
  sol->_Sol[indexCntEdge]->close();

  std::vector < SparseMatrix* >  AIJ(4);

  unsigned muIndex = mlSol->GetIndex("mu1");
  NumericVector  *D = NumericVector::build().release();
  D->init(*sol->_Sol[muIndex]);
  *D = 1.;

  std::vector < SparseMatrix* >  DIJ(4);
  DIJ[0] = SparseMatrix::build().release();
  DIJ[0]->init(nel, nel, nel_loc, nel_loc, 1, 1);
  DIJ[0]->matrix_set_diagonal_values(*D);
  DIJ[1] = NULL;
  DIJ[2] = NULL;
  DIJ[3] = DIJ[0];

  AIJ[0] = SparseMatrix::build().release();
  AIJ[0]->init(2, 2, DIJ);

  AIJ[1] = NULL;
  AIJ[2] = NULL;

  AIJ[2] = SparseMatrix::build().release();
  AIJ[2]->init(2, 2, SIJ);
  AIJ[2]->matrix_RightMatMult(*P);
  AIJ[2]->add(1., *AIJ[0]);

  //MatView ( (static_cast<PetscMatrix*> (AIJ[2]))->mat(), PETSC_VIEWER_STDOUT_SELF);

  AIJ[1] = SparseMatrix::build().release();
  AIJ[2]->get_transpose(*AIJ[1]);


  //AIJ[3] = NULL;

  *D = 0.01;
  DIJ[0]->matrix_set_diagonal_values(*D);
  DIJ[1] = NULL;
  DIJ[2] = NULL;
  DIJ[3] = DIJ[0];

  AIJ[3] = SparseMatrix::build().release();
  AIJ[3]->init(2, 2, DIJ);

  KK->init(2, 2, AIJ);



  for(unsigned k = 0; k < 4; k++) {
    delete PIJ[k];
    delete SIJ[k];
    if(AIJ[k] != NULL) delete AIJ[k];
  }

  delete DIJ[0];
  delete P;
  delete D;


}

void AssembleResMu(MultiLevelProblem& ml_prob) {

  BuildMuMatrix(ml_prob);

  LinearImplicitSystem* mlPdeSys   = &ml_prob.get_system< LinearImplicitSystem> ("mu");   // pointer to the linear implicit system named "Poisson"

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  // Pointers to the mesh (level) object and elem object in mesh (level).
  Mesh *msh = ml_prob._ml_msh->GetLevel(level);
  elem *el = msh->el;

  // Pointers to the multilevel solution, solution (level) and equation (level).
  MultiLevelSolution *mlSol = ml_prob._ml_sol;

  Solution* sol = mlSol->GetSolutionLevel(level);

  LinearEquationSolver *pdeSys = mlPdeSys->_LinSolver[level];

  NumericVector *RES = pdeSys->_RES;
  SparseMatrix *KK = pdeSys->_KK;

  RES->zero();
  //KK->draw();



  //unsigned  dim = msh->GetDimension();
  unsigned dim = 2;
  unsigned DIM = (parameter.surface) ? 3 : 2;

  std::vector < unsigned > indexDx(DIM);
  indexDx[0] = mlSol->GetIndex("Dx1");
  indexDx[1] = mlSol->GetIndex("Dx2");
  if(parameter.surface) indexDx[2] = mlSol->GetIndex("Dx3");
  unsigned solTypeDx = mlSol->GetSolutionType(indexDx[0]);


  unsigned indexW1 = mlSol->GetIndex("weight1");
  unsigned solType1 = mlSol->GetSolutionType(indexW1);

  std::vector< double > dof1;

  sol->_Sol[indexW1]->zero();

  double weight; // gauss point weight

  unsigned iproc = msh->processor_id();
  unsigned nprocs = msh->n_processors();

  std::vector < std::vector < double > > cX(2);

  unsigned vAngleIndex = mlSol->GetIndex("vAngle");
  unsigned vAngleType = mlSol->GetSolutionType(vAngleIndex);
  std::vector <double> vAngle;

  std::vector< double > phi;
  std::vector< double > dphidu;

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofs1  = msh->GetElementDofNumber(iel, solType1);
    unsigned nDofsDx  = msh->GetElementDofNumber(iel, solTypeDx);

    dof1.resize(nDofs1);

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofs1; i++) {
      dof1[i] = msh->GetSolutionDof(i, iel, solType1);
    }

    unsigned nvAngle = msh->GetElementDofNumber(iel, vAngleType);
    vAngle.resize(nvAngle);
    for(unsigned i = 0; i < nvAngle; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, vAngleType);
      vAngle[i] = (*sol->_Sol[vAngleIndex])(idof);
    }

    GetConformalCoordinates(msh, conformalType, iel, solTypeDx, vAngle, cX);

// *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solTypeDx]->GetGaussPointNumber(); ig++) {

      double weight; // gauss point weight
      msh->_finiteElement[ielGeom][solTypeDx]->Jacobian(cX, ig, weight, phi, dphidu);
      const double *phi1 = msh->_finiteElement[ielGeom][solType1]->GetPhi(ig);  // local test function

      for(unsigned i = 0; i < nDofs1; i++) {
        sol->_Sol[indexW1]->add(dof1[i], phi1[i] * weight);
      } // end phi_i loop

    } // end gauss point loop
  } //end element loop for each process*/
  sol->_Sol[indexW1]->close();



  ////////////////////////////////////////////////////

  std::vector < unsigned > indexMu(2);
  indexMu[0] = mlSol->GetIndex("mu1");
  indexMu[1] = mlSol->GetIndex("mu2");

  std::vector < unsigned > indexlMu(2);
  indexlMu[0] = mlSol->GetIndex("lmu1");
  indexlMu[1] = mlSol->GetIndex("lmu2");

  sol->_Sol[indexMu[0]]->zero();
  sol->_Sol[indexMu[1]]->zero();
  sol->_Sol[indexlMu[0]]->zero();
  sol->_Sol[indexlMu[1]]->zero();

  std::vector < unsigned > solPdeIndexMu(2);
  solPdeIndexMu[0] = mlPdeSys->GetSolPdeIndex("mu1");
  solPdeIndexMu[1] = mlPdeSys->GetSolPdeIndex("mu2");
  std::vector < unsigned > localToGlobalPDEMu[2];

  std::vector < std::vector < double > > xhat(DIM);
  std::vector < std::vector < double > > solx(DIM);

  for(unsigned k = 0; k < dim; k++) {
    sol->_Sol[indexMu[k]]->zero();
  }

  std::vector<double> phi_uv0;
  std::vector<double> phi_uv1;

  std::vector<double> area1;

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofs1  = msh->GetElementDofNumber(iel, solType1);
    unsigned nDofsDx  = msh->GetElementDofNumber(iel, solTypeDx);

    dof1.resize(nDofs1);
    area1.resize(nDofs1);
    localToGlobalPDEMu[0].resize(nDofs1);
    localToGlobalPDEMu[1].resize(nDofs1);

    for(int K = 0; K < DIM; K++) {
      xhat[K].resize(nDofsDx);
      solx[K].resize(nDofsDx);
    }

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofs1; i++) {
      dof1[i] = msh->GetSolutionDof(i, iel, solType1);
      area1[i] = (*sol->_Sol[indexW1])(dof1[i]);
      //std::cout<<area1[i]<<" ";
      localToGlobalPDEMu[0][i] = pdeSys->GetSystemDof(indexMu[0], solPdeIndexMu[0], i, iel);
      localToGlobalPDEMu[1][i] = pdeSys->GetSystemDof(indexMu[1], solPdeIndexMu[1], i, iel);
      std::cout << localToGlobalPDEMu[0][i] << " ";
      std::cout << localToGlobalPDEMu[1][i] << " ";
    }
    // local storage of coordinates
    for(unsigned i = 0; i < nDofsDx; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, solTypeDx);
      unsigned xDof  = msh->GetSolutionDof(i, iel, 2);
      for(unsigned K = 0; K < DIM; K++) {
        xhat[K][i] = (*msh->_topology->_Sol[K])(xDof) + (*sol->_SolOld[indexDx[K]])(idof);
        solx[K][i] = (*msh->_topology->_Sol[K])(xDof) + (*sol->_Sol[indexDx[K]])(idof);
      }
    }

    unsigned nvAngle = msh->GetElementDofNumber(iel, vAngleType);
    vAngle.resize(nvAngle);
    for(unsigned i = 0; i < nvAngle; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, vAngleType);
      vAngle[i] = (*sol->_Sol[vAngleIndex])(idof);
    }

    GetConformalCoordinates(msh, conformalType, iel, solTypeDx, vAngle, cX);

// *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solTypeDx]->GetGaussPointNumber(); ig++) {

      double weight; // gauss point weight
      msh->_finiteElement[ielGeom][solTypeDx]->Jacobian(cX, ig, weight, phi, dphidu);
      const double *phi1 = msh->_finiteElement[ielGeom][solType1]->GetPhi(ig);  // local test function

      // Initialize and compute fields at the Gauss points.
      double xhat_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      double solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};

      for(unsigned K = 0; K < DIM; K++) {
        for(int j = 0; j < dim; j++) {
          for(unsigned i = 0; i < nDofsDx; i++) {
            xhat_uv[K][j] += dphidu[i * dim + j] * xhat[K][i];
            solx_uv[K][j] += dphidu[i * dim + j] * solx[K][i];
          }
        }
      }

      // Compute the metric, metric determinant, and area element.
      double g[2][2] = {{0., 0.}, {0., 0.}};
      for(unsigned i = 0; i < dim; i++) {
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned K = 0; K < DIM; K++) {
            g[i][j] += xhat_uv[K][i] * xhat_uv[K][j];
          }
        }
      }
      double detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];

      double normal[3] = {0., 0., 1.};

      if(parameter.surface) {
        normal[0] = (xhat_uv[1][0] * xhat_uv[2][1] - xhat_uv[2][0] * xhat_uv[1][1]) / sqrt(detg);
        normal[1] = (xhat_uv[2][0] * xhat_uv[0][1] - xhat_uv[0][0] * xhat_uv[2][1]) / sqrt(detg);
        normal[2] = (xhat_uv[0][0] * xhat_uv[1][1] - xhat_uv[1][0] * xhat_uv[0][1]) / sqrt(detg);
      }

      boost::math::quaternion <double> N(0, normal[0], normal[1], normal[2]);

      boost::math::quaternion <double> du(0, solx_uv[0][0], solx_uv[1][0], solx_uv[2][0]);
      boost::math::quaternion <double> dv(0, solx_uv[0][1], solx_uv[1][1], solx_uv[2][1]);

      boost::math::quaternion <double> dup = du - N * dv;
      boost::math::quaternion <double> dum = du + N * dv;

      boost::math::quaternion <double> MU = (dum * conj(dup)) / norm(dup);

      double mu[2];

      mu[0] = MU.R_component_1();
      mu[1] = (MU * conj(N)).R_component_1();

      for(unsigned i = 0; i < nDofs1; i++) {
        for(unsigned k = 0; k < dim; k++) {
          //sol->_Sol[indexMu[k]]->add(dof1[i], mu[k] * phi1[i] * weight / area1[i]);
          RES->add(localToGlobalPDEMu[k][i], mu[k] * phi1[i] * weight / area1[i]);
        }
      } // end phi_i loop

    } // end gauss point loop
  } //end element loop for each process*/

  RES->close();

//   for(unsigned k = 0; k < dim; k++) {
//     sol->_Sol[indexMu[k]]->close();
//   }



}

double GetFinalMu(MultiLevelSolution & mlSol) {


  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1u;

  Solution* sol = mlSol.GetSolutionLevel(level);
  Mesh* msh = mlSol._mlMesh->GetLevel(level);

  unsigned dim = 2;

  std::vector < unsigned > indexMu(dim);
  indexMu[0] = mlSol.GetIndex("mu1");
  indexMu[1] = mlSol.GetIndex("mu2");
  unsigned solType1 = mlSol.GetSolutionType(indexMu[0]);

  unsigned iproc = msh->processor_id();
  unsigned nprocs = msh->n_processors();

  double MuNormAverageAfter;
  {
    double MuNormLocalSum = 0.;
    double muNormLocalMax = 0.;
    for(unsigned i = msh->_dofOffset[solType1][iproc]; i < msh->_dofOffset[solType1][iproc + 1]; i++) {
      double muNorm = sqrt(pow((*sol->_Sol[indexMu[0]])(i), 2) + pow((*sol->_Sol[indexMu[1]])(i), 2));
      MuNormLocalSum += muNorm;

      muNormLocalMax = (muNorm > muNormLocalMax) ? muNorm : muNormLocalMax;
    }
    double muNormMax;
    MPI_Allreduce(&muNormLocalMax, &muNormMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    std::cout << "smoothed mu infinity norm = " << muNormMax << std::endl;

    MPI_Allreduce(&MuNormLocalSum, &MuNormAverageAfter, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MuNormAverageAfter /= msh->_dofOffset[solType1][nprocs];

    std::cout << "smoothed mu average norm = " << MuNormAverageAfter << std::endl;
    std::cout << "relative difference = " << (muNormMax - MuNormAverageAfter) / MuNormAverageAfter << std::endl;

  }

  if(parameter.finalSmoothIsOn) {
    for(unsigned i = msh->_dofOffset[solType1][iproc]; i < msh->_dofOffset[solType1][iproc + 1]; i++) {

      double mu[2];
      for(unsigned k = 0; k < 2; k++) {
        mu[k] = (*sol->_Sol[indexMu[k]])(i);
      }

      double norm = sqrt(mu[0] * mu[0] + mu[1] * mu[1]);
      double cosTheta = mu[0] / norm;
      double sinTheta = mu[1] / norm;

      sol->_Sol[indexMu[0]]->set(i, MuNormAverageAfter * cosTheta);
      sol->_Sol[indexMu[1]]->set(i, MuNormAverageAfter * sinTheta);

    }
    for(unsigned k = 0; k < 2; k++) {
      sol->_Sol[indexMu[k]]->close();
    }
  }

}

void BuildPMatrix(MultiLevelSolution & mlSol) {

  double MuNormAverageBefore = EvaluateMu(mlSol);

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1u;

  Solution* sol = mlSol.GetSolutionLevel(level);
  Mesh* msh = mlSol._mlMesh->GetLevel(level);

  unsigned dim = 2;

  std::vector < unsigned > indexMu(dim);
  indexMu[0] = mlSol.GetIndex("mu1");
  indexMu[1] = mlSol.GetIndex("mu2");
  unsigned solType1 = mlSol.GetSolutionType(indexMu[0]);

  unsigned iproc = msh->processor_id();
  unsigned nprocs = msh->n_processors();

  unsigned vAngleIndex = mlSol.GetIndex("vAngle");
  unsigned vAngleType = mlSol.GetSolutionType(vAngleIndex);
  std::vector <double> vAngle;
  std::vector <double> eAngle;

  unsigned indexCntEdge = mlSol.GetIndex("cntEdge");

  unsigned elType = 3;
  unsigned faceType = 2;

  unsigned nel = msh->_dofOffset[elType][nprocs];
  unsigned nel_loc = msh->_dofOffset[elType][iproc + 1] - msh->_dofOffset[elType][iproc];

  unsigned nface = msh->_dofOffset[faceType][nprocs];
  unsigned nface_loc = msh->_dofOffset[faceType][iproc + 1] - msh->_dofOffset[faceType][iproc];


  std::vector < SparseMatrix* > PIJ(4);

  for(unsigned k = 0; k < 4; k++) {
    PIJ[k] = SparseMatrix::build().release();
    PIJ[k]->init(nface, nel, nface_loc, nel_loc, 10, 10);
  }

  std::vector < double > PIJl[4];

  std::vector< unsigned > irow;//loval to global mapping
  std::vector< unsigned > icolumn;//loval to global mapping

  for(int iel = msh->_dofOffset[elType][iproc]; iel < msh->_dofOffset[elType][iproc + 1]; iel++) {

    unsigned nvAngle = msh->GetElementDofNumber(iel, vAngleType);
    vAngle.resize(nvAngle);
    for(unsigned i = 0; i < nvAngle; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, vAngleType);
      vAngle[i] = (*sol->_Sol[vAngleIndex])(idof);
    }

    GetConformalAngles(msh, conformalType, iel, vAngle, eAngle);

    unsigned localDofOffset = msh->GetElementDofNumber(iel, 0);
    unsigned nFaces = msh->GetElementFaceNumber(iel);
    for(unsigned k = 0; k < 4; k++) {
      PIJl[k].resize(nFaces);
    }
    irow.resize(nFaces);
    icolumn.assign(1, iel);
    for(unsigned iface = 0; iface < nFaces; iface++) {

      irow[iface] = msh->GetSolutionDof(localDofOffset + iface, iel, faceType);

      double a = cos(eAngle[iface]);
      double b = sin(eAngle[iface]);

      PIJl[0][iface] = (a * a - b * b);
      PIJl[1][iface] = 2. * a * b;

      PIJl[2][iface] = - 2. * a * b;
      PIJl[3][iface] = (a * a - b * b);

      sol->_Sol[indexCntEdge]->add(irow[iface], 1);
    }
    for(unsigned k = 0; k < 4; k++) {
      PIJ[k]->add_matrix_blocked(PIJl[k], irow, icolumn);
    }
  }
  for(unsigned k = 0; k < 4; k++) {
    PIJ[k]->close();
  }
  sol->_Sol[indexCntEdge]->close();

  for(int iface = msh->_dofOffset[faceType][iproc]; iface < msh->_dofOffset[faceType][iproc + 1]; iface++) {
    double value = (*sol->_Sol[indexCntEdge])(iface);
    if(value > 0.5) {
      sol->_Sol[indexCntEdge]->set(iface, 1. / value);
    }
  }
  sol->_Sol[indexCntEdge]->close();
  for(unsigned k = 0; k < 4; k++) {
    MatDiagonalScale((static_cast<PetscMatrix*>(PIJ[k]))->mat(), (static_cast<PetscVector*>(sol->_Sol[indexCntEdge]))->vec(), NULL);
  }


  std::vector < SparseMatrix* > PIJt(4);
  for(unsigned i = 0; i < 4; i++) {
    PIJt[i] = SparseMatrix::build().release();
    PIJ[i]->get_transpose(*PIJt[i]);
  }

  NumericVector  *D = NumericVector::build().release();
  D->init(*sol->_Sol[indexCntEdge]);
  *D = 1.;

  SparseMatrix* I;
  I = SparseMatrix::build().release();
  I->init(nface, nface, nface_loc, nface_loc, 1, 1);
  I->matrix_set_diagonal_values(*D);

  SparseMatrix* temp;
  temp = SparseMatrix::build().release();



  PtP[0][0] = SparseMatrix::build().release();
  PtP[0][0]->matrix_ABC(*PIJt[0], *I , *PIJ[0], false);
  temp->matrix_ABC(*PIJt[2], *I , *PIJ[2], false);
  PtP[0][0]->add(1., *temp);

  PtP[0][1] = SparseMatrix::build().release();
  PtP[0][1]->matrix_ABC(*PIJt[0], *I , *PIJ[1], false);
  temp->matrix_ABC(*PIJt[2], *I , *PIJ[3], false);
  PtP[0][1]->add(1., *temp);

  PtP[1][0] = SparseMatrix::build().release();
  PtP[1][0]->matrix_ABC(*PIJt[1], *I , *PIJ[0], false);
  temp->matrix_ABC(*PIJt[3], *I , *PIJ[2], false);
  PtP[1][0]->add(1., *temp);

  PtP[1][1] = SparseMatrix::build().release();
  PtP[1][1]->matrix_ABC(*PIJt[1], *I , *PIJ[1], false);
  temp->matrix_ABC(*PIJt[3], *I , *PIJ[3], false);
  PtP[1][1]->add(1., *temp);



  for(unsigned k = 0; k < 4; k++) {
    delete PIJ[k];
    delete PIJt[k];
  }
  delete I;
  delete temp;

//   PtP[0][0]->draw();
//   PtP[0][1]->draw();
//   PtP[1][0]->draw();
//   PtP[1][1]->draw();








}


