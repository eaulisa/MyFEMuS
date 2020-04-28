
double EvaluateMu(MultiLevelSolution& mlSol) {
  
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

// Setting the reference elements to be equilateral triangles.
  std::vector < std::vector < double > > xT(2);
  xT[0].resize(3);
  xT[0] = { -0.5, 0.5, 0., 0., 0.25, -0.25, 0. };

  xT[1].resize(3);
  xT[1] = {0., 0., sqrt(3.) / 2., 0., sqrt(3.) / 4., sqrt(3.) / 4., sqrt(3.) / 6.};

  unsigned solENVNIndex = mlSol.GetIndex("ENVN");
  unsigned solENVNType = mlSol.GetSolutionType(solENVNIndex);

  std::vector<double> phi_uv0;
  std::vector<double> phi_uv1;

  std::vector< double > stdVectorPhi;
  std::vector< double > stdVectorPhi_uv;

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

// *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solTypeDx]->GetGaussPointNumber(); ig++) {

      const double *phix;  // local test function

      const double *phix_uv[dim]; // local test function first order partial derivatives

      double weight; // gauss point weight

      // Get Gauss point weight, test function, and first order derivatives.
      if(ielGeom == QUAD) {
        phix = msh->_finiteElement[ielGeom][solTypeDx]->GetPhi(ig);


        phix_uv[0] = msh->_finiteElement[ielGeom][solTypeDx]->GetDPhiDXi(ig);
        phix_uv[1] = msh->_finiteElement[ielGeom][solTypeDx]->GetDPhiDEta(ig);

        weight = msh->_finiteElement[ielGeom][solTypeDx]->GetGaussWeight(ig);
      }

      // Special adjustments for triangles.
      else {
        msh->_finiteElement[ielGeom][solTypeDx]->Jacobian(xT, ig, weight, stdVectorPhi, stdVectorPhi_uv);

        phix = &stdVectorPhi[0];
        phi_uv0.resize(nDofsDx);
        phi_uv1.resize(nDofsDx);


        for(unsigned i = 0; i < nDofsDx; i++) {
          phi_uv0[i] = stdVectorPhi_uv[i * dim];
          phi_uv1[i] = stdVectorPhi_uv[i * dim + 1];
        }

        phix_uv[0] = &phi_uv0[0];
        phix_uv[1] = &phi_uv1[0];

      }
      const double *phi1 = msh->_finiteElement[ielGeom][solType1]->GetPhi(ig);  // local test function

      // Initialize and compute values of x, Dx, NDx, x_uv at the Gauss points.
      double xhat_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      double solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};

      for(unsigned K = 0; K < DIM; K++) {
        for(int j = 0; j < dim; j++) {
          for(unsigned i = 0; i < nDofsDx; i++) {
            xhat_uv[K][j] += phix_uv[j][i] * xhat[K][i] ;
            solx_uv[K][j] += phix_uv[j][i] * solx[K][i];
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

      double normal[DIM] = {0., 0., 1.};

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


void UpdateMu(MultiLevelSolution& mlSol) {

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

  double angles[2][4] = {
    {0., 0.5 * M_PI, M_PI, 1.5 * M_PI}, // for square
    {0., 2. / 3. * M_PI, 4. / 3 * M_PI} // for equilateral triangle
  };

  std::vector < unsigned > indexMuEdge(dim);
  indexMuEdge[0] = mlSol.GetIndex("mu1Edge");
  indexMuEdge[1] = mlSol.GetIndex("mu2Edge");
  unsigned indexCntEdge = mlSol.GetIndex("cntEdge");
  unsigned solType2 = mlSol.GetSolutionType(indexMuEdge[0]);


  std::cout << "\nNumber of Smoothing steps = " << parameter.numberOfSmoothingSteps << std::endl;

  for(unsigned ismooth = 0; ismooth < parameter.numberOfSmoothingSteps; ismooth++) {

    for(unsigned k = 0; k < dim; k++) {
      sol->_Sol[indexMuEdge[k]]->zero();
    }
    sol->_Sol[indexCntEdge]->zero();

    for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

      short unsigned ielGeom = msh->GetElementType(iel);
      unsigned idx = (ielGeom == QUAD) ? 0 : 1;

      double mu[2];
      for(unsigned k = 0; k < 2; k++) {
        mu[k] = (*sol->_Sol[indexMu[k]])(iel);
      }

      unsigned nDofs0  = msh->GetElementDofNumber(iel, 0);
      for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {

        unsigned idof = msh->GetSolutionDof(nDofs0 + iface, iel, solType2);

        double a = cos(angles[idx][iface]);
        double b = sin(angles[idx][iface]);

        double mu0s = (a * a - b * b) * mu[0] + 2. * a * b * mu[1];
        double mu1s = (a * a - b * b) * mu[1] - 2. * a * b * mu[0];

        sol->_Sol[indexMuEdge[0]]->add(idof, mu0s);
        sol->_Sol[indexMuEdge[1]]->add(idof, mu1s);

        sol->_Sol[indexCntEdge]->add(idof, 1);
      }
    }
    for(unsigned k = 0; k < 2; k++) {
      sol->_Sol[indexMuEdge[k]]->close();
    }
    sol->_Sol[indexCntEdge]->close();



    for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

      short unsigned ielGeom = msh->GetElementType(iel);
      unsigned idx = (ielGeom == QUAD) ? 0 : 1;

      double mu[2] = {0., 0.};
      double cnt = 0.;

      unsigned nDofs0  = msh->GetElementDofNumber(iel, 0);
      for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {

        unsigned idof = msh->GetSolutionDof(nDofs0 + iface, iel, solType2);

        double mu0s = (*sol->_Sol[indexMuEdge[0]])(idof);
        double mu1s = (*sol->_Sol[indexMuEdge[1]])(idof);


        double a = cos(angles[idx][iface]);
        double b = sin(angles[idx][iface]);

        mu[0] += (a * a - b * b) * mu0s - 2. * a * b * mu1s;
        mu[1] += (a * a - b * b) * mu1s + 2. * a * b * mu0s;

        cnt += (*sol->_Sol[indexCntEdge])(idof);
      }

      for(unsigned k = 0; k < 2; k++) {
        sol->_Sol[indexMu[k]]->set(iel, mu[k] / cnt);
      }

    }
    for(unsigned k = 0; k < 2; k++) {
      sol->_Sol[indexMu[k]]->close();
    }
  }


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

}
