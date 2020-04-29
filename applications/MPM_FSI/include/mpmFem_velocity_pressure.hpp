#include "MultiLevelSolution.hpp"


using namespace femus;

double beta = 0.25;
double Gamma = 0.5;
double gravity[3] = {0., -9.81, 0.};
double scalingFactor1 = 1.e-6;
double scalingFactor2 = 1.e-10;
Line* linea;

void AssembleMPMSys (MultiLevelProblem& ml_prob) {

  clock_t AssemblyTime = 0;
  clock_t start_time, end_time;

  TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("MPM_FSI");
  const unsigned  level = my_nnlin_impl_sys.GetLevelToAssemble();
  MultiLevelSolution* ml_sol = ml_prob._ml_sol;
  Solution* mysolution = ml_sol->GetSolutionLevel (level);
  LinearEquationSolver* myLinEqSolver = my_nnlin_impl_sys._LinSolver[level];

  Mesh* mymsh = ml_prob._ml_msh->GetLevel (level);
  elem* myel = mymsh->el;
  SparseMatrix* myKK = myLinEqSolver->_KK;
  NumericVector* myRES =  myLinEqSolver->_RES;
  bool assembleMatrix = my_nnlin_impl_sys.GetAssembleMatrix();

  adept::Stack& s = FemusInit::_adeptStack;
  if (assembleMatrix) s.continue_recording();
  else s.pause_recording();

  const unsigned dim = mymsh->GetDimension();

  const unsigned max_size = static_cast< unsigned > (ceil (pow (3, dim)));   // conservative: based on line3, quad9, hex27

  unsigned iproc  = mymsh->processor_id();

  vector < double > phi;
  vector < double > phi_hat;
  vector < double > phi_tilde;
  vector < adept::adouble> gradphi;
  vector < double > gradphi_hat;
  vector < double > gradphi_tilde;

  phi.reserve (max_size);
  phi_hat.reserve (max_size);

  gradphi.reserve (max_size * dim);
  gradphi_hat.reserve (max_size * dim);

  vector <vector < adept::adouble> > vx (dim);
  vector <vector < double> > vx_hat (dim);

  vector< vector< adept::adouble > > SolVd (dim);
  vector< vector< double > > SolDd (dim);

  vector< vector< int > > dofsVAR (dim);

  vector< vector< double > > Rhs (dim);
  vector< vector< adept::adouble > > aRhs (dim);

  vector < int > dofsAll;

  vector < double > Jac;

  adept::adouble weight;
  double weight_hat = 0.;
  double weight_tilde = 0.;

  double density_MPM = ml_prob.parameters.get<Solid> ("SolidMPM").get_density();
  double E_MPM = ml_prob.parameters.get<Solid> ("SolidMPM").get_young_module();
  double mu_MPM = ml_prob.parameters.get<Solid> ("SolidMPM").get_lame_shear_modulus();
  double nu_MPM = ml_prob.parameters.get<Solid> ("SolidMPM").get_poisson_coeff();
  double lambda_MPM = ml_prob.parameters.get<Solid> ("SolidMPM").get_lame_lambda();
  double K_MPM = E_MPM / (3.* (1. - 2. * nu_MPM)); //bulk modulus

  double dt =  my_nnlin_impl_sys.GetIntervalTime();

  std::cout.precision (10);
  std::cout << "Scaling Factor 1 = " << scalingFactor1 << " Scaling Factor 2 = " << scalingFactor2 << std::endl;

  const char varname[7][5] = {"VX", "VY", "VZ", "DX", "DY", "DZ", "Mat"};
  vector <unsigned> indexSolV (dim);
  vector <unsigned> indexSolD (dim);
  vector <unsigned> indexPdeV (dim);
  unsigned solType = ml_sol->GetSolutionType (&varname[0][0]);

  unsigned indexSolM = ml_sol->GetIndex ("M");
  vector < bool > solidFlag;

  for (unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolV[ivar] = ml_sol->GetIndex (&varname[ivar][0]);
    indexSolD[ivar] = ml_sol->GetIndex (&varname[ivar + 3][0]);
    indexPdeV[ivar] = my_nnlin_impl_sys.GetSolPdeIndex (&varname[ivar][0]);
  }

  unsigned indexSolMat = ml_sol->GetIndex (&varname[6][0]);
  unsigned solTypeMat = ml_sol->GetSolutionType (&varname[6][0]);

  start_time = clock();

  if (assembleMatrix) myKK->zero();

  //line instances
  std::vector<unsigned> markerOffset = linea->GetMarkerOffset();
  unsigned markerOffset1 = markerOffset[iproc];
  unsigned markerOffset2 = markerOffset[iproc + 1];
  std::vector<Marker*> particles = linea->GetParticles();

  //BEGIN loop on elements (for soft stiffness and mass matrix)
  for (int iel = mymsh->_elementOffset[iproc]; iel < mymsh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielt = mymsh->GetElementType (iel);

    unsigned nDofsV = mymsh->GetElementDofNumber (iel, solType);
    unsigned nDofs = dim * nDofsV ;

    std::vector <int> sysDof (nDofs);

    solidFlag.resize (nDofsV);
    for (unsigned  k = 0; k < dim; k++) {
      SolVd[k].resize (nDofsV);
      SolDd[k].resize (nDofsV);
      vx[k].resize (nDofsV);
      vx_hat[k].resize (nDofsV);
    }

    for (unsigned  k = 0; k < dim; k++) {
      aRhs[k].resize (nDofsV);   //resize
      std::fill (aRhs[k].begin(), aRhs[k].end(), 0);
    }

    for (unsigned i = 0; i < nDofsV; i++) {
      unsigned idof = mymsh->GetSolutionDof (i, iel, solType);
      unsigned idofX = mymsh->GetSolutionDof (i, iel, 2);

      solidFlag[i] = ( (*mysolution->_Sol[indexSolM]) (idof) > 0.5 || mymsh->GetSolidMark (idof)) ? true : false;

      for (unsigned  k = 0; k < dim; k++) {
        SolVd[k][i] = (*mysolution->_Sol[indexSolV[k]]) (idof);
        SolDd[k][i] = (*mysolution->_Sol[indexSolD[k]]) (idof);
        sysDof[i + k * nDofsV] = myLinEqSolver->GetSystemDof (indexSolV[k], indexPdeV[k], i, iel);
        vx_hat[k][i] = (*mymsh->_topology->_Sol[k]) (idofX);
      }
    }

    if (assembleMatrix) s.new_recording();

    for (unsigned ig = 0; ig < mymsh->_finiteElement[ielt][solType]->GetGaussPointNumber(); ig++) {

      mymsh->_finiteElement[ielt][solType]->Jacobian (vx_hat, ig, weight_hat, phi_hat, gradphi_hat);

      vector < double > Xg (dim, 0);
      vector < vector < adept::adouble > > GradSolVgssHat (dim);
      vector < vector < double > > GradSolDgssHat (dim);

      for (unsigned  k = 0; k < dim; k++) {
        GradSolVgssHat[k].resize (dim);
        std::fill (GradSolVgssHat[k].begin(), GradSolVgssHat[k].end(), 0);
        GradSolDgssHat[k].assign (dim, 0.);
      }

      for (unsigned i = 0; i < nDofsV; i++) {
        for (unsigned j = 0; j < dim; j++) {
          Xg[j] += phi[i] * vx_hat[j][i];
          for (unsigned  k = 0; k < dim; k++) {
            GradSolVgssHat[k][j] += gradphi_hat[i * dim + j] * SolVd[k][i];
            GradSolDgssHat[k][j] += gradphi_hat[i * dim + j] * SolDd[k][i];
          }
        }
      }

      unsigned idofMat = mymsh->GetSolutionDof (0, iel, solTypeMat);
      double  MPMmaterial = (*mysolution->_Sol[indexSolMat]) (idofMat);
      double scalingFactor = 0;
      if (MPMmaterial < 5) scalingFactor = scalingFactor1;
      else if (MPMmaterial < 9) scalingFactor = scalingFactor2;

      for (unsigned i = 0; i < nDofsV; i++) {
        vector < adept::adouble > softStiffness (dim, 0.);

        for (unsigned k = 0; k < dim; k++) {
          for (unsigned  j = 0; j < dim; j++) {
            softStiffness[k]   +=  mu_MPM * gradphi_hat[i * dim + j] * (GradSolVgssHat[k][j] + GradSolVgssHat[j][k]) * dt; //TODO
//               softStiffness[k]   +=  mu_MPM * gradphi_hat[i * dim + j] * (GradSolDgssHat[k][j] + GradSolDgssHat[j][k]); //note GradSolDgssHat is note an Adept variable as GradSolVgssHat
          }
        }
        for (unsigned  k = 0; k < dim; k++) {
          if (MPMmaterial >= 2) { //soft stiffness contribution
            aRhs[k][i] += - softStiffness[k] * weight_hat * scalingFactor;
          }
          else if (!solidFlag[i]) { // mass matrix contribution
            aRhs[k][i] += phi_hat[i] * SolVd[k][i] * weight_hat;
          }
        }
      }
    }

    std::vector<double> Rhs (nDofs); //resize

    for (int i = 0; i < nDofsV; i++) {
      for (unsigned  k = 0; k < dim; k++) {
        Rhs[ i +  k * nDofsV ] = -aRhs[k][i].value();
      }
    }

    myRES->add_vector_blocked (Rhs, sysDof);

    if (assembleMatrix) {
      Jac.resize (nDofs * nDofs);

      for (unsigned  k = 0; k < dim; k++) {
        s.dependent (&aRhs[k][0], nDofsV);
      }

      for (unsigned  k = 0; k < dim; k++) {
        s.independent (&SolVd[k][0], nDofsV);
      }

      // get and store jacobian matrix (row-major)
      s.jacobian (&Jac[0], true);
      myKK->add_matrix_blocked (Jac, sysDof, sysDof);

      s.clear_independents();
      s.clear_dependents();
    }
  }
  //END loop on elements (for soft stiffness and mass matrix)


  unsigned ielOld = UINT_MAX;

  //BEGIN loop on particles (used as Gauss points)

  for (unsigned iMarker = markerOffset1; iMarker < markerOffset2; iMarker++) {
    //element of particle iMarker
    unsigned iel = particles[iMarker]->GetMarkerElement();
    if (iel != UINT_MAX) {
      short unsigned ielt;
      unsigned nDofsV;

      if (iel != ielOld) {

        ielt = mymsh->GetElementType (iel);
        nDofsV = mymsh->GetElementDofNumber (iel, solType);

        for (int i = 0; i < dim; i++) {
          dofsVAR[i].resize (nDofsV);
          SolVd[i].resize (nDofsV);
          SolDd[i].resize (nDofsV);
          aRhs[i].resize (nDofsV);
          vx[i].resize (nDofsV);
          vx_hat[i].resize (nDofsV);
        }
        dofsAll.resize (0);

        for (unsigned i = 0; i < nDofsV; i++) {
          unsigned idof = mymsh->GetSolutionDof (i, iel, solType);
          unsigned idofX = mymsh->GetSolutionDof (i, iel, 2);
          for (int j = 0; j < dim; j++) {
            SolVd[j][i] = (*mysolution->_Sol[indexSolV[j]]) (idof);
            SolDd[j][i] = (*mysolution->_Sol[indexSolD[j]]) (idof);

            dofsVAR[j][i] = myLinEqSolver->GetSystemDof (indexSolV[j], indexPdeV[j], i, iel);
            aRhs[j][i] = 0.;

            vx_hat[j][i] = (*mymsh->_topology->_Sol[j]) (idofX);
          }
        }

        for (int idim = 0; idim < dim; idim++) {
          dofsAll.insert (dofsAll.end(), dofsVAR[idim].begin(), dofsVAR[idim].end());
        }

        if (assembleMatrix) s.new_recording();

      }

      std::vector <double> xi = particles[iMarker]->GetMarkerLocalCoordinates();
      mymsh->_finiteElement[ielt][solType]->Jacobian (vx_hat, xi, weight_hat, phi_hat, gradphi_hat);

      std::vector <double> SolVpOld (dim);
      particles[iMarker]->GetMarkerVelocity (SolVpOld);

      std::vector <double> SolApOld (dim);
      particles[iMarker]->GetMarkerAcceleration (SolApOld);

      double mass = particles[iMarker]->GetMarkerMass();

      for (int j = 0; j < dim; j++) {
        for (unsigned inode = 0; inode < nDofsV; inode++) {
          vx[j][inode] = vx_hat[j][inode] + SolVd[j][inode] * dt; //TODO
//           vx[j][inode] = vx_hat[j][inode] + SolDd[j][inode];


// std::cout <<"SolDd["<<j<<"]["<<inode<<"]= " << SolDd[j][inode] << " , " << "SolVd["<<j<<"]["<<inode<<"] * dt = " << SolVd[j][inode] * dt << std::endl;

        }
      }

      mymsh->_finiteElement[ielt][solType]->Jacobian (vx, xi, weight, phi, gradphi);

      //BEGIN evaluates SolVp and GradSolVpHat at the particle iMarker
      vector<adept::adouble> SolVp (dim, 0.);
      vector<vector < adept::adouble > > GradSolVpHat (dim);
      vector<vector < double > > GradSolDpHat (dim);
      for (int i = 0; i < dim; i++) {
        GradSolVpHat[i].assign (dim, 0.);
        GradSolDpHat[i].assign (dim, 0.);
      }

      for (int i = 0; i < dim; i++) {
        for (unsigned inode = 0; inode < nDofsV; inode++) {
          SolVp[i] += phi[inode] * SolVd[i][inode];
          for (int j = 0; j < dim; j++) {
            GradSolVpHat[i][j] +=  gradphi_hat[inode * dim + j] * SolVd[i][inode];
            GradSolDpHat[i][j] +=  gradphi_hat[inode * dim + j] * SolDd[i][inode];
          }
        }
      }
      //END evaluates SolVp and GradSolVpHat at the particle iMarker


      //BEGIN computation of the Cauchy Stress
      std::vector < std::vector < double > > FpOld;
      FpOld = particles[iMarker]->GetDeformationGradient();

      adept::adouble FpNew[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
      adept::adouble F[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      adept::adouble B[3][3];
      adept::adouble Id2th[3][3] = {{ 1., 0., 0.}, { 0., 1., 0.}, { 0., 0., 1.}};
      adept::adouble Cauchy[3][3];

      for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
          FpNew[i][j] += GradSolVpHat[i][j] * dt; //TODO
//           FpNew[i][j] += GradSolDpHat[i][j]; // note: GradSolDpHat cannot be an Adept variable like GradSolVpHat
        }
      }

      for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
          for (int k = 0; k < dim; k++) {
            F[i][j] += FpNew[i][k] * FpOld[k][j];
          }
        }
      }

      if (dim == 2) F[2][2] = 1.;

      adept::adouble J_hat =  F[0][0] * F[1][1] * F[2][2] + F[0][1] * F[1][2] * F[2][0] + F[0][2] * F[1][0] * F[2][1]
                              - F[2][0] * F[1][1] * F[0][2] - F[2][1] * F[1][2] * F[0][0] - F[2][2] * F[1][0] * F[0][1];

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          B[i][j] = 0.;

          for (int k = 0; k < 3; k++) {
            //left Cauchy-Green deformation tensor or Finger tensor (B = F*F^T)
            B[i][j] += F[i][k] * F[j][k];
          }
        }
      }

      adept::adouble I1_B = B[0][0] + B[1][1] + B[2][2];

      double mu = mu_MPM;
      double lambda = lambda_MPM;

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          //Cauchy[i][j] = mu_MPM * (B[i][j] - I1_B * Id2th[i][j] / 3.) / pow(J_hat, 5. / 3.)
          //             + K_MPM * (J_hat - 1.) * Id2th[i][j];  //Generalized Neo-Hookean solid, in Allan-Bower's book, for rubbers with very limited compressibility and K >> mu

          Cauchy[i][j] = lambda * log (J_hat) / J_hat * Id2th[i][j] + mu / J_hat * (B[i][j] - Id2th[i][j]); //alternative formulation


        }
      }
      //END computation of the Cauchy Stress

      //BEGIN MPM contributions
      for (unsigned i = 0; i < nDofsV; i++) {
        adept::adouble CauchyDIR[3] = {0., 0., 0.};

        for (int idim = 0.; idim < dim; idim++) {
          for (int jdim = 0.; jdim < dim; jdim++) {
            CauchyDIR[idim] += gradphi[i * dim + jdim] * Cauchy[idim][jdim];
          }
        }

        for (int idim = 0; idim < dim; idim++) {
          aRhs[idim][i] += (phi[i] * gravity[idim] - J_hat * CauchyDIR[idim] / density_MPM
                            - phi[i] * ( (SolVp[idim] - SolVpOld[idim]) / (Gamma * dt) + (1. - 1. / Gamma) * SolApOld[idim])
                           ) * mass;
        }
      }
      //END MPM contributions


      if (iMarker == markerOffset2 - 1 || iel != particles[iMarker + 1]->GetMarkerElement()) {

        for (unsigned i = 0; i < dim; i++) {
          Rhs[i].resize (nDofsV);

          for (int j = 0; j < nDofsV; j++) {
            Rhs[i][j] = -aRhs[i][j].value();
          }
        }

        for (int i = 0; i < dim; i++) {
          myRES->add_vector_blocked (Rhs[i], dofsVAR[i]);
        }

        if (assembleMatrix) {
          for (int i = 0; i < dim; i++) {
            s.dependent (&aRhs[i][0], nDofsV);
            s.independent (&SolVd[i][0], nDofsV);
          }

          Jac.resize ( (dim * nDofsV) * (dim * nDofsV));

          s.jacobian (&Jac[0], true);

          myKK->add_matrix_blocked (Jac, dofsAll, dofsAll);
          s.clear_independents();
          s.clear_dependents();
        }
      }

      ielOld = iel;
    }
    else {
      break;
    }
  }
  //END loop on particles

  myRES->close();
  mysolution->_Sol[indexSolMat]->close();

  if (assembleMatrix) {
    myKK->close();
  }

  // *************************************
  end_time = clock();
  AssemblyTime += (end_time - start_time);
  // ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

}

void GridToParticlesProjection (MultiLevelProblem & ml_prob, Line & linea) {

  clock_t AssemblyTime = 0;
  clock_t start_time, end_time;

  TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("MPM_FSI");
  const unsigned  level = my_nnlin_impl_sys.GetLevelToAssemble();
  MultiLevelSolution* ml_sol = ml_prob._ml_sol;
  Solution* mysolution = ml_sol->GetSolutionLevel (level);

  Mesh* mymsh = ml_prob._ml_msh->GetLevel (level);
  elem* myel = mymsh->el;

  double dt =  my_nnlin_impl_sys.GetIntervalTime();
  const unsigned dim = mymsh->GetDimension();

  unsigned iproc  = mymsh->processor_id();

  vector< vector < double > > SolVd (dim);
  vector< vector < double > > GradSolVpHat (dim);

  for (int i = 0; i < dim; i++) {
    GradSolVpHat[i].resize (dim);
  }

  vector < double > phi_hat;
  vector < double > gradphi_hat;
  vector < double > nablaphi_hat;

  vector <vector < double> > vx_hat (dim);

  double weight;

  const char varname[4][3] = {"VX", "VY", "VZ", "gM"};
  vector <unsigned> indexSolV (dim);
  unsigned solType = ml_sol->GetSolutionType (&varname[0][0]);

  unsigned indexGridMass = ml_sol->GetIndex (&varname[3][0]);
  for (unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolV[ivar] = ml_sol->GetIndex (&varname[ivar][0]);
  }

  std::vector<unsigned> markerOffset = linea.GetMarkerOffset();
  unsigned markerOffset1 = markerOffset[iproc];
  unsigned markerOffset2 = markerOffset[iproc + 1];
  std::vector<Marker*> particles = linea.GetParticles();

  unsigned ielOld = UINT_MAX;

  //BEGIN loop on particles
  for (unsigned iMarker = markerOffset1; iMarker < markerOffset2; iMarker++) {

    unsigned iel = particles[iMarker]->GetMarkerElement();

    if (iel != UINT_MAX) {

      short unsigned ielt;
      unsigned nve;

      if (iel != ielOld) {

        ielt = mymsh->GetElementType (iel);
        nve = mymsh->GetElementDofNumber (iel, solType);

        for (int i = 0; i < dim; i++) {
          SolVd[i].resize (nve);
          vx_hat[i].resize (nve);
        }

        for (unsigned inode = 0; inode < nve; inode++) {
          unsigned idof = mymsh->GetSolutionDof (inode, iel, solType);
          unsigned idofX = mymsh->GetSolutionDof (inode, iel, 2);

          for (int i = 0; i < dim; i++) {
            SolVd[i][inode] = (*mysolution->_Sol[indexSolV[i]]) (idof);
            vx_hat[i][inode] = (*mymsh->_topology->_Sol[i]) (idofX);
          }
        }


      }
      std::vector <double> xi = particles[iMarker]->GetMarkerLocalCoordinates();

      mymsh->_finiteElement[ielt][solType]->Jacobian (vx_hat, xi, weight, phi_hat, gradphi_hat, nablaphi_hat);

      std::vector <double> particleVelOld (dim);
      particles[iMarker]->GetMarkerVelocity (particleVelOld);

      std::vector <double> particleAccOld (dim);
      particles[iMarker]->GetMarkerAcceleration (particleAccOld);

      std::vector <double> particleDispOld (dim);
      particles[iMarker]->GetMarkerDisplacement (particleDispOld);

      std::vector <double> particleVel (dim, 0.);
      for (int i = 0; i < dim; i++) {
        for (unsigned inode = 0; inode < nve; inode++) {
          particleVel[i] += phi_hat[inode] * SolVd[i][inode];
        }
      }

      particles[iMarker]->SetMarkerVelocity (particleVel);

      std::vector <double> particleAcc (dim);
      std::vector <double> particleDisp (dim);
      for (unsigned i = 0; i < dim; i++) {
        particleAcc[i] = (particleVel[i] - particleVelOld[i]) / (Gamma * dt) + (1. - 1. / Gamma) * particleAccOld[i];
        particleDisp[i] = particleDispOld[i] + dt * particleVelOld[i] + 0.5 * dt * dt * (1. - 2. * beta) * particleAccOld[i] + dt * dt * beta * particleAcc[i];
      }

      particles[iMarker]->SetMarkerAcceleration (particleAcc);
      particles[iMarker]->SetMarkerDisplacement (particleDisp);
      particles[iMarker]->UpdateParticleCoordinates();

      for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
          GradSolVpHat[i][j] = 0.;
          for (unsigned inode = 0; inode < nve; inode++) {
            GradSolVpHat[i][j] +=  gradphi_hat[inode * dim + j] * SolVd[i][inode];
          }
        }
      }

      std::vector < std::vector < double > > FpOld;
      FpOld = particles[iMarker]->GetDeformationGradient();

      double FpNew[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
      std::vector < std::vector < double > > Fp (dim);

      for (unsigned i = 0; i < dim; i++) {
        for (unsigned j = 0; j < dim; j++) {
          FpNew[i][j] += GradSolVpHat[i][j] * dt;
        }
      }

      for (unsigned i = 0; i < dim; i++) {
        Fp[i].resize (dim);
        for (unsigned j = 0; j < dim; j++) {
          Fp[i][j] = 0.;
          for (unsigned k = 0; k < dim; k++) {
            Fp[i][j] += FpNew[i][k] * FpOld[k][j];
          }
        }
      }

      particles[iMarker]->SetDeformationGradient (Fp);

      ielOld = iel;
    }
    else {
      break;
    }
  }
  //END loop on particles

  linea.UpdateLineMPM();

  linea.GetParticlesToGridMaterial();
}

void ParticlesToGridProjection (MultiLevelProblem & ml_prob, Line & linea) {

  clock_t AssemblyTime = 0;
  clock_t start_time, end_time;

  TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("MPM_FSI");
  const unsigned  level = my_nnlin_impl_sys.GetLevelToAssemble();
  MultiLevelSolution* ml_sol = ml_prob._ml_sol;
  Solution* mysolution = ml_sol->GetSolutionLevel (level);

  Mesh* mymsh = ml_prob._ml_msh->GetLevel (level);
  elem* myel = mymsh->el;

  double dt =  my_nnlin_impl_sys.GetIntervalTime();
  const unsigned dim = mymsh->GetDimension();

  unsigned iproc  = mymsh->processor_id();

  vector < double > phi_hat;
  vector < double > gradphi_hat;
  vector < double > nablaphi_hat;

  vector <vector < double> > vx_hat (dim);

  double weight;

  const char varname[7][3] = {"VX", "VY", "VZ", "DX", "DY", "DZ", "gM"};
  vector <unsigned> indexSolV (dim);
  vector <unsigned> indexSolD (dim);
  unsigned solType = ml_sol->GetSolutionType (&varname[0][0]);

  unsigned indexGridMass = ml_sol->GetIndex (&varname[6][0]);
  for (unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolV[ivar] = ml_sol->GetIndex (&varname[ivar][0]);
    indexSolD[ivar] = ml_sol->GetIndex (&varname[ivar + 3][0]);
  }

  std::vector<unsigned> markerOffset = linea.GetMarkerOffset();
  unsigned markerOffset1 = markerOffset[iproc];
  unsigned markerOffset2 = markerOffset[iproc + 1];
  std::vector<Marker*> particles = linea.GetParticles();

  unsigned ielOld = UINT_MAX;

  //BEGIN loop on elements to set grid mass and displacement to zero
  for (unsigned idof = mymsh->_dofOffset[solType][iproc]; idof < mymsh->_dofOffset[solType][iproc + 1]; idof++) {
    mysolution->_Sol[indexGridMass]->set (idof, 0.);
    for (int i = 0; i < dim; i++) {
      mysolution->_Sol[indexSolD[i]]->set (idof, 0.);
    }
  }
  mysolution->_Sol[indexGridMass]->close();
  for (int i = 0; i < dim; i++) {
    mysolution->_Sol[indexSolD[i]]->close();
  }
  //END loop on elements to set grid mass and displacement to zero


  //BEGIN loop on particles to build grid mass
  for (unsigned iMarker = markerOffset1; iMarker < markerOffset2; iMarker++) {

    unsigned iel = particles[iMarker]->GetMarkerElement();

    if (iel != UINT_MAX) {

      short unsigned ielt;
      unsigned nve;

      if (iel != ielOld) {

        ielt = mymsh->GetElementType (iel);
        nve = mymsh->GetElementDofNumber (iel, solType);

        for (int i = 0; i < dim; i++) {
          vx_hat[i].resize (nve);
        }

        for (unsigned inode = 0; inode < nve; inode++) {
          unsigned idofX = mymsh->GetSolutionDof (inode, iel, 2);

          for (int i = 0; i < dim; i++) {
            vx_hat[i][inode] = (*mymsh->_topology->_Sol[i]) (idofX);
          }
        }

      }

      std::vector <double> xi = particles[iMarker]->GetMarkerLocalCoordinates();

      mymsh->_finiteElement[ielt][solType]->Jacobian (vx_hat, xi, weight, phi_hat, gradphi_hat, nablaphi_hat);

      double particleMass = particles[iMarker]->GetMarkerMass();

      for (unsigned inode = 0; inode < nve; inode++) {
        unsigned idof = mymsh->GetSolutionDof (inode, iel, solType);
        double massGridLocal = phi_hat[inode] * particleMass;
        mysolution->_Sol[indexGridMass]->add (idof, massGridLocal);
      }

      ielOld = iel;

    }

    else {
      break;
    }
  }

  mysolution->_Sol[indexGridMass]->close();
//END loop on particles to build grid mass

//BEGIN loop on particles to build grid displacement
  for (unsigned iMarker = markerOffset1; iMarker < markerOffset2; iMarker++) {

    unsigned iel = particles[iMarker]->GetMarkerElement();

    if (iel != UINT_MAX) {

      short unsigned ielt;
      unsigned nve;

      if (iel != ielOld) {

        ielt = mymsh->GetElementType (iel);
        nve = mymsh->GetElementDofNumber (iel, solType);

        for (int i = 0; i < dim; i++) {
          vx_hat[i].resize (nve);
        }

        for (unsigned inode = 0; inode < nve; inode++) {
          unsigned idofX = mymsh->GetSolutionDof (inode, iel, 2);

          for (int i = 0; i < dim; i++) {
            vx_hat[i][inode] = (*mymsh->_topology->_Sol[i]) (idofX);
          }
        }

      }

      std::vector <double> xi = particles[iMarker]->GetMarkerLocalCoordinates();

      mymsh->_finiteElement[ielt][solType]->Jacobian (vx_hat, xi, weight, phi_hat, gradphi_hat, nablaphi_hat);

      double particleMass = particles[iMarker]->GetMarkerMass();

      std::vector <double> particleDisp (dim);
      particles[iMarker]->GetMarkerDisplacement (particleDisp);

      for (unsigned inode = 0; inode < nve; inode++) {
        unsigned idof = mymsh->GetSolutionDof (inode, iel, solType);
        double inodeGridMass = (*mysolution->_Sol[indexGridMass]) (idof);
        for (int i = 0; i < dim; i++) {
          double dispLocal = (fabs (inodeGridMass) > 1.e-7) ? (phi_hat[inode] * particleMass * particleDisp[i]) / inodeGridMass : 0.;
          mysolution->_Sol[indexSolD[i]]->add (idof, dispLocal);
        }
      }

      ielOld = iel;

    }

    else {
      break;
    }
  }

  for (int i = 0; i < dim; i++) {
    mysolution->_Sol[indexSolD[i]]->close();
  }
//END loop on particles to build grid displacement

}


unsigned getNumberOfLayers (const double &a, const double &fac) {
  double da = 1. / fac;
  double b =  da;
  unsigned n = 1;

  while (b < a) {
    da /= fac;
    b += da;
    n++;
    if (n >= 100) {
      std::cout << "Error: number of layer is unbounded, try with a smaller factor\n";
      abort();
    }
  }
  return n;
}


