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

  vector< vector< adept::adouble > > SolDd (dim);
  vector< vector< adept::adouble > > SolVd (dim);
  vector< vector< double > > SolVdOld (dim);
  vector< vector< double > > SolAdOld (dim);

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

  const char varname[13][6] = {"DX", "DY", "DZ", "VX", "VY", "VZ", "AXOld", "AYOld", "AZOld", "Mat", "VXOld", "VYOld", "VZOld"};
  vector <unsigned> indexSolD (dim);
  vector <unsigned> indexSolV (dim);
  vector <unsigned> indexSolVOld (dim);
  vector <unsigned> indexSolAOld (dim);
  vector <unsigned> indexPdeV (dim);
  unsigned solType = ml_sol->GetSolutionType (&varname[3][0]);

  unsigned indexSolM = ml_sol->GetIndex ("M");
  vector < bool > solidFlag;

  for (unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolD[ivar] = ml_sol->GetIndex (&varname[ivar][0]);
    indexSolV[ivar] = ml_sol->GetIndex (&varname[ivar + 3][0]);
    indexSolVOld[ivar] = ml_sol->GetIndex (&varname[ivar + 10][0]);
    indexSolAOld[ivar] = ml_sol->GetIndex (&varname[ivar + 6][0]);
    indexPdeV[ivar] = my_nnlin_impl_sys.GetSolPdeIndex (&varname[ivar + 3][0]);
  }

  unsigned indexSolMat = ml_sol->GetIndex (&varname[9][0]);
  unsigned solTypeMat = ml_sol->GetSolutionType (&varname[9][0]);

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
      SolDd[k].resize (nDofsV);
      SolVd[k].resize (nDofsV);
      SolVdOld[k].resize (nDofsV);
      SolAdOld[k].resize (nDofsV);
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
        SolVdOld[k][i] = (*mysolution->_Sol[indexSolVOld[k]]) (idof);
        SolAdOld[k][i] = (*mysolution->_Sol[indexSolAOld[k]]) (idof);
        sysDof[i + k * nDofsV] = myLinEqSolver->GetSystemDof (indexSolV[k], indexPdeV[k], i, iel);
        vx_hat[k][i] = (*mymsh->_topology->_Sol[k]) (idofX);
      }
    }

    if (assembleMatrix) s.new_recording();

    for (unsigned i = 0; i < nDofsV; i++) {
      for (unsigned  k = 0; k < dim; k++) {
        //NOTE
//         SolDd[k][i] = SolVdOld[k][i] * dt + (SolVd[k][i] - SolVdOld[k][i]) * beta * dt / Gamma + SolAdOld[k][i] * dt * dt * (0.5 - beta / Gamma);
        SolDd[k][i] = SolVd[k][i] * dt;
      }
    }

    for (unsigned ig = 0; ig < mymsh->_finiteElement[ielt][solType]->GetGaussPointNumber(); ig++) {

      mymsh->_finiteElement[ielt][solType]->Jacobian (vx_hat, ig, weight_hat, phi_hat, gradphi_hat);

      vector < vector < adept::adouble > > GradSolDgssHat (dim);

      for (unsigned  k = 0; k < dim; k++) {
        GradSolDgssHat[k].resize (dim);
        std::fill (GradSolDgssHat[k].begin(), GradSolDgssHat[k].end(), 0);
      }

      for (unsigned i = 0; i < nDofsV; i++) {
        for (unsigned j = 0; j < dim; j++) {
          for (unsigned  k = 0; k < dim; k++) {
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
            softStiffness[k]   +=  mu_MPM * gradphi_hat[i * dim + j] * (GradSolDgssHat[k][j] + GradSolDgssHat[j][k]);
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
          SolDd[i].resize (nDofsV);
          SolVd[i].resize (nDofsV);
          SolVdOld[i].resize (nDofsV);
          SolAdOld[i].resize (nDofsV);
          aRhs[i].resize (nDofsV);
          vx[i].resize (nDofsV);
          vx_hat[i].resize (nDofsV);
        }
        dofsAll.resize (0);

        for (unsigned i = 0; i < nDofsV; i++) {
          unsigned idof = mymsh->GetSolutionDof (i, iel, solType);
          unsigned idofX = mymsh->GetSolutionDof (i, iel, 2);
          for (int k = 0; k < dim; k++) {
            SolVd[k][i] = (*mysolution->_Sol[indexSolV[k]]) (idof);
            SolVdOld[k][i] = (*mysolution->_Sol[indexSolVOld[k]]) (idof);
            SolAdOld[k][i] = (*mysolution->_Sol[indexSolAOld[k]]) (idof);
            dofsVAR[k][i] = myLinEqSolver->GetSystemDof (indexSolV[k], indexPdeV[k], i, iel);
            aRhs[k][i] = 0.;
            vx_hat[k][i] = (*mymsh->_topology->_Sol[k]) (idofX);
          }
        }

        for (int idim = 0; idim < dim; idim++) {
          dofsAll.insert (dofsAll.end(), dofsVAR[idim].begin(), dofsVAR[idim].end());
        }

        if (assembleMatrix) s.new_recording();

        for (unsigned i = 0; i < nDofsV; i++) {
          for (unsigned  k = 0; k < dim; k++) {
            //NOTE
//             SolDd[k][i] = SolVdOld[k][i] * dt + (SolVd[k][i] - SolVdOld[k][i]) * beta * dt / Gamma + SolAdOld[k][i] * dt * dt * (0.5 - beta / Gamma);
            SolDd[k][i] = SolVd[k][i] * dt;
          }
        }

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
          vx[j][inode] = vx_hat[j][inode] + SolDd[j][inode];
        }
      }

      mymsh->_finiteElement[ielt][solType]->Jacobian (vx, xi, weight, phi, gradphi);

      //BEGIN evaluates SolVp, SolDp and GradSolDpHat at the particle iMarker
      vector<adept::adouble> SolDp (dim, 0.);
      vector<adept::adouble> SolVp (dim, 0.);
      vector<vector < adept::adouble > > GradSolDpHat (dim);
      for (int i = 0; i < dim; i++) {
        GradSolDpHat[i].assign (dim, 0.);
      }

      for (int i = 0; i < dim; i++) {
        for (unsigned inode = 0; inode < nDofsV; inode++) {
          SolDp[i] += phi[inode] * SolDd[i][inode];
          SolVp[i] += phi[inode] * SolVd[i][inode];
          for (int j = 0; j < dim; j++) {
            GradSolDpHat[i][j] +=  gradphi_hat[inode * dim + j] * SolDd[i][inode];
          }
        }
      }
      //END evaluates SolVp, SolDp and GradSolDpHat at the particle iMarker


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
          FpNew[i][j] += GradSolDpHat[i][j];
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
          //NOTE
//           aRhs[idim][i] += (phi[i] * gravity[idim] - J_hat * CauchyDIR[idim] / density_MPM
//                             - phi[i] * ( (SolVp[idim] - SolVpOld[idim]) / (Gamma * dt) + (1. - 1. / Gamma) * SolApOld[idim])
//                            ) * mass;
          aRhs[idim][i] += (phi[i] * gravity[idim] - J_hat * CauchyDIR[idim] / density_MPM
                            - phi[i] * (1. / (beta * dt * dt) * SolDp[idim] - 1. / (beta * dt) * SolVpOld[idim] - (1. - 2.* beta) / (2. * beta) * SolApOld[idim])
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

  TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("MPM_FSI");
  const unsigned  level = my_nnlin_impl_sys.GetLevelToAssemble();
  MultiLevelSolution* ml_sol = ml_prob._ml_sol;
  Solution* mysolution = ml_sol->GetSolutionLevel (level);

  Mesh* mymsh = ml_prob._ml_msh->GetLevel (level);
  elem* myel = mymsh->el;

  double dt =  my_nnlin_impl_sys.GetIntervalTime();
  const unsigned dim = mymsh->GetDimension();

  unsigned iproc  = mymsh->processor_id();

  vector< vector < double > > SolDd (dim);
  vector< vector < double > > SolVd (dim);
  vector< vector < double > > SolVdOld (dim);
  vector< vector < double > > SolAdOld (dim);
  vector< vector < double > > GradSolDpHat (dim);

  for (int i = 0; i < dim; i++) {
    GradSolDpHat[i].resize (dim);
  }

  vector < double > phi_hat;
  vector < double > gradphi_hat;
  vector < double > nablaphi_hat;

  vector <vector < double> > vx_hat (dim);

  double weight;

  const char varname[12][6] = {"DX", "DY", "DZ", "VX", "VY", "VZ", "AXOld", "AYOld", "AZOld", "VXOld", "VYOld", "VZOld"};
  vector <unsigned> indexSolD (dim);
  vector <unsigned> indexSolV (dim);
  vector <unsigned> indexSolVOld (dim);
  vector <unsigned> indexSolAOld (dim);
  unsigned solType = ml_sol->GetSolutionType (&varname[3][0]);

  for (unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolD[ivar] = ml_sol->GetIndex (&varname[ivar][0]);
    indexSolV[ivar] = ml_sol->GetIndex (&varname[ivar + 3][0]);
    indexSolVOld[ivar] = ml_sol->GetIndex (&varname[ivar + 9][0]);
    indexSolAOld[ivar] = ml_sol->GetIndex (&varname[ivar + 6][0]);
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
          SolDd[i].resize (nve);
          SolVd[i].resize (nve);
          SolVdOld[i].resize (nve);
          SolAdOld[i].resize (nve);
          vx_hat[i].resize (nve);
        }

        for (unsigned inode = 0; inode < nve; inode++) {
          unsigned idof = mymsh->GetSolutionDof (inode, iel, solType);
          unsigned idofX = mymsh->GetSolutionDof (inode, iel, 2);

          for (int i = 0; i < dim; i++) {
            SolVd[i][inode] = (*mysolution->_Sol[indexSolV[i]]) (idof);
            SolVdOld[i][inode] = (*mysolution->_Sol[indexSolVOld[i]]) (idof);
            SolAdOld[i][inode] = (*mysolution->_Sol[indexSolAOld[i]]) (idof);
            vx_hat[i][inode] = (*mymsh->_topology->_Sol[i]) (idofX);
          }
        }

        for (unsigned inode = 0; inode < nve; inode++) {
          for (unsigned  i = 0; i < dim; i++) {
            //NOTE
//             SolDd[i][inode] = SolVdOld[i][inode] * dt + (SolVd[i][inode] - SolVdOld[i][inode]) * beta * dt / Gamma + SolAdOld[i][inode] * dt * dt * (0.5 - beta / Gamma);
            SolDd[i][inode] = SolVd[i][inode] * dt;
          }
        }


      }
      std::vector <double> xi = particles[iMarker]->GetMarkerLocalCoordinates();

      mymsh->_finiteElement[ielt][solType]->Jacobian (vx_hat, xi, weight, phi_hat, gradphi_hat, nablaphi_hat);

      //NOTE
      //BEGIN

//       std::vector <double> particleVelOld (dim);
//       particles[iMarker]->GetMarkerVelocity (particleVelOld);
// 
//       std::vector <double> particleAccOld (dim);
//       particles[iMarker]->GetMarkerAcceleration (particleAccOld);
// 
//       std::vector <double> particleDispOld (dim);
//       particles[iMarker]->GetMarkerDisplacement (particleDispOld);
// 
//       std::vector <double> particleVel (dim, 0.);
//       for (int i = 0; i < dim; i++) {
//         for (unsigned inode = 0; inode < nve; inode++) {
//           particleVel[i] += phi_hat[inode] * SolVd[i][inode];
//         }
//       }
// 
//       particles[iMarker]->SetMarkerVelocity (particleVel);
// 
//       std::vector <double> particleAcc (dim);
//       std::vector <double> particleDisp (dim);
//       for (unsigned i = 0; i < dim; i++) {
//         particleAcc[i] = (particleVel[i] - particleVelOld[i]) / (Gamma * dt) + (1. - 1. / Gamma) * particleAccOld[i];
//         particleDisp[i] = particleDispOld[i] + dt * particleVelOld[i] + 0.5 * dt * dt * (1. - 2. * beta) * particleAccOld[i] + dt * dt * beta * particleAcc[i];
//       }
// 
//       particles[iMarker]->SetMarkerAcceleration (particleAcc);
//       particles[iMarker]->SetMarkerDisplacement (particleDisp);
//       particles[iMarker]->UpdateParticleCoordinates();

      //END

      //NOTE
      //BEGIN will be removed eventually

      std::vector <double> particleVelOld (dim);
      particles[iMarker]->GetMarkerVelocity (particleVelOld);

      std::vector <double> particleAccOld (dim);
      particles[iMarker]->GetMarkerAcceleration (particleAccOld);

      std::vector <double> particleDisp (dim, 0.);
      //update displacement and acceleration
      for (int i = 0; i < dim; i++) {
        for (unsigned inode = 0; inode < nve; inode++) {
          particleDisp[i] += phi_hat[inode] * SolVd[i][inode] * dt;
        }
      }

      particles[iMarker]->SetMarkerDisplacement (particleDisp);
      particles[iMarker]->UpdateParticleCoordinates();

      std::vector <double> particleAcc (dim);
      std::vector <double> particleVel (dim);
      for (unsigned i = 0; i < dim; i++) {
        particleAcc[i] = 1. / (beta * dt * dt) * particleDisp[i] - 1. / (beta * dt) * particleVelOld[i] - (1. - 2.* beta) / (2. * beta) * particleAccOld[i];
        particleVel[i] = particleVelOld[i] + dt * ( (1. - Gamma) * particleAccOld[i] + Gamma * particleAcc[i]);
      }

      particles[iMarker]->SetMarkerVelocity (particleVel);
      particles[iMarker]->SetMarkerAcceleration (particleAcc);

      //END

      for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
          GradSolDpHat[i][j] = 0.;
          for (unsigned inode = 0; inode < nve; inode++) {
            GradSolDpHat[i][j] +=  gradphi_hat[inode * dim + j] * SolDd[i][inode];
          }
        }
      }

      std::vector < std::vector < double > > FpOld;
      FpOld = particles[iMarker]->GetDeformationGradient();

      double FpNew[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
      std::vector < std::vector < double > > Fp (dim);

      for (unsigned i = 0; i < dim; i++) {
        for (unsigned j = 0; j < dim; j++) {
          FpNew[i][j] += GradSolDpHat[i][j];
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

void ParticlesToGridProjection (MultiLevelSolution & mlSol, Line & linea) {

  const unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;
  Mesh* msh = mlSol._mlMesh->GetLevel (level);
  Solution* mysolution = mlSol.GetSolutionLevel (level);

  Mesh* mymsh = mlSol._mlMesh->GetLevel (level);
  elem* myel = mymsh->el;

  const unsigned dim = mymsh->GetDimension();

  unsigned iproc  = mymsh->processor_id();

  vector < double > phi_hat;
  vector < double > gradphi_hat;
  vector < double > nablaphi_hat;

  vector <vector < double> > vx_hat (dim);

  double weight;

  const char varname[10][6] = {"DX", "DY", "DZ", "VXOld", "VYOld", "VZOld", "AXOld", "AYOld", "AZOld", "gM"};
  vector <unsigned> indexSolD (dim);
  vector <unsigned> indexSolVOld (dim);
  vector <unsigned> indexSolAOld (dim);
  unsigned solType = mlSol.GetSolutionType (&varname[3][0]);

  unsigned indexGridMass = mlSol.GetIndex (&varname[9][0]);
  for (unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolD[ivar] = mlSol.GetIndex (&varname[ivar][0]);
    indexSolVOld[ivar] = mlSol.GetIndex (&varname[ivar + 3][0]);
    indexSolAOld[ivar] = mlSol.GetIndex (&varname[ivar + 6][0]);
  }

  std::vector<unsigned> markerOffset = linea.GetMarkerOffset();
  unsigned markerOffset1 = markerOffset[iproc];
  unsigned markerOffset2 = markerOffset[iproc + 1];
  std::vector<Marker*> particles = linea.GetParticles();

  unsigned ielOld = UINT_MAX;

  //BEGIN set grid mass, old grid velocity and old grid acceleration to zero
  mysolution->_Sol[indexGridMass]->zero();
  for (int i = 0; i < dim; i++) {
    mysolution->_Sol[indexSolVOld[i]]->zero();
    mysolution->_Sol[indexSolAOld[i]]->zero();
  }
  mysolution->_Sol[indexGridMass]->close();
  for (int i = 0; i < dim; i++) {
    mysolution->_Sol[indexSolVOld[i]]->close();
    mysolution->_Sol[indexSolAOld[i]]->close();
  }
  //END


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
//END

//BEGIN loop on particles to build old grid velocity and old grid acceleration
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

      std::vector <double> particleVel (dim);
      particles[iMarker]->GetMarkerVelocity (particleVel);

      std::vector <double> particleAcc (dim);
      particles[iMarker]->GetMarkerAcceleration (particleAcc);

      for (unsigned inode = 0; inode < nve; inode++) {
        unsigned idof = mymsh->GetSolutionDof (inode, iel, solType);
        double inodeGridMass = (*mysolution->_Sol[indexGridMass]) (idof);
        for (int i = 0; i < dim; i++) {
          double velLocal = (fabs (inodeGridMass) > 1.e-7) ? (phi_hat[inode] * particleMass * particleVel[i]) / inodeGridMass : 0.;
          mysolution->_Sol[indexSolVOld[i]]->add (idof, velLocal);
          double accLocal = (fabs (inodeGridMass) > 1.e-7) ? (phi_hat[inode] * particleMass * particleAcc[i]) / inodeGridMass : 0.;
          mysolution->_Sol[indexSolAOld[i]]->add (idof, accLocal);
        }
      }

      ielOld = iel;

    }

    else {
      break;
    }
  }

  for (int i = 0; i < dim; i++) {
    mysolution->_Sol[indexSolVOld[i]]->close();
    mysolution->_Sol[indexSolAOld[i]]->close();
  }
//END

}

void ProjectVelAcc (MultiLevelProblem &ml_prob) {

  TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("MPM_FSI");
  const unsigned  level = my_nnlin_impl_sys.GetLevelToAssemble();
  MultiLevelSolution* ml_sol = ml_prob._ml_sol;
  Solution* sol = ml_sol->GetSolutionLevel (level);

  Mesh* msh = ml_prob._ml_msh->GetLevel (level);
  elem* myel = msh->el;

  double dt =  my_nnlin_impl_sys.GetIntervalTime();
  const unsigned dim = msh->GetDimension();

  unsigned iproc  = msh->processor_id();
  unsigned nprocs = msh->n_processors();

  vector< vector< double > > solV (dim);
  vector< vector< double > > solA (dim);

  vector< bool > nodeFlag;
  vector< unsigned > idof;

  vector < double > phi;
  vector <vector < double> > vx (dim);  //vx is coordX in assembly of ex30
  vector <vector < double> > xp;

  //variable-name handling
  const char varname[15][6] = {"DX", "DY", "DZ", "VX", "VY", "VZ", "VXOld", "VYOld", "VZOld", "AX", "AY", "AZ",  "AXOld", "AYOld", "AZOld"};

  vector <unsigned> indexSolD (dim);
  vector <unsigned> indexSolV (dim);
  vector <unsigned> indexSolVOld (dim);
  vector <unsigned> indexSolA (dim);
  vector <unsigned> indexSolAOld (dim);

  for (unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolD[ivar] = ml_sol->GetIndex (&varname[ivar][0]);
    indexSolV[ivar] = ml_sol->GetIndex (&varname[ivar + 3][0]);
    indexSolVOld[ivar] = ml_sol->GetIndex (&varname[ivar + 6][0]);
    indexSolA[ivar] = ml_sol->GetIndex (&varname[ivar + 9][0]);
    indexSolAOld[ivar] = ml_sol->GetIndex (&varname[ivar + 12][0]);
  }
  unsigned indexNodeFlag =  ml_sol->GetIndex ("NodeFlag");

  unsigned solType = ml_sol->GetSolutionType (&varname[3][0]);

  sol->_Sol[indexNodeFlag]->zero();

  unsigned counter = 0;

  std::vector < std::vector < std::vector <double > > > aP (3);

  //BEGIN evaluate the displacement and acceleration with Newmark
  for (unsigned idof = msh->_dofOffset[solType][iproc]; idof < msh->_dofOffset[solType][iproc + 1]; idof++) {
    for (unsigned k = 0; k < dim; k++) {
      double v = (*sol->_Sol[indexSolV[k]]) (idof);
      double vOld = (*sol->_Sol[indexSolVOld[k]]) (idof);
      double aOld = (*sol->_Sol[indexSolAOld[k]]) (idof);
      double d = vOld * dt + (v - vOld) * beta * dt / Gamma + aOld * dt * dt * (0.5 - beta / Gamma);
      double a = (v - vOld) / (Gamma * dt) + (1. - 1. / Gamma) * aOld;
      sol->_Sol[indexSolD[k]]->set (idof, d);
      sol->_Sol[indexSolA[k]]->set (idof, a);
    }
  }
  for (unsigned k = 0; k < dim; k++) {
    sol->_Sol[indexSolD[k]]->close();
    sol->_Sol[indexSolA[k]]->close();
  }
  //END

  //BEGIN COPY SOLUTION TO OLD SOLUTION
  for (unsigned k = 0; k < dim; k++) {
    (*sol->_Sol[indexSolVOld[k]]) = (*sol->_Sol[indexSolV[k]]);
    sol->_Sol[indexSolV[k]]->zero();
    (*sol->_Sol[indexSolAOld[k]]) = (*sol->_Sol[indexSolA[k]]);
    sol->_Sol[indexSolA[k]]->zero();
  }
  //END

  //BEGIN loop on elements
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielType = msh->GetElementType (iel);

    unsigned nDofs = msh->GetElementDofNumber (iel, solType);   // number of solution element dofs

    for (unsigned  k = 0; k < dim; k++) {
      solV[k].resize (nDofs);
      solA[k].resize (nDofs);
      vx[k].resize (nDofs);
    }
    nodeFlag.resize (nDofs);
    idof.resize (nDofs);
    xp.resize (nDofs);
    for (unsigned i = 0; i < nDofs; i++) {
      xp[i].resize (dim);
    }

    for (unsigned i = 0; i < nDofs; i++) {
      idof[i] = msh->GetSolutionDof (i, iel, solType);
      unsigned idofX = msh->GetSolutionDof (i, iel, 2);

      for (unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[indexSolVOld[k]]) (idof[i]);
        solA[k][i] = (*sol->_Sol[indexSolAOld[k]]) (idof[i]);
        xp[i][k] = (*msh->_topology->_Sol[k]) (idofX);    // coordinates of the reference configuration;
        vx[k][i] = xp[i][k] + (*sol->_Sol[indexSolD[k]]) (idof[i]);    // coordinates of the deformed configuration
      }
      nodeFlag[i] = ( (*sol->_Sol[indexNodeFlag]) (idof[i]) > 0.5) ? true : false;
    }

    bool aPIsInitialized = false;

    double r;
    std::vector <double> xc;
    GetConvexHullSphere (vx, xc, r, 0.0001); // get the ball that circumscribe the element
    double r2 = r * r;

    std::vector < std::vector< double > > xe; // get the box that encloses the element
    GetBoundingBox (vx, xe, 0.0001);

    for (unsigned i = 0; i < nDofs; i++) { // loop on the nodes of the reference elements now considered as independent points
      if (!nodeFlag[i]) {
        double d2 = 0.;
        for (int k = 0; k < dim; k++) {
          d2 += (xp[i][k] - xc[k]) * (xp[i][k] - xc[k]);
        }
        bool insideHull = true;
        if (d2 > r2) {
          insideHull = false;
        }
        for (unsigned k = 0; k < dim; k++) {
          if (xp[i][k] < xe[k][0] || xp[i][k] > xe[k][1]) {
            insideHull = false;
          }
        }

        if (insideHull) { //rough test
          if (!aPIsInitialized) {
            aPIsInitialized = true;
            std::vector < std::vector <double> > x1 (dim);
            for (unsigned jtype = 0; jtype < solType + 1; jtype++) {
              ProjectNodalToPolynomialCoefficients (aP[jtype], vx, ielType, jtype) ;
            }
          }
          std::vector <double> xi;
          GetClosestPointInReferenceElement (vx, xp[i], ielType, xi);

          bool inverseMapping = GetInverseMapping (solType, ielType, aP, xp[i], xi, 100);
          if (!inverseMapping) {
            std::cout << "InverseMapping failed at " << iel << " " << idof[i] << std::endl;
          }


          bool insideDomain = CheckIfPointIsInsideReferenceDomain (xi, ielType, 1.e-3); // fine testing

          if (inverseMapping && insideDomain) {
            sol->_Sol[indexNodeFlag]->add (idof[i], 1.);
            msh->_finiteElement[ielType][solType]->GetPhi (phi, xi);
            //std::cout << iel << " " << i << "  ";
            counter++;
            for (unsigned k = 0; k < dim; k++) {
              double solVk = 0.;
              double solAk = 0.;
              for (unsigned j = 0; j < nDofs; j++)    {
                solVk += phi[j] * solV[k][j];
                solAk += phi[j] * solA[k][j];
              }
              sol->_Sol[indexSolV[k]]->add (idof[i], solVk);
              sol->_Sol[indexSolA[k]]->add (idof[i], solAk);
            }
          }
        }
      }
    }
  }

  sol->_Sol[indexNodeFlag]->close();

  for (unsigned k = 0; k < dim; k++) {
    sol->_Sol[indexSolV[k]]->close();
    sol->_Sol[indexSolA[k]]->close();
  }

  unsigned c0 = 0;
  for (unsigned i = msh->_dofOffset[solType][iproc]; i < msh->_dofOffset[solType][iproc + 1]; i++) {
    unsigned cnt = static_cast < unsigned > (floor ( (*sol->_Sol[indexNodeFlag]) (i) + 0.5));
    if (cnt == 0) {
      c0++;
    }
    else if (cnt > 1) {
      counter -= (cnt - 1);
      for (unsigned k = 0; k < dim; k++) {
        double velk = (*sol->_Sol[indexSolV[k]]) (i) / cnt;
        sol->_Sol[indexSolV[k]]->set (i, velk);
        double acck = (*sol->_Sol[indexSolA[k]]) (i) / cnt;
        sol->_Sol[indexSolA[k]]->set (i, acck);
      }
    }
  }

  unsigned counterAll;
  MPI_Reduce (&counter, &counterAll, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  std::cout << "COUNTER = " << counterAll << " " << msh->GetTotalNumberOfDofs (solType) << std::endl;

  idof.resize (c0);
  vector< double > xp0 (c0 * dim);

  unsigned c1 = 0;
  for (unsigned i = msh->_dofOffset[solType][iproc]; i < msh->_dofOffset[solType][iproc + 1]; i++) {
    if (static_cast < unsigned > (floor ( (*sol->_Sol[indexNodeFlag]) (i) + 0.5)) == 0) {
      idof[c1] = i;
      for (unsigned k = 0; k < dim; k++) {
        xp0[c1 * dim + k] = (*msh->_topology->_Sol[k]) (i);
      }
      c1++;
      if (c1 == c0) break;
    }
  }

  vector< double > xp1;
  for (unsigned jproc = 0; jproc < nprocs; jproc++) {
    c1 = c0;
    MPI_Bcast (&c1, 1, MPI_UNSIGNED, jproc, PETSC_COMM_WORLD);
    if (c1) {
      xp1.resize (c1 * dim);
      if (iproc == jproc) {
        for (unsigned i = 0; i < c1; i++) {
          for (unsigned k = 0; k < dim; k++) {
            xp1[i * dim + k] = xp0[i * dim + k];
          }
        }
      }
      MPI_Bcast (&xp1[0], c1 * dim, MPI_DOUBLE, jproc, PETSC_COMM_WORLD);

      for (unsigned i = 0; i < c1; i++) {
        std::vector < double > xp (dim);
        for (unsigned k = 0; k < dim; k++) {
          xp[k] = xp1[i * dim + k];
        }
        Marker p (xp, 1, VOLUME, sol, solType, true, 1.);
        unsigned mproc = p.GetMarkerProc (sol);
        if (iproc == mproc) {
          unsigned jel = p.GetMarkerElement();
          short unsigned jelType = msh->GetElementType (jel);
          unsigned nDofs = msh->GetElementDofNumber (jel, solType);   // number of solution element dofs

          for (unsigned  k = 0; k < dim; k++) {
            solV[k].resize (nDofs);
            solA[k].resize (nDofs);
            vx[k].resize (nDofs);
          }

          for (unsigned j = 0; j < nDofs; j++) {
            unsigned jdof = msh->GetSolutionDof (j, jel, solType);
            unsigned jdofX = msh->GetSolutionDof (j, jel, 2);
            for (unsigned  k = 0; k < dim; k++) {
              solV[k][j] = (*sol->_Sol[indexSolVOld[k]]) (jdof);    //velocity to be projected
              solA[k][j] = (*sol->_Sol[indexSolAOld[k]]) (jdof);    //acceleration to be projected
              vx[k][j] = (*msh->_topology->_Sol[k]) (jdofX) + (*sol->_Sol[indexSolD[k]]) (jdof);       // coordinates of the deformed configuration
            }
          }
          std::vector <double> xi = p.GetMarkerLocalCoordinates();
          msh->_finiteElement[jelType][solType]->GetPhi (phi, xi);

          std::vector < double > vel (dim, 0.);
          std::vector < double > acc (dim, 0.);
          for (unsigned k = 0; k < dim; k++) {
            for (unsigned j = 0; j < nDofs; j++)    {
              vel[k] += phi[j] * solV[k][j];
              acc[k] += phi[j] * solA[k][j];
            }
          }
          MPI_Send (&vel[0], dim, MPI_DOUBLE, jproc, 1, PETSC_COMM_WORLD);
          MPI_Send (&acc[0], dim, MPI_DOUBLE, jproc, 1, PETSC_COMM_WORLD);

        }
        if (iproc == jproc) {
          std::vector < double > vel (dim);
          std::vector < double > acc (dim);
          MPI_Recv (&vel[0], dim, MPI_DOUBLE, mproc, 1, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
          MPI_Recv (&acc[0], dim, MPI_DOUBLE, mproc, 1, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
          for (unsigned k = 0; k < dim; k++) {
            sol->_Sol[indexSolV[k]]->set (idof[i], vel[k]);
            sol->_Sol[indexSolA[k]]->set (idof[i], acc[k]);
          }
          counter++;
        }
      }
    }
  }

  for (unsigned k = 0; k < dim; k++) {
    sol->_Sol[indexSolV[k]]->close();
    sol->_Sol[indexSolA[k]]->close();
  }

  MPI_Reduce (&counter, &counterAll, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  std::cout << "COUNTER = " << counterAll << " " << msh->GetTotalNumberOfDofs (solType) << std::endl;
  
}

void CopySolutionToSolutionOld (MultiLevelSolution & mlSol) {

  const unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;
  Mesh* msh = mlSol._mlMesh->GetLevel (level);
  Solution* mysolution = mlSol.GetSolutionLevel (level);

  Mesh* mymsh = mlSol._mlMesh->GetLevel (level);

  const unsigned dim = mymsh->GetDimension();
  
    const char varname[12][6] = {"VX", "VY", "VZ", "VXOld", "VYOld", "VZOld", "AX", "AY", "AZ",  "AXOld", "AYOld", "AZOld"};

  vector <unsigned> indexSolV (dim);
  vector <unsigned> indexSolVOld (dim);
  vector <unsigned> indexSolA (dim);
  vector <unsigned> indexSolAOld (dim);

  for (unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolV[ivar] = mlSol.GetIndex (&varname[ivar][0]);
    indexSolVOld[ivar] = mlSol.GetIndex (&varname[ivar + 3][0]);
    indexSolA[ivar] = mlSol.GetIndex (&varname[ivar + 6][0]);
    indexSolAOld[ivar] = mlSol.GetIndex (&varname[ivar + 9][0]);
  }
  
    //BEGIN COPY SOLUTION TO OLD SOLUTION
  //NOTE
  for (unsigned k = 0; k < dim; k++) {
//     (*mysolution->_Sol[indexSolVOld[k]]) = (*mysolution->_Sol[indexSolV[k]]);
//     (*mysolution->_Sol[indexSolAOld[k]]) = (*mysolution->_Sol[indexSolA[k]]);
      mysolution->_Sol[indexSolV[k]]->zero();
      mysolution->_Sol[indexSolA[k]]->zero();
  }
  //END
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


