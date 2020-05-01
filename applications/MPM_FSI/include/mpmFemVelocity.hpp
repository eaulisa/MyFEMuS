#include "MultiLevelSolution.hpp"


using namespace femus;

double beta = 0.25;
double Gamma = 0.5;
double gravity[3] = {0., -9.81, 0.};
double scalingFactor1 = 1.e-6;
double scalingFactor2 = 1.e-10;
double NeumannFactor = .0;
Line* linea;

void AssembleMPMSys(MultiLevelProblem& ml_prob) {

  clock_t AssemblyTime = 0;
  clock_t start_time, end_time;

  //pointers and references

  TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("MPM_FEM");
  const unsigned  level = my_nnlin_impl_sys.GetLevelToAssemble();
  MultiLevelSolution* ml_sol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* mysolution = ml_sol->GetSolutionLevel(level);     // pointer to the solution (level) object
  LinearEquationSolver* myLinEqSolver = my_nnlin_impl_sys._LinSolver[level];  // pointer to the equation (level) object

  Mesh* mymsh = ml_prob._ml_msh->GetLevel(level);     // pointer to the mesh (level) object
  elem* myel = mymsh->el;   // pointer to the elem object in msh (level)
  SparseMatrix* myKK = myLinEqSolver->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* myRES =  myLinEqSolver->_RES;  // pointer to the global residual vector object in pdeSys (level)
  bool assembleMatrix = my_nnlin_impl_sys.GetAssembleMatrix();

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  const unsigned dim = mymsh->GetDimension();

  // reserve memory for the local standar vectors
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  // data
  unsigned iproc  = mymsh->processor_id();

  vector < double > phi;
  vector < double > phi_hat;

  vector < adept::adouble> gradphi;
  vector < double > gradphi_hat;

  phi.reserve(max_size);
  phi_hat.reserve(max_size);

  gradphi.reserve(max_size * dim);
  gradphi_hat.reserve(max_size * dim);

  vector <vector < adept::adouble> > vx(dim); //vx is coordX in assembly of ex30
  vector <vector < double> > vx_hat(dim);


  vector< vector< adept::adouble > > SolVd(dim);      // local solution (displacement)
  vector< vector< adept::adouble > > SolDd(dim);      // local solution (displacement)
  vector< vector< int > > dofsVAR(dim);

  vector< vector< double > > Rhs(dim);     // local redidual vector
  vector< vector< adept::adouble > > aRhs(dim);     // local redidual vector

  vector < int > dofsAll;
  vector < double > Jac;

  adept::adouble weight;
  double weight_hat = 0.;

  //reading parameters for MPM body
  double density_MPM = ml_prob.parameters.get<Solid> ("SolidMPM").get_density();
  double E_MPM = ml_prob.parameters.get<Solid> ("SolidMPM").get_young_module();
  double mu_MPM = ml_prob.parameters.get<Solid> ("SolidMPM").get_lame_shear_modulus();
  double nu_MPM = ml_prob.parameters.get<Solid> ("SolidMPM").get_poisson_coeff();
  double lambda_MPM = ml_prob.parameters.get<Solid> ("SolidMPM").get_lame_lambda();
  double K_MPM = E_MPM / (3.*(1. - 2. * nu_MPM)); //bulk modulus

  double dt =  my_nnlin_impl_sys.GetIntervalTime();

  std::cout.precision(10);
  std::cout << "SF1 = " << scalingFactor1 << " SF2 = " << scalingFactor2 << " NF = " << NeumannFactor << std::endl;

  //variable-name handling
  const char varname[10][5] = {"VX", "VY", "VZ", "Mat"};
  vector <unsigned> indexSolV(dim);
  vector <unsigned> indexPdeV(dim);
  unsigned solType = ml_sol->GetSolutionType(&varname[0][0]);
  for(unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolV[ivar] = ml_sol->GetIndex(&varname[ivar][0]);
    indexPdeV[ivar] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar][0]);
  }
  unsigned indexSolM = ml_sol->GetIndex("M");
  unsigned indexSolMat = ml_sol->GetIndex(&varname[3][0]);
  unsigned solTypeMat = ml_sol->GetSolutionType(&varname[3][0]);

  vector < bool > solidFlag;

  start_time = clock();

  myKK->zero();
  myRES->zero();

  //line instances
  std::vector<unsigned> markerOffset = linea->GetMarkerOffset();
  unsigned markerOffset1 = markerOffset[iproc];
  unsigned markerOffset2 = markerOffset[iproc + 1];
  std::vector<Marker*> particles = linea->GetParticles();

  //BEGIN loop on elements (to initialize the "soft" stiffness matrix)
  for(int iel = mymsh->_elementOffset[iproc]; iel < mymsh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielt = mymsh->GetElementType(iel);

    unsigned nDofs = mymsh->GetElementDofNumber(iel, solType);    // number of solution element dofs
    unsigned nDofsAll = dim * nDofs ;
    // resize local arrays
    std::vector <int> sysDof(nDofsAll);

    solidFlag.resize(nDofs);
    for(unsigned  k = 0; k < dim; k++) {
      SolVd[k].resize(nDofs);
      SolDd[k].resize(nDofs);
      vx[k].resize(nDofs);
      vx_hat[k].resize(nDofs);
    }

    for(unsigned  k = 0; k < dim; k++) {
      aRhs[k].resize(nDofs);    //resize
      std::fill(aRhs[k].begin(), aRhs[k].end(), 0);    //set aRes to zero
    }

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned idof = mymsh->GetSolutionDof(i, iel, solType);    // global to global mapping between solution node and solution dof
      unsigned idofX = mymsh->GetSolutionDof(i, iel, 2);    // global to global mapping between solution node and solution dof

      solidFlag[i] = ((*mysolution->_Sol[indexSolM])(idof) > 0.5 || mymsh->GetSolidMark(idof)) ? true : false;

      for(unsigned  k = 0; k < dim; k++) {
        SolVd[k][i] = (*mysolution->_Sol[indexSolV[k]])(idof);
        sysDof[i + k * nDofs] = myLinEqSolver->GetSystemDof(indexSolV[k], indexPdeV[k], i, iel);    // global to global mapping between solution node and pdeSys dof
        vx_hat[k][i] = (*mymsh->_topology->_Sol[k])(idofX);
      }
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    for(unsigned  k = 0; k < dim; k++) {
      for(unsigned i = 0; i < nDofs; i++) {
        SolDd[k][i] = SolVd[k][i] * dt;
      }
    }

    for(unsigned i = 0; i < nDofs; i++) {
      for(unsigned  k = 0; k < dim; k++) {
        vx[k][i] = vx_hat[k][i] + SolDd[k][i];
      }
    }

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < mymsh->_finiteElement[ielt][solType]->GetGaussPointNumber(); ig++) {

      mymsh->_finiteElement[ielt][solType]->Jacobian(vx_hat, ig, weight_hat, phi_hat, gradphi_hat);

      vector < double > Xg(dim, 0);
      vector < vector < adept::adouble > > GradSolDgssHat(dim);
      vector < vector < adept::adouble > > GradSolDgssTilde(dim);

      for(unsigned  k = 0; k < dim; k++) {
        GradSolDgssHat[k].resize(dim);
        std::fill(GradSolDgssHat[k].begin(), GradSolDgssHat[k].end(), 0);
      }

      for(unsigned i = 0; i < nDofs; i++) {
        for(unsigned j = 0; j < dim; j++) {
          Xg[j] += phi[i] * vx_hat[j][i];
          for(unsigned  k = 0; k < dim; k++) {
            GradSolDgssHat[k][j] += gradphi_hat[i * dim + j] * SolDd[k][i];
          }
        }
      }

      unsigned idofMat = mymsh->GetSolutionDof(0, iel, solTypeMat);
      double  MPMmaterial = (*mysolution->_Sol[indexSolMat])(idofMat);
      double scalingFactor = 0;// / (1. + 100. * distance);
      if(MPMmaterial < 5) scalingFactor = scalingFactor1;
      else if(MPMmaterial < 9) scalingFactor = scalingFactor2;

      double mu = mu_MPM;

      for(unsigned i = 0; i < nDofs; i++) {
        vector < adept::adouble > softStiffness(dim, 0.);

        for(unsigned k = 0; k < dim; k++) {
          for(unsigned  j = 0; j < dim; j++) {
            softStiffness[k]   +=  mu * gradphi_hat[i * dim + j] * (GradSolDgssHat[k][j] + GradSolDgssHat[j][k]);
          }
        }
        for(unsigned  k = 0; k < dim; k++) {
          if(MPMmaterial >= 2) {
            aRhs[k][i] += - softStiffness[k] * weight_hat * scalingFactor;
          }
          else if(!solidFlag[i]) {
            aRhs[k][i] += phi_hat[i] * SolDd[k][i] * weight_hat;
          }
        }
      }

    } // end gauss point loop

    //copy the value of the adept::adoube aRes in double Res and store them in RES
    std::vector<double> Rhs(nDofsAll);  //resize

    for(int i = 0; i < nDofs; i++) {
      for(unsigned  k = 0; k < dim; k++) {
        Rhs[ i +  k * nDofs ] = -aRhs[k][i].value();
      }
    }

    myRES->add_vector_blocked(Rhs, sysDof);

    if(assembleMatrix) {
      Jac.resize(nDofsAll * nDofsAll);
      // define the dependent variables

      for(unsigned  k = 0; k < dim; k++) {
        s.dependent(&aRhs[k][0], nDofs);
      }

      // define the independent variables
      for(unsigned  k = 0; k < dim; k++) {
        s.independent(&SolVd[k][0], nDofs);
      }

      // get the and store jacobian matrix (row-major)
      s.jacobian(&Jac[0] , true);
      myKK->add_matrix_blocked(Jac, sysDof, sysDof);

      s.clear_independents();
      s.clear_dependents();
    }
  }
  //END building "soft" stiffness matrix

  //initialization of iel
  unsigned ielOld = UINT_MAX;

  //BEGIN loop on particles (used as Gauss points)

  for(unsigned iMarker = markerOffset1; iMarker < markerOffset2; iMarker++) {
    //element of particle iMarker
    unsigned iel = particles[iMarker]->GetMarkerElement();
    if(iel != UINT_MAX) {
      short unsigned ielt;
      unsigned nDofs;

      //update element related quantities only if we are in a different element
      if(iel != ielOld) {

        ielt = mymsh->GetElementType(iel);
        nDofs = mymsh->GetElementDofNumber(iel, solType);
        //initialization of everything is in common fluid and solid

        //Rhs
        for(int i = 0; i < dim; i++) {
          dofsVAR[i].resize(nDofs);
          SolVd[i].resize(nDofs);
          SolDd[i].resize(nDofs);
          aRhs[i].resize(nDofs);
          vx[i].resize(nDofs);
          vx_hat[i].resize(nDofs);
        }
        dofsAll.resize(0);

        //BEGIN copy of the value of Sol at the dofs idof of the element iel
        for(unsigned i = 0; i < nDofs; i++) {
          unsigned idof = mymsh->GetSolutionDof(i, iel, solType); //local 2 global solution
          unsigned idofX = mymsh->GetSolutionDof(i, iel, 2); //local 2 global coordinates
          for(int j = 0; j < dim; j++) {
            SolVd[j][i] = (*mysolution->_Sol[indexSolV[j]])(idof);
            dofsVAR[j][i] = myLinEqSolver->GetSystemDof(indexSolV[j], indexPdeV[j], i, iel); //local 2 global Pde
            aRhs[j][i] = 0.;
            vx_hat[j][i] = (*mymsh->_topology->_Sol[j])(idofX);
          }
        }

        // build dof composition
        for(int idim = 0; idim < dim; idim++) {
          dofsAll.insert(dofsAll.end(), dofsVAR[idim].begin(), dofsVAR[idim].end());
        }
        // start a new recording of all the operations involving adept::adouble variables
        s.new_recording();

        for(unsigned  k = 0; k < dim; k++) {
          for(unsigned i = 0; i < nDofs; i++) {
            SolDd[k][i] = SolVd[k][i] * dt;
          }
        }
      }

      // the local coordinates of the particles are the Gauss points in this context
      std::vector <double> xi = particles[iMarker]->GetMarkerLocalCoordinates();
      mymsh->_finiteElement[ielt][solType]->Jacobian(vx_hat, xi, weight_hat, phi_hat, gradphi_hat);

      std::vector <double> SolVpOld(dim);
      particles[iMarker]->GetMarkerVelocity(SolVpOld);

      std::vector <double> SolApOld(dim);
      particles[iMarker]->GetMarkerAcceleration(SolApOld);

      double mass = particles[iMarker]->GetMarkerMass();

      for(int j = 0; j < dim; j++) {
        for(unsigned inode = 0; inode < nDofs; inode++) {
          vx[j][inode] = vx_hat[j][inode] + SolDd[j][inode];
        }
      }

      mymsh->_finiteElement[ielt][solType]->Jacobian(vx, xi, weight, phi, gradphi); //function to evaluate at the particles

      // displacement and velocity
      //BEGIN evaluates SolDp at the particle iMarker
      vector<adept::adouble> SolDp(dim, 0.);
      vector<vector < adept::adouble > > GradSolDpHat(dim);
      for(int i = 0; i < dim; i++) {
        GradSolDpHat[i].assign(dim, 0.);
      }

      for(int i = 0; i < dim; i++) {
        for(unsigned inode = 0; inode < nDofs; inode++) {
          SolDp[i] += phi[inode] * SolDd[i][inode];
          for(int j = 0; j < dim; j++) {
            GradSolDpHat[i][j] +=  gradphi_hat[inode * dim + j] * SolDd[i][inode];
          }
        }
      }
      //END evaluates SolDp at the particle iMarker


      //BEGIN computation of the Cauchy Stress
      std::vector < std::vector < double > > FpOld;
      FpOld = particles[iMarker]->GetDeformationGradient(); //extraction of the deformation gradient

      adept::adouble FpNew[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
      adept::adouble F[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      adept::adouble B[3][3];
      adept::adouble Id2th[3][3] = {{ 1., 0., 0.}, { 0., 1., 0.}, { 0., 0., 1.}};
      adept::adouble Cauchy[3][3];

      for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
          FpNew[i][j] += GradSolDpHat[i][j];
        }
      }

      for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
          for(int k = 0; k < dim; k++) {
            F[i][j] += FpNew[i][k] * FpOld[k][j];
          }
        }
      }

      if(dim == 2) F[2][2] = 1.;

      adept::adouble J_hat =  F[0][0] * F[1][1] * F[2][2] + F[0][1] * F[1][2] * F[2][0] + F[0][2] * F[1][0] * F[2][1]
                              - F[2][0] * F[1][1] * F[0][2] - F[2][1] * F[1][2] * F[0][0] - F[2][2] * F[1][0] * F[0][1];

      for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
          B[i][j] = 0.;

          for(int k = 0; k < 3; k++) {
            //left Cauchy-Green deformation tensor or Finger tensor (B = F*F^T)
            B[i][j] += F[i][k] * F[j][k];
          }
        }
      }

      adept::adouble I1_B = B[0][0] + B[1][1] + B[2][2];

      double mu = mu_MPM;
      double lambda = lambda_MPM;

      for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
          Cauchy[i][j] = lambda * log(J_hat) / J_hat * Id2th[i][j] + mu / J_hat * (B[i][j] - Id2th[i][j]); //alternative formulation
        }
      }
      //END computation of the Cauchy Stress

      //BEGIN redidual Solid Momentum in moving domain
      for(unsigned i = 0; i < nDofs; i++) {
        adept::adouble CauchyDIR[3] = {0., 0., 0.};

        for(int idim = 0.; idim < dim; idim++) {
          for(int jdim = 0.; jdim < dim; jdim++) {
            CauchyDIR[idim] += gradphi[i * dim + jdim] * Cauchy[idim][jdim];
          }
        }

        for(int idim = 0; idim < dim; idim++) {
          aRhs[idim][i] += (phi[i] * gravity[idim] - J_hat * CauchyDIR[idim] / density_MPM
                            - phi[i] * (1. / (beta * dt * dt) * SolDp[idim] - 1. / (beta * dt) * SolVpOld[idim] - (1. - 2.* beta) / (2. * beta) * SolApOld[idim])
                           ) * mass;
        }
      }
      //END redidual Solid Momentum in moving domain


      if(iMarker == markerOffset2 - 1 || iel != particles[iMarker + 1]->GetMarkerElement()) {

        //copy adouble aRhs into double Rhs
        for(unsigned i = 0; i < dim; i++) {
          Rhs[i].resize(nDofs);

          for(int j = 0; j < nDofs; j++) {
            Rhs[i][j] = -aRhs[i][j].value();
          }
        }

        for(int i = 0; i < dim; i++) {
          myRES->add_vector_blocked(Rhs[i], dofsVAR[i]);
        }

        if(assembleMatrix) {
          //Store equations
          for(int i = 0; i < dim; i++) {
            s.dependent(&aRhs[i][0], nDofs);
            s.independent(&SolVd[i][0], nDofs);
          }

          Jac.resize((dim * nDofs) * (dim * nDofs));

          s.jacobian(&Jac[0], true);

          myKK->add_matrix_blocked(Jac, dofsAll, dofsAll);
          s.clear_independents();
          s.clear_dependents();
        }
      }
      //END local to global assembly

      ielOld = iel;
    }
    else {
      break;
    }
  }
  //END loop on particles

  myRES->close();
  mysolution->_Sol[indexSolMat]->close();

  if(assembleMatrix) {
    myKK->close();
  }

  // *************************************
  end_time = clock();
  AssemblyTime += (end_time - start_time);
  // ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

}

void GridToParticlesProjection(MultiLevelProblem & ml_prob, Line & linea) {

  // ml_prob is the global object from/to where get/set all the data
  // level is the level of the PDE system to be assembled
  // levelMax is the Maximum level of the MultiLevelProblem
  // assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  clock_t AssemblyTime = 0;
  clock_t start_time, end_time;

  //pointers and references

  TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("MPM_FEM");
  //NonLinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<NonLinearImplicitSystem> ("MPM_FEM");
  const unsigned  level = my_nnlin_impl_sys.GetLevelToAssemble();
  MultiLevelSolution* ml_sol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* mysolution = ml_sol->GetSolutionLevel(level);     // pointer to the solution (level) object

  Mesh* mymsh = ml_prob._ml_msh->GetLevel(level);     // pointer to the mesh (level) object
  elem* myel = mymsh->el;   // pointer to the elem object in msh (level)

  double dt =  my_nnlin_impl_sys.GetIntervalTime();
  const unsigned dim = mymsh->GetDimension();

  // data
  unsigned iproc  = mymsh->processor_id();

  // local objects
  vector< vector < double > > SolDd(dim);
  vector< vector < double > > SolVd(dim);
  vector< vector < double > > GradSolDpHat(dim);

  for(int i = 0; i < dim; i++) {
    GradSolDpHat[i].resize(dim);
  }

  vector < double > phi_hat;
  vector < double > gradphi_hat;
  vector < double > nablaphi_hat;

  vector <vector < double> > vx_hat(dim); //vx is coordX in assembly of ex30

  double weight;

  //variable-name handling
  const char varname[3][3] = {"VX", "VY", "VZ"};
  vector <unsigned> indexSolV(dim);
  unsigned solType = ml_sol->GetSolutionType(&varname[0][0]);

  for(unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolV[ivar] = ml_sol->GetIndex(&varname[ivar][0]);
  }

  //line instances
  std::vector<unsigned> markerOffset = linea.GetMarkerOffset();
  unsigned markerOffset1 = markerOffset[iproc];
  unsigned markerOffset2 = markerOffset[iproc + 1];
  std::vector<Marker*> particles = linea.GetParticles();

  //initialization of iel
  unsigned ielOld = UINT_MAX;

  //BEGIN loop on particles
  for(unsigned iMarker = markerOffset1; iMarker < markerOffset2; iMarker++) {

    //element of particle iMarker
    unsigned iel = particles[iMarker]->GetMarkerElement();


    if(iel != UINT_MAX) {

      short unsigned ielt;
      unsigned nve;

      //update element related quantities only if we are in a different element
      if(iel != ielOld) {

        ielt = mymsh->GetElementType(iel);
        nve = mymsh->GetElementDofNumber(iel, solType);

        for(int i = 0; i < dim; i++) {
          SolVd[i].resize(nve);
          SolDd[i].resize(nve);
          vx_hat[i].resize(nve);
        }

        //BEGIN copy of the value of Sol at the dofs idof of the element iel
        for(unsigned inode = 0; inode < nve; inode++) {
          unsigned idof = mymsh->GetSolutionDof(inode, iel, solType); //local 2 global solution
          unsigned idofX = mymsh->GetSolutionDof(inode, iel, 2); //local 2 global solution

          for(int i = 0; i < dim; i++) {
            SolVd[i][inode] = (*mysolution->_Sol[indexSolV[i]])(idof);
            SolDd[i][inode] = SolVd[i][inode] * dt;
            //moving domain
            vx_hat[i][inode] = (*mymsh->_topology->_Sol[i])(idofX);
          }
        }
        //END
      }
      std::vector <double> xi = particles[iMarker]->GetMarkerLocalCoordinates();

      mymsh->_finiteElement[ielt][solType]->Jacobian(vx_hat, xi, weight, phi_hat, gradphi_hat, nablaphi_hat); //function to evaluate at the particles

      std::vector <double> particleVelOld(dim);
      particles[iMarker]->GetMarkerVelocity(particleVelOld);

      std::vector <double> particleAccOld(dim);
      particles[iMarker]->GetMarkerAcceleration(particleAccOld);

      std::vector <double> particleDisp(dim, 0.);
      //update displacement and acceleration
      for(int i = 0; i < dim; i++) {
        for(unsigned inode = 0; inode < nve; inode++) {
          particleDisp[i] += phi_hat[inode] * SolDd[i][inode];
        }
      }

      particles[iMarker]->SetMarkerDisplacement(particleDisp);
      particles[iMarker]->UpdateParticleCoordinates();

      std::vector <double> particleAcc(dim);
      std::vector <double> particleVel(dim);
      for(unsigned i = 0; i < dim; i++) {
        particleAcc[i] = 1. / (beta * dt * dt) * particleDisp[i] - 1. / (beta * dt) * particleVelOld[i] - (1. - 2.* beta) / (2. * beta) * particleAccOld[i];
        particleVel[i] = particleVelOld[i] + dt * ((1. - Gamma) * particleAccOld[i] + Gamma * particleAcc[i]);
      }

      particles[iMarker]->SetMarkerVelocity(particleVel);
      particles[iMarker]->SetMarkerAcceleration(particleAcc);

      //   update the deformation gradient

      for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
          GradSolDpHat[i][j] = 0.;
          for(unsigned inode = 0; inode < nve; inode++) {
            GradSolDpHat[i][j] +=  gradphi_hat[inode * dim + j] * SolDd[i][inode];
          }
        }
      }

      std::vector < std::vector < double > > FpOld;
      FpOld = particles[iMarker]->GetDeformationGradient(); //extraction of the deformation gradient

      double FpNew[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
      std::vector < std::vector < double > > Fp(dim);

      for(unsigned i = 0; i < dim; i++) {
        for(unsigned j = 0; j < dim; j++) {
          FpNew[i][j] += GradSolDpHat[i][j];
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
  //END loop on particles


  //BEGIN loop on elements to update grid velocity and acceleration
  for(unsigned idof = mymsh->_dofOffset[solType][iproc]; idof < mymsh->_dofOffset[solType][iproc + 1]; idof++) {
    for(int i = 0; i < dim; i++) {
      mysolution->_Sol[indexSolV[i]]->set(idof, 0.);
    }
  }

  for(int i = 0; i < dim; i++) {
    mysolution->_Sol[indexSolV[i]]->close();
  }
  //END loop on elements to update grid velocity and acceleration

  linea.UpdateLineMPM();

  linea.GetParticlesToGridMaterial();
}


unsigned getNumberOfLayers(const double & a, const double & fac) {
  double da = 1. / fac;
  double b =  da;
  unsigned n = 1;

  while(b < a) {
    da /= fac;
    b += da;
    n++;
    if(n >= 100) {
      std::cout << "Error: number of layer is unbounded, try with a smaller factor\n";
      abort();
    }
  }
  return n;
}




