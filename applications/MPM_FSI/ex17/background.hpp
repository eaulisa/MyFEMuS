
bool SetBoundaryConditionB(const std::vector < double >&x, const char name[], double &value, const int facename, const double t);
double SetVariableTimeStepB(const double time);

void InitializeBackgroundVariables(MultiLevelSolution &mlSol) {

  unsigned dim = mlSol._mlMesh->GetDimension();

  FEOrder femOrder = FIRST;

  //add variables to mlSol
  mlSol.AddSolution("DX", LAGRANGE, femOrder, 2);
  mlSol.AddSolution("DY", LAGRANGE, femOrder, 2);
  if(dim == 3) mlSol.AddSolution("DZ", LAGRANGE, femOrder, 2);
  mlSol.AddSolution("VX", LAGRANGE, femOrder, 2);
  mlSol.AddSolution("VY", LAGRANGE, femOrder, 2);
  if(dim == 3) mlSol.AddSolution("VZ", LAGRANGE, femOrder, 2);
  
   mlSol.AddSolution("P", LAGRANGE, FIRST, 2);

  mlSol.AddSolution("eflag", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
  mlSol.AddSolution("nflag", LAGRANGE, femOrder, 0, false);
  
  mlSol.Initialize("All");
  //mlSol.Initialize("DX", InitVariableDX);
  //mlSol.Initialize("DY", InitVariableDY);

  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryConditionB);

}


void AssembleMPMSys(MultiLevelProblem& ml_prob) {

  clock_t AssemblyTime = 0;
  clock_t start_time, end_time;

  //pointers and references

  TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("MPM_FSI");
  const unsigned  level = my_nnlin_impl_sys.GetLevelToAssemble();
  MultiLevelSolution* mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* mysolution = mlSol->GetSolutionLevel(level);     // pointer to the solution (level) object
  
  LinearEquationSolver* myLinEqSolver = my_nnlin_impl_sys._LinSolver[level];  // pointer to the equation (level) object

  Mesh* msh = ml_prob._ml_msh->GetLevel(level);     // pointer to the mesh (level) object
  elem* el = msh->el;   // pointer to the elem object in msh (level)
  SparseMatrix* myKK = myLinEqSolver->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* myRES =  myLinEqSolver->_RES;  // pointer to the global residual vector object in pdeSys (level)

  myKK->zero();
  myRES->zero();

  //AssembleGhostPenaltyP(ml_prob, true);
  //AssembleGhostPenaltyP(ml_prob, false);


  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  const unsigned dim = msh->GetDimension();

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  // data
  unsigned iproc  = msh->processor_id();

  vector< vector< adept::adouble > > solDb(dim);   // local background solution (displacement)
  vector< adept::adouble > solP;
  
  vector< vector< double > > solVOld(dim);
  vector< vector< adept::adouble > > solVNew(dim);      // local solution (velocity)
  
  
  vector< double > rhs;    // local redidual vector
  vector< vector< adept::adouble > > aResD(dim);     // local redidual vector
  vector< vector< adept::adouble > > aResV(dim);     // local redidual vector
  vector< adept::adouble > aResP;    // local redidual vector
  vector < double > Jac;

  std::vector <unsigned> sysDofsAll;
  
  vector < double > phi;
  
  vector < double > phiP;
  vector < adept::adouble> gradPhiP;  // phiP_x
  
  vector < double > gradPhiOld;          
  vector < adept::adouble> gradPhi;     // phi_x
  vector < adept::adouble> gradPhiNew;  // phi_x
  
   vector < adept::adouble> nablaphi; //phi_xx
  
  unsigned dim2 = 3 * (dim - 1);
  
  vector <vector < double> > vxOld(dim);
  vector <vector < adept::adouble> > vx(dim);
  vector <vector < adept::adouble> > vxNew(dim);

  double weightOld;
  adept::adouble weight;
  adept::adouble weightNew;


  //reading parameters for MPM body
  double rhoMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_density();
  double EMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_young_module();
  double muMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_lame_shear_modulus();
  double nuMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_poisson_coeff();
  double lambdaMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_lame_lambda();
  
  //reading parameters for fluid FEM domain
  double rhoFluid = ml_prob.parameters.get<Fluid> ("FluidFEM").get_density();
  double muFluid = ml_prob.parameters.get<Fluid> ("FluidFEM").get_viscosity();

  std::cout << muMpm << " " << lambdaMpm << " " << nuMpm << std::endl;
  std::cout << rhoMpm << " " << EMpm << " " << std::endl;
  std::cout << rhoFluid << " " << muFluid << " " << std::endl;

  double dt =  my_nnlin_impl_sys.GetIntervalTime();

  std::cout.precision(10);

  //variable-name handling
  const char varname[10][5] = {"DX", "DY", "DZ", "VX", "VY", "VZ"};

  vector <unsigned> indexSolD(dim);
  vector <unsigned> indexSolV(dim);
  vector <unsigned> indexPdeD(dim);
  vector <unsigned> indexPdeV(dim);
  for(unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolD[ivar] = mlSol->GetIndex(&varname[ivar][0]);
    indexSolV[ivar] = mlSol->GetIndex(&varname[ivar + 3][0]);
    indexPdeD[ivar] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar][0]);       // DX, DY, DZ
    indexPdeV[ivar] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar + 3][0]);   //VX, VY, VZ
  }
  unsigned solType = mlSol->GetSolutionType(&varname[0][0]);

  unsigned indexSolP = mlSol->GetIndex("P");
  unsigned indexPdeP = my_nnlin_impl_sys.GetSolPdeIndex("P");
  unsigned solTypeP = mlSol->GetSolutionType("P");

  unsigned eflagIndex = mlSol->GetIndex("eflag");
  unsigned nflagIndex = mlSol->GetIndex("nflag");
  std::vector < unsigned >  nodeFlag; // local solution

  start_time = clock();


  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielt = msh->GetElementType(iel);

    unsigned nDofs = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solTypeP);  // number of pressure dofs

    unsigned nDofsAll = 2 * dim * nDofs + nDofsP;


    // resize local arrays
    sysDofsAll.resize(nDofsAll);

    nodeFlag.resize(nDofs);
    double tempEflag = (*mysolution->_Sol[eflagIndex])(iel);
    unsigned eFlag = static_cast <unsigned>(floor(tempEflag + 0.25));
    unsigned eFlag1 = 2; // interface or solid
    if(eFlag == 0) {
      eFlag1 = (tempEflag < 0.25) ?  0 : 1; //fluid or fluid-particle
    }


    for(unsigned  k = 0; k < dim; k++) {
      solDb[k].resize(nDofs);
      
      solVOld[k].resize(nDofs);
      solVNew[k].resize(nDofs);

      aResD[k].assign(nDofs, 0.);
      aResV[k].assign(nDofs, 0.);

      vxOld[k].resize(nDofs);
      vx[k].resize(nDofs);
      vxNew[k].resize(nDofs);
    }
    solP.resize(nDofsP);
    aResP.assign(nDofsP, 0.);


    for(unsigned i = 0; i < nDofs; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, solType);

      nodeFlag[i] = (*mysolution->_Sol[nflagIndex])(idof);

      for(unsigned  k = 0; k < dim; k++) {
        solDb[k][i] = (*mysolution->_Sol[indexSolD[k]])(idof); //t_{n+1} -t_n

        solVNew[k][i] = (*mysolution->_Sol[indexSolV[k]])(idof);//t_{n+1}
        solVOld[k][i] = (*mysolution->_SolOld[indexSolV[k]])(idof);//t_n

        sysDofsAll[i + k * nDofs] = myLinEqSolver->GetSystemDof(indexSolD[k], indexPdeD[k], i, iel);
        sysDofsAll[i + (k + dim) * nDofs] = myLinEqSolver->GetSystemDof(indexSolV[k], indexPdeV[k], i, iel);
      }
    }

    for(unsigned i = 0; i < nDofsP; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, solTypeP);
      solP[i] = (*mysolution->_Sol[indexSolP])(idof);
      sysDofsAll[i + (2 * dim) * nDofs] = myLinEqSolver->GetSystemDof(indexSolP, indexPdeP, i, iel);
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    for(unsigned i = 0; i < nDofs; i++) {
      unsigned idofX = msh->GetSolutionDof(i, iel, 2);
      for(unsigned  k = 0; k < dim; k++) {
        vxOld[k][i] = (*msh->_topology->_Sol[k])(idofX); // undeformed background configuration 
        vx[k][i]  = vxOld[k][i] + (1. - af) * solDb[k][i]; // deformed background configuration at alpha_f/theta
        vxNew[k][i]  = vxOld[k][i] + solDb[k][i]; // deformed background configuration at t_{n+1}
      }
    }

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielt][solType]->GetGaussPointNumber(); ig++) {

      msh->_finiteElement[ielt][solTypeP]->Jacobian(vx, ig, weight, phiP, gradPhiP);

      
      msh->_finiteElement[ielt][solType]->Jacobian(vxOld, ig, weightOld, phi, gradPhiOld);
      msh->_finiteElement[ielt][solType]->Jacobian(vx,    ig, weight,    phi, gradPhi, nablaphi);
      msh->_finiteElement[ielt][solType]->Jacobian(vxNew, ig, weightNew, phi, gradPhiNew);
      
      vector < adept::adouble > solVgNew(dim, 0.);
      vector < adept::adouble > solVg(dim, 0.);
      vector < adept::adouble > solVgOld(dim, 0.);

      vector < adept::adouble > solDg(dim, 0.);

      vector < vector < adept::adouble > > gradSolVg(dim);
      vector<vector<adept::adouble> > DeltaSolVg(dim);

      vector < vector < adept::adouble > > gradSolDgOld(dim);
      vector < vector < adept::adouble > > gradSolVgNew(dim);

      for(unsigned  k = 0; k < dim; k++) {
        gradSolDgOld[k].assign(dim, 0);
        gradSolVgNew[k].assign(dim, 0);
        gradSolVg[k].assign(dim, 0.);
        DeltaSolVg[k].resize(dim2, 0.);
      }

      for(unsigned i = 0; i < nDofs; i++) {
        for(unsigned j = 0; j < dim; j++) {
          solVgNew[j] += phi[i] * solVNew[j][i]; // new velocity of background grid
          solVgOld[j] += phi[i] * solVOld[j][i]; // velocity in the undeformed reference configuration
          solDg[j] += phi[i] * solDb[j][i]; // new displacement
          for(unsigned  k = 0; k < dim; k++) {
            gradSolDgOld[k][j] += gradPhiOld[i * dim + j] * solDb[k][i]; //gradient of new solution with respect to deformed reference configuration
            gradSolVgNew[k][j] += gradPhiNew[i * dim + j] * solVNew[k][i]; // gradient of the new velocity with respect to the theta domain
            gradSolVg[k][j] += gradPhi[i * dim + j] * (theta * solVNew[k][i] + (1. - theta) * solVOld[k][i]); // gradient of the theta velocity with respect to the theta domain
          }
        }
        for(unsigned j = 0; j < dim; j++) {
          solVg[j] += theta * solVgNew[j] + (1. - theta) * solVgOld[j]; //velocity in the deformed theta configuration  
        }

        for(unsigned j = 0; j < dim2; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            DeltaSolVg[k][j]   += nablaphi[i * dim2 + j] * (theta * solVNew[k][i] + (1. - theta) * solVOld[k][i]) ; // laplace of the theta velocity with respect to the theta domain
          }
        }
      }

      adept::adouble solPg = 0.;
      for(unsigned i = 0; i < nDofsP; i++) {
        solPg += phiP[i] * solP[i];
      }

      if(eFlag == 0) { //bulk fluid: all this equation is centered at a_f, There are no time derivatives

        adept::adouble F[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned k = 0; k < dim; k++) {
            F[j][k] += gradSolDgOld[j][k]; //tilde F in the ALE equation
          }
        }

        adept::adouble J_Old =  F[0][0] * F[1][1] * F[2][2] + F[0][1] * F[1][2] * F[2][0] + F[0][2] * F[1][0] * F[2][1]
                                - F[2][0] * F[1][1] * F[0][2] - F[2][1] * F[1][2] * F[0][0] - F[2][2] * F[1][0] * F[0][1];

        adept::adouble B[3][3];
        for(unsigned i = 0; i < 3; i++) {
          for(int j = 0; j < 3; j++) {
            B[i][j] = 0.;
            for(unsigned k = 0; k < 3; k++) {
              //left Cauchy-Green deformation tensor or Finger tensor (B = F*F^T)
              B[i][j] += F[i][k] * F[j][k];
            }
          }
        }

        adept::adouble I1_B = B[0][0] + B[1][1] + B[2][2];
        adept::adouble Id2th[3][3] = {{ 1., 0., 0.}, { 0., 1., 0.}, { 0., 0., 1.}};
        adept::adouble sigma[3][3];

        double E = pow(10., eFlag1);
        double nu = 0.4;

        double mu = E / (2. * (1. + nu));
        double lambda = (E * nu) / ((1. + nu) * (1. - 2.*nu));

        for(unsigned j = 0; j < 3; j++) {
          for(unsigned k = 0; k < 3; k++) {
            sigma[j][k] = lambda * log(J_Old) / J_Old * Id2th[j][k] + mu / J_Old * (B[j][k] - Id2th[j][k]);    // alternative formulation
          }
        }
        //END computation of the Cauchy Stress
        for(unsigned i = 0; i < nDofs; i++) {//Auxiliary Equations
          if(nodeFlag[i] == 0) {
            for(unsigned k = 0.; k < dim; k++) {
              adept::adouble cauchy = 0.;
              for(unsigned j = 0.; j < dim; j++) {
                //cauchy += sigma[k][j] * gradPhi[i * dim + j] ;
                cauchy += gradSolDgOld[k][j] * gradPhi[i * dim + j] ;
              }
              aResD[k][i] += cauchy * weight;
            }
          }
        }
      }


//               adept::adouble wlaplaceD  = 0.;
//               for(unsigned  j = 0; j < dim; j++) {
//                 wlaplaceD +=  gradPhiHat[i * dim + j] * (gradSolDgHat[k][j] + gradSolDgHat[j][k]);
//               }
//               aResD[k][i] +=  wlaplaceD * weightHat;



      if(eFlag == 0) {   // only fluid cells

        //start SUPG paramters, tauM, tauC, G to get tauM_SupgPhi
        std::vector <std::vector <adept::adouble> > JacMatrix;
        std::vector <std::vector <double> > JacMatrixHat;
        msh->_finiteElement[ielt][solType]->GetJacobian(vx, ig, weight, JacMatrix); //centered at theta


        if(solTypeP == 4) { //discontinuous pressure <1,\xi,\eta> bases centered at theta
          for(unsigned j = 0; j < dim; j++) {
            gradPhiP[0 * dim + j]  = 0.;
          }
          for(unsigned i = 0; i < dim; i++) {
            for(unsigned j = 0; j < dim; j++) {
              gradPhiP[(i + 1) * dim + j] = JacMatrix[i][j];
            }
          }
        }

        vector<adept::adouble> gradSolPg(dim, 0.); //centered at theta
        for(unsigned i = 0; i < nDofsP; i++) {
          for(unsigned k = 0; k < dim; k++) {
            gradSolPg[k] += solP[i] * gradPhiP[i * dim + k];
          }
        }


        std::vector <std::vector <adept::adouble> > G(dim); // J^T . J //centered at theta
        for(unsigned i = 0; i < dim; i++) {
          G[i].assign(dim, 0.);
          for(unsigned j = 0; j < dim; j++) {
            for(unsigned k = 0; k < dim; k++) {
              G[i][j] += JacMatrix[k][i] * JacMatrix[k][j];
            }
          }
        }

        adept::adouble tauM = 0.;
        double CI = 36.;
        adept::adouble denom = pow(2 * rhoFluid / dt, 2.);
        for(unsigned i = 0; i < dim; i++) {
          for(unsigned j = 0; j < dim; j++) {
            denom += rhoFluid * (solVg[i] - (solDg[i] - 0.) / dt) * G[i][j] * rhoFluid * (solVg[j] - (solDg[j] - 0.) / dt) // we can get the mesh velocity at af
                     + CI * muFluid * muFluid * G[i][j] * G[i][j]; //this could be improved if we had the old mesh velocity
          }
        }
        tauM += 1. / sqrt(denom);

        adept::adouble tauMtrG = 0.;
        for(unsigned k = 0; k < dim; k++) {
          tauMtrG += G[k][k];
        }
        tauMtrG *= tauM;
        adept::adouble tauC = 1. / tauMtrG;

        //end SUPG parameters
        std::vector < adept::adouble > tauM_SupgPhi(nDofs, 0.);
        for(unsigned i = 0; i < nDofs; i++) {
          for(unsigned j = 0; j < dim; j++) {
            tauM_SupgPhi[i] += tauM * (rhoFluid * (solVg[j] - (solDg[j] - 0.) / dt) * gradPhi[i * dim + j]);
          }
        }

        adept::adouble divVg = 0.;
        for(unsigned k = 0; k < dim; k++) {
          divVg += gradSolVg[k][k];
        }

        for(unsigned i = 0; i < nDofs; i++) {
          for(unsigned k = 0; k < dim; k++) {
            adept::adouble wlaplace = 0.;
            adept::adouble SupgLaplace = 0.;
            adept::adouble advection = 0.;

            for(unsigned j = 0; j < dim; j++) {
              wlaplace += muFluid * gradPhi[i * dim + j] * (gradSolVg[k][j] + gradSolVg[j][k]);

              unsigned kdim;
              if(k == j) kdim = j;
              else if(1 == k + j) kdim = dim;       // xy
              else if(2 == k + j) kdim = dim + 2;    // xz
              else if(3 == k + j) kdim = dim + 1;    // yz
              SupgLaplace += (-muFluid * (DeltaSolVg[k][j] + DeltaSolVg[j][kdim])) * tauM_SupgPhi[i];  // SUPG laplace

              advection  +=  rhoFluid * (solVg[j] - (solDg[j] - 0.) / dt) * gradSolVg[k][j] * (phi[i] + tauM_SupgPhi[i]);  // ALE advection + SUPG Advection

            }

            adept::adouble SupgDiv = tauC * divVg * gradPhi[i * dim + k];
            adept::adouble SupgPressure = gradSolPg[k] * tauM_SupgPhi[i];

            aResV[k][i] += rhoFluid * (solVgNew[k] * weightNew - solVgOld[k] * weightOld) / dt * phi[i]
                           + (rhoFluid * (solVgNew[k] - solVgOld[k]) / dt * tauM_SupgPhi[i]
                              + advection
                              + wlaplace + SupgLaplace
                              - weakP * gradPhi[i * dim + k] * solPg
                              + !weakP * phi[i] * gradSolPg[k]
                              + SupgPressure
                              + SupgDiv
                             ) * weight;
          }
        }


        //continuity block
        for(unsigned i = 0; i < nDofsP; i++) {
          for(unsigned k = 0; k < dim; k++) {

            adept::adouble sLaplace = 0.;
            adept::adouble advection = 0.;

            for(unsigned j = 0; j < dim; j++) {
              unsigned kdim;

              if(k == j) kdim = j;
              else if(1 == k + j) kdim = dim;       // xy
              else if(2 == k + j) kdim = dim + 2;   // xz
              else if(3 == k + j) kdim = dim + 1;   // yz

              sLaplace += (- muFluid * (DeltaSolVg[k][j] + DeltaSolVg[j][kdim]));
              advection += rhoFluid * (solVg[j] - (solDg[j] - 0.) / dt) * gradSolVg[k][j];

            }

            aResP[i] +=  phiP[i] * gradSolVgNew[k][k] * weightNew
                         + ((rhoFluid * (solVgNew[k] - solVgOld[k]) / dt + advection +
                             sLaplace +  gradSolPg[k]) * tauM * gradPhiP[i * dim + k]
                           ) * weight;

          }
        }
      }
      else if(eFlag == 2) {   // only solid cells: fake pressure //TODO
        for(unsigned i = 0; i < nDofsP; i++) {
          aResP[i] += 1.0e-10 * phiP[i] * solP[i] * weight;
        }
      }
    } // end gauss point loop




//     //BEGIN BULK PARTICLE
//     if(eFlag > 0) {   // interface and solid
//       while(iBmarker < markerOffsetBulk[iproc + 1] && iel > particlesBulk[iBmarker]->GetMarkerElement()) {
//         iBmarker++;
//       }
//       while(iBmarker < markerOffsetBulk[iproc + 1] && iel == particlesBulk[iBmarker]->GetMarkerElement()) {
// 
//         // the local coordinates of the particles are the Gauss points in this context
//         std::vector <double> xi = particlesBulk[iBmarker]->GetMarkerLocalCoordinates();
//         double dg1 = particlesBulk[iBmarker]->GetMarkerDistance();
// 
//         double U = GetSmoothStepFunction(dg1, eps);
// 
//         double area = particlesBulk[iBmarker]->GetMarkerMass();
// 
//         std::vector <std::vector <adept::adouble> > JacMatrix1;
//         std::vector <std::vector <adept::adouble> > JacMatrix;
//         msh->_finiteElement[ielt][solType]->GetJacobianMatrix(vx, xi, JacMatrix1, JacMatrix); //centered at theta
// 
//         msh->_finiteElement[ielt][solTypeP]->Jacobian(vx, xi, weight, phiP, gradPhiP);
// 
//         msh->_finiteElement[ielt][solType]->Jacobian(vxNew, xi, weightNew, phi, gradPhiNew);
//         msh->_finiteElement[ielt][solType]->Jacobian(vx, xi, weight, phi, gradPhi, nablaphi);
//         msh->_finiteElement[ielt][solType]->Jacobian(vxHat, xi, weightHat, phiHat, gradPhiHat);
// 
// 
// 
// 
//         // BEGIN EVALUATION Quantities at the particles
//         std::vector <double> solVpSOld(dim);
//         particlesBulk[iBmarker]->GetMarkerVelocity(solVpSOld);
// 
//         std::vector <double> solApOld(dim);
//         particlesBulk[iBmarker]->GetMarkerAcceleration(solApOld);
// 
//         vector<adept::adouble> solDp(dim, 0.);
//         vector<adept::adouble> solVp(dim, 0.);
//         vector<adept::adouble> solVpOld(dim, 0.);
//         vector<adept::adouble> solVpTheta(dim, 0.);
// 
//         vector<vector < adept::adouble > > gradSolVpNew(dim);
//         vector<vector < adept::adouble > > gradSolVpTheta(dim);
//         vector<vector < adept::adouble > > DeltaSolVpTheta(dim);
//         vector<vector < adept::adouble > > gradSolDpHat(dim);
//         vector<vector < adept::adouble > > gradSolDpHatNew(dim);
// 
//         for(int j = 0; j < dim; j++) {
//           gradSolVpNew[j].assign(dim, 0.);
//           gradSolVpTheta[j].assign(dim, 0.);
//           gradSolDpHat[j].assign(dim, 0.);
//           gradSolDpHatNew[j].assign(dim, 0.);
//           DeltaSolVpTheta[j].assign(dim2, 0.);
//         }
// 
//         for(int j = 0; j < dim; j++) {
//           for(unsigned i = 0; i < nDofs; i++) {
//             solDp[j] += phi[i] * solDb[j][i];
//             solVp[j] += phi[i] * solV[j][i];
//             solVpOld[j] += phi[i] * solVOld[j][i];
//             for(int k = 0; k < dim; k++) {
//               gradSolVpNew[j][k] +=  gradPhiNew[i * dim + k] * solV[j][i];
//               gradSolVpTheta[j][k] +=  gradPhi[i * dim + k] * (theta * solV[j][i] + (1. - theta) * solVOld[j][i]);
//               gradSolDpHat[k][j] += (1. - af) * solDb[k][i] * gradPhiHat[i * dim + j];
//               gradSolDpHatNew[k][j] += solDb[k][i] * gradPhiHat[i * dim + j];
//             }
//           }
//         }
//         for(unsigned i = 0; i < nDofs; i++) {
//           for(unsigned j = 0; j < dim2; j++) {
//             for(unsigned  k = 0; k < dim; k++) {
//               DeltaSolVpTheta[k][j]   += nablaphi[i * dim2 + j] * (theta * solV[k][i] + (1. - theta) * solVOld[k][i]) ; // laplace of the theta velocity with respect to the theta domain
//             }
//           }
//         }
// 
// 
//         std::vector<adept::adouble> solAp(dim);
//         std::vector < adept::adouble > solApAm(dim);
//         std::vector < adept::adouble > solVpS(dim);
// 
//         for(int j = 0; j < dim; j++) {
//           solVpTheta[j] = theta * solVp[j] + (1. - theta) * solVpOld[j]; //TODO this comes from the particle, try the old mesh velocity
//           solAp[j] = (solDp[j] - 0.) / (beta * dt * dt) - solVpSOld[j] / (beta * dt) + solApOld[j] * (beta - 0.5) / beta ;   //NEWMARK ACCELERATION
//           solVpS[j] = solVpSOld[j] + dt * (Gamma * solAp[j] + (1. - Gamma) * solApOld[j]);   //velocity from the solid at xp, gamma configuration
//           solApAm[j] = (1. - am) * solAp[j] + am * solApOld[j];
//         }
// 
//         //Here we missed the if for piecewise  linear discontinuous //TODO
//         adept::adouble solPg = 0.;
//         vector<adept::adouble> gradSolPg(dim, 0.);
//         for(unsigned i = 0; i < nDofsP; i++) {
//           solPg += phiP[i] * solP[i];
//           for(unsigned k = 0; k < dim; k++) {
//             gradSolPg[k] += solP[i] * gradPhiP[i * dim + k];
//           }
//         }
// 
// 
// 
// 
//         std::vector <std::vector <adept::adouble> > G(dim); // J^T . J //centered at theta
//         for(unsigned i = 0; i < dim; i++) {
//           G[i].assign(dim, 0.);
//           for(unsigned j = 0; j < dim; j++) {
//             for(unsigned k = 0; k < dim; k++) {
//               G[i][j] += JacMatrix[k][i] * JacMatrix[k][j];
//             }
//           }
//         }
// 
//         double CI = 36.;
//         adept::adouble denom = pow(2 * rhoFluid / dtMin, 2.);
//         for(unsigned i = 0; i < dim; i++) {
//           for(unsigned j = 0; j < dim; j++) {
//             denom += rhoFluid * (solVpTheta[i] - (solDp[i] - 0.) / dt) * G[i][j] * rhoFluid * (solVpTheta[j] - (solDp[j] - 0.) / dt) // we can get the mesh velocity at af
//                      + CI * muFluid * muFluid * G[i][j] * G[i][j]; //this could be improved if we had the old mesh velocity
//           }
//         }
//         adept::adouble tauM = 1. / sqrt(denom);
// 
// 
// 
// 
// 
//         //BEGIN computation of the Cauchy Stress
//         std::vector < std::vector < double > > FpOld;
//         FpOld = particlesBulk[iBmarker]->GetDeformationGradient(); //extraction of the deformation gradient
// 
//         adept::adouble FpNew[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
//         adept::adouble FpNewNew[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
//         adept::adouble F[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
//         adept::adouble FNew[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
//         adept::adouble Id2th[3][3] = {{ 1., 0., 0.}, { 0., 1., 0.}, { 0., 0., 1.}};
//         adept::adouble Cauchy[3][3];
// 
//         for(unsigned j = 0; j < dim; j++) {
//           for(unsigned k = 0; k < dim; k++) {
//             FpNew[j][k] += gradSolDpHat[j][k]; // with respect to the reference deformed configuration at alpha_f
//             FpNewNew[j][k] += gradSolDpHatNew[j][k]; // with respect to the reference deformed configuration at alpha_f
//           }
//         }
// 
//         for(unsigned i = 0; i < dim; i++) {
//           for(unsigned j = 0; j < dim; j++) {
//             for(unsigned k = 0; k < dim; k++) {
//               F[i][j] += FpNew[i][k] * FpOld[k][j];
//               FNew[i][j] += FpNewNew[i][k] * FpOld[k][j];
//             }
//           }
//         }
// 
//         if(dim == 2) {
//           F[2][2] = 1.;
//           FNew[2][2] = 1.;
//         }
// 
//         adept::adouble J_hatOld;
//         if(dim == 2) {
//           J_hatOld = FpOld[0][0] * FpOld[1][1] - FpOld[0][1] + FpOld[1][0];
//         }
//         else {
//           J_hatOld = FpOld[0][0] * FpOld[1][1] * FpOld[2][2] + FpOld[0][1] * FpOld[1][2] * FpOld[2][0] + FpOld[0][2] * FpOld[1][0] * FpOld[2][1]
//                      - FpOld[2][0] * FpOld[1][1] * FpOld[0][2] - FpOld[2][1] * FpOld[1][2] * FpOld[0][0] - FpOld[2][2] * FpOld[1][0] * FpOld[0][1];
//         }
// 
//         adept::adouble J_hat =  F[0][0] * F[1][1] * F[2][2] + F[0][1] * F[1][2] * F[2][0] + F[0][2] * F[1][0] * F[2][1]
//                                 - F[2][0] * F[1][1] * F[0][2] - F[2][1] * F[1][2] * F[0][0] - F[2][2] * F[1][0] * F[0][1];
// 
//         adept::adouble J_hatNew =  FNew[0][0] * FNew[1][1] * FNew[2][2] + FNew[0][1] * FNew[1][2] * FNew[2][0] + FNew[0][2] * FNew[1][0] * FNew[2][1]
//                                    - FNew[2][0] * FNew[1][1] * FNew[0][2] - FNew[2][1] * FNew[1][2] * FNew[0][0] - FNew[2][2] * FNew[1][0] * FNew[0][1];
// 
// 
//         if(NeoHookean) {
//           adept::adouble B[3][3];
//           for(unsigned i = 0; i < 3; i++) {
//             for(int j = 0; j < 3; j++) {
//               B[i][j] = 0.;
//               for(unsigned k = 0; k < 3; k++) {
//                 //left Cauchy-Green deformation tensor or Finger tensor (B = F*F^T)
//                 B[i][j] += F[i][k] * F[j][k];
//               }
//             }
//           }
// 
//           adept::adouble I1_B = B[0][0] + B[1][1] + B[2][2];
// 
//           for(unsigned j = 0; j < 3; j++) {
//             for(unsigned k = 0; k < 3; k++) {
//               Cauchy[j][k] = lambdaMpm * log(J_hat) / J_hat * Id2th[j][k] + muMpm / J_hat * (B[j][k] - Id2th[j][k]);     //alternative formulation
//             }
//           }
//           //END computation of the Cauchy Stress
//         }
// 
//         else {
//           adept::adouble E[3][3];
//           adept::adouble S[3][3];
// 
//           for(unsigned i = 0; i < 3; i++) { //E = 0.5(F^T F - I)
//             for(unsigned j = 0; j < 3; j++) {
//               E[i][j] = 0.;
//               for(unsigned k = 0; k < 3; k++) {
//                 E[i][j] += F[k][i] * F[k][j];
//               }
//               E[i][j] = 0.5 * (E[i][j] - Id2th[i][j]);
//             }
//           }
// 
//           adept::adouble traceE = E[0][0] + E[1][1] + E[2][2];
// 
//           for(unsigned i = 0; i < 3; i++) { // S = lambda Tr(E) +  2 mu E
//             for(unsigned j = 0; j < 3; j++) {
//               S[i][j] = lambdaMpm * traceE * Id2th[i][j] + 2. * muMpm * E[i][j];     //alternative formulation
//             }
//           }
// 
//           adept::adouble SFt[3][3];
//           for(unsigned i = 0; i < 3; i++) { // S F^t
//             for(unsigned j = 0; j < 3; j++) {
//               SFt[i][j] = 0.;
//               for(unsigned k = 0; k < 3; k++) {
//                 SFt[i][j] += S[i][k] * F[j][k];
//               }
//             }
//           }
// 
//           for(unsigned i = 0; i < 3; i++) { // 1./J F S F^t
//             for(unsigned j = 0; j < 3; j++) {
//               Cauchy[i][j] = 0.;
//               for(unsigned k = 0; k < 3; k++) {
//                 Cauchy[i][j] += F[i][k] * SFt[k][j] / J_hat;
//               }
//             }
//           }
//         }
// 
// 
// 
//         //BEGIN Navier-Stokes in the bulk interface cells (integration is on the particles in \Omega_f)
//         if((1. - U) > 0 && eFlag == 1) {
// 
//           adept::adouble weightOld = (1 - U) * area * J_hatOld; // we need a * J_hat, to add also in the paper
//           adept::adouble weight = (1 - U) * area * J_hat; // we need a * J_hat, to add also in the paper
//           adept::adouble weightNew = (1 - U) * area * J_hatNew; // we need a * J_hat, to add also in the paper
// 
//           for(unsigned i = 0; i < nDofs; i++) {
//             for(unsigned k = 0; k < dim; k++) {
//               adept::adouble Vlaplace = 0.;
//               adept::adouble advection = 0.;
//               for(unsigned j = 0; j < dim; j++) {
//                 Vlaplace  +=  gradPhi[i * dim + j] * (gradSolVpTheta[k][j] + gradSolVpTheta[j][k]);
//                 advection +=  phi[i] * (solVpTheta[j] - (solDp[j] - 0.) / dt) * gradSolVpTheta[k][j]; //ALE
//               }
// 
//               aResV[k][i] += rhoFluid * phi[i] * (solVp[k] * weightNew - solVpOld[k] * weightOld) / dt +
//                              (rhoFluid * advection
//                               +
//                               muFluid * Vlaplace
//                               - weakP * gradPhi[i * dim + k] * solPg
//                               + !weakP * phi[i] * gradSolPg[k]
//                              ) * weight;
//             }
//           }
// 
//           for(unsigned i = 0; i < nDofsP; i++) {
//             if(eFlag == 1) {
// 
//               for(unsigned  k = 0; k < dim; k++) {
// 
//                 adept::adouble sLaplace = 0.;
//                 adept::adouble advection = 0.;
// 
//                 for(unsigned j = 0; j < dim; j++) {
//                   unsigned kdim;
// 
//                   if(k == j) kdim = j;
//                   else if(1 == k + j) kdim = dim;       // xy
//                   else if(2 == k + j) kdim = dim + 2;   // xz
//                   else if(3 == k + j) kdim = dim + 1;   // yz
// 
//                   sLaplace += (- muFluid * (DeltaSolVpTheta[k][j] + DeltaSolVpTheta[j][kdim]));
//                   advection += rhoFluid * (solVpTheta[j] - (solDp[j] - 0.) / dt) * gradSolVpTheta[k][j];
// 
//                 }
// 
//                 aResP[i] += phiP[i] *  gradSolVpNew[k][k] * weightNew;
// //                             + (rhoFluid * (solVp[k] - solVpOld[k]) / dt + advection +
// //                                sLaplace +  gradSolPg[k]) * tauM * gradPhiP[i * dim + k]
// //                             * weight;
// 
// 
// 
//               }
//             }
//           }
//         }
// 
//         //BEGIN solid Momentum in the bulk interface and solid cells (integration is on the particles in \Omega_s)
//         if(U > 0) {
//           double dM = (eFlag == 1) ? U * area * rhoMpm : area * rhoMpm;
//           for(unsigned i = 0; i < nDofs; i++) {
//             adept::adouble CauchyDIR[3] = {0., 0., 0.};
//             for(unsigned j = 0.; j < dim; j++) {
//               for(unsigned k = 0.; k < dim; k++) {
//                 CauchyDIR[j] += gradPhi[i * dim + k] * Cauchy[j][k];
//               }
//             }
// 
//             for(unsigned k = 0; k < dim; k++) {
// 
//               aResD[k][i] += (phi[i] * solApAm[k] + J_hat * CauchyDIR[k] / rhoMpm - gravity[k] * phi[i])  * dM;
// 
//               if(nodeFlag[i] == 0) { //bulk solid nodes: kinematic: v - dD/dt = 0
//                 aResV[k][i] += -phiHat[i] * (solVp[k] - solVpS[k]) * area; //TODO
//               }
// 
//             }
//           }
//         }
//         iBmarker++;
//       }
//     }
//     //END BULK PARTICLE
// 
//     if(true) {
// 
//       //BEGIN INTERFACE PARTICLE
// 
//       if(eFlag == 1) {  //interface markers
// 
//         double h = sqrt(dim) * sqrt((vxHat[0][0] - vxHat[0][1]) * (vxHat[0][0] - vxHat[0][1]) +
//                                     (vxHat[1][0] - vxHat[1][1]) * (vxHat[1][0] - vxHat[1][1])) ;
// 
//         double thetaI = 1.; // t_{n + 1}
//         double afI = 0.; // t_{n + 1}
//         double afN = af; // t_{n + alpha_f}
// 
// 
// 
//         while(imarkerI < markerOffsetI[iproc + 1] && iel != particleI[imarkerI]->GetMarkerElement()) {
//           imarkerI++;
//         }
// 
//         while(imarkerI < markerOffsetI[iproc + 1] && iel == particleI[imarkerI]->GetMarkerElement()) {
// 
//           // the local coordinates of the particles are the Gauss points in this context
//           std::vector <double> xi = particleI[imarkerI]->GetMarkerLocalCoordinates();
//           msh->_finiteElement[ielt][solType]->Jacobian(vxHat, xi, weightHat, phiHat, gradPhiHat);
//           msh->_finiteElement[ielt][solType]->Jacobian(vx, xi, weight, phi, gradPhi);
// 
// 
//           vector<vector < adept::adouble > > gradSolDpHat(dim);
//           for(int j = 0; j < dim; j++) {
//             gradSolDpHat[j].assign(dim, 0.);
//           }
//           for(int j = 0; j < dim; j++) {
//             for(unsigned i = 0; i < nDofs; i++) {
//               for(int k = 0; k < dim; k++) {
//                 gradSolDpHat[k][j] += ((1. - afN) * solDb[k][i]) * gradPhiHat[i * dim + j];
//               }
//             }
//           }
// 
//           std::vector <std::vector < double > > THat;
//           particleI[imarkerI]->GetMarkerTangent(THat);
// 
//           std::vector < std::vector < double > > FpOld;
//           FpOld = particleI[imarkerI]->GetDeformationGradient(); //extraction of the deformation gradient
//           adept::adouble FpNew[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
//           adept::adouble F[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
// 
//           for(unsigned j = 0; j < dim; j++) {
//             for(unsigned k = 0; k < dim; k++) {
//               FpNew[j][k] += gradSolDpHat[j][k]; // with respect to the reference deformed configuration at alpha_f
//             }
//           }
// 
//           for(unsigned i = 0; i < dim; i++) {
//             for(unsigned j = 0; j < dim; j++) {
//               for(unsigned k = 0; k < dim; k++) {
//                 F[i][j] += FpNew[i][k] * FpOld[k][j];
//               }
//             }
//           }
// 
//           std::vector <std::vector < adept::adouble > > T;
//           T.resize(THat.size());
// 
//           for(unsigned k = 0; k < T.size(); k++) {
//             T[k].assign(dim, 0.);
//             for(unsigned i = 0; i < dim; i++) {
//               for(unsigned j = 0; j < dim; j++) {
//                 T[k][i] += F[i][j] * THat[k][j]; // can be improved
//               }
//             }
//           }
// 
//           adept::adouble weight;
//           std::vector < adept::adouble > N(dim);
//           if(dim == 2) {
//             N[0] =  T[0][1];
//             N[1] = -T[0][0];
//             weight = sqrt(N[0] * N[0] + N[1] * N[1]);
//             N[0] /= weight;
//             N[1] /= weight;
//           }
//           else {
//             N[0] = T[0][1] * T[1][2] - T[0][2] * T[1][1];
//             N[1] = T[0][2] * T[1][0] - T[0][0] * T[1][2];
//             N[2] = T[0][0] * T[1][1] - T[0][1] * T[1][0];
//             weight = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
//             N[0] /= weight;
//             N[1] /= weight;
//             N[2] /= weight;
//           }
// 
//           for(unsigned k = 0; k < dim; k++) {
//             N[k] *= -1.; // fluid to solid normal
//           }
// 
//           msh->_finiteElement[ielt][solTypeP]->GetPhi(phiP, xi);
//           adept::adouble solPp = 0.;
//           for(unsigned i = 0; i < nDofsP; i++) {
//             solPp += phiP[i] * solP[i];
//           }
// 
//           std::vector < adept::adouble > v1(dim, 0.); //fluid velocity in thetaI configuration
//           std::vector < adept::adouble > v2(dim, 0.); //solid velocity in afI configuration
// 
//           std::vector < adept::adouble > tau(dim, 0.);
// 
//           std::vector <double> solVpOld(dim);
//           particleI[imarkerI]->GetMarkerVelocity(solVpOld);
// 
//           std::vector <double> solApOld(dim);
//           particleI[imarkerI]->GetMarkerAcceleration(solApOld);
// 
//           std::vector <adept::adouble> solDp(dim, 0.);
//           std::vector <adept::adouble> solAp(dim);
// 
//           //update displacement and acceleration
//           for(int k = 0; k < dim; k++) {
//             for(unsigned i = 0; i < nDofs; i++) {
//               v1[k] += phi[i] * (thetaI * solV[k][i] + (1. - thetaI) * solVOld[k][i]);
//               solDp[k] += phi[i] * solDb[k][i];
//             }
//             solAp[k] = (solDp[k] - 0.) / (beta * dt * dt) - solVpOld[k] / (beta * dt) + (beta - 0.5) / beta * solApOld[k]; // Newmark acceleration
//             v2[k] = solVpOld[k] + (1. - afI) * (dt * (Gamma * solAp[k] + (1. - Gamma) * solApOld[k]));
//           }
// 
//           for(unsigned k = 0; k < dim; k++) {
//             tau[k] += solPp * N[k];
//             for(unsigned i = 0; i < nDofs; i++) {
//               for(unsigned j = 0; j < dim; j++) {
//                 tau[k] += -muFluid * ((theta * solV[k][i] + (1. - theta) * solVOld[k][i]) * gradPhi[i * dim + j] +
//                                       (theta * solV[j][i] + (1. - theta) * solVOld[j][i]) * gradPhi[i * dim + k]) * N[j];
//               }
//             }
//           }
// 
//           double c = 0.;
//           for(unsigned k = 0; k < dim; k++) {
//             c += (v1[k].value() - v2[k].value()) * (v1[k].value() - v2[k].value());
//           }
//           c = sqrt(c);
// 
//           double thetaM = GAMMA * muFluid / h; //  [rho L/ T]
//           double thetaL = GAMMA * rhoFluid * ((c / 6.) + h / (12. * theta * dtMin)) + thetaM;  // [rho L/T]
// 
//           // *** phi_i loop ***
//           for(unsigned i = 0; i < nDofs; i++) {
//             for(unsigned k = 0; k < dim; k++) {
//               aResV[k][i] += (tau[k] - !weakP * solPp * N[k]) * phi[i] * weight;  // correct sign due to the normal
//               aResV[k][i] += thetaM * (v1[k] - v2[k]) * phi[i] * weight;
// 
//               aResD[k][i] += -tau[k] * phi[i] * weight; // correct sign due to the normal
//               aResD[k][i] +=  thetaM * (v1[k] - v2[k]) * (-phi[i]) * weight;
// 
//               for(unsigned j = 0; j < dim; j++) {
//                 aResV[k][i] += - muFluid * gradPhi[i * dim + j] * N[j] * (v1[k] - v2[k]) * weight;
//                 aResV[k][i] += - muFluid * gradPhi[i * dim + j] * N[k] * (v1[j] - v2[j]) * weight;
// 
//                 aResV[k][i] += thetaL * (v1[j] - v2[j]) * N[j] * phi[i] * N[k] * weight;
//                 aResD[k][i] += thetaL * (v1[j] - v2[j]) * N[j] * (-phi[i]) * N[k] * weight;
// 
//               }
//             }
//           } // end phi_i loop
// 
// 
//           for(unsigned i = 0; i < nDofsP; i++) {
//             for(unsigned k = 0; k < dim; k++) {
//               aResP[i] += - (phiP[i]  * (v1[k] - v2[k]) * N[k]) * weight;
//             }
//           }
//           imarkerI++;
//         }
//       }
//     }
//     //END INTERFACE PARTICLES

    //copy the value of the adept::adoube aRes in double Res and store them in RES
    rhs.resize(nDofsAll);   //resize

    for(int i = 0; i < nDofs; i++) {
      for(unsigned  k = 0; k < dim; k++) {
        rhs[ i +  k * nDofs ] = -aResD[k][i].value();
        rhs[ i + (k + dim) * nDofs ] = -aResV[k][i].value();
      }
    }
    for(int i = 0; i < nDofsP; i++) {
      rhs[ i + (2 * dim) * nDofs ] = -aResP[i].value();
    }

    myRES->add_vector_blocked(rhs, sysDofsAll);


    Jac.resize(nDofsAll * nDofsAll);
    // define the dependent variables

    for(unsigned  k = 0; k < dim; k++) {
      s.dependent(&aResD[k][0], nDofs);
    }
    for(unsigned  k = 0; k < dim; k++) {
      s.dependent(&aResV[k][0], nDofs);
    }
    s.dependent(&aResP[0], nDofsP);

    // define the independent variables
    for(unsigned  k = 0; k < dim; k++) {
      s.independent(&solDb[k][0], nDofs);
    }
    for(unsigned  k = 0; k < dim; k++) {
      s.independent(&solVNew[k][0], nDofs);
    }
    s.independent(&solP[0], nDofsP);

    // get the and store jacobian matrix (row-major)
    s.jacobian(&Jac[0], true);
    myKK->add_matrix_blocked(Jac, sysDofsAll, sysDofsAll);

    s.clear_independents();
    s.clear_dependents();

  }

  myRES->close();
  myKK->close();


//   PetscViewer    viewer1;
//   PetscViewerDrawOpen (PETSC_COMM_WORLD, NULL, NULL, 0, 0, 1800, 1800, &viewer1);
//   PetscObjectSetName ( (PetscObject) viewer1, "FSI matrix");
//   PetscViewerPushFormat (viewer1, PETSC_VIEWER_DRAW_LG);
//   MatView ( (static_cast< PetscMatrix* > (myKK))->mat(), viewer1);
// 
//   double a;
//   std::cin >> a;

// *************************************
  end_time = clock();
  AssemblyTime += (end_time - start_time);
// ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

}


