
#include "./include/GhostPenalty.hpp"

bool SetBoundaryConditionB(const std::vector < double >&x, const char name[], double &value, const int facename, const double t);
double SetVariableTimeStepB(const double time);

double InitalValue0(const std::vector < double >& x) {
  return 0.;
}


void InitializeBackgroundVariables(MultiLevelSolution &mlSol) {

  unsigned dim = mlSol._mlMesh->GetDimension();

  FEOrder femOrder = SECOND;

  //add variables to mlSol
  mlSol.AddSolution("DX", LAGRANGE, femOrder, 2);
  mlSol.AddSolution("DY", LAGRANGE, femOrder, 2);
  if(dim == 3) mlSol.AddSolution("DZ", LAGRANGE, femOrder, 2);
  mlSol.AddSolution("VX", LAGRANGE, femOrder, 2);
  mlSol.AddSolution("VY", LAGRANGE, femOrder, 2);
  if(dim == 3) mlSol.AddSolution("VZ", LAGRANGE, femOrder, 2);

  mlSol.AddSolution("P", LAGRANGE, FIRST, 2);

  mlSol.AddSolution("eflag", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
  mlSol.AddSolution("fldCnt", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
  mlSol.AddSolution("sldCnt", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
  mlSol.AddSolution("nflag", LAGRANGE, SECOND, 0, false);

  //mlSol.Initialize("All");
  
  
  mlSol.Initialize("DX", InitalValue0);
  mlSol.Initialize("DY", InitalValue0);
  
  mlSol.Initialize("VX", InitalValue0);
  mlSol.Initialize("VY", InitalValue0);
  
  mlSol.Initialize("P", InitalValue0);
  mlSol.Initialize("eflag", InitalValue0);
  mlSol.Initialize("fldCnt", InitalValue0);
  mlSol.Initialize("sldCnt", InitalValue0);
  mlSol.Initialize("nflag", InitalValue0);
  

  mlSol.AttachSetBoundaryConditionFunction(par->_bdcFunction);

  mlSol.SetIfFSI(true);
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

  MatResetPreallocation((static_cast< PetscMatrix* >(myKK))->mat());
  MatSetOption((static_cast< PetscMatrix* >(myKK))->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  
  myKK->zero();
  myRES->zero();

  AssembleGhostPenaltyP(ml_prob, true);
  AssembleGhostPenaltyP(ml_prob, false);


  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  const unsigned dim = msh->GetDimension();

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  // data
  unsigned iproc  = msh->processor_id();
  unsigned nprocs  = msh->n_processors();

  vector< vector< adept::adouble > > solD(dim);   // background displacement at n+1
  vector< adept::adouble > solP; // fluid pressure at theta

  vector< vector< double > > solVOld(dim); // fluid velocity at n
  vector< vector< adept::adouble > > solVNew(dim); // fluid velocity at n+1


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

  double gaussWeight;
  adept::adouble agaussWeight;


  //reading parameters for MPM body
  double rhoMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_density();
  double EMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_young_module();
  double muMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_lame_shear_modulus();
  double nuMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_poisson_coeff();
  double lambdaMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_lame_lambda();

  //reading parameters for fluid FEM domain
  double rhoFluid = ml_prob.parameters.get<Fluid> ("FluidFEM").get_density();
  double muFluid = ml_prob.parameters.get<Fluid> ("FluidFEM").get_viscosity();

  std::cout << "mu_s = " << muMpm    << " lambda_s= " << lambdaMpm << " nu_s = " << nuMpm << std::endl;
  std::cout << "rho_s = "<< rhoMpm   << " E_s = " << EMpm << std::endl;
  std::cout << "rho_f = "<< rhoFluid << " mu_f = " << muFluid << std::endl;

  std::cout << "a_f = " << par->_af << " a_m = " << par->_am << std::endl;
  std::cout << "gamma = " << par->_gamma << " beta = " << par->_beta << std::endl;
  
  double dt =  my_nnlin_impl_sys.GetIntervalTime();
  double time =  my_nnlin_impl_sys.GetTime();

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
  
  unsigned meshType = 2;
  std::vector < unsigned >  nodeFlag; // local solution

  start_time = clock();

  const std::vector<std::vector<unsigned> > & ielp = projection->GetIel();
  const std::vector < std::vector < std::vector <double > > > & Xip = projection->GetXi();
  const std::vector<std::vector<unsigned> > & mtypep = projection->GetMtype();
  const std::vector<std::vector<double> > & weightpOld = projection->GetWeight();
  const std::vector < std::vector < std::vector <double > > > & VpOld = projection->GetV();
  const std::vector < std::vector < std::vector <double > > > & ApOld = projection->GetA();
  const std::vector < std::vector<std::vector<double> > > & NpOld = projection->GetN();
  const std::vector < std::vector < std::vector < std::vector <double > > > > & gradDpOld = projection->GetGradD();
  const std::vector < std::vector < std::vector < std::vector <double > > > > & FFpOld = projection->GetF();

  std::vector < unsigned > im(nprocs, 0);
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielt = msh->GetElementType(iel);

    unsigned nDofs = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solTypeP);  // number of pressure dofs

    unsigned nDofsMesh = msh->GetElementDofNumber(iel, meshType);


    unsigned nDofsAll = 2 * dim * nDofs + nDofsP;

    // resize local arrays
    sysDofsAll.resize(nDofsAll);

    nodeFlag.resize(nDofsMesh);
    double tempEflag = (*mysolution->_Sol[eflagIndex])(iel);
    unsigned eFlag = static_cast <unsigned>(floor(tempEflag + 0.2));
    unsigned eFlag1 = 2; // interface or solid
    if(eFlag == 0) {
      eFlag1 = (tempEflag < 0.25) ?  0 : 1; //fluid or fluid-particle
    }


    for(unsigned  k = 0; k < dim; k++) {
      solD[k].resize(nDofs);

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

    for(unsigned i = 0; i < nDofsMesh; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, meshType);
      nodeFlag[i] = (*mysolution->_Sol[nflagIndex])(idof);
    }


    for(unsigned i = 0; i < nDofs; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, solType);

      for(unsigned  k = 0; k < dim; k++) {
        solD[k][i] = (*mysolution->_Sol[indexSolD[k]])(idof); //t_{n+1} -t_n

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
      unsigned idofX = msh->GetSolutionDof(i, iel, meshType);
      for(unsigned  k = 0; k < dim; k++) {
        vxOld[k][i] = (*msh->_topology->_Sol[k])(idofX); // undeformed background configuration at t_{n}
        vx[k][i]  = vxOld[k][i] + (1. - par->_af) * solD[k][i]; // deformed background configuration at alpha_f/theta
        vxNew[k][i]  = vxOld[k][i] + solD[k][i]; // deformed background configuration at t_{n+1}
      }
    }

    double elementArea = 0.;

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielt][solType]->GetGaussPointNumber(); ig++) {

      msh->_finiteElement[ielt][solTypeP]->Jacobian(vx, ig, weight, phiP, gradPhiP);

      msh->_finiteElement[ielt][solType]->Jacobian(vxOld, ig, weightOld, phi, gradPhiOld);
      msh->_finiteElement[ielt][solType]->Jacobian(vx,    ig, weight,    phi, gradPhi, nablaphi);
      msh->_finiteElement[ielt][solType]->Jacobian(vxNew, ig, weightNew, phi, gradPhiNew);

      elementArea += weightOld;

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
          solDg[j] += phi[i] * solD[j][i]; // new displacement
          for(unsigned  k = 0; k < dim; k++) {
            gradSolDgOld[k][j] += gradPhiOld[i * dim + j] * solD[k][i]; //gradient of new solution with respect to deformed reference configuration
            gradSolVgNew[k][j] += gradPhiNew[i * dim + j] * solVNew[k][i]; // gradient of the new velocity with respect to the theta domain
            gradSolVg[k][j] += gradPhi[i * dim + j] * (par->_theta * solVNew[k][i] + (1. - par->_theta) * solVOld[k][i]); // gradient of the theta velocity with respect to the theta domain
          }
        }
        for(unsigned j = 0; j < dim; j++) {
          solVg[j] += par->_theta * solVgNew[j] + (1. - par->_theta) * solVgOld[j]; //velocity in the deformed theta configuration
        }

        for(unsigned j = 0; j < dim2; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            DeltaSolVg[k][j]   += nablaphi[i * dim2 + j] * (par->_theta * solVNew[k][i] + (1. - par->_theta) * solVOld[k][i]) ; // laplace of the theta velocity with respect to the theta domain
          }
        }
      }

      adept::adouble solPg = 0.;
      for(unsigned i = 0; i < nDofsP; i++) {
        solPg += phiP[i] * solP[i];
      }

      if(eFlag == 0) { //bulk fluid: all this equation is centered at a_f, There are no time derivatives

//         adept::adouble F[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
//         for(unsigned j = 0; j < dim; j++) {
//           for(unsigned k = 0; k < dim; k++) {
//             F[j][k] += gradSolDgOld[j][k]; //tilde F in the ALE equation
//           }
//         }
//
//         adept::adouble J_Old =  F[0][0] * F[1][1] * F[2][2] + F[0][1] * F[1][2] * F[2][0] + F[0][2] * F[1][0] * F[2][1]
//                                 - F[2][0] * F[1][1] * F[0][2] - F[2][1] * F[1][2] * F[0][0] - F[2][2] * F[1][0] * F[0][1];
//
//         adept::adouble B[3][3];
//         for(unsigned i = 0; i < 3; i++) {
//           for(int j = 0; j < 3; j++) {
//             B[i][j] = 0.;
//             for(unsigned k = 0; k < 3; k++) {
//               //left Cauchy-Green deformation tensor or Finger tensor (B = F*F^T)
//               B[i][j] += F[i][k] * F[j][k];
//             }
//           }
//         }
//
//         adept::adouble I1_B = B[0][0] + B[1][1] + B[2][2];
//         adept::adouble Id2th[3][3] = {{ 1., 0., 0.}, { 0., 1., 0.}, { 0., 0., 1.}};
//         adept::adouble sigma[3][3];
//
//         double E = pow(10., eFlag1);
//         double nu = 0.4;
//
//         double mu = E / (2. * (1. + nu));
//         double lambda = (E * nu) / ((1. + nu) * (1. - 2.*nu));
//
//         for(unsigned j = 0; j < 3; j++) {
//           for(unsigned k = 0; k < 3; k++) {
//             sigma[j][k] = lambda * log(J_Old) / J_Old * Id2th[j][k] + mu / J_Old * (B[j][k] - Id2th[j][k]);    // alternative formulation
//           }
//         }
//         //END computation of the Cauchy Stress
//         for(unsigned i = 0; i < nDofs; i++) {//Auxiliary Equations
//           if(nodeFlag[i] == 0) {
//             for(unsigned k = 0.; k < dim; k++) {
//               adept::adouble cauchy = 0.;
//               for(unsigned j = 0.; j < dim; j++) {
//                 cauchy += sigma[k][j] * gradPhi[i * dim + j] ;
//                 //cauchy += gradSolDgOld[k][j] * gradPhi[i * dim + j] ;
//               }
//               aResD[k][i] += cauchy * weight;
//             }
//           }
//         }

        double E = pow(10., eFlag1);
        for(unsigned i = 0; i < nDofs; i++) {
          if(nodeFlag[i] == 0) {
            for(unsigned k = 0.; k < dim; k++) {
              adept::adouble wlaplaceD  = 0.;
              for(unsigned  j = 0; j < dim; j++) {
                wlaplaceD +=  E * gradPhiOld[i * dim + j] * (gradSolDgOld[k][j] + gradSolDgOld[j][k]);
              }
              aResD[k][i] +=  wlaplaceD * weightOld;
            }
          }
        }
      }


      if(eFlag == 0) {   // only fluid cells

        //start SUPG paramters, tauM, tauC, G to get tauM_SupgPhi
        std::vector <std::vector <adept::adouble> > Jac;
        std::vector <std::vector <adept::adouble> > JacI;
        msh->_finiteElement[ielt][solType]->GetJacobianMatrix(vx, ig, weight, Jac, JacI); //centered at theta


        if(solTypeP == 4) { //discontinuous pressure <1,\xi,\eta> bases centered at theta
          for(unsigned j = 0; j < dim; j++) {
            gradPhiP[0 * dim + j]  = 0.;
          }
          for(unsigned i = 0; i < dim; i++) {
            for(unsigned j = 0; j < dim; j++) {
              gradPhiP[(i + 1) * dim + j] = JacI[i][j];
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
              G[i][j] += JacI[k][i] * JacI[k][j];
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

        //tauM = 0.;
        //tauC = 0.;

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
                              - /*par->_weakP * */ gradPhi[i * dim + k] * solPg
//                               + !par->_weakP * phi[i] * gradSolPg[k]
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
                         + tauM * ((rhoFluid * (solVgNew[k] - solVgOld[k]) / dt + advection +
                                    sLaplace +  gradSolPg[k]) * gradPhiP[i * dim + k]
                                  ) * weight;

          }
        }
      }
      else if(eFlag == 2) {   // only solid cells: fake pressure //TODO
        if(solTypeP >= 3) {
          for(unsigned i = 0; i < nDofsP; i++) {
            aResP[i] += phiP[i] * solP[i] * weight;
          }
        }
        else {
          for(unsigned i = 0; i < nDofsP; i++) {
            if(nodeFlag[i] == 0) {
              aResP[i] += phiP[i] * solP[i] * weight;
            }
          }
        }


        for(unsigned i = 0; i < nDofs; i++) {
          if(nodeFlag[i] == 0) {
            for(unsigned k = 0.; k < dim; k++) {
              adept::adouble wlaplaceV  = 0.;
              for(unsigned  j = 0; j < dim; j++) {
                wlaplaceV +=  gradPhi[i * dim + j] * gradSolVg[k][j];
              }
              aResV[k][i] +=  1.0e-10 * wlaplaceV * weight;
            }
          }
        }


      }
    } // end gauss point loop



    if(eFlag > 0) {   //BEGIN BULK PARTICLE

      double particleArea = 0.;

      for(unsigned kp = 0; kp < nprocs; kp++) {
        //im[kp] = 0;
        while(im[kp] < ielp[kp].size() && ielp[kp][im[kp]] < iel) {
          im[kp]++;
        }
        unsigned im0 = im[kp];
        while(im[kp] < ielp[kp].size() && iel == ielp[kp][im[kp]]) {
          particleArea += weightpOld[kp][im[kp]];
          im[kp]++;
        }
        im[kp] = im0;
      }
      double scale = 1;//elementArea / particleArea;


      for(unsigned kp = 0; kp < nprocs; kp++) {
        while(im[kp] < ielp[kp].size() && iel == ielp[kp][im[kp]]) {

          // the local coordinates of the particles are the integration points in this context
          std::vector <double> xi(dim);
          for(unsigned  k = 0; k < dim; k++) xi[k] = Xip[kp][k][im[kp]];

          double U = mtypep[kp][im[kp]] / 2.; // U = 0 fluid, U = 0.5 interface, U = 1 solid
          double areaOld = scale * weightpOld[kp][im[kp]];

          msh->_finiteElement[ielt][solType]->Jacobian(vx, xi, agaussWeight, phi, gradPhi);
          msh->_finiteElement[ielt][solType]->Jacobian(vxOld, xi, gaussWeight, phi, gradPhiOld);

          // BEGIN EVALUATION Quantities at the particles need by both fluid and solid particles
          vector<adept::adouble> solDb(dim, 0.);
          vector<vector < adept::adouble > > gradSolDb(dim);
          vector<vector < adept::adouble > > gradSolDbNew(dim);

          for(unsigned j = 0; j < dim; j++) {
            gradSolDb[j].assign(dim, 0.);
            gradSolDbNew[j].assign(dim, 0.);
            for(unsigned i = 0; i < nDofs; i++) {
              solDb[j] += phi[i] * solD[j][i];
              for(unsigned k = 0; k < dim; k++) {
                gradSolDb[j][k] += gradPhiOld[i * dim + k] * (1. - par->_af) * solD[j][i]; // gradient centered at a_f with respect to the old background mesh
                gradSolDbNew[j][k] += gradPhiOld[i * dim + k] * solD[j][i]; // gradient centered at n+1 with respect to the old background mesh
              }
            }
          }

          adept::adouble Fb[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
          adept::adouble FbNew[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};

          for(unsigned j = 0; j < dim; j++) {
            for(unsigned k = 0; k < dim; k++) {
              Fb[j][k] += gradSolDb[j][k]; // with respect to the reference deformed configuration at alpha_f
              FbNew[j][k] += gradSolDbNew[j][k]; // with respect to the reference deformed configuration at n + 1
            }
          }

          adept::adouble Jb;
          adept::adouble JbNew;

          if(dim == 2) {
            Jb = Fb[0][0] * Fb[1][1] - Fb[0][1] * Fb[1][0];
            JbNew = FbNew[0][0] * FbNew[1][1] - FbNew[0][1] * FbNew[1][0];
          }
          else {
            Jb = Fb[0][0] * Fb[1][1] * Fb[2][2] + Fb[0][1] * Fb[1][2] * Fb[2][0] + Fb[0][2] * Fb[1][0] * Fb[2][1] -
                 Fb[2][0] * Fb[1][1] * Fb[0][2] - Fb[2][1] * Fb[1][2] * Fb[0][0] - Fb[2][2] * Fb[1][0] * Fb[0][1];
            JbNew = FbNew[0][0] * FbNew[1][1] * FbNew[2][2] + FbNew[0][1] * FbNew[1][2] * FbNew[2][0] + FbNew[0][2] * FbNew[1][0] * FbNew[2][1] -
                    FbNew[2][0] * FbNew[1][1] * FbNew[0][2] - FbNew[2][1] * FbNew[1][2] * FbNew[0][0] - FbNew[2][2] * FbNew[1][0] * FbNew[0][1];
          }

          vector<double> solVFpOld(dim, 0.); //old fluid velocity
          vector<adept::adouble> solVFp(dim, 0.); // fluid velocity
          vector<adept::adouble> solVFpNew(dim, 0.);  //new fluid velocity
          for(int j = 0; j < dim; j++) {
            for(unsigned i = 0; i < nDofs; i++) {
              solVFpOld[j] += phi[i] * solVOld[j][i];
              solVFpNew[j] += phi[i] * solVNew[j][i];
            }
          }
          for(int j = 0; j < dim; j++) {
            solVFp[j] = par->_theta * solVFpNew[j] + (1. - par->_theta) * solVFpOld[j];
          }


          //BEGIN Navier-Stokes in the bulk interface cells (integration is on the particles in \Omega_f)
          vector<vector < adept::adouble > > gradSolVp(dim);
          adept::adouble solPp = 0.;

          if(eFlag == 1 && (1. - U) > 0.1) {

            msh->_finiteElement[ielt][solType]->Jacobian(vxNew, xi, agaussWeight, phi, gradPhiNew);

            vector<vector < adept::adouble > > gradSolVpNew(dim);
            vector<adept::adouble> gradSolPp(dim, 0.);

            for(unsigned j = 0; j < dim; j++) {
              gradSolVp[j].assign(dim, 0.);
              gradSolVpNew[j].assign(dim, 0.);
              for(unsigned i = 0; i < nDofs; i++) {
                for(unsigned k = 0; k < dim; k++) {
                  gradSolVpNew[j][k] +=  gradPhiNew[i * dim + k] * solVNew[j][i];
                  gradSolVp[j][k] +=  gradPhi[i * dim + k] * (par->_theta * solVNew[j][i] + (1. - par->_theta) * solVOld[j][i]);
                }
              }
            }

            // Here we missed the if for piecewise  linear discontinuous //TODO
            // std::vector <std::vector <adept::adouble> > JacI;
            // std::vector <std::vector <adept::adouble> > Jac;
            // msh->_finiteElement[ielt][solType]->GetJacobianMatrix(vx, xi, Jac, JacI); //centered at theta

            msh->_finiteElement[ielt][solTypeP]->Jacobian(vx, xi, agaussWeight, phiP, gradPhiP);
            for(unsigned i = 0; i < nDofsP; i++) {
              solPp += phiP[i] * solP[i];
              for(unsigned k = 0; k < dim; k++) {
                gradSolPp[k] += solP[i] * gradPhiP[i * dim + k];
              }
            }

            double dAOld = (1. - U) * areaOld; // areaOld = area at n
            adept::adouble dA = (1. - U) * areaOld * Jb; // areaOld * Jb = area at theta
            adept::adouble dANew = (1. - U) * areaOld * JbNew; // areaOld * JbNew = area at n + 1

            for(unsigned i = 0; i < nDofs; i++) {
              for(unsigned k = 0; k < dim; k++) {
                adept::adouble laplace = 0.;
                adept::adouble advection = 0.;
                for(unsigned j = 0; j < dim; j++) {
                  laplace  +=  gradPhi[i * dim + j] * (gradSolVp[k][j] + gradSolVp[j][k]);
                  advection +=  phi[i] * (solVFp[j] - (solDb[j] - 0.) / dt) * gradSolVp[k][j]; //ALE
                }

                aResV[k][i] += rhoFluid * phi[i] * (solVFpNew[k] * dANew - solVFpOld[k] * dAOld) / dt +
                               (rhoFluid * advection
                                + muFluid * laplace
                                - /*par->_weakP **/ gradPhi[i * dim + k] * solPp
//                                 + !par->_weakP * phi[i] * gradSolPp[k]
                               ) * dA;
              }
            }

            for(unsigned i = 0; i < nDofsP; i++) {
              for(unsigned  k = 0; k < dim; k++) {
                aResP[i] += phiP[i] *  gradSolVpNew[k][k] * dANew;
              }
            }
          }


          //BEGIN solid Momentum in the bulk interface and solid cells (integration is on the particles in \Omega_s)
          std::vector < adept::adouble > solVSpNew(dim); // solid particle velocity
          std::vector <double> solVSpOld(dim); //old solid particle velocity
          std::vector <double> solApOld(dim);
          std::vector <adept::adouble> solApNew(dim);


          if(U > 0.1) {

            std::vector < adept::adouble > solAp(dim); //centered at am

            for(unsigned k = 0; k < dim; k++) {
              solVSpOld[k] = VpOld[kp][k][im[kp]];
              solApOld[k] = ApOld[kp][k][im[kp]];
            }

            for(int j = 0; j < dim; j++) {
              solApNew[j] = (solDb[j] - 0.) / (par->_beta * dt * dt) - solVSpOld[j] / (par->_beta * dt) + solApOld[j] * (par->_beta - 0.5) / par->_beta ;   //NEWMARK ACCELERATION
              solVSpNew[j] = solVSpOld[j] + dt * (par->_gamma * solApNew[j] + (1. - par->_gamma) * solApOld[j]);   //velocity from the solid at xp, gamma configuration
              solAp[j] = (1. - par->_am) * solApNew[j] + par->_am * solApOld[j];
              //solAp[j] = solApNew[j];
              //solVSp[j] = solAp[j] = 0.;
            }

            double FpOld[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
            for(unsigned j = 0; j < dim; j++) {
              for(unsigned k = 0; k < dim; k++) {
                FpOld[j][k] += gradDpOld[kp][j][k][im[kp]];
                //FpOld[j][k] = FFpOld[kp][j][k][im[kp]];
              }
            }

            adept::adouble Fp[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
            for(unsigned i = 0; i < dim; i++) {
              for(unsigned j = 0; j < dim; j++) {
                for(unsigned k = 0; k < dim; k++) {
                  Fp[i][j] += Fb[i][k] * FpOld[k][j];
                }
              }
            }

            if(dim == 2) {
              Fp[2][2] = 1.;
            }

            double JpOld;
            adept::adouble Jp;
            if(dim == 2) {
              JpOld = FpOld[0][0] * FpOld[1][1] - FpOld[0][1] * FpOld[1][0];
              Jp = Fp[0][0] * Fp[1][1] - Fp[0][1] * Fp[1][0];
            }
            else {
              JpOld = FpOld[0][0] * FpOld[1][1] * FpOld[2][2] + FpOld[0][1] * FpOld[1][2] * FpOld[2][0] + FpOld[0][2] * FpOld[1][0] * FpOld[2][1] -
                      FpOld[2][0] * FpOld[1][1] * FpOld[0][2] - FpOld[2][1] * FpOld[1][2] * FpOld[0][0] - FpOld[2][2] * FpOld[1][0] * FpOld[0][1];
              Jp = Fp[0][0] * Fp[1][1] * Fp[2][2] + Fp[0][1] * Fp[1][2] * Fp[2][0] + Fp[0][2] * Fp[1][0] * Fp[2][1] -
                   Fp[2][0] * Fp[1][1] * Fp[0][2] - Fp[2][1] * Fp[1][2] * Fp[0][0] - Fp[2][2] * Fp[1][0] * Fp[0][1];
            }

            //BEGIN computation of the Cauchy Stress
            adept::adouble Id2th[3][3] = {{ 1., 0., 0.}, { 0., 1., 0.}, { 0., 0., 1.}};
            adept::adouble Cauchy[3][3];
            if(par->_NeoHookean) {
              adept::adouble B[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
              for(unsigned i = 0; i < 3; i++) {
                for(unsigned j = 0; j < 3; j++) {
                  for(unsigned k = 0; k < 3; k++) {
                    //left Cauchy-Green deformation tensor or Finger tensor (B = F*F^T)
                    B[i][j] += Fp[i][k] * Fp[j][k];
                  }
                }
              }
              for(unsigned j = 0; j < dim; j++) {
                for(unsigned k = 0; k < dim; k++) {
                  Cauchy[j][k] += lambdaMpm * log(Jp) / Jp * Id2th[j][k] + muMpm / Jp * (B[j][k] - Id2th[j][k]);     //alternative formulation
                }
              }
            }
            else { //Saint Venant
              adept::adouble E[3][3];
              adept::adouble S[3][3];
              for(unsigned i = 0; i < 3; i++) { //E = 0.5(F^T F - I)
                for(unsigned j = 0; j < 3; j++) {
                  E[i][j] = 0.;
                  for(unsigned k = 0; k < 3; k++) {
                    E[i][j] += Fp[k][i] * Fp[k][j];
                  }
                  E[i][j] = 0.5 * (E[i][j] - Id2th[i][j]);
                }
              }
              adept::adouble traceE = E[0][0] + E[1][1] + E[2][2];
              for(unsigned i = 0; i < 3; i++) { // S = lambda Tr(E) +  2 mu E
                for(unsigned j = 0; j < 3; j++) {
                  S[i][j] = lambdaMpm * traceE * Id2th[i][j] + 2. * muMpm * E[i][j];     //alternative formulation
                }
              }
              adept::adouble SFt[3][3];
              for(unsigned i = 0; i < 3; i++) { // S F^t
                for(unsigned j = 0; j < 3; j++) {
                  SFt[i][j] = 0.;
                  for(unsigned k = 0; k < 3; k++) {
                    SFt[i][j] += S[i][k] * Fp[j][k];
                  }
                }
              }
              for(unsigned i = 0; i < 3; i++) { // 1./J F S F^t
                for(unsigned j = 0; j < 3; j++) {
                  Cauchy[i][j] = 0.;
                  for(unsigned k = 0; k < 3; k++) {
                    Cauchy[i][j] += Fp[i][k] * SFt[k][j] / Jp;
                  }
                }
              }
            }
            //END computation of the Cauchy Stress

            double dM = U * (areaOld / JpOld) * rhoMpm; // (areaOld / JpOld) = areaHat
            adept::adouble dA = U * areaOld * Jb;
            for(unsigned i = 0; i < nDofs; i++) {
              adept::adouble CauchyDIR[3] = {0., 0., 0.};
              for(unsigned j = 0.; j < dim; j++) {
                for(unsigned k = 0.; k < dim; k++) {
                  CauchyDIR[j] += gradPhi[i * dim + k] * Cauchy[j][k];
                }
              }
              //par->_gravity[1] = (time < 2) ? 5 * sin(M_PI / 2.*time) : 0.;

              for(unsigned k = 0; k < dim; k++) {
                //aResD[k][i] += (phi[i] * solAp[k] + Jp * CauchyDIR[k] / rhoMpm - par->_gravity[k] * phi[i])  * dM;
                aResD[k][i] += phi[i] * solAp[k] * dM + CauchyDIR[k] * dA - phi[i] * par->_gravity[k] * dM;
                //if(nodeFlag[i] == 0) { //bulk solid nodes: kinematic: v - dD/dt = 0
                if(eFlag == 2) {
                  aResV[k][i] += -1.0e10 * phi[i] * (solVFpNew[k] - solVSpNew[k]) * areaOld; //TODO
                }
              }
            }
          }

          //BEGIN Nietsche coupling
          if(fabs(U - 0.5) < 0.1) {
            double h = sqrt(dim) * sqrt((vxOld[0][0] - vxOld[0][1]) * (vxOld[0][0] - vxOld[0][1]) +
                                        (vxOld[1][0] - vxOld[1][1]) * (vxOld[1][0] - vxOld[1][1])) ;

            //double thetaI = 1.; // t_{n + 1}
            //double afI = 0.; // t_{n + 1}
            //double afN = par->_af; // t_{n + alpha_f}

            std::vector < adept::adouble > N(dim);
            for(unsigned k = 0; k < dim; k++) N[k] = -NpOld[kp][k][im[kp]]; // from the fluid to the solid
            if(dim == 2) {
              weight = sqrt(N[0] * N[0] + N[1] * N[1]);
              N[0] /= weight;
              N[1] /= weight;
            }
            else {
              weight = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
              N[0] /= weight;
              N[1] /= weight;
              N[2] /= weight;
            }

            std::vector < adept::adouble > tau(dim, 0.);

            for(unsigned k = 0; k < dim; k++) {
              tau[k] += solPp * N[k];
              for(unsigned j = 0; j < dim; j++) {
                tau[k] += -muFluid * (gradSolVp[k][j] + gradSolVp[j][k]) * N[j];
              }
            }

            double c = 0.;
            for(unsigned k = 0; k < dim; k++) {
              c += (solVFpNew[k].value() - solVSpNew[k].value()) * (solVFpNew[k].value() - solVSpNew[k].value());
            }
            c = sqrt(c);

            double thetaM = par->_GAMMA * muFluid / h; //  [rho L/ T]
            double thetaL = par->_GAMMA * rhoFluid * ((c / 6.) + h / (12. * par->_theta * dt)) + thetaM;  // [rho L/T]

            // *** phi_i loop ***
            for(unsigned i = 0; i < nDofs; i++) {
              for(unsigned k = 0; k < dim; k++) {
                aResV[k][i] += (tau[k] /*- !par->_weakP * solPp * N[k]*/) * phi[i] * weight;  // correct sign due to the normal
                aResV[k][i] += thetaM * (solVFpNew[k] - solVSpNew[k]) * phi[i] * weight;

                aResD[k][i] += -tau[k] * phi[i] * weight; // correct sign due to the normal
                aResD[k][i] +=  thetaM * (solVFpNew[k] - solVSpNew[k]) * (-phi[i]) * weight;

                for(unsigned j = 0; j < dim; j++) {
                  aResV[k][i] += - muFluid * gradPhi[i * dim + j] * N[j] * (solVFpNew[k] - solVSpNew[k]) * weight;
                  aResV[k][i] += - muFluid * gradPhi[i * dim + j] * N[k] * (solVFpNew[j] - solVSpNew[j]) * weight;

                  aResV[k][i] += thetaL * (solVFpNew[j] - solVSpNew[j]) * N[j] * phi[i] * N[k] * weight;
                  aResD[k][i] += thetaL * (solVFpNew[j] - solVSpNew[j]) * N[j] * (-phi[i]) * N[k] * weight;

                }
              }
            } // end phi_i loop

            for(unsigned i = 0; i < nDofsP; i++) {
              for(unsigned k = 0; k < dim; k++) {
                aResP[i] += - (phiP[i]  * (solVFpNew[k] - solVSpNew[k]) * N[k]) * weight;
              }
            }
          }
          im[kp]++;
        }
      }
    }
    //END PARTICLE

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
      s.independent(&solD[k][0], nDofs);
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
//   PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL, NULL, 0, 0, 1800, 1800, &viewer1);
//   PetscObjectSetName((PetscObject) viewer1, "FSI matrix");
//   PetscViewerPushFormat(viewer1, PETSC_VIEWER_DRAW_LG);
//   MatView((static_cast< PetscMatrix* >(myKK))->mat(), viewer1);

//   double a;
//   std::cin >> a;

// *************************************
  end_time = clock();
  AssemblyTime += (end_time - start_time);
// ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

}











