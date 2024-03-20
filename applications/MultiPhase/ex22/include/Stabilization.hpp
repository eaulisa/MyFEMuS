
void AssembleStabilizationTerms(MultiLevelProblem& ml_prob) {

  //this function works both for fluid and solid ghost penalty, the boolean fluid switches between the two

  clock_t start_time;

  //pointers and references

  LinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<LinearImplicitSystem> ("NS");
  const unsigned  level = my_nnlin_impl_sys.GetLevelToAssemble();
  MultiLevelSolution* mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* mysolution = mlSol->GetSolutionLevel(level);     // pointer to the solution (level) object

  LinearEquationSolver* myLinEqSolver = my_nnlin_impl_sys._LinSolver[level];  // pointer to the equation (level) object

  Mesh* msh = ml_prob._ml_msh->GetLevel(level);     // pointer to the mesh (level) object
  elem* el = msh->el;   // pointer to the elem object in msh (level)

  SparseMatrix* myKK = myLinEqSolver->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* myRES =  myLinEqSolver->_RES;  // pointer to the global residual vector object in pdeSys (level)

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  const unsigned dim = msh->GetDimension();
  const unsigned dim2 = 3 * (dim - 1);

  // data
  unsigned iproc  = msh->processor_id();
  unsigned nprocs  = msh->n_processors();

  //quantities for iel will have index1
  //quantities for jel will have index2

  vector< vector< double > > solVOld(dim);
  vector< vector< adept::adouble > > solV(dim);
  vector< vector< adept::adouble > > aResV(dim);     // local redidual vector

  vector< adept::adouble > solP;
  vector< adept::adouble > aResP;     // local redidual vector

  vector< double > rhs; // local redidual vector
  vector < double > Jac;

  std::vector <unsigned> sysDofsAll;

  double weight;
  std::vector < double > phi;
  std::vector < double> gradPhi;
  std::vector < double> nablaPhi;

  vector <vector < double> > vx(dim);

  std::cout.precision(10);

  //variable-name handling
  const char varname[3][5] = {"U", "V", "W"};

  vector <unsigned> indexSolV(dim);
  vector <unsigned> indexPdeV(dim);
  for(unsigned k = 0; k < dim; k++) {
    indexSolV[k] = mlSol->GetIndex(&varname[k][0]);
    indexPdeV[k] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[k][0]);
  }
  unsigned solTypeV = mlSol->GetSolutionType(&varname[0][0]);

  unsigned indexSolP1 = mlSol->GetIndex("P1");
  unsigned indexPdeP1 = my_nnlin_impl_sys.GetSolPdeIndex("P1");

  unsigned indexSolP2 = mlSol->GetIndex("P2");
  unsigned indexPdeP2 = my_nnlin_impl_sys.GetSolPdeIndex("P2");

  unsigned indexSolC = mlSol->GetIndex("C");

  unsigned solTypeP = mlSol->GetSolutionType("P1");

  start_time = clock();


  //flagmark
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    double C = (*mysolution->_Sol[indexSolC])(iel);
    if(C == 0. || C == 1.) {

      short unsigned ielt = msh->GetElementType(iel);

      unsigned nDofsV = msh->GetElementDofNumber(iel, solTypeV);
      unsigned nDofsP = msh->GetElementDofNumber(iel, solTypeP);// number of solution element dofs
      unsigned nDofsAll = dim * nDofsV + nDofsP;

      unsigned indexSolP = (C == 1) ? indexSolP1 : indexSolP2;
      unsigned indexPdeP = (C == 1) ? indexPdeP1 : indexPdeP2;
      double rho = rho1 * C + rho2 * (1 - C);
      double mu = mu1 * C + mu2 * (1 - C);

      // resize local arrays
      sysDofsAll.resize(nDofsAll);

      for(unsigned  k = 0; k < dim; k++) {
        solV[k].resize(nDofsV);
        solVOld[k].resize(nDofsV);
        vx[k].resize(nDofsV);
        aResV[k].assign(nDofsV, 0.);
      }
      solP.resize(nDofsP);
      aResP.assign(nDofsP, 0.);

      for(unsigned i = 0; i < nDofsV; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, solTypeV);

        for(unsigned  k = 0; k < dim; k++) {
          solV[k][i] = (*mysolution->_Sol[indexSolV[k]])(idof);
          solVOld[k][i] = (*mysolution->_SolOld[indexSolV[k]])(idof);
          sysDofsAll[k * nDofsV + i] = myLinEqSolver->GetSystemDof(indexSolV[k], indexPdeV[k], i, iel);
        }
      }


      for(unsigned i = 0; i < nDofsP; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, solTypeP);

        solP[i] = (*mysolution->_Sol[indexSolP])(idof);
        sysDofsAll[dim * nDofsV + i] = myLinEqSolver->GetSystemDof(indexSolP, indexPdeP, i, iel);

      }


      for(unsigned i = 0; i < nDofsV; i++) {
        unsigned idofX = msh->GetSolutionDof(i, iel, 2);
        for(unsigned  k = 0; k < dim; k++) {
          vx[k][i] = (*msh->_topology->_Sol[k])(idofX);
        }
      }
      s.new_recording();

      //

      double elementArea = 0.;

      // *** Gauss point loop ***
      for(unsigned ig = 0; ig < msh->_finiteElement[ielt][solTypeV]->GetGaussPointNumber(); ig++) {

        msh->_finiteElement[ielt][solTypeV]->Jacobian(vx, ig, weight, phi, gradPhi, nablaPhi);

        std::vector <std::vector <double> > Jac;
        std::vector <std::vector <double> > JacI;
        msh->_finiteElement[ielt][solTypeV]->GetJacobianMatrix(vx, ig, weight, Jac, JacI); //centered at theta

        vector < adept::adouble > solVg(dim, 0.);
        vector < vector < adept::adouble > > gradSolVg(dim, vector<adept::adouble>(dim,0.));
        vector < vector < adept::adouble > > DeltaSolVg(dim, vector<adept::adouble>(dim2,0.));
        
        
        vector < double > solVgOld(dim, 0.);
        for(unsigned i = 0; i < nDofsV; i++) {
          for(unsigned j = 0; j < dim; j++) {
            solVg[j] += phi[i] * solV[j][i]; // new velocity of background grid
            solVgOld[j] += phi[i] * solVOld[j][i]; // velocity in the undeformed reference configuration
            for(unsigned  k = 0; k < dim; k++) {
              gradSolVg[k][j] += gradPhi[i * dim + j] * solV[k][i]; // gradient of the new velocity with respect to the theta domain
            }

          }
          for(unsigned j = 0; j < dim2; j++) {
            for(unsigned  k = 0; k < dim; k++) {
              DeltaSolVg[k][j] += nablaPhi[i * dim2 + j] * solV[k][i]; // laplace of the theta velocity with respect to the theta domain
            }
          }

        }




//         if(solTypeP == 4) { //discontinuous pressure <1,\xi,\eta> bases centered at theta
//           for(unsigned j = 0; j < dim; j++) {
//             gradPhiP[0 * dim + j]  = 0.;
//           }
//           for(unsigned i = 0; i < dim; i++) {
//             for(unsigned j = 0; j < dim; j++) {
//               gradPhiP[(i + 1) * dim + j] = JacI[i][j];
//             }
//           }
//         }
//
//         vector<adept::adouble> gradSolPg(dim, 0.); //centered at theta
//         for(unsigned i = 0; i < nDofsP; i++) {
//           for(unsigned k = 0; k < dim; k++) {
//             gradSolPg[k] += solP[i] * gradPhiP[i * dim + k];
//           }
//         }


        std::vector <std::vector <double> > G(dim); // J^(-T) . J^(-1) //centered at theta
        for(unsigned i = 0; i < dim; i++) {
          G[i].assign(dim, 0.);
          for(unsigned j = 0; j < dim; j++) {
            for(unsigned k = 0; k < dim; k++) {
              G[i][j] += JacI[k][i] * JacI[k][j];
            }
          }
        }

        double CI = 36.;
        adept::adouble denom = pow(2 * rho / dt, 2.);
        for(unsigned i = 0; i < dim; i++) {
          for(unsigned j = 0; j < dim; j++) {
            denom += rho * rho * solVg[i] * G[i][j] * solVg[j] + CI * mu * mu * G[i][j] * G[i][j];
          }
        }
        adept::adouble tauM = 1. / sqrt(denom);

        adept::adouble tauMtrG = 0.;
        for(unsigned k = 0; k < dim; k++) {
          tauMtrG += G[k][k];
        }
        tauMtrG *= tauM;
        adept::adouble tauC = 1. / tauMtrG;

        //end SUPG parameters
        std::vector < adept::adouble > tauMsupgPhi(nDofsV, 0.);
        for(unsigned i = 0; i < nDofsV; i++) {
          for(unsigned j = 0; j < dim; j++) {
            tauMsupgPhi[i] += tauM * rho * solVg[j] * gradPhi[i * dim + j];
          }
        }

        adept::adouble divVg = 0.;
        for(unsigned k = 0; k < dim; k++) {
          divVg += gradSolVg[k][k];
        }

        for(unsigned i = 0; i < nDofsV; i++) {
          
          for(unsigned k = 0; k < dim; k++) {
            adept::adouble diffusion = 0.;
            adept::adouble advection = 0.;

            for(unsigned j = 0; j < dim; j++) {
              advection +=  rho * solVg[j] * gradSolVg[k][j]; 
              
              unsigned kdim;
              if(k == j) kdim = j;
              else if(1 == k + j) kdim = dim;        // xy
              else if(2 == k + j) kdim = dim + 2;    // xz
              else if(3 == k + j) kdim = dim + 1;    // yz
              diffusion += (-mu * (DeltaSolVg[k][j] + DeltaSolVg[j][kdim])); 
              

            }

            double f = - rho * phi[i] * g[k];
            adept::adouble pressureGradient = 0.;
            adept::adouble rM = rho * (solVg[k] - solVgOld[k]) / dt + advection + diffusion + pressureGradient - f;   
            
            adept::adouble supgDiv = tauC * divVg * gradPhi[i * dim + k];
            
            
            aResV[k][i] += ( rM * tauMsupgPhi[i] + supgDiv) * weight;
          }
        }
      }


      //copy the value of the adept::adoube aRes in double Res and store them in RES
      rhs.resize(nDofsAll);   //resize
      for(unsigned  k = 0; k < dim; k++) {
        for(int i = 0; i < nDofsV; i++) {
          rhs[ k * nDofsV + i ] = -aResV[k][i].value();
        }
      }
      for(int i = 0; i < nDofsP; i++) {
        rhs[ dim * nDofsV + i ] = -aResP[i].value();
      }
      myRES->add_vector_blocked(rhs, sysDofsAll);


      // define the dependent variables J11 and J12
      for(unsigned  k = 0; k < dim; k++) {
        s.dependent(&aResV[k][0], nDofsV);


      }
      s.dependent(&aResP[0], nDofsP);


      // define the independent variables J11
      for(unsigned  k = 0; k < dim; k++) {
        s.independent(&solV[k][0], nDofsV);
      }

      s.independent(&solP[0], nDofsP);

      Jac.resize(nDofsAll * nDofsAll);
      // get the and store jacobian matrix (row-major)
      s.jacobian(&Jac[0], true);
      myKK->add_matrix_blocked(Jac, sysDofsAll, sysDofsAll);

      s.clear_independents();

      s.clear_dependents(); // for J21 and J22
    }
  }

  // *************************************
  std::cout << "Stabilization Assembly time = " << static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC << std::endl;

}
