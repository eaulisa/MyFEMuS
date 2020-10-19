#include "MultiLevelSolution.hpp"
#include "PetscMatrix.hpp"

using namespace femus;




void GetParticlesToNodeFlag(MultiLevelSolution &mlSol, Line & solidLine, Line & fluidLine);
void GetPressureNeighbor(MultiLevelSolution &mlSol, Line & solidLine, Line & fluidLine);

void ProjectGridVelocity(MultiLevelSolution &mlSol);

double GetSmoothStepFunction(const double &dg1, const double &eps) {

  double a0 = 0.5; // 128./256.;
  double a1 = pow(eps, -1.) * 1.23046875; // 315/256.;
  double a3 = -pow(eps, -3.) * 1.640625; //420./256.;
  double a5 = pow(eps, -5.) * 1.4765625; // 378./256.;
  double a7 = -pow(eps, -7.) * 0.703125; // 180./256.;
  double a9 = pow(eps, -9.) * 0.13671875; // 35./256.;

  double dg2 = dg1 * dg1;
  double xi;
  if(dg1 < -eps)
    xi = 0.;
  else if(dg1 > eps) {
    xi = 1.;
  }
  else {
    xi = (a0 + dg1 * (a1 + dg2 * (a3 + dg2 * (a5 + dg2 * (a7 + dg2 * a9)))));
  }

  return xi;
}


double clamp(double x, double lowerlimit, double upperlimit) {
  if(x < lowerlimit)
    x = lowerlimit;
  if(x > upperlimit)
    x = upperlimit;
  return x;
}

double smoothstep(double edge0, double edge1, double x) {
  // Scale, bias and saturate x to 0..1 range
  x = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
  // Evaluate polynomial
  return x * x * (3 - 2 * x);
}

void AssembleSolid(MultiLevelProblem& ml_prob) {

  clock_t AssemblyTime = 0;
  clock_t start_time, end_time;

  TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("MPM_FSI");
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

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  // data
  unsigned iproc  = msh->processor_id();

  vector< vector< adept::adouble > > solD(dim);      // local solution (displacement)
  vector< adept::adouble > solP;

  vector< vector< double > > solDOld(dim);      // local solution (displacement)
  vector< vector< double > > solVOld(dim);
  vector< vector< double > > solAOld(dim);

  vector< double > rhs;    // local redidual vector
  vector< vector< adept::adouble > > aRhsD(dim);     // local redidual vector
  vector< adept::adouble > aRhsP;    // local redidual vector
  vector < double > Jac;

  std::vector <unsigned> sysDofsAll;

  vector < double > phiD;
  vector < double > phiP;

  vector < double > gradPhiP;
  vector < double > gradPhiV;
  vector < adept::adouble > gradPhiD;
  vector < double > gradPhiDHat;

  vector <vector < adept::adouble> > vx(dim);   //vx is coordX in assembly of ex30
  vector <vector < double> > vxHat(dim);

  double weightDHat;
  adept::adouble weightD;


  //reading parameters for MPM body
  double rhoMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_density();
  double EMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_young_module();
  double muMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_lame_shear_modulus();
  double nuMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_poisson_coeff();
  double lambdaMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_lame_lambda();


  double dt =  my_nnlin_impl_sys.GetIntervalTime();

  std::cout.precision(10);

  //variable-name handling
  const char varname[9][5] = {"UX", "UY", "UZ", "VX", "VY", "VZ", "AX", "AY", "AZ"};

  vector <unsigned> indexSolD(dim);
  vector <unsigned> indexSolAOld(dim);
  vector <unsigned> indexSolVOld(dim);

  vector <unsigned> indexPdeD(dim);

  for(unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolD[ivar] = mlSol->GetIndex(&varname[ivar][0]);
    indexSolVOld[ivar] = mlSol->GetIndex(&varname[ivar + 3][0]); // For Newmark in Solid
    indexSolAOld[ivar] = mlSol->GetIndex(&varname[ivar + 6][0]); // For Newmark in Solid
    indexPdeD[ivar] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar][0]);     // DX, DY, DZ

  }
  unsigned solType = mlSol->GetSolutionType(&varname[0][0]);

  unsigned indexSolP = mlSol->GetIndex("P");
  unsigned indexPdeP = my_nnlin_impl_sys.GetSolPdeIndex("P");
  unsigned solTypeP = mlSol->GetSolutionType("P");

  unsigned eflagIndex = mlSol->GetIndex("eflag");
  unsigned nflagIndex = mlSol->GetIndex("nflag");
  std::vector < unsigned >  nodeFlag; // local solution

  start_time = clock();

  //BEGIN loop on elements (to initialize the "soft" stiffness matrix)
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielt = msh->GetElementType(iel);

    unsigned nDofs = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solTypeP);    // number of solution element dofs
    unsigned nDofsAll = dim * nDofs + nDofsP;

    // resize local arrays
    sysDofsAll.resize(nDofsAll);
    nodeFlag.resize(nDofs);

    unsigned eFlag = static_cast <unsigned>(floor((*mysolution->_Sol[eflagIndex])(iel) + 0.5));

    if(eFlag == 4) { //solid

      for(unsigned  k = 0; k < dim; k++) {
        solD[k].resize(nDofs);
        solDOld[k].resize(nDofs);
        solVOld[k].resize(nDofs);
        solAOld[k].resize(nDofs);
        aRhsD[k].assign(nDofs, 0.);
        vx[k].resize(nDofs);
        vxHat[k].resize(nDofs);
      }
      solP.resize(nDofsP);
      aRhsP.assign(nDofsP, 0.);

      for(unsigned i = 0; i < nDofs; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, solType);

        nodeFlag[i] = (*mysolution->_Sol[nflagIndex])(idof); // set it to 0 for no-marker

        for(unsigned  k = 0; k < dim; k++) {
          solD[k][i] = (*mysolution->_Sol[indexSolD[k]])(idof);
          solDOld[k][i] = (*mysolution->_SolOld[indexSolD[k]])(idof);
          solVOld[k][i] = (*mysolution->_Sol[indexSolVOld[k]])(idof);
          solAOld[k][i] = (*mysolution->_Sol[indexSolAOld[k]])(idof);
          sysDofsAll[k * nDofs + i] = myLinEqSolver->GetSystemDof(indexSolD[k], indexPdeD[k], i, iel);
        }
      }

      for(unsigned i = 0; i < nDofsP; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, solTypeP);
        solP[i] = (*mysolution->_Sol[indexSolP])(idof);
        sysDofsAll[dim * nDofs + i] = myLinEqSolver->GetSystemDof(indexSolP, indexPdeP, i, iel);
      }

      // start a new recording of all the operations involving adept::adouble variables
      s.new_recording();

      for(unsigned i = 0; i < nDofs; i++) {
        unsigned idofX = msh->GetSolutionDof(i, iel, 2);
        for(unsigned  k = 0; k < dim; k++) {
          vxHat[k][i] = (*msh->_topology->_Sol[k])(idofX);
          vx[k][i] = (*msh->_topology->_Sol[k])(idofX) + (1. - af) * solDOld[k][i] + af * solD[k][i];
        }
      }



      // *** Gauss point loop ***
      for(unsigned ig = 0; ig < msh->_finiteElement[ielt][solType]->GetGaussPointNumber(); ig++) {

        msh->_finiteElement[ielt][solType]->Jacobian(vxHat, ig, weightDHat, phiD, gradPhiDHat);
        msh->_finiteElement[ielt][solType]->Jacobian(vx, ig, weightD, phiD, gradPhiD);


        std::vector < adept::adouble > solDg(dim, 0.);
        std::vector < adept::adouble > solDgOld(dim, 0.);
        std::vector < double > solVgOld(dim, 0.);
        std::vector < double > solAgOld(dim, 0.);

        std::vector < std::vector < adept::adouble > > gradSolDgHat(dim);

        for(unsigned  k = 0; k < dim; k++) {
          gradSolDgHat[k].assign(dim, 0.);
        }




        for(unsigned i = 0; i < nDofs; i++) {
          for(unsigned j = 0; j < dim; j++) {
            solDg[j] += phiD[i] * solD[j][i];
            solDgOld[j] += phiD[i] * solDOld[j][i];
            solVgOld[j] += phiD[i] * solVOld[j][i];
            solAgOld[j] += phiD[i] * solAOld[j][i];

            for(unsigned  k = 0; k < dim; k++) {
              gradSolDgHat[k][j] += ((1. - af) * solD[k][i] + af * solDOld[k][i]) * gradPhiDHat[i * dim + j];
            }
          }
        }



        std::vector < adept::adouble > solAgAm(dim);
        for(unsigned  k = 0; k < dim; k++) {
          adept::adouble solAgk = (solDg[k] - solDgOld[k]) / (beta * dt * dt) - solVgOld[k] / (beta * dt) - solAgOld[k] * (1. - 2.* beta) / (2. * beta); //NEWMARK ACCELERATION
          solAgAm[k] = ((1. - am) * solAgk + am * solAgOld[k]);
        }



        double *phiP = msh->_finiteElement[ielt][solTypeP]->GetPhi(ig);
        adept::adouble solPg = 0.;




        for(unsigned i = 0; i < nDofsP; i++) {
          solPg += phiP[i] * solP[i];
        }



        adept::adouble F[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned k = 0; k < dim; k++) {
            F[j][k] += gradSolDgHat[j][k];
          }
        }

        adept::adouble J_hat =  F[0][0] * F[1][1] * F[2][2] + F[0][1] * F[1][2] * F[2][0] + F[0][2] * F[1][0] * F[2][1]
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
        adept::adouble Cauchy[3][3];


        for(unsigned j = 0; j < 3; j++) {
          for(unsigned k = 0; k < 3; k++) {
            Cauchy[j][k] = lambdaMpm * log(J_hat) / J_hat * Id2th[j][k] + muMpm / J_hat * (B[j][k] - Id2th[j][k]);    // alternative formulation
          }
        }

        //END computation of the Cauchy Stress
        for(unsigned i = 0; i < nDofs; i++) {
          adept::adouble CauchyDIR[3] = {0., 0., 0.};
          for(unsigned j = 0.; j < dim; j++) {
            for(unsigned k = 0.; k < dim; k++) {
              CauchyDIR[j] += gradPhiD[i * dim + k] * Cauchy[j][k];
            }
          }
          for(unsigned k = 0; k < dim; k++) {
            aRhsD[k][i] -= (rhoMpm * phiD[i] * solAgAm[k] + CauchyDIR[k]) * weightD;
            if(k == 0) {
              aRhsD[k][i] -= (rhoMpm * phiD[i]) * weightD;
            }
          }
        }
        //continuity block
        for(unsigned i = 0; i < nDofsP; i++) {
          aRhsP[i] -= (phiP[i] * solPg) * weightD;
        }
      } // end gauss point loop

      //copy the value of the adept::adoube aRes in double Res and store them in RES
      rhs.resize(nDofsAll);   //resize

      for(int i = 0; i < nDofs; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          rhs[k * nDofs + i] = -aRhsD[k][i].value();
        }
      }
      for(int i = 0; i < nDofsP; i++) {
        rhs[ dim * nDofs  + i] = -aRhsP[i].value();
      }

      myRES->add_vector_blocked(rhs, sysDofsAll);


      Jac.resize(nDofsAll * nDofsAll);
      // define the dependent variables

      for(unsigned  k = 0; k < dim; k++) {
        s.dependent(&aRhsD[k][0], nDofs);
      }
      s.dependent(&aRhsP[0], nDofsP);

      for(unsigned  k = 0; k < dim; k++) {
        s.independent(&solD[k][0], nDofs);
      }
      s.independent(&solP[0], nDofsP);

      // get the and store jacobian matrix (row-major)
      s.jacobian(&Jac[0], true);
      myKK->add_matrix_blocked(Jac, sysDofsAll, sysDofsAll);

      s.clear_independents();
      s.clear_dependents();

    }
  }


// *************************************
  end_time = clock();
  AssemblyTime += (end_time - start_time);
// ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

}


void AssembleSolidInterface(MultiLevelProblem& ml_prob) {

  clock_t AssemblyTime = 0;
  clock_t start_time, end_time;

  TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("MPM_FSI");
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

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  // data
  unsigned iproc  = msh->processor_id();
  unsigned nprocs = msh->n_processors();

  vector< vector< adept::adouble > > solD(dim);      // local solution (displacement)
  vector< vector< adept::adouble > > solV(dim);      // local solution (velocity)
  vector< adept::adouble > solP;


  vector< vector< double > > solV2(dim);      // local solution (velocity)
  vector< vector< double > > solV2Old(dim);
  vector< double > solP2;
  vector < vector< unsigned > > solV2dofs(dim);
  vector< unsigned > solP2dofs;


  vector< vector< double > > solD1(dim);      // local solution (displacement)
  vector< vector< double > > solD1Old(dim);      // local solution (displacement)
  vector< vector< unsigned > > solD1dofs(dim);      // local solution (displacement)

  vector <vector < adept::adouble> > vx1(dim);   //vx1 are goind
  vector <vector < double> > vx1Hat(dim); //1 solid 2 is fluid
  vector <vector < double> > vx2Hat(dim); //1 solid 2 is fluid

  vector< double > rhs;    // local redidual vector
  vector< vector< adept::adouble > > aRhsD(dim);     // local redidual vector
  vector< vector< adept::adouble > > aRhsV(dim);     // local redidual vector
  vector< adept::adouble > aRhsP;    // local redidual vector
  vector < double > Jac;

  std::vector <unsigned> sysDofsAll;

  vector < double > phiV;
  vector < double > phiD;
  vector < double > phiP;

  vector < double > gradPhiP;
  vector < double > gradPhiV;
  vector < adept::adouble > gradPhiD;
  vector < double > gradPhiDHat;

  unsigned dim2 = 3 * (dim - 1);



  double weightP;
  double weightV;
  double weightDHat;
  adept::adouble weightD;


  //reading parameters for MPM body
  double rhoMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_density();
  double EMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_young_module();
  double muMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_lame_shear_modulus();
  double nuMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_poisson_coeff();
  double lambdaMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_lame_lambda();


  //reading parameters for fluid FEM domain
  double rhoFluid = ml_prob.parameters.get<Fluid> ("FluidFEM").get_density();
  double muFluid = ml_prob.parameters.get<Fluid> ("FluidFEM").get_viscosity();

  double dt =  my_nnlin_impl_sys.GetIntervalTime();

  std::cout.precision(10);

  //variable-name handling
  const char varname[12][5] = {"UX", "UY", "UZ", "UX", "UY", "UZ", "VX", "VY", "VZ", "AX", "AY", "AZ"};

  vector <unsigned> indexSolD(dim);
  vector <unsigned> indexSolV(dim);

  vector <unsigned> indexSolAOld(dim);
  vector <unsigned> indexSolVOld(dim);

  vector <unsigned> indexPdeD(dim);
  vector <unsigned> indexPdeV(dim);
  for(unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolD[ivar] = mlSol->GetIndex(&varname[ivar][0]);
    indexSolV[ivar] = mlSol->GetIndex(&varname[ivar + 3][0]);
    indexSolVOld[ivar] = mlSol->GetIndex(&varname[ivar + 6][0]); // For Newmark in Solid
    indexSolAOld[ivar] = mlSol->GetIndex(&varname[ivar + 9][0]); // For Newmark in Solid
    indexPdeD[ivar] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar][0]);     // DX, DY, DZ
    indexPdeV[ivar] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar + 3][0]); //VX, VY, VZ
  }
  unsigned solType = mlSol->GetSolutionType(&varname[0][0]);

  unsigned indexSolP = mlSol->GetIndex("P");
  unsigned indexPdeP = my_nnlin_impl_sys.GetSolPdeIndex("P");
  unsigned solTypeP = mlSol->GetSolutionType("P");

  unsigned eflagIndex = mlSol->GetIndex("eflag");
  unsigned pElemIndex = mlSol->GetIndex("pElem");
  unsigned nflagIndex = mlSol->GetIndex("nflag");

  std::vector <unsigned> nodeFlag1;
  std::vector <unsigned> nodeFlag2;

  start_time = clock();

  unsigned meshOffset = msh->_elementOffset[iproc];
  unsigned meshOffsetp1 = msh->_elementOffset[iproc + 1];

  for(int iel = meshOffset; iel < meshOffsetp1; iel++) {
    unsigned eFlag = static_cast <unsigned>(floor((*mysolution->_Sol[eflagIndex])(iel) + 0.5));
    if(eFlag < 2) { //solid
      mysolution->_Sol[eflagIndex]->set(iel, 2.);
    }
  }
  mysolution->_Sol[eflagIndex]->close();

  //BEGIN loop on elements (to initialize the "soft" stiffness matrix)
  for(unsigned kproc = 0; kproc < nprocs; kproc++) {
    for(int iel1 = msh->_elementOffset[kproc]; iel1 < msh->_elementOffset[kproc + 1]; iel1++) {

      unsigned interface = 0;
      short unsigned ielt1;
      unsigned nDofsD1;
      unsigned pElem;
      if(iproc == kproc) {
        unsigned eFlag = static_cast <unsigned>(floor((*mysolution->_Sol[eflagIndex])(iel1) + 0.5));
        if(eFlag == 4) { //solid
          ielt1 = msh->GetElementType(iel1);
          pElem = static_cast <unsigned>(floor((*mysolution->_Sol[pElemIndex])(iel1) + 0.5));

          for(unsigned i = el->GetFaceRangeStart(ielt1); i < el->GetFaceRangeEnd(ielt1); i++) {
            unsigned idof = msh->GetSolutionDof(i, iel1, solType);
            if((*mysolution->_Sol[nflagIndex])(idof) == 5) {
              interface = 1;

              nDofsD1 = msh->GetElementDofNumber(iel1, solType);
              for(unsigned  k = 0; k < dim; k++) {
                solD1[k].resize(nDofsD1);
                solD1Old[k].resize(nDofsD1);
                vx1Hat[k].resize(nDofsD1);
                aRhsD[k].assign(nDofsD1, 0.);
                solD1dofs[k].resize(nDofsD1);
              }
              nodeFlag1.resize(nDofsD1);
              for(unsigned i = 0; i < nDofsD1; i++) {
                unsigned idofD = msh->GetSolutionDof(i, iel1, solType); //global dof for solution D
                unsigned idofX = msh->GetSolutionDof(i, iel1, 2); //global dof for mesh coordinates

                nodeFlag1[i] = (*mysolution->_Sol[nflagIndex])(idofD); // set it to 0 for no-marker

                for(unsigned  k = 0; k < dim; k++) {
                  solD1[k][i] = (*mysolution->_Sol[indexSolD[k]])(idofD);
                  solD1Old[k][i] = (*mysolution->_SolOld[indexSolD[k]])(idofD);
                  solD1dofs[k][i] = myLinEqSolver->GetSystemDof(indexSolD[k], indexPdeD[k], i, iel1);
                  vx1Hat[k][i] = (*msh->_topology->_Sol[k])(idofX);
                }
              }
              break;
            }
          }
        }
      }
      MPI_Bcast(&interface, 1, MPI_UNSIGNED, kproc, PETSC_COMM_WORLD);

      if(interface) {
        MPI_Bcast(&ielt1, 1, MPI_UNSIGNED, kproc, PETSC_COMM_WORLD);
        MPI_Bcast(&nDofsD1, 1, MPI_UNSIGNED, kproc, PETSC_COMM_WORLD);
        MPI_Bcast(&pElem, 1, MPI_UNSIGNED, kproc, PETSC_COMM_WORLD);

        if(iproc != kproc) {
          for(unsigned  k = 0; k < dim; k++) {
            solD1[k].resize(nDofsD1);
            solD1Old[k].resize(nDofsD1);
            vx1Hat[k].resize(nDofsD1);
            solD1dofs[k].resize(nDofsD1);
          }
          nodeFlag1.resize(nDofsD1);
        }

        for(unsigned  k = 0; k < dim; k++) {
          MPI_Bcast(solD1[k].data(), solD1[k].size(), MPI_DOUBLE, kproc, PETSC_COMM_WORLD);
          MPI_Bcast(solD1Old[k].data(), solD1Old[k].size(), MPI_DOUBLE, kproc, PETSC_COMM_WORLD);
          MPI_Bcast(vx1Hat[k].data(), vx1Hat[k].size(), MPI_DOUBLE, kproc, PETSC_COMM_WORLD);
          MPI_Bcast(solD1dofs[k].data(), solD1dofs[k].size(), MPI_UNSIGNED, kproc, PETSC_COMM_WORLD);
        }
        MPI_Bcast(nodeFlag1.data(), nodeFlag1.size(), MPI_UNSIGNED, kproc, PETSC_COMM_WORLD);

        for(unsigned j1 = el->GetFaceRangeStart(ielt1); j1 < el->GetFaceRangeEnd(ielt1); j1++) {
          if(nodeFlag1[j1] == 5) {

            unsigned jface1 = j1 - el->GetFaceRangeStart(ielt1); //local face number

            unsigned jfaceType1 = el->GetFaceType(ielt1, jface1);
            unsigned jfaceDofs1 = el->GetNFACENODES(ielt1, jface1, solType);

            std::vector  < std::vector  <  double > > xf1(dim);    // physical coordinates of the face in the af configuration
            for(int k = 0; k < dim; k++) {
              xf1[k].resize(jfaceDofs1);
            }

            for(unsigned i = 0; i < jfaceDofs1; i++) {
              unsigned inode = el->GetIG(ielt1, jface1, i); // local mapping from face to element
              for(unsigned k = 0; k < dim; k++) {
                xf1[k][i] =  vx1Hat[k][inode] + af * solD1[k][inode] + (1 - af) * solD1Old[k][inode];
              }
            }

            unsigned ng1 = msh->_finiteElement[jfaceType1][solType]->GetGaussPointNumber();
            std::vector < Marker *> gp(ng1);
            std::vector < unsigned> jel(ng1);
            std::map <unsigned, unsigned> jelCounter;
            std::map <unsigned, unsigned>::iterator it;

            for(unsigned ig = 0; ig < ng1; ig++) {

              std::vector < double> normal;
              const double* phiD1 = msh->_finiteElement[jfaceType1][solType]->GetPhi(ig);

              std::vector< double > xg(dim, 0.);
              for(unsigned i = 0; i < jfaceDofs1; i++) {
                for(unsigned k = 0; k < dim; k++) {
                  xg[k] += phiD1[i] * xf1[k][i];
                }
              }
              gp[ig] = new Marker(xg, VOLUME, mysolution, solType, pElem);
              jel[ig] = gp[ig]->GetMarkerElement();

              it = jelCounter.find(jel[ig]);
              if(it != jelCounter.end()) {
                jelCounter[jel[ig]] += 1;
              }
              else {
                jelCounter[jel[ig]] = 1;
              }

            }
            std::vector < unsigned > IELT2(jelCounter.size());
            std::vector < std::vector < unsigned > > IG(jelCounter.size());

            for(it == jelCounter.begin(); it != jelCounter.end(); it++) {
              std::cout << it->first <<" "<< it ->second << std::endl;
            }
           

//             {
//               unsigned i = 0;
//               for(i = 0, it == jelCounter.begin(); it != jelCounter.end(); it++, i++) {
//                 IELT2[i] = it->first;
//                 IG[i].reserve(it->second);
//               }
// 
//               std::cout << IELT2[i] << " ";
//             }
// 
//             std::cout << std::endl;

//             for(unsigned ig = 0; ig < gp.size(); ig++) {
//               unsigned i = 0;
//               while(jel[ig] != IELT2[i]) i++;
//               unsigned k = IG[i].size();
//               IG[i].resize(k + 1);
//               IG[i][k] = ig;
//             }
//
//             for(unsigned i = 0; i < IELT2.size(); i++){
//               std::cout << "On element " << IELT2[i] << " we have gauss points\n ";
//               for(unsigned k = 0; k < IG[i].size(); k++){
//                 std::cout<< IG[i][k] << " ";
//               }
//               std::cout << std::endl;
//             }


            for(unsigned ig = 0; ig < gp.size(); ig++) {

              unsigned iel2 = gp[ig]->GetMarkerElement();



              if(iel2 != UINT_MAX) {
                unsigned jproc = msh->IsdomBisectionSearch(iel2 , 3); // return  jproc for piece-wise constant discontinuous type (3)


                unsigned ielt2;
                unsigned nDofsV2;
                unsigned nDofsP2;

                std::vector < double > xi =  gp[ig]->GetMarkerLocalCoordinates();

                if(iproc == jproc) {
                  mysolution->_Sol[eflagIndex]->set(iel2, 1.);

                  ielt2 = msh->GetElementType(iel2);
                  nDofsV2 = msh->GetElementDofNumber(iel2, solType);
                  nDofsP2 = msh->GetElementDofNumber(iel2, solTypeP);
                  for(unsigned  k = 0; k < dim; k++) {
                    solV2[k].resize(nDofsV2);
                    solV2Old[k].resize(nDofsV2);
                    vx2Hat[k].resize(nDofsV2);
                    aRhsV[k].assign(nDofsV2, 0.);
                    solV2dofs[k].resize(nDofsV2);
                  }
                  solP2.resize(nDofsP2);
                  aRhsP.assign(nDofsP2, 0.);
                  solP2dofs.resize(nDofsP2);

                  for(unsigned i = 0; i < nDofsV2; i++) {
                    unsigned idofV = msh->GetSolutionDof(i, iel2, solType); //global dof for solution D
                    unsigned idofX = msh->GetSolutionDof(i, iel2, 2); //global dof for mesh coordinates
                    for(unsigned  k = 0; k < dim; k++) {
                      solV2[k][i] = (*mysolution->_Sol[indexSolV[k]])(idofV);
                      solV2Old[k][i] = (*mysolution->_SolOld[indexSolV[k]])(idofV);
                      solV2dofs[k][i] = myLinEqSolver->GetSystemDof(indexSolV[k], indexPdeV[k], i, iel2);
                      vx2Hat[k][i] = (*msh->_topology->_Sol[k])(idofX);
                    }
                  }
                  for(unsigned i = 0; i < nDofsP2; i++) {
                    unsigned idofP = msh->GetSolutionDof(i, iel2, solTypeP); //global dof for solution D
                    solP2[i] = (*mysolution->_Sol[indexSolP])(idofP);
                    solP2dofs[i] = myLinEqSolver->GetSystemDof(indexSolP, indexPdeP, i, iel2);
                  }

                  MPI_Send(&ielt2, 1, MPI_UNSIGNED, kproc, 0, PETSC_COMM_WORLD);
                  MPI_Send(&nDofsV2, 1, MPI_UNSIGNED, kproc, 1, PETSC_COMM_WORLD);
                  MPI_Send(&nDofsP2, 1, MPI_UNSIGNED, kproc, 2, PETSC_COMM_WORLD);

                  for(unsigned  k = 0; k < dim; k++) {
                    MPI_Send(solV2[k].data(), solV2[k].size(), MPI_DOUBLE, kproc, 3 + (k * 4), PETSC_COMM_WORLD);
                    MPI_Send(solV2Old[k].data(), solV2Old[k].size(), MPI_DOUBLE, kproc, 4 + (k * 4), PETSC_COMM_WORLD);
                    MPI_Send(vx2Hat[k].data(), vx2Hat[k].size(), MPI_DOUBLE, kproc, 5 + (k * 4), PETSC_COMM_WORLD);
                    MPI_Send(solV2dofs[k].data(), solV2dofs[k].size(), MPI_UNSIGNED, kproc, 6 + (k * 4), PETSC_COMM_WORLD);
                  }
                  MPI_Send(solP2.data(), solP2.size(), MPI_DOUBLE, kproc, 3 + (dim * 4), PETSC_COMM_WORLD);
                  MPI_Send(solP2dofs.data(), solP2dofs.size(), MPI_UNSIGNED, kproc, 4 + (dim * 4), PETSC_COMM_WORLD);

                  MPI_Send(xi.data(), xi.size(), MPI_DOUBLE, kproc, 3 + (dim * 4), PETSC_COMM_WORLD);

                }

                if(iproc == kproc) {
                  MPI_Recv(&ielt2, 1, MPI_UNSIGNED, jproc, 0, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                  MPI_Recv(&nDofsV2, 1, MPI_UNSIGNED, jproc, 1, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                  MPI_Recv(&nDofsP2, 1, MPI_UNSIGNED, jproc, 2, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);

                  for(unsigned  k = 0; k < dim; k++) {
                    solV2[k].resize(nDofsV2);
                    solV2Old[k].resize(nDofsV2);
                    vx2Hat[k].resize(nDofsV2);
                    solV2dofs[k].resize(nDofsV2);
                  }
                  solP2.resize(nDofsP2);
                  solP2dofs.resize(nDofsP2);

                  for(unsigned  k = 0; k < dim; k++) {
                    MPI_Recv(solV2[k].data(), solV2[k].size(), MPI_DOUBLE, jproc, 3 + (k * 4), PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(solV2Old[k].data(), solV2Old[k].size(), MPI_DOUBLE, jproc, 4 + (k * 4), PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(vx2Hat[k].data(), vx2Hat[k].size(), MPI_DOUBLE, jproc, 5 + (k * 4), PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(solV2dofs[k].data(), solV2dofs[k].size(), MPI_UNSIGNED, jproc, 6 + (k * 4), PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                  }
                  MPI_Recv(solP2.data(), solP2.size(), MPI_DOUBLE, jproc, 3 + (dim * 4), PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                  MPI_Recv(solP2dofs.data(), solP2dofs.size(), MPI_UNSIGNED, jproc, 4 + (dim * 4), PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                  MPI_Recv(xi.data(), xi.size(), MPI_DOUBLE, jproc, 3 + (dim * 4), PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                }



              }


              pElem = (iel2 != UINT_MAX) ? iel2 : gp[ig]->GetIprocMarkerPreviousElement();
              if(iproc == kproc) {
                mysolution->_Sol[pElemIndex]->set(iel1, pElem);
              }

            }

            for(unsigned ig = 0; ig < gp.size(); ig++) {
              delete gp[ig];
            }
          }
        }
      }
    }

  }

  mysolution->_Sol[eflagIndex]->close();
  mysolution->_Sol[pElemIndex]->close();



  /*


  // resize local arrays
  sysDofsAll.resize(nDofsAll);
  nodeFlag.resize(nDofs);





  for(unsigned  k = 0; k < dim; k++) {
  solD[k].resize(nDofs);
  solDOld[k].resize(nDofs);
  solVOld[k].resize(nDofs);
  solAOld[k].resize(nDofs);
  aRhsD[k].assign(nDofs, 0.);
  vx[k].resize(nDofs);
  vxHat[k].resize(nDofs);
  }
  solP.resize(nDofsP);
  aRhsP.assign(nDofsP, 0.);

  for(unsigned i = 0; i < nDofs; i++) {
  unsigned idof = msh->GetSolutionDof(i, iel, solType);

  nodeFlag[i] = (*mysolution->_Sol[nflagIndex])(idof); // set it to 0 for no-marker

  for(unsigned  k = 0; k < dim; k++) {
  solD[k][i] = (*mysolution->_Sol[indexSolD[k]])(idof);
  solDOld[k][i] = (*mysolution->_SolOld[indexSolD[k]])(idof);
  solVOld[k][i] = (*mysolution->_Sol[indexSolVOld[k]])(idof);
  solAOld[k][i] = (*mysolution->_Sol[indexSolAOld[k]])(idof);
  sysDofsAll[k * nDofs + i] = myLinEqSolver->GetSystemDof(indexSolD[k], indexPdeD[k], i, iel);
  }
  }

  for(unsigned i = 0; i < nDofsP; i++) {
  unsigned idof = msh->GetSolutionDof(i, iel, solTypeP);
  solP[i] = (*mysolution->_Sol[indexSolP])(idof);
  sysDofsAll[dim * nDofs + i] = myLinEqSolver->GetSystemDof(indexSolP, indexPdeP, i, iel);
  }

  // start a new recording of all the operations involving adept::adouble variables
  s.new_recording();

  for(unsigned i = 0; i < nDofs; i++) {
  unsigned idofX = msh->GetSolutionDof(i, iel, 2);
  for(unsigned  k = 0; k < dim; k++) {
  vxHat[k][i] = (*msh->_topology->_Sol[k])(idofX);
  vx[k][i] = (*msh->_topology->_Sol[k])(idofX) + (1. - af) * solDOld[k][i] + af * solD[k][i];
  }
  }



  // *** Gauss point loop ***
  for(unsigned ig = 0; ig < msh->_finiteElement[ielt][solType]->GetGaussPointNumber(); ig++) {

  msh->_finiteElement[ielt][solType]->Jacobian(vxHat, ig, weightDHat, phiD, gradPhiDHat);
  msh->_finiteElement[ielt][solType]->Jacobian(vx, ig, weightD, phiD, gradPhiD);


  std::vector < adept::adouble > solDg(dim, 0.);
  std::vector < adept::adouble > solDgOld(dim, 0.);
  std::vector < double > solVgOld(dim, 0.);
  std::vector < double > solAgOld(dim, 0.);

  std::vector < std::vector < adept::adouble > > gradSolDgHat(dim);

  for(unsigned  k = 0; k < dim; k++) {
  gradSolDgHat[k].assign(dim, 0.);
  }




  for(unsigned i = 0; i < nDofs; i++) {
  for(unsigned j = 0; j < dim; j++) {
    solDg[j] += phiD[i] * solD[j][i];
    solDgOld[j] += phiD[i] * solDOld[j][i];
    solVgOld[j] += phiD[i] * solVOld[j][i];
    solAgOld[j] += phiD[i] * solAOld[j][i];

    for(unsigned  k = 0; k < dim; k++) {
      gradSolDgHat[k][j] += ((1. - af) * solD[k][i] + af * solDOld[k][i]) * gradPhiDHat[i * dim + j];
    }
  }
  }



  std::vector < adept::adouble > solAgAm(dim);
  for(unsigned  k = 0; k < dim; k++) {
  adept::adouble solAgk = (solDg[k] - solDgOld[k]) / (beta * dt * dt) - solVgOld[k] / (beta * dt) - solAgOld[k] * (1. - 2.* beta) / (2. * beta); //NEWMARK ACCELERATION
  solAgAm[k] = ((1. - am) * solAgk + am * solAgOld[k]);
  }



  double *phiP = msh->_finiteElement[ielt][solTypeP]->GetPhi(ig);
  adept::adouble solPg = 0.;




  for(unsigned i = 0; i < nDofsP; i++) {
  solPg += phiP[i] * solP[i];
  }



  adept::adouble F[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
  for(unsigned j = 0; j < dim; j++) {
  for(unsigned k = 0; k < dim; k++) {
    F[j][k] += gradSolDgHat[j][k];
  }
  }

  adept::adouble J_hat =  F[0][0] * F[1][1] * F[2][2] + F[0][1] * F[1][2] * F[2][0] + F[0][2] * F[1][0] * F[2][1]
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
  adept::adouble Cauchy[3][3];


  for(unsigned j = 0; j < 3; j++) {
  for(unsigned k = 0; k < 3; k++) {
    Cauchy[j][k] = lambdaMpm * log(J_hat) / J_hat * Id2th[j][k] + muMpm / J_hat * (B[j][k] - Id2th[j][k]);    // alternative formulation
  }
  }

  //END computation of the Cauchy Stress
  for(unsigned i = 0; i < nDofs; i++) {
  adept::adouble CauchyDIR[3] = {0., 0., 0.};
  for(unsigned j = 0.; j < dim; j++) {
    for(unsigned k = 0.; k < dim; k++) {
      CauchyDIR[j] += gradPhiD[i * dim + k] * Cauchy[j][k];
    }
  }
  for(unsigned k = 0; k < dim; k++) {
    aRhsD[k][i] -= (rhoMpm * phiD[i] * solAgAm[k] + CauchyDIR[k]) * weightD;
    if(k == 0) {
      aRhsD[k][i] -= (rhoMpm * phiD[i]) * weightD;
    }
  }
  }
  //continuity block
  for(unsigned i = 0; i < nDofsP; i++) {
  aRhsP[i] -= (phiP[i] * solPg) * weightD;
  }
  } // end gauss point loop

  //copy the value of the adept::adoube aRes in double Res and store them in RES
  rhs.resize(nDofsAll);   //resize

  for(int i = 0; i < nDofs; i++) {
  for(unsigned  k = 0; k < dim; k++) {
  rhs[k * nDofs + i] = -aRhsD[k][i].value();
  }
  }
  for(int i = 0; i < nDofsP; i++) {
  rhs[ dim * nDofs  + i] = -aRhsP[i].value();
  }

  myRES->add_vector_blocked(rhs, sysDofsAll);


  Jac.resize(nDofsAll * nDofsAll);
  // define the dependent variables

  for(unsigned  k = 0; k < dim; k++) {
  s.dependent(&aRhsD[k][0], nDofs);
  }
  s.dependent(&aRhsP[0], nDofsP);

  for(unsigned  k = 0; k < dim; k++) {
  s.independent(&solD[k][0], nDofs);
  }
  s.independent(&solP[0], nDofsP);

  // get the and store jacobian matrix (row-major)
  s.jacobian(&Jac[0], true);
  myKK->add_matrix_blocked(Jac, sysDofsAll, sysDofsAll);

  s.clear_independents();
  s.clear_dependents();

  }
  }
  }
  */

// *************************************
  end_time = clock();
  AssemblyTime += (end_time - start_time);
// ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

}























void GridToParticlesProjection(MultiLevelProblem & ml_prob,
//                                Line & solidLine, Line & fluidLine, Line & interfaceLine,
                               Line & bulk, Line & lineI) {

  //pointers and references

  TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("MPM_FSI");

  const unsigned  level = my_nnlin_impl_sys.GetLevelToAssemble();
  MultiLevelSolution* mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* mysolution = mlSol->GetSolutionLevel(level);     // pointer to the solution (level) object

  Mesh* msh = ml_prob._ml_msh->GetLevel(level);     // pointer to the mesh (level) object
  elem* el = msh->el;   // pointer to the elem object in msh (level)

  double dt =  my_nnlin_impl_sys.GetIntervalTime();
  const unsigned dim = msh->GetDimension();

  // data
  unsigned iproc  = msh->processor_id();

  // local objects
  vector< vector < double > > solD(dim);
  vector< vector < double > > solDOld(dim);
  vector< vector < double > > gradSolDHat(dim);

  for(int k = 0; k < dim; k++) {
    gradSolDHat[k].resize(dim);
  }

  vector < double > phiHat;
  vector < double > gradPhiHat;

  vector <vector < double> > vxHat(dim);   //vx is coordX in assembly of ex30

  double weightHat;

  //variable-name handling
  const char varname[9][3] = {"DX", "DY", "DZ"};
  vector <unsigned> indexSolD(dim);
  unsigned solType = mlSol->GetSolutionType(&varname[0][0]);

  for(unsigned k = 0; k < dim; k++) {
    indexSolD[k] = mlSol->GetIndex(&varname[k][0]);
  }

  //BEGIN loop on bulk particles

  std::vector<Marker*> particles =  bulk.GetParticles();
  std::vector<unsigned> markerOffset = bulk.GetMarkerOffset();

  unsigned ielOld = UINT_MAX;
  for(unsigned iMarker = markerOffset[iproc]; iMarker < markerOffset[iproc + 1]; iMarker++) {
    unsigned iel = particles[iMarker]->GetMarkerElement();
    if(iel != UINT_MAX) {
      short unsigned ielt;
      unsigned nDofs;
      //update element related quantities only if we are in a different element
      if(iel != ielOld) {


        ielt = msh->GetElementType(iel);
        nDofs = msh->GetElementDofNumber(iel, solType);
        for(int i = 0; i < dim; i++) {
          solD[i].resize(nDofs);
          solDOld[i].resize(nDofs);
          vxHat[i].resize(nDofs);
        }

        for(unsigned inode = 0; inode < nDofs; inode++) {
          unsigned idof = msh->GetSolutionDof(inode, iel, solType);   //local 2 global solution
          unsigned idofX = msh->GetSolutionDof(inode, iel, 2);   //local 2 global solution
          for(int i = 0; i < dim; i++) {
            solDOld[i][inode] = (*mysolution->_SolOld[indexSolD[i]])(idof);
            solD[i][inode] = (*mysolution->_Sol[indexSolD[i]])(idof) - solDOld[i][inode];
            //moving domain
            vxHat[i][inode] = (*msh->_topology->_Sol[i])(idofX) + solDOld[i][inode];
          }
        }

      }
      std::vector <double> xi = particles[iMarker]->GetMarkerLocalCoordinates();

      msh->_finiteElement[ielt][solType]->Jacobian(vxHat, xi, weightHat, phiHat, gradPhiHat);

      std::vector <double> solVpOld(dim);
      particles[iMarker]->GetMarkerVelocity(solVpOld);

      std::vector <double> solApOld(dim);
      particles[iMarker]->GetMarkerAcceleration(solApOld);

      std::vector <double> solDp(dim, 0.);
      //update displacement and acceleration
      for(int i = 0; i < dim; i++) {
        for(unsigned inode = 0; inode < nDofs; inode++) {
          solDp[i] += phiHat[inode] * solD[i][inode];
        }
      }

      particles[iMarker]->SetMarkerDisplacement(solDp);
      particles[iMarker]->UpdateParticleCoordinates();

      std::vector <double> solAp(dim);
      std::vector <double> solVp(dim);
      for(unsigned i = 0; i < dim; i++) {
        solAp[i] = 1. / (beta * dt * dt) * solDp[i] - 1. / (beta * dt) * solVpOld[i] - (1. - 2.* beta) / (2. * beta) * solApOld[i];
        solVp[i] = solVpOld[i] + dt * ((1. - Gamma) * solApOld[i] + Gamma * solAp[i]);
      }

      particles[iMarker]->SetMarkerVelocity(solVp);
      particles[iMarker]->SetMarkerAcceleration(solAp);

      //   update the deformation gradient
      for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
          gradSolDHat[i][j] = 0.;
          for(unsigned inode = 0; inode < nDofs; inode++) {
            gradSolDHat[i][j] +=  gradPhiHat[inode * dim + j] * solD[i][inode];
          }
        }
      }
      std::vector < std::vector < double > > FpOld;
      FpOld = particles[iMarker]->GetDeformationGradient(); //extraction of the deformation gradient
      double FpNew[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
      std::vector < std::vector < double > > Fp(dim);
      for(unsigned i = 0; i < dim; i++) {
        for(unsigned j = 0; j < dim; j++) {
          FpNew[i][j] +=  gradSolDHat[i][j];
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
  //END loop on bulk particles

  //BEGIN loop on interface particles
  particles = lineI.GetParticles();
  markerOffset = lineI.GetMarkerOffset();
  ielOld = UINT_MAX;
  for(unsigned iMarker = markerOffset[iproc]; iMarker < markerOffset[iproc + 1]; iMarker++) {
    unsigned iel = particles[iMarker]->GetMarkerElement();
    if(iel != UINT_MAX) {
      short unsigned ielt;
      unsigned nDofs;
      //update element related quantities only if we are in a different element
      if(iel != ielOld) {
        ielt = msh->GetElementType(iel);
        nDofs = msh->GetElementDofNumber(iel, solType);
        for(int i = 0; i < dim; i++) {
          solD[i].resize(nDofs);
          solDOld[i].resize(nDofs);
          vxHat[i].resize(nDofs);
        }

        for(unsigned inode = 0; inode < nDofs; inode++) {
          unsigned idof = msh->GetSolutionDof(inode, iel, solType);   //local 2 global solution
          unsigned idofX = msh->GetSolutionDof(inode, iel, 2);   //local 2 global solution
          for(int i = 0; i < dim; i++) {
            solDOld[i][inode] = (*mysolution->_SolOld[indexSolD[i]])(idof);
            solD[i][inode] = (*mysolution->_Sol[indexSolD[i]])(idof) - solDOld[i][inode];
            //moving domain
            vxHat[i][inode] = (*msh->_topology->_Sol[i])(idofX) + solDOld[i][inode];
          }
        }

      }
      std::vector <double> xi = particles[iMarker]->GetMarkerLocalCoordinates();

      msh->_finiteElement[ielt][solType]->Jacobian(vxHat, xi, weightHat, phiHat, gradPhiHat);

      std::vector <double> solVpOld(dim);
      particles[iMarker]->GetMarkerVelocity(solVpOld);

      std::vector <double> solApOld(dim);
      particles[iMarker]->GetMarkerAcceleration(solApOld);

      std::vector <double> solDp(dim, 0.);
      //update displacement and acceleration
      for(int i = 0; i < dim; i++) {
        for(unsigned inode = 0; inode < nDofs; inode++) {
          solDp[i] += phiHat[inode] * solD[i][inode];
        }
      }

      particles[iMarker]->SetMarkerDisplacement(solDp);
      particles[iMarker]->UpdateParticleCoordinates();

      std::vector <double> solAp(dim);
      std::vector <double> solVp(dim);
      for(unsigned i = 0; i < dim; i++) {
        solAp[i] = 1. / (beta * dt * dt) * solDp[i] - 1. / (beta * dt) * solVpOld[i] - (1. - 2.* beta) / (2. * beta) * solApOld[i];
        solVp[i] = solVpOld[i] + dt * ((1. - Gamma) * solApOld[i] + Gamma * solAp[i]);
      }

      particles[iMarker]->SetMarkerVelocity(solVp);
      particles[iMarker]->SetMarkerAcceleration(solAp);

      //   update the deformation gradient
      for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
          gradSolDHat[i][j] = 0.;
          for(unsigned inode = 0; inode < nDofs; inode++) {
            gradSolDHat[i][j] +=  gradPhiHat[inode * dim + j] * solD[i][inode];
          }
        }
      }
      std::vector < std::vector < double > > FpOld;
      FpOld = particles[iMarker]->GetDeformationGradient(); //extraction of the deformation gradient
      double FpNew[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
      std::vector < std::vector < double > > Fp(dim);
      for(unsigned i = 0; i < dim; i++) {
        for(unsigned j = 0; j < dim; j++) {
          FpNew[i][j] += gradSolDHat[i][j];
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
  //END loop on interface particle

  std::cout << "Projecting velocity" << std::endl;
  ProjectGridVelocity(*mlSol);

  //BEGIN loop on elements to update grid velocity and acceleration
  for(unsigned idof = msh->_dofOffset[solType][iproc]; idof < msh->_dofOffset[solType][iproc + 1]; idof++) {
    for(int i = 0; i < dim; i++) {
      mysolution->_Sol[indexSolD[i]]->set(idof, 0.);
    }
  }

  for(int i = 0; i < dim; i++) {
    mysolution->_Sol[indexSolD[i]]->close();
  }
  //END loop on elements to update grid velocity and acceleration

  std::cout << "updating interface line" << std::endl;
  lineI.UpdateLineMPM();
  std::cout << "updating bulk" << std::endl;
  bulk.UpdateLineMPM();

  bool updateMat = false;
  std::cout << "GetParticlesToGridMaterial Interfcae" << std::endl;
  lineI.GetParticlesToGridMaterial(updateMat);
  std::cout << "GetParticlesToGridMaterial bulk" << std::endl;
  bulk.GetParticlesToGridMaterial(updateMat);

  std::cout << "Building element flag" << std::endl;
  BuildFlag(*mlSol);

  //GetParticleWeights(*mlSol, &bulk);

}

unsigned getNumberOfLayers(const double & a, const double & fac, const bool inverse = true) {

  double fac1  = (inverse) ? fac : 1. / fac;
  double da = 1. / fac1;
  double b =  da;
  unsigned n = 1;

  while(b < a) {
    da /= fac1;
    b += da;
    n++;
    if(n >= 100) {
      std::cout << "Error: number of layer is unbounded, try with a smaller factor\n";
      abort();
    }
  }
  return n;
}

void ProjectGridVelocity(MultiLevelSolution & mlSol) {

  const unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;
  Mesh* msh = mlSol._mlMesh->GetLevel(level);
  Solution* sol  = mlSol.GetSolutionLevel(level);
  const unsigned dim = msh->GetDimension();

  unsigned iproc  = msh->processor_id();
  unsigned nprocs = msh->n_processors();

  vector< vector< double > > solV(dim);      // local solution (velocity)

  vector< bool > nodeFlag;     // local solution (velocity)
  vector< unsigned > idof;     // local solution (velocity)

  vector < double > phi;
  vector <vector < double> > vx(dim);   //vx is coordX in assembly of ex30
  vector <vector < double> > xp;


  //variable-name handling
  const char varname[10][5] = {"DX", "DY", "DZ", "VX", "VY", "VZ"};


  vector <unsigned> indexSolD(dim);
  vector <unsigned> indexSolV(dim);

  for(unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolD[ivar] = mlSol.GetIndex(&varname[ivar][0]);
    indexSolV[ivar] = mlSol.GetIndex(&varname[ivar + 3][0]);
  }
//   unsigned indexNodeFlag =  mlSol.GetIndex("NodeFlag");
  unsigned indexNodeFlag =  mlSol.GetIndex("nflag");

  unsigned solType = mlSol.GetSolutionType(&varname[0][0]);

  sol->_Sol[indexNodeFlag]->zero();

  for(unsigned k = 0; k < dim; k++) {
    (*sol->_SolOld[indexSolV[k]]) = (*sol->_Sol[indexSolV[k]]);
    sol->_Sol[indexSolV[k]]->zero();
  }

  unsigned counter = 0;

  std::vector < std::vector < std::vector <double > > > aP(3);

  //BEGIN loop on elements
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielType = msh->GetElementType(iel);

    unsigned nDofs = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs

    for(unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofs);
      vx[k].resize(nDofs);
    }
    nodeFlag.resize(nDofs);
    idof.resize(nDofs);
    xp.resize(nDofs);
    for(unsigned i = 0; i < nDofs; i++) {
      xp[i].resize(dim);
    }

    for(unsigned i = 0; i < nDofs; i++) {
      idof[i] = msh->GetSolutionDof(i, iel, solType);
      unsigned idofX = msh->GetSolutionDof(i, iel, 2);

      for(unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_SolOld[indexSolV[k]])(idof[i]);
        xp[i][k] = (*msh->_topology->_Sol[k])(idofX);     // coordinates of the reference configuration;
        vx[k][i] = xp[i][k] + (*sol->_Sol[indexSolD[k]])(idof[i]);     // coordinates of the deformed configuration
      }
      nodeFlag[i] = ((*sol->_Sol[indexNodeFlag])(idof[i]) > 0.5) ? true : false;
    }

    bool aPIsInitialized = false;

    double r;
    std::vector <double> xc;
    GetConvexHullSphere(vx, xc, r, 0.0001);  // get the ball that circumscribe the element
    double r2 = r * r;

    std::vector < std::vector< double > > xe; // get the box that encloses the element
    GetBoundingBox(vx, xe, 0.0001);

    for(unsigned i = 0; i < nDofs; i++) {  // loop on the nodes of the reference elements now considered as independent points
      if(!nodeFlag[i]) {
        double d2 = 0.;
        for(int k = 0; k < dim; k++) {
          d2 += (xp[i][k] - xc[k]) * (xp[i][k] - xc[k]);
        }
        bool insideHull = true;
        if(d2 > r2) {
          insideHull = false;
        }
        for(unsigned k = 0; k < dim; k++) {
          if(xp[i][k] < xe[k][0] || xp[i][k] > xe[k][1]) {
            insideHull = false;
          }
        }

        if(insideHull) {  //rough test
          if(!aPIsInitialized) {
            aPIsInitialized = true;
            std::vector < std::vector <double> > x1(dim);
            for(unsigned jtype = 0; jtype < solType + 1; jtype++) {
              ProjectNodalToPolynomialCoefficients(aP[jtype], vx, ielType, jtype) ;
            }
          }
          std::vector <double> xi;
          GetClosestPointInReferenceElement(vx, xp[i], ielType, xi);

          bool inverseMapping = GetInverseMapping(solType, ielType, aP, xp[i], xi, 100);
          if(!inverseMapping) {
            std::cout << "InverseMapping failed at " << iel << " " << idof[i] << std::endl;
          }


          bool insideDomain = CheckIfPointIsInsideReferenceDomain(xi, ielType, 1.e-3);  // fine testing

          if(inverseMapping && insideDomain) {
            sol->_Sol[indexNodeFlag]->add(idof[i], 1.);
            msh->_finiteElement[ielType][solType]->GetPhi(phi, xi);
            //std::cout << iel << " " << i << "  ";
            counter++;
            for(unsigned k = 0; k < dim; k++) {
              double solVk = 0.;
              for(unsigned j = 0; j < nDofs; j++)    {
                solVk += phi[j] * solV[k][j];
              }
              sol->_Sol[indexSolV[k]]->add(idof[i], solVk);
            }
          }
        }
      }
    }
  }

  sol->_Sol[indexNodeFlag]->close();

  for(unsigned k = 0; k < dim; k++) {
    sol->_Sol[indexSolV[k]]->close();
  }

  unsigned c0 = 0;
  for(unsigned i = msh->_dofOffset[solType][iproc]; i < msh->_dofOffset[solType][iproc + 1]; i++) {
    unsigned cnt = static_cast < unsigned >(floor((*sol->_Sol[indexNodeFlag])(i) + 0.5));
    if(cnt == 0) {
      c0++;
    }
    else if(cnt > 1) {
      counter -= (cnt - 1);
      for(unsigned k = 0; k < dim; k++) {
        double velk = (*sol->_Sol[indexSolV[k]])(i) / cnt;
        sol->_Sol[indexSolV[k]]->set(i, velk);
      }
    }
  }

  unsigned counterAll;
  MPI_Reduce(&counter, &counterAll, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  std::cout << "COUNTER = " << counterAll << " " << msh->GetTotalNumberOfDofs(solType) << std::endl;

  idof.resize(c0);
  vector< double > xp0(c0 * dim);

  unsigned c1 = 0;
  for(unsigned i = msh->_dofOffset[solType][iproc]; i < msh->_dofOffset[solType][iproc + 1]; i++) {
    if(static_cast < unsigned >(floor((*sol->_Sol[indexNodeFlag])(i) + 0.5)) == 0) {
      idof[c1] = i;
      for(unsigned k = 0; k < dim; k++) {
        xp0[c1 * dim + k] = (*msh->_topology->_Sol[k])(i);
      }
      c1++;
      if(c1 == c0) break;
    }
  }

  vector< double > xp1;
  for(unsigned jproc = 0; jproc < nprocs; jproc++) {
    c1 = c0;
    MPI_Bcast(&c1, 1, MPI_UNSIGNED, jproc, PETSC_COMM_WORLD);
    if(c1) {
      xp1.resize(c1 * dim);
      if(iproc == jproc) {
        for(unsigned i = 0; i < c1; i++) {
          for(unsigned k = 0; k < dim; k++) {
            xp1[i * dim + k] = xp0[i * dim + k];
          }
        }
      }
      MPI_Bcast(&xp1[0], c1 * dim, MPI_DOUBLE, jproc, PETSC_COMM_WORLD);

      for(unsigned i = 0; i < c1; i++) {
        std::vector < double > xp(dim);
        for(unsigned k = 0; k < dim; k++) {
          xp[k] = xp1[i * dim + k];
        }
        Marker p(xp, 1, VOLUME, sol, solType, UINT_MAX, 1.);
        unsigned mproc = p.GetMarkerProc(sol);
        if(iproc == mproc) {
          unsigned jel = p.GetMarkerElement();
          short unsigned jelType = msh->GetElementType(jel);
          unsigned nDofs = msh->GetElementDofNumber(jel, solType);    // number of solution element dofs

          for(unsigned  k = 0; k < dim; k++) {
            solV[k].resize(nDofs);
            vx[k].resize(nDofs);
          }

          for(unsigned j = 0; j < nDofs; j++) {
            unsigned jdof = msh->GetSolutionDof(j, jel, solType);
            unsigned jdofX = msh->GetSolutionDof(j, jel, 2);
            for(unsigned  k = 0; k < dim; k++) {
              solV[k][j] = (*sol->_SolOld[indexSolV[k]])(jdof);     //velocity to be projected
              vx[k][j] = (*msh->_topology->_Sol[k])(jdofX) + (*sol->_Sol[indexSolD[k]])(jdof);         // coordinates of the deformed configuration
            }
          }
          std::vector <double> xi = p.GetMarkerLocalCoordinates();
          msh->_finiteElement[jelType][solType]->GetPhi(phi, xi);

          std::vector < double > vel(dim, 0.);
          for(unsigned k = 0; k < dim; k++) {
            for(unsigned j = 0; j < nDofs; j++)    {
              vel[k] += phi[j] * solV[k][j];
            }
          }
          MPI_Send(&vel[0], dim, MPI_DOUBLE, jproc, 1, PETSC_COMM_WORLD);

        }
        if(iproc == jproc) {
          std::vector < double > vel(dim);
          MPI_Recv(&vel[0], dim, MPI_DOUBLE, mproc, 1, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
          for(unsigned k = 0; k < dim; k++) {
            sol->_Sol[indexSolV[k]]->set(idof[i], vel[k]);
          }
          counter++;
        }
      }
    }
  }

  for(unsigned k = 0; k < dim; k++) {
    sol->_Sol[indexSolV[k]]->close();
  }

  MPI_Reduce(&counter, &counterAll, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  std::cout << "COUNTER = " << counterAll << " " << msh->GetTotalNumberOfDofs(solType) << std::endl;

}


















