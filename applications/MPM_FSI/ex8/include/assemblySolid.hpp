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
  vector< vector< adept::adouble > > aResD(dim);     // local redidual vector
  vector< adept::adouble > aResP;    // local redidual vector
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
  double muMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_lame_shear_modulus();
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
        aResD[k].assign(nDofs, 0.);
        vx[k].resize(nDofs);
        vxHat[k].resize(nDofs);
      }
      solP.resize(nDofsP);
      aResP.assign(nDofsP, 0.);

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
          vx[k][i] = (*msh->_topology->_Sol[k])(idofX) +  (1. - af) * solD[k][i] + af * solDOld[k][i];
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
              gradSolDgHat[k][j] += ( (1. - af) * solD[k][i] + af * solDOld[k][i]) * gradPhiDHat[i * dim + j];
            }
          }
        }



        std::vector < adept::adouble > solAgAm(dim);
        for(unsigned  k = 0; k < dim; k++) {
          adept::adouble solAgk = (solDg[k] - solDgOld[k]) / (beta * dt * dt) - solVgOld[k] / (beta * dt) + solAgOld[k] * (beta - 0.5) / beta; //NEWMARK ACCELERATION
          solAgAm[k] = (1. - am) * solAgk + am * solAgOld[k];
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
        adept::adouble sigma[3][3];


        for(unsigned j = 0; j < 3; j++) {
          for(unsigned k = 0; k < 3; k++) {
            sigma[j][k] = lambdaMpm * log(J_hat) / J_hat * Id2th[j][k] + muMpm / J_hat * (B[j][k] - Id2th[j][k]);    // alternative formulation
          }
        }

        //END computation of the Cauchy Stress
        for(unsigned i = 0; i < nDofs; i++) {
          for(unsigned k = 0.; k < dim; k++) {
            adept::adouble cauchy = 0;
            for(unsigned j = 0.; j < dim; j++) {
              cauchy += gradPhiD[i * dim + j] * sigma[k][j];
            }
            aResD[k][i] += (rhoMpm * phiD[i] * solAgAm[k] + cauchy + 0. * rhoMpm * phiD[i] * (k == 1)) * weightDHat;
          }
        }
        //continuity block
        for(unsigned i = 0; i < nDofsP; i++) {
          aResP[i] += (phiP[i] * solPg) * weightD;
        }
      } // end gauss point loop

      //copy the value of the adept::adoube aRes in double Res and store them in RES
      rhs.resize(nDofsAll);   //resize

      for(int i = 0; i < nDofs; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          rhs[k * nDofs + i] = -aResD[k][i].value();
        }
      }
      for(int i = 0; i < nDofsP; i++) {
        rhs[ dim * nDofs  + i] = -aResP[i].value();
      }

      myRES->add_vector_blocked(rhs, sysDofsAll);


      Jac.resize(nDofsAll * nDofsAll);
      // define the dependent variables

      for(unsigned  k = 0; k < dim; k++) {
        s.dependent(&aResD[k][0], nDofs);
      }
      s.dependent(&aResP[0], nDofsP);

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
  vector< unsigned > solV2dofs;
  vector< unsigned > solP2dofs;


  vector< vector< double > > solD1(dim);      // local solution (displacement)
  vector< vector< double > > solD1Old(dim);      // local solution (displacement)
  vector< vector< double > > solV1Old(dim);      // local solution (displacement)
  vector< vector< double > > solA1Old(dim);      // local solution (displacement)
  vector< unsigned > solD1dofs;      // local solution (displacement)

  vector <vector < adept::adouble> > vx1(dim);   //vx1 are goind
  vector <vector < double> > vx1Hat(dim); //1 solid 2 is fluid
  vector <vector < double> > vx2Hat(dim); //1 solid 2 is fluid

  vector< double > rhs;    // local redidual vector
  vector< vector< adept::adouble > > aResD(dim);     // local redidual vector
  vector< vector< adept::adouble > > aResV(dim);     // local redidual vector
  vector< adept::adouble > aResP;    // local redidual vector
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


  //reading parameters for fluid FEM domain
  double rhoFluid = ml_prob.parameters.get<Fluid> ("FluidFEM").get_density();
  double muFluid = ml_prob.parameters.get<Fluid> ("FluidFEM").get_viscosity();

  double dt =  my_nnlin_impl_sys.GetIntervalTime();

  std::cout.precision(10);

  //variable-name handling
  const char varname[12][5] = {"UX", "UY", "UZ", "UX", "UY", "UZ", "VX", "VY", "VZ", "AX", "AY", "AZ"};

  vector <unsigned> indexSolD(dim);
  vector <unsigned> indexSolV(dim);

  vector <unsigned> indexSolVOld(dim);
  vector <unsigned> indexSolAOld(dim);

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
                solV1Old[k].resize(nDofsD1);
                solA1Old[k].resize(nDofsD1);
                vx1Hat[k].resize(nDofsD1);
              }
              solD1dofs.resize(dim * nDofsD1);
              nodeFlag1.resize(nDofsD1);
              for(unsigned i = 0; i < nDofsD1; i++) {
                unsigned idofD = msh->GetSolutionDof(i, iel1, solType); //global dof for solution D
                unsigned idofX = msh->GetSolutionDof(i, iel1, 2); //global dof for mesh coordinates

                nodeFlag1[i] = (*mysolution->_Sol[nflagIndex])(idofD); // set it to 0 for no-marker

                for(unsigned  k = 0; k < dim; k++) {
                  solD1[k][i] = (*mysolution->_Sol[indexSolD[k]])(idofD);
                  solD1Old[k][i] = (*mysolution->_SolOld[indexSolD[k]])(idofD);

                  solV1Old[k][i] = (*mysolution->_Sol[indexSolVOld[k]])(idofD);
                  solA1Old[k][i] = (*mysolution->_Sol[indexSolAOld[k]])(idofD);

                  solD1dofs[k * nDofsD1 + i] = myLinEqSolver->GetSystemDof(indexSolD[k], indexPdeD[k], i, iel1);
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
            solV1Old[k].resize(nDofsD1);
            solA1Old[k].resize(nDofsD1);
            vx1Hat[k].resize(nDofsD1);
          }
          solD1dofs.resize(dim * nDofsD1);
          nodeFlag1.resize(nDofsD1);
        }

        for(unsigned  k = 0; k < dim; k++) {
          MPI_Bcast(solD1[k].data(), solD1[k].size(), MPI_DOUBLE, kproc, PETSC_COMM_WORLD);
          MPI_Bcast(solD1Old[k].data(), solD1Old[k].size(), MPI_DOUBLE, kproc, PETSC_COMM_WORLD);
          MPI_Bcast(solV1Old[k].data(), solV1Old[k].size(), MPI_UNSIGNED, kproc, PETSC_COMM_WORLD);
          MPI_Bcast(solA1Old[k].data(), solA1Old[k].size(), MPI_UNSIGNED, kproc, PETSC_COMM_WORLD);
          MPI_Bcast(vx1Hat[k].data(), vx1Hat[k].size(), MPI_DOUBLE, kproc, PETSC_COMM_WORLD);
        }
        MPI_Bcast(solD1dofs.data(), solD1dofs.size(), MPI_UNSIGNED, kproc, PETSC_COMM_WORLD);
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
              unsigned if2e = el->GetIG(ielt1, jface1, i); // local mapping from face to element
              for(unsigned k = 0; k < dim; k++) {
                xf1[k][i] =  vx1Hat[k][if2e] + (1. - af) * solD1[k][if2e] + af * solD1Old[k][if2e];
              }
            }
            unsigned coarse = 0;
            unsigned fine = 1;
            unsigned qType = coarse;
            unsigned ng1;
            std::vector < Marker *> gp;
            std::vector < unsigned> jel;

          testGaussPoints:
            ng1 = fem[qType][jfaceType1][solType]->GetGaussPointNumber();
            gp.resize(ng1);
            jel.resize(ng1);

            std::map <unsigned, unsigned> jelCounter;
            std::map <unsigned, unsigned>::iterator it;

            for(unsigned ig = 0; ig < ng1; ig++) {

              const double* phiD1 = fem[qType][jfaceType1][solType]->GetPhi(ig);

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

            if(qType == coarse && jelCounter.size() > 1) {
              qType = fine;
              for(unsigned ig = 0; ig < gp.size(); ig++) {
                delete gp[ig];
              }
              goto testGaussPoints;
            }

            if(qType == fine) {
              for(unsigned ig = 0; ig < gp.size(); ig++) {
                std::cout << jel[ig] << " ";
              }
              std::cout << std::endl;
            }

            std::vector < unsigned > iel2All(jelCounter.size());
            std::vector < std::vector < unsigned > > igAll(jelCounter.size());

            {
              unsigned i = 0;
              for(i = 0, it = jelCounter.begin(); it != jelCounter.end(); it++, i++) {
                iel2All[i] = it->first;
                igAll[i].reserve(it->second);
              }
            }

            for(unsigned ig = 0; ig < gp.size(); ig++) {
              unsigned i = 0;
              while(jel[ig] != iel2All[i]) i++;
              unsigned k = igAll[i].size();
              igAll[i].resize(k + 1);
              igAll[i][k] = ig;
            }


            if( iel2All.size() > 1 ) {


              for(unsigned ii = 0; ii < iel2All.size(); ii++) {

                unsigned iel2 = iel2All[ii];

                //std::cout << ii << " " << iel2 << std::endl;
                std::vector < std::vector < double > > xi2(igAll[ii].size());
                for(unsigned ig = 0; ig < igAll[ii].size(); ig++) {
                  xi2[ig] = gp[igAll[ii][ig]]->GetMarkerLocalCoordinates();
                  //std::cout <<  igAll[ii][ig] << " " << xi2[ig][0] << " " << xi2[ig][1] << std::endl;
                }


                std::vector < std::vector < double > > xp(igAll[ii].size());
                for(unsigned ig = 0; ig < igAll[ii].size(); ig++) {
                  xp[ig] = gp[igAll[ii][ig]]->GetIprocMarkerCoordinates();
                  //std::cout << xp[ig][0] << " " << xp[ig][1] << std::endl;
                }

                unsigned nDofsV2 = msh->GetElementDofNumber(iel2, solType);
                for(unsigned  k = 0; k < dim; k++) {
                  vx2Hat[k].resize(nDofsV2);
                }
                for(unsigned i = 0; i < nDofsV2; i++) {
                  unsigned idofX = msh->GetSolutionDof(i, iel2, 2); //global dof for mesh coordinates
                  for(unsigned  k = 0; k < dim; k++) {
                    vx2Hat[k][i] = (*msh->_topology->_Sol[k])(idofX);
                  }
                  //std::cout << vx2Hat[0][i] << " " << vx2Hat[1][i] << std::endl;
                }


              }

              //exit(0);
            }


            for(unsigned ii = 0; ii < iel2All.size(); ii++) {
              unsigned iel2 = iel2All[ii];
              if(iel2 != UINT_MAX) {
                unsigned jproc = msh->IsdomBisectionSearch(iel2 , 3); // return  jproc for piece-wise constant discontinuous type (3)

                unsigned ielt2;
                unsigned nDofsV2;
                unsigned nDofsP2;

                std::vector < std::vector < double > > xi2(igAll[ii].size());
                if(iproc == jproc) {

                  for(unsigned ig = 0; ig < igAll[ii].size(); ig++) {
                    xi2[ig] = gp[igAll[ii][ig]]->GetMarkerLocalCoordinates();
                  }

                  // std::cout << xi2[0] <<" "<<" "xi2[1] << std::endl;

                  mysolution->_Sol[eflagIndex]->set(iel2, 1.);

                  ielt2 = msh->GetElementType(iel2);
                  nDofsV2 = msh->GetElementDofNumber(iel2, solType);
                  nDofsP2 = msh->GetElementDofNumber(iel2, solTypeP);
                  for(unsigned  k = 0; k < dim; k++) {
                    solV2[k].resize(nDofsV2);
                    solV2Old[k].resize(nDofsV2);
                    vx2Hat[k].resize(nDofsV2);
                    aResV[k].assign(nDofsV2, 0.);
                  }
                  solV2dofs.resize(dim * nDofsV2);
                  solP2.resize(nDofsP2);
                  aResP.assign(nDofsP2, 0.);
                  solP2dofs.resize(nDofsP2);

                  for(unsigned i = 0; i < nDofsV2; i++) {
                    unsigned idofV = msh->GetSolutionDof(i, iel2, solType); //global dof for solution D
                    unsigned idofX = msh->GetSolutionDof(i, iel2, 2); //global dof for mesh coordinates
                    for(unsigned  k = 0; k < dim; k++) {
                      solV2[k][i] = (*mysolution->_Sol[indexSolV[k]])(idofV);
                      solV2Old[k][i] = (*mysolution->_SolOld[indexSolV[k]])(idofV);
                      solV2dofs[k * nDofsV2 + i] = myLinEqSolver->GetSystemDof(indexSolV[k], indexPdeV[k], i, iel2);
                      vx2Hat[k][i] = (*msh->_topology->_Sol[k])(idofX);
                    }
                  }
                  for(unsigned i = 0; i < nDofsP2; i++) {
                    unsigned idofP = msh->GetSolutionDof(i, iel2, solTypeP); //global dof for solution D
                    solP2[i] = (*mysolution->_Sol[indexSolP])(idofP);
                    solP2dofs[i] = myLinEqSolver->GetSystemDof(indexSolP, indexPdeP, i, iel2);
                  }

                  if(kproc != jproc) {
                    MPI_Send(&ielt2, 1, MPI_UNSIGNED, kproc, 0, PETSC_COMM_WORLD);
                    MPI_Send(&nDofsV2, 1, MPI_UNSIGNED, kproc, 1, PETSC_COMM_WORLD);
                    MPI_Send(&nDofsP2, 1, MPI_UNSIGNED, kproc, 2, PETSC_COMM_WORLD);

                    for(unsigned  k = 0; k < dim; k++) {
                      MPI_Send(solV2[k].data(), solV2[k].size(), MPI_DOUBLE, kproc, 3 + (k * 3), PETSC_COMM_WORLD);
                      MPI_Send(solV2Old[k].data(), solV2Old[k].size(), MPI_DOUBLE, kproc, 4 + (k * 3), PETSC_COMM_WORLD);
                      MPI_Send(vx2Hat[k].data(), vx2Hat[k].size(), MPI_DOUBLE, kproc, 5 + (k * 3), PETSC_COMM_WORLD);
                    }
                    MPI_Send(solV2dofs.data(), solV2dofs.size(), MPI_UNSIGNED, kproc, 3 + (dim * 3), PETSC_COMM_WORLD);
                    MPI_Send(solP2.data(), solP2.size(), MPI_DOUBLE, kproc, 4 + (dim * 3), PETSC_COMM_WORLD);
                    MPI_Send(solP2dofs.data(), solP2dofs.size(), MPI_UNSIGNED, kproc, 5 + (dim * 3), PETSC_COMM_WORLD);

                    for(unsigned ig = 0; ig < xi2.size(); ig++) {
                      MPI_Send(xi2[ig].data(), xi2[ig].size(), MPI_DOUBLE, kproc, 6 + (dim * 3) + ig, PETSC_COMM_WORLD);
                    }
                  }
                }

                if(kproc != jproc && iproc == kproc) {
                  MPI_Recv(&ielt2, 1, MPI_UNSIGNED, jproc, 0, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                  MPI_Recv(&nDofsV2, 1, MPI_UNSIGNED, jproc, 1, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                  MPI_Recv(&nDofsP2, 1, MPI_UNSIGNED, jproc, 2, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);

                  for(unsigned  k = 0; k < dim; k++) {
                    solV2[k].resize(nDofsV2);
                    solV2Old[k].resize(nDofsV2);
                    vx2Hat[k].resize(nDofsV2);
                  }
                  solV2dofs.resize(dim * nDofsV2);
                  solP2.resize(nDofsP2);
                  solP2dofs.resize(nDofsP2);

                  for(unsigned  k = 0; k < dim; k++) {
                    MPI_Recv(solV2[k].data(), solV2[k].size(), MPI_DOUBLE, jproc, 3 + (k * 3), PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(solV2Old[k].data(), solV2Old[k].size(), MPI_DOUBLE, jproc, 4 + (k * 3), PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(vx2Hat[k].data(), vx2Hat[k].size(), MPI_DOUBLE, jproc, 5 + (k * 3), PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                  }
                  MPI_Recv(solV2dofs.data(), solV2dofs.size(), MPI_UNSIGNED, jproc, 3 + (dim * 3), PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                  MPI_Recv(solP2.data(), solP2.size(), MPI_DOUBLE, jproc, 4 + (dim * 3), PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                  MPI_Recv(solP2dofs.data(), solP2dofs.size(), MPI_UNSIGNED, jproc, 5 + (dim * 3), PETSC_COMM_WORLD, MPI_STATUS_IGNORE);

                  for(unsigned ig = 0; ig < xi2.size(); ig++) {
                    xi2[ig].resize(dim);
                    MPI_Recv(xi2[ig].data(), xi2[ig].size(), MPI_DOUBLE, jproc, 6 + (dim * 3) + ig, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                  }
                }
                // we finished to exchange all needed information, we are ready to assemble!
                if(iproc == kproc || iproc == jproc) {

                  if(iproc == kproc) {
                    for(unsigned k = 0; k < dim; k++) {
                      aResD[k].assign(nDofsD1, 0.);
                    }
                  }


                  //initialize adept varibles
                  for(unsigned k = 0; k < dim; k++) {
                    solD[k].resize(nDofsD1) ;
                    for(unsigned i = 0; i < nDofsD1; i++) {
                      solD[k][i] = solD1[k][i];
                    }
                    solV[k].resize(nDofsV2);
                    for(unsigned i = 0; i < nDofsV2; i++) {
                      solV[k][i] = solV2[k][i];
                    }
                  }
                  solP.resize(nDofsP2);
                  for(unsigned i = 0; i < nDofsP2; i++) {
                    solP[i] = solP2[i];
                  }

                  s.new_recording();

                  //solid
                  std::vector  < std::vector  <  adept::adouble > > xf1(dim);    // physical coordinates of the face in the af configuration
                  for(int k = 0; k < dim; k++) {
                    xf1[k].resize(jfaceDofs1);
                  }

                  for(unsigned i = 0; i < jfaceDofs1; i++) {
                    unsigned if2e = el->GetIG(ielt1, jface1, i); // local mapping from face to element
                    for(unsigned k = 0; k < dim; k++) {
                      xf1[k][i] =  vx1Hat[k][if2e] + (1. - af) * solD[k][if2e] + af * solD1Old[k][if2e];
                    }
                  }

                  for(unsigned ig = 0; ig < igAll[ii].size(); ig++) { // loop on the interface integration points both on the solid and the fluid
                    unsigned ig1 = igAll[ii][ig]; //solid gauss point id
                    unsigned ig2 = ig; //integration point id in the ii fluid element

                    std::vector < adept::adouble > Nsf;
                    fem[qType][jfaceType1][solType]->JacobianSur(xf1, ig1, weightD, phiD, gradPhiD, Nsf); // boundary solid integration

                    std::vector < adept::adouble > N(dim); // normal from the fluid to the solid domain
                    for(unsigned k = 0; k < dim; k++) {
                      N[k] = -Nsf[k];
                    }


                    //std::cout << ii << " " << ig2 << " " << ig1 << " " << xi2[ig2][0] << " " << xi2[ig2][1] << std::endl;

                    fem[qType][ielt2][solType]->Jacobian(vx2Hat, xi2[ig2], weightV, phiV, gradPhiV); //cut fem fluid cell
                    fem[qType][ielt2][solTypeP]->GetPhi(phiP, xi2[ig2]); //cut fem fluid cell

                    adept::adouble solPg = 0.;
                    for(unsigned i = 0; i < nDofsP2; i++) {
                      solPg += phiP[i] * solP[i];
                    }


                    std::vector < adept::adouble > xg1(dim, 0.);
                    std::vector < adept::adouble > xg2(dim, 0.);

                    std::vector < adept::adouble > vs(dim, 0.); // this is the velocity of the solid
                    std::vector < adept::adouble > vf(dim, 0.); // this is the velocity of the fluid


                    std::vector <adept::adouble> solDg(dim, 0.);
                    std::vector <double> solD1gOld(dim, 0.);
                    std::vector <double> solV1gOld(dim, 0.);
                    std::vector <double> solA1gOld(dim, 0.);
                    std::vector <adept::adouble> solA1g(dim);
                    std::vector <adept::adouble> solV1g(dim);

                    //update displacement and acceleration
                    for(int k = 0; k < dim; k++) {
                      for(unsigned i = 0; i < nDofsV2; i++) {
                        vf[k] += phiV[i] * (theta * solV[k][i] + (1. - theta) * solV2Old[k][i]);
                        xg2[k] += phiV[i] * vx2Hat[k][i];
                      }

                      for(unsigned i = 0; i < jfaceDofs1; i++) {
                        unsigned if2e = el->GetIG(ielt1, jface1, i); // local mapping from face to element
                        solDg[k] += phiD[i] * solD[k][if2e];
                        solD1gOld[k] += phiD[i] * solD1Old[k][if2e];
                        solV1gOld[k] += phiD[i] * solV1Old[k][if2e];
                        solA1gOld[k] += phiD[i] * solA1Old[k][if2e];
                        xg1[k] += phiD[i] * xf1[k][i];
                      }

                      solA1g[k] = (solDg[k] - solD1gOld[k]) / (beta * dt * dt)
                                  - solV1gOld[k] / (beta * dt)
                                  + solA1gOld[k] * (beta - 0.5) / beta;
                      solV1g[k] = solV1gOld[k] + dt * ((1. - Gamma) * solA1gOld[k] + Gamma * solA1g[k]);

                      vs[k] = (1. - af) * solV1g[k] +  af * solV1gOld[k];

                    }

                    if( (xg1[0] - xg2[0]) * (xg1[0] - xg2[0]) + (xg1[1] - xg2[1]) * (xg1[1] - xg2[1]) > 1.0e-14 ) {
                      std::cout << ii << iel1 << " " << iel2 << " " << ig1 << " " << ig2 << " " << xg1[0] - xg2[0] << " " << xg1[1] - xg2[1] << std::endl;
                    }


                    std::vector < adept::adouble > tau(dim, 0.);
                    for(unsigned k = 0; k < dim; k++) {
                      tau[k] += solPg * N[k];
                      for(unsigned i = 0; i < nDofsV2; i++) {
                        for(unsigned j = 0; j < dim; j++) {
                          tau[k] += -muFluid * ((theta * solV[k][i] + (1. - theta) * solV2Old[k][i]) * gradPhiV[i * dim + j] +
                                                (theta * solV[j][i] + (1. - theta) * solV2Old[j][i]) * gradPhiV[i * dim + k]) * N[j];
                        }
                      }
                    }

//                     std::cout<<ig << " " << xi2[ig2][0] << " " << xi2[ig2][1] <<" " << tau[0] << " " << tau[1] <<std::endl;
//                     std::cout<<ig << " " << vf[0] << " " << vs[0] <<" " << vf[1] << " " << vs[1] <<std::endl;

                    double h = sqrt((vx2Hat[0][0] - vx2Hat[0][2]) * (vx2Hat[0][0] - vx2Hat[0][2]) +
                                    (vx2Hat[1][0] - vx2Hat[1][2]) * (vx2Hat[1][0] - vx2Hat[1][2])) ;

                    double thetaM = 1 * muFluid / h;
                    double thetaL = 1 * (rhoFluid * h / (theta * dt) + muFluid / h);
                    //thetaM = thetaL;

                    //std::cout << thetaM << " " << thetaL << std::endl;

                    if(iproc == kproc) {
                      for(unsigned i = 0; i < jfaceDofs1; i++) {
                        unsigned if2e = el->GetIG(ielt1, jface1, i); // local mapping from face to element
                        for(unsigned k = 0; k < dim; k++) {
                          aResD[k][if2e] += tau[k] * (-phiD[i]) * weightD;
                          aResD[k][if2e] += thetaM * (vf[k] - vs[k]) * (-phiD[i]) * weightD;
                          for(unsigned j = 0; j < dim; j++) {
                            aResD[k][if2e] +=  thetaL * (vf[j] - vs[j]) * N[j] * (-phiD[i]) * N[k] * weightD;
                          }
                        }
                      }
                    }

                    if(iproc == jproc) {
                      for(unsigned i = 0; i < nDofsV2; i++) {
                        for(unsigned k = 0; k < dim; k++) {
                          aResV[k][i] +=  (tau[k] - solPg * N[k]) * phiV[i] * weightD;
                          aResV[k][i] +=  thetaM * (vf[k] - vs[k]) * phiV[i] * weightD;

                          for(unsigned j = 0; j < dim; j++) {
                            aResV[k][i] += (muFluid * gradPhiV[i * dim + j] * N[j] * (vf[k] - vs[k])) * weightD;
                            aResV[k][i] += (muFluid * gradPhiV[i * dim + j] * N[k] * (vf[j] - vs[j])) * weightD;
                            aResV[k][i] +=  thetaL * (vf[j] - vs[j]) * N[j] * phiV[i] * N[k] * weightD;
                          }
                        }
                      } // end phi_i loop

                      for(unsigned i = 0; i < nDofsP2; i++) {
                        for(unsigned k = 0; k < dim; k++) {
                          aResP[i] += - (phiP[i]  * (vf[k] - vs[k]) * N[k]) * weightD;
                        }
                      }
                    }
                  }//close ig loop

                  //call adept

                  sysDofsAll = solD1dofs;
                  sysDofsAll.insert(sysDofsAll.end(), solV2dofs.begin(), solV2dofs.end());
                  sysDofsAll.insert(sysDofsAll.end(), solP2dofs.begin(), solP2dofs.end());

                  if(iproc == kproc) {
                    unsigned nDofsD1All = dim * nDofsD1;

                    rhs.resize(nDofsD1All);   //resize

                    for(int i = 0; i < nDofsD1; i++) {
                      for(unsigned  k = 0; k < dim; k++) {
                        rhs[ i +  k * nDofsD1 ] = -aResD[k][i].value();
                      }
                    }

                    myRES->add_vector_blocked(rhs, solD1dofs);
                    unsigned nDofsAll = dim * (nDofsD1 + nDofsV2) + nDofsP2;


                    Jac.resize(nDofsD1All * nDofsAll);
                    // define the dependent variables

                    for(unsigned  k = 0; k < dim; k++) {
                      s.dependent(&aResD[k][0], nDofsD1);
                    }

                    // define the independent variables
                    for(unsigned  k = 0; k < dim; k++) {
                      s.independent(&solD[k][0], nDofsD1);
                    }
                    for(unsigned  k = 0; k < dim; k++) {
                      s.independent(&solV[k][0], nDofsV2);
                    }
                    s.independent(&solP[0], nDofsP2);

                    // get the and store jacobian matrix (row-major)
                    s.jacobian(&Jac[0], true);
                    myKK->add_matrix_blocked(Jac, solD1dofs, sysDofsAll);

                    s.clear_independents();
                    s.clear_dependents();
                  }

                  if(iproc == jproc) {

                    std::vector < unsigned >solVP2dofs = solV2dofs;
                    solVP2dofs.insert(solVP2dofs.end(), solP2dofs.begin(), solP2dofs.end());

                    unsigned nDofsVP2All = dim * nDofsV2 + nDofsP2; // number of rows in the fluid block

                    rhs.resize(nDofsVP2All);   //resize

                    for(int i = 0; i < nDofsV2; i++) {
                      for(unsigned  k = 0; k < dim; k++) {
                        rhs[ i +  k * nDofsV2 ] = -aResV[k][i].value();
                      }
                    }
                    for(int i = 0; i < nDofsP2; i++) {
                      rhs[ i +  dim * nDofsV2 ] = -aResP[i].value();
                    }
                    myRES->add_vector_blocked(rhs, solVP2dofs);

                    unsigned nDofsAll = dim * (nDofsD1 + nDofsV2) + nDofsP2; // number of columns

                    Jac.resize(nDofsVP2All * nDofsAll);
                    // define the dependent variables

                    for(unsigned  k = 0; k < dim; k++) {
                      s.dependent(&aResV[k][0], nDofsV2);
                    }
                    s.dependent(&aResP[0], nDofsP2);

                    // define the independent variables
                    for(unsigned  k = 0; k < dim; k++) {
                      s.independent(&solD[k][0], nDofsD1);
                    }
                    for(unsigned  k = 0; k < dim; k++) {
                      s.independent(&solV[k][0], nDofsV2);
                    }
                    s.independent(&solP[0], nDofsP2);

                    // get the and store jacobian matrix (row-major)
                    s.jacobian(&Jac[0], true);
                    myKK->add_matrix_blocked(Jac, solVP2dofs, sysDofsAll);

                    s.clear_independents();
                    s.clear_dependents();
                  }
                }
              }




            }

            //if(iel2All.size() > 1) exit(0);

            for(unsigned ig = 0; ig < gp.size(); ig++) {
              unsigned iel2 = gp[ig]->GetMarkerElement();
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

  // *************************************
  end_time = clock();
  AssemblyTime += (end_time - start_time);
  // ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

}




