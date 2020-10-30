#include "MultiLevelSolution.hpp"
#include "PetscMatrix.hpp"

using namespace femus;


void AssembleBoundaryLayer(MultiLevelProblem& ml_prob) {

  // ml_prob is the global object from/to where get/set all the data
  // level is the level of the PDE system to be assembled
  // levelMax is the Maximum level of the MultiLevelProblem
  // assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  clock_t AssemblyTime = 0;
  clock_t start_time, end_time;

  //pointers and references

  TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("MPM_FSI");
  const unsigned  level = my_nnlin_impl_sys.GetLevelToAssemble();
  MultiLevelSolution* mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* mysolution = mlSol->GetSolutionLevel(level);     // pointer to the solution (level) object

  //Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

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

  vector < double > gradPhiV;
  vector < adept::adouble > gradPhiD;
  vector < double > gradPhiDHat;


  vector <vector < adept::adouble> > vx(dim);   //vx is coordX in assembly of ex30
  vector <vector < double> > vxHat(dim);

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

    if(eFlag == 3) { // porous media

      for(unsigned  k = 0; k < dim; k++) {
        solD[k].resize(nDofs);
        solDOld[k].resize(nDofs);
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
          sysDofsAll[k * nDofs + i] = myLinEqSolver->GetSystemDof(indexSolD[k], indexPdeD[k], i, iel);
        }
      }

      for(unsigned i = 0; i < nDofsP; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, solTypeP);
        solP[i] = (*mysolution->_Sol[indexSolP])(idof);
        sysDofsAll[dim * nDofs + i] = myLinEqSolver->GetSystemDof(indexSolP, indexPdeP, i, iel);
      }

      s.new_recording();

      for(unsigned i = 0; i < nDofs; i++) {
        unsigned idofX = msh->GetSolutionDof(i, iel, 2);
        for(unsigned  k = 0; k < dim; k++) {
          vxHat[k][i] = (*msh->_topology->_Sol[k])(idofX);
          vx[k][i] = (*msh->_topology->_Sol[k])(idofX) + af * solD[k][i] + (1. - af) * solDOld[k][i];
        }
      }

      // *** Gauss point loop ***
      for(unsigned ig = 0; ig < msh->_finiteElement[ielt][solType]->GetGaussPointNumber(); ig++) {

        msh->_finiteElement[ielt][solType]->Jacobian(vxHat, ig, weightDHat, phiD, gradPhiDHat);
        msh->_finiteElement[ielt][solType]->Jacobian(vx, ig, weightD, phiD, gradPhiD);

        std::vector < std::vector < adept::adouble > > gradSolDgHat(dim);

        for(unsigned  k = 0; k < dim; k++) {
          gradSolDgHat[k].assign(dim, 0.);
        }

        for(unsigned i = 0; i < nDofs; i++) {
          for(unsigned j = 0; j < dim; j++) {

            for(unsigned  k = 0; k < dim; k++) {
              gradSolDgHat[k][j] += (af * solD[k][i] + (1. - af) * solDOld[k][i]) * gradPhiDHat[i * dim + j];
            }
          }
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


        double E = 1;
        double nu = 0.4;

        double mu = E / (2. * (1. + nu));
        double lambda = (E * nu) / ((1. + nu) * (1. - 2.*nu));

        for(unsigned j = 0; j < 3; j++) {
          for(unsigned k = 0; k < 3; k++) {
            Cauchy[j][k] = lambda * log(J_hat) / J_hat * Id2th[j][k] + mu / J_hat * (B[j][k] - Id2th[j][k]);    // alternative formulation
          }
        }
        //END computation of the Cauchy Stress
        for(unsigned i = 0; i < nDofs; i++) {
          if(nodeFlag[i] == 4) {
            adept::adouble CauchyDIR[3] = {0., 0., 0.};
            for(unsigned j = 0.; j < dim; j++) {
              for(unsigned k = 0.; k < dim; k++) {
                CauchyDIR[j] += gradPhiD[i * dim + k] * Cauchy[j][k];
              }
            }
            for(unsigned k = 0; k < dim; k++) {
              aResD[k][i] +=   CauchyDIR[k] * weightD;
            }
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



void AssembleBoundaryLayerProjection(MultiLevelProblem& ml_prob) {

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
//   vector< vector< double > > solV1Old(dim);      // local solution (displacement)
//   vector< vector< double > > solA1Old(dim);      // local solution (displacement)
  vector< unsigned > solD1dofs;      // local solution (displacement)

  vector <vector < adept::adouble> > vx1(dim);   //vx1 are goind
  vector <vector < double> > vx1Hat(dim); //1 solid 2 is fluid
  vector <vector < double> > vx2Hat(dim); //1 solid 2 is fluid

  vector< double > rhs;    // local redidual vector
  //vector< vector< adept::adouble > > aRhsD(dim);     // local redidual vector
  vector< vector< adept::adouble > > aResV(dim);     // local redidual vector
  vector< adept::adouble > aResP;    // local redidual vector
  vector < double > Jac;

  std::vector <unsigned> sysDofsAll;

  vector < double > phiV;
  vector < double > phiD;
  vector < double > phiP;

// vector < double > gradPhiP;
  vector < double > gradPhiV;

  vector < adept::adouble > gradPhiD;
  //vector < double > gradPhiDHat;

  unsigned dim2 = 3 * (dim - 1);



  //double weightP;
  double weightV;
// double weightDHat;
  adept::adouble weightD;


  //reading parameters for MPM body
//   double rhoMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_density();
//   double EMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_young_module();
//   double muMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_lame_shear_modulus();
//   double nuMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_poisson_coeff();
//   double lambdaMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_lame_lambda();


  //reading parameters for fluid FEM domain
  double rhoFluid = ml_prob.parameters.get<Fluid> ("FluidFEM").get_density();
  double muFluid = ml_prob.parameters.get<Fluid> ("FluidFEM").get_viscosity();

  double dt =  my_nnlin_impl_sys.GetIntervalTime();

  std::cout.precision(10);

  //variable-name handling
  const char varname[12][5] = {"UX", "UY", "UZ", "UX", "UY", "UZ", "VX", "VY", "VZ", "AX", "AY", "AZ"};

  vector <unsigned> indexSolD(dim);
  vector <unsigned> indexSolV(dim);

//   vector <unsigned> indexSolAOld(dim);
//   vector <unsigned> indexSolVOld(dim);

  vector <unsigned> indexPdeD(dim);
  vector <unsigned> indexPdeV(dim);
  for(unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolD[ivar] = mlSol->GetIndex(&varname[ivar][0]);
    indexSolV[ivar] = mlSol->GetIndex(&varname[ivar + 3][0]);
//     indexSolVOld[ivar] = mlSol->GetIndex(&varname[ivar + 6][0]); // For Newmark in Solid
//     indexSolAOld[ivar] = mlSol->GetIndex(&varname[ivar + 9][0]); // For Newmark in Solid
    indexPdeD[ivar] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar][0]);     // DX, DY, DZ
    indexPdeV[ivar] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar + 3][0]); //VX, VY, VZ
  }
  unsigned solType = mlSol->GetSolutionType(&varname[0][0]);

  unsigned indexSolP = mlSol->GetIndex("P");
  unsigned indexPdeP = my_nnlin_impl_sys.GetSolPdeIndex("P");
  unsigned solTypeP = mlSol->GetSolutionType("P");

  unsigned eflagIndex = mlSol->GetIndex("eflag");
  unsigned pElemIndex = mlSol->GetIndex("pElem");

  start_time = clock();

  unsigned meshOffset = msh->_elementOffset[iproc];
  unsigned meshOffsetp1 = msh->_elementOffset[iproc + 1];

  //BEGIN loop on elements (to initialize the "soft" stiffness matrix)
  for(unsigned kproc = 0; kproc < nprocs; kproc++) {
    for(int iel1 = msh->_elementOffset[kproc]; iel1 < msh->_elementOffset[kproc + 1]; iel1++) {

      short unsigned ielt1;
      unsigned nDofsD1;
      unsigned pElem;
      unsigned eFlag1 = 0;

      if(iproc == kproc) {
        eFlag1 = static_cast <unsigned>(floor((*mysolution->_Sol[eflagIndex])(iel1) + 0.5));
        if(eFlag1 == 3) { // boundary Layer
          ielt1 = msh->GetElementType(iel1);
          nDofsD1 = msh->GetElementDofNumber(iel1, solType);
          pElem = static_cast <unsigned>(floor((*mysolution->_Sol[pElemIndex])(iel1) + 0.5));
          for(unsigned  k = 0; k < dim; k++) {
            solD1[k].resize(nDofsD1);
            solD1Old[k].resize(nDofsD1);
            vx1Hat[k].resize(nDofsD1);
          }
          solD1dofs.resize(dim * nDofsD1);
          for(unsigned i = 0; i < nDofsD1; i++) {
            unsigned idofD = msh->GetSolutionDof(i, iel1, solType); //global dof for solution D
            unsigned idofX = msh->GetSolutionDof(i, iel1, 2); //global dof for mesh coordinates

            for(unsigned  k = 0; k < dim; k++) {
              solD1[k][i] = (*mysolution->_Sol[indexSolD[k]])(idofD);
              solD1Old[k][i] = (*mysolution->_SolOld[indexSolD[k]])(idofD);
              vx1Hat[k][i] = (*msh->_topology->_Sol[k])(idofX);
              solD1dofs[k * nDofsD1 + i] = myLinEqSolver->GetSystemDof(indexSolD[k], indexPdeD[k], i, iel1);
            }
          }
        }
      }


      MPI_Bcast(&eFlag1, 1, MPI_UNSIGNED, kproc, PETSC_COMM_WORLD);

      if(eFlag1 == 3) {

        MPI_Bcast(&ielt1, 1, MPI_UNSIGNED, kproc, PETSC_COMM_WORLD);
        MPI_Bcast(&nDofsD1, 1, MPI_UNSIGNED, kproc, PETSC_COMM_WORLD);
        MPI_Bcast(&pElem, 1, MPI_UNSIGNED, kproc, PETSC_COMM_WORLD);

        if(iproc != kproc) {
          for(unsigned  k = 0; k < dim; k++) {
            solD1[k].resize(nDofsD1);
            solD1Old[k].resize(nDofsD1);
            vx1Hat[k].resize(nDofsD1);
          }
          solD1dofs.resize(dim * nDofsD1);
        }


        for(unsigned  k = 0; k < dim; k++) {
          MPI_Bcast(solD1[k].data(), solD1[k].size(), MPI_DOUBLE, kproc, PETSC_COMM_WORLD);
          MPI_Bcast(solD1Old[k].data(), solD1Old[k].size(), MPI_DOUBLE, kproc, PETSC_COMM_WORLD);
          MPI_Bcast(vx1Hat[k].data(), vx1Hat[k].size(), MPI_DOUBLE, kproc, PETSC_COMM_WORLD);
        }
        MPI_Bcast(solD1dofs.data(), solD1dofs.size(), MPI_UNSIGNED, kproc, PETSC_COMM_WORLD);

        unsigned coarse = 0;
        unsigned fine = 1;
        unsigned qType = fine;
        unsigned ng1;
        std::vector < Marker *> gp;
        std::vector < unsigned> jel;

      testGaussPoints:
        ng1 = fem[qType][ielt1][solType]->GetGaussPointNumber();
        gp.resize(ng1);
        jel.resize(ng1);

        std::map <unsigned, unsigned> jelCounter;
        std::map <unsigned, unsigned>::iterator it;

        for(unsigned ig = 0; ig < ng1; ig++) {

          const double* phiD1 = fem[qType][ielt1][solType]->GetPhi(ig);
          std::vector< double > xg(dim, 0.);
          for(unsigned i = 0; i < nDofsD1; i++) {
            for(unsigned k = 0; k < dim; k++) {
              xg[k] += phiD1[i] * (vx1Hat[k][i] + af * solD1[k][i] + (1. - af) * solD1Old[k][i]);
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

//         if(qType == fine) {
//           for(unsigned ig = 0; ig < gp.size(); ig++) {
//             std::cout << jel[ig] << " ";
//           }
//           std::cout << std::endl;
//         }

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


        for(unsigned ii = 0; ii < iel2All.size(); ii++) {
          unsigned iel2 = iel2All[ii];
          if(iel2 != UINT_MAX) {
            unsigned jproc = msh->IsdomBisectionSearch(iel2 , 3); // return  jproc for piece-wise constant discontinuous type (3)

            if(iproc == jproc) {

              unsigned eFlag2 = static_cast <unsigned>(floor((*mysolution->_Sol[eflagIndex])(iel2) + 0.5));

              if(eFlag2 == 1) {
                unsigned ielt2;
                unsigned nDofsV2;
                unsigned nDofsP2;

                std::vector < std::vector < double > > xi2(igAll[ii].size());
                for(unsigned ig = 0; ig < igAll[ii].size(); ig++) {
                  xi2[ig] = gp[igAll[ii][ig]]->GetMarkerLocalCoordinates();
                }
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

                //solid coordinates in the af domain
                for(unsigned k = 0; k < dim; k++) {
                  vx1[k]. resize(nDofsD1);
                  for(unsigned i = 0; i < nDofsD1; i++) {
                    vx1[k][i] =  vx1Hat[k][i] + af * solD[k][i] + (1. - af) * solD1Old[k][i];
                  }
                }

                for(unsigned ig = 0; ig < igAll[ii].size(); ig++) { // loop on the interface integration points both on the solid and the fluid
                  unsigned ig1 = igAll[ii][ig]; //solid gauss point id
                  unsigned ig2 = ig; //integration point id in the ii fluid element


                  fem[qType][ielt1][solType]->Jacobian(vx1, ig1, weightD, phiD, gradPhiD); // boundary solid integration

                  fem[qType][ielt2][solType]->Jacobian(vx2Hat, xi2[ig2], weightV, phiV, gradPhiV); //cut fem fluid cell
                  fem[qType][ielt2][solTypeP]->GetPhi(phiP, xi2[ig2]); //cut fem fluid cell

                  adept::adouble solPg = 0.;
                  for(unsigned i = 0; i < nDofsP2; i++) {
                    solPg += phiP[i] * solP[i];
                  }

                  std::vector < double > solV2gOld(dim, 0.);
                  std::vector < adept::adouble > solV2gTheta(dim, 0.);
                  std::vector < adept::adouble > solV2g(dim, 0.);
                  std::vector < std::vector < adept::adouble > > gradSolV2g(dim);

                  for(int k = 0; k < dim; k++) {
                    gradSolV2g[k].assign(dim, 0.);
                    for(unsigned i = 0; i < nDofsV2; i++) {
                      solV2gOld[k] += phiV[i] * solV2Old[k][i];
                      solV2gTheta[k] += phiV[i] * (theta * solV[k][i] + (1. - theta) * solV2Old[k][i]);
                      solV2g[k] += phiV[i] * solV[k][i];
                      for(unsigned j = 0; j < dim; j++) {
                        gradSolV2g[k][j] += gradPhiV[i * dim + j] * (theta * solV[k][i] + (1. - theta) * solV2Old[k][i]);
                      }
                    }
                  }

                  std::vector < adept::adouble > xg1(dim, 0.);
                  std::vector < adept::adouble > xg2(dim, 0.);
                 
                  for(int k = 0; k < dim; k++) {
                    for(unsigned i = 0; i < nDofsV2; i++) {
                      xg2[k] += vx2Hat[k][i] *  phiV[i];
                    }
                    for(unsigned i = 0; i < nDofsD1; i++) {
                      xg1[k] += vx1[k][i] *  phiD[i];
                    }
                  }

                  if( (xg1[0] - xg2[0]) * (xg1[0] - xg2[0]) + (xg1[1] - xg2[1]) * (xg1[1] - xg2[1]) > 1.0e-14 ) {
                    std::cout << ii <<" "<< iel1 << " " << iel2 << " " << ig1 << " " << ig2 << " " << xg1[0] - xg2[0] << " " << xg1[1] - xg2[1] << std::endl;
                    std::cout << ii <<" "<< xg1[0]  << " " <<  xg2[0] << " "<< xg1[1] << " " << xg2[1] << std::endl;
                  }



                  adept::adouble divVg = 0.;
                  for(unsigned k = 0; k < dim; k++) {
                    divVg += gradSolV2g[k][k];
                  }

                  for(unsigned i = 0; i < nDofsV2; i++) {
                    for(unsigned k = 0; k < dim; k++) {

                      adept::adouble wlaplace = 0.;
                      adept::adouble advection = 0.;

                      for(unsigned j = 0; j < dim; j++) {
                        wlaplace += muFluid * gradPhiV[i * dim + j] * (gradSolV2g[k][j] + gradSolV2g[j][k]);
                        advection += rhoFluid * solV2gTheta[j] * gradSolV2g[k][j] * phiV[i];

                      }

                      aResV[k][i] += (rhoFluid * (solV2g[k] - solV2gOld[k]) / dt * phiV[i]
                                      + advection
                                      + wlaplace
                                      - gradPhiV[i * dim + k] * solPg) * weightD;
                    }
                  }


                  //continuity block
                  for(unsigned i = 0; i < nDofsP2; i++) {
                    aResP[i] += phiP[i] * divVg * weightD;
                  }



                }//close ig loop

                std::vector < unsigned > solVP2dofs = solV2dofs;
                solVP2dofs.insert(solVP2dofs.end(), solP2dofs.begin(), solP2dofs.end());

                sysDofsAll = solD1dofs;
                sysDofsAll.insert(sysDofsAll.end(), solVP2dofs.begin(), solVP2dofs.end());

                rhs.resize(solVP2dofs.size());   //resize

                for(int i = 0; i < nDofsV2; i++) {
                  for(unsigned  k = 0; k < dim; k++) {
                    rhs[ i +  k * nDofsV2 ] = -aResV[k][i].value();
                  }
                }
                for(int i = 0; i < nDofsP2; i++) {
                  rhs[ i +  dim * nDofsV2 ] = -aResP[i].value();
                }
                myRES->add_vector_blocked(rhs, solVP2dofs);

                Jac.resize(solVP2dofs.size() * sysDofsAll.size());
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
                // std::cout << Jac.size() << " " << solVP2dofs.size() << " " << sysDofsAll.size() << " " << solVP2dofs.size() * sysDofsAll.size() << std::endl << std::flush;

                myKK->add_matrix_blocked(Jac, solVP2dofs, sysDofsAll);

                s.clear_independents();
                s.clear_dependents();
              }
            }
          }
        }


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

  mysolution->_Sol[eflagIndex]->close();
  mysolution->_Sol[pElemIndex]->close();

// *************************************
  end_time = clock();
  AssemblyTime += (end_time - start_time);
// ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

}

