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
  vector< vector< adept::adouble > > solV(dim);      // local solution (velocity)
  vector< adept::adouble > solP;

  vector< vector< double > > solDOld(dim);      // local solution (displacement)
  vector< vector< double > > solVOld(dim);
  vector< vector< double > > solAOld(dim);

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

  vector < double> nablaphiP;
  vector < double> nablaphiV;

  unsigned dim2 = 3 * (dim - 1);

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
          vx[k][i] = (*msh->_topology->_Sol[k])(idofX) + (1. - af) * solDOld[k][i] + af * solD[k][i];
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
              gradSolDgHat[k][j] += ((1. - af) * solD[k][i] + af * solDOld[k][i]) * gradPhiDHat[i * dim + j];
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
              aRhsD[k][i] -=   CauchyDIR[k] * weightD;
            }
          }
        }

//         unsigned group = msh->GetElementGroup(iel);
//         double E = (group == 13) ? 2 : 1;
//         double nu = 0.35;
//
//         double mu = E / (2. * (1. + nu));
//         double lambda = (E * nu) / ((1. + nu) * (1. - 2.*nu));
//
//         adept::adouble divD = 0.;
//         for(int I = 0; I < dim; I++) {
//           divD += gradSolDgHat[I][I];
//         }
//
//         // *** phiA_i loop ***
//         for(unsigned i = 0; i < nDofs; i++) {
//           if(nodeFlag[i] == 4) {
//             for(unsigned  I = 0; I < dim; I++) {  //momentum equation in k
//               adept::adouble term = 0.;
//               for(unsigned J = 0; J < dim; J++) {  // second index j in each equation
//                 term +=  mu * gradPhiDHat[i * dim + J] * (gradSolDgHat[I][J] + gradSolDgHat[J][I]); // diffusion
//               }
//               term +=  gradPhiDHat[i * dim + I]  * lambda * divD;
//
//               aRhsD[I][i] -=   term * weightDHat;
//             }
//           }
//         } // end phiA_i loop




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


