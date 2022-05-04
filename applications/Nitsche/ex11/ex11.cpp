/** tutorial/Ex1
 * This example shows how to:
 * initialize a femus application;
 * define the multilevel-mesh object mlMsh;
 * read from the file ./input/square.neu the coarse-level mesh and associate it to mlMsh;
 * add in mlMsh uniform refined level-meshes;
 * define the multilevel-solution object mlSol associated to mlMsh;
 * add in mlSol different types of finite element solution variables;
 * initialize the solution varables;
 * define vtk and gmv writer objects associated to mlSol;
 * print vtk and gmv binary-format files in ./output directory.
 **/

#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "PetscMatrix.hpp"
#include "LinearImplicitSystem.hpp"
#include "CutFemWeight.hpp"
#include "Fem.hpp"
#include "PetscMatrix.hpp"
#include "slepceps.h"



const double alpha1 = .5;
const double alpha2 = 3.;


#include "GhostPenalty.hpp"

using namespace femus;
void AssembleNitscheProblem_AD(MultiLevelProblem& mlProb);

void BuildFlag(MultiLevelSolution& mlSol, const std::vector<double> &a, const double &d); // the normal a points to Omega1, a[0] x +a[1] y +a[2] z + d = 0

void getNormalInReferenceSystem(const std::vector < std::vector<double> > &xv, const std::vector<double> &aIn, std::vector<double> &aOut);

void GetPlaneInTheParentElement(const std::vector < std::vector<double> > &xv,
                                const std::vector<double> &aIn, const double &dIn, unsigned & efla,
                                std::vector<double> &aOut, double &dOut);





unsigned DIM = 2;

const std::vector<double> N = {-1., 0.33};
const double D = 0.;

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = false; //dirichlet
  // solution does not depend on y and z, only on x. Then dirichlet boundary conditions are fixed only on the left-right side of the domain. Everyelse it is Neumann.
  // u1 has 0 Dirichlet on the left
  // u2 has 0 Dirichlet on the left

  if(DIM == 1 || DIM == 3) {
    if(facename == 1  && !strcmp(SolName, "u1")) dirichlet = true;
    if(facename == 2  && !strcmp(SolName, "u2")) dirichlet = true;
  }
  else if(DIM == 2) {
    if(facename == 4  && !strcmp(SolName, "u1")) dirichlet = true;
    if(facename == 2  && !strcmp(SolName, "u2")) dirichlet = true;
  }
  value = 0.;
  return dirichlet;
}


int main(int argc, char** args) {

  SlepcInitialize(&argc, &args, PETSC_NULL, PETSC_NULL);
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);


  // define multilevel mesh
  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;

  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;

  unsigned nx = 211; // this should always be a odd number

  double lengthX = 3.;
  double length = 1.;

  if(DIM == 1) {
    mlMsh.GenerateCoarseBoxMesh(nx, 0, 0, -lengthX / 2, lengthX / 2, 0., 0., 0., 0., EDGE3, "seventh");
  }
  else if(DIM == 2) {
    mlMsh.GenerateCoarseBoxMesh(nx, 80, 0, -lengthX / 2, lengthX / 2, 0, 1, 0., 0., QUAD9, "seventh");
  }
  else if(DIM == 3) {
    mlMsh.ReadCoarseMesh("./input/cube.neu", "seventh", scalingFactor);
    mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels, NULL);
//     mlMsh.EraseCoarseLevels (numberOfUniformLevels - 1);
//     numberOfUniformLevels = 1;
  }

  //mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , NULL);
  mlMsh.PrintInfo();

  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol(&mlMsh);

  mlSol.AddSolution("u1", LAGRANGE, SECOND);
  mlSol.AddSolution("u2", LAGRANGE, SECOND);

  mlSol.AddSolution("eflag", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
  mlSol.AddSolution("nflag", LAGRANGE, SECOND, 0, false);

  mlSol.Initialize("All");

  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("All");  //??

  BuildFlag(mlSol, N, D);

  MultiLevelProblem ml_prob(&mlSol);

  // ******* Add FEM system to the MultiLevel problem *******
  LinearImplicitSystem& system = ml_prob.add_system < LinearImplicitSystem > ("Nitsche");

  // add solution "u" to system
  system.AddSolutionToSystemPDE("u1");
  system.AddSolutionToSystemPDE("u2");

  // attach the assembling function to system
  system.SetAssembleFunction(AssembleNitscheProblem_AD);

  // time loop parameter
  system.SetMaxNumberOfLinearIterations(1);

  system.SetSparsityPatternMinimumSize(200u);
  
  system.init();

  // ******* Print solution *******
  mlSol.SetWriter(VTK);
  mlSol.GetWriter()->SetDebugOutput(true);

  std::vector<std::string> print_vars;
  print_vars.push_back("All");

  system.MGsolve();

  mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, 0);

  ml_prob.clear();

  return 0;
}

void AssembleNitscheProblem_AD(MultiLevelProblem& ml_prob) {
  //  ml_prob is a big guy. it has multilevel problem, solution, mesh
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
//automatic differentiation//automatic differentiation//automatic differentiation  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled
  //
  // call the adept stack object


  CutFemWeight <double, double> quad  = CutFemWeight<double, double>(QUAD, 7, "legendre");
  Fem fem(quad.GetGaussQuadratureOrder(), quad.GetDimension());

  adept::Stack& s = FemusInit::_adeptStack;//due to automatic differentiation

  //  extract pointers to the several objects that we are going to use

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("Nitsche");   // pointer to the linear implicit system named "Poisson"

  const unsigned level = mlPdeSys->GetLevelToAssemble(); // We have different level of meshes. we assemble the problem on the specified one.

  MultiLevelSolution* mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution *sol  = mlSol->GetSolutionLevel(level); //extracting top level solution where we are going look for solution
  Mesh     *msh   = mlSol->_mlMesh->GetLevel(level); //extracting the mesh relating to that level
  elem *el = msh->el;  // pointer to the elem object in msh (level)
  unsigned iproc  = msh->processor_id();//extracting the process id


  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*            JAC = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*           RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  JAC->zero();
  RES->zero();

  AssembleGhostPenaltyP(ml_prob, true);
  AssembleGhostPenaltyP(ml_prob, false);


  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  //solution variable
  unsigned solu1Index = mlSol->GetIndex("u1");    // get the position of "u" in the ml_sol object
  unsigned solu2Index = mlSol->GetIndex("u2");    // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType(solu1Index);    // get the finite element type for "u"

  unsigned solu1PdeIndex;
  solu1PdeIndex = mlPdeSys->GetSolPdeIndex("u1");    // get the position of "u" in the pdeSys object
  unsigned solu2PdeIndex;
  solu2PdeIndex = mlPdeSys->GetSolPdeIndex("u2");    // get the position of "u" in the pdeSys object

  vector < adept::adouble >  solu1; // local solution
  vector < adept::adouble >  solu2; // local solution

  unsigned eflagIndex = mlSol->GetIndex("eflag");
  unsigned nflagIndex = mlSol->GetIndex("nflag");

  vector < unsigned >  nodeFlag; // local solution

  vector < vector < double > > x(dim);    // local coordinates. x is now dim x m matrix.
  unsigned biquadraticType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  double weight; // gauss point weight

  vector< adept::adouble > aResu1; // local residual vector
  vector< adept::adouble > aResu2; // local residual vector

  vector< unsigned > l2GMap; // local to global mapping
  vector< double > Res; // local redidual vector
  vector < double > Jac;

//   JAC->zero(); // Set to zero all the entries of the Global Matrix
//   RES->zero(); // Set to zero all the entries of the Global Residual

  std::vector < std::vector < std::vector <double > > > aP(3);

  std::vector < double > Av;
  std::vector < std::vector < double > > Bv(2);

  Mat A, B;
  EPS eps;


  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);  // number of solution element dofs


    unsigned eFlag = static_cast <unsigned>(floor((*sol->_Sol[eflagIndex])(iel) + 0.5));

    // resize local arrays
    l2GMap.resize(2 * nDofu); //at the end we have both u1 and u2. therefore size should be double
    solu1.resize(nDofu);
    solu2.resize(nDofu);
    nodeFlag.resize(nDofu);

    for(int k = 0; k < dim; k++) {
      x[k].resize(nDofu);
    }

    aResu1.assign(nDofu, 0.);    //assign
    aResu2.assign(nDofu, 0.);    //assign

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofu; i++) {
      unsigned iDof = msh->GetSolutionDof(i, iel, soluType);    // local to global mapping between solution node and solution dof
      solu1[i] = (*sol->_Sol[solu1Index])(iDof);                  // global extraction and local storage for the solution
      solu2[i] = (*sol->_Sol[solu2Index])(iDof);                  // global extraction and local storage for the solution
      nodeFlag[i] = (*sol->_Sol[nflagIndex])(iDof);

      l2GMap[i] = pdeSys->GetSystemDof(solu1Index, solu1PdeIndex, i, iel);    // local to global mapping between u1 and pdeSys dof
      l2GMap[nDofu + i] = pdeSys->GetSystemDof(solu2Index, solu2PdeIndex, i, iel);    // local to global mapping between u2 and pdeSys dof

    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofu; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, biquadraticType);    // global to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k])(xDof);      // global extraction and local storage for the element coordinates
      }
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();


    if(eFlag == 0 || eFlag == 2) {
      // *** Element Gauss point loop ***
      for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
        // *** get gauss point weight, test function and test function partial derivatives (inputs: coordinates and gause point) ***
        msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x);

        // evaluate the solution, the solution derivatives and the coordinates in the gauss point
        //why gradient of solution adouble?, because solution is adouble variable
        vector < adept::adouble > gradSolu1g(dim, 0.);
        vector < adept::adouble > gradSolu2g(dim, 0.);
        //gradient of the solution at gauss point // 'k' is for dimension and 'i' is for nodes of the mesh
        for(unsigned i = 0; i < nDofu; i++) {
          for(unsigned k = 0; k < dim; k++) {
            gradSolu1g[k] += phi_x[i * dim + k] * solu1[i];
            gradSolu2g[k] += phi_x[i * dim + k] * solu2[i];
          }
        }

        // *** phi_i loop ***
        for(unsigned i = 0; i < nDofu; i++) {

          adept::adouble graduGradphi1 = 0.;
          adept::adouble graduGradphi2 = 0.;
          //right hand side a(u,v)
          for(unsigned k = 0; k < dim; k++) {
            graduGradphi1 += gradSolu1g[k] * phi_x[i * dim + k];
            graduGradphi2 += gradSolu2g[k] * phi_x[i * dim + k];
          }


          if(eFlag == 0) { //we are on the left of the domain
            aResu1[i] += (- phi[i] + alpha1 * graduGradphi1) * weight;
            if(nodeFlag[i] != 1) {
              aResu2[i] += (graduGradphi2) * weight;//fake equations for u2 in Omega1 open
            }
          }
          else if(eFlag == 2) {
            aResu2[i] += (- phi[i] + alpha2 * graduGradphi2) * weight;
            if(nodeFlag[i] != 1) {
              aResu1[i] += (graduGradphi1) * weight; //fake equations for u1 in Omega2 open
            }
          }

        } // end phi_i loop
      } // end gauss point loop
    }

    else {



      std::vector<double> a;
      double d = 0;
      unsigned elemFlag;

      GetPlaneInTheParentElement(x, N, D, elemFlag, a, d);

      //getNormalInReferenceSystem(x, N, a);

     // std::cout << iel << " " << a[0] << " " << a[1] << " " << d << std::endl;

      //plane passsing through center *************

      std::vector<double> eqPoly1;
      quad.GetWeightWithMap(0, a, d, eqPoly1);//s=0 (volume integral)

      std::vector<double> eqPoly2;
      quad.GetWeightWithMap(0, {-a[0], -a[1]}, -d, eqPoly2);

      std::vector<double> eqPolyI;
      quad.GetWeightWithMap(-1, a, d, eqPolyI);//s=-1 (boundary integral)

      const elem_type *thisfem = fem.GetFiniteElement(ielGeom, soluType);

//       //solve eigenvalu problems to get theta, gamma1, gamma2
//       Av.assign(nDofu * nDofu, 0.);
//       Bv[0].assign(nDofu * nDofu, 0.);
//       Bv[1].assign(nDofu * nDofu, 0.);
// 
//       for(unsigned ig = 0; ig < thisfem->GetGaussPointNumber(); ig++) {
//         // *** get gauss point weight, test function and test function partial derivativese (inputs: coordinates and gause point) ***
// 
//         std::vector<std::vector<double>> Jac, JacI;
//         thisfem->GetJacobianMatrix(x, ig, weight, Jac, JacI);
// 
//         thisfem->Jacobian(x, ig, weight, phi, phi_x);
// 
//         double dsN = 0;
//         for(unsigned i = 0; i < dim; i++) {
//           double dsi = 0.;
//           for(unsigned j = 0; j < dim; j++) {
//             dsi += JacI[j][i] * a[j];
//           }
//           dsN += dsi * dsi;
//         }
//         dsN = sqrt(dsN);
// 
// 
//         for(int i = 0; i < nDofu; i++) {
//           for(int j = 0; j < nDofu; j++) {
//             for(unsigned k = 0; k < dim; k++) {
//               Bv[0][ i * nDofu + j] += phi_x[i * dim + k] * phi_x[j * dim + k] * weight * eqPoly1[ig];
//               Bv[1][ i * nDofu + j] += phi_x[i * dim + k] * phi_x[j * dim + k] * weight * eqPoly2[ig];
//             }
//           }
// 
//           double gradPhiiDotN = 0.;
//           for(unsigned k = 0; k < dim; k++) {
//             gradPhiiDotN += phi_x[i * dim + k] * N[k];
//           }
//           for(int j = 0; j < nDofu; j++) {
//             double gradPhijDotN = 0.;
//             for(unsigned k = 0; k < dim; k++) {
//               gradPhijDotN += phi_x[j * dim + k] * N[k];
//             }
//             Av[ i * nDofu + j] += gradPhiiDotN * gradPhijDotN  * weight * eqPolyI[ig] * dsN;
//           }
//         }
//       }
// 
//       /* Careful B has one zero eigenvalue (it is singular) with nullspace x = [1,1,1,...]^T,
//       * thus the generalized eigenvalue problem
//       * $$A u = \lambda B u$$
//       * (or $B^{-1} A x = \lambda x$) has one indetermined eigenvalue, that makes the SLEPC solve very unstable,
//       * even using its built-in deflation method.
//       * Fortunately, A has at least one zero eigenvalue, with the same x = [1,1,1,...]^T being an element
//       * of its nullspace. Then, it is possible to deflate A and B simultaneously:
//       * $Ad = A - x^T. a1$ and $Bd = B - x^T.b1$, where a1 and b1 are the first rows of A and B, respectively.
//       * The generalized eigenvalue problem $Ab u = \lambda Bb u$, with matrices Ab and Bb,
//       * obtained as block matrices from Ad and Bd, removing the first row and the first column,
//       * has the same eigenvalues of the original one, except the indetermine one.
//       * Note that Bb is now invertible, and SLEPC has no problem in solving the deflated
//       * generalized eigenvalue problem.
//       */
// 
// 
//       double C[2];
//       for(unsigned k = 0; k < 2; k++) {
//         MatCreateSeqDense(PETSC_COMM_SELF, nDofu - 1, nDofu - 1, NULL, &A);
//         MatCreateSeqDense(PETSC_COMM_SELF, nDofu - 1, nDofu - 1, NULL, &B);
// 
//         for(int i = 0; i < nDofu - 1; i++) {
//           for(int j = 0; j < nDofu - 1; j++) {
//             double value;
//             value = Bv[k][(i + 1) * nDofu + (j + 1)] - Bv[k][j + 1];
//             MatSetValues(B, 1, &i, 1, &j, &value, INSERT_VALUES);
//             value = Av[(i + 1) * nDofu + (j + 1)] - Av[j + 1];
//             MatSetValues(A, 1, &i, 1, &j, &value, INSERT_VALUES);
//           }
//         }
// 
//         MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
//         MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
// 
//         MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
//         MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
// 
//         EPSCreate(PETSC_COMM_SELF, &eps);
//         EPSSetOperators(eps, A, B);
//         EPSSetFromOptions(eps);
//         EPSSetWhichEigenpairs(eps, EPS_LARGEST_MAGNITUDE);
//         EPSSolve(eps);
// 
//         double imaginary;
//         EPSGetEigenpair(eps, 0, &C[k], &imaginary, NULL, NULL);
// 
//         EPSDestroy(&eps);
//         MatDestroy(&A);
//         MatDestroy(&B);
//       }
// 
// 
//       double ia1C1 = 1. / (alpha1 * C[0]);
//       double ia2C2 = 1. / (alpha2 * C[1]);
// 
//       double den = ia1C1 + ia2C2;
// 
//       double gamma1 = ia1C1 / den;
//       double gamma2 = ia2C2 / den;
// 
//       double theta = 2. / den;


      double h = sqrt(dim) * sqrt((x[0][0] - x[0][1]) * (x[0][0] - x[0][1]) +
                                  (x[1][0] - x[1][1]) * (x[1][0] - x[1][1])) ;
      double gamma1 = 0.5;
      double gamma2 = 0.5;
      double theta = 10. * 0.5 * (alpha1 + alpha2) / h;

      for(unsigned ig = 0; ig < thisfem->GetGaussPointNumber(); ig++) {
        // *** get gauss point weight, test function and test function partial derivativese (inputs: coordinates and gause point) ***

        std::vector<std::vector<double>> Jac, JacI;
        thisfem->GetJacobianMatrix(x, ig, weight, Jac, JacI);

        thisfem->Jacobian(x, ig, weight, phi, phi_x);

        double dsN = 0;
        for(unsigned i = 0; i < dim; i++) {
          double dsi = 0.;
          for(unsigned j = 0; j < dim; j++) {
            dsi += JacI[j][i] * a[j];
          }
          dsN += dsi * dsi;
        }
        dsN = sqrt(dsN);

        vector < adept::adouble > gradSolu1g(dim, 0.);
        vector < adept::adouble > gradSolu2g(dim, 0.);

        for(unsigned i = 0; i < nDofu; i++) {
          for(unsigned k = 0; k < dim; k++) {
            gradSolu1g[k] += phi_x[i * dim + k] * solu1[i];
            gradSolu2g[k] += phi_x[i * dim + k] * solu2[i];
          }
        }

        adept::adouble solu1g  = 0.;
        adept::adouble solu2g  = 0.;
        adept::adouble alphaGradSoluDotN = 0.;


        for(unsigned i = 0; i < nDofu; i++) {
          solu1g += phi[i] * solu1[i];
          solu2g += phi[i] * solu2[i];
          for(unsigned k = 0; k < dim; k++) {
            alphaGradSoluDotN += (alpha1 * gamma1  * solu1[i] + alpha2 * gamma2 * solu2[i]) * phi_x[i * dim + k] * N[k];
          }
        }

        // *** phi_i loop ***
        for(unsigned i = 0; i < nDofu; i++) {

          adept::adouble gradu1Gradphi = 0.;
          adept::adouble gradu2Gradphi = 0.;
          //right hand side a(u,v)
          for(unsigned k = 0; k < dim; k++) {
            gradu1Gradphi += gradSolu1g[k] * phi_x[i * dim + k];
            gradu2Gradphi += gradSolu2g[k] * phi_x[i * dim + k];
          }
          aResu1[i] += (- phi[i] + alpha1 * gradu1Gradphi) * weight * eqPoly1[ig];
          aResu2[i] += (- phi[i] + alpha2 * gradu2Gradphi) * weight * eqPoly2[ig];

          adept::adouble gradPhiDotN = 0.;
          for(unsigned k = 0; k < dim; k++) {
            gradPhiDotN += phi_x[i * dim + k] * N[k];
          }

          aResu1[i] += (solu2g - solu1g) * (-phi[i] * theta - alpha1 * gamma1 * gradPhiDotN) * weight * eqPolyI[ig] * dsN;
          aResu1[i] += alphaGradSoluDotN * (+phi[i]) * weight * eqPolyI[ig] * dsN;
          aResu2[i] += (solu2g - solu1g) * (+phi[i] * theta - alpha2 * gamma2 * gradPhiDotN) * weight * eqPolyI[ig] * dsN;
          aResu2[i] += alphaGradSoluDotN * (-phi[i]) * weight * eqPolyI[ig] * dsN;


        } // end phi_i loop
      } // end gauss point loop
    }

    //copy the value of the adept::adoube aRes in double Res and store
    Res.resize(2 * nDofu);    //resize

    for(int i = 0; i < nDofu; i++) {
      Res[i] = - aResu1[i].value();
      Res[nDofu + i] = - aResu2[i].value();
    }

    RES->add_vector_blocked(Res, l2GMap);

    // define the dependent variables
    s.dependent(&aResu1[0], nDofu);
    s.dependent(&aResu2[0], nDofu);

    // define the independent variables
    s.independent(&solu1[0], nDofu);
    s.independent(&solu2[0], nDofu);

    // get the jacobian matrix (ordered by row major )
    Jac.resize(2 * nDofu * 2 * nDofu);    //resize
    s.jacobian(&Jac[0], true);

    //store K in the global matrix KK
    JAC->add_matrix_blocked(Jac, l2GMap, l2GMap);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();
  JAC->close();

//   PetscViewer    viewer;
//   //PetscViewerDrawOpen (PETSC_COMM_WORLD, NULL, NULL, 0, 0, 900, 900, &viewer);
//   //PetscObjectSetName ( (PetscObject) viewer, "FSI matrix");
//   //PetscViewerPushFormat (viewer, PETSC_VIEWER_DRAW_LG);
//   MatView ( (static_cast<PetscMatrix*> (KK))->mat(), viewer);
//   VecView ( (static_cast<PetscVector*> (RES))->vec(), viewer);

  //double a;
  //std::cin >> a;



  //***************** END ASSEMBLY *******************
}

void BuildFlag(MultiLevelSolution & mlSol, const std::vector<double> &a, const double & d) {

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1; //why -1? c/c++ notation
  //object we have is multilevel Solution. we need to extract highest level of the solution.

  Solution *sol  = mlSol.GetSolutionLevel(level); //extracting top level solution where we are going look for solution
  Mesh     *msh   = mlSol._mlMesh->GetLevel(level); //extracting the mesh relating to that level
  unsigned iproc  = msh->processor_id();//extracting the process id

  unsigned eflagIndex = mlSol.GetIndex("eflag"); //identifier where eflag is stored
  unsigned nflagIndex = mlSol.GetIndex("nflag");

  unsigned nflagType = mlSol.GetSolutionType(nflagIndex); //giving us the type order of finite element (lagrange second)
  //we already know eflag type as it is piecewise constant

  unsigned biquadraticType = 2; //2 is the indicator of biquadratic type

  sol->_Sol[eflagIndex]->zero();//setting elagindex to zero
  sol->_Sol[nflagIndex]->zero();

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  std::vector < std::vector<double> >  x(dim); //standard vector which will store coordinates of mesh
  //each element has several nodes. each nodes has dim number of coodinates.

  //loop on the element
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel); //since we are paticular element we can extract the geomery of the element. it can be square, triangle..1D 2D
    unsigned nDofs  = msh->GetElementDofNumber(iel, nflagType);  // number of solution element dofs //number of ndofs depend on two parameters, element you are looking at and type of solution

    for(int k = 0; k < dim; k++) {
      x[k].resize(nDofs);  // rezising to storing coodinates related to nDofs
    }
    //extracting coodinates of the ndofs
    double pcnt = 0;
    double mcnt = 0;
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, biquadraticType);    // local to global mapping between coordinates node and coordinate dofs
      double dist = d;
      for(unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k])(xDof);
        dist += a[k] * x[k][i];

      }

      if(dist > 0) pcnt++;
      else if(dist < 0) mcnt++;
    }

    if(pcnt == 0) {
      sol->_Sol[eflagIndex]->set(iel, 2.);
    }
    else if(mcnt != 0) {
      sol->_Sol[eflagIndex]->set(iel, 1.);
      for(unsigned i = 0; i < nDofs; i++) {
        unsigned iDof = msh->GetSolutionDof(i, iel, nflagType);
        sol->_Sol[nflagIndex]->set(iDof, 1.);
      }
    }
  }
  sol->_Sol[eflagIndex]->close();
  sol->_Sol[nflagIndex]->close();
}


void getNormalInReferenceSystem(const std::vector < std::vector<double> > &xv, const std::vector<double> &aIn, std::vector<double> &aOut) {
//aIn is a input
  unsigned dim = 2;

  const double& x1 = xv[0][0];
  const double& x2 = xv[0][1];
  const double& x3 = xv[0][2];
  const double& x4 = xv[0][3];
  const double& y1 = xv[1][0];
  const double& y2 = xv[1][1];
  const double& y3 = xv[1][2];
  const double& y4 = xv[1][3];

  double J[2][2];

  double dx12 = x1 - x2;
  double dx34 = x3 - x4;
  double dy12 = y1 - y2;
  double dy34 = y3 - y4;

  double dx14 = (x1 - x4);
  double dy23 = (y2 - y3);
  double dx23 = (x2 - x3);
  double dy14 = (y1 - y4);

  double u = 0;
  double v = 0;

  J[0][0] = 0.25 * ((-1. + v) * dx12 + (1. + v) * dx34);
  J[0][1] = 0.25 * ((-1. + u) * dx14 - (1. + u) * dx23);
  J[1][0] = 0.25 * ((-1. + v) * dy12 + (1. + v) * dy34);
  J[1][1] = 0.25 * ((-1. + u) * dy14 - (1. + u) * dy23);


  aOut.assign(dim, 0);
  for(unsigned k = 0; k < dim; k++) {
    for(unsigned j = 0; j < dim; j++) {
      aOut[k] += J[j][k] * aIn[j]; ///aIn=physical normal , aOut=fem normal/normal in the parent system, J[j][k] is jacobian transpose
    }
  }
  double aOutNorm = sqrt(aOut[0] * aOut[0] + aOut[1] * aOut[1]);
  aOut[0] /= aOutNorm;
  aOut[1] /= aOutNorm;
}

void GetPlaneInTheParentElement(const std::vector < std::vector<double> > &xv,
                                const std::vector<double> &aIn, const double & dIn, unsigned & efla,
                                std::vector<double> &aOut, double & dOut) {

  unsigned dim = 2;

  const double& x1 = xv[0][0];
  const double& x2 = xv[0][1];
  const double& x3 = xv[0][2];
  const double& x4 = xv[0][3];
  const double& y1 = xv[1][0];
  const double& y2 = xv[1][1];
  const double& y3 = xv[1][2];
  const double& y4 = xv[1][3];

  std::vector<double> xm;

  //setting characteristic epsilon
  double h1 = (fabs(x2 - x1) + fabs(x3 - x2) + fabs(x4 - x3) + fabs(x4 - x1)) / 4.;
  double h2 = (fabs(y2 - y1) + fabs(y3 - y2) + fabs(y4 - y3) + fabs(y4 - y1)) / 4.;

  double h = sqrt(h1 * h1 + h2 * h2);
  double eps = 1.0e-10 * h;

  //extracting number of vertices (4)
  unsigned l = 4;

  //a vector to store distances
  std::vector<double> dist(l, dIn);
  std::vector<double> distf(l);

  double den = 0.;
  for(unsigned k = 0 ; k < dim ; k++) {
    den += aIn[k] * aIn[k];
  }
  den = 1. / sqrt(den);

  unsigned t = 0;
  //calculating perpendicular distances
  for(unsigned i = 0 ; i < l ; i++) {
    for(unsigned k = 0 ; k < dim ; k++) {
      dist[i] += aIn[k] * xv[k][i];
    }
    dist[i] = dist[i] / sqrt(den);

    ////taking care of extream cases (updating dist and filling distf)
    if(fabs(dist[i]) < eps) {  //interface is very close to a vertex
      distf[i] = (dist[i] < 0.) ? -eps : eps;
      dist[i] = 0. ;
      t++;
    }
    else {
      distf[i] = dist[i];
    }
  }

  //if we have a interface very close to vertex
  if(t > 0) {
    unsigned pcnt = 0;
    double mcnt = 0;
    for(unsigned i = 0; i < l; i++) {
      if(dist[i] > 0.) {
        pcnt++;
      }
      else if(dist[i] < 0.) {
        mcnt++;
      }
      dist[i] = distf[i];
    }

    if(pcnt == 0) {
      efla = 2;
      return;
    }
    else if(mcnt == 0) {
      efla = 0;
      return;
    }
  }

  //calculating s
  std::vector < std::vector <double> > xe(dim, std::vector<double>(2, 0));
  unsigned j = 0;
  for(unsigned i = 0; i < l; i++) {
    unsigned i1 = (i + 1) % l;
    if(dist[i] * dist[i1] < 0) {
      double s = fabs(dist[i]) / (fabs(dist[i]) + fabs(dist[i1]));
      for(unsigned k = 0; k < dim; k++) {
        xe[k][j] = (1 - s) * xv[k][i] + s * xv[k][i1];
      }
      j++;
    }
  }

  if(j == 0) {
    efla = (dist[0] < 0) ? 2 : 0;
    return;
  }
  else {
    efla = 1;
    xm.resize(dim);
    for(unsigned k = 0; k < dim; k++) {
      xm[k] = (xe[k][0] + xe[k][1]) / 2.;
    }

    std::vector<double> xi(dim);
    double &u = xi[0];
    double &v = xi[1];

    std::vector < std::vector < double > > J(2, std::vector<double>(2));

    double dx12 = x1 - x2;
    double dx34 = x3 - x4;
    double dy12 = y1 - y2;
    double dy34 = y3 - y4;
    double hu = dx34 * dy12 - dx12 * dy34;

    double dx14 = (x1 - x4);
    double dy23 = (y2 - y3);
    double dx23 = (x2 - x3);
    double dy14 = (y1 - y4);
    double hv = dx14 * dy23 - dx23 * dy14;

    double eps2 = 1.0e-10 * h * h;

    if(fabs(hu) > eps2) {//edges 1 and 3 are not parallel
      double gu = -x4 * y1 + x3 * y2 - x2 * y3 + x1 * y4;
      double f = xm[0] * (dy12 + dy34) - xm[1] * (dx12 + dx34);
      double fpgu = f + gu;

      double det = sqrt(hu * (- 2. * xm[0] * (dy14 + dy23)
                              + (2. * xm[1] - y3 - y4) * (x1 + x2)
                              - (2. * xm[1] - y1 - y2) * (x3 + x4))
                        + fpgu * fpgu);
      u = (fpgu + det) / hu;

      if(fabs(hv) > eps2) { //edges 2 and 4 are not parallel
        double gv = -x4 * y3 + x3 * y4 - x2 * y1 + x1 * y2;
        v = (f + gv - det) / hv;

        J[0][0] = 0.25 * ((-1. + v) * dx12 + (1. + v) * dx34);
        J[0][1] = 0.25 * ((-1. + u) * dx14 - (1. + u) * dx23);
        J[1][0] = 0.25 * ((-1. + v) * dy12 + (1. + v) * dy34);
        J[1][1] = 0.25 * ((-1. + u) * dy14 - (1. + u) * dy23);

      }
      else { //edges 2 and 4 are parallel
        J[0][1] = 0.25 * ((-1. + u) * dx14 - (1. + u) * dx23);
        J[1][1] = 0.25 * ((-1. + u) * dy14 - (1. + u) * dy23);

        v = (J[0][1] > eps) ?
            (0.25 * ((-1. + u) * (x1 + x4) - (1. + u) * (x3 + x2)) + xm[0]) / J[0][1] :
            (0.25 * ((-1. + u) * (y1 + y4) - (1. + u) * (y3 + y2)) + xm[1]) / J[1][1];

        J[0][0] = 0.25 * ((-1. + v) * dx12 + (1. + v) * dx34);
        J[1][0] = 0.25 * ((-1. + v) * dy12 + (1. + v) * dy34);

      }
    }
    else if(fabs(hv) > eps2) {  //edges 1 and 3 are parallel, but edges 2 and 4 are not
      // std::cout << "1 and 3 are parallel\n";
      double f = xm[0] * (dy12 + dy34) - xm[1] * (dx12 + dx34);
      double gv = -x4 * y3 + x3 * y4 - x2 * y1 + x1 * y2;
      double fpgv = f + gv;

      double det = sqrt(hv * (- 2. * xm[0] * (dy12 - dy34)
                              + (2. * xm[1] - y2 - y3) * (x1 + x4)
                              - (2. * xm[1] - y1 - y4) * (x2 + x3))
                        +  fpgv * fpgv);

      v = (fpgv - det) / hv;

      J[0][0] = 0.25 * ((-1. + v) * dx12 + (1. + v) * dx34);
      J[1][0] = 0.25 * ((-1. + v) * dy12 + (1. + v) * dy34);

      u = (fabs(J[0][0]) > eps) ?
          (0.25 * ((-1. + v) * (x1 + x2) - (1. + v) * (x3 + x4)) + xm[0]) / J[0][0] :
          (0.25 * ((-1. + v) * (y1 + y2) - (1. + v) * (y3 + y4)) + xm[1]) / J[1][0];

      J[0][1] = 0.25 * ((-1. + u) * dx14 - (1. + u) * dx23);
      J[1][1] = 0.25 * ((-1. + u) * dy14 - (1. + u) * dy23);
    }
    else { //edges 1 and 3, and  edges 2 and 4 are parallel
      //   std::cout << "Romboid\n";
      std::vector<std::vector<unsigned> > idx = {{3, 1}, {0, 2}};

      double A[2][2] = {{-dy14, dy23}, {dy12, dy34}};
      double B[2][2] = {{dx14, -dx23}, {-dx12, -dx34}};

      for(unsigned k = 0; k < 2; k++) {
        double d[2];
        for(unsigned j = 0 ; j < 2; j++) {
          double Ckj = - A[k][j] * xv[0][idx[k][j]] - B[k][j] * xv[1][idx[k][j]];
          d[j] = (A[k][j] * xm[0] + B[k][j] * xm[1] + Ckj) / sqrt(A[k][j] * A[k][j] + B[k][j] * B[k][j]);
        }
        xi[k] = -1. + 2. * d[0] / (d[0] + d[1]);
      }

      J[0][0] = 0.25 * ((-1. + v) * dx12 + (1. + v) * dx34);
      J[0][1] = 0.25 * ((-1. + u) * dx14 - (1. + u) * dx23);
      J[1][0] = 0.25 * ((-1. + v) * dy12 + (1. + v) * dy34);
      J[1][1] = 0.25 * ((-1. + u) * dy14 - (1. + u) * dy23);
    }

    aOut.assign(dim, 0);
    for(unsigned k = 0; k < dim; k++) {
      for(unsigned j = 0; j < dim; j++) {
        aOut[k] += J[j][k] * aIn[j];
      }
    }
    double aOutNorm = sqrt(aOut[0] * aOut[0] + aOut[1] * aOut[1]);
    aOut[0] /= aOutNorm;
    aOut[1] /= aOutNorm;
    dOut = - aOut[0] * u - aOut[1] * v;

  }



}
