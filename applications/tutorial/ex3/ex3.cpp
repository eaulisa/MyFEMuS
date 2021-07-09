/** tutorial/Ex3
 * This example shows how to set and solve the weak form of the multidimensional Chafee–Infante steady state nonlinear equation
 *                     -\Delta u + lambda (u^3-u) = f(x) \text{ on }\Omega,
 *            u=0 \text{ on } \Gamma,
 * on a box domain $\Omega$ with boundary $\Gamma$;
 * all the coarse-level meshes are removed;
 * a multilevel problem and an equation system are initialized;
 * a direct solver is used to solve the problem.
 **/

#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "adept.h"


using namespace femus;

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0;

  if(facename == 2)
    dirichlet = false;

  return dirichlet;
}

void AssembleNonlinearProblem(MultiLevelProblem& ml_prob);

void AssembleNonlinearProblem_AD(MultiLevelProblem& ml_prob);

std::pair < double, double > GetErrorNorm(MultiLevelSolution* mlSol);

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);


  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  mlMsh.ReadCoarseMesh("./input/square_quad.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/square_tri.neu","seventh",scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/square_mixed.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/cube_hex.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/cube_wedge.neu","seventh",scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/cube_tet.neu","seventh",scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/cube_mixed.neu","seventh",scalingFactor);

  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
    probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();
  unsigned maxNumberOfMeshes;

  if(dim == 2) {
    maxNumberOfMeshes = 5;
  }
  else {
    maxNumberOfMeshes = 4;
  }

  std::vector < std::vector < double > > l2Norm;
  l2Norm.resize(maxNumberOfMeshes);

  std::vector < std::vector < double > > semiNorm;
  semiNorm.resize(maxNumberOfMeshes);

  for(unsigned i = 1; i < maxNumberOfMeshes; i++) {    // loop on the mesh level

    unsigned numberOfUniformLevels = i ;
    unsigned numberOfSelectiveLevels = 0;
    mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

    // erase all the coarse mesh levels
    mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

    // print mesh info
    mlMsh.PrintInfo();

    FEOrder feOrder[3] = {FIRST, SERENDIPITY, SECOND};
    l2Norm[i].resize(3);
    semiNorm[i].resize(3);

    for(unsigned j = 0; j < 3; j++) {    // loop on the FE Order
      // define the multilevel solution and attach the mlMsh object to it
      MultiLevelSolution mlSol(&mlMsh);

      // add variables to mlSol
      mlSol.AddSolution("u", LAGRANGE, feOrder[j]);
      mlSol.Initialize("All");

      // attach the boundary condition function and generate boundary data
      mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
      mlSol.GenerateBdc("u");

      // define the multilevel problem attach the mlSol object to it
      MultiLevelProblem mlProb(&mlSol);

      // add system Chafee–Infante in mlProb as a Linear Implicit System
      NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("Chafee–Infante");

      // add solution "u" to system
      system.AddSolutionToSystemPDE("u");

      // attach the assembling function to system
      //system.SetAssembleFunction(AssembleNonlinearProblem_AD);
      system.SetAssembleFunction(AssembleNonlinearProblem);

      // initilaize and solve the system
      system.init();
      system.SetOuterSolver(PREONLY);
      system.MGsolve();

      std::pair< double , double > norm = GetErrorNorm(&mlSol);
      l2Norm[i][j]  = norm.first;
      semiNorm[i][j] = norm.second;
      // print solutions
      std::vector < std::string > variablesToBePrinted;
      variablesToBePrinted.push_back("All");

      VTKWriter vtkIO(&mlSol);
      vtkIO.SetDebugOutput(true);
      vtkIO.Write(DEFAULT_OUTPUTDIR, "linear", variablesToBePrinted, i + j * 10);
      vtkIO.Write(DEFAULT_OUTPUTDIR, "quadratic", variablesToBePrinted, i + j * 10);
      vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, i + j * 10);

    }
  }

  // print the seminorm of the error and the order of convergence between different levels
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "l2 ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";

  for(unsigned i = 1; i < maxNumberOfMeshes; i++) {
    std::cout << i + 1 << "\t";
    std::cout.precision(14);

    for(unsigned j = 0; j < 3; j++) {
      std::cout << l2Norm[i][j] << "\t";
    }

    std::cout << std::endl;

    if(i < maxNumberOfMeshes - 1) {
      std::cout.precision(3);
      std::cout << "\t\t";

      for(unsigned j = 0; j < 3; j++) {
        std::cout << log(l2Norm[i][j] / l2Norm[i + 1][j]) / log(2.) << "\t\t\t";
      }

      std::cout << std::endl;
    }

  }

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "SEMINORM ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";

  for(unsigned i = 1; i < maxNumberOfMeshes; i++) {
    std::cout << i + 1 << "\t";
    std::cout.precision(14);

    for(unsigned j = 0; j < 3; j++) {
      std::cout << semiNorm[i][j] << "\t";
    }

    std::cout << std::endl;

    if(i < maxNumberOfMeshes - 1) {
      std::cout.precision(3);
      std::cout << "\t\t";

      for(unsigned j = 0; j < 3; j++) {
        std::cout << log(semiNorm[i][j] / semiNorm[i + 1][j]) / log(2.) << "\t\t\t";
      }

      std::cout << std::endl;
    }

  }
}


double GetExactSolutionValue(const std::vector < double >& x) {
  double pi = acos(-1.);
  return cos(pi * x[0]) * cos(pi * x[1]);
};

void GetExactSolutionGradient(const std::vector < double >& x, std::vector < double >& solGrad) {
  double pi = acos(-1.);
  solGrad[0]  = -pi * sin(pi * x[0]) * cos(pi * x[1]);
  solGrad[1] = -pi * cos(pi * x[0]) * sin(pi * x[1]);
};

double GetExactSolutionLaplace(const std::vector < double >& x) {
  double pi = acos(-1.);
  return -2.*pi * pi * cos(pi * x[0]) * cos(pi * x[1]); 
};

/**
 * Given the non linear problem
 *
 *        - \Delta u + lambda ( u^3 -u ) = f(x),
 *
 * in the unit box centered in the origin with
 *
 *        f(x) = - \Delta u_e + lambda ( u_e^3 - u_e),
 *        u_e = \cos ( \pi * x ) * \cos( \pi * y ),
 *
 *the following function assembles the residual std::vector Res(u^i) and the exact Jacobian matrix J(u^i)
 * for the Newton iteration, i.e.
 *
 *        J(u^i) w = Res(u^i) = f(x) - ( - \Delta u^i + lambda ( u_i^3 - u^i ),
 *        u^{i+1} = u^i + w,
 *        where
 *        J(u^i) w = - \Delta w  + lambda ( 3 u_i ^2 w - w ).
 *
 **/

void AssembleNonlinearProblem(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system< NonLinearImplicitSystem > ("Chafee–Infante");   // pointer to the linear implicit system named "Chafee–Infante"

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    K         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  const unsigned lambda = 1;
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex;
  soluIndex = mlSol->GetIndex("u");    // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex("u");    // get the position of "u" in the pdeSys object

  std::vector < double >  solu; // local solution

  std::vector < std::vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector < unsigned > sysDof; // local to global pdeSys dofs
  std::vector < double > phi;  // local test function
  std::vector < double > phi_x; // local test function first order partial derivatives
  double weight; // gauss point weight

  std::vector < double > Res; // local residual vector
  std::vector < double > J; // local Jacobian matrix

  K->zero(); // Set to zero all the entries of the Global Matrix
  RES->zero(); // Set to zero all the entries of the Global Residual

  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofs  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs

    // resize local arrays
    sysDof.resize(nDofs);
    solu.resize(nDofs);

    for(unsigned k = 0; k < dim; k++) {
      x[k].resize(nDofs);
    }

    Res.assign(nDofs, 0.);
    J.assign(nDofs * nDofs, 0.);

    // local storage of global mapping
    for(unsigned i = 0; i < nDofs; i++) {
      //solution
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
      sysDof[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof

      // coordinates
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k])(xDof);      // global extraction and local storage for the element coordinates
      }
    }

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      double soluGauss = 0;
      std::vector < double > soluGauss_x(dim, 0.);
      std::vector < double > xGauss(dim, 0.);

      for(unsigned i = 0; i < nDofs; i++) {
        soluGauss += phi[i] * solu[i];

        for(unsigned k = 0; k < dim; k++) {
          soluGauss_x[k] += phi_x[i * dim + k] * solu[i];
          xGauss[k] += x[k][i] * phi[i];
        }
      }

      // *** phi_i loop ***
      for(unsigned i = 0; i < nDofs; i++) {

        double nonLinearTerm = phi[i] * lambda * (soluGauss * soluGauss * soluGauss  - soluGauss);

        double mLaplace = 0.;
        for(unsigned k = 0; k < dim; k++) {
          mLaplace   +=  phi_x[i * dim + k] * soluGauss_x[k];
        }

        double exactSolValue = GetExactSolutionValue(xGauss);
        double exactSolLaplace = GetExactSolutionLaplace(xGauss);

        double f = phi[i] * (- exactSolLaplace + lambda * (exactSolValue * exactSolValue * exactSolValue - exactSolValue));
        Res[i] += (f - (mLaplace + nonLinearTerm)) * weight;

        // *** phi_j loop ***
        for(unsigned j = 0; j < nDofs; j++) {
          nonLinearTerm = phi[i] * lambda * (3. * soluGauss * soluGauss - 1.) * phi[j];
          mLaplace = 0.;

          for(unsigned k = 0; k < dim; k++) {
            mLaplace += (phi_x[i * dim + k] * phi_x[j * dim + k]);
          }
          J[i * nDofs + j] += (mLaplace + nonLinearTerm) * weight;
        } // end phi_j loop
      } // end phi_i loop
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    RES->add_vector_blocked(Res, sysDof);

    K->add_matrix_blocked(J, sysDof, sysDof);
  } //end element loop for each process

  RES->close();
  K->close();
 
  // ***************** END ASSEMBLY *******************
}


/**
 * Given the non linear problem
 *
 *        - \Delta u + lambda ( u^3 -u ) = f(x),
 *
 * in the unit box centered in the origin with
 *
 *        f(x) = - \Delta u_e + lambda ( u_e^3 - u_e),
 *        u_e = \cos ( \pi * x ) * \cos( \pi * y ),
 *
 *the following function assembles the residual std::vector Res(u^i) and using automatic differentiation gets
 *the exact Jacobian matrix J(u^i) for the Newton iteration, i.e.
 *
 *        J(u^i) w = Res(u^i) = f(x) - ( - \Delta u^i + lambda ( u_i^3 - u^i ),
 *        u^{i+1} = u^i + w,
 *        where
 *        J(u^i) w = - \Delta w  + lambda ( 3 u_i ^2 w - w ).
 *
 **/

void AssembleNonlinearProblem_AD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use

  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system< NonLinearImplicitSystem > ("Chafee–Infante");   // pointer to the linear implicit system named "Chafee–Infante"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    K         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual std::vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  const unsigned lambda = 1;
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex;
  soluIndex = mlSol->GetIndex("u");    // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex("u");    // get the position of "u" in the pdeSys object

  std::vector < adept::adouble >  solu; // local solution

  std::vector < std::vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector< unsigned > sysDof; // local to global pdeSys dofs
  std::vector <double> phi;  // local test function
  std::vector <double> phi_x; // local test function first order partial derivatives
  double weight; // gauss point weight

  std::vector< double > Res; // local redidual std::vector
  std::vector< adept::adouble > mRes; // local redidual std::vector

  std::vector < double > Jac; // local Jacobian matrix (ordered by column, adept)

  K->zero(); // Set to zero all the entries of the Global Matrix
  RES->zero(); // Set to zero all the entries of the Global residual

  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofs  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs

    // resize local arrays
    sysDof.resize(nDofs);
    solu.resize(nDofs);

    for(unsigned i = 0; i < dim; i++) {
      x[i].resize(nDofs);
    }

    mRes.assign(nDofs, 0.);

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
      sysDof[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof

      // local storage of coordinates
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k])(xDof);      // global extraction and local storage for the element coordinates
      }
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      adept::adouble soluGauss = 0;
      std::vector < adept::adouble > soluGauss_x(dim, 0.);
      std::vector < double > xGauss(dim, 0.);

      for(unsigned i = 0; i < nDofs; i++) {
        soluGauss += phi[i] * solu[i];

        for(unsigned k = 0; k < dim; k++) {
          soluGauss_x[k] += phi_x[i * dim + k] * solu[i];
          xGauss[k] += x[k][i] * phi[i];
        }
      }

      // *** phi_i loop ***
      for(unsigned i = 0; i < nDofs; i++) {

        adept::adouble nonLinearTerm = phi[i] * lambda * (soluGauss * soluGauss * soluGauss - soluGauss);
        adept::adouble mLaplace = 0.;

        for(unsigned k = 0; k < dim; k++) {
          mLaplace   +=  phi_x[i * dim + k] * soluGauss_x[k];
        }

        double exactSolValue = GetExactSolutionValue(xGauss);
        std::vector < double > exactSolGrad(dim);
        GetExactSolutionGradient(xGauss , exactSolGrad);
        double exactSolLaplace = GetExactSolutionLaplace(xGauss);


        double f = phi[i] * (- exactSolLaplace + lambda * (exactSolValue * exactSolValue * exactSolValue - exactSolValue)) ;
        mRes[i] -= (f - (mLaplace + nonLinearTerm)) * weight;

      } // end phi_i loop
    } // end gauss point loop


    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    Res.resize(nDofs);
    //copy the value of the adept::adoube mRes in double Res and store
    for(int i = 0; i < nDofs; i++) {
      Res[i] = -mRes[i].value();
    }
    RES->add_vector_blocked(Res, sysDof);

    // define the dependent variables
    s.dependent(&mRes[0], nDofs);

    // define the independent variables
    s.independent(&solu[0], nDofs);

    Jac.resize(nDofs * nDofs);
    // get the jacobian matrix (ordered by column)
    s.jacobian(&Jac[0], true);

    //store Jact in the global matrix KK
    K->add_matrix_blocked(Jac, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();

  K->close();

  // ***************** END ASSEMBLY *******************
}


std::pair < double, double > GetErrorNorm(MultiLevelSolution * mlSol) {

  unsigned level = mlSol->_mlMesh->GetNumberOfLevels() - 1u;
  //  extract pointers to the several objects that we are going to use
  Mesh*     msh = mlSol->_mlMesh->GetLevel(level);    // pointer to the mesh (level) object
  elem*     el  = msh->el;  // pointer to the elem object in msh (level)
  Solution* sol = mlSol->GetSolutionLevel(level);    // pointer to the solution (level) object

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned iproc = msh->processor_id(); // get the process_id (for parallel computation)
  //I believe you covered this, but how do you specify the partition of the mesh for parrelel computation, I remember you showing this but I don't remember how you know this partition

  //solution variable
  unsigned soluIndex = mlSol->GetIndex("u");    // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  std::vector < double >  solu; // local solution
  std::vector < vector < double > > x(dim);  // local coordinates

  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE BIQUADRATIC)

  std::vector <double> phi;  // local test function
  std::vector <double> phi_x; // local test function first order partial derivatives
  double weight; // gauss point weight

  double norm[2] = {0., 0.};

  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {


    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs

    // resize local arrays
    solu.resize(nDofu);
    for(unsigned i = 0; i < dim; i++) {
      x[i].resize(nDofu);
    }

    // local storage of global solution and coordinates
    for(unsigned i = 0; i < nDofu; i++) {
      unsigned uDof = msh->GetSolutionDof(i, iel, soluType);    // global to local mapping of soluType variables
      solu[i] = (*sol->_Sol[soluIndex])(uDof);      // global to local u solution extraction

      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to local mapping of xType variables
      for(unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k])(xDof);  // global to local coordinate extraction
      }
    }


    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      double solu_gss = 0;
      std::vector < double > gradSolu_gss(dim, 0.);
      std::vector < double > x_gss(dim, 0.);

      for(unsigned i = 0; i < nDofu; i++) {
        solu_gss += phi[i] * solu[i];

        for(unsigned k = 0; k < dim; k++) {
          gradSolu_gss[k] += phi_x[i * dim + k] * solu[i];
          x_gss[k] += x[k][i] * phi[i];
        }
      }

      double exactSol = GetExactSolutionValue(x_gss);
      norm[0] += (exactSol - solu_gss) * (exactSol - solu_gss) * weight;

      std::vector <double> exactGradSol(dim);
      GetExactSolutionGradient(x_gss, exactGradSol);

      for(unsigned k = 0; k < dim ; k++) {
        norm[1]   += ((gradSolu_gss[k] - exactGradSol[k]) * (gradSolu_gss[k] - exactGradSol[k])) * weight;
      }

    } // end gauss point loop
  } //end element loop for each process

  double normSum[2] = {0., 0.};
  MPI_Reduce(&norm, &normSum, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);  //MPI_library

  return std::pair < double, double > (sqrt(normSum[0]), sqrt(normSum[1]));//use of std pair with l2 and semi norms

}
