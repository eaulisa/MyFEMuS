/** tutorial/Ex2
 * This example shows how to set and solve the weak form of the Poisson problem
 *                    $$ \Delta u = \Delta u_exact \text{ on }\Omega, $$
 *                    $$ u=0 \text{ on } \Gamma, $$
 * on a square domain $\Omega$ with boundary $\Gamma$;
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
#include "LinearImplicitSystem.hpp"
#include "adept.h"


using namespace femus;

bool SetBoundaryCondition(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time) {

  bool dirichlet = true; //dirichlet
  value = 0;

  return dirichlet;
}

void AssemblePoissonProblem(MultiLevelProblem& ml_prob);

void AssemblePoissonProblem_AD(MultiLevelProblem& ml_prob);

std::pair < double, double > GetErrorNorm(MultiLevelSolution* mlSol);

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  mlMsh.ReadCoarseMesh("./input/square_quad.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/cube_tet.neu", "seventh", scalingFactor);
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
    probably in furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();
  unsigned maxNumberOfMeshes;

  if(dim == 2) {
    maxNumberOfMeshes = 7;
  }
  else {
    maxNumberOfMeshes = 4;
  }

  std::vector < std::vector < double > > l2Norm; //get familiar l2Norm[i][j] std::pair, std::vector, std::map,  dynamic memory allocation
  l2Norm.resize(maxNumberOfMeshes);
  // std vector allows for dynamic memory allocation in a more efficient way than restructuring an array, the std vector uses more memory upon initialization, but it restructures in a more efficient way.
  // std pair stores two objects and is a special case of the std tuple
  //std map is a container object 

  vector < vector < double > > semiNorm;
  semiNorm.resize(maxNumberOfMeshes);

  for(unsigned i = 0; i < maxNumberOfMeshes; i++) {    // loop on the mesh level

    unsigned numberOfUniformLevels = i + 1;
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

      // add system Poisson in mlProb as a Linear Implicit System
      LinearImplicitSystem& system = mlProb.add_system < LinearImplicitSystem > ("Poisson");

      // add solution "u" to system
      system.AddSolutionToSystemPDE("u");

      // attach the assembling function to system
      system.SetAssembleFunction(AssemblePoissonProblem);

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
      vtkIO.SetGraphVariable("u");
      //vtkIO.SetDebugOutput(true);
      vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, i * 10 + j);
    }
  }

  // print the seminorm of the error and the order of convergence between different levels
  //Loops for l2norm of the error
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "l2 ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";

  for(unsigned i = 0; i < maxNumberOfMeshes; i++) {
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
        std::cout << log(l2Norm[i][j] / l2Norm[i + 1][j]) / log(2.) << "\t\t\t";//Definiton of k
      }

      std::cout << std::endl;
    }

  }

  //Loops for the seminorm of the error
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "SEMINORM ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";

  for(unsigned i = 0; i < maxNumberOfMeshes; i++) {
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
        std::cout << log(semiNorm[i][j] / semiNorm[i + 1][j]) / log(2.) << "\t\t\t";//Definiton of k
      }

      std::cout << std::endl;
    }

  }

  return 0;
}

//Exact solution
double GetExactSolutionValue(const std::vector < double >& x) {
  double pi = acos(-1.);
  return cos(pi * x[0]) * cos(pi * x[1]);
};

//Grad of exact solution    
void GetExactSolutionGradient(const std::vector < double >& x, vector < double >& solGrad) {
  double pi = acos(-1.);
  solGrad[0]  = -pi * sin(pi * x[0]) * cos(pi * x[1]);
  solGrad[1] = -pi * cos(pi * x[0]) * sin(pi * x[1]);
};

//Laplace of exact solution
double GetExactSolutionLaplace(const std::vector < double >& x) {
  double pi = acos(-1.);
  return -pi * pi * cos(pi * x[0]) * cos(pi * x[1]) - pi * pi * cos(pi * x[0]) * cos(pi * x[1]);
};

/**
 * This function assemble the stiffness matrix K and the residual vector RES
 * such that
 *                  K w = RES = F - K u0,
 * and consequently
 *        u = u0 + w satisfies K u = F
 **/
void AssemblePoissonProblem(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("Poisson");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*              K = pdeSys->_KK;  // pointer to the global stiffness matrix object in pdeSys (level)
  NumericVector*           RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex;
  soluIndex = mlSol->GetIndex("u");    // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex("u");    // get the position of "u" in the pdeSys object


  std::vector < double >  solu; // local solution

  std::vector < std::vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector <double> phi;  // local test function
  std::vector <double> phi_x; // local test function first order partial derivatives
  double weight; // gauss point weight


  std::vector< unsigned > l2GMap; // local to global mapping
  std::vector< double > Res; // local residual
  std::vector < double > Jac;

  K->zero(); // Set to zero all the entries of the Global Matrix
  RES->zero();
  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs

    // resize local arrays
    l2GMap.resize(nDofu);
    solu.resize(nDofu);

    for(unsigned i = 0; i < dim; i++) {
      x[i].resize(nDofu);
    }

    Res.assign(nDofu, 0.);
    Jac.assign(nDofu * nDofu, 0.);

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofu; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);  // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex])(solDof);  // global extraction and local storage for the solution
      l2GMap[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel); // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofu; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType); // global to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k])(xDof); // global extraction and local storage for the element coordinates
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
        solu_gss += phi[i] * solu[i]; // We dont use this one for this problem.

        for(unsigned k = 0; k < dim; k++) {
          gradSolu_gss[k] += phi_x[i * dim + k] * solu[i];
          x_gss[k] += x[k][i] * phi[i];
        }
      }

      // *** phi_i loop ***
      for(unsigned i = 0; i < nDofu; i++) {

        double gradPhiGradu = 0.;

        for(unsigned k = 0; k < dim; k++) {
          gradPhiGradu += phi_x[i * dim + k] * gradSolu_gss[k];
        }

        double srcTerm = - GetExactSolutionLaplace(x_gss);
        Res[i] += (srcTerm * phi[i] - gradPhiGradu) * weight;

        for(unsigned j = 0; j < nDofu; j++) {
          double gradPhigradPhi = 0.;

          for(unsigned k = 0; k < dim; k++) {
            gradPhigradPhi += phi_x[i * dim + k] * phi_x[j * dim + k];
          }
          Jac[i * nDofu + j] += gradPhigradPhi * weight;
        }
      } // end phi_i loop
    } // end gauss point loop
    //The above loops uses the weak formulation and Gauss integration scheme and then compares to the exact solution with Gauss Quadrature, I believe...

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //store local res in global RES
    RES->add_vector_blocked(Res, l2GMap);

    //store K in the global matrix KK
    K->add_matrix_blocked(Jac, l2GMap, l2GMap);

  } //end element loop for each process

  RES->close();
  K->close();

  // ***************** END ASSEMBLY *******************
}

/**
 * This function assemble the stiffnes matrix K and the residual vector RES
 * such that
 *                  K w = RES = F - K u0,
 * and consequently
 *        u = u0 + w satisfies K u = F
 * The matrix K is not assembled but it is obtained using automatic differentiation K = grad_u (-RES(u0))
 **/
void AssemblePoissonProblem_AD(MultiLevelProblem & ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object


  adept::Stack& s = FemusInit::_adeptStack; //Question on this?

  //  extract pointers to the several objects that we are going to use

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("Poisson");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*              K = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*           RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)

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

  std::vector <double> phi;  // local test function
  std::vector <double> phi_x; // local test function first order partial derivatives
  double weight; // gauss point weight

  std::vector< adept::adouble > mRes; // local redidual std::vector
  std::vector< unsigned > l2GMap; // local to global mapping
  std::vector< double > Res; // local redidual std::vector
  std::vector < double > Jac;

  K->zero(); // Set to zero all the entries of the Global Matrix
  RES->zero();
  // element loop: each process loops only on the elements that owns

  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs

    // resize local arrays
    l2GMap.resize(nDofu);
    solu.resize(nDofu);

    for(unsigned i = 0; i < dim; i++) {
      x[i].resize(nDofu);
    }

    mRes.assign(nDofu, 0.);
    Jac.assign(nDofu * nDofu, 0.);

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofu; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
      l2GMap[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofu; i++) {
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
      adept::adouble solu_gss = 0;
      std::vector < adept::adouble > gradSolu_gss(dim, 0.);
      std::vector < double > x_gss(dim, 0.);

      for(unsigned i = 0; i < nDofu; i++) {
        solu_gss += phi[i] * solu[i]; // We dont use this one for this problem.

        for(unsigned k = 0; k < dim; k++) {
          gradSolu_gss[k] += phi_x[i * dim + k] * solu[i];
          x_gss[k] += x[k][i] * phi[i];
        }
      }

      // *** phi_i loop ***
      for(unsigned i = 0; i < nDofu; i++) {

        adept::adouble gradPhiGradu = 0.;

        for(unsigned k = 0; k < dim; k++) {
          gradPhiGradu += phi_x[i * dim + k] * gradSolu_gss[k];
        }

        double srcTerm = - GetExactSolutionLaplace(x_gss);
        mRes[i] -= (srcTerm * phi[i] - gradPhiGradu) * weight;

      } // end phi_i loop

    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    Res.resize(nDofu);    //resize

    for(unsigned i = 0; i < nDofu; i++) {
      Res[i] = -mRes[i].value();
    }
    RES->add_vector_blocked(Res, l2GMap);

    // define the dependent variables
    s.dependent(&mRes[0], nDofu);

    // define the independent variables
    s.independent(&solu[0], nDofu);

    // get the jacobian matrix (ordered by row major )
    s.jacobian(&Jac[0], true);

    //store K in the global matrix KK
    K->add_matrix_blocked(Jac, l2GMap, l2GMap);

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

  double seminorm = 0.;
  double l2norm = 0.;

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
      l2norm += (exactSol - solu_gss) * (exactSol - solu_gss) * weight;

      std::vector <double> exactGradSol(dim);
      GetExactSolutionGradient(x_gss, exactGradSol);

      for(unsigned k = 0; k < dim ; k++) {
        seminorm   += ((gradSolu_gss[k] - exactGradSol[k]) * (gradSolu_gss[k] - exactGradSol[k])) * weight;
      }

    } // end gauss point loop
  } //end element loop for each process

  double l2normAll;
  MPI_Reduce(&l2norm, &l2normAll, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);  //MPI_library

  double seminormAll;
  MPI_Reduce(&seminorm, &seminormAll, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);  //MPI_library

  return std::pair < double, double > (sqrt(l2normAll), sqrt(seminormAll));//use of std pair with l2 and semi norms
  //I could not find the use of std map, when do you use this?

}
