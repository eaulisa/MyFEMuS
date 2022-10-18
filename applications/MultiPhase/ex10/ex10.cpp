/** \file Ex6.cpp
 *  \brief This example shows how to set and solve the weak form
 *   of the Navier-Stokes Equation
 *
 *  \f{eqnarray*}
 *  && \mathbf{V} \cdot \nabla \mathbf{V} - \nabla \cdot \nu (\nabla \mathbf{V} +(\nabla \mathbf{V})^T)
 *  +\nabla P = 0 \\
 *  && \nabla \cdot \mathbf{V} = 0
 *  \f}
 *  in a unit box domain (in 2D and 3D) with given vertical velocity 1 on
 *  the left boundary and walls elsewhere.
 *  \author Eugenio Aulisa
 */

#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "adept.h"


#include "petsc.h"
#include "petscmat.h"
#include "PetscMatrix.hpp"


#include "CutFemWeight.hpp"
#include "CDWeights.hpp"
#include "Fem.hpp"

typedef double TypeIO;
typedef cpp_bin_float_oct TypeA;
typedef cpp_bin_float_oct oct;

// CutFemWeight <double, double> quad = CutFemWeight<double, double>(QUAD, 5, "legendre");
CutFemWeight <TypeIO, TypeA> quad  = CutFemWeight<TypeIO, TypeA >(QUAD, 1, "legendre");
CutFemWeight <TypeIO, TypeA> tri  = CutFemWeight<TypeIO, TypeA >(TRI, 1, "legendre");
Fem fem = Fem(quad.GetGaussQuadratureOrder(), quad.GetDimension());


#include "../include/MyMarker/MyMarker.hpp"
#include "../include/MyMarker/MyMarker.cpp"
#include "../include/Cloud.hpp"
#include "MyEigenFunctions.hpp"

#include <fstream>
#include <iostream>

using namespace femus;

void SetVelocity(Solution *sol, const std::vector<std::string> &U, const double &time, const double &T);

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = false; //dirichlet

  if(!strcmp(SolName, "U")) {  // strcmp compares two string in lexiographic sense.
    value = 0;
  }
  else if(!strcmp(SolName, "V")) {
    value = 0;
  }
  else if(!strcmp(SolName, "W")) {
    value = 0.;
  }
  else if(!strcmp(SolName, "P")) {
    dirichlet = false;
    value = 0.;
  }
  return dirichlet;
}

double SetInitialCondition(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {

  double value = 0.;

  if(!strcmp(name, "U")) {
    value = -x[1];
    //value = -0.1;
  }
  else if(!strcmp(name, "V")) {
    value = x[0];
    //value = 0;
  }
  else if(!strcmp(name, "W")) {
    value = 0.;
  }
  else if(!strcmp(name, "P")) {
    value = 0.;
  }

  return value;
}


void AssembleBoussinesqAppoximation_AD(MultiLevelProblem& ml_prob);    //, unsigned level, const unsigned &levelMax, const bool &assembleMatrix );
void AssembleBoussinesqAppoximation(MultiLevelProblem& ml_prob);
void TestMarkersAndCloud(MultiLevelProblem & ml_prob);

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  //mlMsh.ReadCoarseMesh("./input/cube_hex.neu", "seventh", scalingFactor);
//   mlMsh.ReadCoarseMesh("./input/square_quad.neu", "seventh", scalingFactor);
  
  mlMsh.GenerateCoarseBoxMesh(2, 2, 0, -0.5, 0.5, -0.5, 0.5, 0., 0., TRI6, "seventh"); 

  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
     probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();

  unsigned numberOfUniformLevels = 7;
  unsigned nMax = 4 * pow(2,6);
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels, numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // erase all the coarse mesh levels
  mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("U", LAGRANGE, SECOND, 2);
  mlSol.AddSolution("V", LAGRANGE, SECOND, 2);
  if(dim == 3) mlSol.AddSolution("W", LAGRANGE, SECOND, 2);
  mlSol.AddSolution("P",  DISCONTINUOUS_POLYNOMIAL, FIRST);
  
  mlSol.AddSolution("C", DISCONTINUOUS_POLYNOMIAL, ZERO, false);
  mlSol.AddSolution("Cn", LAGRANGE, SECOND, false);

  std::vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol.GetIndex("U");
  solVIndex[1] = mlSol.GetIndex("V");
  if(dim == 3) solVIndex[2] = mlSol.GetIndex("W");

//    //Taylor-hood
//    mlSol.AddSolution("U", LAGRANGE, SERENDIPITY);
//    mlSol.AddSolution("V", LAGRANGE, SERENDIPITY);
//    if (dim == 3) mlSol.AddSolution("W", LAGRANGE, SERENDIPITY);
//    mlSol.AddSolution("P", LAGRANGE, FIRST);

//    // Bad FEM pair - no LBB condition
//    mlSol.AddSolution("U", LAGRANGE, FIRST);
//    mlSol.AddSolution("V", LAGRANGE, FIRST);
//    if (dim == 3) mlSol.AddSolution("W", LAGRANGE, FIRST);
//    mlSol.AddSolution("P", LAGRANGE, FIRST);

// define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);

  mlSol.Initialize("All");
  mlSol.Initialize("U", SetInitialCondition, &mlProb);
  mlSol.Initialize("V", SetInitialCondition, &mlProb);

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.FixSolutionAtOnePoint("P");

  mlSol.GenerateBdc("All");

  // add system Poisson in mlProb as a Linear Implicit System
  NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("NS");

  // add solution "u" to system
  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");

  if(dim == 3) system.AddSolutionToSystemPDE("W");

  system.AddSolutionToSystemPDE("P");

  // attach the assembling function to system
//   system.SetAssembleFunction(TestMarkersAndCloud);

  const unsigned level = mlMsh.GetNumberOfLevels() - 1;
  Mesh*          msh          = mlProb._ml_msh->GetLevel(level);    // pointer to the mesh (level) object

  Solution* sol = mlSol.GetSolutionLevel(level);

  unsigned iproc = sol->processor_id();
  unsigned nprocs = sol->n_processors();

  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(true);

  // BEGIN Testing the class Cloud
  Cloud cld;
  Cloud cldint;
  std::vector<std::string> velocity = {"U", "V"};
  std::cout << "Testing the class Cloud \n";

  double period = 4;
  unsigned nIterations = 320;

  double time = 0.;
//   cld.AddQuadric({1.,0.,1.,0.,-0.5,0.04}, 8, sol);
//   cld.AddQuadric({1.,0.,1.,0.,+0.5,0.04}, 8, sol);
//   cldint.AddInteriorQuadric({1.,0.,1.,0.,-0.5,0.04}, 8, sol);
//   cldint.AddInteriorQuadric({1.,0.,1.,0.,+0.5,0.04}, 8, sol);
  
  cld.AddQuadric({0.,0.,0.,0.,1.,0.01}, 8, sol);
  cldint.AddInteriorQuadric({0.,0.,0.,0.,1.,0.01}, 8, sol);

  cldint.RebuildInteriorMarkers(cld, "C","Cn");
  SetVelocity(sol, velocity, time, period );
  cld.PrintCSV("markerBefore",0);
  cld.PrintCSV("marker",0);
  cldint.PrintCSV("markerInternalBefore",0);
  cldint.PrintCSV("markerInternal",0);
  
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);


  double dt = period / nIterations;

  for(unsigned it = 1; it <= nIterations; it++) {
    std::cout << "ITERATION " << it << "\n";
    
    sol->CopySolutionToOldSolution();

    time += dt;
    SetVelocity(sol, velocity, time, period);
    
    cld.RKAdvection(4, velocity, dt);
    cldint.RKAdvection(4, velocity, dt);
    cldint.PrintCSV("markerInternalBefore",it);
    
    cld.PrintCSV("markerBefore",it);
    cld.ComputeQuadraticBestFit();
    
    cld.RebuildMarkers(8, 12, 8);
    
    cldint.RebuildInteriorMarkers(cld, "C", "Cn");
    cldint.PrintCSV("markerInternal",it);
    cld.PrintCSV("marker",it);
    vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, it);

//     for(unsigned kp = 0; kp < nprocs; kp++) {
//       if(msh->processor_id() == kp) {
//         for(unsigned iel = msh->_elementOffset[kp]; iel < msh->_elementOffset[kp + 1]; iel++) {
//           std::cerr << "iel = " << iel << "   ";
//           const std::vector<double> &a = cld.GetQuadraticBestFitCoefficients(iel);
//           for(unsigned i = 0; i < a.size(); i++) std::cerr << a[i] << "  ";
//           std::cerr << "\n" << std::flush;
//         }
//       }
//       MPI_Barrier(MPI_COMM_WORLD);
//     }
//     std::cerr << std::endl;

  }
  
  // END Testing the class Cloud

  // initilaize and solve the system
  //   system.init();
  //   system.SetOuterSolver(PREONLY);
  //   system.MGsolve();


  return 0;
}

void SetVelocity(Solution *sol, const std::vector<std::string> &U, const double &time, const double &T) {

  Mesh* msh = sol->GetMesh();    // pointer to the mesh (level) object
  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  std::vector < unsigned > uIndex(dim);
  for(unsigned k = 0; k < dim; k++) {
    uIndex[k] = sol->GetIndex(U[k].c_str());
  }
  unsigned uType = sol->GetSolutionType(uIndex[0]);

  std::vector < double > xv(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned nDofsU = msh->GetElementDofNumber(iel, uType);
    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsU; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // local to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        xv[k] = (*msh->_topology->_Sol[k])(xDof);
      }
      unsigned uDof = msh->GetSolutionDof(i, iel, uType);    // local to global mapping between solution node and solution dof
      //rotation;
//       sol->_Sol[uIndex[0]]->set(uDof, -xv[1]);
//       sol->_Sol[uIndex[1]]->set(uDof, xv[0]);
      
      //single vortex;
      double x = xv[0] + 0.5;
      double y = xv[1] + 0.5;
      double u = -2. * sin(M_PI * x) * sin(M_PI * x) * sin(M_PI * y) * cos(M_PI * y) * cos(M_PI * time / T);
      double v =  2. * sin(M_PI * x) * cos(M_PI * x) * sin(M_PI * y) * sin(M_PI * y) * cos(M_PI * time / T);
      
      //double x = xv[0] + 0.25;
      //double y = xv[1] /*+ 0.5*/;
      
//       double u = - cos(M_PI * 2 * x) * cos(M_PI * 2 * y);
//       double v = - sin(M_PI * 2 * x) * sin(M_PI * 2 * y);
      sol->_Sol[uIndex[0]]->set(uDof, u);
      sol->_Sol[uIndex[1]]->set(uDof, v);
    }
  }
  for(unsigned  k = 0; k < dim; k++) {
    sol->_Sol[uIndex[k]]->close();
  }

}




void AssembleBoussinesqAppoximation_AD(MultiLevelProblem & ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<NonLinearImplicitSystem> ("NS");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual std::vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  std::vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("V");    // get the position of "V" in the ml_sol object
  if(dim == 3) solVIndex[2] = mlSol->GetIndex("W");       // get the position of "V" in the ml_sol object

  unsigned solVType = mlSol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

  unsigned solPIndex;
  solPIndex = mlSol->GetIndex("P");    // get the position of "P" in the ml_sol object
  unsigned solPType = mlSol->GetSolutionType(solPIndex);    // get the finite element type for "u"

  std::vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("U");    // get the position of "U" in the pdeSys object
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("V");    // get the position of "V" in the pdeSys object
  if(dim == 3) solVPdeIndex[2] = mlPdeSys->GetSolPdeIndex("W");

  unsigned solPPdeIndex;
  solPPdeIndex = mlPdeSys->GetSolPdeIndex("P");    // get the position of "P" in the pdeSys object

  std::vector < std::vector < adept::adouble > >  solV(dim);    // local solution
  std::vector < adept::adouble >  solP; // local solution

  std::vector< std::vector < adept::adouble > > mResV(dim);    // local redidual std::vector
  std::vector< adept::adouble > mResP; // local redidual std::vector

  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector <double> phiV;  // local test function for velocity
  std::vector <double> phiV_x; // local test function first order partial derivatives

  double* phiP; // local test function for the pressure
  double weight; // gauss point weight

  std::vector< unsigned > sysDof; // local to global pdeSys dofs
  std::vector< double > Res; // local redidual std::vector
  std::vector < double > Jac;

  RES->zero(); // Set to zero all the entries of the Global Residual std::vector
  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);

    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);    // number of solution element dofs

    unsigned nDofsVP = dim * nDofsV + nDofsP;
    // resize local arrays
    sysDof.resize(nDofsVP);

    for(unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofsV);
      coordX[k].resize(nDofsV);
    }
    solP.resize(nDofsP);

    for(unsigned  k = 0; k < dim; k++) {
      mResV[k].assign(nDofsV, 0.);
    }
    mResP.assign(nDofsP, 0.);

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // local to global mapping between solution node and solution dof

      for(unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);      // global extraction and local storage for the solution
        sysDof[k * nDofsV + i] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }

    for(unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);    // local to global mapping between solution node and solution dof
      solP[i] = (*sol->_Sol[solPIndex])(solPDof);      // global extraction and local storage for the solution
      sysDof[dim * nDofsV + i ] = pdeSys->GetSystemDof(solPIndex, solPPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // local to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solVType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solVType]->Jacobian(coordX, ig, weight, phiV, phiV_x);
      phiP = msh->_finiteElement[ielGeom][solPType]->GetPhi(ig);

      std::vector < adept::adouble > solV_gss(dim, 0);
      std::vector < std::vector < adept::adouble > > gradSolV_gss(dim);

      for(unsigned  k = 0; k < dim; k++) {
        gradSolV_gss[k].assign(dim, 0.);
      }

      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          solV_gss[k] += solV[k][i] * phiV[i];
        }
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            gradSolV_gss[k][j] += solV[k][i] * phiV_x[i * dim + j];
          }
        }
      }

      adept::adouble solP_gss = 0;
      for(unsigned i = 0; i < nDofsP; i++) {
        solP_gss += phiP[i] * solP[i];
      }

      double nu = 1. / 250.;

      // *** phiV_i loop ***
      for(unsigned i = 0; i < nDofsV; i++) {
        std::vector < adept::adouble > NSV(dim, 0.);

        for(unsigned  k = 0; k < dim; k++) {  //momentum equation in k
          for(unsigned j = 0; j < dim; j++) {  // second index j in each equation
            NSV[k]   +=  nu * phiV_x[i * dim + j] * (gradSolV_gss[k][j] + gradSolV_gss[j][k]); // laplace
            NSV[k]   +=  phiV[i] * (solV_gss[j] * gradSolV_gss[k][j]); // non-linear term
          }
          NSV[k] += -solP_gss * phiV_x[i * dim + k]; // pressure gradient
        }
        for(unsigned  k = 0; k < dim; k++) {
          mResV[k][i] += - NSV[k] * weight;
        }
      } // end phiV_i loop

      // *** phiP_i loop ***
      for(unsigned i = 0; i < nDofsP; i++) {
        for(int k = 0; k < dim; k++) {
          mResP[i] += - (-gradSolV_gss[k][k]) * phiP[i]  * weight;
        }
      } // end phiP_i loop

    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube mRes in double Res and store them in RES
    Res.resize(nDofsVP);    //resize

    for(int i = 0; i < nDofsV; i++) {
      for(unsigned  k = 0; k < dim; k++) {
        Res[ k * nDofsV + i ] = -mResV[k][i].value();
      }
    }

    for(int i = 0; i < nDofsP; i++) {
      Res[ dim * nDofsV + i ] = -mResP[i].value();
    }

    RES->add_vector_blocked(Res, sysDof);

    //Extract and store the Jacobian

    Jac.resize(nDofsVP * nDofsVP);
    // define the dependent variables

    for(unsigned  k = 0; k < dim; k++) {
      s.dependent(&mResV[k][0], nDofsV);
    }
    s.dependent(&mResP[0], nDofsP);

    // define the independent variables
    for(unsigned  k = 0; k < dim; k++) {
      s.independent(&solV[k][0], nDofsV);
    }
    s.independent(&solP[0], nDofsP);

    // get the and store jacobian matrix (row-major)
    s.jacobian(&Jac[0], true);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();
  KK->close();
  // ***************** END ASSEMBLY *******************

}


//Attempting to create J by hand
void AssembleBoussinesqAppoximation(MultiLevelProblem & ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled


  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<NonLinearImplicitSystem> ("NS");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual std::vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  std::vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("V");    // get the position of "V" in the ml_sol object
  if(dim == 3) solVIndex[2] = mlSol->GetIndex("W");       // get the position of "V" in the ml_sol object

  unsigned solVType = mlSol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

  unsigned solPIndex;
  solPIndex = mlSol->GetIndex("P");    // get the position of "P" in the ml_sol object
  unsigned solPType = mlSol->GetSolutionType(solPIndex);    // get the finite element type for "u"

  std::vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("U");    // get the position of "U" in the pdeSys object
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("V");    // get the position of "V" in the pdeSys object
  if(dim == 3) solVPdeIndex[2] = mlPdeSys->GetSolPdeIndex("W");

  unsigned solPPdeIndex;
  solPPdeIndex = mlPdeSys->GetSolPdeIndex("P");    // get the position of "P" in the pdeSys object

  std::vector < std::vector < double > >  solV(dim);    // local solution
  std::vector < double >  solP; // local solution

  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector <double> phiV;  // local test function for velocity
  std::vector <double> phiV_x; // local test function first order partial derivatives

  double* phiP; // local test function for the pressure
  double weight; // gauss point weight

  std::vector< unsigned > sysDof; // local to global pdeSys dofs
  std::vector< double > Res; // local redidual std::vector
  std::vector < double > Jac;

  RES->zero(); // Set to zero all the entries of the Global Residual std::vector
  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {


    short unsigned ielGeom = msh->GetElementType(iel);

    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);    // number of solution element dofs

    Jac.assign((3 * nDofsV * nDofsV + nDofsP) * (3 * nDofsV * nDofsV + nDofsP), 0.);

    unsigned nDofsVP = dim * nDofsV + nDofsP;
    // resize local arrays
    sysDof.resize(nDofsVP);
    Res.assign(nDofsVP, 0.);
    Jac.assign(nDofsVP * nDofsVP, 0.);

    for(unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofsV);
      coordX[k].resize(nDofsV);
    }
    solP.resize(nDofsP);


    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // local to global mapping between solution node and solution dof

      for(unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);      // global extraction and local storage for the solution
        sysDof[k * nDofsV + i] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }

    for(unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);    // local to global mapping between solution node and solution dof
      solP[i] = (*sol->_Sol[solPIndex])(solPDof);      // global extraction and local storage for the solution
      sysDof[dim * nDofsV + i ] = pdeSys->GetSystemDof(solPIndex, solPPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // local to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }



    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solVType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solVType]->Jacobian(coordX, ig, weight, phiV, phiV_x);
      phiP = msh->_finiteElement[ielGeom][solPType]->GetPhi(ig);

      std::vector < double > solV_gss(dim, 0);
      std::vector < std::vector < double > > gradSolV_gss(dim);

      for(unsigned  k = 0; k < dim; k++) {
        gradSolV_gss[k].assign(dim, 0.);
      }

      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          solV_gss[k] += solV[k][i] * phiV[i];
        }
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            gradSolV_gss[k][j] += solV[k][i] * phiV_x[i * dim + j];
          }
        }
      }

      double solP_gss = 0;
      for(unsigned i = 0; i < nDofsP; i++) {
        solP_gss += phiP[i] * solP[i];
      }

      double nu = 1. / 500.;

      // *** phiV_i loop ***
      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned  I = 0; I < dim; I++) {  //momentum equation in k
          double NSV = 0.;
          for(unsigned J = 0; J < dim; J++) {  // second index j in each equation
            NSV   +=  nu * phiV_x[i * dim + J] * (gradSolV_gss[I][J] + gradSolV_gss[J][I]); // diffusion
            NSV   +=  phiV[i] * (solV_gss[J] * gradSolV_gss[I][J]); // nonlinear term
          }
          NSV += - phiV_x[i * dim + I] * solP_gss; // pressure gradient
          Res[I * nDofsV + i] -=  NSV * weight;
        }
      } // end phiV_i loop

      // *** phiP_i loop ***
      for(unsigned i = 0; i < nDofsP; i++) {
        for(int I = 0; I < dim; I++) {
          Res[dim * nDofsV + i] -= -gradSolV_gss[I][I] * phiP[i]  * weight; //continuity
        }
      } // end phiP_i loop
      // end gauss point loop


      //--------------------------------------------------------------------------------------------------------
      // Add the local Matrix/Vector into the global Matrix/Vector

      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned I = 0; I < dim; I++) { //row velocity blocks or dimension
          unsigned VIrow = I * nDofsV + i;
          for(unsigned j = 0; j < nDofsV; j++) {
            unsigned VIcolumn = I * nDofsV + j;
            for(unsigned J = 0; J < dim ; J++) { //column velocity blocks or dimension
              unsigned VJcolumn = J * nDofsV + j;
              Jac[ VIrow * nDofsVP + VIcolumn] += nu * phiV_x[i * dim + J] * phiV_x[j * dim + J] * weight; //diagonal diffusion
              Jac[ VIrow * nDofsVP + VJcolumn] += nu * phiV_x[i * dim + J] * phiV_x[j * dim + I] * weight; //off-diagonal diffusion

              Jac[ VIrow * nDofsVP + VIcolumn] += phiV[i] * solV_gss[J] * phiV_x[j * dim + J] * weight; //diagonal nonlinear
              Jac[ VIrow * nDofsVP + VJcolumn] += phiV[i] * phiV[j] * gradSolV_gss[I][J] * weight; //off-diagonal nonlinear
            }
          }

          for(unsigned j = 0; j < nDofsP; j++) {
            unsigned Pcolumn = dim * nDofsV + j;
            Jac[VIrow * nDofsVP + Pcolumn] += - phiV_x[i * dim + I] * phiP[j] * weight; //pressure gradient
            Jac[Pcolumn * nDofsVP + VIrow] += - phiV_x[i * dim + I] * phiP[j] * weight; //continuity
          }

        }
      }

//       //--------------------------------------------------------------------------------------------------------
//       // Add the local Matrix/Vector into the global Matrix/Vector
//
//       double nonLinear = 0.;
//       double nonLinear_v = 0.;
//       double nonLinear_w = 0.;
//       double laplce = 0.;
//
//
//       for(unsigned I = 0; I < dim; I++) {
//         for(unsigned i = 0; i < nDofsV; i++) {
//           for(unsigned  J = 0; J < dim; J++) {
//             for(unsigned j = 0; j < nDofsV; j++) {
//               for(unsigned k = 0; k < dim; k++) {
//
//                 laplce = 0.;
//                 nonLinear = 0.;
//
//
//                 if(i % 3  == 0) {
//                   if(k == 0) {
//                     laplce += 2 * phiV_x[i * dim + k] * phiV_x[j * dim + k] + phiV_x[i * dim + k + 1] * phiV_x[j * dim + k + 1] + phiV_x[i * dim + k + 2] * phiV_x[j * dim + k + 2];
//                     nonLinear += (solV_gss[0] * phiV_x[j * dim + k] + solV_gss[1] * phiV_x[j * dim + k + 1] + solV_gss[2] * phiV_x[j * dim + k + 2] + gradSolV_gss[k][k]) * phiV[i];
//                   }
//
//                   else if(k == 1) {
//                     laplce += phiV_x[i * dim + k] * phiV_x[j * dim + k];
//                     nonLinear += phiV[i] * phiV[j] * gradSolV_gss[0][k];
//
//                   }
//
//                   else if(k == 2) {
//                     laplce += phiV_x[i * dim + k] * phiV_x[j * dim + k];
//                     nonLinear += phiV[i] * phiV[j] * gradSolV_gss[0][k];
//
//                   }
//                 }
//
//                 if(i % 3  == 1) {
//                   if(k == 1) {
//                     laplce += 2 * phiV_x[i * dim + k] * phiV_x[j * dim + k] + phiV_x[i * dim + k + 1] * phiV_x[j * dim + k + 1] + phiV_x[i * dim + k - 1] * phiV_x[j * dim + k - 1];
//                     nonLinear += (solV_gss[0] * phiV_x[j * dim + k] + solV_gss[1] * phiV_x[j * dim + k + 1] + solV_gss[2] * phiV_x[j * dim + k + 2] + gradSolV_gss[k][k]) * phiV[i];
//                   }
//
//                   else if(k == 0) {
//                     laplce += phiV_x[i * dim + k] * phiV_x[j * dim + k + 1];
//                     nonLinear += phiV[i] * phiV[j] * gradSolV_gss[1][k];
//
//                   }
//
//                   else if(k == 2) {
//                     laplce += phiV_x[i * dim + k] * phiV_x[j * dim + k - 1];
//                     nonLinear += phiV[i] * phiV[j] * gradSolV_gss[1][k];
//
//                   }
//                 }
//
//                 if(i % 3  == 2) {
//                   if(k == 2) {
//                     laplce += 2 * phiV_x[i * dim + k] * phiV_x[j * dim + k] + phiV_x[i * dim + k - 1] * phiV_x[j * dim + k - 1] + phiV_x[i * dim + k - 2] * phiV_x[j * dim + k - 2];
//                     nonLinear += (solV_gss[0] * phiV_x[j * dim + k] + solV_gss[1] * phiV_x[j * dim + k + 1] + solV_gss[2] * phiV_x[j * dim + k + 2] + gradSolV_gss[k][k]) * phiV[i];
//                   }
//
//                   else if(k == 0) {
//                     laplce += phiV_x[i * dim + k] * phiV_x[j * dim + k + 2];
//                     nonLinear += phiV[i] * phiV[j] * gradSolV_gss[2][k];
//
//                   }
//
//                   else if(k == 1) {
//                     laplce += phiV_x[i * dim + k] * phiV_x[j * dim + k + 1];
//                     nonLinear += phiV[i] * phiV[j] * gradSolV_gss[2][k];
//
//                   }
//                 }
//
//               }
//
//               Jac[3 * i * nDofsV +  i * nDofsP + j] = nu * (laplce + nonLinear) * weight;
//               Jac[3 * i * nDofsV +  i * nDofsP + j] = nu * laplce * weight;
//               Jac[(3 * i + 1) * nDofsV +  i * nDofsP + j] = nu * (laplce + nonLinear) * weight;
//               Jac[(3 * i + 1) * nDofsV +  i * nDofsP + j] = nu * laplce * weight;
//               Jac[(3 * i + 2) * nDofsV +  i * nDofsP + j] = nu * (laplce + nonLinear) * weight;
//               Jac[(3 * i + 2) * nDofsV +  i * nDofsP + j] = nu * laplce * weight;
//
//               Jac[3 * nDofsV * (nDofsV + i) +  nDofsV * nDofsP + i * nDofsP + j] = nu * (laplce + nonLinear) * weight;
//               Jac[3 * nDofsV * (nDofsV + i) +  nDofsV * nDofsP + i * nDofsP + j] = nu * laplce * weight;
//               Jac[3 * nDofsV * (nDofsV + i) + nDofsV + nDofsV * nDofsP + i * nDofsP + j] = nu * (laplce + nonLinear) * weight;
//               Jac[3 * nDofsV * (nDofsV + i) + nDofsV + nDofsV * nDofsP + i * nDofsP + j] = nu * laplce * weight;
//               Jac[3 * nDofsV * (nDofsV + i) + 2 * nDofsV + nDofsV * nDofsP + i * nDofsP + j] = (laplce + nonLinear) * weight;
//               Jac[3 * nDofsV * (nDofsV + i) + 2 * nDofsV + nDofsV * nDofsP + i * nDofsP + j] = laplce;
//
//               Jac[6 * nDofsV * (nDofsV + i) + 2 * nDofsV * nDofsP + i * nDofsP + j] = nu * (laplce + nonLinear) * weight;
//               Jac[6 * nDofsV * (nDofsV + i) + 2 * nDofsV * nDofsP + i * nDofsP + j] = nu * laplce * weight;
//               Jac[6 * nDofsV * (nDofsV + i) + nDofsV + 2 * nDofsV * nDofsP + i * nDofsP + j] = nu * (laplce + nonLinear) * weight;
//               Jac[6 * nDofsV * (nDofsV + i) + nDofsV + 2 * nDofsV * nDofsP + i * nDofsP + j] = nu * laplce * weight;
//               Jac[6 * nDofsV * (nDofsV + i) + 2 * nDofsV + 2 * nDofsV * nDofsP + i * nDofsP + j] = nu * (laplce + nonLinear) * weight;
//               Jac[6 * nDofsV * (nDofsV + i) + 2 * nDofsV + 2 * nDofsV * nDofsP + i * nDofsP + j] = nu * laplce * weight;
//
//
//             }
//           }
//         }
//       }
//       double minusP = 0.;
//       double minusV = 0.;
//
//       for(unsigned i = 0; i < nDofsV; i++) {
//         for(unsigned j = 0; j < nDofsP; j++) {
//           minusP = 0;
//           for(unsigned  J = 0; J < dim; J++) {
//
//             minusP += phiV_x[i * dim + J];
//           }
//
//           Jac[3 * (i + 1) * nDofsV +  i * nDofsP + j] = -phiP[j] * minusP * weight;
//           Jac[9 * (i + 1) * nDofsV + 3 * nDofsP * nDofsV + i * nDofsP + j] = -phiP[j] * minusV * weight;
//
//         }
//
//       }
//     }
//
//   }

    }

    RES->add_vector_blocked(Res, sysDof);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);


  } //end element loop for each process


  RES->close();
  KK->close();

//  VecView ( (static_cast<PetscVector*> (RES))->vec(),  PETSC_VIEWER_STDOUT_SELF);
//MatView ( (static_cast<PetscMatrix*> (KK))->mat(), PETSC_VIEWER_STDOUT_SELF);

//   PetscViewer    viewer;
//   PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL, NULL, 0, 0, 900, 900, &viewer);
//   PetscObjectSetName((PetscObject) viewer, "PWilmore matrix");
//   PetscViewerPushFormat(viewer, PETSC_VIEWER_DRAW_LG);
//   MatView((static_cast<PetscMatrix*>(KK))->mat(), viewer);
//   double a;
//   std::cin >> a;


}

/*

        for(unsigned k = 0; k < dim; k++) {

          laplce = 0.;
          nonLinear = 0.;


          if(i % 3  == 0) {
            if(k == 0) {
              laplce += 2 * phiV_x[i * dim + k] * phiV_x[j * dim + k] + phiV_x[i * dim + k + 1] * phiV_x[j * dim + k + 1] + phiV_x[i * dim + k + 2] * phiV_x[j * dim + k + 2];
              nonLinear += (solV_gss[0] * phiV_x[j * dim + k] + solV_gss[1] * phiV_x[j * dim + k + 1] + solV_gss[2] * phiV_x[j * dim + k + 2] + gradSolV_gss[k][k]) * phiV[i];
            }

            else if(k == 1) {
              laplce += phiV_x[i * dim + k] * phiV_x[j * dim + k];
              nonLinear += phiV[i] * phiV[j] * gradSolV_gss[0][k];

            }

            else if(k == 2) {
              laplce += phiV_x[i * dim + k] * phiV_x[j * dim + k];
              nonLinear += phiV[i] * phiV[j] * gradSolV_gss[0][k];

            }
          }

          if(i % 3  == 1) {
            if(k == 1) {
              laplce += 2 * phiV_x[i * dim + k] * phiV_x[j * dim + k] + phiV_x[i * dim + k + 1] * phiV_x[j * dim + k + 1] + phiV_x[i * dim + k - 2] * phiV_x[j * dim + k - 2];
              nonLinear += (solV_gss[0] * phiV_x[j * dim + k] + solV_gss[1] * phiV_x[j * dim + k + 1] + solV_gss[2] * phiV_x[j * dim + k + 2] + gradSolV_gss[k][k]) * phiV[i];
            }

            else if(k == 0) {
              laplce += phiV_x[i * dim + k] * phiV_x[j * dim + k + 1];
              nonLinear += phiV[i] * phiV[j] * gradSolV_gss[1][k];

            }

            else if(k == 2) {
              laplce += phiV_x[i * dim + k] * phiV_x[j * dim + k - 1];
              nonLinear += phiV[i] * phiV[j] * gradSolV_gss[1][k];

            }
          }

          if(i % 3  == 2) {
            if(k == 2) {
              laplce += 2 * phiV_x[i * dim + k] * phiV_x[j * dim + k] + phiV_x[i * dim + k - 1] * phiV_x[j * dim + k - 1] + phiV_x[i * dim + k - 2] * phiV_x[j * dim + k - 2];
              nonLinear += (solV_gss[0] * phiV_x[j * dim + k] + solV_gss[1] * phiV_x[j * dim + k + 1] + solV_gss[2] * phiV_x[j * dim + k + 2] + gradSolV_gss[k][k]) * phiV[i];
            }

            else if(k == 0) {
              laplce += phiV_x[i * dim + k] * phiV_x[j * dim + k + 2];
              nonLinear += phiV[i] * phiV[j] * gradSolV_gss[2][k];

            }

            else if(k == 1) {
              laplce += phiV_x[i * dim + k] * phiV_x[j * dim + k + 1];
              nonLinear += phiV[i] * phiV[j] * gradSolV_gss[2][k];

            }
          }

        }

        Jac[3 * i * nDofsV +  i * nDofsP + j] = nu * laplce + nonLinear;
        Jac[(3 * i + 1) * nDofsV +  i * nDofsP + j] = nu * laplce + nonLinear;
        Jac[(3 * i + 2) * nDofsV +  i * nDofsP + j] = nu * laplce + nonLinear;

        Jac[3 * nDofsV * (nDofsV + i) +  nDofsV * nDofsP + i * nDofsP + j] = nu * laplce + nonLinear;
        Jac[3 * nDofsV * (nDofsV + i) + nDofsV + nDofsV * nDofsP + i * nDofsP + j] = nu * laplce + nonLinear;
        Jac[3 * nDofsV * (nDofsV + i) + 2 * nDofsV + nDofsV * nDofsP + i * nDofsP + j] = laplce + nonLinear;

        Jac[6 * nDofsV * (nDofsV + i) + 2 * nDofsV * nDofsP + i * nDofsP + j] = nu * laplce + nonLinear;
        Jac[6 * nDofsV * (nDofsV + i) + nDofsV + 2 * nDofsV * nDofsP + i * nDofsP + j] = nu * laplce + nonLinear;
        Jac[6 * nDofsV * (nDofsV + i) + 2 * nDofsV + 2 * nDofsV * nDofsP + i * nDofsP + j] = nu * laplce + nonLinear;


      }
    }*/
//   }
// }
// double minusP = 0.;
// double minusV = 0.;
//
// for(unsigned i = 0; i < nDofsV; i++) {
//   for(unsigned j = 0; j < nDofsP; j++) {
//     minusP = 0;
//     for(unsigned  J = 0; J < dim; J++) {
//
//       minusP += phiV_x[i * dim + J];
//     }
//
//     Jac[3 * (i + 1) * nDofsV +  i * nDofsP + j] = -phiP[j] * minusP;
//     Jac[9 * (i + 1) * nDofsV + 3 * nDofsP * nDofsV + i * nDofsP + j] = -phiP[j] * minusV;
//
//   }
//
// }
// }
// }
// } //end element loop for each process







