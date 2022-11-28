

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
//CutFemWeight <TypeIO, TypeA> quad  = CutFemWeight<TypeIO, TypeA >(QUAD, 5, "legendre");
//CutFemWeight <TypeIO, TypeA> tri  = CutFemWeight<TypeIO, TypeA >(TRI, 5, "legendre");
//Fem fem = Fem(quad.GetGaussQuadratureOrder(), quad.GetDimension());


#include "../include/MyMarker/MyMarker.hpp"
#include "../include/MyMarker/MyMarker.cpp"
#include "../include/Cloud.hpp"
#include "MyEigenFunctions.hpp"
#include "../include/AdaptiveSplit.hpp"

const double R = 0.24;
const double XC = 0.;
const double YC = 0.;

Cloud *cld;
Cloud *cldint;

#include <fstream>
#include <iostream>

using namespace femus;

AdaptiveSplit *asplit;

void FakeAssembly(MultiLevelSolution* mlSol);

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

  mlMsh.GenerateCoarseBoxMesh(2, 2, 0, -0.5, 0.5, -0.5, 0.5, 0., 0., QUAD9, "seventh");

  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
     probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();

  unsigned numberOfUniformLevels = 2;
  unsigned nMax = 4 * pow(2, 6);
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

  mlSol.AddSolution("Q", DISCONTINUOUS_POLYNOMIAL, ZERO, false);
  mlSol.AddSolution("DIC", DISCONTINUOUS_POLYNOMIAL, ZERO, false); //double intersection conuter

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


  const unsigned level = mlMsh.GetNumberOfLevels() - 1;
  Mesh*          msh          = mlProb._ml_msh->GetLevel(level);    // pointer to the mesh (level) object

  Solution* sol = mlSol.GetSolutionLevel(level);


  unsigned iproc = sol->processor_id();
  unsigned nprocs = sol->n_processors();

  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(true);




  //return 0;

  // BEGIN Testing the class Cloud
  cld = new Cloud(sol);
  cldint = new Cloud(sol);

  std::vector<std::string> velocity = {"U", "V"};
  std::cout << "Testing the class Cloud \n";

  double period = 4;
  unsigned nIterations = 320;

  double time = 0.;
  cld->AddEllipses({{XC, YC}}, {{R, R}}, {9});
  cldint->AddInteriorEllipses({{0., 0.}}, {{R, R}});

  cldint->RebuildInteriorMarkers(*cld, "C", "Cn");
  SetVelocity(sol, velocity, time, period);
  cld->PrintCSV("markerBefore", 0);
  cld->PrintCSV("marker", 0);
  cldint->PrintCSV("markerInternalBefore", 0);
  cldint->PrintCSV("markerInternal", 0);

  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);

  unsigned qorder = 5;
  unsigned lmax = 3;
  asplit = new AdaptiveSplit(qorder, lmax);

  FakeAssembly(&mlSol);

  cld->PrintCSV("markerBefore", 1);
  cld->PrintCSV("marker", 1);
  cldint->PrintCSV("markerInternalBefore", 1);
  cldint->PrintCSV("markerInternal", 1);

  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 1);


  delete asplit;
  delete cld;
  delete cldint;

  return 0;
}


void FakeAssembly(MultiLevelSolution* mlSol) {


  const unsigned level = mlSol->_mlMesh->GetNumberOfLevels() - 1;
  Solution*  sol = mlSol->GetSolutionLevel(level);    // pointer to the solution (level) object

  Mesh*  msh          = sol->GetMesh();    // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  unsigned dim = msh->GetDimension();
  unsigned iproc = msh->processor_id(); // get the process_id (for parallel computation)

  unsigned solCIndex = mlSol->GetIndex("C");

  std::vector < std::vector < double > > xv(dim);    // local coordinates

  std::vector <double> phiV;  // local test function for velocity
  std::vector <double> phiV_x; // local test function first order partial derivatives

  double* phiP; // local test function for the pressure
  double weight; // gauss point weight

  unsigned solVType = 2;
  unsigned solPType = 0;
  unsigned solXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  double int1 = 0.;
  double intS1 = 0.;
  double volume1 = 0.;
  double volume2 = 0.;
  double surface = 0.;

  double int1b = 0.;
  double intS1b = 0.;
  double volume1b = 0.;
  double volume2b = 0.;
  double surfaceb = 0.;

  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    unsigned ielType = msh->GetElementType(iel);
    double C = (*sol->_Sol[solCIndex])(iel);
    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);    // number of solution element dofs
    unsigned nDofsX = msh->GetElementDofNumber(iel, solXType);    // number of solution element dofs

    xv.assign(dim, std::vector<double>(nDofsX));
    // local storage of coordinates
    for(unsigned i = 0; i < nDofsX; i++) {
      unsigned xvDof  = msh->GetSolutionDof(i, iel, solXType);    // local to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        xv[k][i] = (*msh->_topology->_Sol[k])(xvDof);      // global extraction and local storage for the element coordinates
      }
    }

    unsigned cut = (cld->GetNumberOfMarker(iel) > 0) ? 1 : 0;

    if(cut == 1) {

      std::vector<double> a;
      double d;
      std::vector<std::vector<double>> Jac, JacI;

      std::cout << iel << std::endl;

      const elem_type *femV = asplit->GetFiniteElementFine(ielType, solVType);
      const elem_type *femP = asplit->GetFiniteElementFine(ielType, solPType);
            
      femV->GetJacobianMatrix(xv, {1. / 3, 1. / 3}, weight, Jac, JacI); //TODO in all ex
      cld->GetElementQuantities(iel, Jac, asplit->GetXiFather(), asplit->GetDsFather(),  asplit->GetNiFather());

      
       {
        asplit->Split(xv, ielType, false);

        const std::vector <double> &weight1 = asplit->GetWeight1();
        const std::vector <double> &weight2 = asplit->GetWeight2();
        const std::vector <double> &weightI = asplit->GetWeightI();
        const std::vector<std::vector <double>> &xg = asplit->GetXg();
        const std::vector<std::vector <double>> &xig = asplit->GetXig();

        for(unsigned i = 0; i < weight1.size(); i++) {

          femV->Jacobian(xv, i, weight, phiV, phiV_x);// cheap

          int1b += (xg[i][0] * xg[i][0] + xg[i][1] * xg[i][1]) * weight1[i];
          intS1b += (xg[i][0] * xg[i][0] + xg[i][1] * xg[i][1]) * weightI[i];
        }

        volume1b += asplit->GetVolume1();
        volume2b += asplit->GetVolume2();
        surfaceb += asplit->GetSurfaceArea();
      }
           

      {
        asplit->Split(xv, ielType, true);

        const std::vector <double> &weight1 = asplit->GetWeight1();
        const std::vector <double> &weight2 = asplit->GetWeight2();
        const std::vector <double> &weightI = asplit->GetWeightI();
        const std::vector<std::vector <double>> &xg = asplit->GetXg();
        const std::vector<std::vector <double>> &xig = asplit->GetXig();

        for(unsigned i = 0; i < weight1.size(); i++) {

          femV->Jacobian(xv, xig[i], weight, phiV, phiV_x);// expensive for real adaptivity

          int1 += (xg[i][0] * xg[i][0] + xg[i][1] * xg[i][1]) * weight1[i];
          intS1 += (xg[i][0] * xg[i][0] + xg[i][1] * xg[i][1]) * weightI[i];
        }

        volume1 += asplit->GetVolume1();
        volume2 += asplit->GetVolume2();
        surface += asplit->GetSurfaceArea();

        double sum1 = 0;
        double sum2 = 0;
        for(unsigned i = 0; i < weight1.size(); i++) {
          sum1 +=  weight1[i];
          sum2 +=  weight2[i];
        }
        sol->_Sol[solCIndex]->set(iel, asplit->GetC());

      }

     
    }
    else {
      const elem_type *femV = asplit->GetFiniteElementCoarse(ielType, solVType);
      const elem_type *femP = asplit->GetFiniteElementCoarse(ielType, solPType);  
        
      for(unsigned ig = 0; ig < femV->GetGaussPointNumber(); ig++) {
        // *** get gauss point weight, test function and test function partial derivatives ***
        femV->Jacobian(xv, ig, weight, phiV, phiV_x);
        phiP = femP->GetPhi(ig);

        std::vector<double> xg(dim, 0);
        for(unsigned k = 0; k < dim; k++) {
          for(unsigned i = 0; i < xv[k].size(); i++) {
            xg[k] += phiV[i] * xv[k][i];
          }
        }


        if(C == 1) {
          int1 += (xg[0] * xg[0] + xg[1] * xg[1]) * weight;
          int1b += (xg[0] * xg[0] + xg[1] * xg[1]) * weight;
          volume1 += weight;
          volume1b += weight;
        }
        else {
          volume2 += weight;
          volume2b += weight;
        }
      }
    } //end element loop for each process
  }

  std::cout << volume1 << " " << volume2 << " " << surface << " " << volume1 + volume2 << std::endl;
  std::cout << volume1b << " " << volume2b << " " << surfaceb << " " << volume1b + volume2b << std::endl;
  std::cout << M_PI * R * R << " " << 1 - M_PI * R * R << " " << 2 * M_PI * R << " " << 1. << std::endl;


  std::cout << int1 << " " << int1b << " " << M_PI / 2 * pow(R, 4) << std::endl;
  std::cout << intS1 << " " << intS1b << " " << 2 * M_PI * pow(R, 3) << std::endl;

  sol->_Sol[solCIndex]->close();
}













void SetVelocity(Solution * sol, const std::vector<std::string> &U, const double & time, const double & T) {

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


//       double x = xv[0];
//       double y = xv[1];
//       double u = 0.;
//       double v =  (y<0.)?0.1:-0.1;
//       if(y==0) v= 0;



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





