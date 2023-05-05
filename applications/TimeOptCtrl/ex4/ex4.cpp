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
#include "LinearImplicitSystem.hpp"
#include "TransientSystem.hpp"
#include "adept.h"

#include "petsc.h"
#include "petscmat.h"
#include "PetscMatrix.hpp"

#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define BLUE    "\033[34m"      /* Blue */
#define RESET   "\033[0m"       /* Reset */

const double betaU = .00001;
const double alphaU = 0.00001;
const double gammaU = .1;
const double betaV = 0.00001;
const double alphaV = 0.00001;
const double t0 = 1.;
const double H0 = 0.5;
const double Re = 50.;
const double E0 = 0.;

const bool oneDimDisp = !false;

bool cleanFile = true;

const unsigned numberOfIterations = 4;
unsigned iext;

using namespace femus;

double flc4hs(double const &x, double const &eps);
double dflc4hs(double const &x, double const &eps);

double SetVariableTimeStep(const double time) {
  return 0.05;
}

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0;

  double c = flc4hs(time - t0, t0);
  double dc = dflc4hs(time - t0, t0);
  if(facename == 1 || facename == 3) { //outflow
    if(!strcmp(SolName, "U1c")  || !strcmp(SolName, "V1c") || !strcmp(SolName, "V2c") ||
        !strcmp(SolName, "bU1") || !strcmp(SolName, "bV1") || !strcmp(SolName, "bV2") ||
        !strcmp(SolName, "lU1") || !strcmp(SolName, "lV1") || !strcmp(SolName, "lV2") || //TO CHECK
        !strcmp(SolName, "U1i") || !strcmp(SolName, "V1i") || !strcmp(SolName, "V2i")
      ) {
      dirichlet = false;
    }
//     else if(!strcmp(SolName, "U1c")  ||
//             !strcmp(SolName, "bU1")  ||
//             !strcmp(SolName, "U1i")
//            ) {
//       value = H0 * (-dc * time - c) * M_PI / 4. * cos(M_PI / 4. * (x[1] - c * time)) * x[0];
//     }
  }
//   else if(facename == 5) {
//     if(!strcmp(SolName, "U1c")  ||
//         !strcmp(SolName, "bU1")  ||
//         !strcmp(SolName, "U1i") ||
//         !strcmp(SolName, "V1c")  ||
//         !strcmp(SolName, "bV1")  ||
//         !strcmp(SolName, "V1i")
//       ) {
//       value = H0 * (-dc * time - c) * M_PI / 4. * cos(M_PI / 4. * (x[1] - c * time));
//     }
//   }


  else if(facename == 4) { //slip condition on the axis
    if(!strcmp(SolName, "U2c")  || !strcmp(SolName, "V2c")  ||
        !strcmp(SolName, "bU2")  || !strcmp(SolName, "bV2") ||
        !strcmp(SolName, "lU2")  || !strcmp(SolName, "lV2") || //TO CHECK
        !strcmp(SolName, "U2i")  || !strcmp(SolName, "V2i")
      ) {
      dirichlet = false;
    }
  }
  else if(facename == 2 || facename == 5) {
    if(!strcmp(SolName, "U1c")) {
      value = -H0 * (-dc * time - c) * M_PI / 4. * sin(M_PI / 4. * (x[1] - c * time));
      //value = H0 * cos(time) * flc4hs(x[1] - 3, 1) * flc4hs(5 - x[1], 1);

    }
    else if(!strcmp(SolName, "V1c") || !strcmp(SolName, "V2c") || // V1c = U1c && V2c = U2c
            !strcmp(SolName, "bU1") || !strcmp(SolName, "bU2") || // the shape velocity is free to move
            !strcmp(SolName, "bV1") || !strcmp(SolName, "bV2") || // bV1 = bU1 && bV2 = bU2
            !strcmp(SolName, "U1i") || !strcmp(SolName, "U2i") || // U1i = bU1 && U2i = bU2
            !strcmp(SolName, "V1i") || !strcmp(SolName, "V2i")    // V1i = bV1 && V2i = bV2
           ) {
      dirichlet = false;
    }
  }


  if(!strcmp(SolName, "Pc")  || !strcmp(SolName, "bP") || !strcmp(SolName, "lP") || !strcmp(SolName, "Pi")) {
    dirichlet = false;
  }



//   if(!strcmp(SolName, "bV2") || !strcmp(SolName, "V2i")) {
//     if(x[0] < 0. && x[1] < 0.5 && x[1] > -0.5 && x[2] < 0.5 && x[2] > -0.5) dirichlet = false;
//   }
//   if(!strcmp(SolName, "V2c")) {
//     if(x[0] < 0. && x[1] < 0.5 && x[1] > -0.5 && x[2] < 0.5 && x[2] > -0.5) value = flc4hs(time - t0, t0) - 2 * flc4hs(time - 6 * t0, t0) + 2 * flc4hs(time - 11 * t0, t0) - 2 * flc4hs(time - 16 * t0, t0) + 2 * flc4hs(time - 21 * t0, t0);
//   }
//   else if(!strcmp(SolName, "bP") || !strcmp(SolName, "lP") || !strcmp(SolName, "Pc")) {
//     dirichlet = false;
//     value = 0.;
//   }

  return dirichlet;
}

void AssembleSteadyStateControl(MultiLevelProblem& ml_prob);
void AssembleSystemZi(MultiLevelProblem& ml_prob);
void AssembleManifactureSolution(MultiLevelProblem& ml_prob);
void GetError(MultiLevelProblem& ml_prob);

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  //mlMsh.ReadCoarseMesh("./input/cube_hex.neu", "seventh", scalingFactor);
  mlMsh.ReadCoarseMesh("./input/channel.neu", "seventh", scalingFactor);
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
     probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();

  unsigned numberOfUniformLevels = 3;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels, numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // erase all the coarse mesh levels
  mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("U1c", LAGRANGE, SECOND, 2);
  mlSol.AddSolution("U2c", LAGRANGE, SECOND, 2);
  mlSol.AddSolution("V1c", LAGRANGE, SECOND, 2);
  mlSol.AddSolution("V2c", LAGRANGE, SECOND, 2);
  mlSol.AddSolution("Pc",  DISCONTINUOUS_POLYNOMIAL, FIRST);

  // add variables to mlSol
  mlSol.AddSolution("bU1", LAGRANGE, SECOND);
  mlSol.AddSolution("bU2", LAGRANGE, SECOND);
  mlSol.AddSolution("bV1", LAGRANGE, SECOND);
  mlSol.AddSolution("bV2", LAGRANGE, SECOND);
  mlSol.AddSolution("bP",  DISCONTINUOUS_POLYNOMIAL, FIRST);

  mlSol.AddSolution("lU1", LAGRANGE, SECOND);
  mlSol.AddSolution("lU2", LAGRANGE, SECOND);
  mlSol.AddSolution("lV1", LAGRANGE, SECOND);
  mlSol.AddSolution("lV2", LAGRANGE, SECOND);
  mlSol.AddSolution("lP",  DISCONTINUOUS_POLYNOMIAL, FIRST);

  mlSol.AddSolution("U1i", LAGRANGE, SECOND, 2);
  mlSol.AddSolution("U2i", LAGRANGE, SECOND, 2);
  mlSol.AddSolution("V1i", LAGRANGE, SECOND, 2);
  mlSol.AddSolution("V2i", LAGRANGE, SECOND, 2);
  mlSol.AddSolution("Pi",  DISCONTINUOUS_POLYNOMIAL, FIRST);

  char u1Name[10];
  char u2Name[10];
  char v1Name[10];
  char v2Name[10];
  char pName[10];
  const unsigned level = 0;
  Solution* sol = mlSol.GetSolutionLevel(level);
  Mesh* msh = sol->GetMesh();
  unsigned    iproc = msh->processor_id();

  for(unsigned i = 0; i < numberOfIterations; i++) {
    sprintf(u1Name, "U1%d", i); //cascade solution
    mlSol.AddSolution(u1Name, LAGRANGE, SECOND, 2); //
    sprintf(u2Name, "U2%d", i); //cascade solution
    mlSol.AddSolution(u2Name, LAGRANGE, SECOND, 2); //
    sprintf(v1Name, "V1%d", i); //cascade solution
    mlSol.AddSolution(v1Name, LAGRANGE, SECOND, 2); //
    sprintf(v2Name, "V2%d", i); //cascade solution
    mlSol.AddSolution(v2Name, LAGRANGE, SECOND, 2); //
    sprintf(pName, "P%d", i); //cascade solution
    mlSol.AddSolution(pName,  DISCONTINUOUS_POLYNOMIAL, FIRST);
  }
  mlSol.AddSolution("U10Older", LAGRANGE, SECOND, false);
  mlSol.AddSolution("U20Older", LAGRANGE, SECOND, false);
  mlSol.AddSolution("V10Older", LAGRANGE, SECOND, false);
  mlSol.AddSolution("V20Older", LAGRANGE, SECOND, false);

  mlSol.AddSolution("DXc", LAGRANGE, SECOND, false);
  mlSol.AddSolution("DYc", LAGRANGE, SECOND, false);

  mlSol.AddSolution("DXi", LAGRANGE, SECOND, false);
  mlSol.AddSolution("DYi", LAGRANGE, SECOND, false);



  mlSol.Initialize("All");

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

  mlSol.GenerateBdc("All");
  mlSol.GenerateBdc("U1c", "Time_dependent");
  mlSol.GenerateBdc("bU1", "Time_dependent");
  mlSol.GenerateBdc("U1i", "Time_dependent");
  mlSol.GenerateBdc("V1c", "Time_dependent");
  mlSol.GenerateBdc("bV1", "Time_dependent");
  mlSol.GenerateBdc("V1i", "Time_dependent");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);

  // add system Poisson in mlProb as a Linear Implicit System
  TransientNonlinearImplicitSystem& systemC = mlProb.add_system < TransientNonlinearImplicitSystem > ("ManSol");

  systemC.AddSolutionToSystemPDE("U1c");
  systemC.AddSolutionToSystemPDE("U2c");
  systemC.AddSolutionToSystemPDE("V1c");
  systemC.AddSolutionToSystemPDE("V2c");
  systemC.AddSolutionToSystemPDE("Pc");

  // attach the assembling function to system
  systemC.SetAssembleFunction(AssembleManifactureSolution);
  systemC.AttachGetTimeIntervalFunction(SetVariableTimeStep);

  // initilaize and solve the system
  systemC.init();
  systemC.SetOuterSolver(PREONLY);

  TransientLinearImplicitSystem& system = mlProb.add_system < TransientLinearImplicitSystem > ("systembZ");

  // add solution "u" to system

  system.AddSolutionToSystemPDE("bU1");
  system.AddSolutionToSystemPDE("bU2");
  system.AddSolutionToSystemPDE("bV1");
  system.AddSolutionToSystemPDE("bV2");
  system.AddSolutionToSystemPDE("bP");

  system.AddSolutionToSystemPDE("lU1");
  system.AddSolutionToSystemPDE("lU2");
  system.AddSolutionToSystemPDE("lV1");
  system.AddSolutionToSystemPDE("lV2");
  system.AddSolutionToSystemPDE("lP");

  // attach the assembling function to system
  system.SetAssembleFunction(AssembleSteadyStateControl);
  system.AttachGetTimeIntervalFunction(SetVariableTimeStep);


  // initilaize and solve the system
  system.init();
  system.SetOuterSolver(PREONLY);


  TransientNonlinearImplicitSystem& systemi = mlProb.add_system < TransientNonlinearImplicitSystem > ("systemZi");

  systemi.AddSolutionToSystemPDE("U1i");
  systemi.AddSolutionToSystemPDE("U2i");
  systemi.AddSolutionToSystemPDE("V1i");
  systemi.AddSolutionToSystemPDE("V2i");
  systemi.AddSolutionToSystemPDE("Pi");

  // attach the assembling function to system
  systemi.SetAssembleFunction(AssembleSystemZi);
  systemi.AttachGetTimeIntervalFunction(SetVariableTimeStep);

  // initilaize and solve the system
  systemi.init();
  systemi.SetOuterSolver(PREONLY);

  std::vector<unsigned> solUiIndex(dim);
  solUiIndex[0] = mlSol.GetIndex("U1i");
  solUiIndex[1] = mlSol.GetIndex("U2i");

  std::vector<unsigned> solUcIndex(dim);
  solUcIndex[0] = mlSol.GetIndex("U1c");
  solUcIndex[1] = mlSol.GetIndex("U2c");

  std::vector<unsigned> solDXiIndex(dim);
  solDXiIndex[0] = mlSol.GetIndex("DXi");
  solDXiIndex[1] = mlSol.GetIndex("DYi");

  std::vector<unsigned> solDXcIndex(dim);
  solDXcIndex[0] = mlSol.GetIndex("DXc");
  solDXcIndex[1] = mlSol.GetIndex("DYc");


  for(unsigned i = msh->_dofOffset[2][iproc]; i < msh->_dofOffset[2][iproc + 1]; i++) {
    double x = (*msh->_topology->_Sol[0])(i);
    double y = (*msh->_topology->_Sol[1])(i);
    double u = H0 * x * cos(M_PI / 4. * y);

    sol->_Sol[solDXiIndex[0]]->set(i, u);
    sol->_Sol[solDXcIndex[0]]->set(i, u);

    //msh->_topology->_Sol[0]->set(i, x + u);
  }

  sol->_Sol[solDXiIndex[0]]->close();
  sol->_Sol[solDXcIndex[0]]->close();

  //msh->_topology->_Sol[0]->close();


  // print solutions
  std::vector < std::string > variablesToBePrinted = {"U1i", "U2i", "V1i", "V2i", "Pi"};
  std::vector<std::string> mov_vars = {"DXi", "DYi"};

  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(false);
  vtkIO.SetMovingMesh(mov_vars);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);


  std::vector < std::string > variablesToBePrintedc = {"U1c", "U2c", "V1c", "V2c", "Pc"};
  std::vector<std::string> mov_varsc = {"DXc", "DYc"};

  VTKWriter vtkIOc(&mlSol);
  vtkIOc.SetDebugOutput(false);
  vtkIOc.SetMovingMesh(mov_varsc);
  vtkIOc.Write("outputc", "biquadratic", variablesToBePrintedc, 0);


  std::vector < std::string > variablesToBePrinted1 = {"bU1", "bU2", "bV1", "bV2", "bP", "lU1", "lU2", "lV1", "lV2", "lP"};
  VTKWriter vtkIO1(&mlSol);
  vtkIO1.SetDebugOutput(true);
  vtkIO1.SetMovingMesh(mov_vars);
  vtkIO1.Write("output1", "biquadratic", variablesToBePrinted1, 0);

  *(sol->_SolOld[mlProb._ml_sol->GetIndex("U10")]) = *(sol->_Sol[mlProb._ml_sol->GetIndex("U10")]);
  *(sol->_SolOld[mlProb._ml_sol->GetIndex("U20")]) = *(sol->_Sol[mlProb._ml_sol->GetIndex("U20")]);
  *(sol->_SolOld[mlProb._ml_sol->GetIndex("V10")]) = *(sol->_Sol[mlProb._ml_sol->GetIndex("V10")]);
  *(sol->_SolOld[mlProb._ml_sol->GetIndex("V20")]) = *(sol->_Sol[mlProb._ml_sol->GetIndex("V20")]);

  for(unsigned t = 0; t < 500; t++) {
    double dt =  system.GetIntervalTime();
    double time =  systemC.GetTime();

    *(sol->_Sol[mlProb._ml_sol->GetIndex("U10Older")]) = *(sol->_SolOld[mlProb._ml_sol->GetIndex("U10")]);
    *(sol->_Sol[mlProb._ml_sol->GetIndex("U20Older")]) = *(sol->_SolOld[mlProb._ml_sol->GetIndex("U20")]);
    *(sol->_Sol[mlProb._ml_sol->GetIndex("V10Older")]) = *(sol->_SolOld[mlProb._ml_sol->GetIndex("V10")]);
    *(sol->_Sol[mlProb._ml_sol->GetIndex("V20Older")]) = *(sol->_SolOld[mlProb._ml_sol->GetIndex("V20")]);

    mlSol.CopySolutionToOldSolution();
    systemC.MGsolve();

    for(iext = 0; iext < numberOfIterations; iext++) {
      sprintf(u1Name, "U1%d", iext);
      sprintf(u2Name, "U2%d", iext);
      sprintf(v1Name, "V1%d", iext);
      sprintf(v2Name, "V2%d", iext);
      sprintf(pName, "P%d", iext);


      *(sol->_Sol[mlProb._ml_sol->GetIndex("U1i")]) = *(sol->_Sol[mlProb._ml_sol->GetIndex(u1Name)]);
      *(sol->_SolOld[mlProb._ml_sol->GetIndex("U1i")]) = *(sol->_SolOld[mlProb._ml_sol->GetIndex(u1Name)]);
      *(sol->_Sol[mlProb._ml_sol->GetIndex("U2i")]) = *(sol->_Sol[mlProb._ml_sol->GetIndex(u2Name)]);
      *(sol->_SolOld[mlProb._ml_sol->GetIndex("U2i")]) = *(sol->_SolOld[mlProb._ml_sol->GetIndex(u2Name)]);
      *(sol->_Sol[mlProb._ml_sol->GetIndex("V1i")]) = *(sol->_Sol[mlProb._ml_sol->GetIndex(v1Name)]);
      *(sol->_SolOld[mlProb._ml_sol->GetIndex("V1i")]) = *(sol->_SolOld[mlProb._ml_sol->GetIndex(v1Name)]);
      *(sol->_Sol[mlProb._ml_sol->GetIndex("V2i")]) = *(sol->_Sol[mlProb._ml_sol->GetIndex(v2Name)]);
      *(sol->_SolOld[mlProb._ml_sol->GetIndex("V2i")]) = *(sol->_SolOld[mlProb._ml_sol->GetIndex(v2Name)]);
      *(sol->_Sol[mlProb._ml_sol->GetIndex("Pi")]) = *(sol->_Sol[mlProb._ml_sol->GetIndex(pName)]);

      system.SetTime(time);
      systemi.SetTime(time);

      system.MGsolve();
      systemi.MGsolve();

      *(sol->_Sol[mlProb._ml_sol->GetIndex(u1Name)]) = *(sol->_Sol[mlProb._ml_sol->GetIndex("U1i")]);
      *(sol->_SolOld[mlProb._ml_sol->GetIndex(u1Name)]) = *(sol->_SolOld[mlProb._ml_sol->GetIndex("U1i")]);
      *(sol->_Sol[mlProb._ml_sol->GetIndex(u2Name)]) = *(sol->_Sol[mlProb._ml_sol->GetIndex("U2i")]);
      *(sol->_SolOld[mlProb._ml_sol->GetIndex(u2Name)]) = *(sol->_SolOld[mlProb._ml_sol->GetIndex("U2i")]);
      *(sol->_Sol[mlProb._ml_sol->GetIndex(v1Name)]) = *(sol->_Sol[mlProb._ml_sol->GetIndex("V1i")]);
      *(sol->_SolOld[mlProb._ml_sol->GetIndex(v1Name)]) = *(sol->_SolOld[mlProb._ml_sol->GetIndex("V1i")]);
      *(sol->_Sol[mlProb._ml_sol->GetIndex(v2Name)]) = *(sol->_Sol[mlProb._ml_sol->GetIndex("V2i")]);
      *(sol->_SolOld[mlProb._ml_sol->GetIndex(v2Name)]) = *(sol->_SolOld[mlProb._ml_sol->GetIndex("V2i")]);
      *(sol->_Sol[mlProb._ml_sol->GetIndex(pName)]) = *(sol->_Sol[mlProb._ml_sol->GetIndex("Pi")]);


      //GetError(mlProb);

    }

    for(unsigned k = 0; k < dim; k++)  {
      for(unsigned i = msh->_dofOffset[2][iproc]; i < msh->_dofOffset[2][iproc + 1]; i++) {

        double u = (*sol->_Sol[solUiIndex[k]])(i);
        double dx = (*sol->_Sol[solDXiIndex[k]])(i);
        sol->_Sol[solDXiIndex[k]]->set(i, dx + u * dt);

        u = (*sol->_Sol[solUcIndex[k]])(i);
        dx = (*sol->_Sol[solDXcIndex[k]])(i);
        sol->_Sol[solDXcIndex[k]]->set(i, dx + u * dt);

      }
      sol->_Sol[solDXiIndex[k]]->close();
      sol->_Sol[solDXcIndex[k]]->close();
    }
    std::cout << BLUE << " /***********************************************************/\n" << RESET;
    std::cout << BLUE << " /***********************************************************/\n" << RESET;
    std::cout << BLUE << " /*********************** REAL OUTPUT ***********************/\n" << RESET;
    std::cout << BLUE << " /***********************************************************/\n" << RESET;
    std::cout << BLUE << " /***********************************************************/\n" << RESET;
    vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, t + 1);
    vtkIOc.Write("outputc", "biquadratic", variablesToBePrintedc, t + 1);
    vtkIO1.Write("output1", "biquadratic", variablesToBePrinted1, t + 1);
  }

  return 0;
}


void AssembleSteadyStateControl(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  TransientLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<TransientLinearImplicitSystem> ("systembZ");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh    = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el     = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual std::vector object in pdeSys (level)

  double dt =  mlPdeSys->GetIntervalTime();
  double time =  mlPdeSys->GetTime();

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  std::vector < unsigned > solVcIndex(dim);
  solVcIndex[0] = mlSol->GetIndex("V1c");
  solVcIndex[1] = mlSol->GetIndex("V2c");

  char U1i[10], U2i[10], V1i[10], V2i[10];
  sprintf(U1i, "U1%d", iext);
  sprintf(U2i, "U2%d", iext);
  sprintf(V1i, "V1%d", iext);
  sprintf(V2i, "V2%d", iext);


  std::vector < unsigned > solUiIndex(dim);
  std::vector < unsigned > solViIndex(dim);
  solUiIndex[0] = mlSol->GetIndex(U1i);
  solUiIndex[1] = mlSol->GetIndex(U2i);
  solViIndex[0] = mlSol->GetIndex(V1i);
  solViIndex[1] = mlSol->GetIndex(V2i);

  char U1im1[10], U2im1[10], V1im1[10], V2im1[10];
  if(iext > 0) {
    sprintf(U1im1, "U1%d", iext - 1);
    sprintf(U2im1, "U2%d", iext - 1);
    sprintf(V1im1, "V1%d", iext - 1);
    sprintf(V2im1, "V2%d", iext - 1);
  }
  else {
    sprintf(U1im1, "U1%d", 0);
    sprintf(U2im1, "U2%d", 0);
    sprintf(V1im1, "V1%d", 0);
    sprintf(V2im1, "V2%d", 0);
  }

  std::vector < unsigned > solUim1Index(dim);
  std::vector < unsigned > solVim1Index(dim);
  solUim1Index[0] = mlSol->GetIndex(U1im1);
  solUim1Index[1] = mlSol->GetIndex(U2im1);
  solVim1Index[0] = mlSol->GetIndex(V1im1);
  solVim1Index[1] = mlSol->GetIndex(V2im1);

  std::vector < unsigned > solU0OlderIndex(dim);
  std::vector < unsigned > solV0OlderIndex(dim);
  solU0OlderIndex[0] = mlSol->GetIndex("U10Older");
  solU0OlderIndex[1] = mlSol->GetIndex("U20Older");
  solV0OlderIndex[0] = mlSol->GetIndex("V10Older");
  solV0OlderIndex[1] = mlSol->GetIndex("V20Older");

  //solution variable
  std::vector < unsigned > solbUIndex(dim);
  std::vector < unsigned > solbVIndex(dim);
  std::vector < unsigned > sollUIndex(dim);
  std::vector < unsigned > sollVIndex(dim);


  solbUIndex[0] = mlSol->GetIndex("bU1");
  solbUIndex[1] = mlSol->GetIndex("bU2");
  solbVIndex[0] = mlSol->GetIndex("bV1");
  solbVIndex[1] = mlSol->GetIndex("bV2");

  sollUIndex[0] = mlSol->GetIndex("lU1");
  sollUIndex[1] = mlSol->GetIndex("lU2");
  sollVIndex[0] = mlSol->GetIndex("lV1");
  sollVIndex[1] = mlSol->GetIndex("lV2");

  unsigned solType = mlSol->GetSolutionType(solbUIndex[0]);

  unsigned solbPIndex;
  unsigned sollPIndex;
  solbPIndex = mlSol->GetIndex("bP");
  sollPIndex = mlSol->GetIndex("lP");
  unsigned solPType = mlSol->GetSolutionType(solbPIndex);

  std::vector < unsigned > solbUPdeIndex(dim);
  std::vector < unsigned > solbVPdeIndex(dim);
  std::vector < unsigned > sollUPdeIndex(dim);
  std::vector < unsigned > sollVPdeIndex(dim);
  solbUPdeIndex[0] = mlPdeSys->GetSolPdeIndex("bU1");
  solbUPdeIndex[1] = mlPdeSys->GetSolPdeIndex("bU2");
  solbVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("bV1");
  solbVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("bV2");

  sollUPdeIndex[0] = mlPdeSys->GetSolPdeIndex("lU1");
  sollUPdeIndex[1] = mlPdeSys->GetSolPdeIndex("lU2");
  sollVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("lV1");
  sollVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("lV2");


  unsigned solbPPdeIndex;
  unsigned sollPPdeIndex;
  solbPPdeIndex = mlPdeSys->GetSolPdeIndex("bP");
  sollPPdeIndex = mlPdeSys->GetSolPdeIndex("lP");


  std::vector < std::vector < double > >  solVc(dim);

  std::vector < std::vector < double > >  solUim1(dim);
  std::vector < std::vector < double > >  solVim1(dim);
  std::vector < std::vector < double > >  solVim1Old(dim);

  std::vector < std::vector < double > >  solUiOld(dim);
  std::vector < std::vector < double > >  solViOld(dim);

  std::vector < std::vector < adept::adouble > >  solbU(dim);
  std::vector < std::vector < adept::adouble > >  solbV(dim);
  std::vector < std::vector < adept::adouble > >  sollU(dim);
  std::vector < std::vector < adept::adouble > >  sollV(dim);

  std::vector < adept::adouble >  solbP;
  std::vector < adept::adouble >  sollP;

  std::vector< std::vector < adept::adouble > > mResbU(dim);
  std::vector< std::vector < adept::adouble > > mResbV(dim);
  std::vector< std::vector < adept::adouble > > mReslU(dim);
  std::vector< std::vector < adept::adouble > > mReslV(dim);
  std::vector< adept::adouble > mResbP;
  std::vector< adept::adouble > mReslP;

  std::vector < std::vector < adept::adouble > > x(dim);
  std::vector < std::vector < double > > xOld(dim);
  std::vector < std::vector < double > > xm1(dim);

  std::vector < std::vector < adept::adouble > > DX(dim);


  unsigned xType = 2;

  std::vector <double> phi;
  std::vector <double> phix;

  double* phiP;
  double weight;
  double weightOld;

  std::vector < unsigned > solDXIndex(dim);
  solDXIndex[0] = mlSol->GetIndex("DXi");
  solDXIndex[1] = mlSol->GetIndex("DYi");
  std::vector < std::vector < double > > xHat(dim);
  std::vector <double> phixHat;
  std::vector <double> phiHat;
  std::vector <double> phixOld;
  std::vector <double> phiOld;
  double weightHat;

  std::vector< unsigned > sysDof;
  std::vector< double > Res;
  std::vector < double > Jac;

  RES->zero();
  KK->zero();

  std::vector < std::vector < std::vector <double > > > aP(3);

  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielType = msh->GetElementType(iel);
    short unsigned ielGroup = msh->GetElementGroup(iel);

    unsigned nDofs = msh->GetElementDofNumber(iel, solType);
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);

    unsigned nDofsAll = 2 * dim * nDofs + nDofsP;
    sysDof.resize(2 * nDofsAll);

    for(unsigned  k = 0; k < dim; k++) {
      solVc[k].resize(nDofs);

      solUim1[k].resize(nDofs);
      solVim1[k].resize(nDofs);
      solVim1Old[k].resize(nDofs);

      solUiOld[k].resize(nDofs);
      solViOld[k].resize(nDofs);

      solbU[k].resize(nDofs);
      solbV[k].resize(nDofs);

      sollU[k].resize(nDofs);
      sollV[k].resize(nDofs);

      xOld[k].resize(nDofs);
      x[k].resize(nDofs);
      xm1[k].resize(nDofs);
      xHat[k].resize(nDofs);

      DX[k].resize(nDofs);
    }
    solbP.resize(nDofsP);
    sollP.resize(nDofsP);

    for(unsigned  k = 0; k < dim; k++) {
      mResbU[k].assign(nDofs, 0.);
      mResbV[k].assign(nDofs, 0.);
      mReslU[k].assign(nDofs, 0.);
      mReslV[k].assign(nDofs, 0.);
    }
    mResbP.assign(nDofsP, 0.);
    mReslP.assign(nDofsP, 0.);

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, solType);

      for(unsigned  k = 0; k < dim; k++) {

        solVc[k][i] = (*sol->_Sol[solVcIndex[k]])(solDof);

        if(iext != 0) { // used to build the approximate time derivative coming form the previous iteration
          solUim1[k][i] = (*sol->_Sol[solUim1Index[k]])(solDof);
          solVim1[k][i] = (*sol->_Sol[solVim1Index[k]])(solDof);
          solVim1Old[k][i] = (*sol->_SolOld[solVim1Index[k]])(solDof);
        }
        else { // used at the zero iteration we use to build a approximate/delaied time derivative
          solUim1[k][i] = (*sol->_SolOld[solUim1Index[k]])(solDof);
          solVim1[k][i] = (*sol->_SolOld[solVim1Index[k]])(solDof);
          solVim1Old[k][i] = (*sol->_Sol[solV0OlderIndex[k]])(solDof);
        }

        solUiOld[k][i] = (*sol->_SolOld[solUiIndex[k]])(solDof);
        solViOld[k][i] = (*sol->_SolOld[solViIndex[k]])(solDof);

        solbU[k][i] = (*sol->_Sol[solbUIndex[k]])(solDof);
        solbV[k][i] = (*sol->_Sol[solbVIndex[k]])(solDof);
        sollU[k][i] = (*sol->_Sol[sollUIndex[k]])(solDof);
        sollV[k][i] = (*sol->_Sol[sollVIndex[k]])(solDof);

        sysDof[k * nDofs + i] = pdeSys->GetSystemDof(solbUIndex[k], solbUPdeIndex[k], i, iel);
        sysDof[(dim + k) * nDofs + i] = pdeSys->GetSystemDof(solbVIndex[k], solbVPdeIndex[k], i, iel);
        sysDof[nDofsAll + k * nDofs + i] = pdeSys->GetSystemDof(sollUIndex[k], sollUPdeIndex[k], i, iel);
        sysDof[nDofsAll + (dim + k) * nDofs + i] = pdeSys->GetSystemDof(sollVIndex[k], sollVPdeIndex[k], i, iel);
      }
    }

    for(unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);
      solbP[i] = (*sol->_Sol[solbPIndex])(solPDof);
      sollP[i] = (*sol->_Sol[sollPIndex])(solPDof);
      sysDof[2 * dim * nDofs + i ] = pdeSys->GetSystemDof(solbPIndex, solbPPdeIndex, i, iel);
      sysDof[nDofsAll + 2 * dim * nDofs + i ] = pdeSys->GetSystemDof(sollPIndex, sollPPdeIndex, i, iel);
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // local storage of coordinates
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // local to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        xHat[k][i] = (*msh->_topology->_Sol[k])(xDof);
        xOld[k][i] = xHat[k][i] + (*sol->_Sol[solDXIndex[k]])(xDof);
        xm1[k][i] = xOld[k][i] + solUim1[k][i] * dt;    // current deformed mesh at the previous iteration
        DX[k][i] = (*sol->_Sol[solDXIndex[k]])(xDof) + solbU[k][i] * dt;     // current deformed bmesh
        x[k][i] = xHat[k][i] + DX[k][i];     // current deformed bmesh
      }
    }


    double iRe = 1. / Re;


    bool elementIs2 = false;
    bool elementIs1or3 = false;
    std::vector<bool> nodeIsControlBoundary(nDofs, false);
    std::vector<bool> nodeIs1Or3(nDofs, false);
    for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
      int facename = -(msh->el->GetFaceElementIndex(iel, iface) + 1);

      if(facename == 2 || facename == 5) {
        elementIs2 = true;
        unsigned nve = msh->GetElementFaceDofNumber(iel, iface, solType);

        for(unsigned i = 0; i < nve; i++) {
          nodeIsControlBoundary[ msh->GetLocalFaceVertexIndex(iel, iface, i)] = true;
        }

        const unsigned faceType = msh->GetElementFaceType(iel, iface);
        unsigned faceDofs = msh->GetElementFaceDofNumber(iel, iface, solType);
        std::vector  < std::vector  <  double > > faceXHat(dim);
        for(int k = 0; k < dim; k++) {
          faceXHat[k].resize(faceDofs);
        }
        for(unsigned i = 0; i < faceDofs; i++) {
          unsigned inode = msh->GetLocalFaceVertexIndex(iel, iface, i);    // face-to-element local node mapping.
          for(unsigned k = 0; k < dim; k++) {
            faceXHat[k][i] =  xHat[k][inode];
          }
        }


        for(unsigned itype = 0; itype < solType + 1; itype++) {
          ProjectNodalToPolynomialCoefficients(aP[itype], xHat, ielType, itype);
        }

        for(unsigned ig = 0; ig  <  msh->_finiteElement[faceType][solType]->GetGaussPointNumber(); ig++) {

          std::vector < double> normalHat;
          std::vector < double> tauHat(dim);
          msh->_finiteElement[faceType][solType]->JacobianSur(faceXHat, ig, weightHat, phiHat, phixHat, normalHat);
          tauHat[0] = -normalHat[1];
          tauHat[1] = normalHat[0];

          std::vector< double > xgHat(dim, 0.); // physical coordinates of the face Gauss point
          for(unsigned i = 0; i < faceDofs; i++) {
            for(unsigned k = 0; k < dim; k++) {
              xgHat[k] += phiHat[i] * faceXHat[k][i];
            }
          }

          std::vector <double> xi;//local coordinates of the face gauss point with respect to iel
          GetClosestPointInReferenceElement(faceXHat, xgHat, ielType, xi);
          bool inverseMapping = GetInverseMapping(solType, ielType, aP, xgHat, xi, 100);

          double weightFake;
          msh->_finiteElement[ielType][solType]->Jacobian(xHat, xi, weightFake, phiHat, phixHat);
          std::vector < double > phixDotTauHat(nDofs, 0.);
          std::vector < adept::adouble > solDX_xHatgDotTauHatg(dim, 0.);

          for(unsigned i = 0; i < nDofs; i++) {
            for(unsigned  k = 0; k < dim; k++) {
              for(unsigned j = 0; j < dim; j++) {
                solDX_xHatgDotTauHatg[k] += DX[k][i] * phixHat[i * dim + j] * tauHat[j];
              }
            }
            for(unsigned j = 0; j < dim; j++) {
              phixDotTauHat[i] += phixHat[i * dim + j] * tauHat[j];
            }
          }

          for(unsigned i = 0; i < nDofs; i++) {
            for(unsigned  k = 0; k < dim; k++) {  //momentum equation in k
              mReslV[k][i] += -gammaU * phixDotTauHat[i] * solDX_xHatgDotTauHatg[k] * weightHat;
            }
          }
        }
      }
    }
    for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
      
      int facename = -(msh->el->GetFaceElementIndex(iel, iface) + 1);  
      if(facename == 1 || facename == 3) {
        elementIs1or3 = true;
        unsigned nve = msh->GetElementFaceDofNumber(iel, iface, solType);
        const unsigned felt = msh->GetElementFaceType(iel, iface);
        for(unsigned i = 0; i < nve; i++) {
          nodeIs1Or3[ msh->GetLocalFaceVertexIndex(iel, iface, i)] = true;

          unsigned inode = msh->GetLocalFaceVertexIndex(iel, iface, i);

          if(nodeIsControlBoundary[inode]) {
            double c = flc4hs(time - t0, t0);
            if(facename == 1) {
              //std::cout << "I am down\n";  
              mReslV[0][inode] += gammaU * H0 * M_PI / 4. * sin(M_PI / 4. * (- c * time));
            }
            else {
              //std::cout << "I am up\n";    
              mReslV[0][inode] += -gammaU * H0 * M_PI / 4. * sin(M_PI / 4. * (- c * time));
            }
          }
        }
      }
    }

    bool elementIsCorner = elementIs2 * elementIs1or3;
    //if(elementIsCorner) std::cout << BLUE << iel << " is corner cell\n" << BLACK;






    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielType][solType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielType][solType]->Jacobian(xm1, ig, weight, phi, phix);
      msh->_finiteElement[ielType][solType]->Jacobian(xOld, ig, weightOld, phiHat, phixHat);
      msh->_finiteElement[ielType][solType]->Jacobian(xHat, ig, weightHat, phiHat, phixHat);
      phiP = msh->_finiteElement[ielType][solPType]->GetPhi(ig);

      std::vector < double > solVcg(dim, 0);

      std::vector < double > solVim1g(dim, 0);
      std::vector < double > solVim1Oldg(dim, 0);

      std::vector < double > solUiOldg(dim, 0);
      std::vector < double > solViOldg(dim, 0);

      std::vector < adept::adouble > solbUg(dim, 0);
      std::vector < adept::adouble > solbVg(dim, 0);
      std::vector < adept::adouble > sollUg(dim, 0);
      std::vector < adept::adouble > sollVg(dim, 0);
      std::vector < double > xg(dim, 0);

      std::vector < std::vector < adept::adouble > > solbUxg(dim, std::vector < adept::adouble >(dim, 0.));
      std::vector < std::vector < adept::adouble > > solbVxg(dim, std::vector < adept::adouble >(dim, 0.));
      std::vector < std::vector < adept::adouble > > sollUxg(dim, std::vector < adept::adouble >(dim, 0.));
      std::vector < std::vector < adept::adouble > > sollVxg(dim, std::vector < adept::adouble >(dim, 0.));

      std::vector < std::vector < adept::adouble > > solbU_xHatg(dim, std::vector < adept::adouble >(dim, 0.));
      std::vector < std::vector < adept::adouble > > sollU_xHatg(dim, std::vector < adept::adouble >(dim, 0.));
      std::vector < std::vector < adept::adouble > > DX_xHatg(dim, std::vector < adept::adouble >(dim, 0.));


      for(unsigned i = 0; i < nDofs; i++) {
        for(unsigned  k = 0; k < dim; k++) {

          solVcg[k] += solVc[k][i] * phi[i];

          solVim1g[k] += solVim1[k][i] * phi[i];
          solVim1Oldg[k] += solVim1Old[k][i] * phi[i];

          solUiOldg[k] += solUiOld[k][i] * phi[i];
          solViOldg[k] += solViOld[k][i] * phi[i];

          solbUg[k] += solbU[k][i] * phi[i];
          solbVg[k] += solbV[k][i] * phi[i];
          sollUg[k] += sollU[k][i] * phi[i];
          sollVg[k] += sollV[k][i] * phi[i];

          xg[k] += x[k][i].value() * phi[i];
        }
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            solbUxg[k][j] += solbU[k][i] * phix[i * dim + j];
            solbVxg[k][j] += solbV[k][i] * phix[i * dim + j];
            sollUxg[k][j] += sollU[k][i] * phix[i * dim + j];
            sollVxg[k][j] += sollV[k][i] * phix[i * dim + j];

            solbU_xHatg[k][j] += solbU[k][i] * phixHat[i * dim + j];
            sollU_xHatg[k][j] += sollU[k][i] * phixHat[i * dim + j];
            DX_xHatg[k][j] += DX[k][i] * phixHat[i * dim + j];
          }
        }
      }



      double yg = xg[1] / 4.;

      adept::adouble solbPg = 0;
      adept::adouble sollPg = 0;
      for(unsigned i = 0; i < nDofsP; i++) {
        solbPg += phiP[i] * solbP[i];
        sollPg += phiP[i] * sollP[i];
      }

      std::vector < adept::adouble > lvNSVb(dim, 0.); // this is the Navier Stokes Equation mulitplied by the lagrange multiplier lv used for the variation of the jacobian


      for(unsigned  k = 0; k < dim; k++) {  //momentum equation in k
        for(unsigned j = 0; j < dim; j++) {  // second index j in each equation
          lvNSVb[k]   +=  iRe * sollVxg[k][j] * (solbVxg[k][j] + 0 * solbVxg[j][k]);
          lvNSVb[k]   +=  sollVg[k] * (solViOldg[j] - solUiOldg[j]) * solbVxg[k][j];
        }
        lvNSVb[k] += (1 - 0.1 * (iext == 0)) * (solVim1g[k] - solVim1Oldg[k]) / dt * sollVg[k] - solbPg * sollVxg[k][k];
      }

      adept::adouble divVlp = 0.;
      for(unsigned k = 0; k < dim; k++) {
        divVlp += solbVxg[k][k];
      }
      divVlp *= solbPg;


      // *** phiV_i loop ***
      for(unsigned i = 0; i < nDofs; i++) {
        std::vector < adept::adouble > ALEb(dim, 0.);
        std::vector < adept::adouble > ALEl(dim, 0.);

        std::vector < adept::adouble > NSVb(dim, 0.);
        std::vector < adept::adouble > NSVl(dim, 0.);

        for(unsigned  k = 0; k < dim; k++) {  //momentum equation in k
          for(unsigned j = 0; j < dim; j++) {  // second index j in each equation

            ALEb[k] += ((j == 0) || !oneDimDisp) * phixHat[i * dim + j] * (DX_xHatg[k][j] + 0 * DX_xHatg[j][k] /*- 1 * (j == k)*/);
            ALEl[k] += ((j == 0) || !oneDimDisp) * (phixHat[i * dim + j] * dt * (sollU_xHatg[k][j] + 0 * sollU_xHatg[j][k]) + betaU * phixHat[i * dim + j] * (solbUxg[k][j] + 0 * solbUxg[j][k]));

            NSVb[k]   +=  iRe * phix[i * dim + j] * (solbVxg[k][j] + 0 * solbVxg[j][k]);
            NSVl[k]   +=  iRe * phix[i * dim + j] * (sollVxg[k][j] + 0 * sollVxg[j][k]) + betaV * phix[i * dim + j] * (solbVxg[k][j] + 0 * solbVxg[j][k]);
            NSVb[k]   +=  phi[i] * (solViOldg[j] - solUiOldg[j]) * solbVxg[k][j];
            NSVl[k]   +=  sollVg[k] * (solViOldg[j] - solUiOldg[j]) * phix[i * dim + j];
          }
          NSVb[k] += (1 - 0.1 * (iext == 0)) * (solVim1g[k] - solVim1Oldg[k]) / dt * phi[i] - solbPg * phix[i * dim + k];
          ALEl[k] += alphaU * solbUg[k] * phi[i];
          NSVl[k] += -sollPg * phix[i * dim + k]  + alphaV * solbVg[k] * phi[i] + /*(ielGroup == 6)* */ (solbVg[k] - solVcg[k]) * phi[i];
        }

        for(unsigned  k = 0; k < dim; k++) {

          if(k == 0 || !oneDimDisp) {
            mResbU[k][i] += - ALEb[k] * weightHat;
          }
          else {
            mResbU[k][i] += solbU[k][i];
          }
          mResbV[k][i] += - NSVb[k] * weight;

          if(!nodeIsControlBoundary[i]) {
            if(k == 0 || !oneDimDisp) {
              mReslU[k][i] += - ALEl[k] * weightHat;
            }
            else {
              mReslU[k][i] += sollU[k][i];
            }
            mReslV[k][i] += - NSVl[k] * weight;

            //if(k == 0) mReslU[k][i] += - (lvNSVb[k] - divVlp) * (phix[i * dim + 0] * dt * (-1. + solbUxg[1][1] * dt) - phix[i * dim + 1] * dt * solbUxg[1][0] * dt) * weightOld;
            //else mReslU[k][i] += - (lvNSVb[k] - divVlp) * (phix[i * dim + 1] * dt * (-1. + solbUxg[0][0] * dt) - phix[i * dim + 0] * dt * solbUxg[0][1] * dt) * weightOld;
          }
          else {
            if(k == 0 || !oneDimDisp) {
              mReslU[k][i] += solbV[k][i] - solbU[k][i];
              mReslV[k][i] += - (NSVl[k] + ALEl[k]) * weight;
            }
            else {
              mReslU[k][i] += solbU[k][i];
              mReslV[k][i] += solbV[k][i];
            }
          }
        }
      } // end phiV_i loop

      // *** phiP_i loop ***
      for(unsigned i = 0; i < nDofsP; i++) {
        for(int k = 0; k < dim; k++) {
          mResbP[i] += - (-solbVxg[k][k]) * phiP[i]  * weight;
          mReslP[i] += - (-sollVxg[k][k]) * phiP[i]  * weight;
        }
      } // end phiP_i loop
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube mRes in double Res and store them in RES
    Res.resize(2 * nDofsAll);    //resize

    for(int i = 0; i < nDofs; i++) {
      for(unsigned  k = 0; k < dim; k++) {
        Res[ k * nDofs + i ] = -mReslU[k][i].value();
        Res[(k + dim) * nDofs + i ] = -mReslV[k][i].value();
        Res[nDofsAll +  k * nDofs + i ] = -mResbU[k][i].value();
        Res[nDofsAll + (k + dim) * nDofs + i ] = -mResbV[k][i].value();
      }
    }

    for(int i = 0; i < nDofsP; i++) {
      Res[ 2 * dim * nDofs + i ] = -mReslP[i].value();
      Res[nDofsAll + 2 * dim * nDofs + i ] = -mResbP[i].value();
    }

    RES->add_vector_blocked(Res, sysDof);

    //Extract and store the Jacobian

    Jac.resize(2 * nDofsAll * 2 * nDofsAll);
    // define the dependent variables

    for(unsigned  k = 0; k < dim; k++) {
      s.dependent(&mReslU[k][0], nDofs);
    }
    for(unsigned  k = 0; k < dim; k++) {
      s.dependent(&mReslV[k][0], nDofs);
    }
    s.dependent(&mReslP[0], nDofsP);

    for(unsigned  k = 0; k < dim; k++) {
      s.dependent(&mResbU[k][0], nDofs);
    }
    for(unsigned  k = 0; k < dim; k++) {
      s.dependent(&mResbV[k][0], nDofs);
    }
    s.dependent(&mResbP[0], nDofsP);






    // define the independent variables
    for(unsigned  k = 0; k < dim; k++) {
      s.independent(&solbU[k][0], nDofs);
    }
    for(unsigned  k = 0; k < dim; k++) {
      s.independent(&solbV[k][0], nDofs);
    }
    s.independent(&solbP[0], nDofsP);

    for(unsigned  k = 0; k < dim; k++) {
      s.independent(&sollU[k][0], nDofs);
    }
    for(unsigned  k = 0; k < dim; k++) {
      s.independent(&sollV[k][0], nDofs);
    }
    s.independent(&sollP[0], nDofsP);


    // get the and store jacobian matrix (row-major)
    s.jacobian(&Jac[0], true);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();
  KK->close();

  //KK->draw();
  // ***************** END ASSEMBLY *******************

}


void AssembleSystemZi(MultiLevelProblem & ml_prob) {
  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  TransientNonlinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<TransientNonlinearImplicitSystem> ("systemZi");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  double dt =  mlPdeSys->GetIntervalTime();
  double time =  mlPdeSys->GetTime();


  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual std::vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)


  //solution variable
  std::vector < unsigned > solbUIndex(dim);
  solbUIndex[0] = mlSol->GetIndex("bU1");
  solbUIndex[1] = mlSol->GetIndex("bU2");

  std::vector < unsigned > solbVIndex(dim);
  solbVIndex[0] = mlSol->GetIndex("bV1");
  solbVIndex[1] = mlSol->GetIndex("bV2");

  /*
    //solution variable
    std::vector < unsigned > solbUIndex(dim);
    solbUIndex[0] = mlSol->GetIndex("U1c");
    solbUIndex[1] = mlSol->GetIndex("U2c");

    std::vector < unsigned > solbVIndex(dim);
    solbVIndex[0] = mlSol->GetIndex("V1c");
    solbVIndex[1] = mlSol->GetIndex("V2c");*/



  //solution variable
  std::vector < unsigned > solUIndex(dim);
  solUIndex[0] = mlSol->GetIndex("U1i");
  solUIndex[1] = mlSol->GetIndex("U2i");


  std::vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("V1i");
  solVIndex[1] = mlSol->GetIndex("V2i");

  unsigned solType = mlSol->GetSolutionType(solVIndex[0]);

  unsigned solPIndex;
  solPIndex = mlSol->GetIndex("Pi");
  unsigned solPType = mlSol->GetSolutionType(solPIndex);

  std::vector < unsigned > solUPdeIndex(dim);
  solUPdeIndex[0] = mlPdeSys->GetSolPdeIndex("U1i");
  solUPdeIndex[1] = mlPdeSys->GetSolPdeIndex("U2i");


  std::vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("V1i");
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("V2i");


  unsigned solPPdeIndex;
  solPPdeIndex = mlPdeSys->GetSolPdeIndex("Pi");

  std::vector < std::vector < double > >  solbU(dim);
  std::vector < std::vector < double > >  solbV(dim);

  std::vector < std::vector < adept::adouble > >  solU(dim);
  std::vector < std::vector < double > >  solUOld(dim);

  std::vector < std::vector < adept::adouble > >  solV(dim);
  std::vector < std::vector < double > >  solVOld(dim);

  std::vector < adept::adouble >  solP;

  std::vector< std::vector < adept::adouble > > mResU(dim);
  std::vector< std::vector < adept::adouble > > mResV(dim);
  std::vector< adept::adouble > mResP;


  std::vector < unsigned > solDXIndex(dim);
  solDXIndex[0] = mlSol->GetIndex("DXi");
  solDXIndex[1] = mlSol->GetIndex("DYi");

//   std::vector < unsigned > solxHatIndex(dim);
//   solxHatIndex[0] = mlSol->GetIndex("xHat");
//   solxHatIndex[1] = mlSol->GetIndex("yHat");
  std::vector < std::vector < double > > xHat(dim);
  std::vector <double> phixHat;
  std::vector <double> phiHat;
  double weightHat;

  std::vector < std::vector < adept::adouble > > x(dim);
  unsigned xType = 2;

  std::vector <double> phi;
  std::vector <adept::adouble> phix;

  double* phiP;
  adept::adouble weight;

  std::vector< unsigned > sysDof;
  std::vector< double > Res;
  std::vector < double > Jac;

  RES->zero();
  KK->zero();

  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielType = msh->GetElementType(iel);

    unsigned nDofs = msh->GetElementDofNumber(iel, solType);
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);

    unsigned nDofsAll = 2 * dim * nDofs + nDofsP;
    sysDof.resize(nDofsAll);

    for(unsigned  k = 0; k < dim; k++) {
      solbU[k].resize(nDofs);
      solbV[k].resize(nDofs);
      solU[k].resize(nDofs);
      solV[k].resize(nDofs);
      solUOld[k].resize(nDofs);
      solVOld[k].resize(nDofs);
      x[k].resize(nDofs);
      xHat[k].resize(nDofs);
    }
    solP.resize(nDofsP);

    for(unsigned  k = 0; k < dim; k++) {
      mResU[k].assign(nDofs, 0.);
      mResV[k].assign(nDofs, 0.);
    }
    mResP.assign(nDofsP, 0.);

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, solType);

      for(unsigned  k = 0; k < dim; k++) {
        solbU[k][i] = (*sol->_Sol[solbUIndex[k]])(solDof);
        solbV[k][i] = (*sol->_Sol[solbVIndex[k]])(solDof);
        solU[k][i] = (*sol->_Sol[solUIndex[k]])(solDof);
        solUOld[k][i] = (*sol->_SolOld[solUIndex[k]])(solDof);
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solDof);
        solVOld[k][i] = (*sol->_SolOld[solVIndex[k]])(solDof);

        sysDof[k * nDofs + i] = pdeSys->GetSystemDof(solUIndex[k], solUPdeIndex[k], i, iel);
        sysDof[(dim + k) * nDofs + i] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);

      }
    }

    for(unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);
      solP[i] = (*sol->_Sol[solPIndex])(solPDof);
      sysDof[2 * dim * nDofs + i ] = pdeSys->GetSystemDof(solPIndex, solPPdeIndex, i, iel);
    }

    bool elementIs2 = false;
    bool elementIs1or3 = false;
    std::vector<bool> nodeIsControlBoundary(nDofs, false);
    std::vector<bool> nodeIs1Or3(nDofs, false);
    for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
      int facename = -(msh->el->GetFaceElementIndex(iel, iface) + 1);

      if(facename == 2 || facename == 5) {
        elementIs2 = true;
        unsigned nve = msh->GetElementFaceDofNumber(iel, iface, solType);
        const unsigned felt = msh->GetElementFaceType(iel, iface);
        for(unsigned i = 0; i < nve; i++) {
          nodeIsControlBoundary[ msh->GetLocalFaceVertexIndex(iel, iface, i)] = true;
        }
      }
      if(facename == 1 || facename == 3) {
        elementIs1or3 = true;
        unsigned nve = msh->GetElementFaceDofNumber(iel, iface, solType);
        const unsigned felt = msh->GetElementFaceType(iel, iface);
        for(unsigned i = 0; i < nve; i++) {
          nodeIs1Or3[ msh->GetLocalFaceVertexIndex(iel, iface, i)] = true;
        }
      }
    }

    bool elementIsCorner = elementIs2 * elementIs1or3;

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();


    // local storage of coordinates
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // local to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        xHat[k][i] = (*msh->_topology->_Sol[k])(xDof);
        x[k][i] = xHat[k][i] + (*sol->_Sol[solDXIndex[k]])(xDof) + solU[k][i] * dt;      // global extraction and local storage for the element coordinates
      }
    }


    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielType][solType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielType][solType]->Jacobian(x, ig, weight, phi, phix);
      msh->_finiteElement[ielType][solType]->Jacobian(xHat, ig, weightHat, phi, phixHat);
      phiP = msh->_finiteElement[ielType][solPType]->GetPhi(ig);

      std::vector < adept::adouble > solUg(dim, 0);
      std::vector < double > solUOldg(dim, 0);
      std::vector < adept::adouble > solVg(dim, 0);
      std::vector < double > solVOldg(dim, 0);

      std::vector < std::vector < adept::adouble > > solUxg(dim, std::vector < adept::adouble >(dim, 0.));
      std::vector < std::vector < adept::adouble > > solVxg(dim, std::vector < adept::adouble >(dim, 0.));
      std::vector < std::vector < adept::adouble > > x_xHatg(dim, std::vector < adept::adouble >(dim, 0.));


      for(unsigned i = 0; i < nDofs; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          solUg[k] += solU[k][i] * phi[i];
          solUOldg[k] += solUOld[k][i] * phi[i];
          solVg[k] += solV[k][i] * phi[i];
          solVOldg[k] += solVOld[k][i] * phi[i];
        }
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            solUxg[k][j] += solU[k][i] * phix[i * dim + j];
            solVxg[k][j] += solV[k][i] * phix[i * dim + j];
            x_xHatg[k][j] += x[k][i] * phixHat[i * dim + j];
          }
        }
      }

      adept::adouble solPg = 0;
      for(unsigned i = 0; i < nDofsP; i++) {
        solPg += phiP[i] * solP[i];
      }
      double iRe = 1. / Re;

      // *** phiV_i loop ***
      for(unsigned i = 0; i < nDofs; i++) {
        std::vector < adept::adouble > ALE(dim, 0.);
        std::vector < adept::adouble > NSV(dim, 0.);

        for(unsigned  k = 0; k < dim; k++) {  //momentum equation in k
          for(unsigned j = 0; j < dim; j++) {  // second index j in each equation

            ALE[k] += ((j == 0) || !oneDimDisp) * phixHat[i * dim + j] * (x_xHatg[k][j] + 0 * x_xHatg[j][k]  - 1 * (j == k));

            NSV[k] += iRe * phix[i * dim + j] * (solVxg[k][j] + 0 * solVxg[j][k]);
            NSV[k] +=  phi[i] * (solVOldg[j] - solUOldg[j]) * solVxg[k][j];
          }
          NSV[k] += (solVg[k] - solVOldg[k]) / dt * phi[i] - solPg * phix[i * dim + k];
        }
        for(unsigned  k = 0; k < dim; k++) {
          if(!nodeIsControlBoundary[i]) {
            if(k == 0 || !oneDimDisp) {
              mResU[k][i] += - ALE[k] * weightHat;
            }
            else {
              mResU[k][i] += solU[k][i];
            }
            mResV[k][i] += - NSV[k] * weight;
          }
          else {
            if(k == 0 || !oneDimDisp) {
              mResU[k][i] += solU[k][i] - solbU[k][i];
            }
            else {
              mResU[k][i] += solU[k][i];
            }
            mResV[k][i] += solV[k][i] - solbV[k][i];
          }
        }
      } // end phiV_i loop

      // *** phiP_i loop ***
      for(unsigned i = 0; i < nDofsP; i++) {
        for(int k = 0; k < dim; k++) {
          mResP[i] += - (-solVxg[k][k]) * phiP[i]  * weight;
        }
      } // end phiP_i loop
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube mRes in double Res and store them in RES
    Res.resize(nDofsAll);    //resize

    for(int i = 0; i < nDofs; i++) {
      for(unsigned  k = 0; k < dim; k++) {
        Res[ k * nDofs + i ] = -mResU[k][i].value();
      }
      for(unsigned  k = 0; k < dim; k++) {
        Res[(dim + k) * nDofs + i ] = -mResV[k][i].value();
      }
    }

    for(int i = 0; i < nDofsP; i++) {
      Res[ 2 * dim * nDofs + i ] = -mResP[i].value();
    }

    RES->add_vector_blocked(Res, sysDof);

    //Extract and store the Jacobian

    Jac.resize(nDofsAll * nDofsAll);
    // define the dependent variables


    for(unsigned  k = 0; k < dim; k++) {
      s.dependent(&mResU[k][0], nDofs);
    }
    for(unsigned  k = 0; k < dim; k++) {
      s.dependent(&mResV[k][0], nDofs);
    }
    s.dependent(&mResP[0], nDofsP);

    // define the independent variables
    for(unsigned  k = 0; k < dim; k++) {
      s.independent(&solU[k][0], nDofs);
    }
    for(unsigned  k = 0; k < dim; k++) {
      s.independent(&solV[k][0], nDofs);
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

// KK->draw();
  // ***************** END ASSEMBLY *******************

}



void AssembleManifactureSolution(MultiLevelProblem & ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  TransientNonlinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<TransientNonlinearImplicitSystem> ("ManSol");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  double dt =  mlPdeSys->GetIntervalTime();
  double time =  mlPdeSys->GetTime();


  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual std::vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  std::vector < unsigned > solUIndex(dim);
  solUIndex[0] = mlSol->GetIndex("U1c");
  solUIndex[1] = mlSol->GetIndex("U2c");


  std::vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("V1c");
  solVIndex[1] = mlSol->GetIndex("V2c");

  unsigned solType = mlSol->GetSolutionType(solVIndex[0]);

  unsigned solPIndex;
  solPIndex = mlSol->GetIndex("Pc");
  unsigned solPType = mlSol->GetSolutionType(solPIndex);

  std::vector < unsigned > solUPdeIndex(dim);
  solUPdeIndex[0] = mlPdeSys->GetSolPdeIndex("U1c");
  solUPdeIndex[1] = mlPdeSys->GetSolPdeIndex("U2c");


  std::vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("V1c");
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("V2c");


  unsigned solPPdeIndex;
  solPPdeIndex = mlPdeSys->GetSolPdeIndex("Pc");

  std::vector < std::vector < adept::adouble > >  solU(dim);
  std::vector < std::vector < double > >  solUOld(dim);

  std::vector < std::vector < adept::adouble > >  solV(dim);
  std::vector < std::vector < double > >  solVOld(dim);

  std::vector < adept::adouble >  solP;

  std::vector< std::vector < adept::adouble > > mResU(dim);
  std::vector< std::vector < adept::adouble > > mResV(dim);
  std::vector< adept::adouble > mResP;

  std::vector < std::vector < double > > xHat(dim);
  std::vector <double> phixHat;
  std::vector <double> phiHat;
  double weightHat;

  std::vector < unsigned > solDXIndex(dim);
  solDXIndex[0] = mlSol->GetIndex("DXc");
  solDXIndex[1] = mlSol->GetIndex("DYc");

  std::vector < std::vector < adept::adouble > > x(dim);
  unsigned xType = 2;

  std::vector <double> phi;
  std::vector <adept::adouble> phix;

  double* phiP;
  adept::adouble weight;

  std::vector< unsigned > sysDof;
  std::vector< double > Res;
  std::vector < double > Jac;

  RES->zero();
  KK->zero();

  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielType = msh->GetElementType(iel);

    unsigned nDofs = msh->GetElementDofNumber(iel, solType);
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);

    unsigned nDofsAll = 2 * dim * nDofs + nDofsP;
    sysDof.resize(nDofsAll);

    for(unsigned  k = 0; k < dim; k++) {
      solU[k].resize(nDofs);
      solV[k].resize(nDofs);
      solUOld[k].resize(nDofs);
      solVOld[k].resize(nDofs);
      x[k].resize(nDofs);
      xHat[k].resize(nDofs);
    }
    solP.resize(nDofsP);

    for(unsigned  k = 0; k < dim; k++) {
      mResU[k].assign(nDofs, 0.);
      mResV[k].assign(nDofs, 0.);
    }
    mResP.assign(nDofsP, 0.);

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, solType);

      for(unsigned  k = 0; k < dim; k++) {
        solU[k][i] = (*sol->_Sol[solUIndex[k]])(solDof);
        solUOld[k][i] = (*sol->_SolOld[solUIndex[k]])(solDof);
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solDof);
        solVOld[k][i] = (*sol->_SolOld[solVIndex[k]])(solDof);

        sysDof[k * nDofs + i] = pdeSys->GetSystemDof(solUIndex[k], solUPdeIndex[k], i, iel);
        sysDof[(dim + k) * nDofs + i] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);

      }
    }

    for(unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);
      solP[i] = (*sol->_Sol[solPIndex])(solPDof);
      sysDof[2 * dim * nDofs + i ] = pdeSys->GetSystemDof(solPIndex, solPPdeIndex, i, iel);
    }

    bool elementIs2 = false;
    bool elementIs1or3 = false;
    std::vector<bool> nodeIsControlBoundary(nDofs, false);
    std::vector<bool> nodeIs1Or3(nDofs, false);
    for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
      int facename = -(msh->el->GetFaceElementIndex(iel, iface) + 1);

      if(facename == 2 || facename == 5) {
        elementIs2 = true;
        unsigned nve = msh->GetElementFaceDofNumber(iel, iface, solType);
        const unsigned felt = msh->GetElementFaceType(iel, iface);
        for(unsigned i = 0; i < nve; i++) {
          nodeIsControlBoundary[ msh->GetLocalFaceVertexIndex(iel, iface, i)] = true;
        }
      }
      if(facename == 1 || facename == 3) {
        elementIs1or3 = true;
        unsigned nve = msh->GetElementFaceDofNumber(iel, iface, solType);
        const unsigned felt = msh->GetElementFaceType(iel, iface);
        for(unsigned i = 0; i < nve; i++) {
          nodeIs1Or3[ msh->GetLocalFaceVertexIndex(iel, iface, i)] = true;
        }
      }
    }

    bool elementIsCorner = elementIs2 * elementIs1or3;

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();


    // local storage of coordinates
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // local to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        xHat[k][i] = (*msh->_topology->_Sol[k])(xDof);
        x[k][i] = xHat[k][i] + (*sol->_Sol[solDXIndex[k]])(xDof) + solU[k][i] * dt;       // global extraction and local storage for the element coordinates
      }
    }


    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielType][solType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielType][solType]->Jacobian(x, ig, weight, phi, phix);
      msh->_finiteElement[ielType][solType]->Jacobian(xHat, ig, weightHat, phi, phixHat);
      phiP = msh->_finiteElement[ielType][solPType]->GetPhi(ig);

      std::vector < adept::adouble > solUg(dim, 0);
      std::vector < double > solUOldg(dim, 0);
      std::vector < adept::adouble > solVg(dim, 0);
      std::vector < double > solVOldg(dim, 0);

      std::vector < std::vector < adept::adouble > > solUxg(dim, std::vector < adept::adouble >(dim, 0.));
      std::vector < std::vector < adept::adouble > > solVxg(dim, std::vector < adept::adouble >(dim, 0.));
      std::vector < std::vector < adept::adouble > > x_xHatg(dim, std::vector < adept::adouble >(dim, 0.));


      for(unsigned i = 0; i < nDofs; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          solUg[k] += solU[k][i] * phi[i];
          solUOldg[k] += solUOld[k][i] * phi[i];
          solVg[k] += solV[k][i] * phi[i];
          solVOldg[k] += solVOld[k][i] * phi[i];
        }
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            solUxg[k][j] += solU[k][i] * phix[i * dim + j];
            solVxg[k][j] += solV[k][i] * phix[i * dim + j];
            x_xHatg[k][j] += x[k][i] * phixHat[i * dim + j];
          }
        }
      }

      adept::adouble solPg = 0;
      for(unsigned i = 0; i < nDofsP; i++) {
        solPg += phiP[i] * solP[i];
      }
      double iRe = 1. / Re;

      // *** phiV_i loop ***
      for(unsigned i = 0; i < nDofs; i++) {
        std::vector < adept::adouble > ALE(dim, 0.);
        std::vector < adept::adouble > NSV(dim, 0.);

        for(unsigned  k = 0; k < dim; k++) {  //momentum equation in k
          for(unsigned j = 0; j < dim; j++) {  // second index j in each equation

            ALE[k] += ((j == 0) || !oneDimDisp) * phixHat[i * dim + j] * (x_xHatg[k][j] + 0 * x_xHatg[j][k]  - 1 * (j == k));

            NSV[k] += iRe * phix[i * dim + j] * (solVxg[k][j] + 0 * solVxg[j][k]);
            NSV[k] +=  phi[i] * (solVOldg[j] - solUOldg[j]) * solVxg[k][j];
          }
          NSV[k] += (solVg[k] - solVOldg[k]) / dt * phi[i] - solPg * phix[i * dim + k];
        }
        for(unsigned  k = 0; k < dim; k++) {

          if(k == 0 || !oneDimDisp) mResU[k][i] += - ALE[k] * weightHat;
          else mResU[k][i] += solU[k][i];

          if(!nodeIsControlBoundary[i]) {
            mResV[k][i] += - NSV[k] * weight;
          }
          else {
            mResV[k][i] += solV[k][i] - solU[k][i];
          }
        }
      } // end phiV_i loop

      // *** phiP_i loop ***
      for(unsigned i = 0; i < nDofsP; i++) {
        for(int k = 0; k < dim; k++) {
          mResP[i] += - (-solVxg[k][k]) * phiP[i]  * weight;
        }
      } // end phiP_i loop
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube mRes in double Res and store them in RES
    Res.resize(nDofsAll);    //resize

    for(int i = 0; i < nDofs; i++) {
      for(unsigned  k = 0; k < dim; k++) {
        Res[ k * nDofs + i ] = -mResU[k][i].value();
      }
      for(unsigned  k = 0; k < dim; k++) {
        Res[(dim + k) * nDofs + i ] = -mResV[k][i].value();
      }
    }

    for(int i = 0; i < nDofsP; i++) {
      Res[ 2 * dim * nDofs + i ] = -mResP[i].value();
    }

    RES->add_vector_blocked(Res, sysDof);

    //Extract and store the Jacobian

    Jac.resize(nDofsAll * nDofsAll);
    // define the dependent variables


    for(unsigned  k = 0; k < dim; k++) {
      s.dependent(&mResU[k][0], nDofs);
    }
    for(unsigned  k = 0; k < dim; k++) {
      s.dependent(&mResV[k][0], nDofs);
    }
    s.dependent(&mResP[0], nDofsP);

    // define the independent variables
    for(unsigned  k = 0; k < dim; k++) {
      s.independent(&solU[k][0], nDofs);
    }
    for(unsigned  k = 0; k < dim; k++) {
      s.independent(&solV[k][0], nDofs);
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

// KK->draw();
  // ***************** END ASSEMBLY *******************

}



void GetError(MultiLevelProblem & ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled


  char uName[10];
  char vName[10];
  sprintf(uName, "V1%d", iext);
  sprintf(vName, "V2%d", iext);

  //char sName[10];
  //sprintf(sName, "top%d", iext);

  //  extract pointers to the several objects that we are going to use
  //TransientLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<TransientLinearImplicitSystem> (sName);   // pointer to the linear implicit system named "Poisson"
  TransientLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<TransientLinearImplicitSystem> ("systemZi");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable

  double time =  mlPdeSys->GetTime();

  std::vector<unsigned> VcIndex(dim);
  VcIndex[0] = mlSol->GetIndex("V1c");
  VcIndex[1] = mlSol->GetIndex("V2c");

  std::vector<unsigned> VIndex(dim);
  VIndex[0] = mlSol->GetIndex(uName);
  VIndex[1] = mlSol->GetIndex(vName);

  unsigned solType = mlSol->GetSolutionType(VcIndex[0]);

  std::vector<std::vector < double >> solVc(dim);
  std::vector<std::vector < double >> solV(dim);

  std::vector < std::vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector <double> phi, gradPhi;  // local test function for velocity
  double weight;

  double  localIntegral[5] = {0, 0, 0, 0, 0};



  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {



    unsigned group = msh->GetElementGroup(iel);
    if(group == group) {

      short unsigned ielGeom = msh->GetElementType(iel);
      unsigned nDofs = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
      for(unsigned  k = 0; k < dim; k++) {
        x[k].resize(nDofs);
      }
      for(unsigned k = 0; k < dim; k++) {
        solV[k].resize(nDofs);
        solVc[k].resize(nDofs);
      }

      // local storage of global mapping and solution
      for(unsigned i = 0; i < nDofs; i++) {
        unsigned solDof = msh->GetSolutionDof(i, iel, solType);
        for(unsigned k = 0; k < dim; k++) {
          solVc[k][i] = (*sol->_Sol[VcIndex[k]])(solDof);
          solV[k][i] = (*sol->_Sol[VIndex[k]])(solDof);
        }
      }

      // local storage of coordinates
      for(unsigned i = 0; i < nDofs; i++) {
        unsigned xDof  = msh->GetSolutionDof(i, iel, xType);
        for(unsigned k = 0; k < dim; k++) {
          x[k][i] = (*msh->_topology->_Sol[k])(xDof);
        }
      }



      // *** Gauss point loop ***
      for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {
        msh->_finiteElement[ielGeom][solType]->Jacobian(x, ig, weight, phi, gradPhi);

        std::vector<double> Vg(dim, 0.);
        std::vector<double> Vcg(dim, 0.);

        for(unsigned i = 0; i < nDofs; i++) {
          for(unsigned k = 0; k < dim; k++) {
            Vg[k] += solV[k][i] * phi[i];
            Vcg[k] += solVc[k][i] * phi[i];
          }
        }

        localIntegral[0] += Vcg[0] * weight;
        localIntegral[1] += Vcg[1] * weight;
        localIntegral[2] += Vg[0] * weight;
        localIntegral[3] += Vg[1] * weight;
        localIntegral[4] += ((Vg[0] - Vcg[0]) * (Vg[0] - Vcg[0]) + (Vg[1] - Vcg[1]) * (Vg[1] - Vcg[1])) * weight;

      } // end gauss point loop

    }
  } //end element loop for each process

  double globalIntegral[5] = {0., 0., 0., 0., 0.};

  MPI_Reduce(localIntegral, globalIntegral, 5, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


  if(iproc == 0) {

    std::string filename = "SolutionIntegral" + std::to_string(iext) + ".txt";
    std::ofstream fout;
    if(cleanFile) {
      fout.open(filename.c_str(), std::ios::out);
      if(iext == numberOfIterations - 1) cleanFile = false;
    }
    else {
      fout.open(filename.c_str(), std::ios::app);
    }
    fout << time << " " << globalIntegral[0] << " " << globalIntegral[1] << " " << globalIntegral[2] << " " << globalIntegral[3] << " " << globalIntegral[4] << std::endl;
    fout.close();
  }


// ***************** END ASSEMBLY *******************

}



















double flc4hs(double const & x, double const & eps) {

  double r = x / eps;
  if(r < -1) {
    return 0.;
  }
  else if(r < 1.) {
    double r2 = r * r;
    double r3 = r * r2;
    double r5 = r3 * r2;
    double r7 = r5 * r2;
    double r9 = r7 * r2;
    return (128. + 315. * r - 420. * r3 + 378. * r5 - 180. * r7 + 35. * r9) / 256.;
  }
  else {
    return 1.;
  }
}

double dflc4hs(double const & x, double const & eps) {

  double r = x / eps;
  if(r < -1) {
    return 0.;
  }
  else if(r < 1.) {
    double r2 = r * r;
    double r4 = r2 * r2;
    double r6 = r4 * r2;
    double r8 = r6 * r2;
    return (315. - 3. * 420. * r2 + 5 * 378. * r4 - 7. * 180. * r6 + 9.  * 35. * r8) / (256. * eps);
  }
  else {
    return 0.;
  }
}










