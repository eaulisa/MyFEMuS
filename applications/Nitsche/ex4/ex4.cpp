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
#include "slepceps.h"

#include "LinearImplicitSystem.hpp"
#include "Marker.hpp"
#include "Line.hpp"

#include "Fluid.hpp"
#include "Solid.hpp"
#include "Parameter.hpp"

#include "../NewDraft/NewDraft.hpp"
// #include <eigen3/Eigen/Dense>
// #include <eigen3/Eigen/Eigenvalues>
// #include <eigen3/unsupported/Eigen/KroneckerProduct>
// #include </usr/include/eigen3/Eigen/src/Core/util/DisableStupidWarnings.h>
// #include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <algorithm>    // std::minmax_element



using namespace femus;

Line* line3;
Line* lineI;

unsigned DIM = 2;
double eps = 0.5;
double deps;
double a0;
double a1;
double a3;
double a5;
double a7;
double a9;


void PrintMatlabMatrix(Eigen::MatrixXd &A);
void AssembleNitscheProblem_AD(MultiLevelProblem& mlProb);

void BuildFlag(MultiLevelSolution& mlSol);
void GetInterfaceElementEigenvalues(MultiLevelSolution& mlSol);
void GetInterfaceElementEigenvaluesAD(MultiLevelSolution& mlSol);

void GetParticleWeights(MultiLevelSolution& mlSol);



bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = false; //dirichlet
  if(facename == 1)  dirichlet = true;
  value = 0.;
  return dirichlet;
}


int main(int argc, char** args) {

  if(DIM != 2 && DIM != 3) {
    std::cout << "Wrong Dimension!" << std::endl;
    return 0;
  }

  // init Petsc-MPI communicator

  SlepcInitialize(&argc, &args, PETSC_NULL, PETSC_NULL);
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;

  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;

  unsigned nx = 11; // this should always be a odd number
  unsigned ny = 11; // this should always be a odd number
  unsigned nz = 1;

  const double length = 1.;
  const double lengthx = 0.5;


// TRI3-6-7 not working????
  if(DIM == 2) {
    mlMsh.GenerateCoarseBoxMesh(nx, ny, 0, -lengthx / 2, lengthx / 2, -length / 2, length / 2, 0., 0., QUAD9, "seventh");
  }
  else if(DIM == 3) {
    nz = ny;
    mlMsh.GenerateCoarseBoxMesh(nx, ny, nz, -lengthx / 2, lengthx / 2, -length / 2, length / 2, -length / 2, length / 2, HEX27, "seventh");
  }

  double Hx = lengthx / nx;
  double Hy = length / ny;
  double H = sqrt(Hx * Hx + Hy * Hy);

  deps = H * eps; // eps1
  a0 = 0.5; // 128./256.;
  a1 = pow(deps, -1.) * 1.23046875; // 315/256.;
  a3 = -pow(deps, -3.) * 1.640625; //420./256.;
  a5 = pow(deps, -5.) * 1.4765625; // 378./256.;
  a7 = -pow(deps, -7.) * 0.703125; // 180./256.;
  a9 = pow(deps, -9.) * 0.13671875; // 35./256.;



  double Lref = 1.;
  double Uref = 1.;
  double rho1 = 1000;
  double rho2 = 10;
  double mu1 = 1.e-03;
  double mu2 = 1.e-04;


  Parameter par(Lref, Uref);

  // Generate Solid Object
  Fluid fluid1(par, mu1, rho1, "Newtonian");
  Fluid fluid2(par, mu2, rho2, "Newtonian");

  //mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , NULL);
  mlMsh.PrintInfo();

  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol(&mlMsh);

  FEOrder femOrder = SECOND;

  mlSol.AddSolution("VX1", LAGRANGE, femOrder);
  mlSol.AddSolution("VY1", LAGRANGE, femOrder);
  if(DIM == 3) mlSol.AddSolution("VZ1", LAGRANGE, femOrder);
  mlSol.AddSolution("P1", DISCONTINUOUS_POLYNOMIAL, ZERO, 2);


  mlSol.AddSolution("VX2", LAGRANGE, femOrder);
  mlSol.AddSolution("VY2", LAGRANGE, femOrder);
  if(DIM == 3) mlSol.AddSolution("VZ2", LAGRANGE, femOrder);
  mlSol.AddSolution("P2", DISCONTINUOUS_POLYNOMIAL, ZERO, 2);



  mlSol.AddSolution("eflag", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
  mlSol.AddSolution("nflag", LAGRANGE, femOrder, 0, false);

  mlSol.AddSolution("CM1", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
  mlSol.AddSolution("CM2", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);

  mlSol.AddSolution("CL1", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
  mlSol.AddSolution("CL2", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);


  mlSol.Initialize("All");

  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("All");



  MultiLevelProblem ml_prob(&mlSol);

  ml_prob.parameters.set<Fluid> ("Fluid1") = fluid1;
  ml_prob.parameters.set<Fluid> ("Fluid2") = fluid2;


  // ******* Add FEM system to the MultiLevel problem *******
  LinearImplicitSystem& system = ml_prob.add_system < LinearImplicitSystem > ("Nitsche");

  // add solution "u" to system
  system.AddSolutionToSystemPDE("VX1");
  system.AddSolutionToSystemPDE("VY1");
  if(DIM == 3) system.AddSolutionToSystemPDE("VZ1");
// system.AddSolutionToSystemPDE("P1");

  system.AddSolutionToSystemPDE("VX2");
  system.AddSolutionToSystemPDE("VY2");
  if(DIM == 3) system.AddSolutionToSystemPDE("VZ2");
// system.AddSolutionToSystemPDE("P2");


  // attach the assembling function to system
  system.SetAssembleFunction(AssembleNitscheProblem_AD);

  // time loop parameter
  system.SetMaxNumberOfLinearIterations(1);

  system.init();

  //BEGIN init particles

  std::vector < std::vector <double> > xp;
  std::vector <double> wp;
  std::vector <double> dist;
  std::vector < MarkerType > markerType;

  

//// INIT Rectange Particle Initilization 

  // inner bulk solid markers + outer shell fluid markers
//   double Ls = 0.2;
//   double Hs = 3 * Ls;
//   double Lf = 0.4; 
//   double Hf = 4 * Hs;
//   unsigned rows = 30;
//   std::vector < double> xc = {-0.1, -0.5};
//   unsigned cols = (ceil((Ls / Hs) * (rows + 0.5) - 1)) ; // ensures dx~=dy in rectangle.
//   unsigned nbl = 3; // odd number
//  
//   double dL = Hs / rows;
//   double DB = 0.1 * dL;
//  
//   InitRectangleParticle(DIM, Ls, Hs, Lf, dL, DB, nbl, xc, markerType, xp, wp, dist);
// 
// 
//   unsigned solType = 2;
//   line3 = new Line(xp, wp, dist, markerType, mlSol.GetLevel(numberOfUniformLevels - 1), solType);
// 
//   std::vector < std::vector < std::vector < double > > >  line3Points(1);
//   line3->GetLine(line3Points[0]);
//   PrintLine(DEFAULT_OUTPUTDIR, "SolidMarkers", line3Points, 0);
// 
// 
//   //interface markers
// 
//   unsigned FI = 1;
//   std::vector < std::vector < std::vector < double > > > T;
//   InitRectangleInterface(DIM, Ls, Hs, dL, DB, FI, xc, markerType, xp, T);
//   
// 
// 
//   lineI = new Line(xp, T, markerType, mlSol.GetLevel(numberOfUniformLevels - 1), solType);
// 
//   std::vector < std::vector < std::vector < double > > > lineIPoints(1);
//   lineI->GetLine(lineIPoints[0]);
//   PrintLine(DEFAULT_OUTPUTDIR, "interfaceMarkers", lineIPoints, 0);

 
  
//// END Rectange Particle Initilization 
  
  
  
  //// INIT Ball Marker Initilization
  
  std::vector<double> VxL = { - 0.5 * lengthx, -0.5 * length, -0.5 * length };
  std::vector<double> VxR = {  0.5 * lengthx,  0.5 * length, 0.5 * length };

  double xc = 0.;
  double yc = 0.;
  double zc = 0.;
  double R = 0.125;
  double Rmax = 0.225;
  double DR = H / 10.;
  unsigned nbl = 5;
  unsigned FI = 5;
  std::vector < double> Xc = {xc, yc, zc};

  InitBallVolumeParticles(DIM, VxL, VxR, Xc, markerType, R, Rmax, DR, nbl, FI, xp, wp, dist);


  unsigned solType = 2;
  line3 = new Line(xp, wp, dist, markerType, mlSol.GetLevel(numberOfUniformLevels - 1), solType);

  std::vector < std::vector < std::vector < double > > >  line3Points(1);
  line3->GetLine(line3Points[0]);
  PrintLine(DEFAULT_OUTPUTDIR, "bulk3", line3Points, 0);


  ///interface stuff

  std::vector < std::vector < std::vector < double > > > T;

  InitBallInterfaceParticles(DIM, R, DR, FI, Xc, markerType, xp, T);

  lineI = new Line(xp, T, markerType, mlSol.GetLevel(numberOfUniformLevels - 1), solType);

  std::vector < std::vector < std::vector < double > > > lineIPoints(1);
  lineI->GetLine(lineIPoints[0]);
  PrintLine(DEFAULT_OUTPUTDIR, "interfaceLine", lineIPoints, 0);
//   
  
 //// END Ball Marker Initilization
  
  BuildFlag(mlSol);
  GetParticleWeights(mlSol);
  GetInterfaceElementEigenvalues(mlSol);
  
  

  

  system.MGsolve();

  // ******* Print solution *******
  mlSol.SetWriter(VTK);
  mlSol.GetWriter()->SetDebugOutput(true);

  std::vector<std::string> print_vars;
  print_vars.push_back("All");


  mlSol.GetWriter()->Write("./output", "biquadratic", print_vars, 0);

//   std::vector<std::string> mov_vars1;
//   mov_vars1.push_back("DX1");
//   mov_vars1.push_back("DY1");
//   if(DIM == 3) mov_vars1.push_back("DZ1");
//   mlSol.GetWriter()->SetMovingMesh(mov_vars1);
//   mlSol.GetWriter()->Write("./outputD1", "linear", print_vars, 0);
//
//
//   std::vector<std::string> mov_vars2;
//   mov_vars2.push_back("DX2");
//   mov_vars2.push_back("DY2");
//   if(DIM == 3) mov_vars2.push_back("DZ2");
//   mlSol.GetWriter()->SetMovingMesh(mov_vars2);
//   mlSol.GetWriter()->Write("./outputD2", "linear", print_vars, 0);


  ml_prob.clear();

  return 0;
}


void AssembleNitscheProblem_AD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object


  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("Nitsche");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble(); // We have different level of meshes. we assemble the problem on the specified one.

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*             KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*           RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  double rho1 = ml_prob.parameters.get<Solid> ("Fluid1").get_density();
  double rho2 = ml_prob.parameters.get<Solid> ("Fluid2").get_density();

  double mu1 = ml_prob.parameters.get<Solid> ("Fluid1").get_lame_shear_modulus();

  double mu2 = ml_prob.parameters.get<Solid> ("Fluid2").get_lame_shear_modulus();


  double g[DIM] = {0., -1., 0.};

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  std::vector< unsigned > solV1Index(dim);
  solV1Index[0] = mlSol->GetIndex("VX1");
  solV1Index[1] = mlSol->GetIndex("VY1");
  if(dim == 3)solV1Index[2] = mlSol->GetIndex("VZ1");

  std::vector< unsigned > solV2Index(dim);
  solV2Index[0] = mlSol->GetIndex("VX2");
  solV2Index[1] = mlSol->GetIndex("VY2");
  if(dim == 3)solV2Index[2] = mlSol->GetIndex("VZ2");

  unsigned solVType = mlSol->GetSolutionType(solV1Index[0]);

  std::vector< unsigned > solV1PdeIndex(dim);
  solV1PdeIndex[0] = mlPdeSys->GetSolPdeIndex("VX1");
  solV1PdeIndex[1] = mlPdeSys->GetSolPdeIndex("VY1");
  if(dim == 3) solV1PdeIndex[2] = mlPdeSys->GetSolPdeIndex("VZ1");

  std::vector< unsigned > solV2PdeIndex(dim);
  solV2PdeIndex[0] = mlPdeSys->GetSolPdeIndex("VX2");
  solV2PdeIndex[1] = mlPdeSys->GetSolPdeIndex("VY2");
  if(dim == 3) solV2PdeIndex[2] = mlPdeSys->GetSolPdeIndex("VZ2");

  std::vector < std::vector < adept::adouble > > solV1(dim); // local solution
  std::vector < std::vector < adept::adouble > > solV2(dim); // local solution

  unsigned CMIndex[2];
  unsigned CLIndex[2];

  CMIndex[0] = mlSol->GetIndex("CM1");
  CMIndex[1] = mlSol->GetIndex("CM2");

  CLIndex[0] = mlSol->GetIndex("CL1");
  CLIndex[1] = mlSol->GetIndex("CL2");

  unsigned eflagIndex = mlSol->GetIndex("eflag");
  unsigned nflagIndex = mlSol->GetIndex("nflag");

  vector < unsigned >  nodeFlag; // local solution

  vector < vector < double > > x(dim);    // local coordinates. x is now dim x m matrix.
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  double weight; // gauss point weight

  std::vector< std::vector< adept::adouble > > aResV1(dim); // local redidual vector
  std::vector< std::vector< adept::adouble > > aResV2(dim); // local redidual vector

  vector< unsigned > l2GMap; // local to global mapping
  vector< double > Res; // local redidual vector
  vector < double > Jac;

  KK->zero(); // Set to zero all the entries of the Global Matrix
  RES->zero(); // Set to zero all the entries of the Global Residual

  std::vector < std::vector < std::vector <double > > > aP(3);

//   std::vector<Marker*> particle1 = line1->GetParticles();
//   std::vector<unsigned> markerOffset1 = line1->GetMarkerOffset();
//   unsigned imarker1 = markerOffset1[iproc];
//
//   std::vector<Marker*> particle2 = line2->GetParticles();
//   std::vector<unsigned> markerOffset2 = line2->GetMarkerOffset();
//   unsigned imarker2 = markerOffset2[iproc];

  std::vector<Marker*> particle3 = line3->GetParticles();
  std::vector<unsigned> markerOffset3 = line3->GetMarkerOffset();
  unsigned imarker3 = markerOffset3[iproc];

  std::vector<Marker*> particleI = lineI->GetParticles();
  std::vector<unsigned> markerOffsetI = lineI->GetMarkerOffset();
  unsigned imarkerI = markerOffsetI[iproc];


  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofV  = msh->GetElementDofNumber(iel, solVType);  // number of solution element dofs

    unsigned eFlag = static_cast <unsigned>(floor((*sol->_Sol[eflagIndex])(iel) + 0.5));

    // resize local arrays
    l2GMap.resize(dim * 2 * nDofV);
    for(int k = 0; k < dim; k++) {
      solV1[k].resize(nDofV);
      solV2[k].resize(nDofV);
      aResV1[k].assign(nDofV, 0.);    //resize
      aResV2[k].assign(nDofV, 0.);    //resize
    }
    nodeFlag.resize(nDofV);

    for(int k = 0; k < dim; k++) {
      x[k].resize(nDofV);
    }


    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofV; i++) {
      unsigned iDof = msh->GetSolutionDof(i, iel, solVType);
      nodeFlag[i] = (*sol->_Sol[nflagIndex])(iDof);
      for(unsigned k = 0; k < dim; k++) {
        solV1[k][i] = (*sol->_Sol[solV1Index[k]])(iDof);
        solV2[k][i] = (*sol->_Sol[solV2Index[k]])(iDof);

        l2GMap[k * nDofV + i] = pdeSys->GetSystemDof(solV1Index[k], solV1PdeIndex[k], i, iel);
        l2GMap[(dim + k) * nDofV + i] = pdeSys->GetSystemDof(solV2Index[k], solV2PdeIndex[k], i, iel);

      }
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofV; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, 2);    // global to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k])(xDof);      // global extraction and local storage for the element coordinates
      }
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();
    if(eFlag == 0 || eFlag == 2) {
      // *** Element Gauss point loop ***
      for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solVType]->GetGaussPointNumber(); ig++) {
        // *** get gauss point weight, test function and test function partial derivatives ***
        msh->_finiteElement[ielGeom][solVType]->Jacobian(x, ig, weight, phi, phi_x);

        // evaluate the solution, the solution derivatives and the coordinates in the gauss point
        std::vector < std::vector < adept::adouble > > gradSolV1(dim);
        std::vector < std::vector < adept::adouble > > gradSolV2(dim);

        for(unsigned k = 0; k < dim; k++) {
          gradSolV1[k].assign(nDofV, 0.);
          gradSolV2[k].assign(nDofV, 0.);
        }

        for(unsigned i = 0; i < nDofV; i++) {
          for(unsigned k = 0; k < dim; k++) {
            for(unsigned j = 0; j < dim; j++) {
              gradSolV1[k][j] += phi_x[i * dim + j] * solV1[k][i];
              gradSolV2[k][j] += phi_x[i * dim + j] * solV2[k][i];
            }
          }
        }

        // *** phi_i loop ***
        for(unsigned i = 0; i < nDofV; i++) {

          adept::adouble divV1 = 0.;
          adept::adouble divV2 = 0.;
          for(unsigned k = 0; k < dim; k++) {
            divV1 += gradSolV1[k][k] * (eFlag == 0);
            divV2 += gradSolV2[k][k] * (eFlag == 2);
          }

          for(unsigned k = 0; k < dim; k++) {

            adept::adouble sigma1 = 0.;
            adept::adouble sigma2 = 0.;

            for(unsigned j = 0; j < dim; j++) {
              sigma1 += mu1 * (gradSolV1[k][j] + gradSolV1[j][k] * (eFlag == 0)) * phi_x[i * dim + j];
              sigma2 += mu2 * (gradSolV2[k][j] + gradSolV2[j][k] * (eFlag == 2)) * phi_x[i * dim + j];
            }
            //sigma1 += lambda1 * divD1 * phi_x[i * dim + k];
            //sigma2 += lambda2 * divD2 * phi_x[i * dim + k];

            if(eFlag == 0) {
              aResV1[k][i] += (- rho1 * g[k] * phi[i] + sigma1) * weight;
              if(nodeFlag[i] != 1) { // fake equation for sugular matrix
                aResV2[k][i] += sigma2 * weight;
              }
            }
            else if(eFlag == 2) {
              aResV2[k][i] += (- rho2 * g[k] * phi[i] + sigma2) * weight;
              if(nodeFlag[i] != 1) { // fake equation for sugular matrix
                aResV1[k][i] += sigma1 * weight;
              }
            }
          }

        } // end phi_i loop
      } // end gauss point loop
    }

    else {

      double iM1C1 = 1. / (mu1 * (*sol->_Sol[CMIndex[0]])(iel)); // 1/Ce1: element-wise largest eigenvalue in region1
      double iM2C2 = 1. / (mu2 * (*sol->_Sol[CMIndex[1]])(iel)); // 1/Ce2: element-wise largest eigenvalue in region2

      double denM = iM1C1 + iM2C2;

      double gammaM1 = iM1C1 / denM; // <u>_gamma  = gammaM1 * u1 + gammaM2 * u2
      double gammaM2 = iM2C2 / denM;

      double thetaM = 8. * 10 / denM;   // penalty parameter, sharp version thetaM = 2 / denM

      //std::cout << thetaM <<" ";

      //double iL1C1 = 1. / (lambda1 * (*sol->_Sol[CLIndex[0]])(iel));
      //double iL2C2 = 1. / (lambda2 * (*sol->_Sol[CLIndex[1]])(iel));

      //double denL = iL1C1 + iL2C2;

      // double gammaL1 = iL1C1 / denL;
      // double gammaL2 = iL2C2 / denL;

      // double thetaL = 4. / denL;


//       double gammaL1 = 0.5;
//       double gammaL2 = 0.5;
//
//       double gammaM1 = 0.5;
//       double gammaM2 = 0.5;
//
//
//       double thetaL = lambda1;
//       double thetaM = mu1;


      //bulk3
      while(imarker3 < markerOffset3[iproc + 1] && iel > particle3[imarker3]->GetMarkerElement()) {
        imarker3++;
      }
      while(imarker3 < markerOffset3[iproc + 1] && iel == particle3[imarker3]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particle3[imarker3]->GetMarkerLocalCoordinates();
        msh->_finiteElement[ielGeom][solVType]->Jacobian(x, xi, weight, phi, phi_x);
        double weight = particle3[imarker3]->GetMarkerMass();

        double dg1 = particle3[imarker3]->GetMarkerDistance();

        double dg2 = dg1 * dg1;
        double chi;
        if(dg1 < -deps)
          chi = 0.;
        else if(dg1 > deps) {
          chi = 1.;
        }
        else {
          chi = (a0 + dg1 * (a1 + dg2 * (a3 + dg2 * (a5 + dg2 * (a7 + dg2 * a9)))));
        }

        // evaluate the solution, the solution derivatives and the coordinates in the gauss point
        std::vector < std::vector < adept::adouble > > gradSolV1(dim);
        std::vector < std::vector < adept::adouble > > gradSolV2(dim);
        for(unsigned k = 0; k < dim; k++) {
          gradSolV1[k].assign(nDofV, 0.);
          gradSolV2[k].assign(nDofV, 0.);
        }

        for(unsigned i = 0; i < nDofV; i++) {
          for(unsigned k = 0; k < dim; k++) {
            for(unsigned j = 0; j < dim; j++) {
              gradSolV1[k][j] += phi_x[i * dim + j] * solV1[k][i];
              gradSolV2[k][j] += phi_x[i * dim + j] * solV2[k][i];
            }
          }
        }

        // *** phi_i loop ***
        for(unsigned i = 0; i < nDofV; i++) {

          adept::adouble divV1 = 0.;
          adept::adouble divV2 = 0.;
          for(unsigned k = 0; k < dim; k++) {
            divV1 += gradSolV1[k][k];
            divV2 += gradSolV2[k][k];
          }

          for(unsigned k = 0; k < dim; k++) {
            adept::adouble sigma1 = 0.;
            adept::adouble sigma2 = 0.;
            for(unsigned j = 0; j < dim; j++) {
              sigma1 += mu1 * (gradSolV1[k][j] + gradSolV1[j][k]) * phi_x[i * dim + j];
              sigma2 += mu2 * (gradSolV2[k][j] + gradSolV2[j][k]) * phi_x[i * dim + j];
            }
            //sigma1 += lambda1 * divD1 * phi_x[i * dim + k];

            aResV1[k][i] += (1. - chi) * (- rho1 * g[k] * phi[i] + sigma1) * weight;
            aResV2[k][i] += chi * (- rho2 * g[k] * phi[i] + sigma2) * weight;

          }
        } // end phi_i loop
        imarker3++;
      }

      // interface
      while(imarkerI < markerOffsetI[iproc + 1] && iel > particleI[imarkerI]->GetMarkerElement()) {
        imarkerI++;
      }
      while(imarkerI < markerOffsetI[iproc + 1] && iel == particleI[imarkerI]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particleI[imarkerI]->GetMarkerLocalCoordinates();
        msh->_finiteElement[ielGeom][solVType]->Jacobian(x, xi, weight, phi, phi_x);

        std::vector <std::vector < double > > T;
        particleI[imarkerI]->GetMarkerTangent(T);

        double weight;
        std::vector < double > N(dim);
        if(dim == 2) {
          N[0] =  T[0][1];
          N[1] = -T[0][0];
          weight = sqrt(N[0] * N[0] + N[1] * N[1]);
          N[0] /= weight;
          N[1] /= weight;
        }
        else {
          N[0] = T[0][1] * T[1][2] - T[0][2] * T[1][1];
          N[1] = T[0][2] * T[1][0] - T[0][0] * T[1][2];
          N[2] = T[0][0] * T[1][1] - T[0][1] * T[1][0];
          weight = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
          N[0] /= weight;
          N[1] /= weight;
          N[2] /= weight;
        }

        std::vector < adept::adouble > u1(dim, 0.);
        std::vector < adept::adouble > u2(dim, 0.);

        std::vector < adept::adouble > tau(dim, 0.);

        for(unsigned i = 0; i < nDofV; i++) {
          for(unsigned k = 0; k < dim; k++) {
            u1[k] += phi[i] * solV1[k][i];
            u2[k] += phi[i] * solV2[k][i];

            for(unsigned j = 0; j < dim; j++) {
              //tau[k] += (gammaL1 * lambda1 * solD1[j][i] * phi_x[i * dim + j] +
              //           gammaL2 * lambda2 * solD2[j][i] * phi_x[i * dim + j]) * N[k];

              tau[k] += (gammaM1 * mu1 * solV1[k][i] * phi_x[i * dim + j] +
                         gammaM2 * mu2 * solV2[k][i] * phi_x[i * dim + j]) * N[j];

              tau[k] += (gammaM1 * mu1 * solV1[j][i] * phi_x[i * dim + k] +
                         gammaM2 * mu2 * solV2[j][i] * phi_x[i * dim + k]) * N[j];

            }
          }
        }
        // *** phi_i loop ***
        for(unsigned i = 0; i < nDofV; i++) {
          for(unsigned k = 0; k < dim; k++) {
            aResV1[k][i] += tau[k] * phi[i] * weight;
            aResV1[k][i] += -thetaM * (u2[k] - u1[k]) * phi[i] * weight;

            aResV2[k][i] += -tau[k] * phi[i] * weight;
            aResV2[k][i] +=  thetaM * (u2[k] - u1[k]) * phi[i] * weight;

            for(unsigned j = 0; j < dim; j++) {
              //aResD1[k][i] += -gammaL1 * (lambda1 * phi_x[i * dim + k] * N[j] * (u2[j] - u1[j])) * weight;
              aResV1[k][i] += -gammaM1 * (mu1 * phi_x[i * dim + j] * N[j] * (u2[k] - u1[k])) * weight;
              aResV1[k][i] += -gammaM1 * (mu1 * phi_x[i * dim + j] * N[k] * (u2[j] - u1[j])) * weight;

              //aResD1[k][i] += -thetaL * (u2[j] - u1[j]) * N[j] * phi[i] * N[k] * weight;

              //aResD2[k][i] += -gammaL2 * (lambda2 * phi_x[i * dim + k] * N[j] * (u2[j] - u1[j])) * weight;
              aResV2[k][i] += -gammaM2 * (mu2 * phi_x[i * dim + j] * N[j] * (u2[k] - u1[k])) * weight;
              aResV2[k][i] += -gammaM2 * (mu2 * phi_x[i * dim + j] * N[k] * (u2[j] - u1[j])) * weight;

              //aResD2[k][i] +=  thetaL * (u2[j] - u1[j]) * N[j] * phi[i] * N[k] * weight;
            }


          }
        } // end phi_i loop
        imarkerI++;
      }

    }

    //copy the value of the adept::adoube aRes in double Res and store
    Res.resize(2 * dim * nDofV);    //resize

    for(unsigned k = 0; k < dim; k++) {
      for(int i = 0; i < nDofV; i++) {
        Res[k * nDofV + i] = - aResV1[k][i].value();
        Res[(k + dim) * nDofV + i] = - aResV2[k][i].value();
      }
    }

    RES->add_vector_blocked(Res, l2GMap);

    // define the dependent variables
    for(unsigned k = 0; k < dim; k++) {
      s.dependent(&aResV1[k][0], nDofV);
    }
    for(unsigned k = 0; k < dim; k++) {
      s.dependent(&aResV2[k][0], nDofV);
    }

    // define the independent variables
    for(unsigned k = 0; k < dim; k++) {
      s.independent(&solV1[k][0], nDofV);
    }
    for(unsigned k = 0; k < dim; k++) {
      s.independent(&solV2[k][0], nDofV);
    }

    // get the jacobian matrix (ordered by row major )
    Jac.resize(2 * dim * nDofV * 2 * dim * nDofV);   //resize
    s.jacobian(&Jac[0], true);

    //store K in the global matrix KK
    KK->add_matrix_blocked(Jac, l2GMap, l2GMap);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();
  KK->close();

//   PetscViewer viewer;
//   PetscViewerDrawOpen (PETSC_COMM_WORLD, NULL, NULL, 0, 0, 900, 900, &viewer);
//   PetscObjectSetName ( (PetscObject) viewer, "FSI matrix");
//   PetscViewerPushFormat (viewer, PETSC_VIEWER_DRAW_LG);
//   MatView ( (static_cast<PetscMatrix*> (KK))->mat(), viewer);
//   //VecView ( (static_cast<PetscVector*> (RES))->vec(), viewer);
//
//   double a;
//   std::cin >> a;

  //***************** END ASSEMBLY *******************
}

void BuildFlag(MultiLevelSolution& mlSol) {

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution *sol  = mlSol.GetSolutionLevel(level);
  Mesh     *msh   = mlSol._mlMesh->GetLevel(level);
  unsigned iproc  = msh->processor_id();

  unsigned eflagIndex = mlSol.GetIndex("eflag");
  unsigned nflagIndex = mlSol.GetIndex("nflag");

  unsigned nflagType = mlSol.GetSolutionType(nflagIndex);

  std::vector<Marker*> particle3 = line3->GetParticles();
  std::vector<unsigned> markerOffset3 = line3->GetMarkerOffset();
  unsigned imarker3 = markerOffset3[iproc];

  sol->_Sol[eflagIndex]->zero();
  sol->_Sol[nflagIndex]->zero();

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    bool inside = false;
    bool outside = false;
    bool interface = false;

    while(imarker3 < markerOffset3[iproc + 1] && iel > particle3[imarker3]->GetMarkerElement()) {
      imarker3++;
    }
    while(imarker3 < markerOffset3[iproc + 1] && iel == particle3[imarker3]->GetMarkerElement()) {
      double dg1 = particle3[imarker3]->GetMarkerDistance();
      if(dg1 < -1.0e-10) {
        outside = true;
        if(inside) {
          interface = true;
          break;
        }
      }
      else if(dg1 > 1.0e-10) {
        inside = true;
        if(outside) {
          interface = true;
          break;
        }
      }
      else {
        interface = true;
        break;
      }
      imarker3++;
    }
    if(interface) {
      sol->_Sol[eflagIndex]->set(iel, 1.);
      unsigned nDofu  = msh->GetElementDofNumber(iel, nflagType);  // number of solution element dofs
      for(unsigned i = 0; i < nDofu; i++) {
        unsigned iDof = msh->GetSolutionDof(i, iel, nflagType);
        sol->_Sol[nflagIndex]->set(iDof, 1.);
      }
    }
    else if(inside) {
      sol->_Sol[eflagIndex]->set(iel, 2.);
    }
  }

  sol->_Sol[eflagIndex]->close();
  sol->_Sol[nflagIndex]->close();

}


void GetInterfaceElementEigenvalues(MultiLevelSolution& mlSol) {

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution *sol  = mlSol.GetSolutionLevel(level);
  Mesh     *msh   = mlSol._mlMesh->GetLevel(level);
  unsigned iproc  = msh->processor_id();

  unsigned eflagIndex = mlSol.GetIndex("eflag");

  unsigned CMIndex[2];
  CMIndex[0] = mlSol.GetIndex("CM1");
  CMIndex[1] = mlSol.GetIndex("CM2");

  unsigned CLIndex[2];
  CLIndex[0] = mlSol.GetIndex("CL1");
  CLIndex[1] = mlSol.GetIndex("CL2");

  unsigned solIndex = mlSol.GetIndex("VX1");
  unsigned soluType = mlSol.GetSolutionType(solIndex);

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  std::vector < std::vector<double> >  x(dim);

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  double weight; // gauss point weight

  std::vector < std::vector < std::vector <double > > > aP(3);

  std::vector<Marker*> particle3 = line3->GetParticles();
  std::vector<unsigned> markerOffset3 = line3->GetMarkerOffset();
  unsigned imarker3 = markerOffset3[iproc];


  std::vector<Marker*> particleI = lineI->GetParticles();
  std::vector<unsigned> markerOffsetI = lineI->GetMarkerOffset();
  unsigned imarkerI = markerOffsetI[iproc];

  sol->_Sol[CMIndex[0]]->zero();
  sol->_Sol[CMIndex[1]]->zero();

  sol->_Sol[CLIndex[0]]->zero();
  sol->_Sol[CLIndex[1]]->zero();


  Eigen::MatrixXd AM;
  Eigen::MatrixXd AL;
  Eigen::MatrixXd BM[2];
  Eigen::MatrixXd BL[2];


  clock_t eigenTime = 0;

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    unsigned eFlag = static_cast <unsigned>(floor((*sol->_Sol[eflagIndex])(iel) + 0.5));
    if(eFlag == 1) {

      short unsigned ielGeom = msh->GetElementType(iel);
      unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);  // number of solution element dofs

      unsigned sizeAll = dim * nDofu;

      AM.resize(sizeAll, sizeAll);
      AM.setZero();

      AL.resize(sizeAll, sizeAll);
      AL.setZero();

      for(unsigned k = 0; k < 2; k++) {
        BM[k].resize(sizeAll, sizeAll);
        BM[k].setZero();

        BL[k].resize(sizeAll, sizeAll);
        BL[k].setZero();
      }

//       aM.assign(sizeAll * sizeAll, 0.);
//       bM[0].assign(sizeAll * sizeAll, 0.);
//       bM[1].assign(sizeAll * sizeAll, 0.);
//
//       aL.assign(sizeAll * sizeAll, 0.);
//       bL[0].assign(sizeAll * sizeAll, 0.);
//       bL[1].assign(sizeAll * sizeAll, 0.);

      for(int k = 0; k < dim; k++) {
        x[k].resize(nDofu);
      }

      for(unsigned i = 0; i < nDofu; i++) {
        unsigned xDof  = msh->GetSolutionDof(i, iel, 2);    // global to global mapping between coordinates node and coordinate dof
        for(unsigned k = 0; k < dim; k++) {
          x[k][i] = (*msh->_topology->_Sol[k])(xDof);      // global extraction and local storage for the element coordinates
        }
      }

      while(imarker3 < markerOffset3[iproc + 1] && iel > particle3[imarker3]->GetMarkerElement()) {
        imarker3++;
      }
      while(imarker3 < markerOffset3[iproc + 1] && iel == particle3[imarker3]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particle3[imarker3]->GetMarkerLocalCoordinates();
        msh->_finiteElement[ielGeom][soluType]->Jacobian(x, xi, weight, phi, phi_x);
        double weight = particle3[imarker3]->GetMarkerMass();

        double dg1 = particle3[imarker3]->GetMarkerDistance();

        double dg2 = dg1 * dg1;
        double chi;
        if(dg1 < -deps)
          chi = 0.;
        else if(dg1 > deps) {
          chi = 1.;
        }
        else {
          chi = (a0 + dg1 * (a1 + dg2 * (a3 + dg2 * (a5 + dg2 * (a7 + dg2 * a9)))));
        }

        for(unsigned k = 0; k < dim; k++) {
          for(unsigned i = 0; i < nDofu; i++) {
            for(unsigned l = 0; l < dim; l++) {
              for(unsigned j = 0; j < nDofu; j++) {
                BM[0](nDofu * k + i, k * nDofu + j) += (1. - chi) * 0.5 * phi_x[i * dim + l] * phi_x[j * dim + l] * weight;
                BM[0](nDofu * k + i, l * nDofu + j) += (1. - chi) * 0.5 * phi_x[i * dim + l] * phi_x[j * dim + k] * weight;

                BL[0](nDofu * k + i, l * nDofu + j) += (1. - chi) * phi_x[i * dim + k] * phi_x[j * dim + l] * weight;

                BM[1](nDofu * k + i, k * nDofu + j) += chi * 0.5 * phi_x[i * dim + l] * phi_x[j * dim + l] * weight;
                BM[1](nDofu * k + i, l * nDofu + j) += chi * 0.5 * phi_x[i * dim + l] * phi_x[j * dim + k] * weight;

                BL[1](nDofu * k + i, l * nDofu + j) += chi * phi_x[i * dim + k] * phi_x[j * dim + l] * weight;

              }
            }
          }
        }
        imarker3++;
      }



      // interface
      while(imarkerI < markerOffsetI[iproc + 1] && iel > particleI[imarkerI]->GetMarkerElement()) {
        imarkerI++;
      }
      while(imarkerI < markerOffsetI[iproc + 1] && iel == particleI[imarkerI]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particleI[imarkerI]->GetMarkerLocalCoordinates();
        msh->_finiteElement[ielGeom][soluType]->Jacobian(x, xi, weight, phi, phi_x);

        std::vector <std::vector < double > > T;
        particleI[imarkerI]->GetMarkerTangent(T);

        double weight;
        std::vector < double > N(dim);
        if(dim == 2) {
          N[0] =  T[0][1];
          N[1] = -T[0][0];
          weight = sqrt(N[0] * N[0] + N[1] * N[1]);
          N[0] /= weight;
          N[1] /= weight;
        }
        else {
          N[0] = T[0][1] * T[1][2] - T[0][2] * T[1][1];
          N[1] = T[0][2] * T[1][0] - T[0][0] * T[1][2];
          N[2] = T[0][0] * T[1][1] - T[0][1] * T[1][0];
          weight = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
          N[0] /= weight;
          N[1] /= weight;
          N[2] /= weight;
        }

        // *** phi_i loop ***

        for(unsigned k = 0; k < dim; k++) {
          for(int i = 0; i < nDofu; i++) {

            double gradPhiiDotN = 0.;
            for(unsigned l = 0; l < dim; l++) {
              gradPhiiDotN += phi_x[i * dim + l] * N[l];
            }
            for(int j = 0; j < nDofu; j++) {
              for(unsigned l = 0; l < dim; l++) {

                AM(nDofu * k + i, k * nDofu + j) += 0.5 * gradPhiiDotN * 0.5 * N[l] *  phi_x[j * dim + l]  * weight;
                AM(nDofu * k + i, l * nDofu + j) += 0.5 * gradPhiiDotN * 0.5 * N[l] *  phi_x[j * dim + k]  * weight;

                AL(nDofu * k + i, l * nDofu + j) += phi_x[i * dim + k] * phi_x[j * dim + l] * weight;

              }
              for(unsigned l1 = 0; l1 < dim; l1++) {
                for(unsigned l2 = 0; l2 < dim; l2++) {
                  AM(nDofu * k + i, l1 * nDofu + j) += 0.5 * N[k] * phi_x[i * dim + l1] * 0.5 * N[l2] *  phi_x[j * dim + l2]  * weight;
                  AM(nDofu * k + i, l2 * nDofu + j) += 0.5 * N[k] * phi_x[i * dim + l1] * 0.5 * N[l2] *  phi_x[j * dim + l1]  * weight;
                }
              }

            }
          } // end phi_i loop
        }
        imarkerI++;
      }

      double perturbation = 1e-10;

      std::cout << "======================EIGEN===================================" << std::endl;

      clock_t start = clock();

      Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges;
      double inf = 1e+10;

      for(unsigned k = 0; k < 2; k++) {
        double BM0Lk = BM[k].norm();

        for(unsigned i = 0; i < sizeAll; i++) {
          BM[k](i, i) += perturbation * BM0Lk;
        }

        ges.compute(AM, BM[k], false);
        std::complex < double > temp;
        Eigen::VectorXcd eig;

        eig = ges.eigenvalues();
        double emax0 = 0.;
        for(unsigned i = 0; i < sizeAll; i++) {
          temp = eig(i);
          if(fabs(real(temp)) > emax0 && fabs(real(temp)) < inf) {
            emax0 = fabs(real(temp));
          }
        }
        std::cout << iel << " " << emax0 << std::endl;
        sol->_Sol[CMIndex[k]]->set(iel, emax0);
      }

      for(unsigned k = 0; k < 2; k++) {
        double norm = BL[k].norm();

        for(unsigned i = 0; i < sizeAll; i++) {
          BL[k](i, i) += perturbation * norm;
        }

        ges.compute(AL, BL[k], false);
        std::complex < double > temp;
        Eigen::VectorXcd eig;

        eig = ges.eigenvalues();
        double emax0 = 0.;
        for(unsigned i = 0; i < sizeAll; i++) {
          temp = eig(i);
          if(fabs(real(temp)) > emax0 && fabs(real(temp)) < inf) {
            emax0 = fabs(real(temp));
          }
        }
        std::cout << iel << " " << emax0 << std::endl;
        sol->_Sol[CLIndex[k]]->set(iel, emax0);
      }
      eigenTime += (clock() - start);
    } // end of eflag loop
  } //end of element loop

  sol->_Sol[CMIndex[0]]->close();
  sol->_Sol[CMIndex[1]]->close();

  sol->_Sol[CLIndex[0]]->close();
  sol->_Sol[CLIndex[1]]->close();

  //std::cout << std::endl << "petsc TIME:\t" << static_cast<double>(petscTime) / CLOCKS_PER_SEC << std::endl;
  std::cout << std::endl << "Eigen TIME:\t" << static_cast<double>(eigenTime) / CLOCKS_PER_SEC << std::endl;

}


void GetInterfaceElementEigenvaluesAD(MultiLevelSolution& mlSol) {

  adept::Stack& s = FemusInit::_adeptStack;

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution *sol  = mlSol.GetSolutionLevel(level);
  Mesh     *msh   = mlSol._mlMesh->GetLevel(level);
  unsigned iproc  = msh->processor_id();

  unsigned eflagIndex = mlSol.GetIndex("eflag");

  unsigned CMIndex[2];
  CMIndex[0] = mlSol.GetIndex("CM1");
  CMIndex[1] = mlSol.GetIndex("CM2");

  unsigned CLIndex[2];
  CLIndex[0] = mlSol.GetIndex("CL1");
  CLIndex[1] = mlSol.GetIndex("CL2");

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  std::vector < unsigned > solIndex(dim);
  solIndex[0] = mlSol.GetIndex("VX1");
  solIndex[1] = mlSol.GetIndex("VX2");
  if(dim == 3) {
    solIndex[2] = mlSol.GetIndex("VX3");
  }
  unsigned soluType = mlSol.GetSolutionType(solIndex[0]);


  std::vector < std::vector<adept::adouble> >  solV(dim);
  std::vector < std::vector<double> >  x(dim);

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  double weight; // gauss point weight

  std::vector < std::vector < std::vector <double > > > aP(3);

  std::vector<Marker*> particle3 = line3->GetParticles();
  std::vector<unsigned> markerOffset3 = line3->GetMarkerOffset();
  unsigned imarker3 = markerOffset3[iproc];


  std::vector<Marker*> particleI = lineI->GetParticles();
  std::vector<unsigned> markerOffsetI = lineI->GetMarkerOffset();
  unsigned imarkerI = markerOffsetI[iproc];

  sol->_Sol[CMIndex[0]]->zero();
  sol->_Sol[CMIndex[1]]->zero();

  sol->_Sol[CLIndex[0]]->zero();
  sol->_Sol[CLIndex[1]]->zero();


  Eigen::MatrixXd AM;
  Eigen::MatrixXd AL;
  Eigen::MatrixXd BM[2];
  Eigen::MatrixXd BL[2];


  std::vector< adept::adouble > resAM;
  std::vector< adept::adouble > resAL;;
  std::vector< adept::adouble > resBM[2];
  std::vector< adept::adouble > resBL[2];

  std::vector< double > JacAM;
  std::vector< double > JacAL;;
  std::vector< double > JacBM[2];
  std::vector< double > JacBL[2];


  clock_t eigenTime = 0;

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    unsigned eFlag = static_cast <unsigned>(floor((*sol->_Sol[eflagIndex])(iel) + 0.5));
    if(eFlag == 1) {

      short unsigned ielGeom = msh->GetElementType(iel);
      unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);  // number of solution element dofs

      const unsigned sizeAll = dim * nDofu;

      AM.resize(sizeAll, sizeAll);
      AM.setZero();

      AL.resize(sizeAll, sizeAll);
      AL.setZero();

      for(unsigned k = 0; k < 2; k++) {
        BM[k].resize(sizeAll, sizeAll);
        BM[k].setZero();

        BL[k].resize(sizeAll, sizeAll);
        BL[k].setZero();
      }

      for(int k = 0; k < dim; k++) {
        x[k].resize(nDofu);
        solV[k].resize(nDofu);
      }

      resAM.assign(sizeAll, 0.);
      resAL.assign(sizeAll, 0.);
      resBM[0].assign(sizeAll, 0.);
      resBM[1].assign(sizeAll, 0.);
      resBL[0].assign(sizeAll, 0.);
      resBL[1].assign(sizeAll, 0.);

      for(unsigned i = 0; i < nDofu; i++) {

        unsigned uDof  = msh->GetSolutionDof(i, iel, soluType);
        unsigned xDof  = msh->GetSolutionDof(i, iel, 2);
        for(unsigned k = 0; k < dim; k++) {
          solV[k][i] = (*sol->_Sol[solIndex[k]])(uDof);
          x[k][i] = (*msh->_topology->_Sol[k])(xDof);
        }
      }

      s.new_recording();

      while(imarker3 < markerOffset3[iproc + 1] && iel > particle3[imarker3]->GetMarkerElement()) {
        imarker3++;
      }
      while(imarker3 < markerOffset3[iproc + 1] && iel == particle3[imarker3]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particle3[imarker3]->GetMarkerLocalCoordinates();
        msh->_finiteElement[ielGeom][soluType]->Jacobian(x, xi, weight, phi, phi_x);
        double weight = particle3[imarker3]->GetMarkerMass();

        std::vector<adept::adouble> solVg(dim, 0.);
        std::vector < std::vector<adept::adouble> > gradSolVg(dim);
        for(unsigned k = 0; k < dim; k++) {
          gradSolVg[k].assign(dim, 0.);
          for(unsigned i = 0; i < nDofu; i++) {
            solVg[k] += phi[i] * solV[k][i];
            for(unsigned l = 0; l < dim; l++) {
              gradSolVg[k][l] += phi_x[i * dim + l] * solV[k][i];
            }
          }
        }

        double dg1 = particle3[imarker3]->GetMarkerDistance();

        double dg2 = dg1 * dg1;
        double chi;
        if(dg1 < -deps)
          chi = 0.;
        else if(dg1 > deps) {
          chi = 1.;
        }
        else {
          chi = (a0 + dg1 * (a1 + dg2 * (a3 + dg2 * (a5 + dg2 * (a7 + dg2 * a9)))));
        }

        for(unsigned k = 0; k < dim; k++) {
          for(unsigned i = 0; i < nDofu; i++) {
            for(unsigned l = 0; l < dim; l++) {

              resBM[0][nDofu * k + i] += (1. - chi) * 0.5 * phi_x[i * dim + l] * gradSolVg[k][l] * weight;
              resBM[0][nDofu * k + i] += (1. - chi) * 0.5 * phi_x[i * dim + l] * gradSolVg[l][k] * weight;
              resBL[0][nDofu * k + i] += (1. - chi) * phi_x[i * dim + k] * gradSolVg[l][l] * weight;

              resBM[1][nDofu * k + i] += chi * 0.5 * phi_x[i * dim + l] * gradSolVg[k][l]  * weight;
              resBM[1][nDofu * k + i] += chi * 0.5 * phi_x[i * dim + l]  * gradSolVg[l][k] * weight;
              resBL[1][nDofu * k + i] += chi * phi_x[i * dim + k] * gradSolVg[l][l] * weight;

//               for(unsigned j = 0; j < nDofu; j++) {
//                 BM[0](nDofu * k + i, k * nDofu + j) += (1. - chi) * 0.5 * phi_x[i * dim + l] * phi_x[j * dim + l] * weight;
//                 BM[0](nDofu * k + i, l * nDofu + j) += (1. - chi) * 0.5 * phi_x[i * dim + l] * phi_x[j * dim + k] * weight;
// 
//                 BL[0](nDofu * k + i, l * nDofu + j) += (1. - chi) * phi_x[i * dim + k] * phi_x[j * dim + l] * weight;
// 
//                 BM[1](nDofu * k + i, k * nDofu + j) += chi * 0.5 * phi_x[i * dim + l] * phi_x[j * dim + l] * weight;
//                 BM[1](nDofu * k + i, l * nDofu + j) += chi * 0.5 * phi_x[i * dim + l] * phi_x[j * dim + k] * weight;
// 
//                 BL[1](nDofu * k + i, l * nDofu + j) += chi * phi_x[i * dim + k] * phi_x[j * dim + l] * weight;
// 
//               }
            }
          }
        }
        imarker3++;
      }



      // interface
      while(imarkerI < markerOffsetI[iproc + 1] && iel > particleI[imarkerI]->GetMarkerElement()) {
        imarkerI++;
      }
      while(imarkerI < markerOffsetI[iproc + 1] && iel == particleI[imarkerI]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particleI[imarkerI]->GetMarkerLocalCoordinates();
        msh->_finiteElement[ielGeom][soluType]->Jacobian(x, xi, weight, phi, phi_x);

        std::vector<adept::adouble> solVg(dim, 0.);
        std::vector < std::vector<adept::adouble> > gradSolVg(dim);
        for(unsigned k = 0; k < dim; k++) {
          gradSolVg[k].assign(dim, 0.);
          for(unsigned i = 0; i < nDofu; i++) {
            solVg[k] += phi[i] * solV[k][i];
            for(unsigned l = 0; l < dim; l++) {
              gradSolVg[k][l] += phi_x[i * dim + l] * solV[k][i];
            }
          }
        }

        std::vector <std::vector < double > > T;
        particleI[imarkerI]->GetMarkerTangent(T);

        double weight;
        std::vector < double > N(dim);
        if(dim == 2) {
          N[0] =  T[0][1];
          N[1] = -T[0][0];
          weight = sqrt(N[0] * N[0] + N[1] * N[1]);
          N[0] /= weight;
          N[1] /= weight;
        }
        else {
          N[0] = T[0][1] * T[1][2] - T[0][2] * T[1][1];
          N[1] = T[0][2] * T[1][0] - T[0][0] * T[1][2];
          N[2] = T[0][0] * T[1][1] - T[0][1] * T[1][0];
          weight = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
          N[0] /= weight;
          N[1] /= weight;
          N[2] /= weight;
        }

        // *** phi_i loop ***

        for(unsigned k = 0; k < dim; k++) {
          for(int i = 0; i < nDofu; i++) {

            double gradPhiiDotN = 0.;
            for(unsigned l = 0; l < dim; l++) {
              gradPhiiDotN += phi_x[i * dim + l] * N[l];
            }

            for(unsigned l = 0; l < dim; l++) {

              resAM[nDofu * k + i] += 0.5 * gradPhiiDotN * 0.5 * N[l] *  gradSolVg[k][l]  * weight;
              resAM[nDofu * k + i] += 0.5 * gradPhiiDotN * 0.5 * N[l] *  gradSolVg[l][k]  * weight;
              resAL[nDofu * k + i] += phi_x[i * dim + k] * gradSolVg[l][l] * weight;
            }
            for(unsigned l1 = 0; l1 < dim; l1++) {
              for(unsigned l2 = 0; l2 < dim; l2++) {
                resAM[nDofu * k + i] += 0.5 * N[k] * phi_x[i * dim + l1] * 0.5 * N[l2] *  gradSolVg[l1][l2] * weight;
                resAM[nDofu * k + i] += 0.5 * N[k] * phi_x[i * dim + l1] * 0.5 * N[l2] *  gradSolVg[l2][l1] * weight;
              }
            }

//             for(int j = 0; j < nDofu; j++) {
//               for(unsigned l = 0; l < dim; l++) {
// 
//                 AM(nDofu * k + i, k * nDofu + j) += 0.5 * gradPhiiDotN * 0.5 * N[l] *  phi_x[j * dim + l]  * weight;
//                 AM(nDofu * k + i, l * nDofu + j) += 0.5 * gradPhiiDotN * 0.5 * N[l] *  phi_x[j * dim + k]  * weight;
// 
//                 AL(nDofu * k + i, l * nDofu + j) += phi_x[i * dim + k] * phi_x[j * dim + l] * weight;
// 
//               }
//               for(unsigned l1 = 0; l1 < dim; l1++) {
//                 for(unsigned l2 = 0; l2 < dim; l2++) {
//                   AM(nDofu * k + i, l1 * nDofu + j) += 0.5 * N[k] * phi_x[i * dim + l1] * 0.5 * N[l2] *  phi_x[j * dim + l2]  * weight;
//                   AM(nDofu * k + i, l2 * nDofu + j) += 0.5 * N[k] * phi_x[i * dim + l1] * 0.5 * N[l2] *  phi_x[j * dim + l1]  * weight;
//                 }
//               }
//             }
          } // end phi_i loop
        }
        imarkerI++;
      }



      JacBM[0].resize(sizeAll * sizeAll);
      JacBM[1].resize(sizeAll * sizeAll);
      JacBL[0].resize(sizeAll * sizeAll);
      JacBL[1].resize(sizeAll * sizeAll);


      for(unsigned  k = 0; k < dim; k++) {
        s.independent(&solV[k][0], nDofu);
      }
      s.dependent(&resAM[0], sizeAll);
      JacAM.resize(sizeAll * sizeAll);
      s.jacobian(&JacAM[0] , true);
      s.clear_dependents();

      s.dependent(&resAL[0], sizeAll);
      JacAL.resize(sizeAll * sizeAll);
      s.jacobian(&JacAL[0] , true);
      s.clear_dependents();

      s.dependent(&resBM[0][0], sizeAll);
      JacBM[0].resize(sizeAll * sizeAll);
      s.jacobian(&JacBM[0][0] , true);
      s.clear_dependents();

      s.dependent(&resBM[1][0], sizeAll);
      JacBM[1].resize(sizeAll * sizeAll);
      s.jacobian(&JacBM[1][0] , true);
      s.clear_dependents();

      s.dependent(&resBL[0][0], sizeAll);
      JacBL[0].resize(sizeAll * sizeAll);
      s.jacobian(&JacBL[0][0] , true);
      s.clear_dependents();

      s.dependent(&resBL[1][0], sizeAll);
      JacBL[1].resize(sizeAll * sizeAll);
      s.jacobian(&JacBL[1][0] , true);
      s.clear_dependents();

      s.clear_independents();

      for(unsigned i = 0; i < sizeAll; i++) {
        for(unsigned j = 0; j < sizeAll; j++) {
          AM(i, j) = JacAM[i * sizeAll + j];
          AL(i, j) = JacAL[i * sizeAll + j];
          BM[0](i, j) = JacBM[0][i * sizeAll + j];
          BM[1](i, j) = JacBM[1][i * sizeAll + j];
          BL[0](i, j) = JacBL[0][i * sizeAll + j];
          BL[1](i, j) = JacBL[1][i * sizeAll + j];
        }
      }

      double perturbation = 1e-10;

      std::cout << "======================EIGEN===================================" << std::endl;

      clock_t start = clock();

      Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges;
      double inf = 1e+10;

      for(unsigned k = 0; k < 2; k++) {
        double BM0Lk = BM[k].norm();

        for(unsigned i = 0; i < sizeAll; i++) {
          BM[k](i, i) += perturbation * BM0Lk;
        }

        ges.compute(AM, BM[k], false);
        std::complex < double > temp;
        Eigen::VectorXcd eig;

        eig = ges.eigenvalues();
        double emax0 = 0.;
        for(unsigned i = 0; i < sizeAll; i++) {
          temp = eig(i);
          if(fabs(real(temp)) > emax0 && fabs(real(temp)) < inf) {
            emax0 = fabs(real(temp));
          }
        }
        std::cout << iel << " " << emax0 << std::endl;
        sol->_Sol[CMIndex[k]]->set(iel, emax0);
      }

      for(unsigned k = 0; k < 2; k++) {
        double norm = BL[k].norm();

        for(unsigned i = 0; i < sizeAll; i++) {
          BL[k](i, i) += perturbation * norm;
        }

        ges.compute(AL, BL[k], false);
        std::complex < double > temp;
        Eigen::VectorXcd eig;

        eig = ges.eigenvalues();
        double emax0 = 0.;
        for(unsigned i = 0; i < sizeAll; i++) {
          temp = eig(i);
          if(fabs(real(temp)) > emax0 && fabs(real(temp)) < inf) {
            emax0 = fabs(real(temp));
          }
        }
        std::cout << iel << " " << emax0 << std::endl;
        sol->_Sol[CLIndex[k]]->set(iel, emax0);
      }
      eigenTime += (clock() - start);
    } // end of eflag loop
  } //end of element loop

  sol->_Sol[CMIndex[0]]->close();
  sol->_Sol[CMIndex[1]]->close();

  sol->_Sol[CLIndex[0]]->close();
  sol->_Sol[CLIndex[1]]->close();

  //std::cout << std::endl << "petsc TIME:\t" << static_cast<double>(petscTime) / CLOCKS_PER_SEC << std::endl;
  std::cout << std::endl << "Eigen TIME:\t" << static_cast<double>(eigenTime) / CLOCKS_PER_SEC << std::endl;
}



void GetParticleWeights(MultiLevelSolution & mlSol) {

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution *sol  = mlSol.GetSolutionLevel(level);
  Mesh     *msh   = mlSol._mlMesh->GetLevel(level);
  unsigned iproc  = msh->processor_id();

  unsigned eflagIndex = mlSol.GetIndex("eflag");

  unsigned soluType = 0; //Linear Lagrange

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  std::vector < std::vector<double> >  x(dim);

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  double weight; // gauss point weight

  std::vector<Marker*> particle3 = line3->GetParticles();
  std::vector<unsigned> markerOffset3 = line3->GetMarkerOffset();

  unsigned imarker3 = markerOffset3[iproc];

  unsigned m = 3;  // Chebyshev degree

  unsigned i0, i1;
  if(dim == 3) {
    i0 = 0;
    i1 = 3;
  }
  else if(dim == 2) {
    i0 = 3;
    i1 = 5;
  }
  else {
    i0 = 5;
    i1 = 6;
  }

//grab the gauss points with elemtype and degree
  std::string name[6] = {"hex", "tet", "wedge", "quad", "tri", "line"};

  unsigned ng[6];
  Eigen::MatrixXd xg[6];
  Eigen::VectorXd wg[6];
  std::vector < double > jac[6];
  Eigen::MatrixXd Pg[6];

  for(unsigned i = i0; i < i1; i++) {

    const Gauss *gauss = new  Gauss(name[i].c_str(), "fourth");
    ng[i] = gauss->GetGaussPointsNumber();

    std::vector< const double * > Xg(dim);

    for(unsigned k = 0; k < dim; k++) {
      Xg[k] = gauss->GetGaussCoordinatePointer(k);
    }

    xg[i].resize(dim, ng[i]);
    for(unsigned k = 0; k < dim ; k++) {
      for(unsigned j = 0; j < ng[i]; j++) {
        xg[i](k, j) = Xg[k][j];
      }
    }

    const double *Wg = gauss->GetGaussWeightsPointer();
    wg[i].resize(ng[i]);
    for(unsigned j = 0; j < ng[i]; j++) {
      wg[i](j) = Wg[j];
    }

    jac[i].resize(ng[i]);

    Eigen::Tensor<double, 3, Eigen::RowMajor> PmG;
    GetChebXInfo(m, dim, ng[i], xg[i], PmG);

    GetMultiDimChebMatrix(dim, m, ng[i], PmG, Pg[i]);

    delete gauss;
  }

  Eigen::VectorXd F; // holds P_{n}(x_g) * wg * J(xg), n = 0,1,..m, g = 1,2,..ng  multidimensional
  F.resize(pow(m + 1, dim));
  Eigen::MatrixXd A; // // holds P_{n}(x_p) * wp , n = 0,1,..m, p = 1,2,..np multidimensional

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned eFlag = static_cast <unsigned>(floor((*sol->_Sol[eflagIndex])(iel) + 0.5));

    if(eFlag == 1) {
      short unsigned ielGeom = msh->GetElementType(iel);
      unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);  // number of solution element dofs

      for(int k = 0; k < dim; k++) {
        x[k].resize(nDofu);
      }

      for(unsigned i = 0; i < nDofu; i++) {
        unsigned xDof  = msh->GetSolutionDof(i, iel, 2);    // global to global mapping between coordinates node and coordinate dof
        for(unsigned k = 0; k < dim; k++) {
          x[k][i] = (*msh->_topology->_Sol[k])(xDof);     // global extraction and local storage for the element coordinates
        }
      }


      for(unsigned ig = 0; ig < ng[ielGeom] ; ig++) { // gauss loop to get Jacobians
        std::vector <double> xi(dim);
        for(unsigned k = 0; k < dim; k++) {
          xi[k] = xg[ielGeom](k, ig);
        }
        msh->_finiteElement[ielGeom][soluType]->Jacobian(x, xi, jac[ielGeom][ig], phi, phi_x);
      }

      //Assemble F
      F.setZero();
      for(unsigned i = 0; i < pow(m + 1, dim); i++) {
        for(unsigned j = 0; j < ng[ielGeom] ; j++) {
          F(i) += Pg[ielGeom](i, j) * jac[ielGeom][j] * wg[ielGeom](j);
        }
      }

      // identify the first particle inside iel
      while(imarker3 < markerOffset3[iproc + 1] && iel > particle3[imarker3]->GetMarkerElement()) {
        imarker3++;
      }
      unsigned imarker0 = imarker3;
      // loop on all particles inside iel to find how many particles are in iel
      while(imarker3 < markerOffset3[iproc + 1] && iel == particle3[imarker3]->GetMarkerElement()) {
        imarker3++;
      }
      unsigned nmarker = imarker3 - imarker0;
      imarker3 = imarker0;

      unsigned cnt = 0;
      Eigen::MatrixXd xP(dim, nmarker);
      Eigen::VectorXd wP(nmarker);


      // loop on all particles inside iel
      while(imarker3 < markerOffset3[iproc + 1] && iel == particle3[imarker3]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particle3[imarker3]->GetMarkerLocalCoordinates();
        //msh->_finiteElement[ielGeom][soluType]->Jacobian(x, xi, weight, phi, phi_x);
        double weight = particle3[imarker3]->GetMarkerMass();

        for(unsigned k = 0; k < dim; k++) {
          xP(k, cnt) = xi[k];
        }
        wP[cnt] = weight;

        cnt++;
        imarker3++;
      }

      Eigen::Tensor<double, 3, Eigen::RowMajor> PmX;
      GetChebXInfo(m, dim, nmarker, xP, PmX);

      GetMultiDimChebMatrix(dim, m, nmarker, PmX, A); //multidimensional Chebyshev polynomial evaluation in particle points up to m

      Eigen::VectorXd w_new;
      SolWeightEigen(A, F, wP, w_new); // New weights for iel are avaliable at this point

//       for(unsigned j = 0; j < nmarker; j++) {
//         std::cout << xP(0, j) << " " << xP(1, j) << " " << wP(j) << " " << w_new(j) << std::endl;
//       }

      // loop on all particles inside iel to attach the optimized weights to iel particles.
      imarker3 = imarker0;
      cnt = 0;
      while(imarker3 < markerOffset3[iproc + 1] && iel == particle3[imarker3]->GetMarkerElement()) {

        particle3[imarker3]->SetMarkerMass(w_new[cnt]);

        imarker3++;
        cnt++;
      }

    } // end of interface loop
  } // end of iel loop

}


void PrintMatlabMatrix(Eigen::MatrixXd &A) {

  std::cout << " = [";
  for(unsigned i = 0; i < A.rows(); i++) {
    for(unsigned j = 0; j < A.cols(); j++) {
      std::cout << A(i, j) << " ";
    }
    std::cout << ";" << std::endl;
  }
  std::cout << "];" << std::endl;
}





