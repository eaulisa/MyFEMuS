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

using namespace femus;

Line* line1;
Line* line2;
Line* line3;
Line* lineI;

unsigned DIM = 2;

void AssembleNitscheProblem_AD(MultiLevelProblem& mlProb);

void BuildFlag(MultiLevelSolution& mlSol);
void GetInterfaceElementEigenvalues(MultiLevelSolution& mlSol);

// void InitBallParticles(const unsigned & dim, const unsigned & ng, std::vector<double> &VxL, std::vector<double> &VxR,
//                        const std::vector < double> &xc, const double & R, const double & Rmax, const double & DR,
//                        std::vector < std::vector <double> > &xp, std::vector <double> &wp, std::vector <double> &dist);


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

  if(DIM == 2) {
    mlMsh.GenerateCoarseBoxMesh(nx, ny, 0, -lengthx / 2, lengthx / 2, -length / 2, length / 2, 0., 0., QUAD9, "seventh");
  }
  else if(DIM == 3) {
    nz = ny;
    mlMsh.GenerateCoarseBoxMesh(nx, ny, nz, -lengthx / 2, lengthx / 2, -length / 2, length / 2, -length / 2, length / 2, HEX27, "seventh");
  }

  double Lref = 1.;
  double Uref = 1.;
  double rho1 = 1000;
  double rho2 = 100;
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

  FEOrder femOrder = FIRST;

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

  //init marker

  std::vector<double> VxL = { - 0.5 * lengthx, -0.5 * length,-0.5 * length };
  std::vector<double> VxR = {  0.5 * lengthx,  0.5 * length,0.5 * length };

  
  double xc = 0.;
  double yc = 0.;
  double zc = 0.;
  double R = 0.125;
  double Rmax = 0.225;
  double DR2 = R / 10.;
  std::vector < double> Xc = {xc, yc, zc};
  
  std::vector < std::vector <double> > xp;
  std::vector <double> wp;
  std::vector <double> dist;
  
  InitBallParticles(DIM, VxL, VxR, Xc, R, Rmax, DR2, xp, wp, dist);
  
  std::vector < MarkerType > markerType(xp.size(), VOLUME);

  unsigned solType = 2;
  line3 = new Line(xp, wp, markerType, mlSol.GetLevel(numberOfUniformLevels - 1), solType);

  std::vector < std::vector < std::vector < double > > >  line3Points(1);
  line3->GetLine(line3Points[0]);
  PrintLine(DEFAULT_OUTPUTDIR, "bulk3", line3Points, 0);
  
  
  //BEGIN init particles
  unsigned size = 1;
  std::vector < std::vector < double > > x; // marker

 

  x.resize(size);
  x[0].resize(DIM, 0.);
  x[0][0] = xc;
  x[0][1] = yc;

  
  
  

//   Eigen::VectorXd wP = Eigen::VectorXd::Map(&wp[0], wp.size());
//   Eigen::MatrixXd xP(xp.size(), xp[0].size());
//   for(int i = 0; i < xp.size(); ++i) {
//     xP.row(i) = Eigen::VectorXd::Map(&xp[i][0], xp[0].size());
//   }
// 
//   //std::cout << " NumPoints =  "<< xP.cols() << std::endl;
// 
//   PrintMarkers(DIM, xP, dist, wP, wP, 0, 0);

//  return 1;




  unsigned NTHETA = 600; // number of partitions on the outer circle
  unsigned NR = NTHETA / (2 * M_PI); // number of partitions on the radial direction
  double DR = R / (NR + 0.5);

  std::vector < double > volume(size, M_PI * 0.25 * DR * DR);   // volume of the central particle

  double R2 = R - 0.5 * DR;

  for(unsigned i = 0; i < NR; i++) {
    double  r = R2 - i * DR;
    unsigned Ntheta = static_cast <unsigned>(ceil(NTHETA * r / R)); // number of partitions on the current circle
    double dtheta = 2 * M_PI / Ntheta;
    unsigned sizeOld = x.size();
    x.resize(sizeOld + Ntheta);
    volume.resize(sizeOld + Ntheta, M_PI * ((r + 0.5 * DR) * (r + 0.5 * DR) - (r - 0.5 * DR) * (r - 0.5 * DR)) / Ntheta);

    for(unsigned s = sizeOld; s < x.size(); s++) {
      x[s].resize(DIM);
    }
    for(unsigned j = 0; j < Ntheta; j++) {
      x[sizeOld + j][0] = xc + r * cos(j * dtheta);
      x[sizeOld + j][1] = yc + r * sin(j * dtheta);
    }
  }


  size = x.size();



  markerType.assign(size, VOLUME);
  line2 = new Line(x, volume, markerType, mlSol.GetLevel(numberOfUniformLevels - 1), solType);

  std::vector < std::vector < std::vector < double > > >  line2Points(1);
  line2->GetLine(line2Points[0]);
  PrintLine(DEFAULT_OUTPUTDIR, "bulk2", line2Points, 0);

  x.resize(0);
  volume.resize(0);

  NR = static_cast <unsigned>(ceil(0.1 / DR));

  double R1 = R + 0.5 * DR;

  for(unsigned i = 0; i < NR; i++) {
    double  r = R1 + i * DR;
    unsigned Ntheta = static_cast <unsigned>(ceil(NTHETA * r / R)); // number of partitions on the current circle
    double dtheta = 2 * M_PI / Ntheta;
    unsigned sizeOld = x.size();
    x.resize(sizeOld + Ntheta);
    volume.resize(sizeOld + Ntheta, M_PI * ((r + 0.5 * DR) * (r + 0.5 * DR) - (r - 0.5 * DR) * (r - 0.5 * DR)) / Ntheta);

    for(unsigned s = sizeOld; s < x.size(); s++) {
      x[s].resize(DIM);
    }
    for(unsigned j = 0; j < Ntheta; j++) {
      x[sizeOld + j][0] = xc + r * cos(j * dtheta);
      x[sizeOld + j][1] = yc + r * sin(j * dtheta);
    }
  }

  size = x.size();

  markerType.assign(size, VOLUME);

  line1 = new Line(x, volume, markerType, mlSol.GetLevel(numberOfUniformLevels - 1), solType);

  std::vector < std::vector < std::vector < double > > >  line1Points(1);
  line1->GetLine(line1Points[0]);
  PrintLine(DEFAULT_OUTPUTDIR, "bulk1", line1Points, 0);


  //interface marker initialization
  double Ntheta = NTHETA * 2;
  x.resize(Ntheta);
  markerType.assign(Ntheta, INTERFACE);

  for(unsigned i = 0; i < Ntheta; i++) {
    x[i].assign(DIM, 0.);
  }
  //std::vector < double > area(Ntheta, 2 * M_PI * R / Ntheta);  // uniform marker volume

  std::vector < std::vector < std::vector < double > > > T;
  T.resize(Ntheta);
  for(unsigned i = 0; i < Ntheta; i++) {
    T[i].resize(DIM - 1);
    for(unsigned k = 0; k < DIM - 1; k++) {
      T[i][k].resize(DIM, 0.);
    }
  }

  double arcLenght = 2. * M_PI * R / Ntheta;
  double dtheta = 2 * M_PI / Ntheta;

  //BEGIN initialization
  for(unsigned i = 0; i < Ntheta; i++) {

    x[i][0] = xc + R * cos(i * dtheta);
    x[i][1] = yc + R * sin(i * dtheta);

    T[i][0][0] = -arcLenght * sin(i * dtheta);
    T[i][0][1] = arcLenght * cos(i * dtheta);

  }

  lineI = new Line(x, T, markerType, mlSol.GetLevel(numberOfUniformLevels - 1), solType);

  std::vector < std::vector < std::vector < double > > > lineIPoints(1);
  lineI->GetLine(lineIPoints[0]);
  PrintLine(DEFAULT_OUTPUTDIR, "interfaceLine", lineIPoints, 0);
  //END interface markers


  BuildFlag(mlSol);

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

  std::vector<Marker*> particle1 = line1->GetParticles();
  std::vector<unsigned> markerOffset1 = line1->GetMarkerOffset();
  unsigned imarker1 = markerOffset1[iproc];


  std::vector<Marker*> particle2 = line2->GetParticles();
  std::vector<unsigned> markerOffset2 = line2->GetMarkerOffset();
  unsigned imarker2 = markerOffset2[iproc];

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

      double iM1C1 = 1. / (mu1 * (*sol->_Sol[CMIndex[0]])(iel));
      double iM2C2 = 1. / (mu2 * (*sol->_Sol[CMIndex[1]])(iel));

      double denM = iM1C1 + iM2C2;

      double gammaM1 = iM1C1 / denM;
      double gammaM2 = iM2C2 / denM;

      double thetaM = 8. / denM;

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



      //bulk1
      while(imarker1 < markerOffset1[iproc + 1] && iel > particle1[imarker1]->GetMarkerElement()) {
        imarker1++;
      }
      while(imarker1 < markerOffset1[iproc + 1] && iel == particle1[imarker1]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particle1[imarker1]->GetMarkerLocalCoordinates();
        msh->_finiteElement[ielGeom][solVType]->Jacobian(x, xi, weight, phi, phi_x);
        double weight = particle1[imarker1]->GetMarkerMass();


        // evaluate the solution, the solution derivatives and the coordinates in the gauss point
        std::vector < std::vector < adept::adouble > > gradSolV1(dim);
        for(unsigned k = 0; k < dim; k++) {
          gradSolV1[k].assign(nDofV, 0.);
        }

        for(unsigned i = 0; i < nDofV; i++) {
          for(unsigned k = 0; k < dim; k++) {
            for(unsigned j = 0; j < dim; j++) {
              gradSolV1[k][j] += phi_x[i * dim + j] * solV1[k][i];
            }
          }
        }

        // *** phi_i loop ***
        for(unsigned i = 0; i < nDofV; i++) {

          adept::adouble divV1 = 0.;
          for(unsigned k = 0; k < dim; k++) {
            divV1 += gradSolV1[k][k];
          }

          for(unsigned k = 0; k < dim; k++) {
            adept::adouble sigma1 = 0.;
            for(unsigned j = 0; j < dim; j++) {
              sigma1 += mu1 * (gradSolV1[k][j] + gradSolV1[j][k]) * phi_x[i * dim + j];
            }
            //sigma1 += lambda1 * divD1 * phi_x[i * dim + k];
            aResV1[k][i] += (- rho1 * g[k] * phi[i] + sigma1) * weight;
          }
        } // end phi_i loop
        imarker1++;
      }

      //bulk2
      while(imarker2 < markerOffset2[iproc + 1] && iel > particle2[imarker2]->GetMarkerElement()) {
        imarker2++;
      }
      while(imarker2 < markerOffset2[iproc + 1] && iel == particle2[imarker2]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particle2[imarker2]->GetMarkerLocalCoordinates();
        msh->_finiteElement[ielGeom][solVType]->Jacobian(x, xi, weight, phi, phi_x);
        double weight = particle2[imarker2]->GetMarkerMass();

        std::vector < std::vector < adept::adouble > > gradSolV2(dim);

        for(unsigned k = 0; k < dim; k++) {
          gradSolV2[k].assign(nDofV, 0.);
        }

        for(unsigned i = 0; i < nDofV; i++) {
          for(unsigned k = 0; k < dim; k++) {
            for(unsigned j = 0; j < dim; j++) {
              gradSolV2[k][j] += phi_x[i * dim + j] * solV2[k][i];
            }
          }
        }

        // *** phi_i loop ***
        for(unsigned i = 0; i < nDofV; i++) {
          adept::adouble divV2 = 0.;
          for(unsigned k = 0; k < dim; k++) {
            divV2 += gradSolV2[k][k];
          }
          for(unsigned k = 0; k < dim; k++) {
            adept::adouble sigma2 = 0.;
            for(unsigned j = 0; j < dim; j++) {
              sigma2 += mu2 * (gradSolV2[k][j] + gradSolV2[j][k]) * phi_x[i * dim + j];
            }
            //sigma2 += lambda2 * divD2 * phi_x[i * dim + k];
            aResV2[k][i] += (- rho2 * g[k] * phi[i] + sigma2) * weight;
          }
        }
        imarker2++;
      }

      // interface
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

  sol->_Sol[eflagIndex]->zero();
  sol->_Sol[nflagIndex]->zero();

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  std::vector < std::vector<double> >  x(dim);


  std::vector<Marker*> particle2 = line2->GetParticles();
  std::vector<unsigned> markerOffset2 = line2->GetMarkerOffset();

  for(unsigned imarker = markerOffset2[iproc]; imarker < markerOffset2[iproc + 1]; imarker++) {
    unsigned iel = particle2[imarker]->GetMarkerElement();
    sol->_Sol[eflagIndex]->set(iel, 2.);
  }


  std::vector<Marker*> particleI = lineI->GetParticles();
  std::vector<unsigned> markerOffsetI = lineI->GetMarkerOffset();

  for(unsigned imarker = markerOffsetI[iproc]; imarker < markerOffsetI[iproc + 1]; imarker++) {
    unsigned iel = particleI[imarker]->GetMarkerElement();
    sol->_Sol[eflagIndex]->set(iel, 1.);
    unsigned nDofu  = msh->GetElementDofNumber(iel, nflagType);  // number of solution element dofs
    for(unsigned i = 0; i < nDofu; i++) {
      unsigned iDof = msh->GetSolutionDof(i, iel, nflagType);
      sol->_Sol[nflagIndex]->set(iDof, 1.);
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

  std::vector<Marker*> particle1 = line1->GetParticles();
  std::vector<unsigned> markerOffset1 = line1->GetMarkerOffset();
  unsigned imarker1 = markerOffset1[iproc];

  std::vector<Marker*> particle2 = line2->GetParticles();
  std::vector<unsigned> markerOffset2 = line2->GetMarkerOffset();
  unsigned imarker2 = markerOffset2[iproc];

  std::vector<Marker*> particleI = lineI->GetParticles();
  std::vector<unsigned> markerOffsetI = lineI->GetMarkerOffset();
  unsigned imarkerI = markerOffsetI[iproc];

  sol->_Sol[CMIndex[0]]->zero();
  sol->_Sol[CMIndex[1]]->zero();

  sol->_Sol[CLIndex[0]]->zero();
  sol->_Sol[CLIndex[1]]->zero();

  std::vector < double > aM;
  std::vector < double > aM0;
  std::vector < std::vector < double > > bM(2);

  std::vector < double > aM1;
  std::vector < double > bM1;


  std::vector < double > aL;
  std::vector < double > aL0;
  std::vector < std::vector < double > > bL(2);

  std::vector < double > aL1;
  std::vector < double > bL1;


//   std::vector < double > a2;
//   std::vector < std::vector < double > > b2(2);

  Mat A, B;
  Vec v;
  EPS eps;

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned eFlag = static_cast <unsigned>(floor((*sol->_Sol[eflagIndex])(iel) + 0.5));
    if(eFlag == 1) {

      short unsigned ielGeom = msh->GetElementType(iel);
      unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);  // number of solution element dofs

      unsigned sizeAll = dim * nDofu;

      aM.assign(sizeAll * sizeAll, 0.);
      bM[0].assign(sizeAll * sizeAll, 0.);
      bM[1].assign(sizeAll * sizeAll, 0.);

      aL.assign(sizeAll * sizeAll, 0.);
      bL[0].assign(sizeAll * sizeAll, 0.);
      bL[1].assign(sizeAll * sizeAll, 0.);

      for(int k = 0; k < dim; k++) {
        x[k].resize(nDofu);  // Now we
      }

      for(unsigned i = 0; i < nDofu; i++) {
        unsigned xDof  = msh->GetSolutionDof(i, iel, 2);    // global to global mapping between coordinates node and coordinate dof
        for(unsigned k = 0; k < dim; k++) {
          x[k][i] = (*msh->_topology->_Sol[k])(xDof);      // global extraction and local storage for the element coordinates
        }
      }

      //bulk1
      while(imarker1 < markerOffset1[iproc + 1] && iel > particle1[imarker1]->GetMarkerElement()) {
        imarker1++;
      }
      while(imarker1 < markerOffset1[iproc + 1] && iel == particle1[imarker1]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particle1[imarker1]->GetMarkerLocalCoordinates();
        msh->_finiteElement[ielGeom][soluType]->Jacobian(x, xi, weight, phi, phi_x);
        double weight = particle1[imarker1]->GetMarkerMass();

        // *** phi_i loop ***

        for(unsigned k = 0; k < dim; k++) {
          for(unsigned i = 0; i < nDofu; i++) {
            for(unsigned l = 0; l < dim; l++) {
              for(unsigned j = 0; j < nDofu; j++) {
                bM[0][((nDofu * k) + i) * sizeAll + (k * nDofu + j)] += 0.5 * phi_x[i * dim + l] * phi_x[j * dim + l] * weight;
                bM[0][((nDofu * k) + i) * sizeAll + (l * nDofu + j)] += 0.5 * phi_x[i * dim + l] * phi_x[j * dim + k] * weight;

                bL[0][((nDofu * k) + i) * sizeAll + (l * nDofu + j)] += phi_x[i * dim + k] * phi_x[j * dim + l] * weight;

              }
            }
          }
        }
        imarker1++;
      }

      //bulk2
      while(imarker2 < markerOffset2[iproc + 1] && iel > particle2[imarker2]->GetMarkerElement()) {
        imarker2++;
      }
      while(imarker2 < markerOffset2[iproc + 1] && iel == particle2[imarker2]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particle2[imarker2]->GetMarkerLocalCoordinates();
        msh->_finiteElement[ielGeom][soluType]->Jacobian(x, xi, weight, phi, phi_x);
        double weight = particle2[imarker2]->GetMarkerMass();

        // *** phi_i loop ***

        for(unsigned k = 0; k < dim; k++) {
          for(unsigned i = 0; i < nDofu; i++) {
            for(unsigned l = 0; l < dim; l++) {
              for(unsigned j = 0; j < nDofu; j++) {
                bM[1][((nDofu * k) + i) * sizeAll + (k * nDofu + j)] += 0.5 * phi_x[i * dim + l] * phi_x[j * dim + l] * weight;
                bM[1][((nDofu * k) + i) * sizeAll + (l * nDofu + j)] += 0.5 * phi_x[i * dim + l] * phi_x[j * dim + k] * weight;

                bL[1][((nDofu * k) + i) * sizeAll + (l * nDofu + j)] += phi_x[i * dim + k] * phi_x[j * dim + l] * weight;
              }
            }
          }
        }
        imarker2++;
      }

      // interface
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

                aM[((nDofu * k) + i) * sizeAll + (k * nDofu + j) ] += 0.5 * gradPhiiDotN * 0.5 * N[l] *  phi_x[j * dim + l]  * weight;
                aM[((nDofu * k) + i) * sizeAll + (l * nDofu + j) ] += 0.5 * gradPhiiDotN * 0.5 * N[l] *  phi_x[j * dim + k]  * weight;

                aL[((nDofu * k) + i) * sizeAll + (l * nDofu + j)] += phi_x[i * dim + k] * phi_x[j * dim + l] * weight;

              }
              for(unsigned l1 = 0; l1 < dim; l1++) {
                for(unsigned l2 = 0; l2 < dim; l2++) {
                  aM[((nDofu * k) + i) * sizeAll + (l1 * nDofu + j)] += 0.5 * N[k] * phi_x[i * dim + l1] * 0.5 * N[l2] *  phi_x[j * dim + l2]  * weight;
                  aM[((nDofu * k) + i) * sizeAll + (l2 * nDofu + j)] += 0.5 * N[k] * phi_x[i * dim + l1] * 0.5 * N[l2] *  phi_x[j * dim + l1]  * weight;
                }
              }

            }
          } // end phi_i loop
        }
        imarkerI++;
      }

      unsigned sizeAll0 = sizeAll;
      aM0 = aM;
      aL0 = aL;

      for(unsigned s = 0; s < 2; s++) {

        sizeAll = sizeAll0;

        //BEGIN DEFLATION

        unsigned sizeAll1 = dim * (nDofu - 1);
        aM1.resize(sizeAll1 * sizeAll1);
        bM1.resize(sizeAll1 * sizeAll1);

        MatCreateSeqDense(PETSC_COMM_SELF, sizeAll1, sizeAll1, NULL, &B);

        for(int k = 0; k < dim; k++) {
          for(int i = 0; i < nDofu - 1; i++) {

            int ip = i + 1;
            int i1 = (nDofu - 1) * k + i;
            for(int l = 0; l < dim; l++) {
              for(int j = 0; j < nDofu - 1; j++) {
                int jp = j + 1;
                int j1 = (nDofu - 1) * l + j;
                double value;
                value = aM0[((nDofu * k) + ip) * sizeAll0 + (nDofu * l + jp)] - aM0[(nDofu * k) * sizeAll0 + (nDofu * l + jp)];
                aM1[i1 * sizeAll1 + j1] = value;

                value = bM[s][((nDofu * k) + ip) * sizeAll0 + (nDofu * l + jp)] - bM[s][(nDofu * k) * sizeAll0 + (nDofu * l + jp)];
                bM1[i1 * sizeAll1 + j1] = value;
                MatSetValues(B, 1, &i1, 1, &j1, &value, INSERT_VALUES);

              }
            }
          }
        }

        MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);

        sizeAll = sizeAll1;
        aM.swap(aM1);
        bM[s].swap(bM1);

        double real = 0.;
        while(fabs(real) < 1.0e-12) {

          MatCreateVecs(B, &v, NULL);

          EPSCreate(PETSC_COMM_SELF, &eps);
          EPSSetOperators(eps, B, NULL);
          EPSSetFromOptions(eps);
          EPSSetWhichEigenpairs(eps, EPS_SMALLEST_MAGNITUDE);
          EPSSolve(eps);

          double imaginary;
          EPSGetEigenpair(eps, 0, &real, &imaginary, v, NULL);
          EPSDestroy(&eps);

          if(fabs(real) < 1.0e-12) {
            PetscScalar *pv;
            VecGetArray(v, &pv);
            unsigned ii = 0;
            for(unsigned i = 1; i < sizeAll; i++) {
              if(fabs(pv[i]) > fabs(pv[ii])) ii = i;
            }

            unsigned sizeAll1 = sizeAll - 1;

            aM1.resize(sizeAll1 * sizeAll1);
            bM1.resize(sizeAll1 * sizeAll1);

            MatDestroy(&B);

            MatCreateSeqDense(PETSC_COMM_SELF, sizeAll1, sizeAll1, NULL, &B);

            for(unsigned i = 0; i < sizeAll; i++) {
              if(i != ii) {
                int i1 = (i < ii) ? i : i - 1;
                for(unsigned j = 0; j < sizeAll; j++) {
                  if(j != ii) {
                    int j1 = (j < ii) ? j : j - 1;
                    double value;
                    value = aM[i * sizeAll + j] - 1. / pv[ii] * pv[i] * aM[ii * sizeAll + j];
                    aM1[i1 * sizeAll1 + j1] = value;

                    value = bM[s][i * sizeAll + j] - 1. / pv[ii] * pv[i] * bM[s][ii * sizeAll + j];
                    bM1[i1 * sizeAll1 + j1] = value;
                    MatSetValues(B, 1, &i1, 1, &j1, &value, INSERT_VALUES);
                  }
                }
              }
            }
            MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
            VecRestoreArray(v, &pv);

            sizeAll = sizeAll1;
            aM.swap(aM1);
            bM[s].swap(bM1);
          }
          else {
            MatCreateSeqDense(PETSC_COMM_SELF, sizeAll, sizeAll, NULL, &A);
            for(int i = 0; i < sizeAll; i++) {
              for(int j = 0; j < sizeAll; j++) {
                MatSetValues(A, 1, &i, 1, &j, &aM[i * sizeAll + j], INSERT_VALUES);
              }
            }
            MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
          }

          VecDestroy(&v);
        }

        //END DEFLATION

        EPSCreate(PETSC_COMM_SELF, &eps);
        EPSSetOperators(eps, A, B);
        EPSSetFromOptions(eps);
        EPSSetWhichEigenpairs(eps, EPS_LARGEST_MAGNITUDE);
        EPSSolve(eps);

        EPSGetEigenpair(eps, 0, &real, NULL, NULL, NULL);
        std::cout << real << " " << std::endl;

        sol->_Sol[CMIndex[s]]->set(iel, real);

        EPSDestroy(&eps);
        MatDestroy(&A);
        MatDestroy(&B);

      }

      for(unsigned s = 0; s < 2; s++) {

        sizeAll = sizeAll0;

        //BEGIN DEFLATION

        unsigned sizeAll1 = dim * (nDofu - 1);
        aL1.resize(sizeAll1 * sizeAll1);
        bL1.resize(sizeAll1 * sizeAll1);

        MatCreateSeqDense(PETSC_COMM_SELF, sizeAll1, sizeAll1, NULL, &B);

        for(int k = 0; k < dim; k++) {
          for(int i = 0; i < nDofu - 1; i++) {

            int ip = i + 1;
            int i1 = (nDofu - 1) * k + i;
            for(int l = 0; l < dim; l++) {
              for(int j = 0; j < nDofu - 1; j++) {
                int jp = j + 1;
                int j1 = (nDofu - 1) * l + j;
                double value;
                value = aL0[((nDofu * k) + ip) * sizeAll0 + (nDofu * l + jp)] - aL0[(nDofu * k) * sizeAll0 + (nDofu * l + jp)];
                aL1[i1 * sizeAll1 + j1] = value;

                value = bL[s][((nDofu * k) + ip) * sizeAll0 + (nDofu * l + jp)] - bL[s][(nDofu * k) * sizeAll0 + (nDofu * l + jp)];
                bL1[i1 * sizeAll1 + j1] = value;
                MatSetValues(B, 1, &i1, 1, &j1, &value, INSERT_VALUES);

              }
            }
          }
        }

        MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);

        sizeAll = sizeAll1;
        aL.swap(aL1);
        bL[s].swap(bL1);

        double real = 0.;
        while(fabs(real) < 1.0e-10) {

          MatCreateVecs(B, &v, NULL);

          EPSCreate(PETSC_COMM_SELF, &eps);
          EPSSetOperators(eps, B, NULL);
          EPSSetFromOptions(eps);
          EPSSetWhichEigenpairs(eps, EPS_SMALLEST_MAGNITUDE);
          EPSSolve(eps);

          double imaginary;
          EPSGetEigenpair(eps, 0, &real, &imaginary, v, NULL);
          EPSDestroy(&eps);

          if(fabs(real) < 1.0e-10) {
            PetscScalar *pv;
            VecGetArray(v, &pv);
            unsigned ii = 0;
            for(unsigned i = 1; i < sizeAll; i++) {
              if(fabs(pv[i]) > fabs(pv[ii])) ii = i;
            }

            unsigned sizeAll1 = sizeAll - 1;

            aL1.resize(sizeAll1 * sizeAll1);
            bL1.resize(sizeAll1 * sizeAll1);

            MatDestroy(&B);

            MatCreateSeqDense(PETSC_COMM_SELF, sizeAll1, sizeAll1, NULL, &B);

            for(unsigned i = 0; i < sizeAll; i++) {
              if(i != ii) {
                int i1 = (i < ii) ? i : i - 1;
                for(unsigned j = 0; j < sizeAll; j++) {
                  if(j != ii) {
                    int j1 = (j < ii) ? j : j - 1;
                    double value;
                    value = aL[i * sizeAll + j] - 1. / pv[ii] * pv[i] * aL[ii * sizeAll + j];
                    aL1[i1 * sizeAll1 + j1] = value;

                    value = bL[s][i * sizeAll + j] - 1. / pv[ii] * pv[i] * bL[s][ii * sizeAll + j];
                    bL1[i1 * sizeAll1 + j1] = value;
                    MatSetValues(B, 1, &i1, 1, &j1, &value, INSERT_VALUES);
                  }
                }
              }
            }
            MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
            VecRestoreArray(v, &pv);

            sizeAll = sizeAll1;
            aL.swap(aL1);
            bL[s].swap(bL1);
          }
          else {
            MatCreateSeqDense(PETSC_COMM_SELF, sizeAll, sizeAll, NULL, &A);
            for(int i = 0; i < sizeAll; i++) {
              for(int j = 0; j < sizeAll; j++) {
                MatSetValues(A, 1, &i, 1, &j, &aL[i * sizeAll + j], INSERT_VALUES);
              }
            }
            MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
          }

          VecDestroy(&v);
        }

        //END DEFLATION

        EPSCreate(PETSC_COMM_SELF, &eps);
        EPSSetOperators(eps, A, B);
        EPSSetFromOptions(eps);
        EPSSetWhichEigenpairs(eps, EPS_LARGEST_MAGNITUDE);
        EPSSolve(eps);

        EPSGetEigenpair(eps, 0, &real, NULL, NULL, NULL);
        std::cout << real << " " << std::endl;

        sol->_Sol[CLIndex[s]]->set(iel, real);

        EPSDestroy(&eps);
        MatDestroy(&A);
        MatDestroy(&B);

      }
    }
  }

  sol->_Sol[CMIndex[0]]->close();
  sol->_Sol[CMIndex[1]]->close();

  sol->_Sol[CLIndex[0]]->close();
  sol->_Sol[CLIndex[1]]->close();
}





