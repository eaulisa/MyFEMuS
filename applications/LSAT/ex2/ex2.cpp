/** \file Ex13.cpp
 *  \brief This example shows how to set and solve the weak form
 *   of the time dependent Stress-Strain Equation, for a beam with pressure
 *
 *  \nabla \cdot \sigma + \mass*acceleration = F
 *  in a Beam domain (in 2D and 3D) clamped on one end
 *
 *  \author Eugenio Aulisa
 */


#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "LinearImplicitSystem.hpp"
#include "TransientSystem.hpp"
#include "adept.h"

const unsigned DIM = 2;
const double BETA = 0.25;
const double GAMMA = 0.5;
bool UseNewmarkUpdateWithD = true;
double dt = 0.025;
bool withDisturbance = false;

static double PStar = 0.0;
static std::vector<double> PStarNodes;
std::vector<unsigned> g_controlNodeDofs;


struct WNodeIDs {
  unsigned W0;  // node for W-equation #1
  unsigned X1;  // node for W-equation #2
  unsigned Y1;  // node for W-equation #3
};
WNodeIDs g_specialWnodes;


double SetVariableTimeStep(const double time) {
  return dt;
}


using namespace std;
using namespace femus;

struct RegionBox {
  double xMin, xMax;
  double yMin, yMax;
};


void SetRegions(Solution *sol, /*const RegionBox& boxB,*/ const RegionBox& boxC, const RegionBox* boxBd = nullptr);

void SetPrescribedFields(Solution* sol, const double& time, const std::string& R, const std::string& D = "");

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true;
  if(!strcmp(SolName, "Wi")) dirichlet = false;
  else{
    dirichlet = true;
    value = 0.;
  }



  // if((DIM == 2 && facename == 4) || (DIM == 3 && facename == 5)) {   // left boundary condition.
  //   dirichlet = true;
  // }
  //
  // if((DIM == 2 && facename == 3) || (DIM == 3 && facename == 4)) {   // top boundary condition.
  //   dirichlet = false;
  //   value = 500.;
  // }

  return dirichlet;
}


void NewmarkUpdate(MultiLevelSolution *mlSol);
void NewmarkUpdateWithD(MultiLevelSolution *mlSol);

void AssembleResAD(MultiLevelProblem& ml_prob);
void AssembleResP(MultiLevelProblem& ml_prob);

std::pair<double, std::vector<double>> PrecomputePstarIntegrals(Solution* sol);


std::vector<double> ComputeL2NormCascadeOverC(Solution* sol, unsigned cascadeIterations);


unsigned FindClosestNode(const Mesh* msh, const std::vector<double>& x0);
std::vector<double> MarkControlNodes(const Mesh* msh, const std::vector<std::vector<double>>& points);
std::vector<unsigned> GetControlNodeIndices(const Mesh* msh, const std::vector<std::vector<double>>& points);

// std::vector<unsigned> GetGlobalNodeIDsForW(const Mesh* msh, unsigned elemID, const std::vector<unsigned>& localNodeIDs);
WNodeIDs GetGlobalNodeIDsForW(const Mesh* msh, unsigned elemID);

int main(int argc, char** args) {
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  MultiLevelMesh mlMsh;

  unsigned nx = 10;
  unsigned ny = 2;
  unsigned nz = 1;

  double length = 1;
  double lengthx = M_PI;

  if (DIM == 2) {
    mlMsh.GenerateCoarseBoxMesh(nx, ny, 0, 0., lengthx, 0., length, 0., 0., QUAD9, "seventh");
  } else if (DIM == 3) {
    nz = ny;
    mlMsh.GenerateCoarseBoxMesh(nx, ny, nz, 0., lengthx, 0., length, 0., length, HEX27, "seventh");
  }

  unsigned numberOfUniformLevels = 4;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels, numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);
  mlMsh.PrintInfo();

  MultiLevelSolution mlSol(&mlMsh);

  mlSol.AddSolution("Zi", LAGRANGE, SECOND, 2);
  mlSol.AddSolution("Xi", LAGRANGE, SECOND);
  mlSol.AddSolution("Yi", LAGRANGE, SECOND);
  mlSol.AddSolution("Ei", LAGRANGE, SECOND, false);
  mlSol.AddSolution("Z", LAGRANGE, SECOND, false);

  mlSol.AddSolution("Wi", LAGRANGE, SECOND, 2);

  mlSol.AddSolution("P", LAGRANGE, SECOND, false);

  unsigned cascadeIterations = 3;
  for (unsigned j = 0; j < cascadeIterations; j++) {
    std::string Zj = "Z" + std::to_string(j);
    std::string ZjOld = "Z" + std::to_string(j) + "Old";
    std::string Ej = "E" + std::to_string(j);
    mlSol.AddSolution(Zj.c_str(), LAGRANGE, SECOND, false);
    mlSol.AddSolution(ZjOld.c_str(), LAGRANGE, SECOND, false);
    mlSol.AddSolution(Ej.c_str(), LAGRANGE, SECOND, false);
  }

  mlSol.AddSolution("R", LAGRANGE, SECOND, false);
  // mlSol.AddSolution("B", DISCONTINUOUS_POLYNOMIAL, ZERO, false);
  mlSol.AddSolution("C", DISCONTINUOUS_POLYNOMIAL, ZERO, false);
  if(withDisturbance) {
    mlSol.AddSolution("d", LAGRANGE, SECOND, false);
    mlSol.AddSolution("Bd", DISCONTINUOUS_POLYNOMIAL, ZERO, false);
  }
  mlSol.Initialize("All");
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("All");

  MultiLevelProblem mlProb(&mlSol);

  // New stationary system for P
  LinearImplicitSystem& systemP = mlProb.add_system<LinearImplicitSystem>("LP");
  systemP.AddSolutionToSystemPDE("P");
  systemP.SetAssembleFunction(AssembleResP);
  systemP.init();



  TransientLinearImplicitSystem& system = mlProb.add_system<TransientLinearImplicitSystem>("LSAT");

  system.AddSolutionToSystemPDE("Zi");
  system.AddSolutionToSystemPDE("Xi");
  system.AddSolutionToSystemPDE("Yi");
  system.AddSolutionToSystemPDE("Wi");
  system.SetAssembleFunction(AssembleResAD);
  system.AttachGetTimeIntervalFunction(SetVariableTimeStep);
  system.init();
  system.SetOuterSolver(PREONLY);

  const unsigned level = mlMsh.GetNumberOfLevels() - 1;
  Solution* sol = mlSol.GetSolutionLevel(level);
  Mesh* msh = sol->GetMesh();

  std::vector<std::vector<double>> controlPoints = {
  {0.3, 0.3},
  {1.4, 0.7},
  {2.9, 0.5}
  };

  g_controlNodeDofs = GetControlNodeIndices(msh, controlPoints);


  // RegionBox boxB{M_PI/3., 2*M_PI/3., 0., 1};
  RegionBox boxC{M_PI/3., 2*M_PI/3., 0., 1};
  if(withDisturbance) {
    RegionBox boxBd{2*M_PI/3., M_PI, 0., 1};
    // SetRegions(sol, boxB, boxC, &boxBd);
    SetRegions(sol, boxC, &boxBd);
  }
  else SetRegions(sol, boxC);
  // else SetRegions(sol, boxB, boxC);

  sol->_Sol[mlSol.GetIndex("P")]->zero();

  std::vector<std::string> variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(false);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);


  systemP.MGsolve();

  auto [IntegralP, WeightsP] = PrecomputePstarIntegrals(sol);

  PStar = IntegralP;
  PStarNodes = std::move(WeightsP);

  // Find and number 3 nodes for W solution
  const unsigned targetElem = 0;
  g_specialWnodes = GetGlobalNodeIDsForW(msh, targetElem);

  // ****here we solve the system for P****

  const unsigned n_timesteps = 150;

  for (unsigned j = 0; j < cascadeIterations; j++) {
    std::string ZjOld = "Z" + std::to_string(j) + "Old";
    sol->_Sol[(mlSol.GetIndex(ZjOld.c_str()))]->zero();
  }

  int world_rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  std::ofstream errorFile;
  if (world_rank == 0) {
    errorFile.open("error.dat", std::ios::out | std::ios::trunc);
    if (!errorFile.is_open()) {
      std::cerr << "ERROR: could not open error.dat for writing\n";
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    errorFile.setf(std::ios::scientific);
    errorFile << "# time";
    for (unsigned j = 0; j < cascadeIterations; ++j) errorFile << "   E" << j;
    errorFile << "\n";
    errorFile << std::setprecision(8);
  }

  // Time loop
  for (unsigned t = 1; t <= n_timesteps; t++) {
    if(withDisturbance) SetPrescribedFields(sol, t * dt, "R", "d");
    else SetPrescribedFields(sol, t * dt, "R");

    std::shared_ptr<NumericVector> dOriginal;
    if (withDisturbance) dOriginal = sol->_Sol[mlSol.GetIndex("d")]->clone();

    sol->_Sol[mlSol.GetIndex("Z")]->zero();
    *(sol->_Sol[mlSol.GetIndex("Ei")]) = *(sol->_Sol[(mlSol.GetIndex("R"))]);

    for (unsigned j = 0; j < cascadeIterations; j++) {
      std::string Zj = "Z" + std::to_string(j);
      std::string ZjOld = "Z" + std::to_string(j) + "Old";
      std::string Ej = "E" + std::to_string(j);

      //For withDisturbance, and j = 1,2,.., we set the disturbance to be zero
      if (withDisturbance && j > 0) sol->_Sol[mlSol.GetIndex("d")]->zero();

      *(sol->_Sol[mlSol.GetIndex("Zi")]) = *(sol->_Sol[(mlSol.GetIndex(ZjOld.c_str()))]);
      *(sol->_SolOld[mlSol.GetIndex("Zi")]) = *(sol->_Sol[(mlSol.GetIndex(ZjOld.c_str()))]);
      system.MGsolve();

      *(sol->_Sol[mlSol.GetIndex("Z")]) += *(sol->_Sol[(mlSol.GetIndex("Zi"))]);
      *(sol->_Sol[mlSol.GetIndex(Zj.c_str())]) = *(sol->_Sol[(mlSol.GetIndex("Zi"))]);
      *(sol->_Sol[mlSol.GetIndex(Ej.c_str())]) = *(sol->_Sol[(mlSol.GetIndex("Ei"))]);
      *(sol->_Sol[mlSol.GetIndex(Ej.c_str())]) -= *(sol->_Sol[(mlSol.GetIndex("Zi"))]);
      *(sol->_Sol[mlSol.GetIndex("Ei")]) = *(sol->_Sol[(mlSol.GetIndex(Ej.c_str()))]);
    }

    // restore prescribed disturbance before visualization output
    if (withDisturbance) *(sol->_Sol[mlSol.GetIndex("d")]) = *dOriginal;

    vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, t);

    std::vector<double> L2_E_Cascade = ComputeL2NormCascadeOverC(sol, cascadeIterations);
    if (world_rank == 0) {
      errorFile << std::setw(14) << t * dt;
      for (double val : L2_E_Cascade) errorFile << " " << std::setw(14) << val;
      errorFile << "\n";
      errorFile.flush();
    }

    for (unsigned j = 0; j < cascadeIterations; j++) {
      std::string Zj = "Z" + std::to_string(j);
      std::string ZjOld = "Z" + std::to_string(j) + "Old";
      *(sol->_Sol[(mlSol.GetIndex(ZjOld.c_str()))]) = *(sol->_Sol[mlSol.GetIndex(Zj.c_str())]);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (world_rank == 0) errorFile.close();

  return 0;
}



double flc4hs(double const & x, double const & eps) {

  double r = x / eps;
  if (r < -1) {
    return 0.;
  }
  else if (r < 1.) {
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

/*double GetTargetSolution(const std::vector<double> &xv, const double &time) {
  // return cos(xv[0] * xv[1]) * flc4hs(time - .5, 0.5) ;
  return xv[0] * (M_PI - xv[0]) * xv[1] * (1. - xv[1]) * sin(time) * flc4hs(time - .5, 0.5) ;
}*/

double GetTargetSolution(const std::vector<double> &xv, const double &time) {

  // return cos(xv[0] * xv[1]) * flc4hs(time - .5, 0.5) ;
  return xv[0] * (M_PI - xv[0]) * xv[1] * (1. - xv[1]) * sin(time) * flc4hs(time - .5, 0.5) ;
}


/*double GetDisturbanceSolution(const std::vector<double>& xv, const double& time) {
  return xv[1] * (1. - xv[1]) * sin(2. * time) * flc4hs(time - .5, 0.5);
}*/

double GetDisturbanceSolution(const std::vector<double>& xv, const double& time) {
    return (xv[0] - 2*M_PI/3.0)*(M_PI - xv[0]) * xv[1]*(1 - xv[1]) * sin(2*time) * flc4hs(time - 0.5, 0.5);
}

void SetPrescribedFields(Solution* sol, const double& time, const std::string& R, const std::string& D) {

  Mesh* msh = sol->GetMesh();
  const unsigned dim = msh->GetDimension();
  unsigned iproc = msh->processor_id();

  unsigned rIndex = sol->GetIndex(R.c_str());
  unsigned rType = sol->GetSolutionType(rIndex);

  bool hasDisturbance = false;
  unsigned dIndex = 0;
  unsigned dType = 0;

  // Check if disturbance should be applied and if the field exists
  if (withDisturbance && !D.empty()) {
    if (sol->GetIndex(D.c_str()) != static_cast<unsigned>(-1)) {
      hasDisturbance = true;
      dIndex = sol->GetIndex(D.c_str());
      dType = sol->GetSolutionType(dIndex);
    }
  }

  std::vector<double> xv(dim);
  unsigned xType = 2;  // coordinates always quadratic

  for (unsigned iel = msh->_elementOffset[iproc];
       iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned nDofsR = msh->GetElementDofNumber(iel, rType);

    for (unsigned i = 0; i < nDofsR; i++) {
      unsigned xDof = msh->GetSolutionDof(i, iel, xType);
      for (unsigned k = 0; k < dim; k++)
        xv[k] = (*msh->_topology->_Sol[k])(xDof);

      unsigned uDofR = msh->GetSolutionDof(i, iel, rType);
      double rVal = GetTargetSolution(xv, time);
      sol->_Sol[rIndex]->set(uDofR, rVal);

      if (hasDisturbance) {
        unsigned uDofD = msh->GetSolutionDof(i, iel, dType);
        double dVal = GetDisturbanceSolution(xv, time);
        sol->_Sol[dIndex]->set(uDofD, dVal);
      }
    }
  }

  sol->_Sol[rIndex]->close();
  if (hasDisturbance) sol->_Sol[dIndex]->close();
}

bool CheckIfInside(const std::vector<double>& xv, const RegionBox& box) {
  return (xv[0] > box.xMin && xv[0] < box.xMax &&
          xv[1] > box.yMin && xv[1] < box.yMax);
}

// double CheckIfInsideB(const std::vector<double> &xv) {
//   return (xv[0] > 0.5 && xv[0] < 1. && xv[1] > 0.25 && xv[1] < 0.75);
// }
//
// double CheckIfInsideC(const std::vector<double> &xv) {
//   return (xv[0] > 0.5 && xv[0] < 1. && xv[1] > 0.25 && xv[1] < 0.75);
// }

void SetRegions(Solution* sol,
                  // const RegionBox& boxB,
                  const RegionBox& boxC,
                  const RegionBox* boxBd) {

  Mesh* msh = sol->GetMesh();
  const unsigned dim = msh->GetDimension();
  unsigned iproc = msh->processor_id();

  // unsigned IndexB = sol->GetIndex("B");
  unsigned IndexC = sol->GetIndex("C");
  unsigned IndexBd = 0;
  bool hasBd = (boxBd != nullptr);
  if (hasBd) IndexBd = sol->GetIndex("Bd");

  std::vector<double> xv(dim);
  unsigned xType = 2;

  for (unsigned iel = msh->_elementOffset[iproc];
       iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned nDofs = msh->GetElementDofNumber(iel, 1);
    // bool elementInB = true;
    bool elementInC = true;
    bool elementInBd = true;

    for (unsigned i = 0; i < nDofs; ++i) {
      unsigned xDof = msh->GetSolutionDof(i, iel, xType);
      for (unsigned k = 0; k < dim; ++k)
        xv[k] = (*msh->_topology->_Sol[k])(xDof);

      // elementInB  = elementInB  && CheckIfInside(xv, boxB);
      elementInC  = elementInC  && CheckIfInside(xv, boxC);
      if (hasBd) elementInBd = elementInBd && CheckIfInside(xv, *boxBd);
    }

    // sol->_Sol[IndexB]->set(iel, elementInB ? 1.0 : 0.0);
    sol->_Sol[IndexC]->set(iel, elementInC ? 1.0 : 0.0);
    if (hasBd) sol->_Sol[IndexBd]->set(iel, elementInBd ? 1.0 : 0.0);
  }

  // sol->_Sol[IndexB]->close();
  sol->_Sol[IndexC]->close();
  if (hasBd) sol->_Sol[IndexBd]->close();

}





//Assemble Residual using A to update D amd V
void AssembleResAD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled


  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  TransientLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<TransientLinearImplicitSystem> ("LSAT");   // pointer to the linear implicit system named "Beam"
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

  double dt =  mlPdeSys->GetIntervalTime();
  double time =  mlPdeSys->GetTime();


  //solution variable
  // Note: for P* use PStar = IntegralP;
  // PStarNodes = std::move(WeightsP);
  unsigned solIndexZ = mlSol->GetIndex("Zi");
  unsigned solIndexX = mlSol->GetIndex("Xi");
  unsigned solIndexY = mlSol->GetIndex("Yi");

  unsigned solIndexW = mlSol->GetIndex("Wi");

  unsigned solIndexE = mlSol->GetIndex("Ei");

  unsigned solIndexP = mlSol->GetIndex("P");

  // unsigned solIndexB = mlSol->GetIndex("B");
  unsigned solIndexC = mlSol->GetIndex("C");
  unsigned solIndexBd;
  unsigned solIndexd;
  if(withDisturbance) {
    solIndexBd = mlSol->GetIndex("Bd");
    solIndexd = mlSol->GetIndex("d");
  }



  unsigned solType = mlSol->GetSolutionType(solIndexZ);


  unsigned solPdeIndexZ = mlPdeSys->GetSolPdeIndex("Zi");
  unsigned solPdeIndexX = mlPdeSys->GetSolPdeIndex("Xi");
  unsigned solPdeIndexY = mlPdeSys->GetSolPdeIndex("Yi");

  unsigned solPdeIndexW = mlPdeSys->GetSolPdeIndex("Wi");


  std::vector < double > solZOld;    // local solution

  std::vector < adept::adouble > solZ;    // local solution
  std::vector < adept::adouble > solX;
  std::vector < adept::adouble > solY;

  std::vector < adept::adouble > solW;

  adept::adouble PstarZ = 0;
  adept::adouble PstarX = 0;
  adept::adouble PstarY = 0;
  adept::adouble CstarCPstarZ = 0;
  adept::adouble CstarCPstarY = 0;
  double CstarCPstarP = 0;
  double CstarCPstarE = 0.;
  adept::adouble w = 0;
  adept::adouble x1 = 0;
  adept::adouble y1 = 0;

  std::vector < double > solE;
  std::vector < double > solP;
  std::vector < double > solBd;

  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector <double> phi;  // local test function for velocity
  std::vector <double> gradPhi; // local test function first order partial derivatives
  double weight; // gauss point weight

  std::vector < unsigned > sysDof; // local to global pdeSys dofs
  std::vector < adept::adouble > aRes;
  std::vector < double > res; // local redidual std::vector
  std::vector < double > Jac;

  RES->zero(); // Set to zero all the entries of the Global Residual std::vector
  KK->zero(); // Set to zero all the entries of the Global Matrix

  double wOld = 0.;

  for (unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    double CsC = (*sol->_Sol[solIndexC])(iel);

    unsigned nDofs = msh->GetElementDofNumber(iel, solType);

    solZ.resize(nDofs);
    solX.resize(nDofs);
    solY.resize(nDofs);
    solW.resize(nDofs);
    solP.resize(nDofs);
    solE.resize(nDofs);

    for (unsigned i = 0; i < nDofs; i++) {
      unsigned iDof = msh->GetSolutionDof(i, iel, coordXType);
      solP[i] = (*sol->_Sol[solIndexP])(iDof);
      solE[i] = (*sol->_Sol[solIndexE])(iDof);

      bool isX1node = (iDof == g_specialWnodes.X1);
      if(isX1node) x1 += solW[i];

      bool isW0node = (iDof == g_specialWnodes.W0);
      if(isW0node){
        w += solW[i];
        wOld = (*sol->_SolOld[solIndexW])(iDof);
      }

      bool isY1node = (iDof == g_specialWnodes.Y1);
      if(isY1node) y1 += solW[i];

      PstarZ += PStarNodes[iDof] * solZ[i];
      CstarCPstarZ += (CsC > 0.5) * PStarNodes[iDof] * solZ[i];
      PstarX += PStarNodes[iDof] * solX[i];
      PstarY += PStarNodes[iDof] * solY[i];
      CstarCPstarY += (CsC > 0.5) * PStarNodes[iDof] * solY[i];
      CstarCPstarP += (CsC > 0.5) * PStarNodes[iDof] * solP[i];
      CstarCPstarE += (CsC > 0.5) * PStarNodes[iDof] * solE[i];
    }
  }

  // element loop: each process loops only on the elements that owns
  for (unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);

    // double BBs = (*sol->_Sol[solIndexB])(iel);
    double CsC = (*sol->_Sol[solIndexC])(iel);
    double Bds = 0;
    if(withDisturbance) Bds = (*sol->_Sol[solIndexBd])(iel);;

    unsigned nDofs = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs

    unsigned nUnkn = 4;
    unsigned nDofsAll = nUnkn * nDofs;

    sysDof.resize(nDofsAll);
    aRes.assign(nDofsAll, 0.);

    solZOld.resize(nDofs);
    solZ.resize(nDofs);
    solX.resize(nDofs);
    solY.resize(nDofs);
    solW.resize(nDofs);
    solE.resize(nDofs);
    solP.resize(nDofs);
    if(withDisturbance) solBd.resize(nDofs);

    for (unsigned  k = 0; k < dim; k++) {
      coordX[k].resize(nDofs);
    }

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned iDof = msh->GetSolutionDof(i, iel, solType);

      solZOld[i] = (*sol->_SolOld[solIndexZ])(iDof);
      solZ[i] = (*sol->_Sol[solIndexZ])(iDof);
      solX[i] = (*sol->_Sol[solIndexX])(iDof);
      solY[i] = (*sol->_Sol[solIndexY])(iDof);

      solW[i] = (*sol->_Sol[solIndexW])(iDof);

      solP[i] = (*sol->_Sol[solIndexP])(iDof);

      solE[i] = (*sol->_Sol[solIndexE])(iDof);
      if(withDisturbance) solBd[i] = (*sol->_Sol[solIndexd])(iDof);

      for (unsigned k = 0; k < nUnkn; k++) {
        unsigned solIndex = (k == 0) ? solIndexZ :
                             (k == 1) ? solIndexX :
                             (k == 2) ? solIndexY : solIndexW;
        unsigned solPdeIndex = (k == 0) ? solPdeIndexZ :
                                (k == 1) ? solPdeIndexX :
                                (k == 2) ? solPdeIndexY : solPdeIndexW;
        sysDof[k * nDofs + i] = pdeSys->GetSystemDof(solIndex, solPdeIndex, i, iel);
      }


    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // local to global mapping between coordinates node and coordinate dof
      for (unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }

    s.new_recording();

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solType]->Jacobian(coordX, ig, weight, phi, gradPhi);

      double ZOldg = 0.;
      adept::adouble Zg = 0.;
      adept::adouble Xg = 0.;
      adept::adouble Yg = 0.;

      double r = 0.;
      double d = 0.;

      std::vector < adept::adouble > gradZg(dim, 0.);
      std::vector < adept::adouble > gradXg(dim, 0.);
      std::vector < adept::adouble > gradYg(dim, 0.);

      for (unsigned i = 0; i < nDofs; i++) {

        ZOldg += solZOld[i] * phi[i];
        Zg += solZ[i] * phi[i];
        Xg += solX[i] * phi[i];
        Yg += solY[i] * phi[i];

        r += solE[i] * phi[i];
        if(withDisturbance) d += solBd[i] * phi[i];

        for (unsigned j = 0; j < dim; j++) {
          gradZg[j] += solZ[i] * gradPhi[i * dim + j];
          gradXg[j] += solX[i] * gradPhi[i * dim + j];
          gradYg[j] += solY[i] * gradPhi[i * dim + j];
        }
      }

      double alpha = 1.0e-12;
      double beta = 0.125;
      double h = 1.;
      double a = -5.;
      double b = 5.;

      double mu = 1.;

      double TODOSolW = 0.;

      // *** phiA_i loop ***
      for (unsigned i = 0; i < nDofs; i++) {
        unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);

        adept::adouble aResZ = (Zg - ZOldg) / dt * phi[i] + h * a * solP[i] * w * phi[i] + h * b * solP[i] * (- h * b * PstarX + b * x1) / alpha * phi[i];
        adept::adouble aResX = - (CsC > 0.5) * (r - (1. - beta) * (Zg + h * solP[i] * w) - beta * (Yg + h * solP[i] * y1)) * phi[i];
        adept::adouble aResY =  + h * a * solP[i] * y1 * phi[i] + h * b * solP[i] * (- h * b * PstarX + b * x1) / alpha * phi[i];
        adept::adouble aResW = 0.0;


        for (unsigned j = 0; j < dim; j++) { // second index j in each equation
          aResZ +=  mu * gradPhi[i * dim + j] * gradZg[j]; // diffusion
          aResX +=  mu * gradPhi[i * dim + j] * gradXg[j]; // diffusion
          aResY +=  mu * gradPhi[i * dim + j] * gradYg[j]; // diffusion
        }


        // Looking for special 3 nodes to assemble the W solution
        unsigned iDof = msh->GetSolutionDof(i, iel, solType);

        bool isW0node = (iDof == g_specialWnodes.W0);
        bool isX1node = (iDof == g_specialWnodes.X1);
        bool isY1node = (iDof == g_specialWnodes.Y1);

        // b * TODOSolX1 / alpha //TODO
/*
        if (isW0node) {
          aResW = ( (1 / dt - a) * w - wOld / dt - b * b / alpha * (- h * PstarX + x1)) * phi[i]; // TODO equation for W
        }
        else if (isX1node) {
          aResW = ( h * a * PstarX - a * x1 - h * (CstarCPstarE - (1 - beta) * (CstarCPstarZ + h * CstarCPstarP * w) - beta * ( CstarCPstarY + h * CstarCPstarP * y1)) ) * phi[i]; // TODO equation for X1
        }
        else if (isY1node) {
          aResW = ( - a * y1 - b * b / alpha * (- h * PstarX + x1)) * phi[i]; // TODO equation for Y1
        }
        else {*/
          aResW = solW[i] * phi[i]; // identity
        // }

        aRes[0 * nDofs + i] += aResZ * weight;
        aRes[1 * nDofs + i] += aResX * weight;
        aRes[2 * nDofs + i] += aResY * weight;
        aRes[3 * nDofs + i] += aResW * weight;

      } // end phiA_i loop
    }

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    res.resize(nDofsAll);
    //copy the value of the adept::adoube mRes in double Res and store
    for (int i = 0; i < nDofsAll; i++) {
      res[i] = -aRes[i].value();
    }

    // define the dependent variables
    s.dependent(aRes.data(), nDofsAll);
    s.independent(solZ.data(), nDofs);
    s.independent(solX.data(), nDofs);
    s.independent(solY.data(), nDofs);
    s.independent(solW.data(), nDofs);

    Jac.assign(nDofsAll * nDofsAll, 0.);
    // get the jacobian matrix (ordered by column)
    s.jacobian(Jac.data(), true);

    RES->add_vector_blocked(res, sysDof);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();
  KK->close();
  //KK->draw();

  // double a;
  // std::cin>>a;

}


void AssembleResP(MultiLevelProblem& ml_prob) {
  adept::Stack& s = FemusInit::_adeptStack;

  auto* mlPdeSys = &ml_prob.get_system<LinearImplicitSystem>("LP");
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh* msh = ml_prob._ml_msh->GetLevel(level);
  MultiLevelSolution* mlSol = ml_prob._ml_sol;
  Solution* sol = ml_prob._ml_sol->GetSolutionLevel(level);

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];
  SparseMatrix* KK = pdeSys->_KK;
  NumericVector* RES = pdeSys->_RES;

  RES->zero();
  KK->zero();

  const unsigned pIndex    = mlSol->GetIndex("P");
  const unsigned pType     = mlSol->GetSolutionType(pIndex);
  const unsigned pPdeIndex = mlPdeSys->GetSolPdeIndex("P");

  const unsigned dim        = msh->GetDimension();
  const unsigned coordXType = 2;
  const unsigned iproc      = msh->processor_id();

  // ------------- element loop -------------
  for (unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; ++iel) {

    const short unsigned ielGeom = msh->GetElementType(iel);
    const unsigned nDofs = msh->GetElementDofNumber(iel, pType);

    // coordinates at element nodes
    std::vector<std::vector<double>> X(dim, std::vector<double>(nDofs));
    for (unsigned i = 0; i < nDofs; ++i) {
      const unsigned xd = msh->GetSolutionDof(i, iel, coordXType);
      for (unsigned k = 0; k < dim; ++k) X[k][i] = (*msh->_topology->_Sol[k])(xd);
    }

    // global dof map
    std::vector<unsigned> sysDof(nDofs);
    for (unsigned i = 0; i < nDofs; ++i) {
      sysDof[i] = pdeSys->GetSystemDof(pIndex, pPdeIndex, i, iel);
    }

    // local unknowns as adoubles
    std::vector<adept::adouble> p(nDofs);
    for (unsigned i = 0; i < nDofs; ++i) {
      const unsigned gd = msh->GetSolutionDof(i, iel, pType);
      p[i] = (*sol->_Sol[pIndex])(gd);
    }

    s.new_recording();

    // local residual
    std::vector<adept::adouble> aRes(nDofs, 0.0);

    // Gauss integration: a(u,v) = ∫ grad p · grad v
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][pType]->GetGaussPointNumber(); ++ig) {
      double w;
      std::vector<double> phi, dphi;
      msh->_finiteElement[ielGeom][pType]->Jacobian(X, ig, w, phi, dphi);

      // grad p at GP
      std::vector<adept::adouble> gradPg(dim, 0.0);
      for (unsigned a = 0; a < nDofs; ++a)
        for (unsigned k = 0; k < dim; ++k)
          gradPg[k] += p[a] * dphi[a*dim + k];

      // test with grad v_i
      for (unsigned i = 0; i < nDofs; ++i) {
        adept::adouble gi = 0.0;
        for (unsigned k = 0; k < dim; ++k) gi += dphi[i*dim + k] * gradPg[k];
        aRes[i] += gi * w;
      }
    }

    // assemble local residual and Jacobian
    std::vector<double> Re(nDofs);
    for (unsigned i = 0; i < nDofs; ++i) Re[i] = -aRes[i].value();   // move to RHS

    s.dependent(aRes.data(), nDofs);
    s.independent(p.data(),  nDofs);

    std::vector<double> Je(nDofs * nDofs, 0.0);
    s.jacobian(Je.data(), true);

    RES->add_vector_blocked(Re, sysDof);
    KK->add_matrix_blocked(Je, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();
  }

  // point-loads: δ at selected nodes (B = Σ_i δ_{x_i})
  for (unsigned gdof : g_controlNodeDofs) {
    RES->add(gdof, 1.0);
  }

  RES->close();
  KK->close();
}







std::pair<double, std::vector<double>> PrecomputePstarIntegrals(Solution* sol) {
  Mesh* msh = sol->GetMesh();
  const unsigned dim = msh->GetDimension();
  unsigned iproc = msh->processor_id();

  unsigned pIndex = sol->GetIndex("P");
  unsigned pType  = sol->GetSolutionType(pIndex);
  unsigned coordXType = 2;  // quadratic coordinates

  // global scalar integral ∫Ω P dx
  double localIntegral = 0.0;

  // nodal weights W_i = ∫Ω P φ_i dx
  std::vector<double> nodeWeights(sol->_Sol[pIndex]->size(), 0.0);

  // element loop (each process owns its range)
  for (unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; ++iel) {
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofs = msh->GetElementDofNumber(iel, pType);

    std::vector<std::vector<double>> coordX(dim, std::vector<double>(nDofs));
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned coordXDof = msh->GetSolutionDof(i, iel, coordXType);
      for (unsigned k = 0; k < dim; k++)
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);
    }

    std::vector<double> pVal(nDofs);
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned pDof = msh->GetSolutionDof(i, iel, pType);
      pVal[i] = (*sol->_Sol[pIndex])(pDof);
    }

    // Gauss integration
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][pType]->GetGaussPointNumber(); ig++) {
      double weight;
      std::vector<double> phi, gradPhi;
      msh->_finiteElement[ielGeom][pType]->Jacobian(coordX, ig, weight, phi, gradPhi);

      double Pg = 0.;
      for (unsigned i = 0; i < nDofs; i++) Pg += pVal[i] * phi[i];

      localIntegral += Pg * weight;

      for (unsigned i = 0; i < nDofs; i++) {
        unsigned pDof = msh->GetSolutionDof(i, iel, pType);
        nodeWeights[pDof] += Pg * phi[i] * weight;
      }
    }
  }

  // MPI reductions
  double globalIntegral = 0.0;
  MPI_Allreduce(&localIntegral, &globalIntegral, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  std::vector<double> globalNodeWeights(nodeWeights.size());
  MPI_Allreduce(nodeWeights.data(), globalNodeWeights.data(),
                nodeWeights.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return {globalIntegral, globalNodeWeights};
}




std::vector<double> ComputeL2NormCascadeOverC(Solution* sol, unsigned cascadeIterations) {
  Mesh* msh = sol->GetMesh();
  const unsigned dim = msh->GetDimension();
  unsigned iproc = msh->processor_id();

  unsigned solIndexC = sol->GetIndex("C");
  unsigned coordXType = 2;

  std::vector<unsigned> solIndexE(cascadeIterations);
  std::vector<unsigned> solTypeE(cascadeIterations);
  for (unsigned j = 0; j < cascadeIterations; j++) {
    std::string Ej = "E" + std::to_string(j);
    solIndexE[j] = sol->GetIndex(Ej.c_str());
    solTypeE[j] = sol->GetSolutionType(solIndexE[j]);
  }

  std::vector<double> local_integral(cascadeIterations, 0.0);

  for (unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    double CsC = (*sol->_Sol[solIndexC])(iel);
    if (CsC < 0.5) continue;  // outside C

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofs = msh->GetElementDofNumber(iel, solTypeE[0]);

    std::vector<std::vector<double>> coordX(dim, std::vector<double>(nDofs));
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned coordXDof = msh->GetSolutionDof(i, iel, coordXType);
      for (unsigned k = 0; k < dim; k++)
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);
    }

    // store E_j values for each cascade at element dofs
    std::vector<std::vector<double>> solE(cascadeIterations, std::vector<double>(nDofs));
    for (unsigned j = 0; j < cascadeIterations; j++) {
      for (unsigned i = 0; i < nDofs; i++) {
        unsigned eDof = msh->GetSolutionDof(i, iel, solTypeE[j]);
        solE[j][i] = (*sol->_Sol[solIndexE[j]])(eDof);
      }
    }

    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solTypeE[0]]->GetGaussPointNumber(); ig++) {
      double weight;
      std::vector<double> phi, gradPhi;
      msh->_finiteElement[ielGeom][solTypeE[0]]->Jacobian(coordX, ig, weight, phi, gradPhi);

      for (unsigned j = 0; j < cascadeIterations; j++) {
        double Eg = 0.;
        for (unsigned i = 0; i < nDofs; i++) Eg += solE[j][i] * phi[i];
        local_integral[j] += Eg * Eg * weight;
      }
    }
  }

  std::vector<double> global_integral(cascadeIterations, 0.0);
  MPI_Allreduce(local_integral.data(), global_integral.data(), cascadeIterations,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (double& val : global_integral) val = std::sqrt(val);
  return global_integral;
}


unsigned FindClosestNode(const Mesh* msh, const std::vector<double>& x0) {
  const unsigned dim = msh->GetDimension();
  unsigned coordXType = 2;  // coordinates always quadratic

  double minDist2 = std::numeric_limits<double>::max();
  unsigned closestDof = static_cast<unsigned>(-1);

  unsigned nNodes = msh->_topology->_Sol[0]->size();  // total number of coordinate dofs

  for (unsigned i = 0; i < nNodes; ++i) {
    std::vector<double> xv(dim);
    for (unsigned k = 0; k < dim; ++k)
      xv[k] = (*msh->_topology->_Sol[k])(i);

    double dist2 = 0.0;
    for (unsigned k = 0; k < dim; ++k) {
      double diff = xv[k] - x0[k];
      dist2 += diff * diff;
    }

    if (dist2 < minDist2) {
      minDist2 = dist2;
      closestDof = i;
    }
  }

  return closestDof;
}


std::vector<double> MarkControlNodes(const Mesh* msh, const std::vector<std::vector<double>>& points) {
  const unsigned dim = msh->GetDimension();
  unsigned nNodes = msh->_topology->_Sol[0]->size();

  std::vector<double> controlFlag(nNodes, 0.0);

  // helper lambda to compute squared distance
  auto dist2 = [&](unsigned node, const std::vector<double>& x0) {
    double d2 = 0.0;
    for (unsigned k = 0; k < dim; ++k) {
      double diff = (*msh->_topology->_Sol[k])(node) - x0[k];
      d2 += diff * diff;
    }
    return d2;
  };

  // loop over all input points
  for (const auto& x0 : points) {
    double minD2 = std::numeric_limits<double>::max();
    unsigned closest = 0;

    for (unsigned i = 0; i < nNodes; ++i) {
      double d2 = dist2(i, x0);
      if (d2 < minD2) {
        minD2 = d2;
        closest = i;
      }
    }

    controlFlag[closest] = 1.0;
  }

  return controlFlag;
}


std::vector<unsigned> GetControlNodeIndices(const Mesh* msh,
                                            const std::vector<std::vector<double>>& points) {
  const unsigned dim = msh->GetDimension();
  unsigned nNodes = msh->_topology->_Sol[0]->size();

  std::vector<unsigned> controlIndices;
  controlIndices.reserve(points.size());

  auto dist2 = [&](unsigned node, const std::vector<double>& x0) {
    double d2 = 0.0;
    for (unsigned k = 0; k < dim; ++k) {
      double diff = (*msh->_topology->_Sol[k])(node) - x0[k];
      d2 += diff * diff;
    }
    return d2;
  };

  for (const auto& x0 : points) {
    double minD2 = std::numeric_limits<double>::max();
    unsigned closest = 0;

    for (unsigned i = 0; i < nNodes; ++i) {
      double d2 = dist2(i, x0);
      if (d2 < minD2) {
        minD2 = d2;
        closest = i;
      }
    }

    controlIndices.push_back(closest);
  }

  return controlIndices;
}

WNodeIDs GetGlobalNodeIDsForW(const Mesh* msh, unsigned elemID) {
  const unsigned solType = 2; // quadratic coordinates

  WNodeIDs ids;
  ids.W0 = msh->GetSolutionDof(0, elemID, solType);
  ids.X1 = msh->GetSolutionDof(1, elemID, solType);
  ids.Y1 = msh->GetSolutionDof(2, elemID, solType);

  return ids;
}


