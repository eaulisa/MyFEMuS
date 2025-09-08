
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


// #include "CutFemWeight.hpp"
// #include "CDWeights.hpp"
// #include "Fem.hpp"
//
// typedef double TypeIO;
// typedef cpp_bin_float_oct TypeA;
// typedef cpp_bin_float_oct oct;
/*
// CutFemWeight <double, double> quad = CutFemWeight<double, double>(QUAD, 5, "legendre");
CutFemWeight <TypeIO, TypeA> quad  = CutFemWeight<TypeIO, TypeA >(QUAD, 1, "legendre");
CutFemWeight <TypeIO, TypeA> tri  = CutFemWeight<TypeIO, TypeA >(TRI, 1, "legendre");
Fem fem = Fem(tri.GetGaussQuadratureOrder(), tri.GetDimension());*/


#include "../include/MyMarker/MyMarker.hpp"
#include "../include/MyMarker/MyMarker.cpp"
//#include "../include/Cloud.hpp"
//#include "MyEigenFunctions.hpp"

#include <fstream>
#include <iostream>

using namespace femus;

void SetVelocity(Solution *sol, const std::vector<std::string> &U, const double &time, const double &T) {};

void LevelSetAdvection(const unsigned & stages, const std::vector<std::string> &U, std::string &PSI, const double & dt, Solution *_sol);

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet

  if(!strcmp(SolName, "U")) {  // strcmp compares two string in lexiographic sense.
    value = -x[1];
  }
  else if(!strcmp(SolName, "V")) {
    value = x[0];
  }
  else if(!strcmp(SolName, "W")) {
    value = 0.;
  }
  else if(!strcmp(SolName, "PSI")) {
    value = -1.;
  }
  return dirichlet;
}

double SetInitialCondition(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {

  double value = 0.;

  if(!strcmp(name, "U")) {
    value = -x[1];
  }
  else if(!strcmp(name, "V")) {
    value = x[0];
  }
  else if(!strcmp(name, "W")) {
    value = 0.;
  }
  else if(!strcmp(name, "PSI")) {
    double sigma = 0.18016836131796748;
    double r = 0.15;
    double xc = 0.;
    double yc = 0.25;
    value = exp( (r * r - x[0] * x[0] - (x[1] - yc) * (x[1] - yc) ) / (sigma * sigma) ) - 1.;
  }

  return value;
}


double GetAnalyticPSI(const std::vector < double >& x) {
  return -1.;
}



int main(int argc, char** args) {


  FemusInit mpinit(argc, args, MPI_COMM_WORLD);


  MultiLevelMesh mlMsh;

  double scalingFactor = 1.;

  //mlMsh.GenerateCoarseBoxMesh(256, 256, 0, -0.75, 0.75, -0.75, 0.75, 0., 0., TRI6, "seventh");
  mlMsh.GenerateCoarseBoxMesh(16 * 3 / 2, 16 * 3 / 2, 0, -0.75, 0.75, -0.75, 0.75, 0., 0., QUAD9, "seventh");

  unsigned dim = mlMsh.GetDimension();

  unsigned numberOfUniformLevels = 2;
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

  mlSol.AddSolution("PSI", LAGRANGE, SECOND, 2);


  MultiLevelProblem mlProb(&mlSol);

  mlSol.Initialize("All");
  mlSol.Initialize("U", SetInitialCondition, &mlProb);
  mlSol.Initialize("V", SetInitialCondition, &mlProb);
  if(dim == 3) mlSol.Initialize("W", SetInitialCondition, &mlProb);
  mlSol.Initialize("PSI", SetInitialCondition, &mlProb);

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

  mlSol.GenerateBdc("All");


  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");


  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(true);

  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);

  const unsigned level = mlMsh.GetNumberOfLevels() - 1;
  Solution* sol = mlSol.GetSolutionLevel(level);


  unsigned nIterations = 100;
  double time = 0;
  double period = 2 * M_PI;
  double dt = period/nIterations;

  std::vector<std::string> velocity = {"U", "V"};
  std::string levelSet = {"PSI"};

  unsigned stages = 4;

  for(unsigned it = 1; it <= nIterations; it++) {
    std::cout << "ITERATION " << it << "\n";
    sol->CopySolutionToOldSolution();
    time += dt;
    SetVelocity(sol, velocity, time, period);
    LevelSetAdvection(stages, velocity, levelSet, dt, sol);
    vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, it);
  }
  return 0;
}

double const XI[6][27][3] = {
  {
    //hex
    {-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1}, {-1, -1, 1}, {1, -1, 1}, {1, 1, 1}, {-1, 1, 1}, {0, -1, -1},
    {1, 0, -1}, {0, 1, -1}, {-1, 0, -1}, {0, -1, 1}, {1, 0, 1}, {0, 1, 1}, {-1, 0, 1}, {-1, -1, 0},
    {1, -1, 0}, {1, 1, 0}, {-1, 1, 0}, {0, -1, 0}, {1, 0, 0}, {0, 1, 0}, {-1, 0, 0}, {0, 0, -1},
    {0, 0, 1}, {0, 0, 0}
  },
  {
    // tet
    {0, 0, 0},        {1, 0, 0},       {0, 1, 0},   {0, 0, 1},            //0->4
    {0.5, 0, 0},      {0.5, 0.5, 0},   {0, 0.5, 0},
    {0.,  0, 0.5},    {0.5, 0., 0.5},  {0, 0.5, 0.5},                     //5->9
    {1. / 3., 1. / 3., 0.}, {1. / 3., 0., 1. / 3.}, {1. / 3., 1. / 3., 1. / 3.}, {0., 1. / 3., 1. / 3.}, //external faces of internal tetrahedra
    {0.25, 0.25, 0.25} //34
  },
  {
    //wedge
    {0, 0, -1},      {1, 0, -1},      {0, 1, -1},                     //vertici triangoli
    {0, 0, 1},       {1, 0, 1},       {0, 1, 1},
    {0.5, 0, -1},    {0.5, 0.5, -1},  {0, 0.5, -1},                   //midpoints triangoli
    {0.5, 0, 1},     {0.5, 0.5, 1},   {0, 0.5, 1},
    {0, 0, 0},       {1, 0, 0},      {0, 1, 0},                       //midpoints quadrati
    {0.5, 0, 0},     {0.5, 0.5, 0},    {0, 0.5, 0}, //0->17           //facce quadrati
    {1. / 3., 1. / 3., -1}, {1. / 3., 1. / 3., 1}, {1. / 3., 1. / 3., 0}
  },
  {
    //quad
    { -1, -1}, {1, -1}, {1, 1}, { -1, 1},
    { 0, -1}, {1, 0}, {0, 1}, { -1, 0}, {0, 0}
  },
  {
    //tri
    {0, 0},         {1, 0},         {0, 1},
    {0.5, 0},       {0.5, 0.5},     {0, 0.5},
    {1. / 3., 1. / 3.}
  },
  {
    //line
    {-1,}, {1,}, {0,}
  }
};


void LevelSetAdvection(const unsigned & stages, const std::vector<std::string> &U, std::string &PSI, const double & dt, Solution *_sol) {

  Mesh *msh = _sol->GetMesh();
  unsigned dim = msh->GetDimension();
  unsigned iproc = msh->processor_id();
  elem* el = msh->el;

  map<unsigned, bool> pSearch;

  unsigned cnt = 0;

  std::vector < unsigned > solUIndex(dim);
  for(unsigned k = 0; k < dim; k++) {
    solUIndex[k] = _sol->GetIndex(U[k].c_str());
  }

  unsigned solPSIIndex = _sol->GetIndex(PSI.c_str());

  unsigned solUType = _sol->GetSolutionType(solUIndex[0]);
  unsigned solPSIType = _sol->GetSolutionType(solPSIIndex);

  std::vector<std::vector<std::vector<double>>> solU(stages);
  std::vector<std::vector<std::vector<double>>> solUOld(stages);
  std::vector<std::vector <double>> F(stages);
  std::vector<std::vector <double>> X(stages);
  std::vector<unsigned> iel(stages);
  std::vector<unsigned> ielType(stages);
  std::vector<unsigned> nUDofs(stages);

  for(unsigned j  = 0; j < stages; j++) {
    solU[j].resize(dim);
    solUOld[j].resize(dim);
    F[j].resize(dim);
    X[j].resize(dim);
  }
  std::vector<std::vector<double>> c = {{}, {}, {}, {1., 0.5, 0.5, 0.}};
  std::vector<std::vector<std::vector<double>>> a = {
    {},
    {},
    {},
    {{}, {0.5}, {0., 0.5}, {0., 0., 1.}}
  };
  std::vector<std::vector<double>> b = {{}, {}, {}, {1. / 6., 1. / 3., 1. / 3., 1. / 6.}};

  std::vector<double> phi;

  MyMarker mrk;

  for(unsigned jel = msh->_elementOffset[iproc]; jel < msh->_elementOffset[iproc + 1]; jel++) {

    iel[0] = jel;
    ielType[0] = msh->GetElementType(iel[0]);

    //extract velocity quantities in the iel[0] element
    nUDofs[0] = msh->GetElementDofNumber(iel[0], solUType);
    for(unsigned k = 0; k < dim; k++) {
      solU[0][k].resize(nUDofs[0]);
      solUOld[0][k].resize(nUDofs[0]);
    }
    for(unsigned i = 0; i < nUDofs[0]; i++) {
      unsigned iDof = msh->GetSolutionDof(i, iel[0], solUType);
      for(unsigned k = 0; k < dim; k++) {
        solU[0][k][i] = (*_sol->_Sol[solUIndex[k]])(iDof);
        solUOld[0][k][i] = (*_sol->_SolOld[solUIndex[k]])(iDof);
      }
    }
    unsigned nPSIDofs = msh->GetElementDofNumber(jel, solPSIType);
    for(unsigned j = 0; j < nPSIDofs; j++) {

      unsigned jDof = msh->GetSolutionDof(j, jel, solPSIType);
      std::vector <double> xp(dim, 0.);
      unsigned idof  = msh->GetSolutionDof(j, jel, 2);    // local to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        xp[k] = (*msh->_topology->_Sol[k])(idof);      // global extraction and local storage for the element coordinates
      }
      X[0] = xp;

      std::vector <double> xi(dim, 0.);
      for(unsigned k = 0; k < dim; k++) xi[k] = XI[ielType[0]][j][k];
      msh->_finiteElement[ielType[0]][solUType]->GetPhi(phi, xi);

      F[0].assign(dim, 0);
      for(unsigned i = 0; i < nUDofs[0]; i++) {
        for(unsigned k = 0; k < dim; k++) {
          F[0][k] -= ((1. - c[stages - 1][0]) * solUOld[0][k][i] + c[stages - 1][0] * solU[0][k][i]) * phi[i];
        }
      }

      bool insideLocalDomain = true;
      for(unsigned rk = 1; rk < stages; rk++) {
        X[rk] = X[0];
        for(unsigned jk = 0; jk < rk; jk++) {
          for(unsigned k = 0; k < dim; k++) {
            X[rk][k] += a[stages - 1][rk][jk] * F[jk][k] * dt;
          }
        }
        insideLocalDomain = mrk.SerialElementSearchWithInverseMapping(X[rk], _sol, solUType, iel[rk - 1]);
        if(!insideLocalDomain) {
          unsigned kel = mrk.GetElement();
          if(kel == UINT_MAX) {
            _sol->_Sol[solPSIIndex]->set(jDof, GetAnalyticPSI(X[0]));
          }
          else {
            pSearch[jDof] = true;
          }
          break;
        }
        iel[rk] = mrk.GetElement();
        bool sameElement = false;
        for(unsigned jk = 0; jk < rk; jk++) {
          if(iel[rk] == iel[jk]) {
            sameElement = true;
            ielType[rk] = ielType[jk];
            nUDofs[rk] =  nUDofs[jk];
            solU[rk] = solU[jk];
            solUOld[rk] = solUOld[jk];
            break;
          }
        }
        if(sameElement == false) {
          ielType[rk] = msh->GetElementType(iel[rk]);
          nUDofs[rk] = msh->GetElementDofNumber(iel[rk], solUType);
          for(unsigned k = 0; k < dim; k++) {
            solU[rk][k].resize(nUDofs[rk]);
            solUOld[rk][k].resize(nUDofs[rk]);
          }
          for(unsigned i = 0; i < nUDofs[rk]; i++) {
            unsigned iDof = msh->GetSolutionDof(i, iel[rk], solUType);
            for(unsigned k = 0; k < dim; k++) {
              solU[rk][k][i] = (*_sol->_Sol[solUIndex[k]])(iDof);
              solUOld[rk][k][i] = (*_sol->_SolOld[solUIndex[k]])(iDof);
            }
          }
        }
        msh->_finiteElement[ielType[rk]][solUType]->GetPhi(phi, mrk.GetIprocLocalCoordinates());
        F[rk].assign(dim, 0);
        for(unsigned i = 0; i < nUDofs[rk]; i++) {
          for(unsigned k = 0; k < dim; k++) {
            F[rk][k] -= ((1. - c[stages - 1][rk]) * solUOld[rk][k][i] + c[stages - 1][rk] * solU[rk][k][i]) * phi[i];
          }
        }
      }
      if(insideLocalDomain) {

        for(unsigned k = 0; k < dim; k++) {
          for(unsigned rk = 0; rk < stages; rk++) {
            xp[k] += b[stages - 1][rk] * F[rk][k] * dt;
          }
        }
        insideLocalDomain = mrk.SerialElementSearchWithInverseMapping(xp, _sol, solUType, iel[0]);
        unsigned  kel = mrk.GetElement();
        if(insideLocalDomain) {
          unsigned kelType = msh->GetElementType(kel);
          msh->_finiteElement[kelType][solPSIType]->GetPhi(phi, mrk.GetIprocLocalCoordinates());
          unsigned nDofsPSIKel = msh->GetElementDofNumber(kel, solPSIType);
          double PSI = 0.;
          for(unsigned i = 0; i < nDofsPSIKel; i++) {
            unsigned iDof = msh->GetSolutionDof(i, kel, solPSIType);

            PSI += (*_sol->_SolOld[solPSIIndex])(iDof) * phi[i];

          }
          _sol->_Sol[solPSIIndex]->set(jDof, PSI);

        }
        else {

          if(kel == UINT_MAX) {
            _sol->_Sol[solPSIIndex]->set(jDof,  GetAnalyticPSI(X[0]));
          }
          else {
            pSearch[jDof]  = true;
          }
        }
      }
    }
  }
// _ypNew.resize(cnt);
// _elemNew.resize(cnt);
// _NNew.resize(cnt);
// _kappaNew.resize(cnt);
// _dsNew.resize(cnt);
//
// map<unsigned, bool>::iterator it;
//
// unsigned iproc = _sol->processor_id();
// unsigned nprocs = _sol->n_processors();

// if(nprocs > 1) {
//   for(unsigned kp = 0; kp < nprocs; kp++) {
//
//     unsigned np;
//     if(iproc == kp) {
//       np = pSearch.size();
//     }
//     MPI_Bcast(&np, 1, MPI_UNSIGNED, kp, MPI_COMM_WORLD);
//
//     if(np > 0) {
//       if(iproc == kp) {
//         it =  pSearch.begin();
//       }
//
//       for(unsigned jcnt = 0; jcnt < np; jcnt++) {
//         unsigned j;
//         if(iproc == kp) {
//           j = it->first;
//
//           iel[0] = _elem[_map[j]];
//           ielType[0] = msh->GetElementType(iel[0]);
//           nDofs[0] = msh->GetElementDofNumber(iel[0], solType);
//
//           for(unsigned k = 0; k < dim; k++) {
//             solU[0][k].resize(nDofs[0]);
//             solUOld[0][k].resize(nDofs[0]);
//           }
//
//           for(unsigned i = 0; i < nDofs[0]; i++) {
//             unsigned iDof = msh->GetSolutionDof(i, iel[0], solType);
//             for(unsigned k = 0; k < dim; k++) {
//               solU[0][k][i] = (*_sol->_Sol[solUIndex[k]])(iDof);
//               solUOld[0][k][i] = (*_sol->_SolOld[solUIndex[k]])(iDof);
//             }
//           }
//           msh->_finiteElement[ielType[0]][solType]->GetPhi(phi, _yi[_map[j]]);
//           X[0] = _yp[_map[j]];
//           F[0].assign(dim, 0);
//           for(unsigned i = 0; i < nDofs[0]; i++) {
//             for(unsigned k = 0; k < dim; k++) {
//               F[0][k] += ((1. - c[stages - 1][0]) * solUOld[0][k][i] + c[stages - 1][0] * solU[0][k][i]) * phi[i];
//             }
//           }
//         }
//         MPI_Bcast(&iel[0], 1, MPI_UNSIGNED, kp, MPI_COMM_WORLD);
//         MPI_Bcast(X[0].data(), X[0].size(), MPI_DOUBLE, kp, MPI_COMM_WORLD);
//         MPI_Bcast(F[0].data(), F[0].size(), MPI_DOUBLE, kp, MPI_COMM_WORLD);
//
//         bool insideLocalDomain = true;
//         for(unsigned rk = 1; rk < stages; rk++) {
//           X[rk] = X[0];
//           for(unsigned jk = 0; jk < rk; jk++) {
//             for(unsigned k = 0; k < dim; k++) {
//               X[rk][k] += a[stages - 1][rk][jk] * F[jk][k] * dt;
//             }
//           }
//           insideLocalDomain = _mrk.ParallelElementSearchWithInverseMapping(X[rk], _sol, solType, iel[rk - 1]);
//           if(!insideLocalDomain) {
//             break;
//           }
//
//           iel[rk] = _mrk.GetElement();
//           if(_mrk.GetProc() == iproc) {
//             ielType[rk] = msh->GetElementType(iel[rk]);
//             nDofs[rk] = msh->GetElementDofNumber(iel[rk], solType);
//             for(unsigned k = 0; k < dim; k++) {
//               solU[rk][k].resize(nDofs[rk]);
//               solUOld[rk][k].resize(nDofs[rk]);
//             }
//             for(unsigned i = 0; i < nDofs[rk]; i++) {
//               unsigned iDof = msh->GetSolutionDof(i, iel[rk], solType);
//               for(unsigned k = 0; k < dim; k++) {
//                 solU[rk][k][i] = (*_sol->_Sol[solUIndex[k]])(iDof);
//                 solUOld[rk][k][i] = (*_sol->_SolOld[solUIndex[k]])(iDof);
//               }
//             }
//
//             msh->_finiteElement[ielType[rk]][solType]->GetPhi(phi, _mrk.GetIprocLocalCoordinates());
//             F[rk].assign(dim, 0);
//             for(unsigned i = 0; i < nDofs[rk]; i++) {
//               for(unsigned k = 0; k < dim; k++) {
//                 F[rk][k] += ((1. - c[stages - 1][rk]) * solUOld[rk][k][i] + c[stages - 1][rk] * solU[rk][k][i]) * phi[i];
//               }
//             }
//           }
//           MPI_Bcast(F[rk].data(), F[rk].size(), MPI_DOUBLE, _mrk.GetProc(), MPI_COMM_WORLD);
//         }
//
//         if(insideLocalDomain) {
//           std::vector<double> ypNew = X[0];
//           for(unsigned k = 0; k < dim; k++) {
//             for(unsigned rk = 0; rk < stages; rk++) {
//               ypNew[k] += b[stages - 1][rk] * F[rk][k] * dt;
//             }
//           }
//           insideLocalDomain = _mrk.ParallelElementSearch(ypNew, _sol, solType, iel[0]);
//           if(insideLocalDomain) {
//             if(kp == iproc) {
//               MPI_Send(_N[_map[j]].data(), _N[_map[j]].size(), MPI_DOUBLE, _mrk.GetProc(), 1, MPI_COMM_WORLD);
//               MPI_Send(&_kappa[_map[j]], 1, MPI_DOUBLE, _mrk.GetProc(), 2, MPI_COMM_WORLD);
//               MPI_Send(&_ds[_map[j]], 1, MPI_DOUBLE, _mrk.GetProc(), 3, MPI_COMM_WORLD);
//             }
//
//             if(_mrk.GetProc() == iproc) {
//               _ypNew.resize(cnt + 1);
//               _ypNew[cnt] = ypNew;
//               _elemNew.resize(cnt + 1);
//               _elemNew[cnt] = _mrk.GetElement();
//               _NNew.resize(cnt + 1, std::vector<double> (dim));
//               MPI_Recv(_NNew[cnt].data(), _NNew[cnt].size(), MPI_DOUBLE, kp, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//               _kappaNew.resize(cnt + 1);
//               MPI_Recv(&_kappaNew[cnt], 1, MPI_DOUBLE, kp, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//               _dsNew.resize(cnt + 1);
//               MPI_Recv(&_dsNew[cnt], 1, MPI_DOUBLE, kp, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//               cnt++;
//             }
//           }
//         }
//         if(iproc == kp) {
//           it++;
//         }
//       }
//     }
//   }
// }
//
//
//
//
// _yp.swap(_ypNew);
// _elem.swap(_elemNew);
// _kappa.swap(_kappaNew);
// _ds.swap(_dsNew);
// _N.swap(_NNew);
//
// _yi.assign(cnt, std::vector<double>(dim, 0));
//
// CreateMap();


  _sol->_SolOld[solPSIIndex]->close();
}










