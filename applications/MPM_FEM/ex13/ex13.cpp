#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "adept.h"

#include "PetscMatrix.hpp"

using namespace femus;

bool SetBoundaryCondition (const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = false; 
  value = 0.;

  return dirichlet;
}

double InitalValueU (const std::vector < double >& x) {
  return x[0] * x[0] + x[1] * x[1];
}

void BuidProjection (MultiLevelProblem& ml_prob);

void ProjectSolutionIntoGradient (MultiLevelProblem& ml_prob);

int main (int argc, char** args) {

  FemusInit mpinit (argc, args, MPI_COMM_WORLD);


  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;

  //mlMsh.GenerateCoarseBoxMesh (2., 0, 0, -0.5, 0.5, 0., 0., 0., 0., EDGE3, "seventh");

  //mlMsh.ReadCoarseMesh ("./input/square_quad.neu", "seventh", scalingFactor);
  mlMsh.ReadCoarseMesh ("./input/square_tri.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/square_mixed.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/cube_hex.neu","seventh",scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/cube_wedge.neu","seventh",scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/cube_tet.neu","seventh",scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/cube_mixed.neu","seventh",scalingFactor);

  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
   *    probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();
  unsigned maxNumberOfMeshes;

  if (dim == 1) {
    maxNumberOfMeshes = 10;
  }
  else if (dim == 2) {
    maxNumberOfMeshes = 7;
  }
  else {
    maxNumberOfMeshes = 6;
  }

  vector < vector < double > > l2Norm;
  l2Norm.resize (maxNumberOfMeshes);

  vector < vector < double > > semiNorm;
  semiNorm.resize (maxNumberOfMeshes);

  for (unsigned i = 1; i < maxNumberOfMeshes; i++) {   // loop on the mesh level

    unsigned numberOfUniformLevels = i ;
    unsigned numberOfSelectiveLevels = 0;
    mlMsh.RefineMesh (numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

    // erase all the coarse mesh levels
    mlMsh.EraseCoarseLevels (numberOfUniformLevels - 1);

    // print mesh info
    mlMsh.PrintInfo();

    FEOrder feOrder[3] = {FIRST, SERENDIPITY, SECOND};
    l2Norm[i].resize (3);
    semiNorm[i].resize (3);

    for (unsigned j = 0; j < 3; j++) {   // loop on the FE Order
      // define the multilevel solution and attach the mlMsh object to it
      MultiLevelSolution mlSol (&mlMsh);



      // add variables to mlSol
      mlSol.AddSolution ("u", LAGRANGE, feOrder[j]);

      std::string Uxname[3] = {"ux", "uy", "uz"};
      for (unsigned k = 0; k < dim; k++) {
        mlSol.AddSolution (Uxname[k].c_str(), LAGRANGE, feOrder[j], 0, false);
      }
      
      mlSol.AddSolution ("weight", LAGRANGE, feOrder[j], 0, false);

      mlSol.Initialize ("All");
      mlSol.Initialize("u", InitalValueU);

      // attach the boundary condition function and generate boundary data
      mlSol.AttachSetBoundaryConditionFunction (SetBoundaryCondition);
      mlSol.GenerateBdc ("All");

      // define the multilevel problem attach the mlSol object to it
      MultiLevelProblem mlProb (&mlSol);

      LinearImplicitSystem* systemP[3];
      std::string Pname[3] = {"Px", "Py", "Pz"};
      for (unsigned k = 0; k < dim; k++) {
        systemP[k] = &mlProb.add_system < LinearImplicitSystem > (Pname[k]);
        systemP[k]->AddSolutionToSystemPDE ("u");
        systemP[k]->init();
      }
      
      BuidProjection (mlProb);
      ProjectSolutionIntoGradient (mlProb);

      std::vector < std::string > variablesToBePrinted;
      variablesToBePrinted.push_back ("All");

      VTKWriter vtkIO (&mlSol);
      vtkIO.SetDebugOutput (true);
      vtkIO.Write (DEFAULT_OUTPUTDIR, "linear", variablesToBePrinted, i + j * 10);
      vtkIO.Write (DEFAULT_OUTPUTDIR, "quadratic", variablesToBePrinted, i + j * 10);
      vtkIO.Write (DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, i + j * 10);
    }
  }

  return 1;
}




void BuidProjection (MultiLevelProblem& ml_prob) {

  double scale = 4.;

  adept::Stack& s = FemusInit::_adeptStack;

  MultiLevelSolution*  mlSol = ml_prob._ml_sol;
  unsigned level = mlSol->_mlMesh->GetNumberOfLevels() - 1u;

  Solution* sol = ml_prob._ml_sol->GetSolutionLevel (level);
  Mesh* msh = ml_prob._ml_msh->GetLevel (level);
  elem* el = msh->el;

  unsigned  dim = msh->GetDimension();

  std::vector < LinearImplicitSystem* > mlSysP (dim);
  std::vector < LinearEquationSolver* > sysP (dim);
  std::vector < SparseMatrix*> P (dim);

  std::string Pname[3] = {"Px", "Py", "Pz"};
  std::string Uname = "u";
  unsigned soluIndex = mlSol->GetIndex (Uname.c_str());
  for (unsigned k = 0; k < dim; k++) {
    mlSysP[k] =  &ml_prob.get_system< LinearImplicitSystem > (Pname[k]);
    sysP[k] = mlSysP[k]->_LinSolver[level];
    P[k] = sysP[k]->_KK;
    P[k]->zero();
  }

  vector < vector < double > > x (dim);
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector< int > sysDof;
  vector <double> phi;
  double* phi2;
  vector <double> phi_x;
  double weight;

  unsigned    iproc = msh->processor_id();

  unsigned solwIndex = mlSol->GetIndex ("weight");
  unsigned solType = mlSol->GetSolutionType (solwIndex);
  sol->_Sol[solwIndex]->zero();

  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nDofs  = msh->GetElementDofNumber (iel, solType);
    sysDof.resize (nDofs);
    for (int k = 0; k < dim; k++) {
      x[k].resize (nDofs);
    }
    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof (i, iel, solType);
      sysDof[i] = msh->GetSolutionDof (i, iel, solType);
    }
    // local storage of coordinates
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned xDof  = msh->GetSolutionDof (i, iel, xType);
      for (unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k]) (xDof);
      }
    }

    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {
      msh->_finiteElement[ielGeom][solType]->Jacobian (x, ig, weight, phi, phi_x);
      phi2 = (solType != 1) ? msh->_finiteElement[ielGeom][solType]->GetPhi (ig) : msh->_finiteElement[ielGeom][2]->GetPhi (ig);
      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofs; i++) {
        sol->_Sol[solwIndex]->add (sysDof[i], phi2[i] * weight);
      } // end phi_i loop
    } // end gauss point loop

  } //end element loop for each process*/
  sol->_Sol[solwIndex]->close();


  //solution variable

  
  {
    unsigned soluType = mlSol->GetSolutionType (soluIndex);
    if (soluType != solType) {
      std::cout << "error weight and u should be of the same type\n";
      abort();
    }
  }
  unsigned soluPdeIndex = mlSysP[0]->GetSolPdeIndex ("u");
  std::vector < adept::adouble > solu;

  std::vector <double> Jac;
  std::vector < std::vector< adept::adouble > > aRes (dim); // local redidual vector

  std::vector < double > solw;

  //BEGIN element loop
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nDofs  = msh->GetElementDofNumber (iel, solType);
    solu.resize (nDofs);
    solw.resize (nDofs);
    sysDof.resize (nDofs);
    for (int k = 0; k < dim; k++) {
      x[k].resize (nDofs);
      aRes[k].assign (nDofs, 0.);
    }
    Jac.resize (nDofs * nDofs);
    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof (i, iel, solType);
      solu[i] = (*sol->_Sol[soluIndex]) (solDof);
      solw[i] = (*sol->_Sol[solwIndex]) (solDof);
      sysDof[i] = sysP[0]->GetSystemDof (soluIndex, soluPdeIndex, i, iel);
    }
    // local storage of coordinates
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned xDof  = msh->GetSolutionDof (i, iel, xType);
      for (unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k]) (xDof);
      }
    }

    s.new_recording();
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {

      msh->_finiteElement[ielGeom][solType]->Jacobian (x, ig, weight, phi, phi_x);
      phi2 = (solType != 1) ? msh->_finiteElement[ielGeom][solType]->GetPhi (ig) : msh->_finiteElement[ielGeom][2]->GetPhi (ig);
      std::vector < adept::adouble > solux_g (dim, 0.);
      for (unsigned i = 0; i < nDofs; i++) {
        for (unsigned k = 0; k < dim; k++) {
          solux_g[k] += phi_x[i * dim + k] * solu[i];
        }
      }

      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofs; i++) {
        for (unsigned k = 0; k < dim; k++) {
          aRes[k][i] += solux_g[k] * phi2[i] * weight / solw[i]; //smoothed gradient part.
        }
      } // end phi_i loop
    } // end gauss point loop

    s.independent (&solu[0], nDofs);
    for (unsigned k = 0; k < dim; k++) {
      s.dependent (&aRes[k][0], nDofs);
      s.jacobian (&Jac[0], true);
      P[k]->add_matrix_blocked (Jac, sysDof, sysDof);
      s.clear_dependents();
    }
    s.clear_independents();
  } //end element loop for each process*/

  for (unsigned k = 0; k < dim; k++) {
    P[k]->close();
  }

}

void ProjectSolutionIntoGradient (MultiLevelProblem& ml_prob) {

  MultiLevelSolution*  mlSol = ml_prob._ml_sol;
  unsigned level = mlSol->_mlMesh->GetNumberOfLevels() - 1u;

  Solution* sol = ml_prob._ml_sol->GetSolutionLevel (level);
  Mesh* msh = ml_prob._ml_msh->GetLevel (level);
  elem* el = msh->el;

  unsigned  dim = msh->GetDimension();

  std::vector < LinearImplicitSystem* > mlSysP (dim + 1);
  std::vector < LinearEquationSolver* > sysP (dim + 1);
  std::vector < SparseMatrix*> P (dim + 1);

  std::string Pname[3] = {"Px", "Py", "Pz"};
  std::string Uxname[3] = {"ux", "uy", "uz"};
  std::vector < unsigned > solIndex (dim);
  for (unsigned k = 0; k < dim; k++) {
    mlSysP[k] =  &ml_prob.get_system< LinearImplicitSystem > (Pname[k]);
    sysP[k] = mlSysP[k]->_LinSolver[level];
    solIndex[k] = mlSol->GetIndex (Uxname[k].c_str());
    P[k] = sysP[k]->_KK;
  }

  unsigned soluIndex = mlSol->GetIndex ("u");

  for (unsigned k = 0; k < dim; k++) {
    (*sol->_Sol[solIndex[k]]).matrix_mult ( (*sol->_Sol[soluIndex]), (*P[k]));
  }

}

