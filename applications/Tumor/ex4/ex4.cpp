/*
First two equations of the porous media-system; solve for disp-u and pressure-p
-u_t*gradPhi+(1/mu_w)*(K*p)*gradPhi-F*Phi=0
-[nu*(gradu+gradu^T):gradw+lamda*(divu)*divw-p*divw]-(p*I*n)w_{boundary} = 0
for test space(phi,w)
*/

#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "TransientSystem.hpp"
#include "adept.h"
#include <math.h>


using namespace femus;
unsigned DIM = 2;
double dt = 0.01;
double beta = 2.;
double mu_w = 9.11e-4;
double k0 = 1.82e-15;
double E = 5000;
double nu = 0.4;
double mu = E / (2. * (1. + nu));
double lambda = ((1. + nu) * (1. - 2.*nu)) / (E * nu);

double GetTimeStep(const double time) {
  return dt;
}


bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true;
  value = 0.;

  if(!strcmp(SolName, "P")) {   // left boundary condition.
    value = 400;
  }

  return dirichlet;
}


void AssemblePorousMediaSystem(MultiLevelProblem& ml_prob);    //, unsigned level, const unsigned &levelMax, const bool &assembleMatrix );


int main(int argc, char** args) {



  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;

  unsigned nx = 10;
  unsigned ny = 2;
  unsigned nz = 1;

  double length = .1;
  double lengthx = 1;

  if(DIM == 2) {
    mlMsh.GenerateCoarseBoxMesh(nx, ny, 0, 0., lengthx, 0., length, 0., 0., QUAD9, "seventh");
  }
  else if(DIM == 3) {
    nz = ny;
    mlMsh.GenerateCoarseBoxMesh(nx, ny, nz, 0., lengthx, 0., length, 0., length,  HEX27, "seventh");
  }

  unsigned dim = mlMsh.GetDimension();

  unsigned numberOfUniformLevels = 4;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels, numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // erase all the coarse mesh levels
  //mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("DX", LAGRANGE, SECOND,2);
  mlSol.AddSolution("DY", LAGRANGE, SECOND,2);
  if(dim == 3) mlSol.AddSolution("DZ", LAGRANGE, SECOND,2);
  mlSol.AddSolution("P",  LAGRANGE, FIRST);

  //  Taylor-hood
  //  mlSol.AddSolution("U", LAGRANGE, SERENDIPITY);
  //  mlSol.AddSolution("V", LAGRANGE, SERENDIPITY);
  //  if (dim == 3) mlSol.AddSolution("W", LAGRANGE, SERENDIPITY);
  //  mlSol.AddSolution("P", LAGRANGE, FIRST);


  mlSol.Initialize("All");

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  //mlSol.FixSolutionAtOnePoint("P");

  mlSol.GenerateBdc("All");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);

  // add system Poisson in mlProb as a Linear Implicit System
  TransientLinearImplicitSystem & system = mlProb.add_system< TransientLinearImplicitSystem > ("PorousMedia");

  // add solution "u" to system
  system.AddSolutionToSystemPDE("DX");
  system.AddSolutionToSystemPDE("DY");
  if(dim == 3) system.AddSolutionToSystemPDE("DZ");
  system.AddSolutionToSystemPDE("P");

  // attach the assembling function to system
  system.SetAssembleFunction(AssemblePorousMediaSystem);
  //system.AttachGetTimeIntervalFunction(SetVariableTimeStep);

  // initilaize and solve the system
  system.init();
  //system.SetOuterSolver(PREONLY);

  system.AttachGetTimeIntervalFunction(GetTimeStep);
  const unsigned int n_timesteps = 60;


  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(true);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

   
  

  for(unsigned t = 1; t <= n_timesteps; t++) {

    system.CopySolutionToOldSolution();
    system.MGsolve();
    mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, t);


  }


  return 0;
}


void AssemblePorousMediaSystem(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<NonLinearImplicitSystem> ("PorousMedia");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  //solution variable
  vector < unsigned > solDIndex(dim);
  solDIndex[0] = mlSol->GetIndex("DX");    // get the position of "U" in the ml_sol object
  solDIndex[1] = mlSol->GetIndex("DY");    // get the position of "V" in the ml_sol object

  if(dim == 3) solDIndex[2] = mlSol->GetIndex("DZ");       // get the position of "V" in the ml_sol object

  unsigned solDType = mlSol->GetSolutionType(solDIndex[0]);    // get the finite element type for "u"

  unsigned solPIndex;
  solPIndex = mlSol->GetIndex("P");                         // get the position of "P" in the ml_sol object
  unsigned solPType = mlSol->GetSolutionType(solPIndex);    // get the finite element type for "u"

  vector < unsigned > solDPdeIndex(dim);
  solDPdeIndex[0] = mlPdeSys->GetSolPdeIndex("DX");    // get the position of "DX" in the pdeSys object
  solDPdeIndex[1] = mlPdeSys->GetSolPdeIndex("DY");    // get the position of "DY" in the pdeSys object

  if(dim == 3) solDPdeIndex[2] = mlPdeSys->GetSolPdeIndex("DZ");

  unsigned solPPdeIndex;
  solPPdeIndex = mlPdeSys->GetSolPdeIndex("P");    // get the position of "P" in the pdeSys object



  vector < vector < double > >  solDOld(dim);  // local solution
  vector < vector < adept::adouble > >  solD(dim);      // local solution
  vector < adept::adouble >  solP; // local solution

  vector< vector < adept::adouble > > aResD(dim);    // local redidual vector
  vector< adept::adouble > aResP; // local redidual vector, scalar one

  vector < vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for(unsigned  k = 0; k < dim; k++) {
    solD[k].reserve(maxSize);
    solDOld[k].reserve(maxSize);
    aResD[k].reserve(maxSize);
    coordX[k].reserve(maxSize);
  }

  solP.reserve(maxSize);
  aResP.reserve(maxSize);


  vector <double> phiD;  // local test function for velocity
  vector < double > gradPhiD;

  phiD.reserve(maxSize);
  gradPhiD.reserve(maxSize * dim);

  vector <double> phiP; // local test function for the pressure
  vector < double> gradPhiP;

  double weight; // gauss point weight

  vector< int > sysDof; // local to global pdeSys dofs
  sysDof.reserve((dim + 1) * maxSize);

  vector< double > Res; // local redidual vector
  Res.reserve((dim + 1) * maxSize);

  vector < double > Jac;
  Jac.reserve((dim + 1) * maxSize * (dim + 1) * maxSize);

  RES->zero(); // Set to zero all the entries of the Global Residual vector
  KK->zero(); // Set to zero all the entries of the Global Matrix

  
  
  
  // element loop: each process loops only on the elements that owns
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);

    unsigned nDofsD = msh->GetElementDofNumber(iel, solDType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);    // number of solution element dofs

    unsigned nDofsAll = dim * nDofsD + nDofsP;

    sysDof.resize(nDofsAll);
    Res.assign(nDofsAll, 0.);
    Jac.assign(nDofsAll * nDofsAll, 0.);

    
    
    for(unsigned  k = 0; k < dim; k++) {
      solDOld[k].resize(nDofsD);
      solD[k].resize(nDofsD);
      coordX[k].resize(nDofsD);
    }
    solP.resize(nDofsP);

    std::cout << "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"<< nDofsD  << std::endl;
    
    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsD; i++) {
      unsigned iDof = msh->GetSolutionDof(i, iel, solDType);    // local to global mapping between solution node and solution dof
       
      for(unsigned  k = 0; k < dim; k++) {
          
        solDOld[k][i] = (*sol->_SolOld[solDIndex[k]])(iDof);
         std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" << std::endl;
        solD[k][i] = (*sol->_Sol[solDIndex[k]])(iDof);      // global extraction and local storage for the solution
        sysDof[k * nDofsD + i] = pdeSys->GetSystemDof(solDIndex[k], solDPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
      
    }

    for(unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);    // local to global mapping between solution node and solution dof
      solP[i] = (*sol->_Sol[solPIndex])(solPDof);      // global extraction and local storage for the solution
      sysDof[dim * nDofsD + i ] = pdeSys->GetSystemDof(solPIndex, solPPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

   
    // local storage of coordinates
    for(unsigned i = 0; i < nDofsD; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // local to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }


    
    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Face Gauss point loop (boundary Integral) ***
//     for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
//       int faceIndex = el->GetBoundaryIndex(iel, jface);
//       // look for boundary faces
//       if(faceIndex == 3) {
//         const unsigned faceGeom = msh->GetElementFaceType(iel, jface);
//         unsigned faceDofs = msh->GetElementFaceDofNumber(iel, jface, solDType);
//
//         vector  < vector  <  double> > faceCoordinates(dim);    // A matrix holding the face coordinates rowwise.
//         for(int k = 0; k < dim; k++) {
//           faceCoordinates[k].resize(faceDofs);
//         }
//         for(unsigned i = 0; i < faceDofs; i++) {
//           unsigned inode = msh->GetLocalFaceVertexIndex(iel, jface, i);    // face-to-element local node mapping.
//           for(unsigned k = 0; k < dim; k++) {
//             faceCoordinates[k][i] =  coordX[k][inode]; // We extract the local coordinates on the face from local coordinates on the element.
//           }
//         }
//         for(unsigned ig = 0; ig  <  msh->_finiteElement[faceGeom][solDType]->GetGaussPointNumber(); ig++) {
//           // We call the method GetGaussPointNumber from the object finiteElement in the mesh object msh.
//
//           std::vector < double > normal;
//           //msh->_finiteElement[faceGeom][solDType]->JacobianSur(faceCoordinates, ig, weight, phiD, gradPhiD, normal);
//
//
//           msh->_finiteElement[faceGeom][solPType]->JacobianSur(faceCoordinates, ig, weight, phiP, gradPhiP, normal);
//           //phiP = msh->_finiteElement[faceGeom][solPType]->GetPhi(ig);
//
//
//           vector< double > xg(dim, 0.);
//           for(unsigned i = 0; i < faceDofs; i++) {
//             for(unsigned k = 0; k < dim; k++) {
//               xg[k] += phiP[i] * faceCoordinates[k][i]; // xg(ig)= \sum_{i=0}^faceDofs phi[i](xig) facecoordinates[i]
//             }
//           }
//           double P0;
//           SetBoundaryCondition(xg, "P", P0, faceIndex, 0.);   // return P0
//
//           // *** phi_i loop ***
//           for(unsigned i = 0; i < faceDofs; i++) {
//             unsigned inode = msh->GetLocalFaceVertexIndex(iel, jface, i);
//             for(unsigned k = 0; k < dim; k++) {
//               aResP[inode] +=  -P0 * phiP[i] * normal[k] * weight;
//             }
//           }
//         }
//       }
//     }

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solDType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solDType]->Jacobian(coordX, ig, weight, phiD, gradPhiD);
      msh->_finiteElement[ielGeom][solPType]->Jacobian(coordX, ig, weight, phiP, gradPhiP);

      //phiP = msh->_finiteElement[ielGeom][solPType]->GetPhi(ig);

      vector < adept::adouble > solD_gss(dim, 0);
      vector < double > solDOld_gss(dim, 0);
      vector < vector < adept::adouble > > gradSolD_gss(dim);
      vector < vector < double > > gradSolDOld_gss(dim);

      vector < adept::adouble > gradSolP_gss(dim, 0);
      adept::adouble solP_gss = 0;




      for(unsigned  k = 0; k < dim; k++) {
        gradSolD_gss[k].assign(dim, 0.);
        gradSolDOld_gss[k].assign(dim, 0.);
      }

      for(unsigned i = 0; i < nDofsD; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          solD_gss[k] += solD[k][i] * phiD[i];
          solDOld_gss[k] += solDOld[k][i] * phiD[i];

        }
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            gradSolD_gss[k][j] += solD[k][i] * gradPhiD[i * dim + j];
            gradSolDOld_gss[k][j] += solDOld[k][i] * gradPhiD[i * dim + j];
          }
        }
      }


      for(unsigned i = 0; i < nDofsP; i++) {
        solP_gss += phiP[i] * solP[i];
      }


      for(unsigned i = 0; i < nDofsD; i++) {
        for(unsigned k = 0; k < dim; k++) {
          gradSolP_gss[k] += solP[i] * gradPhiP[i * dim + k];
        }
      }

      adept::adouble div = 0;
      double divOld = 0;
      for(unsigned k = 0; k < dim; k++) {
        div += gradSolD_gss[k][k];
        divOld += gradSolDOld_gss[k][k];
      }

      vector< vector < adept::adouble > > K = {{k0 * exp(beta * div), 0}, {0, k0 * exp(beta * div)}};


    // phiD_i loop
      for(unsigned i = 0; i < nDofsD; i++) {
        std::vector < adept::adouble > sigma(dim, 0.);

        for(unsigned  k = 0; k < dim; k++) {  
          for(unsigned j = 0; j < dim; j++) {  
            sigma[k] +=  mu * gradPhiD[i * dim + j] * (gradSolD_gss[k][j] + gradSolD_gss[j][k]); 

          }
          sigma[k] += (-solP_gss + lambda * div) * gradPhiD[i * dim + k]; 
        }
        for(unsigned  k = 0; k < dim; k++) {
          aResD[k][i] += - sigma[k] * weight;
        }
      } // end phiD_i loop

      
     //phiP_i loop
      for(unsigned i = 0; i < nDofsP; i++) {
        adept::adouble permBlock = 0.;
        for(unsigned k = 0; k < dim; k++) {
          for(unsigned j = 0; j < dim; j++) {
            permBlock += (1 / mu_w) * K[k][j] * gradSolP_gss[j] * gradPhiD[i * dim + k];
          }
        }
        aResP[i] += ((div - divOld) / dt * phiP[i] + permBlock) * weight;

      } // end phiP_i loop
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store them in RES
    Res.resize(nDofsAll);    //resize

    for(int i = 0; i < nDofsD; i++) {
      for(unsigned  k = 0; k < dim; k++) {
        Res[ k * nDofsD + i ] = -aResD[k][i].value();
      }
    }

    for(int i = 0; i < nDofsP; i++) {
      Res[ dim * nDofsD + i ] = -aResP[i].value();
    }

    RES->add_vector_blocked(Res, sysDof);

    //Extract and store the Jacobian

    Jac.resize(nDofsAll * nDofsAll);
    // define the dependent variables

    for(unsigned  k = 0; k < dim; k++) {
      s.dependent(&aResD[k][0], nDofsD);
    }
    s.dependent(&aResP[0], nDofsP);

    // define the independent variables
    for(unsigned  k = 0; k < dim; k++) {
      s.independent(&solD[k][0], nDofsD);
    }
    s.independent(&solP[0], nDofsP);

    // get the and store jacobian matrix (row-major)
    s.jacobian(&Jac[0], true);  // This is rowwise order.
    KK->add_matrix_blocked(Jac, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();
  KK->close();
  // ***************** END ASSEMBLY *******************
}





