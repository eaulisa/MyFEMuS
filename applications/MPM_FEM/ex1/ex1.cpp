// tests the inclusion algorithm for all 2D and 3D elements

#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "XDMFWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "LinearImplicitSystem.hpp"
#include "adept.h"
#include "Marker.hpp"
#include "Line.hpp"
#include "TransientSystem.hpp"
#include "OprtrTypeEnum.hpp"

using namespace femus;


bool SetRefinementFlag(const std::vector < double >& x, const int& elemgroupnumber, const int& level) {

  bool refine = 0;

  if(elemgroupnumber == 6 && level < 4) refine = 1;
  if(elemgroupnumber == 7 && level < 5) refine = 1;
  if(elemgroupnumber == 8 && level < 6) refine = 1;

  return refine;

}

bool SetBoundaryCondition(const std::vector < double >& x, const char name[], double& value, const int facename, const double time) {
  bool test = 1; //dirichlet
  value = 0.;
  return test;
}

void AssembleMPMSys(MultiLevelProblem& ml_prob);
void GridToParticlesProjection(MultiLevelProblem& ml_prob);
void AssembleFEM(MultiLevelProblem& ml_prob);

Line* linea;

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;
  unsigned numberOfUniformLevels = 4; //for refinement in 3D
  //unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;


  /* element types
  0 = HEX
  1 = TET
  2 = WEDGE
  3 = QUAD
  4 = TRI
   */

  unsigned solType = 2;

  std::cout << " --------------------------------------------------     TEST     --------------------------------------------------" << std::endl;

  mlMsh.ReadCoarseMesh("./input/adaptiveRef4.neu", "fifth", scalingFactor);
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , SetRefinementFlag);

  unsigned dim = mlMsh.GetDimension();

//   PrintLine(DEFAULT_OUTPUTDIR, linea.GetLine(), false, 0);

  MultiLevelSolution mlSol(&mlMsh);
  // add variables to mlSol
  mlSol.AddSolution("DX", LAGRANGE, SECOND, 2);
  if(dim > 1) mlSol.AddSolution("DY", LAGRANGE, SECOND, 2);
  if(dim > 2) mlSol.AddSolution("DZ", LAGRANGE, SECOND, 2);

  mlSol.AddSolution("U", LAGRANGE, SECOND, 2);
  if(dim > 1) mlSol.AddSolution("V", LAGRANGE, SECOND, 2);
  if(dim > 2) mlSol.AddSolution("W", LAGRANGE, SECOND, 2);

  mlSol.AddSolution("AU", LAGRANGE, SECOND, 2);
  if(dim > 1) mlSol.AddSolution("AV", LAGRANGE, SECOND, 2);
  if(dim > 2) mlSol.AddSolution("AW", LAGRANGE, SECOND, 2);

  mlSol.Initialize("All");

  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

  // ******* Set boundary conditions *******
  mlSol.GenerateBdc("DX", "Steady");
  if(dim > 1) mlSol.GenerateBdc("DY", "Steady");
  if(dim > 2) mlSol.GenerateBdc("DZ", "Steady");

  MultiLevelProblem ml_prob(&mlSol);

  // ******* Add MPM system to the MultiLevel problem *******
  NonLinearImplicitSystem& system = ml_prob.add_system<NonLinearImplicitSystem> ("MPM_FEM");
  system.AddSolutionToSystemPDE("DX");
  if(dim > 1)system.AddSolutionToSystemPDE("DY");
  if(dim > 2) system.AddSolutionToSystemPDE("DZ");

  // ******* System MPM Assembly *******
  system.SetAssembleFunction(AssembleMPMSys);
  //system.SetAssembleFunction(AssembleFEM);
  // ******* set MG-Solver *******
  system.SetMgType(V_CYCLE);


  system.SetAbsoluteLinearConvergenceTolerance(1.0e-10);
  system.SetMaxNumberOfLinearIterations(10);
  system.SetNonLinearConvergenceTolerance(1.e-9);
  system.SetMaxNumberOfNonLinearIterations(20);


  system.SetNumberPreSmoothingStep(1);
  system.SetNumberPostSmoothingStep(1);

  // ******* Set Preconditioner *******

  system.SetMgSmoother(GMRES_SMOOTHER);

  system.init();

  // ******* Set Smoother *******
  //system.SetSolverFineGrids(RICHARDSON);
  //system.SetRichardsonScaleFactor(.5);
  system.SetSolverFineGrids(GMRES);

  system.SetPreconditionerFineGrids(ILU_PRECOND);

  system.SetTolerances(1.e-12, 1.e-20, 1.e+50, 20, 20);


  unsigned rows = 60;
  unsigned columns = 120;
  unsigned size = rows * columns;

  std::vector < std::vector < double > > x; // marker
  std::vector < MarkerType > markerType;

  x.resize(size);
  markerType.resize(size);

  std::vector < std::vector < std::vector < double > > > line(1);
  std::vector < std::vector < std::vector < double > > > line0(1);

  for(unsigned j = 0; j < size; j++) {
    x[j].assign(dim, 0.);
    markerType[j] = VOLUME;
  }

  // double pi = acos(-1.);

  //BEGIN initialization

  for(unsigned i = 0; i < rows; i++) {
    for(unsigned j = 0; j < columns; j++) {


      x[i * columns + j][0] = -0.5 + 0.0001 + ((0.625 + 0.0001) / (columns - 1)) * j;
      x[i * columns + j][1] = -0.0625 + 0.0001 + ((0.25 + 0.0001) / (rows - 1)) * i;
      if(dim == 3) {
        x[j][2] = 0.;
      }
    }
  }


//   for(unsigned i = 0; i < rows; i++) {
//     for(unsigned j = 0; j < columns; j++) {
//       x[i * columns + j][0] = -0.46875 + (0.0625) * j;
//       if(i == 1) {
//         x[i * columns + j][1] = 0.03125;
//       }
//       else {
//         x[i * columns + j][1] = -0.03125;
//       }
//     }
//   }

  //END

  linea = new Line(x, markerType, mlSol.GetLevel(numberOfUniformLevels - 1), solType);

  linea->GetLine(line0[0]);
  PrintLine(DEFAULT_OUTPUTDIR, line0, false, 0);

  system.MGsolve();

  GridToParticlesProjection(ml_prob);

  linea->UpdateLineMPM();

  linea->GetLine(line[0]);
  PrintLine(DEFAULT_OUTPUTDIR, line, false, 1);

  // ******* Print solution *******
  mlSol.SetWriter(VTK);

  std::vector<std::string> mov_vars;
  mov_vars.push_back("DX");
  mov_vars.push_back("DY");
  mov_vars.push_back("DZ");
  mlSol.GetWriter()->SetMovingMesh(mov_vars);

  std::vector<std::string> print_vars;
  print_vars.push_back("All");

  mlSol.GetWriter()->SetDebugOutput(true);
  mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars);


  delete linea;
  return 0;

} //end main



void AssembleMPMSys(MultiLevelProblem& ml_prob) {

  // ml_prob is the global object from/to where get/set all the data
  // level is the level of the PDE system to be assembled
  // levelMax is the Maximum level of the MultiLevelProblem
  // assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  clock_t AssemblyTime = 0;
  clock_t start_time, end_time;

  //pointers and references

  NonLinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<NonLinearImplicitSystem> ("MPM_FEM");
  const unsigned  level = my_nnlin_impl_sys.GetLevelToAssemble();
  MultiLevelSolution* ml_sol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* mysolution = ml_sol->GetSolutionLevel(level);     // pointer to the solution (level) object
  LinearEquationSolver* myLinEqSolver = my_nnlin_impl_sys._LinSolver[level];  // pointer to the equation (level) object

  Mesh* mymsh = ml_prob._ml_msh->GetLevel(level);     // pointer to the mesh (level) object
  elem* myel = mymsh->el;   // pointer to the elem object in msh (level)
  SparseMatrix* myKK = myLinEqSolver->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* myRES =  myLinEqSolver->_RES;  // pointer to the global residual vector object in pdeSys (level)
  bool assembleMatrix = my_nnlin_impl_sys.GetAssembleMatrix();

// call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;
  if(assembleMatrix) s.continue_recording();
  else s.pause_recording();

  const unsigned dim = mymsh->GetDimension();

  // reserve memory for the local standar vectors
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  // data
  unsigned nel    = mymsh->GetNumberOfElements();
  unsigned iproc  = mymsh->processor_id();

  // local objects
  vector<adept::adouble> SolDp(dim);
  vector<adept::adouble> SolDgss(dim);
  vector<double> SolVp(dim);
  vector<double> SolAp(dim);
  vector < vector < double > > Fp(dim);

  vector<vector < adept::adouble > > GradSolDp(dim);

  for(int i = 0; i < dim; i++) {
    GradSolDp[i].resize(dim);
  }

  vector < double > phi;
  vector < double > phi_hat;
  vector < adept::adouble> gradphi;
  vector < double > gradphi_hat;
  vector < adept::adouble> nablaphi;
  vector < double > nablaphi_hat;

  phi.reserve(max_size);
  phi_hat.reserve(max_size);

  gradphi.reserve(max_size * dim);
  gradphi_hat.reserve(max_size * dim);

  vector <vector < adept::adouble> > vx(dim); //vx is coordX in assembly of ex30
  vector <vector < double> > vx_hat(dim);

  vector< vector< adept::adouble > > SolDd(dim);      // local solution
  vector< vector< double > > SolVd(dim);      // local solution
  vector< vector< double > > SolAd(dim);      // local solution

  vector< vector< int > > dofsVAR(dim);

  vector< vector< double > > Rhs(dim);     // local redidual vector
  vector< vector< adept::adouble > > aRhs(dim);     // local redidual vector

  vector < int > dofsAll;

  vector < double > Jac;

  adept::adouble weight;
  double weight_hat = 0.;
  adept::adouble weightFake; //WARNING remove after testing

  // gravity
  double gravity[3] = {0., -9.81 * 3, 0.}; //TODO use the actual value for gravity
  std::vector <double> gravityP(dim);

  double E = 1000000;
  double nu = 0.3;

  double K = E / (3.*(1. - 2.*nu));
  double lambda = E * nu / ((1. + nu) * (1. - 2 * nu));
  double mu = 1.5 * (K - lambda);

  //variable-name handling
  const char varname[9][3] = {"DX", "DY", "DZ", "U", "V", "W", "AU", "AV", "AW"}; //TODO is there a reason why varname[9][3] and not just varname[9] ?
  vector <unsigned> indexSolD(dim);
  vector <unsigned> indexSolV(dim);
  vector <unsigned> indexSolA(dim);
  vector <unsigned> indexPdeD(dim);
  unsigned solType = ml_sol->GetSolutionType(&varname[0][0]);

  for(unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolD[ivar] = ml_sol->GetIndex(&varname[ivar][0]);
    indexPdeD[ivar] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar][0]);
    indexSolV[ivar] = ml_sol->GetIndex(&varname[3 + ivar][0]);
    indexSolA[ivar] = ml_sol->GetIndex(&varname[6 + ivar][0]);
  }

  start_time = clock();

  if(assembleMatrix) myKK->zero();

  //line instances
  std::vector<unsigned> markerOffset = linea->GetMarkerOffset();
  unsigned markerOffset1 = markerOffset[iproc];
  unsigned markerOffset2 = markerOffset[iproc + 1];
  std::vector<Marker*> particles = linea->GetParticles();
  std::map<unsigned, std::vector < std::vector < std::vector < std::vector < double > > > > > aX;


  //BEGIN loop on elements (to initialize the "soft" stiffness matrix)
  for(int iel = mymsh->_elementOffset[iproc]; iel < mymsh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielt = mymsh->GetElementType(iel);

    unsigned nDofsV = mymsh->GetElementDofNumber(iel, solType);    // number of solution element dofs
    // unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);    // number of solution element dofs
    unsigned nDofsX = mymsh->GetElementDofNumber(iel, 2);    // number of coordinate element dofs

    unsigned nDofsVP = dim * nDofsV ;//+ nDofsP;
    // resize local arrays
    std::vector <int> sysDof(nDofsVP);

    for(unsigned  k = 0; k < dim; k++) {
      SolDd[k].resize(nDofsV);
      vx[k].resize(nDofsX);
      vx_hat[k].resize(nDofsX);
    }

    for(unsigned  k = 0; k < dim; k++) {
      aRhs[k].resize(nDofsV);    //resize
      std::fill(aRhs[k].begin(), aRhs[k].end(), 0);    //set aRes to zero
    }

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = mymsh->GetSolutionDof(i, iel, solType);    // global to global mapping between solution node and solution dof

      for(unsigned  k = 0; k < dim; k++) {
        SolDd[k][i] = (*mysolution->_Sol[indexSolD[k]])(solVDof);      // global extraction and local storage for the solution
        sysDof[i + k * nDofsV] = myLinEqSolver->GetSystemDof(indexSolD[k], indexPdeD[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = mymsh->GetSolutionDof(i, iel, 2);   // global to global mapping between coordinates node and coordinate dof

      for(unsigned k = 0; k < dim; k++) {
        vx_hat[k][i] = (*mymsh->_topology->_Sol[k])(coordXDof);
        vx[k][i] = vx_hat[k][i] /*+ SolDd[k][i]*/;      // TODO should we add also SolDd ???

      }
    }

    // start a new recording of all the operations involving adept::adouble variables
    if(assembleMatrix) s.new_recording();

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < mymsh->_finiteElement[ielt][solType]->GetGaussPointNumber(); ig++) {

      mymsh->_finiteElement[ielt][solType]->Jacobian(vx, ig, weight, phi, gradphi, nablaphi);
      mymsh->_finiteElement[ielt][solType]->Jacobian(vx_hat, ig, weight_hat, phi_hat, gradphi_hat, nablaphi_hat);


      vector < adept::adouble > SolDgss(dim, 0);
      vector < vector < adept::adouble > > GradSolDgss(dim);
      vector < vector < adept::adouble > > GradSolDhatgss(dim);

      for(unsigned  k = 0; k < dim; k++) {
        GradSolDgss[k].resize(dim);
        std::fill(GradSolDgss[k].begin(), GradSolDgss[k].end(), 0);

        GradSolDhatgss[k].resize(dim);
        std::fill(GradSolDhatgss[k].begin(), GradSolDhatgss[k].end(), 0);
      }

      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          SolDgss[k] += phi[i] * SolDd[k][i];
        }

        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            GradSolDgss[k][j] += gradphi[i * dim + j] * SolDd[k][i];
            GradSolDhatgss[k][j] += gradphi_hat[i * dim + j] * SolDd[k][i];
          }
        }
      }

      //BEGIN computation of the Cauchy Stress

      adept::adouble Cauchy[dim][dim];
      std::vector < std::vector <adept::adouble> > Fp(dim);
      adept::adouble J_hat; //det of the deformation gradient
      adept::adouble B[dim][dim]; //left Cauchy-Green deformation tensor
      adept::adouble devB[dim][dim]; //deviatoric part of B
      adept::adouble C[dim][dim]; //right Cauchy-Green deformation tensor
      std::vector < std::vector <double> > I(dim);
      adept::adouble I_1 = 0.; //trace of the right Cauchy-Green deformation tensor
      adept::adouble shear = 6. * 1.e5 ; //shear modulus in Pascal for an hyperelastic material
      adept::adouble bulk = 6. * 1.e6 ; //bulk modulus in Pascal for an hyperelastic material
      adept::adouble C1 = shear * 0.5;
      adept::adouble D1 = bulk * 0.5;


      for(int i = 0; i < dim; i++) {

        I[i].resize(dim);
        Fp[i].resize(dim);

        for(int j = 0; j < dim; j++) {

          if(i == j) {
            I[i][i] = 1.;
            Fp[i][i] = 1.;
          }
          else {
            I[i][j] = 0.;
            Fp[i][j] = 0.;
          }

        }
      }

      for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
          Fp[i][j] += GradSolDhatgss[i][j];
        }
      }

//       std::cout << "-------------------- DEFORMATION GRADIENT --------------------" << std::endl;
//       for(unsigned i = 0 ; i < dim; i++) {
//         for(unsigned j = 0 ; j < dim; j++) {
//           std::cout << Fp[i][j] << " " ;
//         }
//         std::cout << std::endl;
//       }

      for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
          B[i][j] = 0.;
          C[i][j] = 0.;

          for(int k = 0; k < dim; k++) {
            B[i][j] += Fp[i][k] * Fp[j][k];
            C[i][j] += Fp[k][i] * Fp[k][j];
          }
        }
      }

      if(dim == 1) {
        J_hat = Fp[0][0];
      }

      else if(dim == 2) {
        J_hat = (Fp[0][0] * Fp[1][1]) - (Fp[0][1] * Fp[1][0]);
      }

      else {
        J_hat =   Fp[0][0] * Fp[1][1] * Fp[2][2] + Fp[0][1] * Fp[1][2] * Fp[2][0] + Fp[0][2] * Fp[1][0] * Fp[2][1]
                  - Fp[2][0] * Fp[1][1] * Fp[0][2] - Fp[2][1] * Fp[1][2] * Fp[0][0] - Fp[2][2] * Fp[1][0] * Fp[0][1];
      }


      for(unsigned i = 0 ; i < dim; i++) {
        I_1 += C[i][i];
      }

      for(unsigned i = 0 ; i < dim; i++) {
        for(unsigned j = 0 ; j < dim; j++) {
          devB[i][j] = B[i][j]  - (1. / 3.) * I_1 * I[i][j];
          Cauchy[i][j] = (1 / J_hat) * (2. * D1 * J_hat * (J_hat - 1.) * I[i][j] + 2. * C1 * devB[i][j] / (pow(J_hat , (2. / 3.))));
        }
      }

//       std::cout << "-------------------- CAUCHY STRESS TENSOR --------------------" << std::endl;
//       for(unsigned i = 0 ; i < dim; i++) {
//         for(unsigned j = 0 ; j < dim; j++) {
//           std::cout << Cauchy[i][j] << " " ;
//         }
//         std::cout << std::endl;
//       }

      //END computation of the Cauchy Stress


      for(unsigned i = 0; i < nDofsV; i++) {
        vector < adept::adouble > softStiffness(dim, 0.);

        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            softStiffness[k]   +=  0.0001 * mu * gradphi[i * dim + j] * (GradSolDgss[k][j] + GradSolDgss[j][k]);
            //softStiffness[k]   +=  0.0001 * Cauchy[k][j] * gradphi[i * dim + j];
          }
        }


        for(unsigned  k = 0; k < dim; k++) {
          aRhs[k][i] += - (softStiffness[k])  * weight;
        }
      }

    } // end gauss point loop


    //copy the value of the adept::adoube aRes in double Res and store them in RES
    std::vector<double> Rhs(nDofsVP);    //resize

    for(int i = 0; i < nDofsV; i++) {
      for(unsigned  k = 0; k < dim; k++) {
        Rhs[ i +  k * nDofsV ] = -aRhs[k][i].value();
      }
    }

    myRES->add_vector_blocked(Rhs, sysDof);

    if(assembleMatrix) {
      Jac.resize(nDofsVP * nDofsVP);
      // define the dependent variables

      for(unsigned  k = 0; k < dim; k++) {
        s.dependent(&aRhs[k][0], nDofsV);
      }

      // define the independent variables
      for(unsigned  k = 0; k < dim; k++) {
        s.independent(&SolDd[k][0], nDofsV);
      }

      // get the and store jacobian matrix (row-major)
      s.jacobian(&Jac[0] , true);
      myKK->add_matrix_blocked(Jac, sysDof, sysDof);

      s.clear_independents();
      s.clear_dependents();
    }
  }
  //END building "soft" stiffness matrix


  //initialization of iel
  unsigned ielOld = UINT_MAX;


  //BEGIN loop on particles (used as Gauss points)
  for(unsigned iMarker = markerOffset1; iMarker < markerOffset2; iMarker++) {

    //element of particle iMarker
    unsigned iel = particles[iMarker]->GetMarkerElement();
    if(iel != UINT_MAX) {
      short unsigned ielt;
      unsigned nve;
      //update element related quantities only if we are in a different element
      if(iel != ielOld) {

        ielt = mymsh->GetElementType(iel);
        nve = mymsh->GetElementDofNumber(iel, solType);
        //initialization of everything is in common fluid and solid

        //Rhs
        for(int i = 0; i < dim; i++) {
          dofsVAR[i].resize(nve);
          SolDd[indexPdeD[i]].resize(nve);
          SolVd[i].resize(nve);
          SolAd[i].resize(nve);
          aRhs[indexPdeD[i]].resize(nve);
        }

        dofsAll.resize(0);


        for(int i = 0; i < dim; i++) {
          vx[i].resize(nve);
          vx_hat[i].resize(nve);
        }


        //BEGIN copy of the value of Sol at the dofs idof of the element iel
        for(unsigned i = 0; i < nve; i++) {
          unsigned idof = mymsh->GetSolutionDof(i, iel, solType); //local 2 global solution

          for(int j = 0; j < dim; j++) {
            SolDd[indexPdeD[j]][i] = (*mysolution->_Sol[indexSolD[j]])(idof);

            dofsVAR[j][i] = myLinEqSolver->GetSystemDof(indexSolD[j], indexPdeD[j], i, iel); //local 2 global Pde
            aRhs[indexPdeD[j]][i] = 0.;

            //Fixed coordinates (Reference frame)
            vx_hat[j][i] = (*mymsh->_topology->_Sol[j])(idof);

            SolVd[j][i] = (*mysolution->_Sol[indexSolV[j]])(idof);
            SolAd[j][i] = (*mysolution->_Sol[indexSolA[j]])(idof);
          }
        }
        //END

        // build dof composition
        for(int idim = 0; idim < dim; idim++) {
          dofsAll.insert(dofsAll.end(), dofsVAR[idim].begin(), dofsVAR[idim].end());
        }

        if(assembleMatrix) s.new_recording();

        for(unsigned idim = 0; idim < dim; idim++) {
          for(int j = 0; j < nve; j++) {
            vx[idim][j]    = vx_hat[idim][j] + SolDd[indexPdeD[idim]][j];
          }
        }
        // start a new recording of all the operations involving adept::adouble variables
        if(assembleMatrix) s.new_recording();
      }

      bool elementUpdate = (aX.find(iel) != aX.end()) ? false : true;  //TODO to be removed after we include FindLocalCoordinates in the advection
      particles[iMarker]->FindLocalCoordinates(solType, aX[iel], elementUpdate, mysolution, 0);

      // the local coordinates of the particles are the Gauss points in this context
      std::vector <double> xi = particles[iMarker]->GetMarkerLocalCoordinates();

      // the mass of the particle acts as weight
      //double mass = 0.217013888889; //TODO use this with the Cauchy stress formulation
      //double density = 10000.;
      double mass = 0.0217013888889;
      double density = 1000.;

      mymsh->_finiteElement[ielt][solType]->Jacobian(vx, xi, weightFake, phi, gradphi, nablaphi); //function to evaluate at the particles

      // displacement and velocity
      //BEGIN evaluates SolDp at the particle iMarker
      for(int i = 0; i < dim; i++) {
        SolDp[i] = 0.;
        gravityP[i] = 0.;

        for(int j = 0; j < dim; j++) {
          GradSolDp[i][j] = 0.;
        }

        for(unsigned inode = 0; inode < nve; inode++) {

          SolDp[i] += phi[inode] * SolDd[indexPdeD[i]][inode];
          gravityP[i] += phi[inode] * gravity[i];   //TODO added this, the value of gravity at the particle

          for(int j = 0; j < dim; j++) {
            GradSolDp[i][j] +=  gradphi[inode * dim + j] * SolDd[indexPdeD[i]][inode];
          }
        }
      }

      //END evaluates SolDp at the particle iMarker


      //BEGIN computation of the Cauchy Stress

      adept::adouble Cauchy[dim][dim];
      adept::adouble J_hat; //det of the deformation gradient
      adept::adouble B[dim][dim]; //left Cauchy-Green deformation tensor
      adept::adouble devB[dim][dim]; //deviatoric part of B
      adept::adouble C[dim][dim]; //right Cauchy-Green deformation tensor
      std::vector < std::vector <double> > I(dim);
      adept::adouble I_1 = 0.; //trace of the right Cauchy-Green deformation tensor
      adept::adouble shear = 6. * 1.e5 ; //shear modulus in Pascal for an hyperelastic material
      adept::adouble bulk = 6. * 1.e6 ; //bulk modulus in Pascal for an hyperelastic material
      adept::adouble C1 = shear * 0.5;
      adept::adouble D1 = bulk * 0.5;

      Fp = particles[iMarker]->GetDeformationGradient(); //extraction of the deformation gradient

      if(dim == 1) {
        J_hat = Fp[0][0];
      }

      else if(dim == 2) {
        J_hat = (Fp[0][0] * Fp[1][1]) - (Fp[0][1] * Fp[1][0]);
      }

      else {
        J_hat =   Fp[0][0] * Fp[1][1] * Fp[2][2] + Fp[0][1] * Fp[1][2] * Fp[2][0] + Fp[0][2] * Fp[1][0] * Fp[2][1]
                  - Fp[2][0] * Fp[1][1] * Fp[0][2] - Fp[2][1] * Fp[1][2] * Fp[0][0] - Fp[2][2] * Fp[1][0] * Fp[0][1];
      }

      for(int i = 0; i < dim; i++) {

        I[i].resize(dim);

        for(int j = 0; j < dim; j++) {

          if(i == j) {
            I[i][i] = 1.;
          }
          else {
            I[i][j] = 0.;
          }

          B[i][j] = 0.;
          C[i][j] = 0.;

          for(int k = 0; k < dim; k++) {
            B[i][j]     += Fp[i][k] * Fp[j][k];
            C[i][j]     += Fp[k][i] * Fp[k][j];
          }
        }
      }

      for(unsigned i = 0 ; i < dim; i++) {
        I_1 += C[i][i];
      }

      for(unsigned i = 0 ; i < dim; i++) {
        for(unsigned j = 0 ; j < dim; j++) {
          devB[i][j] = B[i][j]  - (1. / 3.) * I_1 * I[i][j];
          Cauchy[i][j] = (1 / J_hat) * (2. * D1 * J_hat * (J_hat - 1.) * I[i][j] + 2. * C1 * devB[i][j] / (pow(J_hat , (2. / 3.))));
        }
      }

//       std::cout << "--------------------  Cauchy Stress Tensor --------------------" << std::endl;
//       for(unsigned i = 0 ; i < dim; i++) {
//         for(unsigned j = 0 ; j < dim; j++) {
//           std::cout << Cauchy[i][j] << " " ;
//         }
//         std::cout << std::endl;
//       }


      //END computation of the Cauchy Stress


      //BEGIN computation of local residual
      adept::adouble divergence = 0.;
      for(unsigned i = 0; i < dim; i++) {
        divergence += GradSolDp[i][i];
      }

      for(unsigned k = 0; k < nve; k++) {
        for(unsigned i = 0; i < dim; i++) {
          adept::adouble weakLaplace = 0.;
          for(unsigned j = 0; j < dim; j++) {
            weakLaplace += 0.5 * (GradSolDp[i][j] + GradSolDp[j][i]) * gradphi[k * dim + j] ;
            //weakLaplace += Cauchy[i][j] * gradphi[k * dim + j] ;
          }
          aRhs[indexPdeD[i]][k] += - ((2. * mu * weakLaplace + lambda * divergence * gradphi[k * dim + i]) / density - gravityP[i] * phi[k]) * mass;
          //aRhs[indexPdeD[i]][k] += - (weakLaplace / density - gravityP[i] * phi[k]) * mass;
        }
      }
      //END

      if(iMarker == markerOffset2 - 1 || iel != particles[iMarker + 1]->GetMarkerElement()) {

        //copy adouble aRhs into double Rhs
        for(unsigned i = 0; i < dim; i++) {
          Rhs[indexPdeD[i]].resize(nve);

          for(int j = 0; j < nve; j++) {
            Rhs[indexPdeD[i]][j] = -aRhs[indexPdeD[i]][j].value();
          }
        }

        for(int i = 0; i < dim; i++) {
          myRES->add_vector_blocked(Rhs[indexPdeD[i]], dofsVAR[i]);
        }

        if(assembleMatrix) {
          //Store equations
          for(int i = 0; i < dim; i++) {
            s.dependent(&aRhs[indexPdeD[i]][0], nve);
            s.independent(&SolDd[indexPdeD[i]][0], nve);
          }

          Jac.resize((dim * nve) * (dim * nve));

          s.jacobian(&Jac[0], true);

          myKK->add_matrix_blocked(Jac, dofsAll, dofsAll);
          s.clear_independents();
          s.clear_dependents();
        }
      }
      //END local to global assembly

      ielOld = iel;
    }
    else {
      break;
    }
  }
  //END loop on particles

  myRES->close();

  if(assembleMatrix) {
    myKK->close();
  }

  // *************************************
  end_time = clock();
  AssemblyTime += (end_time - start_time);
  // ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

}


void GridToParticlesProjection(MultiLevelProblem& ml_prob) {

  // ml_prob is the global object from/to where get/set all the data
  // level is the level of the PDE system to be assembled
  // levelMax is the Maximum level of the MultiLevelProblem
  // assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  clock_t AssemblyTime = 0;
  clock_t start_time, end_time;

  //pointers and references

  NonLinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<NonLinearImplicitSystem> ("MPM_FEM");
  const unsigned  level = my_nnlin_impl_sys.GetLevelToAssemble();
  MultiLevelSolution* ml_sol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* mysolution = ml_sol->GetSolutionLevel(level);     // pointer to the solution (level) object

  Mesh* mymsh = ml_prob._ml_msh->GetLevel(level);     // pointer to the mesh (level) object
  elem* myel = mymsh->el;   // pointer to the elem object in msh (level)

  const unsigned dim = mymsh->GetDimension();

  // data
  unsigned nel    = mymsh->GetNumberOfElements();
  unsigned iproc  = mymsh->processor_id();

  // local objects

  vector< vector < double > > SolDd(dim);
  vector< vector < double > > SolAd(dim);
  vector< vector < double > > GradSolDp(dim);
  vector < vector < double > > Fp(dim); //deformation vector
  vector < vector < double > > FpOld;

  for(int i = 0; i < dim; i++) {
    GradSolDp[i].resize(dim);
    Fp[i].resize(dim);
  }

  vector < double > phi;
  vector < double > gradphi;
  vector < double > nablaphi;

  vector <vector < double> > vx(dim); //vx is coordX in assembly of ex30

  double weight;

  //variable-name handling
  const char varname[6][3] = {"DX", "DY", "DZ", "AU", "AV", "AW"};
  vector <unsigned> indexSolD(dim);
  vector <unsigned> indexSolA(dim);
  unsigned solType = ml_sol->GetSolutionType(&varname[0][0]);

  for(unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolD[ivar] = ml_sol->GetIndex(&varname[ivar][0]);
    indexSolA[ivar] = ml_sol->GetIndex(&varname[ivar][0]);
  }

  //line instances
  std::vector<unsigned> markerOffset = linea->GetMarkerOffset();
  unsigned markerOffset1 = markerOffset[iproc];
  unsigned markerOffset2 = markerOffset[iproc + 1];
  std::vector<Marker*> particles = linea->GetParticles();
  std::map<unsigned, std::vector < std::vector < std::vector < std::vector < double > > > > > aX;

  //initialization of iel
  unsigned ielOld = UINT_MAX;
  //declaration of element instances

  //BEGIN loop on particles (used as Gauss points)
  for(unsigned iMarker = markerOffset1; iMarker < markerOffset2; iMarker++) {

    //element of particle iMarker
    unsigned iel = particles[iMarker]->GetMarkerElement();
    if(iel != UINT_MAX) {

      short unsigned ielt;
      unsigned nve;
      std::vector <double> particleDisp(dim, 0.);
      std::vector <double> particleAcc(dim, 0.);

      //update element related quantities only if we are in a different element
      if(iel != ielOld) {

        ielt = mymsh->GetElementType(iel);
        nve = mymsh->GetElementDofNumber(iel, solType);
        //initialization of everything is in common fluid and solid

        for(int i = 0; i < dim; i++) {
          SolDd[i].resize(nve);
          SolAd[i].resize(nve);
        }

        // ----------------------------------------------------------------------------------------
        // coordinates, solutions, displacement, velocity dofs

        for(int i = 0; i < dim; i++) {
          vx[i].resize(nve);
        }

        //BEGIN copy of the value of Sol at the dofs idof of the element iel
        for(unsigned inode = 0; inode < nve; inode++) {
          unsigned idof = mymsh->GetSolutionDof(inode, iel, solType); //local 2 global solution

          for(int i = 0; i < dim; i++) {
            SolDd[i][inode] = (*mysolution->_Sol[indexSolD[i]])(idof);
            SolAd[i][inode] = (*mysolution->_Sol[indexSolA[i]])(idof);

            //Fixed coordinates (Reference frame)
            vx[i][inode] = (*mymsh->_topology->_Sol[i])(idof) + SolDd[i][inode];
          }
        }
        //END

      }

      bool elementUpdate = (aX.find(iel) != aX.end()) ? false : true;  //TODO to be removed after we include FindLocalCoordinates in the advection
      particles[iMarker]->FindLocalCoordinates(solType, aX[iel], elementUpdate, mysolution, 0);
      std::vector <double> xi = particles[iMarker]->GetMarkerLocalCoordinates();

      mymsh->_finiteElement[ielt][solType]->Jacobian(vx, xi, weight, phi, gradphi, nablaphi); //function to evaluate at the particles

      //update displacement and acceleration
      for(int i = 0; i < dim; i++) {
        for(unsigned inode = 0; inode < nve; inode++) {
          particleDisp[i] += phi[inode] * SolDd[i][inode];
          particleAcc[i] += phi[inode] * SolAd[i][inode];
        }
      }

      particles[iMarker]->SetMarkerDisplacement(particleDisp);
      particles[iMarker]->SetMarkerAcceleration(particleAcc);

      //movement of the particles
      particles[iMarker]->UpdateParticleCoordinates();

      particles[iMarker]->GetElementSerial(iel, mysolution, 0.);
      particles[iMarker]->SetIprocMarkerPreviousElement(iel);

      // update the deformation gradient
      //BEGIN evaluates SolDp at the particle iMarker
      for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
          GradSolDp[i][j] = 0.;
        }
        for(unsigned inode = 0; inode < nve; inode++) {
          for(int j = 0; j < dim; j++) {
            GradSolDp[i][j] +=  gradphi[inode * dim + j] * SolDd[i][inode];
          }
        }
      }

      FpOld = particles[iMarker]->GetDeformationGradient();

      for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
          Fp[i][j] = (GradSolDp[i][j] + 1.) * FpOld[i][j];
        }
      }

      particles[iMarker]->SetDeformationGradient(Fp);

      ielOld = iel;
    }
    else {
      break;
    }
  }
  //END loop on particles


}




void AssembleFEM(MultiLevelProblem & ml_prob) {

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<NonLinearImplicitSystem> ("MPM_FEM");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  //solution variable
  vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("DX");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("DY");    // get the position of "V" in the ml_sol object
  if(dim == 3) solVIndex[2] = mlSol->GetIndex("DZ");       // get the position of "V" in the ml_sol object

  unsigned solVType = mlSol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

  vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("DX");    // get the position of "U" in the pdeSys object
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("DY");    // get the position of "V" in the pdeSys object

  if(dim == 3) solVPdeIndex[2] = mlPdeSys->GetSolPdeIndex("DZ");

  vector < vector < adept::adouble > >  solV(dim);    // local solution

  vector< vector < adept::adouble > > aResV(dim);    // local redidual vector

  vector < vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for(unsigned  k = 0; k < dim; k++) {
    solV[k].reserve(maxSize);
    aResV[k].reserve(maxSize);
    coordX[k].reserve(maxSize);
  }

  vector <double> phiV;  // local test function
  vector <double> phiV_x; // local test function first order partial derivatives
  vector <double> phiV_xx; // local test function second order partial derivatives

  phiV.reserve(maxSize);
  phiV_x.reserve(maxSize * dim);
  phiV_xx.reserve(maxSize * dim2);

  double weight; // gauss point weight

  vector< int > sysDof; // local to global pdeSys dofs
  sysDof.reserve((dim + 1) *maxSize);

  vector< double > Res; // local redidual vector
  Res.reserve((dim + 1) *maxSize);

  vector < double > Jac;
  Jac.reserve((dim + 1) *maxSize * (dim + 1) *maxSize);

  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);

    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
    // unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);    // number of solution element dofs
    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType);    // number of coordinate element dofs

    unsigned nDofsVP = dim * nDofsV ;//+ nDofsP;
    // resize local arrays
    sysDof.resize(nDofsVP);

    for(unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofsV);
      coordX[k].resize(nDofsX);
    }

    for(unsigned  k = 0; k < dim; k++) {
      aResV[k].resize(nDofsV);    //resize
      std::fill(aResV[k].begin(), aResV[k].end(), 0);    //set aRes to zero
    }

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // global to global mapping between solution node and solution dof

      for(unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);      // global extraction and local storage for the solution
        sysDof[i + k * nDofsV] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // global to global mapping between coordinates node and coordinate dof

      for(unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solVType]->GetGaussPointNumber(); ig++) {

      msh->_finiteElement[ielGeom][solVType]->Jacobian(coordX, ig, weight, phiV, phiV_x, phiV_xx);

      vector < adept::adouble > solV_gss(dim, 0);
      vector < vector < adept::adouble > > gradSolV_gss(dim);

      for(unsigned  k = 0; k < dim; k++) {
        gradSolV_gss[k].resize(dim);
        std::fill(gradSolV_gss[k].begin(), gradSolV_gss[k].end(), 0);
      }

      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          solV_gss[k] += phiV[i] * solV[k][i];
        }

        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            gradSolV_gss[k][j] += phiV_x[i * dim + j] * solV[k][i];
          }
        }
      }


      double x = coordX[0][nDofsV - 1];
      double y = coordX[1][nDofsV - 1];

      double nu = (x < 0.0625 * 2 && y < 0.0625 && y > -0.0625) ?  1. : 0;

      double gravity[3] = {0, 0, 0};

      gravity[1] = (x < 0.0625 * 2 && y < 0.0625 && y > -0.0625) ?  -9.81 : 0.0;

      // *** phiV_i loop ***
      for(unsigned i = 0; i < nDofsV; i++) {
        vector < adept::adouble > NSV(dim, 0.);

        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            NSV[k]   +=  nu * phiV_x[i * dim + j] * (gradSolV_gss[k][j] + gradSolV_gss[j][k]);
          }
        }


        for(unsigned  k = 0; k < dim; k++) {
          aResV[k][i] += - (NSV[k] - gravity[k] * phiV[i]) * weight;
        }
      }

    } // end gauss point loop


    //copy the value of the adept::adoube aRes in double Res and store them in RES
    Res.resize(nDofsVP);    //resize

    for(int i = 0; i < nDofsV; i++) {
      for(unsigned  k = 0; k < dim; k++) {
        Res[ i +  k * nDofsV ] = -aResV[k][i].value();
      }
    }

    RES->add_vector_blocked(Res, sysDof);

    //Extarct and store the Jacobian

    Jac.resize(nDofsVP * nDofsVP);
    // define the dependent variables

    for(unsigned  k = 0; k < dim; k++) {
      s.dependent(&aResV[k][0], nDofsV);
    }

    // define the independent variables
    for(unsigned  k = 0; k < dim; k++) {
      s.independent(&solV[k][0], nDofsV);
    }

    // get the and store jacobian matrix (row-major)
    s.jacobian(&Jac[0] , true);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();

  KK->close();

  // ***************** END ASSEMBLY *******************
}





