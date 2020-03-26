#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include "MultiLevelSolution.hpp"

//THIS IS THE 2D ASSEMBLY FOR THE LOCAL FETI PROBLEM WITH 2 SUBDOMAINS

using namespace femus;

double kappa = 1;

const elem_type *fem = new const elem_type_2D ("quad", "linear", "second");   //to use a different quadrature rule in the inner integral

const elem_type *femQuadrature = new const elem_type_2D ("quad", "linear", "eighth");   //to use a different quadrature rule in the inner integral

void AssembleLocalSys (MultiLevelProblem& ml_prob) {
  adept::Stack& s = FemusInit::_adeptStack;

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("Local");
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel (level);
  elem*                     el = msh->el;

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel (level);

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];
  SparseMatrix*             KK = pdeSys->_KK;
  NumericVector*           RES = pdeSys->_RES;

  const unsigned  dim = msh->GetDimension();
  unsigned dim2 = (3 * (dim - 1) + ! (dim - 1));
  const unsigned maxSize = static_cast< unsigned > (ceil (pow (3, dim)));       // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the noumber of processes (for parallel computation)

  unsigned soluIndex = mlSol->GetIndex ("u");   // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType (soluIndex);   // get the finite element type for "u"

  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex ("u");   // get the position of "u" in the pdeSys object

  vector < adept::adouble >  solu; // local solution for the local assembly (it uses adept)
  solu.reserve (maxSize);

  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector < vector < double > > x1 (dim);

  for (unsigned k = 0; k < dim; k++) {
    x1[k].reserve (maxSize);
  }

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight

  phi.reserve (maxSize);
  phi_x.reserve (maxSize * dim);
  phi_xx.reserve (maxSize * dim2);

  vector< adept::adouble > aRes; // local redidual vector
  aRes.reserve (maxSize);

  vector< int > l2GMap1; // local to global mapping
  l2GMap1.reserve (maxSize);

  vector< double > Res1; // local redidual vector
  Res1.reserve (maxSize);

  vector < double > Jac11;
  Jac11.reserve (maxSize * maxSize);

  KK->zero(); // Set to zero all the entries of the Global Matrix

  //BEGIN local assembly

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nDofu  = msh->GetElementDofNumber (iel, soluType);   // number of solution element dofs
    unsigned nDofx = msh->GetElementDofNumber (iel, xType);   // number of coordinate element dofs

    // resize local arrays
    l2GMap1.resize (nDofu);
    solu.resize (nDofu);

    for (int i = 0; i < dim; i++) {
      x1[i].resize (nDofx);
    }

    aRes.resize (nDofu);   //resize
    std::fill (aRes.begin(), aRes.end(), 0);   //set aRes to zero

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofu; i++) {
      unsigned solDof = msh->GetSolutionDof (i, iel, soluType);   // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex]) (solDof);     // global extraction and local storage for the solution
      l2GMap1[i] = pdeSys->GetSystemDof (soluIndex, soluPdeIndex, i, iel);   // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof (i, iel, xType);   // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x1[jdim][i] = (*msh->_topology->_Sol[jdim]) (xDof);     // global extraction and local storage for the element coordinates
      }
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian (x1, ig, weight, phi, phi_x, boost::none);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point

      vector < adept::adouble > gradSolu_gss (dim, 0.);
      vector < double > x_gss (dim, 0.);

      for (unsigned i = 0; i < nDofu; i++) {
        for (unsigned jdim = 0; jdim < dim; jdim++) {
          gradSolu_gss[jdim] += phi_x[i * dim + jdim] * solu[i];
          x_gss[jdim] += x1[jdim][i] * phi[i];
        }
      }

      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofu; i++) {

        adept::adouble laplace = 0.;

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          laplace   +=  kappa * phi_x[i * dim + jdim] * gradSolu_gss[jdim];
        }
        double srcTerm =  - 1. ; // so f = 1
        aRes[i] += (srcTerm * phi[i] + laplace) * weight;

      } // end phi_i loop
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    Res1.resize (nDofu);   //resize

    for (int i = 0; i < nDofu; i++) {
      Res1[i] = - aRes[i].value();
    }

    RES->add_vector_blocked (Res1, l2GMap1);

    // define the dependent variables
    s.dependent (&aRes[0], nDofu);

    // define the independent variables
    s.independent (&solu[0], nDofu);

    // get the jacobian matrix (ordered by row major )
    Jac11.resize (nDofu * nDofu);   //resize
    s.jacobian (&Jac11[0], true);

    //store K in the global matrix KK
    KK->add_matrix_blocked (Jac11, l2GMap1, l2GMap1);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  //END local assembly

  RES->close();

  KK->close();

//     Mat A = ( static_cast<PetscMatrix*> ( KK ) )->mat();
//     MatAssemblyBegin ( A, MAT_FINAL_ASSEMBLY );
//     MatAssemblyEnd ( A, MAT_FINAL_ASSEMBLY );
//     PetscViewer viewer;
//     MatView ( A, viewer );

//     Vec v = ( static_cast< PetscVector* > ( RES ) )->vec();
//     VecView(v,PETSC_VIEWER_STDOUT_WORLD);

  // ***************** END ASSEMBLY *******************
}

void AssembleLocalSysFETI (MultiLevelProblem& ml_prob) {
  adept::Stack& s = FemusInit::_adeptStack;

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("Local_FETI");
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel (level);
  elem*                     el = msh->el;

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel (level);

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];
  SparseMatrix*             KK = pdeSys->_KK;
  NumericVector*           RES = pdeSys->_RES;

  const unsigned  dim = msh->GetDimension();
  unsigned dim2 = (3 * (dim - 1) + ! (dim - 1));
  const unsigned maxSize = static_cast< unsigned > (ceil (pow (3, dim)));       // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the noumber of processes (for parallel computation)

  unsigned solu1Index = mlSol->GetIndex ("u1");   // get the position of "u1" in the ml_sol object
  unsigned solu1Type = mlSol->GetSolutionType (solu1Index);   // get the finite element type for "u1"

  unsigned solu2Index = mlSol->GetIndex ("u2");   // get the position of "u2" in the ml_sol object
  unsigned solu2Type = mlSol->GetSolutionType (solu2Index);   // get the finite element type for "u2"

  unsigned solmuIndex = mlSol->GetIndex ("mu");   // get the position of "mu" in the ml_sol object
  unsigned solmuType = mlSol->GetSolutionType (solmuIndex);   // get the finite element type for "mu"

  unsigned u1FlagIndex = mlSol->GetIndex ("u1Flag");
  unsigned u1FlagType = mlSol->GetSolutionType (u1FlagIndex);

  unsigned u2FlagIndex = mlSol->GetIndex ("u2Flag");
  unsigned u2FlagType = mlSol->GetSolutionType (u2FlagIndex);

  unsigned muFlagIndex = mlSol->GetIndex ("muFlag");
  unsigned muFlagType = mlSol->GetSolutionType (muFlagIndex);

  unsigned solu1PdeIndex;
  solu1PdeIndex = mlPdeSys->GetSolPdeIndex ("u1");   // get the position of "u1" in the pdeSys object

  unsigned solu2PdeIndex;
  solu2PdeIndex = mlPdeSys->GetSolPdeIndex ("u2");   // get the position of "u2" in the pdeSys object

  unsigned solmuPdeIndex;
  solmuPdeIndex = mlPdeSys->GetSolPdeIndex ("mu");   // get the position of "mu" in the pdeSys object

  vector < double >  solu1;
  solu1.reserve (maxSize);

  vector < double >  solu2;
  solu2.reserve (maxSize);

  vector < double >  solmu;
  solmu.reserve (maxSize);

  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector < vector < double > > x1 (dim);

  for (unsigned k = 0; k < dim; k++) {
    x1[k].reserve (maxSize);
  }

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight

  phi.reserve (maxSize);
  phi_x.reserve (maxSize * dim);
  phi_xx.reserve (maxSize * dim2);

  vector< int > l2GMapu1;
  l2GMapu1.reserve (maxSize);

  vector< int > l2GMapu2;
  l2GMapu2.reserve (maxSize);

  vector< int > l2GMapmu; // local to global mapping for mu
  l2GMapmu.reserve (maxSize);

  vector< double > Resu1; // local redidual vector for u1
  Resu1.reserve (maxSize);

  vector< double > Resu2; // local redidual vector for u2
  Resu2.reserve (maxSize);

  vector< double > Resmu; // local redidual vector for mu
  Resmu.reserve (maxSize);

  vector < double > Jacu1;  // stiffness matrix for u1
  Jacu1.reserve (maxSize * maxSize);

  vector < double > Jacu2;  // stiffness matrix for u2
  Jacu2.reserve (maxSize * maxSize);

  vector < double > Jacmu;  // diag block for mu
  Jacmu.reserve (maxSize * maxSize);
  vector < double > M1;
  M1.reserve (maxSize * maxSize);
  vector < double > M2;
  M2.reserve (maxSize * maxSize);

  KK->zero(); // Set to zero all the entries of the Global Matrix

  //BEGIN local assembly

  //BEGIN creation of the flags for the assembly procedure

  //flag = 1 assemble
  //flag = 0 don't assemble

  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType (iel);
    short unsigned ielGroup = msh->GetElementGroup (iel);
    unsigned nDof  = msh->GetElementDofNumber (iel, solu1Type);
    unsigned nDofmu  = msh->GetElementDofNumber (iel, solmuType);

    double epsilon = 1.e-7;
    double rightBound =  epsilon;
    double leftBound =  - epsilon;

    std::vector < double > xCoords (nDof);

    for (unsigned i = 0; i < nDof; i++) {
      unsigned solDof  = msh->GetSolutionDof (i, iel, solu1Type);
      unsigned xDof  = msh->GetSolutionDof (i, iel, xType);
      xCoords[i] = (*msh->_topology->_Sol[0]) (xDof);

      if (xCoords[i] < rightBound) sol->_Sol[u1FlagIndex]->add (solDof, 1.);

      if (xCoords[i] > leftBound) sol->_Sol[u2FlagIndex]->add (solDof, 1.);

      if (xCoords[i] < rightBound && xCoords[i] > leftBound)  sol->_Sol[muFlagIndex]->add (solDof, 1.);

    }

  }

  sol->_Sol[u1FlagIndex]->close();
  sol->_Sol[u2FlagIndex]->close();
  sol->_Sol[muFlagIndex]->close();

  for (unsigned idof = msh->_dofOffset[solu1Type][iproc]; idof < msh->_dofOffset[solu1Type][iproc + 1]; idof++) {

    double u1Flag = (*sol->_Sol[u1FlagIndex]) (idof);
    if (u1Flag > 0) sol->_Sol[u1FlagIndex]->set (idof, 1.);
    else {
      sol->_Bdc[solu1Index]->set (idof, 0.);
      sol->_Sol[solu1Index]->set (idof, 0.);
    }

    double u2Flag = (*sol->_Sol[u2FlagIndex]) (idof);
    if (u2Flag > 0) sol->_Sol[u2FlagIndex]->set (idof, 1.);
    else {
      sol->_Bdc[solu2Index]->set (idof, 0.);
      sol->_Sol[solu2Index]->set (idof, 0.);
    }

    double muFlag = (*sol->_Sol[muFlagIndex]) (idof);
    if (muFlag > 0) sol->_Sol[muFlagIndex]->set (idof, 1.);
    else {
      sol->_Bdc[solmuIndex]->set (idof, 0.);
      sol->_Sol[solmuIndex]->set (idof, 0.);
    }

  }

  sol->_Sol[u1FlagIndex]->close();
  sol->_Sol[u2FlagIndex]->close();
  sol->_Sol[muFlagIndex]->close();

  sol->_Sol[solu1Index]->close();
  sol->_Sol[solu2Index]->close();
  sol->_Sol[solmuIndex]->close();

  sol->_Bdc[solu1Index]->close();
  sol->_Bdc[solu2Index]->close();
  sol->_Bdc[solmuIndex]->close();

  //END creation of the flags for the assembly procedure


  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    bool interfaceElem = false;

    short unsigned ielGeom = msh->GetElementType (iel);
    short unsigned ielGroup = msh->GetElementGroup (iel);
    unsigned nDof1  = msh->GetElementDofNumber (iel, solu1Type);

    l2GMapu1.resize (nDof1);
    l2GMapu2.resize (nDof1);
    l2GMapmu.resize (nDof1);

    solu1.resize (nDof1);
    solu2.resize (nDof1);
    solmu.resize (nDof1);

    Jacu1.assign (nDof1 * nDof1, 0.);
    Resu1.assign (nDof1, 0.);

    Jacu2.assign (nDof1 * nDof1, 0.);
    Resu2.assign (nDof1, 0.);

    Jacmu.assign (nDof1 * nDof1, 0.);
    M1.assign (nDof1 * nDof1, 0.);
    M2.assign (nDof1 * nDof1, 0.);
    Resmu.assign (nDof1, 0.);

    for (int k = 0; k < dim; k++) {
      x1[k].resize (nDof1);
    }

    for (unsigned i = 0; i < nDof1; i++) {
      l2GMapu1[i] = pdeSys->GetSystemDof (solu1Index, solu1PdeIndex, i, iel);
      l2GMapu2[i] = pdeSys->GetSystemDof (solu2Index, solu2PdeIndex, i, iel);
      l2GMapmu[i] = pdeSys->GetSystemDof (solmuIndex, solmuPdeIndex, i, iel);

      unsigned solDofu1 = msh->GetSolutionDof (i, iel, solu1Type);
      unsigned solDofu2 = msh->GetSolutionDof (i, iel, solu2Type);
      unsigned solDofmu = msh->GetSolutionDof (i, iel, solmuType);

      solu1[i] = (*sol->_Sol[solu1Index]) (solDofu1);
      solu2[i] = (*sol->_Sol[solu2Index]) (solDofu2);
      solmu[i] = (*sol->_Sol[solmuIndex]) (solDofmu);

      unsigned xDof  = msh->GetSolutionDof (i, iel, xType);
      for (unsigned k = 0; k < dim; k++) {
        x1[k][i] = (*msh->_topology->_Sol[k]) (xDof);
      }
      if (fabs (x1[0][i]) < 1.e-10) interfaceElem = true;

    }

    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solu1Type]->GetGaussPointNumber(); ig++) {

      msh->_finiteElement[ielGeom][solu1Type]->Jacobian (x1, ig, weight, phi, phi_x, boost::none);

      bool ielU1 = (ielGroup == 7) ? true : false;
      bool ielU2 = (ielGroup == 8) ? true : false;

//       vector < vector < double > > laplace (nDofu);
      double srcTerm =  1. ; // so f = 1
      for (unsigned i = 0; i < nDof1; i++) {
//         laplace[i].assign (nDofu, 0.);
        for (unsigned j = 0; j < nDof1; j++) {
          double laplace = 0.;
          for (unsigned kdim = 0; kdim < dim; kdim++) {
            laplace += (phi_x[i * dim + kdim] * phi_x[j * dim + kdim]) * weight;
          }
          if (ielU1) {
            Jacu1[i * nDof1 + j] -= laplace;
            Resu1[i] +=  laplace * solu1[j] - srcTerm * weight * phi[i];
          }
          if (ielU2) {
            Jacu2[i * nDof1 + j] -= laplace;
            Resu2[i] +=  laplace * solu2[j] - srcTerm * weight * phi[i];
          }

        }

        if (interfaceElem) {
          if (fabs (x1[0][i]) < 1.e-10) {
            //lumped mass matrix
            double Mlumped = phi[i] * weight;
            M1[ i * nDof1 + i ] +=  Mlumped;
            M2[ i * nDof1 + i ] += - Mlumped;
            Resu1[i] -= Mlumped * solmu[i];
            Resu2[i] -= - Mlumped * solmu[i];
            Resmu[i] -= Mlumped * (solu1[i] - solu2[i]);
          }
        }

      }

    } // end gauss point loop

    KK->add_matrix_blocked (Jacu1, l2GMapu1, l2GMapu1);
    RES->add_vector_blocked (Resu1, l2GMapu1);

    KK->add_matrix_blocked (Jacu2, l2GMapu2, l2GMapu2);
    RES->add_vector_blocked (Resu2, l2GMapu2);

    KK->add_matrix_blocked (Jacmu, l2GMapmu, l2GMapmu);

    if (interfaceElem) {
      KK->add_matrix_blocked (M1, l2GMapmu, l2GMapu1); //M1 (mu rows and u1 columns)
      KK->add_matrix_blocked (M2, l2GMapmu, l2GMapu2); //M2 (mu rows and u2 columns)
      KK->add_matrix_blocked (M1, l2GMapu1, l2GMapmu); //M1 transpose (u1 rows and mu columns)
      KK->add_matrix_blocked (M2, l2GMapu2, l2GMapmu); //M2 transpose (u2 rows and mu columns)
      RES->add_vector_blocked (Resmu, l2GMapmu);
    }


  } //end element loop for each process

  //END local assembly

  RES->close();

  KK->close();

//     Mat A = ( static_cast<PetscMatrix*> ( KK ) )->mat();
//     MatAssemblyBegin ( A, MAT_FINAL_ASSEMBLY );
//     MatAssemblyEnd ( A, MAT_FINAL_ASSEMBLY );
//     PetscViewer viewer;
//     MatView ( A, viewer );

//     Vec v = ( static_cast< PetscVector* > ( RES ) )->vec();
//     VecView(v,PETSC_VIEWER_STDOUT_WORLD);

//     PetscViewer    viewer;
//   PetscViewerDrawOpen (PETSC_COMM_WORLD, NULL, NULL, 0, 0, 900, 900, &viewer);
//   PetscObjectSetName ( (PetscObject) viewer, "Nonlocal FETI matrix");
//   PetscViewerPushFormat (viewer, PETSC_VIEWER_DRAW_LG);
//   MatView ( (static_cast<PetscMatrix*> (KK))->mat(), viewer);
//   double a;
//   std::cin >> a;
//   abort();

  // ***************** END ASSEMBLY *******************
}









