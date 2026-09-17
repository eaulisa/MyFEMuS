#pragma once

void AssembleNormalFromField(MultiLevelProblem& ml_prob, const char* psiName);
void AssembleNormalLumped(MultiLevelProblem& ml_prob);

void AssembleCurvature(MultiLevelProblem& ml_prob);
void AssembleCurvatureLumped(MultiLevelProblem& ml_prob);

void AssembleCurvature(MultiLevelProblem& ml_prob);
void AssembleLevelSet(MultiLevelProblem& ml_prob);

void RestrictFineSystem(
    LinearImplicitSystem& system,
    MultiLevelMesh& mlMsh,
    const unsigned levelF,
    const unsigned targetLevel);

void CopyRestrictedSystem(
    LinearImplicitSystem& fineSystem,
    const unsigned targetLevel,
    LinearEquationSolver& coarseSolver);

void CopyAndProlongRestrictedField(
    MultiLevelSolution& fineSol,
    MultiLevelMesh& fineMesh,
    MultiLevelSolution& restrictedSol,
    const std::string& fieldName,
    const unsigned level0,
    const unsigned targetLevel,
    const unsigned levelF);

void AssembleCurvatureRestricted(MultiLevelProblem& mlProbRestricted);
void AssembleLevelSetRestricted(MultiLevelProblem& mlProbRestricted);

void AssembleNormal(MultiLevelProblem& ml_prob) {
  AssembleNormalFromField(ml_prob, "Psi");
}

void AssembleNormalAux(MultiLevelProblem& ml_prob) {
  AssembleNormalFromField(ml_prob, "PsiAux");
}

void RestrictFineSystem(
    LinearImplicitSystem& system,
    MultiLevelMesh& mlMsh,
    const unsigned levelF,
    const unsigned targetLevel) {

  std::vector<SparseMatrix*> PP =
      system.GetProjectionMatrix();

  std::vector<SparseMatrix*> RR =
      system.GetRestrictionMatrix();

  std::vector<SparseMatrix*> PPamr =
      system.GetAMRProjectionMatrix();

  std::vector<SparseMatrix*> RRamr =
      system.GetAMRRestrictionMatrix();

  std::vector<LinearEquationSolver*> LinSolver =
      system.GetLinearSolver();

  for(unsigned level = levelF;
      level > targetLevel;
      --level) {

    if(!mlMsh.GetLevel(level)->GetIfHomogeneous()
       && level == levelF) {

      if(!RRamr[level]) {

        LinSolver[level]->_RESC->matrix_mult_transpose(
            *LinSolver[level]->_RES,
            *PPamr[level]);

        *(LinSolver[level]->_RES) =
            *(LinSolver[level]->_RESC);

        LinSolver[level]->SwapMatrices();

        LinSolver[level]->_KK->matrix_PtAP(
            *PPamr[level],
            *LinSolver[level]->_KKamr,
            false);
      }
      else {

        LinSolver[level]->_RESC->matrix_mult(
            *LinSolver[level]->_RES,
            *RRamr[level]);

        *(LinSolver[level]->_RES) =
            *(LinSolver[level]->_RESC);

        LinSolver[level]->SwapMatrices();

        LinSolver[level]->_KK->matrix_ABC(
            *RRamr[level],
            *LinSolver[level]->_KKamr,
            *PPamr[level],
            false);
      }
    }

    if(!RR[level]) {

      LinSolver[level - 1]->_RES->matrix_mult_transpose(
          *LinSolver[level]->_RES,
          *PP[level]);

      LinSolver[level - 1]->_KK->matrix_PtAP(
          *PP[level],
          *LinSolver[level]->_KK,
          false);
    }
    else {

      LinSolver[level - 1]->_RES->matrix_mult(
          *LinSolver[level]->_RES,
          *RR[level]);

      LinSolver[level - 1]->_KK->matrix_ABC(
          *RR[level],
          *LinSolver[level]->_KK,
          *PP[level],
          false);
    }
  }
}

void CopyRestrictedSystem(
    LinearImplicitSystem& fineSystem,
    const unsigned targetLevel,
    LinearEquationSolver& coarseSolver) {

  SparseMatrix* KKC = coarseSolver._KK;
  NumericVector* RESC = coarseSolver._RES;

  KKC->zero();
  RESC->zero();

  RESC->close();
  KKC->close();

  std::vector<LinearEquationSolver*> LinSolver =
      fineSystem.GetLinearSolver();

  KKC->matrix_add(
      1.,
      *LinSolver[targetLevel]->_KK,
      "different_nonzero_pattern");

  *RESC += *LinSolver[targetLevel]->_RES;

  RESC->close();
  KKC->close();
}

// void AssembleCurvatureRestricted(
//     MultiLevelProblem& mlProbRestricted) {

//   LinearImplicitSystem& restrictedSystem =
//       mlProbRestricted
//           .get_system<LinearImplicitSystem>("K");

//   const unsigned levelRestricted =
//       restrictedSystem.GetLevelToAssemble();

//   MultiphaseParams param =
//       mlProbRestricted.GetMultiphaseParams();

//   MultiLevelProblem* mlProbFine =
//       param.mlProbF;

//   const unsigned levelF =
//       param.levelF;

//   const unsigned level0 =
//       param.level0;

//   const unsigned targetLevel =
//       levelRestricted + level0;

//   LinearImplicitSystem& fineSystem =
//       mlProbFine
//           ->get_system<LinearImplicitSystem>("K");

//   AssembleCurvature(*mlProbFine);

//   RestrictFineSystem(
//       fineSystem,
//       *mlProbFine->_ml_msh,
//       levelF,
//       targetLevel);

//   CopyRestrictedSystem(
//       fineSystem,
//       targetLevel,
//       *restrictedSystem._LinSolver[levelRestricted]);
// }

// void AssembleLevelSetRestricted(
//     MultiLevelProblem& mlProbRestricted) {

//   LinearImplicitSystem& restrictedSystem =
//       mlProbRestricted
//           .get_system<LinearImplicitSystem>("PsiAux");

//   const unsigned levelRestricted =
//       restrictedSystem.GetLevelToAssemble();

//   MultiphaseParams param =
//       mlProbRestricted.GetMultiphaseParams();

//   MultiLevelProblem* mlProbFine =
//       param.mlProbF;

//   const unsigned levelF =
//       param.levelF;

//   const unsigned level0 =
//       param.level0;

//   const unsigned targetLevel =
//       levelRestricted + level0;

//   LinearImplicitSystem& fineSystem =
//       mlProbFine
//           ->get_system<LinearImplicitSystem>("PsiAux");

//   AssembleLevelSet(*mlProbFine);

//   RestrictFineSystem(
//       fineSystem,
//       *mlProbFine->_ml_msh,
//       levelF,
//       targetLevel);

//   CopyRestrictedSystem(
//       fineSystem,
//       targetLevel,
//       *restrictedSystem._LinSolver[levelRestricted]);
// }

void CopyAndProlongRestrictedField(
    MultiLevelSolution& fineSol,
    MultiLevelMesh& fineMesh,
    MultiLevelSolution& restrictedSol,
    const std::string& fieldName,
    const unsigned level0,
    const unsigned targetLevel,
    const unsigned levelF) {

  const unsigned fineIndex = fineSol.GetIndex(fieldName.c_str());
  const unsigned restrictedIndex = restrictedSol.GetIndex(fieldName.c_str());
  const unsigned fieldType = fineSol.GetSolutionType(fineIndex);

  for(unsigned l = level0; l <= targetLevel; ++l) {

    Solution* src = restrictedSol.GetSolutionLevel(l - level0);
    Solution* dst = fineSol.GetSolutionLevel(l);

    *(dst->_Sol[fineIndex]) = *(src->_Sol[restrictedIndex]);
    dst->_Sol[fineIndex]->close();
  }

  for(unsigned l = targetLevel; l < levelF; ++l) {

    Solution* solC = fineSol.GetSolutionLevel(l);
    Solution* solF = fineSol.GetSolutionLevel(l + 1);
    Mesh* mshF = fineMesh.GetLevel(l + 1);

    solF->_Sol[fineIndex]->matrix_mult(
        *(solC->_Sol[fineIndex]),
        *(mshF->GetCoarseToFineProjection(fieldType)));

    solF->_Sol[fineIndex]->close();
  }
}

// void AssembleCurvature(MultiLevelProblem& ml_prob) {

//   MultiphaseParams params = ml_prob.GetMultiphaseParams();

//   LinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<LinearImplicitSystem> ("K");
//   const unsigned level = mlPdeSys->GetLevelToAssemble();

//   Mesh* msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
//   elem* el = msh->el;  // pointer to the elem object in msh (level)

//   MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
//   Solution* sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

//   LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
//   SparseMatrix* KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
//   NumericVector* RES = pdeSys->_RES; // pointer to the global residual std::vector object in pdeSys (level)

//   MatSetOption((static_cast< PetscMatrix* >(KK))->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

//   const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

//   unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

//   std::vector<unsigned> solNIndex(dim);
//   solNIndex[0] = mlSol->GetIndex("NX");
//   solNIndex[1] = mlSol->GetIndex("NY");
//   if(dim == 3) solNIndex[2] = mlSol->GetIndex("NZ");
//   unsigned solNType = mlSol->GetSolutionType("NX");

//   unsigned  solKIndex;
//   solKIndex = mlSol->GetIndex("K");    // get the position of "U" in the ml_sol object

//   unsigned  solKPdeIndex;
//   solKPdeIndex = mlPdeSys->GetSolPdeIndex("K");    // get the position of "U" in the pdeSys object


//   unsigned solKType = mlSol->GetSolutionType(solKIndex);

//   // std::vector < double >  psi; // local solution

//   std::vector < std::vector < double > > coordX(dim);    // local coordinates
//   unsigned solXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

//   std::vector < std::vector < double > > normal(dim);
//   std::vector < double >  K;

//   std::vector <double> phi;  // local test function for velocity
//   std::vector <double> phi_x; // local test function first order partial derivatives
//   std::vector <double> bdphi;  // local test function for velocity
//   std::vector <double> bdphi_x;

//   std::vector <double> phiN;
//   std::vector <double> phiN_x;
//   std::vector <double> bdphiN;  // local test function for velocity
//   std::vector <double> bdphiN_x;

//   std::vector < double> normal_face;
//   std::vector < double> normal_faceN;
//   double weight_face = 0.;
//   double weight_faceN = 0.;


//   double weight; // gauss point weight
//   double weightN;

//   std::vector< unsigned > sysDof; // local to global pdeSys dofs
//   std::vector< double > Res; // local redidual std::vector
//   std::vector < double > Jac;

//   KK->zero();
//   RES->zero();

//   // element loop: each process loops only on the elements that owns
//   for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

//     short unsigned ielGeom = msh->GetElementType(iel);

//     unsigned nDofs = msh->GetElementDofNumber(iel, solKType);
//     unsigned nDofsN = msh->GetElementDofNumber(iel, solNType);
//     unsigned nDofsX = msh->GetElementDofNumber(iel, solXType);

//     // resize local arrays
//     sysDof.resize(nDofs);
//     Res.assign(nDofs, 0.);
//     Jac.assign(nDofs * nDofs, 0.);

//     K.resize(nDofs);
//     for(unsigned  d = 0; d < dim; d++) {
//       normal[d].resize(nDofsN);
//       coordX[d].resize(nDofsX);
//     }


//     for(unsigned i = 0; i < nDofs; i++) {
//       unsigned KDof  = msh->GetSolutionDof(i, iel, solKType);
//       K[i] = (*sol->_Sol[solKIndex])(KDof);
//       sysDof[i] = pdeSys->GetSystemDof(solKIndex, solKPdeIndex, i, iel);
//     }

//     for(unsigned i = 0; i < nDofsN; i++) {
//       unsigned normalDof = msh->GetSolutionDof(i, iel, solNType);
//       for (unsigned d = 0; d < dim; d++) {
//         normal[d][i] = (*sol->_Sol[solNIndex[d]])(normalDof);
//       }
//     }

//     for(unsigned i = 0; i < nDofsX; i++) {
//       unsigned coordXDof  = msh->GetSolutionDof(i, iel, solXType);
//       for(unsigned k = 0; k < dim; k++) {
//         coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);
//       }
//     }

//     const elem_type *femK = msh->_finiteElement[ielGeom][solKType];
//     const elem_type *femN = msh->_finiteElement[ielGeom][solNType];
//     const elem_type* femX = msh->_finiteElement[ielGeom][solXType];

//     double cellMeasure = 0.;

//     std::vector<double> phiX;
//     std::vector<double> phiX_x;
//     double weightX = 0.;

//     for (unsigned ig = 0; ig < femX->GetGaussPointNumber(); ig++) {

//       femX->Jacobian(coordX, ig, weightX, phiX, phiX_x);

//       cellMeasure += weightX;
//     }

//     const double h = std::pow(cellMeasure, 1.0 / static_cast<double>(dim));

//     double alpha = 1.;
//     double hC = params.hC;
//     double epsilon = hC * hC * alpha;

//     // *** Gauss point loop ***
//     for(unsigned ig = 0; ig < femN->GetGaussPointNumber(); ig++) {
//       // *** get gauss point weight, test function and test function partial derivatives ***
//       femK->Jacobian(coordX, ig, weight, phi, phi_x);
//       femN->Jacobian(coordX, ig, weightN, phiN, phiN_x);

//       double K_g = 0.;
//       std::vector<double> gradK_g(dim, 0.);

//       for (unsigned j = 0; j < nDofs; j++) {
//         K_g += K[j] * phi[j];

//         for (unsigned d = 0; d < dim; d++) {
//           gradK_g[d] += K[j] * phi_x[j * dim + d];
//         }
//       }

//       std::vector<double> normal_g(dim, 0.);
//       for(unsigned d = 0; d < dim; d++) {
//         for (unsigned j = 0; j < nDofsN; j++) {
//           normal_g[d] += normal[d][j] * phiN[j];
//         }
//       }
//       double abs_normal_g = 0.;
//       for(unsigned d = 0; d < dim; d++) {
//         abs_normal_g += normal_g[d] * normal_g[d];
//       }
//       abs_normal_g = sqrt(abs_normal_g);

//       for(unsigned d = 0; d < dim; d++) {
//         normal_g[d] /= abs_normal_g;
//       }

//       // *** phiV_i loop ***
//       for(unsigned i = 0; i < nDofs; i++) {
//         double rhs = 0.;
//         for(unsigned  d = 0; d < dim; d++) {  //momentum equation in k
//           rhs -= phi_x[i * dim + d] * normal_g[d];
//           rhs -= epsilon * phi_x[i * dim + d] * gradK_g[d];
//         }
//         rhs -= K_g/*K[i]*/ * phi[i];
//         Res[i] += rhs * weight;
//       } // end phiV_i loop


//       //--------------------------------------------------------------------------------------------------------
//       // Add the local Matrix/Vector into the global Matrix/Vector

//       for(unsigned i = 0; i < nDofs; i++) {
//         // for(unsigned I = 0; I < dim; I++) { //row velocity blocks or dimension
//         unsigned VIrow = i;
//         for(unsigned j = 0; j < nDofs; j++) {
//           unsigned VIcolumn = j;

//           double helmotz_filter = 0.;

//           for(unsigned d = 0; d < dim; d++) {
//             helmotz_filter += phi_x[i * dim + d] *
//                       phi_x[j * dim + d];
//           }

//           // VIcolumn = VIrow;
//           Jac[ VIrow * nDofs + VIcolumn] += (phi[i] * phi[j] + epsilon * helmotz_filter) * weight ; // inertia


//         }
//         // }
//       }
//     }

//     // *** Face Gauss point loop (boundary Integral) ***
//     for ( unsigned jface = 0; jface < msh->GetElementFaceNumber ( iel ); jface++ ) {
//       int faceIndex = el->GetBoundaryIndex(iel, jface);
//       // look for boundary faces

//       if ( faceIndex > 0 ) {
//         const unsigned faceGeom = msh->GetElementFaceType ( iel, jface );
//         unsigned faceDofs = msh->GetElementFaceDofNumber (iel, jface, solKType);
//         unsigned faceDofsN = msh->GetElementFaceDofNumber (iel, jface, solNType);
//         unsigned faceDofsX = msh->GetElementFaceDofNumber (iel, jface, solXType);
//         std::vector  < std::vector  <  double> > faceCoordinates ( dim ); // A matrix holding the face coordinates rowwise.
//         for ( int k = 0; k < dim; k++ ) {
//           faceCoordinates[k].resize (faceDofsX);
//         }
//         for ( unsigned i = 0; i < faceDofsX; i++ ) {
//           unsigned inode = msh->GetLocalFaceVertexIndex ( iel, jface, i ); // face-to-element local node mapping.
//           for ( unsigned k = 0; k < dim; k++ ) {
//             faceCoordinates[k][i] =  coordX[k][inode]; // We extract the local coordinates on the face from local coordinates on the element.
//           }
//         }
//         for ( unsigned ig = 0; ig  <  msh->_finiteElement[faceGeom][solKType]->GetGaussPointNumber(); ig++ ) {
//           // We call the method GetGaussPointNumber from the object finiteElement in the mesh object msh.

//           msh->_finiteElement[faceGeom][solKType]->JacobianSur ( faceCoordinates, ig, weight_face, bdphi, bdphi_x, normal_face );
//           msh->_finiteElement[faceGeom][solNType]->JacobianSur ( faceCoordinates, ig, weight_faceN, bdphiN, bdphiN_x, normal_faceN );

//           std::vector<double> normal_g(dim, 0.);
//           for(unsigned d = 0; d < dim; d++) {
//             for (unsigned j = 0; j < faceDofsN; j++) {
//               unsigned jnode = msh->GetLocalFaceVertexIndex (iel, jface, j );
//               normal_g[d] += normal[d][jnode] * bdphiN[j];
//             }
//           }

//           // *** phi_i loop ***
//           for ( unsigned i = 0; i < faceDofs; i++ ) {
//             double rhs_bd = 0;
//             unsigned inode = msh->GetLocalFaceVertexIndex (iel, jface, i );
//             for( unsigned d = 0; d < dim; d++) {
//               rhs_bd +=  bdphi[i] * normal_face[d] * normal_g[d];
//             }
//             Res[inode] += rhs_bd * weight_face;
//           }
//         }
//       }
//     }

//     RES->add_vector_blocked(Res, sysDof);
//     KK->add_matrix_blocked(Jac, sysDof, sysDof);


//   } //end element loop for each process

//   RES->close();
//   KK->close();


// }

void AssembleCurvature(MultiLevelProblem& ml_prob) {

  LinearImplicitSystem* mlPdeSysOut = &ml_prob.get_system<LinearImplicitSystem>("K");
  const unsigned levelOut = mlPdeSysOut->GetLevelToAssemble();

  MultiphaseParams mParam = ml_prob.GetMultiphaseParams();
  const bool restricted = (mParam.mlProbF != nullptr);

  MultiLevelProblem* mlProbAssembly = restricted ? mParam.mlProbF : &ml_prob;
  LinearImplicitSystem* mlPdeSys = &mlProbAssembly->get_system<LinearImplicitSystem>("K");

  const unsigned level = restricted ? mParam.levelF : levelOut;
  const double hC = mParam.hC;

  Mesh* msh = mlProbAssembly->_ml_msh->GetLevel(level);
  elem* el = msh->el;

  MultiLevelSolution* mlSol = mlProbAssembly->_ml_sol;
  Solution* sol = mlSol->GetSolutionLevel(level);

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];
  SparseMatrix* KK = pdeSys->_KK;
  NumericVector* RES = pdeSys->_RES;

  MatSetOption((static_cast<PetscMatrix*>(KK))->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  const unsigned dim = msh->GetDimension();
  const unsigned iproc = msh->processor_id();

  std::vector<unsigned> solNIndex(dim);
  solNIndex[0] = mlSol->GetIndex("NX");
  solNIndex[1] = mlSol->GetIndex("NY");
  if(dim == 3) solNIndex[2] = mlSol->GetIndex("NZ");

  const unsigned solNType = mlSol->GetSolutionType("NX");
  const unsigned solKIndex = mlSol->GetIndex("K");
  const unsigned solKPdeIndex = mlPdeSys->GetSolPdeIndex("K");
  const unsigned solKType = mlSol->GetSolutionType(solKIndex);
  const unsigned solXType = 2;

  std::vector<std::vector<double>> coordX(dim);
  std::vector<std::vector<double>> normal(dim);
  std::vector<double> K;

  std::vector<double> phi, phi_x;
  std::vector<double> bdphi, bdphi_x;
  std::vector<double> phiN, phiN_x;
  std::vector<double> bdphiN, bdphiN_x;

  std::vector<double> normal_face;
  std::vector<double> normal_faceN;

  double weight = 0.;
  double weightN = 0.;
  double weight_face = 0.;
  double weight_faceN = 0.;

  std::vector<unsigned> sysDof;
  std::vector<double> Res;
  std::vector<double> Jac;

  KK->zero();
  RES->zero();

  const double alpha = 1.e-2;
  const double epsilon = alpha * hC * hC;

  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; ++iel) {

    const short unsigned ielGeom = msh->GetElementType(iel);

    const unsigned nDofs = msh->GetElementDofNumber(iel, solKType);
    const unsigned nDofsN = msh->GetElementDofNumber(iel, solNType);
    const unsigned nDofsX = msh->GetElementDofNumber(iel, solXType);

    sysDof.resize(nDofs);
    Res.assign(nDofs, 0.);
    Jac.assign(nDofs * nDofs, 0.);

    K.resize(nDofs);

    for(unsigned d = 0; d < dim; ++d) {
      normal[d].resize(nDofsN);
      coordX[d].resize(nDofsX);
    }

    for(unsigned i = 0; i < nDofs; ++i) {
      const unsigned KDof = msh->GetSolutionDof(i, iel, solKType);
      K[i] = (*sol->_Sol[solKIndex])(KDof);
      sysDof[i] = pdeSys->GetSystemDof(solKIndex, solKPdeIndex, i, iel);
    }

    for(unsigned i = 0; i < nDofsN; ++i) {
      const unsigned normalDof = msh->GetSolutionDof(i, iel, solNType);
      for(unsigned d = 0; d < dim; ++d)
        normal[d][i] = (*sol->_Sol[solNIndex[d]])(normalDof);
    }

    for(unsigned i = 0; i < nDofsX; ++i) {
      const unsigned coordXDof = msh->GetSolutionDof(i, iel, solXType);
      for(unsigned d = 0; d < dim; ++d)
        coordX[d][i] = (*msh->_topology->_Sol[d])(coordXDof);
    }

    const elem_type* femK = msh->_finiteElement[ielGeom][solKType];
    const elem_type* femN = msh->_finiteElement[ielGeom][solNType];

    for(unsigned ig = 0; ig < femN->GetGaussPointNumber(); ++ig) {

      femK->Jacobian(coordX, ig, weight, phi, phi_x);
      femN->Jacobian(coordX, ig, weightN, phiN, phiN_x);

      double K_g = 0.;
      std::vector<double> gradK_g(dim, 0.);

      for(unsigned j = 0; j < nDofs; ++j) {
        K_g += K[j] * phi[j];

        for(unsigned d = 0; d < dim; ++d)
          gradK_g[d] += K[j] * phi_x[j * dim + d];
      }

      std::vector<double> normal_g(dim, 0.);

      for(unsigned d = 0; d < dim; ++d) {
        for(unsigned j = 0; j < nDofsN; ++j)
          normal_g[d] += normal[d][j] * phiN[j];
      }

      double abs_normal_g = 0.;

      for(unsigned d = 0; d < dim; ++d)
        abs_normal_g += normal_g[d] * normal_g[d];

      abs_normal_g = std::sqrt(abs_normal_g);

      if(abs_normal_g > 1.e-14) {
        for(unsigned d = 0; d < dim; ++d)
          normal_g[d] /= abs_normal_g;
      }

      for(unsigned i = 0; i < nDofs; ++i) {

        double rhs = 0.;

        for(unsigned d = 0; d < dim; ++d) {
          rhs -= phi_x[i * dim + d] * normal_g[d];
          rhs -= epsilon * phi_x[i * dim + d] * gradK_g[d];
        }

        rhs -= K_g * phi[i];
        Res[i] += rhs * weight;

        for(unsigned j = 0; j < nDofs; ++j) {

          double helmholtz_filter = 0.;

          for(unsigned d = 0; d < dim; ++d)
            helmholtz_filter += phi_x[i * dim + d] * phi_x[j * dim + d];

          Jac[i * nDofs + j] +=
              (phi[i] * phi[j] + epsilon * helmholtz_filter) * weight;
        }
      }
    }

    for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); ++jface) {

      const int faceIndex = el->GetBoundaryIndex(iel, jface);

      if(faceIndex > 0) {

        const unsigned faceGeom = msh->GetElementFaceType(iel, jface);
        const unsigned faceDofs = msh->GetElementFaceDofNumber(iel, jface, solKType);
        const unsigned faceDofsN = msh->GetElementFaceDofNumber(iel, jface, solNType);
        const unsigned faceDofsX = msh->GetElementFaceDofNumber(iel, jface, solXType);

        std::vector<std::vector<double>> faceCoordinates(dim);

        for(unsigned d = 0; d < dim; ++d)
          faceCoordinates[d].resize(faceDofsX);

        for(unsigned i = 0; i < faceDofsX; ++i) {
          const unsigned inode = msh->GetLocalFaceVertexIndex(iel, jface, i);

          for(unsigned d = 0; d < dim; ++d)
            faceCoordinates[d][i] = coordX[d][inode];
        }

        for(unsigned ig = 0;
            ig < msh->_finiteElement[faceGeom][solKType]->GetGaussPointNumber();
            ++ig) {

          msh->_finiteElement[faceGeom][solKType]->JacobianSur(
              faceCoordinates, ig, weight_face, bdphi, bdphi_x, normal_face);

          msh->_finiteElement[faceGeom][solNType]->JacobianSur(
              faceCoordinates, ig, weight_faceN, bdphiN, bdphiN_x, normal_faceN);

          std::vector<double> normal_g(dim, 0.);

          for(unsigned d = 0; d < dim; ++d) {
            for(unsigned j = 0; j < faceDofsN; ++j) {
              const unsigned jnode = msh->GetLocalFaceVertexIndex(iel, jface, j);
              normal_g[d] += normal[d][jnode] * bdphiN[j];
            }
          }

          for(unsigned i = 0; i < faceDofs; ++i) {

            double rhs_bd = 0.;
            const unsigned inode = msh->GetLocalFaceVertexIndex(iel, jface, i);

            for(unsigned d = 0; d < dim; ++d)
              rhs_bd += bdphi[i] * normal_face[d] * normal_g[d];

            Res[inode] += rhs_bd * weight_face;
          }
        }
      }
    }

    RES->add_vector_blocked(Res, sysDof);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);
  }

  RES->close();
  KK->close();

  if(restricted) {
    const unsigned targetLevel = levelOut + mParam.level0;

    RestrictFineSystem(
        *mlPdeSys,
        *mlProbAssembly->_ml_msh,
        mParam.levelF,
        targetLevel);

    CopyRestrictedSystem(
        *mlPdeSys,
        targetLevel,
        *mlPdeSysOut->_LinSolver[levelOut]);
  }
}

void AssembleCurvatureLumped(MultiLevelProblem& ml_prob) {

  LinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<LinearImplicitSystem> ("K");
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh* msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem* el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix* KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* RES = pdeSys->_RES; // pointer to the global residual std::vector object in pdeSys (level)

  MatSetOption((static_cast< PetscMatrix* >(KK))->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  std::vector<unsigned> solNIndex(dim);
  solNIndex[0] = mlSol->GetIndex("NX");
  solNIndex[1] = mlSol->GetIndex("NY");
  if(dim == 3) solNIndex[2] = mlSol->GetIndex("NZ");
  unsigned solNType = mlSol->GetSolutionType("NX");

  unsigned  solKIndex;
  solKIndex = mlSol->GetIndex("K");    // get the position of "U" in the ml_sol object

  unsigned  solKPdeIndex;
  solKPdeIndex = mlPdeSys->GetSolPdeIndex("K");    // get the position of "U" in the pdeSys object


  unsigned solKType = mlSol->GetSolutionType(solKIndex);

  // std::vector < double >  psi; // local solution

  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  unsigned solXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector < std::vector < double > > normal(dim);
  std::vector < double >  K;

  std::vector <double> phi;  // local test function for velocity
  std::vector <double> phi_x; // local test function first order partial derivatives
  std::vector <double> bdphi;  // local test function for velocity
  std::vector <double> bdphi_x;

  std::vector <double> phiN;
  std::vector <double> phiN_x;
  std::vector <double> bdphiN;  // local test function for velocity
  std::vector <double> bdphiN_x;

  std::vector < double> normal_face;
  std::vector < double> normal_faceN;
  double weight_face = 0.;
  double weight_faceN = 0.;


  double weight; // gauss point weight
  double weightN;

  std::vector< unsigned > sysDof; // local to global pdeSys dofs
  std::vector< double > Res; // local redidual std::vector
  std::vector < double > Jac;

  KK->zero();
  RES->zero();

  // GATHER NEIGHBOUR LEVELS FOR ALL FINE LEVEL ELEMENTS

  const unsigned nprocs = msh->n_processors();

  const unsigned firstElem = msh->_elementOffset[iproc];
  const unsigned lastElem  = msh->_elementOffset[iproc + 1];
  const unsigned nLocalElem = lastElem - firstElem;

  std::vector<std::vector<int>> neighElem(nLocalElem);
  std::vector<std::vector<int>> neighLevel(nLocalElem);

  std::vector<std::vector<int>> requestedElem(nprocs);
  std::vector<std::vector<std::pair<unsigned, unsigned>>> remoteSlot(nprocs);

  for(unsigned iel = firstElem; iel < lastElem; ++iel) {

      const unsigned localIel = iel - firstElem;

      const int iel_level = msh->el->GetElementLevel(iel);

      if(iel_level != static_cast<int>(level)) {
          continue;
      }

      const unsigned nFaces = msh->GetElementFaceNumber(iel);

      neighElem[localIel].assign(nFaces, -1);
      neighLevel[localIel].assign(nFaces, -2);

      for(unsigned jface = 0; jface < nFaces; ++jface) {

          const int jel = el->GetFaceElementIndex(iel, jface) - 1;

          neighElem[localIel][jface] = jel;

          if(jel < 0) {
              neighLevel[localIel][jface] = -1;
              continue;
          }

          const unsigned jproc = msh->IsdomBisectionSearch(jel, 3);

          if(jproc == iproc) {
              neighLevel[localIel][jface] = msh->el->GetElementLevel(jel);
          }
          else {
              requestedElem[jproc].push_back(jel);
              remoteSlot[jproc].push_back(std::make_pair(localIel, jface));
          }
      }
  }

  std::vector<int> sendCounts(nprocs, 0);
  std::vector<int> recvCounts(nprocs, 0);

  for(unsigned p = 0; p < nprocs; ++p) {
      sendCounts[p] = static_cast<int>(requestedElem[p].size());
  }

  MPI_Alltoall(sendCounts.data(),1,MPI_INT,recvCounts.data(),1,MPI_INT,PETSC_COMM_WORLD);

  std::vector<int> sendDispls(nprocs, 0);
  std::vector<int> recvDispls(nprocs, 0);

  for(unsigned p = 1; p < nprocs; ++p) {
      sendDispls[p] = sendDispls[p - 1] + sendCounts[p - 1];
      recvDispls[p] = recvDispls[p - 1] + recvCounts[p - 1];
  }

  int totalSend = 0;
  int totalRecv = 0;

  if(nprocs > 0) {
      totalSend = sendDispls[nprocs - 1] + sendCounts[nprocs - 1];
      totalRecv = recvDispls[nprocs - 1] + recvCounts[nprocs - 1];
  }

  std::vector<int> sendElem(totalSend);

  for(unsigned p = 0; p < nprocs; ++p) {
      for(unsigned q = 0; q < requestedElem[p].size(); ++q) {
          sendElem[sendDispls[p] + static_cast<int>(q)] = requestedElem[p][q];
      }
  }

  std::vector<int> recvElem(totalRecv);

  MPI_Alltoallv(sendElem.data(),sendCounts.data(),sendDispls.data(),MPI_INT,
      recvElem.data(),recvCounts.data(),recvDispls.data(),MPI_INT,PETSC_COMM_WORLD);

  std::vector<int> sendLevelBack(totalRecv, -1);

  for(int q = 0; q < totalRecv; ++q) {
      const int jel = recvElem[q];
      sendLevelBack[q] =msh->el->GetElementLevel(jel);
  }

  std::vector<int> recvLevelBack(totalSend, -1);

  MPI_Alltoallv(sendLevelBack.data(),recvCounts.data(),recvDispls.data(),MPI_INT,
      recvLevelBack.data(),sendCounts.data(),sendDispls.data(),MPI_INT,PETSC_COMM_WORLD);

  for(unsigned p = 0; p < nprocs; ++p) {
      for(unsigned q = 0; q < remoteSlot[p].size(); ++q) {
          const unsigned localIel = remoteSlot[p][q].first;
          const unsigned jface = remoteSlot[p][q].second;
          const int position = sendDispls[p] + static_cast<int>(q);
          neighLevel[localIel][jface] = recvLevelBack[position];
      }
  }

  // END GATHER

  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    int iel_level = msh->el->GetElementLevel(iel);

    if (iel_level != level)
      continue;

    short unsigned ielGeom = msh->GetElementType(iel);

    unsigned nDofs = msh->GetElementDofNumber(iel, solKType);
    unsigned nDofsN = msh->GetElementDofNumber(iel, solNType);
    unsigned nDofsX = msh->GetElementDofNumber(iel, solXType);

    // resize local arrays
    sysDof.resize(nDofs);
    Res.assign(nDofs, 0.);
    Jac.assign(nDofs * nDofs, 0.);

    K.resize(nDofs);
    for(unsigned  d = 0; d < dim; d++) {
      normal[d].resize(nDofsN);
      coordX[d].resize(nDofsX);
    }


    for(unsigned i = 0; i < nDofs; i++) {
      unsigned KDof  = msh->GetSolutionDof(i, iel, solKType);
      K[i] = (*sol->_Sol[solKIndex])(KDof);
      sysDof[i] = pdeSys->GetSystemDof(solKIndex, solKPdeIndex, i, iel);
    }

    for(unsigned i = 0; i < nDofsN; i++) {
      unsigned normalDof = msh->GetSolutionDof(i, iel, solNType);
      for (unsigned d = 0; d < dim; d++) {
        normal[d][i] = (*sol->_Sol[solNIndex[d]])(normalDof);
      }
    }

    for(unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, solXType);
      for(unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);
      }
    }

    const elem_type *femK = msh->_finiteElement[ielGeom][solKType];
    const elem_type *femN = msh->_finiteElement[ielGeom][solNType];
    const elem_type* femX = msh->_finiteElement[ielGeom][solXType];

    double cellMeasure = 0.;

    std::vector<double> phiX;
    std::vector<double> phiX_x;
    double weightX = 0.;

    for (unsigned ig = 0; ig < femX->GetGaussPointNumber(); ig++) {

      femX->Jacobian(coordX, ig, weightX, phiX, phiX_x);

      cellMeasure += weightX;
    }

    const double h = std::pow(cellMeasure, 1.0 / static_cast<double>(dim));

    const double alpha = 0.;

    const double epsilon = alpha * h * h ;// */ 1.e-6;

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < femN->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      femK->Jacobian(coordX, ig, weight, phi, phi_x);
      femN->Jacobian(coordX, ig, weightN, phiN, phiN_x);

      double K_g = 0.;
      std::vector<double> gradK_g(dim, 0.);

      for (unsigned j = 0; j < nDofs; j++) {
        K_g += K[j] * phi[j];

        for (unsigned d = 0; d < dim; d++) {
          gradK_g[d] += K[j] * phi_x[j * dim + d];
        }
      }

      std::vector<double> normal_g(dim, 0.);
      for(unsigned d = 0; d < dim; d++) {
        for (unsigned j = 0; j < nDofsN; j++) {
          normal_g[d] += normal[d][j] * phiN[j];
        }
      }

      // *** phiV_i loop ***
      for(unsigned i = 0; i < nDofs; i++) {
        double rhs = 0.;
        for(unsigned  d = 0; d < dim; d++) {  //momentum equation in k
          rhs -= phi_x[i * dim + d] * normal_g[d];
          rhs -= epsilon * phi_x[i * dim + d] * gradK_g[d];
        }
        rhs -= /*K_g*/K[i] * phi[i];
        Res[i] += rhs * weight;
      } // end phiV_i loop


      //--------------------------------------------------------------------------------------------------------
      // Add the local Matrix/Vector into the global Matrix/Vector

      for(unsigned i = 0; i < nDofs; i++) {
        // for(unsigned I = 0; I < dim; I++) { //row velocity blocks or dimension
        unsigned VIrow = i;
        for(unsigned j = 0; j < nDofs; j++) {
          unsigned VIcolumn = j;

          double helmotz_filter = 0.;

          for(unsigned d = 0; d < dim; d++) {
            helmotz_filter += phi_x[i * dim + d] *
                      phi_x[j * dim + d];
          }

          VIcolumn = VIrow;
          Jac[ VIrow * nDofs + VIcolumn] += (phi[i] * phi[j] + epsilon * helmotz_filter) * weight ; // inertia


        }
        // }
      }
    }

    // *** Face Gauss point loop (boundary Integral) ***
    for ( unsigned jface = 0; jface < msh->GetElementFaceNumber ( iel ); jface++ ) {
      int faceIndex = el->GetBoundaryIndex(iel, jface);

      // int neigh_level = (faceIndex >= 0) ? msh->el->GetElementLevel(faceIndex) : -1;
      int neigh_level = neighLevel[iel - msh->_elementOffset[iproc]][jface];

      if ( /*faceIndex > 0*/  neigh_level != level) {
        const unsigned faceGeom = msh->GetElementFaceType ( iel, jface );
        unsigned faceDofs = msh->GetElementFaceDofNumber (iel, jface, solKType);
        unsigned faceDofsN = msh->GetElementFaceDofNumber (iel, jface, solNType);
        unsigned faceDofsX = msh->GetElementFaceDofNumber (iel, jface, solXType);
        std::vector  < std::vector  <  double> > faceCoordinates ( dim ); // A matrix holding the face coordinates rowwise.
        for ( int k = 0; k < dim; k++ ) {
          faceCoordinates[k].resize (faceDofsX);
        }
        for ( unsigned i = 0; i < faceDofsX; i++ ) {
          unsigned inode = msh->GetLocalFaceVertexIndex ( iel, jface, i ); // face-to-element local node mapping.
          for ( unsigned k = 0; k < dim; k++ ) {
            faceCoordinates[k][i] =  coordX[k][inode]; // We extract the local coordinates on the face from local coordinates on the element.
          }
        }
        for ( unsigned ig = 0; ig  <  msh->_finiteElement[faceGeom][solKType]->GetGaussPointNumber(); ig++ ) {
          // We call the method GetGaussPointNumber from the object finiteElement in the mesh object msh.

          msh->_finiteElement[faceGeom][solKType]->JacobianSur ( faceCoordinates, ig, weight_face, bdphi, bdphi_x, normal_face );
          msh->_finiteElement[faceGeom][solNType]->JacobianSur ( faceCoordinates, ig, weight_faceN, bdphiN, bdphiN_x, normal_faceN );

          std::vector<double> normal_g(dim, 0.);
          for(unsigned d = 0; d < dim; d++) {
            for (unsigned j = 0; j < faceDofsN; j++) {
              unsigned jnode = msh->GetLocalFaceVertexIndex (iel, jface, j );
              normal_g[d] += normal[d][jnode] * bdphiN[j];
            }
          }

          // *** phi_i loop ***
          for ( unsigned i = 0; i < faceDofs; i++ ) {
            double rhs_bd = 0;
            unsigned inode = msh->GetLocalFaceVertexIndex (iel, jface, i );
            for( unsigned d = 0; d < dim; d++) {
              rhs_bd +=  bdphi[i] * normal_face[d] * normal_g[d];
            }
            Res[inode] += rhs_bd * weight_face;
          }
        }
      }
    }

    RES->add_vector_blocked(Res, sysDof);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);


  } //end element loop for each process

  RES->close();
  KK->close();

  int glob_offset_dof = msh->_dofOffset[solKType][iproc];

  for(int i = 0; i < msh->_ownSize[solKType][iproc]; i++) {
    int global_dof = i + glob_offset_dof;
    if (fabs((*KK)(global_dof, global_dof)) > 1.e-20)
      sol->_Sol[solKIndex]->set(global_dof, (*RES)(global_dof) / ((*KK)(global_dof, global_dof)));
  }

  sol->_Sol[solKIndex]->close();


}

void AssembleNormalFromField(MultiLevelProblem& ml_prob, const char* psiName) {

  MultiphaseParams params = ml_prob.GetMultiphaseParams();

  LinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<LinearImplicitSystem> ("N");
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh* msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem* el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object

  SparseMatrix* KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* RES = pdeSys->_RES; // pointer to the global residual std::vector object in pdeSys (level)

  MatSetOption((static_cast< PetscMatrix* >(KK))->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  unsigned psiIndex = mlSol->GetIndex(psiName);
  unsigned psiType = mlSol->GetSolutionType(psiName);

  std::vector < unsigned > solNIndex(dim);
  solNIndex[0] = mlSol->GetIndex("NX");    // get the position of "U" in the ml_sol object
  solNIndex[1] = mlSol->GetIndex("NY");    // get the position of "V" in the ml_sol object
  if(dim == 3) solNIndex[2] = mlSol->GetIndex("NZ");       // get the position of "V" in the ml_sol object

  std::vector < unsigned > solNPdeIndex(dim);
  solNPdeIndex[0] = mlPdeSys->GetSolPdeIndex("NX");    // get the position of "U" in the pdeSys object
  solNPdeIndex[1] = mlPdeSys->GetSolPdeIndex("NY");    // get the position of "V" in the pdeSys object
  if(dim == 3) solNPdeIndex[2] = mlPdeSys->GetSolPdeIndex("NZ");

  unsigned solNType = mlSol->GetSolutionType(solNIndex[0]);

  std::vector < double >  psi; // local solution

  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  unsigned solXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector <double> phiN;  // local test function for velocity
  std::vector <double> phiN_x; // local test function first order partial derivatives

  std::vector <double> phiPsi;
  std::vector <double>  phiPsi_x;

  std::vector<std::vector<double>> N(dim);
  double weight; // gauss point weight
  double weightPsi;

  std::vector< unsigned > sysDof; // local to global pdeSys dofs
  std::vector< double > Res; // local redidual std::vector
  std::vector < double > Jac;

  KK->zero();
  RES->zero();

  // element loop: each process loops only on the elements that owns
  double alpha = 1.e-2;
  double hC = params.hC;
  double epsilon = hC * hC * alpha;
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    int iel_level = msh->el->GetElementLevel(iel);

    // if (iel_level != level)
    //   continue;

    short unsigned ielGeom = msh->GetElementType(iel);

    unsigned nDofsN = msh->GetElementDofNumber(iel, solNType);
    unsigned nDofsX = msh->GetElementDofNumber(iel, solXType);

    unsigned nDofs =  dim * nDofsN;

    // resize local arrays
    sysDof.resize(nDofs);
    Res.assign(nDofs, 0.);
    Jac.assign(nDofs * nDofs, 0.);

    for(unsigned  k = 0; k < dim; k++) {
      coordX[k].resize(nDofsX);
      N[k].resize(nDofsN);
    }

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsN; i++) {
      unsigned solNDof = msh->GetSolutionDof(i, iel, solNType);
      for(unsigned  d = 0; d < dim; d++) {
        N[d][i] = (*sol->_Sol[solNIndex[d]])(solNDof);
        sysDof[d * nDofsN + i] = pdeSys->GetSystemDof(solNIndex[d], solNPdeIndex[d], i, iel);
      }
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, solXType);
      for(unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);
      }
    }

    unsigned nDofsPsi;

    nDofsPsi = msh->GetElementDofNumber(iel, psiType);
    psi.resize(nDofsPsi);
    for(unsigned i = 0; i < nDofsPsi; i++) {
      unsigned psiDof = msh->GetSolutionDof(i, iel, psiType);
      psi[i] = (*sol->_Sol[psiIndex])(psiDof);
    }

    const elem_type *femPsi = msh->_finiteElement[ielGeom][psiType];
    const elem_type *femN = msh->_finiteElement[ielGeom][solNType];

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < femN->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      femN->Jacobian(coordX, ig, weight, phiN, phiN_x);
      femPsi->Jacobian(coordX, ig, weightPsi, phiPsi, phiPsi_x);


      std::vector<double> NN(dim, 0.);
      for (unsigned i = 0; i < nDofsPsi; i++) {
        for(unsigned d = 0; d < dim; d++) {
          NN[d] -= psi[i] * phiPsi_x[i * dim + d];
        }
      }
      double det = 0;
      for (unsigned d = 0; d < dim; d++) {
        det += NN[d] * NN[d];
      }
      det = std::sqrt(det + 1.e-10);
      for (unsigned d = 0; d < dim; d++) {
        NN[d] /= det;
      }

      std::vector<double> N_g(dim, 0.);
      std::vector<std::vector<double>> gradN_g(dim);
      for (unsigned d = 0; d < dim; d ++)
        gradN_g[d].resize(dim);
      for (unsigned i = 0; i < nDofsN; i++) {
        for(unsigned d = 0; d < dim; d++) {
          N_g[d] += N[d][i] * phiN[i];
          for(unsigned k = 0; k < dim; k++) {
            gradN_g[d][k] += N[d][i] * phiN_x[i * dim + k];
          }
        }
      }
      
      // *** phiV_i loop ***
      for(unsigned i = 0; i < nDofsN; i++) {
        for(unsigned  d = 0; d < dim; d++) {  //momentum equation in k
          double rhs = 0.;
          rhs += phiN[i] * (NN[d] - N_g[d]/*N[d][i]*/);
          for(unsigned  k = 0; k < dim; k++) {
            rhs -= epsilon * phiN_x[i * dim + k] * gradN_g[d][k];
          }
          Res[d * nDofsN + i] +=  rhs * weight;
        }
      } // end phiV_i loop


      //--------------------------------------------------------------------------------------------------------
      // Add the local Matrix/Vector into the global Matrix/Vector

      for(unsigned i = 0; i < nDofsN; i++) {
        for(unsigned d = 0; d < dim; d++) { //row velocity blocks or dimension
          unsigned VIrow = d * nDofsN + i;
          for(unsigned j = 0; j < nDofsN; j++) {
            unsigned VIcolumn = d * nDofsN + j;

            double helmotz_filter = 0;
            for(unsigned k = 0; k < dim; k++) {
              helmotz_filter += phiN_x[i * dim + k] *
                        phiN_x[j * dim + k];
            }

            // VIcolumn = VIrow;
            Jac[ VIrow * nDofs + VIcolumn] += (phiN[i] * phiN[j] + epsilon * helmotz_filter)* weight ; // inertia


          }
        }
      }
    }

    RES->add_vector_blocked(Res, sysDof);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);


  } //end element loop for each process

  RES->close();
  KK->close();



}

void AssembleNormalLumped(MultiLevelProblem& ml_prob) {
  LinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<LinearImplicitSystem> ("N");
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh* msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem* el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object

  SparseMatrix* KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* RES = pdeSys->_RES; // pointer to the global residual std::vector object in pdeSys (level)

  MatSetOption((static_cast< PetscMatrix* >(KK))->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  unsigned psiIndex = mlSol->GetIndex("Psi");
  unsigned psiType = mlSol->GetSolutionType("Psi");

  std::vector < unsigned > solNIndex(dim);
  solNIndex[0] = mlSol->GetIndex("NX");    // get the position of "U" in the ml_sol object
  solNIndex[1] = mlSol->GetIndex("NY");    // get the position of "V" in the ml_sol object
  if(dim == 3) solNIndex[2] = mlSol->GetIndex("NZ");       // get the position of "V" in the ml_sol object

  std::vector < unsigned > solNPdeIndex(dim);
  solNPdeIndex[0] = mlPdeSys->GetSolPdeIndex("NX");    // get the position of "U" in the pdeSys object
  solNPdeIndex[1] = mlPdeSys->GetSolPdeIndex("NY");    // get the position of "V" in the pdeSys object
  if(dim == 3) solNPdeIndex[2] = mlPdeSys->GetSolPdeIndex("NZ");

  unsigned solNType = mlSol->GetSolutionType(solNIndex[0]);

  std::vector < double >  psi; // local solution

  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  unsigned solXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector <double> phiN;  // local test function for velocity
  std::vector <double> phiN_x; // local test function first order partial derivatives

  std::vector <double> phiPsi;
  std::vector <double>  phiPsi_x;

  std::vector<std::vector<double>> N(dim);
  double weight; // gauss point weight
  double weightPsi;

  std::vector< unsigned > sysDof; // local to global pdeSys dofs
  std::vector< double > Res; // local redidual std::vector
  std::vector < double > Jac;

  KK->zero();
  RES->zero();

  // element loop: each process loops only on the elements that owns
  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    int iel_level = msh->el->GetElementLevel(iel);

    // if (iel_level != level)
    //   continue;

    short unsigned ielGeom = msh->GetElementType(iel);

    unsigned nDofsN = msh->GetElementDofNumber(iel, solNType);
    unsigned nDofsX = msh->GetElementDofNumber(iel, solXType);

    unsigned nDofs =  dim * nDofsN;

    // resize local arrays
    sysDof.resize(nDofs);
    Res.assign(nDofs, 0.);
    Jac.assign(nDofs * nDofs, 0.);

    for(unsigned  k = 0; k < dim; k++) {
      coordX[k].resize(nDofsX);
      N[k].resize(nDofsN);
    }

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsN; i++) {
      unsigned solNDof = msh->GetSolutionDof(i, iel, solNType);
      for(unsigned  d = 0; d < dim; d++) {
        N[d][i] = (*sol->_Sol[solNIndex[d]])(solNDof);
        sysDof[d * nDofsN + i] = pdeSys->GetSystemDof(solNIndex[d], solNPdeIndex[d], i, iel);
      }
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, solXType);
      for(unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);
      }
    }

    unsigned nDofsPsi;

    nDofsPsi = msh->GetElementDofNumber(iel, psiType);
    psi.resize(nDofsPsi);
    for(unsigned i = 0; i < nDofsPsi; i++) {
      unsigned psiDof = msh->GetSolutionDof(i, iel, psiType);
      psi[i] = (*sol->_Sol[psiIndex])(psiDof);
    }

    const elem_type *femPsi = msh->_finiteElement[ielGeom][psiType];
    const elem_type *femN = msh->_finiteElement[ielGeom][solNType];

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < femN->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      femN->Jacobian(coordX, ig, weight, phiN, phiN_x);
      femPsi->Jacobian(coordX, ig, weightPsi, phiPsi, phiPsi_x);


      std::vector<double> NN(dim, 0.);
      for (unsigned i = 0; i < nDofsPsi; i++) {
        for(unsigned d = 0; d < dim; d++) {
          NN[d] -= psi[i] * phiPsi_x[i * dim + d];
        }
      }
      double det = 0;
      for (unsigned d = 0; d < dim; d++) {
        det += NN[d] * NN[d];
      }
      det = std::sqrt(det + 1.e-10);
      for (unsigned d = 0; d < dim; d++) {
        NN[d] /= det;
      }

      std::vector<double> N_g(dim, 0.);
      for (unsigned i = 0; i < nDofsN; i++) {
        for(unsigned d = 0; d < dim; d++) {
          N_g[d] += N[d][i] * phiN[i];
        }
      }

      // *** phiV_i loop ***
      for(unsigned i = 0; i < nDofsN; i++) {
        for(unsigned  d = 0; d < dim; d++) {  //momentum equation in k
          double rhs = 0.;
          rhs += phiN[i] * (NN[d] - /*N_g[d]*/N[d][i]);
          Res[d * nDofsN + i] +=  rhs * weight;
        }
      } // end phiV_i loop


      //--------------------------------------------------------------------------------------------------------
      // Add the local Matrix/Vector into the global Matrix/Vector

      for(unsigned i = 0; i < nDofsN; i++) {
        for(unsigned d = 0; d < dim; d++) { //row velocity blocks or dimension
          unsigned VIrow = d * nDofsN + i;
          for(unsigned j = 0; j < nDofsN; j++) {
            unsigned VIcolumn = d * nDofsN + j;

            VIcolumn = VIrow;
            Jac[ VIrow * nDofs + VIcolumn] += phiN[i] * phiN[j] * weight ; // inertia


          }
        }
      }
    }

    RES->add_vector_blocked(Res, sysDof);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);


  } //end element loop for each process

  RES->close();
  KK->close();

  int sol_offset = msh->_dofOffset[solNType][iproc];

  for(unsigned k = 0; k < solNPdeIndex.size(); k++) {

    unsigned indexSol = solNIndex[k];
    int sys_offset = mlPdeSys->_LinSolver[level]->KKoffset[k][iproc];

    for(int i = 0; i < msh->_ownSize[solNType][iproc]; i++) {

      int sol_dof = sol_offset + i;
      int sys_dof = sys_offset + i;

      if (fabs((*KK)(sys_dof, sys_dof)) > 1.e-20)
        sol->_Sol[indexSol]->set( sol_dof, (*RES)(sys_dof) / (*KK)(sys_dof, sys_dof));
    }

    sol->_Sol[indexSol]->close();
  }

}

void AssembleLevelSet(MultiLevelProblem& ml_prob) {

  LinearImplicitSystem* mlPdeSysRestricted =
      &ml_prob.get_system<LinearImplicitSystem>("PsiAux");

  const unsigned levelRestricted =
      mlPdeSysRestricted->GetLevelToAssemble();

  MultiphaseParams mParam =
      ml_prob.GetMultiphaseParams();

  MultiLevelProblem* mlProbFine =
      mParam.mlProbF;

  const unsigned levelF =
      mParam.levelF;

  const unsigned level0 =
      mParam.level0;

  const unsigned targetLevel =
      levelRestricted + level0;

  LinearImplicitSystem* mlPdeSys =
      &mlProbFine->get_system<LinearImplicitSystem>("PsiAux");

  Mesh* msh =
      mlProbFine->_ml_msh->GetLevel(levelF);

  MultiLevelSolution* mlSol =
      mlProbFine->_ml_sol;

  Solution* sol =
      mlSol->GetSolutionLevel(levelF);

  LinearEquationSolver* pdeSys =
      mlPdeSys->_LinSolver[levelF];

  SparseMatrix* KK =
      pdeSys->_KK;

  NumericVector* RES =
      pdeSys->_RES;

  MatSetOption(
      (static_cast<PetscMatrix*>(KK))->mat(),
      MAT_NEW_NONZERO_ALLOCATION_ERR,
      PETSC_FALSE);

  const unsigned dim =
      msh->GetDimension();

  const unsigned iproc =
      msh->processor_id();

  const unsigned psiIndex =
      mlSol->GetIndex("Psi");

  const unsigned psiAuxIndex =
      mlSol->GetIndex("PsiAux");

  const unsigned psiType =
      mlSol->GetSolutionType(psiIndex);

  const unsigned psiAuxType =
      mlSol->GetSolutionType(psiAuxIndex);

  const unsigned psiAuxPdeIndex =
      mlPdeSys->GetSolPdeIndex("PsiAux");

  const unsigned solXType = 2;

  std::vector<std::vector<double>> coordX(dim);

  std::vector<double> psi;
  std::vector<double> psiAux;

  std::vector<double> phi;
  std::vector<double> phi_x;

  std::vector<double> phiPsi;
  std::vector<double> phiPsi_x;

  std::vector<unsigned> sysDof;
  std::vector<double> Res;
  std::vector<double> Jac;

  KK->zero();
  RES->zero();

  for(unsigned iel = msh->_elementOffset[iproc];
      iel < msh->_elementOffset[iproc + 1];
      ++iel) {

    const unsigned ielGeom =
        msh->GetElementType(iel);

    const unsigned nDofs =
        msh->GetElementDofNumber(iel, psiAuxType);

    const unsigned nDofsPsi =
        msh->GetElementDofNumber(iel, psiType);

    const unsigned nDofsX =
        msh->GetElementDofNumber(iel, solXType);

    psi.resize(nDofsPsi);
    psiAux.resize(nDofs);

    sysDof.resize(nDofs);
    Res.assign(nDofs, 0.);
    Jac.assign(nDofs * nDofs, 0.);

    for(unsigned d = 0; d < dim; ++d)
      coordX[d].resize(nDofsX);

    for(unsigned i = 0; i < nDofsPsi; ++i) {

      const unsigned dof =
          msh->GetSolutionDof(i, iel, psiType);

      psi[i] =
          (*sol->_Sol[psiIndex])(dof);
    }

    for(unsigned i = 0; i < nDofs; ++i) {

      const unsigned dof =
          msh->GetSolutionDof(i, iel, psiAuxType);

      psiAux[i] =
          (*sol->_Sol[psiAuxIndex])(dof);

      sysDof[i] =
          pdeSys->GetSystemDof(
              psiAuxIndex,
              psiAuxPdeIndex,
              i,
              iel);
    }

    for(unsigned i = 0; i < nDofsX; ++i) {

      const unsigned dof =
          msh->GetSolutionDof(i, iel, solXType);

      for(unsigned d = 0; d < dim; ++d)
        coordX[d][i] =
            (*msh->_topology->_Sol[d])(dof);
    }

    const elem_type* fem =
        msh->_finiteElement[ielGeom][psiAuxType];

    const elem_type* femPsi =
        msh->_finiteElement[ielGeom][psiType];

    for(unsigned ig = 0;
        ig < fem->GetGaussPointNumber();
        ++ig) {

      double weight = 0.;
      double weightPsi = 0.;

      fem->Jacobian(
          coordX,
          ig,
          weight,
          phi,
          phi_x);

      femPsi->Jacobian(
          coordX,
          ig,
          weightPsi,
          phiPsi,
          phiPsi_x);

      double psi_g = 0.;
      double psiAux_g = 0.;

      for(unsigned j = 0; j < nDofsPsi; ++j)
        psi_g += psi[j] * phiPsi[j];

      for(unsigned j = 0; j < nDofs; ++j)
        psiAux_g += psiAux[j] * phi[j];

      for(unsigned i = 0; i < nDofs; ++i) {

        Res[i] +=
            phi[i] *
            (psi_g - psiAux_g) *
            weight;

        for(unsigned j = 0; j < nDofs; ++j)
          Jac[i * nDofs + j] +=
              phi[i] *
              phi[j] *
              weight;
      }
    }

    RES->add_vector_blocked(
        Res,
        sysDof);

    KK->add_matrix_blocked(
        Jac,
        sysDof,
        sysDof);
  }

  RES->close();
  KK->close();

  RestrictFineSystem(
      *mlPdeSys,
      *mlProbFine->_ml_msh,
      levelF,
      targetLevel);

  CopyRestrictedSystem(
      *mlPdeSys,
      targetLevel,
      *mlPdeSysRestricted->_LinSolver[levelRestricted]);
}