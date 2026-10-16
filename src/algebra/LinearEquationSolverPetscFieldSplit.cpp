/*=========================================================================

  Program: FEMUS
  Module: LinearEquationSolverPetscFieldSplit
  Authors: Eugenio Aulisa, Guoyi Ke

  Copyright (c) FEMTTU
  All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

  =========================================================================*/

#include "FemusConfig.hpp"

#ifdef HAVE_PETSC

#include "LinearEquationSolverPetscFieldSplit.hpp"

namespace femus {

  void LinearEquationSolverPetscFieldSplit::SetFieldSplitTree(FieldSplitTree* fieldSplitTree) {
    _fieldSplitTree = fieldSplitTree;
  }

  void LinearEquationSolverPetscFieldSplit::BuildBdcIndex(const vector <unsigned>& variable_to_be_solved) {
    if(_fieldSplitTree != NULL) _fieldSplitTree->BuildIndexSet(KKoffset, _iproc, _nprocs, _msh->GetLevel(), this);
    else FielSlipTreeIsNotDefined();
    LinearEquationSolverPetsc::BuildBdcIndex(variable_to_be_solved);
  }

  void LinearEquationSolverPetscFieldSplit::SetPreconditioner(KSP& subksp, PC& subpc) {
    if(_fieldSplitTree != NULL) _fieldSplitTree->SetPC(subksp, _msh->GetLevel());
    else FielSlipTreeIsNotDefined();
  }

  void LinearEquationSolverPetscFieldSplit::FielSlipTreeIsNotDefined() {
    std::cout << "Error! No FieldSplitTree object has been passed to the FEMuS_FIELDSPLIT system" << std::endl;
    std::cout << "Define a FieldSplitTree object FS and pass it to the FEMuS_FIELDSPLIT system with" << std::endl;
    std::cout << "system.SetFieldSplitTree(&FS);" << std::endl;
    abort();
  }

  void LinearEquationSolverPetscFieldSplit::MGSolve (const bool ksp_clean) {

    PetscLogDouble t1;
    PetscLogDouble t2;
    PetscTime (&t1);
    if (ksp_clean) {
      Mat KK = (static_cast< PetscMatrix* > (_KK))->mat();
      KSPSetOperators (_ksp, KK, KK);
      KSPSetTolerances (_ksp, _rtol, _abstol, _dtol, _maxits);

      if (_mgSolverType != PREONLY) {
        KSPSetInitialGuessKnoll (_ksp, PETSC_TRUE);
      }
      else {
        KSPSetInitialGuessKnoll (_ksp, PETSC_FALSE);
        KSPSetNormType (_ksp, KSP_NORM_NONE);
      }
      KSPSetFromOptions (_ksp);
      KSPSetNormType(_ksp, KSP_NORM_UNPRECONDITIONED);
      KSPGMRESSetRestart (_ksp, _restart);
      KSPOrthogonalizationSetCGSRefinementType(_ksp, KSP_ORTHOGONALIZATION_CGS_REFINE_IFNEEDED);
      KSPSetUp (_ksp);

      //BEGIN setup for the multigrid inside the field-split

      unsigned levelMax = _msh->GetLevel() + 1;

      KSP* kspMG = this->GetKSP();
      PC pcMG;
      KSPGetPC (*kspMG, &pcMG);

      for(unsigned level = 0; level < levelMax; level++ ) {
        KSP kspLevel;
        if (level == 0) {
          PCMGGetCoarseSolve (pcMG, &kspLevel);
        }
        else {
          PCMGGetSmoother (pcMG, level, &kspLevel);
        }

        _fieldSplitTree->MGSolveSetup(kspLevel, level);
      }

      //END setup for the multigrid inside the field-split

    }

    ZerosBoundaryResiduals();

    KSPSolve (_ksp, (static_cast< PetscVector* > (_RES))->vec(), (static_cast< PetscVector* > (_EPSC))->vec());

    _RESC->matrix_mult (*_EPSC, *_KK);
    *_RES -= *_RESC;
    *_EPS += *_EPSC;

    if (_printSolverInfo) {
      int its;
      KSPGetIterationNumber (_ksp, &its);

      KSPConvergedReason reason;
      KSPGetConvergedReason (_ksp, &reason);

      PetscReal rnorm;
      KSPGetResidualNorm (_ksp, &rnorm);

      PetscTime (&t2);
      PetscPrintf (PETSC_COMM_WORLD, "       *************** MG linear solver time: %e \n", t2 - t1);
      PetscPrintf (PETSC_COMM_WORLD, "       *************** Number of outer ksp solver iterations = %i \n", its);
      PetscPrintf (PETSC_COMM_WORLD, "       *************** Convergence reason = %i \n", reason);
      PetscPrintf (PETSC_COMM_WORLD, "       *************** Residual norm = %10.8g \n", rnorm);
    }
  }

  void LinearEquationSolverPetscFieldSplit::MGSetLevel (LinearEquationSolver *LinSolver, const unsigned &maxlevel,
      const vector <unsigned> &variable_to_be_solved,
      SparseMatrix* PP, SparseMatrix* RR,
      const unsigned &npre, const unsigned &npost) {
    LinearEquationSolverPetsc::MGSetLevel (LinSolver, maxlevel, variable_to_be_solved, PP, RR, npre, npost);

    KSP* kspMG = LinSolver->GetKSP();
    PC pcMG;
    KSPGetPC (*kspMG, &pcMG);

    KSP kspLevel;
    unsigned level = _msh->GetLevel();

    if (level == 0) {
      PCMGGetCoarseSolve (pcMG, &kspLevel);
    }
    else {
      PCMGGetSmoother (pcMG, level, &kspLevel);
    }

    Mat P = (level != 0) ? (static_cast<PetscMatrix*>(PP))->mat() : nullptr;
    Mat R = (level != 0) ? ((PP == RR) ? P : (static_cast<PetscMatrix*>(RR))->mat()) : nullptr;
    _fieldSplitTree->MGSetLevel (kspLevel, level, P, R);

    // if (level == 0) {
    //   PCMGGetCoarseSolve (pcMG, &kspLevel);
    //   KSPSetTolerances (kspLevel, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, npre);
    //   this->SetSolver (kspLevel, _levelSolverType);
    // }

  }

} //end namespace femus

#endif
