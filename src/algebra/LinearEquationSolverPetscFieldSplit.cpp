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
        KSPSetInitialGuessKnoll (_ksp, PETSC_FALSE);
      }
      else {
        KSPSetInitialGuessKnoll (_ksp, PETSC_FALSE);
        //KSPSetNormType (_ksp, KSP_NORM_NONE);
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
        //PetscCallAbort(PETSC_COMM_WORLD, KSPSetInitialGuessKnoll(kspLevel, PETSC_FALSE));
        //PetscCallAbort(PETSC_COMM_WORLD, KSPSetInitialGuessNonzero(kspLevel, PETSC_FALSE));

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

  void LinearEquationSolverPetscFieldSplit::Solve(const bool ksp_clean) {

    PetscLogDouble t1;
    PetscLogDouble t2;
    PetscCallAbort(PETSC_COMM_WORLD, PetscTime(&t1));

    if (ksp_clean) {
      Mat KK = (static_cast<PetscMatrix*>(_KK))->mat();
      PetscCallAbort(PETSC_COMM_WORLD, KSPSetOperators(_ksp, KK, KK));
      PetscCallAbort(PETSC_COMM_WORLD, KSPSetTolerances(_ksp, _rtol, _abstol, _dtol, _maxits));
      PetscCallAbort(PETSC_COMM_WORLD, KSPSetInitialGuessKnoll(_ksp, PETSC_FALSE));
      PetscCallAbort(PETSC_COMM_WORLD, KSPSetFromOptions(_ksp));
      PetscCallAbort(PETSC_COMM_WORLD, KSPSetNormType(_ksp, KSP_NORM_UNPRECONDITIONED));
      PetscCallAbort(PETSC_COMM_WORLD, KSPGMRESSetRestart(_ksp, _restart));
      PetscCallAbort(PETSC_COMM_WORLD, KSPOrthogonalizationSetCGSRefinementType(_ksp, KSP_ORTHOGONALIZATION_CGS_REFINE_IFNEEDED));
      PetscCallAbort(PETSC_COMM_WORLD, KSPSetUp(_ksp));
    }

    PC pcMG;
    PetscCallAbort(PETSC_COMM_WORLD, KSPGetPC(_ksp, &pcMG));

    PetscInt numberOfLevels;
    PetscCallAbort(PETSC_COMM_WORLD, PCMGGetLevels(pcMG, &numberOfLevels));

    if (ksp_clean) {
      for (unsigned level = 0; level < static_cast<unsigned>(numberOfLevels); ++level) {
        KSP kspLevel;

        if (level == 0) {
          PetscCallAbort(PETSC_COMM_WORLD, PCMGGetCoarseSolve(pcMG, &kspLevel));
        }
        else {
          PetscCallAbort(PETSC_COMM_WORLD, PCMGGetSmoother(pcMG, level, &kspLevel));
        }

        _fieldSplitTree->MGSolveSetup(kspLevel, level);
      }
    }

    KSP kspFine;
    PetscCallAbort(PETSC_COMM_WORLD, PCMGGetSmoother(pcMG, numberOfLevels - 1, &kspFine));

    KSPNormType previousNormType;
    PetscBool previousInitialGuessNonzero;
    PetscBool previousInitialGuessKnoll;
    PetscCallAbort(PETSC_COMM_WORLD, KSPGetNormType(kspFine, &previousNormType));
    PetscCallAbort(PETSC_COMM_WORLD, KSPGetInitialGuessNonzero(kspFine, &previousInitialGuessNonzero));
    PetscCallAbort(PETSC_COMM_WORLD, KSPGetInitialGuessKnoll(kspFine, &previousInitialGuessKnoll));

    PetscErrorCode (*previousConvergenceTest)(KSP, PetscInt, PetscReal, KSPConvergedReason*, void*) = nullptr;
    PetscErrorCode (*previousConvergenceDestroy)(void*) = nullptr;
    void* previousConvergenceContext = nullptr;
    PetscCallAbort(PETSC_COMM_WORLD, KSPGetAndClearConvergenceTest(kspFine, &previousConvergenceTest, &previousConvergenceContext, &previousConvergenceDestroy));

    void* solveConvergenceContext = nullptr;
    PetscCallAbort(PETSC_COMM_WORLD, KSPConvergedDefaultCreate(&solveConvergenceContext));
    PetscCallAbort(PETSC_COMM_WORLD, KSPSetConvergenceTest(kspFine, KSPConvergedDefault, solveConvergenceContext, KSPConvergedDefaultDestroy));
    PetscCallAbort(PETSC_COMM_WORLD, KSPSetNormType(kspFine, KSP_NORM_UNPRECONDITIONED));
    PetscCallAbort(PETSC_COMM_WORLD, KSPOrthogonalizationSetCGSRefinementType(kspFine, KSP_ORTHOGONALIZATION_CGS_REFINE_IFNEEDED));

    PetscCallAbort(PETSC_COMM_WORLD, KSPSetInitialGuessKnoll(kspFine, PETSC_TRUE));
    PetscCallAbort(PETSC_COMM_WORLD, KSPSetInitialGuessNonzero(kspFine, PETSC_FALSE));

    ZerosBoundaryResiduals();

    PetscCallAbort(PETSC_COMM_WORLD, KSPSolve(kspFine, (static_cast<PetscVector*>(_RES))->vec(), (static_cast<PetscVector*>(_EPSC))->vec()));

    PetscCallAbort(PETSC_COMM_WORLD, KSPSetConvergenceTest(kspFine, previousConvergenceTest, previousConvergenceContext, previousConvergenceDestroy));
    //PetscCallAbort(PETSC_COMM_WORLD, KSPSetTolerances(kspFine, previousRtol, previousAtol, previousDtol, previousMaxIts));
    PetscCallAbort(PETSC_COMM_WORLD, KSPSetNormType(kspFine, previousNormType));
    PetscCallAbort(PETSC_COMM_WORLD, KSPSetInitialGuessKnoll(kspFine, previousInitialGuessKnoll));
    PetscCallAbort(PETSC_COMM_WORLD, KSPSetInitialGuessNonzero(kspFine, previousInitialGuessNonzero));

    _RESC->matrix_mult(*_EPSC, *_KK);
    *_RES -= *_RESC;
    *_EPS += *_EPSC;

    if (_printSolverInfo) {
      PetscInt its;
      KSPConvergedReason reason;
      PetscReal rnorm;

      PetscCallAbort(PETSC_COMM_WORLD, KSPGetIterationNumber(kspFine, &its));
      PetscCallAbort(PETSC_COMM_WORLD, KSPGetConvergedReason(kspFine, &reason));
      PetscCallAbort(PETSC_COMM_WORLD, KSPGetResidualNorm(kspFine, &rnorm));
      PetscCallAbort(PETSC_COMM_WORLD, PetscTime(&t2));

      PetscCallAbort(PETSC_COMM_WORLD, PetscPrintf(PETSC_COMM_WORLD, "       *************** Fine-level linear solver time: %e\n", static_cast<double>(t2 - t1)));
      PetscCallAbort(PETSC_COMM_WORLD, PetscPrintf(PETSC_COMM_WORLD, "       *************** Number of fine-level KSP iterations = %" PetscInt_FMT "\n", its));
      PetscCallAbort(PETSC_COMM_WORLD, PetscPrintf(PETSC_COMM_WORLD, "       *************** Convergence reason = %i\n", static_cast<int>(reason)));
      PetscCallAbort(PETSC_COMM_WORLD, PetscPrintf(PETSC_COMM_WORLD, "       *************** Residual norm = %10.8g\n", static_cast<double>(rnorm)));
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

  }

} //end namespace femus

#endif
