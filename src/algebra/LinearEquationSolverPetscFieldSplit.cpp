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

        PC pcLevel;
        KSPGetPC (kspLevel, &pcLevel);

        PetscBool isFieldSplit = PETSC_FALSE;

        PetscObjectTypeCompare((PetscObject) pcLevel, PCFIELDSPLIT, &isFieldSplit);

        if (isFieldSplit) {
          PetscInt nsplit = 0;
          KSP* kspLevelSplit = NULL;
          PCSetUp(pcLevel);
          PCFieldSplitGetSubKSP(pcLevel, &nsplit, &kspLevelSplit);// one for each split

          if(nsplit > 1) {
            for (unsigned i = 0; i < nsplit; ++i) {
              if(_fieldSplitTree->GetChild(i)->GetPreconditioner() == MG_PRECOND ||
                  _fieldSplitTree->GetChild(i)->GetPreconditioner() == MULTIGRID_PRECOND) {

                PC pcMGSplit;
                KSPGetPC(kspLevelSplit[i], &pcMGSplit);

                PetscBool isGMG = PETSC_FALSE;
                PetscBool isAMG = PETSC_FALSE;

                PetscObjectTypeCompare((PetscObject)pcMGSplit, PCMG, &isGMG);
                PetscObjectTypeCompare((PetscObject)pcMGSplit, PCHMG, &isAMG);

                if(_fieldSplitTree->GetChild(i)->GetPreconditioner() == MULTIGRID_PRECOND) {

                  if(isAMG) {
                    continue;
                  }

                  if(!isGMG) {
                    throw std::runtime_error(
                      "MULTIGRID_PRECOND requires pc_type hmg or mg");
                  }
                }

                if(_fieldSplitTree->GetChild(i)->GetPreconditioner() == MG_PRECOND && !isGMG) {
                  throw std::runtime_error(
                    "MG_PRECOND requires PCMG");
                }

                PCMGSetLevels (pcMGSplit, level + 1, NULL);
                PCMGSetType (pcMGSplit, PC_MG_MULTIPLICATIVE);

                for(unsigned l = 0; l <= level; l++) {

                  KSP subksp;
                  int npre = 1, npost = 1;

                  if (l == 0) {
                    PCMGGetCoarseSolve(pcMGSplit, &subksp);
                    KSPSetType(subksp, KSPPREONLY);
                    KSPSetTolerances(subksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1);
                  }
                  else {
                    PCMGGetSmoother(pcMGSplit, l, &subksp);
                    //DEFAULT
                    KSPSetType(subksp, KSPCHEBYSHEV); // set in options file
                    // KSPSetType(subksp, KSPRICHARDSON);
                    // bool selfScale = false;  // true = automatic, false = fixed
                    // double scale = 0.1;
                    //
                    // if (selfScale) {
                    //   KSPRichardsonSetSelfScale(subksp, PETSC_TRUE);
                    // }
                    // else {
                    //   KSPRichardsonSetSelfScale(subksp, PETSC_FALSE);
                    //   KSPRichardsonSetScale(subksp, scale);
                    // }

                    KSPSetTolerances(subksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 2);
                    KSPSetNormType(subksp, KSP_NORM_NONE);
                  }

                  KSPSetOperators (subksp, _fieldSplitTree->GetChild(i)->GetKgmg()[l], _fieldSplitTree->GetChild(i)->GetKgmg()[l]);
                  PCMGSetGalerkin(pcMGSplit, PC_MG_GALERKIN_NONE);

                  if (l > 0) {
                    PCMGSetInterpolation (pcMGSplit, l, _fieldSplitTree->GetChild(i)->GetPgmg()[l]);
                    PCMGSetRestriction (pcMGSplit, l, _fieldSplitTree->GetChild(i)->GetRgmg()[l]);
                  }

                  PC subpc;
                  KSPGetPC (subksp, &subpc);

                  int parallelOverlapping = 0;

                  if (l == 0) {
                    PCSetType(subpc, PCBJACOBI);
                    //PCSetType(subpc, PCLU);
                    //PCFactorSetMatSolverType(subpc, MATSOLVERMUMPS);
                    // PCSetType(subpc, PCHMG);
                    // PCHMGSetInnerPCType(subpc, PCGAMG);
                    // PCHMGSetUseSubspaceCoarsening(subpc, PETSC_FALSE);
                    // PCHMGUseMatMAIJ(subpc, PETSC_FALSE);
                  }
                  else {
                    //DEFAULT
                    PCSetType(subpc, PCBJACOBI); // set in options file
                  }

                }

                KSPSetFromOptions(kspLevelSplit[i]);
                KSPSetUp(kspLevelSplit[i]);

              }

            }

          }

          PetscFree(kspLevelSplit);
        }

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

    PC pcLevel;
    KSPGetPC (kspLevel, &pcLevel);

    PetscBool isFieldSplit = PETSC_FALSE;

    PetscObjectTypeCompare((PetscObject) pcLevel, PCFIELDSPLIT, &isFieldSplit);

    if (isFieldSplit) {
      PetscInt nsplit = 0;
      KSP* kspLevelSplit = NULL;
      PCSetUp(pcLevel);
      PCFieldSplitGetSubKSP(pcLevel, &nsplit, &kspLevelSplit);

      if(nsplit > 1) {
        for (unsigned i = 0; i < nsplit; ++i) {
          if(_fieldSplitTree->GetChild(i)->GetPreconditioner() == MG_PRECOND ||
              _fieldSplitTree->GetChild(i)->GetPreconditioner() == MULTIGRID_PRECOND) {

            std::vector <Mat> &Kgmg = _fieldSplitTree->GetChild(i)->GetKgmg();

            if (level + 1 > Kgmg.size())   Kgmg.resize(level + 1, NULL);

            Mat Pmat;
            KSPGetOperators(kspLevelSplit[i], &Kgmg[level], &Pmat);

            if(level > 0) {
              std::vector <Mat> &Pgmg = _fieldSplitTree->GetChild(i)->GetPgmg();
              std::vector <Mat> &Rgmg = _fieldSplitTree->GetChild(i)->GetRgmg();
              if (level + 1 > Pgmg.size())   Pgmg.resize(level + 1, NULL);
              if (level + 1 > Rgmg.size())   Rgmg.resize(level + 1, NULL);
              if (Pgmg[level] != NULL) {
                MatDestroy(&Pgmg[level]);
              }
              Pgmg[level] = NULL;
              if (Rgmg[level] != NULL) {
                MatDestroy(&Rgmg[level]);
              }
              Rgmg[level] = NULL;

              std::vector < std::vector < IS > > &isSplit = _fieldSplitTree->GetISSplit();
              Mat P = static_cast<PetscMatrix*>(PP)->mat();
              Mat R = static_cast<PetscMatrix*>(RR)->mat();
              bool sameObject = (PP == RR);

              MatCreateSubMatrix(P, isSplit[level][i], isSplit[level - 1][i], MAT_INITIAL_MATRIX, &Pgmg[level]);
              if(sameObject) {
                Rgmg[level] = Pgmg[level];
                PetscObjectReference((PetscObject)Rgmg[level]);
              }
              else {
                MatCreateSubMatrix(R, isSplit[level - 1][i], isSplit[level][i], MAT_INITIAL_MATRIX, &Rgmg[level]);
              }
            }
          }
        }
      }
      PetscFree(kspLevelSplit);
    }
  }

} //end namespace femus

#endif
