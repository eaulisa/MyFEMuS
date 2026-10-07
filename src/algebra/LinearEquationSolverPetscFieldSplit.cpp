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

  // void LinearEquationSolverPetscFieldSplit::MGSolve (const bool ksp_clean) {
  //
  //   unsigned levelMax = _msh->GetLevel() + 1;
  //
  //   KSP* kspMG = this->GetKSP();
  //   PC pcMG;
  //   KSPGetPC (*kspMG, &pcMG);
  //
  //   for(unsigned level = 0; level < levelMax; level++ ) {
  //
  //     KSP kspLevel;
  //     if (level == 0) {
  //       PCMGGetCoarseSolve (pcMG, &kspLevel);
  //     }
  //     else {
  //       PCMGGetSmoother (pcMG, level, &kspLevel);
  //     }
  //
  //     PC pcLevel;
  //     KSPGetPC (kspLevel, &pcLevel);
  //
  //     PetscBool isFieldSplit = PETSC_FALSE;
  //
  //     PetscObjectTypeCompare((PetscObject) pcLevel, PCFIELDSPLIT, &isFieldSplit);
  //
  //     if (isFieldSplit) {
  //       PetscInt nsplit = 0;
  //       KSP* kspLevelSplit = NULL;
  //       PCSetUp(pcLevel);
  //       PCFieldSplitGetSubKSP(pcLevel, &nsplit, &kspLevelSplit);// one for each split
  //
  //       if(nsplit > 1) {
  //         for (unsigned i = 0; i < nsplit; ++i) {
  //           if(_fieldSplitTree->GetChild(i)->GetPreconditioner() == AMG_PRECOND) {
  //             std::cout << "BBBBBBBBBB " << level << " " << i << "\n";
  //
  //             PC pcMGSplit;
  //             KSPGetPC (kspLevelSplit[i], &pcMGSplit); // The multigrid preconditioner of the FS[level][i]
  //
  //             PCSetType (pcMGSplit, PCMG);
  //             PCMGSetLevels (pcMGSplit, level + 1, NULL);
  //             PCMGSetType (pcMGSplit, PC_MG_MULTIPLICATIVE);
  //
  //             for(unsigned l = 0; l <= level; l++) {
  //               KSP subksp;
  //               int npre = 1, npost = 1;
  //               if (l == 0) {
  //                 PCMGGetCoarseSolve (pcMGSplit, &subksp);
  //                 KSPSetTolerances (subksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, npre);
  //               }
  //               else {
  //                 PCMGGetSmoother (pcMGSplit, l, &subksp);
  //                 KSPSetTolerances (subksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, npre);
  //               }
  //
  //               SolverType levelSolverType = PREONLY;
  //               double richardsonScaleFactor = 1.;
  //               if (l != 0 && levelSolverType == PREONLY) {
  //                 levelSolverType = RICHARDSON;
  //               }
  //               SetPetscSolverType (subksp, levelSolverType, &richardsonScaleFactor);
  //               KSPSetFromOptions (subksp);
  //
  //               std::cout << "CCCCCCCCCC " << l << " " << i << std::endl << std::flush;
  //
  //               KSPSetOperators (subksp, _fieldSplitTree->GetChild(i)->GetKgmg()[l], _fieldSplitTree->GetChild(i)->GetKgmg()[l]);
  //
  //
  //               std::cout << "CCCC1 " << l << " " << i << std::endl << std::flush;
  //               PC subpc;
  //               KSPGetPC (subksp, &subpc);
  //               std::cout << "CCCC2 " << l << " " << i << std::endl << std::flush;
  //
  //               int parallelOverlapping = 0;
  //               if (l == 0 && i == 0) PetscPreconditioner::set_petsc_preconditioner_type (MLU_PRECOND, subpc, parallelOverlapping);
  //               else PetscPreconditioner::set_petsc_preconditioner_type (ILU_PRECOND, subpc, parallelOverlapping);
  //               PetscReal zero = 1.e-16;
  //               PCFactorSetZeroPivot (subpc, zero);
  //               PCFactorSetShiftType (subpc, MAT_SHIFT_NONZERO);
  //
  //
  //
  //               //SetPreconditioner (subksp, subpc);
  //
  //               std::cout << "CCCC3 " << l << " " << i << std::endl << std::flush;
  //
  //               if (l < level) {
  //                 PCMGSetX (pcMGSplit, l, _fieldSplitTree->GetChild(i)->GetXgmg()[l]);
  //                 std::cout << "X " << l << " " << i << std::endl << std::flush;
  //                 PCMGSetRhs (pcMGSplit, l, _fieldSplitTree->GetChild(i)->GetRHSgmg()[l]);
  //                 std::cout << "Rhs " << l << " " << i << std::endl << std::flush;
  //               }
  //
  //               // KSPSetUp (subksp);
  //               if (l > 0) {
  //                 PCMGSetR (pcMGSplit, l, _fieldSplitTree->GetChild(i)->GetRESgmg()[l]);
  //                 std::cout << "RES " << level << " " << i << std::endl << std::flush;
  //                 PCMGSetInterpolation (pcMGSplit, l, _fieldSplitTree->GetChild(i)->GetPgmg()[l]);
  //                 std::cout << "P " << level << " " << i << std::endl << std::flush;
  //                 PCMGSetRestriction (pcMGSplit, l, _fieldSplitTree->GetChild(i)->GetRgmg()[l]);
  //                 std::cout << "R " << level << " " << i << std::endl << std::flush;
  //
  //                 if (npre != npost) {
  //                   KSP subkspUp;
  //                   PCMGGetSmootherUp (pcMGSplit, l, &subkspUp);
  //                   KSPSetTolerances (subkspUp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, npost);
  //                   this->SetSolver (subkspUp, levelSolverType);
  //                   KSPSetPC (subkspUp, subpc);
  //                   PC subpcUp;
  //                   KSPGetPC (subkspUp, &subpcUp);
  //                   KSPSetUp (subkspUp);
  //                 }
  //               }
  //             }
  //           }
  //         }
  //       }
  //       PetscFree(kspLevelSplit);
  //     }
  //
  //   }
  //
  //
  //   std::cout << "BBBBBBBBBB " << std::endl << std::flush;
  //   LinearEquationSolverPetsc::MGSolve (ksp_clean);
// }

  void LinearEquationSolverPetscFieldSplit::MGSolve (const bool ksp_clean) {

    PetscLogDouble t1;
    PetscLogDouble t2;
    PetscTime (&t1);
    // std::cout << "AAAAAA\n";
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
      std::cout << "AA1\n";
      KSPSetFromOptions (_ksp);
      KSPSetNormType(_ksp, KSP_NORM_UNPRECONDITIONED);

      std::cout << "AA11\n";
      KSPGMRESSetRestart (_ksp, _restart);
      std::cout << "AA12\n";
      KSPOrthogonalizationSetCGSRefinementType(_ksp, KSP_ORTHOGONALIZATION_CGS_REFINE_IFNEEDED);
      std::cout << "AA13\n";
      KSPSetUp (_ksp);
      std::cout << "AA2\n";

      if(true) {

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
              for (unsigned i = 0; i < 1 + 0 * nsplit; ++i) {
                if(_fieldSplitTree->GetChild(i)->GetPreconditioner() == AMG_PRECOND) {
                  std::cout << "BBBBBBBBBB " << level << " " << i << "\n";

                  PC pcMGSplit;
                  KSPGetPC (kspLevelSplit[i], &pcMGSplit); // The multigrid preconditioner of the FS[level][i]

                  PCSetType (pcMGSplit, PCMG);
                  PCMGSetLevels (pcMGSplit, level + 1, NULL);
                  PCMGSetType (pcMGSplit, PC_MG_MULTIPLICATIVE);

                  for(unsigned l = 0; l <= level; l++) {

                    KSP subksp;
                    int npre = 1, npost = 1;

                    if (l == 0) {

                      PCMGGetCoarseSolve (pcMGSplit, &subksp);

                      KSPSetType (subksp, KSPPREONLY);

                      KSPSetTolerances (subksp,
                                        PETSC_DEFAULT,
                                        PETSC_DEFAULT,
                                        PETSC_DEFAULT,
                                        1);
                    }
                    else {

                      npre = 4;
                      npost = 4;

                      PCMGGetSmoother (pcMGSplit, l, &subksp);

                      KSPSetType (subksp, KSPGMRES);

                      KSPSetTolerances (subksp,
                                        PETSC_DEFAULT,
                                        PETSC_DEFAULT,
                                        PETSC_DEFAULT,
                                        4);

                      KSPSetNormType (subksp, KSP_NORM_NONE);
                    }

                    // KSP subksp;
                    // int npre = 1, npost = 1;
                    // if (l == 0) {
                    //   PCMGGetCoarseSolve (pcMGSplit, &subksp);
                    //   KSPSetTolerances (subksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, npre);
                    // }
                    // else {
                    //   PCMGGetSmoother (pcMGSplit, l, &subksp);
                    //   KSPSetTolerances (subksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, npre);
                    // }
                    //
                    // SolverType levelSolverType = PREONLY;
                    // double richardsonScaleFactor = 1.;
                    // if (l != 0 && levelSolverType == PREONLY) {
                    //   levelSolverType = RICHARDSON;
                    // }
                    // SetPetscSolverType (subksp, levelSolverType, &richardsonScaleFactor);
                    // //KSPSetFromOptions (subksp);

                    std::cout << "CCCCCCCCCC " << l << " " << i << std::endl << std::flush;

                    if(l == level) KSPSetOperators (subksp, _fieldSplitTree->GetChild(i)->GetKgmg()[l], _fieldSplitTree->GetChild(i)->GetKgmg()[l]);

                    PCMGSetGalerkin(pcMGSplit, PC_MG_GALERKIN_BOTH);

                    if (l > 0) {
                      // PCMGSetR (pcMGSplit, l, _fieldSplitTree->GetChild(i)->GetRESgmg()[l]);
                      // std::cout << "RES " << level << " " << i << std::endl << std::flush;
                      PCMGSetInterpolation (pcMGSplit, l, _fieldSplitTree->GetChild(i)->GetPgmg()[l]);
                      std::cout << "P " << level << " " << i << std::endl << std::flush;
                      // PCMGSetRestriction (pcMGSplit, l, _fieldSplitTree->GetChild(i)->GetRgmg()[l]);
                      // std::cout << "R " << level << " " << i << std::endl << std::flush;

                      // if (npre != npost) {
                      //   KSP subkspUp;
                      //   PCMGGetSmootherUp (pcMGSplit, l, &subkspUp);
                      //   KSPSetTolerances (subkspUp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, npost);
                      //   this->SetSolver (subkspUp, levelSolverType);
                      //   KSPSetPC (subkspUp, subpc);
                      //   PC subpcUp;
                      //   KSPGetPC (subkspUp, &subpcUp);
                      //   KSPSetUp (subkspUp);
                      // }
                    }

                    std::cout << "CCCC1 " << l << " " << i << std::endl << std::flush;
                    PC subpc;
                    KSPGetPC (subksp, &subpc);
                    std::cout << "CCCC2 " << l << " " << i << std::endl << std::flush;
                    //
                    // int parallelOverlapping = 0;
                    // if (l == 0 && i == 0) PetscPreconditioner::set_petsc_preconditioner_type (MLU_PRECOND, subpc, parallelOverlapping);
                    // else PetscPreconditioner::set_petsc_preconditioner_type (ILU_PRECOND, subpc, parallelOverlapping);
                    // PetscReal zero = 1.e-16;
                    // PCFactorSetZeroPivot (subpc, zero);
                    // PCFactorSetShiftType (subpc, MAT_SHIFT_NONZERO);

                    int parallelOverlapping = 0;

                    if (l == 0) {

                      KSPSetType(subksp, KSPPREONLY);

                      PC subpc;
                      KSPGetPC(subksp, &subpc);

                      PCSetType(subpc, PCLU);
                      PCFactorSetMatSolverType(subpc, MATSOLVERMUMPS);

                      PetscReal zero = 1.e-16;
                      PCFactorSetZeroPivot(subpc, zero);
                      PCFactorSetShiftType(subpc, MAT_SHIFT_NONZERO);

                      // PetscPreconditioner::set_petsc_preconditioner_type (
                      //   MLU_PRECOND,
                      //   subpc,
                      //   parallelOverlapping);
                      //
                      // PetscReal zero = 1.e-16;
                      //
                      // PCFactorSetZeroPivot (subpc, zero);
                      // PCFactorSetShiftType (subpc, MAT_SHIFT_NONZERO);
                    }
                    else {

                      PCSetType (subpc, PCILU);
                    }

                    //SetPreconditioner (subksp, subpc);

                    std::cout << "CCCC3 " << l << " " << i << std::endl << std::flush;

                    // if (l < level) {
                    //   PCMGSetX (pcMGSplit, l, _fieldSplitTree->GetChild(i)->GetXgmg()[l]);
                    //   std::cout << "X " << l << " " << i << std::endl << std::flush;
                    //   PCMGSetRhs (pcMGSplit, l, _fieldSplitTree->GetChild(i)->GetRHSgmg()[l]);
                    //   std::cout << "Rhs " << l << " " << i << std::endl << std::flush;
                    // }

                    // KSPSetUp (subksp);
                    // PCMGSetGalerkin(pcMGSplit, PC_MG_GALERKIN_BOTH);
                    //
                    // if (l > 0) {
                    //   // PCMGSetR (pcMGSplit, l, _fieldSplitTree->GetChild(i)->GetRESgmg()[l]);
                    //   // std::cout << "RES " << level << " " << i << std::endl << std::flush;
                    //   PCMGSetInterpolation (pcMGSplit, l, _fieldSplitTree->GetChild(i)->GetPgmg()[l]);
                    //   std::cout << "P " << level << " " << i << std::endl << std::flush;
                    //   // PCMGSetRestriction (pcMGSplit, l, _fieldSplitTree->GetChild(i)->GetRgmg()[l]);
                    //   // std::cout << "R " << level << " " << i << std::endl << std::flush;
                    //
                    //   // if (npre != npost) {
                    //   //   KSP subkspUp;
                    //   //   PCMGGetSmootherUp (pcMGSplit, l, &subkspUp);
                    //   //   KSPSetTolerances (subkspUp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, npost);
                    //   //   this->SetSolver (subkspUp, levelSolverType);
                    //   //   KSPSetPC (subkspUp, subpc);
                    //   //   PC subpcUp;
                    //   //   KSPGetPC (subkspUp, &subpcUp);
                    //   //   KSPSetUp (subkspUp);
                    //   // }
                    // }
                  }
                }
              }
            }
            PetscFree(kspLevelSplit);
          }

        }

        std::cout << "BBBBBBBBBB " << std::endl << std::flush;

      }

//       PetscViewer    viewer;
//       PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,NULL,0,0,1800,1800,&viewer);
//       PetscObjectSetName((PetscObject)viewer,"FSI matrix");
//       PetscViewerPushFormat(viewer,PETSC_VIEWER_DRAW_LG);
//       MatView(KK,viewer);
//
//       VecView((static_cast< PetscVector* >(_RES))->vec(),viewer);
//       double a;
//       std::cin>>a;

    }

    std::cout << "BBBBB\n";

    ZerosBoundaryResiduals();
    std::cout << "BB1\n";
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
          if(_fieldSplitTree->GetChild(i)->GetPreconditioner() == AMG_PRECOND) {

            std::vector <Mat> &Kgmg = _fieldSplitTree->GetChild(i)->GetKgmg();
            std::vector <Vec> &Xgmg = _fieldSplitTree->GetChild(i)->GetXgmg();
            std::vector <Vec> &RESgmg = _fieldSplitTree->GetChild(i)->GetRESgmg();
            std::vector <Vec> &RHSgmg = _fieldSplitTree->GetChild(i)->GetRHSgmg();

            if (level + 1 > Kgmg.size())   Kgmg.resize(level + 1, NULL);
            if (level + 1 > Xgmg.size())   Xgmg.resize(level + 1, NULL);
            if (level + 1 > RESgmg.size()) RESgmg.resize(level + 1, NULL);
            if (level + 1 > RHSgmg.size()) RHSgmg.resize(level + 1, NULL);

            Mat Pmat;
            KSPGetOperators(kspLevelSplit[i], &Kgmg[level], &Pmat);

            if (Xgmg[level] != NULL)   VecDestroy(&Xgmg[level]);
            if (RESgmg[level] != NULL) VecDestroy(&RESgmg[level]);
            if (RHSgmg[level] != NULL) VecDestroy(&RHSgmg[level]);

            MatCreateVecs(Kgmg[level], &Xgmg[level], NULL);
            VecDuplicate(Xgmg[level], &RESgmg[level]);
            VecDuplicate(Xgmg[level], &RHSgmg[level]);

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
