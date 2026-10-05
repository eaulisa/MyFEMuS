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

  void LinearEquationSolverPetscFieldSplit::MGSetLevel (LinearEquationSolver *LinSolver, const unsigned &maxlevel,
      const vector <unsigned> &variable_to_be_solved,
      SparseMatrix* PP, SparseMatrix* RR,
      const unsigned &npre, const unsigned &npost) {
    LinearEquationSolverPetsc::MGSetLevel (LinSolver, maxlevel, variable_to_be_solved, PP, RR, npre, npost);

    KSP* kspMG = LinSolver->GetKSP();
    PC pcMG;
    KSPGetPC (*kspMG, &pcMG);

    KSP subksp;
    unsigned level = _msh->GetLevel();

    if (level == 0) {
      PCMGGetCoarseSolve (pcMG, &subksp);
    }
    else {
      PCMGGetSmoother (pcMG, level, &subksp);
    }

    std::cout << "AAAAAAAA " << level << std::endl;

    PC subpc;
    KSPGetPC (subksp, &subpc);

    PetscBool isFieldSplit = PETSC_FALSE;

    PetscObjectTypeCompare((PetscObject) subpc, PCFIELDSPLIT, &isFieldSplit);

    if (isFieldSplit) {

      PetscInt nsplit = 0;
      KSP* subksps = NULL;

      PCSetUp(subpc);

      PCFieldSplitGetSubKSP(subpc, &nsplit, &subksps);

      std::vector<std::vector <Mat>> Ks = _fieldSplitTree->GetKSplit();
      std::vector<std::vector <Mat>> Ps = _fieldSplitTree->GetPSplit();
      std::vector<std::vector <Mat>> Is = _fieldSplitTree->GetISplit();
      std::vector < std::vector < IS > > IS = _fieldSplitTree->GetISSplit();

      Mat I = static_cast<PetscMatrix*>(PP)->mat();

      if(level + 1 > Ks.size()) {
        Ks.resize(level + 1);
        Ps.resize(level + 1);
        Is.resize(level + 1);
      }

      Ks[level].resize(nsplit);
      Ps[level].resize(nsplit);
      if(level > 0) {
        Is[level].resize(nsplit);
        for (unsigned i = 0; i < Is[level].size(); ++i) {
          Is[level][i] = NULL;
        }
      }

      for (unsigned i = 0; i < nsplit; ++i) {

        std::cout << "CCCCCCCCC" << i << std::endl;

        PC  pc_split = NULL;

        KSPGetOperators(subksps[i], &Ks[level][i], &Ps[level][i]);
        KSPGetPC(subksps[i], &pc_split);

        // A_split  = operator solved by this split
        // P_split  = matrix used to build its preconditioner
        // pc_split = actual PC associated with this split

        if(level > 0) {
          MatCreateSubMatrix(I, IS[level][i], IS[level - 1][i], MAT_INITIAL_MATRIX, &Is[level][i]);
        }
      }

      std::cout << "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB\n";
    }

  }

} //end namespace femus

#endif
