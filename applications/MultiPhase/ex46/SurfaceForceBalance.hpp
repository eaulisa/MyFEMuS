#ifndef EX46_SURFACE_FORCE_BALANCE_HPP
#define EX46_SURFACE_FORCE_BALANCE_HPP

// Diagnostic for ex46: is the discrete Laplace--Beltrami force in range(B^T)?
// Standalone PETSc matrices are assembled before boundary elimination.
// All arrays use the original Navier--Stokes global DOF numbering.
#include <algorithm>
#include <cmath>
#include <vector>
#include "petscksp.h"

namespace ex46_balance {

class SurfaceForceRepresentability {
public:
  void initialize(PetscInt nLocal, PetscInt nGlobal) {
    PetscCallAbort(PETSC_COMM_WORLD,
        MatCreateAIJ(PETSC_COMM_WORLD, nLocal, nLocal, nGlobal, nGlobal,
                     32, nullptr, 32, nullptr, &gradient_));
    PetscCallAbort(PETSC_COMM_WORLD,
        MatCreateAIJ(PETSC_COMM_WORLD, nLocal, nLocal, nGlobal, nGlobal,
                     64, nullptr, 64, nullptr, &mass_));

    // Diagnostic only: tolerate rows with more neighbors than preallocated.
    PetscCallAbort(PETSC_COMM_WORLD,
        MatSetOption(gradient_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
    PetscCallAbort(PETSC_COMM_WORLD,
        MatSetOption(mass_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));

    PetscCallAbort(PETSC_COMM_WORLD,
        VecCreateMPI(PETSC_COMM_WORLD, nLocal, nGlobal, &force_));
    PetscCallAbort(PETSC_COMM_WORLD, VecSet(force_, 0.0));
  }

  // localJac is the ORIGINAL ex46 Jacobian (pressure columns have sign -B^T).
  // localForce contains +f_sigma (opposite to its contribution to Res).
  // localMass is the component-independent consistent velocity L2 mass matrix.
  // freeV and freeP1/freeP2 are the FEMuS _Bdc > 0.5 flags.
  void addElement(const std::vector<unsigned>& sysDof,
                  unsigned dim, unsigned nV, unsigned nP,
                  const std::vector<char>& freeV,
                  const std::vector<char>& freeP1,
                  const std::vector<char>& freeP2,
                  const std::vector<double>& localForce,
                  const std::vector<double>& localMass,
                  const std::vector<double>& localJac) {

    const unsigned nVP = dim * nV + 2 * nP;

    for (unsigned I = 0; I < dim; ++I) {
      std::vector<PetscInt> rows;
      std::vector<unsigned> localRows;
      std::vector<PetscScalar> load;

      for (unsigned i = 0; i < nV; ++i) {
        if (!freeV[I * nV + i]) continue;
        rows.push_back(static_cast<PetscInt>(sysDof[I * nV + i]));
        localRows.push_back(i);
        load.push_back(localForce[I * nV + i]);
      }

      const PetscInt nr = static_cast<PetscInt>(rows.size());
      if (!nr) continue;

      PetscCallAbort(PETSC_COMM_WORLD,
          VecSetValues(force_, nr, rows.data(), load.data(), ADD_VALUES));

      std::vector<PetscScalar> massBlock(nr * nr);
      for (PetscInt i = 0; i < nr; ++i)
        for (PetscInt j = 0; j < nr; ++j)
          massBlock[i * nr + j] =
              localMass[localRows[i] * nV + localRows[j]];

      PetscCallAbort(PETSC_COMM_WORLD,
          MatSetValues(mass_, nr, rows.data(), nr, rows.data(),
                       massBlock.data(), ADD_VALUES));

      for (unsigned phase = 0; phase < 2; ++phase) {
        const std::vector<char>& freeP = phase == 0 ? freeP1 : freeP2;
        std::vector<PetscInt> cols;
        std::vector<unsigned> localCols;

        for (unsigned j = 0; j < nP; ++j) {
          if (!freeP[j]) continue;
          cols.push_back(static_cast<PetscInt>(
              sysDof[dim * nV + phase * nP + j]));
          localCols.push_back(j);
        }

        const PetscInt nc = static_cast<PetscInt>(cols.size());
        if (!nc) continue;

        std::vector<PetscScalar> gradBlock(nr * nc);
        for (PetscInt i = 0; i < nr; ++i) {
          const unsigned iRow = I * nV + localRows[i];
          for (PetscInt j = 0; j < nc; ++j) {
            const unsigned jCol = dim * nV + phase * nP + localCols[j];
            gradBlock[i * nc + j] = -localJac[iRow * nVP + jCol];
          }
        }

        PetscCallAbort(PETSC_COMM_WORLD,
            MatSetValues(gradient_, nr, rows.data(),
                         nc, cols.data(), gradBlock.data(), ADD_VALUES));
      }
    }
  }

  PetscErrorCode report() {
    const MPI_Comm comm = PETSC_COMM_WORLD;
    PetscCallAbort(comm, MatAssemblyBegin(gradient_, MAT_FINAL_ASSEMBLY));
    PetscCallAbort(comm, MatAssemblyEnd(gradient_, MAT_FINAL_ASSEMBLY));
    PetscCallAbort(comm, MatAssemblyBegin(mass_, MAT_FINAL_ASSEMBLY));
    PetscCallAbort(comm, MatAssemblyEnd(mass_, MAT_FINAL_ASSEMBLY));
    PetscCallAbort(comm, VecAssemblyBegin(force_));
    PetscCallAbort(comm, VecAssemblyEnd(force_));

    // Embed the velocity mass in the full NS layout.
    // Pressure and constrained/inactive velocity rows have zero mass;
    // put identity there. Their force and gradient rows are zero.
    // This makes the full matrix SPD without altering any active row.
    Vec diagonal = nullptr;
    PetscCallAbort(comm, VecDuplicate(force_, &diagonal));
    PetscCallAbort(comm, MatGetDiagonal(mass_, diagonal));
    PetscInt nLocal = 0;
    PetscCallAbort(comm, VecGetLocalSize(diagonal, &nLocal));
    PetscScalar* diagArray = nullptr;
    PetscCallAbort(comm, VecGetArray(diagonal, &diagArray));
    for (PetscInt i = 0; i < nLocal; ++i)
      diagArray[i] = PetscAbsScalar(diagArray[i]) < 1.e-30 ? 1.0 : 0.0;
    PetscCallAbort(comm, VecRestoreArray(diagonal, &diagArray));
    PetscCallAbort(comm, MatDiagonalSet(mass_, diagonal, ADD_VALUES));
    PetscCallAbort(comm, VecDestroy(&diagonal));
    PetscCallAbort(comm, MatAssemblyBegin(mass_, MAT_FINAL_ASSEMBLY));
    PetscCallAbort(comm, MatAssemblyEnd(mass_, MAT_FINAL_ASSEMBLY));

    PetscReal normF = 0.0;
    PetscCallAbort(comm, VecNorm(force_, NORM_2, &normF));
    if (normF == 0.0) {
      PetscCallAbort(comm, PetscPrintf(comm,
          "[ex46 balance] Zero surface force; projection skipped.\n"));
      return PETSC_SUCCESS;
    }

    // Consistent mass solve used for all M^{-1} actions.
    // KSP massSolver = nullptr;
    // PetscCallAbort(comm, KSPCreate(comm, &massSolver));
    // PetscCallAbort(comm, KSPSetOptionsPrefix(massSolver, "ex46_mass_"));
    // PetscCallAbort(comm, KSPSetOperators(massSolver, mass_, mass_));
    // PetscCallAbort(comm, KSPSetType(massSolver, KSPCG));
    // PC massPc = nullptr;
    // PetscCallAbort(comm, KSPGetPC(massSolver, &massPc));
    // PetscCallAbort(comm, PCSetType(massPc, PCJACOBI));
    // PetscCallAbort(comm, KSPSetTolerances(
    //     massSolver, 1.e-12, 1.e-14, PETSC_DEFAULT, 2000));
    // PetscCallAbort(comm, KSPSetFromOptions(massSolver));


    // Fixed, accurate M^{-1} operator.
    // Direct factorization is important for the outer CG solver.
    KSP massSolver = nullptr;

    PetscCall(KSPCreate(comm, &massSolver));
    PetscCall(KSPSetOptionsPrefix(
        massSolver, "ex46_mass_"));

    PetscCall(KSPSetOperators(
        massSolver, mass_, mass_));

    PetscCall(KSPSetType(massSolver, KSPPREONLY));

    PC massPc = nullptr;
    PetscCall(KSPGetPC(massSolver, &massPc));
    PetscCall(PCSetType(massPc, PCLU));

    PetscCall(PCFactorSetMatSolverType(
        massPc, MATSOLVERMUMPS));

    PetscCall(KSPSetErrorIfNotConverged(
        massSolver, PETSC_TRUE));

    PetscCall(KSPSetFromOptions(massSolver));
    PetscCall(KSPSetUp(massSolver));


    Vec pressureE = nullptr, projected = nullptr, residualE = nullptr;
    Vec pressureM = nullptr, residualM = nullptr;
    Vec tmp = nullptr, rhs = nullptr, orth = nullptr;
    PetscCallAbort(comm, MatCreateVecs(gradient_, &pressureE, &projected));
    PetscCallAbort(comm, VecDuplicate(projected, &residualE));
    PetscCallAbort(comm, VecDuplicate(pressureE, &pressureM));
    PetscCallAbort(comm, VecDuplicate(projected, &residualM));
    PetscCallAbort(comm, VecDuplicate(projected, &tmp));
    PetscCallAbort(comm, VecDuplicate(pressureE, &rhs));
    PetscCallAbort(comm, VecDuplicate(pressureE, &orth));

    // First projection: Euclidean least squares, LSQR on G p = f.
    KSP lsqr = nullptr;
    PetscCallAbort(comm, KSPCreate(comm, &lsqr));
    PetscCallAbort(comm, KSPSetOptionsPrefix(lsqr, "ex46_euclid_"));
    PetscCallAbort(comm, KSPSetOperators(lsqr, gradient_, gradient_));
    PetscCallAbort(comm, KSPSetType(lsqr, KSPLSQR));
    PC lsqrPc = nullptr;
    PetscCallAbort(comm, KSPGetPC(lsqr, &lsqrPc));
    PetscCallAbort(comm, PCSetType(lsqrPc, PCNONE));
    PetscCall(KSPSetTolerances(
    lsqr,
    1.e-16,
    1.e-16,
    PETSC_DEFAULT,
    30000));
    PetscCallAbort(comm, KSPSetFromOptions(lsqr));
    PetscCallAbort(comm, VecSet(pressureE, 0.0));
    PetscCallAbort(comm, KSPSolve(lsqr, force_, pressureE));
    KSPConvergedReason reasonE;
    PetscCallAbort(comm, KSPGetConvergedReason(lsqr, &reasonE));
    PetscCallAbort(comm, MatMult(gradient_, pressureE, projected));
    PetscCallAbort(comm, VecWAXPY(residualE, -1.0, projected, force_));
    PetscReal normE = 0.0, orthE = 0.0;
    PetscCallAbort(comm, VecNorm(residualE, NORM_2, &normE));
    PetscCallAbort(comm, MatMultTranspose(gradient_, residualE, orth));
    PetscCallAbort(comm, VecNorm(orth, NORM_2, &orthE));

    // Mass-weighted projection:
    //   (G^T M^{-1} G) p = G^T M^{-1} f.
    // Use a shell matrix, never form M^{-1} or the normal equations explicitly.
    PetscCallAbort(comm, KSPSolve(massSolver, force_, tmp));
    PetscScalar dotF;
    PetscCallAbort(comm, VecDot(force_, tmp, &dotF));
    const double normFM = std::sqrt(
        std::max(0.0, static_cast<double>(PetscRealPart(dotF))));
    PetscCallAbort(comm, MatMultTranspose(gradient_, tmp, rhs));

    ShellContext context;
    context.gradient = gradient_;
    context.massSolver = massSolver;
    PetscCallAbort(comm, VecDuplicate(force_, &context.work));
    PetscCallAbort(comm, VecDuplicate(force_, &context.riesz));

    PetscInt localSize = 0, globalSize = 0;
    PetscCallAbort(comm, VecGetLocalSize(force_, &localSize));
    PetscCallAbort(comm, VecGetSize(force_, &globalSize));
    Mat schur = nullptr;
    PetscCallAbort(comm, MatCreateShell(
        comm, localSize, localSize, globalSize, globalSize, &context, &schur));
    PetscCallAbort(comm, MatShellSetOperation(
        schur, MATOP_MULT,
        reinterpret_cast<void (*)(void)>(&applyWeightedOperator)));

    
    // Approximate weighted pressure operator:
    //
    // P = G^T diag(M)^{-1} G
    //
    // Used ONLY as a preconditioner.
    // The true operator remains G^T M^{-1} G.

    Vec invMassDiagonal = nullptr;
    PetscCall(VecDuplicate(force_, &invMassDiagonal));

    PetscCall(MatGetDiagonal(mass_, invMassDiagonal));
    PetscCall(VecReciprocal(invMassDiagonal));

    Mat scaledGradient = nullptr;
    PetscCall(MatDuplicate(
        gradient_, MAT_COPY_VALUES, &scaledGradient));

    PetscCall(MatDiagonalScale(
        scaledGradient, invMassDiagonal, nullptr));

    Mat weightedPmat = nullptr;

    PetscCall(MatTransposeMatMult(
        gradient_,
        scaledGradient,
        MAT_INITIAL_MATRIX,
        PETSC_DEFAULT,
        &weightedPmat));

    // Inactive DOFs and the velocity columns of the full
    // NS layout have zero pressure-operator diagonal.
    // Give these rows a unit diagonal in the preconditioner
    // only. Do not modify the true projection operator.

    PetscCall(MatSetOption(
        weightedPmat,
        MAT_NEW_NONZERO_ALLOCATION_ERR,
        PETSC_FALSE));

    Vec precondDiagonal = nullptr;
    PetscCall(MatCreateVecs(
        weightedPmat, &precondDiagonal, nullptr));

    PetscCall(MatGetDiagonal(
        weightedPmat, precondDiagonal));

    PetscInt nDiagonalLocal = 0;
    PetscCall(VecGetLocalSize(
        precondDiagonal, &nDiagonalLocal));

    PetscScalar* values = nullptr;
    PetscCall(VecGetArray(precondDiagonal, &values));

    for (PetscInt i = 0; i < nDiagonalLocal; ++i) {
      values[i] =
          PetscAbsScalar(values[i]) < 1.e-30 ? 1.0 : 0.0;
    }

    PetscCall(VecRestoreArray(precondDiagonal, &values));

    PetscCall(MatDiagonalSet(
        weightedPmat, precondDiagonal, ADD_VALUES));

    PetscCall(MatAssemblyBegin(
        weightedPmat, MAT_FINAL_ASSEMBLY));

    PetscCall(MatAssemblyEnd(
        weightedPmat, MAT_FINAL_ASSEMBLY));

    PetscCall(VecDestroy(&precondDiagonal));
    PetscCall(VecDestroy(&invMassDiagonal));
    PetscCall(MatDestroy(&scaledGradient));


    // KSP weighted = nullptr;
    // PetscCallAbort(comm, KSPCreate(comm, &weighted));
    // PetscCallAbort(comm, KSPSetOptionsPrefix(weighted, "ex46_weighted_"));
    // PetscCallAbort(comm, KSPSetOperators(weighted, schur, schur));
    // PetscCallAbort(comm, KSPSetType(weighted, KSPCG));
    // PC weightedPc = nullptr;
    // PetscCallAbort(comm, KSPGetPC(weighted, &weightedPc));
    // PetscCallAbort(comm, PCSetType(weightedPc, PCNONE));
    // PetscCallAbort(comm, KSPSetTolerances(
    //     weighted, 1.e-10, 1.e-14, PETSC_DEFAULT, 10000));
    // PetscCallAbort(comm, KSPSetFromOptions(weighted));
    // PetscCallAbort(comm, VecSet(pressureM, 0.0));
    // PetscCallAbort(comm, KSPSolve(weighted, rhs, pressureM));

    KSP weighted = nullptr;

    PetscCall(KSPCreate(comm, &weighted));
    PetscCall(KSPSetOptionsPrefix(
        weighted, "ex46_weighted_"));

    // Exact matrix-free operator + approximate preconditioner.
    PetscCall(KSPSetOperators(
        weighted, schur, weightedPmat));

    PetscCall(KSPSetType(weighted, KSPCG));

    PC weightedPc = nullptr;
    PetscCall(KSPGetPC(weighted, &weightedPc));
    PetscCall(PCSetType(weightedPc, PCJACOBI));

    PetscCall(KSPSetTolerances(
        weighted,
        1.e-16,
        1.e-16,
        PETSC_DEFAULT,
        30000));

    // Start from the Euclidean LSQR pressure.
    PetscCall(VecCopy(pressureE, pressureM));
    PetscCall(KSPSetInitialGuessNonzero(
        weighted, PETSC_TRUE));

    PetscCall(KSPSetFromOptions(weighted));

    PetscCall(KSPSolve(weighted, rhs, pressureM));

    KSPConvergedReason reasonM;
    PetscCallAbort(comm, KSPGetConvergedReason(weighted, &reasonM));

    PetscCallAbort(comm, MatMult(gradient_, pressureM, projected));
    PetscCallAbort(comm, VecWAXPY(residualM, -1.0, projected, force_));
    PetscCallAbort(comm, KSPSolve(massSolver, residualM, tmp));
    PetscScalar dotR;
    PetscCallAbort(comm, VecDot(residualM, tmp, &dotR));
    const double normRM = std::sqrt(
        std::max(0.0, static_cast<double>(PetscRealPart(dotR))));
    PetscCallAbort(comm, MatMultTranspose(gradient_, tmp, orth));
    PetscReal orthM = 0.0;
    PetscCallAbort(comm, VecNorm(orth, NORM_2, &orthM));


    // Euclidean norm of the MASS-WEIGHTED solution residual.
    PetscReal euclideanAtWeighted = 0.0;

    PetscCall(VecNorm(
        residualM, NORM_2, &euclideanAtWeighted));

    // Mass-weighted norm of the EUCLIDEAN solution residual.
    PetscCall(KSPSolve(
        massSolver, residualE, tmp));

    PetscScalar dotEuclideanWeighted = 0.0;

    PetscCall(VecDot(
        residualE, tmp, &dotEuclideanWeighted));

    const PetscReal massAtEuclidean =
        PetscSqrtReal(PetscMax(
            0.0, PetscRealPart(dotEuclideanWeighted)));

    PetscCall(PetscPrintf(comm,
        "\n--- Cross-projection verification ---\n"
        "Euclidean residual (LSQR)       = %.12e\n"
        "Euclidean residual (weighted)   = %.12e\n"
        "Mass residual (LSQR)            = %.12e\n"
        "Mass residual (weighted)        = %.12e\n",
        (double)normE,
        (double)euclideanAtWeighted,
        (double)massAtEuclidean,
        (double)normRM));



    // Relative optimality of the weighted projection.
    PetscReal rhsNorm = 0.0;
    PetscCall(VecNorm(rhs, NORM_2, &rhsNorm));

    const double relativeOptimality =
        rhsNorm > 0.0 ? orthM / rhsNorm : 0.0;

    PetscInt weightedIterations = 0;
    PetscCall(KSPGetIterationNumber(
        weighted, &weightedIterations));

    PetscReal weightedKspResidual = 0.0;
    PetscCall(KSPGetResidualNorm(
        weighted, &weightedKspResidual));

    PetscCall(PetscPrintf(comm,
        "\n--- Weighted projection convergence ---\n"
        "Iterations                = %d\n"
        "KSP residual              = %.12e\n"
        "Relative optimality       = %.12e\n"
        "Convergence reason        = %d\n",
        (int)weightedIterations,
        (double)weightedKspResidual,
        relativeOptimality,
        (int)reasonM));


    PetscCallAbort(comm, PetscPrintf(comm,
        "\n========== EX46 SURFACE FORCE REPRESENTABILITY ==========\n"
        "||f_sigma||_2                          = %.12e\n"
        "Euclidean min ||Gp-f||_2                = %.12e\n"
        "Euclidean relative error               = %.12e\n"
        "Euclidean ||G^T (f-Gp)||_2              = %.12e\n"
        "||f_sigma||_(M^-1)                     = %.12e\n"
        "Mass-weighted min ||Gp-f||_(M^-1)       = %.12e\n"
        "Mass-weighted relative error           = %.12e\n"
        "Mass-weighted ||G^T M^-1 (f-Gp)||_2    = %.12e\n"
        "LSQR reason / weighted CG reason       = %d / %d\n"
        "=========================================================\n",
        static_cast<double>(normF), static_cast<double>(normE),
        normF > 0.0 ? static_cast<double>(normE / normF) : 0.0,
        static_cast<double>(orthE), normFM, normRM,
        normFM > 0.0 ? normRM / normFM : 0.0,
        static_cast<double>(orthM), static_cast<int>(reasonE),
        static_cast<int>(reasonM)));

    if (reasonE <= 0 || reasonM <= 0)
      PetscCallAbort(comm, PetscPrintf(comm,
          "[ex46 balance] WARNING: projection solver did not converge; "
          "reported errors may overestimate the true minima.\n"));

    PetscCallAbort(comm, KSPDestroy(&weighted));
    PetscCallAbort(comm, MatDestroy(&schur));
    PetscCallAbort(comm, VecDestroy(&context.work));
    PetscCallAbort(comm, VecDestroy(&context.riesz));
    PetscCallAbort(comm, KSPDestroy(&lsqr));
    PetscCallAbort(comm, KSPDestroy(&massSolver));
    PetscCallAbort(comm, VecDestroy(&pressureE));
    PetscCallAbort(comm, VecDestroy(&pressureM));
    PetscCallAbort(comm, VecDestroy(&projected));
    PetscCallAbort(comm, VecDestroy(&residualE));
    PetscCallAbort(comm, VecDestroy(&residualM));
    PetscCallAbort(comm, VecDestroy(&tmp));
    PetscCallAbort(comm, VecDestroy(&rhs));
    PetscCallAbort(comm, VecDestroy(&orth));
  }

  void destroy() {
    if (gradient_) PetscCallAbort(PETSC_COMM_WORLD, MatDestroy(&gradient_));
    if (mass_) PetscCallAbort(PETSC_COMM_WORLD, MatDestroy(&mass_));
    if (force_) PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&force_));
  }

private:
  struct ShellContext {
    Mat gradient = nullptr;
    KSP massSolver = nullptr;
    Vec work = nullptr;
    Vec riesz = nullptr;
  };

  static PetscErrorCode applyWeightedOperator(Mat mat, Vec x, Vec y) {
    ShellContext* ctx = nullptr;
    PetscCall(MatShellGetContext(mat, &ctx));
    PetscCall(MatMult(ctx->gradient, x, ctx->work));
    PetscCall(KSPSolve(ctx->massSolver, ctx->work, ctx->riesz));
    PetscCall(MatMultTranspose(ctx->gradient, ctx->riesz, y));
    return 0;
  }

  Mat gradient_ = nullptr;
  Mat mass_ = nullptr;
  Vec force_ = nullptr;
};

} // namespace ex46_balance
#endif