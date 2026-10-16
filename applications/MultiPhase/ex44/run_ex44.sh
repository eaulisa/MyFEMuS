#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# USER CONFIGURATION
# ============================================================
# ------------------------------------------------------------
# MPI
# ------------------------------------------------------------

NPROC=8

# ------------------------------------------------------------
# Problem
# ------------------------------------------------------------

DIM=2
SIMULATION="bb1"

# ------------------------------------------------------------
# Grid configurations
# ------------------------------------------------------------

UNIFORM_LEVELS=(2)
ADAPTIVE_LEVELS=(10)

LEVEL_OFFSET=2

# ------------------------------------------------------------
# Time discretization
# ------------------------------------------------------------

NSTEPS=(1200)

# ------------------------------------------------------------
# Velocity multigrid
# ------------------------------------------------------------

MULTIGRID_TYPE="mg"

# ------------------------------------------------------------
# Executable / options
# ------------------------------------------------------------

EXECUTABLE="./MultiPhase_ex44"
OPTIONS_FILE="./input/FS_solver.options"

# ------------------------------------------------------------
# PETSc configuration
# ------------------------------------------------------------

KSP_MONITOR=true
KSP_VIEW=false
OPTIONS_LEFT=true


# ============================================================
# Basic checks
# ============================================================

if [[ "${MULTIGRID_TYPE}" != "hmg" &&
      "${MULTIGRID_TYPE}" != "mg" ]]; then
  echo "ERROR: MULTIGRID_TYPE must be 'hmg' or 'mg'."
  exit 1
fi

if (( NPROC < 1 )); then
  echo "ERROR: NPROC must be >= 1."
  exit 1
fi

if (( LEVEL_OFFSET < 0 )); then
  echo "ERROR: LEVEL_OFFSET must be >= 0."
  exit 1
fi


# ============================================================
# Generate FS_solver.options
# ============================================================

generate_options_file() {

  local uniform_levels="$1"
  local adaptive_levels="$2"
  local level_offset="$3"

  local nlevels

  nlevels=$((uniform_levels + adaptive_levels - level_offset))

  if (( nlevels < 1 )); then

    echo "ERROR: Number of solver levels must be >= 1." >&2
    echo "       N = uniform_levels + adaptive_levels - level_offset" >&2
    echo "       N = ${uniform_levels} + ${adaptive_levels} - ${level_offset}" >&2

    return 1
  fi

  mkdir -p "$(dirname "${OPTIONS_FILE}")"

  : > "${OPTIONS_FILE}"


  # ==========================================================
  # Level 0
  # ==========================================================

  cat >> "${OPTIONS_FILE}" <<EOF_OPTIONS
# ============================================================
# Level 0
# ============================================================

# Velocity multigrid type:
#   hmg -> AMG
#   mg  -> GMG

-ksp_type fgmres
-ksp_rtol 1.e-8
-ksp_atol 1.e-12
-ksp_divtol 1.e+50
-ksp_max_it 50

-level-0ksp_max_it 10
-level-0fieldsplit_0_pc_type ${MULTIGRID_TYPE}

# -level-0fieldsplit_1_pc_type ${MULTIGRID_TYPE}

-level-0fieldsplit_1_ksp_type preonly
-level-0fieldsplit_1_pc_type lu
-level-0fieldsplit_1_pc_factor_mat_solver_type mumps
EOF_OPTIONS


  # ==========================================================
  # Outer multigrid levels >= 1
  # ==========================================================

  for ((i=1; i<nlevels; ++i)); do

    cat >> "${OPTIONS_FILE}" <<EOF_OPTIONS
# ============================================================
# Level ${i}
# ============================================================

# ------------------------------------------------------------
# FieldSplit level solver
# ------------------------------------------------------------


-level-${i}ksp_type fgmres
#-level-${i}ksp_pc_side right
-level-${i}ksp_norm_type unpreconditioned
#-level-${i}ksp_monitor
-level-${i}ksp_view_final_residual
-level-${i}ksp_rtol 1.e-8
-level-${i}ksp_atol 1.e-12
-level-${i}ksp_divtol 1.e+50
-level-${i}ksp_max_it 10

-level-${i}pc_type fieldsplit
-level-${i}pc_fieldsplit_type schur
-level-${i}pc_fieldsplit_schur_fact_type upper
-level-${i}pc_fieldsplit_schur_precondition selfp


# -level-${i}ksp_norm_type unpreconditioned
#-level-${i}ksp_monitor

# -level-${i}ksp_type gmres
# -level-${i}ksp_rtol 1.e-8
# -level-${i}ksp_atol 1.e-12
# -level-${i}ksp_divtol 1.e+50
# -level-${i}ksp_max_it 10
#
# -level-${i}pc_fieldsplit_schur_fact_type upper
# -level-${i}pc_fieldsplit_schur_precondition selfp

# ------------------------------------------------------------
# Velocity
# ------------------------------------------------------------

-level-${i}fieldsplit_0_ksp_type preonly
-level-${i}fieldsplit_0_pc_type ${MULTIGRID_TYPE}
-level-${i}fieldsplit_0_mg_levels_ksp_type chebyshev
-level-${i}fieldsplit_0_mg_levels_ksp_max_it 4
-level-${i}fieldsplit_0_mg_levels_pc_type jacobi

# AMG-specific options

# -level-${i}fieldsplit_0_pc_hmg_reuse_interpolation true
# -level-${i}fieldsplit_0_pc_hmg_use_subspace_coarsening false
# -level-${i}fieldsplit_0_pc_hmg_use_matmaij false
# -level-${i}fieldsplit_0_pc_hmg_coarsening_component 0
# -level-${i}fieldsplit_0_hmg_inner_pc_type gamg
# -level-${i}fieldsplit_0_hmg_inner_pc_gamg_aggressive_square_graph false

# ------------------------------------------------------------
# Pressure
# ------------------------------------------------------------

# -level-${i}fieldsplit_1_ksp_type preonly
# -level-${i}fieldsplit_1_pc_type ${MULTIGRID_TYPE}
# -level-${i}fieldsplit_1_mg_levels_ksp_type chebyshev
# -level-${i}fieldsplit_1_mg_levels_ksp_max_it 4
# -level-${i}fieldsplit_1_mg_levels_ksp_norm_type none
# -level-${i}fieldsplit_1_mg_levels_pc_type jacobi

# -level-${i}fieldsplit_1_ksp_norm_type unpreconditioned
#-level-${i}fieldsplit_1_ksp_monitor


# -level-${i}fieldsplit_1_ksp_type preonly
# -level-${i}fieldsplit_1_pc_type lu
# -level-${i}fieldsplit_1_pc_factor_mat_solver_type mumps

# -level-${i}fieldsplit_1_ksp_type richardson
# -level-${i}fieldsplit_1_ksp_richardson_scale 1
# -level-${i}fieldsplit_1_ksp_max_it 1




-level-${i}fieldsplit_1_ksp_type gmres
#-level-${i}fieldsplit_1_ksp_pc_side right
-level-${i}fieldsplit_1_ksp_norm_type unpreconditioned
#-level-${i}fieldsplit_1_ksp_monitor
-level-${i}fieldsplit_1_ksp_initial_guess_nonzero false
-level-${i}fieldsplit_1_ksp_knoll false
-level-${i}fieldsplit_1_ksp_rtol 1.e-2
-level-${i}fieldsplit_1_ksp_atol 1.e-12
-level-${i}fieldsplit_1_ksp_divtol 1.e+50
-level-${i}fieldsplit_1_ksp_max_it 3
-level-${i}fieldsplit_1_ksp_gmres_restart 3

-level-${i}fieldsplit_1_pc_type asm
-level-${i}fieldsplit_1_pc_use_amat false
-level-${i}fieldsplit_1_pc_asm_type restrict
-level-${i}fieldsplit_1_pc_asm_overlap 1
-level-${i}fieldsplit_1_pc_asm_local_type multiplicative

-level-${i}fieldsplit_1_sub_ksp_type preonly
-level-${i}fieldsplit_1_sub_pc_type ilu
-level-${i}fieldsplit_1_sub_pc_factor_mat_solver_type petsc
-level-${i}fieldsplit_1_sub_pc_factor_levels 1
-level-${i}fieldsplit_1_sub_pc_factor_shift_type nonzero


# -level-${i}fieldsplit_1_ksp_type preonly
# -level-${i}fieldsplit_1_pc_type asm
# -level-${i}fieldsplit_1_pc_use_amat false
# -level-${i}fieldsplit_1_pc_asm_type restrict
# -level-${i}fieldsplit_1_pc_asm_overlap 1
# -level-${i}fieldsplit_1_pc_asm_local_type multiplicative
#
# -level-${i}fieldsplit_1_sub_ksp_type preonly
# -level-${i}fieldsplit_1_sub_pc_type ilu
# -level-${i}fieldsplit_1_sub_pc_factor_mat_solver_type petsc
# -level-${i}fieldsplit_1_sub_pc_factor_levels 1
# -level-${i}fieldsplit_1_sub_pc_factor_shift_type nonzero

# AMG-specific options

# -level-${i}fieldsplit_1_pc_hmg_reuse_interpolation true
# -level-${i}fieldsplit_1_pc_hmg_use_subspace_coarsening false
# -level-${i}fieldsplit_1_pc_hmg_use_matmaij false
# -level-${i}fieldsplit_1_pc_hmg_coarsening_component 0
# -level-${i}fieldsplit_1_hmg_inner_pc_type gamg
# -level-${i}fieldsplit_1_hmg_inner_pc_gamg_aggressive_square_graph false











EOF_OPTIONS

  done

  printf '%s\n' "${nlevels}"
}


# ============================================================
# PETSc command-line options
# ============================================================

PETSC_OPTIONS=()

if [[ "${KSP_MONITOR}" == true ]]; then
  PETSC_OPTIONS+=(
    -ksp_monitor_true_residual
  )
fi

if [[ "${KSP_VIEW}" == true ]]; then
  PETSC_OPTIONS+=(
    -ksp_view
  )
fi

if [[ "${OPTIONS_LEFT}" == true ]]; then
  PETSC_OPTIONS+=(
    -options_left
  )
fi

if (( NPROC > 1 )); then
  PETSC_OPTIONS+=(
    -matptap_via hypre
  )
fi


# ============================================================
# Number of runs
# ============================================================

TOTAL_RUNS=$(( \
  ${#UNIFORM_LEVELS[@]} * \
  ${#ADAPTIVE_LEVELS[@]} * \
  ${#NSTEPS[@]} \
))

RUN_INDEX=0


# ============================================================
# Run all requested configurations
# ============================================================

for uniform_levels in "${UNIFORM_LEVELS[@]}"; do

  for adaptive_levels in "${ADAPTIVE_LEVELS[@]}"; do

    # --------------------------------------------------------
    # Generate PETSc options for this grid
    # --------------------------------------------------------

    NLEVELS=$(
      generate_options_file \
        "${uniform_levels}" \
        "${adaptive_levels}" \
        "${LEVEL_OFFSET}"
    )


    # --------------------------------------------------------
    # Run all requested time resolutions on this grid
    # --------------------------------------------------------

    for nsteps in "${NSTEPS[@]}"; do

      RUN_INDEX=$((RUN_INDEX + 1))

      echo
      echo "============================================================"
      echo "ex43 run ${RUN_INDEX}/${TOTAL_RUNS}"
      echo "============================================================"
      echo "MPI ranks        : ${NPROC}"
      echo "Uniform levels   : ${uniform_levels}"
      echo "Adaptive levels  : ${adaptive_levels}"
      echo "Level offset     : ${LEVEL_OFFSET}"
      echo "Solver levels    : ${NLEVELS}"
      echo "Number of steps  : ${nsteps}"
      echo "Dimension        : ${DIM}"
      echo "Simulation       : ${SIMULATION}"
      echo "Velocity MG      : ${MULTIGRID_TYPE}"
      echo "Options file     : ${OPTIONS_FILE}"
      echo "============================================================"
      echo

      mpirun -n "${NPROC}" \
        "${EXECUTABLE}" \
        "${PETSC_OPTIONS[@]}" \
        -ksp_knoll true \
        -fp_trap\
        -ksp_error_if_not_converged true\
        -ksp_converged_reason\
        -options_file "${OPTIONS_FILE}" \
        --uniform-levels "${uniform_levels}" \
        --adaptive-levels "${adaptive_levels}" \
        --level-offset "${LEVEL_OFFSET}" \
        --nsteps "${nsteps}" \
        --dim "${DIM}" \
        --simulation "${SIMULATION}"

      echo

    done
  done
done
