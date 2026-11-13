#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# USER CONFIGURATION
# ============================================================
# ------------------------------------------------------------
# MPI
# ------------------------------------------------------------

NPROC=4
NPRINT=1000
NREINIT=1

# ------------------------------------------------------------
# Problem
# ------------------------------------------------------------

DIM=2
SIMULATION="par"
SYMMETRIC=1
AXISYMMETRIC=0

# ------------------------------------------------------------
# Grid configurations
# ------------------------------------------------------------

UNIFORM_LEVELS=(1)
ADAPTIVE_LEVELS=(4)
LEVEL_OFFSET_PSI=0
LEVEL_OFFSET=0

# ------------------------------------------------------------
# Time discretization
# ------------------------------------------------------------

NSTEPS=(100000)

# ------------------------------------------------------------
# Cutfem integration
# ------------------------------------------------------------

CUTFEM_INTEGRATION="quadratic"     # linear | quadratic
SURFACE_MOMENT_ORDER=10

# ------------------------------------------------------------
# Velocity multigrid
# ------------------------------------------------------------

MULTIGRID_TYPE="mg"

# ------------------------------------------------------------
# Executable / options
# ------------------------------------------------------------

EXECUTABLE="./MultiPhase_ex46"
OPTIONS_FILE="./input/FS_solver.options"

# ------------------------------------------------------------
# PETSc configuration
# ------------------------------------------------------------

KSP_MONITOR=false
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

if (( SYMMETRIC == 1 && AXISYMMETRIC == 1 )); then
  echo "ERROR: SYMMETRIC and AXISYMMETRIC cannot both be 1."
  exit 1
fi

if [[ "${CUTFEM_INTEGRATION}" != "linear" &&
      "${CUTFEM_INTEGRATION}" != "quadratic" ]]; then
  echo "ERROR: CUTFEM_INTEGRATION must be linear or quadratic."
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
-ksp_atol 1.e-10
-ksp_divtol 1.e+50
-ksp_max_it 200

-level-0ksp_max_it 200

-level-0ksp_norm_type unpreconditioned
# -level-0ksp_monitor

-level-0fieldsplit_0_pc_type ${MULTIGRID_TYPE}
-level-0fieldsplit_0_sub_ksp_type preonly
-level-0fieldsplit_0_sub_pc_type lu
-level-0fieldsplit_0_sub_factor_mat_solver_type mumps


-level-0fieldsplit_1_pc_type ${MULTIGRID_TYPE}
-level-0fieldsplit_1_sub_ksp_type preonly
-level-0fieldsplit_1_sub_pc_type lu
-level-0fieldsplit_1_sub_factor_mat_solver_type mumps

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

-level-${i}ksp_norm_type unpreconditioned
# -level-${i}ksp_view_final_residual
# -level-${i}ksp_monitor

-level-${i}ksp_type gmres
-level-${i}ksp_rtol 1.e-8
-level-${i}ksp_atol 1.e-10
-level-${i}ksp_divtol 1.e+50
-level-${i}ksp_max_it 200

-level-${i}pc_fieldsplit_schur_fact_type upper
-level-${i}pc_fieldsplit_schur_precondition selfp

# ------------------------------------------------------------
# Velocity
# ------------------------------------------------------------

-level-${i}fieldsplit_0_ksp_type preonly
-level-${i}fieldsplit_0_pc_type ${MULTIGRID_TYPE}
-level-${i}fieldsplit_0_mg_levels_ksp_type chebyshev
-level-${i}fieldsplit_0_mg_levels_ksp_max_it 4
-level-${i}fieldsplit_0_mg_levels_pc_type jacobi

# ------------------------------------------------------------
# Pressure
# ------------------------------------------------------------

-level-${i}fieldsplit_1_ksp_type preonly
-level-${i}fieldsplit_1_pc_type ${MULTIGRID_TYPE}
-level-${i}fieldsplit_1_mg_levels_ksp_type chebyshev
-level-${i}fieldsplit_1_mg_levels_ksp_max_it 4
-level-${i}fieldsplit_1_mg_levels_pc_type jacobi

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
      echo "MPI ranks          : ${NPROC}"
      echo "Uniform levels     : ${uniform_levels}"
      echo "Adaptive levels    : ${adaptive_levels}"
      echo "Level offset       : ${LEVEL_OFFSET}"
      echo "Level offset psi   : ${LEVEL_OFFSET_PSI}"
      echo "CutFEM integration : ${CUTFEM_INTEGRATION}"
      echo "Solver levels      : ${NLEVELS}"
      echo "Number of steps    : ${nsteps}"
      echo "Dimension          : ${DIM}"
      echo "Simulation         : ${SIMULATION}"
      echo "Velocity MG        : ${MULTIGRID_TYPE}"
      echo "Options file       : ${OPTIONS_FILE}"
      echo "============================================================"
      echo

      mpirun -n "${NPROC}" \
        "${EXECUTABLE}" \
        "${PETSC_OPTIONS[@]}" \
        -options_file "${OPTIONS_FILE}" \
        --uniform-levels "${uniform_levels}" \
        --adaptive-levels "${adaptive_levels}" \
        --level-offset "${LEVEL_OFFSET}" \
        --level-offset-psi "${LEVEL_OFFSET_PSI}" \
        --cutfem-integration "${CUTFEM_INTEGRATION}" \
        --surface-moment-order "${SURFACE_MOMENT_ORDER}" \
        --nsteps "${nsteps}" \
        --dim "${DIM}" \
        --simulation "${SIMULATION}" \
        --symmetric "${SYMMETRIC}" \
        --axisymmetric "${AXISYMMETRIC}" \
        --nprint "${NPRINT}" \
        --nreinit "${NREINIT}"  

      echo

    done
  done
done
