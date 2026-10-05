#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Default run parameters
# ============================================================
NPROC=2
DIM=2
SIMULATION="rb2"

# Default run parameters
UNIFORM_LEVELS=1
LEVEL_OFFSET=1

# Parameters that may contain multiple values
ADAPTIVE_LEVELS=(2)
NSTEPS=(1200)

EXECUTABLE="./MultiPhase_ex43"
OPTIONS_FILE="./input/FS_solver.options"

# ============================================================
# Usage
# ============================================================
usage() {
  cat <<EOF_USAGE
Usage:
  $0 [--uniform-levels n] \
     [--adaptive-levels n1 n2 ...] \
     [--level-offset n] \
     [--nsteps n1 n2 ...]

Examples:
  $0 --uniform-levels 2 --adaptive-levels 2 3 4
  $0 --level-offset 0 --nsteps 600 1200 2400
  $0 --uniform-levels 2 --adaptive-levels 2 3 4 --level-offset 1 --nsteps 600 1200

Only --adaptive-levels and --nsteps accept multiple values.
The script runs the Cartesian product of those two lists.

--uniform-levels and --level-offset each accept exactly one value.
If an option is omitted, its default value defined at the top of the
script is used.
EOF_USAGE
}

# ============================================================
# Read a list of non-negative integers following an option.
#
# Arguments:
#   $1 = option name (for error messages)
#   $2 = minimum allowed value
#   remaining arguments = command-line tokens
#
# Parsing is performed directly below to keep array assignment simple.
# ============================================================
while (( $# > 0 )); do
  case "$1" in
    --uniform-levels)
      if (( $# < 2 )) || [[ "$2" == --* ]]; then
        echo "ERROR: --uniform-levels requires exactly one value."
        exit 1
      fi
      if [[ ! "$2" =~ ^[1-9][0-9]*$ ]]; then
        echo "ERROR: uniform levels must be a positive integer."
        echo "       Invalid value: $2"
        exit 1
      fi
      UNIFORM_LEVELS="$2"
      shift 2
      ;;

    --adaptive-levels)
      shift
      ADAPTIVE_LEVELS=()

      while (( $# > 0 )) && [[ "$1" != --* ]]; do
        if [[ ! "$1" =~ ^[0-9]+$ ]]; then
          echo "ERROR: adaptive levels must be non-negative integers."
          echo "       Invalid value: $1"
          exit 1
        fi
        ADAPTIVE_LEVELS+=("$1")
        shift
      done

      if (( ${#ADAPTIVE_LEVELS[@]} == 0 )); then
        echo "ERROR: --adaptive-levels requires at least one value."
        exit 1
      fi
      ;;

    --level-offset)
      if (( $# < 2 )) || [[ "$2" == --* ]]; then
        echo "ERROR: --level-offset requires exactly one value."
        exit 1
      fi
      if [[ ! "$2" =~ ^[0-9]+$ ]]; then
        echo "ERROR: level offset must be a non-negative integer."
        echo "       Invalid value: $2"
        exit 1
      fi
      LEVEL_OFFSET="$2"
      shift 2
      ;;

    --nsteps)
      shift
      NSTEPS=()

      while (( $# > 0 )) && [[ "$1" != --* ]]; do
        if [[ ! "$1" =~ ^[1-9][0-9]*$ ]]; then
          echo "ERROR: nsteps must be positive integers."
          echo "       Invalid value: $1"
          exit 1
        fi
        NSTEPS+=("$1")
        shift
      done

      if (( ${#NSTEPS[@]} == 0 )); then
        echo "ERROR: --nsteps requires at least one value."
        exit 1
      fi
      ;;

    -h|--help)
      usage
      exit 0
      ;;

    *)
      echo "ERROR: unknown argument '$1'."
      echo
      usage
      exit 1
      ;;
  esac
done

# ============================================================
# Generate FS_solver.options for one grid configuration
# ============================================================
generate_options_file() {
  local uniform_levels="$1"
  local adaptive_levels="$2"
  local level_offset="$3"
  local nlevels

  nlevels=$((uniform_levels + adaptive_levels - level_offset))

  if (( nlevels < 1 )); then
    echo "ERROR: Number of solver levels must be >= 1." >&2
    echo "       N = uniform-levels + adaptive-levels - level-offset" >&2
    echo "       N = ${uniform_levels} + ${adaptive_levels} - ${level_offset} = ${nlevels}" >&2
    return 1
  fi

  mkdir -p "$(dirname "${OPTIONS_FILE}")"
  : > "${OPTIONS_FILE}"

  for ((i=1; i<nlevels; ++i)); do
    cat >> "${OPTIONS_FILE}" <<EOF_OPTIONS
# ============================================================
# Level ${i}
# ============================================================

# Fieldsplit level solver/preconditioner
-level-${i}ksp_richardson_scale .4

-level-${i}ksp_rtol 1.e-8
-level-${i}ksp_atol 1.e-12
-level-${i}ksp_divtol 1.e+50
-level-${i}ksp_max_it 2
-level-${i}ksp_norm_type none

-level-${i}pc_fieldsplit_schur_fact_type upper
-level-${i}pc_fieldsplit_schur_precondition selfp

# ------------------------------------------------------------
# Velocity - PETSc HMG with GAMG hierarchy
# ------------------------------------------------------------
-level-${i}fieldsplit_0_ksp_type preonly

-level-${i}fieldsplit_0_pc_type hmg
-level-${i}fieldsplit_0_pc_hmg_reuse_interpolation true
-level-${i}fieldsplit_0_pc_hmg_use_subspace_coarsening false
-level-${i}fieldsplit_0_pc_hmg_use_matmaij false
-level-${i}fieldsplit_0_pc_hmg_coarsening_component 0

-level-${i}fieldsplit_0_hmg_inner_pc_type gamg
-level-${i}fieldsplit_0_hmg_inner_pc_gamg_aggressive_square_graph false

-level-${i}fieldsplit_0_mg_levels_ksp_type chebyshev
-level-${i}fieldsplit_0_mg_levels_ksp_max_it 4
-level-${i}fieldsplit_0_mg_levels_ksp_norm_type none
-level-${i}fieldsplit_0_mg_levels_pc_type jacobi

# ------------------------------------------------------------
# Pressure - PETSc HMG with GAMG hierarchy
# ------------------------------------------------------------
-level-${i}fieldsplit_1_ksp_type preonly

-level-${i}fieldsplit_1_pc_type hmg
-level-${i}fieldsplit_1_pc_hmg_reuse_interpolation true
-level-${i}fieldsplit_1_pc_hmg_use_subspace_coarsening false
-level-${i}fieldsplit_1_pc_hmg_use_matmaij false
-level-${i}fieldsplit_1_pc_hmg_coarsening_component 0

-level-${i}fieldsplit_1_hmg_inner_pc_type gamg
-level-${i}fieldsplit_1_hmg_inner_pc_gamg_aggressive_square_graph false

-level-${i}fieldsplit_1_mg_levels_ksp_type chebyshev
-level-${i}fieldsplit_1_mg_levels_ksp_max_it 2
-level-${i}fieldsplit_1_mg_levels_ksp_norm_type none
-level-${i}fieldsplit_1_mg_levels_pc_type sor
-level-${i}fieldsplit_1_mg_levels_pc_sor_local_symmetric

EOF_OPTIONS
  done

  printf '%s\n' "${nlevels}"
}

# ============================================================
# Run all requested combinations
# ============================================================
TOTAL_RUNS=$(( ${#ADAPTIVE_LEVELS[@]} * ${#NSTEPS[@]} ))
RUN_INDEX=0

for adaptive_levels in "${ADAPTIVE_LEVELS[@]}"; do

  N=$(generate_options_file \
    "${UNIFORM_LEVELS}" \
    "${adaptive_levels}" \
    "${LEVEL_OFFSET}")

  for nsteps in "${NSTEPS[@]}"; do
    RUN_INDEX=$((RUN_INDEX + 1))

    echo "============================================================"
    echo "ex43 run ${RUN_INDEX}/${TOTAL_RUNS}"
    echo "============================================================"
    echo "MPI ranks       : ${NPROC}"
    echo "Uniform levels  : ${UNIFORM_LEVELS}"
    echo "Adaptive levels : ${adaptive_levels}"
    echo "Level offset    : ${LEVEL_OFFSET}"
    echo "Solver levels N : ${N}"
    echo "Number of steps : ${nsteps}"
    echo "Dimension       : ${DIM}"
    echo "Simulation      : ${SIMULATION}"
    echo "Options file    : ${OPTIONS_FILE}"
    echo "============================================================"

    mpirun -n "${NPROC}" "${EXECUTABLE}" \
      -options_left \
      -options_file "${OPTIONS_FILE}" \
      --uniform-levels "${UNIFORM_LEVELS}" \
      --adaptive-levels "${adaptive_levels}" \
      --level-offset "${LEVEL_OFFSET}" \
      --nsteps "${nsteps}" \
      --dim "${DIM}" \
      --simulation "${SIMULATION}"

    echo
  done
done


      #-matptap_via allatonce \
