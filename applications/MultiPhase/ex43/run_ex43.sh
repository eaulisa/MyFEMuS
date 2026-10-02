# #!/usr/bin/env bash
#
# set -e
#
# # ============================================================
# # Default run parameters
# # ============================================================
# NPROC=8
# UNIFORM_LEVELS=1
# ADAPTIVE_LEVELS=2
# LEVEL_OFFSET=1
# NSTEPS=1200
# DIM=3
# SIMULATION="rb1"
#
# EXECUTABLE="./MultiPhase_ex43"
# OPTIONS_FILE="./input/FS_solver.options"
#
# # ============================================================
# # Optional command-line overrides
# #
# # Usage:
# #   ./run_ex43.sh [uniform-levels] [adaptive-levels] [level-offset]
# #
# # Example:
# #   ./run_ex43.sh 2 3 1
# # ============================================================
# if [ $# -ge 1 ]; then
#   UNIFORM_LEVELS="$1"
# fi
#
# if [ $# -ge 2 ]; then
#   ADAPTIVE_LEVELS="$2"
# fi
#
# if [ $# -ge 3 ]; then
#   LEVEL_OFFSET="$3"
# fi
#
# # ============================================================
# # Number of solver levels
# # ============================================================
# N=$((UNIFORM_LEVELS + ADAPTIVE_LEVELS - LEVEL_OFFSET))
#
# # if [ "$N" -lt 1 ]; then
# #   echo "ERROR: Number of solver levels must be >= 1."
# #   echo "       N = uniform-levels + adaptive-levels - level-offset"
# #   echo "       N = ${UNIFORM_LEVELS} + ${ADAPTIVE_LEVELS} - ${LEVEL_OFFSET} = ${N}"
# #   exit 1
# # fi
#
# echo "============================================================"
# echo "ex43 run configuration"
# echo "============================================================"
# echo "MPI ranks       : ${NPROC}"
# echo "Uniform levels  : ${UNIFORM_LEVELS}"
# echo "Adaptive levels : ${ADAPTIVE_LEVELS}"
# echo "Level offset    : ${LEVEL_OFFSET}"
# echo "Solver levels N : ${N}"
# echo "Options file    : ${OPTIONS_FILE}"
# echo "============================================================"
#
# # ============================================================
# # Generate FS_solver.options
# # ============================================================
# mkdir -p "$(dirname "${OPTIONS_FILE}")"
#
# : > "${OPTIONS_FILE}"
#
# for ((i=1; i<N; ++i)); do
#
#   cat >> "${OPTIONS_FILE}" <<EOF
# # ============================================================
# # Level ${i}
# # ============================================================
#
# # Fieldsplit level solver/preconditioner
# -level-${i}ksp_richardson_scale .4
#
# -level-${i}ksp_rtol 1.e-8
# -level-${i}ksp_atol 1.e-12
# -level-${i}ksp_divtol 1.e+50
# -level-${i}ksp_max_it 2
# -level-${i}ksp_norm_type none
#
# -level-${i}pc_fieldsplit_schur_fact_type upper
# -level-${i}pc_fieldsplit_schur_precondition selfp
#
# # Velocity
#
# # -level-${i}fieldsplit_0_ksp_type preonly
#
# # -level-${i}fieldsplit_0_pc_hmg_reuse_interpolation true
# # -level-${i}fieldsplit_0_pc_hmg_use_subspace_coarsening false
# # -level-${i}fieldsplit_0_pc_hmg_use_matmaij false
# # -level-${i}fieldsplit_0_pc_hmg_coarsening_component 0
# # -level-${i}fieldsplit_0_hmg_inner_pc_type gamg
# # -level-${i}fieldsplit_0_hmg_inner_pc_gamg_aggressive_square_graph false
# #
# # -level-${i}fieldsplit_0_mg_levels_ksp_type chebyshev
# # -level-${i}fieldsplit_0_mg_levels_ksp_max_it 4
# # -level-${i}fieldsplit_0_mg_levels_ksp_norm_type none
# # -level-${i}fieldsplit_0_mg_levels_pc_type jacobi
#
# -level-${i}fieldsplit_0_ksp_type preonly
#
# -level-${i}fieldsplit_0_pc_type hypre
# -level-${i}fieldsplit_0_pc_hypre_type boomeramg
#
# -level-${i}fieldsplit_0_pc_hypre_boomeramg_max_iter 1
# -level-${i}fieldsplit_0_pc_hypre_boomeramg_tol 0.0
#
# -level-${i}fieldsplit_0_pc_hypre_boomeramg_grid_sweeps_down 1
# -level-${i}fieldsplit_0_pc_hypre_boomeramg_grid_sweeps_up 1
# -level-${i}fieldsplit_0_pc_hypre_boomeramg_grid_sweeps_coarse 1
#
# -level-${i}fieldsplit_0_pc_hypre_boomeramg_relax_type_all SOR/Jacobi
#
# -level-${i}fieldsplit_0_pc_hypre_boomeramg_coarsen_type HMIS
# -level-${i}fieldsplit_0_pc_hypre_boomeramg_interp_type ext+i
#
# # Pressure
# # -level-${i}fieldsplit_1_pc_hmg_reuse_interpolation true
# # -level-${i}fieldsplit_1_pc_hmg_use_subspace_coarsening false
# # -level-${i}fieldsplit_1_pc_hmg_use_matmaij false
# # -level-${i}fieldsplit_1_pc_hmg_coarsening_component 0
# # -level-${i}fieldsplit_1_hmg_inner_pc_type gamg
# #  -level-${i}fieldsplit_1_hmg_inner_pc_gamg_aggressive_square_graph false
#
# # -level-${i}fieldsplit_1_mg_levels_ksp_type chebyshev
# # -level-${i}fieldsplit_1_mg_levels_ksp_max_it 2
# # -level-${i}fieldsplit_1_mg_levels_ksp_norm_type none
# # -level-${i}fieldsplit_1_mg_levels_pc_type sor
# # -level-${i}fieldsplit_1_mg_levels_pc_sor_local_symmetric
#
# -level-${i}fieldsplit_1_ksp_type preonly
#
# -level-${i}fieldsplit_1_pc_type hypre
# -level-${i}fieldsplit_1_pc_hypre_type boomeramg
#
# -level-${i}fieldsplit_1_pc_hypre_boomeramg_max_iter 1
# -level-${i}fieldsplit_1_pc_hypre_boomeramg_tol 0.0
#
# -level-${i}fieldsplit_1_pc_hypre_boomeramg_grid_sweeps_down 1
# -level-${i}fieldsplit_1_pc_hypre_boomeramg_grid_sweeps_up 1
# -level-${i}fieldsplit_1_pc_hypre_boomeramg_grid_sweeps_coarse 1
#
# -level-${i}fieldsplit_1_pc_hypre_boomeramg_relax_type_all SOR/Jacobi
#
# EOF
#
# done
#
# echo "Generated ${OPTIONS_FILE} for levels 1 through ${N}."
# echo
#
# # ============================================================
# # Run ex43
# # ============================================================
# mpirun -n "${NPROC}" "${EXECUTABLE}" \
#   -matptap_via allatonce \
#   -ksp_monitor_true_residual \
#   -log_view_memory \
#   -ksp_view \
#   -options_left \
#   -options_file "${OPTIONS_FILE}" \
#   --uniform-levels "${UNIFORM_LEVELS}" \
#   --adaptive-levels "${ADAPTIVE_LEVELS}" \
#   --level-offset "${LEVEL_OFFSET}" \
#   --nsteps "${NSTEPS}" \
#   --dim "${DIM}" \
#   --simulation "${SIMULATION}"

#!/usr/bin/env bash

set -e

# ============================================================
# Default run parameters
# ============================================================
NPROC=8
UNIFORM_LEVELS=1
ADAPTIVE_LEVELS=2
LEVEL_OFFSET=1
NSTEPS=1200
DIM=3
SIMULATION="rb1"

# ============================================================
# AMG backend
#
# Available:
#   gamg   -> PETSc HMG with GAMG-generated hierarchy
#   hypre  -> HYPRE BoomerAMG
# ============================================================
AMG_BACKEND="gamg"

EXECUTABLE="./MultiPhase_ex43"
OPTIONS_FILE="./input/FS_solver.options"

# ============================================================
# Optional command-line overrides
#
# Usage:
#   ./run_ex43.sh [uniform-levels] [adaptive-levels] [level-offset]
#
# Example:
#   ./run_ex43.sh 2 3 1
# ============================================================
if [ $# -ge 1 ]; then
  UNIFORM_LEVELS="$1"
fi

if [ $# -ge 2 ]; then
  ADAPTIVE_LEVELS="$2"
fi

if [ $# -ge 3 ]; then
  LEVEL_OFFSET="$3"
fi

# ============================================================
# Check AMG backend
# ============================================================
case "${AMG_BACKEND}" in
  gamg|hypre)
    ;;
  *)
    echo "ERROR: AMG_BACKEND must be either 'gamg' or 'hypre'."
    echo "       Current value: ${AMG_BACKEND}"
    exit 1
    ;;
esac

# ============================================================
# Number of solver levels
# ============================================================
N=$((UNIFORM_LEVELS + ADAPTIVE_LEVELS - LEVEL_OFFSET))

echo "============================================================"
echo "ex43 run configuration"
echo "============================================================"
echo "MPI ranks       : ${NPROC}"
echo "Uniform levels  : ${UNIFORM_LEVELS}"
echo "Adaptive levels : ${ADAPTIVE_LEVELS}"
echo "Level offset    : ${LEVEL_OFFSET}"
echo "Solver levels N : ${N}"
echo "AMG backend     : ${AMG_BACKEND}"
echo "Options file    : ${OPTIONS_FILE}"
echo "============================================================"

# ============================================================
# Generate FS_solver.options
# ============================================================
mkdir -p "$(dirname "${OPTIONS_FILE}")"

: > "${OPTIONS_FILE}"

for ((i=1; i<N; ++i)); do

  cat >> "${OPTIONS_FILE}" <<EOF
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

EOF

  # ==========================================================
  # GAMG through PETSc HMG
  # ==========================================================
  if [ "${AMG_BACKEND}" = "gamg" ]; then

    cat >> "${OPTIONS_FILE}" <<EOF
# ------------------------------------------------------------
# Velocity - PETSc HMG / GAMG
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
# Pressure - PETSc HMG / GAMG
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

EOF

  # ==========================================================
  # HYPRE BoomerAMG
  # ==========================================================
  elif [ "${AMG_BACKEND}" = "hypre" ]; then

    cat >> "${OPTIONS_FILE}" <<EOF
# ------------------------------------------------------------
# Velocity - HYPRE BoomerAMG
# ------------------------------------------------------------
-level-${i}fieldsplit_0_ksp_type preonly

-level-${i}fieldsplit_0_pc_type hypre
-level-${i}fieldsplit_0_pc_hypre_type boomeramg

-level-${i}fieldsplit_0_pc_hypre_boomeramg_max_iter 1
-level-${i}fieldsplit_0_pc_hypre_boomeramg_tol 0.0

-level-${i}fieldsplit_0_pc_hypre_boomeramg_grid_sweeps_down 1
-level-${i}fieldsplit_0_pc_hypre_boomeramg_grid_sweeps_up 1
-level-${i}fieldsplit_0_pc_hypre_boomeramg_grid_sweeps_coarse 1

-level-${i}fieldsplit_0_pc_hypre_boomeramg_relax_type_all SOR/Jacobi
-level-${i}fieldsplit_0_pc_hypre_boomeramg_coarsen_type HMIS
-level-${i}fieldsplit_0_pc_hypre_boomeramg_interp_type ext+i

# ------------------------------------------------------------
# Pressure - HYPRE BoomerAMG
# ------------------------------------------------------------
-level-${i}fieldsplit_1_ksp_type preonly

-level-${i}fieldsplit_1_pc_type hypre
-level-${i}fieldsplit_1_pc_hypre_type boomeramg

-level-${i}fieldsplit_1_pc_hypre_boomeramg_max_iter 1
-level-${i}fieldsplit_1_pc_hypre_boomeramg_tol 0.0

-level-${i}fieldsplit_1_pc_hypre_boomeramg_grid_sweeps_down 1
-level-${i}fieldsplit_1_pc_hypre_boomeramg_grid_sweeps_up 1
-level-${i}fieldsplit_1_pc_hypre_boomeramg_grid_sweeps_coarse 1

-level-${i}fieldsplit_1_pc_hypre_boomeramg_relax_type_all SOR/Jacobi

EOF

  fi

done

echo "Generated ${OPTIONS_FILE} for levels 1 through $((N - 1))."
echo "AMG backend: ${AMG_BACKEND}"
echo

# ============================================================
# Run ex43
# ============================================================
mpirun -n "${NPROC}" "${EXECUTABLE}" \
  -matptap_via allatonce \
  -ksp_monitor_true_residual \
  -log_view_memory \
  -ksp_view \
  -options_left \
  -options_file "${OPTIONS_FILE}" \
  --uniform-levels "${UNIFORM_LEVELS}" \
  --adaptive-levels "${ADAPTIVE_LEVELS}" \
  --level-offset "${LEVEL_OFFSET}" \
  --nsteps "${NSTEPS}" \
  --dim "${DIM}" \
  --simulation "${SIMULATION}"

