#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# PARASITIC CURRENTS: EDIT ONLY THIS CONFIGURATION SECTION
# ============================================================
NPROC=8
DIM=3                               # 2 | 3
MESHES=(imported)                  # generated | imported; may run both

# Imported mesh: Q1-projected QUAD9/HEX27 .neu files.
IMPORTED_MESH_2D="./input/unstructured_PC_coarse_Q1.neu"
IMPORTED_MESH_3D="./input/unstructured_cube_hex27_Q1.neu"
MESH_SCALE=100                      # .neu coordinates / MESH_SCALE

# Generated mesh: number of *coarse* elements in each direction.
NX=8
NY=8
NZ=8

UNIFORM_LEVELS=(1)
NSTEPS=(100)                       # dt = 0.03 / NSTEPS; not actual solves
SOLVE_STEPS=1                       # actual NS steps, with fixed interface
INITIAL_PRESSURE="young-laplace"   # young-laplace | zero

# Keep uploaded ex46 defaults initially; vary one switch at a time.
PRESSURE_GP=1
VELOCITY_GP=0
STABILIZATION=0

SURFACE_MOMENT_ORDER=17
SURFACE_MOMENT_TOL=1.e-12
SURFACE_WEIGHTS=1                   # 3D: three signed-normal moments
SURFACE_GK="7/15"
SURFACE_EXTENDED_BUDGET=1
SURFACE_PRESSURE_DIAGNOSTIC=1

MULTIGRID_TYPE="mg"               # mg | hmg
EXECUTABLE="./MultiPhase_ex46_par"
OPTIONS_FILE="./input/FS_solver.options"

# ============================================================
# No need to edit below
# ============================================================
SOURCE_SCRIPT_DIR="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"
[[ "$DIM" == 2 || "$DIM" == 3 ]] || { echo "DIM must be 2 or 3" >&2; exit 1; }
[[ "$PRESSURE_GP" == 0 || "$PRESSURE_GP" == 1 ]] || exit 1
[[ "$VELOCITY_GP" == 0 || "$VELOCITY_GP" == 1 ]] || exit 1
[[ "$STABILIZATION" == 0 || "$STABILIZATION" == 1 ]] || exit 1
[[ "$SURFACE_PRESSURE_DIAGNOSTIC" == 0 || "$SURFACE_PRESSURE_DIAGNOSTIC" == 1 ]] || exit 1
[[ "$INITIAL_PRESSURE" == "zero" || "$INITIAL_PRESSURE" == "young-laplace" ]] || exit 1
[[ "$MULTIGRID_TYPE" == "mg" || "$MULTIGRID_TYPE" == "hmg" ]] || exit 1
(( NPROC > 0 && SOLVE_STEPS > 0 )) || exit 1

if [[ "$DIM" == 2 ]]; then IMPORTED_MESH="$IMPORTED_MESH_2D";
else IMPORTED_MESH="$IMPORTED_MESH_3D"; fi

# Preserve the original field-split solver profile for comparisons.
write_solver_options() {
  local levels="$1"
  mkdir -p "$(dirname "$OPTIONS_FILE")"
  cat > "$OPTIONS_FILE" <<EOF
-ksp_type fgmres
-ksp_rtol 1.e-20
-ksp_atol 1.e-20
-ksp_divtol 1.e+50
-ksp_max_it 200
-level-0ksp_max_it 200
-level-0ksp_norm_type unpreconditioned
-level-0fieldsplit_0_pc_type $MULTIGRID_TYPE
-level-0fieldsplit_0_sub_ksp_type preonly
-level-0fieldsplit_0_sub_pc_type lu
-level-0fieldsplit_0_sub_factor_mat_solver_type mumps
-level-0fieldsplit_1_pc_type $MULTIGRID_TYPE
-level-0fieldsplit_1_sub_ksp_type preonly
-level-0fieldsplit_1_sub_pc_type lu
-level-0fieldsplit_1_sub_factor_mat_solver_type mumps
EOF
  local i
  for ((i=1; i<levels; ++i)); do
    cat >> "$OPTIONS_FILE" <<EOF
-level-${i}ksp_norm_type unpreconditioned
-level-${i}ksp_type gmres
-level-${i}ksp_rtol 1.e-20
-level-${i}ksp_atol 1.e-20
-level-${i}ksp_max_it 200
-level-${i}pc_fieldsplit_schur_fact_type upper
-level-${i}pc_fieldsplit_schur_precondition selfp
-level-${i}fieldsplit_0_ksp_type preonly
-level-${i}fieldsplit_0_pc_type $MULTIGRID_TYPE
-level-${i}fieldsplit_0_mg_levels_ksp_type chebyshev
-level-${i}fieldsplit_0_mg_levels_ksp_max_it 4
-level-${i}fieldsplit_0_mg_levels_pc_type jacobi
-level-${i}fieldsplit_1_ksp_type preonly
-level-${i}fieldsplit_1_pc_type $MULTIGRID_TYPE
-level-${i}fieldsplit_1_mg_levels_ksp_type chebyshev
-level-${i}fieldsplit_1_mg_levels_ksp_max_it 4
-level-${i}fieldsplit_1_mg_levels_pc_type jacobi
EOF
  done
}

if [[ "$SURFACE_EXTENDED_BUDGET" == 1 ]]; then
  export FEMUS_CUTFEM_UNLIMITED=1
else
  unset FEMUS_CUTFEM_UNLIMITED || true
fi

if [[ "$SURFACE_PRESSURE_DIAGNOSTIC" == 1 ]]; then
  [[ -f "$SOURCE_SCRIPT_DIR/diagnostic/solve_surface_pressure_range.py" ]] || {
    echo "Missing Python pressure-range postprocessor in diagnostic/" >&2; exit 1;
  }
fi

mkdir -p levelsetdiagnostic
run=0
for mesh in "${MESHES[@]}"; do
  [[ "$mesh" == generated || "$mesh" == imported ]] || {
    echo "Unknown mesh mode: $mesh" >&2; exit 1;
  }
  if [[ "$mesh" == imported && ! -f "$IMPORTED_MESH" ]]; then
    echo "Missing imported mesh: $IMPORTED_MESH" >&2; exit 1
  fi
  for level in "${UNIFORM_LEVELS[@]}"; do
    (( level >= 1 )) || { echo "UNIFORM_LEVELS must be >= 1" >&2; exit 1; }
    write_solver_options "$level"
    for nsteps in "${NSTEPS[@]}"; do
      run=$((run + 1))
      prefix="levelsetdiagnostic/par_${mesh}_D${DIM}_U${level}_NS${nsteps}_run${run}"
      printf '\n=== PAR run %d: mesh=%s dim=%d levels=%d dt=0.03/%d solves=%d ===\n' \
        "$run" "$mesh" "$DIM" "$level" "$nsteps" "$SOLVE_STEPS"
      args=(
        --mesh "$mesh" --mesh-file "$IMPORTED_MESH" --mesh-scale "$MESH_SCALE"
        --dim "$DIM" --uniform-levels "$level" --nsteps "$nsteps"
        --solve-steps "$SOLVE_STEPS" --nx "$NX" --ny "$NY" --nz "$NZ"
        --initial-pressure "$INITIAL_PRESSURE"
        --pressure-gp "$PRESSURE_GP" --velocity-gp "$VELOCITY_GP"
        --stabilization "$STABILIZATION"
        --surface-moment-order "$SURFACE_MOMENT_ORDER"
        --surface-moment-tolerance "$SURFACE_MOMENT_TOL"
        --surface-weights "$SURFACE_WEIGHTS"
        --gk-regular-2d "$SURFACE_GK" --gk-tangent-2d "$SURFACE_GK"
        --gk-regular-3d "$SURFACE_GK" --gk-tangent-3d "$SURFACE_GK"
        --surface-pressure-diagnostic "$SURFACE_PRESSURE_DIAGNOSTIC"
        --surface-pressure-prefix "${prefix}_surface_pressure"
        --output-prefix "$prefix"
      )
      petsc=(-options_file "$OPTIONS_FILE")
      if (( NPROC > 1 )); then petsc+=(-matptap_via hypre); fi
      mpirun -n "$NPROC" "$EXECUTABLE" "${petsc[@]}" "${args[@]}"
      if [[ "$SURFACE_PRESSURE_DIAGNOSTIC" == 1 ]]; then
        OPENBLAS_NUM_THREADS=1 python3 \
          "$SOURCE_SCRIPT_DIR/diagnostic/solve_surface_pressure_range.py" \
          "${prefix}_surface_pressure"
      fi
      printf 'Parasitic-current metrics: %s_currents.dat\n' "$prefix"
    done
  done
done
