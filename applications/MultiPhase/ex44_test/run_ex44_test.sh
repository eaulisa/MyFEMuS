
#!/usr/bin/env bash
set -euo pipefail

# ------------------ EDIT ONLY THIS SECTION ------------------
ELEMENTS=(0)            # 0=HEX27, 3=QUAD9

MODE="moment"
MOMENT_BACKEND="native"

DEGREE=1
ALGOIM_ORDER=20

PSI="torus"
CURVED=0
LEVELS=5
CURVATURE=0.50

# Nonconvex parameters
NONCONVEX_QUAD_A=8
NONCONVEX_QUAD_B=8

NONCONVEX_HEX_A=4
NONCONVEX_HEX_B=4

# Sweep over C (same C for QUAD and HEX)
C_VALUES=(16)

CUTFEM_ABS_TOL=1e-6
CUTFEM_REL_TOL=1e-6
ALGOIM_SUBDIVISIONS=0
VTK=1

GK_REGULAR_2D="${GK_REGULAR_2D:-3/7}"
GK_TANGENT_2D="${GK_TANGENT_2D:-7/15}"
GK_REGULAR_3D="${GK_REGULAR_3D:-3/7}"
GK_TANGENT_3D="${GK_TANGENT_3D:-7/15}"

EXECUTABLE="./MultiPhase_ex44_test"
OUTPUT_DIR="results_C"
# ------------------------------------------------------------

cd "$(dirname "${BASH_SOURCE[0]}")"
mkdir -p "$OUTPUT_DIR"

element_names=()
for element in "${ELEMENTS[@]}"; do
  case "$element" in
    0) element_names+=(hex) ;;
    1) element_names+=(tet) ;;
    2) element_names+=(wedge) ;;
    3) element_names+=(quad) ;;
    4) element_names+=(tri) ;;
    *) echo "Unsupported element ID: $element" >&2; exit 1 ;;
  esac
done

if ((${#element_names[@]} == 0)); then
  echo "ELEMENTS must not be empty" >&2
  exit 1
fi

elements_csv="$(IFS=,; echo "${element_names[*]}")"

for C in "${C_VALUES[@]}"; do

  echo "=========================================="
  echo " Running C = $C"
  echo "=========================================="

  cmd=("$EXECUTABLE"
    --gk-regular-2d "$GK_REGULAR_2D"
    --gk-tangent-2d "$GK_TANGENT_2D"
    --gk-regular-3d "$GK_REGULAR_3D"
    --gk-tangent-3d "$GK_TANGENT_3D"
    --elements "$elements_csv"
    --psi "$PSI"
    --curved "$CURVED"
    --levels "$LEVELS"
    --curvature "$CURVATURE"
    --mode "$MODE"
    --moment-backend "$MOMENT_BACKEND"
    --degree "$DEGREE"
    --cutfem-abs-tol "$CUTFEM_ABS_TOL"
    --cutfem-rel-tol "$CUTFEM_REL_TOL"
    --algoim-order "$ALGOIM_ORDER"
    --algoim-subdivisions "$ALGOIM_SUBDIVISIONS"
    --vtk "$VTK")

  if [[ "$PSI" == "nonconvex" ]]; then
    cmd+=(
      --nonconvex-quad-a "$NONCONVEX_QUAD_A"
      --nonconvex-quad-b "$NONCONVEX_QUAD_B"
      --nonconvex-quad-c "$C"
      --nonconvex-hex-a "$NONCONVEX_HEX_A"
      --nonconvex-hex-b "$NONCONVEX_HEX_B"
      --nonconvex-hex-c "$C"
    )
  fi

  "${cmd[@]}" 2>&1 | tee "$OUTPUT_DIR/C_${C}.prn"

done

echo "All C values completed."
