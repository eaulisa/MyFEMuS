#!/usr/bin/env bash
set -euo pipefail

# ------------------ EDIT ONLY THIS SECTION ------------------
ELEMENTS=(0 1 3 4)    # 0=HEX, 1=TET, 3=QUAD, 4=TRI
PSI="ball"            # ball or saye
CURVED=1              # 0=straight, 1=curved
LEVELS=3              # AMR refinements (test levels 0..LEVELS)
EPSILON=0.20          # manufactured curved BALL deformation
MODE="moment"         # moment, algoim, compare
DEGREE=7              # moment-fitting degree
ALGOIM_ORDER=5
ALGOIM_SUBDIVISIONS=0
VTK=1                # 0=off, 1=on
EXECUTABLE="./MultiPhase_ex44_test"
# ------------------------------------------------------------

cd "$(dirname "${BASH_SOURCE[0]}")"

element_names=()
for element in "${ELEMENTS[@]}"; do
  case "$element" in
    0) element_names+=(hex) ;;
    1) element_names+=(tet) ;;
    3) element_names+=(quad) ;;
    4) element_names+=(tri) ;;
    *) echo "Unsupported element ID: $element" >&2; exit 1 ;;
  esac
done

elements_csv="$(IFS=,; echo "${element_names[*]}")"

"$EXECUTABLE" \
  --elements "$elements_csv" \
  --psi "$PSI" \
  --curved "$CURVED" \
  --levels "$LEVELS" \
  --epsilon "$EPSILON" \
  --mode "$MODE" \
  --degree "$DEGREE" \
  --algoim-order "$ALGOIM_ORDER" \
  --algoim-subdivisions "$ALGOIM_SUBDIVISIONS" \
  --vtk "$VTK"