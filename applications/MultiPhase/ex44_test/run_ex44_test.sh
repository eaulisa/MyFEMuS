#!/usr/bin/env bash
set -euo pipefail

# ------------------ EDIT ONLY THIS SECTION ------------------
ELEMENTS=(0 1 3 4)     # 0=HEX, 1=TET, 3=QUAD, 4=TRI
PSI="ball"             # ball, saye, stella (Stella only on straight QUAD9)
CURVED=1               # 0=straight, 1=analytical curved BALL
LEVELS=4               # 0 through LEVELS
CURVATURE=0.40         # only used for curved BALL
STELLA_A=4             # >4 recommended to keep oval strictly inside [-0.5,0.5]^2
STELLA_B=4             # =4 means touching the corresponding cell sides
STELLA_C=160           # x^2*y^2 coefficient
MODE="moment"          # moment, algoim, compare
DEGREE=8
ALGOIM_ORDER=9
ALGOIM_SUBDIVISIONS=0
VTK=1
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
if ((${#element_names[@]} == 0)); then
  echo "ELEMENTS must contain at least one element ID" >&2
  exit 1
fi
elements_csv="$(IFS=,; echo "${element_names[*]}")"

cmd=("$EXECUTABLE"
  --elements "$elements_csv"
  --psi "$PSI"
  --curved "$CURVED"
  --levels "$LEVELS"
  --curvature "$CURVATURE"
  --mode "$MODE"
  --degree "$DEGREE"
  --algoim-order "$ALGOIM_ORDER"
  --algoim-subdivisions "$ALGOIM_SUBDIVISIONS"
  --vtk "$VTK")

if [[ ",${PSI}," == *",stella,"* ]]; then
  cmd+=(--stella-a "$STELLA_A"
        --stella-b "$STELLA_B"
        --stella-c "$STELLA_C")
fi

"${cmd[@]}"
