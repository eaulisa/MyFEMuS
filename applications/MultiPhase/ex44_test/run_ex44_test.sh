#!/usr/bin/env bash
set -euo pipefail

# ------------------ EDIT ONLY THIS SECTION ------------------
ELEMENTS=(0 3)           # 0=HEX, 3=QUAD; TRI/TET remain available for BALL
PSI="nonconvex"          # ball or nonconvex (QUAD9/HEX27 only)
CURVED=0                 # 0=straight; 1=analytical curved BALL only
LEVELS=4                 # levels 0..LEVELS, preserving the uploaded default
CURVATURE=0.40           # used only for curved BALL

# nonconvex: QUAD9 uses the original Stella Q2 oval.
NONCONVEX_QUAD_A=4        # a>=4; a>4 for a strictly interior oval
NONCONVEX_QUAD_B=4        # b>=4; b>4 for a strictly interior oval
NONCONVEX_QUAD_C=160

# nonconvex: HEX27 uses the Lauricella F_D interface reference.
NONCONVEX_HEX_A=4        # a,b>=1 in [-1,1]^3 coordinates
NONCONVEX_HEX_B=4
NONCONVEX_HEX_C=160000      # chi=c/(a*b)>3 for nonconvex geometry

MODE="moment"            # moment, algoim, compare
DEGREE=8
CUTFEM_ABS_TOL=1e-8      # CutFEM moment-fitting absolute tolerance
CUTFEM_REL_TOL=1e-8      # CutFEM moment-fitting relative tolerance
ALGOIM_ORDER=8
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
  --cutfem-abs-tol "$CUTFEM_ABS_TOL"
  --cutfem-rel-tol "$CUTFEM_REL_TOL"
  --algoim-order "$ALGOIM_ORDER"
  --algoim-subdivisions "$ALGOIM_SUBDIVISIONS"
  --vtk "$VTK")

if [[ ",$PSI," == *",nonconvex,"* ]]; then
  cmd+=(--nonconvex-quad-a "$NONCONVEX_QUAD_A"
        --nonconvex-quad-b "$NONCONVEX_QUAD_B"
        --nonconvex-quad-c "$NONCONVEX_QUAD_C"
        --nonconvex-hex-a "$NONCONVEX_HEX_A"
        --nonconvex-hex-b "$NONCONVEX_HEX_B"
        --nonconvex-hex-c "$NONCONVEX_HEX_C")
fi

"${cmd[@]}"
