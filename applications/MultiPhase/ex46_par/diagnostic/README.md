# Surface–pressure analysis (work in progress)

The Python analysis is implemented, but ex46 does not yet export its input or
invoke it. It therefore cannot yet measure the production assembly.

The proposed input directory contains `metadata.json` with an integer
`mpi_ranks` and, for each rank `N`, the following CSV files with headers:

* `dofs_N.csv`: `system_dof,kind,free,prescribed,fixed`. The kind is `U`,
  `V`, `W`, `P1`, or `P2`; free is 0 or 1. Pressure values are the
  gauge-adjusted Young–Laplace value and the actual fixed value. Records must
  agree if a system DOF occurs on more than one rank.
* `matrix_N.csv`: `row,col,value`, the pressure contribution to the ex46
  velocity residual (the negative of its Jacobian pressure block).
* `forces_N.csv`: `row,surface,pressure_action`. The latter is independently
  accumulated at quadrature points using the prescribed pressure.

Entries may repeat: the analysis sums rank and element contributions, filters
velocity rows using the exported constraints, and minimizes over free pressure
columns while retaining fixed pressure contributions. The matrix includes
constrained columns for the prescribed-action check. The exporter must obtain
constraints and fixed values from the actual FEMuS solution on its owning rank.

Run `python3 solve_surface_pressure_range.py OUTPUT_DIRECTORY`. NumPy and SciPy
are required. The program writes `result.json` and returns a failure status for
inconsistent pressure action, failed optimality, failed negative control, or
LSMR nonconvergence. It reconstructs residuals in long-double arithmetic and
checks scaled normal equations independently of LSMR's stopping code.
