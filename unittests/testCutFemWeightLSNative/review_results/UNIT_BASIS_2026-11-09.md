# Direct unit-cell CutFEM reconstruction and moment basis

LINE, QUAD, and HEX now reconstruct the LS polynomial directly from FEMuS-ordered nodal values at unit-cell coordinates 0, 1/2, and 1. The tensor Lagrange polynomials are 1-3u+2u^2, 4u-4u^2, and -u+2u^2. The gradient is differentiated from that unit-coordinate polynomial. No completed moment table is reconstructed into parent-coordinate monomials. TRI and TET retain their existing standard simplex coordinates and nodal reconstruction.

`CutFemWeight_LS` integrates unit monomials and uses that same basis in the weight-fitting system. The tensor Gram entries are products of 1/(a+b+1). A shifted Legendre factorization evaluates the same inverse Gram action in multiprecision during setup, preserving monomial moments and avoiding a poorly conditioned explicit dense inverse. The moment-to-weight matrix is retained in TypeA and both output sums use compensated accumulation. Long-double unit moments reach fitting without an intermediate double conversion.

The unit HEX path uses the existing HEX -> QUAD -> LINE hierarchy, graph partitioning, per-moment refinement masks, and cached graph branch. Unit-cube geometry is selected explicitly; symmetry assumptions about the origin of the old cube are disabled for this path. QUAD uses the same shared LINE kernel on its unit reconstruction. Direction diagnostics and boundary-face half weights are retained.

FEMuS Gaussian coordinates and the application API remain unchanged. Their existing volume weights supply the factor 2^dimension. Reference-surface moments are scaled by 1/2 before fitting, and physical-surface moments reconstructed with the unit-coordinate geometry cofactor are scaled by 2^-dimension. These factors preserve reference and physical surface measures. Delta moments use the directly reconstructed unit-coordinate derivative.

The benchmark compares against e2a9c1056812fe568a54b412e39acd3671e07bdd. Ten fixed physical graph families per element are tested at degrees 2, 5, and 10, with three alternating paired runs, identical -O3 -DNDEBUG -march=x86-64 -mtune=generic builds, and CPU 0 affinity. Each case uses calibrated batches and the median of five measurements. The table averages the ten cases and three runs. Times measure moment integration, excluding one-time construction of the fitting matrix and full ex44 assembly or solves.

The requested absolute and relative tolerances are 1e-10 in the selected basis: parent-coordinate moments for the baseline, unit-cell moments and unit-cell measures for the new tensor path. This intentionally changes the error target with the requested basis; it is not a guarantee of the same absolute error for every old parent monomial. End-to-end fitted weights are checked separately against analytic moments and the established independent regression references. TRI/TET are unchanged timing controls; small timing differences there are noise/code-layout effects.

| Degree | Element | Parent basis ms/cell | Unit basis ms/cell | Speedup | Outer slices old/new |
|---:|:---|---:|---:|---:|---:|
| 2 | QUAD | 0.1363 | 0.1061 | 1.28x | 0.0/0.0 |
| 2 | TRI | 0.1083 | 0.1015 | 1.07x | 0.0/0.0 |
| 2 | HEX | 4.5626 | 3.5508 | 1.28x | 13.6/13.6 |
| 2 | TET | 7.4123 | 7.0467 | 1.05x | 40.5/40.5 |
| 5 | QUAD | 0.1876 | 0.1183 | 1.59x | 0.0/0.0 |
| 5 | TRI | 0.1179 | 0.1176 | 1.00x | 0.0/0.0 |
| 5 | HEX | 9.4090 | 4.4323 | 2.12x | 19.6/13.6 |
| 5 | TET | 9.3850 | 9.2592 | 1.01x | 46.5/46.5 |
| 10 | QUAD | 0.4402 | 0.2342 | 1.88x | 0.0/0.0 |
| 10 | TRI | 0.1543 | 0.1494 | 1.03x | 0.0/0.0 |
| 10 | HEX | 51.3170 | 13.5015 | 3.80x | 69.1/24.1 |
| 10 | TET | 15.7142 | 15.3506 | 1.02x | 58.5/58.5 |

All 386 native cases pass (36,860 accuracy assertions), as do 8 supplemental cases (3,176 assertions), 19 focused executions, and 147 comparisons against planar rules (16,258 weight and 4,552 moment comparisons; maximum normalized moment difference 1.57e-13). New tests pass 4,212 checks for direct nodal reconstruction, gradients against FEMuS shape functions, unit moments through degree 10, LS scaling, coordinate-face traces with independent direction validation, and fitted Gaussian weights. The maximum plane unit-moment error is 5.42e-20; the maximum fitted-moment error is 1.20e-13. ASan/UBSan pass the same 4,212 new checks; leak detection is disabled because the runtime uses ptrace. Actual FEMuS FE and Gaussian-quadrature sources are compiled in these checks. Full MPI/PETSc ex44 execution was not run in this workspace.

Reproduce the new focused check with `python3 unittests/testCutFemWeightLSNative/build.py --tests unit_element_moments` and `python3 unittests/testCutFemWeightLSNative/run_focused.py --tests unit_element_moments`. Reproduce timings with `python3 unittests/testCutFemWeightLSNative/run_unit_basis_benchmark.py`. Dependency include directories can be supplied through the existing build options/environment variables. The benchmark saves individual case CSVs and source hashes; the committed summary and validation JSONs contain aggregate results and the tested header hashes.
