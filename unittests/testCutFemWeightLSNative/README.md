# Native reference CutFemWeight_LS verification

This standalone suite uses the actual FEMuS finite elements and Gauss tables. It does not require PETSc or build the rest of FEMuS. Requirements are a C++14 compiler, Boost headers, Eigen 3 headers, Python 3, NumPy, and mpmath.

## Build and run

From this directory:

```sh
python build.py --boost /path/to/boost/include --eigen /path/to/eigen3
python run.py
python run.py --selection supplemental
python run_focused.py
python compare_planar.py
python run_focused.py --extended
```

The Boost option is unnecessary when Boost is on the compiler's default include path. Eigen defaults to `/usr/include/eigen3`. The build also accepts `--repo`, `--cxx`, `--jobs`, and `--target regression|contracts|geometry`. Environment variables `BOOST_INCLUDEDIR`, `EIGEN3_INCLUDE_DIR`, and `CXX` supply the corresponding defaults.

For a targeted rerun, `build.py --tests native_suite tet_standard_moments` rebuilds named translation units, `run.py --selection original --element TET` selects one element, and `run_focused.py --tests tet_standard_moments simplex_high_degree` selects focused probes.

Every executable is built from these test sources and FEMuS's `src/fe` and `src/quadrature` sources. Build commands, logs, dependency lists, and implementation hashes are written to `build/`. Per-case numerical results are written to `results/candidate/`; aggregate results are JSON files beside the scripts. The runners return a nonzero exit status for failed checks, integration exceptions, missing executables, or timeouts.

Direct parent-coordinate moments are checked by `direct_parent_moments` and compared with the pinned previous implementation by `python run_direct_moment_benchmark.py`; the latter writes `review_results/direct_moment_summary.json`.

The certified single-root comparison is reproduced with `python run_direct_moment_benchmark.py --baseline 188a0034139361a1772fc0d8a5c0151c2cba8469 --folder single_root --summary single_root_summary.json --rounds 3`; `single_graph_root` verifies both root branches and certificate forwarding.

## Coverage and independent references

| Coverage | Cases |
| --- | --- |
| Elements and fitting degrees | LINE3, QUAD9, TRI7, HEX27, TET15; degrees 0 through 5 |
| Level sets | Affine axis and oblique cuts; quadratic monotone fields; coupled tensor and enriched simplex polynomials; closed circles and spheres; two line roots |
| Edge cases | Full and empty cells; interface on a cell face; near-affine fields; tiny caps; positive scaling by powers of two and sign reversal |
| Geometry | Reference cells; affine anisotropic scaling; curved shear; variable-Jacobian stretching |
| Integrals | Positive volume, raw delta, reference surface, physical surface, and weighted native monomial moments |
| Formulations | Norm included inside the moments; raw delta correction multiplied by the pointwise FEMuS gradient or cofactor norm |
| Additional checks | FE field/gradient and geometry/Jacobian evaluation; repeated cache output; API rejection and independence; multiprecision; both include orders; original planar class and all tetrahedral bases |

The independent Python references use exact interval, simplex, ball, and sphere formulas where available. Curved graphs use an independently derived root parameterization and Gauss-Legendre integration at two orders. Those references do not call the implementation's integration kernels. They require the two reference orders to agree before assessing a result. The C++ driver independently evaluates G, grad G, and the geometry Jacobian with the unchanged FEMuS finite-element basis.

The fitting space is the full total-degree polynomial space P_qM. Consequently, multiplying a raw delta correction by a nonpolynomial gradient norm after fitting does not generally give the same finite-degree result as including that norm inside the moments. The suite asserts equality only when the complete polynomial integrand is within the degree budget. Otherwise it records the outside formulation's error against the independent surface reference. Variable-Jacobian volume tests use the same degree-budget rule. Physical surface tests use reference Gauss weights and the cofactor metric without a second determinant factor.

## Strict small-cut regression

The default original matrix contains 386 cases, including the TET15 degree-5 cap `x > 1 - 1/256` at requested absolute and relative tolerances of 1e-13. This is a required regression. Its 45-second process budget covers raw delta, reference surface, physical surface, and a cached repeat. A failure or timeout is reported with a nonzero exit status.

The strict cap's normalized weighted-moment limit is 2e-13, including all representable inside/outside formulations. It can be selected alone or as part of the complete matrix:

```sh
python run.py --stress
python run.py --selection original
```

`--selection supplemental` runs eight additional cases: six degree-5 cap variants at 1e-13 covering level-set scaling by 2^-20 and 2^20, sign reversal, and affine physical geometry; the separately named `_relaxed_1e-10` cap request; and the scaled degree-5 plane at 1e-13 alongside its original-matrix counterpart at 1e-12. The relaxed request never substitutes for a failed strict request. `regular` is retained as an alias for the original strict matrix.

`test_tet_slice_regression.cpp` reconstructs the diagnosed first slice from actual FEMuS nodes and calls the existing line/triangle kernels. It rejects invalid corrected intervals, preserves valid thin intervals, checks every degree-5 moment, and requires nonzero accurate interface mass from every completed chart. An ill-conditioned chart may explicitly reject; regular charts and chart selection must succeed at the original final tolerance. The same test can be compiled against the preserved pre-fix headers. `test_root_candidate_regression.cpp` separately checks consolidation of duplicate candidates for the same simple root while preserving genuinely distinct nearby roots, independent polynomial identities, and the positive strip next to a boundary.

## Focused checks

`run_focused.py` runs API contracts, output range rejection, standalone multiprecision, native coefficient agreement with actual FEMuS values and gradients, and direct QUAD/TRI/TET/HEX moment probes. QUAD and HEX probes include degree 10. `--extended` also runs the degree-8/10 TRI/TET checks, with a five-minute process budget. The build checks both planar/LS include orders and rejects accidental dependencies on the planar or Algoim headers from the standalone LS tests.

`compare_planar.py` checks 147 affine rules against the unchanged original `CutFemWeight` class, including tetrahedral base selections 0, 1, and 2. Its tiny-cap comparisons require the original strict result.

The exact original-commit snapshots, original/candidate comparison files, and optional direct Algoim application check are external verification artifacts rather than prerequisites for this candidate-only suite.

## QUAD/TRI refinement review

`python3 run_review_validation.py` checks moments through total degree 10, 144 circle convergence cases, 48 translated ellipse cases, selective Taylor/full-GK fallback, and independent 60-digit graph references. AddressSanitizer and UndefinedBehaviorSanitizer run with recovery disabled. The runner also checks 30 mesh-refinement configurations: QUAD and TRI, 4 through 64 subdivisions per coordinate, and degrees 2, 6, and 10. The mesh probe reconstructs the same circle from nodal FE values and checks global area and interface length against analytic values.

Timing uses 108 coefficient, nodal, and reconstructed-input cases plus 240 deterministic stress cases. Identical benchmark sources and compiler flags run against pinned `rifat`, the earlier reviewed commit, and the candidate. Runs are pinned to one allowed CPU and interleaved in opposite orders. The comparison rejects failed integrations, missing moments, nonfinite values, and momentwise disagreements. Results, source hashes, environment information, and acceptance are written to `review_results/` and summarized in the final-validation workflow logs.

Existing pinned worktrees can be supplied with `--baseline-repo` and `--reviewed-repo`; their commit identities and CutFEM header contents are verified before use. `BOOST_INCLUDEDIR`, `EIGEN3_INCLUDE_DIR`, and `CXX` configure dependencies consistently with `build.py`. The historical TRI source-comparison workflow treats old revisions as reports; the candidate probe and final-validation workflow enforce failures.

The [2026-11-09 recheck report](review_results/RECHECK_2026-11-09.md) records fixes, accuracy checks, mesh-refinement results, measured speedups, and remaining timing regressions.

## HEX/TET graph refinement

`python3 run_hex_tet_validation.py` tests the default HEX/TET local-graph paths at degrees 2, 6, and 10. HEX uses rectangular subregions and the QUAD integrator. TET keeps convex regions clipped from the standard tetrahedron and triangulates each slice for the TRI integrator. Both retain the parent FE polynomial, certify a nonzero graph derivative, and select directions per region. Critical-point splits isolate uncut cores before integrating closed surfaces. Exact TET zero faces retain the endpoint trace rule and its half-weight convention.

The checks cover analytic plane and sphere moments, tetrahedral clipping and parent-polynomial identities, original-face restriction, strict caps, scaling and sign reversal, reference and physical surface weights, 18 sphere mesh-refinement configurations up to 512 HEX / 3,072 TET elements, tolerances 1e-8 through 1e-12, and small enclosed spheres. ASan/UBSan run with recovery disabled. Timings use identical sources and flags against pinned pre-3D and rifat revisions; the summary reports common successful cases and any baseline failures. No requested monomial degree exceeds 10.

The existing CutFemWeight_LS API selects the new paths automatically. `--skip-benchmarks` and `--skip-sanitizers` support focused reruns; the final workflow runs both. Raw local outputs go to `results/hex_tet/`; compact acceptance and source hashes go to `review_results/hex_tet_summary.json`. The separate report records the numerical results and timing limitations.

`python3 validate_enriched_3d.py` adds closed TET15 interior-bubble and coupled non-graph HEX tube fields, with both signs and both raw-delta and reference-area modes at degrees 2, 6, and 10. Every requested moment is compared with independent transformed Gauss quadrature at orders 96 and 160. These tests also require phase certification by Bernstein coefficient refinement when a direct control-net bound is inconclusive.

## Unit-coordinate baseline validation

After building, run `python3 run_wide_unit_validation.py`, `python3 run_wide_unit_mesh.py`, `python3 run_wide_double_output.py`, and `python3 run_wide_unit_sanitizers.py`. The separate unit-baseline workflow also reruns the established native, supplemental, focused, and planar checks. Pass `--leak-detection` to enable LeakSanitizer outside a ptrace environment.

The wide matrix checks QUAD/TRI/HEX/TET at degrees 2, 5, and 10 against independent analytic or angular references. It covers full/empty cells, boundary interfaces, axis and oblique graphs, two sheets, enclosed and tangent balls, ellipses, and simplex bubbles; tolerance sweeps, sign reversal, scaling, direction validation, and prescribed directions add further cases. Every requested unit monomial is checked for raw volume, raw delta, and reference surface moments. Fitted corrections and nonlinear physical-surface cofactors are checked separately. A failed or missing case, unresolved graph region, nonfinite result, or changed production source rejects the run.

The mesh driver reconstructs a translated circle or sphere from actual nodal values and uses the complete fitted physical-surface API on cut cells. It checks physical volume, surface measure, and their x-squared moments through 64 subdivisions per coordinate in 2D and 8 in 3D, including 3,072 tetrahedra. CSV output records actual cell splits and time at every level. Compact results and production hashes are saved in `review_results/wide_unit*.json`; full numerical outputs are in `results/wide_unit*/`.
