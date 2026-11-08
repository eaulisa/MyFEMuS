# Native reference CutFemWeight_LS verification

This standalone suite uses the actual FEMuS finite elements and Gauss tables. It does not require PETSc or build the rest of FEMuS. Requirements are a C++14 compiler, Boost headers, Eigen 3 headers, Python 3, and NumPy.

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
