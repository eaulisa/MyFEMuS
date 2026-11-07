# Native reference CutFemWeight_LS verification

This standalone suite uses the actual FEMuS finite elements and Gauss tables. It does not require PETSc or build the rest of FEMuS. Requirements are a C++14 compiler, Boost headers, Eigen 3 headers, Python 3, and NumPy.

## Build and run

From this directory:

```sh
python build.py --boost /path/to/boost/include --eigen /path/to/eigen3
python run.py
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

## Strict small-cut stress case

The original comparison matrix contains 386 cases. Its TET15 degree-5 cap `x > 1 - 1/256` requests absolute and relative tolerances of 1e-13. That case exceeded the 45-second process budget in both the original and native-reference implementations during the initial comparison. A run includes raw delta, reference surface, and physical surface in that budget.

The default portable suite keeps the other 385 cases and runs this same cap separately at a requested tolerance of 1e-10, giving 386 regular cases. The relaxed case has a distinct `_relaxed_1e-10` identifier. The strict case is retained and can be run explicitly:

```sh
python run.py --stress
python run.py --selection original
```

`--stress` runs the strict case alone. `--selection original` runs the unchanged original 386-case matrix. A timeout remains a reported failure and produces a nonzero exit status; it is never counted as a passing test. The initial relaxed run completed in approximately 0.31 seconds for each implementation, with native volume error below 1.5e-21 and surface-moment error below 1.2e-18.

## Focused checks

`run_focused.py` runs API contracts, output range rejection, standalone multiprecision, native coefficient agreement with actual FEMuS values and gradients, and direct QUAD/TRI/TET/HEX moment probes. QUAD and HEX probes include degree 10. `--extended` also runs the slower degree-8/10 TRI/TET checks, with a five-minute process budget. The build checks both planar/LS include orders and rejects accidental dependencies on the planar or Algoim headers from the standalone LS tests.

`compare_planar.py` checks 147 affine rules against the unchanged original `CutFemWeight` class, including tetrahedral base selections 0, 1, and 2. The tiny-cap comparison uses the explicitly named relaxed result when the strict result is unavailable.

The exact original-commit snapshots, original/candidate comparison files, and optional direct Algoim application check are external verification artifacts rather than prerequisites for this candidate-only suite.
