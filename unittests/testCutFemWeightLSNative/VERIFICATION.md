# CutFemWeight_LS native reference verification

## Outcome

The final native-reference implementation passes all 386 cases in the explicitly documented regular set, with 36,860 numerical, finite-element, and cache checks. The unchanged original comparison matrix contains 386 cases: 385 pass with 36,463 checks, and its strict degree-5 tetrahedral tiny-cap case exceeds the 45-second process budget in both the original and final implementations. The same cap at a requested tolerance of 1e-10 passes 397 additional checks in 0.31 seconds. The strict timeout is retained as a failed stress result and is not counted as a pass.

Seven scaling failures found in the first candidate sweep were corrected and all pass in the final sweep. Standalone compilation, API contracts, multiprecision, native coefficient extraction, focused element moments, extended degree-8/10 simplex checks, and the separate original planar-class comparisons pass.

The exact original source is commit `cb19f37a14a7601b103da176bc44e30195d6f667`. The final header SHA-256 manifest is `build/sources_candidate.json`; the same manifest is saved with the portable suite. The delivered headers match the complete validated source snapshot byte for byte. The full 14-translation-unit build, original numerical matrix, focused checks, and planar comparisons reported here all refer to that matching implementation. The saved source manifest and the `before_tet_budget` snapshot retain its exact provenance; discarded experimental results are excluded from the final numerical results.

## Native coordinates and monomial formulation

The tests evaluate monomials, G, and grad G directly on the standard FEMuS cells: [-1,1] for LINE, [-1,1]^2 for QUAD, [-1,1]^3 for HEX, and the standard nonnegative unit simplexes for TRI and TET. The fitting space remains the full total-degree space P_qM. The unchanged FEMuS `Basis.hpp` and actual element definitions provide the independent double-precision interpolation checks.

The verified correction has the native-reference formulation `c(x_q) = m(x_q)^T M^-1 f`, with `M_ab = integral_K m_a m_b`. Volume moments use `f_a = integral_K m_a H(G)`. Raw-delta moments use `integral_K m_a delta(G)`. Reference-surface moments include `|grad G|` inside that integral; physical-surface moments include `|cof(J) grad G|`. Ordinary reference Gauss weights multiply the returned correction. The original/final and planar comparisons verify that this coordinate and basis change preserves the fitted integrals without retaining the old planar integration dependency.

## Main numerical matrix

| Set | Cases | Checks | Result |
| --- | ---: | ---: | --- |
| Original implementation, unchanged comparison matrix | 386 | 36,463 | 385 pass; one strict tiny-cap timeout |
| Final implementation, unchanged comparison matrix | 386 | 36,463 | 385 pass; the same strict tiny-cap timeout |
| Explicit relaxed tiny-cap follow-up | 1 | 397 | Pass |
| Documented regular set: 385 original cases plus the named relaxed cap | 386 | 36,860 | All pass |

The matrix covers LINE3, QUAD9, TRI7, HEX27, and TET15, fitting degrees 0 through 5, and raw delta, reference-surface, and physical-surface modes. Cases include affine axis and oblique cuts, monotone quadratics, coupled tensor polynomials, enriched simplex fields, circles and spheres, two roots on a line, full/empty cells, boundary interfaces, near-affine fields, tiny cuts, sign reversal, and level-set scaling by powers of two. Geometry maps include identity, affine anisotropic scaling, curved shear, and nonconstant-Jacobian stretching.

The main sweep instantiates `CutFemWeight_LS<double,long double>`. Separate tests instantiate `cpp_bin_float_oct` as the arithmetic type and as the output type.

For each case the C++ driver evaluates the actual unchanged FEMuS shape functions and their derivatives to independently reconstruct G, grad G, the physical geometry, and its Jacobian. Every fitted monomial through degree qM is integrated. Cache repeats must return exactly the same vectors.

Independent Python references use exact interval, simplex, ball, and sphere moments where available. Other curved cuts use independently derived graph roots and Gauss-Legendre quadrature with 40 and 72 points per parameter direction. The largest discrepancy between those two reference quadratures was 7.11e-15. The largest asserted normalized error of the final implementation against these references was 1.55e-12, in the TET15 degree-2 sphere case, whose requested integration tolerance was 1e-9. The largest ratio of any asserted error to its test limit was 0.00199.

The test normalization is `max(1, |zeroth integral|, |target moment|)`. This keeps raw delta results comparable when scaling G changes delta by factors near one million. Final Gauss weights and coordinates are the actual FEMuS double tables.

## Original/final numerical equivalence

For the 385 successful original/final pairs, the comparison checks 70,288 correction values and 38,598 integrated moments. Every pair also passes its independent reference checks.

| Quantity | Maximum normalized original/final difference |
| --- | ---: |
| Positive-volume moment | 6.19e-14 |
| Raw-delta moment | 1.32e-13 |
| Reference-surface moment, norm inside | 1.12e-13 |
| Reference-surface estimate, norm outside | 1.12e-13 |
| Physical-surface moment, cofactor norm inside | 1.47e-13 |
| Physical-surface estimate, cofactor norm outside | 1.50e-13 |
| Physical-volume moment | 1.02e-13 |

The largest normalized difference in an individual correction value is 1.22e-10, in the degree-5 TET15 zero-face case. Integrated moments are the primary equivalence measure because a change of polynomial basis can alter cancellation in individual correction values.

The largest absolute raw-delta moment difference is 1.61e-7 for the scaled HEX plane. Its zeroth raw-delta integral is approximately 4,194,304, so that difference is approximately 3.85e-14 relative. This is not an unscaled surface error.

The separate original planar class comparison passes 147 affine rules, including tetrahedral bases 0, 1, and 2: 16,258 correction comparisons and 4,552 moment comparisons. Its maximum normalized moment difference is 1.57e-13. The three tiny-cap planar comparisons explicitly use the successful relaxed native request; their records include the actual candidate case identifier and requested tolerance.

## The two surface formulations

The suite checks both the raw-delta correction multiplied by a pointwise FEMuS gradient/cofactor norm and the correction obtained with that norm inside the moment integral. It demands equality only when the complete integrand is in the fitted polynomial space. Nonpolynomial norms multiplied outside the projection are recorded as approximations against the independent surface reference.

For example, on the degree-5 coupled QUAD graph, the reference interface length is approximately 2.049993397776071. The norm-inside rule has maximum weighted surface-moment error 1.64e-14. The outside formulation's constant surface estimate differs by approximately 2.55e-5. For the degree-5 circle, the outside formulation has about 4.05% relative surface error while the norm-inside surface moments agree to about 1.13e-13. These finite-degree differences are expected projection behavior; the tests do not hide them by asserting an inappropriate equivalence.

Physical surface tests use the cofactor metric and reference Gauss weights, without multiplying a second geometry determinant. Variable-Jacobian volume checks account for the degree of the determinant in the integrand.

## Additional completed checks

| Check | Result |
| --- | --- |
| Standalone portable compilation | 14 test translation units; 19 actual FE/Gauss object sources; all pass |
| LS dependency isolation | No planar, old cut-FEM, or Algoim headers in standalone LS dependency lists |
| Both include orders | Compile successfully; LS has no planar inheritance/conversion |
| API contracts | 157 runtime checks; 29 removed signatures inaccessible |
| Native G/grad G extraction | 15,360 comparisons against actual FEMuS evaluations; max error 7.14e-15 |
| QUAD/TRI direct native moments | 2,654 checks in 128 cases; max error 4.11e-16 |
| Degree-10 QUAD quadrant integration | 8,448 checks; max error 7.41e-17 |
| Scaled QUAD plane and curved graph | 216 checks; error below 0.00184 of requested per-moment budget |
| TET native moments/directions/metrics | 640 checks in 28 cases; max error 1.05e-16 |
| Degree-8/10 TRI/TET probes | 3,372 checks; max error 5.66e-17 |
| HEX degree-5/10 planes and curved geometry | All probes pass; reported plane/geometry errors below 3.81e-16 |
| Multiprecision and output range | Standalone cpp_bin_float_oct arithmetic, wider output type, atomic overflow rejection, recovery, and cache checks pass |
| Optional direct Algoim application caller | Eight analytic volume/delta checks across QUAD/TRI/HEX/TET; max error 4.45e-16 |

The optional Algoim caller check was built in C++17 separately; it is not a dependency of the portable LS suite. The standalone suite does not claim a full PETSc-based FEMuS application build.

The degree-8/10 probes above exercise the direct moment kernels. The broad end-to-end class tests, including Gram inversion, correction construction, and actual FEMuS Gauss evaluation, cover degrees 0 through 5. The direct high-degree results are not presented as a degree-10 end-to-end fitting guarantee.

## Remaining limits and timing observation

The strict stress case is TET15 degree 5 with `G = x - (1 - 1/256)`, identity geometry, and requested absolute and relative tolerances 1e-13. The 45-second budget covers raw delta, reference surface, physical identity surface, and a cached repeat. Both original and final implementations exceed that budget. At tolerance 1e-10, the final result completes in 0.31 seconds with volume error below 1.49e-21 and interface-moment error below 1.12e-18. The strict case remains available through `run.py --stress`, and `--selection original` reproduces the exact original matrix with a nonzero exit status when it times out.

A separate accuracy-passing performance regression remains for TET15 degree 5 with `G = 2^-20 * (x - 0.25)`, identity geometry, auto direction, validation 0, and tolerances 1e-12. The measured complete request sequence takes 27.45 seconds in the final implementation versus 0.54 seconds in the original. The final unscaled case takes 0.55 seconds. These timings were observed during concurrent verification and are diagnostic measurements rather than controlled benchmarks; individual mode times were not recorded. The exact executable arguments are `4 5 plane identity 9.5367431640625e-07 1e-12 0`.

## Reproduction and artifacts

The portable suite is `MyFEMuS_LS_native/unittests/testCutFemWeightLSNative`. Its README provides compiler, Boost, Eigen, and NumPy requirements and full commands. It builds actual FE/Gauss sources without PETSc, keeps numerical assertions separate from recorded outside-norm approximations, and returns nonzero for unexpected failures. The slower degree-8/10 simplex probe is available with `run_focused.py --extended`.

Key machine-readable artifacts in this verification directory are `baseline_all_summary.json`, `candidate_all_summary.json`, `candidate_regular_summary.json`, `baseline_candidate_comparison.json`, `planar_comparison.json`, `extended_results.json`, `build_results.json`, and `build/sources_candidate.json`. Per-case JSON retains the command, timing, generated correction arrays, weighted moments, independent FE errors, and all failure details. The first failing candidate results remain separately under `results/candidate_initial`; final results are under `results/candidate`.
