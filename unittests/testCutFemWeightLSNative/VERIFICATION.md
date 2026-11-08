# CutFemWeight_LS native reference verification

## Outcome

The final implementation passes the unchanged 386-case numerical matrix, all eight required supplemental cases, all 13 focused executable runs, and all 147 comparisons with the unchanged planar class. The main and supplemental matrices contain 40,036 numerical, finite-element, and cache assertions. There are no failed cases, integration exceptions, or timeouts in these final results.

The original degree-5 tetrahedral tiny-cap request now passes at its original absolute and relative tolerances of 1e-13. The scaled degree-5 plane passes at both 1e-12 and 1e-13. The separately named relaxed cap remains supplementary evidence and never replaces a strict request.

The final complete build contains 16 test translation units and uses 19 actual FEMuS finite-element and Gauss-table object sources. Every build, final numerical result, and isolated timing reported below uses the same frozen implementation. The complete SHA-256 manifest is `build/source_manifest.json`; the two headers changed during the defect correction have these final hashes:

| Header | SHA-256 |
| --- | --- |
| `AdaptiveLineSliceMoments.hpp` | `cdc7b416804a68d6ba0ec43bc473418a2b45c45879511f8574fccbf1bba554f8` |
| `Tet15LevelSetMoments.hpp` | `20e407c6601f3ce8f497b31dd3446d15db7566494e38464e331c7abc9fcff7a4` |

The original comparison snapshot is commit `cb19f37a14a7601b103da176bc44e30195d6f667`. The previously delivered native implementation was also preserved from the git index before any defect edits. Earlier failed attempts, reports, and timings are retained in explicitly labeled historical directories; they are not counted as final passes.

## Regressions for the demonstrated defects

`test_tet_slice_regression.cpp` reconstructs the diagnosed first triangle slice from actual FEMuS TET15 nodes. It retains the small enriched coefficients, checks every degree-5 triangle moment against an independent geometric formula, and tests corrected intervals, final output rounding, and atomic failure behavior.

The exact same test source was compiled against the preserved pre-fix headers and the final headers. Its SHA-256 is `6940528cb8c5e66b4a272ae572bbe92396fb697546f6f6b1224fefd821b7f23a`.

| Same-source direct regression | Before correction | Final implementation |
| --- | --- | --- |
| Overall result | 612 checks, 33 failures | 530 checks, no failures |
| Invalid corrected intervals | Reversed and zero-width intervals accepted | Explicit rejection |
| Ill-conditioned first chart | False zero interface mass returned | Explicit rejection before returning moments |
| Maximum observed moment error | 1.95304e-3 | 9.25638e-18 |

The check counts differ because the final implementation rejects the invalid chart before returning the moment table. A completed chart must return nonzero, accurate interface mass; the test does not accept a false zero as an alternative to an exception. Valid one-ULP intervals remain positive. Regular charts and chart selection pass every requested degree-5 moment at the original final tolerance.

The separate `test_root_candidate_regression.cpp` passes 11 checks. It verifies that duplicate candidates for the same simple root are consolidated while retaining a genuine boundary strip of width 1.83173e-20, distinct roots separated by 1e-24 and 1.86265e-9, and separate polynomial provenance. The identical probe fails on the saved pre-consolidation header and passes on the final header. Consolidation does not weaken corrected-width or provenance checks.

The degree-0 output-cast regression uses a full tetrahedron with `G=x+1` and zero relative tolerance. At an absolute tolerance of 8e-18 the double output is rejected without modifying either output table. At 2e-17 it succeeds, and the returned error estimate includes the actual double rounding of 1/6. This protects the final requested budget after the native moment transformation and output cast.

The nested TET correction distinguishes working allocations from the final requested tolerance. Bounded values, their error estimates, and cut state stay together; an unattainable inner working allocation does not discard a usable bounded sample. The outer integrator still checks the full requested tolerance and explicitly rejects an irreducible inner error that alone exceeds that tolerance.

Evidence: `defect_fix_before/test_tet_slice_regression_before_results.json`, `defect_fix_after/test_tet_slice_regression_after_results.json`, `defect_fix_after/test_root_candidate_regression_after_results.json`, and the preserved root-candidate diagnosis records.

## Main numerical matrix

| Set | Cases | Assertions | Result |
| --- | ---: | ---: | --- |
| Unchanged original matrix | 386 | 36,860 | All pass |
| Required supplemental matrix | 8 | 3,176 | All pass |
| Combined final class tests | 394 | 40,036 | All pass |
| Preserved original implementation | 386 | 36,463 | 385 pass; strict tiny cap timed out |

The original case definitions and requested tolerances are unchanged. The eight supplemental requests comprise six strict tiny-cap variants covering level-set scale, sign, and affine physical geometry; the separately named relaxed cap; and the scaled degree-5 plane at 1e-13 alongside its original-matrix counterpart at 1e-12.

The class matrix covers LINE3, QUAD9, TRI7, HEX27, and TET15 with fitting degrees 0 through 5. Fields include affine axis and oblique cuts, monotone quadratics, coupled tensor polynomials, enriched simplex fields, circles and spheres, two roots on a line, full and empty cells, boundary interfaces, near-affine fields, tiny cuts, sign reversal, and scaling by powers of two. Geometry maps include identity, affine anisotropic scaling, curved shear, and nonconstant-Jacobian stretching.

Each case records positive volume, raw delta, reference surface with the norm inside and outside the moments, physical surface with the cofactor norm inside and outside, and physical volume. Every native monomial through the fitting degree is evaluated. Repeated cache requests must return exactly the same vectors.

The main matrix instantiates `CutFemWeight_LS<double,long double>`. Separate checks instantiate `cpp_bin_float_oct` as the arithmetic type and as the output type.

Evidence: `candidate_original_summary.json`, `candidate_supplemental_summary.json`, and `results/candidate/`.

## Native coordinates and fitting space

Monomials, G, and grad G are evaluated on the standard FEMuS reference cells: [-1,1] for LINE, [-1,1]^2 for QUAD, [-1,1]^3 for HEX, and the standard nonnegative unit simplexes for TRI and TET. The fitting space remains the full total-degree space P_qM. FEMuS `Basis.hpp`, the original finite-element definitions, the planar class, and the Gauss tables remain unchanged.

The verified correction is

\[
c(x_q)=m(x_q)^T M^{-1}f,\qquad M_{ab}=\int_K m_a m_b.
\]

Native tensor-cell monomial integrals vanish for an odd exponent and otherwise factor into products of `2/(exponent+1)`. Standard-simplex monomial integrals are `product(exponent!)/(dimension+sum(exponent))!`. These determine the native Gram matrix without a change to the double-valued FEMuS shape-function API.

Volume moments use `f_a = integral_K m_a H(G)`. Raw-delta moments use `integral_K m_a delta(G)`. Reference-surface moments include `|grad G|`; physical-surface moments include `|cof(J) grad G|`. The returned corrections are multiplied by the ordinary reference Gauss weights. Physical surface integration does not introduce a second determinant factor.

The independent driver reconstructs G, grad G, geometry, and its Jacobian directly with the actual FEMuS finite-element basis. The separate coefficient test makes 15,360 such comparisons and has maximum absolute error 7.13492e-15.

## Analytic references and error policy

Independent Python references use exact interval, simplex, ball, and sphere formulas where available. Other curved cuts use independently derived graph roots and Gauss-Legendre quadrature at 40 and 72 points per parameter direction. Those references do not call the implementation's integration kernels. The largest discrepancy between the two reference orders is 7.10543e-15.

The normalized error is the absolute moment error divided by `max(1, |zeroth integral|, |target moment|)`. This distinguishes a raw delta integral magnified by level-set scaling from an unscaled geometric surface integral. The general class-test allowance is `20*requested_tolerance + 3e-11`, allowing for double Gauss tables and the fitting calculation. The strict cap and its strict scaled variants use the tighter limit `2*requested_tolerance`, namely 2e-13. Direct kernel tests additionally assess their reported moment error estimates against the requested budgets.

| Accuracy measurement | Final value |
| --- | ---: |
| Largest asserted normalized error in the original 386 cases | 1.18794e-13 |
| Largest asserted normalized error in all 394 cases | 1.45661e-13 |
| Largest fraction of a class-test allowance used | 0.72831 |
| Unscaled strict-cap maximum volume error | 1.48562e-21 |
| Unscaled strict-cap maximum interface error across both formulations | 1.11385e-18 |

The largest normalized error across all cases occurs in the strict cap scaled by 2^-20. Its raw delta mass is approximately 8; the largest absolute raw delta moment error is 1.16529e-12. The corresponding geometric surface errors remain about 1.11e-18. These are reconstructed class moments using double Gauss data; the table does not assert exact arithmetic or a universal equality between an internal estimate and final floating-point reconstruction error.

### The two surface formulations

A raw delta correction multiplied by the pointwise FEMuS gradient or cofactor norm is asserted equivalent to including that norm inside the moments only when the complete integrand belongs to the fitting space. For a nonpolynomial norm, the outside formulation is recorded as a finite-degree approximation against the independent surface reference. Variable-Jacobian volume tests apply the corresponding degree rule to the determinant.

For the degree-5 coupled QUAD graph, the reference interface length is 2.049993397776071. Including the norm inside gives maximum weighted surface-moment error 1.63203e-14. Multiplying the norm outside gives a constant surface estimate differing by 2.54880e-5. For the degree-5 circle, the outside estimate has approximately 4.054% relative error, while the inside surface moments agree within 1.12799e-13. These approximation differences are recorded openly and are not turned into inappropriate exactness assertions.

Evidence: `defect_fix_after/final_fixed_results.json` and the individual case assessments.

## Original and final numerical equivalence

The comparison was regenerated from the final fixed records and the preserved original records. It contains 385 successful pairs, 70,288 correction values, and 38,598 integrated moments. The original strict tiny-cap timeout is excluded from paired comparison because it has no successful original result; the final case passes independently.

| Quantity | Maximum normalized original/final moment difference |
| --- | ---: |
| Positive volume | 6.18522e-14 |
| Raw delta | 1.52611e-12 |
| Reference surface, norm inside | 3.81528e-13 |
| Reference surface estimate, norm outside | 4.34874e-13 |
| Physical surface, cofactor norm inside | 3.81500e-13 |
| Physical surface estimate, cofactor norm outside | 4.34847e-13 |
| Physical volume | 1.01030e-13 |

The largest normalized difference in an individual correction is 1.21152e-10, for the degree-5 TET zero-face volume rule. Integrated moments are the primary equivalence measurement because the polynomial basis change and floating-point cancellation can change individual corrections more than their integrated moments.

The largest absolute raw delta moment difference is 1.61119e-7 for the scaled HEX plane. Its zeroth raw delta integral is approximately 4,194,304, so the difference is about 3.85e-14 relative. Both implementations also pass their independent reference checks on every successful pair.

The separate planar comparison passes all 147 affine rules, including tetrahedral bases 0, 1, and 2. It checks 16,258 correction values and 4,552 moments; the largest normalized moment difference is 1.56783e-13. All three tiny-cap tetrahedral comparisons use the successful native request at the original tolerance of 1e-13.

Evidence: `baseline_candidate_comparison.json` and `planar_comparison_results.json`.

## Complete-request timings

These final timings were collected one request at a time after compilation and preceding checks had finished. Each includes the complete class sequence: raw delta, reference surface, physical identity surface, correction construction, and a cached repeat. Each request passes 397 assertions.

| Request | Requested absolute and relative tolerance | Final isolated time |
| --- | ---: | ---: |
| TET degree 5, `G=x-(1-1/256)` | 1e-13 | 0.320 s |
| TET degree 5, `G=2^-20*(x-0.25)` | 1e-12 | 0.493 s |
| Same scaled plane | 1e-13 | 0.620 s |

The preserved original and previous native implementations exceeded the 45-second limit for the strict cap. The previous native scaled-plane request at 1e-12 took 27.45 seconds during concurrent verification. That older timing is diagnostic rather than a controlled benchmark, so it is not used to claim a precise speedup factor. The final timings are complete class requests, distinct from the much smaller direct-kernel timings in the diagnostic records.

The intermediate 1e-13 success and 1e-12 failure were different requested tolerances. A later preflight also exposed an explicit 1e-13 root-bracket rejection. Both exact requests remain required final cases and now pass; neither was removed or relaxed.

Evidence: `defect_fix_after/strict_cap_full_class_isolated.json`, `scaled_plane_full_class_isolated.json`, and `scaled_plane_1e13_full_class_isolated.json`.

## Additional completed checks and scope

| Check | Final result |
| --- | --- |
| Standalone compilation | All 16 translation units pass; 19 actual FE/Gauss object sources |
| LS dependency isolation | No planar, old cut-FEM, or Algoim headers in standalone LS dependencies |
| Include orders | Both orders compile; no planar inheritance or conversion |
| API contracts | 157 runtime checks; 29 removed signatures inaccessible |
| Native G and grad G | 15,360 comparisons; maximum absolute error 7.13492e-15 |
| QUAD/TRI direct moments | 2,654 checks in 128 cases; maximum error 4.10262e-16 |
| Degree-10 QUAD quadrants | 8,448 checks; maximum error 7.40510e-17 |
| Scaled QUAD plane and graph | 216 checks; maximum error 0.00183529 of the requested budget |
| TET directions and metrics | 640 checks in 28 cases; maximum error 1.04409e-16 |
| Degree-8/10 TRI/TET kernels | 3,372 checks; maximum error 9.25638e-18 |
| HEX degree-5/10 kernels and geometry | All probes pass; plane and curved-geometry errors below 3.81e-16 |
| First-slice and root-candidate regressions | 530 and 11 checks pass |
| Multiprecision and output range | Arithmetic/output instantiation, atomic rejection, recovery, and cache checks pass |

The direct degree-8/10 results exercise moment kernels. The broad class matrix, including Gram inversion, correction construction, and actual Gauss evaluation, covers degrees 0 through 5. Kernel results are not a degree-10 end-to-end fitting guarantee.

An earlier separate C++17 direct Algoim caller check passed eight analytic volume/delta comparisons across QUAD, TRI, HEX, and TET, with maximum error below 4.45e-16. That caller has no LS dependency and was unchanged by these numerical fixes; this is retained prior evidence, not an extra final LS build.

The suite uses actual FEMuS sources without PETSc. It does not claim a full PETSc-based application build, performance guarantees on other machines, or universal accuracy for nonpolynomial outside-norm estimates.

## Reproduction and evidence

The accompanying `README.md` gives portable compiler, Boost, Eigen, Python, and NumPy requirements and commands. The final required sequence is:

```sh
python build.py --boost /path/to/boost/include --eigen /path/to/eigen3
python run.py --selection original
python run.py --selection supplemental
python run_focused.py --extended
python compare_planar.py
```

Runners exit nonzero for failed assertions, integration exceptions, nonfinite results, missing executables, or timeouts. The original strict cap is selectable with `python run.py --stress`. Its relaxed counterpart is always separately named.

Current aggregate evidence is `build_results.json`, `candidate_original_summary.json`, `candidate_supplemental_summary.json`, `extended_results.json`, `planar_comparison_results.json`, and the regenerated `baseline_candidate_comparison.json`. `build/source_manifest.json` and `build/test_source_manifest.json` identify the implementation and test sources. Per-case JSON retains executable arguments, timings, correction arrays, weighted moments, independent finite-element checks, approximation records, and failure details. The external `defect_fix_after/final_fixed_results.json` cross-checks all final gates and ties the same-source regressions and isolated timings to the final manifest.
