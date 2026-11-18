# ex47 validation

The prepared-topology/chart revision starts from `2b2dc6cc95ba8787d7c1ffa43c3d5fc460a55ce6` on `barbi-ex47-levelset-data`, after the selective-measure revision and the user formatting commit `26cbb141a881550597e0106ba9d0691eddbfad04`. The original prototype derives from `barbi` commit `ffa844b5c0742a15962dd0cb06ccb30d2946e13b`. All executable checks are in `ex47.cpp` and run from its main. The production CutFEM implementation is unchanged.

Environment: GCC 14.2.0, Boost 1.86.0, Release builds plus a Debug AddressSanitizer/UndefinedBehaviorSanitizer build; the class and main pass explicit C++14 syntax checks. Leak detection was disabled for the sanitizer run because of the execution environment.

| Run | Result | Maximum scaled error | Elapsed time |
| --- | ---: | ---: | ---: |
| Default, 64 cases and 32 points | 804,228 checks passed | 7.47e-15 | 0.330 s |
| Stress, 256 cases and 64 points, with benchmark | 4,939,947 checks passed | 1.44e-14 | 1.355 s |
| ASan/UBSan, default cases and points | 804,228 checks passed; no sanitizer errors | 7.47e-15 | 6.531 s |

Scaled errors are `abs(actual-expected)/(1+abs(expected))`. Timings are measurements on this runtime, not production assembly benchmarks. The full MPI/PETSc application build was not run; this revision checks the standalone target. The unchanged CMake target was also checked through a parent-project build in the preceding revision.

Coverage includes:

- Complete fixed-size LINE topology, stable signed slopes, multiplicities, endpoints, positive intervals, explicit zero fibers, and compatibility/root-record surface methods.
- All supported axis permutations and all 6 TRI/WEDGE or 24 TET barycentric permutations, root-component ordering, recursive restrictions, elevated swapped WEDGE slices, and original-gradient omission.
- Exact nodal zero certificates and restriction on every boundary, including sloping simplex faces and scales 1e-200 through 1e200; nonzero values down to 1e-300 remain nonzero.
- Full-degree external polynomial chart composition, increased coordinate degree, and explicit rejection rather than truncation above degree 32.
- Cached projection polynomials, compensated/shifted discriminants, boundary quartics, unit exports, tighter derivative bounds, strict/weak directions, zero-fiber exclusions, dependence/parity, and cache detachment on copies.
- LINE3, QUAD9, TRI7, HEX27, TET15, and WEDGE21 reconstruction and first/second derivatives against the actual FEMuS basis implementations and node ordering.
- Bernstein/power agreement, bounds, normalization, nodal-first detection, configurable depth, both phases, small interior cuts with uniform nodal signs, tangencies, and inconclusive control nets.
- Every slice direction, recursive 3D-to-2D-to-1D restriction, exact reference maps and measure factors, endpoint/collapsed slices, direct Bernstein contraction, and retained original gradient components.
- Arbitrary external polynomial components through coordinate degree 6, including preservation of three components through dimension reduction.
- Owned nodal/slice snapshots, independent copies, reset/context/depth changes, derivative/Hessian/graph invalidation, actual parent destruction, and reuse of already allocated child storage without parent coefficient reconstruction.
- Input validation, nonfinite values, invalid layouts/coordinates/depths, aliasing, and one-dimensional termination.
- Multiprecision root references, asymmetric roots, near tangencies, cancellation, inherited slice roundoff margins, and input scales from 1e-200 to 1e200.
- Nonuniform GK15 slice sweeps and nested volume/interface integration against analytic half-domain and planar-interface measures for all six elements.
- Physical circle/sphere pullbacks on affine and non-affine bilinear/trilinear QUAD/HEX geometry, including analytic gradients and slices.

- Fixed-output gradients, line roots and scalar evaluation, with compatibility-interface comparisons.
- All four scalar and all signed normal surface densities against analytic Jacobians, with prepared cofactor-gradient polynomial restriction and matrix-input comparisons; affine/non-affine geometry on every element and bi/triquadratic QUAD/HEX geometry.
- Zero steady-state allocations in the prepared quadrature and selective-surface paths under the dedicated counting build.
- Exact equality between each selected surface density and explicit batch evaluation, including cofactor-matrix variants and individual/shared normal components.
- Omitted original-gradient context through recursive slices, explicit missing-context errors, restoration from a retained ancestor, identical roots/classification/maps/local derivatives, and snapshot ownership.
- Multiprecision delta/inverse-metric comparisons and scale invariance for small cuts with radii down to 1e-12 across input scales 1e-200 through 1e200.

The double/long-double evaluation comparison uses a quadruple-precision reference (`cpp_bin_float_quad`) for the retained polynomial. It measures the arithmetic kernel, not an independent reconstruction of the supplied nodal data.

| Element | Double maximum error | Long-double maximum error | Double time | Long-double time |
| --- | ---: | ---: | ---: | ---: |
| LINE3 | 1.89e-16 | 5.33e-20 | 0.000445 s | 0.00175 s |
| QUAD9 | 1.11e-16 | 5.8e-20 | 0.000484 s | 0.00376 s |
| TRI7 | 5.11e-16 | 2.26e-19 | 0.000442 s | 0.00381 s |
| HEX27 | 2.65e-16 | 1.07e-19 | 0.00115 s | 0.00932 s |
| TET15 | 1.55e-15 | 4.8e-19 | 0.00114 s | 0.00919 s |
| WEDGE21 | 7.52e-16 | 4.64e-19 | 0.00114 s | 0.00972 s |

Each timing covers 128,000 evaluations. Routine Horner kernels use long-double multiply/add; the sensitive root discriminant retains compensated fused arithmetic. These results do not justify changing production precision or predict production integration cost.

The new main also passes `-Wall -Wextra -Wpedantic -Werror` syntax checking with existing FEMuS/Boost headers treated as system headers. The runner passes Bash syntax checking and its standalone launch; CLI help and four invalid-argument cases were checked. `git diff --check` passes.

Source SHA-256:

- `include/CutFemLevelSetData.hpp`: `a4dc3228ecad5c58e16192983c889ead7a71f0ccd55df18b3feafaccedfaf638`
- `include/CutFemLevelSetDataAnalysis.hpp`: `ba4ed2f569aba86d9d82e6976b23b77d28237a84950e47b5bb986e5a47807244`
- `include/CutFemLevelSetDataCharts.hpp`: `2c432a1a799283af69cf0c065a2f36a0cb1af1861dfb0df37b48b5d3e815c252`
- `ex47.cpp`: `c6bdb004996d9ea633dad9bf2ff4becd4b8e68cbff083ebf961a1bdc48d0f137`

See [README](README.md) for reproducible build/run commands, [PERFORMANCE](PERFORMANCE.md) for the five-run kernel comparison, and [TODO](TODO.md) for the remaining production integration work.

## Sensitive-root regressions

The new references use `cpp_bin_float_100` and the retained parent polynomial, not the intended unrounded analytic circle/sphere. Errors below are relative density errors, separate from the scaled-error statistic above.

| Regression | Maximum relative error |
| --- | ---: |
| Density from the previous rounded-root derivative, same test roots | 1.318e-11 |
| Density from retained discriminant slopes | 1.085e-19 |
| Last-ULP fold restrictions in QUAD, HEX, TET, and WEDGE | 1.882e-17 |

The tests also cover centered folds, vanished quadratic coefficients with a retained negative constant of order 1e-200, ancestor-based refinement of 3D slice discriminants, zero/repeated/endpoint roots, and translated LINE reference/physical measures whose exact gradient-to-slope ratio is one. Root-surface gradients preserve the directional chain rule, including oblique barycentric charts. These targeted regressions supplement the prior broad accuracy tests, which did not expose the rounded-root loss.

All executable tests remain in ex47 main. Allocation-tracking builds pass 27 data-reuse, 45 prepared-algebra, 61 existing-quadrature, and 130 selected-surface timing rows. Fixed-output warmed paths allocate zero C++ heap objects; the existing vector-returning wrappers retain their expected allocations. First preparation, copied-cache detachment after a state change, and the sloping-boundary convenience constructor are outside that steady-state claim. LeakSanitizer cannot run under this execution environment's tracing; AddressSanitizer and UndefinedBehaviorSanitizer run with leak detection disabled.


## Local weight integration

The production class remains unchanged. `Ex47::CutFemWeight_LS` and its local recursive drivers use prepared data from LINE through QUAD/TRI to HEX/TET/WEDGE. Run `--test-weights` to execute those stages in order; `--weight-dimension 1|2|3` selects a stage.

The completed Release suite passes 32,513 weight checks. The original data suite passes 4,939,941 checks with `--cases 256 --points 64`. The direct compensated contraction preserves the translated-root density error of 1.085e-19 and gives a maximum last-ULP fold density error of 1.920e-17. No source polynomial is approximated or truncated by this optimization.

Errors in the following table are `abs(actual-expected)/(1+abs(expected))`, maximized over the analytic volume and four scalar surface measures and every tested moment. LINE uses a 100-digit reference built from the supplied nodal doubles; ball cases use closed-form volume moments and independent angular surface integration for an explicit non-affine geometry.

| Element | Fitted degrees compared with production | Maximum prepared raw-moment error | Raw degree-13 error |
| --- | --- | ---: | ---: |
| LINE3 | 2, 5, 10 | 2.587e-19 | — |
| QUAD9 | 2, 5, 10 | 9.616e-18 | — |
| TRI7 | 2, 5, 10 | 1.478e-16 | — |
| HEX27 | 2, 5, 10 | 1.821e-17 | 1.821e-17 |
| TET15 | 2, 5 | 2.029e-16 | 2.029e-16 |
| WEDGE21 | 2, 5 | 7.062e-17 | 7.062e-17 |

Every direction number is exercised with enriched nodal fields and every scalar surface mode. Additional checks cover GK 3/7 and 15/31 rules, validation modes 0/1/2, geometry-cache invalidation, all three signed normals, the TET barycentric normal orientation, centered-sphere parity, and QUAD circles down to radius 1e-12. The preceding sensitive-root/data tests cover further cancellation, scale, chart, and boundary cases.

The fixed fitting operator and supplied Gauss tables impose a separate residual floor. Feeding independent exact moments through the production fitter reproduces the measured residual: for example, the largest LINE degree-10 residual is 1.185e-11, while its raw integration error is 2.587e-19. Other fitted ball cases have residuals between about 6e-15 and 1.1e-13. Four requests per LINE degree fail in the original class on a tiny two-root cut that the prepared path retains; those failures are counted separately and excluded from the old error maximum.

High-order simplex fitting is deliberately excluded from the default timed validation loop because constructing its unchanged multiprecision fitting matrix is expensive. Raw degree-13 integration is checked independently for every 3D element. All executable checks remain in the ex47 main; none changes production assembly.

The final AddressSanitizer/UndefinedBehaviorSanitizer build passes all 32,513 weight checks and 2,551,441 data checks (`--cases 128 --points 64`) with leak detection disabled for the tracing environment. The Release comparison uses GCC with `-O2 -DNDEBUG`; sanitizers use `-O1 -g0 -DNDEBUG`.
