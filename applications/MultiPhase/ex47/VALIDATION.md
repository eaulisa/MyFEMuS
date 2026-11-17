# ex47 validation

The selective-measure revision starts from `67f00f5edd5e2de1ca61cba0a54d83eb7b3e75ea` on `barbi-ex47-levelset-data`, which includes the user formatting commit `26cbb141a881550597e0106ba9d0691eddbfad04`. The original prototype derives from `barbi` commit `ffa844b5c0742a15962dd0cb06ccb30d2946e13b`. All executable checks are in `ex47.cpp` and run from its main. The production CutFEM implementation is unchanged.

Environment: GCC 14.2.0, Boost 1.86.0, Release builds plus a Debug AddressSanitizer/UndefinedBehaviorSanitizer build; the class and main pass explicit C++14 syntax checks. Leak detection was disabled for the sanitizer run because of the execution environment.

| Run | Result | Maximum scaled error | Elapsed time |
| --- | ---: | ---: | ---: |
| Default, 64 cases and 32 points | 637,146 checks passed | 7.47e-15 | 0.172 s |
| Stress, 256 cases and 64 points, with benchmark | 4,752,864 checks passed | 1.44e-14 | 1.179 s |
| ASan/UBSan, default cases and points | 637,146 checks passed; no sanitizer errors | 7.47e-15 | 2.205 s |
| CMake parent-project `add_subdirectory`, 2 cases and 3 points | 6,694 checks passed | 1.09e-15 | 0.025 s |

Scaled errors are `abs(actual-expected)/(1+abs(expected))`. Timings are measurements on this runtime, not production assembly benchmarks. The full MPI/PETSc application build was not run; standalone and parent-project builds of the new target were checked.

Coverage includes:

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

- `include/CutFemLevelSetData.hpp`: `b7b597e76befbad88d75a8b9a46b2362a7a3c5d2ab1b4b8d112e4387bf1c3f6e`
- `ex47.cpp`: `90ef7489c330f86d6573fd01a44015b1d6f260ec0eb1ed6be91615ca3c156c6d`

See [README](README.md) for reproducible build/run commands, [PERFORMANCE](PERFORMANCE.md) for the five-run kernel comparison, and [TODO](TODO.md) for the remaining production integration work.
