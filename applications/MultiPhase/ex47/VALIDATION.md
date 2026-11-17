# ex47 validation

Based on the latest `barbi` commit `ffa844b5c0742a15962dd0cb06ccb30d2946e13b`, on the separate branch `barbi-ex47-levelset-data`. All executable checks are in `ex47.cpp` and run from its main. The production CutFEM implementation is unchanged.

Environment: GCC 14.2.0, C++14, Boost 1.86.0, Release builds plus a Debug AddressSanitizer/UndefinedBehaviorSanitizer build. Leak detection was disabled for the sanitizer run because of the execution environment.

| Run | Result | Maximum scaled error | Elapsed time |
| --- | ---: | ---: | ---: |
| Default, 64 cases and 32 points | 588,975 checks passed | 7.47e-15 | 0.19 s |
| Stress, 256 cases and 64 points, with benchmark | 4,561,461 checks passed | 1.44e-14 | 1.15 s |
| ASan/UBSan, default cases and points | 588,975 checks passed; no sanitizer errors | 7.47e-15 | 2.57 s |
| CMake parent-project `add_subdirectory`, 2 cases and 3 points | 4,775 checks passed | 1.09e-15 | 0.013 s |

Scaled errors are `abs(actual-expected)/(1+abs(expected))`. Timings are measurements on this runtime, not production assembly benchmarks. The full MPI/PETSc application build was not run; standalone and parent-project builds of the new target were checked.

Coverage includes:

- LINE3, QUAD9, TRI7, HEX27, TET15, and WEDGE21 reconstruction and first/second derivatives against the actual FEMuS basis implementations and node ordering.
- Bernstein/power agreement, bounds, normalization, nodal-first detection, configurable depth, both phases, small interior cuts with uniform nodal signs, tangencies, and inconclusive control nets.
- Every slice direction, recursive 3D-to-2D-to-1D restriction, exact reference maps and measure factors, endpoint/collapsed slices, direct Bernstein contraction, and retained original gradient components.
- Arbitrary external polynomial components through coordinate degree 6, including preservation of three components through dimension reduction.
- Owned nodal/slice snapshots, independent copies, reset/context/depth changes, derivative/Hessian invalidation, and reuse of already allocated child storage without parent coefficient reconstruction.
- Input validation, nonfinite values, invalid layouts/coordinates/depths, aliasing, and one-dimensional termination.
- Multiprecision root references, asymmetric roots, near tangencies, cancellation, inherited slice roundoff margins, and input scales from 1e-200 to 1e200.
- Nonuniform GK15 slice sweeps and nested volume/interface integration against analytic half-domain and planar-interface measures for all six elements.
- Physical circle/sphere pullbacks on affine and non-affine bilinear/trilinear QUAD/HEX geometry, including analytic gradients and slices.

The double/long-double evaluation comparison uses a quadruple-precision reference (`cpp_bin_float_quad`) for the retained polynomial. It measures the arithmetic kernel, not an independent reconstruction of the supplied nodal data.

| Element | Double maximum error | Long-double maximum error | Double time | Long-double time |
| --- | ---: | ---: | ---: | ---: |
| LINE3 | 1.89e-16 | 5.33e-20 | 0.00042 s | 0.00115 s |
| QUAD9 | 1.11e-16 | 5.80e-20 | 0.00054 s | 0.00208 s |
| TRI7 | 5.11e-16 | 2.26e-19 | 0.00043 s | 0.00219 s |
| HEX27 | 2.65e-16 | 1.07e-19 | 0.00102 s | 0.00491 s |
| TET15 | 1.55e-15 | 4.80e-19 | 0.00074 s | 0.00462 s |
| WEDGE21 | 7.52e-16 | 4.64e-19 | 0.00085 s | 0.00461 s |

Each timing covers 128,000 evaluations. Routine Horner kernels use long-double multiply/add; the sensitive root discriminant retains compensated fused arithmetic. These results do not justify changing production precision or predict production integration cost.

The new main also passes `-Wall -Wextra -Wpedantic -Werror` syntax checking with existing FEMuS/Boost headers treated as system headers. The runner passes Bash syntax checking and its standalone launch; CLI help and four invalid-argument cases were checked. `git diff --check` passes.

Source SHA-256:

- `include/CutFemLevelSetData.hpp`: `1716025d0f7cee9644543269ae6019b18022abf783f4666a960c3adab00ad3b5`
- `ex47.cpp`: `fa6f290bb49de35b6837a04fbbde41e1a23cecfa602872101caaf74a9a3c2bf3`

See [README](README.md) for reproducible build/run commands and [TODO](TODO.md) for the remaining production integration work.
