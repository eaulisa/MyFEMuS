# ex47 validation

The speed revision starts from the exact user formatting commit `26cbb141a881550597e0106ba9d0691eddbfad04` on `barbi-ex47-levelset-data`; the original prototype derives from `barbi` commit `ffa844b5c0742a15962dd0cb06ccb30d2946e13b`. All executable checks are in `ex47.cpp` and run from its main. The production CutFEM implementation is unchanged.

Environment: GCC 14.2.0, Boost 1.86.0, Release builds plus a Debug AddressSanitizer/UndefinedBehaviorSanitizer build; the class and main pass explicit C++14 syntax checks. Leak detection was disabled for the sanitizer run because of the execution environment.

| Run | Result | Maximum scaled error | Elapsed time |
| --- | ---: | ---: | ---: |
| Default, 64 cases and 32 points | 612,534 checks passed | 7.47e-15 | 0.175 s |
| Stress, 256 cases and 64 points, with benchmark | 4,655,676 checks passed | 1.44e-14 | 0.949 s |
| ASan/UBSan, default cases and points | 612,534 checks passed; no sanitizer errors | 7.47e-15 | 2.210 s |
| CMake parent-project `add_subdirectory`, 2 cases and 3 points | 5,518 checks passed | 1.09e-15 | 0.021 s |

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
- Zero steady-state allocations in the prepared quadrature paths under the dedicated counting build.

The double/long-double evaluation comparison uses a quadruple-precision reference (`cpp_bin_float_quad`) for the retained polynomial. It measures the arithmetic kernel, not an independent reconstruction of the supplied nodal data.

| Element | Double maximum error | Long-double maximum error | Double time | Long-double time |
| --- | ---: | ---: | ---: | ---: |
| LINE3 | 1.89e-16 | 5.33e-20 | 0.00055 s | 0.00133 s |
| QUAD9 | 1.11e-16 | 5.8e-20 | 0.000425 s | 0.00215 s |
| TRI7 | 5.11e-16 | 2.26e-19 | 0.000416 s | 0.00196 s |
| HEX27 | 2.65e-16 | 1.07e-19 | 0.000798 s | 0.00503 s |
| TET15 | 1.55e-15 | 4.8e-19 | 0.00073 s | 0.0049 s |
| WEDGE21 | 7.52e-16 | 4.64e-19 | 0.00074 s | 0.00489 s |

Each timing covers 128,000 evaluations. Routine Horner kernels use long-double multiply/add; the sensitive root discriminant retains compensated fused arithmetic. These results do not justify changing production precision or predict production integration cost.

The new main also passes `-Wall -Wextra -Wpedantic -Werror` syntax checking with existing FEMuS/Boost headers treated as system headers. The runner passes Bash syntax checking and its standalone launch; CLI help and four invalid-argument cases were checked. `git diff --check` passes.

Source SHA-256:

- `include/CutFemLevelSetData.hpp`: `aa7b58626ca1f92a90e6ac040945aadb61926a2d376c3e6a3f40918a35521236`
- `ex47.cpp`: `6573a6b26f36b0a0cff4f5bbe621026ea75e1b9bbb3412e1e28f7577d8358bf4`

See [README](README.md) for reproducible build/run commands, [PERFORMANCE](PERFORMANCE.md) for the five-run kernel comparison, and [TODO](TODO.md) for the remaining production integration work.
