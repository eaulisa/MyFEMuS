# Prepared quadrature kernel timings

Baseline: the exact class from the user formatting commit `26cbb141a881550597e0106ba9d0691eddbfad04`. Both executables use the same added benchmark driver and deterministic inputs; the baseline has no functional class changes. GCC 14.2.0, `-O3 -DNDEBUG`, Boost 1.86.0; double has 53 significand bits and long double has 64 on this runtime. No fast-math flags or precision downgrades are used.

Each number is the median of five interleaved baseline/optimized process runs, with 20,000 calls per kernel after 64 warm-up calls. Final runs had no concurrent builds. Results describe these local kernels, not complete production integration or assembly. Checksums agree at the printed precision.

## Existing interfaces

`restrict` includes coefficient contraction, all original gradient components, maps, Bernstein controls, and classification. `nested_restrict` includes both 3D-to-2D and 2D-to-LINE restrictions. `restrict_graph` additionally certifies a graph direction on the new slice. `components_degree4` restricts three arbitrary polynomial components of coordinate degree 4. These representative restriction timings use the last reference axis; correctness tests exercise every axis.

| Element | Kernel | Before (ns) | After (ns) | Speedup |
| --- | --- | ---: | ---: | ---: |
| LINE3 | `value` | 37.3 | 29.4 | 1.27x |
| LINE3 | `root_gradient` | 44.1 | 45.2 | 0.98x |
| LINE3 | `bernstein` | 91.5 | 77.1 | 1.19x |
| LINE3 | `line_roots` | 895.2 | 868.1 | 1.03x |
| QUAD9 | `value` | 82.4 | 39.2 | 2.10x |
| QUAD9 | `root_gradient` | 172.7 | 85.9 | 2.01x |
| QUAD9 | `bernstein` | 185.3 | 174.8 | 1.06x |
| QUAD9 | `restrict` | 1164.4 | 556.5 | 2.09x |
| QUAD9 | `components_degree4` | 1000.7 | 621.2 | 1.61x |
| QUAD9 | `restrict_graph` | 1473.7 | 619.4 | 2.38x |
| TRI7 | `value` | 60.2 | 40.4 | 1.49x |
| TRI7 | `root_gradient` | 136.8 | 85.7 | 1.60x |
| TRI7 | `bernstein` | 253.8 | 159.5 | 1.59x |
| TRI7 | `restrict` | 1102.8 | 675.3 | 1.63x |
| TRI7 | `components_degree4` | 1642.8 | 1087.7 | 1.51x |
| TRI7 | `restrict_graph` | 1467.4 | 713.5 | 2.06x |
| HEX27 | `value` | 110.1 | 63.3 | 1.74x |
| HEX27 | `root_gradient` | 294.1 | 179.1 | 1.64x |
| HEX27 | `bernstein` | 465.9 | 410.1 | 1.14x |
| HEX27 | `restrict` | 3050.1 | 1388.2 | 2.20x |
| HEX27 | `nested_restrict` | 4218.0 | 1956.2 | 2.16x |
| HEX27 | `components_degree4` | 4991.4 | 3081.2 | 1.62x |
| HEX27 | `restrict_graph` | 6601.8 | 1592.3 | 4.15x |
| TET15 | `value` | 109.6 | 67.1 | 1.63x |
| TET15 | `root_gradient` | 281.3 | 183.3 | 1.54x |
| TET15 | `bernstein` | 846.1 | 476.7 | 1.77x |
| TET15 | `restrict` | 4884.9 | 1831.5 | 2.67x |
| TET15 | `nested_restrict` | 6476.0 | 2919.9 | 2.22x |
| TET15 | `components_degree4` | 7810.2 | 3965.3 | 1.97x |
| TET15 | `restrict_graph` | 7465.1 | 2237.6 | 3.34x |
| WEDGE21 | `value` | 113.8 | 68.9 | 1.65x |
| WEDGE21 | `root_gradient` | 320.2 | 189.0 | 1.69x |
| WEDGE21 | `bernstein` | 690.3 | 454.5 | 1.52x |
| WEDGE21 | `restrict` | 2319.2 | 1337.4 | 1.73x |
| WEDGE21 | `nested_restrict` | 3648.4 | 2136.6 | 1.71x |
| WEDGE21 | `components_degree4` | 4848.0 | 3235.2 | 1.50x |
| WEDGE21 | `restrict_graph` | 4216.9 | 1563.3 | 2.70x |

The vector-returning LINE gradient wrapper can be a few nanoseconds slower in this sample. The fixed-output gradient interface below avoids that allocation and is the intended quadrature interface.

## Fixed-output interfaces

| Element | Kernel | Time (ns) |
| --- | --- | ---: |
| LINE3 | `root_gradient_array` | 41.5 |
| LINE3 | `line_value_scalar` | 20.8 |
| LINE3 | `line_roots_array` | 780.3 |
| LINE3 | `line_surface_all` | 97.9 |
| QUAD9 | `root_gradient_array` | 74.3 |
| QUAD9 | `line_value_scalar` | 14.1 |
| QUAD9 | `line_roots_array` | 757.6 |
| QUAD9 | `line_surface_all` | 181.9 |
| TRI7 | `root_gradient_array` | 86.7 |
| TRI7 | `line_value_scalar` | 15.7 |
| TRI7 | `line_roots_array` | 782.1 |
| TRI7 | `line_surface_all` | 180.4 |
| HEX27 | `root_gradient_array` | 170.6 |
| HEX27 | `line_value_scalar` | 15.7 |
| HEX27 | `line_roots_array` | 769.1 |
| HEX27 | `line_surface_all` | 247.7 |
| TET15 | `root_gradient_array` | 167.2 |
| TET15 | `line_value_scalar` | 15.1 |
| TET15 | `line_roots_array` | 753.6 |
| TET15 | `line_surface_all` | 284.7 |
| WEDGE21 | `root_gradient_array` | 173.7 |
| WEDGE21 | `line_value_scalar` | 15.1 |
| WEDGE21 | `line_roots_array` | 742.7 |
| WEDGE21 | `line_surface_all` | 252.5 |

`line_surface_all` computes all four scalar densities and all original normal components together, starting from a prepared line record and supplied cofactor-gradient values. It excludes cofactor-polynomial preparation/evaluation, moment accumulation, outer quadrature, and final moment-to-weight fitting. LINE records without an intersection use a regular chart point for timing; correctness tests evaluate actual regular roots. The matrix-input surface interface is separately checked against the prepared-polynomial interface.

## Allocation checks

A separate Release build with `EX47_TRACK_ALLOCATIONS=ON` instruments global C++ allocation in the main and asserts zero allocations across the 20,000 warmed calls for restriction, nested restriction, restriction plus graph certification, component restriction, scalar evaluation, Bernstein evaluation, and every fixed-output kernel. It passes for all six element families. Timings above come from the build without allocation instrumentation. Vector compatibility wrappers retain their expected output allocation.

## Changes responsible for reuse

- Immutable reference indices, simplex multipliers, and Bernstein contraction plans are constructed once and shared.
- Strided Horner contraction avoids decoding/rebuilding polynomial exponents for every coefficient.
- Slice maps are composed sparsely in retained storage; certification subdivision uses fixed-size stack buffers.
- Bounds and graph directions are cached. Derivative bounds use Bernstein differences without building or converting derivative polynomials.
- Point validation is shared across gradient components; native quadratic evaluation uses fixed-degree Horner kernels.
- Line roots and gradients have fixed-output interfaces; the compensated long-double discriminant is retained.
- All surface density families share the line derivative and metric evaluation. Affine elements can supply one constant cofactor matrix; non-affine elements can supply restricted cofactor-gradient polynomials.

See [README](README.md) for interfaces and commands, and [VALIDATION](VALIDATION.md) for accuracy and sanitizer results.
