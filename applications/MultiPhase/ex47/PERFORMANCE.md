# Prepared quadrature kernel timings

The first comparison below records commit `67f00f5ed`. The latest selective-measure revision is documented in the final section.

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

## Selective surface evaluation

The selective revision starts from `67f00f5edd5e2de1ca61cba0a54d83eb7b3e75ea`. Dedicated functions compute one requested measure, and only explicitly named `GetLineAllSurfaceFactors` methods compute every measure. Normal-family methods share the line derivative and physical metric; single-component methods produce only the requested component while retaining the full metric norm.

These measurements compare explicit all-measure evaluation with its selected replacement in the same executable. A switch selects the benchmark's loop once, outside all point evaluations. Every pair uses identical curved-level-set data and its checksums agree exactly. Each number is the median of five isolated process runs, with 20,000 calls after 128 warm-up calls. No concurrent builds were running. The 64 evaluation points are regular line-chart points near an actual sphere/circle intersection; correctness tests separately cover exact roots, analytic geometry, and multiprecision references over scales 1e-200 through 1e200 and small-cut radii down to 1e-12.

`physical`, `tangential`, `normal0`, and `normals` consume prepared cofactor-gradient values and use slices without original-gradient context. `matrix_*` methods compute the cofactor-gradient from a supplied matrix and retain original-gradient context. `restrict_selected` compares the same recursive restriction with original-gradient contraction omitted; it does not change roots, bounds, maps, or classification. Root solving, geometry-polynomial preparation, outer quadrature, moment accumulation, and final weight fitting are outside the density timings.

| Element | Selected kernel | All measures/full context (ns) | Selected (ns) | Speedup |
| --- | --- | ---: | ---: | ---: |
| LINE3 | `delta` | 75.9 | 15.9 | 4.78x |
| LINE3 | `reference` | 81.1 | 42.2 | 1.92x |
| LINE3 | `physical` | 74.8 | 21.6 | 3.46x |
| LINE3 | `tangential` | 74.3 | 21.3 | 3.48x |
| LINE3 | `normal0` | 74.5 | 23.8 | 3.12x |
| LINE3 | `normals` | 75.9 | 41.1 | 1.85x |
| LINE3 | `matrix_physical` | 151.7 | 53.5 | 2.84x |
| LINE3 | `matrix_tangential` | 154.2 | 54.4 | 2.83x |
| LINE3 | `matrix_normal0` | 149.5 | 62.0 | 2.41x |
| LINE3 | `matrix_normals` | 153.9 | 89.0 | 1.73x |
| QUAD9 | `delta` | 176.2 | 17.5 | 10.04x |
| QUAD9 | `reference` | 165.9 | 88.3 | 1.88x |
| QUAD9 | `physical` | 170.6 | 62.5 | 2.73x |
| QUAD9 | `tangential` | 169.3 | 61.7 | 2.74x |
| QUAD9 | `normal0` | 169.3 | 64.4 | 2.63x |
| QUAD9 | `normals` | 165.3 | 91.9 | 1.80x |
| QUAD9 | `matrix_physical` | 242.7 | 109.6 | 2.21x |
| QUAD9 | `matrix_tangential` | 239.4 | 109.4 | 2.19x |
| QUAD9 | `matrix_normal0` | 253.6 | 130.9 | 1.94x |
| QUAD9 | `matrix_normals` | 237.6 | 143.6 | 1.65x |
| QUAD9 | `restrict_selected` | 571.2 | 464.6 | 1.23x |
| TRI7 | `delta` | 158.1 | 14.9 | 10.62x |
| TRI7 | `reference` | 157.6 | 83.8 | 1.88x |
| TRI7 | `physical` | 158.6 | 60.7 | 2.61x |
| TRI7 | `tangential` | 163.4 | 61.8 | 2.65x |
| TRI7 | `normal0` | 164.3 | 62.9 | 2.61x |
| TRI7 | `normals` | 158.9 | 90.6 | 1.75x |
| TRI7 | `matrix_physical` | 241.6 | 112.1 | 2.15x |
| TRI7 | `matrix_tangential` | 247.3 | 112.7 | 2.20x |
| TRI7 | `matrix_normal0` | 265.8 | 114.7 | 2.32x |
| TRI7 | `matrix_normals` | 243.2 | 156.4 | 1.55x |
| TRI7 | `restrict_selected` | 509.3 | 399.6 | 1.27x |
| HEX27 | `delta` | 252.7 | 15.5 | 16.29x |
| HEX27 | `reference` | 230.8 | 141.5 | 1.63x |
| HEX27 | `physical` | 247.3 | 107.8 | 2.29x |
| HEX27 | `tangential` | 243.4 | 107.5 | 2.26x |
| HEX27 | `normal0` | 249.6 | 110.2 | 2.27x |
| HEX27 | `normals` | 233.1 | 136.4 | 1.71x |
| HEX27 | `matrix_physical` | 340.4 | 174.9 | 1.95x |
| HEX27 | `matrix_tangential` | 336.4 | 177.8 | 1.89x |
| HEX27 | `matrix_normal0` | 339.2 | 194.6 | 1.74x |
| HEX27 | `matrix_normals` | 360.0 | 202.7 | 1.78x |
| HEX27 | `restrict_selected` | 1935.8 | 1598.9 | 1.21x |
| TET15 | `delta` | 244.6 | 15.5 | 15.82x |
| TET15 | `reference` | 257.8 | 144.0 | 1.79x |
| TET15 | `physical` | 242.4 | 110.5 | 2.19x |
| TET15 | `tangential` | 231.0 | 107.1 | 2.16x |
| TET15 | `normal0` | 234.8 | 112.2 | 2.09x |
| TET15 | `normals` | 250.9 | 138.6 | 1.81x |
| TET15 | `matrix_physical` | 343.1 | 186.8 | 1.84x |
| TET15 | `matrix_tangential` | 352.8 | 179.4 | 1.97x |
| TET15 | `matrix_normal0` | 339.1 | 195.9 | 1.73x |
| TET15 | `matrix_normals` | 346.2 | 213.4 | 1.62x |
| TET15 | `restrict_selected` | 2683.8 | 1888.2 | 1.42x |
| WEDGE21 | `delta` | 243.9 | 15.4 | 15.85x |
| WEDGE21 | `reference` | 241.8 | 147.3 | 1.64x |
| WEDGE21 | `physical` | 249.4 | 107.5 | 2.32x |
| WEDGE21 | `tangential` | 239.4 | 107.6 | 2.23x |
| WEDGE21 | `normal0` | 239.5 | 111.0 | 2.16x |
| WEDGE21 | `normals` | 230.7 | 134.3 | 1.72x |
| WEDGE21 | `matrix_physical` | 349.4 | 184.2 | 1.90x |
| WEDGE21 | `matrix_tangential` | 348.1 | 182.1 | 1.91x |
| WEDGE21 | `matrix_normal0` | 348.2 | 200.1 | 1.74x |
| WEDGE21 | `matrix_normals` | 337.5 | 201.8 | 1.67x |
| WEDGE21 | `restrict_selected` | 1456.9 | 941.5 | 1.55x |

The allocation-counting build passes both `--benchmark-surface` and `--benchmark-quadrature`, asserting zero C++ heap allocations in every warmed selected density and restriction kernel. The existing all-measure helper remains available only as an explicit batch operation; no precision or roundoff protection was reduced.

## Prepared topology, charts, and projection revision

Baseline: `2b2dc6cc95ba8787d7c1ffa43c3d5fc460a55ce6`, with the exact preceding class and the same common ex47 benchmark functions. Both use GCC 14.2.0, Boost 1.86.0, and `-O3 -DNDEBUG`. Numbers are medians of five interleaved process runs with no concurrent builds. Comparable checksums agree at the recorded precision; graph-direction checksums are excluded from equality checks because tighter bounds can legitimately select another direction. These are prepared-data kernels, not full production assembly timings.

### Fresh slicing and roots

Each call creates a new LINE state in retained storage, solves its roots, and returns the root count and first root. The parent and any first 3D-to-2D slice are prepared outside the loop. Original-gradient context is omitted in both versions.

| Element | Before (ns) | After (ns) | Speedup |
| --- | ---: | ---: | ---: |
| QUAD9 | 1318.9 | 702.1 | 1.88x |
| TRI7 | 1190.6 | 651.3 | 1.83x |
| HEX27 | 1375.7 | 798.5 | 1.72x |
| TET15 | 1421.7 | 646.8 | 2.20x |
| WEDGE21 | 1229.9 | 559.9 | 2.20x |

The same workload with the selected physical density, using supplied cofactor-gradient values:

| Element | Before (ns) | After (ns) | Speedup |
| --- | ---: | ---: | ---: |
| QUAD9 | 1461.2 | 856.7 | 1.71x |
| TRI7 | 1298.9 | 625.3 | 2.08x |
| HEX27 | 1518.2 | 957.1 | 1.59x |
| TET15 | 1439.7 | 824.1 | 1.75x |
| WEDGE21 | 1288.4 | 721.5 | 1.79x |

### Reusing roots in the same slice

The old fixed-output interface solves the same quadratic on every call. The new interface returns the retained topology. This reuse is useful when volume and selected surface families consume the same slice; it does not represent the cost of solving a new slice.

| Element | Before (ns) | After (ns) | Speedup |
| --- | ---: | ---: | ---: |
| LINE3 | 649.9 | 23.3 | 27.87x |
| QUAD9 | 846.2 | 22.9 | 36.88x |
| TRI7 | 798.1 | 21.7 | 36.82x |
| HEX27 | 827.8 | 26.4 | 31.35x |
| TET15 | 850.3 | 21.4 | 39.71x |
| WEDGE21 | 808.3 | 20.9 | 38.62x |

### Existing restriction and certification costs

These calls retain every original gradient component, as in the preceding benchmark. Additional state ownership and cache invalidation have a measurable cost. Default direction selection uses strict bounds only; weak-graph/zero-fiber preparation is performed when explicitly requested. QUAD corner signs reject impossible strict directions before tighter extrema are considered.

| Element | Kernel | Before (ns) | After (ns) | Before/after |
| --- | --- | ---: | ---: | ---: |
| QUAD9 | `restrict` | 600.1 | 641.0 | 0.94x |
| QUAD9 | `restrict_graph` | 639.6 | 686.8 | 0.93x |
| TRI7 | `restrict` | 637.3 | 652.8 | 0.98x |
| TRI7 | `restrict_graph` | 727.2 | 718.8 | 1.01x |
| HEX27 | `restrict` | 1279.6 | 1399.5 | 0.91x |
| HEX27 | `nested_restrict` | 1934.3 | 2412.3 | 0.80x |
| HEX27 | `restrict_graph` | 1626.0 | 1692.9 | 0.96x |
| TET15 | `restrict` | 1681.8 | 1716.2 | 0.98x |
| TET15 | `nested_restrict` | 2562.3 | 2571.6 | 1.00x |
| TET15 | `restrict_graph` | 2084.3 | 2212.2 | 0.94x |
| WEDGE21 | `restrict` | 1215.0 | 1183.4 | 1.03x |
| WEDGE21 | `nested_restrict` | 1930.0 | 1964.9 | 0.98x |
| WEDGE21 | `restrict_graph` | 1476.4 | 1710.6 | 0.86x |

The restriction/certification benchmark therefore should not be described as a universal speedup. The measured slice-plus-root workload benefits from retained analysis and a cheap well-conditioned discriminant path. Sensitive contractions keep compensation and ancestor recovery rather than reducing arithmetic precision.

### New preparation methods

Chart construction and external component composition are setup operations. Prepare a selected chart outside the quadrature-point loop. These timings use reused output storage; root-record delta evaluation consumes an already prepared regular root.

| Element | Axis chart (ns) | Barycentric chart (ns) | External degree 3 (ns) | Root-record delta (ns) |
| --- | ---: | ---: | ---: | ---: |
| QUAD9 | 1021.0 | — | 1378.8 | 19.5 |
| TRI7 | 451.6 | 2243.2 | 3800.3 | 18.3 |
| HEX27 | 3155.8 | — | 7083.5 | 18.8 |
| TET15 | 1630.6 | 11121.5 | 23399.0 | 18.5 |
| WEDGE21 | 4906.5 | 18771.8 | 18615.2 | 18.5 |

The record itself is 1200 bytes versus 864 bytes before this revision on this ABI; these figures exclude owned coefficient buffers and optional cache allocations. LINE/graph/projection buffers are created lazily. An immutable ancestor snapshot is retained once when restriction needs it; descendants share it. Copy-on-write cache preparation preserves copied records. The class owns no geometry, GK objects, moment matrices, or final fitted weights.

Allocation tracking passes 27 data-reuse, 45 prepared-algebra, 61 quadrature, and 130 selected-surface benchmark rows. Warmed fixed-output paths allocate zero C++ heap objects; vector-returning compatibility wrappers allocate their outputs as before. First preparation and detachment after copying are not included in this steady-state claim. Sloping-boundary convenience restriction constructs a temporary chart; repeated traversal should retain that chart explicitly.

Accuracy includes 4,939,947 stress checks, 100-digit retained-parent references, last-ULP folds in every 3D family, hidden tiny cuts, all supported charts/boundaries, and every selected surface mode. Delta error in the translated near-double-root test decreases from 1.318e-11 to 1.085e-19; maximum last-ULP fold density error is 1.882e-17. See [VALIDATION](VALIDATION.md) for exact scope and sanitizer results.
