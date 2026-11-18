# Prepared quadrature kernel timings

Historical record from before promotion into the native library; see [MIGRATION.md](MIGRATION.md) for the current layout and validation.

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


## Local weight-class comparison

This compares the production `::CutFemWeight_LS` with the local `Ex47::CutFemWeight_LS`. GCC Release flags are `-O2 -DNDEBUG`, and fitting arithmetic is long double. Both objects retain their Gauss points and unchanged multiprecision-built fitting operator outside the timed loops. Every timed call clears the final-weight cache. Two different prepared element records alternate; the geometry/metric must refresh when changing fields. The five interleaved runs report medians. No compilation or sanitizer process runs concurrently with these measurements.

The 2D/3D fixtures contain closed reference circles/spheres strictly inside each element and an explicitly non-affine physical map. These timings describe that demanding closed-surface fixture. LINE uses two distinct quadratic nodal fields. All benchmark checksums agree with production within the stated comparison tolerances.

Entries below are production time divided by prepared time, over the tested degrees: values above 1 mean faster; values below 1 mean slower.

| Element | Delta | Reference surface | Physical surface | Tangential | Three normals |
| --- | ---: | ---: | ---: | ---: | ---: |
| LINE3 | 5.79–11.50× | 5.56–9.14× | 5.29–7.30× | 4.92–9.19× | — |
| QUAD9 | 0.89–0.95× | 1.28–1.29× | 1.91–2.05× | 1.91–2.04× | — |
| TRI7 | 1.02–1.06× | 1.42–1.42× | 1.81–1.83× | 1.92–1.95× | — |
| HEX27 | 0.84–0.94× | 1.24–1.30× | 2.98–4.55× | 3.09–3.24× | 3.61–3.92× |
| TET15 | 0.96–0.96× | 1.26–1.84× | 4.09–4.58× | 4.36–4.73× | 5.45–7.19× |
| WEDGE21 | 1.07–1.25× | 1.86–1.87× | 4.11–4.15× | 3.69–4.07× | 4.49–4.83× |

Physical/tangential and normal-family requests improve substantially. Raw delta remains 4–19% slower in some QUAD/HEX/TET cases; TRI/WEDGE are neutral or faster. The native-coordinate path and retained sensitive-root information preserve small cuts, but they do not produce a universal runtime improvement. This prototype is therefore not an automatic production switch. A shared traversal for all three normal components remains pending; the measured normal-family method performs three selected integrations and reuses its metric/positive weights.

The local integration prepares all full-strip Gauss rules and shared GK rule objects before weight integration. It reuses prepared data, derivative/projection caches, LINE topology, and native physical metric polynomials. Sensitive ancestor fibers are contracted directly in compensated arithmetic, without approximation or three repeated whole-polynomial evaluations. Full adaptive integration still uses vector workspaces and makes allocations; the zero-allocation claims in earlier kernel sections do not apply to the complete weight integration.

Absolute median times in microseconds are given below. LINE includes degrees 2, 5, and 10; all other elements include degrees 2 and 5.

| Element | Degree | Measure | Production µs | Prepared µs | Speed ratio |
| --- | ---: | --- | ---: | ---: | ---: |
| LINE3 | 2 | delta | 5.456 | 0.474 | 11.505× |
| LINE3 | 2 | reference | 5.697 | 0.623 | 9.140× |
| LINE3 | 2 | physical | 6.285 | 0.875 | 7.184× |
| LINE3 | 2 | tangential | 5.629 | 0.612 | 9.193× |
| LINE3 | 5 | delta | 5.871 | 0.613 | 9.580× |
| LINE3 | 5 | reference | 6.185 | 0.852 | 7.259× |
| LINE3 | 5 | physical | 6.092 | 0.834 | 7.303× |
| LINE3 | 5 | tangential | 6.158 | 0.860 | 7.157× |
| LINE3 | 10 | delta | 7.495 | 1.295 | 5.787× |
| LINE3 | 10 | reference | 7.389 | 1.330 | 5.556× |
| LINE3 | 10 | physical | 7.491 | 1.416 | 5.291× |
| LINE3 | 10 | tangential | 7.484 | 1.520 | 4.923× |
| QUAD9 | 2 | delta | 440.421 | 495.564 | 0.889× |
| QUAD9 | 2 | reference | 724.346 | 565.275 | 1.281× |
| QUAD9 | 2 | physical | 1125.108 | 549.992 | 2.046× |
| QUAD9 | 2 | tangential | 1137.749 | 557.879 | 2.039× |
| QUAD9 | 5 | delta | 516.450 | 542.836 | 0.951× |
| QUAD9 | 5 | reference | 751.381 | 581.444 | 1.292× |
| QUAD9 | 5 | physical | 1139.707 | 596.066 | 1.912× |
| QUAD9 | 5 | tangential | 1175.590 | 614.758 | 1.912× |
| TRI7 | 2 | delta | 1296.684 | 1266.932 | 1.023× |
| TRI7 | 2 | reference | 1982.807 | 1400.326 | 1.416× |
| TRI7 | 2 | physical | 5887.096 | 3216.398 | 1.830× |
| TRI7 | 2 | tangential | 3090.777 | 1612.943 | 1.916× |
| TRI7 | 5 | delta | 1439.914 | 1364.717 | 1.055× |
| TRI7 | 5 | reference | 2112.825 | 1487.732 | 1.420× |
| TRI7 | 5 | physical | 2967.841 | 1635.936 | 1.814× |
| TRI7 | 5 | tangential | 6077.175 | 3119.516 | 1.948× |
| HEX27 | 2 | delta | 25152.567 | 29949.017 | 0.840× |
| HEX27 | 2 | reference | 38305.100 | 29442.036 | 1.301× |
| HEX27 | 2 | physical | 137801.930 | 30316.367 | 4.545× |
| HEX27 | 2 | tangential | 98213.886 | 30352.890 | 3.236× |
| HEX27 | 2 | normal family | 346939.680 | 88607.343 | 3.915× |
| HEX27 | 5 | delta | 27206.716 | 28801.361 | 0.945× |
| HEX27 | 5 | reference | 41050.532 | 33102.370 | 1.240× |
| HEX27 | 5 | physical | 100228.210 | 33668.516 | 2.977× |
| HEX27 | 5 | tangential | 102229.120 | 33098.897 | 3.089× |
| HEX27 | 5 | normal family | 306261.250 | 84874.187 | 3.608× |
| TET15 | 2 | delta | 72416.759 | 75345.529 | 0.961× |
| TET15 | 2 | reference | 146141.110 | 79315.215 | 1.843× |
| TET15 | 2 | physical | 828222.470 | 180694.390 | 4.584× |
| TET15 | 2 | tangential | 405096.780 | 85716.464 | 4.726× |
| TET15 | 2 | normal family | 1698011.100 | 236297.370 | 7.186× |
| TET15 | 5 | delta | 79188.933 | 82534.198 | 0.959× |
| TET15 | 5 | reference | 163616.360 | 130088.920 | 1.258× |
| TET15 | 5 | physical | 428402.470 | 104624.560 | 4.095× |
| TET15 | 5 | tangential | 426234.700 | 97689.477 | 4.363× |
| TET15 | 5 | normal family | 1367812.800 | 250971.640 | 5.450× |
| WEDGE21 | 2 | delta | 86240.694 | 69027.144 | 1.249× |
| WEDGE21 | 2 | reference | 140906.200 | 75503.012 | 1.866× |
| WEDGE21 | 2 | physical | 323080.650 | 78603.611 | 4.110× |
| WEDGE21 | 2 | tangential | 570376.210 | 139992.770 | 4.074× |
| WEDGE21 | 2 | normal family | 1609862.000 | 358150.360 | 4.495× |
| WEDGE21 | 5 | delta | 73677.757 | 68781.466 | 1.071× |
| WEDGE21 | 5 | reference | 136847.390 | 73721.110 | 1.856× |
| WEDGE21 | 5 | physical | 589133.510 | 141839.650 | 4.154× |
| WEDGE21 | 5 | tangential | 367380.700 | 99666.227 | 3.686× |
| WEDGE21 | 5 | normal family | 1050922.700 | 217566.820 | 4.830× |
