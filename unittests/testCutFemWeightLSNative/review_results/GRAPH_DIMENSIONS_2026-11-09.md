# Graph integration timing: 2D versus 3D, 2026-11-09

Production source: `73190d0b7013d574e9509be21b327447058b22ba`. Ten globally monotone graph families per element and degree; all families are quadratic polynomials represented exactly by the FE spaces, up to nodal rounding. The 3D family extends its 2D counterpart with a second independent coordinate. The named saddle family reduces to a parabola in 2D.

Direct native nodal integration computes all positive-volume and reference-surface moments of total degree at most q. Absolute and relative tolerances are both 1e-10. Nodal inputs are prepared outside timing; nodal coefficient extraction, surface gradients, adaptive integration, result allocation, conversion, and a checksum of every requested moment are timed. Gauss-weight fitting and ex44 assembly are outside this measurement.

Compiled with g++ 14.2, C++14, -O3 -DNDEBUG -march=x86-64 -mtune=generic; one process pinned to CPU 0. Each case uses a warm-up, calibrated batches of at least approximately 8 ms, and the median of five batches. The complete experiment was repeated twice. The table averages the two medians per case, then averages the ten cases. Ratios divide summed 3D times by summed 2D times. Wall-clock variation remains; rounded ratios are appropriate.

| Degree | Comparison | Mean 2D (ms/cell) | Mean 3D (ms/cell) | Ratio | Ratio, run 1 / run 2 | Mean outer slices |
|---:|:---|---:|---:|---:|:---|---:|
| 2 | HEX/QUAD | 0.166 | 5.021 | 30.3x | 25.9x / 35.5x | 13.6 |
| 2 | TET/TRI | 0.127 | 9.549 | 74.9x | 80.0x / 70.5x | 40.5 |
| 5 | HEX/QUAD | 0.206 | 9.864 | 48.0x | 48.9x / 47.0x | 19.6 |
| 5 | TET/TRI | 0.162 | 11.316 | 70.0x | 70.6x / 69.4x | 46.5 |
| 10 | HEX/QUAD | 0.540 | 80.267 | 148.6x | 156.2x / 141.0x | 69.1 |
| 10 | TET/TRI | 0.235 | 22.217 | 94.7x | 93.8x / 95.6x | 58.5 |

All 120 element/degree/case combinations completed in both passes, with positive finite zeroth volume and surface moments, finite checksums, and no 3D geometric splitting or fallback. These are performance probes using the kernels' internal tolerance checks; this experiment does not add independent reference solutions for the curved cases.

| Graph family | q2 HEX/QUAD | q2 TET/TRI | q5 HEX/QUAD | q5 TET/TRI | q10 HEX/QUAD | q10 TET/TRI |
|:---|---:|---:|---:|---:|---:|---:|
| axis_plane | 2.6x | 53.6x | 2.4x | 52.7x | 3.1x | 50.7x |
| tilted_plane | 31.7x | 45.9x | 30.1x | 40.9x | 68.6x | 35.3x |
| near_plane | 28.2x | 37.7x | 27.5x | 38.6x | 63.7x | 32.2x |
| gentle_convex | 29.2x | 41.8x | 29.5x | 41.4x | 151.8x | 37.1x |
| gentle_concave | 31.9x | 128.0x | 57.8x | 175.1x | 110.4x | 163.1x |
| convex | 19.1x | 70.9x | 64.4x | 15.9x | 156.7x | 73.6x |
| saddle | 43.1x | 47.7x | 68.3x | 43.9x | 175.6x | 41.4x |
| mixed_curvature | 33.0x | 272.0x | 43.3x | 342.1x | 197.9x | 371.1x |
| implicit_quadratic | 29.1x | 38.0x | 65.6x | 36.3x | 104.0x | 93.0x |
| variable_slope | 54.5x | 41.8x | 29.2x | 40.9x | 186.2x | 47.3x |

At degree 2 the HEX axis-plane shortcut uses one slice and all other HEX cases use 15. TET cases use 30 to 105 slices, averaging 40.5. At degree 10 the maxima are 150 for HEX and 165 for TET. Several degree-10 HEX cases revisit the same graph region to tighten the combined error budget; the graph-region diagnostic therefore counts visits, including repeats.

The 20x estimate understates the total integration cost for most tested cases. Besides outer sample count, each 3D call pays for tighter inner integration budgets, slice construction, gradient preparation, moment transforms, topology checks, and final error accumulation. The timings alone do not identify which of these is the largest optimization target.

Benchmark source: `../benchmark_graph_dimensions.cpp`. Reuse the existing standalone FE objects from `build.py`. Load `build/hex_tet_validation_command.json`, replace its source with `benchmark_graph_dimensions.cpp`, replace -O2 with -O3, add -DNDEBUG -march=x86-64 -mtune=generic, and change its output executable to `build/benchmark_graph_dimensions`. Run `taskset -c 0 build/benchmark_graph_dimensions 2`, then repeat with 5 and 10.
