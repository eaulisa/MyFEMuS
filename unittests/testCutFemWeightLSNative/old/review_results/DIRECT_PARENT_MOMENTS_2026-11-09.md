# Direct parent-coordinate CutFEM moments, 2026-11-09

The regular default paths now accumulate monomials in the original FEMuS parent reference coordinates inside the shared LINE evaluator. Local coefficient, root-finding, and quadrature coordinate maps remain available. Their scale and measure are applied while sampling; completed moment tables are no longer reconstructed for regular QUAD/HEX integration or normalized TET slices. The default TRI path uses its standard parent triangle directly.

The requested 3D moment mask is projected onto the required 2D monomials and passed through QUAD/TRI into LINE. The zeroth volume and surface moments are retained when needed for cut classification. The certified 3D graph direction is passed into the slice solver, allowing it to skip a repeated derivative certification and exclude unnecessary fiber-discriminant event searches. Actual slice boundary roots still depend on slice height and are evaluated there.

The 16/32-point graph trial is retained for mapped endpoints and degrees above five; certified regular lower-degree graphs use adaptive GK7/15 directly. Taylor and full-GK fallbacks remain. Explicit heavy validation, forced diagnostic paths, and established exceptional boundary traces retain their legacy routes, including their existing transformations. No application API changes or physical Jacobian/metric changes are required.

Comparison uses the previous result commit `bd1ddd478433e40da0e09f57a1f11f5d2ce0f07b`, the same ten quadratic graph families per element, degrees 2/5/10, positive-volume and reference-surface moments, and absolute/relative tolerances 1e-10. Both versions are rebuilt with identical C++14 -O3 -DNDEBUG -march=x86-64 -mtune=generic flags and pinned to one CPU. Each case uses calibrated batches and the median of five measurements; two complete runs reverse the version order. Times below are averages over the ten cases and two runs. Speedup is old time divided by new time.

| Degree | Element | Previous ms/cell | Direct ms/cell | Speedup | Mean outer slices |
|---:|:---|---:|---:|---:|---:|
| 2 | QUAD | 0.137 | 0.141 | 0.97x | 0.0 |
| 2 | TRI | 0.118 | 0.114 | 1.04x | 0.0 |
| 2 | HEX | 4.549 | 4.616 | 0.99x | 13.6 |
| 2 | TET | 8.236 | 7.377 | 1.12x | 40.5 |
| 5 | QUAD | 0.199 | 0.214 | 0.93x | 0.0 |
| 5 | TRI | 0.160 | 0.127 | 1.26x | 0.0 |
| 5 | HEX | 9.472 | 10.137 | 0.93x | 19.6 |
| 5 | TET | 10.987 | 9.907 | 1.11x | 46.5 |
| 10 | QUAD | 0.507 | 0.471 | 1.08x | 0.0 |
| 10 | TRI | 0.227 | 0.162 | 1.40x | 0.0 |
| 10 | HEX | 75.326 | 53.326 | 1.41x | 69.1 |
| 10 | TET | 20.931 | 15.792 | 1.33x | 58.5 |

Degree-10 HEX time falls from 75.326 to 53.326 ms/cell (29% less time, 1.41x speedup), and TET falls from 20.931 to 15.792 ms/cell (25% less time, 1.33x speedup). Degree-2 HEX is essentially unchanged. Degree-5 QUAD and HEX take about 7% more time. TRI improves most at degrees 5 and 10. These results support a measured tradeoff rather than a universal speedup.

Outer slice counts are unchanged in this matrix. Eliminating the 4*3^q moment-reconstruction budget therefore removes part of the cost but does not explain the large degree-dependent increase by itself; substantial outer adaptive work remains. These timings cover native moment integration rather than full ex44 assembly or solves.

Validation completed without numerical failures: 17 focused checks, including 408 signed/reflected/sparse/Taylor parent-coordinate checks; 386 original native cases with 36,860 assertions; eight supplemental cases with 3,176 assertions; 147 planar comparisons; 12,345 independent 2D regression checks; 34,116 analytic 3D checks in 132 cases; and 6,080 enriched 3D comparisons in 24 cases. The complete direct-coordinate mesh/tolerance matrix passed; the final degree-dependent quadrature policy was rechecked on the 2D oracle, all three analytic 3D degrees, a coarse degree-10 TET sphere mesh, and the complete enriched 3D matrix. ASan/UBSan checks passed for the direct LINE and 3D analytic paths. The final enriched maximum absolute discrepancy is 1.5632e-13. The full GitHub workflow repeats its complete matrix on the pushed source.

Reproduce the timing comparison with `python unittests/testCutFemWeightLSNative/run_direct_moment_benchmark.py` after `build.py`; use `run_focused.py --extended`, `run_hex_tet_validation.py`, and `validate_enriched_3d.py` for validation. Curated exact timing and source hashes are in `direct_moment_summary.json`; acceptance details are in `direct_moment_validation.json`.
