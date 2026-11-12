# Certified single-root timing, 2026-11-09

Baseline: `188a0034139361a1772fc0d8a5c0151c2cba8469`, whose full GitHub validation completed successfully. The candidate changes only the shared LINE root evaluation and certificate forwarding.

For each certified regular graph interval, the derivative sign is determined once from the already available midpoint root and forwarded through GK7/15, the 16/32-point graph trial, and adaptive children. At each fast quadrature sample only the corresponding stable quadratic branch is evaluated. Its discriminant square root is reused for the interface derivative magnitude. A known boundary root is reused directly. Exactly linear fibers keep their exact linear solve. Weak/mapped graph endpoints, uncertified samples, and full multiple-root fallbacks keep the previous general solver. The level-set polynomial, requested moments, tolerance tests, and quadrature rules are preserved. Initial topology and interval classification still use the existing general root/event logic.

The complete integration benchmark rebuilds both versions with identical C++14 -O3 -DNDEBUG -march=x86-64 -mtune=generic flags. Ten graph families per element, degrees 2/5/10, positive-volume and reference-surface moments, absolute and relative tolerances 1e-10, and one pinned CPU are used. Each case uses calibrated batches and the median of five measurements. Three complete runs alternate version order. Times below average the ten cases and three runs; speedup is baseline divided by candidate.

| Degree | Element | Previous ms/cell | Single root ms/cell | Speedup | Speedup range across runs |
|---:|:---|---:|---:|---:|:---|
| 2 | QUAD | 0.143 | 0.129 | 1.103x | 1.076x–1.135x |
| 2 | TRI | 0.115 | 0.104 | 1.111x | 1.079x–1.130x |
| 2 | HEX | 4.539 | 4.497 | 1.009x | 0.992x–1.022x |
| 2 | TET | 7.314 | 7.215 | 1.014x | 1.006x–1.018x |
| 5 | QUAD | 0.202 | 0.185 | 1.092x | 1.042x–1.156x |
| 5 | TRI | 0.126 | 0.115 | 1.102x | 1.058x–1.162x |
| 5 | HEX | 9.832 | 9.249 | 1.063x | 1.047x–1.075x |
| 5 | TET | 9.447 | 9.247 | 1.022x | 1.004x–1.042x |
| 10 | QUAD | 0.461 | 0.425 | 1.086x | 1.076x–1.098x |
| 10 | TRI | 0.159 | 0.147 | 1.082x | 1.069x–1.094x |
| 10 | HEX | 53.252 | 49.693 | 1.072x | 1.038x–1.113x |
| 10 | TET | 15.525 | 14.874 | 1.044x | 1.029x–1.072x |

QUAD/TRI take roughly 8–10% less time. HEX takes about 6–7% less time at degrees 5 and 10. TET takes about 1–4% less time. Low-degree 3D differences are small and should be treated cautiously. Every individual timing case retained its outer slice count: True. These are native moment-integration timings, not full ex44 assembly/solver timings.

A separate microbenchmark compares the old and single-root branches inside the same O3 binary, evaluating a curved quadratic fiber, reference-surface metric, and all requested fiber moments. It uses three repetitions of median-of-five calibrated measurements.

| Degree | Previous microseconds/fiber | Single root microseconds/fiber | Speedup |
|---:|---:|---:|---:|
| 2 | 4.691 | 3.983 | 1.178x |
| 5 | 5.054 | 4.256 | 1.187x |
| 10 | 5.999 | 5.323 | 1.127x |

Validation passed: 2,508 focused root/moment comparisons covering positive and negative slopes, both local quadratic branches, boundary roots, nearly linear equations, and certificate forwarding through both graph quadrature rules; 18 focused checks; 386 original native cases (36,860 assertions); eight supplemental cases (3,176 assertions); 12,345 independent 2D checks; 34,116 analytic 3D checks; and 6,080 enriched 3D checks. The maximum normalized fiber difference from the old solver was 8.794e-19; the largest enriched 3D discrepancy was 1.5632e-13. ASan/UBSan passed the root regression.

Reproduce the full comparison with `python unittests/testCutFemWeightLSNative/run_direct_moment_benchmark.py --baseline 188a0034139361a1772fc0d8a5c0151c2cba8469 --folder single_root --summary single_root_summary.json --rounds 3`. The fiber test accepts an extra argument to run its microbenchmark; compile it with the O3 flags above for comparable timings. Exact values and source hashes are in `single_root_summary.json` and `single_root_validation.json`.
