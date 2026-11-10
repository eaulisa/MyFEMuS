# 2D tangency-strip experiment — 2026-11-10

Branch: `rifat-cutfem-2d-tangent-strips`. Baseline: `1a394fdc97bedc5a6cff9b1e6c68f16261ba57e9` on `rifat-cutfem-hex-tet-result`. The `rifat` benchmark remains unchanged.

The experiment passes the numerical checks, including fitted weights, but does not provide a general speed improvement. Two-sheet cuts benefit greatly; circles and the simplex bubble generally regress. Keep this branch as an experiment; do not replace the baseline with it.

## Implementation

Existing certified and weak graph integration bypasses the new method. No LINE, QUAD graph, TRI graph, HEX or TET integration kernel was changed. The dispatch additions are in `UnitElementMoments.hpp` and `Tri7DirectionSelection.hpp`; the new implementation is `TangencyStripIntegration.hpp`.

For a non-graph cell, reconstruct the same nodal level set as before. In each candidate direction write `G(x,y)=A(x)y²+B(x)y+C(x)`. Locate regular tangencies from the discriminant `B²−4AC=0`, check the repeated root lies inside the square or triangle, and include all boundary crossings and zeros of A in the ordered outer events. Classify root order, clipping and positive-phase segments once per interval. Choose the direction using the number of active intervals, tangencies and representative slopes; try the alternate direction if needed.

Near a simple fold `(xt,yt)`, use the explicit parabola `Q=Gx(xt,yt)(x−xt)+Gyy(xt,yt)(y−yt)²/2`. Start with delta equal to 5% of its adjacent outer interval. Set `x=xt±delta*u²`, so the two parabolic roots are linear in u. Integrate the polynomial volume and raw-delta terms by analytic coefficient integration. Surface norms use a convergent analytic series with a truncation estimate; they are not evaluated by GK inside the strips, and this is not a general elementary closed-form surface primitive. The model uses its own consistent gradient/cofactor surface metric.

Estimate the model error against the regularized analytic expansion of the original G, including truncation and endpoint-location roundoff. Halve delta only for moments that have not passed; retain completed components. The regular complement uses selective 7/15 GK, then 15/31 GK, and subdivision only for failing components. Quadrature nodes reuse the classified topology and share one discriminant evaluation and powers across moments.

A unit-moment-only prototype failed two fitted QUAD circle checks. Independent per-moment approximation biases combined in the signed FEMuS basis. The final QUAD endpoint approximation reserves a factor `4*3^q`: the coefficient 1-norm of `(2x−1)^i(2y−1)^j` is `3^(i+j)` and the parent-cell area factor is 4. This reserve applies only to the new endpoint model acceptance; it does not change graph integration, bulk GK tolerances, or the stored unit monomial basis. It is conservative and contributes strongly to the poor QUAD result.

If a direction is degenerate, the strip expansion fails its estimates, or the budget is exhausted, discard its local result and retain the existing integrator. Prescribed-direction requests, the strongest direction-validation mode and degrees above 10 retain the existing paths.

## Numerical validation

- 405 independent-reference cases at q=2,5,10; 25,110 checked unit moments. Raw delta, reference surface and curved physical geometry are covered. The independent oracle uses 80-digit root calculations and separate 128/256-point integration.
- 180 graph cases have serialized moments identical to baseline and never enter the new method.
- 173 wider 2D cases; 31,068 checks, including level-set scaling, small cuts, forced directions, fitted weights and boundary cases.
- 36 physical-mesh configurations, through 4,096 QUAD or 8,192 TRI cells and q=10.
- 8 degree-10 double-output fitted-weight cases; 394 original/supplemental native cases (40,036 assertions); 19 focused executions; 147 planar comparisons.
- 135 cases pass AddressSanitizer and UndefinedBehaviorSanitizer. Local leak detection is disabled because the runtime is traced; CI enables it.
- 16 unchanged HEX/TET cases at q=2,10; 7,104 checks and serialized moment equality.
- All 30 standalone translation units build. Reports bind to SHA-256 hashes of the implementation headers. Full MPI/PETSc ex44 was not built or run in this environment.

Of 207 eligible calls, 133 complete using the new method and 74 fall back. These counts cover the three modes and three degrees; the 405 total cases must not be presented as 405 successful uses of the new path.

## Absolute timings

Times below are milliseconds per full-cell moment integration at absolute/relative tolerance 1e-10. Each includes nodal reconstruction, both moment families and direction/topology work, but excludes fitted-weight setup and output formatting. Identical C++14 O3 generic-x86 flags and fixed CPU affinity are used. Five batches of at least 8 ms are measured; neighboring level sets alternate so the geometry cache cannot hide reconstruction costs. Medians are compared against the pinned baseline. Hardware: AMD EPYC 9V74, CPU 0.

Surface-mode examples at q=10:

| Element | Curve | Baseline ms | Experiment ms | Baseline/experiment | Result |
|---|---|---:|---:|---:|---|
| QUAD | circle | 4.0195 | 37.4575 | 0.11x | legacy fallback |
| TRI | circle | 2.8333 | 28.6694 | 0.10x | new path |
| QUAD | two_sheets | 0.6026 | 0.1079 | 5.59x | new path |
| TRI | two_sheets | 6.0470 | 0.4218 | 14.34x | new path |
| TRI | rotated_ellipse | 14.9579 | 16.5848 | 0.90x | new path |
| TRI | bubble | 5.6150 | 30.3925 | 0.18x | new path |
| QUAD | tensor_lobes | 21.1136 | 63.4968 | 0.33x | legacy fallback |

Mean of the 13 QUAD / 12 TRI non-graph-labeled surface-mode examples (includes fallbacks):

| Element | q | Baseline ms | Experiment ms | Experiment/baseline |
|---|---:|---:|---:|---:|
| QUAD | 2 | 4.6535 | 13.3232 | 2.86x |
| QUAD | 5 | 4.9543 | 12.5680 | 2.54x |
| QUAD | 10 | 5.9009 | 28.9320 | 4.90x |
| TRI | 2 | 4.0020 | 10.6559 | 2.66x |
| TRI | 5 | 4.3638 | 9.2935 | 2.13x |
| TRI | 10 | 5.5202 | 12.8795 | 2.33x |

The QUAD q=10 circle spends 64 model refinements / 128 analytic panels across the two attempted directions, then falls back. The TRI q=10 circle completes with 40 model refinements / 84 analytic panels and 922 bulk samples; its smallest strip delta is about 1.4e-8. Thus locating tangencies and eliminating geometric cell splits does not remove the main cost: certifying a sufficiently accurate parabola makes the strips very thin, and the remaining GK interval still approaches a steep endpoint. The two-sheet gains come from cached interval topology and analytic inner integration; they do not establish a speed benefit from tangent strips.

Separate sustained timing on the original ten graphs per element (three interleaved rounds):

| Element | q | Baseline ms | Experiment ms | Experiment/baseline |
|---|---:|---:|---:|---:|
| QUAD | 2 | 0.1206 | 0.1250 | 1.036x |
| TRI | 2 | 0.1422 | 0.1439 | 1.012x |
| QUAD | 5 | 0.1554 | 0.1293 | 0.832x |
| TRI | 5 | 0.1555 | 0.1558 | 1.002x |
| QUAD | 10 | 0.2503 | 0.2804 | 1.120x |
| TRI | 10 | 0.1872 | 0.2029 | 1.084x |

Graph outputs are unchanged. Treat timing fluctuations as measurement/compiler-layout effects, not an algorithmic graph speedup.

## Reproduction

The dedicated `.github/workflows/cutfem-2d-tangency-strips.yml` builds and runs the complete experiment. From `unittests/testCutFemWeightLSNative`, run `python3 build.py --jobs 2` with Eigen/Boost paths if necessary; then the original/supplemental runners, focused/planar checks, the 2D wide-moment/mesh/double-output runners, `run_tangency_strip_experiment.py --rounds 1 --measure`, `run_tangency_preservation.py`, `run_tangency_sanitizers.py`, and finally `run_tangency_acceptance.py`. The workflow contains the exact commands and uploads raw numerical evidence.

Committed JSON reports contain all cases, timings, fallback counts, errors and implementation hashes. The experiment report records all three modes at q=2,5,10.
