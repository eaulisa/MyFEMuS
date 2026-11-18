# Unified outer quadrature

The baseline is `rifat-cutfem-endpoint-only` at `f5f81a8a088cb990780f326c4a8692ac98265536` (`all example worked so far`). The final design follows the November 10 correction: G7/K15 is the primary rule on every outer interval. `rifat` is unchanged.

| Endpoint classification | Physical map | Primary rule |
|---|---|---|
| Regular / regular | Original coordinate over the whole interval | G7/K15 |
| Tangent / regular | `x = xL + L s^2`, `s in [0,1]` | G7/K15 |
| Regular / tangent | `x = xR - L s^2`, `s in [0,1]` | G7/K15 |
| Tangent / tangent | Two halves, with `delta = L/2` from the respective endpoint | G7/K15 on each half |

`OuterIntervalQuadrature.hpp` implements the common plan, exact map and Jacobian, embedded error estimate, independent volume/surface tolerances, selective masks, and adaptive subdivision. Regular children stay in the original coordinate. Mapped children retain the original endpoint map and subdivide its parameter. Accepted components are excluded from child evaluations. Certified higher contact exponents remain 4 or 6 and are counted; they are not replaced with 2.

There is no L/6 strip, central bulk, normal K31 escalation, or Gauss16/32 escalation in the native outer hierarchy. The historical `EndpointStripQuadrature` name remains for the cached general-curve evaluator, whose intervals now use the common full-interval planner. HEX calls the updated QUAD implementation; TET calls the updated full TRI implementation. `LineIntegration.hpp` and the exact LINE3 integration are unchanged.

Before evaluating a mapped K15 panel, the closest node is mapped with `fma` and checked against the endpoint. A collapsed node triggers K7, then K3. Each emergency rule must pass the embedded estimate, propagated inner errors, and an adjacent-coordinate sensitivity check. The probe goes toward the endpoint when that neighbor is distinct; otherwise it goes one representable point away. It never evaluates the singular endpoint. Its geometry parameter is recomputed while its quadrature Jacobian stays fixed, so rounding uncertainty is not mistaken for a change of integration measure. Only unresolved components continue to the next emergency rule. A precision failure never triggers repeated bisection toward the endpoint and never silently discards a contribution: the integration throws an explicit coordinate-precision error when no rule establishes tolerance.

Projection classification now distinguishes regular events from true fold/contact events, filters irrelevant folds outside the slice domain, and recognizes simple boundary fold onsets. A coefficient-scaled stationary-contact check prevents a rounded event at middle-coordinate zero from being mistaken for a regular endpoint. Isolated definite face contacts with nonzero parent normal derivative remain regular. Existing compensated 2D discriminant arithmetic and the nonzero discriminant constant for near-boundary folds are retained.

## Reproduction

Build the native targets with the standalone CMake project and run CTest as described in `README.md`. `cutfem_outer_interval_contract` tests all four interval plans, fourth/sixth-order maps, selective refinement, independent H/D budgets, emergency K7/K3 acceptance, inaccurate emergency rejection, and explicit below-resolution failure in both dimensions.

The common degree-10 driver uses identical nodal fields and physical geometry in both checkouts. Build `outer_interval_benchmark.cpp` against the f5f81a8 headers and the native FE support library for the baseline; against the new headers, CMake supplies `cutfem_outer_interval_benchmark`. To count the old panels without changing numerical formulas, apply `instrument_outer_baseline.py BASELINE_CHECKOUT` and compile that driver with `-DFEMUS_BASELINE_OUTER_COUNTERS`. The script only accepts the pinned baseline, adds guarded integer counters to three headers, and refuses unrelated modifications.

```sh
python3 unittests/testCutFemWeightLSNative/instrument_outer_baseline.py /path/to/baseline
python3 unittests/testCutFemWeightLSNative/run_outer_interval_benchmark.py --baseline /path/to/baseline-benchmark --candidate build-cutfem/cutfem_outer_interval_benchmark --output unittests/testCutFemWeightLSNative/unified_outer_results
```

The benchmark pins both sequential processes to one available CPU and records the median of three integrations. The raw JSON includes every H/D moment through total degree 10. The comparison CSV records maximum absolute/relative errors, worst powers, root-sample counts, regular/mapped K15 panels, emergency panels, subdivisions, collapses, and 3D slice calls. Legacy K31 and Gauss16/32 panels are recorded separately. Empty baseline mask/plan counters mean unavailable, not zero work.

The independent references use exact simplex/plane formulas, shifted ball moments, Wick expansions for rotated ellipsoid volume/delta moments, and independent angular or graph parameterizations for surface metrics and curved geometry. Two orders of the independent numerical reference must agree before a result is accepted. The 180 cases per checkout cover QUAD, TRI, HEX and TET; planes and oblique/slab cuts; near planes and curved graphs; balls, ellipses/ellipsoids and rotations; near tangencies and tiny cuts; physical curved metrics; forced directions; and tolerance sequences from 1e-6 through 1e-12. No physical monomial above total degree 10 is tested.

## Changed implementation

The production changes are `OuterIntervalQuadrature.hpp` (new shared controller), `LineDomainIntegration.hpp`, `EndpointStripQuadrature.hpp`, `SurfaceEndpointIntegration.hpp`, and `EndpointProjection.hpp`. Public FEMuS entry points, exact LINE3 kernels, FE reconstruction, physical metrics and moment indexing are preserved. Test/build/documentation changes and validation artifacts are listed by the accompanying commit.

## Removed code

Removed `StripWidth`, `MappedStrip`, the superseded `NearTangent` helper, `MappedPoint`, `MappedDiscriminant`, `MappedDerivativeMagnitude`, `GaussKronrodPanel`, `GraphGaussPairPanel`, duplicate `Panel`/`ResolvedPanel` implementations and their old adaptive loops; sine endpoint map; duplicate outer panel/adaptive loops; normal K31 and Gauss16/32 escalation; dead counters for those paths. Removed `endpoint_strip_circle.cpp`, `endpoint_strip_benchmark.cpp`, `run_endpoint_strip_circle.py`, `run_endpoint_strip_comparison.py`, `run_endpoint_strip_sanitizers.py`, and `run_endpoint_strip_validation.py`. Geometry coverage is provided by the common degree-10 driver, retained circle/ellipse regressions, and the new contract tests. The retained mapped-discriminant and single-root tests now exercise the current implementation. Historical results remain explicitly historical.

## Results

Both revisions completed and passed all 180 degree-10 cases. All monomial components through degree 10 were checked. Maximum errors across all fields, metrics and directions are:

| Element | Baseline H error | New H error | Baseline D/Gamma error | New D/Gamma error |
|---|---:|---:|---:|---:|
| QUAD | 2.318e-15 | 6.384e-16 | 4.936e-13 | 4.858e-15 |
| TRI | 2.345e-15 | 1.804e-16 | 4.935e-13 | 8.882e-15 |
| HEX | 6.245e-16 | 6.245e-16 | 1.167e-14 | 6.799e-15 |
| TET | 1.362e-16 | 6.939e-17 | 1.167e-14 | 6.800e-15 |

`unified_outer_results/comparison.csv` contains absolute/relative maxima, worst powers, all work counts, and every per-case runtime. Large relative errors of tiny moments are judged with the requested combined absolute/relative tolerance, not a relative-only criterion. The compressed `moment_records.jsonl.gz` retains every unrounded numerical component emitted by both drivers, together with the counters. `metadata.json` pins the compiler, CPU and header hashes.

The table below sums each case's median integration time and work, excluding full/empty phases. It is a workload comparison, not an application wall-clock prediction.

| Element | Baseline summed medians (ms) | New summed medians (ms) | Baseline root samples | New root samples | Baseline 2D slices | New 2D slices |
|---|---:|---:|---:|---:|---:|---:|
| QUAD | 62.854 | 70.834 | 2721 | 3060 | 0 | 0 |
| TRI | 55.881 | 60.196 | 2751 | 2820 | 0 | 0 |
| HEX | 10877.290 | 14208.011 | 175952 | 243300 | 4452 | 2586 |
| TET | 12778.621 | 17765.514 | 188630 | 256560 | 5304 | 3225 |

Representative degree-10 cases at atol=rtol=1e-10, automatic direction:

| Element/case/metric | Baseline (ms) | New (ms) | Baseline/new root samples | Baseline/new slices |
|---|---:|---:|---:|---:|
| QUAD plane, metric 0 | 0.088 | 0.125 | 15/15 | 0/0 |
| TRI graph, metric 2 | 1.040 | 0.346 | 63/15 | 0/0 |
| HEX graph, metric 2 | 38.254 | 29.590 | 945/675 | 15/15 |
| HEX ball, metric 1 | 179.136 | 181.159 | 3420/4500 | 135/60 |
| HEX tiny_ball, metric 1 | 180.848 | 138.181 | 3420/2700 | 135/60 |
| HEX ellipse, metric 2 | 807.165 | 1485.906 | 10488/22500 | 228/180 |
| TET tiny_ball, metric 1 | 174.804 | 137.603 | 3420/2700 | 135/60 |
| TET ellipse, metric 2 | 1060.481 | 2187.236 | 10488/22500 | 228/180 |

Metric 0 is raw delta; 1 is reference surface; 2 is the physical curved surface metric. The sampled HEX and TET workloads use approximately 42% and 39% fewer 2D slices, respectively. Total root work can increase because transformed degree-10 moments have higher degree in the map parameter and the embedded G7 error can require subdivision even when K15 is accurate. The old implementation could instead accept K31 or a Gauss16/32 pair. Removing that escalation, as requested, makes some curved ellipsoid cases approximately 1.8-2.1 times slower. That regression is retained and reported; no tolerance was relaxed and no higher normal rule was restored. The summed medians above expose the resulting workload regression; the reduction in slice calls alone is not a claim of universal speedup.

There were no coordinate-collapse or emergency-rule events in the 180 ordinary benchmark cases. Constructed resolution tests exercise both dimensions separately: K7 is accepted at delta=4096 long-double epsilons and K3 at 128 epsilons for a constant physical integrand. At one epsilon all candidates collapse and an explicit precision error is required. A degree-10 offset-scaled moment and a physically sampled singular kernel at 4096 epsilons force inaccurate emergency K7/K3 rejection. At 64 epsilons K3 accepts a constant physical integrand and rejects the singular kernel when tolerance is not established. No contribution is silently zeroed. Fourth- and sixth-order analytic maps and eight actual HEX/TET higher-contact comparisons also pass through moment degree 10.

Additional validation:

- All 25 CTests pass, including all 22 original baseline tests, the shared interval contract, retained circle/ellipse/metric regressions, and cached-root tests.
- The degree-10 surface regression driver passes 24,364 checks across 44 plane, scale, sphere and metric cases; its maximum observed error is 1.78e-15.
- The eight degree-10 higher-contact comparisons pass with maximum H difference 4.65e-15 and D/Gamma difference 2.54e-14 against a tighter independent projection.
- AddressSanitizer and UndefinedBehaviorSanitizer pass the interval and exponent tests. LeakSanitizer is disabled because this execution environment uses ptrace.
- Exact LINE3 integration and `LineIntegration.hpp` are byte-identical to f5f81a8 (SHA-256 `877452cf6e15ec3b4f429cd9a6a1e7960a7d70e2a45e1c3053a0f9d8b4b6d741`).

The actual `MultiPhase_ex44_test` target was compiled with the FEMuS library, MPI and PETSc. The checked-in `run_ex44_test.sh` completed for all four elements through level 4 with its degree-8 setting, first unchanged for the curved ball and then with only `CURVED=0` for the straight ball. The checked-in application's local absolute and relative moment tolerances are both 1e-5; the standalone degree-10 validations use the stricter tolerances stated above. Each application run produced 80 reference comparisons; all 320 baseline/new, straight/curved rows are in `ex44_refinement.csv`.

| Application geometry | Baseline maximum volume error | New maximum volume error | Baseline maximum surface error | New maximum surface error |
|---|---:|---:|---:|---:|
| Straight | 2.379e-10 | 1.719e-11 | 4.917e-09 | 1.122e-08 |
| Curved | 1.843e-10 | 2.144e-11 | 1.265e-08 | 1.264e-08 |

The application surface errors include physical mapping and moment fitting and are global sums; they are separate from the native per-moment integration errors. No refinement run aborted.

Limitations: tolerances below coordinate precision are explicitly unsupported; genuinely non-transverse/repeated level sets still require a mathematically supported topology and can be rejected. Higher contacts retain their certified exponent rather than being represented as simple folds. Timings are CPU/compiler-specific, use three repetitions, and describe these inputs only. The remaining curved-ellipsoid refinement cost is visible in the per-case table.

An inherited limitation also remains in the special TET barycentric reparameterization (`integrationAxis=1`), which tightens the internal budget by `4^q`. At q=10 and a physical curved metric, forced graph, sphere and rotated-ellipsoid probes fail to establish the requested tolerance in both f5f81a8 and this revision. They are listed separately in `tet_permutation_limits.json`; they are not counted as passing cases in the 180-case suite. Automatic direction selection and the ordinary Cartesian directions used by that suite pass.

The final commit is the branch commit containing this report and its source-hash manifest; its SHA is supplied in the accompanying completion message.
