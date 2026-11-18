# Deferred Bernstein restriction experiment

Deferring child Bernstein restriction accelerates monomial-only consumers, but gives no convincing improvement to the existing complete weight computation. The median ratio across 58 weight cases is 1.004x. This experiment is retained separately on barbi-ex47-lazy-bernstein; the integration branch remains at its eager implementation.

## Implementation

A parent captures its Bernstein coefficients and error padding once in a shared immutable snapshot. Each child retains that snapshot and at most two coordinate restrictions. Restrict() computes power coefficients immediately and defers Bernstein coefficients until bounds, classification, boundary analysis, or explicit Bernstein access requests them. The resulting child coefficients are cached. Parent reset/destruction cannot alter surviving children. Root snapshot preparation still materializes the root Bernstein representation once, as the eager implementation did; this trial defers child restriction, not root conversion.

The fixed-size chain avoids allocating a new descriptor for every slice. Intermediate Bernstein restriction is reconstructed along the chain when requested; it is not independently cached as a shared node. Known zero-boundary cases retain the existing eager pinning and error padding. Quadrature rules, moments, tolerances, and surface-density definitions are unchanged.

## Restriction benchmark

Ratios are eager / deferred; above one favors deferred. Five interleaved runs of 20,000 operations after warmup, pinned to CPU 1. All 22 checksum rows match exactly. Tensor and simplex coordinates were chosen so root cases intersect the level set. The isolated driver uses the benchmark functions stored in ex47.cpp and the same FEM basis objects. LINE is an unchanged evaluation control and has no restriction optimization.

| Element | Restrict + power access | Restrict + roots | Restrict + bounds | Nested + power | Nested + bounds |
|---|---:|---:|---:|---:|---:|
| QUAD9 | 1.575x | 1.065x | 1.163x | - | - |
| TRI7 | 2.138x | 1.058x | 1.059x | - | - |
| HEX27 | 1.824x | 1.076x | 1.025x | 1.789x | 1.131x |
| TET15 | 3.694x | 1.009x | 0.992x | 3.008x | 0.994x |
| WEDGE21 | 2.457x | 0.931x | 1.072x | 2.170x | 1.104x |

## Complete weight computation

The existing ex47 benchmark alternates two retained elements, clears final weights, checks checksums against the production class, and reports five-batch medians. Both eager and deferred binaries use Release -O3 builds. These are separate process runs pinned to CPU 0, rather than five interleaved full-suite runs; process/environment variation limits small percentage conclusions. The unchanged production class runs inside both processes and supplies a control.

| Element | Cases | Eager / deferred median | Production control ratio |
|---|---:|---:|---:|
| LINE3 | 12 | 0.979x | 1.007x |
| QUAD9 | 8 | 0.958x | 0.964x |
| TRI7 | 8 | 1.049x | 1.009x |
| HEX27 | 10 | 1.059x | 1.076x |
| TET15 | 10 | 0.993x | 1.026x |
| WEDGE21 | 10 | 1.000x | 1.019x |

The apparent HEX improvement also occurs in the unchanged production control. These measurements do not establish an integration speedup. QUAD/TRI results also vary, although their integration path does not rely on the changed recursive slicing.

## Why the integration gain disappears

PlanarEndpointIntegration::Integrate calls data.GetBounds() immediately. A 3D slice therefore requests the Bernstein coefficients just deferred by Restrict(). Graph certification and restricted-boundary analysis can also request them. Root analysis of restricted LINE data requests boundary information, which can materialize Bernstein coefficients as well. Moving this work later does not eliminate it for these consumers.

The monomial-only improvement remains useful for callers that evaluate or restrict a slice without asking for Bernstein-derived metadata. Changing the current integration consumers to avoid those requests would be separate work and was not included in this trial.

## Accuracy and lifetime checks

- 5,902,997 main checks pass, including new deferred-state, copy, parent reset/destruction, nested ancestor reset, and cache-reuse checks.
- Maximum scaled data error remains 2.3646164e-14, matching the eager baseline.
- 32,513 weight/moment checks pass across all six elements, including degree-13 HEX/TET/WEDGE moments.
- Address and undefined-behavior sanitizers pass 157,059 checks; leak detection was disabled.
- Existing allocation-instrumented data benchmarks report zero warmed heap allocations.

## Reproduction

Build this branch using the standalone ex47 CMake target in Release mode. For an eager comparison, create a worktree at 4ec74de89 and copy this branch's ex47.cpp into that worktree before building: the benchmark calls use the same existing public API. The deferred-ownership main assertions are intended for the trial build; run the benchmark options on the eager comparison.

```sh
python applications/MultiPhase/ex47/native_tests/benchmark_deferred_bernstein.py --eager /path/to/eager/MultiPhase_ex47 --lazy /path/to/lazy/MultiPhase_ex47 --output /tmp/bernstein-results
```

The recorded micro results used a smaller extracted driver for faster compilation; the same benchmark is available as --benchmark-bernstein in the main. The runner also accepts --eager-micro, --lazy-micro, and --micro-cpu. Raw logs, per-case timings, compiler flags, binary hashes, and source hashes accompany this report.
