# Regularized endpoint quadrature: full QUAD/TRI comparison

10 November 2026. All validation passed. Maximum requested total monomial degree: 10.

## Implementation

The endpoint experiment is now called by the existing unit QUAD and standard TRI interfaces only for eligible non-graph requests. Certified graph inputs, prescribed directions, independent-direction validation, and unsupported boundary or singular configurations retain their existing paths. Integration proceeds with the actual FE level-set polynomial; no parabolic interface approximation or model-refinement loop is used.

Tangencies and cell-boundary events define the outer intervals. For each regular fold, a strip uses x = a + delta*u^2 or x = b - delta*u^2, with algebraic cancellation of the vanishing Jacobian/discriminant factors. Automatic delta is one sixth of its event-interval width, capped to prevent overlap; this gives delta approximately 0.05 for the radius-0.15 circle. The regular bulk keeps the original coordinate. Cached roots and topology are shared by all moments. GK15 is followed by GK31 on remaining components, then bisection only for components still failing.

New non-graph unit-basis absolute and relative budgets are 0.01 times the caller budget, protecting the conversion to reference moments. The comparisons request absolute and relative tolerance 1e-10. The retained baseline uses its original caller budget. Tolerances were not relaxed.

## Coverage and accuracy

| Check | Coverage | Result |
|---|---|---|
| Pinned baseline comparison | 405 cases; 37,944 unit and signed-reference moment checks; q=2,5,10; raw delta, reference surface, and curved physical surface | PASS |
| New method completion | 207 eligible calls, 207 completed, zero fallback | PASS |
| Graph preservation | 180 graph comparisons, bit-for-bit identical moments | PASS |
| Wider fitted-weight matrix | 173 QUAD/TRI cases; 31,068 checks; sign reversal, scaling, directions, tolerances 1e-8 through 1e-12 | PASS |
| Physical mesh refinement | 36 configurations; 144 global volume/surface/x-squared checks; up to 4,096 QUAD or 8,192 TRI cells | PASS |
| Degree-ten double-output weights | 8 cases; 3,186 checks | PASS |
| Existing native QUAD/TRI cases | 155 cases; 10,814 checks | PASS |
| ASan and UBSan | All 405 comparison configurations | PASS |
| API, precision, output range, include order, LINE and standard QUAD/TRI moments | Seven builds; five runtime checks, including 2,654 standard-moment checks | PASS |

The numeric matrices report 83,156 checks in total, with API checks and sanitizer repetitions additional. The largest absolute error in the comparison is 8.38e-13 in the unit basis and 1.54e-12 in the signed FEMuS reference basis. The largest curved-physical-surface error is 3.94e-13. The largest observed error divided by its complete tolerance/oracle limit is 0.01461.

The 45 input cells are exactly the prior 23 QUAD and 22 TRI fixtures: ten graph families, circles with several boundary placements and sizes, rotated ellipse, two sheets, hyperbola, bubble, perturbed circle, tensor lobes, hyperbolic sheets, and nearly parallel sheets. The tensor-lobe polynomial contains x^2*y^2 and is not in the TRI7 space, so the prior TRI exclusion is preserved. TRI clipped-circle and hyperbola requests remain on their existing triangle graph paths; these account for the 18 additional calls outside the 207 new-method calls and 180 primary graph comparisons.

Independent references locate events at 80-digit precision, then use fixed 128- and 256-point Gauss rules with a different, whole-interval sine-squared map. Signed QUAD reference moments are integrated directly in the signed basis for the oracle. The test checks raw moment conversion and physical factors separately. The comparison limit is 1e-10*(1+abs(reference)) + 8*oracle_discrepancy + 3e-12. Wider fitted-weight checks retain their established oracle limits.

All mesh levels pass at all three degrees. Maximum QUAD mesh errors are 2.69e-14 for volume and 6.51e-14 for surface; maximum TRI errors are 1.87e-14 and 8.32e-14. These circle fields are represented exactly by the FE space, so errors are near rounding and a monotone convergence rate is not inferred.

## Timing method

Baseline commit: 1a394fdc97bedc5a6cff9b1e6c68f16261ba57e9, before either strip experiment. Both binaries were built from the same fixture with identical C++14/O3/NDEBUG/generic-x86-64 flags and actual FEMuS FE objects. Baseline integration headers come from a git archive of that pinned commit.

Host CPU: INTEL(R) XEON(R) PLATINUM 8573C; both processes pinned to CPU 0. Compilation and other numerical runs were complete before this comparison. For each cell, three outer rounds each take the median of five batches lasting at least 8 ms; the table reports the median of those three round values. Nearby nodal constants alternate in each batch to force geometry reconstruction and prevent cached-topology timing.

Times cover nodal polynomial reconstruction, direction/event work, metrics/cofactors, and volume plus surface moments. They exclude nodal input generation, all weight fitting, fitting-matrix setup, and output formatting. These are same-session kernel comparisons, not the full MPI/PETSc ex44 application.

## Degree 10: volume plus reference-surface moments

| Curve | QUAD baseline (ms) | QUAD mapped (ms) | Speedup | TRI baseline (ms) | TRI mapped (ms) | Speedup |
|---|---:|---:|---:|---:|---:|---:|
| Circle | 5.2668 | 2.5570 | 2.06x | 4.5315 | 2.1634 | 2.09x |
| Clipped circle | 5.0646 | 2.8121 | 1.80x | 5.9966 | 3.0763 | 1.95x |
| Edge circle | 1.6857 | 0.6420 | 2.63x | 6.9806 | 3.0029 | 2.32x |
| Rotated ellipse | 11.2921 | 3.3209 | 3.40x | 17.2691 | 2.6367 | 6.55x |
| Two sheets | 0.9827 | 0.1293 | 7.60x | 7.6274 | 0.5268 | 14.48x |
| Hyperbola | 4.3614 | 7.9244 | 0.55x | 3.2719 | 2.5995 | 1.26x |
| Bubble | 16.7532 | 4.6843 | 3.58x | 9.8170 | 5.9812 | 1.64x |
| Perturbed circle | 22.3128 | 2.2262 | 10.02x | 30.6611 | 2.7082 | 11.32x |
| Tensor lobes | 39.5854 | 18.1594 | 2.18x | — | — | — |
| Hyperbolic sheets | 5.5864 | 1.3437 | 4.16x | 4.9225 | 1.6454 | 2.99x |
| Nearly parallel sheets | 2.0302 | 0.5096 | 3.98x | 4.9134 | 0.9369 | 5.24x |
| Corner circle | 3.7104 | 1.9342 | 1.92x | 4.9045 | 2.3562 | 2.08x |
| Tiny circle | 5.2558 | 2.3216 | 2.26x | 2.3031 | 2.3747 | 0.97x |

The main remaining regression is the QUAD hyperbola: 4.36 to 7.92 ms, about 1.82 times slower. It has no mapped tangent strip in the chosen projection and uses 506 bulk samples with four bisections. The comparison retains this case; no per-geometry exclusion was added after timing. The tiny TRI circle is about 3% slower in this run. Differences on paths retained from the baseline, including the TRI clipped circle and hyperbola, should not be attributed to the new endpoint map.

## Aggregate reference-surface time

Sum of the per-cell median times over the 13 QUAD and 12 TRI non-primary-graph fixtures; each row uses the same cells for its baseline and candidate.

| q | QUAD baseline (ms) | QUAD mapped (ms) | Speedup | TRI baseline (ms) | TRI mapped (ms) | Speedup |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 68.453 | 34.713 | 1.97x | 63.739 | 22.271 | 2.86x |
| 5 | 80.938 | 43.151 | 1.88x | 65.681 | 23.793 | 2.76x |
| 10 | 123.888 | 48.565 | 2.55x | 103.199 | 30.008 | 3.44x |

## Degree 2: reference surface

| Curve | QUAD baseline (ms) | QUAD mapped (ms) | Speedup | TRI baseline (ms) | TRI mapped (ms) | Speedup |
|---|---:|---:|---:|---:|---:|---:|
| Circle | 3.9122 | 1.7856 | 2.19x | 3.0014 | 1.6954 | 1.77x |
| Clipped circle | 3.1417 | 2.0106 | 1.56x | 2.4162 | 2.2540 | 1.07x |
| Edge circle | 0.9275 | 0.3797 | 2.44x | 3.4699 | 2.2534 | 1.54x |
| Rotated ellipse | 7.2964 | 2.7664 | 2.64x | 11.0938 | 2.1320 | 5.20x |
| Two sheets | 0.4788 | 0.0614 | 7.80x | 10.6852 | 0.4101 | 26.06x |
| Hyperbola | 2.3884 | 4.4969 | 0.53x | 2.8471 | 1.9959 | 1.43x |
| Bubble | 9.1633 | 4.1500 | 2.21x | 5.7077 | 4.6251 | 1.23x |
| Perturbed circle | 10.2297 | 1.8405 | 5.56x | 11.7586 | 1.6304 | 7.21x |
| Tensor lobes | 19.6611 | 12.3514 | 1.59x | — | — | — |
| Hyperbolic sheets | 3.9420 | 1.0953 | 3.60x | 3.3890 | 1.1234 | 3.02x |
| Nearly parallel sheets | 1.0576 | 0.3818 | 2.77x | 3.2714 | 0.7964 | 4.11x |
| Corner circle | 2.1815 | 1.5056 | 1.45x | 4.2320 | 1.6843 | 2.51x |
| Tiny circle | 4.0730 | 1.8880 | 2.16x | 1.8663 | 1.6702 | 1.12x |

## Degree 5: reference surface

| Curve | QUAD baseline (ms) | QUAD mapped (ms) | Speedup | TRI baseline (ms) | TRI mapped (ms) | Speedup |
|---|---:|---:|---:|---:|---:|---:|
| Circle | 4.7160 | 2.0366 | 2.32x | 3.1585 | 1.8517 | 1.71x |
| Clipped circle | 3.4913 | 2.2463 | 1.55x | 2.3710 | 2.7027 | 0.88x |
| Edge circle | 1.1271 | 0.4406 | 2.56x | 3.2770 | 2.4214 | 1.35x |
| Rotated ellipse | 8.2991 | 5.5508 | 1.50x | 11.1020 | 2.4149 | 4.60x |
| Two sheets | 0.5891 | 0.0713 | 8.26x | 7.9956 | 0.4192 | 19.07x |
| Hyperbola | 2.6731 | 4.7485 | 0.56x | 3.0435 | 2.3786 | 1.28x |
| Bubble | 12.6035 | 4.5143 | 2.79x | 6.4661 | 4.1567 | 1.56x |
| Perturbed circle | 11.2068 | 1.8923 | 5.92x | 14.1583 | 1.7522 | 8.08x |
| Tensor lobes | 23.8329 | 16.4748 | 1.45x | — | — | — |
| Hyperbolic sheets | 4.5551 | 1.1248 | 4.05x | 3.7286 | 1.2370 | 3.01x |
| Nearly parallel sheets | 1.2393 | 0.4115 | 3.01x | 4.0538 | 0.8300 | 4.88x |
| Corner circle | 2.2108 | 1.6405 | 1.35x | 4.1846 | 1.7565 | 2.38x |
| Tiny circle | 4.3939 | 1.9986 | 2.20x | 2.1415 | 1.8725 | 1.14x |

## Degree 10: curved physical surface

| Curve | QUAD baseline (ms) | QUAD mapped (ms) | Speedup | TRI baseline (ms) | TRI mapped (ms) | Speedup |
|---|---:|---:|---:|---:|---:|---:|
| Circle | 6.3807 | 2.9211 | 2.18x | 4.4466 | 2.6608 | 1.67x |
| Clipped circle | 6.3497 | 3.3845 | 1.88x | 4.7611 | 4.0264 | 1.18x |
| Edge circle | 2.6047 | 1.2565 | 2.07x | 7.6376 | 3.8594 | 1.98x |
| Rotated ellipse | 16.3711 | 4.3264 | 3.78x | 23.5565 | 3.6778 | 6.41x |
| Two sheets | 1.3139 | 0.2081 | 6.31x | 9.0684 | 0.6545 | 13.86x |
| Hyperbola | 4.4991 | 8.6246 | 0.52x | 3.0495 | 3.4316 | 0.89x |
| Bubble | 22.6440 | 7.4208 | 3.05x | 10.0543 | 6.5701 | 1.53x |
| Perturbed circle | 32.1709 | 3.1452 | 10.23x | 24.1906 | 2.5555 | 9.47x |
| Tensor lobes | 37.3476 | 27.5523 | 1.36x | — | — | — |
| Hyperbolic sheets | 8.5914 | 4.3337 | 1.98x | 5.2750 | 1.9590 | 2.69x |
| Nearly parallel sheets | 2.1722 | 0.6580 | 3.30x | 6.2073 | 1.1312 | 5.49x |
| Corner circle | 4.4359 | 2.7365 | 1.62x | 6.1819 | 3.2102 | 1.93x |
| Tiny circle | 7.3930 | 3.0488 | 2.42x | 3.5725 | 2.8659 | 1.25x |

## Reproduction and evidence

The complete per-case timings, sample counts, split counts, errors, source hashes, and CPU settings are in endpoint_strip_comparison.json. The acceptance file endpoint_strip_acceptance.json names all supporting reports. The raw per-case CSV outputs are reproducible in results/endpoint_strips.

With Boost, Eigen, a C++14 compiler, Python NumPy, and mpmath available:

```sh
python3 run_endpoint_strip_validation.py --boost /path/to/boost/include --eigen /path/to/eigen3
```

The driver builds the required executables, compares both variants, repeats the wide fitted-weight and mesh cases, runs sanitizers and native/API tests, and checks final source hashes. ASan and UBSan passed here; leak detection is disabled because LeakSanitizer cannot run in this ptrace environment. No degree above ten was requested. Full MPI/PETSc ex44 was not built or run.
