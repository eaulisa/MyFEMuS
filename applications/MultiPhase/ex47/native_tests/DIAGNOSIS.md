# Accuracy and speed findings

The local ex47 prototype still has numerical regressions in 3D. The original production implementation and its tests are unchanged. The complete comparison is in `results/REPORT.md`; individual invocations and timings are in the CSV files beside it.

## Corrected certified-boundary roots

The first 386-case comparison flagged eight TRI/TET zero-edge/face cases. A root known to lie on a zero boundary could round one ULP inward and receive full delta weight instead of half weight. One TRI result was approximately `0.503556714` instead of `0.5`.

`GetLineAnalysis()` now pins a root only when exact boundary metadata proves its location. The prepared graph sampler also preserves that certified boundary coordinate. Arbitrary strip bounds are checked against their own boundary polynomial. All eight cases now pass, and all 386 public-class cases pass 36,860 original assertions per backend. `boundary_fix.json` retains the before/after assessments.

## Remaining prototype differences

| Fixture | Result | Interpretation |
|---|---|---|
| Curved HEX27, degree 10, first reflected case, raw delta | Prototype mass `1.5707910280271792`; reference `1.5707963267948966` at requested tolerance `1e-10` | About `5.30e-6` error; adaptive retries exhaust the evaluation budget and retain an uncertified result. |
| Same nodal fixture through actual public weight APIs, degree 2 | Original passes; prototype rejects an endpoint at coordinate resolution | The failure also exists through the class interface; it is not confined to a test adapter. |
| Enriched HEX bubble, degrees 2/5/10 | Delta mass `34.04504412723863` versus `34.04504408913909` | About `3.81e-8` error against the original unchanged oracle limit of about `7.01e-9`; fitted weights show the same underlying discrepancy. |
| Degree-10 nested cubic HEX contacts | Original completes; prototype rejects two cubic-contact cases | The higher-order endpoint path is incomplete for these projections. Default degree-2 nested-contact checks pass. |
| Outer-interval HEX graph, curved physical metric, forced order 1, degrees 2/10 | Original passes; prototype reports no complete finite outer rule at coordinate resolution | Automatic selection and the other tested orders succeed for the same field/metric. |
| Unattainable-tolerance warning contract | Prototype LINE weights retain the correct mass but omit the expected aggregated warning | The LINE integration routine produces error estimates; its public wrapper discards them without a tolerance check. This is a diagnostic contract failure. |

A control run of the original native-coordinate HEX kernel also fails on the perturbed curved fixture, while the original unit-coordinate path passes. This supports a chart/endpoint-conditioning explanation. Native coordinates place some small offsets next to `-1` or `1`, where representing an absolute coordinate can lose information that remains representable next to zero. The exact offending event in each failing case has not yet been isolated.

One concrete code path to inspect is `GetFiberProjection()`: cached discriminant coefficients are formed with ordinary long-double products and sums. Compensated LINE evaluation occurs later and cannot restore precision already lost when forming a projection event. This is a candidate contributor, not a verified complete root cause. Preserve event offsets and parameter derivatives through both slicing levels, and validate the fixes against these exact fixtures before production integration.

## Shared limitations and detection configuration

The forced nearly parallel TRI plane direction is rejected by both implementations. HEX sphere refinement at division 4, degree 6/10 and tolerance `1e-11` also fails the strict reference in both implementations. These are recorded as shared failures. The exact 32 CTest invocations pass in the original; the two CTest sphere-mesh requests use tolerance `1e-9` and pass in both implementations.

The old detector uses Bernstein depth 12. The prototype default is depth 1, as requested earlier. A strictly positive polynomial can therefore remain uncertified and be treated as a possible cut. The legacy strict-certification test fails under that default and passes with depth 12 selected. This is a conservative configuration difference, not a missed small cut or a wrong integral.

## Measured performance

All binaries use GNU C++14 with `-O2 -DNDEBUG`. Public API batches are pinned to one CPU, interleaved three times, and clear final-weight caches before each call. Prepared polynomial/geometry data and the fixed rule are retained. Constructor times and sample ranges are exported separately. Drivers with their own repeated timers retain those timers; public API timing excludes independent-oracle work.

The median public-case speedups are approximately 5.7× for LINE weights, 1.1–1.8× for 2D/3D raw/reference calls, and 2.1–5.1× for 3D physical-surface calls. Across the 58 retained-element/surface-family rows, 3D physical surfaces are 3.1–4.7× faster and the three-normal family is 3.7–5.3× faster. Raw HEX ball delta is about 8% slower at degree 2 and 2% slower at degree 5.

The clearest regressions occur on inexpensive certified boundary cuts. For the degree-0 QUAD lower-boundary case, a raw/reference call takes about 15.2 µs originally and 32.4 µs in the prototype. Both use 15 quadrature samples; the prototype additionally restricts 15 LINE records and analyzes roots whose boundary location is already known. The upper-boundary and degree-2/5 cases show the same pattern. The relevant path is `FastGraphSample()` → `PreparedPlane::Fiber()` → `Restrict()` → `GetLineAnalysis()`.

Use a certified-boundary/known-root path that avoids unnecessary LINE topology and restriction work, and defer lower-dimensional metadata that the requested density does not need. Retain the general and compensated paths for difficult cuts. This optimization remains pending; the comparison does not claim that the regressions have been removed.

The failing 3D projections also incur selective retries, very small local tolerance allocations, and evaluation/coordinate-resolution exhaustion. Their wall time is retained in the report, but no accepted speedup is assigned to an inaccurate or unsupported result.

A fixture that aborts on its first failed assertion does not complete its remaining internal subcases. The report records source-level fixture outcomes, the completed parameterized requests, and these failures explicitly.
