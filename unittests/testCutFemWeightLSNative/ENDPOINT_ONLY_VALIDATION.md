The current outer integration update and acceptance evidence are in `UNIFIED_OUTER_QUADRATURE.md`; the results below describe the earlier endpoint-only revision.

# Native CutFEM endpoint-only validation

## Source and scope

New branch: `rifat-cutfem-endpoint-only` in `eaulisa/MyFEMuS`.

The starting tip was verified before implementation and again before publication: **`8bad259c5a51060e4c2674db49a7ae34b75a0367`**, on `rifat-cutfem-2d-endpoint-quadrature`. Its parent is the previously accepted `aea2c13efd1ac883103461d3218790cfbe70e4f0`. The unchanged `rifat` benchmark is **`0ff0669494e603be248d02c44c0a2e9b81f72b4d`**. This commit is based directly on the verified tip; the comparison JSON records that base commit plus the exact candidate header SHA256 manifest.

The native hierarchy is **HEX → QUAD → LINE** and **TET → full TRI → LINE**. The original parent finite-element polynomial is restricted algebraically. No polynomial is fitted from slice samples, and quartic projected constraints are retained in dynamically sized polynomial arrays. FEMuS TET15 face and interior enrichment is retained. Barycentric direction changes permute the nodal basis and transform moments exactly.

Retained paths are the LINE quadratic kernel and graph kernel, planar graph integration, planar endpoint strips with actual-integrand quadrature, and full parent 3D slices with projection events and endpoint quadrature. Topology intervals, endpoint strips, and adaptive quadrature panels remain. There is no spatial partition tree, geometric subdivision, Taylor dispatcher, Full GK fallback dispatcher, or Algoim dependency in the native path. GK rules within the retained algorithms remain.

Selective masks stop work for moments whose accumulated error passes. The necessary slice dependencies remain active. The graph dispatcher also checks accumulated errors after cancellation and retries only failed outputs. Known trace roots and zero-face/edge metadata, discriminants, branch signs, and interval topology are reused where valid; roots at a new quadrature coordinate are evaluated there. Slice-coordinate, endpoint-map, and surface/cofactor weights are applied once. API-output rounding and strict rejection remain checked.

## Transformation exponent

See [ENDPOINT_ONLY_MATHEMATICS.md](ENDPOINT_ONLY_MATHEMATICS.md) for the derivation and analytical examples. The definition is `s=a+delta*t^p`, `0<=t<=1`, with Jacobian `p*delta*t^(p-1)`; the right strip reverses the limits.

For a density `s^alpha`, boundedness requires `p*(alpha+1)-1>=0`. An ordinary fold has `alpha=-1/2`, so the minimum integer p is 2. For the **transverse** projected model `D(u,s)=c*s+d*u^m+...`, with nonzero c,d and nonzero quadratic leading coefficient, the integrated fold density has `alpha=1/m-1/2`.

| Contact multiplicity m | Minimum integer p for boundedness of that slice density | Sufficient power clearing both algebraic branches |
|---|---:|---:|
| 2 | 1 | 2 |
| 3 | 2 | 6 |
| 4 | 2 | 4 |

The code uses the latter powers and their least common multiple at a common event; **12 is a sufficient common power for these stated local models**. It is not a universal theorem for every supported finite-element polynomial. HEX27 has a full tensor Q2 polynomial: its coordinate discriminant has degree four in each projected coordinate. The actual enriched TET15 space can also produce cubic and quartic projected contacts. Contact multiplicity determines these powers; steepness affects conditioning and constants. Arbitrarily larger powers can work, so there is no maximum admissible p. No finite p makes every possible singular or logarithmic input analytic or every divergent delta integral finite.

Vanishing A with nonzero B is handled as a linear fiber. An entire zero fiber needs a valid transverse chart. Persistent repeated factors and detected non-transverse general projected endpoints are rejected. Certified regular graph boundary traces are treated separately. A regular saddle can have a logarithmic slice density; power transformations can bound it but cannot make it analytic. No fallback conceals a rejection.

## Executed checks

All moment requests are at total degree ten or below.

| Check | Executed result |
|---|---|
| Standalone actual-FEMuS-source builds | 34 original targets passed; subsequent changed test utilities also rebuilt |
| Maintained CMake/CTest | **18/18 passed**, including LINE, QUAD/TRI, enriched TET, HEX, analytical p checks, nested endpoints, and sphere mesh regressions |
| Original native/fitting matrix | **386 cases, 36,860 assertions, zero failures** |
| Accepted planar geometry matrix | **405 cases, 37,944 checks, zero failures**; degrees 2,5,10; raw delta, reference surface, physical surface |
| Focused tests | All 19 invocations passed, including degree-8/10 simplex, directions, traces, scaling, strict rounding, physical geometry, and fitted weights |
| Actual FEMuS basis/gradient comparison | 15,360 checks; maximum absolute difference `7.135e-15` |
| Circle/ellipse convergence regressions | 144 circle cases at tolerances 1e-8/1e-10/1e-12, 48 translated ellipses; 12,116 combined checks passed |
| ASan + UBSan | QUAD/TRI direction/metric suite, nested HEX/TET contacts, and API contracts passed; no sanitizer reports; leak detection disabled |
| Removed API contracts | 29 obsolete signatures inaccessible; 157 contract checks passed |
| Original LINE kernel identity | Five retained kernel bodies are byte-identical to the accepted start; hashes in `line_kernel_identity.json` |

The 405-case matrix retains closed and clipped circles, ellipses, two sheets, hyperbolas, bubbles, perturbed circles, tensor lobes on QUAD, nearly parallel sheets, corner tangencies, tiny circles, and graph cases. The full TRI restriction is evaluated in its enriched basis. Native `[-1,1]` QUAD moments are requested directly and checked independently, including the degree-10 hyperbola. Postprocessing rounded unit monomials with a degree-10 binomial transformation is ill-conditioned and is not substituted for the direct native-coordinate check.

Maximum planar errors across the 405-case matrix were `2.970e-15` for unit volume and `8.380e-13` for interface moments across the three measures. Direct FEMuS-coordinate maxima were `9.407e-15` and `9.599e-13`. The independent references retain their doubled-order discrepancy in each record; the largest error-to-test-budget ratio was `0.009143`.

### Closed-interface tolerance checks

Degree 10; maximum absolute error over all requested monomials, against 60-digit ball/sphere formulas. HEX/QUAD/LINE radius is 0.6; TRI/TET radius is 0.125 and center coordinates are 0.25. Surface means reference surface measure, rather than raw delta.

| Element | Tolerance | Volume maximum error | Surface maximum error |
|---|---:|---:|---:|
| LINE | 1e-8 | 4.000e-18 | 1.000e-17 |
| LINE | 1e-10 | 4.000e-18 | 1.000e-17 |
| LINE | 1e-12 | 4.000e-18 | 1.000e-17 |
| QUAD | 1e-8 | 9.804e-16 | 1.441e-15 |
| QUAD | 1e-10 | 3.415e-17 | 1.862e-16 |
| QUAD | 1e-12 | 3.415e-17 | 1.862e-16 |
| TRI | 1e-8 | 1.604e-15 | 4.202e-15 |
| TRI | 1e-10 | 2.480e-17 | 9.038e-17 |
| TRI | 1e-12 | 1.065e-17 | 9.038e-17 |
| HEX | 1e-8 | 5.268e-17 | 6.263e-15 |
| HEX | 1e-10 | 5.268e-17 | 6.263e-15 |
| HEX | 1e-12 | 5.268e-17 | 6.263e-15 |
| TET | 1e-8 | 4.148e-18 | 2.260e-17 |
| TET | 1e-10 | 1.082e-19 | 2.260e-17 |
| TET | 1e-12 | 1.082e-19 | 2.260e-17 |

Errors already at input/output roundoff do not keep decreasing with tolerance. No asymptotic convergence order is inferred from that plateau.

### Nested projected endpoints

Forced non-graph cubic/quartic contacts are compared with a different, regular graph direction using a reference tolerance 100 times tighter. Both raw delta and reference surface are tested. At degree 10 and tolerance 1e-12:

| Element | Contact | p | Volume maximum error | Reference surface maximum error |
|---|---:|---:|---:|---:|
| HEX | 3 | 6 | 2.711e-19 | 3.795e-19 |
| HEX | 4 | 4 | 2.050e-19 | 1.247e-18 |
| TET | 3 | 6 | 4.470e-16 | 4.904e-19 |
| TET | 4 | 4 | 2.154e-18 | 8.586e-18 |

The degree-6/tolerance-1e-11 and degree-2/tolerance-1e-9 checks also pass. The model tests distinguish bounded p=2 substitutions from powers that clear the cubic or quartic branches and require rejection of a non-transverse model. An event-rounding regression checks that stationary points expose a fold even when the approximate event leaves a small polynomial residual.

### Mesh refinement

LINE: five meshes with 4–64 cells and degrees through 10; maximum volume error `5.421e-18`, surface error `1.129e-16`.

QUAD/TRI: 30 circle-mesh configurations, 4–64 divisions per axis, degrees 2,6,10. All pass; maximum area error `1.713e-17`, perimeter error `7.351e-17`.

HEX/TET: sphere radius 0.613 with translated center (0.071,-0.043,0.029), degree 2, tolerance 1e-10. The volume and surface x-squared moments are checked in addition to those shown:

| Element | Divisions | Cells | Volume error | Surface error |
|---|---:|---:|---:|---:|
| HEX | 4 | 64 | 1.143e-15 | 4.857e-17 |
| HEX | 8 | 512 | 6.451e-18 | 2.602e-17 |
| HEX | 16 | 4096 | 2.168e-18 | 7.373e-18 |
| TET | 4 | 384 | 5.660e-15 | 1.531e-16 |
| TET | 8 | 3072 | 9.575e-16 | 3.209e-17 |
| TET | 16 | 24576 | 1.528e-16 | 4.337e-19 |

These geometries are exactly represented by the finite-element spaces, so fine-grid errors approach roundoff rather than a geometry approximation rate. Mesh tests do not introduce spatial subdivision into the integrator.

## Same-session comparisons and timings

All **60 common cases pass on all three implementations**: degrees 2 and 10, all five element types, planar, nearly planar, closed, two-root/two-sheet, small-cap, and tiny closed interfaces. Maximum candidate/accepted differences are `2.359e-16` for volume and `6.052e-14` for surface; candidate/rifat maxima are `1.721e-15` and `6.217e-15`.

Compiler: GNU g++ 14.2.0, C++14, `-O3 -DNDEBUG -march=x86-64 -mtune=generic`; same Boost, sandbox x86_64 hardware, and pinned logical CPU 0. The host CPU model is not exposed; these are not measurements on the user's i9. Both absolute and relative tolerances are 1e-10. A warm-up integration precedes each timed process. Seven paired batches alternate revision order; all revisions use the same repetition count for a case. JSON formatting and process startup are excluded. Full timing samples, min/max/median, commands, binary hashes, and candidate header hashes are in `comparison_final.json`. These are repeated nodal integrations with a fixed input; fitted-weight setup is excluded. Shared-host scheduling causes timing variation, particularly for microsecond cases; small differences should not be interpreted as robust improvements.

Representative degree-10 median times, milliseconds per integration:

| Element/case | Candidate | Accepted start | rifat |
|---|---:|---:|---:|
| LINE / plane | 0.00644202 | 0.00708478 | 0.00630668 |
| LINE / two roots/sheets | 0.00742769 | 0.00605586 | 0.00606619 |
| QUAD / closed | 1.68047 | 4.23535 | 6.94269 |
| QUAD / two roots/sheets | 0.139403 | 1.87342 | 2.42402 |
| TRI / closed | 1.3217 | 1.26801 | 5.89497 |
| TRI / two roots/sheets | 0.487367 | 0.365772 | 13.6965 |
| HEX / plane | 0.227186 | 0.279558 | 0.318889 |
| HEX / closed | 226.967 | 2010.52 | 6284.68 |
| HEX / two roots/sheets | 0.254463 | 0.974206 | 8.20312 |
| TET / plane | 6.82781 | 7.2313 | 107.738 |
| TET / closed | 192.512 | 934.526 | 946.745 |
| TET / two roots/sheets | 510.754 | 305.3 | 14860.1 |

Across the twelve common cases per element, geometric-mean speedup (baseline time / candidate time):

| Element | Against accepted | Against rifat |
|---|---:|---:|
| LINE | 0.894x | 0.895x |
| QUAD | 1.928x | 1.774x |
| TRI | 0.933x | 1.800x |
| HEX | 2.212x | 2.808x |
| TET | 1.346x | 4.370x |

Regressions are retained in the data. In particular, the degree-10 TET two-sheet case is about 1.67 times slower than the accepted start (while much faster than rifat). TRI closed/two-sheet cases and some tiny cuts or planar requests regress against one baseline. The five-element harness reports a roughly 11% aggregate LINE slowdown despite identical numerical kernel bodies; direct-kernel measurements are provided separately to distinguish harness/layout variation from kernel work. The cleanup does not establish a universal speed improvement.

The separate direct-kernel benchmark uses only `FELevelSetMoments.hpp`, so the other element dispatchers do not affect the measured function. Its candidate/accepted executable `.text` sections are **byte-identical** (SHA256 `bb26227693f9a1f006a5e4343dec5cfe0844e60628d09cb9498b34aad1f4c074`). Two-root medians are 4.830/4.436 microseconds; observed ranges overlap (3.826–5.256 and 3.657–5.188 microseconds). Fast-graph medians are 62.102/66.006 microseconds. Identical machine code establishes that this cleanup introduced no LINE kernel work; timing variation remains visible rather than being reported as an algorithmic speedup. The corresponding rifat medians are 4.297 and 60.855 microseconds. See `line_comparison.json` and `endpoint_line_comparison.cpp`.

## Removed implementation and test paths

The removed native files are `AdaptiveLineSliceMoments.hpp` and `GraphCellIntegration.hpp`. The reusable polynomial utilities moved to `LineSliceUtilities.hpp`; retained geometry/range kernels now live with `SurfaceEndpointIntegration.hpp`. New projection and chart helpers are `EndpointProjection.hpp`, `EndpointSimplexChart.hpp`, and `PlanarEndpointIntegration.hpp`.

Removed method families include `Series*`, `SliceTaylor`, `IntegrateTaylorInterval`, `LineDomainIntegration::TaylorPanel`, `IntegrateCell`, the general sample/Full GK slice dispatcher, `IntegrateNormalizedByLineSlices`, `IntegrateNormalizedCell`, and their `useTaylor` options. Removed spatial methods include `Clip`, `Slab`, `Partition`, `Split`, `CriticalSplit`, `BoundSmallRegion`, `RefinedPhase`, TRI polygon/derivative/longest-edge subdivision, `Hex27SplitBernstein`, `Tet15SplitHeight`, and `Tet15SplitTriangle`, plus obsolete slice/slab event and integration wrappers. Uncalled line-order, full-moment, derivative, and old surface-slice helpers were removed after a repository caller audit.

Deleted test/tool files: `graph_cell_splitting.cpp` and `tri_domain_splitting.cpp` were replaced by `endpoint_surface_regressions.cpp` and `tri_endpoint_regressions.cpp`; their useful accuracy, metric, boundary, and difficult geometry cases remain. `apply_review_fixes.py`, `finish_review_fixes.py`, `run_direct_moment_benchmark.py`, and `run_unit_basis_benchmark.py` were removed. Three historical native-validation workflows were replaced by `cutfem-endpoint-only.yml`. Obsolete spatial-counter assertions and labels were removed. The old planar CutFEM library, `cutfemOld`, independently used components, Algoim as an independent library, and unrelated tests were left intact.

Source searches and compiler dependency files show no active native reference to the deleted methods/files or to the old CutFEM/Algoim implementation. All native validation Python files pass syntax compilation. No native code contains a Taylor, Full GK fallback, geometric split dispatcher, or calls into legacy implementations.

## Remaining limitations and incomplete requirement

- General projected inputs outside the stated transverse contact models are not universally certified. Persistent repeated factors, detected non-transverse endpoints, exhausted quadrature budgets, or coordinate-resolution loss return an explicit error. A different admissible direction may use the same endpoint algorithm; no alternate integration algorithm is retained.
- Repeated level-set factors can make raw delta moments undefined/divergent. There is no endpoint power that repairs them. The branch makes no general analytic-convergence claim for logarithmic saddle projections.
- The unchanged coefficient-only LINE kernel retains its old endpoint snap tolerance. A root at `1e-24` beside a root at zero is outside its accuracy claim and may receive the old half-boundary weight; the regression labels that case unsupported. It cannot be distinguished in the ordinary double nodal encoding of that example. Resolvable close interior roots are tested. This limitation was not hidden by rewriting the preserved LINE algorithm.
- **Full FEMuS/ex44 linking is incomplete in this environment.** Configuration with MPI/PETSc and actual dependencies succeeds, and the ex44 source object compiles against the final native headers. Building the full target stops in unchanged `LinearEquationSolverPetsc.cpp:359` and `LinearEquationSolverPetscFieldSplit.cpp:66,136,180`: installed PETSc 3.22.5 lacks `KSPOrthogonalizationSetCGSRefinementType` and `KSP_ORTHOGONALIZATION_CGS_REFINE_IFNEEDED`. Those solver files are unchanged from the starting commit. No solver compatibility workaround was introduced. Consequently an ex44 executable/MPI run and the complete unrelated FEMuS test suite are not claimed. Logs are included.

## Reproduction and evidence

The README gives the standalone CMake and broader suite commands. The comparison can be rebuilt from independent checkouts with:

```sh
python3 unittests/testCutFemWeightLSNative/run_endpoint_only_comparison.py --accepted-root /path/to/accepted --rifat-root /path/to/rifat --boost /path/to/boost
python3 unittests/testCutFemWeightLSNative/run_outer_interval_benchmark.py --baseline /path/to/baseline-benchmark --candidate build-cutfem/cutfem_outer_interval_benchmark --output unittests/testCutFemWeightLSNative/unified_outer_results
python3 unittests/testCutFemWeightLSNative/run_endpoint_analytic_convergence.py --binary /path/to/comparison_candidate
```

The last two reference checks need Python numpy and mpmath. See `endpoint_only_results/` for source manifests, the original suite summary, focused stdout and checks, CTest output, sanitizer results, 405-case reference results, separate volume/surface errors, tolerance and mesh tables, nested-contact results, all paired timing samples, and ex44 build evidence. Historical `review_results` files are not acceptance evidence for this commit.
