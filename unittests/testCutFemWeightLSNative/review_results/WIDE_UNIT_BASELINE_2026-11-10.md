# Wide unit-coordinate integration baseline — 2026-11-10

The corrected implementation passes the complete integration matrix below. The starting commit was `e6db9f8eefa5eb905fa984866c3ebacfe42a5c7e` on `rifat-cutfem-hex-tet-result`; the production source digest tested is `e8c92b89e49d94a0e4b7a060f1189b28026d818e3350155bf206db4ff18770be`. The `rifat` benchmark branch is unchanged.

## Corrections established by the wider tests

- Keep the simplex inverse-Gram product in oct precision until each final moment-to-weight matrix entry is formed. Casting the inverse first gave a TRI degree-10 plane surface moment of about `0.7499999762946` instead of `0.75`. The corrected plane errors are about `1e-14` or smaller. Reusing each Gram-matrix row product also removes repeated work in orthonormalization.
- Bound TRI polynomials on the triangle with Bernstein coefficients and split along the actual affine derivative-zero line, including affine factors of the simplex bubble derivatives. The enclosed TRI bubble previously exhausted the evaluation budget after over 1,500 splits; the corrected bubble takes 2–3 geometric splits in the focused tolerance checks.
- Apply the tangency endpoint mapping to full/multiple-root integration intervals as well as certified graph intervals, retaining the stable factored discriminant.
- Share the absolute error budget across strips and topology intervals. Higher-degree slab surface integrals previously overspent the accumulated budget. A circle crossing a coarse TRI mesh also produced a legitimate interval only about `1e-16` wide and an impractically small width-proportional budget; the corrected allocation keeps the original final tolerances.
- If a prescribed TET slice direction fails, retry with the locally certified graph directions. The regression is the plane with prescribed direction 4, whose outer slices are nearly parallel to the surface. Successful prescribed slices retain their original path.

## Completed checks

| Check | Coverage | Result |
| --- | --- | --- |
| New unit-coordinate matrix | 321 cases, 129,630 moment checks; QUAD/TRI/HEX/TET; q=2,5,10 | PASS |
| Double-output fitting | 16 degree-10 cases, 16,822 assertions; plane, curved graph, sphere, variable Jacobian | PASS |
| Established native matrix | 386 original + 8 supplemental cases, 40,036 assertions | PASS |
| Physical circle/sphere meshes | 54 configurations, 216 global volume/surface/x-squared checks | PASS |
| Focused APIs, basis reconstruction, roots and reuse | 19 executions, including 4,212 unit-element checks | PASS |
| Original planar-rule comparisons | 147 cases; 16,258 weight and 4,552 moment comparisons | PASS |
| 3D partition identities and analytic moments | 132 cases, 34,116 checks through degree 10 | PASS |
| Partition/endpoint/fallback/convergence review | 12,345 checks, including 144 circles and 48 translated ellipses | PASS |
| Enriched non-graph 3D fields | 24 cases, 6,080 checks, both signs and raw/surface modes | PASS |
| ASan + UBSan | 15 cases, 1,590 complete moment rows, through degree 10 | PASS |

The new matrix includes full/empty cells, boundary interfaces, axis and oblique planes, two sheets, enclosed/tiny/tangent balls, ellipses, and simplex bubbles. It checks raw volume, raw delta, reference surface, fitted weights, and nonlinear physical cofactors. Sign reversal, scaling by powers of two, requested directions, and tolerances from `1e-8` to `1e-12` are included. Physical-volume checks respect the polynomial degree needed by the variable determinant. No requested monomial degree exceeds 10.

Independent references use exact ball/simplex/plane formulas or independently parameterized angular quadrature at orders 80 and 144. Each moment is checked against `2*tol*(1+abs(reference)) + 8*reference_discrepancy + 2e-12`; the largest observed error divided by this complete limit is `0.016883`. Missing moments, integration exceptions, unresolved 3D graph cells, nonfinite results, and changed production headers reject the matrix. The complete per-case assessments are in [wide_unit_cases.json](wide_unit_cases.json).

Across the wide matrix the recorded geometric split totals are 1,045 QUAD splits, 8,485 TRI splits, and 8,996 3D cell splits; the TRI total includes the intentionally unsuitable prescribed TET slice test before its graph retry. Accepted regions recorded 204,630 QUAD and 466,560 TRI slice calls. These totals cover repeated raw/surface/physical calls and are not averages per element.

## Physical mesh accuracy

All degrees 2, 5, and 10 pass at every tested level. There are six levels per 2D element (2 through 64 subdivisions per coordinate) and three per 3D element (2 through 8). A translated circle/sphere is reconstructed from the actual nodal FE values; cut cells use the complete fitted physical-surface API. Errors are already close to rounding, so a monotone convergence rate is not inferred.

| Element | Configurations | Largest mesh | Maximum volume error | Maximum surface error |
| --- | ---: | ---: | ---: | ---: |
| QUAD | 18 | 4,096 cells | 3.508e-14 | 6.513e-14 |
| TRI | 18 | 8,192 cells | 2.020e-14 | 8.312e-14 |
| HEX | 9 | 512 cells | 2.164e-14 | 1.221e-13 |
| TET | 9 | 3,072 cells | 2.530e-13 | 8.117e-14 |

The maximum volume x-squared error is `1.345e-14`; the maximum surface x-squared error is `2.747e-14`. All 61,404 cell visits are represented in [wide_unit_mesh_summary.json](wide_unit_mesh_summary.json), with split counts and wall times.

## Absolute graph timings

Mean milliseconds per element for the same ten graph inputs, three interleaved rounds, pinned to one CPU. These measure positive-volume and reference-surface moment integration; they exclude weight-matrix setup and the physical geometry/cofactor path. The host is an AMD EPYC 9V74, with the recorded generic x86-64 O3 flags.

| q | QUAD ms | TRI ms | HEX ms | TET ms |
| ---: | ---: | ---: | ---: | ---: |
| 2 | 0.104724 | 0.125574 | 3.572934 | 7.597828 |
| 5 | 0.116500 | 0.137720 | 4.308785 | 9.824062 |
| 10 | 0.225294 | 0.168940 | 13.216808 | 17.346505 |

The robustness changes have a measurable TRI/TET overhead in this graph benchmark: against the pinned pre-unit implementation `e2a9c1056812fe568a54b412e39acd3671e07bdd`, TRI is about 15–23% slower and TET about 5–16% slower, depending on degree. QUAD/HEX retain their unit-basis gains. The comparison includes the change in tensor moment basis and tolerance measure, as recorded in [wide_unit_timing.json](wide_unit_timing.json).

The complete physical-surface path on coarse TET sphere meshes remains expensive: the 48-tetrahedron level took about 152–305 seconds depending on q. The 3,072-tetrahedron level took about 49–82 seconds. These mesh runs were concurrent, so their wall times are diagnostic rather than pinned performance comparisons. Accuracy passes; reducing this cost remains further work.

## Scope and reproduction

These are standalone integration tests linked with the actual FEMuS FE bases and Gaussian quadrature sources. The full MPI/PETSc ex44 application was not built or run here. ASan and UBSan completed with recovery disabled; leak detection was disabled because LeakSanitizer cannot execute in this ptrace environment. The CI job enables leak detection on its native runner.

```sh
python3 build.py --jobs 2 --eigen /path/to/eigen3 --boost /path/to/boost/include
python3 run.py --selection original
python3 run.py --selection supplemental
python3 run_focused.py --extended
python3 compare_planar.py
python3 run_wide_double_output.py
python3 run_wide_unit_validation.py
python3 run_wide_unit_mesh.py
python3 run_wide_unit_sanitizers.py
./build/graph_cell_splitting 2
./build/graph_cell_splitting 6
./build/graph_cell_splitting 10
python3 validate_enriched_3d.py
```

The dedicated `cutfem-unit-baseline.yml` workflow reruns the main matrix, fitted mesh checks, non-graph fields, existing regressions, and sanitizers on pushes to this branch. [wide_unit_acceptance.json](wide_unit_acceptance.json) records counts, source hashes, and scope for this accepted integration baseline.
