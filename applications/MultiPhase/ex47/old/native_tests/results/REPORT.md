# Complete native comparison

Compared all 50 C++ sources in `unittests/testCutFemWeightLSNative` and one additional public-API curved HEX27 reproduction. Original fixtures, tolerances, and independent references are retained.

The 386 public-class cases pass 36,860 original assertions per backend after the certified-boundary-root fix. The broader fixtures still expose prototype regressions; see `failures.json` and the table below.

Timing ratios are original / prepared; values above one favor prepared. Public API batches clear final weights before every call and retain rule/element preparation. Constructor timing is separate. Benchmark drivers use their internal timers; other fixture timings include startup and setup. Failed cases have no accepted speedup. Primitive adapter timings include coordinate/format conversion; unchanged mathematical helpers and the legacy planar reference remain shared infrastructure.

## Public API timing

| Element | Raw/reference weights | Physical surface |
|---|---:|---:|
| LINE3 | 5.725× | 6.175× |
| QUAD9 | 1.071× | 1.405× |
| TRI7 | 1.535× | 2.773× |
| HEX27 | 1.785× | 2.111× |
| TET15 | 1.842× | 5.057× |

The exact CTest requests pass 32/32 originally and 28/32 in the prototype. One prototype difference is the requested shallower, conservative detection default; its depth-12 control passes. Shared failures in the extended parameterized runs are listed separately in `summary.json`.

The accuracy-only versions of the timing drivers pass 59,520 additional paired moment comparisons: graph degrees 2/10 and all five degree-13 volume/tangential/normal moment families. Their reconstruction work is excluded from timed drivers.

## Every native source

| Test | Cases | Original pass | Prepared pass | Median speed ratio | Timing |
|---|---:|---:|---:|---:|---|
| `benchmark_graph_dimensions` | 2 | 2 | 2 | 3.656× | sum of benchmark case medians |
| `benchmark_quad_tri` | 1 | 1 | 1 | 0.648× | sum of benchmark case medians |
| `benchmark_surface_reuse` | 1 | 1 | 1 | 2.747× | same/different element calls; setup excluded |
| `direct_parent_moments` | 1 | 1 | 1 | 0.706× | whole fixture, including startup/setup |
| `endpoint_comparison` | 70 | 70 | 70 | 1.621× | internal seconds per call |
| `endpoint_coordinate_resolution` | 1 | 1 | 1 | 0.635× | whole fixture, including startup/setup |
| `endpoint_exponents` | 1 | 1 | 1 | 0.627× | whole fixture, including startup/setup |
| `endpoint_line_comparison` | 2 | 2 | 2 | 41.166× | internal seconds per call |
| `endpoint_nested_contacts` | 2 | 2 | 1 | 1.671× | whole fixture, including startup/setup |
| `endpoint_surface_regressions` | 3 | 3 | 3 | 0.968× | whole fixture, including startup/setup |
| `enriched_3d` | 24 | 24 | 24 | 1.198× | whole fixture, including startup/setup |
| `hex_curved_endpoint_regressions` | 1 | 1 | 0 | — | whole fixture, including startup/setup |
| `hex_sphere_boundary_tangencies` | 1 | 1 | 1 | 0.719× | whole fixture, including startup/setup |
| `hex_tet_validation` | 60 | 58 | 58 | 1.110× | driver integration seconds |
| `line_integration` | 1 | 1 | 1 | 4.217× | whole fixture, including startup/setup |
| `line_mesh_refinement` | 1 | 1 | 1 | 1.030× | whole fixture, including startup/setup |
| `mesh_refinement` | 1 | 1 | 1 | 1.176× | whole fixture, including startup/setup |
| `native_suite` | 386 | 386 | 386 | 2.364× | median per-case API speedup |
| `nodal_cut_detection` | 2 | 2 | 1 | 0.981× | whole fixture, including startup/setup |
| `normal_surface` | 1 | 1 | 1 | 4.753× | public API calls; construction/oracles excluded |
| `outer_interval_benchmark` | 8 | 8 | 6 | 2.493× | sum of internal case medians |
| `outer_interval_contract` | 1 | 1 | 1 | 0.892× | whole fixture, including startup/setup |
| `planar_suite` | 147 | 147 | 147 | 0.894× | whole fixture, including startup/setup |
| `quad_curved_endpoint_regressions` | 1 | 1 | 1 | 0.963× | whole fixture, including startup/setup |
| `quad_native_quadrants` | 1 | 1 | 1 | 1.020× | whole fixture, including startup/setup |
| `quad_scaled_plane` | 1 | 1 | 1 | 1.173× | whole fixture, including startup/setup |
| `quad_tri_standard_moments` | 1 | 1 | 1 | 0.743× | whole fixture, including startup/setup |
| `quadrature_rule_selection` | 1 | 1 | 0 | — | public API calls; construction/oracles excluded |
| `review_regressions` | 2 | 2 | 2 | 1.009× | whole fixture, including startup/setup |
| `simplex_high_degree` | 1 | 1 | 1 | 0.764× | whole fixture, including startup/setup |
| `single_graph_root` | 2 | 2 | 2 | 0.932× | whole fixture, including startup/setup |
| `surface_rule_reuse` | 1 | 1 | 1 | 10.896× | public API calls; construction/oracles excluded |
| `tangential_surface` | 1 | 1 | 1 | 3.384× | public API calls; construction/oracles excluded |
| `test_contract` | 1 | 1 | 1 | 2.198× | public API calls; construction/oracles excluded |
| `test_hex_native` | 2 | 2 | 2 | 1.028× | whole fixture, including startup/setup |
| `test_ls_multiprecision` | 1 | 1 | 1 | 0.492× | public API calls; construction/oracles excluded |
| `test_ls_output_range` | 1 | 1 | 1 | 0.966× | public API calls; construction/oracles excluded |
| `test_order_ls_first` | 1 | 1 | 1 | 0.636× | public API calls; construction/oracles excluded |
| `test_order_planar_first` | 1 | 1 | 1 | 0.516× | public API calls; construction/oracles excluded |
| `test_reference_coefficients` | 1 | 1 | 1 | 0.974× | whole fixture, including startup/setup |
| `test_root_candidate_regression` | 1 | 1 | 1 | 0.688× | whole fixture, including startup/setup |
| `test_tet_slice_regression` | 1 | 1 | 1 | 0.913× | whole fixture, including startup/setup |
| `tet_standard_moments` | 1 | 1 | 1 | 1.380× | whole fixture, including startup/setup |
| `tolerance_warnings` | 1 | 1 | 0 | — | public API calls; construction/oracles excluded |
| `tri_circle_probe` | 1 | 1 | 1 | 1.040× | whole fixture, including startup/setup |
| `tri_endpoint_regressions` | 1 | 1 | 1 | 0.934× | whole fixture, including startup/setup |
| `unit_element_moments` | 1 | 1 | 1 | 1.098× | public API calls; construction/oracles excluded |
| `wedge21_moments` | 1 | 1 | 1 | 1.642× | public API calls; construction/oracles excluded |
| `wide_unit_mesh` | 12 | 12 | 12 | 2.069× | mesh seconds; rule setup excluded |
| `wide_unit_validation` | 12 | 11 | 8 | 1.763× | case seconds; fitting setup excluded |

## Work behind the slow cases

41 of 772 public-case/API-family rows are over 10% slower; 20 have separated three-batch timing ranges. The largest low-order QUAD boundary regressions use the same 15 samples in both implementations. The prepared path additionally creates/restricts 15 LINE records and performs root analysis even though the boundary root is certified. This adds metadata and topology work to an otherwise inexpensive analytic-root path. The fixed quadrature and fitting matrix are reused.

| Case | API | Original µs | Prepared µs | Original samples | Prepared samples | LINE restrictions |
|---|---|---:|---:|---:|---:|---:|
| `e1_q0_boundary_low_identity_s1_v0` | weights | 15.20 | 32.36 | 15 | 15 | 15 |
| `e1_q0_boundary_high_identity_s1_v0` | weights | 14.31 | 29.49 | 15 | 15 | 15 |
| `e1_q0_boundary_low_identity_s1_v0` | physical | 14.78 | 29.01 | 15 | 15 | 15 |
| `e1_q0_boundary_high_identity_s1_v0` | physical | 15.76 | 29.63 | 15 | 15 | 15 |
| `e1_q2_boundary_low_identity_s1_v0` | weights | 20.33 | 35.26 | 15 | 15 | 15 |
| `e2_q1_diagonal_identity_s1_v0` | weights | 33.47 | 55.15 | 15 | 15 | 16 |
| `e1_q2_boundary_low_identity_s1_v0` | physical | 20.61 | 33.21 | 15 | 15 | 15 |
| `e1_q2_boundary_high_identity_s1_v0` | weights | 23.06 | 36.01 | 15 | 15 | 15 |
| `e1_q2_boundary_high_identity_s1_v0` | physical | 22.58 | 34.51 | 15 | 15 | 15 |
| `e1_q5_zero_face_identity_s1_v0` | weights | 46.77 | 69.10 | 15 | 15 | 15 |

## Two retained elements and all surface families

Five-batch medians from the existing ex47 benchmark; final weights are cleared and the two element records alternate. Both checksum and accuracy checks remain active. The table gives median speed ratios across degrees for each family.

| Element | Raw delta | Reference | Physical | Tangential | Three normals |
|---|---:|---:|---:|---:|---:|
| LINE3 | 9.424× | 8.549× | 6.176× | 7.868× | — |
| QUAD9 | 0.910× | 1.274× | 1.850× | 1.871× | — |
| TRI7 | 1.053× | 1.435× | 1.931× | 1.981× | — |
| HEX27 | 0.950× | 1.269× | 3.135× | 3.094× | 3.827× |
| TET15 | 1.009× | 1.999× | 4.097× | 4.584× | 5.041× |
| WEDGE21 | 1.050× | 2.008× | 4.098× | 4.081× | 5.034× |

## Files

`native_cases.csv` gives each public fixture and API family, independent error, median/range, and work counts. `construction.csv` separates constructor cost. `fixtures.csv` gives every executable invocation; `ctest.csv` identifies the exact CTest requests. `benchmark_cases.csv` expands driver subcases, and `retained_modes.csv` records all six elements and requested surface families. `regressions.csv` retains slower rows and timing-range checks. `sources.csv` gives the source-level table. `failures.json` retains failed assertions, shared limitations, and prototype-only failures. `summary.json` records source/binary checksums and coverage.

See `../DIAGNOSIS.md` for failure analysis and `../README.md` for reproduction commands.
