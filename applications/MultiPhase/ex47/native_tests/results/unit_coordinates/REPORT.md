# Unit-coordinate recovery comparison

All runtime changes remain local to ex47. Tensor cells use [0,1], simplices retain their standard coordinates, and WEDGE uses a unit-height triangle product. FEMuS node ordering and the original quadrature strategy are retained.

All 50 native C++ sources were compiled and compared, alongside the curved-HEX reproduction. The 386 public cases pass 36,860 assertions per backend. The additional graph and degree-13 surface drivers pass 59,520 paired moment checks.

Main data, fitted-weight, curved-HEX, warning, and sanitizer validation details are recorded in validation.json. Leak detection was disabled for the address/undefined-behavior sanitizer run.

Both implementations pass all 32 CTest commands at matched detection depth 12. The requested prepared default depth 1 conservatively classifies one additional uncertain element as cut; that default-depth fixture difference is retained separately. Extended fixture failures are listed below without exclusion.

## Public API timing

Ratios are original / prepared; above one favors prepared. Each case uses three interleaved batches pinned to one CPU. Every timed result is checked. Final weights are cleared before each call; rule and element preparation are retained. Other fixture timings are single runs and include driver/adapter overhead, so use the repeated public API results for performance conclusions.

| Element | Raw/reference weights | Physical surface |
|---|---:|---:|
| LINE3 | 1.583× | 4.830× |
| QUAD9 | 1.004× | 1.312× |
| TRI7 | 1.138× | 2.187× |
| HEX27 | 2.161× | 2.154× |
| TET15 | 1.967× | 3.087× |

53 of 772 public API rows are over 10% slower; 17 have separated batch ranges. Performance is not uniformly better.

The inner planar quadrature constructs no LINE records and restricts no metric per sample. Slower inexpensive cases retain prepared fiber bookkeeping, coefficient-resolution checks, endpoint metadata, and complete interval analysis where the original graph path can use a single certified root. LINE calls also retain prepared-data dispatch and error checks. Identical evaluation counts in the rows below identify algebra/metadata overhead rather than extra adaptive subdivision. This is a source-level diagnosis supported by work counters, rather than an instruction-level profile. These costs remain follow-up optimization work; quadrature rules are unchanged.

| Case | API | Original µs | Prepared µs | Original samples | Prepared samples |
|---|---|---:|---:|---:|---:|
| `e0_q0_plane_identity_s1_v0` | weights | 1.86 | 3.33 | 0 | 0 |
| `e0_q2_empty_identity_s1_v0` | weights | 1.06 | 1.77 | 0 | 0 |
| `e1_q0_monotone_identity_s1_v0` | weights | 29.49 | 49.34 | 15 | 15 |
| `e0_q2_full_identity_s1_v0` | weights | 0.99 | 1.61 | 0 | 0 |
| `e0_q0_empty_identity_s1_v0` | weights | 0.80 | 1.29 | 0 | 0 |
| `e0_q0_full_identity_s1_v0` | weights | 0.79 | 1.27 | 0 | 0 |
| `e2_q4_diagonal_identity_s1_v0` | weights | 41.33 | 65.22 | 15 | 15 |
| `e0_q5_empty_identity_s1_v0` | weights | 1.34 | 2.12 | 0 | 0 |
| `e1_q1_monotone_identity_s1_v0` | weights | 30.72 | 46.51 | 15 | 15 |
| `e1_q2_boundary_low_identity_s1_v0` | physical | 20.72 | 31.26 | 15 | 15 |
| `e0_q0_diagonal_identity_s1_v0` | weights | 2.70 | 4.01 | 0 | 0 |
| `e0_q0_boundary_high_identity_s1_v0` | weights | 1.46 | 2.14 | 0 | 0 |

## Retained elements and selected surface measures

Five-batch medians, alternating two retained elements and clearing final weights. Ratios are medians over the tested degrees; signed normals comprise three component calls.

| Element | Raw delta | Reference | Physical | Tangential | Three normals |
|---|---:|---:|---:|---:|---:|
| LINE3 | 7.015× | 6.505× | 6.184× | 6.106× | — |
| QUAD9 | 0.887× | 1.368× | 1.961× | 1.930× | — |
| TRI7 | 0.951× | 1.340× | 1.799× | 1.924× | — |
| HEX27 | 0.951× | 1.550× | 3.763× | 3.633× | 4.453× |
| TET15 | 1.025× | 2.117× | 4.533× | 4.789× | 5.428× |
| WEDGE21 | 1.068× | 1.937× | 4.097× | 4.054× | 4.685× |

## Every native source

| Test | Cases | Original pass | Prepared pass |
|---|---:|---:|---:|
| `benchmark_graph_dimensions` | 2 | 2 | 2 |
| `benchmark_quad_tri` | 1 | 1 | 1 |
| `benchmark_surface_reuse` | 1 | 1 | 1 |
| `direct_parent_moments` | 1 | 1 | 1 |
| `endpoint_comparison` | 70 | 70 | 70 |
| `endpoint_coordinate_resolution` | 1 | 1 | 1 |
| `endpoint_exponents` | 1 | 1 | 1 |
| `endpoint_line_comparison` | 2 | 2 | 2 |
| `endpoint_nested_contacts` | 2 | 2 | 2 |
| `endpoint_surface_regressions` | 3 | 3 | 3 |
| `enriched_3d` | 24 | 24 | 24 |
| `hex_curved_endpoint_regressions` | 1 | 1 | 1 |
| `hex_sphere_boundary_tangencies` | 1 | 1 | 1 |
| `hex_tet_validation` | 60 | 58 | 60 |
| `line_integration` | 1 | 1 | 1 |
| `line_mesh_refinement` | 1 | 1 | 1 |
| `mesh_refinement` | 1 | 1 | 1 |
| `native_suite` | 386 | 386 | 386 |
| `nodal_cut_detection` | 2 | 2 | 1 |
| `normal_surface` | 1 | 1 | 1 |
| `outer_interval_benchmark` | 8 | 8 | 8 |
| `outer_interval_contract` | 1 | 1 | 1 |
| `planar_suite` | 147 | 147 | 147 |
| `quad_curved_endpoint_regressions` | 1 | 1 | 1 |
| `quad_native_quadrants` | 1 | 1 | 1 |
| `quad_scaled_plane` | 1 | 1 | 1 |
| `quad_tri_standard_moments` | 1 | 1 | 1 |
| `quadrature_rule_selection` | 1 | 1 | 1 |
| `review_regressions` | 2 | 2 | 2 |
| `simplex_high_degree` | 1 | 1 | 1 |
| `single_graph_root` | 2 | 2 | 2 |
| `surface_rule_reuse` | 1 | 1 | 1 |
| `tangential_surface` | 1 | 1 | 1 |
| `test_contract` | 1 | 1 | 1 |
| `test_hex_native` | 2 | 2 | 2 |
| `test_ls_multiprecision` | 1 | 1 | 1 |
| `test_ls_output_range` | 1 | 1 | 1 |
| `test_order_ls_first` | 1 | 1 | 1 |
| `test_order_planar_first` | 1 | 1 | 1 |
| `test_reference_coefficients` | 1 | 1 | 1 |
| `test_root_candidate_regression` | 1 | 1 | 1 |
| `test_tet_slice_regression` | 1 | 1 | 1 |
| `tet_standard_moments` | 1 | 1 | 1 |
| `tolerance_warnings` | 1 | 1 | 1 |
| `tri_circle_probe` | 1 | 1 | 1 |
| `tri_endpoint_regressions` | 1 | 1 | 1 |
| `unit_element_moments` | 1 | 1 | 1 |
| `wedge21_moments` | 1 | 1 | 1 |
| `wide_unit_mesh` | 12 | 12 | 12 |
| `wide_unit_validation` | 12 | 11 | 11 |

## Extended fixture exceptions

| Case | Original | Prepared |
|---|---|---|
| `hex_tet_validation_mesh_HEX_10_4_1e-11` | FAIL | PASS |
| `hex_tet_validation_mesh_HEX_6_4_1e-11` | FAIL | PASS |
| `nodal_cut_detection` | PASS | FAIL |
| `wide_unit_validation_1_2` | FAIL | FAIL |

The depth-1 detection difference is intentional; depth 12 passes. Shared nearly-parallel TRI limitations and strict HEX mesh differences remain visible in failures.json. No fixture assertion or public integration tolerance was relaxed. One last-ULP internal WEDGE discriminant check changed from 3e-16 to 4e-16 after remapping; the observed relative error is 3.0385e-16 and retained-slope density error is 1.4626e-16.

See native_cases.csv for all public errors, medians, ranges, and work counts; retained_modes.csv for all six elements and surface families; fixtures.csv and failures.json for every extended request; ctest.csv for matched-depth commands; summary.json for coverage and source/binary checksums. The runtime coordinate audit is ../../../UNIT_COORDINATE_AUDIT.md.
