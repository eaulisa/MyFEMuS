# Curved ex44_test endpoint fix

Baseline: `7c8cbb894ef3314b0c8d68272de7d8d5d8967bef` (`testing curved elements`) on `rifat-cutfem-endpoint-only`. The checked-in runner has CURVED=1, CURVATURE=0.40, LEVELS=4, degree 8, moment mode, VTK enabled, and HEX/TET/QUAD/TRI.

## Reported HEX failure

The exact element-3 nodal values and coordinates from the report reproduce `Unsupported endpoint topology: An endpoint map requires a positive simple-fold discriminant`. The full curved application also exposed an earlier orientation with `LINE domain integration did not reach the requested tolerance`.

A fold nearly coincides with an element-face crossing. The two roots were refined using different rounded evaluations of the same slice polynomial. Cancellation in B*B-4*A*C produced inconsistent event positions and an artificial interval about 7.59e-19 wide; its mapped reduced discriminant was negative (about -4.34e-20). This is a numerical topology error, not a singular physical sphere.

`LineDomainIntegration.hpp` now refines both discriminant events and boundary crossings against the original slice polynomial using compensated long-double sums and products. This retains the small coefficient contributions from the nodal field and makes the two event calculations consistent. The fundamental LINE integration kernels and numerical quadrature routines are unchanged.

## Related near-boundary fold

The broader curved runner exposed a QUAD cell at level 2 whose fold lies just outside the cell. Its positive endpoint discriminant is extremely small, so bulk GK treats an almost singular density as a smooth endpoint and exhausts refinement.

`EndpointStripQuadrature.hpp` now also applies the retained quadratic map when the positive endpoint discriminant is small compared with its inward change across the strip. The endpoint constant is preserved in `D(a+h)=D(a)+h*R(h)`. Thus the transformed surface ratio is `2*delta/sqrt(D(a)/t^2 +/- delta*R(h))`, with h=+/-delta*t^2. The strip width, selective moment adaptivity, full polynomial, and actual-integrand quadrature are retained. No coefficient is refitted and no negative discriminant is clamped away.

## Full application runs

Both runs use the final compiled MultiPhase_ex44_test, one MPI rank, degree 8, default absolute and relative tolerances 1e-8, levels 0 through 4, all four geometries, and VTK. The curved run uses the checked-in run_ex44_test.sh unchanged. The straight regression uses identical options with --curved 0 in a separate output directory. Both exit 0. All 160 result rows are finite and meet the absolute-plus-relative analytical threshold.

Final curved level, constant integrand:

| Element | Elements | Cut elements | Volume/area absolute error | Surface/curve absolute error |
|---|---:|---:|---:|---:|
| HEX | 2136 | 296 | 7.77156117e-16 | 8.10462808e-15 |
| TET | 2339 | 248 | 2.67928041e-15 | 1.73697834e-10 |
| QUAD | 184 | 36 | 1.66533454e-16 | 1.73194792e-14 |
| TRI | 184 | 32 | 1.24900090e-16 | 1.55431223e-15 |

The eight CSV files in `curved_amr_results/` preserve every level, both integrands, exact references, and separate volume and surface errors for the curved and straight runs.

## Maintained tests

- 22/22 native CTests pass.
- New HEX regression: exact reported data, eight chart reflections, raw delta/reference surface/physical surface moments, total degree through 10, tolerances 1e-10 and 1e-12; 48 cases and 27,742 checks. Maximum reference-volume error 2.6319008e-17 and physical-surface error 3.701264e-17.
- New QUAD regression: exact failing application cell, four chart reflections, the same three metrics and two tolerances; 24 cases and 3,234 checks. Maximum reference-area error 1.3859865e-15 and physical-curve error 2.0618815e-16.
- References use analytical ball/disk monomial integrals and an independent angular parametrization of each curved interface, converged separately at Gauss orders 64 and 128.
- AddressSanitizer and UndefinedBehaviorSanitizer pass both new regression executables; LeakSanitizer is disabled because process inspection is unavailable in this environment.
- Existing LINE, graph, endpoint-contact, coordinate-resolution, and mesh tests remain active and pass.

Build environment and PETSc 3.22 compatibility compiler settings match `EX44_AMR_SPHERE_FIX.md`. The complete curved runner took 199.94 seconds and the straight regression took 161.11 seconds in this environment; these include diagnostic printing, VTK, and concurrent checks.

No geometric subdivision or competing integration path was added. Existing explicit failures for unsupported singular input or unmet error estimates remain active.
