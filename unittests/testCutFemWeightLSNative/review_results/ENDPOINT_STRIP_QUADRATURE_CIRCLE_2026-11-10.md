# Numerical quadrature on regularized endpoint strips

Circle experiment, 10 November 2026.

This report describes the circle-only snapshot at commit aea2c13efd1ac883103461d3218790cfbe70e4f0; use that commit to reproduce its original timing controls. The later general QUAD/TRI comparison is recorded in ENDPOINT_STRIP_QUADRATURE_2D_2026-11-10.md.

The prototype evaluates the original FE level-set polynomial and applies GK quadrature after the substitution x = a + delta*u^2 at a left fold, or x = b - delta*u^2 at a right fold. It keeps the original integration direction in the strips and the bulk. It does not approximate the curve.

For the surface integrand, factor D(x) = (x-a)*Q(x-a). The mapped Jacobian divided by the root derivative is evaluated as 2*delta/sqrt(delta*Q(delta*u^2)); the factor u cancels algebraically. This avoids subtracting nearly equal discriminant terms at the endpoint and avoids a numerical 0/0. The right-end expression includes the corresponding orientation sign.

The event list, topology classification, and direction handling reuse the published 5840de91e9d3d6564b562975ece93e010e6c9550 experiment. Accepted moment components are frozen; GK15 is tried first, then GK31 on the remaining components, then bisection. All moment components share the roots at each node.

## Scope

This is a standalone circle prototype on branch rifat-cutfem-2d-endpoint-quadrature, based on baseline 1a394fdc97bedc5a6cff9b1e6c68f16261ba57e9. Existing production graph, QUAD, TRI, HEX, and TET kernels are unchanged. The new header is not dispatched by the production integrator. General non-graph geometry, curved physical mappings, fitted weights, and ex44 have not been validated for this prototype.

## Accuracy and timing

G = 0.15^2 - (x-0.29)^2 - (y-0.23)^2. The circle lies entirely inside both the unit square and the reference triangle. Degree q is the maximum total monomial degree. Requested reference accuracy is absolute 1e-10 plus relative 1e-10. The two strip prototypes use absolute and relative 1e-12 in the unit basis; baseline uses 1e-10. This is the same tightened bulk budget used in the earlier small-strip diagnostic.

All 39 cases pass at q = 2, 5, 10: 30 mapped cases across five widths, six baseline cases, and three switched-direction QUAD comparisons. There are 2,418 volume/surface moment checks against 80-digit analytic circle/disk formulas. QUAD values are also transformed to the signed FEMuS reference basis, including factors four for volume and two for surface. No degree above ten was tested.

Times include nodal polynomial reconstruction, metrics, event construction, and volume plus surface moments. They exclude weight fitting. They are medians of nine batches of 16 calls, after two warm-up batches. Nodal constants alternate by 2.910383045673370361328125e-11, forcing geometry reconstruction instead of timing cached topology. All comparisons below were measured serially in the same session.

| QUAD, q=10: strip width | Median time (ms) | Strip samples | Bulk samples | Strip bisections | Bulk bisections |
|---:|---:|---:|---:|---:|---:|
| 0.015 | 5.338 | 30 | 260 | 0 | 3 |
| 0.03 | 3.336 | 30 | 138 | 0 | 1 |
| 0.043933983 | 2.088 | 30 | 46 | 0 | 0 |
| 0.05 | 2.047 | 30 | 46 | 0 | 0 |
| 0.075 | 3.030 | 92 | 46 | 0 | 0 |

| q=10 comparison | QUAD (ms) | TRI (ms) |
|---|---:|---:|
| Existing integration | 4.368 | 3.094 |
| Direction switch inside strips | 3.336 | not tested |
| Square-root map, delta=0.05 | 2.047 | 1.764 |

At delta=0.05, QUAD maximum signed-reference errors are 8.692e-15 for volume and 7.504e-17 for surface. Each endpoint strip accepts GK15 (15 nodes), while the middle accepts the GK15/31 sequence (46 nodes). Total: 76 samples, zero bisections. The bulk interval is 0.2 wide. No model-error test, strip shrinking, Taylor series, or fallback is used.

A preliminary narrow-strip control at delta=1e-8 required 4,124 bulk samples and 45 bulk bisections. It missed the reference volume tolerance slightly (maximum absolute error about 1.22e-10), so it is not included among the passing cases above. This prototype still budgets errors in the unit basis; selecting a general reference-basis error budget remains future work.

## Verification

AddressSanitizer and UndefinedBehaviorSanitizer pass for degree-ten QUAD at delta approximately 0.043934 and 0.075 and TRI at delta=0.05. LeakSanitizer is unavailable because this environment runs the process under ptrace.

Implementation and harness SHA-256 hashes are recorded in endpoint_strip_circle.json; sanitizer results are in endpoint_strip_circle_sanitizers.json.

## Reproduce

From unittests/testCutFemWeightLSNative with Boost, Eigen, a C++14 compiler, and Python mpmath installed:

```sh
python3 build.py --tests endpoint_strip_circle --jobs 2
python3 run_endpoint_strip_circle.py
```

Use build.py --boost and --eigen when headers are not in the usual system locations.
