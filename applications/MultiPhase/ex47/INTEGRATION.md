# Native prepared weight integration

`src/quadrature/cutfem_ls/CutFemWeight_LS.hpp` provides the native global class. The implementation was moved from ex47; the former native code is archived under `src/quadrature/cutfem_ls/old` and is not included or built.

```cpp
CutFemWeight_LS<long double,long double> weights(HEX,5,"legendre");
CutFemLevelSetData element(HEX,nodalValues,1);
std::vector<long double> positive;
std::vector<long double> surface;
weights.GetPhysicalSurfaceWeightsFromLevelSetData(element,coordinates,positive,surface,1e-10,1e-10);
```

Construct the weight object once for a fixed element type, fitting degree, and Gauss family. Retain an element data object across classification and integration. Geometry is a separate vector indexed by physical component, then FEM node. Only physical, tangential, and signed-normal requests need it.

| Request | Prepared-data interface | Nodal-input interface |
| --- | --- | --- |
| Positive volume and delta | `GetWeightsFromLevelSetData` | `GetWeightsFromNodalValues` |
| Positive volume and reference surface | Same method with `referenceSurface=true` | Same flag on nodal method |
| Physical surface | `GetPhysicalSurfaceWeightsFromLevelSetData` | `GetPhysicalSurfaceWeightsFromNodalGeometry` |
| Inverse-metric tangential measure | `GetTangentialSurfaceWeightsFromLevelSetData` | `GetTangentialSurfaceWeightsFromNodalGeometry` |
| Three signed normal measures, 3D only | `GetNormalSurfaceWeightsFromLevelSetData` | `GetNormalSurfaceWeightsFromNodalGeometry` |
| Positive volume, tangential and three signed normals, 3D | `GetTangentialAndNormalSurfaceWeightsFromLevelSetData` | `GetTangentialAndNormalSurfaceWeightsFromNodalGeometry` |
| Classification | `element.GetPhase()` | `GetLevelSetPhaseFromNodalValues` |

Prepared overloads require a nodal root record matching the weight object's element. A restricted child retains its ancestor's derivative context and is used internally by the recursive integration. Nodal wrappers retain an owned record, exposed by `GetPreparedLevelSetData()`. Snapshot changes invalidate the relevant preparation.

`SetIntegrationDirection`, `SetDirectionValidation`, `SetOuterQuadratureRules`, `ClearCache`, Gauss accessors, and external-moment fitting retain the ex47 interface. `SetMaximumCertificationDepth` controls the owned nodal record; an external record uses its own depth. Direction numbering includes the TRI/TET barycentric option. Validation mode 2 compares alternative prepared projections. GK rules are shared immutable objects selected before integration. Applications select the rules on the weight object before calling the prepared-data overloads. The former native physical-surface overload's extra Boolean argument is removed.

The implementation uses prepared coefficient, graph, boundary, projection, LINE root/interval, gradient, and restriction methods. Integration drivers in the native directory retain the ex47 endpoint maps and adaptive moment policy. Unit coordinates are used throughout coefficient preparation and recursive integration. Tensor polynomials are reconstructed directly on `[0,1]`; simplex coordinates retain their standard unit domains. Final moments use the unit-monomial fitting basis. The public boundary applies the reference/physical density factors. The moment-to-weight matrix and Gauss points are unchanged and are built only in the weight constructor.

Physical cofactor-gradient polynomials are prepared once per matching field/geometry snapshot, trimmed only where high powers are exactly zero, and restricted without lowering their degree. This cache survives rule, tolerance, surface-mode, and final-weight cache changes. Delta/reference requests do not prepare geometry. Reference requests preserve all original gradient components through slices; prepared physical requests restrict only their external cofactor-gradient components. Signed normal components retain their original orientation in barycentric charts.

The combined 3D call shares the restricted field, physical cofactor-gradient and its norm at common quadrature nodes across the tangential and signed normal weights. Each density selects its projection and refines its own unresolved moments. Shared sample storage lasts for one element call.

Known endpoint roots use the parameterized map's slope when the absolute coordinate has rounded to the endpoint. `GetRegularRootDeltaFactor` permits density evaluation directly from this regular root record. Other LINE evaluations reuse the complete prepared topology and stable signed slopes. Sensitive ancestor fibers are contracted directly in compensated algebra, skipping zero polynomial terms, rather than recovered from three full ancestor evaluations. This preserves the same polynomial and introduces no approximation. Slice classification is lazy, so LINE root-only work does not pay for unused Bernstein sign certification.

## Quadrature selection on HEX27, TET15, and WEDGE21

Let `t` denote the outer slicing coordinate and `u,v` the transverse coordinates. A vertical tangent inside a slice and a saddle of the outer height are different events. They can occur at different levels of the same integral.

```text
integrate(field, requested moments and weight)
    Classify the element and account for any face coincident with the interface.
    If the phase is certified and the interface contributes nothing, use exact full or empty moments.
    Look for a transverse direction with a derivative bounded strictly away from zero.
    If one exists, use it as the inner graph direction and order the remaining outer axes by derivative bound.
    Otherwise, try outer axes in order of increasing maximum absolute derivative bound.
        For each candidate t, find its regular interior critical points once for this field and chart.
        Prefer an axis without a relevant saddle. Honor an explicitly requested outer axis.
        If every candidate has a relevant saddle, use the first in the original order.
    If the field and density are independent of an outer product coordinate, integrate its powers exactly.
    Otherwise, get the event heights for the chosen outer axis and inner chart.
    For each interval between consecutive event heights, choose the panel rule from the table below.
        At each outer node, restrict the original field and requested density to the actual QUAD or TRI slice.
        On that slice, use a certified graph when available. Otherwise locate folds and boundary intersections.
        On each LINE, integrate positive intervals analytically and sum the density at each simple root.
    Accumulate the outer rule difference and the weighted inner error for each requested moment.
    Refine the one-dimensional panels or retry eligible directions for unresolved moments within budget.
    Return the fitted weights and report any unmet tolerance.
```

| Event at a panel endpoint | Coordinate and quadrature | Reason |
| --- | --- | --- |
| Regular interval | Ordinary configured Gauss–Kronrod pair, default 7/15 | The slice density varies smoothly. |
| Ordinary vertical tangent or fold | Use `t = t* ± L s²`, include the Jacobian, and use the configured tangent Gauss–Kronrod pair, default 7/15 | The square substitution removes the square-root behavior. For a higher supported projected contact the code uses the power determined by its multiplicity. |
| Regular interior saddle of the outer height | Use normalized distance `s = abs(t − t*) / L` and the generalized 8- and 12-point Gaussian rules on the adjacent panel | The slice density can contain `A(s) + B(s) log(s)`. The N-point rule integrates both `s^j` and `s^j log(s)` for `0 ≤ j < N`. The two rules give the local error estimate. |
| Regular interior height minimum or maximum | Insert the critical height as an interval boundary and use the ordinary or power-mapped rule selected by the contact test | The transverse Hessian is definite. A logarithmic rule is specific to an indefinite Hessian. |
| Singularities at both ends | Bisect the one-dimensional interval and apply each endpoint's rule on its own half | Each map or logarithmic rule is centered at one endpoint. With one singular endpoint its panel spans the interval. |

For a regular interior critical point, the code solves `G = G_u = G_v = 0` and checks `G_t ≠ 0`. The sign of `G_uu G_vv − G_uv²` distinguishes a saddle (negative) from an extremum (positive). The candidate point must lie inside the reference cube, the tetrahedron `x + y + z < 1`, or the prism `x + y < 1, z < 1`, respectively, with positive coordinates. The derivative bound orders the candidate axes cheaply. It cannot distinguish the two Hessian signs. For a nearly bilinear transverse field, the bilinear stationary equations supply extra Newton seeds. Every accepted point is checked on the full interpolant. A singular level set or degenerate Hessian requires separate analysis.

A scalar weight demotes a direction with any detected saddle. For a signed normal weight the ranking uses the physical cofactor-gradient at each saddle and demotes the direction when the requested component is dominant, including ties. Once an axis is selected, its saddle heights receive the logarithmic rule. After adaptive bisection only the child touching the saddle uses that rule. The other child uses the configured tangent Gauss–Kronrod pair in the inherited distance coordinate. All subdivisions are intervals of quadrature on the original element.

If child accuracy or coordinate resolution limits a panel, the integrator can try lower Gauss–Kronrod pairs. The final error check includes the error returned by the slices.

Critical points and projection resultants are shared across weights and retries for the same reference element, normalized level-set polynomial, and outer axis. The endpoint list is also shared for a matching inner direction and graph/general mode. The dominant physical normal components are evaluated once for the current geometry and chart and recomputed when the metric changes. Each density retains its own moment integration and error estimates.

## Tests and timing

All executable checks remain in `ex47.cpp` main. Run `--test-weights --weight-dimension 1`, then 2, then 3, or omit the dimension to run those stages in order. `--benchmark-weights` uses the same selector. The original data tests remain the default main path.

Accuracy checks compare fitted weights, raw prepared moments, and independent references: 100-digit LINE interpolation/integration, analytic ball volume moments, and separate angular surface integrals for an explicit non-affine geometry. Nodal and prepared interfaces are also compared. Fitted orders 2 and 5 cover every element; order 10 also covers LINE, QUAD, TRI, and HEX. Raw order-13 moments cover every 3D element without repeatedly constructing the expensive high-order simplex fitting matrix. Tests also exercise all direction numbers, enriched nodal fields, alternate GK rules, geometry invalidation, signed normals, and tiny circles.

The supplied Gauss tables use finite decimal precision. The tests therefore also fit independent exact moments through the unchanged production fitting operator and measure its residual. Raw integration error and this fitting floor are reported separately. A fitted residual must be judged against that measured floor.

Timing constructs the Gauss/fitting objects outside the measured loop and clears final fitted weights on every iteration. Two distinct prepared element records alternate, while reusable workspace storage remains available. Geometry-dependent metrics are rebuilt when switching fields and reused within a requested normal family. Five interleaved runs report medians. These measure repeated integration and fitting rather than final-weight cache lookup. The benchmark reports `nodal_us`, `prepared_us`, and `nodal_over_prepared` for the two input interfaces of the same native implementation. Runtime ratios are machine- and field-dependent; the benchmark prints both absolute times.

The unit-coordinate adapter analyzes local fibers directly and evaluates only the requested surface density. Topology, discriminants, endpoint slopes, and moment coordinates use the same current-slice polynomial. Ancestor-aware algebra remains available separately in the data class. No inner quadrature sample constructs a LINE object.
