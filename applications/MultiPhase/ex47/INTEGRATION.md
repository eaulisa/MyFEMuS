# Local prepared weight integration

`include/CutFemWeight_LS.hpp` duplicates the production interface inside `namespace Ex47`. All new runtime code remains in ex47; including the local header also makes the original global class available. No production assembly is switched over.

```cpp
Ex47::CutFemWeight_LS<long double,long double> weights(HEX,5,"legendre");
CutFemLevelSetData element(HEX,nodalValues,1);
std::vector<long double> positive;
std::vector<long double> surface;
weights.GetPhysicalSurfaceWeightsFromLevelSetData(element,coordinates,positive,surface,1e-10,1e-10);
```

Construct the weight object once for a fixed element type, fitting degree, and Gauss family. Retain an element data object across classification and integration. Geometry is a separate vector indexed by physical component, then FEM node. Only physical, tangential, and signed-normal requests need it.

| Request | Prepared-data interface | Nodal compatibility interface |
| --- | --- | --- |
| Positive volume and delta | `GetWeightsFromLevelSetData` | `GetWeightsFromNodalValues` |
| Positive volume and reference surface | Same method with `referenceSurface=true` | Same flag on nodal method |
| Physical surface | `GetPhysicalSurfaceWeightsFromLevelSetData` | `GetPhysicalSurfaceWeightsFromNodalGeometry` |
| Inverse-metric tangential measure | `GetTangentialSurfaceWeightsFromLevelSetData` | `GetTangentialSurfaceWeightsFromNodalGeometry` |
| Three signed normal measures, 3D only | `GetNormalSurfaceWeightsFromLevelSetData` | `GetNormalSurfaceWeightsFromNodalGeometry` |
| Classification | `element.GetPhase()` | `GetLevelSetPhaseFromNodalValues` |

Prepared overloads require a nodal root record matching the weight object's element. A restricted child retains its ancestor's derivative context and is used internally by the recursive integration. Nodal wrappers retain an owned record, exposed by `GetPreparedLevelSetData()`. Snapshot changes invalidate the relevant preparation.

`SetIntegrationDirection`, `SetDirectionValidation`, `SetOuterQuadratureRules`, `ClearCache`, Gauss accessors, and external-moment fitting retain the original interface. `SetMaximumCertificationDepth` controls the owned nodal record; an external record uses its own depth. Direction numbering remains the production numbering, including the TRI/TET barycentric option. Validation mode 2 compares alternative prepared projections. GK rules are shared immutable objects selected before integration.

The duplicate uses prepared coefficient, graph, boundary, projection, LINE root/interval, gradient, and restriction methods. Integration drivers under `include/prepared` retain the existing endpoint maps and adaptive moment policy. Unit coordinates are used throughout coefficient preparation and recursive integration. Tensor polynomials are reconstructed directly on `[0,1]`; simplex coordinates retain their standard unit domains. Final moments already use the original unit-monomial fitting basis. The public boundary applies the same reference/physical density factors as production. The original moment-to-weight matrix and Gauss points are unchanged and are built only in the weight constructor.

Physical cofactor-gradient polynomials are prepared once per matching field/geometry snapshot, trimmed only where high powers are exactly zero, and restricted without lowering their degree. This cache survives rule, tolerance, surface-mode, and final-weight cache changes. Delta/reference requests do not prepare geometry. Reference requests preserve all original gradient components through slices; prepared physical requests restrict only their external cofactor-gradient components. Signed normal components retain their original orientation in barycentric charts.

Only the requested scalar density is evaluated. The three-normal public wrapper currently performs three selected integrations, reusing its positive weights and geometry; a shared traversal of the whole normal family remains a separate production optimization.

Known endpoint roots use the parameterized map's slope when the absolute coordinate has rounded to the endpoint. `GetRegularRootDeltaFactor` permits density evaluation directly from this regular root record. Other LINE evaluations reuse the complete prepared topology and stable signed slopes. Sensitive ancestor fibers are contracted directly in compensated algebra, skipping zero polynomial terms, rather than recovered from three full ancestor evaluations. This preserves the same polynomial and introduces no approximation. Slice classification is lazy, so LINE root-only work does not pay for unused Bernstein sign certification.

## Tests and timing

All executable checks remain in `ex47.cpp` main. Run `--test-weights --weight-dimension 1`, then 2, then 3, or omit the dimension to run those stages in order. `--benchmark-weights` uses the same selector. The original data tests remain the default main path.

Accuracy checks compare both fitted classes, raw prepared moments, and independent references: 100-digit LINE interpolation/integration, analytic ball volume moments, and separate angular surface integrals for an explicit non-affine geometry. Fitted orders 2 and 5 cover every element; order 10 also covers LINE, QUAD, TRI, and HEX. Raw order-13 moments cover every 3D element without repeatedly constructing the expensive high-order simplex fitting matrix. Tests also exercise all direction numbers, enriched nodal fields, alternate GK rules, geometry invalidation, signed normals, and tiny circles.

The supplied Gauss tables use finite decimal precision. The tests therefore also fit independent exact moments through the unchanged production fitting operator and measure its residual. Raw integration error and this fitting floor are reported separately. A fitted residual must be judged against that measured floor.

Timing constructs the Gauss/fitting objects outside the measured loop and clears final fitted weights on every iteration. Two distinct prepared element records alternate, while reusable workspace storage remains available. Geometry-dependent metrics are rebuilt when switching fields and reused within a requested normal family. Five interleaved runs report medians. These measure repeated integration and fitting rather than final-weight cache lookup. `ClearCache` in the production class also clears its geometry preparation, and the new class also refreshes its metric when changing fields. Runtime ratios are machine- and field-dependent; the benchmark prints both absolute times.

The unit-coordinate adapter analyzes local fibers directly and evaluates only the requested surface density. Topology, discriminants, endpoint slopes, and moment coordinates use the same current-slice polynomial. Ancestor-aware algebra remains available separately in the data class. No inner quadrature sample constructs a LINE object.
