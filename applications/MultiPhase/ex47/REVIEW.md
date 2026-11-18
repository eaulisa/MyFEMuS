# Review of the prepared level-set and weight classes

Baseline: `4ec74de89` on `barbi-ex47-levelset-data`. The deferred-Bernstein experiment is excluded. This review covers `CutFemLevelSetData`, its analysis/chart implementation, the local `CutFemWeight_LS`, prepared contexts/moments, and the recursive endpoint integration adapters. Production `src/` is unchanged.

## Findings and changes

| Finding | Consequence | Resolution |
|---|---|---|
| Imported power polynomials were re-interpolated on simplex FEM nodes during vertex permutation | Monomial and Bernstein representations could describe different fields | Compose imported coefficients algebraically; reject charts whose coordinate degree exceeds the supported layout |
| Native chart changes discarded the nodal source and retained only rounded powers | A small cut could change from two simple roots into one double root after an axis swap | Preserve a chart-local nodal source, including compensated low parts for charts of slices; detach reusable storage when an older copy still owns it |
| External-moment validation wrote earlier entries before checking later entries | An exception could leave a partially overwritten moment pair | Validate the entire input before mutating the fitting state |
| External fitting retained the previous integrated-moment validity flag | `GetLastIntegratedMoments()` depended on call history and could report supplied moments as integrated | Invalidate that flag after successful external input validation; preserve it when validation fails |
| WEDGE surface bookkeeping scaled the third component a second time | The auxiliary independence field disagreed with the metric actually used for density evaluation | Remove the redundant factor; retain the necessary reference-metric height conversion |
| `_nodalSurfaceGradient` had no readers or assignments | Dead cache storage and redundant clears remained from the duplicated class | Remove the member and its clears |
| Local shifted-discriminant and compensated-root validation occurred late or incompletely | Invalid requests could modify outputs before throwing | Check endpoint finiteness, planar dimension, and ordered finite root bounds before calculation |

The polynomial inconsistency has a direct reproduction: import `G(x,y)=x²y²` on TRI and swap the two non-origin vertices. At `(0.2,0.3)`, the original implementation returned `0.0025` from monomials while Bernstein evaluation and the original mapped field returned `0.0036`. Both representations now return `0.0036`.

A second independent reproduction uses QUAD nodal values equal to one except for the center value `-1e-24`. Direct restriction through the center finds two simple roots. Before the fix, swapping axes first produced one double root; it now preserves both simple roots. The tests cover tensor and simplex/wedge charts, charts of slices, parent/output reuse, and center values down to `-1e-200` where representable tensor nodal information can be retained.

The chart fix distinguishes two contracts. Native nodal fields retain their known FEM-space closure, including supported restrictions and vertex permutations; that existing path remains unchanged. Imported powers are the polynomial supplied, not permission to project it onto a different nodal space. They are transformed coefficient by coefficient. If a barycentric substitution introduces a genuine coordinate power above two, the data layout cannot represent it and the operation throws with an empty output object. No nonzero higher power is silently discarded. Higher-degree external components retain the separate `TransformParentComponents()` interface.

The WEDGE cleanup does not change the surface measure: the extra factor was applied only to the auxiliary field examined for coordinate dependence. Multiplying a component by a nonzero constant did not change that dependence test. The selected reference/physical metric and evaluated density now also agree in their stored coefficients.

## Representation and coordinate audit

| Operation | Representation and contract |
|---|---|
| Nodal reconstruction | Actual LINE3, QUAD9, TRI7, HEX27, TET15, and WEDGE21 bases and FEMuS ordering |
| Polynomial coordinates | Unit tensor domains, standard simplices, unit triangle times unit height |
| Initial cut selection | Nodal crossing/contact first; Bernstein sign certification otherwise |
| Uncertain sign | Retain as a possible cut; maximum subdivision depth defaults to one |
| Bounds and graph tests | Padded Bernstein bounds; analytic QUAD derivative extrema and weak-graph checks where applicable |
| Values, derivatives, restrictions | Retained unit power coefficients |
| Restricted context | Immediate-parent map, original-root map, original gradient components, and reference embedding Jacobian are separate quantities |
| Boundary certificates | Exact zero traces of the appropriate nodal/control representation, with retained-source checks for cancellation-sensitive slices |
| Geometry | Required only by physical/inverse-metric/normal evaluation, not by level-set classification |
| Fixed quadrature | Shared GK/full-strip tables; moment-to-weight construction occurs in the weight constructor |
| Final fitting | Existing unit-monomial moments and the original FEMuS Gauss rule; reference measure factors are applied at the public boundary |

Native-space closure is structural. TRI7 is quadratic plus its cubic barycentric bubble. TET15 contains quadratic, face-bubble, and volume-bubble terms; fixing a coordinate leaves a TRI7 field. WEDGE21 is TRI7 times a quadratic height polynomial; its coordinate slices are TRI7 or QUAD9. This justifies retaining the native interpolation path. Arbitrary imported coordinate-degree-two powers need not have this closure, which is why they now use the explicit coefficient path.

Tensor reference-surface derivatives require the uniform native/unit factor; WEDGE requires a height-only factor. Physical cofactors are constructed from the supplied geometry in unit coordinates. The public physical, inverse-metric, and signed-normal factors are distinct and are retained. Their scaling, level-set sign reversal, boundary half-weights, and degree-13 raw moments are covered by the existing checks.

The new `HasNodalPolynomialSource()` query distinguishes a nodal polynomial in the active chart from nodal data belonging only to an ancestor. Local projection methods can safely use the former for compensated evaluation; they must not mix the latter into a rounded slice projection. Canonical simplex thirds are evaluated with a compensated remainder when forming nodal snapshots of slices. This is reference-basis arithmetic, not a test-coordinate adjustment. Optional nodal tails are stored only when needed; source storage is reused when unshared, and warmed allocation checks cover the chart operations.

## Numerical special cases retained deliberately

| Mechanism | General reason for retaining it |
|---|---|
| Compensated source contraction | Recover cancellation in evaluating the polynomial determined by the supplied double nodal values |
| Factorized tensor Lagrange evaluation | Preserve a very small nodal value that can disappear in shifted power coefficients |
| Midpoint discriminant identity | For any quadratic on `[0,1]`, `D=(G(1)-G(0))²-4aG(1/2)`; this avoids subtracting nearly equal shifted terms |
| Stable quadratic roots and retained slopes | Use the stable quadratic formula and the discriminant slope instead of recomputing a small derivative at an already rounded root |
| Root-aware gradient reconstruction | Enforce the known directional derivative through the retained chart Jacobian before forming surface norms |
| Separate local and ancestor-aware fiber methods | Projection topology, local roots, and endpoint slopes must refer to the same polynomial; original-nodal refinement is a separate algebraic operation |
| Compensated near-degenerate projection events | Refine quadratic-dominated cancellation-sensitive event polynomials without changing tolerances; this optional path is not a general quartic certificate |
| Legacy weak-graph qualification | Preserve the original generic derivative-range/no-zero-fiber policy when the stronger prepared certificate is inconclusive |
| Endpoint maps and selective retries | Smooth regular folds and refine only failed moments using the existing quadrature policy |
| Warning/exception paths | Return uncertified results with their warnings where supported, or reject unsupported topology; do not force a passing result |

The constants selecting compensated evaluation are conditioning/performance heuristics. They select a more accurate evaluation of the same algebra; they are not mesh IDs, sphere radii, fitted target moments, or replacements for requested quadrature tolerances. No runtime branch matching a test example's coordinates, element identifier, resolution, or expected answer was found.

Bernstein subdivision used for classification is not reused as an integration partition. The integration partition comes from projection events and endpoint maps. Parent coefficients are retained, but these are different kinds of subdivision information.

The rounded-root versus retained-slope tests explain a necessary distinction: a root location can be accurate in absolute terms while its recomputed small derivative is inaccurate relatively. Removing the retained-slope path because an example passes at ordinary resolution would damage small-cut accuracy.

## Remaining precision finding

The extra probe with QUAD nodal values equal to one except for a center value `-1e-24` also exposes a limit in full 2D integration. Cut detection retains it and direct fiber analysis resolves its roots, but the projection event path rejects it with `A tangency strip requires a regular simple fold`. The discriminant is quartic; its rounded power representation has lost the tiny interval, and the current optional compensated-event gate is limited to quadratic-dominated cases. This is a numerically unresolved regular cut, not a mathematical singularity of the original nodal field.

This review does not force that case to pass by changing tolerances, clipping roots, or replacing the field with a circle. Generalizing compensated projection construction and endpoint certification is separate work. The current explicit rejection must remain visible; the successful standard regression matrix is not a claim that all sub-ULP or cancellation-sensitive quartic cuts are supported.

## Cache and state audit

Nodal polynomial and geometry caches compare their actual input values. Coefficients remain valid when only mesh labels change. Restrictions invalidate child analyses; shared analysis/projection records use copy-on-write before replacing cached polynomial-dependent results. Retained source records outlive parent reset. Geometry cofactors are prepared once for matching nodal/geometry inputs and restricted as polynomial components.

Final weight cache keys include field values, geometry values, tolerances, and GK configuration. Direction/validation setters invalidate final weights. The positive moment pair is reused across 3D surface modes for matching inputs. Clearing final weights does not rebuild quadrature rules or discard valid polynomial preparation.

`GetLastIntegratedMoments()` means the last pair actually integrated. A cache hit does not replace it. External fitting now invalidates it because that operation overwrites the fitting pair with supplied values. Invalid external inputs are rejected before state changes. These distinctions are explicitly tested.

## Limits that the review does not conceal

- Cut classification is conservative; an uncertified element is not proof of an intersection. The requested depth-one default can classify more elements as candidates than the original depth-twelve default.
- Root records can retain slopes after a root location has rounded, but this does not make arbitrary sub-ULP multidimensional projection intervals or composite metric evaluations exactly resolvable.
- Adaptive error estimates and weak-graph qualifications are floating-point procedures, not a formal interval-arithmetic proof for every polynomial.
- Repeated/non-transverse projection factors and singular delta measures are not universally supported. The existing nearly-parallel TRI exception remains visible in the extended results.
- Imported polynomials can require a larger coordinate degree after chart changes. The current fixed layout rejects such transformations; the nodal FEM spaces are the supported closed family.
- Output conversion warnings do not constitute an exhaustive certificate of every fitted Gauss moment. High-order fitting can have a precision floor even when raw moments are much more accurate. The LINE degree-ten floor is independently reproduced by fitting exact reference moments through the original fitter; it is not hidden by modifying root integration.
- The near-fold `4e-16` assertion belongs to a multiprecision comparison in the test main. It is not a runtime integration tolerance or a uniform relative-error guarantee arbitrarily close to a multiple root. No assertion or integration tolerance was relaxed during this review.
- The inherited endpoint-strip root-admission guard includes a fixed `1e-12` slack. It is a topology guard, not a requested-accuracy certificate; it was not introduced or enlarged here. Replacing that policy would need a separate error-analysis change rather than an example-specific threshold adjustment.

## Validation

Fresh main, weight, sanitizer, complete native-source, public-case, high-degree moment, and timing results are recorded in `native_tests/results/deep_review/`. All newly added tests remain in the ex47 main. They include imported chart identities, unsupported-degree rejection, repeated native charts and slices, external-fitting transaction state, local-discriminant output preservation, and WEDGE metric consistency. Production files and original native fixture assertions are unchanged.

| Final check | Result |
|---|---|
| Data main, all six elements | 5,903,365 checks; maximum scaled error `2.365e-14` |
| Weight main, all six elements | 32,525 checks; degree-13 raw 3D errors at most `2.031e-16` |
| Public native regression matrix | 386/386 for each class; 36,860 assertions per class; maximum class difference `1.333e-15` |
| Complete native-source matrix | 50 fixture sources, 385 requests; prepared 383 pass, two documented exceptions |
| Matched-depth CTest control | 32/32 for each class, detection depth 12 |
| Graph/surface moment comparison | 59,520 checks pass; largest scaled class difference `1.771e-12` |
| Additional fixture dump pairing | 26,260 checks against unchanged production binaries pass |
| ASAN/UBSAN | 157,427 data and 3,657 LINE weight checks pass; leak detection disabled under ptrace |
| Warmed reuse/allocation checks | Zero instrumented allocations in the recorded algebra and data operations |
| Interleaved timing | Three batches, 386 cases and 772 method rows; all assertions pass for both classes |

The complete matrix retains two exceptions without changing their assertions: `nodal_cut_detection` expects the original depth-twelve certification whereas the requested prepared default is one; `wide_unit_validation_1_2` reaches an unsupported nearly parallel TRI projection. The matched-depth control passes. Detailed failures, per-case results, paired comparisons, source hashes, and timing ranges are in the accompanying results directory.

Timing compares the final prepared implementation with production. Median production/prepared time ratios are 2.83 for LINE, 1.05 for QUAD, 1.59 for TRI, 2.14 for HEX, and 2.53 for TET. These are medians of per-method ratios, not total application speedups. Individual slower cases remain, particularly QUAD and inexpensive LINE calls where preparation/state overhead matters. The CSV preserves every row and sample range. This measurement is not a controlled before/after attribution for each review fix; it does not establish that every operation is faster. WEDGE is covered by the main accuracy/allocation checks rather than this five-element public matrix.
