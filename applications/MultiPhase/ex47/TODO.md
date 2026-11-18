# CutFemLevelSetData: local ex47 status

Historical record from before promotion into the native library; see [MIGRATION.md](MIGRATION.md) for the current layout and validation.

Implemented and exercised only in ex47:

- [x] Check nodal signs before constructing Bernstein coefficients; retain conservative contact and uncertified classifications.
- [x] Make maximum certification depth configurable, with default 1 instead of 12 and depth 0 for a single Bernstein test.
- [x] Retain parent power/Bernstein coefficients and normalization, with lazy preparation and owned nodal snapshots.
- [x] Support the actual LINE3, QUAD9, TRI7, HEX27, TET15, and WEDGE21 bases in long-double coefficient algebra.
- [x] Provide classification reasons, bounds, lazy derivative/Hessian caches, graph certification, context/version metadata, and explicit reference-coordinate conventions.
- [x] Reuse prepared element records and prove ownership/invalidation behavior in main, including copies and parent resets.
- [x] Add recursive `Restrict(axis, coordinate, child)` using the same data class for 3D, 2D, and 1D objects.
- [x] Restrict power coefficients directly with Horner evaluation and Bernstein coefficients directly from the parent control net.
- [x] Preserve simplex/wedge rescaling, normalization, degenerate slices, root-coordinate maps, and every original derivative component.
- [x] Carry parent roundoff margins through slice classification and derivative bounds, with cancellation regressions.
- [x] Reuse child storage and specialize identity/zero-origin chart transformations.
- [x] Keep geometry separate and support exact restriction of external polynomial components without reducing their degree or component count.
- [x] Preserve long-double arithmetic through coefficient conversion and restriction; compare double and extended arithmetic against multiprecision references in main.
- [x] Keep all executable tests inside ex47 main, including nested volume/surface and nonuniform GK slice checks.
- [x] Keep fixed quadrature objects, moment matrices, and final weight arrays outside the new class.
- [x] Share immutable reference indices, simplex multipliers, and exact Bernstein restriction plans; use strided Horner kernels and fixed-size subdivision buffers.
- [x] Cache bounds and graph certificates, with derivative bounds obtained directly from Bernstein controls.
- [x] Add allocation-free gradient, scalar-line, and root interfaces while preserving vector compatibility wrappers.
- [x] Evaluate all existing regular-root surface densities together from prepared cofactor-gradient values or a prepared cofactor matrix.
- [x] Test geometry/cofactor preparation and recursive restriction for every element and all surface modes, including bi/triquadratic QUAD/HEX geometry.
- [x] Benchmark the actual quadrature kernels against formatting commit `26cbb141` and assert zero steady-state heap allocations in the fixed-output paths.
- [x] Select one surface measure with a dedicated evaluator; reserve all-measure evaluation for explicitly named batch methods.
- [x] Share only the derivative and metric required by a requested normal family, with a separate single-component interface.
- [x] Omit original-gradient preparation/contraction for delta and prepared-cofactor-gradient modes; propagate this choice through recursive restriction.
- [x] Check selected results against batch/analytic/multiprecision references, including small cuts and input scales, and benchmark selection outside point loops.
- [x] Retain complete fixed-size LINE topology, stable signed root slopes, multiplicities, endpoint tags, positive intervals, and identically-zero status.
- [x] Add prepared-root overloads for every selected and batch surface measure, retaining compensated parent restriction near tangencies and cancellation.
- [x] Detect exact zero nodal traces and restrict every boundary, including sloping simplex faces, with original gradients and reference embedding maps.
- [x] Add supported axis/barycentric charts, elevated swapped WEDGE slice layouts, native-to-unit coefficient exports, and full-degree external component transforms.
- [x] Add strict/weak graph information, tighter QUAD derivative extrema, conservative zero-fiber status, and exact coefficient dependence/reflection queries.
- [x] Retain normalized projection polynomials lazily and expose compensated discriminants, planar boundary polynomials, and shifted endpoint expansions.
- [x] Keep LINE/graph/projection caches optional and preserve copied records when reusable workspaces detach from shared prepared state.
- [x] Exercise every chart/boundary and near-tangent/tiny-cut case in main; benchmark fresh and retained paths against the preceding class.

Local weight integration in this prototype:

- [x] Duplicate the weight class as `Ex47::CutFemWeight_LS` under ex47/include and preserve the production class for side-by-side tests.
- [x] Route nodal wrappers through owned prepared data and provide overloads accepting retained root-element data.
- [x] Use prepared classification, graphs, zero boundaries, projections, LINE topology, original gradients, and component restrictions throughout the recursive integration.
- [x] Preserve endpoint parameter derivatives and offset information; avoid rebuilding a LINE for a delta density when its slope is already known.
- [x] Integrate in native coordinates and transform final moments once into the unchanged fitting basis.
- [x] Retain physical cofactor-gradient polynomials across surface-mode, tolerance, rule, and final-weight cache changes; invalidate on changed field or geometry.
- [x] Keep Gauss points, GK objects, and the fitting matrix outside the element and quadrature-point loops.
- [x] Exercise LINE first, then QUAD/TRI, then HEX/TET/WEDGE in main, with independent moments and original/prepared timing comparisons.

Pending production integration after acceptance of this prototype:

- [ ] Move the accepted class into the shared library and retain its records across actual mesh-level classification and assembly loops.
- [ ] Add prepared-data overloads to `CutFemWeight_LS` while preserving its nodal wrappers; reuse records in every volume/surface weight mode.
- [ ] Reuse known endpoint-graph roots and signed chart derivatives in prepared root evaluators, preserving offset/parameter information when an absolute sample coordinate rounds to an endpoint.
- [ ] Integrate all requested surface moment families in one traversal using the shared densities, roots, and prepared geometry.
- [ ] Connect record lifecycle and geometry-dependent cache invalidation to actual level-set updates, refinement, geometry changes, and MPI redistribution.
- [ ] Remove unnecessary intermediate double conversions from the existing TRI7/WEDGE21 reconstruction and TRI/TET moment-output paths, with integration regressions.
- [ ] Validate the resulting complete integration paths, high moment orders, and production memory/runtime behavior before selecting lower precision in well-conditioned kernels.

Incremental updates between neighboring slice coordinates remain deferred; restrictions are evaluated directly at GK points. Existing production quadrature and its one-time multiprecision matrix construction are unchanged.

Native-suite comparison follow-up (local ex47 only):

- [ ] Preserve endpoint offsets/parameter derivatives through native-coordinate HEX slicing; `native_tests/curved_public_regression.cpp` and the curved HEX27 native fixture currently fail while the original unit-coordinate path passes.
- [ ] Resolve the enriched HEX bubble surface-moment error and the degree-10 cubic-contact HEX rejection before production integration.
- [ ] Restore forced physical-surface HEX projection order 1 in the outer-interval graph benchmark; automatic directions and the other tested directions succeed.
- [ ] Skip LINE restriction/root analysis for a certified boundary root; the low-order QUAD boundary benchmarks currently use the same 15 samples but take about twice as long.
- [ ] Preserve and check analytic LINE moment error estimates at the public weight interface so unattainable tolerances retain the original aggregated warning behavior.
- [ ] Keep depth-1 conservative detection distinct from depth-12 strict-certification comparisons; the shallower default intentionally retains more possible cuts.
