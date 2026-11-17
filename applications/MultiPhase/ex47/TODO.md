# CutFemLevelSetData: local ex47 status

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

Pending production integration after acceptance of this prototype:

- [ ] Move the accepted class into the shared library and retain its records across actual mesh-level classification and assembly loops.
- [ ] Add prepared-data overloads to `CutFemWeight_LS` while preserving its nodal wrappers; reuse records in every volume/surface weight mode.
- [ ] Connect record lifecycle and geometry-dependent cache invalidation to actual level-set updates, refinement, geometry changes, and MPI redistribution.
- [ ] Remove unnecessary intermediate double conversions from the existing TRI7/WEDGE21 reconstruction and TRI/TET moment-output paths, with integration regressions.
- [ ] Validate the resulting complete integration paths, high moment orders, and production memory/runtime behavior before selecting lower precision in well-conditioned kernels.

Incremental updates between neighboring slice coordinates remain deferred; restrictions are evaluated directly at GK points. Existing production quadrature and its one-time multiprecision matrix construction are unchanged.
