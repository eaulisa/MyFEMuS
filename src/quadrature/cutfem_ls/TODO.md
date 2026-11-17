# CutFemWeight_LS tasks

Checked items are implemented and tested only in the local [ex47 prototype](../../../applications/MultiPhase/ex47/README.md). Production assembly and integration remain unchanged; partially completed integration tasks stay open. See [local status](../../../applications/MultiPhase/ex47/TODO.md) and [validation](../../../applications/MultiPhase/ex47/VALIDATION.md).

- [x] Reverse explicit cut-detection checks: test nodal level-set signs first, then compute Bernstein coefficients only when all nodal values have the same strict sign.
- [x] Expose the maximum Bernstein subdivision depth as a configurable parameter instead of the hard-coded depth of 12; allow depth 0 for a single Bernstein test without subdivision.
- [x] Set a substantially lower default subdivision depth; the ex47 prototype uses depth 1.
- [ ] Retain and share the parent power coefficients and parent Bernstein coefficients between cut detection and integration, including the different surface-weight computations; rebuild them only when the nodal level set or reference-element context changes, independently of moment order, integration tolerances, and Gauss-Kronrod settings.
- [x] Introduce a lightweight per-element `CutFemLevelSetData` class whose constructor receives the element type, nodal level-set values, and detection settings; own the nodal snapshot and classify automatically, constructing coefficients only when first needed.
- [x] Keep parent power/Bernstein coefficients, classification reason, explicit reference-coordinate convention, normalization scale, and lazily computed parent derivatives in the element data; distinguish certified phases, nodal crossings, contacts, and unresolved cases.
- [ ] Retain the same element-data record from classification through assembly and reuse it for volume, scalar surface, and all three normal-weight computations; associate records with the mesh level, element, and nodal-data version, and rebuild or transfer them after refinement or redistribution.
- [ ] Keep geometry data separate and invalidate geometry-dependent derived quantities when either the geometry or level-set state changes; obtain slice, subregion, and chart coefficients from retained parent data.
- [x] Keep fixed Gauss rules and moment-to-weight matrices in the reusable quadrature objects; keep final weight arrays in integration workspace by default rather than permanently storing every weight family in every element record.
- [ ] Add integration interfaces that consume prepared element data and migrate detection and assembly to use them; preserve the existing nodal interfaces as compatibility wrappers.
- [x] Add recursive `Restrict(axis, coordinate, child)` to `CutFemLevelSetData`, producing a lower-dimensional object of the same class: HEX to QUAD, TET to TRI, WEDGE to TRI for height slices or QUAD for triangle-coordinate slices, and QUAD/TRI to LINE.
- [x] Initialize slice data directly from parent coefficients using Horner contraction for power coefficients and Bernstein contraction/de Casteljau evaluation when Bernstein coefficients are needed; avoid rebuilding nodal values or invoking the nodal constructor for slices.
- [x] Carry the slice domain, active coordinates, reference-domain rescaling, normalization, classification, and derivative context through each restriction; keep integration dimension separate from physical dimension and retain all three physical surface-factor components for slices of a 3D element.
- [x] Allocate one reusable child workspace per integration level; overwrite slice-dependent data and invalidate derived caches when the slicing coordinate changes while retaining allocated storage and the shared fixed quadrature objects.
- [x] Specialize identity and zero-origin scaling chart transformations to avoid unnecessary generic polynomial transformations during slicing.
- [x] Keep parent power/Bernstein coefficients, derivatives, normalization, and all recursive restrictions in `long double`, with explicit precision boundaries for double nodal inputs and final outputs.
- [ ] Audit and remove unnecessary intermediate conversions to `double`, particularly TRI7 coefficients, their use in WEDGE21, and TRI/TET moment outputs before moment-to-weight accumulation; preserve existing public input/output compatibility.
- [ ] Preserve extended precision in sensitive root/discriminant calculations, geometry/cofactor polynomial algebra, moment accumulation, and compensated weight construction; retain the one-time multiprecision construction of the fixed moment-to-weight matrices.
- [ ] Compare accuracy and runtime of double and extended-precision paths on all supported elements, including small cuts, near tangencies, cancellation, and high moment orders, before considering lower precision in frequently called well-conditioned kernels.

Preserve the roundoff margins and classify opposite nodal signs, zero nodal values, and unresolved Bernstein cases as possible cuts.

Incremental updates between neighboring slice coordinates are deferred; the current plan evaluates restrictions directly at the nonuniform Gauss-Kronrod points.

The ex47 speed revision starts from formatting commit `26cbb141`, shares immutable reference tables, uses allocation-free prepared restriction/certification, and supplies every regular-root surface density from prepared cofactor-gradient data or a constant affine cofactor matrix. See [kernel timings](../../../applications/MultiPhase/ex47/PERFORMANCE.md); migration into the production integration traversal remains pending.
