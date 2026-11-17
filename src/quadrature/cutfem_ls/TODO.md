# CutFemWeight_LS tasks

- [ ] Reverse explicit cut-detection checks: test nodal level-set signs first, then compute Bernstein coefficients only when all nodal values have the same strict sign.
- [ ] Expose the maximum Bernstein subdivision depth as a configurable parameter instead of the hard-coded depth of 12; allow depth 0 for a single Bernstein test without subdivision.
- [ ] Set a substantially lower default subdivision depth; the exact default remains to be chosen.
- [ ] Retain and share the parent power coefficients and parent Bernstein coefficients between cut detection and integration, including the different surface-weight computations; rebuild them only when the nodal level set or reference-element context changes, independently of moment order, integration tolerances, and Gauss-Kronrod settings.
- [ ] Introduce a lightweight per-element `CutFemLevelSetData` class whose constructor receives the element type, nodal level-set values, and detection settings; own the nodal snapshot and classify automatically, constructing coefficients only when first needed.
- [ ] Keep parent power/Bernstein coefficients, classification reason, explicit reference-coordinate convention, normalization scale, and lazily computed parent derivatives in the element data; distinguish certified phases, nodal crossings, contacts, and unresolved cases.
- [ ] Retain the same element-data record from classification through assembly and reuse it for volume, scalar surface, and all three normal-weight computations; associate records with the mesh level, element, and nodal-data version, and rebuild or transfer them after refinement or redistribution.
- [ ] Keep geometry data separate and invalidate geometry-dependent derived quantities when either the geometry or level-set state changes; obtain slice, subregion, and chart coefficients from retained parent data.
- [ ] Keep fixed Gauss rules and moment-to-weight matrices in the reusable quadrature objects; keep final weight arrays in integration workspace by default rather than permanently storing every weight family in every element record.
- [ ] Add integration interfaces that consume prepared element data and migrate detection and assembly to use them; preserve the existing nodal interfaces as compatibility wrappers.

Preserve the roundoff margins and classify opposite nodal signs, zero nodal values, and unresolved Bernstein cases as possible cuts.
